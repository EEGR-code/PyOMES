# -*- coding: utf-8 -*-
"""Thermodynamic parameter configuration for consistent acid-base chemistry.

Provides :class:`ThermodynamicConfig` — a single source of truth for
pKa reference values and van 't Hoff enthalpy corrections (dH) used
across the speciation engine, gas-liquid transfer corrections, and
context functions.

Usage
-----
>>> from PyOMES.chemistry.thermo_params import ThermodynamicConfig
>>> tc = ThermodynamicConfig.bsm2_default()
>>> tc.apply_to_cv(cv)
>>> ctx = tc.context_kwargs_at_current_T()

Each acid species is defined with:
- Reference pKa value(s) at a specified reference temperature
- A temperature correction method (``"none"`` or ``"van_t_hoff"``)
- Enthalpy values (dH) when using van 't Hoff correction

Different species may have different reference temperatures, allowing
thermodynamic data from multiple literature sources to be combined.
"""

from __future__ import annotations

import math
import warnings
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple

from ..units import R_J_PER_MOL_K as _R_J
_VALID_CORRECTIONS = ("none", "van_t_hoff")

# ════════════════════════════════════════════════════════════════════════
#  Unit conversion helpers
# ════════════════════════════════════════════════════════════════════════
# All internal storage uses: pKa (dimensionless), J/mol, K.
# These helpers convert user-supplied values to internal units at the
# point of entry so no unit tracking is needed downstream.

_DH_FACTORS = {
    "J/mol":    1.0,
    "j/mol":    1.0,
    "kJ/mol":   1000.0,
    "kj/mol":   1000.0,
    "cal/mol":  4.184,
    "Cal/mol":  4.184,
    "kcal/mol": 4184.0,
}

_VALID_K_FORMATS = ("pKa", "pka", "Ka", "ka", "lnKa", "lnka")
_VALID_DH_UNITS = tuple(_DH_FACTORS.keys())
_VALID_T_UNITS = ("K", "k", "C", "c", "°C")


def _convert_dH_to_J(value: float, unit: str) -> float:
    """Convert an enthalpy value to J/mol."""
    factor = _DH_FACTORS.get(unit)
    if factor is None:
        # Try case-insensitive match
        for key, f in _DH_FACTORS.items():
            if key.lower() == unit.lower():
                factor = f
                break
    if factor is None:
        raise ValueError(
            f"Unknown dH unit {unit!r}. Options: {list(_DH_FACTORS.keys())}")
    return float(value) * factor


def _convert_dH_tuple(values, unit: str) -> Tuple[float, ...]:
    """Convert a tuple of enthalpy values to J/mol."""
    if values is None:
        return None
    if not isinstance(values, (tuple, list)):
        values = (values,)
    return tuple(_convert_dH_to_J(v, unit) for v in values)


def _convert_T_to_K(value: float, unit: str) -> float:
    """Convert temperature to Kelvin."""
    u = unit.strip()
    if u in ("K", "k"):
        return float(value)
    if u in ("C", "c", "°C"):
        return float(value) + 273.15
    raise ValueError(
        f"Unknown temperature unit {unit!r}. Options: 'K', 'C', '°C'")


def _Ka_to_pKa(Ka: float) -> float:
    """Convert Ka (mol/L) to pKa."""
    if Ka <= 0:
        raise ValueError(f"Ka must be positive, got {Ka}")
    return -math.log10(Ka)


def _lnKa_to_pKa(lnKa: float) -> float:
    """Convert ln(Ka) to pKa."""
    return -lnKa / math.log(10)


@dataclass(frozen=True)
class AcidDefinition:
    """One acid species with its reference pKa and correction method.

    Parameters
    ----------
    pKas : tuple of float
        Reference pKa value(s) at ``T_ref_K``.
    T_ref_K : float
        Reference temperature (K) for these pKa values.
    correction : str
        Temperature correction method: ``"none"`` or ``"van_t_hoff"``.
    dH_J_per_mol : tuple of float
        Enthalpy of dissociation (J/mol) per pKa step.
    """
    pKas: Tuple[float, ...]
    T_ref_K: float = 298.15
    correction: str = "none"
    dH_J_per_mol: Tuple[float, ...] = ()

    def __post_init__(self):
        if not isinstance(self.pKas, tuple):
            object.__setattr__(self, "pKas", tuple(self.pKas))
        if not isinstance(self.dH_J_per_mol, tuple):
            object.__setattr__(self, "dH_J_per_mol", tuple(self.dH_J_per_mol))
        if len(self.dH_J_per_mol) == 0:
            object.__setattr__(self, "dH_J_per_mol", (0.0,) * len(self.pKas))

    def pKas_at_T(self, T_K: float) -> Tuple[float, ...]:
        """Compute pKa(s) at temperature T_K."""
        if self.correction == "none":
            return self.pKas
        if abs(T_K - self.T_ref_K) < 0.01:
            return self.pKas
        result = []
        for pKa_ref, dH in zip(self.pKas, self.dH_J_per_mol):
            if abs(dH) < 1e-10:
                result.append(pKa_ref)
            else:
                Ka_ref = 10.0 ** (-pKa_ref)
                Ka_T = Ka_ref * math.exp(
                    -dH / _R_J * (1.0 / T_K - 1.0 / self.T_ref_K))
                result.append(-math.log10(max(Ka_T, 1e-30)))
        return tuple(result)


@dataclass(frozen=True)
class WaterDefinition:
    """Water autoionisation constant definition."""
    pKw: float = 14.0
    T_ref_K: float = 298.15
    correction: str = "none"
    dH_J_per_mol: float = 0.0

    def Kw_at_T(self, T_K: float) -> float:
        Kw_ref = 10.0 ** (-self.pKw)
        if self.correction == "none" or abs(self.dH_J_per_mol) < 1e-10:
            return Kw_ref
        if abs(T_K - self.T_ref_K) < 0.01:
            return Kw_ref
        return Kw_ref * math.exp(
            -self.dH_J_per_mol / _R_J * (1.0 / T_K - 1.0 / self.T_ref_K))

    def pKw_at_T(self, T_K: float) -> float:
        return -math.log10(max(self.Kw_at_T(T_K), 1e-30))


@dataclass
class ThermodynamicConfig:
    """Single source of truth for acid-base thermodynamic parameters.

    .. deprecated::
        chemistry-unification-3b (2026-06-02).  Activity-model
        configuration migrates to :class:`~PyOMES.thermo.ThermoFramework`;
        pKa data migrates to :class:`~PyOMES.chemistry.ChemistryDatabase`
        stock modules.  ``ThermodynamicConfig`` remains operational for
        legacy code paths; it will be removed when ``ChemistryDatabase``
        covers all callers.  New code should use ``ThermoFramework`` and
        ``ChemistryDatabase`` instead.

    Internally holds an :class:`~fermenter.chemistry.equilibria.EquilibriumSet`
    that stores all acid/water definitions.  The ``add_acid`` and
    ``set_water`` methods delegate to the set, adding ``category`` and
    ``n_active`` parameters for the charge balance.

    Build using :meth:`add_acid` and :meth:`set_water`, or use a
    factory method like :meth:`bsm2_default`.

    Parameters
    ----------
    T_ref_K : float
        Default reference temperature (K).  Used as fallback when
        :meth:`add_acid` is called without an explicit ``T_ref_K``.
    """
    T_ref_K: float = 298.15
    _equilibria: Any = field(default=None, repr=False)
    _cv_ref: Any = field(default=None, repr=False)

    def __post_init__(self):
        if self._equilibria is None:
            from .equilibria import EquilibriumSet
            self._equilibria = EquilibriumSet(T_ref_K=self.T_ref_K)

    @property
    def equilibria(self):
        """The underlying :class:`EquilibriumSet`."""
        return self._equilibria

    # ── Builder methods (delegate to EquilibriumSet) ──────────────────

    def add_acid(self, species, *, pKas=None, Ka=None, lnKa=None,
                 correction="none", dH=None, dH_unit="J/mol",
                 dH_J_per_mol=None, T_ref=None, T_unit="K",
                 T_ref_K=None, category="acid", n_active=None,
                 total_key=None):
        """Add or replace an acid species definition.

        Delegates to :meth:`EquilibriumSet.add` with full unit
        conversion support.  See that method for parameter details.

        Parameters
        ----------
        species : str
            Species identifier (e.g. ``"CO2"``, ``"S_ac"``).
        pKas, Ka, lnKa : tuple of float
            Dissociation constant(s) — provide exactly one format.
        correction : str
            ``"none"`` or ``"van_t_hoff"``.
        dH, dH_unit : tuple + str
            Enthalpy with unit (``"J/mol"``, ``"kJ/mol"``, etc.).
        dH_J_per_mol : tuple of float
            Enthalpy in J/mol (takes precedence over dH + dH_unit).
        T_ref, T_unit : float + str
            Reference temperature with unit (``"K"`` or ``"C"``).
        T_ref_K : float
            Reference temperature in K (takes precedence).
        category : str
            Charge balance role: ``"acid"`` (default), ``"cation_acid"``,
            ``"inorganic_acid"``, or ``"strong_ion"``.
        n_active : int, optional
            Number of dissociation steps active in the charge balance.
            Defaults to ``len(pKas)``.
        total_key : str, optional
            Context dict key for total concentration.  Defaults to
            ``"CT_{species}"``.

        Returns
        -------
        self
        """
        self._equilibria.add(
            species, category=category,
            pKas=pKas, Ka=Ka, lnKa=lnKa,
            n_active=n_active, correction=correction,
            dH=dH, dH_unit=dH_unit, dH_J_per_mol=dH_J_per_mol,
            T_ref=T_ref, T_unit=T_unit, T_ref_K=T_ref_K,
            total_key=total_key)
        return self

    def set_water(self, *, pKw=None, Kw=None, correction="none",
                  dH=None, dH_unit="J/mol", dH_J_per_mol=None,
                  T_ref=None, T_unit="K", T_ref_K=None):
        """Set the water autoionisation definition.

        Delegates to :meth:`EquilibriumSet.set_water` with full unit
        conversion support.

        Returns self for method chaining.
        """
        self._equilibria.set_water(
            pKw=pKw, Kw=Kw, correction=correction,
            dH=dH, dH_unit=dH_unit, dH_J_per_mol=dH_J_per_mol,
            T_ref=T_ref, T_unit=T_unit, T_ref_K=T_ref_K)
        return self

    # ── Factory methods ───────────────────────────────────────────────

    @staticmethod
    def bsm2_default():
        """BSM2-canonical parameters (Rosen & Jeppsson 2006, Table 3).

        CO₂ is diprotic but only the first dissociation is active in
        the charge balance (``n_active=1``), matching BSM2's monoprotic
        treatment.  VFA pKas have no temperature correction.
        """
        from .equilibria import EquilibriumSet
        tc = ThermodynamicConfig(T_ref_K=298.15)
        tc._equilibria = EquilibriumSet.bsm2_default()
        return tc

    # ── Live temperature from attached zone ───────────────────────────

    @property
    def T_K_live(self):
        """Current operating temperature from the attached zone."""
        if self._cv_ref is not None:
            phases = getattr(self._cv_ref, "phases", None)
            liq = phases.get("liquid") if phases is not None else None
            if liq is not None:
                return float(getattr(liq, "T_K", self.T_ref_K))
        return self.T_ref_K

    # ── Apply to CV ───────────────────────────────────────────────────

    def apply_to_cv(self, cv):
        """Attach this config to a CV as a back-reference.

        Stores ``self`` on ``cv._thermo_config`` and sets ``self._cv_ref
        = cv`` so :attr:`T_K_live` (used by ``chem_env_fn`` to read the
        zone temperature dynamically) can find the live phase state.

        As of ``chemistry-unification-2`` this is the only remaining
        responsibility.  Pre-Phase-2 the call also populated the
        gas-liquid link's ``speciation_corrections`` (pKa/dH/n_active
        snapshots driving inline ``f_molecular(pH, T_K)`` reads).
        Those corrections are gone — the link now reads molecular
        fractions from :attr:`PropertyResult.alphas` populated by the
        speciation engine, whose pKas come from declared equilibrium
        reactions.  The engine-side hook was already a no-op after
        Phase 1.  Phase 3 will absorb this config into
        ``ChemistryDatabase`` and the back-reference can move there.
        """
        self._cv_ref = cv
        cv._thermo_config = self

    # ── Temperature-corrected accessors ───────────────────────────────

    def pKa_at_T(self, species, T_K):
        """Compute pKa(s) at temperature T_K using the species' own T_ref."""
        if not self._equilibria.has(species):
            return ()
        return self._equilibria.get(species).pKas_at_T(T_K)

    def Kw_at_T(self, T_K):
        return self._equilibria.water.Kw_at_T(T_K)

    def pKw_at_T(self, T_K):
        return self._equilibria.water.pKw_at_T(T_K)

    def context_kwargs_at_current_T(self):
        """Speciation kwargs pre-corrected to the zone's current temperature.

        Returns pKa and Kw values already corrected to operating T.
        The dH values are zero so the speciation engine does NOT apply
        a second correction.
        """
        T = self.T_K_live
        co2_pKas = self.pKa_at_T("CO2", T)
        nh4_pKas = self.pKa_at_T("NH4", T)

        acid_pKas = {}
        for eq_def in self._equilibria:
            if eq_def.name not in ("CO2", "NH4") and eq_def.category in ("acid",):
                corrected = eq_def.pKas_at_T(T)
                if len(corrected) == 1:
                    acid_pKas[eq_def.name] = corrected[0]

        result = {
            "acid_pKas": acid_pKas,
            "Kw": self.Kw_at_T(T),
            "dH_Kw_J_per_mol": 0.0,
            "dH_TIC1_J_per_mol": 0.0,
            "dH_TIC2_J_per_mol": 0.0,
            "dH_NH_J_per_mol": 0.0,
        }
        if len(co2_pKas) >= 1:
            result["pKa1_TIC"] = co2_pKas[0]
        if len(co2_pKas) >= 2:
            result["pKa2_TIC"] = co2_pKas[1]
        if len(nh4_pKas) >= 1:
            result["pKa_NH"] = nh4_pKas[0]
        return result

    # ── Legacy-compatible property accessors ──────────────────────────
    # These provide backward compatibility with code that reads the old
    # separate dicts.  They construct views from the EquilibriumSet.

    @property
    def acids(self):
        """Legacy dict view: ``{species: AcidDefinition}``."""
        result = {}
        for eq_def in self._equilibria:
            result[eq_def.name] = AcidDefinition(
                pKas=eq_def.pKas,
                T_ref_K=eq_def.T_ref_K,
                correction=eq_def.correction,
                dH_J_per_mol=eq_def.dH_J_per_mol)
        return result

    @property
    def water(self):
        """Legacy water accessor."""
        w = self._equilibria.water
        return WaterDefinition(
            pKw=w.pKw, T_ref_K=w.T_ref_K,
            correction=w.correction, dH_J_per_mol=w.dH_J_per_mol)

    @property
    def pKa_ref(self):
        """Reference pKa values (backward-compatible dict view)."""
        return {eq_def.name: eq_def.pKas for eq_def in self._equilibria}

    @property
    def dH_J_per_mol(self):
        """dH values (backward-compatible dict view)."""
        return {eq_def.name: eq_def.dH_J_per_mol for eq_def in self._equilibria}

    @property
    def pKw_ref(self):
        return self._equilibria.water.pKw

    @property
    def dH_Kw_J_per_mol(self):
        return self._equilibria.water.dH_J_per_mol

    def context_dH_kwargs(self):
        """Legacy: return dH keys for context dict."""
        co2_dH = self.dH_J_per_mol.get("CO2", (0.0, 0.0))
        nh4_dH = self.dH_J_per_mol.get("NH4", (0.0,))
        return {
            "dH_Kw_J_per_mol": self.dH_Kw_J_per_mol,
            "dH_TIC1_J_per_mol": co2_dH[0] if len(co2_dH) > 0 else 0.0,
            "dH_TIC2_J_per_mol": co2_dH[1] if len(co2_dH) > 1 else 0.0,
            "dH_NH_J_per_mol": nh4_dH[0] if len(nh4_dH) > 0 else 0.0,
        }

    def context_acid_pKas(self):
        """Legacy: return VFA pKa dict at reference T."""
        result = {}
        for eq_def in self._equilibria:
            if eq_def.name not in ("CO2", "NH4") and len(eq_def.pKas) == 1:
                result[eq_def.name] = eq_def.pKas[0]
        return result

    # ── Charge balance helper ─────────────────────────────────────────

    def compute_CT_cation_from_charge_balance(self, pH, acid_totals,
                                               CT_anion=0.0, T_K=None):
        """Compute inert strong cation concentration to close the charge balance.

        Iterates over all acid-base systems registered in this
        :class:`ThermodynamicConfig`, using the stored ``n_active`` and
        temperature-corrected pKa values for each.  This guarantees
        consistency with the runtime speciation solver, which reads
        from the same equilibrium definitions.

        Parameters
        ----------
        pH : float
            Target pH.
        acid_totals : dict
            ``{acid_name: total_concentration_mol_L}`` using the same
            names as registered equilibrium definitions.  Example::

                {"CO2": 0.0951, "NH4": 0.0945,
                 "S_ac": 0.00140, "S_pro": 0.000157, ...}
        CT_anion : float
            Total inert strong anion concentration (mol_charge/L).
            Equivalent to BSM2's ``S_anion``.
        T_K : float, optional
            Temperature (K).  Defaults to the config's reference T.

        Returns
        -------
        float
            Monovalent-equivalent inert cation concentration (mol/L)
            needed to close the charge balance at the target pH.
            Clamped to >= 0 since a physical cation concentration
            cannot be negative.

        Notes
        -----
        ``CT_cation`` and ``CT_anion`` are monovalent-equivalent inert
        charge sinks for use when specific ion identities are unknown
        (e.g. BSM2's ``S_cation`` / ``S_anion``).  They contribute to
        the charge balance and ionic strength (as z=1 ions) but do not
        participate in complexation, precipitation, or other
        ion-specific chemistry.

        When the contributing ions are known (Na⁺, Ca²⁺, Cl⁻, etc.),
        use the specific strong ion keys instead — these carry the
        correct valence for ionic strength and enable ion-specific
        chemistry at higher speciation levels.

        ``CT_cation`` assumes +1 charge per mole.  If the real unknown
        ions are polyvalent, the ionic strength contribution will be
        underestimated (by a factor of z² for an ion of charge z).
        """
        T = T_K if T_K is not None else self.T_K_live
        H = 10.0 ** (-pH)
        Kw = self.Kw_at_T(T)

        cations = H
        anions = Kw / H + CT_anion

        for eq_def in self._equilibria:
            name = eq_def.name
            CT = acid_totals.get(name, 0.0)
            if CT <= 0.0:
                continue

            # Temperature-corrected pKas, limited to n_active
            active_pKas = eq_def.active_pKas_at_T(T)
            active_Kas = [10.0 ** (-pk) for pk in active_pKas]

            if not active_Kas:
                continue

            if eq_def.category == "cation_acid":
                # BH⁺ → B + H⁺ (e.g. NH4⁺/NH3): protonated form is cation
                Ka = active_Kas[0]
                denom = Ka + H
                if denom > 0:
                    cations += CT * H / denom
            elif eq_def.category in ("acid", "inorganic_acid"):
                # Polyprotic acid: anion charge from active dissociations
                n = len(active_Kas)
                prods = [1.0]
                for Ka in active_Kas:
                    prods.append(prods[-1] * Ka)
                denom = sum(prods[i] * H ** (n - i) for i in range(n + 1))
                if denom > 0:
                    for i in range(1, n + 1):
                        alpha_i = prods[i] * H ** (n - i) / denom
                        anions += i * CT * alpha_i

        return max(0.0, anions - cations)

    # ── Deprecated wrapper ────────────────────────────────────────────

    def compute_CT_Na_from_charge_balance(self, pH, S_IC, S_IN,
                                           vfa_conc_mol_L, CT_Cl,
                                           T_K=None):
        """Compute CT_Na required to achieve a given pH.

        .. deprecated::
            Use :meth:`compute_CT_cation_from_charge_balance` instead,
            which respects ``n_active`` for each acid system and uses
            a generalised ``acid_totals`` dict.

        This wrapper converts the legacy arguments to an ``acid_totals``
        dict and delegates to
        :meth:`compute_CT_cation_from_charge_balance`.
        """
        acid_totals = {"CO2": S_IC, "NH4": S_IN}
        acid_totals.update(vfa_conc_mol_L)
        return self.compute_CT_cation_from_charge_balance(
            pH=pH, acid_totals=acid_totals, CT_anion=CT_Cl, T_K=T_K,
        )

    # ── Repr ──────────────────────────────────────────────────────────

    def __repr__(self):
        species = sorted(eq.name for eq in self._equilibria)
        corrections = {eq.name: eq.correction for eq in self._equilibria}
        return (f"ThermodynamicConfig(species={species}, "
                f"corrections={corrections}, "
                f"pKw={self._equilibria.water.pKw}, "
                f"T_default={self.T_ref_K} K)")


# ════════════════════════════════════════════════════════════════════════
#  Thermodynamic parameter validation
# ════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class ThermoSnapshot:
    """Frozen record of thermodynamic parameters found in one subsystem."""
    source: str
    species: str
    pKas: Tuple[float, ...]
    dH_pKas: Tuple[float, ...]
    n_active: Optional[int] = None


def _zone_label(obj):
    return getattr(obj, "label", None) or type(obj).__name__


def _find_gl_link(cv):
    """Return the KineticGasLiquidLink on a fermenter CV, or None.

    Works for both canonical (transfer_models kwarg) and legacy
    (internal_interfaces=[KineticGasLiquidLink(...)]) construction
    paths, since in both cases the link ends up in internal_interfaces.
    """
    from ..core.gas_liquid_link import KineticGasLiquidLink
    return next(
        (i for i in getattr(cv, "internal_interfaces", [])
         if isinstance(i, KineticGasLiquidLink)),
        None,
    )


def _collect_from_cv(cv, prefix=""):
    label = prefix or _zone_label(cv)
    snapshots = []
    tc = getattr(cv, "_thermo_config", None)
    if tc is not None:
        for eq_def in tc.equilibria:
            snapshots.append(ThermoSnapshot(
                source=f"{label}/thermo_config", species=eq_def.name,
                pKas=eq_def.pKas, dH_pKas=eq_def.dH_J_per_mol,
                n_active=eq_def.n_active))
    # Note: the transfer-link snapshot path (chemistry-unification-2)
    # was dropped along with link.speciation_corrections. The
    # property_solvers fallback path was dropped in state-unification
    # C4d -- pKas now live on declared EquilibriumReactions in
    # cv.reaction_system, collected by the eq_set walk above.
    # Phase 3b's ChemistryDatabase will rewrite this collection
    # entirely.
    return snapshots


def collect_thermo_params(model):
    """Extract thermodynamic parameters from all zones in a model.

    state-unification C4e: ``chem_env_fns`` kwarg removed. pKa /
    dH values flow exclusively from declared
    :class:`~PyOMES.reactions.equilibrium.EquilibriumReaction` instances
    in ``cv.reaction_system`` and from the engine's
    ``EquilibriumSet`` (collected by :func:`_collect_from_cv`).
    """
    snapshots = []
    cv_pairs = []
    model_type = type(model).__name__
    if model_type == "MultiCVSystem":
        for key, cv in model.cvs.items():
            if _find_gl_link(cv) is not None:
                cv_pairs.append((str(key), cv))
    elif isinstance(model, (list, tuple)):
        for i, item in enumerate(model):
            cv_pairs.append((_zone_label(item) or f"zone_{i}", item))
    elif _find_gl_link(model) is not None:
        cv_pairs.append((_zone_label(model), model))
    else:
        warnings.warn(f"validate_thermodynamics: unrecognised model type {model_type!r}", UserWarning, stacklevel=2)
    for label, cv in cv_pairs:
        snapshots.extend(_collect_from_cv(cv, prefix=label))
    return snapshots


def _check_within_zone(snapshots):
    issues = []
    by_zs = defaultdict(list)
    for s in snapshots:
        by_zs[(s.source.split("/")[0], s.species)].append(s)
    for (zone, species), group in by_zs.items():
        pg = [s for s in group if len(s.pKas) > 0]
        if len(pg) > 1 and len(set(s.pKas for s in pg)) > 1:
            issues.append(f"Within-zone pKa mismatch for '{species}' in '{zone}': " + "; ".join(f"{s.source}={s.pKas}" for s in pg))
        # ── dH mismatch (skip trivial zeros from pre-corrected contexts) ──
        dg = [s for s in group if any(abs(v) > 0.0 for v in s.dH_pKas)]
        if len(dg) > 1 and len(set(s.dH_pKas for s in dg)) > 1:
            issues.append(f"Within-zone dH mismatch for '{species}' in '{zone}': " + "; ".join(f"{s.source}={s.dH_pKas}" for s in dg))
        # ── n_active mismatch ─────────────────────────────────────────
        # If n_active is set on one subsystem but not another for the
        # same species, or set to different values, the speciation solver
        # and transfer link will compute different f_molecular fractions.
        na_group = [s for s in group if s.n_active is not None]
        if na_group:
            # Check: do all subsystems with n_active agree?
            na_vals = set(s.n_active for s in na_group)
            if len(na_vals) > 1:
                issues.append(
                    f"Within-zone n_active mismatch for '{species}' in "
                    f"'{zone}': "
                    + "; ".join(f"{s.source} n_active={s.n_active}"
                                for s in na_group))
            # Check: is any subsystem using more pKas than n_active allows?
            for s in group:
                if s.n_active is None and len(s.pKas) > 1:
                    # This subsystem uses all pKas, but another has n_active set
                    n_ref = na_group[0].n_active
                    if len(s.pKas) > n_ref:
                        issues.append(
                            f"n_active inconsistency for '{species}' in "
                            f"'{zone}': {s.source} uses all {len(s.pKas)} "
                            f"pKas but {na_group[0].source} limits to "
                            f"n_active={n_ref}. The transfer link will "
                            f"compute a different f_molecular than the "
                            f"speciation solver, causing systematic gas "
                            f"transfer bias.")
    return issues


def _check_across_zones(snapshots):
    issues = []
    fundamental = {"CO2", "NH4", "Kw"}
    by_sp = defaultdict(list)
    for s in snapshots:
        if s.species in fundamental:
            by_sp[s.species].append(s)
    for species, group in by_sp.items():
        zones = set(s.source.split("/")[0] for s in group)
        if len(zones) < 2:
            continue
        pg = [s for s in group if len(s.pKas) > 0]
        if len(pg) > 1 and len(set(s.pKas for s in pg)) > 1:
            issues.append(f"Across-zone pKa mismatch for '{species}': " + "; ".join(f"{s.source}={s.pKas}" for s in pg))
        dg = [s for s in group if any(abs(v) > 0.0 for v in s.dH_pKas)]
        if len(dg) > 1 and len(set(s.dH_pKas for s in dg)) > 1:
            issues.append(f"Across-zone dH mismatch for '{species}': " + "; ".join(f"{s.source}={s.dH_pKas}" for s in dg))
    return issues


def validate_thermodynamics(model, raise_on_error=False):
    """Check thermodynamic parameter consistency across all zones.

    Checks performed:

    1. **Within-zone pKa consistency**: all sources of pKas in a
       zone (the ``ThermodynamicConfig`` snapshot on the CV plus
       any declared ``EquilibriumReaction`` instances on
       ``cv.reaction_system``) must agree on each species' pKa
       values.

    2. **Within-zone dH consistency**: van 't Hoff enthalpies must
       match across sources. Catches the common error of different
       literature sources using different dH values (e.g. 7600 vs
       7646 for CO₂).

    3. **Within-zone n_active consistency**: if the speciation
       engine uses ``n_active=1`` for CO₂ (first dissociation
       only), every subsystem reading pKas for CO₂ must agree.

    4. **Across-zone consistency**: fundamental equilibria (CO₂,
       NH₄⁺, Kw) should use the same parameters in all reactor
       zones.

    Parameters
    ----------
    model : ControlVolume, MultiCVSystem, or list
        The model(s) to validate.
    raise_on_error : bool
        If True, raise ``ValueError`` on the first inconsistency.

    Returns
    -------
    list of str
        Warning messages (empty if fully consistent).

    Notes
    -----
    state-unification C4e removed the ``chem_env_fns`` parameter --
    chem_env is gone, and pKas now flow exclusively through declared
    ``EquilibriumReaction`` instances + the engine's
    ``EquilibriumSet``. The pre-Phase-2
    ``_check_f_molecular_consistency`` runtime cross-check is also
    removed -- the gas-liquid link no longer carries its own pKas,
    so there is no second source to drift against.
    """
    snapshots = collect_thermo_params(model)
    issues = []
    issues.extend(_check_within_zone(snapshots))
    issues.extend(_check_across_zones(snapshots))
    for issue in issues:
        if raise_on_error:
            raise ValueError(issue)
        warnings.warn(issue, UserWarning, stacklevel=2)
    return issues
