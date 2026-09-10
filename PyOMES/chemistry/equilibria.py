# -*- coding: utf-8 -*-
"""Equilibrium set: a declarative collection of acid-base equilibria.

An :class:`EquilibriumSet` formally describes which acid-base systems
participate in a speciation charge balance.  Each system is represented
by an :class:`EquilibriumDef` that carries its thermodynamic parameters
(pKa, dH, reference temperature, correction method) and charge balance
role (category, number of active dissociation steps).

Usage
-----
Build from scratch (approach 1 — blank + add):

>>> eq = EquilibriumSet()
>>> eq.set_water(pKw=14.0, correction="van_t_hoff", dH_J_per_mol=55900.0)
>>> eq.add("CO2", category="inorganic_acid",
...        pKas=(6.35, 10.33), n_active=1,
...        correction="van_t_hoff", dH_J_per_mol=(7646.0, 14900.0))
>>> eq.add("S_ac", category="acid", pKas=(4.76,))

Load a preset and modify (approach 2 — defaults + tweak):

>>> eq = EquilibriumSet.bsm2_default()
>>> eq.remove("S_va")
>>> eq.add("H2S", category="acid", pKas=(7.0,),
...        correction="van_t_hoff", dH=(20.0,), dH_unit="kJ/mol")

Categories
----------
Each equilibrium falls into one of these categories, which determines
how it enters the charge balance:

``"acid"``
    Weak acid (HA → A⁻ + H⁺).  Anionic species contribute negative
    charge.  Handles monoprotic and polyprotic.  Examples: VFAs, H₂S.

``"cation_acid"``
    Cationic acid (BH⁺ → B + H⁺).  The protonated form is a cation.
    Example: NH₄⁺/NH₃.

``"inorganic_acid"``
    Polyprotic inorganic acid.  Same algebra as ``"acid"`` but
    semantically distinct.  Examples: CO₂/carbonate, phosphate.

``"strong_ion"``
    Fully dissociated, no equilibrium.  Fixed charge.  Included for
    completeness; typically handled via the context dict.
    Examples: Na⁺, Cl⁻.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

# Reuse unit conversion helpers from thermo_params
from .thermo_params import (
    _R_J, _VALID_CORRECTIONS,
    _convert_dH_to_J, _convert_dH_tuple,
    _convert_T_to_K, _Ka_to_pKa, _lnKa_to_pKa,
)

_VALID_CATEGORIES = ("acid", "cation_acid", "inorganic_acid", "strong_ion")


# ════════════════════════════════════════════════════════════════════════
#  EquilibriumDef
# ════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class EquilibriumDef:
    """One equilibrium system in the charge balance.

    Parameters
    ----------
    name : str
        Unique identifier (e.g. ``"CO2"``, ``"S_ac"``, ``"phosphate"``).
    category : str
        Charge balance role: ``"acid"``, ``"cation_acid"``,
        ``"inorganic_acid"``, or ``"strong_ion"``.
    pKas : tuple of float
        Reference pKa value(s) at ``T_ref_K``.  Length equals the
        number of dissociation steps defined (``n_protons``).
    n_active : int
        Number of dissociation steps used in the charge balance.
        Must be ``<= len(pKas)``.  Defaults to ``len(pKas)``.

        This allows storing full thermodynamic data while using only
        a subset in the charge balance.  For example, BSM2 defines
        CO₂ as diprotic ``pKas=(6.35, 10.33)`` but sets ``n_active=1``
        so only the first dissociation (HCO₃⁻) participates in the
        charge balance.  The second pKa is still available for other
        calculations (e.g. gas transfer speciation corrections).
    T_ref_K : float
        Reference temperature (K) for the pKa values.
    correction : str
        Temperature correction method: ``"none"`` or ``"van_t_hoff"``.
    dH_J_per_mol : tuple of float
        Van 't Hoff enthalpy (J/mol) per dissociation step.
        Length must match ``pKas``.
    total_key : str
        Key in the context dict for the total concentration of this
        species.  E.g. ``"CT_TIC"`` for CO₂, ``"CT_S_ac"`` for acetate.
    """
    name: str
    category: str
    pKas: Tuple[float, ...]
    n_active: int
    T_ref_K: float = 298.15
    correction: str = "none"
    dH_J_per_mol: Tuple[float, ...] = ()
    total_key: str = ""
    # Ordered ladder of Species objects (most protonated → least protonated).
    # When non-empty, _compute_species_eq emits using Species.id directly
    # instead of synthesised {name}_HA / {name}_A- keys.
    # Set by EquilibriumSet.add(species_refs=...) or engine._add_acid.
    species_refs: Tuple[Any, ...] = ()

    def __post_init__(self):
        if not isinstance(self.pKas, tuple):
            object.__setattr__(self, "pKas", tuple(self.pKas))
        if not isinstance(self.dH_J_per_mol, tuple):
            object.__setattr__(self, "dH_J_per_mol", tuple(self.dH_J_per_mol))
        if len(self.dH_J_per_mol) == 0:
            object.__setattr__(self, "dH_J_per_mol", (0.0,) * len(self.pKas))
        if not isinstance(self.species_refs, tuple):
            object.__setattr__(self, "species_refs", tuple(self.species_refs))

    @property
    def n_protons(self) -> int:
        """Total number of dissociation steps defined."""
        return len(self.pKas)

    def pKas_at_T(self, T_K: float) -> Tuple[float, ...]:
        """Compute all pKa(s) at temperature T_K."""
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

    def active_pKas_at_T(self, T_K: float) -> Tuple[float, ...]:
        """Compute only the active pKa(s) at temperature T_K.

        Returns the first ``n_active`` entries from :meth:`pKas_at_T`.
        These are the values used in the charge balance.
        """
        return self.pKas_at_T(T_K)[:self.n_active]


# ════════════════════════════════════════════════════════════════════════
#  WaterDef (reused within EquilibriumSet)
# ════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class WaterDef:
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


# ════════════════════════════════════════════════════════════════════════
#  EquilibriumSet
# ════════════════════════════════════════════════════════════════════════

class EquilibriumSet:
    """Ordered collection of equilibrium definitions for a charge balance.

    Build from scratch or load a preset, then add/remove systems.
    """

    def __init__(self, *, T_ref_K: float = 298.15):
        """Create an empty set with a default reference temperature.

        Parameters
        ----------
        T_ref_K : float
            Default reference temperature (K) used when ``add()`` is
            called without specifying ``T_ref_K``.
        """
        self.T_ref_K = T_ref_K
        self._equilibria: Dict[str, EquilibriumDef] = {}
        self._water: WaterDef = WaterDef(pKw=14.0, T_ref_K=T_ref_K)
        self._insertion_order: List[str] = []

    # ── Water ─────────────────────────────────────────────────────────

    def set_water(self, *, pKw=None, Kw=None, correction="none",
                  dH=None, dH_unit="J/mol", dH_J_per_mol=None,
                  T_ref=None, T_unit="K", T_ref_K=None):
        """Define water autoionisation.

        Accepts the same unit flexibility as :meth:`add`.

        Returns self for method chaining.
        """
        import warnings

        if pKw is not None and Kw is not None:
            raise ValueError("set_water: provide only one of pKw or Kw")
        if Kw is not None:
            if Kw <= 0:
                raise ValueError(f"Kw must be positive, got {Kw}")
            pKw = -math.log10(Kw)
        elif pKw is None:
            pKw = 14.0

        if T_ref_K is not None:
            t_ref = float(T_ref_K)
        elif T_ref is not None:
            t_ref = _convert_T_to_K(T_ref, T_unit)
        else:
            t_ref = self.T_ref_K

        if dH_J_per_mol is not None:
            dH_internal = float(dH_J_per_mol)
        elif dH is not None:
            dH_internal = _convert_dH_to_J(dH, dH_unit)
        else:
            dH_internal = 0.0

        correction = correction.lower().strip()
        if correction not in _VALID_CORRECTIONS:
            raise ValueError(
                f"Unknown correction {correction!r} for water. "
                f"Options: {_VALID_CORRECTIONS}")
        if correction == "van_t_hoff" and abs(dH_internal) < 1e-10:
            raise ValueError(
                "correction='van_t_hoff' for water requires nonzero dH")

        self._water = WaterDef(
            pKw=pKw, T_ref_K=t_ref,
            correction=correction, dH_J_per_mol=dH_internal)
        return self

    @property
    def water(self) -> WaterDef:
        """The water autoionisation definition."""
        return self._water

    # ── Add / remove / query ──────────────────────────────────────────

    def add(self, name, *, category, pKas=None, Ka=None, lnKa=None,
            n_active=None, correction="none",
            dH=None, dH_unit="J/mol", dH_J_per_mol=None,
            T_ref=None, T_unit="K", T_ref_K=None,
            total_key=None, species_refs=()):
        """Add or replace an equilibrium definition.

        Parameters
        ----------
        name : str
            Unique identifier for this equilibrium system.
        category : str
            ``"acid"``, ``"cation_acid"``, ``"inorganic_acid"``,
            or ``"strong_ion"``.
        pKas, Ka, lnKa : tuple of float
            Dissociation constant(s).  Provide exactly one format.
        n_active : int, optional
            Number of dissociation steps active in the charge balance.
            Defaults to ``len(pKas)``.  Set lower to store inactive
            pKas as metadata (e.g. ``n_active=1`` for BSM2 monoprotic
            CO₂ while storing both pKa₁ and pKa₂).
        correction : str
            ``"none"`` or ``"van_t_hoff"``.
        dH, dH_unit : tuple + str
            Enthalpy values with unit.  See :meth:`ThermodynamicConfig.add_acid`.
        dH_J_per_mol : tuple of float
            Enthalpy in J/mol (takes precedence over dH + dH_unit).
        T_ref, T_unit : float + str
            Reference temperature with unit.
        T_ref_K : float
            Reference temperature in K (takes precedence).
        total_key : str, optional
            Context dict key for total concentration.  Defaults to
            ``"CT_{name}"`` (e.g. ``"CT_CO2"``, ``"CT_S_ac"``).

        Returns
        -------
        self

        Examples
        --------
        >>> eq.add("CO2", category="inorganic_acid",
        ...        pKas=(6.35, 10.33), n_active=1,
        ...        correction="van_t_hoff",
        ...        dH=(7.646, 14.9), dH_unit="kJ/mol")
        >>> eq.add("H2S", category="acid",
        ...        Ka=(1.07e-7,), correction="van_t_hoff",
        ...        dH=(20.0,), dH_unit="kJ/mol",
        ...        T_ref=20, T_unit="C")
        >>> eq.add("S_ac", category="acid", pKas=(4.76,))
        """
        import warnings

        # ── Validate category ─────────────────────────────────────
        category = category.lower().strip()
        if category not in _VALID_CATEGORIES:
            raise ValueError(
                f"Unknown category {category!r} for '{name}'. "
                f"Options: {_VALID_CATEGORIES}")

        # ── Convert dissociation constants to pKa ─────────────────
        if category == "strong_ion":
            pKas_internal = ()
        else:
            n_provided = sum(x is not None for x in (pKas, Ka, lnKa))
            if n_provided == 0:
                raise ValueError(
                    f"add('{name}'): must provide one of pKas, Ka, or lnKa")
            if n_provided > 1:
                raise ValueError(
                    f"add('{name}'): provide only one of pKas, Ka, or lnKa")

            if Ka is not None:
                if not isinstance(Ka, (tuple, list)):
                    Ka = (Ka,)
                pKas_internal = tuple(_Ka_to_pKa(k) for k in Ka)
            elif lnKa is not None:
                if not isinstance(lnKa, (tuple, list)):
                    lnKa = (lnKa,)
                pKas_internal = tuple(_lnKa_to_pKa(lk) for lk in lnKa)
            else:
                if not isinstance(pKas, (tuple, list)):
                    pKas = (pKas,)
                pKas_internal = tuple(float(p) for p in pKas)

        # ── n_active ──────────────────────────────────────────────
        if n_active is None:
            n_active = len(pKas_internal)
        else:
            n_active = int(n_active)
            if n_active < 0:
                raise ValueError(
                    f"n_active must be >= 0, got {n_active} for '{name}'")
            if n_active > len(pKas_internal):
                raise ValueError(
                    f"n_active ({n_active}) cannot exceed number of pKas "
                    f"({len(pKas_internal)}) for '{name}'")

        # ── Convert temperature ───────────────────────────────────
        if T_ref_K is not None:
            t_ref = float(T_ref_K)
        elif T_ref is not None:
            t_ref = _convert_T_to_K(T_ref, T_unit)
        else:
            t_ref = self.T_ref_K

        # ── Convert enthalpy ──────────────────────────────────────
        if category == "strong_ion":
            dH_internal = ()
        elif dH_J_per_mol is not None:
            if not isinstance(dH_J_per_mol, (tuple, list)):
                dH_J_per_mol = (dH_J_per_mol,)
            dH_internal = tuple(float(d) for d in dH_J_per_mol)
        elif dH is not None:
            dH_internal = _convert_dH_tuple(dH, dH_unit)
        else:
            dH_internal = None

        # ── Validate correction + dH ─────────────────────────────
        correction = correction.lower().strip()
        if correction not in _VALID_CORRECTIONS:
            raise ValueError(
                f"Unknown correction {correction!r} for '{name}'. "
                f"Options: {_VALID_CORRECTIONS}")

        if category != "strong_ion":
            if correction == "van_t_hoff":
                if dH_internal is None:
                    raise ValueError(
                        f"correction='van_t_hoff' for '{name}' requires dH")
                if len(dH_internal) != len(pKas_internal):
                    raise ValueError(
                        f"dH length ({len(dH_internal)}) must match pKas "
                        f"length ({len(pKas_internal)}) for '{name}'")
            elif correction == "none":
                if dH_internal is not None and any(
                        abs(d) > 1e-10 for d in dH_internal):
                    warnings.warn(
                        f"add('{name}'): dH provided but correction='none' "
                        f"— values will be ignored.",
                        UserWarning, stacklevel=2)
                dH_internal = (0.0,) * len(pKas_internal)

        # ── Default total_key ─────────────────────────────────────
        if total_key is None:
            total_key = f"CT_{name}"

        # ── Store ─────────────────────────────────────────────────
        is_new = name not in self._equilibria
        self._equilibria[name] = EquilibriumDef(
            name=name, category=category,
            pKas=pKas_internal, n_active=n_active,
            T_ref_K=t_ref, correction=correction,
            dH_J_per_mol=dH_internal if dH_internal is not None else (),
            total_key=total_key,
            species_refs=tuple(species_refs))
        if is_new:
            self._insertion_order.append(name)
        return self

    def remove(self, name: str) -> "EquilibriumSet":
        """Remove an equilibrium by name.

        Raises
        ------
        KeyError
            If ``name`` is not in the set.
        """
        if name not in self._equilibria:
            raise KeyError(
                f"Equilibrium '{name}' not found. "
                f"Available: {self.names}")
        del self._equilibria[name]
        self._insertion_order.remove(name)
        return self

    def has(self, name: str) -> bool:
        """Check if an equilibrium is registered."""
        return name in self._equilibria

    def get(self, name: str) -> EquilibriumDef:
        """Get an equilibrium definition by name.

        Raises
        ------
        KeyError
            If ``name`` is not in the set.
        """
        if name not in self._equilibria:
            raise KeyError(
                f"Equilibrium '{name}' not found. "
                f"Available: {self.names}")
        return self._equilibria[name]

    @property
    def names(self) -> List[str]:
        """List of registered equilibrium names (insertion order)."""
        return list(self._insertion_order)

    def __len__(self) -> int:
        return len(self._equilibria)

    def __contains__(self, name: str) -> bool:
        return name in self._equilibria

    def __iter__(self):
        """Iterate over EquilibriumDef objects in insertion order."""
        for name in self._insertion_order:
            yield self._equilibria[name]

    def __getitem__(self, name: str) -> EquilibriumDef:
        return self.get(name)

    # ── Factory presets ───────────────────────────────────────────────

    @staticmethod
    def bsm2_default() -> "EquilibriumSet":
        """BSM2-canonical ADM1 equilibria (Rosen & Jeppsson 2006).

        CO₂ is diprotic but only the first dissociation is active in
        the charge balance (``n_active=1``), matching BSM2's monoprotic
        treatment.  The second pKa is stored for use by gas transfer
        speciation corrections.

        VFA pKas have no temperature correction.
        """
        from .common_species import (
            CO2, HCO3_minus, CO3_2minus, NH4_plus, NH3,
        )
        eq = EquilibriumSet(T_ref_K=298.15)
        eq.set_water(pKw=14.0, correction="van_t_hoff",
                     dH_J_per_mol=55900.0)
        eq.add("CO2", category="inorganic_acid",
               pKas=(6.35, 10.33), n_active=1,
               correction="van_t_hoff",
               dH_J_per_mol=(7646.0, 14900.0),
               total_key="CT_TIC",
               species_refs=(CO2, HCO3_minus, CO3_2minus))
        eq.add("NH4", category="cation_acid",
               pKas=(9.25,), correction="van_t_hoff",
               dH_J_per_mol=(51965.0,),
               total_key="CT_NH_T",
               species_refs=(NH4_plus, NH3))
        # DEPRECATED: string-based VFA entries without species_refs.
        # The _HA/_A- legacy fallback is removed in the PARTITION_MODEL phase.
        eq.add("S_ac",  category="acid", pKas=(4.76,))
        eq.add("S_pro", category="acid", pKas=(4.88,))
        eq.add("S_bu",  category="acid", pKas=(4.82,))
        eq.add("S_va",  category="acid", pKas=(4.86,))
        return eq

    @staticmethod
    def bsm2_diprotic_co2() -> "EquilibriumSet":
        """BSM2 with full diprotic CO₂ (both dissociations active).

        Identical to :meth:`bsm2_default` except ``n_active=2`` for CO₂,
        so CO₃²⁻ participates in the charge balance.
        """
        eq = EquilibriumSet.bsm2_default()
        # Replace CO2 with n_active=2, preserving species_refs.
        co2 = eq.get("CO2")
        eq.add("CO2", category=co2.category,
               pKas=co2.pKas, n_active=2,
               correction=co2.correction,
               dH_J_per_mol=co2.dH_J_per_mol,
               T_ref_K=co2.T_ref_K,
               total_key=co2.total_key,
               species_refs=co2.species_refs)
        return eq

    @staticmethod
    def bsm2_with_sulfide() -> "EquilibriumSet":
        """BSM2 + H₂S/HS⁻ equilibrium."""
        eq = EquilibriumSet.bsm2_default()
        eq.add("H2S", category="acid", pKas=(7.0,),
               correction="van_t_hoff",
               dH_J_per_mol=(20000.0,))
        return eq

    @staticmethod
    def adm1_full() -> "EquilibriumSet":
        """Full ADM1: BSM2 (diprotic CO₂) + phosphate + bisulfate."""
        from .common_species import (
            H3PO4, H2PO4_minus, HPO4_2minus, PO4_3minus,
            HSO4_minus, SO4_2minus,
        )
        eq = EquilibriumSet.bsm2_diprotic_co2()
        eq.add("phosphate", category="inorganic_acid",
               pKas=(2.15, 7.20, 12.35), correction="none",
               total_key="CT_P",
               species_refs=(H3PO4, H2PO4_minus, HPO4_2minus, PO4_3minus))
        eq.add("bisulfate", category="inorganic_acid",
               pKas=(1.99,), correction="none",
               total_key="CT_SO4",
               species_refs=(HSO4_minus, SO4_2minus))
        return eq

    # ── Repr ──────────────────────────────────────────────────────────

    def __repr__(self) -> str:
        entries = []
        for eq_def in self:
            active_str = (f"{eq_def.n_active}/{eq_def.n_protons}"
                          if eq_def.n_active < eq_def.n_protons
                          else str(eq_def.n_protons))
            entries.append(
                f"{eq_def.name}({eq_def.category}, "
                f"pKas={eq_def.pKas}, active={active_str}, "
                f"correction={eq_def.correction!r})")
        water_str = (f"water(pKw={self._water.pKw}, "
                     f"correction={self._water.correction!r})")
        return (f"EquilibriumSet([{water_str}, "
                + ", ".join(entries) + "])")

    def summary(self) -> str:
        """Human-readable summary of all registered equilibria."""
        lines = [f"EquilibriumSet ({len(self)} equilibria, "
                 f"T_ref_default={self.T_ref_K} K):"]
        lines.append(f"  Water: pKw={self._water.pKw} at "
                     f"{self._water.T_ref_K} K, "
                     f"correction={self._water.correction!r}")
        for eq_def in self:
            active_str = (f" (active: {eq_def.n_active}/{eq_def.n_protons})"
                          if eq_def.n_active < eq_def.n_protons else "")
            lines.append(
                f"  {eq_def.name:12s}  {eq_def.category:18s}  "
                f"pKas={eq_def.pKas}{active_str}  "
                f"T_ref={eq_def.T_ref_K} K  "
                f"correction={eq_def.correction!r}  "
                f"key={eq_def.total_key!r}")
        return "\n".join(lines)
