# -*- coding: utf-8 -*-
"""Newton-Raphson speciation engine.

:class:`NRChemicalEquilibriumEngine` is a drop-in alternative to
:class:`~PyOMES.chemical_equilibrium.engine.BisectionChemicalEquilibriumEngine` that solves aqueous
equilibrium chemistry via a full Newton-Raphson system in log-activity
space rather than the 1-D charge-balance bisection used by the existing
engine.

The two engines expose the same ``.solve(phases=...)`` interface and
emit the same output-dict format.  The NR engine is selected by
passing ``solver="newton_raphson"`` to :class:`~PyOMES.reactions.reaction_system.ReactionSystem`.

Key differences vs the existing engine
---------------------------------------
- The NR engine handles **arbitrary reaction networks** (cross-component
  species, metal complexation) once the chemistry is declared.
- H⁺, OH⁻, and H₂O are written to ``phase.n_mol`` as **engine-owned
  derived species** after every solve (the existing engine only writes
  H⁺ and OH⁻).
- Activity correction is an outer fixed-point loop (identical to the
  existing engine) wrapping the inner NR loop.
- Warmstarting uses the full log-activity vector, not just log[H⁺].

Scope (this implementation)
----------------------------
- Single liquid phase; acid-base and metal-complexation networks.
- Precipitation/phase detection: deferred (see design doc).
- Redox (pe as a master): deferred.
"""
from __future__ import annotations

import logging
import warnings
from typing import Any, Dict, FrozenSet, List, Optional

import numpy as np
from scipy.sparse import csr_matrix

from .nr_tableau import NRTableau, _vant_hoff_log_K, build_tableau
from .nr_solver import (
    NRSolverCache, _gamma_safe, solve_nr,
    _build_gammas, _compute_concentrations, _residual_and_jacobian,
)
from .activity_models import make_activity_model
from .protocols import EquilibriumResult, SparseJacobian, SpeciationJacobian

logger = logging.getLogger(__name__)

# Water molar mass and density for H₂O writeback
_M_WATER_G_MOL = 18.015
_RHO_WATER_G_L = 1000.0   # g/L at 25 °C
_C_WATER_MOL_L = _RHO_WATER_G_L / _M_WATER_G_MOL   # ≈ 55.51 mol/L

# Canonical strong-ion species IDs → CT_* keys used for reading from phases
_STRONG_ION_SPECIES_TO_KEY: Dict[str, str] = {
    "K+":           "CT_K",
    "Na+":          "CT_Na",
    "Cl-":          "CT_Cl",
    "NO3-":         "CT_NO3",
    "Mg++":         "CT_Mg",
    "Ca++":         "CT_Ca",
    "Zn++":         "CT_Zn",
    "Mn++":         "CT_Mn",
    "Cu++":         "CT_Cu",
    "Co++":         "CT_Co",
    # Ferrous iron is a dissolved, divalent spectator until an Fe(II)
    # complexation/hydrolysis ladder is declared.  It must still enter
    # electroneutrality while Fe(III) is handled by its own component.
    "Fe2+":         "CT_Fe2",
    "Mo7O24------": "CT_Mo7O24",
    "MoO4--":       "CT_MoO4",
    "S_cat":        "CT_cation",
    "S_an":         "CT_anion",
}


_JACOBIAN_DISABLED_MSG = (
    "Jacobian retention is disabled. "
    "Construct NRChemicalEquilibriumEngine with retain_jacobian=True to enable white-box methods."
)

_LN10 = float(np.log(10.0))

# Strong-ion charge dictionary (mirrors the one in nr_solver.py)
_STRONG_CHARGES: Dict[str, int] = {
    "CT_K": +1, "CT_Na": +1, "CT_cation": +1,
    "CT_Cl": -1, "CT_NO3": -1, "CT_anion": -1,
    "CT_Mg": +2, "CT_Ca": +2, "CT_Zn": +2, "CT_Mn": +2, "CT_Cu": +2, "CT_Co": +2,
    "CT_Fe2": +2,
    "CT_Mo7O24": -6, "CT_MoO4": -2,
}


class NRChemicalEquilibriumEngine:
    """Newton-Raphson aqueous speciation engine.

    Construct via :meth:`from_reactions`; subsequent ``solve(**kwargs)``
    calls use the pre-built :class:`~PyOMES.chemical_equilibrium.nr_tableau.NRTableau`.

    Parameters
    ----------
    tableau : NRTableau
        Pre-built log-linear tableau.
    use_activity : bool
        Apply Davies activity-coefficient corrections.  Default ``False``.
    activity_model : str
        ``"davies"`` (default), ``"ideal"``.
    T_C : float
        Operating temperature (°C).  Default ``25.0``.
    use_warmstart : bool
        Reuse the last converged state as initial guess.  Default ``True``.
    retain_jacobian : bool
        Cache the converged NR Jacobian matrix after each ``solve()`` call,
        unlocking the white-box methods ``jacobian_dg_dz()``,
        ``jacobian_dg_dy()``, ``jacobian_dz_dy()``, and ``residual()``.
        When ``False`` (default) those methods raise ``RuntimeError``.
    """

    def __init__(
        self,
        tableau: NRTableau,
        *,
        use_activity: bool = False,
        activity_model: str = "davies",
        T_C: float = 25.0,
        use_warmstart: bool = True,
        max_log_activity: float = 50.0,
        min_component_total: float = 1e-20,
        retain_jacobian: bool = False,
        thermo=None,
    ):
        self._tableau = tableau
        if retain_jacobian and any(sec.phase == "gas" for sec in tableau.secondaries):
            raise NotImplementedError(
                "NRChemicalEquilibriumEngine: retain_jacobian=True is not yet "
                "supported for a tableau with folded gas-liquid secondaries "
                "(LAYER1_GAP_CLOSURE CP1/CP2). The white-box Jacobian "
                "machinery (residual(), jacobian_dg_dz(), jacobian_dz_dy()) "
                "keys its internal caches by bare species_id, which "
                "collides whenever a gas secondary shares an id with its "
                "liquid parent (the common Henry case, e.g. both 'CO2') — "
                "CP3 confirmed algebraic_species() itself needed no change "
                "(it already unions bare ids across phases correctly for "
                "set-membership use), but did not thread a phase-aware "
                "distinction through the white-box methods, so this guard "
                "stays in place. The main solve() path is unaffected."
            )
        if thermo is not None:
            self._liquid_activity = thermo.liquid_activity
            use_activity = thermo.use_activity
            activity_model = thermo.activity_model
        else:
            self._liquid_activity = None
        self.use_activity = bool(use_activity)
        self.activity_model = str(activity_model)
        self.T_C = float(T_C)
        self.use_warmstart = bool(use_warmstart)
        self.max_log_activity = float(max_log_activity)
        self.min_component_total = float(min_component_total)
        self.retain_jacobian = bool(retain_jacobian)

        self._cache = NRSolverCache() if use_warmstart else None
        self.n_solve_calls = 0
        self._precipitation_reactions: List = []  # set by from_reactions()

        # White-box cached state (populated when retain_jacobian=True)
        self._cached_jacobian: Optional[np.ndarray] = None
        self._cached_concentrations: Optional[Dict[str, float]] = None
        self._cached_x: Optional[np.ndarray] = None
        self._cached_I: Optional[float] = None

    # ── Factory ───────────────────────────────────────────────────────

    @classmethod
    def from_reactions(
        cls,
        equilibrium_reactions,
        *,
        precipitation_reactions=None,
        use_activity: bool = False,
        activity_model: str = "davies",
        T_K: float = 298.15,
        use_warmstart: bool = True,
        max_log_activity: float = 50.0,
        min_component_total: float = 1e-20,
        retain_jacobian: bool = False,
    ) -> "NRChemicalEquilibriumEngine":
        """Build an :class:`NRChemicalEquilibriumEngine` from declared equilibrium reactions.

        Parameters
        ----------
        equilibrium_reactions : iterable of EquilibriumConstraint
            One flat list of all equilibrium constraints for the system —
            typically ``EquilibriumReaction``, but any
            ``EquilibriumConstraint``-conforming item is accepted (e.g.
            ``KspEquilibrium``). Each item is classified via
            :func:`~PyOMES.reactions.equilibrium.classify_equilibrium_constraint`:
            single-phase (acid-base) items feed the NR tableau; solid-liquid
            items (one ``StoichiometryEntry(phase="solid")`` — the mineral —
            plus one or more ``StoichiometryEntry(phase="liquid")`` dissolved
            ionic products) are auto-detected as precipitation reactions,
            no separate kwarg required. Must include a water-dissociation
            reaction. The ``log_K`` on a precipitation item is the Ksp
            (solid activity = 1 by convention).  When any are present,
            :meth:`solve` runs an outer active-set loop so that IAP = Ksp
            for every active mineral.  The ``"minerals"`` key in the output
            dict reports ``xi_mol_L`` and ``SI`` for each mineral.

            .. note::
                For CV/SolidPhase integration (ξ writeback to n_mol, lazy phase
                creation) see NR_PRECIPITATION_CV_INTEGRATION.md (Phase 2).
        precipitation_reactions : list of EquilibriumReaction, optional
            Deprecated — solid-liquid items are now auto-classified from
            ``equilibrium_reactions``. Kept for one phase
            (``EQUILIBRIUM_CONSTRAINT_UNIFICATION`` CP2); emits
            ``DeprecationWarning`` and is merged with the auto-detected set.
        use_activity : bool
            Enable Davies activity corrections.  Default ``False``.
        activity_model : str
            Activity model name (``"davies"`` or ``"ideal"``).
        T_K : float
            Operating temperature (K).  Van't Hoff correction is applied
            to all log_K values at this temperature.
        use_warmstart : bool
            Keep a warmstart cache across ``solve()`` calls.
        max_log_activity : float
            Passed to :func:`~PyOMES.chemical_equilibrium.nr_solver.solve_nr`.
            See :meth:`__init__` for details.
        min_component_total : float
            Passed to :func:`~PyOMES.chemical_equilibrium.nr_solver.solve_nr`.
            See :meth:`__init__` for details.

        Returns
        -------
        NRChemicalEquilibriumEngine
        """
        from ..reactions.equilibrium import (
            EquilibriumConstraint, classify_equilibrium_constraint,
        )

        equilibrium_reactions = list(equilibrium_reactions)
        auto_precipitation: List = []
        for item in equilibrium_reactions:
            if not isinstance(item, EquilibriumConstraint):
                continue
            try:
                kind = classify_equilibrium_constraint(item)
            except ValueError:
                continue
            if kind == "solid_liquid":
                auto_precipitation.append(item)

        if precipitation_reactions:
            warnings.warn(
                "NRChemicalEquilibriumEngine.from_reactions(precipitation_reactions=...) "
                "is deprecated; solid-liquid items are now auto-classified "
                "from equilibrium_reactions via StoichiometryEntry(phase="
                "'solid') tags (EQUILIBRIUM_CONSTRAINT_UNIFICATION CP2). "
                "Pass them directly in equilibrium_reactions instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            auto_precipitation = auto_precipitation + list(precipitation_reactions)

        tableau = build_tableau(equilibrium_reactions, T_K=T_K)
        T_C = float(T_K) - 273.15
        engine = cls(
            tableau,
            use_activity=use_activity,
            activity_model=activity_model,
            T_C=T_C,
            use_warmstart=use_warmstart,
            max_log_activity=max_log_activity,
            min_component_total=min_component_total,
            retain_jacobian=retain_jacobian,
        )
        if auto_precipitation:
            engine._precipitation_reactions = auto_precipitation
        return engine

    # ── Public solve interface ────────────────────────────────────────

    def solve(self, **kwargs) -> EquilibriumResult:
        """Solve the NR speciation system.

        Accepts either of two call patterns:

        **Phase-based** (primary path within ``cv.advance()``)::

            result = engine.solve(phases={"liquid": liq_phase})
            result.apply_to_phases({"liquid": liq_phase})  # explicit commit

        **Direct** (testing and scripting)::

            result = engine.solve(
                totals={"CO2": 0.05, "NH3": 0.04},   # mol/L, keyed by master id
                strong_ions={"CT_Na": 0.01},          # mol/L
                T_K=308.15,                            # optional override
            )

        In the phase-based path, ``totals`` are read by summing ``n_mol``
        across all species in each connected component.  ``strong_ions`` are
        read from ``n_mol`` via the canonical :data:`_STRONG_ION_SPECIES_TO_KEY`
        mapping. ``solve()`` makes no side-effecting writes to phases —
        call ``result.apply_to_phases(phases)`` explicitly to commit.

        Parameters
        ----------
        phases : dict, optional
            ``{"liquid": phase_obj, "gas": phase_obj}`` where each
            ``phase_obj`` has ``.n_mol`` (``dict[str, float]``) and
            ``.V_L`` (float, litres). ``"gas"`` is required only if the
            tableau folds any gas-liquid equilibrium (CP1/CP2 of
            ``LAYER1_GAP_CLOSURE``); its totals and volume feed the
            volume-aware mass balance for those components.
        totals : dict, optional
            Direct override: ``{master_id: C_total (mol/L)}``. For a
            folded component this must already be the total *across both
            phases*, divided by ``V_liq_L`` — see ``V_liq_L``/``V_gas_L``.
        strong_ions : dict, optional
            Direct override: ``{CT_Na: mol/L, ...}``.
        T_K : float, optional
            Temperature override (K).
        V_liq_L, V_gas_L : float, optional
            Direct-pattern only: liquid/gas volumes (litres), required
            (both, > 0) when the tableau folds any gas-liquid equilibrium.
            Ignored in the phase-based pattern (read from ``phases``
            instead).

        Returns
        -------
        EquilibriumResult
            ``pH``, ``pH_conc``, ``logH``, ``aH``, ``gamma_H``, ``gamma_OH``,
            ``ionic_strength``, ``charge_residual`` are the common meta
            fields. ``species_mol_L`` carries one entry per liquid-phase
            species in the tableau (plus H2O); ``partial_pressures_atm``
            carries one entry per folded gas-phase secondary (atm, not
            mol/L — kept out of ``species_mol_L`` to avoid mislabeling a
            pressure as a concentration). When precipitation reactions are
            declared, ``saturation_indices`` carries per-mineral SI and
            ``extra["minerals_xi_mol_L"]`` carries per-mineral ξ (mol/L
            precipitated).

        .. note::
            For CV/SolidPhase integration (writing ξ back to ``n_mol``,
            lazy SolidPhase creation, ODE state update) see
            ``NR_PRECIPITATION_CV_INTEGRATION.md`` (Phase 2).
        """
        phases = kwargs.pop("phases", None)
        T_K_override = kwargs.pop("T_K", None)
        T_K = float(T_K_override) if T_K_override is not None else 273.15 + self.T_C

        # -- Extract totals and strong ions ----------------------------
        if phases is not None:
            totals, strong_ions, V_liq_L, V_gas_L = self._read_from_phases(phases)
        else:
            totals = dict(kwargs.pop("totals", {}))
            strong_ions = dict(kwargs.pop("strong_ions", {}))
            V_liq_L = kwargs.pop("V_liq_L", None)
            V_gas_L = kwargs.pop("V_gas_L", None)

        am = (
            self._liquid_activity
            if self._liquid_activity is not None
            else make_activity_model(self.use_activity, self.activity_model)
        )

        self.n_solve_calls += 1

        if self._precipitation_reactions:
            out = self._solve_with_precipitation(
                totals, strong_ions, T_K, am, V_liq_L=V_liq_L, V_gas_L=V_gas_L,
            )
        else:
            out = solve_nr(
                self._tableau,
                totals,
                strong_ions,
                T_K=T_K,
                activity_model=am,
                cache=self._cache if self.use_warmstart else None,
                max_log_activity=self.max_log_activity,
                min_component_total=self.min_component_total,
                retain_jacobian=self.retain_jacobian,
                V_liq_L=V_liq_L,
                V_gas_L=V_gas_L,
            )

        if self.retain_jacobian:
            self._cached_jacobian = out.pop("_jacobian_matrix", None)
            self._cached_x = out.pop("_log_activities", None)
            self._cached_I = out.pop("_ionic_strength_final", None)
            # Snapshot all species concentrations for jacobian_dz_dy().
            # Bare species_id is safe here (never sec.c_key) because
            # __init__ already refuses retain_jacobian=True for a tableau
            # with any gas secondary — see that guard's docstring.
            self._cached_concentrations = {
                sp_id: float(out[sp_id])
                for sp_id in (
                    list(self._tableau.masters)
                    + [sec.species_id for sec in self._tableau.secondaries]
                )
                if sp_id in out
            }

        # -- Build EquilibriumResult -------------------------------------
        # species_mol_L is the full masters + liquid secondaries + H2O set
        # — mirrors what this engine has always written back (no canonical
        # restriction, unlike BisectionChemicalEquilibriumEngine, since the NR tableau owns
        # its full species set). Gas-phase secondaries (CP1/CP2 of
        # LAYER1_GAP_CLOSURE) go to partial_pressures_atm instead — their
        # value is a pressure (atm), not a concentration (mol/L), and
        # sec.c_key (not sec.species_id) is required to read them out of
        # `out` without colliding with a same-named liquid entry.
        species_mol_L: Dict[str, float] = {}
        partial_pressures_atm: Dict[str, float] = {}
        for sp_id in self._tableau.masters:
            if sp_id in out and out[sp_id] is not None:
                try:
                    species_mol_L[sp_id] = float(out[sp_id])
                except (TypeError, ValueError):
                    pass
        for sec in self._tableau.secondaries:
            if sec.c_key not in out or out[sec.c_key] is None:
                continue
            try:
                value = float(out[sec.c_key])
            except (TypeError, ValueError):
                continue
            if sec.phase == "gas":
                partial_pressures_atm[sec.species_id] = value
            else:
                species_mol_L[sec.species_id] = value
        species_mol_L["H2O"] = _C_WATER_MOL_L

        minerals = out.get("minerals", {}) or {}
        saturation_indices = {
            label: float(info["SI"]) for label, info in minerals.items()
        }
        extra: Dict[str, Any] = {}
        if minerals:
            extra["minerals_xi_mol_L"] = {
                label: float(info["xi_mol_L"]) for label, info in minerals.items()
            }

        pH_val = out.get("pH")
        result = EquilibriumResult(
            pH=float(pH_val) if pH_val is not None else float("nan"),
            pH_conc=out.get("pH_conc"),
            logH=out.get("logH"),
            aH=out.get("aH"),
            gamma_H=out.get("gamma_H"),
            gamma_OH=out.get("gamma_OH"),
            ionic_strength=out.get("IonicStrength"),
            charge_residual=out.get("charge_residual"),
            species_mol_L=species_mol_L,
            partial_pressures_atm=partial_pressures_atm,
            saturation_indices=saturation_indices,
            extra=extra,
        )
        return result

    # ── BisectionChemicalEquilibriumEngine-compatible helpers ───────────────────────────

    def algebraic_species(self) -> FrozenSet[str]:
        """Species IDs written to ``phase.n_mol`` by this engine.

        These form the algebraic state z (complement of the differential
        state y) in any DAE formulation. Already includes folded gas-
        liquid species (CP1/CP2 of ``LAYER1_GAP_CLOSURE``) alongside
        ordinary acid-base and solid-liquid species — the underlying set
        union is over bare ``species_id`` regardless of
        ``SecondaryEntry.phase``, so no change was needed here to pick
        them up (verified by CP3; see :meth:`gas_liquid_species` for the
        narrower subset that actually has a folded gas-liquid row).
        """
        ids = set(self._tableau.masters)
        ids.update(sec.species_id for sec in self._tableau.secondaries)
        ids.add("H2O")
        return frozenset(ids)

    def gas_liquid_species(self) -> FrozenSet[str]:
        """Species with a folded gas-liquid row in this engine's tableau.

        A strict subset of :meth:`algebraic_species` — deliberately
        *not* every species this engine resolves (that would also
        include ordinary acid-base species like liquid CO2/HCO3-/CO3--
        that have no gas-liquid coupling declared at all). Used by
        :meth:`~PyOMES.core.control_volume.ControlVolume.step_internal_transfer`
        (CP3 of ``LAYER1_GAP_CLOSURE``) to skip a separately-declared
        ``transfer_models`` entry only for species this engine already
        resolves simultaneously with acid-base — for a species with
        ordinary acid-base chemistry but no folded Henry/Raoult row, a
        ``transfer_models`` entry is still the correct (and only)
        mechanism, so it must not be filtered.

        Not part of :class:`~PyOMES.chemical_equilibrium.protocols.ChemicalEquilibriumEngineProtocol`
        — only :class:`NRChemicalEquilibriumEngine` currently supports gas-liquid
        folding; callers duck-type via ``getattr(engine,
        "gas_liquid_species", None)``.
        """
        return frozenset(
            sec.species_id for sec in self._tableau.secondaries
            if sec.phase == "gas"
        )

    def reset_cache(self) -> None:
        """Reset warmstart cache and any retained white-box state."""
        if self._cache is not None:
            self._cache.log_x = None
            self._cache.I_last = None
            self._cache.jacobian = None
        self._cached_jacobian = None
        self._cached_concentrations = None
        self._cached_x = None
        self._cached_I = None

    def reset_counters(self) -> None:
        """Reset solve-call counter."""
        self.n_solve_calls = 0

    @property
    def tableau(self) -> NRTableau:
        return self._tableau

    # ── Private helpers ───────────────────────────────────────────────

    def _read_from_phases(self, phases):
        """Extract totals, strong ions, and volumes from a phases dict.

        Returns
        -------
        (totals, strong_ions, V_liq_L, V_gas_L)
            ``V_gas_L`` is ``None`` when no ``"gas"`` entry is present in
            *phases* (fine unless the tableau folds a gas-liquid row, in
            which case ``solve_nr`` raises a clear error downstream).

        For a component with folded gas-liquid secondaries
        (``comp.gas_species_ids``, CP1/CP2 of ``LAYER1_GAP_CLOSURE``), the
        total now sums moles across *both* phases before dividing by
        ``V_liq_L`` — matching the existing convention already used by
        ``PartitionModel.equilibrium_a_moles`` (``n_total`` spans both
        phases there too), not the liquid-only sum this used before CP2.
        """
        liq = phases.get("liquid") if hasattr(phases, "get") else None
        if liq is None:
            return {}, {}, None, None

        V_L = float(getattr(liq, "V_L", 1.0))
        if V_L <= 0.0:
            return {}, {}, None, None

        n_mol = getattr(liq, "n_mol", {})

        gas = phases.get("gas") if hasattr(phases, "get") else None
        gas_n_mol = getattr(gas, "n_mol", {}) if gas is not None else {}
        V_gas_L = float(getattr(gas, "V_L", 0.0)) if gas is not None else None

        totals: Dict[str, float] = {}
        for comp in self._tableau.components:
            # Sum n_mol for all liquid-phase species in the component
            CT_mol = sum(float(n_mol.get(sp_id, 0.0))
                         for sp_id in comp.species_ids)
            if comp.gas_species_ids:
                CT_mol += sum(float(gas_n_mol.get(sp_id, 0.0))
                              for sp_id in comp.gas_species_ids)
            if CT_mol <= 0.0:
                # Fallback: the master might be the only species in n_mol
                # (first solve, before derived species have been written back)
                CT_mol = float(n_mol.get(comp.master_id, 0.0))
            totals[comp.master_id] = CT_mol / V_L

        # Strong ions
        strong_ions: Dict[str, float] = {}
        for sp_id, ct_key in _STRONG_ION_SPECIES_TO_KEY.items():
            if sp_id in n_mol:
                val = float(n_mol[sp_id]) / V_L
                if val > 0.0:
                    strong_ions[ct_key] = val

        return totals, strong_ions, V_L, V_gas_L

    # ── Precipitation active-set loop ────────────────────────────────

    def _solve_with_precipitation(
        self,
        totals: Dict[str, float],
        strong_ions: Dict[str, float],
        T_K: float,
        am,
        *,
        tol_prec: float = 1e-8,
        max_outer_prec: int = 30,
        max_inner_prec: int = 20,
        fd_eps: float = 1e-7,
        V_liq_L: Optional[float] = None,
        V_gas_L: Optional[float] = None,
    ) -> Dict[str, Any]:
        """Outer active-set loop that enforces IAP = Ksp for each active mineral.

        Algorithm (MINTEQ/PHREEQC-style active-set Newton):

        1. Start with no active minerals (ξ = 0 for all).
        2. Outer loop: compute effective totals/strong-ions adjusted by ξ,
           run inner solve_nr(), compute IAP and SI for every mineral.
        3. Inner ξ-convergence loop (fixed active set): finite-difference
           Newton on ξ to drive f_j(ξ) = log10(IAP_j) − log10Ksp_j → 0
           for each active mineral j.  Clamp ξ_j ≥ 0.
        4. Active-set update: remove minerals with ξ_j < 0 (set to 0);
           add most-supersaturated inactive mineral if SI > 0.
        5. Repeat until active set is unchanged.

        .. note::
            Ca²⁺ (and similar precipitating ions with no dissolved
            complexation) stays a strong ion.  The outer loop adjusts its
            effective contribution via ``strong_ions["CT_Ca"] -= ξ_calcite``.
            If Ca²⁺ complexation is declared in future, it enters the NR
            tableau naturally via the BFS step; no new apparatus is required.
        """
        rxns = self._precipitation_reactions

        # Build a lookup: species_id → component master_id (for NR secondary species)
        _secondary_to_master: Dict[str, str] = {}
        for comp in self._tableau.components:
            for sp_id in comp.species_ids:
                _secondary_to_master[sp_id] = comp.master_id
        # Masters map to themselves
        for m_id in self._tableau.masters:
            _secondary_to_master.setdefault(m_id, m_id)

        # Build per-reaction metadata.
        # Each entry in dissolved_info[j] is a tuple:
        #   (species_id, charge, coeff, kind, iap_key, total_adj_key)
        # where:
        #   species_id  — id of the dissolved species in the reaction
        #   charge      — ionic charge (for activity coefficient)
        #   coeff       — stoichiometric coefficient (positive = product of dissolution)
        #   kind        — "strong" or "nr"
        #   iap_key     — key for IAP: CT_* for strong ions (actual dissolved conc),
        #                 species_id for NR species (from solve_nr output)
        #   total_adj_key — key to decrement when precipitation occurs:
        #                   CT_* for strong ions; master_id for NR species
        labels = []
        log_Ksps = []
        dissolved_info: List[List] = []

        for rxn in rxns:
            label = rxn.label or rxn.stoichiometry[0].species.id
            labels.append(label)
            log_ksp = _vant_hoff_log_K(
                float(rxn.log_K),
                rxn.dH_J_per_mol,
                T_K,
                float(rxn.T_ref_K),
            )
            log_Ksps.append(log_ksp)

            d_info = []
            for e in rxn.stoichiometry:
                if e.phase == "solid":
                    continue
                sp_id = e.species.id
                coeff = float(e.coefficient)
                charge = int(e.species.charge)
                strong_key = _STRONG_ION_SPECIES_TO_KEY.get(sp_id)
                if strong_key is not None:
                    d_info.append((sp_id, charge, coeff, "strong", strong_key, strong_key))
                else:
                    # NR species: adjust the MASTER component total
                    master_id = _secondary_to_master.get(sp_id, sp_id)
                    d_info.append((sp_id, charge, coeff, "nr", sp_id, master_id))
            dissolved_info.append(d_info)

        n_rxns = len(rxns)
        xi = np.zeros(n_rxns)   # mol/L precipitated for each mineral
        active = [False] * n_rxns

        def _effective_inputs(xi_arr):
            """Compute totals and strong_ions adjusted by current xi.

            When xi_j > 0 (precipitation), each dissolved product of mineral j
            has its component total reduced by coeff * xi_j.  For strong ions
            (e.g. Ca²⁺), the CT_* key is decremented.  For NR species (e.g.
            CO₃²⁻), the master component total (e.g. CT_CO₂) is decremented.
            """
            eff_totals = dict(totals)
            eff_strong = dict(strong_ions)
            for j in range(n_rxns):
                xj = xi_arr[j]
                if xj == 0.0:
                    continue
                for (sp_id, charge, coeff, kind, iap_key, total_adj_key) in dissolved_info[j]:
                    if kind == "strong":
                        eff_strong[total_adj_key] = max(
                            eff_strong.get(total_adj_key, 0.0) - coeff * xj, 0.0
                        )
                    else:
                        eff_totals[total_adj_key] = max(
                            eff_totals.get(total_adj_key, 0.0) - coeff * xj, 0.0
                        )
            return eff_totals, eff_strong

        def _run_inner(xi_arr, *, _retain_jac: bool = False):
            """Run solve_nr with adjusted inputs; return output dict."""
            et, es = _effective_inputs(xi_arr)
            return solve_nr(
                self._tableau,
                et,
                es,
                T_K=T_K,
                activity_model=am,
                cache=None,
                max_log_activity=self.max_log_activity,
                min_component_total=self.min_component_total,
                retain_jacobian=_retain_jac,
                V_liq_L=V_liq_L,
                V_gas_L=V_gas_L,
            )

        def _compute_iap_and_si(out, xi_arr, j):
            """Compute log10(IAP) and SI for mineral j given inner solve output.

            For strong ions: use effective dissolved concentration (CT_* − xi).
            For NR species: use the concentration from the inner solve output
            (solve_nr already accounts for the adjusted component total).
            """
            I = float(out.get("IonicStrength", 0.0))
            log_iap = 0.0
            et, es = _effective_inputs(xi_arr)
            for (sp_id, charge, coeff, kind, iap_key, total_adj_key) in dissolved_info[j]:
                if kind == "strong":
                    c_sp = es.get(iap_key, 0.0)
                else:
                    c_sp = float(out.get(iap_key, 0.0))
                g = _gamma_safe(am, abs(charge), I, T_K)
                log_a = np.log10(max(g * c_sp, 1e-300))
                log_iap += coeff * log_a
            si = log_iap - log_Ksps[j]
            return log_iap, si

        # ── Outer active-set loop ────────────────────────────────────────
        prev_active = None
        out = _run_inner(xi)

        for _outer in range(max_outer_prec):
            out = _run_inner(xi)

            # ── Inner ξ-convergence loop (fixed active set) ──────────
            active_indices = [j for j, a in enumerate(active) if a]
            if active_indices:
                for _inner in range(max_inner_prec):
                    # Compute residuals f_j = log10(IAP_j) - log10Ksp_j
                    f = np.zeros(len(active_indices))
                    for idx_i, j in enumerate(active_indices):
                        _, si = _compute_iap_and_si(out, xi, j)
                        f[idx_i] = si   # = log10(IAP) - log10(Ksp)

                    if np.max(np.abs(f)) < tol_prec:
                        break

                    # Finite-difference Jacobian: df_j / d(xi_j) for active minerals
                    # (we assume off-diagonal terms are small — diagonal FD only)
                    J_fd = np.zeros(len(active_indices))
                    for idx_i, j in enumerate(active_indices):
                        xi_p = xi.copy()
                        xi_p[j] += fd_eps
                        out_p = _run_inner(xi_p)
                        _, si_p = _compute_iap_and_si(out_p, xi_p, j)
                        J_fd[idx_i] = (si_p - f[idx_i]) / fd_eps

                    # Newton step with physical step-size constraint.
                    # The step Δxi must not exceed the minimum available
                    # component total (can't precipitate more than exists in
                    # solution).  Use a 0.9 safety factor.
                    et_cur, es_cur = _effective_inputs(xi)
                    for idx_i, j in enumerate(active_indices):
                        if abs(J_fd[idx_i]) < 1e-30:
                            continue
                        delta = -f[idx_i] / J_fd[idx_i]
                        if delta > 0:
                            # Compute maximum feasible step for mineral j
                            xi_max_step = float("inf")
                            for (sp_id, charge, coeff, kind, iap_key, total_adj_key) \
                                    in dissolved_info[j]:
                                if coeff <= 0:
                                    continue
                                if kind == "strong":
                                    avail = es_cur.get(total_adj_key, 0.0)
                                else:
                                    avail = et_cur.get(total_adj_key, 0.0)
                                if avail <= 0:
                                    xi_max_step = 0.0
                                    break
                                xi_max_step = min(xi_max_step, avail / coeff * 0.9)
                            delta = min(delta, xi_max_step)
                        xi[j] = max(xi[j] + delta, 0.0)

                    out = _run_inner(xi)

            # ── Active-set update ────────────────────────────────────
            # 1. Remove minerals with xi < 0 (already clamped to 0 above;
            #    but if xi hit exactly 0 during Newton, deactivate)
            for j in active_indices:
                if xi[j] <= 0.0:
                    active[j] = False
                    xi[j] = 0.0

            # 2. Compute SI for all inactive minerals; add the most
            #    supersaturated one if SI > 0
            best_j, best_si = None, 0.0
            for j in range(n_rxns):
                if active[j]:
                    continue
                _, si = _compute_iap_and_si(out, xi, j)
                if si > best_si:
                    best_si = si
                    best_j = j
            if best_j is not None:
                active[best_j] = True

            # 3. Check convergence of active set
            current_active = tuple(active)
            if current_active == prev_active:
                break
            prev_active = current_active

        # ── Final solve with converged xi ────────────────────────────
        out = _run_inner(xi, _retain_jac=self.retain_jacobian)

        # ── Build "minerals" output dict ──────────────────────────────
        minerals: Dict[str, Dict[str, float]] = {}
        for j, label in enumerate(labels):
            _, si = _compute_iap_and_si(out, xi, j)
            minerals[label] = {"xi_mol_L": float(xi[j]), "SI": float(si)}

        out["minerals"] = minerals
        return out

    # ── White-box methods (WhiteBoxEngineProtocol) ────────────────────

    def jacobian_dg_dz(self, **kwargs) -> SparseJacobian:
        """Converged NR Jacobian ∂g/∂z (in log-activity space).

        Shape: ``(m, m)`` where ``m = len(tableau.masters)``.
        Rows are residual equations (mass-balance rows first, charge balance
        last); columns are master log-activities in ``tableau.masters`` order.

        Requires ``retain_jacobian=True`` at construction and a prior
        ``solve()`` call.

        Raises
        ------
        RuntimeError
            If ``retain_jacobian=False`` or ``solve()`` has not been called.
        """
        if not self.retain_jacobian:
            raise RuntimeError(_JACOBIAN_DISABLED_MSG)
        if self._cached_jacobian is None:
            raise RuntimeError(
                "No Jacobian cached. Call solve() first."
            )
        m = len(self._tableau.masters)
        row_ids = tuple(
            comp.master_id for comp in self._tableau.components
        ) + ("H+",)
        col_ids = tuple(self._tableau.masters)
        return SparseJacobian(
            matrix=csr_matrix(self._cached_jacobian),
            row_ids=row_ids,
            col_ids=col_ids,
        )

    def jacobian_dg_dy(self, **kwargs) -> SparseJacobian:
        """Structural Jacobian ∂g/∂y (residuals vs. component totals).

        Shape: ``(m, m-1)`` where the top ``(m-1, m-1)`` block is ``-I``
        (each mass-balance row depends on its own component total with
        coefficient ``-1``) and the bottom row is zero (charge balance is
        independent of component totals).

        This matrix is structural and never changes after construction.

        Requires ``retain_jacobian=True`` at construction.
        """
        if not self.retain_jacobian:
            raise RuntimeError(_JACOBIAN_DISABLED_MSG)
        m = len(self._tableau.masters)
        n_comp = m - 1
        # Build as lil for efficient entry, convert to csr for output
        from scipy.sparse import lil_matrix
        dg_dy = lil_matrix((m, n_comp))
        for i in range(n_comp):
            dg_dy[i, i] = -1.0
        row_ids = tuple(
            comp.master_id for comp in self._tableau.components
        ) + ("H+",)
        col_ids = tuple(comp.master_id for comp in self._tableau.components)
        return SparseJacobian(
            matrix=dg_dy.tocsr(),
            row_ids=row_ids,
            col_ids=col_ids,
        )

    def residual(self, **kwargs) -> np.ndarray:
        """Evaluate the NR residual ``g(y, z)`` without updating internal state.

        The returned vector has ``m = len(tableau.masters)`` entries:
        mass-balance rows first (one per non-H⁺ component), charge-balance
        row last.

        Parameters
        ----------
        totals : dict
            ``{master_id: C_total_mol_L}`` — component totals ``y``.
        log_activities : array-like, optional
            Master log10-activities ``z = x``.  If not supplied, uses the
            cached value from the most recent ``solve()`` call.
        strong_ions : dict, optional
            ``{CT_Na: mol_L, ...}`` — strong-ion concentrations.
        T_K : float, optional
            Temperature (K); defaults to engine's ``T_C``.

        Raises
        ------
        RuntimeError
            If ``retain_jacobian=False`` or no cached state is available.
        ValueError
            If ``totals=`` is not provided.
        """
        if not self.retain_jacobian:
            raise RuntimeError(_JACOBIAN_DISABLED_MSG)

        totals = kwargs.pop("totals", None)
        if totals is None:
            raise ValueError("residual() requires totals= kwarg.")
        log_activities = kwargs.pop("log_activities", None)
        strong_ions = kwargs.pop("strong_ions", {})
        T_K_override = kwargs.pop("T_K", None)
        T_K = float(T_K_override) if T_K_override is not None else 273.15 + self.T_C

        if log_activities is not None:
            x = np.asarray(log_activities, dtype=float)
        elif self._cached_x is not None:
            x = self._cached_x.copy()
        else:
            raise RuntimeError(
                "No cached log-activities. Call solve() first, or supply log_activities=."
            )

        I = self._cached_I if self._cached_I is not None else 0.01

        strong_charge = sum(
            _STRONG_CHARGES.get(key, 0) * float(val)
            for key, val in strong_ions.items()
        )

        all_charges: Dict[str, int] = dict(self._tableau.master_charges)
        for sec in self._tableau.secondaries:
            all_charges[sec.species_id] = sec.charge

        am = (
            self._liquid_activity
            if self._liquid_activity is not None
            else make_activity_model(self.use_activity, self.activity_model)
        )
        gammas = _build_gammas(self._tableau, all_charges, I, T_K, am)
        c = _compute_concentrations(
            self._tableau, x, gammas, max_log_activity=self.max_log_activity
        )

        pin_val = -float(self.max_log_activity)
        pin_specs = []
        for row_idx, comp in enumerate(self._tableau.components):
            CT = float(totals.get(comp.master_id, 0.0))
            x_col = self._tableau.masters.index(comp.master_id)
            if CT < float(self.min_component_total):
                pin_specs.append((row_idx, x_col, pin_val))

        R, _ = _residual_and_jacobian(
            self._tableau, x, c, totals, strong_charge,
            pin_specs=pin_specs if pin_specs else None,
        )
        return R

    def jacobian_dz_dy(self, **kwargs) -> SpeciationJacobian:
        """Total sensitivity ∂c/∂y via the implicit function theorem.

        Computes ``∂c_j / ∂T_i`` for all algebraic species ``j`` and all
        component totals ``i`` analytically from the retained NR Jacobian and
        the chain rule through the log-linear formula:

        ``∂c_j/∂T_i = Σ_k (ln10 · ν_{jk} · c_j) · (∂x_k/∂T_i)``

        where ``∂x/∂T = -(∂g/∂z)^{-1} · (∂g/∂y)`` by the implicit function
        theorem.

        Returns a :class:`~PyOMES.chemical_equilibrium.protocols.SpeciationJacobian`
        with rows in sorted-alphabetical algebraic-species order (matching
        :class:`~PyOMES.chemical_equilibrium.numerical_gradient.NumericalGradientEquilibriumEngine`)
        and columns in ``tableau.components`` order.

        Raises
        ------
        RuntimeError
            If ``retain_jacobian=False``, or ``solve()`` has not been called,
            or the converged Jacobian is singular.
        """
        if not self.retain_jacobian:
            raise RuntimeError(_JACOBIAN_DISABLED_MSG)
        if self._cached_jacobian is None or self._cached_concentrations is None:
            raise RuntimeError(
                "No cached state. Call solve() first."
            )

        m = len(self._tableau.masters)
        n_comp = m - 1

        # ∂g/∂y: (m, n_comp) dense; -I in top block, 0 in charge-balance row
        dg_dy = np.zeros((m, n_comp))
        for i in range(n_comp):
            dg_dy[i, i] = -1.0

        # ∂x/∂T = -(J_NR)^{-1} * dg_dy  →  shape (m, n_comp)
        J = self._cached_jacobian
        try:
            dx_dy = np.linalg.solve(J, -dg_dy)
        except np.linalg.LinAlgError:
            raise RuntimeError(
                "NR Jacobian is singular; cannot compute jacobian_dz_dy()."
            )

        component_ids = tuple(comp.master_id for comp in self._tableau.components)
        alg_ids = tuple(sorted(self.algebraic_species()))
        n_alg = len(alg_ids)
        alg_id_to_row = {sp: i for i, sp in enumerate(alg_ids)}
        master_to_xcol = {m_id: k for k, m_id in enumerate(self._tableau.masters)}

        dz_dy = np.zeros((n_alg, n_comp))

        # Masters: ∂c_k/∂T_i = ln10 * c_k * dx_dy[k_col, i]
        for k, m_id in enumerate(self._tableau.masters):
            row = alg_id_to_row.get(m_id)
            if row is None:
                continue
            c_k = self._cached_concentrations.get(m_id, 0.0)
            dz_dy[row, :] = _LN10 * c_k * dx_dy[k, :]

        # Secondaries: ∂c_j/∂T_i = Σ_k ln10 * ν_jk * c_j * dx_dy[k, i]
        for sec in self._tableau.secondaries:
            row = alg_id_to_row.get(sec.species_id)
            if row is None:
                continue
            c_j = self._cached_concentrations.get(sec.species_id, 0.0)
            for nu_k_id, nu_jk in sec.nu.items():
                k_col = master_to_xcol.get(nu_k_id)
                if k_col is None:
                    continue
                dz_dy[row, :] += nu_jk * _LN10 * c_j * dx_dy[k_col, :]

        # H2O is in algebraic_species() but has no log-activity variable;
        # its sensitivity is negligible (dilute solution) — row stays at 0.

        return SpeciationJacobian(
            dz_dy=dz_dy,
            algebraic_ids=alg_ids,
            component_ids=component_ids,
        )
