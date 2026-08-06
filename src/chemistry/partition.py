# -*- coding: utf-8 -*-
"""PartitionModel protocol and the Henry/Raoult/Ksp *Equilibrium implementations.

A PartitionModel describes how a species distributes between two phases at
equilibrium.  HenryEquilibrium covers linear (Henry's law) gas-liquid
systems.  Future non-linear models (Langmuir, Freundlich) implement the
same protocol.

Capacity parameters are phase-agnostic floats so the protocol is usable for
gas-liquid (volumes in L), solid-liquid (solid mass in kg), etc.

The partition_ratio() return value selects the solver path in KineticGasLiquidLink:
  - float  → analytical exponential step solution (exact for linear models)
  - None   → ODE instantaneous-rate path (required for non-linear models)

``HenryEquilibrium``, ``RaoultEquilibrium``, and single-ion ``KspEquilibrium``
each also satisfy :class:`~PyOMES.reactions.equilibrium.EquilibriumConstraint`
(``stoichiometry``/``log_K``/``dH_J_per_mol``/``T_ref_K``) — the same
constructed instance can be passed both to ``KineticGasLiquidLink``'s
``partition_models=`` dict and into the reaction list feeding
``NRChemicalEquilibriumEngine``/``ChemicalEquilibriumEngine``, per
``EQUILIBRIUM_CONSTRAINT_UNIFICATION`` CP1.

Activity correction in HenryEquilibrium
----------------------------------------
Henry's law relates gas-phase partial pressure to liquid-phase fugacity:
    p_i = f_i^liq = γ_i × C_i / kH

so the effective Henry constant becomes kH / γ_i, reducing the corrected
partition ratio by a factor of γ_i.  A Davies activity correction (γ < 1
for ionic species) reduces the effective partition toward liquid, keeping
more solute dissolved.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import (
    TYPE_CHECKING, Any, Dict, Optional, Protocol, Sequence, Tuple, Union,
    runtime_checkable,
)

from ..units import R_J_PER_MOL_K as _R_J_MOL
from .species import Species

if TYPE_CHECKING:
    from PyOMES.thermo import ThermoFramework
    from PyOMES.reactions.stoichiometry import StoichiometryEntry

_R_L_ATM_MOL_K = 0.0820574  # L·atm/(mol·K) — must match core.phases.R_L_ATM_MOL_K


def _resolve_species(species: Union[str, Species, None]) -> Optional[Species]:
    """Resolve a species id string or ``Species`` object to a ``Species``.

    Deferred import of ``PyOMES.reactions.stoichiometry`` avoids a
    chemistry↔reactions circular import: this module is imported during
    ``PyOMES.chemistry`` package initialization, and
    ``PyOMES.reactions.stoichiometry``/``PyOMES.reactions.equilibrium``
    import back from ``PyOMES.chemistry.species``.
    """
    if species is None or isinstance(species, Species):
        return species
    from ..reactions.stoichiometry import _get_common_species
    lookup = _get_common_species()
    if species not in lookup:
        raise ValueError(
            f"Unknown species id {species!r} — pass a Species object "
            "directly, or use an id from PyOMES.chemistry.common_species."
        )
    return lookup[species]


@runtime_checkable
class PartitionModel(Protocol):
    """Protocol for phase-partition relationships."""

    def equilibrium_a_moles(
        self,
        n_total: float,
        capacity_a: float,
        capacity_b: float,
        T_K: float,
        alpha: float = 1.0,
    ) -> float:
        """Moles in phase 'a' at equilibrium given n_total across both phases."""
        ...

    def partition_ratio(
        self,
        capacity_a: float,
        capacity_b: float,
        T_K: float,
        alpha: float = 1.0,
    ) -> Optional[float]:
        """Dimensionless partition ratio for the analytical step solution.

        Returns the ratio ``n_a_eq / n_b_eq`` at equilibrium for linear
        models.  Return ``None`` for non-linear models (Langmuir, etc.)
        that require the ODE instantaneous-rate path instead.
        """
        ...


@dataclass(frozen=True)
class HenryEquilibrium:
    """Henry's law gas-liquid partition (Sander solubility convention).

    Parameters
    ----------
    H_ref : float
        Henry solubility at T_ref in mol m⁻³ Pa⁻¹ (Sander convention).
    dlnH : float
        d(ln kH)/d(1/T) in K (van 't Hoff temperature sensitivity).
        Positive for exothermic dissolution (more soluble at lower T).
    T_ref : float
        Reference temperature for H_ref (K).
    thermo : ThermoFramework or None
        When provided, applies the liquid-phase activity correction
        ``γ_i`` via ``thermo.liquid_activity``.  Callers must also
        pass ``species_id``, ``x_mol``, and ``charge`` to the
        ``partition_ratio`` / ``equilibrium_a_moles`` methods for the
        correction to take effect.  Without ``thermo``, behaviour is
        unchanged (γ_i = 1).
    gas_species, liquid_species : str or Species, optional
        The gas-phase and liquid-phase forms this constant relates
        (e.g. ``gas_species="CO2", liquid_species="CO2"``). Optional —
        only needed to also satisfy
        :class:`~PyOMES.reactions.equilibrium.EquilibriumConstraint`
        (i.e. to appear in the reaction list fed to
        ``NRChemicalEquilibriumEngine``/``ChemicalEquilibriumEngine``); the
        ``PartitionModel`` role (``partition_ratio``/
        ``equilibrium_a_moles``) does not need them. Accepts a species
        id string (resolved against ``PyOMES.chemistry.common_species``)
        or a ``Species`` object directly.
    label : str
        Human-readable label (optional, for diagnostics/logging — e.g.
        ``ReactionSystem.show_reactions()``).
    """

    H_ref: float
    dlnH: float
    T_ref: float = 298.15
    thermo: Optional[object] = field(default=None, compare=False)  # ThermoFramework | None
    gas_species: Union[str, Species, None] = None
    liquid_species: Union[str, Species, None] = None
    label: str = ""

    def partition_ratio(
        self,
        capacity_a: float,
        capacity_b: float,
        T_K: float,
        alpha: float = 1.0,
        *,
        species_id: str = "",
        x_mol: Optional[Dict[str, float]] = None,
        charge: Optional[Dict[str, int]] = None,
    ) -> float:
        """Return (kH / (alpha × γ_i)) × R × T × V_liq / V_gas.

        ratio >> 1 → mostly in liquid; ratio << 1 → mostly in gas.

        Parameters
        ----------
        capacity_a, capacity_b : float
            Phase capacities (volumes in L for gas-liquid).
        T_K : float
            Temperature (K).
        alpha : float
            Ionisation fraction (gas-transferable fraction).
        species_id : str
            Species key used to look up γ_i from the activity model result.
            Ignored when ``thermo`` is None.
        x_mol : dict, optional
            Full liquid-phase composition in mol/L (used by activity model).
        charge : dict, optional
            Charge per species (used by activity model).
        """
        kH = self._kH_mol_L_atm(T_K)
        gamma = self._gamma(species_id, x_mol, charge, T_K)
        return (kH / (max(1e-12, alpha) * gamma)) * _R_L_ATM_MOL_K * T_K * capacity_a / capacity_b

    def equilibrium_a_moles(
        self,
        n_total: float,
        capacity_a: float,
        capacity_b: float,
        T_K: float,
        alpha: float = 1.0,
        *,
        species_id: str = "",
        x_mol: Optional[Dict[str, float]] = None,
        charge: Optional[Dict[str, int]] = None,
    ) -> float:
        """Moles in the liquid (phase 'a') at gas-liquid equilibrium."""
        r = self.partition_ratio(
            capacity_a, capacity_b, T_K, alpha,
            species_id=species_id, x_mol=x_mol, charge=charge,
        )
        return r * n_total / (1.0 + r)

    def _kH_mol_L_atm(self, T_K: float) -> float:
        H_ref_mol_L_atm = (self.H_ref / 1000.0) * 101325.0
        return H_ref_mol_L_atm * math.exp(self.dlnH * (1.0 / T_K - 1.0 / self.T_ref))

    def _gamma(
        self,
        species_id: str,
        x_mol: Optional[Dict[str, float]],
        charge: Optional[Dict[str, int]],
        T_K: float,
    ) -> float:
        """Return γ_i from thermo.liquid_activity, or 1.0 if not available."""
        if self.thermo is None or not x_mol or not charge:
            return 1.0
        gammas = self.thermo.liquid_activity.gamma_all(x_mol, T_K, charge=charge)
        return gammas.get(species_id, 1.0)

    # ── EquilibriumConstraint conformance ───────────────────────────────────

    @property
    def T_ref_K(self) -> float:
        return self.T_ref

    @property
    def log_K(self) -> float:
        """log10(kH) at T_ref — reference mass-action constant for gas ⇌ liquid."""
        return math.log10(self._kH_mol_L_atm(self.T_ref))

    @property
    def dH_J_per_mol(self) -> float:
        """Van 't Hoff enthalpy derived from dlnH (§14.1: dH = -dlnH × R)."""
        return -self.dlnH * _R_J_MOL

    @property
    def stoichiometry(self) -> Tuple["StoichiometryEntry", ...]:
        """gas ⇌ liquid stoichiometry, when gas_species/liquid_species are set.

        Returns an empty tuple when either is unset (e.g. an instance
        constructed for ``PartitionModel``-only use) rather than
        raising, so ``isinstance(x, EquilibriumConstraint)`` stays safe.
        """
        if self.gas_species is None or self.liquid_species is None:
            return ()
        from ..reactions.stoichiometry import StoichiometryEntry
        gas = _resolve_species(self.gas_species)
        liq = _resolve_species(self.liquid_species)
        return (
            StoichiometryEntry(species=gas, phase="gas", coefficient=-1.0),
            StoichiometryEntry(species=liq, phase="liquid", coefficient=+1.0),
        )


def HenryPartition(*args: Any, **kwargs: Any) -> HenryEquilibrium:
    """Deprecated alias for :class:`HenryEquilibrium`.

    Kept for one phase (``EQUILIBRIUM_CONSTRAINT_UNIFICATION`` CP1) so
    existing call sites/tests continue to work while call sites
    migrate. Emits ``DeprecationWarning`` and constructs-and-returns a
    real ``HenryEquilibrium`` (not a subclass), so ``type(x) is
    HenryEquilibrium`` and equality both hold.
    """
    warnings.warn(
        "HenryPartition is deprecated; use HenryEquilibrium instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return HenryEquilibrium(*args, **kwargs)


# ── Raoult's law partition for solvent (H₂O) ───────────────────────────────

_M_WATER = 18.015   # g/mol
_P_SAT_REF = 0.03169  # atm — saturation pressure of pure water at 298.15 K
_T_REF_WATER = 298.15  # K


@dataclass(frozen=True)
class RaoultEquilibrium:
    """Raoult's law vapour pressure partition for the solvent (H₂O).

    Uses Clausius-Clapeyron for temperature correction:
        P_sat(T) = P_sat_ref × exp(−ΔH_vap/R × (1/T − 1/T_ref))

    The partial pressure of water vapour in equilibrium with a liquid
    of mole-fraction x_w is:
        p_w = x_w × P_sat(T)

    where ``n_total`` mol spans liquid (phase 'a') and gas (phase 'b'),
    and the returned ratio n_liq/n_gas is derived from the condition
    p_w = x_w × P_sat(T) at equilibrium (see ``partition_ratio``).

    Satisfies ``PartitionModel``:
    - ``partition_ratio`` returns a float (analytical exponential step)
    - ``equilibrium_a_moles`` solves for n_liq at equilibrium

    Also satisfies
    :class:`~PyOMES.reactions.equilibrium.EquilibriumConstraint` via the
    ``gas ⇌ liquid`` stoichiometry named by ``gas_species``/
    ``liquid_species`` (default ``"H2O"``/``"H2O"``).

    Parameters
    ----------
    P_sat_ref : float
        Saturation pressure of pure water at T_ref (atm).
        Default 0.03169 atm (25 °C).
    dH_vap : float
        Enthalpy of vaporisation (J/mol).  Default 44 011 J/mol (H₂O at 25 °C).
    T_ref : float
        Reference temperature for P_sat_ref (K).  Default 298.15 K.
    C_water_mol_L : float
        Molar concentration of pure water (mol/L).  Default 55.51 mol/L.
        Used to convert liquid volume (capacity_a in L) to moles of water
        when computing the mole fraction x_w.
    gas_species, liquid_species : str or Species, optional
        The gas-phase and liquid-phase forms this constant relates.
        Default ``"H2O"``/``"H2O"`` (resolved against
        ``PyOMES.chemistry.common_species``); accepts a ``Species``
        object directly instead.
    label : str
        Human-readable label (optional, for diagnostics/logging — e.g.
        ``ReactionSystem.show_reactions()``).
    """

    P_sat_ref: float = _P_SAT_REF
    dH_vap: float = 44011.0   # J/mol at 25 °C
    T_ref: float = _T_REF_WATER
    C_water_mol_L: float = 55.51  # mol/L pure water at 25 °C
    gas_species: Union[str, Species, None] = "H2O"
    liquid_species: Union[str, Species, None] = "H2O"
    label: str = ""

    def P_sat(self, T_K: float) -> float:
        """Saturation pressure of pure water (atm) at T_K via Clausius-Clapeyron."""
        return self.P_sat_ref * math.exp(
            -self.dH_vap / _R_J_MOL * (1.0 / T_K - 1.0 / self.T_ref)
        )

    def partition_ratio(
        self,
        capacity_a: float,
        capacity_b: float,
        T_K: float,
        alpha: float = 1.0,
    ) -> float:
        """Dimensionless ratio n_liq/n_gas for H₂O at equilibrium.

        Derived from p_w = x_w × P_sat(T) with x_w ≈ n_liq_water / (C_w × V_liq),
        and the ideal-gas condition n_gas = p_w × V_gas / (R T).

        For a single-species system where n_total = n_liq + n_gas:
            n_liq/n_gas = C_w × V_liq × R T / (P_sat × V_gas)

        ``alpha`` is accepted for protocol compatibility but is meaningless
        for Raoult's law (water is not ionised); it is silently ignored.
        """
        Ps = self.P_sat(T_K)
        if Ps <= 0.0:
            return float("inf")
        return (self.C_water_mol_L * capacity_a * _R_L_ATM_MOL_K * T_K) / (Ps * capacity_b)

    def equilibrium_a_moles(
        self,
        n_total: float,
        capacity_a: float,
        capacity_b: float,
        T_K: float,
        alpha: float = 1.0,
    ) -> float:
        """Moles of H₂O in the liquid phase at vapour-liquid equilibrium."""
        r = self.partition_ratio(capacity_a, capacity_b, T_K, alpha)
        return r * n_total / (1.0 + r)

    # ── EquilibriumConstraint conformance ───────────────────────────────────

    @property
    def T_ref_K(self) -> float:
        return self.T_ref

    @property
    def log_K(self) -> float:
        """log10(1/P_sat_ref) — reference mass-action constant for gas ⇌ liquid."""
        return -math.log10(self.P_sat_ref)

    @property
    def dH_J_per_mol(self) -> float:
        """Van 't Hoff enthalpy: dH = -dH_vap (Clausius-Clapeyron sign flip, §14.1)."""
        return -self.dH_vap

    @property
    def stoichiometry(self) -> Tuple["StoichiometryEntry", ...]:
        """gas ⇌ liquid stoichiometry, when gas_species/liquid_species are set.

        Returns an empty tuple when either is unset rather than
        raising, so ``isinstance(x, EquilibriumConstraint)`` stays safe.
        """
        if self.gas_species is None or self.liquid_species is None:
            return ()
        from ..reactions.stoichiometry import StoichiometryEntry
        gas = _resolve_species(self.gas_species)
        liq = _resolve_species(self.liquid_species)
        return (
            StoichiometryEntry(species=gas, phase="gas", coefficient=-1.0),
            StoichiometryEntry(species=liq, phase="liquid", coefficient=+1.0),
        )


def RaoultPartition(*args: Any, **kwargs: Any) -> RaoultEquilibrium:
    """Deprecated alias for :class:`RaoultEquilibrium`. See :func:`HenryPartition`."""
    warnings.warn(
        "RaoultPartition is deprecated; use RaoultEquilibrium instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return RaoultEquilibrium(*args, **kwargs)


# ── Ksp (solubility-product) solid-liquid equilibrium ───────────────────────

class KspEquilibrium:
    """Solubility-product (Ksp) solid-liquid equilibrium.

    Formalizes the shipped ``NR_PRECIPITATION_SPECIATION`` convention
    (``log_K = Ksp``) as a named constructor satisfying both
    ``PartitionModel`` and
    :class:`~PyOMES.reactions.equilibrium.EquilibriumConstraint`, instead
    of requiring hand-built
    :class:`~PyOMES.reactions.equilibrium.EquilibriumReaction` objects.

    The ``PartitionModel`` role (``partition_ratio``/
    ``equilibrium_a_moles``) is well-defined only for the single-ion
    case — exactly one solid entry and one dissolved entry in
    ``stoichiometry`` — per the wrapping litmus test (a genuinely
    multi-ion solubility equilibrium requires the active-set NR solver,
    not a closed-form partition ratio). ``partition_ratio()`` returns
    ``None`` for the single-ion case (per ``PartitionModel``'s own
    contract for non-linear/capacity-limited models — solubility is a
    concentration cap, not a fixed ratio); both methods raise
    ``NotImplementedError`` for multi-ion stoichiometry, though the
    type still satisfies ``EquilibriumConstraint`` (its reaction-list
    role is unaffected).

    Parameters
    ----------
    stoichiometry : list of StoichiometryEntry or str
        Solid ⇌ dissolved-ion(s) stoichiometry, either as
        ``StoichiometryEntry`` objects or a string (see
        :class:`~PyOMES.reactions.equilibrium.EquilibriumReaction` for
        the string format).
    Ksp : float
        Solubility product (mass-action units matching the
        stoichiometry's dissolved-species exponents).
    species : dict[str, Species], optional
        Caller-supplied species for locally declared minerals not in
        ``common_species``. Only used when *stoichiometry* is a string.
    dH_J_per_mol : float, optional
        Van 't Hoff reaction enthalpy (J/mol) at ``T_ref_K``.
    T_ref_K : float
        Reference temperature for ``Ksp``/``dH_J_per_mol`` (default
        ``298.15``).
    label : str
        Human-readable label (optional, for diagnostics/logging).
    """

    def __init__(
        self,
        stoichiometry: Union[Sequence["StoichiometryEntry"], str],
        *,
        Ksp: float,
        species: Optional[Dict[str, Species]] = None,
        dH_J_per_mol: Optional[float] = None,
        T_ref_K: float = 298.15,
        label: str = "",
    ):
        if isinstance(stoichiometry, str):
            from ..reactions.stoichiometry import _parse_stoichiometry
            stoichiometry = _parse_stoichiometry(
                stoichiometry, species, reaction_type="equilibrium"
            )
        self.stoichiometry: Tuple["StoichiometryEntry", ...] = tuple(stoichiometry)
        self.Ksp: float = float(Ksp)
        self.dH_J_per_mol: Optional[float] = (
            float(dH_J_per_mol) if dH_J_per_mol is not None else None
        )
        self.T_ref_K: float = float(T_ref_K)
        self.label: str = str(label)

    @property
    def log_K(self) -> float:
        return math.log10(self.Ksp)

    def _is_single_ion(self) -> bool:
        solid = [e for e in self.stoichiometry if e.phase == "solid"]
        dissolved = [e for e in self.stoichiometry if e.phase != "solid"]
        return len(solid) == 1 and len(dissolved) == 1

    def partition_ratio(
        self,
        capacity_a: float,
        capacity_b: float,
        T_K: float,
        alpha: float = 1.0,
    ) -> Optional[float]:
        """``None`` for single-ion (solubility is a cap, not a ratio).

        Raises ``NotImplementedError`` for multi-ion stoichiometry.
        """
        if not self._is_single_ion():
            raise NotImplementedError(
                "KspEquilibrium.partition_ratio is defined only for "
                "single-ion solubility equilibria (one solid + one "
                "dissolved species); multi-ion Ksp requires the "
                "active-set NR solver."
            )
        return None

    def equilibrium_a_moles(
        self,
        n_total: float,
        capacity_a: float,
        capacity_b: float,
        T_K: float,
        alpha: float = 1.0,
    ) -> float:
        """``min(n_total, Ksp(T) × capacity_a)`` — single-ion solubility cap.

        Raises ``NotImplementedError`` for multi-ion stoichiometry.
        """
        if not self._is_single_ion():
            raise NotImplementedError(
                "KspEquilibrium.equilibrium_a_moles is defined only for "
                "single-ion solubility equilibria (one solid + one "
                "dissolved species); multi-ion Ksp requires the "
                "active-set NR solver."
            )
        from ..reactions.equilibrium import vant_hoff_log_K
        Ksp_T = 10.0 ** vant_hoff_log_K(self, T_K)
        return min(n_total, Ksp_T * capacity_a)

    def __repr__(self) -> str:
        lbl = f" [{self.label}]" if self.label else ""
        return f"KspEquilibrium(Ksp={self.Ksp:.3g}{lbl})"


# ── Multispecies partition model protocol and VLE implementation ────────────

@runtime_checkable
class MultispeciesPartitionModel(Protocol):
    """Phase-partition protocol for models where equilibrium of one species
    depends on the concentrations of all other species.

    Use cases:
    - Coupled VLE with non-ideal gas EOS (Peng-Robinson mixing rules)
    - Competitive adsorption (Langmuir multicomponent)
    - NRTL liquid–liquid partitioning

    Unlike ``PartitionModel`` (per-species, independent), this protocol
    solves all species simultaneously in a single call.
    """

    def equilibrium_all_a_moles(
        self,
        n_total_all: Dict[str, float],
        capacity_a: float,
        capacity_b: float,
        T_K: float,
    ) -> Dict[str, float]:
        """Moles of each species in phase 'a' at simultaneous equilibrium.

        Parameters
        ----------
        n_total_all : dict
            ``{species_id: n_total_i}`` — total moles of each species
            summed across both phases.
        capacity_a : float
            Phase 'a' capacity (volume in L for gas-liquid; mass in kg
            for solid-liquid).
        capacity_b : float
            Phase 'b' capacity.
        T_K : float
            Temperature (K).

        Returns
        -------
        dict
            ``{species_id: n_a_i}`` — moles in phase 'a' at equilibrium.
            Species absent from the return dict are treated as n_a = 0.
        """
        ...


def _kH_mol_L_atm_from_ref(H_ref: float, dlnH: float, T_K: float, T_ref: float) -> float:
    """Convert Sander kH (mol m⁻³ Pa⁻¹) to mol/(L·atm) and apply van 't Hoff."""
    H_mol_L_atm = (H_ref / 1000.0) * 101325.0
    return H_mol_L_atm * math.exp(dlnH * (1.0 / T_K - 1.0 / T_ref))


@dataclass(frozen=True)
class MultispeciesVLEPartition:
    """Coupled gas-liquid VLE for multiple volatile species.

    With ``IdealGasEOS``, each species decouples and the solution is
    analytical (same formula as ``HenryEquilibrium`` without the ``alpha``
    correction).  Non-ideal EOS (Peng-Robinson) coupling is deferred —
    passing a non-ideal EOS raises ``NotImplementedError``.

    Satisfies ``MultispeciesPartitionModel``:
        ``equilibrium_all_a_moles(n_total_all, V_liq, V_gas, T_K)``

    Parameters
    ----------
    kH_ref : dict
        Henry solubility at T_ref per species (mol m⁻³ Pa⁻¹, Sander
        convention).  Only species with kH_ref entries are solved.
    dlnH : dict
        d(ln kH)/d(1/T) per species (K).  Missing entries default to 0.
    T_ref : float
        Reference temperature for kH_ref (K).  Default 298.15 K.

    Examples
    --------
    >>> from PyOMES.chemistry import MultispeciesVLEPartition
    >>> # H2S and CO2 in a 1 L liquid / 0.1 L headspace system
    >>> H_H2S = 0.10 * 1000 / 101325       # 0.10 mol/L/atm → Sander units
    >>> H_CO2 = 3.4e-4 * 1000 / 101325     # 3.4e-4 mol/L/atm → Sander units
    >>> vle = MultispeciesVLEPartition(
    ...     kH_ref={"H2S": H_H2S, "CO2": H_CO2},
    ...     dlnH={"H2S": 2100.0, "CO2": 2400.0},
    ... )
    >>> n_liq = vle.equilibrium_all_a_moles(
    ...     {"H2S": 0.01, "CO2": 0.02}, V_liq=1.0, V_gas=0.1, T_K=308.15
    ... )
    """

    kH_ref: Dict[str, float]
    dlnH: Dict[str, float] = field(default_factory=dict)
    T_ref: float = 298.15

    def equilibrium_all_a_moles(
        self,
        n_total_all: Dict[str, float],
        capacity_a: float,
        capacity_b: float,
        T_K: float,
    ) -> Dict[str, float]:
        """Solve VLE for all species using IdealGasEOS (analytical, decoupled).

        For each species i with Henry constant kH_i(T):
            r_i = kH_i × R × T × V_liq / V_gas
            n_liq_i = r_i × n_total_i / (1 + r_i)

        Species not in ``kH_ref`` are ignored (returned as absent).
        """
        from PyOMES.equilibria.vle import IdealGasEOS  # local import avoids circular
        # Only IdealGasEOS path is implemented; non-ideal EOS raises NotImplementedError
        # at point of use (deferred: see THERMODYNAMIC_MODEL_ARCHITECTURE §CP7).
        result = {}
        for species_id, n_total in n_total_all.items():
            if species_id not in self.kH_ref:
                continue
            H_ref = self.kH_ref[species_id]
            dlnH  = self.dlnH.get(species_id, 0.0)
            kH = _kH_mol_L_atm_from_ref(H_ref, dlnH, T_K, self.T_ref)
            r = kH * _R_L_ATM_MOL_K * T_K * capacity_a / capacity_b
            result[species_id] = r * float(n_total) / (1.0 + r)
        return result
