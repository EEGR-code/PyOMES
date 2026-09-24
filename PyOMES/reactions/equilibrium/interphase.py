# -*- coding: utf-8 -*-
"""Henry, Raoult and Ksp phase-equilibrium constraints.

Each class here is declared once and plays two roles:

- As a :class:`~PyOMES.chemistry.partition.PartitionModel` it describes how a
  species distributes between two phases at equilibrium
  (``partition_ratio``/``equilibrium_a_moles``).
- As an :class:`~PyOMES.reactions.equilibrium.EquilibriumConstraint`
  (``stoichiometry``/``log_K``/``dH_J_per_mol``/``T_ref_K``) it can go into the
  reaction list feeding ``NRChemicalEquilibriumEngine``/
  ``ChemicalEquilibriumEngine``.

The same constructed instance can therefore be passed both to
``KineticGasLiquidLink``'s ``partition_models=`` dict and into a
``ReactionSystem``.

- :class:`HenryEquilibrium` — linear (Henry's law) gas-liquid partition.
- :class:`RaoultEquilibrium` — Raoult's law vapour-liquid partition for the
  solvent (H₂O).
- :class:`KspEquilibrium` — solubility-product solid-liquid equilibrium.

The ``partition_ratio()`` return value selects the solver path in
``KineticGasLiquidLink``:
  - float  → analytical exponential step solution (exact for linear models)
  - None   → ODE instantaneous-rate path (required for non-linear models)

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
from dataclasses import dataclass, field
from typing import Dict, Optional, Sequence, Tuple, Union

from PyOMES.chemistry import common_species
from PyOMES.chemistry.partition import _kH_mol_L_atm_from_ref
from PyOMES.chemistry.species import Species
from PyOMES.units import R_J_PER_MOL_K as _R_J_MOL
from PyOMES.units import R_L_ATM_PER_MOL_K
from .constraint import vant_hoff_log_K
from PyOMES.reactions.stoichiometry import StoichiometryEntry, _parse_stoichiometry


def _resolve_species(species: Union[str, Species, None]) -> Optional[Species]:
    """Resolve a species id string or ``Species`` object to a ``Species``.

    Looks the id up in ``PyOMES.chemistry.common_species`` directly — the
    same source :func:`PyOMES.reactions.stoichiometry._get_common_species`
    builds its own lookup from.
    """
    if species is None or isinstance(species, Species):
        return species
    lookup = {
        v.id: v for v in vars(common_species).values() if isinstance(v, Species)
    }
    if species not in lookup:
        raise ValueError(
            f"Unknown species id {species!r} — pass a Species object "
            "directly, or use an id from PyOMES.chemistry.common_species."
        )
    return lookup[species]


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
        return (kH / (max(1e-12, alpha) * gamma)) * R_L_ATM_PER_MOL_K * T_K * capacity_a / capacity_b

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
        return _kH_mol_L_atm_from_ref(self.H_ref, self.dlnH, T_K, self.T_ref)

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
        """Van 't Hoff enthalpy derived from dlnH: dH = -dlnH × R."""
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
        gas = _resolve_species(self.gas_species)
        liq = _resolve_species(self.liquid_species)
        return (
            StoichiometryEntry(species=gas, phase="gas", coefficient=-1.0),
            StoichiometryEntry(species=liq, phase="liquid", coefficient=+1.0),
        )


# ── Raoult's law partition for solvent (H₂O) ───────────────────────────────

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
        return (self.C_water_mol_L * capacity_a * R_L_ATM_PER_MOL_K * T_K) / (Ps * capacity_b)

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
        """Van 't Hoff enthalpy: dH = -dH_vap (Clausius-Clapeyron sign flip)."""
        return -self.dH_vap

    @property
    def stoichiometry(self) -> Tuple["StoichiometryEntry", ...]:
        """gas ⇌ liquid stoichiometry, when gas_species/liquid_species are set.

        Returns an empty tuple when either is unset rather than
        raising, so ``isinstance(x, EquilibriumConstraint)`` stays safe.
        """
        if self.gas_species is None or self.liquid_species is None:
            return ()
        gas = _resolve_species(self.gas_species)
        liq = _resolve_species(self.liquid_species)
        return (
            StoichiometryEntry(species=gas, phase="gas", coefficient=-1.0),
            StoichiometryEntry(species=liq, phase="liquid", coefficient=+1.0),
        )


# ── Ksp (solubility-product) solid-liquid equilibrium ───────────────────────

class KspEquilibrium:
    """Solubility-product (Ksp) solid-liquid equilibrium.

    A named constructor for the solubility-product convention
    (``log_K = log10(Ksp)``) satisfying both ``PartitionModel`` and
    :class:`~PyOMES.reactions.equilibrium.EquilibriumConstraint`, instead
    of requiring hand-built
    :class:`~PyOMES.reactions.equilibrium.EquilibriumReaction` objects.

    The ``PartitionModel`` role (``partition_ratio``/
    ``equilibrium_a_moles``) is well-defined only for the single-ion
    case — exactly one solid entry and one dissolved entry in
    ``stoichiometry`` — because a genuinely multi-ion solubility
    equilibrium requires the active-set NR solver, not a closed-form
    partition ratio. ``partition_ratio()`` returns
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
        Ksp_T = 10.0 ** vant_hoff_log_K(self, T_K)
        return min(n_total, Ksp_T * capacity_a)

    def __repr__(self) -> str:
        lbl = f" [{self.label}]" if self.label else ""
        return f"KspEquilibrium(Ksp={self.Ksp:.3g}{lbl})"
