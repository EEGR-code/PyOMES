# -*- coding: utf-8 -*-
"""External boundary objects for ControlVolume.

An :class:`ExternalBoundary` sits at the edge of a ControlVolume and
computes molar fluxes entering or leaving a single phase.  Unlike
internal :class:`PhaseInterface` objects (which transfer material between
two phases *within* a CV), external boundaries represent flows across the
system boundary — feeds, vents, drains, and connections to other CVs.

The distinction is physically meaningful: internal transfers conserve the
CV's total inventory, while external fluxes intentionally change it.

Concrete implementations provided here:

* :class:`GasFeed` — continuous volumetric gas supply (vvm-based)
* :class:`PressureReliefVent` — headspace pressure relief (instant or smooth)

Usage
-----
>>> feed = GasFeed(vvm_min=1.0, y={"O2": 0.21, "N2": 0.79})
>>> vent = PressureReliefVent(P_set_atm=1.0, mode="instant")
>>> rec = apply_boundary(feed, cv, dt_h=0.01)
>>> print(rec.mol_applied)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Protocol, runtime_checkable

from .phases import GasPhase, LiquidPhase, Phase, R_L_ATM_MOL_K
from ..control.descriptors import MutableDict as _MutableDict, MutableScalar as _MS


# ════════════════════════════════════════════════════════════════════════
#  Protocol
# ════════════════════════════════════════════════════════════════════════

@runtime_checkable
class ExternalBoundary(Protocol):
    """Contract for an external material boundary on a ControlVolume.

    Implementations must define:
      - ``phase_key``: the CV phase this boundary acts on.
      - ``compute_flux(cv, dt_h)``: returns species molar fluxes
        (mol/h) *into* the target phase.  Positive = entering the CV,
        negative = leaving the CV.

    The boundary receives the full CV so it can inspect any phase
    (e.g. a gas feed may need the liquid volume for vvm calculation).
    """

    @property
    def phase_key(self) -> str: ...

    @property
    def label(self) -> str: ...

    def compute_flux(
        self,
        cv: Any,  # ControlVolume (avoid circular import at protocol level)
        dt_h: float,
        instantaneous: bool = False,
    ) -> Dict[str, float]:
        """Compute species molar fluxes (mol/h) into the target phase.

        Parameters
        ----------
        cv : ControlVolume
            The control volume this boundary is attached to.
        dt_h : float
            Timestep duration (hours).
        instantaneous : bool
            If True, return the instantaneous rate at the current state
            (dn/dt).  Use this when calling from within an ODE derivative
            function, where the integrator handles time-stepping.

            If False (default), return a step-averaged rate suitable for
            discrete Euler-style integration over the full ``dt_h``.

        Returns
        -------
        dict
            ``{species_id: flux_mol_per_h}`` — positive = into the CV.
        """
        ...


# ════════════════════════════════════════════════════════════════════════
#  Diagnostics record
# ════════════════════════════════════════════════════════════════════════

@dataclass
class ExternalFluxRecord:
    """Auditable record of a single external-flux application.

    Returned by :func:`apply_boundary` for each boundary evaluation.

    Attributes
    ----------
    boundary_label : str
        Human-readable label for the boundary (e.g. ``"gas_feed"``).
    phase_key : str
        Which CV phase was targeted.
    flux_mol_per_h : dict
        Species fluxes as computed (mol/h).
    dt_h : float
        Timestep duration (hours).
    mol_applied : dict
        Actual moles added/removed: ``flux × dt`` for each species.
    """
    boundary_label: str
    phase_key: str
    flux_mol_per_h: Dict[str, float] = field(default_factory=dict)
    dt_h: float = 0.0
    mol_applied: Dict[str, float] = field(default_factory=dict)

    @property
    def total_mol_applied(self) -> float:
        """Net moles added to (positive) or removed from (negative) the CV."""
        return sum(self.mol_applied.values())

    def summary(self) -> str:
        lines = [f"ExternalFluxRecord({self.boundary_label!r} → {self.phase_key!r}):"]
        lines.append(f"  dt_h = {self.dt_h:.6g}")
        for sp in sorted(self.mol_applied.keys()):
            mol = self.mol_applied[sp]
            rate = self.flux_mol_per_h.get(sp, 0.0)
            lines.append(f"  {sp:12s}  {rate:+.6e} mol/h  →  {mol:+.6e} mol")
        lines.append(f"  net = {self.total_mol_applied:+.6e} mol")
        return "\n".join(lines)


# ════════════════════════════════════════════════════════════════════════
#  apply_boundary — orchestrator helper
# ════════════════════════════════════════════════════════════════════════

def apply_boundary(
    boundary: ExternalBoundary,
    cv: Any,  # ControlVolume
    dt_h: float,
) -> ExternalFluxRecord:
    """Evaluate a boundary, apply its flux to the CV, and return a record.

    This is the recommended way to apply external boundaries:

    >>> rec = apply_boundary(feed, cv, dt_h=0.01)
    >>> assert rec.mol_applied["O2"] > 0  # O2 was fed

    Parameters
    ----------
    boundary : ExternalBoundary
        The boundary object to evaluate.
    cv : ControlVolume
        The control volume to modify.
    dt_h : float
        Timestep duration (hours).

    Returns
    -------
    ExternalFluxRecord
    """
    flux = boundary.compute_flux(cv, dt_h)
    cv.apply_external_flux(boundary.phase_key, flux, dt_h)
    mol_applied = {sp: rate * dt_h for sp, rate in flux.items()}
    return ExternalFluxRecord(
        boundary_label=boundary.label,
        phase_key=boundary.phase_key,
        flux_mol_per_h=dict(flux),
        dt_h=float(dt_h),
        mol_applied=mol_applied,
    )


# ════════════════════════════════════════════════════════════════════════
#  GasFeed — continuous volumetric gas supply
# ════════════════════════════════════════════════════════════════════════

class GasFeed:
    """Continuous gas feed into the headspace, specified by vvm and composition.

    The molar feed rate is computed from the volumetric gas flow rate
    (vvm × V_liquid) and the ideal gas law at the inlet conditions.

    Parameters
    ----------
    vvm_min : float
        Gas volume per liquid volume per minute (L_gas / L_liq / min).
    y : dict
        Inlet gas composition as mole fractions, e.g.
        ``{"O2": 0.21, "N2": 0.79}``.  Normalised internally.
    P_inlet_atm : float
        Inlet gas pressure (atm).  Default 1.0.
    phase_key : str
        CV phase key for the gas phase (default ``"gas"``).
    liquid_phase_key : str
        CV phase key for the liquid phase (default ``"liquid"``).
        Used to read liquid volume for vvm calculation.
    label : str
        Human-readable label.
    V_liq_override : float or None
        If set, use this fixed liquid volume instead of reading from
        the CV's liquid phase.  Useful when no liquid phase exists.
    T_override_K : float or None
        If set, use this temperature instead of reading from the gas
        phase.  Useful for inlet-temperature corrections.
    """

    # vvm_min and y declared as descriptors so the reflective walker can
    # write to them without needing separate _set_*_unchecked methods.
    vvm_min: float = _MS(float)
    y: Dict[str, float] = _MutableDict(value_type=float)

    def __init__(
        self,
        vvm_min: float = 0.0,
        y: Optional[Dict[str, float]] = None,
        P_inlet_atm: float = 1.0,
        phase_key: str = "gas",
        liquid_phase_key: str = "liquid",
        label: str = "gas_feed",
        V_liq_override: Optional[float] = None,
        T_override_K: Optional[float] = None,
    ):
        self._phase_key = str(phase_key)
        self._liquid_phase_key = str(liquid_phase_key)
        self._label = str(label)

        self.P_inlet_atm = float(P_inlet_atm)
        self.V_liq_override = float(V_liq_override) if V_liq_override is not None else None
        self.T_override_K = float(T_override_K) if T_override_K is not None else None

        # RunContext back-reference — wired by Simulation._wire_cv_context
        # when this boundary is owned by a running Simulation.
        self._context = None

        # Descriptor initialisation (_context=None so raise_if_running passes).
        self.vvm_min = float(vvm_min)

        # Normalise composition and store directly into the backing slot.
        raw = dict(y or {"O2": 0.21, "N2": 0.79})
        y_sum = sum(max(0.0, float(v)) for v in raw.values())
        if y_sum > 0.0:
            normalised = {k: max(0.0, float(v)) / y_sum for k, v in raw.items()}
        else:
            normalised = {k: 0.0 for k in raw}
        self.__dict__["_y"] = normalised

    # ── Protocol properties ─────────────────────────────────────────

    @property
    def phase_key(self) -> str:
        return self._phase_key

    @property
    def label(self) -> str:
        return self._label


    # ── Flux computation ────────────────────────────────────────────

    def compute_flux(self, cv: Any, dt_h: float,
                     instantaneous: bool = False) -> Dict[str, float]:
        """Compute species feed rates (mol/h) into the gas phase.

        The calculation:
          Q_gas = vvm × V_liq  (L/min)
          n_dot_total = P_inlet × Q_gas × 60 / (R × T)  (mol/h)
          n_dot_i = y_i × n_dot_total

        Returns
        -------
        dict
            ``{species: mol_per_h}`` — always non-negative.
        """
        if self.vvm_min <= 0.0:
            return {k: 0.0 for k in self._y}

        # Liquid volume
        if self.V_liq_override is not None:
            V_liq = self.V_liq_override
        elif self._liquid_phase_key in cv.phases:
            V_liq = float(cv.phases[self._liquid_phase_key].V_L)
        else:
            return {k: 0.0 for k in self._y}

        if V_liq <= 0.0:
            return {k: 0.0 for k in self._y}

        # Temperature
        if self.T_override_K is not None:
            T_K = self.T_override_K
        elif self._phase_key in cv.phases:
            T_K = float(cv.phases[self._phase_key].T_K)
        else:
            T_K = 298.15

        T_K = max(T_K, 1.0)

        Q_gas_L_per_min = self.vvm_min * V_liq
        n_dot_total = (self.P_inlet_atm * Q_gas_L_per_min * 60.0) / (R_L_ATM_MOL_K * T_K)

        return {sp: yi * n_dot_total for sp, yi in self._y.items()}

    def __repr__(self):
        species = ", ".join(f"{k}={v:.3f}" for k, v in self._y.items())
        return f"GasFeed(vvm={self.vvm_min}, P_in={self.P_inlet_atm} atm, y=[{species}])"


# ════════════════════════════════════════════════════════════════════════
#  PressureReliefVent — headspace pressure relief
# ════════════════════════════════════════════════════════════════════════

def _smooth_vent_fraction(
    excess_atm: float,
    dt_h: float,
    k_vent_per_h: float,
    smooth_width_atm: float,
) -> float:
    """Smooth venting fraction (C¹ approximation to step relief).

    Returns the fraction of total headspace moles to vent over the
    timestep.  Uses a smooth-ReLU activation:
      smooth_excess = 0.5 × (x + √(x² + w²))
    where x = P - P_set and w = smooth_width_atm.
    """
    excess = float(excess_atm)
    dt = max(0.0, float(dt_h))
    k = max(0.0, float(k_vent_per_h))
    w = max(1e-12, float(smooth_width_atm))

    smooth_excess = 0.5 * (excess + math.sqrt(excess * excess + w * w))
    if smooth_excess <= 0.0 or dt <= 0.0 or k <= 0.0:
        return 0.0

    rate_per_h = k * (smooth_excess / max(w, 1e-12))
    frac = 1.0 - math.exp(-rate_per_h * dt)
    return min(max(frac, 0.0), 0.999999999)


class PressureReliefVent:
    """Headspace pressure relief vent.

    Removes gas proportionally across species when headspace pressure
    exceeds the setpoint.  Two modes are available:

    ``"instant"``
        Vent exactly the moles needed to reach ``P_set_atm`` in one step.
        This is a hard clamp: if ``P > P_set``, remove ``n_total - n_target``
        moles, split by current mole fractions.

    ``"smooth"``
        Vent a fraction of headspace moles that increases smoothly as
        pressure rises above ``P_set_atm``.  Controlled by ``k_vent_per_h``
        (rate constant) and ``smooth_width_atm`` (activation sharpness).
        More physical and better-conditioned for stiff integrators.

    Flux convention: returns *negative* fluxes (material leaving the CV).

    Parameters
    ----------
    P_set_atm : float
        Pressure setpoint (atm).  Venting activates when P > P_set.
    mode : str
        ``"instant"`` or ``"smooth"``.
    k_vent_per_h : float
        Rate constant for smooth mode (1/h).  Default 500.
    smooth_width_atm : float
        Smooth-ReLU width for smooth mode (atm).  Default 0.01.
    phase_key : str
        CV phase key for the gas phase.  Default ``"gas"``.
    label : str
        Human-readable label.
    """

    def __init__(
        self,
        P_set_atm: float = 1.0,
        mode: str = "instant",
        k_vent_per_h: float = 500.0,
        smooth_width_atm: float = 0.01,
        phase_key: str = "gas",
        label: str = "pressure_relief_vent",
    ):
        if mode not in ("instant", "smooth"):
            raise ValueError(f"mode must be 'instant' or 'smooth', got {mode!r}")

        self._phase_key = str(phase_key)
        self._label = str(label)
        self.P_set_atm = float(P_set_atm)
        self.mode = str(mode)
        self.k_vent_per_h = float(k_vent_per_h)
        self.smooth_width_atm = float(smooth_width_atm)

    # ── Protocol properties ─────────────────────────────────────────

    @property
    def phase_key(self) -> str:
        return self._phase_key

    @property
    def label(self) -> str:
        return self._label

    # ── Flux computation ────────────────────────────────────────────

    def compute_flux(self, cv: Any, dt_h: float,
                     instantaneous: bool = False) -> Dict[str, float]:
        """Compute vent fluxes (mol/h) — negative values leaving the gas phase.

        Parameters
        ----------
        cv : ControlVolume
        dt_h : float
            Timestep (hours).

        Returns
        -------
        dict
            ``{species: flux_mol_per_h}``.  Negative = leaving the CV.
        """
        dt_h = max(float(dt_h), 1e-30)
        gas = cv.phases[self._phase_key]

        n_mol = gas.n_mol
        n_total = sum(max(0.0, float(v)) for v in n_mol.values())
        if n_total <= 0.0:
            return {k: 0.0 for k in n_mol}

        P = gas.P_atm
        if P <= self.P_set_atm:
            return {k: 0.0 for k in n_mol}

        if self.mode == "instant":
            vent_mol = self._instant_vent_mol(gas, n_total)
        else:
            vent_mol = self._smooth_vent_mol(P, n_total, dt_h)

        if vent_mol <= 0.0:
            return {k: 0.0 for k in n_mol}

        # Clamp to available moles
        vent_mol = min(vent_mol, n_total * 0.999999999)

        # Split proportionally by mole fraction
        y = {k: max(0.0, float(v)) / n_total for k, v in n_mol.items()} if n_total > 0 else {}

        # Return as negative flux (leaving the CV), expressed as rate
        return {sp: -(yi * vent_mol) / dt_h for sp, yi in y.items() if yi > 0.0}

    def _instant_vent_mol(self, gas: GasPhase, n_total: float) -> float:
        """Moles to vent for instant (hard-clamp) pressure relief."""
        T_K = max(float(gas.T_K), 1.0)
        V_L = max(float(gas.V_L), 1e-30)
        n_target = (self.P_set_atm * V_L) / (R_L_ATM_MOL_K * T_K)
        return max(0.0, n_total - max(0.0, n_target))

    def _smooth_vent_mol(self, P_atm: float, n_total: float, dt_h: float) -> float:
        """Moles to vent for smooth (C¹) pressure relief."""
        excess = P_atm - self.P_set_atm
        frac = _smooth_vent_fraction(excess, dt_h, self.k_vent_per_h, self.smooth_width_atm)
        return frac * n_total

    def __repr__(self):
        return (f"PressureReliefVent(P_set={self.P_set_atm} atm, "
                f"mode={self.mode!r})")


# ════════════════════════════════════════════════════════════════════════
#  MembraneGasBoundary — gas-permeable membrane to external atmosphere
# ════════════════════════════════════════════════════════════════════════

class MembraneGasBoundary:
    """Gas-permeable membrane connecting headspace to an external atmosphere.

    Models a membrane (e.g. silicone film, breathable well-plate seal)
    through which gas species permeate driven by partial-pressure
    differences between the headspace and the external atmosphere.

    The molar flux for each species is::

        flux_i = permeability_i × area × (p_external_i − p_headspace_i)

    Positive flux = gas enters the headspace (permeation inward, e.g.
    O₂ replenishment).  Negative flux = gas leaves the headspace
    (permeation outward, e.g. CO₂ venting).

    The external atmosphere composition is treated as a fixed boundary
    condition (infinite reservoir) — it does not change regardless of
    the flux.

    Parameters
    ----------
    permeability : dict
        ``{species_id: permeability_mol_per_h_per_atm_per_m2}`` —
        membrane permeability for each transferable species.  Species
        not listed are impermeable.
    area_m2 : float
        Membrane area (m²).
    external_atmosphere : dict
        ``{species_id: partial_pressure_atm}`` — partial pressures in
        the external atmosphere.  Defaults to standard dry air.
    phase_key : str
        CV phase key for the headspace gas (default ``"gas"``).
    label : str
        Human-readable label.

    Notes
    -----
    Typical permeability values for silicone membranes (rough order of
    magnitude, in mol/h/atm/m²):

    - O₂: 0.01 – 0.1
    - CO₂: 0.05 – 0.5  (CO₂ permeability is typically 3–10× that of O₂)
    - N₂: 0.002 – 0.02 (N₂ permeability is typically 0.2× that of O₂)

    These values vary widely with membrane material, thickness, and
    temperature.  Always check against manufacturer data sheets.

    Example
    -------
    >>> membrane = MembraneGasBoundary(
    ...     permeability={"O2": 0.01, "CO2": 0.05, "N2": 0.002},
    ...     area_m2=3.14e-5,   # ~6 mm diameter well
    ... )
    >>> flux = membrane.compute_flux(well_plate_glv, dt_h=0.01)
    """

    # Standard dry air partial pressures at 1 atm (approximate)
    _DEFAULT_ATMOSPHERE = {
        "O2": 0.2095,
        "CO2": 0.0004,
        "N2": 0.7808,
    }

    def __init__(
        self,
        permeability: Optional[Dict[str, float]] = None,
        area_m2: float = 1e-4,
        external_atmosphere: Optional[Dict[str, float]] = None,
        phase_key: str = "gas",
        label: str = "membrane",
    ):
        self._phase_key = str(phase_key)
        self._label = str(label)
        self.area_m2 = float(area_m2)

        # Permeability: mol/h/atm/m² for each species
        self.permeability = {
            k: float(v) for k, v in (permeability or {}).items()
        }

        # External atmosphere partial pressures
        if external_atmosphere is not None:
            self.external_atmosphere = {
                k: float(v) for k, v in external_atmosphere.items()
            }
        else:
            self.external_atmosphere = dict(self._DEFAULT_ATMOSPHERE)

    # ── Protocol properties ─────────────────────────────────────────

    @property
    def phase_key(self) -> str:
        return self._phase_key

    @property
    def label(self) -> str:
        return self._label

    # ── Flux computation ────────────────────────────────────────────

    def compute_flux(self, cv: Any, dt_h: float,
                     instantaneous: bool = False) -> Dict[str, float]:
        """Compute membrane permeation fluxes (mol/h).

        For each permeable species::

            flux_i = permeability_i × area × (p_ext_i − p_hs_i)

        Positive = entering the headspace (e.g. O₂ from air).
        Negative = leaving the headspace (e.g. CO₂ to air).

        Parameters
        ----------
        cv : ControlVolume
            Must have ``cv.phases["gas"]`` with a ``p_atm`` property.
        dt_h : float
            Timestep duration (hours).

        Returns
        -------
        dict
            ``{species_id: flux_mol_per_h}``
        """
        gas = cv.phases.get(self._phase_key)
        if gas is None:
            return {}

        p_headspace = gas.p_atm if hasattr(gas, "p_atm") else {}
        A = self.area_m2

        flux: Dict[str, float] = {}
        for species, perm in self.permeability.items():
            perm = float(perm)
            if perm <= 0.0:
                continue

            p_ext = float(self.external_atmosphere.get(species, 0.0))
            p_hs = float(p_headspace.get(species, 0.0))

            rate = perm * A * (p_ext - p_hs)
            if abs(rate) > 0.0:
                flux[species] = rate

        return flux

    # ── Mutators ────────────────────────────────────────────────────

    def set_external_atmosphere(self, atmosphere: Dict[str, float]) -> None:
        """Update the external atmosphere composition.

        Useful for modelling sealed incubators where the atmosphere
        changes, or for simulating altitude or hypoxic conditions.
        """
        self.external_atmosphere = {
            k: float(v) for k, v in atmosphere.items()
        }

    def set_permeability(self, species: str, value: float) -> None:
        """Update the permeability for a single species."""
        self.permeability[species] = float(value)

    # ── Representation ──────────────────────────────────────────────

    def __repr__(self):
        species = sorted(self.permeability.keys())
        perm_str = ", ".join(f"{sp}={self.permeability[sp]:.3g}" for sp in species)
        return (f"MembraneGasBoundary(area={self.area_m2:.4g} m², "
                f"perm=[{perm_str}])")


# ════════════════════════════════════════════════════════════════════════
#  LiquidFeed — continuous liquid feed into the reactor
# ════════════════════════════════════════════════════════════════════════

class LiquidFeed:
    """Continuous liquid feed into a reactor's liquid phase.

    Models the inlet stream of a CSTR, fed-batch, or perfusion reactor.
    Each species enters at a rate determined by the feed volumetric flow
    rate and the feed composition::

        flux_i = C_feed_i × Q   (mol/h, positive = into reactor)

    Parameters
    ----------
    Q_L_per_h : float
        Volumetric feed flow rate (L/h).  Must be >= 0.
    feed_conc_mol_L : dict
        ``{species_id: concentration_mol_L}`` — composition of the
        feed stream.  Species not listed have zero feed concentration.
    phase_key : str
        CV phase key (default ``"liquid"``).
    label : str
        Human-readable label.

    Example
    -------
    >>> feed = LiquidFeed(
    ...     Q_L_per_h=320.0,
    ...     feed_conc_mol_L={"AceticAcid": 0.0167, "Yeast": 0.0001},
    ... )
    """

    def __init__(
        self,
        Q_L_per_h: float = 0.0,
        feed_conc_mol_L: Optional[Dict[str, float]] = None,
        phase_key: str = "liquid",
        label: str = "liquid_feed",
    ):
        self._phase_key = str(phase_key)
        self._label = str(label)
        self.Q_L_per_h = float(Q_L_per_h)
        self.feed_conc_mol_L = {
            k: float(v) for k, v in (feed_conc_mol_L or {}).items()
        }

    @property
    def phase_key(self) -> str:
        return self._phase_key

    @property
    def label(self) -> str:
        return self._label

    def compute_flux(self, cv: Any, dt_h: float,
                     instantaneous: bool = False) -> Dict[str, float]:
        """Compute feed fluxes (mol/h) into the liquid phase.

        flux_i = C_feed_i × Q  (always non-negative).

        The feed rate is independent of ``dt_h`` and ``instantaneous``
        — it is always an instantaneous rate.
        """
        Q = self.Q_L_per_h
        if Q <= 0.0:
            return {}
        return {sp: C * Q for sp, C in self.feed_conc_mol_L.items() if C > 0}

    def __repr__(self):
        n_sp = len(self.feed_conc_mol_L)
        return f"LiquidFeed(Q={self.Q_L_per_h:.4g} L/h, {n_sp} species)"


# ════════════════════════════════════════════════════════════════════════
#  LiquidDrain — continuous liquid drain from a perfectly mixed reactor
# ════════════════════════════════════════════════════════════════════════

class LiquidDrain:
    """Continuous liquid drain from a perfectly mixed reactor.

    Models the outlet stream of a CSTR or perfusion reactor.  All
    species leave at a rate proportional to their current reactor
    concentration (perfect mixing assumption)::

        flux_i = −C_reactor_i × Q   (mol/h, negative = out of reactor)

    Uses an exponential formulation for unconditional stability::

        effective_flux_i = −n_i × (1 − exp(−D × dt)) / dt

    where D = Q / V_liq is the dilution rate.

    Parameters
    ----------
    Q_L_per_h : float
        Volumetric drain flow rate (L/h).  Must be >= 0.
    species_filter : set or None
        If provided, only these species are drained.  Use this to model
        cell-recycle or membrane retention (e.g. exclude biomass).
    phase_key : str
        CV phase key (default ``"liquid"``).
    label : str
        Human-readable label.

    Example
    -------
    >>> drain = LiquidDrain(Q_L_per_h=320.0)
    >>> # With cell recycle (retain biomass):
    >>> drain = LiquidDrain(Q_L_per_h=320.0,
    ...     species_filter={"AceticAcid", "CO2", "O2", "N2"})
    """

    def __init__(
        self,
        Q_L_per_h: float = 0.0,
        species_filter: Optional[set] = None,
        phase_key: str = "liquid",
        label: str = "liquid_drain",
    ):
        self._phase_key = str(phase_key)
        self._label = str(label)
        self.Q_L_per_h = float(Q_L_per_h)
        self.species_filter = set(species_filter) if species_filter is not None else None

    @property
    def phase_key(self) -> str:
        return self._phase_key

    @property
    def label(self) -> str:
        return self._label

    def compute_flux(self, cv: Any, dt_h: float,
                     instantaneous: bool = False) -> Dict[str, float]:
        """Compute drain fluxes (mol/h) out of the liquid phase.

        When ``instantaneous=False`` (default), uses an exponential
        formulation for unconditional stability at any D × dt product.
        This gives the step-averaged removal rate over ``dt_h``.

        When ``instantaneous=True``, returns the simple linear rate
        ``−D × n`` (the true derivative dn/dt at the current state).
        Use this when calling from an ODE derivative function.
        """
        Q = self.Q_L_per_h
        if Q <= 0.0:
            return {}

        liq = cv.phases.get(self._phase_key)
        if liq is None:
            return {}

        V = float(liq.V_L)
        D = Q / V  # dilution rate (1/h)

        if instantaneous:
            # Instantaneous rate: dn/dt = -D × n
            flux: Dict[str, float] = {}
            for sp, n in liq.n_mol.items():
                if n <= 0.0:
                    continue
                if self.species_filter is not None and sp not in self.species_filter:
                    continue
                flux[sp] = -(D * n)
            return flux

        # Step-averaged exponential drain: n(dt) = n₀ × exp(-D×dt)
        dt = max(float(dt_h), 1e-30)
        exponent = D * dt
        if exponent > 50.0:
            frac = 1.0
        else:
            frac = 1.0 - math.exp(-exponent)

        flux: Dict[str, float] = {}
        for sp, n in liq.n_mol.items():
            if n <= 0.0:
                continue
            if self.species_filter is not None and sp not in self.species_filter:
                continue
            flux[sp] = -(n * frac) / dt
        return flux

    def __repr__(self):
        filt = f", filter={sorted(self.species_filter)}" if self.species_filter else ""
        return f"LiquidDrain(Q={self.Q_L_per_h:.4g} L/h{filt})"


# ════════════════════════════════════════════════════════════════════════
#  ProportionalGasOutlet — BSM2-style continuous gas outlet
# ════════════════════════════════════════════════════════════════════════

class ProportionalGasOutlet:
    """BSM2-style proportional gas outlet.

    Gas flow rate is proportional to overpressure::

        q_gas = k_p × (P_gas − P_atm)    [L/h]

    Each tracked species is carried out proportionally to its partial
    pressure (equivalently, its mole fraction in the vent stream).

    When ``include_water_vapour=True``, the total headspace pressure
    used for the driving force includes the analytical water vapour
    saturation pressure (Raoult's law, Clausius-Clapeyron)::

        P_gas = P_dry_species + P_water_sat(T)

    This matches BSM2's treatment where water vapour contributes to
    total pressure and therefore to the gas outlet flow, but water
    is not tracked as a gas-phase state variable.  The water vapour
    molar loss rate is stored in :attr:`last_water_loss_mol_per_h`
    for the liquid-phase water sink (stage 3).

    Flux convention: returns *negative* fluxes (material leaving the CV).

    Parameters
    ----------
    k_p_L_per_h_per_atm : float
        Pipe friction coefficient (L/h/atm).
        BSM2 default: 5e4 m³/d/bar ≈ 2.111e6 L/h/atm.
    P_atm : float
        Downstream (atmospheric) pressure (atm).  Default 1.0.
    include_water_vapour : bool
        If True, add the analytical water vapour saturation pressure
        to the total headspace pressure for the overpressure driving
        force.  Default False (backward compatible).
    phase_key : str
        CV phase key for the gas phase.  Default ``"gas"``.
    label : str
        Human-readable label.

    Notes
    -----
    Conversion from BSM2 units::

        k_p_BSM2 = 5e4 m³/d/bar
        k_p_L_per_h_per_atm = k_p_BSM2 × 1000 / 24 × 1.01325
                             ≈ 2.111e6 L/h/atm
    """

    def __init__(
        self,
        k_p_L_per_h_per_atm: float = 2.0e6,
        P_atm: float = 1.0,
        include_water_vapour: bool = False,
        phase_key: str = "gas",
        label: str = "proportional_gas_outlet",
    ):
        self._phase_key = str(phase_key)
        self._label = str(label)
        self.k_p = float(k_p_L_per_h_per_atm)
        self.P_atm = float(P_atm)
        self.include_water_vapour = bool(include_water_vapour)

        # Water vapour loss rate from the most recent compute_flux call.
        # Read by stage 3 (liquid H₂O sink boundary) to remove the
        # corresponding water moles from the liquid phase.
        self.last_water_loss_mol_per_h: float = 0.0

    # ── Protocol properties ─────────────────────────────────────────

    @property
    def phase_key(self) -> str:
        return self._phase_key

    @property
    def label(self) -> str:
        return self._label

    # ── Flux computation ────────────────────────────────────────────

    def compute_flux(self, cv: Any, dt_h: float,
                     instantaneous: bool = False) -> Dict[str, float]:
        """Compute continuous gas outlet fluxes (mol/h).

        The gas outlet rate is always instantaneous (proportional to
        current overpressure).  The ``instantaneous`` parameter is
        accepted for protocol consistency but does not change behaviour.

        Parameters
        ----------
        cv : ControlVolume or snapshot wrapper
            Must have ``cv.phases["gas"]`` with ``P_atm``, ``T_K``,
            ``V_L``, and ``n_mol``.
        dt_h : float
            Timestep (hours).  Not used directly (flux is a rate).
        instantaneous : bool
            Accepted for protocol consistency.  No effect — the gas
            outlet always returns an instantaneous rate.

        Returns
        -------
        dict
            ``{species: flux_mol_per_h}``.  Negative = leaving the CV.
            Only tracked (dry) species appear.  Water vapour loss is
            stored in :attr:`last_water_loss_mol_per_h`.
        """
        gas = cv.phases.get(self._phase_key)
        if gas is None:
            self.last_water_loss_mol_per_h = 0.0
            return {}

        T_K = max(float(gas.T_K), 1.0)

        # Total headspace pressure for the driving force
        P_dry = gas.P_atm
        if self.include_water_vapour:
            P_water = gas.P_water_sat_atm
        else:
            P_water = 0.0
        P_gas = P_dry + P_water

        overpressure = P_gas - self.P_atm
        if overpressure <= 0.0:
            self.last_water_loss_mol_per_h = 0.0
            return {k: 0.0 for k in gas.n_mol}

        # Volumetric gas flow rate (L/h) — total (wet) gas
        q_gas_L_h = self.k_p * overpressure

        V_gas = max(float(gas.V_L), 1e-30)

        # Each tracked species i has partial pressure p_i = n_i RT / V.
        # Its molar removal rate is:
        #   n_dot_i = p_i × q_gas / (R × T) = n_i × q_gas / V_gas
        # This is independent of whether water vapour is included —
        # the water fraction of the vent stream does not dilute the
        # tracked species' concentrations, because water vapour is
        # analytically imposed (not tracked as moles in the gas phase).
        flux: Dict[str, float] = {}
        for sp, n in gas.n_mol.items():
            n = max(0.0, float(n))
            if n > 0.0:
                flux[sp] = -(n * q_gas_L_h / V_gas)

        # Water vapour loss: the vent stream also carries water vapour
        # at partial pressure P_water_sat.
        #   n_dot_water = P_water × q_gas / (R × T)
        if self.include_water_vapour and P_water > 0.0:
            self.last_water_loss_mol_per_h = (
                P_water * q_gas_L_h / (R_L_ATM_MOL_K * T_K)
            )
        else:
            self.last_water_loss_mol_per_h = 0.0

        return flux

    def __repr__(self):
        wv = ", water_vapour=True" if self.include_water_vapour else ""
        return (f"ProportionalGasOutlet(k_p={self.k_p:.3e} L/h/atm, "
                f"P_atm={self.P_atm} atm{wv})")


