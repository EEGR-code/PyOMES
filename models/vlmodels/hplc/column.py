# -*- coding: utf-8 -*-
"""HPLC column model: 1D advection-dispersion with Langmuir adsorption.

Models the elution of solute species through a packed chromatography
column by discretising the column into N finite-volume cells, each
containing a liquid mobile phase and a solid stationary phase with
Langmuir adsorption kinetics.

The model solves the coupled system simultaneously using an adaptive
ODE integrator (``scipy.integrate.solve_ivp``), giving accurate
elution profiles even for sharp fronts and closely eluting peaks.

Mathematical formulation
------------------------
For each cell *i* (of *N*) and each species *s*:

**Mobile phase (liquid):**

.. math::

    \\varepsilon \\frac{\\partial C_i^s}{\\partial t} =
        \\frac{Q}{V_{cell}} (C_{i-1}^s - C_i^s)
        + \\frac{D_{ax}}{\\Delta z^2} (C_{i-1}^s - 2C_i^s + C_{i+1}^s)
        - (1-\\varepsilon) \\rho_b \\frac{\\partial q_i^s}{\\partial t}

**Solid phase (adsorbed), linear driving force:**

.. math::

    \\frac{\\partial q_i^s}{\\partial t} = k_f^s (q_{eq}^s(C_i^s) - q_i^s)

**Langmuir equilibrium isotherm:**

.. math::

    q_{eq}^s = \\frac{q_{max}^s K^s C}{1 + \\sum_j K^j C_j}

where the denominator sums over all competing species (competitive
Langmuir).

Boundary conditions
-------------------
- **Inlet (cell 0):** Danckwerts BC with time-varying inlet concentration.
- **Outlet (cell N-1):** Zero-gradient (free outflow).

Usage
-----
>>> from PyOMES.models.hplc_column import HPLCColumn, LangmuirSpecies
>>> col = HPLCColumn(
...     length_cm=25.0, diameter_cm=0.46,
...     void_fraction=0.4, bulk_density_g_mL=0.6,
...     flow_rate_mL_min=1.0, n_cells=100,
...     species=[
...         LangmuirSpecies("glucose",  q_max=0.05, K_L=2.0, kf=10.0),
...         LangmuirSpecies("fructose", q_max=0.05, K_L=5.0, kf=10.0),
...         LangmuirSpecies("sucrose",  q_max=0.04, K_L=12.0, kf=8.0),
...     ],
... )
>>> result = col.simulate(t_end_min=30.0, inlet_fn=my_pulse)
>>> result.plot()

See Also
--------
:class:`~fermenter.core.phases.LiquidPhase` — mobile phase per cell.
:class:`~fermenter.core.phases.SolidPhase` — adsorbed inventory per cell.
:class:`~fermenter.core.links.AdvectiveLink` — inter-cell advective flow.
"""

from __future__ import annotations

import math
import time as _time
import warnings
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np


# ════════════════════════════════════════════════════════════════════════
#  Acid-base helpers
# ════════════════════════════════════════════════════════════════════════

def _f_protonated(pH, pKa):
    """Fraction of a monoprotic acid in its protonated (neutral) form.

    For an acid HA ⇌ H⁺ + A⁻ with dissociation constant Ka = 10^(-pKa):

        f_HA = [H⁺] / ([H⁺] + Ka)

    Parameters
    ----------
    pH : float or ndarray
        Solution pH.
    pKa : float
        Acid dissociation constant (negative log).

    Returns
    -------
    float or ndarray
        Fraction protonated (0–1).  At pH = pKa, f_HA = 0.5.
        At pH << pKa, f_HA → 1.  At pH >> pKa, f_HA → 0.
    """
    H = 10.0 ** (-np.asarray(pH, dtype=float))
    Ka = 10.0 ** (-float(pKa))
    return H / (H + Ka)


# ════════════════════════════════════════════════════════════════════════
#  Species definition
# ════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class LangmuirSpecies:
    """Langmuir isotherm parameters for one solute species.

    Parameters
    ----------
    name : str
        Species identifier (e.g. ``"glucose"``).
    q_max : float
        Maximum adsorption capacity (mol / g_solid).
    K_L : float
        Langmuir equilibrium constant (L / mol).  Higher K means
        stronger binding and longer retention.
    kf : float
        Linear driving force mass transfer coefficient (1/min).
        Controls how fast the system approaches equilibrium within
        each cell.  Typical HPLC values: 1–100 /min.
    D_ax : float or None
        Axial dispersion coefficient (cm²/min) for this species.
        If None, uses the column-level default.
    MW : float
        Molecular weight (g/mol).  Used for mass-based reporting.
    """
    name: str
    q_max: float = 0.05      # mol/g
    K_L: float = 5.0         # L/mol
    kf: float = 10.0         # 1/min
    D_ax: Optional[float] = None
    MW: float = 180.16       # g/mol (glucose default)


@dataclass(frozen=True)
class PartitionSpecies:
    """Linear partition/exclusion parameters for one solute species.

    Appropriate for ion-exclusion chromatography (e.g. Aminex HPX-87H)
    where retention is governed by partitioning of solute into stagnant
    pore liquid inside resin beads, not by adsorption to discrete sites.

    The equilibrium is a simple linear partition::

        C_internal = K_D × C_external

    where ``K_D`` is the distribution coefficient (0–1, dimensionless).
    There is no saturation, no competition between species, and no
    Langmuir denominator.

    For ion-exclusion columns with acidic mobile phase, ``K_D``
    encapsulates three effects:

    - **Size exclusion** — larger molecules access fewer pores
      (lower K_D).
    - **Donnan exclusion** — anionic species are repelled by fixed
      sulfonate charges on the resin (lower K_D for dissociated acids).
    - **Hydrophobic partitioning** — non-polar solutes interact with
      the polystyrene backbone (higher K_D).

    Parameters
    ----------
    name : str
        Species identifier.
    K_D : float
        Distribution coefficient (dimensionless, 0–1 for pure
        exclusion, >1 possible with hydrophobic partitioning).
        The fraction of the internal pore volume accessible to
        this species.  Higher K_D → longer retention.
    kf : float
        Linear driving force mass transfer coefficient (1/min).
    intraparticle_porosity : float or None
        Porosity within the resin bead (ε_p).  If None, uses the
        column-level default.  Typical: 0.3–0.6 for gel-type resins.
    D_ax : float or None
        Axial dispersion coefficient (cm²/min) for this species.
        If None, uses the column-level default.
    MW : float
        Molecular weight (g/mol).
    pKa : float or None
        Acid dissociation constant.  If provided, allows Donnan
        exclusion effects to be computed from mobile phase pH.
        Currently informational only.
    """
    name: str
    K_D: float = 0.5         # dimensionless
    kf: float = 30.0         # 1/min
    intraparticle_porosity: Optional[float] = None
    D_ax: Optional[float] = None
    MW: float = 180.16       # g/mol
    pKa: Optional[float] = None


# ════════════════════════════════════════════════════════════════════════
#  Simulation result
# ════════════════════════════════════════════════════════════════════════

@dataclass
class HPLCResult:
    """Results from an HPLC column simulation.

    Attributes
    ----------
    t_min : ndarray, shape (n_times,)
        Time points (minutes).
    C_outlet : dict
        ``{species_name: ndarray}`` — outlet (detector) concentration
        timeseries (mol/L) for each species.
    C_all : dict
        ``{species_name: ndarray shape (n_times, n_cells)}`` — full
        spatiotemporal concentration field.
    q_all : dict
        ``{species_name: ndarray shape (n_times, n_cells)}`` — adsorbed
        concentration (mol/g) in each cell over time.
    column : HPLCColumn
        Reference to the column that produced this result.
    runtime_s : float
        Wall-clock simulation time (seconds).
    solver_info : dict
        Solver diagnostics (n_eval, success, message).
    """
    t_min: np.ndarray
    C_outlet: Dict[str, np.ndarray]
    C_all: Dict[str, np.ndarray]
    q_all: Dict[str, np.ndarray]
    column: Any
    runtime_s: float = 0.0
    solver_info: Dict[str, Any] = field(default_factory=dict)

    def plot(self, figsize=(12, 5), species=None):
        """Plot the chromatogram (outlet concentration vs time).

        Parameters
        ----------
        figsize : tuple
            Figure size (width, height) in inches.
        species : list of str, optional
            Species to plot.  If None, plots all.

        Returns
        -------
        matplotlib.figure.Figure
        """
        import matplotlib.pyplot as plt

        names = species or list(self.C_outlet.keys())
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        for name in names:
            if name in self.C_outlet:
                ax.plot(self.t_min, self.C_outlet[name] * 1000,
                        label=name, lw=1.5)
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Concentration (mM)")
        ax.set_title("HPLC Chromatogram (detector at column outlet)")
        ax.legend()
        ax.set_xlim(0, self.t_min[-1])
        ax.set_ylim(bottom=0)
        return fig

    def plot_spatial(self, t_min_snapshot, figsize=(12, 5), species=None):
        """Plot concentration profile along the column at a given time.

        Parameters
        ----------
        t_min_snapshot : float
            Time (minutes) at which to take the snapshot.
        """
        import matplotlib.pyplot as plt

        names = species or list(self.C_all.keys())
        idx = np.argmin(np.abs(self.t_min - t_min_snapshot))

        z = np.linspace(0, self.column.length_cm, self.column.n_cells)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

        for name in names:
            ax1.plot(z, self.C_all[name][idx] * 1000, label=name, lw=1.5)
            ax2.plot(z, self.q_all[name][idx] * 1000, label=name, lw=1.5)

        ax1.set_xlabel("Position (cm)"); ax1.set_ylabel("C (mM)")
        ax1.set_title(f"Mobile phase at t = {self.t_min[idx]:.2f} min")
        ax1.legend()

        ax2.set_xlabel("Position (cm)"); ax2.set_ylabel("q (mmol/g)")
        ax2.set_title(f"Adsorbed phase at t = {self.t_min[idx]:.2f} min")
        ax2.legend()
        return fig


# ════════════════════════════════════════════════════════════════════════
#  Inlet functions (convenience)
# ════════════════════════════════════════════════════════════════════════

def pulse_inlet(
    species_conc: Dict[str, float],
    t_start_min: float = 0.0,
    duration_min: float = 0.5,
) -> Callable:
    """Create a rectangular pulse inlet function.

    Parameters
    ----------
    species_conc : dict
        ``{species_name: concentration_mol_L}`` during the pulse.
    t_start_min : float
        Start time of the pulse (minutes).
    duration_min : float
        Duration of the pulse (minutes).

    Returns
    -------
    callable
        ``inlet_fn(t_min) -> {species: C_mol_L}``
    """
    t_end = t_start_min + duration_min

    def fn(t):
        if t_start_min <= t <= t_end:
            return dict(species_conc)
        return {sp: 0.0 for sp in species_conc}

    return fn


def step_inlet(
    species_conc: Dict[str, float],
    t_start_min: float = 0.0,
) -> Callable:
    """Create a step-change inlet function.

    Returns constant concentration from t_start_min onwards.
    """
    def fn(t):
        if t >= t_start_min:
            return dict(species_conc)
        return {sp: 0.0 for sp in species_conc}

    return fn


def injection_inlet(
    sample_conc: Dict[str, float],
    mobile_phase_conc: Optional[Dict[str, float]] = None,
    t_start_min: float = 0.5,
    injection_volume_uL: float = 20.0,
    flow_rate_mL_min: float = 0.6,
    tau_mix_s: float = 0.5,
) -> Callable:
    """Create a realistic injection inlet with extra-column mixing.

    Models the injection as a rectangular sample plug convolved with
    an exponential decay that represents extra-column band broadening
    from the injection valve, connecting tubing, and column inlet frit.

    During the injection window, sample replaces mobile phase (constant
    total flow — the pump maintains volumetric flow rate, so the sample
    displaces mobile phase rather than adding to it).

    Outside the injection window, mobile phase flows continuously.
    This is important for species present in the mobile phase (e.g.
    sulfate from H₂SO₄) — they are absent from the sample plug,
    creating depletion zones.

    The exponential mixing model assumes a well-mixed dead volume
    V_mix at the column inlet.  When the injection valve switches,
    the concentration in V_mix evolves as::

        dC/dt = (Q / V_mix) × (C_in − C)

    giving exponential rise (at injection start) and decay (at
    injection end) with time constant τ = V_mix / Q.

    Parameters
    ----------
    sample_conc : dict
        ``{species_name: concentration_mol_L}`` of the sample.
        Species not listed are assumed absent from the sample (0.0).
    mobile_phase_conc : dict or None
        ``{species_name: concentration_mol_L}`` of the mobile phase.
        Species not listed are assumed absent (0.0).
        If None, all species are 0.0 in the mobile phase.
    t_start_min : float
        Time the injection valve switches (minutes).
    injection_volume_uL : float
        Sample loop volume (µL).
    flow_rate_mL_min : float
        Pump flow rate (mL/min).  Used to compute injection duration.
    tau_mix_s : float
        Extra-column mixing time constant (seconds).  Typical values:
        0.3–1.0 s for standard analytical HPLC with 0.17 mm tubing,
        0.1–0.3 s for UHPLC with low-dispersion fittings.
        Set to 0.0 for a perfect rectangular pulse (no mixing).

    Returns
    -------
    callable
        ``inlet_fn(t_min) -> {species: C_mol_L}``

    Notes
    -----
    The returned function evaluates the analytical solution of the
    exponential mixing model, avoiding any numerical ODE integration.
    For a rectangular input switching from C_mobile to C_sample at
    t_start and back at t_end::

        C(t) = C_mobile + (C_sample − C_mobile) × h(t)

    where h(t) is the convolution of the rectangular pulse with the
    exponential kernel::

        h(t) = 0                                      for t < t_start
        h(t) = 1 − exp(−(t − t_start)/τ)             for t_start ≤ t < t_end
        h(t) = [1 − exp(−(t − t_start)/τ)]
             × exp(−(t − t_end)/τ)                    simplified below

    This is the standard exponentially-modified rectangle (EMR).
    """
    if flow_rate_mL_min <= 0:
        raise ValueError(
            f"flow_rate_mL_min must be positive, got {flow_rate_mL_min}")
    if tau_mix_s < 0:
        raise ValueError(
            f"tau_mix_s must be non-negative, got {tau_mix_s}")
    if injection_volume_uL <= 0:
        raise ValueError(
            f"injection_volume_uL must be positive, got {injection_volume_uL}")

    mobile = dict(mobile_phase_conc) if mobile_phase_conc else {}
    sample = dict(sample_conc)

    # Ensure both dicts have the same keys
    all_species = set(mobile.keys()) | set(sample.keys())
    for sp in all_species:
        mobile.setdefault(sp, 0.0)
        sample.setdefault(sp, 0.0)

    duration_min = (injection_volume_uL / 1000.0) / flow_rate_mL_min
    t_end_min = t_start_min + duration_min
    tau_min = tau_mix_s / 60.0  # convert to minutes

    if tau_min < 1e-10:
        # No mixing: perfect rectangular pulse
        def fn(t):
            if t_start_min <= t <= t_end_min:
                return dict(sample)
            return dict(mobile)
        return fn

    inv_tau = 1.0 / tau_min

    def fn(t):
        if t < t_start_min:
            return dict(mobile)

        # Exponential rise from mobile → sample at t_start
        dt_start = t - t_start_min
        rise = 1.0 - math.exp(-dt_start * inv_tau)

        if t <= t_end_min:
            # Still injecting: exponential approach to sample conc
            frac_sample = rise
        else:
            # Injection ended: exponential decay back to mobile
            dt_end = t - t_end_min
            # Concentration at the moment injection ended
            rise_at_end = 1.0 - math.exp(-duration_min * inv_tau)
            frac_sample = rise_at_end * math.exp(-dt_end * inv_tau)

        # Linear blend between mobile and sample
        result = {}
        for sp in all_species:
            result[sp] = mobile[sp] + frac_sample * (sample[sp] - mobile[sp])
        return result

    return fn


# ════════════════════════════════════════════════════════════════════════
#  Divert / peak-shaving functions (for recycle chromatography)
# ════════════════════════════════════════════════════════════════════════

def make_divert_windows(
    windows: list,
) -> Callable:
    """Divert outlet to collection during specified time windows.

    During a divert window the column inlet receives pure mobile phase
    (all concentrations zero) rather than recycled outlet.  Material
    eluting during these windows is effectively collected and removed
    from the recycle loop.

    Parameters
    ----------
    windows : list of (t_start, t_end)
        Time windows (minutes) during which outlet is diverted.

    Returns
    -------
    callable
        ``divert_fn(t, C_outlet) -> C_inlet`` compatible with
        :meth:`HPLCColumn.simulate`.

    Examples
    --------
    Collect sucrose between 7 and 8 min on each pass (period 7.5 min):

    >>> windows = [(7.5 * n - 0.5, 7.5 * n + 0.5) for n in range(1, 15)]
    >>> divert = make_divert_windows(windows)
    """
    sorted_windows = sorted(windows)

    def fn(t, C_out):
        for t_start, t_end in sorted_windows:
            if t_start <= t <= t_end:
                return np.zeros_like(C_out)
            if t_start > t:
                break
        return C_out

    return fn


def make_divert_threshold(
    threshold_mol_L: float,
) -> Callable:
    """Divert outlet when total concentration exceeds a threshold.

    This mimics a detector-triggered collection valve: when the RI
    signal exceeds the threshold, the outlet is sent to collection
    rather than recycled.

    Parameters
    ----------
    threshold_mol_L : float
        Total outlet concentration (mol/L) above which diversion
        is triggered.

    Returns
    -------
    callable
        ``divert_fn(t, C_outlet) -> C_inlet``
    """
    def fn(t, C_out):
        if C_out.sum() > threshold_mol_L:
            return np.zeros_like(C_out)
        return C_out

    return fn


def make_peak_shave_schedule(
    column: 'HPLCColumn',
    species_to_collect: list,
    n_passes: int = 10,
    margin_factor: float = 2.5,
) -> Callable:
    """Auto-generate divert windows from predicted retention times.

    For each species in ``species_to_collect``, computes the predicted
    retention time and creates divert windows centred on each pass of
    that species.  The window half-width is set by ``margin_factor``
    times an estimated peak sigma (from plate count).

    This is the standard "heart-cutting" approach: well-resolved peaks
    are diverted to collection on each pass so they don't collide with
    poorly-resolved peaks on subsequent passes.

    Parameters
    ----------
    column : HPLCColumn
        The column (must have species registered).
    species_to_collect : list of str
        Species names to divert (e.g. ``["sucrose"]``).
    n_passes : int
        Number of recycle passes to generate windows for.
    margin_factor : float
        Window half-width = margin_factor * sigma_peak, where sigma_peak
        is estimated from the plate count.  Default 2.5 captures >99%
        of the peak area.

    Returns
    -------
    callable
        ``divert_fn(t, C_outlet) -> C_inlet``
    """
    windows = []
    for sp_name in species_to_collect:
        t_R = column.predicted_retention_min(sp_name)
        N_plates = max(column.effective_plates, 10.0)
        sigma = t_R / math.sqrt(N_plates)
        half_w = margin_factor * sigma
        for n in range(1, n_passes + 1):
            centre = n * t_R
            windows.append((centre - half_w, centre + half_w))

    return make_divert_windows(windows)


# ════════════════════════════════════════════════════════════════════════
#  Spatial discretisation schemes (from shared numerics module)
# ════════════════════════════════════════════════════════════════════════

from PyOMES.numerics.spatial import (                     # noqa: E402
    available_advection_schemes,
    available_dispersion_schemes,
    get_advection_scheme,
    get_dispersion_scheme,
    register_advection_scheme,      # re-exported for convenience
    register_dispersion_scheme,     # re-exported for convenience
)


# ════════════════════════════════════════════════════════════════════════
#  HPLC Column
# ════════════════════════════════════════════════════════════════════════

class HPLCColumn:
    """1D packed chromatography column with Langmuir adsorption.

    Parameters
    ----------
    length_cm : float
        Column length (cm).
    diameter_cm : float
        Column inner diameter (cm).
    void_fraction : float
        Inter-particle void fraction ε (0–1).  Typical: 0.35–0.45.
        Liquid occupies ε × V_column; solid occupies (1−ε) × V_column.
    bulk_density_g_mL : float
        Packing bulk density (g solid / mL bed volume).
        Related to void fraction and particle density:
        ρ_bulk = (1 − ε) × ρ_particle.  Typical: 0.5–0.8 g/mL.
    flow_rate_mL_min : float
        Volumetric flow rate through the column (mL/min).
    n_cells : int
        Number of finite-volume cells.  More cells → more accurate
        but slower.  50–200 is typical.
    species : list of LangmuirSpecies
        Solute species with their Langmuir parameters.
    D_ax_default : float
        Default axial dispersion coefficient (cm²/min).
        Applied to species where ``D_ax`` is not set.
        Set to 0.0 to disable dispersion (pure plug flow + adsorption).
    T_K : float
        Column temperature (K).  Default 298.15 (25°C).
    label : str
        Human-readable label.
    """

    _FROZEN_ATTRS = frozenset({
        "length_cm", "diameter_cm", "void_fraction", "flow_rate_mL_min",
        "bulk_density_g_mL", "D_ax_default", "T_K", "n_cells",
        "mobile_phase_pH", "K_D_reference_pH", "advection_scheme",
        "dispersion_scheme", "intraparticle_porosity", "species",
    })

    def __setattr__(self, name, value):
        if getattr(self, "_frozen", False) and name in self._FROZEN_ATTRS:
            raise AttributeError(
                f"Cannot modify '{name}' after construction. "
                f"Create a new HPLCColumn with the desired parameters.")
        super().__setattr__(name, value)

    def __init__(
        self,
        length_cm: float = 25.0,
        diameter_cm: float = 0.46,
        void_fraction: float = 0.40,
        bulk_density_g_mL: float = 0.60,
        flow_rate_mL_min: float = 1.0,
        n_cells: Optional[int] = None,
        species: Optional[List] = None,
        D_ax_default: float = 0.001,
        intraparticle_porosity: float = 0.50,
        T_K: float = 298.15,
        label: str = "HPLC_column",
        target_plates: Optional[int] = None,
        advection_scheme: str = "upwind",
        dispersion_scheme: str = "central_2nd",
        mobile_phase_pH: Optional[float] = None,
        K_D_reference_pH: Optional[float] = None,
    ):
        self.length_cm = float(length_cm)
        self.diameter_cm = float(diameter_cm)
        self.void_fraction = float(void_fraction)
        self.bulk_density_g_mL = float(bulk_density_g_mL)
        self.flow_rate_mL_min = float(flow_rate_mL_min)
        self.species = tuple(species or ())
        self.D_ax_default = float(D_ax_default)
        self.intraparticle_porosity = float(intraparticle_porosity)
        self.T_K = float(T_K)
        self.label = str(label)
        self.mobile_phase_pH = float(mobile_phase_pH) if mobile_phase_pH is not None else None

        # ── Physical parameter validation ─────────────────────────
        if self.length_cm <= 0:
            raise ValueError(f"length_cm must be positive, got {self.length_cm}")
        if self.diameter_cm <= 0:
            raise ValueError(f"diameter_cm must be positive, got {self.diameter_cm}")
        if not (0 < self.void_fraction < 1):
            raise ValueError(
                f"void_fraction must be in (0, 1), got {self.void_fraction}")
        if self.flow_rate_mL_min <= 0:
            raise ValueError(
                f"flow_rate_mL_min must be positive, got {self.flow_rate_mL_min}")
        if self.bulk_density_g_mL < 0:
            raise ValueError(
                f"bulk_density_g_mL must be non-negative, got {self.bulk_density_g_mL}")
        if self.D_ax_default < 0:
            raise ValueError(
                f"D_ax_default must be non-negative, got {self.D_ax_default}")
        if self.T_K <= 0:
            raise ValueError(f"T_K must be positive, got {self.T_K}")
        if not self.species:
            raise ValueError("At least one species must be provided.")

        # ── Species-level validation ──────────────────────────────
        for sp in self.species:
            if isinstance(sp, PartitionSpecies):
                if sp.K_D < 0:
                    raise ValueError(
                        f"K_D must be non-negative for '{sp.name}', got {sp.K_D}")
            elif isinstance(sp, LangmuirSpecies):
                if sp.q_max <= 0:
                    raise ValueError(
                        f"q_max must be positive for '{sp.name}', got {sp.q_max}")
                if sp.K_L < 0:
                    raise ValueError(
                        f"K_L must be non-negative for '{sp.name}', got {sp.K_L}")
            if sp.kf < 0:
                raise ValueError(
                    f"kf must be non-negative for '{sp.name}', got {sp.kf}")

        # ── Soft warnings for unusual configurations ──────────────
        if self.void_fraction > 0.95:
            warnings.warn(
                f"void_fraction={self.void_fraction:.3f} is unusually high "
                f"(typical range 0.25–0.50 for packed HPLC columns). "
                f"Phase ratio (1-ε)/ε = "
                f"{(1 - self.void_fraction) / self.void_fraction:.4f} "
                f"may give negligible retention.",
                UserWarning, stacklevel=2,
            )

        # K_D_reference_pH: the pH at which the supplied K_D values
        # were experimentally determined.  If set, K_D values for
        # species with pKa are back-corrected to K_D_neutral (the
        # intrinsic partition coefficient of the fully protonated
        # form), then re-corrected to the operating mobile_phase_pH.
        # If not set but mobile_phase_pH is set, assumes K_D values
        # were measured at mobile_phase_pH (so no net correction at
        # that pH, but corrections apply at other pH values for
        # Phase B future use).
        if K_D_reference_pH is not None:
            self.K_D_reference_pH = float(K_D_reference_pH)
        elif mobile_phase_pH is not None:
            self.K_D_reference_pH = self.mobile_phase_pH
            # Check if any species have pKa — if so, the pH correction
            # will silently cancel out at this pH.
            has_pKa = any(
                isinstance(sp, PartitionSpecies) and sp.pKa is not None
                for sp in self.species
            )
            if has_pKa:
                warnings.warn(
                    f"mobile_phase_pH={mobile_phase_pH} set without explicit "
                    f"K_D_reference_pH. Assuming K_D values were calibrated at "
                    f"pH {mobile_phase_pH} — pH correction will have no effect "
                    f"at this pH. Set K_D_reference_pH to the pH at which K_D "
                    f"values were experimentally determined to enable "
                    f"pH-dependent retention shifts.",
                    UserWarning, stacklevel=2,
                )
        else:
            self.K_D_reference_pH = None

        # Validate advection scheme
        advection_scheme = advection_scheme.lower().strip()
        if advection_scheme not in available_advection_schemes():
            raise ValueError(
                f"Unknown advection scheme {advection_scheme!r}. "
                f"Available: {available_advection_schemes()}"
            )
        self.advection_scheme = advection_scheme

        # Validate dispersion scheme
        dispersion_scheme = dispersion_scheme.lower().strip()
        if dispersion_scheme not in available_dispersion_schemes():
            raise ValueError(
                f"Unknown dispersion scheme {dispersion_scheme!r}. "
                f"Available: {available_dispersion_schemes()}"
            )
        self.dispersion_scheme = dispersion_scheme

        # Derived geometry (needed for n_cells calculation)
        self._A_cm2 = math.pi / 4 * self.diameter_cm ** 2
        self._V_col_mL = self._A_cm2 * self.length_cm
        _u = self.flow_rate_mL_min / (self._A_cm2 * self.void_fraction)

        # Determine n_cells: either from target_plates or directly.
        # The relationship between n_cells and plate count depends on
        # the advection scheme:
        #   upwind:     N_plates ≈ L×u / (2×(D_ax + u×Δz/2))
        #   tvd_vanleer: N_plates ≈ L×u / (2×D_ax)  (numerical diffusion negligible)
        if target_plates is not None and n_cells is not None:
            raise ValueError(
                "Specify either n_cells or target_plates, not both."
            )
        if target_plates is not None:
            if advection_scheme == "upwind":
                HETP_target = self.length_cm / float(target_plates)
                D_total_needed = HETP_target * _u / 2.0
                D_num_needed = D_total_needed - self.D_ax_default
                if D_num_needed <= 0:
                    self.n_cells = max(1000, int(self.length_cm / 0.01))
                else:
                    dz_needed = 2.0 * D_num_needed / _u
                    self.n_cells = max(10, int(math.ceil(self.length_cm / dz_needed)))
            else:
                # TVD and higher-order schemes: numerical diffusion is
                # negligible, so plate count is set by physical D_ax.
                # We just need enough cells to resolve the peaks spatially.
                # Rule of thumb: ~20 cells per peak standard deviation.
                # σ_z = L / sqrt(N_plates), so cells_per_sigma = L / (N × σ_z)
                N_target = float(target_plates)
                sigma_z = self.length_cm / math.sqrt(max(N_target, 1))
                self.n_cells = max(50, int(math.ceil(
                    20 * self.length_cm / sigma_z
                )))
                # Check: can the physical D_ax actually deliver the
                # requested plates?
                N_physical = self.length_cm * _u / (2 * max(self.D_ax_default, 1e-30))
                if N_physical < N_target:
                    warnings.warn(
                        f"target_plates={int(N_target)} exceeds the physical "
                        f"plate count of {N_physical:.0f} set by D_ax="
                        f"{self.D_ax_default:.4f} cm²/min. The TVD scheme "
                        f"cannot exceed the physical limit. Reduce D_ax or "
                        f"increase column length to achieve higher plate counts.",
                        UserWarning, stacklevel=2,
                    )
            self._target_plates = int(target_plates)
        elif n_cells is not None:
            self.n_cells = int(n_cells)
            self._target_plates = None
        else:
            self.n_cells = 100
            self._target_plates = None

        # Derived geometry (cell-level)
        self._dz = self.length_cm / self.n_cells
        self._V_cell_mL = self._A_cm2 * self._dz
        self._V_liq_cell_mL = self._V_cell_mL * self.void_fraction
        self._V_solid_cell_mL = self._V_cell_mL * (1 - self.void_fraction)
        self._m_solid_cell_g = self._V_cell_mL * self.bulk_density_g_mL
        self._u_cm_min = _u

        # Species indexing
        self._sp_names = [sp.name for sp in self.species]
        self._n_sp = len(self.species)
        self._sp_idx = {sp.name: i for i, sp in enumerate(self.species)}

        if len(self._sp_idx) != self._n_sp:
            dupes = [n for n in self._sp_names if self._sp_names.count(n) > 1]
            raise ValueError(
                f"Duplicate species names: {sorted(set(dupes))}. "
                f"Each species in the column must have a unique name.")

        # ── pH-dependent partition: back-calculate K_D_neutral ────
        # The user-supplied K_D values were measured at K_D_reference_pH.
        # We recover K_D_neutral (intrinsic, fully protonated form):
        #
        #   K_D_neutral = K_D_observed / f_HA(reference_pH, pKa)
        #
        # Then in the ODE, K_D_eff = K_D_neutral × f_HA(mobile_phase_pH, pKa).
        # When mobile_phase_pH == reference_pH, K_D_eff == K_D_observed.
        # When mobile_phase_pH > reference_pH, acids are more dissociated,
        # K_D_eff drops, and they elute earlier.
        self._K_D_neutral = {}
        if self.K_D_reference_pH is not None:
            for sp in self.species:
                if isinstance(sp, PartitionSpecies):
                    if sp.pKa is not None:
                        f_HA_ref = _f_protonated(self.K_D_reference_pH, sp.pKa)
                        self._K_D_neutral[sp.name] = sp.K_D / max(f_HA_ref, 1e-30)
                    else:
                        self._K_D_neutral[sp.name] = sp.K_D

        self._frozen = True

    # ── Properties ───────────────────────────────────────────────────

    @property
    def column_volume_mL(self) -> float:
        return self._V_col_mL

    @property
    def void_volume_mL(self) -> float:
        return self._V_col_mL * self.void_fraction

    @property
    def dead_time_min(self) -> float:
        """Time for an unretained species to traverse the column."""
        return self.void_volume_mL / self.flow_rate_mL_min

    @property
    def interstitial_velocity_cm_min(self) -> float:
        return self._u_cm_min

    @property
    def cell_length_cm(self) -> float:
        return self._dz

    def predicted_retention_min(self, species_name: str) -> float:
        """Predicted retention time for a species based on column parameters.

        For partition species::

            t_R = (ε + (1-ε) × ε_p × K_D_eff) × V_col / Q

        For Langmuir species (linear limit, low concentration)::

            t_R = (ε + (1-ε) × ρ_b × q_max × K_L) × V_col / Q

        Parameters
        ----------
        species_name : str
            Name of the species (must be in the column's species list).

        Returns
        -------
        float
            Predicted retention time in minutes.
        """
        if species_name not in self._sp_idx:
            raise ValueError(
                f"Species '{species_name}' not in column. "
                f"Available: {self._sp_names}")
        idx = self._sp_idx[species_name]
        sp = self.species[idx]
        eps = self.void_fraction
        V_col = self._V_col_mL
        Q = self.flow_rate_mL_min

        if isinstance(sp, PartitionSpecies):
            ep = (sp.intraparticle_porosity
                  if sp.intraparticle_porosity is not None
                  else self.intraparticle_porosity)
            # Use effective K_D if pH correction is active
            if (self.mobile_phase_pH is not None
                    and sp.name in self._K_D_neutral):
                K_D = (self._K_D_neutral[sp.name]
                       * (_f_protonated(self.mobile_phase_pH, sp.pKa)
                          if sp.pKa is not None else 1.0))
            else:
                K_D = sp.K_D
            return (eps + (1 - eps) * ep * K_D) * V_col / Q
        elif isinstance(sp, LangmuirSpecies):
            # Linear limit: k' = (1-ε) × ρ_b × q_max × K_L / ε
            rho_b = self.bulk_density_g_mL
            return (eps + (1 - eps) * rho_b * sp.q_max * sp.K_L) * V_col / Q
        else:
            # Fallback: unretained
            return self.dead_time_min

    def predicted_last_elution_min(self) -> float:
        """Predicted retention time of the most retained species.

        Useful as a default ``recycle_after_min`` value — switching to
        recycle mode after the slowest species has had one pass ensures
        all injected material is within the column.

        Returns
        -------
        float
            Maximum predicted retention time across all species.
        """
        return max(self.predicted_retention_min(sp.name)
                   for sp in self.species)

    @property
    def peclet_number(self) -> float:
        """Column Peclet number: Pe = u × L / D_ax."""
        D_ax = max(self.D_ax_default, 1e-30)
        return self._u_cm_min * self.length_cm / D_ax

    @property
    def effective_plates(self) -> float:
        """Effective number of theoretical plates.

        For the ``upwind`` scheme, includes numerical dispersion:

            D_total = D_ax + u × Δz / 2
            HETP = 2 × D_total / u
            N_plates = L / HETP

        For ``tvd_vanleer`` and other higher-order schemes, numerical
        dispersion is negligible and plates are set by physical D_ax:

            N_plates ≈ L × u / (2 × D_ax)
        """
        if self.advection_scheme == "upwind":
            D_num = self._u_cm_min * self._dz / 2.0
            D_total = self.D_ax_default + D_num
        else:
            # TVD and higher-order: numerical diffusion negligible
            D_total = max(self.D_ax_default, 1e-30)
        HETP = 2.0 * D_total / self._u_cm_min
        return self.length_cm / HETP

    def summary(self) -> str:
        lines = [
            f"HPLCColumn: {self.label}",
            f"  Geometry:   L={self.length_cm:.1f} cm, "
            f"d={self.diameter_cm:.2f} cm, "
            f"V_col={self._V_col_mL:.3f} mL",
            f"  Packing:    ε={self.void_fraction:.2f}, "
            f"ρ_bulk={self.bulk_density_g_mL:.2f} g/mL",
            f"  Flow:       Q={self.flow_rate_mL_min:.3f} mL/min, "
            f"u={self._u_cm_min:.3f} cm/min",
            f"  Void time:  t₀={self.dead_time_min:.2f} min",
            f"  Dispersion: D_ax={self.D_ax_default:.4f} cm²/min, "
            f"Pe={self.peclet_number:.0f}, scheme={self.dispersion_scheme}",
            f"  Advection:  {self.advection_scheme}",
            f"  Grid:       {self.n_cells} cells, "
            f"Δz={self._dz:.4f} cm",
            f"  Plates:     N_eff={self.effective_plates:.0f}"
            + (f" (target={self._target_plates})" if self._target_plates else ""),
        ]
        if self.mobile_phase_pH is not None:
            lines.append(f"  Mobile pH:  {self.mobile_phase_pH:.2f} "
                         f"(Donnan correction active)")
        lines.append(f"  Species ({self._n_sp}):")
        for sp in self.species:
            if isinstance(sp, LangmuirSpecies):
                D = sp.D_ax if sp.D_ax is not None else self.D_ax_default
                lines.append(
                    f"    {sp.name:15s} [Langmuir] q_max={sp.q_max:.4f} mol/g, "
                    f"K_L={sp.K_L:.2f} L/mol, kf={sp.kf:.1f} /min"
                )
            elif isinstance(sp, PartitionSpecies):
                ep = sp.intraparticle_porosity if sp.intraparticle_porosity is not None else self.intraparticle_porosity
                extra = f", ε_p={ep:.2f}, kf={sp.kf:.1f} /min"
                if sp.pKa is not None:
                    extra += f", pKa={sp.pKa:.2f}"
                    if sp.name in self._K_D_neutral:
                        extra += f", K_D_neutral={self._K_D_neutral[sp.name]:.3f}"
                lines.append(
                    f"    {sp.name:15s} [Partition] K_D={sp.K_D:.3f}" + extra
                )
        return "\n".join(lines)

    # ── State vector packing ─────────────────────────────────────────

    def _state_size(self) -> int:
        """Total ODE state variables: 2 × n_species × n_cells."""
        return 2 * self._n_sp * self.n_cells

    def _pack_y0(
        self,
        C_init: Optional[Dict[str, float]] = None,
        q_init: Optional[Dict[str, float]] = None,
    ) -> np.ndarray:
        """Pack initial conditions into a flat state vector.

        Layout: [C_s0_cell0, ..., C_sN_cell0, q_s0_cell0, ..., q_sN_cell0,
                 C_s0_cell1, ..., q_sN_cellN-1]
        i.e., for each cell: all C values then all q values.
        """
        y0 = np.zeros(self._state_size())
        C_init = C_init or {}
        q_init = q_init or {}

        sp_name_set = set(self._sp_names)
        unknown_C = set(C_init.keys()) - sp_name_set
        if unknown_C:
            warnings.warn(
                f"C_init contains species not in the column: {sorted(unknown_C)}. "
                f"These will be ignored.",
                UserWarning, stacklevel=2,
            )
        unknown_q = set(q_init.keys()) - sp_name_set
        if unknown_q:
            warnings.warn(
                f"q_init contains species not in the column: {sorted(unknown_q)}. "
                f"These will be ignored.",
                UserWarning, stacklevel=2,
            )

        for cell in range(self.n_cells):
            base = cell * 2 * self._n_sp
            for s, sp in enumerate(self.species):
                y0[base + s] = C_init.get(sp.name, 0.0)
                y0[base + self._n_sp + s] = q_init.get(sp.name, 0.0)
        return y0

    def _unpack(self, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Unpack flat state vector into C and q arrays.

        Returns
        -------
        C : ndarray, shape (n_cells, n_species) — mol/L
        q : ndarray, shape (n_cells, n_species) — mol/g
        """
        C = np.zeros((self.n_cells, self._n_sp))
        q = np.zeros((self.n_cells, self._n_sp))
        for cell in range(self.n_cells):
            base = cell * 2 * self._n_sp
            C[cell, :] = y[base:base + self._n_sp]
            q[cell, :] = y[base + self._n_sp:base + 2 * self._n_sp]
        return C, q

    # ── ODE right-hand side ──────────────────────────────────────────

    def _build_ode(self, inlet_fn: Callable, recycle_after_min: Optional[float] = None,
                   divert_fn: Optional[Callable] = None):
        """Build the ODE function for solve_ivp.

        Handles two isotherm types simultaneously:

        - **LangmuirSpecies**: competitive Langmuir adsorption.
          q is mol/g_solid.  Phase ratio: (1−ε) × ρ_b / ε.

        - **PartitionSpecies**: linear partition into stagnant pore
          liquid (ion-exclusion / size-exclusion).  q is mol/L_pore.
          Phase ratio: (1−ε) × ε_p / ε.  No competition.

        Fully vectorised over cells using numpy.

        Parameters
        ----------
        inlet_fn : callable
            ``inlet_fn(t_min) -> {species: C_mol_L}``
        recycle_after_min : float or None
            If set, after this time the column outlet (last cell)
            is fed back as the column inlet, creating a closed-loop
            recycle system.  Before this time, ``inlet_fn`` is used
            normally.
        divert_fn : callable or None
            If set, applied during recycle to modify the outlet before
            it is fed back to the inlet.  Signature:
            ``divert_fn(t_min, C_outlet_array) -> C_inlet_array``.
            Returning zeros diverts all material to collection;
            returning C_outlet unchanged gives full recycle.
            Only called when recycle is active (t >= recycle_after_min).
            Use :func:`make_divert_windows`, :func:`make_divert_threshold`,
            or :func:`make_peak_shave_schedule` to create common patterns.
        """
        N = self.n_cells
        n_sp = self._n_sp
        eps = self.void_fraction
        rho_b = self.bulk_density_g_mL
        eps_p_default = self.intraparticle_porosity
        dz = self._dz

        # Classify species
        is_langmuir = np.array([isinstance(sp, LangmuirSpecies)
                                for sp in self.species])
        lang_idx = np.where(is_langmuir)[0]
        part_idx = np.where(~is_langmuir)[0]
        n_lang = len(lang_idx)

        # Common per-species arrays
        kf = np.array([sp.kf for sp in self.species])
        D_ax = np.array([
            sp.D_ax if sp.D_ax is not None else self.D_ax_default
            for sp in self.species
        ])

        # Langmuir parameters (indexed into lang_idx positions)
        if n_lang > 0:
            q_max_L = np.array([self.species[i].q_max for i in lang_idx])
            K_L = np.array([self.species[i].K_L for i in lang_idx])
        else:
            q_max_L = K_L = np.array([])

        # Partition parameters (indexed into part_idx positions)
        # When mobile_phase_pH is set and species have pKa values,
        # K_D_eff = K_D_neutral × f_HA(pH, pKa) accounts for Donnan
        # exclusion of the dissociated (anionic) form.
        if len(part_idx) > 0:
            if self.mobile_phase_pH is not None and self._K_D_neutral:
                K_D_eff = np.array([
                    self._K_D_neutral[self.species[i].name]
                    * (_f_protonated(self.mobile_phase_pH, self.species[i].pKa)
                       if self.species[i].pKa is not None else 1.0)
                    for i in part_idx
                ])
            else:
                K_D_eff = np.array([self.species[i].K_D for i in part_idx])
        else:
            K_D_eff = np.array([])

        # Per-species phase ratio
        phase_ratio = np.zeros(n_sp)
        for i, sp in enumerate(self.species):
            if isinstance(sp, LangmuirSpecies):
                phase_ratio[i] = (1 - eps) * rho_b / eps
            else:
                ep = (sp.intraparticle_porosity
                      if sp.intraparticle_porosity is not None
                      else eps_p_default)
                phase_ratio[i] = (1 - eps) * ep / eps

        # Transport coefficients
        u = self.flow_rate_mL_min / (self._A_cm2 * eps)

        # Build the advection callable from the selected scheme
        advect = get_advection_scheme(self.advection_scheme, u, dz, N, n_sp)

        # Build the dispersion callable from the selected scheme
        disperse = get_dispersion_scheme(self.dispersion_scheme, D_ax, dz, N, n_sp)

        sp_names = self._sp_names
        stride = 2 * n_sp
        _recycle_t = recycle_after_min  # capture in closure
        _divert = divert_fn            # capture in closure

        def f(t, y):
            Y = y.reshape(N, stride)
            C = np.maximum(Y[:, :n_sp], 0.0)
            q = np.maximum(Y[:, n_sp:], 0.0)

            # ── Equilibrium for each species type ────────────────
            q_eq = np.zeros_like(C)

            if n_lang > 0:
                C_L = C[:, lang_idx]
                denom = 1.0 + C_L @ K_L
                q_eq[:, lang_idx] = (
                    q_max_L[None, :] * K_L[None, :] * C_L
                    / denom[:, None]
                )

            if len(part_idx) > 0:
                q_eq[:, part_idx] = K_D_eff[None, :] * C[:, part_idx]

            # ── Mass transfer rate ───────────────────────────────
            dqdt = kf[None, :] * (q_eq - q)

            # ── Inlet concentration ──────────────────────────────
            if _recycle_t is not None and t >= _recycle_t:
                # Recycle: outlet of last cell feeds back to inlet
                C_outlet = C[-1, :]  # shape (n_sp,)
                if _divert is not None:
                    # Peak shaving: divert_fn decides what gets recycled
                    C_recycled = _divert(t, C_outlet)
                    C_in = np.maximum(C_recycled, 0.0).reshape(1, -1)
                else:
                    C_in = C_outlet.reshape(1, -1)
            else:
                C_in_dict = inlet_fn(t)
                C_in = np.array([[C_in_dict.get(sp, 0.0) for sp in sp_names]])
            adv = advect(C, C_in)

            # ── Axial dispersion (delegated to scheme) ───────────
            disp = disperse(C, C_in)

            # ── Assemble dC/dt (per-species phase ratio) ─────────
            dCdt = adv + disp - phase_ratio[None, :] * dqdt

            dY = np.empty_like(Y)
            dY[:, :n_sp] = dCdt
            dY[:, n_sp:] = dqdt
            return dY.ravel()

        return f

    # ── Simulation ───────────────────────────────────────────────────

    def simulate(
        self,
        t_end_min: float,
        inlet_fn: Callable,
        *,
        C_init: Optional[Dict[str, float]] = None,
        q_init: Optional[Dict[str, float]] = None,
        n_output: int = 500,
        method: str = "LSODA",
        rtol: float = 1e-6,
        atol: float = 1e-9,
        max_step: Optional[float] = None,
        recycle_after_min: Optional[float] = None,
        divert_fn: Optional[Callable] = None,
    ) -> HPLCResult:
        """Simulate the column and return the elution profile.

        Parameters
        ----------
        t_end_min : float
            Simulation end time (minutes).
        inlet_fn : callable
            ``inlet_fn(t_min) -> {species_name: C_mol_L}``.
            Returns inlet concentrations at time ``t``.
            Use :func:`pulse_inlet` or :func:`step_inlet` for common
            injection modes, or provide a custom function.
        C_init : dict, optional
            Initial liquid concentration in each cell (mol/L).
            Default: all zeros (clean column).
        q_init : dict, optional
            Initial adsorbed concentration (mol/g).
            Default: all zeros (clean column).
        n_output : int
            Number of output timepoints for the result.
        method : str
            ODE solver method.  ``"LSODA"`` (auto stiff/nonstiff) is
            recommended.  ``"RK45"`` also works well.  Avoid ``"Radau"``
            for large N — the implicit Jacobian is O(N²) per step.
        rtol, atol : float
            Solver tolerances.
        max_step : float, optional
            Maximum solver step size (minutes).  If None, auto-selected
            as the minimum of: cell transit time (CFL), 1% of simulation
            time, and the output interval.  For narrow pulse injections,
            set this to less than the pulse duration to ensure the solver
            doesn't step over the pulse.
        recycle_after_min : float, optional
            If set, switches to closed-loop recycle mode at this time.
            Before ``recycle_after_min``, the column inlet is fed by
            ``inlet_fn`` (normal operation with injection).  After
            ``recycle_after_min``, the column outlet (last cell, after
            the detector) is fed directly back to the inlet.  This
            models recycle chromatography where the eluent passes
            through the column multiple times for improved resolution.

            A typical value is :meth:`dead_time_min` (one void volume)
            or :meth:`predicted_retention_min` for the species of
            interest.
        divert_fn : callable, optional
            Peak-shaving function applied during recycle mode.
            Signature: ``divert_fn(t_min, C_outlet_array) -> C_inlet_array``.
            When the function returns zeros for some or all species,
            that material is diverted to collection rather than recycled.
            Only active when ``recycle_after_min`` is set and
            ``t >= recycle_after_min``.

            Use the convenience factories to build common patterns:

            - :func:`make_divert_windows` — divert during fixed time windows
            - :func:`make_divert_threshold` — divert when signal exceeds threshold
            - :func:`make_peak_shave_schedule` — auto-generate windows from
              predicted retention times (heart-cutting)

        Returns
        -------
        HPLCResult
        """
        from scipy.integrate import solve_ivp

        y0 = self._pack_y0(C_init, q_init)
        ode_fn = self._build_ode(inlet_fn, recycle_after_min=recycle_after_min,
                                 divert_fn=divert_fn)
        t_eval = np.linspace(0.0, t_end_min, n_output)

        # Auto max_step: ensure the solver resolves the cell transit
        # time and doesn't step over narrow injection pulses.
        if max_step is None:
            cfl_dt = self._dz / max(self._u_cm_min, 1e-30)
            output_dt = t_end_min / max(n_output, 1)
            max_step = min(cfl_dt, output_dt, t_end_min / 100)

            # Detect narrow inlet features by probing the inlet function
            # at fine resolution over the first 10% of the simulation.
            # If the inlet has sharp transitions (pulse injection), the
            # max_step must be smaller than the pulse duration.
            probe_dt = max_step / 10
            t_probe = np.arange(0, min(t_end_min * 0.1, 5.0), probe_dt)
            if len(t_probe) > 1:
                sp0 = self._sp_names[0] if self._sp_names else None
                if sp0:
                    vals = np.array([inlet_fn(t).get(sp0, 0.0) for t in t_probe])
                    nonzero = np.where(vals > 0)[0]
                    if len(nonzero) >= 2:
                        # Pulse width ≈ time span of nonzero inlet
                        pulse_width = t_probe[nonzero[-1]] - t_probe[nonzero[0]]
                        if pulse_width > 0 and pulse_width < max_step * 2:
                            max_step = min(max_step, pulse_width / 3)

        kwargs = dict(
            method=method, rtol=rtol, atol=atol,
            dense_output=False, t_eval=t_eval,
            max_step=max_step,
        )

        t0 = _time.time()
        sol = solve_ivp(ode_fn, [0.0, t_end_min], y0, **kwargs)
        runtime = _time.time() - t0

        if not sol.success:
            warnings.warn(
                f"HPLC solver failed: {sol.message}. "
                f"Results may be incomplete.",
                RuntimeWarning,
            )

        # Unpack results at each output time
        n_times = len(sol.t)
        C_all = {sp.name: np.zeros((n_times, self.n_cells))
                 for sp in self.species}
        q_all = {sp.name: np.zeros((n_times, self.n_cells))
                 for sp in self.species}
        C_outlet = {sp.name: np.zeros(n_times) for sp in self.species}

        for ti in range(n_times):
            C, q = self._unpack(sol.y[:, ti])
            for s, sp in enumerate(self.species):
                C_all[sp.name][ti, :] = C[:, s]
                q_all[sp.name][ti, :] = q[:, s]
                C_outlet[sp.name][ti] = max(0.0, C[-1, s])

        return HPLCResult(
            t_min=sol.t,
            C_outlet=C_outlet,
            C_all=C_all,
            q_all=q_all,
            column=self,
            runtime_s=runtime,
            solver_info={
                "success": sol.success,
                "message": sol.message,
                "n_eval": sol.nfev,
                "method": method,
            },
        )

    # ── Fermenter integration helpers ────────────────────────────────

    def build_fermenter_cells(self):
        """Create fermenter LiquidPhase + SolidPhase objects for each cell.

        Returns
        -------
        list of dict
            ``[{"liquid": LiquidPhase, "solid": SolidPhase}, ...]``
            One per cell.  These can be attached to ControlVolumes for
            integration with the fermenter's MultiCVSystem.
        """
        from PyOMES.core.phases import LiquidPhase, SolidPhase

        cells = []
        for i in range(self.n_cells):
            liq_n_mol = {sp.name: 0.0 for sp in self.species}
            sol_n_mol = {f"{sp.name}_ads": 0.0 for sp in self.species}

            liq = LiquidPhase(
                n_mol=liq_n_mol,
                V_L=self._V_liq_cell_mL * 1e-3,  # mL → L
                T_K=self.T_K,
            )
            sol = SolidPhase(
                n_mol=sol_n_mol,
                V_L=self._V_solid_cell_mL * 1e-3,  # mL → L
                T_K=self.T_K,
            )
            cells.append({"liquid": liq, "solid": sol})
        return cells

    def __repr__(self):
        return (f"HPLCColumn(L={self.length_cm} cm, d={self.diameter_cm} cm, "
                f"ε={self.void_fraction}, {self._n_sp} species, "
                f"{self.n_cells} cells)")
