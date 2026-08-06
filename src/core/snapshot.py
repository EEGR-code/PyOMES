# -*- coding: utf-8 -*-
"""Read-only snapshots of CV / Simulation state for controllers.

A :class:`CVSnapshot` is the typed view that the controller protocol
consumes (decision 3). It is a frozen dataclass populated once per
step by the orchestrator from the post-:meth:`ControlVolume.advance`
state.

The :func:`build_cv_snapshot` builder reads fields off the CV's
phases directly. Post-STATE_UNIFICATION the
:attr:`AdvanceResult.properties` field was removed; pH and
derived chemistry-state live on the phases via
``phase.pH`` (raises if no ``"H+"`` species) and
``liquid.properties`` (PropertyCalculator outputs) /
``liquid.speciation`` (legacy ``"IonicStrength"`` key until
``chemistry-unification-3b``).

This module lands C2 of SIMULATION_CLASS. The
:class:`Controller` protocol that consumes these snapshots lands at
C7; orchestrator wiring lands at C9.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


# ════════════════════════════════════════════════════════════════════════
#  Snapshot dataclasses
# ════════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class CVSnapshot:
    """Read-only view of one ControlVolume at one timestep.

    Constructed by :func:`build_cv_snapshot` from a CV's current
    state. Given to controllers and recorders without granting them
    write access. Frozen at the attribute level — attempting to
    reassign a field raises :class:`dataclasses.FrozenInstanceError`.

    Attributes
    ----------
    cv_key : str
        The CV's key within its owning :class:`Simulation`.
    t_h : float
        Wall-clock time at which the snapshot was taken (hours).
    pH : float or None
        Liquid pH if a speciation engine has populated
        ``n_mol["H+"]`` on the liquid phase; otherwise ``None``.
        Read via ``cv.phases["liquid"].pH`` with the
        ``ValueError`` caught for the no-speciation case.
    ionic_strength : float or None
        Liquid ionic strength (mol/L) if available. Prefers
        ``liquid.properties["ionic_strength"]`` (PropertyCalculator
        output, post-STATE_UNIFICATION C5); falls back to
        ``liquid.speciation["IonicStrength"]`` (legacy until
        chemistry-unification-3b).
    T_K : float
        Liquid temperature if a liquid phase is present, else gas
        temperature, else ``298.15`` as a final default.
    V_liq_L : float
        Liquid volume (L); ``0.0`` if no liquid phase.
    V_gas_L : float
        Gas volume (L); ``0.0`` if no gas phase.
    P_gas_atm : float
        Total gas pressure (atm) from ideal-gas law on the gas
        phase; ``0.0`` if no gas phase.
    n_gas_mol : dict
        Per-species gas moles ``{species_id: mol}``.
    n_liq_mol : dict
        Per-species liquid moles ``{species_id: mol}``.
    y_gas : dict
        Per-species gas mole fractions.
    species_properties : dict
        Liquid phase ``properties`` dict (PropertyCalculator
        outputs and any forwarded scalar derived properties).
        Alias for ``properties.get("liquid", {})``.
    sensors : dict
        Derived quantities convenient for controllers
        (``DO_mol_L``, ``T_C``, etc.). Extension point; the exact
        membership grows as controllers are ported in C8.
    phase_keys : list
        Sorted list of phase keys present on the CV (e.g.
        ``["gas", "liquid"]`` for a fermenter, or
        ``["mobile", "stationary"]`` for HPLC).
    n_mol : dict
        ``{phase_key: {species_id: mol}}`` for every phase.
        General form of ``n_liq_mol`` / ``n_gas_mol``.
    V_L : dict
        ``{phase_key: volume_L}`` for every phase.
        General form of ``V_liq_L`` / ``V_gas_L``.
    properties : dict
        ``{phase_key: {property_key: value}}`` for every phase.
        General form of ``species_properties``.
    """

    cv_key: str
    t_h: float
    pH: Optional[float]
    ionic_strength: Optional[float]
    T_K: float
    V_liq_L: float
    V_gas_L: float
    P_gas_atm: float
    n_gas_mol: Dict[str, float] = field(default_factory=dict)
    n_liq_mol: Dict[str, float] = field(default_factory=dict)
    y_gas: Dict[str, float] = field(default_factory=dict)
    species_properties: Dict[str, float] = field(default_factory=dict)
    sensors: Dict[str, Any] = field(default_factory=dict)
    # General phase-agnostic storage — keyed by phase_key.
    # New domains (HPLC, membrane, cell culture) write controllers
    # against these instead of the liquid/gas convenience aliases above.
    phase_keys: List[str] = field(default_factory=list)
    n_mol: Dict[str, Dict[str, float]] = field(default_factory=dict)
    V_L: Dict[str, float] = field(default_factory=dict)
    properties: Dict[str, Dict[str, float]] = field(default_factory=dict)


@dataclass(frozen=True)
class SimulationSnapshot:
    """Multi-CV snapshot across an entire Simulation at one timestep.

    Given to multi-CV-aware controllers (rare). Single-CV controllers
    receive a :class:`CVSnapshot` directly.

    Attributes
    ----------
    t_h : float
        Wall-clock time at which the snapshot was taken (hours).
    cvs : dict
        ``{cv_key: CVSnapshot}`` for every CV in the Simulation.
    """

    t_h: float
    cvs: Dict[str, CVSnapshot] = field(default_factory=dict)


# ════════════════════════════════════════════════════════════════════════
#  Builders
# ════════════════════════════════════════════════════════════════════════


def build_cv_snapshot(
    cv: Any,
    cv_key: str,
    t_h: float,
    advance_result: Any = None,
) -> CVSnapshot:
    """Build a :class:`CVSnapshot` from a CV's current state.

    Parameters
    ----------
    cv : ControlVolume
        The CV to snapshot. Phases keyed ``"gas"`` and ``"liquid"``
        are pulled when present; missing phases yield sensible
        defaults.
    cv_key : str
        The CV's key in the owning Simulation (echoed onto the
        snapshot).
    t_h : float
        Wall-clock time at which the snapshot is being taken.
    advance_result : AdvanceResult, optional
        Result of the CV's last :meth:`advance` call. Today unused
        by the snapshot builder — pH and derived chemistry-state
        live on phases directly post-STATE_UNIFICATION. Kept on
        the signature as a forward-compat hook for fields that
        controllers might want later (e.g. ``transfer_record``
        diagnostics surfaced through ``sensors``).

    Returns
    -------
    CVSnapshot
        A frozen view of the CV's current state.

    Notes
    -----
    The builder consults phase attributes directly:

    * ``cv.phases["liquid"].pH`` (catches the no-``H+`` ``ValueError``
      and surfaces ``None``)
    * ``cv.phases["liquid"].properties["ionic_strength"]`` (preferred,
      from ``PropertyCalculator``) with fallback to
      ``cv.phases["liquid"].speciation["IonicStrength"]`` (legacy)
    * ``cv.phases["liquid"].properties`` forwarded onto
      ``species_properties``
    * ``cv.phases["gas"].P_atm`` / ``cv.phases["gas"].y`` /
      ``cv.phases["gas"].n_mol``
    """
    liquid = cv.phases.get("liquid")
    gas = cv.phases.get("gas")

    # Temperature: prefer liquid, fall back to gas, then default.
    if liquid is not None:
        T_K = float(liquid.T_K)
    elif gas is not None:
        T_K = float(gas.T_K)
    else:
        T_K = 298.15

    # General phase-agnostic dicts — cover every phase the CV exposes.
    phase_keys: List[str] = sorted(cv.phases.keys())
    n_mol: Dict[str, Dict[str, float]] = {
        pk: dict(ph.n_mol) for pk, ph in cv.phases.items()
    }
    V_L_map: Dict[str, float] = {
        pk: float(ph.V_L) for pk, ph in cv.phases.items()
    }
    properties_map: Dict[str, Dict[str, float]] = {
        pk: dict(ph.properties)
        for pk, ph in cv.phases.items()
        if hasattr(ph, "properties")
    }

    # Liquid convenience aliases.
    if liquid is not None:
        V_liq_L = float(liquid.V_L)
        n_liq_mol: Dict[str, float] = dict(liquid.n_mol)
        # pH: surface None when H+ absent or computation fails.
        if "H+" in liquid.n_mol:
            try:
                pH: Optional[float] = float(liquid.pH)
            except ValueError:
                pH = None
        else:
            pH = None
        # Ionic strength: PropertyCalculator output wins; fall back to
        # the legacy speciation dict key used until chemistry-unification-3b.
        is_val = liquid.properties.get("ionic_strength")
        if is_val is None:
            is_val = liquid.speciation.get("IonicStrength")
        ionic_strength: Optional[float] = (
            float(is_val) if is_val is not None else None
        )
        species_properties: Dict[str, float] = dict(liquid.properties)
    else:
        V_liq_L = 0.0
        n_liq_mol = {}
        pH = None
        ionic_strength = None
        species_properties = {}

    # Gas convenience aliases.
    if gas is not None:
        V_gas_L = float(gas.V_L)
        n_gas_mol: Dict[str, float] = dict(gas.n_mol)
        y_gas: Dict[str, float] = dict(gas.y)
        P_gas_atm = float(gas.P_atm)
    else:
        V_gas_L = 0.0
        n_gas_mol = {}
        y_gas = {}
        P_gas_atm = 0.0

    # Sensors: derived quantities convenient for controllers.
    sensors: Dict[str, Any] = {"T_C": T_K - 273.15}
    if V_liq_L > 0.0 and "O2" in n_liq_mol:
        sensors["DO_mol_L"] = n_liq_mol["O2"] / V_liq_L

    return CVSnapshot(
        cv_key=str(cv_key),
        t_h=float(t_h),
        pH=pH,
        ionic_strength=ionic_strength,
        T_K=T_K,
        V_liq_L=V_liq_L,
        V_gas_L=V_gas_L,
        P_gas_atm=P_gas_atm,
        n_gas_mol=n_gas_mol,
        n_liq_mol=n_liq_mol,
        y_gas=y_gas,
        species_properties=species_properties,
        sensors=sensors,
        phase_keys=phase_keys,
        n_mol=n_mol,
        V_L=V_L_map,
        properties=properties_map,
    )


def build_simulation_snapshot(
    sim: Any,
    t_h: float,
    results: Optional[Dict[str, Any]] = None,
) -> SimulationSnapshot:
    """Build a :class:`SimulationSnapshot` aggregating every CV.

    Parameters
    ----------
    sim : Simulation
        The Simulation whose CVs to snapshot. Only ``sim.cvs`` is
        consulted (duck-typed to avoid a circular import).
    t_h : float
        Wall-clock time at which the snapshot is being taken.
    results : dict, optional
        ``{cv_key: AdvanceResult}`` from the most recent step;
        forwarded into each per-CV builder. ``None`` is equivalent
        to an empty dict.

    Returns
    -------
    SimulationSnapshot
    """
    results = dict(results or {})
    return SimulationSnapshot(
        t_h=float(t_h),
        cvs={
            key: build_cv_snapshot(cv, key, t_h, results.get(key))
            for key, cv in sim.cvs.items()
        },
    )
