# -*- coding: utf-8 -*-
"""Solver protocol and implementations for ControlVolume time-stepping.

A :class:`StepSolver` encapsulates the integration strategy used by
:meth:`ControlVolume.advance` when called with ``solver=...``.  When
``solver`` is omitted, ``advance()`` runs its own inline sequential
operator-split body — that default is not one of the classes below.

:class:`SimultaneousEulerSolver` uses a snapshot-based compute-then-apply
pattern with proportional clamping (explicit Euler, unsplit/simultaneous
composition).

For higher-order adaptive integration, see :class:`SimultaneousAdaptiveSolver`
(same unsplit/simultaneous composition, adaptive continuous integration
instead of a single Euler batch).

Usage
-----
>>> from PyOMES.core.solvers import SimultaneousEulerSolver
>>> solver = SimultaneousEulerSolver()
>>> result = cv.advance(dt_h=0.01, t_h=t_h, solver=solver)
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Protocol, runtime_checkable

from .phases import Phase, R_L_ATM_MOL_K
from .interfaces import AdvanceResult
from .links import LinkFlowRecord
from .boundaries import ExternalFluxRecord
from .clamping import proportional_clamp, floor_clamp, floor_nonnegative


# ════════════════════════════════════════════════════════════════════════
#  Protocol
# ════════════════════════════════════════════════════════════════════════

@runtime_checkable
class StepSolver(Protocol):
    """Protocol for time-stepping strategies.

    Implementations read the current CV state, compute all sub-system
    contributions, and update the CV in-place.
    """

    def solve_step(
        self,
        cv: "ControlVolume",
        dt_h: float,
        t_h: float,
        external_source_terms: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> "AdvanceResult":
        """Advance the CV by one timestep.

        Parameters
        ----------
        cv : ControlVolume
            The control volume to advance (mutated in-place).
        dt_h : float
            Timestep duration (hours).
        t_h : float
            Current simulation time (hours). Forwarded to the
            reaction environment for time-dependent rate laws.
            State-unification C4 replaced the legacy ``chem_env``
            dict with this single argument; strong ions now live in
            ``phase.n_mol`` and totals are derived inline by the
            speciation engine.
        external_source_terms : dict or None
            Additional source terms ``{phase_key: {species: mol/h}}``.

        Returns
        -------
        AdvanceResult
        """
        ...


# ════════════════════════════════════════════════════════════════════════
#  SequentialAdvanceSolver
# ════════════════════════════════════════════════════════════════════════

class SequentialAdvanceSolver:
    """The operator-split default, promoted to a real :class:`StepSolver`.

    Runs the same Lie-Trotter operator-split body ``ControlVolume.advance()``
    always ran inline: speciation solve, property calculators, external
    source terms, boundary fluxes, reactions (against the speciation
    pinned at step 1), internal transfer, conservation check. Moved here
    verbatim so the default is reachable, wrappable, and placeable in a
    config-driven solver dispatch dict — it was previously ~70 lines of
    inline logic on :class:`~PyOMES.core.control_volume.ControlVolume` with
    no object identity of its own.

    ``cv.advance(dt_h, t_h)`` and
    ``cv.advance(dt_h, t_h, solver=SequentialAdvanceSolver())`` are
    byte-identical calls — ``solver=None`` is pure sugar for this class,
    constructed with its default ``clamp_fn``.

    Parameters
    ----------
    clamp_fn : callable or None
        ``clamp_fn(deltas, current_mol, dt_h) -> deltas``, applied
        independently to each sub-step's delta dict (feeds, each
        boundary, reactions) as it is computed — this solver is
        sequential/Lie-Trotter, not simultaneous, so there is no single
        combined delta dict to clamp once the way
        :class:`SimultaneousEulerSolver` does. Default
        :func:`~PyOMES.core.clamping.proportional_clamp`. ``None``
        disables clamping — every sub-step is applied with
        ``clamp=False``, so a large enough ``dt_h``/rate combination
        can genuinely drive a species negative (see
        :class:`~PyOMES.monitoring.accuracy.AccuracyMonitor`'s
        ``negative_mole`` check).
    """

    def __init__(self, clamp_fn: Optional[Any] = proportional_clamp):
        self.clamp_fn = clamp_fn

    def solve_step(
        self,
        cv: "ControlVolume",
        dt_h: float,
        t_h: float = 0.0,
        external_source_terms: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> "AdvanceResult":
        # 1. Speciation solve on the pre-step state. The engine
        #    writes derived species back to phase.n_mol via
        #    _refresh_derived; subsequent sub-steps read them from
        #    n_mol directly.
        # ``reaction_system`` here is a ReactionSystem proper when
        # going through the canonical builder path; some test
        # fixtures and the factory's single-substrate convenience
        # path attach a bare KineticReaction / BlackBoxReactionModel
        # instead, which has no ``.engine`` lazy property — guard
        # accordingly.
        if cv.reaction_system is not None and hasattr(cv.reaction_system, "engine"):
            engine = cv.reaction_system.engine
            if engine is not None:
                liq = cv.phases.get("liquid")
                if liq is not None:
                    result = engine.solve(phases=cv.phases, T_K=float(liq.T_K))
                    result.apply_to_phases(cv.phases)

        # 1b. Property calculators (viscosity, density, ...). Run on
        #     the post-speciation state so derived species are
        #     visible if a calculator needs them. Results land on
        #     phase.properties[calc.key] for kinetic rate laws to
        #     read via env.prop(key). See
        #     PyOMES/core/property_calculator.py for the protocol.
        if cv.property_calculators:
            cv._run_property_calculators()

        # 2. External source terms (feeds, dosing, inter-CV transport).
        # clamp_fn applied per phase_key's delta dict as it is computed
        # (this solver is sequential, not simultaneous — there is no
        # single combined dict to clamp once). clamp=False on apply_flux:
        # clamping (or its deliberate absence, clamp_fn=None) already
        # happened above; apply_flux's own floor must not override it.
        if external_source_terms:
            for phase_key, species_rates in external_source_terms.items():
                if phase_key in cv.phases:
                    rates = species_rates
                    if self.clamp_fn is not None:
                        rates = self.clamp_fn(
                            rates, cv.phases[phase_key].n_mol, dt_h
                        )
                        _monitor_clamp_invoked(cv, species_rates, rates, dt_h)
                    cv.phases[phase_key].apply_flux(rates, dt_h, clamp=False)

        # 2b. Apply state-dependent boundary fluxes.
        for boundary in cv.boundaries:
            flux = boundary.compute_flux(cv, dt_h)
            phase_key = boundary.phase_key
            if phase_key in cv.phases:
                if self.clamp_fn is not None:
                    clamped_flux = self.clamp_fn(flux, cv.phases[phase_key].n_mol, dt_h)
                    _monitor_clamp_invoked(cv, flux, clamped_flux, dt_h)
                    flux = clamped_flux
                cv.phases[phase_key].apply_flux(flux, dt_h, clamp=False)

        # 3. Reactions (post-feed state). compute_reaction_rates() reads the
        #    speciation state pinned by step 1 (operator-splitting contract:
        #    feeds in step 2 do not shift the same-step pH). compute_rhs()
        #    would re-run speciation on the post-feed totals and break that
        #    invariant, so advance() uses compute_reaction_rates directly.
        rxn_sources = None
        if cv.reaction_system is not None:
            rhs = cv.compute_reaction_rates(t_h)
            rxn_sources = {}
            for pk, sp_rates in rhs.items():
                rxn_sources[pk] = {}
                phase = cv.phases[pk]
                clamped_rates = sp_rates
                if self.clamp_fn is not None:
                    clamped_rates = self.clamp_fn(sp_rates, phase.n_mol, dt_h)
                    _monitor_clamp_invoked(cv, sp_rates, clamped_rates, dt_h)
                for sp, rate in sp_rates.items():
                    applied_rate = clamped_rates.get(sp, rate)
                    phase.n_mol[sp] = phase.n_mol.get(sp, 0.0) + applied_rate * dt_h
                    rxn_sources[pk][sp] = rate

        # 4. Internal transfer (gas-liquid Henry, kinetic kLa, …).
        transfer_diag = cv.step_internal_transfer(dt_h)

        # 5. ConservationMonitor check (state-unification C6).
        #    Element + charge accounting on the post-step state.
        if cv._conservation_monitor is not None:
            cv._conservation_monitor.check_step(cv.phases)

        # clamp_fn=None: nothing above guaranteed non-negativity: check
        # whether that actually produced negative inventory anywhere.
        if self.clamp_fn is None:
            for phase in cv.phases.values():
                _monitor_negative_mole(cv, phase.n_mol)

        return AdvanceResult(
            transfer=transfer_diag,
            reaction_sources=rxn_sources,
        )


# ════════════════════════════════════════════════════════════════════════
#  Snapshot helpers
# ════════════════════════════════════════════════════════════════════════

class _SnapshotCV:
    """Read-only wrapper around snapshot phases for boundaries/links.

    Boundaries call ``cv.phases[key]`` to read state; that is the only
    surface this wrapper needs to provide.
    """

    def __init__(self, **phases):
        self.phases = dict(phases)

    def __contains__(self, key: str) -> bool:
        return key in self.phases

    def __getitem__(self, key: str) -> Phase:
        return self.phases[key]


def _monitor_pH_post_step(cv) -> None:
    """Read the final pH off ``cv.phases["liquid"].pH`` and feed
    it to ``cv._accuracy_monitor``.

    No-op if the CV pre-dates Phase 4 (the monitor attribute is set
    in ``ControlVolume.__init__`` from chemistry-unification-4
    onwards; the ``getattr`` guard keeps test fixtures that
    hand-construct phase-and-solver scenarios without going through
    ``ControlVolume.__init__`` working). Also no-op if the CV has
    no speciation populated ("H+" not in n_mol).
    """
    monitor = getattr(cv, "_accuracy_monitor", None)
    if monitor is None:
        return
    liq = cv.phases.get("liquid") if hasattr(cv, "phases") else None
    if liq is None or "H+" not in liq.n_mol:
        return
    try:
        monitor.check_pH_jump(float(liq.pH))
    except ValueError:
        return


def _monitor_clamp_invoked(
    cv, deltas_before: Dict[str, float], deltas_after: Dict[str, float], dt_h: float,
) -> None:
    """Feed a before/after clamp_fn delta-dict pair to
    ``cv._accuracy_monitor`` (STEP_SOLVER_INTERFACE_REFINEMENT.md
    item 8b). No-op if the CV has no monitor attached."""
    monitor = getattr(cv, "_accuracy_monitor", None)
    if monitor is None:
        return
    monitor.check_clamp_invoked(deltas_before, deltas_after, dt_h)


def _monitor_negative_mole(cv, n_mol: Dict[str, float]) -> None:
    """Feed a post-step n_mol dict to ``cv._accuracy_monitor`` — only
    meaningful when clamp_fn=None (STEP_SOLVER_INTERFACE_REFINEMENT.md
    item 8a). No-op if the CV has no monitor attached."""
    monitor = getattr(cv, "_accuracy_monitor", None)
    if monitor is None:
        return
    monitor.check_negative_mole(n_mol)


# ════════════════════════════════════════════════════════════════════════
#  SimultaneousEulerSolver
# ════════════════════════════════════════════════════════════════════════

class SimultaneousEulerSolver:
    """Explicit Euler with snapshot-based computation and swappable clamping.

    All sub-systems compute from a frozen snapshot of the phase state
    and deltas are applied in one batch (unsplit/simultaneous composition).

    Produces identical results to the pre-refactoring sequential advance
    for small timesteps, with improved consistency for large timesteps
    (no ordering dependence). Not the default — ``ControlVolume.advance()``'s
    inline sequential body runs when ``solver`` is omitted.

    Parameters
    ----------
    clamp_fn : callable or None
        ``clamp_fn(deltas, current_mol, dt_h) -> deltas`` — restores
        non-negativity on the combined delta dict before it's applied.
        Default :func:`~PyOMES.core.clamping.proportional_clamp`
        (preserves relative stoichiometry). Pass
        :func:`~PyOMES.core.clamping.floor_clamp` for a per-species
        floor, a bespoke composite (see ``PyOMES/core/clamping.py``'s
        module docstring), or ``None`` to disable clamping entirely —
        useful for diagnosing whether a given ``dt_h``/kinetics
        combination is aggressive enough to need it (see
        :class:`~PyOMES.monitoring.accuracy.AccuracyMonitor`'s
        ``negative_mole`` check).
    """

    def __init__(self, clamp_fn: Optional[Any] = proportional_clamp):
        self.clamp_fn = clamp_fn

    def solve_step(
        self,
        cv: Any,
        dt_h: float,
        t_h: float,
        external_source_terms: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> AdvanceResult:
        """Advance the CV by one Euler step using the snapshot pattern.

        Phase-generic (STEP_SOLVER_INTERFACE_REFINEMENT.md item 6):
        snapshots and deltas cover whatever phases ``cv`` has (gas,
        liquid, solid, or any subset) — only a ``"liquid"`` phase is
        required, since speciation and reactions are liquid-scoped.
        Every registered ``cv.internal_interfaces`` entry is processed
        generically via its ``phase_a_key``/``phase_b_key``, not just
        a single ``KineticGasLiquidLink`` found by isinstance search.

        state-unification C4d: speciation routes through
        cv.reaction_system.engine; PropertyResult / property_solvers
        gone.
        """
        # Avoid circular imports
        from .control_volume import ControlVolume

        # ── 0. VALIDATE phase keys ──────────────────────────────────
        if "liquid" not in cv.phases:
            raise ValueError(
                f"SimultaneousEulerSolver requires a 'liquid' phase "
                f"(speciation and reactions are liquid-scoped); got "
                f"{list(cv.phases.keys())}"
            )

        # ── 1. SNAPSHOT (every phase the CV has) ──────────────────────
        snap_phases: Dict[str, Phase] = {
            k: p.snapshot() for k, p in cv.phases.items()
        }
        snap_cv = _SnapshotCV(**snap_phases)
        snap_liq = snap_phases["liquid"]

        # ── 2. COMPUTE all deltas from snapshot ──────────────────────
        deltas: Dict[str, Dict[str, float]] = {k: {} for k in cv.phases}

        # 2a. Boundaries
        boundary_records: List[ExternalFluxRecord] = []
        for boundary in cv.boundaries:
            flux = boundary.compute_flux(snap_cv, dt_h)
            mol_applied = {sp: rate * dt_h for sp, rate in flux.items()}
            boundary_records.append(ExternalFluxRecord(
                boundary_label=boundary.label,
                phase_key=boundary.phase_key,
                flux_mol_per_h=dict(flux),
                dt_h=float(dt_h),
                mol_applied=mol_applied,
            ))
            target = deltas.setdefault(boundary.phase_key, {})
            for sp, rate in flux.items():
                target[sp] = target.get(sp, 0.0) + rate

        # 2b. Speciation (on the full snapshot phase set). state-unification
        # C4d: the engine runs directly via cv.reaction_system.engine;
        # apply_to_phases() commits derived species (H+, OH-, CO2aq,
        # HCO3-, CO3--, NH3, NH4+, strong ions) back to snap_liq.n_mol.
        # Internal interfaces then read alpha values inline from
        # snap_liq.n_mol; the reaction rate evaluation also reads pH
        # from snap_liq.pH directly.
        if cv.reaction_system is not None:
            engine = getattr(cv.reaction_system, "engine", None)
            if engine is not None:
                result = engine.solve(
                    phases=snap_phases,
                    T_K=float(snap_liq.T_K),
                )
                result.apply_to_phases(snap_phases)

        # 2c. Internal transfer — every registered PhaseInterface
        # (gas-liquid Henry, liquid-solid Ksp, or any future pair),
        # not just one gas-liquid link. Each interface's own
        # phase_a_key/phase_b_key select which snapshot phases it
        # reads/writes. Positive flux = a→b (PhaseInterface convention).
        transfer_flows: List[tuple] = []
        for iface in cv.internal_interfaces:
            state_a = snap_phases[iface.phase_a_key]
            state_b = snap_phases[iface.phase_b_key]
            flow = iface.compute_flux(state_a, state_b, dt_h)
            transfer_flows.append((iface, flow))
            a_target = deltas.setdefault(iface.phase_a_key, {})
            b_target = deltas.setdefault(iface.phase_b_key, {})
            for sp, rate in flow.items():
                a_target[sp] = a_target.get(sp, 0.0) - rate
                b_target[sp] = b_target.get(sp, 0.0) + rate

        # 2d. Reactions (from snapshot liquid with injected speciation).
        # A read-only CV is built inline so compute_reaction_rates can
        # build a ReactionEnvironment from the frozen state.
        rxn_sources: Optional[Dict[str, Dict[str, float]]] = None
        if cv.reaction_system is not None:
            snap_liq_cv = ControlVolume(
                phases={"liquid": snap_liq},
                reaction_system=cv.reaction_system,
                label="snapshot_liquid",
            )
            rxn_sources = snap_liq_cv.compute_reaction_rates(t_h)

        # ── 3. COLLECT reactions + external source terms into deltas ──
        if rxn_sources is not None:
            for phase_key, sp_rates in rxn_sources.items():
                if phase_key not in cv.phases:
                    continue
                target = deltas.setdefault(phase_key, {})
                for sp, rate in sp_rates.items():
                    target[sp] = target.get(sp, 0.0) + rate

        if external_source_terms:
            for phase_key, species_rates in external_source_terms.items():
                if phase_key not in cv.phases:
                    continue
                target = deltas.setdefault(phase_key, {})
                for sp, rate in species_rates.items():
                    target[sp] = target.get(sp, 0.0) + rate

        # ── 4. CLAMP (per phase) ───────────────────────────────────────
        # clamp_fn=None deliberately skips this step so raw (possibly
        # negative-inducing) deltas reach APPLY below with clamp=False
        # — apply_flux's own floor must not silently re-clamp them, or
        # clamp_fn=None's diagnostic purpose (see AccuracyMonitor's
        # negative_mole check) would be defeated.
        if self.clamp_fn is not None:
            for phase_key in list(deltas):
                clamped = self.clamp_fn(
                    deltas[phase_key], snap_phases[phase_key].n_mol, dt_h
                )
                _monitor_clamp_invoked(cv, deltas[phase_key], clamped, dt_h)
                deltas[phase_key] = clamped

        # ── 5. APPLY to live state ───────────────────────────────────
        # clamp=False: clamping (if any) already happened in step 4
        # over the *combined* multi-source delta dict; apply_flux's
        # own per-call floor would be redundant when clamp_fn ran, and
        # would wrongly mask clamp_fn=None's negative values otherwise.
        for phase_key, phase_deltas in deltas.items():
            cv.phases[phase_key].apply_flux(phase_deltas, dt_h, clamp=False)

        # clamp_fn=None: nothing above guaranteed non-negativity; check
        # whether that actually produced negative inventory anywhere.
        if self.clamp_fn is None:
            for phase_key in deltas:
                _monitor_negative_mole(cv, cv.phases[phase_key].n_mol)

        # ── 6. FINAL SPECIATION on live state ────────────────────────
        # state-unification C4d: re-solve on the post-step phase so
        # n_mol's derived species reflect the new totals. apply_to_phases()
        # commits n_mol["H+"] etc. so subsequent reads of
        # cv.phases["liquid"].pH see the post-step value.
        if cv.reaction_system is not None:
            engine = getattr(cv.reaction_system, "engine", None)
            if engine is not None:
                result = engine.solve(
                    phases=cv.phases,
                    T_K=float(cv.phases["liquid"].T_K),
                )
                result.apply_to_phases(cv.phases)

        # Accuracy monitor: pH-jump check across the step-closing
        # speciation re-solve. Defers silently when the CV pre-dates
        # the chemistry-unification-4 monitor wiring.
        _monitor_pH_post_step(cv)

        # ── 7. Build result ──────────────────────────────────────────
        # AdvanceResult.transfer_record is a single optional slot (not
        # a list); if more than one interface produced nonzero flux,
        # only the first is captured here — step_internal_transfer's
        # TransferDiagnostics (used by the sequential path) remains the
        # complete multi-interface audit trail.
        transfer_record = None
        for iface, flow in transfer_flows:
            if flow:
                total_transferred = sum(
                    abs(v) * float(dt_h) for v in flow.values()
                )
                transfer_record = LinkFlowRecord(
                    link_label=iface.label,
                    source=iface.phase_a_key,
                    sink=iface.phase_b_key,
                    flow_mol_per_h=dict(flow),
                    total_mol_transferred=total_transferred,
                )
                break

        return AdvanceResult(
            reaction_sources=rxn_sources,
            transfer_record=transfer_record,
            boundary_records=boundary_records,
        )

    def __repr__(self):
        return "SimultaneousEulerSolver()"


# ════════════════════════════════════════════════════════════════════════
#  SimultaneousAdaptiveSolver
# ════════════════════════════════════════════════════════════════════════

class SimultaneousAdaptiveSolver:
    """Adaptive ODE solver via ``scipy.integrate.solve_ivp``.

    Unsplit/simultaneous composition, same as :class:`SimultaneousEulerSolver`:
    handles all sub-systems (boundaries, transfer, reactions, speciation)
    within a single ODE derivative function — but integrates continuously
    with adaptive step control instead of a single Euler batch. Speciation
    is solved algebraically at every internal sub-step (DAE approach,
    matching BSM2).

    Parameters
    ----------
    rtol : float
        Relative tolerance for solve_ivp.  Default 1e-6.
    atol : float
        Absolute tolerance for solve_ivp.  Default 1e-9.
    max_step : float
        Maximum internal step size (hours).  Default 1.0.
    method : str
        scipy solve_ivp method.  Default ``"DOP853"`` (explicit,
        high-order).  For stiff systems (e.g. AD with H₂ kinetics),
        use ``"Radau"`` or ``"BDF"`` (implicit).
    freeze_speciation : bool
        If True, solve speciation once at the start of each macro
        step and reuse that result for all internal ODE sub-steps.
        This mimics BSM2/PyADM1's approach where ion states are
        frozen during integration and updated algebraically after
        the step.  Default False (solve speciation at every
        derivative evaluation — more physically correct but gives
        different gas transfer dynamics than BSM2).
    """

    def __init__(self, rtol: float = 1e-6, atol: float = 1e-9,
                 max_step: float = 1.0, method: str = "DOP853",
                 freeze_speciation: bool = False,
                 use_engine_jacobian: bool = False):
        self.rtol = rtol
        self.atol = atol
        self.max_step = max_step
        self.method = method
        self.freeze_speciation = freeze_speciation
        self.use_engine_jacobian = use_engine_jacobian
        if freeze_speciation and use_engine_jacobian:
            raise ValueError(
                "freeze_speciation=True and use_engine_jacobian=True are incompatible: "
                "the Jacobian path requires live speciation solves inside the ODE."
            )
        if freeze_speciation:
            import warnings
            warnings.warn(
                "freeze_speciation=True: speciation will be computed once per "
                "ODE step and reused for all internal substeps. This improves "
                "performance but may reduce accuracy for pH-sensitive systems.",
                UserWarning, stacklevel=2,
            )
    def solve_step(
        self,
        cv: Any,
        dt_h: float,
        t_h: float,
        external_source_terms: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> AdvanceResult:
        """Advance by dt_h using adaptive scipy ODE solver.

        Phase-generic (STEP_SOLVER_INTERFACE_REFINEMENT.md item 6,
        checkpoint 7b): the ODE state vector, RHS closure, and Jacobian
        cover whatever phases ``cv`` has — only a ``"liquid"`` phase is
        required, since speciation and reactions are liquid-scoped
        (matching ``SequentialAdvanceSolver``/``SimultaneousEulerSolver``).
        Every registered ``cv.internal_interfaces`` entry is processed
        generically via its ``phase_a_key``/``phase_b_key``, not just a
        single gas-liquid link found by isinstance search.

        state-unification C4d: speciation routes through
        cv.reaction_system.engine; PropertyResult / property_solvers
        gone.
        """
        import numpy as np
        from scipy.integrate import solve_ivp
        from .state_vector import StateVector

        # ── 0. VALIDATE phase keys ──────────────────────────────────
        if "liquid" not in cv.phases:
            raise ValueError(
                f"SimultaneousAdaptiveSolver requires a 'liquid' phase "
                f"(speciation and reactions are liquid-scoped); got "
                f"{list(cv.phases.keys())}"
            )

        # _SV_KEY is an internal placeholder for this solver's single-CV
        # StateVector scope; never surfaced to the caller.
        _SV_KEY = "_cv"
        sv = StateVector({_SV_KEY: cv}, exclude_species=frozenset({"H+"}))
        y0 = sv.pack()
        phase_layout = sv.phase_species(_SV_KEY)   # {phase_key: [species, ...]}
        phase_offset = sv.phase_offset(_SV_KEY)     # {phase_key: starting index}

        # Save original state for restoration on failure
        y0_copy = y0.copy()

        # Capture CV internals for the derivative function
        boundaries = cv.boundaries
        liq_phase = cv.phases["liquid"]
        liq_V_L = liq_phase.V_L
        liq_T_K = liq_phase.T_K

        # Pre-assemble external source rate arrays, per phase
        ext_rates: Dict[str, Dict[str, float]] = {pk: {} for pk in phase_layout}
        if external_source_terms:
            for pk, sp_rates in external_source_terms.items():
                if pk not in ext_rates:
                    continue
                target = ext_rates[pk]
                for sp, rate in sp_rates.items():
                    target[sp] = target.get(sp, 0.0) + rate

        # Pre-cache references for hot loop. state-unification C4d:
        # speciation routes through cv.reaction_system.engine
        # directly. PropertyResult / SpeciationPropertySolver are
        # gone; the engine writes derived species back to the
        # snap_liq phase's n_mol each call.
        rxn_model = cv.reaction_system
        has_rxn = rxn_model is not None
        engine = getattr(rxn_model, "engine", None) if has_rxn else None
        has_speciation = engine is not None

        # Phase 5 — detect gray-box capability for analytical Jacobian
        import warnings as _warnings
        from ..chemical_equilibrium.protocols import GrayBoxEngineProtocol
        _JAC_METHODS = frozenset({"BDF", "Radau", "LSODA"})
        _engine_is_graybox = (
            has_speciation
            and isinstance(engine, GrayBoxEngineProtocol)
            and self.method in _JAC_METHODS
        )
        if _engine_is_graybox and not self.use_engine_jacobian:
            _warnings.warn(
                "engine satisfies GrayBoxEngineProtocol but use_engine_jacobian=False. "
                "Pass use_engine_jacobian=True to SimultaneousAdaptiveSolver to enable analytical Jacobian.",
                UserWarning,
                stacklevel=2,
            )

        # Mutable slot used by jac_callable to skip speciation re-solve for
        # component-column perturbations; None = run normal speciation
        _alg_override = [None]

        # Pre-build reusable phase objects (mutated in-place each call),
        # one per phase the CV actually has. Every phase is topped up
        # with its *current* full n_mol (not just ODE-tracked species)
        # so excluded/algebraic species (H+ on "liquid") have a home for
        # speciation context and the reaction environment — a no-op for
        # phases with no excluded species (gas, solid).
        _temp_phases: Dict[str, Any] = {}
        for pk, species in phase_layout.items():
            live_phase = cv.phases[pk]
            seed = dict.fromkeys(species, 0.0)
            for sp, val in live_phase.n_mol.items():
                if sp not in seed:
                    seed[sp] = val
            _temp_phases[pk] = type(live_phase)(
                n_mol=seed, V_L=live_phase.V_L, T_K=live_phase.T_K,
            )
        _temp_liq = _temp_phases["liquid"]
        _snap_cv = _SnapshotCV(**_temp_phases)

        # Import ReactionEnvironment once
        if has_rxn:
            from ..reactions.environment import ReactionEnvironment

        # Frozen-speciation cache (when freeze_speciation=True, the
        # engine solve runs once at the first f() call and the
        # resulting n_mol H+ / OH- / ladder species are pinned for
        # subsequent f() calls within the same solve_ivp).
        _frozen_pH = [None]
        _do_freeze = self.freeze_speciation and has_speciation

        def f(t, y):
            # 1+2. Unpack into reusable phase objects (in-place mutation)
            #      and write to live CV (kept for any side-effects;
            #      chem_env_fn no longer reads from this since strong
            #      ions live in phase.n_mol post-C4). floor_nonnegative
            #      guards against adaptive-integrator overshoot on a
            #      trial state -- a distinct concern from clamp_fn
            #      (discrete-step solvers only; see clamping.py).
            y_safe = floor_nonnegative(y)
            for pk, species in phase_layout.items():
                offset = phase_offset[pk]
                temp_phase = _temp_phases[pk]
                live_phase = cv.phases[pk]
                for i, sp in enumerate(species):
                    val = float(y_safe[offset + i])
                    temp_phase.n_mol[sp] = val
                    live_phase.n_mol[sp] = val

            # 3. Speciation (algebraic — the DAE part). apply_to_phases()
            #    commits derived species back to _temp_liq.n_mol;
            #    subsequent reaction RHS evaluation reads pH from
            #    _temp_liq.pH directly. Committed unconditionally on every
            #    RHS evaluation (including rejected trial steps) — safe
            #    because _temp_phases are scratch snapshots rebuilt each
            #    call and discarded except when the accepted final state
            #    is persisted via sv.unpack(..., floor=True).
            pH = None
            spec_result = None  # engine writes to n_mol directly; no result object
            if has_speciation:
                if _do_freeze and _frozen_pH[0] is not None:
                    # Reuse cached H+ from start of step
                    _temp_liq.n_mol["H+"] = 10.0 ** (-_frozen_pH[0]) * liq_V_L
                    pH = _frozen_pH[0]
                elif _alg_override[0] is not None:
                    # Jacobian path: use analytically predicted algebraic state
                    for sp, val in _alg_override[0].items():
                        _temp_liq.n_mol[sp] = float(val)
                    try:
                        pH = float(_temp_liq.pH)
                    except (ValueError, AttributeError):
                        pH = None
                else:
                    result = engine.solve(
                        phases=_temp_phases,
                        T_K=float(_temp_liq.T_K),
                    )
                    result.apply_to_phases(_temp_phases)
                    try:
                        pH = float(_temp_liq.pH)
                    except ValueError:
                        pH = None
                    if _do_freeze and pH is not None:
                        _frozen_pH[0] = pH

            # 5. Boundary fluxes (reuse snap objects — already point to _temp phases)
            boundary_rates: Dict[str, Dict[str, float]] = {pk: {} for pk in phase_layout}
            for b in boundaries:
                flux = b.compute_flux(_snap_cv, dt_h, instantaneous=True)
                if b.phase_key not in boundary_rates:
                    continue
                target = boundary_rates[b.phase_key]
                for sp, rate in flux.items():
                    target[sp] = target.get(sp, 0.0) + rate

            # 6. Internal transfer — every registered PhaseInterface
            # (gas-liquid Henry, liquid-solid Ksp, or any future pair),
            # not just one gas-liquid link. dn/dt needed here, not a
            # step-averaged rate, hence instantaneous=True.
            transfer_rates: Dict[str, Dict[str, float]] = {pk: {} for pk in phase_layout}
            for iface in cv.internal_interfaces:
                state_a = _temp_phases[iface.phase_a_key]
                state_b = _temp_phases[iface.phase_b_key]
                flow = iface.compute_flux(state_a, state_b, dt_h, instantaneous=True)
                a_target = transfer_rates.setdefault(iface.phase_a_key, {})
                b_target = transfer_rates.setdefault(iface.phase_b_key, {})
                for sp, rate in flow.items():
                    a_target[sp] = a_target.get(sp, 0.0) - rate
                    b_target[sp] = b_target.get(sp, 0.0) + rate

            # 7. Reaction rates (direct — avoid ControlVolume construction)
            rxn_rates: Dict[str, Dict[str, float]] = {pk: {} for pk in phase_layout}
            if has_rxn:
                conc = {sp: float(n) / liq_V_L
                        for sp, n in _temp_liq.n_mol.items()}
                properties = {}
                if spec_result is not None:
                    ionic = getattr(spec_result, "ionic_strength", None)
                    if ionic is not None:
                        properties["ionic_strength"] = float(ionic)
                    for k, v in getattr(spec_result, "species", {}).items():
                        properties[k] = float(v)
                env = ReactionEnvironment(
                    T_K=liq_T_K, V_L=liq_V_L, pH=pH,
                    concentrations=conc, properties=properties,
                    t_h=t,
                )
                rxn_sources = rxn_model.compute_rates(env)
                for phase_key, sp_rates in rxn_sources.items():
                    if phase_key not in rxn_rates:
                        continue
                    target = rxn_rates[phase_key]
                    for sp, val in sp_rates.items():
                        target[sp] = target.get(sp, 0.0) + val

            # 8. Assemble dy/dt
            dydt = np.zeros_like(y)
            for pk, species in phase_layout.items():
                offset = phase_offset[pk]
                b_rates = boundary_rates.get(pk, {})
                t_rates = transfer_rates.get(pk, {})
                r_rates = rxn_rates.get(pk, {})
                e_rates = ext_rates.get(pk, {})
                for i, sp in enumerate(species):
                    rate = b_rates.get(sp, 0.0)
                    rate += t_rates.get(sp, 0.0)
                    rate += r_rates.get(sp, 0.0)
                    rate += e_rates.get(sp, 0.0)
                    dydt[offset + i] = rate

            return dydt

        # Build analytical Jacobian callable (Phase 5) for implicit methods.
        # For each liquid-component column, skips the speciation re-solve by
        # using the analytical dz_dy prediction; other columns use forward FD.
        jac_callable = None
        if self.use_engine_jacobian and _engine_is_graybox:
            liq_species = phase_layout.get("liquid", [])
            liq_offset = phase_offset.get("liquid", 0)
            n_liq = len(liq_species)
            _diff_species = set(liq_species)  # ODE differential species; excluded from override

            def jac_callable(t, y):
                n = sv.n_total
                J = np.zeros((n, n))

                # Baseline: full evaluation with speciation (populates engine cache)
                f0 = f(t, y)

                # Component totals (mol/L) for dz_dy lookup
                y_safe = floor_nonnegative(y)
                totals = {
                    sp: float(y_safe[liq_offset + i]) / liq_V_L
                    for i, sp in enumerate(liq_species)
                }
                try:
                    spec_jac = engine.jacobian_dz_dy(totals=totals)
                    comp_id_list = list(spec_jac.component_ids)
                    comp_id_set = set(comp_id_list)
                    # Moles of secondary algebraic species at baseline
                    z_base = {
                        sp: float(_temp_liq.n_mol.get(sp, 0.0))
                        for sp in spec_jac.algebraic_ids
                        if sp not in _diff_species
                    }
                except (RuntimeError, ValueError):
                    spec_jac = None
                    comp_id_set = set()
                    z_base = {}

                for j in range(n):
                    eps = max(1e-8, 1e-5 * abs(float(y[j])))
                    y_pert = y.copy()
                    y_pert[j] += eps

                    # For liquid species that are speciation components: use
                    # analytical prediction for secondary algebraic species
                    liq_j = j - liq_offset
                    if spec_jac is not None and 0 <= liq_j < n_liq:
                        sp_j = liq_species[liq_j]
                        if sp_j in comp_id_set:
                            comp_j_idx = comp_id_list.index(sp_j)
                            # dz_dy[i,j] = ∂c_i/∂C_j (mol/L per mol/L);
                            # eps in moles → Δn_alg = dz_dy * eps
                            _alg_override[0] = {
                                sp: z_base[sp] + spec_jac.dz_dy[i, comp_j_idx] * eps
                                for i, sp in enumerate(spec_jac.algebraic_ids)
                                if sp not in _diff_species
                            }
                            J[:, j] = (f(t, y_pert) - f0) / eps
                            _alg_override[0] = None
                            continue

                    # Standard forward difference for gas and non-component columns
                    J[:, j] = (f(t, y_pert) - f0) / eps

                return J

        # Run the adaptive solver
        result = solve_ivp(
            f, [0.0, dt_h], y0,
            method=self.method,
            rtol=self.rtol,
            atol=self.atol,
            max_step=self.max_step,
            dense_output=False,
            jac=jac_callable,
        )

        if not result.success:
            # Restore original state and warn
            sv.unpack(y0_copy, floor=True)
            import warnings
            warnings.warn(
                f"SimultaneousAdaptiveSolver ({self.method}) failed: {result.message}. "
                f"State restored to pre-step values.",
                RuntimeWarning,
            )
        else:
            # Write final state back to CV
            sv.unpack(result.y[:, -1], floor=True)

        # Accuracy monitor: scipy step-rejection rate. result.nfev is
        # the function-evaluation count (proxy for attempted steps);
        # len(result.t) is the accepted-step count. The ratio spikes
        # when the integrator is fighting stiffness or a tight
        # tolerance — a regime warning, not a solver-failure warning
        # (which is handled separately above with RuntimeWarning).
        monitor = getattr(cv, "_accuracy_monitor", None)
        if monitor is not None and result.success:
            n_accepted = int(getattr(result.t, "size", len(result.t)))
            n_attempted = int(getattr(result, "nfev", n_accepted))
            monitor.check_scipy_rejections(
                n_accepted=n_accepted,
                n_attempted=n_attempted,
            )

        # Final speciation on live state. state-unification C4d:
        # apply_to_phases() commits derived species (H+ included) back to
        # the live liquid phase's n_mol.
        if has_speciation:
            eq_result = engine.solve(phases=cv.phases, T_K=float(liq_phase.T_K))
            eq_result.apply_to_phases(cv.phases)

        # Accuracy monitor: pH-jump check across the step-closing
        # speciation re-solve.
        _monitor_pH_post_step(cv)

        # Build result (lightweight — no per-boundary records for adaptive solver)
        return AdvanceResult()

    def __repr__(self):
        return (f"SimultaneousAdaptiveSolver(method={self.method!r}, rtol={self.rtol}, "
                f"atol={self.atol}, max_step={self.max_step}, "
                f"use_engine_jacobian={self.use_engine_jacobian})")


def _write_temp_state(cv, gas_mol, liq_mol):
    """Temporarily write state dicts to CV phases.

    Used by SimultaneousAdaptiveSolver so context functions that read
    ``cv.phases["liquid"].n_mol`` see the current ODE state.
    """
    for sp, val in gas_mol.items():
        cv.phases["gas"].n_mol[sp] = val
    for sp, val in liq_mol.items():
        cv.phases["liquid"].n_mol[sp] = val
