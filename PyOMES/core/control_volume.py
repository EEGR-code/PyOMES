# -*- coding: utf-8 -*-
"""ControlVolume — container for phases connected by internal interfaces.

A ControlVolume owns one or more :class:`~PyOMES.core.phases.Phase`
objects (gas, liquid, solid) and zero or more internal
:class:`~PyOMES.core.interfaces.PhaseInterface` objects that compute
transfer between them.

External fluxes (feeds, vents, connections to other CVs) are applied
directly to individual phases via ``cv.phases[key].apply_flux()``.
The ControlVolume itself does not manage external connectivity — that
is the orchestrator's responsibility.

Usage
-----
>>> from PyOMES.core import GasPhase, LiquidPhase, ControlVolume
>>> gas = GasPhase({"O2": 1.0, "CO2": 0.1, "N2": 3.0}, V_L=200.0, T_K=305.15)
>>> liquid = LiquidPhase({"CO2": 0.01}, V_L=800.0, T_K=305.15)
>>> cv = ControlVolume(phases={"gas": gas, "liquid": liquid},
...                    internal_interfaces=[my_gl_interface])
>>> diag = cv.step_internal_transfer(dt_h=0.01)
>>> assert diag.is_conserved()
"""

from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional

import numpy as np

# Strong-corrector map: user-facing name → ion added to n_mol.
# Only the charge-carrying ion is added; the counter-ion (OH- for bases,
# H+ for acids) is produced automatically by the charge-balance solver.
_STRONG_CORRECTOR_ION: Dict[str, str] = {
    "NaOH": "Na+",
    "KOH":  "K+",
}

from .phases import Phase, GasPhase, LiquidPhase
from .interfaces import PhaseInterface, TransferDiagnostics, AdvanceResult
from .lifecycle import _LockableList, raise_if_running
from .solvers import SequentialAdvanceSolver


def _build_transfer_link(phases, transfer_models, phase_pair=None, label=""):
    """Construct a KineticGasLiquidLink from a transfer_models dict.

    Detects the GasPhase/LiquidPhase pair automatically unless
    ``phase_pair=("gas_key", "liquid_key")`` is provided explicitly.

    Returns the constructed link (with DeprecationWarning suppressed
    since this is an internal construction path, not a user call).
    """
    import warnings
    from .transfer_models import KineticTransferModel, EquilibriumTransferModel
    from .gas_liquid_link import KineticGasLiquidLink

    if phase_pair is not None:
        gas_key, liq_key = phase_pair
    else:
        gas_key = next(
            (k for k, p in phases.items() if isinstance(p, GasPhase)), None
        )
        liq_key = next(
            (k for k, p in phases.items() if isinstance(p, LiquidPhase)), None
        )
        if gas_key is None or liq_key is None:
            raise ValueError(
                "transfer_models requires a GasPhase and a LiquidPhase in "
                "cv.phases. Found: "
                + str({k: type(v).__name__ for k, v in phases.items()})
                + ". Pass phase_pair=('gas_key', 'liquid_key') to override "
                "auto-detection, or use internal_interfaces directly for "
                "non-standard phase pairs."
            )

    cv_key = label or "cv"
    partition_models = {sp: m.partition_model for sp, m in transfer_models.items()}
    kLa = {
        sp: m.k_transfer
        for sp, m in transfer_models.items()
        if isinstance(m, KineticTransferModel)
    }
    equilibrium_species = {
        sp
        for sp, m in transfer_models.items()
        if isinstance(m, EquilibriumTransferModel)
    }
    molecular_driving_force = {
        sp
        for sp, m in transfer_models.items()
        if isinstance(m, KineticTransferModel) and m.transfer_basis == "molecular"
    }

    link_label = f"{label}_gl_transfer" if label else "gl_transfer"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        link = KineticGasLiquidLink(
            gas_cv_key=cv_key,
            gas_phase_key=gas_key,
            liquid_cv_key=cv_key,
            liquid_phase_key=liq_key,
            partition_models=partition_models,
            kLa=kLa,
            equilibrium_species=equilibrium_species,
            molecular_driving_force=molecular_driving_force,
            _label=link_label,
        )
    return link


class OrchestrationWarning(UserWarning):
    """A CV owned by a :class:`~PyOMES.core.simulation.Simulation` was
    advanced directly via :meth:`ControlVolume.advance` instead of
    through ``sim.run(...)``.

    Not an error — deliberately stepping one CV in isolation as a
    diagnostic (e.g. "what would this CV's kinetics alone look like")
    is legitimate. But it silently skips inter-CV links, controllers,
    and profiles for that step, which produces a complete-looking but
    physically wrong answer if the caller didn't intend it.
    """


class ControlVolume:
    """Container for phases and their internal interfaces.

    Parameters
    ----------
    phases : dict
        Mapping of ``{key: Phase}`` objects, e.g.
        ``{"gas": GasPhase(...), "liquid": LiquidPhase(...)}``.
    internal_interfaces : list of PhaseInterface
        Interfaces that transfer material between phases within this CV.
        Each interface's ``phase_a_key`` and ``phase_b_key`` must match
        keys in ``phases``.
    boundaries : list, optional
        :class:`~PyOMES.core.boundaries.ExternalBoundary` objects
        (gas feed, vent, membrane, liquid feed/drain, …) applied each
        timestep in ``advance()``.  Stored as a mutable list, so
        boundaries can be appended after construction.  The default
        sequential ``advance()`` body computes each boundary's flux
        from the current CV state and applies it via
        :meth:`Phase.apply_flux`; the snapshot solver (when used via
        ``advance(solver=...)``) computes from a frozen snapshot
        instead.
    reaction_system : object or None
        A :class:`~PyOMES.reactions.reaction_system.ReactionSystem`
        instance (or, for tests, any object implementing the
        :class:`~PyOMES.reactions.protocols.ReactionModel` protocol).
        If set, ``advance()`` integrates kinetic and black-box
        reactions after applying any external source terms.  The
        sequential ``advance()`` body applies a single forward Euler
        evaluation over ``dt_h`` for the reaction sub-step; callers
        needing higher accuracy should pass
        ``solver=SimultaneousAdaptiveSolver(...)`` to :meth:`advance`, which
        short-circuits the sequential body and uses an adaptive
        monolithic integrator. Equilibrium reactions held in the
        system are silently skipped during integration (they have
        no rate) and routed elsewhere by the system's internal
        pre-bucketing.
    label : str
        Human-readable label for this CV (e.g. ``"fermenter"``, ``"ring_1"``).
    transfer_models : dict, optional
        ``{species_id: TransferModel}`` — the canonical way to declare
        intra-CV phase mass transfer.  Each entry is either a
        :class:`~PyOMES.core.transfer_models.KineticTransferModel` (rate-limited,
        kLa-style) or an
        :class:`~PyOMES.core.transfer_models.EquilibriumTransferModel`
        (instantaneous partition).  A :class:`~PyOMES.core.gas_liquid_link.KineticGasLiquidLink`
        is constructed automatically from the dict and appended to
        ``internal_interfaces``.  Requires a :class:`~PyOMES.core.phases.GasPhase`
        and a :class:`~PyOMES.core.phases.LiquidPhase` to be present in
        ``phases`` (or provide ``phase_pair`` to override detection).
    phase_pair : tuple of (str, str), optional
        ``(gas_phase_key, liquid_phase_key)`` — overrides automatic
        phase-pair detection when the CV has more than two phases or
        non-standard phase keys.  Only used when ``transfer_models`` is
        provided.
    """

    def __init__(
        self,
        phases: Dict[str, Phase],
        internal_interfaces: Optional[List[PhaseInterface]] = None,
        boundaries: Optional[List] = None,
        reaction_system: Optional[Any] = None,
        property_calculators: Optional[List[Any]] = None,
        label: str = "",
        chemistry_db: Optional[Any] = None,
        transfer_models: Optional[Dict[str, Any]] = None,
        phase_pair: Optional[tuple] = None,
    ):
        self.phases = dict(phases)
        self.label = str(label)
        # transfer_models is the canonical construction-time interface for
        # intra-CV phase transfer. The factory builds a KineticGasLiquidLink
        # and appends it to internal_interfaces.
        # TODO Option B: make cv.transfer_models a live dict that propagates
        # k_transfer changes to the internal link at runtime (needed for
        # controller write paths via ParamPath).
        self.transfer_models: Dict[str, Any] = dict(transfer_models or {})
        # _explicit_interfaces tracks interfaces passed directly by the caller
        # (not auto-generated from transfer_models) so snapshot() can
        # reconstruct without duplicating the auto-generated link.
        self._explicit_interfaces = list(internal_interfaces or [])
        ifaces = list(self._explicit_interfaces)
        if self.transfer_models:
            ifaces.append(
                _build_transfer_link(
                    self.phases, self.transfer_models,
                    phase_pair=phase_pair, label=self.label,
                )
            )
        self.internal_interfaces = ifaces
        # chemistry-unification-3b C9: optional ChemistryDatabase attaches
        # the ThermoFramework + species + reactions bundle to this CV.
        self.chemistry_db = chemistry_db
        # C6 lockable wrappers — refuse mutation while owning
        # Simulation is running. _context wired by Simulation.__init__.
        self.boundaries = _LockableList(
            boundaries or [],
            label=f"ControlVolume({self.label!r}).boundaries",
        )
        # state-unification C5: scalar derived-property calculators
        # (viscosity, density, ...) — runs once before kinetic
        # reactions in advance(); writes results to
        # phase.properties[calc.key]. See PropertyCalculator
        # protocol in PyOMES/core/property_calculator.py.
        self.property_calculators = _LockableList(
            property_calculators or [],
            label=f"ControlVolume({self.label!r}).property_calculators",
        )
        # C6: reaction_system reassignment is gated via a property.
        # Underscore storage is the actual reference; the property
        # getter / setter below provides the gate.
        self._reaction_system = reaction_system

        # Validate that every interface references existing phases
        for iface in self.internal_interfaces:
            a = iface.phase_a_key
            b = iface.phase_b_key
            if a not in self.phases:
                raise KeyError(
                    f"Interface references phase_a_key={a!r}, "
                    f"but phases only contains: {list(self.phases.keys())}"
                )
            if b not in self.phases:
                raise KeyError(
                    f"Interface references phase_b_key={b!r}, "
                    f"but phases only contains: {list(self.phases.keys())}"
                )

        # Validate cross-reaction Species consistency. Raise on hard
        # conflicts, warn on redeclaration soft conflicts.  Skip
        # silently for reaction models that don't expose declared
        # stoichiometry (e.g. a bare BlackBoxReactionModel attached
        # directly without wrapping in a ReactionSystem).
        if self.reaction_system is not None:
            from ..chemistry.species_check import check_species_consistency
            rxns = getattr(self.reaction_system, "reactions", None)
            if rxns is None and hasattr(self.reaction_system, "stoichiometry"):
                rxns = [self.reaction_system]
            if rxns is not None:
                check_species_consistency(rxns, soft_conflicts="ignore")

                # Drive gas-liquid link speciation_keys from declared
                # cross-phase equilibrium reactions
                # (chemistry-unification-3). Each link reads the
                # ReactionSystem's pre-bucketed cross_phase_equilibria
                # list and records its own gas → liquid mapping based
                # on its phase keys, so multi-link CVs wire correctly.
                # The method is additive — pre-existing speciation_keys
                # entries survive when no cross-phase reaction overrides
                # them; speciation_ladders is also built from the
                # single-phase equilibria graph.
                from .gas_liquid_link import KineticGasLiquidLink
                for iface in self.internal_interfaces:
                    if isinstance(iface, KineticGasLiquidLink):
                        iface.derive_speciation_keys(self.reaction_system)

        # AccuracyMonitor: per-CV check runner for cheap numerical
        # accuracy heuristics (pH change, Newton iters, scipy step
        # rejections, ionic-strength regime). Thresholds + throttle
        # are read off PyOMES.config.warnings at every emission, so
        # configuration changes between simulations propagate without
        # re-instantiation. state-unification C4d: monitor attaches
        # via cv.reaction_system.attach_monitor(monitor), which
        # propagates to the speciation engine on first build.
        from ..monitoring import AccuracyMonitor, ConservationMonitor
        self._accuracy_monitor = AccuracyMonitor()
        if self.reaction_system is not None and hasattr(
            self.reaction_system, "attach_monitor"
        ):
            self.reaction_system.attach_monitor(self._accuracy_monitor)

        # ConservationMonitor (state-unification C6): element +
        # charge accounting per advance() step. Auto-attached
        # when a ReactionSystem is present; the species registry
        # is harvested from the reaction stoichiometries so the
        # monitor can compute atom/charge contributions from
        # n_mol entries.
        self._conservation_monitor = ConservationMonitor()
        if self.reaction_system is not None and hasattr(
            self.reaction_system, "attach_conservation_monitor"
        ):
            self._conservation_monitor.set_species_registry(
                self._collect_species_registry()
            )
            self.reaction_system.attach_conservation_monitor(
                self._conservation_monitor
            )

        # τ_min estimate is left as None for now — a future phase can
        # plumb a Jacobian-free 1/max_rate estimate from the kinetic
        # reactions. check_dt_vs_tau_min no-ops on None, so the hook
        # is in place but the check is silently deferred.
        self._tau_min_h_estimate: Optional[float] = None

        # simulation-class C1: RunContext back-reference. Set by
        # Simulation.__init__ when this CV is added; consulted by
        # lockable mutators at C6. None when the CV is unowned
        # (freshly constructed, snapshotted, or already removed
        # from a Simulation). Annotated as Any to avoid a circular
        # import — the real type is PyOMES.core.simulation.RunContext.
        self._context: Optional[Any] = None

        # P3 species-vector-map cache (framework-polish).  Species
        # membership is invariant for the life of a ReactionSystem
        # (lifecycle gating blocks mid-run mutation), so the result of
        # _build_species_vector_map can be memoised keyed by
        # id(_reaction_system).  Cleared by _set_reaction_system_unchecked.
        self._species_vector_map_cache: Optional[tuple] = None
        self._species_vector_map_cache_rs_id: Optional[int] = None

    # ── Phase access helpers ───────────────────────────────────────────

    def __getitem__(self, key: str) -> Phase:
        """Shorthand: ``cv["gas"]`` is equivalent to ``cv.phases["gas"]``."""
        return self.phases[key]

    def __contains__(self, key: str) -> bool:
        return key in self.phases

    @property
    def phase_keys(self) -> list:
        return list(self.phases.keys())

    # ── Gated reaction_system access (C6) ──────────────────────────────

    @property
    def reaction_system(self) -> Optional[Any]:
        return self._reaction_system

    @reaction_system.setter
    def reaction_system(self, value: Optional[Any]) -> None:
        raise_if_running(self, "reaction_system")
        self._set_reaction_system_unchecked(value)

    def _set_reaction_system_unchecked(self, value: Optional[Any]) -> None:
        """Orchestrator-mediated unchecked path (Pattern B). Bypasses
        the lifecycle gate."""
        self._reaction_system = value
        self._species_vector_map_cache = None
        self._species_vector_map_cache_rs_id = None

    # ── System-wide mole inventory ─────────────────────────────────────

    def total_mol(self) -> Dict[str, float]:
        """Sum moles of each species across all phases.

        Returns a dict ``{species: total_mol}`` aggregated over every
        phase in the CV.
        """
        totals: Dict[str, float] = {}
        for phase in self.phases.values():
            for species, n in phase.total_mol().items():
                totals[species] = totals.get(species, 0.0) + float(n)
        return totals

    # ── Internal transfer step ─────────────────────────────────────────

    def _gas_liquid_engine_owned_species(self) -> frozenset:
        """Species this CV's speciation engine already resolves via a
        folded gas-liquid row (CP1/CP2 of ``LAYER1_GAP_CLOSURE``).

        Used by :meth:`step_internal_transfer` to skip a separately-
        declared ``transfer_models`` entry for the same species — without
        this, both the engine's own simultaneous fold and a
        ``KineticTransferModel``/``EquilibriumTransferModel`` would move
        the same species between phases, double-applying the transfer.

        Deliberately narrower than ``engine.algebraic_species()`` (which
        also includes ordinary acid-base species with no gas-liquid
        coupling at all — filtering against the full set would silently
        disable legitimate, still-needed ``transfer_models`` entries for
        those). Empty when there is no reaction-system engine, or the
        engine doesn't support gas-liquid folding (only
        ``NRChemicalEquilibriumEngine`` does, via ``gas_liquid_species()``,
        duck-typed rather than a formal protocol requirement since no
        other engine implements it yet).
        """
        if self.reaction_system is None or not hasattr(self.reaction_system, "engine"):
            return frozenset()
        engine = self.reaction_system.engine
        if engine is None:
            return frozenset()
        gas_liquid_species = getattr(engine, "gas_liquid_species", None)
        if gas_liquid_species is None:
            return frozenset()
        return frozenset(gas_liquid_species())

    def step_internal_transfer(
        self,
        dt_h: float,
    ) -> TransferDiagnostics:
        """Evaluate all internal interfaces and apply fluxes.

        For each interface:
          1. Call ``interface.compute_flux(phase_a, phase_b, dt_h)``
          2. Drop any species already resolved by a folded gas-liquid row
             in this CV's speciation engine (CP3 of
             ``LAYER1_GAP_CLOSURE`` — see
             :meth:`_gas_liquid_engine_owned_species`), so the engine's
             own simultaneous solve and this transfer model don't both
             move the same species between phases.
          3. Apply ``-flux`` to phase_a (source loses material)
          4. Apply ``+flux`` to phase_b (sink gains material)

        Returns :class:`TransferDiagnostics` with per-species mass balance
        residuals.  For purely conservative transfer (no reactions),
        ``diag.is_conserved()`` should be True.

        State-unification C4 dropped the ``property_results`` thread:
        kinetic interfaces (e.g. CO2 partitioning needing alpha) now
        read derived species inline from ``phase.n_mol`` populated
        by the speciation engine on the current step.

        Parameters
        ----------
        dt_h : float
            Timestep duration (hours).

        Returns
        -------
        TransferDiagnostics
        """
        # Snapshot moles before
        mol_before = self.total_mol()

        owned = self._gas_liquid_engine_owned_species()

        # Evaluate and apply each interface
        all_fluxes = []
        for iface in self.internal_interfaces:
            phase_a = self.phases[iface.phase_a_key]
            phase_b = self.phases[iface.phase_b_key]

            flux = iface.compute_flux(phase_a, phase_b, dt_h)
            if owned:
                flux = {sp: rate for sp, rate in flux.items() if sp not in owned}
            all_fluxes.append(dict(flux))

            # Apply: a loses, b gains
            neg_flux = {sp: -rate for sp, rate in flux.items()}
            phase_a.apply_flux(neg_flux, dt_h)
            phase_b.apply_flux(flux, dt_h)

        # Snapshot moles after
        mol_after = self.total_mol()

        # Compute residuals
        all_species = set(mol_before.keys()) | set(mol_after.keys())
        residual = {}
        for sp in sorted(all_species):
            before = mol_before.get(sp, 0.0)
            after = mol_after.get(sp, 0.0)
            residual[sp] = after - before

        total_before = sum(mol_before.values())
        total_after = sum(mol_after.values())

        return TransferDiagnostics(
            species_before=mol_before,
            species_after=mol_after,
            residual=residual,
            fluxes=all_fluxes,
            total_mol_before=total_before,
            total_mol_after=total_after,
        )

    # ── External flux application ──────────────────────────────────────

    def apply_external_flux(
        self,
        phase_key: str,
        flux_mol_per_h: Dict[str, float],
        dt_h: float,
    ) -> None:
        """Apply an external flux to a specific phase.

        This is used by the orchestrator for feeds, vents, and inter-CV
        transport.  It is *not* tracked by ``step_internal_transfer``
        diagnostics (because external fluxes intentionally change the
        CV's total inventory).

        Strong-corrector aliases (``"NaOH"`` → ``"Na+"``,
        ``"KOH"`` → ``"K+"``) are resolved here so that controller
        dosing via :class:`~PyOMES.control.cv_loops.PHController`
        is consistent with :meth:`equilibrate_to_pH`, which uses
        the same mapping.

        Parameters
        ----------
        phase_key : str
            Key into ``self.phases``.
        flux_mol_per_h : dict
            Species fluxes (mol/h).  Positive = into the phase.
        dt_h : float
            Timestep duration (hours).
        """
        resolved = {
            _STRONG_CORRECTOR_ION.get(sp, sp): mol_h
            for sp, mol_h in flux_mol_per_h.items()
        }
        self.phases[phase_key].apply_flux(resolved, dt_h)

    # ── Advance (reactions + equilibrium + properties) ─────────────────

    def advance(
        self,
        dt_h: float,
        t_h: float = 0.0,
        *,
        external_source_terms: Optional[Dict[str, Dict[str, float]]] = None,
        solver: Optional[Any] = None,
    ) -> AdvanceResult:
        """Advance the CV by one timestep.

        Default sequential body, implemented by
        :class:`~PyOMES.core.solvers.SequentialAdvanceSolver` (see
        docs/dev/implementation/shipped/ORDERING.md for the design rationale):

          1. Speciation solve on the pre-step state. The
             ``cv.reaction_system.engine`` writes derived species
             (H+, OH-, CO2, HCO3-, CO3--, NH3, NH4+, strong ions)
             to ``phase.n_mol`` so subsequent sub-steps read current
             speciation from the phase.
          2. Apply external source terms (orchestrator-supplied fluxes).
          2b. Apply state-dependent external boundaries (feeds, vents,
              membranes, …).  Each boundary's flux is computed from the
              current CV state and applied via
              :meth:`Phase.apply_flux`.  No diagnostic record is
              collected here — that is the snapshot solver's job.
          3. Integrate reactions (operating on the post-feed state,
             with speciation pinned at the pre-step values written to
             ``n_mol`` in step 1).
          4. Internal equilibrium / kinetic transfer (gas-liquid, etc.).

        Speciation runs first so subsequent sub-steps read pH /
        derived-species concentrations from ``phase.n_mol``
        directly. No more :class:`PropertyResult` threading — the
        phase IS the source of truth for derived state, populated
        by the engine's privileged ``_refresh_derived`` writeback.

        Feed runs *before* reactions because, in a continuously fed CV,
        material that arrives in the timestep is immediately available
        for reaction.  This is one valid Lie-Trotter operator splitting
        with O(dt) error; it is **not** equivalent to the simultaneous
        explicit-Euler step that ``SimultaneousEulerSolver`` performs (see
        ORDERING_CRITIQUE.md, Review 2).

        Solver dispatch
        ---------------
        ``solver.solve_step(self, dt_h, t_h, external_source_terms=...)``
        is always called — ``solver=None`` (the default) is pure sugar
        for ``SequentialAdvanceSolver()``, so this and every other
        :class:`StepSolver` (e.g. :class:`SimultaneousEulerSolver`,
        :class:`SimultaneousAdaptiveSolver`) are dispatched identically.

        Parameters
        ----------
        dt_h : float
            Timestep duration (hours).
        t_h : float
            Current simulation time (hours). Forwarded to the
            reaction environment so time-dependent rate laws (lag
            phases, scheduled feeds, induced kinetics) can read
            ``env.t_h``. Owner is the orchestrator
            (:class:`~PyOMES.core.simulation.Simulation` orchestrator).
        external_source_terms : dict, optional
            ``{phase_key: {species_id: mol_per_h}}`` — feed/dose rates
            applied to individual phases before reactions run.
        solver : StepSolver, optional
            Time-stepping strategy.  When ``None`` (default), dispatches
            to :class:`~PyOMES.core.solvers.SequentialAdvanceSolver`,
            which implements the sequential body documented above.

        Returns
        -------
        AdvanceResult
            Reaction sources and transfer diagnostics for this step.
            Derived species (pH, ionic strength, individual species
            concentrations) live on ``cv.phases[...]`` — read them
            from there.

        Ownership guard
        ----------------
        If this CV is owned by a :class:`~PyOMES.core.simulation.Simulation`
        (``self._context is not None``), calling this method directly
        emits :class:`OrchestrationWarning` — inter-CV links,
        controllers, and profiles will **not** be applied for this
        step. This remains a legitimate way to advance an *unowned*
        CV (unit tests, notebooks exploring a bare CV's behaviour) and
        is the trusted internal building block every ``SystemSolver``
        calls (via :meth:`_advance_unchecked`, which bypasses this
        guard). For an owned CV, prefer ``sim.run(...)`` — even for a
        single CV with zero links — since that is where ``t_h``
        bookkeeping, the recorder, and lifecycle safety actually live.
        """
        if self._context is not None:
            warnings.warn(
                "cv.advance() called directly on a CV owned by a "
                "Simulation; inter-CV links/controllers/profiles will "
                "NOT be applied this step. Use sim.run(...) unless "
                "this is intentional.",
                OrchestrationWarning,
                stacklevel=2,
            )
        return self._advance_unchecked(
            dt_h, t_h,
            external_source_terms=external_source_terms,
            solver=solver,
        )

    def _advance_unchecked(
        self,
        dt_h: float,
        t_h: float = 0.0,
        *,
        external_source_terms: Optional[Dict[str, Dict[str, float]]] = None,
        solver: Optional[Any] = None,
    ) -> AdvanceResult:
        """Trusted internal entry point for :meth:`advance`.

        Identical semantics, minus the :class:`OrchestrationWarning`
        ownership guard. Every ``SystemSolver.advance_system()``
        implementation calls this directly — they are the trusted
        orchestrator, so the guard must not fire on that path.
        """
        # solver=None is sugar for SequentialAdvanceSolver() — the
        # sequential body used to live here inline; it is now a real
        # StepSolver (PyOMES/core/solvers.py) so it can be passed
        # explicitly, wrapped, or placed in a config-driven dispatch
        # dict like any other StepSolver.
        solver = solver if solver is not None else SequentialAdvanceSolver()
        return solver.solve_step(
            self, dt_h, t_h,
            external_source_terms=external_source_terms,
        )

    # ── CV compute interface (Phase A — CV_COMPUTE_INTERFACE) ─────────────

    def compute_rhs(self, t_h: float) -> Dict[str, Dict[str, float]]:
        """ODE right-hand side: speciation write-back then kinetic rates.

        Runs the speciation engine (idempotent for a given totals state),
        then evaluates kinetic reaction rates from the current phase state.
        Returns ``{phase_key: {species_id: rate_mol_per_h}}``.  Does not
        mutate CV state beyond the speciation write-back.

        Used by ``advance()`` for its Euler step and by system-level
        solvers (Phase C+) as the per-CV contribution to the joint RHS.
        """
        if self.reaction_system is None:
            return {}
        if hasattr(self.reaction_system, "engine"):
            engine = self.reaction_system.engine
            if engine is not None:
                liq = self.phases.get("liquid")
                if liq is not None:
                    result = engine.solve(phases=self.phases, T_K=float(liq.T_K))
                    result.apply_to_phases(self.phases)
        return self.compute_reaction_rates(t_h)

    def compute_differential_rhs(
        self,
        t_h: float,
        y_state: Dict[str, Dict[str, float]],
        z_state: Dict[str, Dict[str, float]],
    ) -> Dict[str, Dict[str, float]]:
        """DAE interface: kinetic rates at explicit (y_state, z_state).

        Evaluates reaction rates without calling the speciation engine and
        without touching ``phase.n_mol``.  Structural metadata (V_L, T_K,
        properties) is read from ``self.phases``; species concentrations
        come from the merged ``y_state``/``z_state`` dicts.

        Required by ``DAEStepSolver`` (Phase F) and ``DAESystemSolver``
        (Phase G).
        """
        if self.reaction_system is None:
            return {}
        from ..reactions.environment import ReactionEnvironment

        liq = self.phases.get("liquid")
        if liq is not None:
            pk = "liquid"
            ref_phase = liq
        else:
            pk = next(iter(self.phases), None)
            if pk is None:
                return {}
            ref_phase = self.phases[pk]

        V_L = float(ref_phase.V_L)
        T_K = float(ref_phase.T_K)

        merged = {}
        merged.update(y_state.get(pk, {}))
        merged.update(z_state.get(pk, {}))

        concentrations = {sp: float(n) / V_L for sp, n in merged.items()}

        pH = None
        h_mol = merged.get("H+", 0.0)
        if h_mol > 0.0:
            pH = float(-np.log10(h_mol / V_L))

        properties: Dict[str, float] = {}
        properties.update(getattr(ref_phase, "properties", {}))
        ionic = None
        if hasattr(ref_phase, "speciation"):
            ionic = ref_phase.speciation.get("IonicStrength")
        if ionic is not None and "ionic_strength" not in properties:
            properties["ionic_strength"] = float(ionic)

        env = ReactionEnvironment(
            T_K=T_K,
            V_L=V_L,
            pH=pH,
            concentrations=concentrations,
            properties=properties,
            t_h=float(t_h),
        )
        rates = self.reaction_system.compute_rates(env)
        return {p: dict(sr) for p, sr in rates.items()}

    def compute_algebraic_residual(
        self,
        t_h: float,
        y_state: Dict[str, Dict[str, float]],
        z_state: Dict[str, Dict[str, float]],
    ) -> Dict[str, Dict[str, float]]:
        """DAE interface stub: speciation charge-balance residual.

        Activated in Phase F when ``DAEStepSolver`` is implemented.
        """
        raise NotImplementedError(
            "compute_algebraic_residual is a Phase F stub; "
            "activate when DAEStepSolver is implemented."
        )

    def compute_jacobian(self, t_h: float) -> Optional["np.ndarray"]:
        """Return Jacobian of compute_rhs w.r.t. the species vector, or None.

        Returns ``None`` until a model-specific override provides the
        analytical Jacobian — a stub extension point, not currently
        consumed by any shipped solver. The one shipped analytical-Jacobian
        path, ``SimultaneousAdaptiveSolver(use_engine_jacobian=True)``, uses
        ``engine.jacobian_dz_dy()`` (the speciation engine's own gray-box
        Jacobian, a distinct object from this CV-level RHS Jacobian) —
        not this method. Neither does ``MonolithicODESolver``, which
        never passes a ``jac=`` callable to its ``solve_ivp`` call.
        """
        return None

    def snapshot_state(self) -> Dict[str, Dict[str, float]]:
        """Return ``{phase_key: dict(phase.n_mol)}`` — deep copy of mole inventory.

        Use with :meth:`restore_state` for reversible trial-state evaluation::

            snap = cv.snapshot_state()
            # ... mutate cv state ...
            cv.restore_state(snap)
        """
        return {pk: dict(phase.n_mol) for pk, phase in self.phases.items()}

    def restore_state(self, snapshot: Dict[str, Dict[str, float]]) -> None:
        """Overwrite ``phase.n_mol`` from a snapshot taken by :meth:`snapshot_state`."""
        for pk, n_mol_dict in snapshot.items():
            if pk in self.phases:
                phase = self.phases[pk]
                phase.n_mol.clear()
                phase.n_mol.update(n_mol_dict)

    # ── Read-only reaction rate computation ──────────────────────────────

    def compute_reaction_rates(
        self,
        t_h: float = 0.0,
    ) -> Dict[str, Dict[str, float]]:
        """Compute reaction source terms from current state (read-only).

        Returns ``{phase_key: {species: mol/h}}`` without modifying any
        phase state.  Used by the snapshot-based advance pattern (see
        ``SimultaneousEulerSolver``).

        The :class:`ReactionEnvironment` is built from current
        ``phase.pH`` (derived from ``n_mol["H+"]``) and
        ``phase.n_mol``. Callers wanting a *pre-speciation* read must
        ensure the engine has already solved on the phase (the
        sequential ``cv.advance()`` body handles this by running
        speciation in step 1 before any rate evaluation).

        Parameters
        ----------
        t_h : float
            Current simulation time (hours). Forwarded into the
            :class:`ReactionEnvironment` for time-dependent rate
            laws. Default ``0.0``.

        Returns
        -------
        dict of dict
            ``{phase_key: {species_id: rate_mol_per_h}}``.
        """
        if self.reaction_system is None:
            return {}
        if not hasattr(self.reaction_system, "compute_rates"):
            return {}
        env = self._build_reaction_environment(t_h)
        rates = self.reaction_system.compute_rates(env)
        # Return a deep copy to ensure no aliasing
        return {pk: dict(sp_rates) for pk, sp_rates in rates.items()}

    # ── ConservationMonitor wiring (state-unification C6) ──────────────

    def _collect_species_registry(self) -> Dict[str, Any]:
        """Build ``{species_id: Species}`` from the reaction
        stoichiometries.

        Used by :meth:`__init__` to seed the
        :class:`ConservationMonitor`'s registry without requiring
        the user to declare species twice. Species not appearing
        in any reaction stoichiometry (e.g. unnamed strong-ion
        lumps ``S_cat`` / ``S_an``) are absent from the registry
        and silently skipped during conservation accounting.
        """
        registry: Dict[str, Any] = {}
        if self.reaction_system is None:
            return registry
        rxns = getattr(self.reaction_system, "reactions", None)
        if rxns is None:
            # Bare KineticReaction attached directly
            if hasattr(self.reaction_system, "stoichiometry"):
                rxns = [self.reaction_system]
            else:
                return registry
        for rxn in rxns:
            stoich = getattr(rxn, "stoichiometry", None)
            if stoich is None:
                continue
            for entry in stoich:
                sp = entry.species
                registry[sp.id] = sp
        # Supplement with any species from common_species found in phase.n_mol
        # (e.g. Cl-, Na+, K+ which never appear in reaction stoichiometry but
        # must be included for charge conservation accounting).
        _catalog = self._common_species_catalog()
        for phase in self.phases.values():
            for sp_id in phase.n_mol:
                if sp_id not in registry:
                    known = _catalog.get(sp_id)
                    if known is not None:
                        registry[sp_id] = known
        return registry

    @staticmethod
    def _common_species_catalog() -> Dict[str, Any]:
        """Return ``{id: Species}`` for every Species declared in
        ``PyOMES.chemistry.common_species``.  Lazy import avoids
        circular-import risk at module load time."""
        from ..chemistry import common_species as _cs
        from ..chemistry.species import Species
        return {obj.id: obj for obj in vars(_cs).values() if isinstance(obj, Species)}

    # ── Property calculators (state-unification C5) ────────────────────

    def _run_property_calculators(self) -> None:
        """Invoke each attached :class:`PropertyCalculator` on its
        target phase and store the result on ``phase.properties``.

        Runs once per ``advance()`` step, after speciation and
        before reactions. Each calculator declares its target via
        ``phase_key``; calculators whose target phase is absent are
        silently skipped. Pressure passed to ``compute`` is taken
        from the phase itself when it exposes ``P_atm`` (gas phases),
        otherwise from the gas phase if present, else 1.0 atm.

        Calculators are independent — order of evaluation should not
        matter. If a calculator raises, the exception propagates
        (no silent suppression; the user must register only
        calculators they trust).
        """
        gas = self.phases.get("gas")
        for calc in self.property_calculators:
            phase = self.phases.get(calc.phase_key)
            if phase is None:
                continue
            T_K = float(phase.T_K)
            if hasattr(phase, "P_atm"):
                P_atm = float(phase.P_atm)
            elif gas is not None:
                P_atm = float(gas.P_atm)
            else:
                P_atm = 1.0
            value = calc.compute(phase, T_K, P_atm)
            phase.properties[calc.key] = float(value)

    # ── Reaction integration helpers ───────────────────────────────────

    def _build_reaction_environment(
        self,
        t_h: float = 0.0,
    ) -> "ReactionEnvironment":
        """Build a ReactionEnvironment from current phase state.

        Reads concentrations from the liquid phase's ``n_mol`` (mol/L
        via division by ``V_L``). pH comes from ``phase.pH``, which
        derives from ``n_mol["H+"]`` (the speciation engine writes
        it on each solve). When the liquid phase has no ``H+`` entry
        (a CV without speciation), ``env.pH`` is ``None``; reaction
        models that require pH must handle this themselves via
        :attr:`ReactionEnvironment.has_pH`.

        Parameters
        ----------
        t_h : float
            Current simulation time (hours). Forwarded to the
            environment for time-dependent rate laws.

        Returns
        -------
        ReactionEnvironment
        """
        from ..reactions.environment import ReactionEnvironment

        # Liquid phase concentrations (mol/L)
        concentrations: Dict[str, float] = {}
        liq = self.phases.get("liquid")
        if liq is not None:
            V_L = float(liq.V_L)
            for sp, n in liq.n_mol.items():
                concentrations[sp] = float(n) / V_L
        else:
            gas = self.phases.get("gas")
            if gas is not None:
                V_L = float(gas.V_L)
            else:
                raise ValueError(
                    f"CV {self.label!r}: _build_reaction_environment requires "
                    f"at least one phase with a volume; neither 'liquid' nor "
                    f"'gas' phase is present."
                )

        # Temperature (prefer liquid, fall back to gas)
        T_K = 298.15
        for phase_key in ("liquid", "gas"):
            ph = self.phases.get(phase_key)
            if ph is not None:
                T_K = float(getattr(ph, "T_K", T_K))
                break

        # pH from phase.pH (raises if n_mol["H+"] missing). For CVs
        # without speciation this is the explicit no-pH case.
        pH = None
        if liq is not None and "H+" in liq.n_mol:
            try:
                pH = liq.pH
            except ValueError:
                pH = None

        # Forward scalar derived properties. C5: PropertyCalculator
        # outputs land on phase.properties; legacy IonicStrength is
        # still on phase.speciation until the speciation dict is
        # retired alongside chemistry-unification-3b.
        properties: Dict[str, float] = {}
        if liq is not None:
            properties.update(liq.properties)
            ionic = liq.speciation.get("IonicStrength")
            if ionic is not None and "ionic_strength" not in properties:
                properties["ionic_strength"] = float(ionic)

        return ReactionEnvironment(
            T_K=T_K,
            V_L=V_L,
            pH=pH,
            concentrations=concentrations,
            properties=properties,
            t_h=float(t_h),
        )

    def _integrate_reactions(
        self,
        dt_h: float,
        t_h: float = 0.0,
    ) -> Dict[str, Dict[str, float]]:
        """Thin wrapper around compute_reaction_rates + Euler apply.

        Kept for backward compatibility; advance() now inlines this
        logic directly.  Returns the net source terms (mol/h) applied.
        """
        rhs = self.compute_reaction_rates(t_h)
        sources: Dict[str, Dict[str, float]] = {}
        for pk, sp_rates in rhs.items():
            sources[pk] = {}
            for sp, rate in sp_rates.items():
                phase = self.phases[pk]
                phase.n_mol[sp] = max(
                    0.0, phase.n_mol.get(sp, 0.0) + rate * dt_h
                )
                sources[pk][sp] = rate
        return sources

    def _build_species_vector_map(self):
        """Identify all species the reaction model touches.

        Evaluates the reaction model once with zero rates to discover
        which (phase_key, species_id) pairs it produces source terms for.
        Also includes any species already present in those phases.

        The result is cached keyed by ``id(self._reaction_system)`` and
        cleared when :meth:`_set_reaction_system_unchecked` is called.
        Species membership is invariant for the life of a
        :class:`~PyOMES.reactions.ReactionSystem` (lifecycle gating blocks
        mid-run mutation), so the cache is valid for the full run.

        Returns
        -------
        species_order : list of (phase_key, species_id)
            Ordered list defining the state vector layout.
        phase_map : dict
            ``{(phase_key, species_id): index}`` for quick lookup.
        """
        if self.reaction_system is None:
            return [], {}

        cache_key = id(self._reaction_system)
        if (self._species_vector_map_cache is not None and
                self._species_vector_map_cache_rs_id == cache_key):
            return self._species_vector_map_cache

        # Discover species from a dummy evaluation
        from ..reactions.environment import ReactionEnvironment
        dummy_env = ReactionEnvironment(T_K=298.15, V_L=1.0)
        try:
            dummy_sources = self.reaction_system.compute_rates(dummy_env)
        except (TypeError, ValueError, AttributeError, KeyError):
            dummy_sources = {}

        # Collect all (phase, species) pairs
        pairs = set()
        for pk, sp_dict in dummy_sources.items():
            if pk in self.phases:
                for sp in sp_dict:
                    pairs.add((pk, sp))
        # Also include species already in those phases
        for pk in {p for p, _ in pairs}:
            for sp in self.phases[pk].n_mol:
                pairs.add((pk, sp))

        species_order = sorted(pairs)
        phase_map = {pair: idx for idx, pair in enumerate(species_order)}
        self._species_vector_map_cache = (species_order, phase_map)
        self._species_vector_map_cache_rs_id = cache_key
        return species_order, phase_map

    # ── Convenience accessors ──────────────────────────────────────────

    def species_concentration(
        self, phase_key: str, species_id: str, default: float = 0.0
    ) -> float:
        """Return the concentration (mol/L) of a species in a phase.

        Parameters
        ----------
        phase_key : str
            Phase key (e.g. ``"liquid"``).
        species_id : str
            Species identifier.
        default : float
            Value returned if the species or phase is not found.

        Returns
        -------
        float
            Concentration in mol/L.
        """
        phase = self.phases.get(phase_key)
        if phase is None:
            return default
        V = float(phase.V_L)
        return float(phase.n_mol.get(species_id, 0.0)) / V

    # ── Snapshot ───────────────────────────────────────────────────────

    def snapshot(self) -> "ControlVolume":
        """Return an independent deep copy of this CV and all its phases."""
        return ControlVolume(
            phases={k: p.snapshot() for k, p in self.phases.items()},
            internal_interfaces=list(self._explicit_interfaces),  # shared (stateless)
            transfer_models=self.transfer_models,  # factory rebuilds link fresh
            boundaries=list(self.boundaries),  # shared (stateless)
            reaction_system=self.reaction_system,  # shared (stateless)
            property_calculators=list(self.property_calculators),  # shared (stateless)
            label=self.label,
        )

    # ── Initialization utilities ─────────────────────────────────────────────

    def equilibrate_to_pH(
        self,
        corrector_id: str,
        ph_target: float,
        *,
        liquid_key: str = "liquid",
        max_add_mol_per_L: float = 1.0,
        tol_pH: float = 1e-5,
    ) -> float:
        """Add *corrector_id* to the liquid until *ph_target* is reached.

        The search is evaluated through :meth:`advance` at ``dt_h=0`` so
        all multi-phase equilibria (gas–liquid partitioning etc.) are
        accounted for alongside liquid speciation.

        Strong-corrector shorthand: ``"NaOH"`` and ``"KOH"`` are recognised
        without requiring dissolution reactions in the database — the method
        adds the corresponding strong cation (``Na+`` / ``K+``) to
        ``n_mol``, and the charge-balance solver automatically adjusts
        ``OH-`` / ``H+`` to compensate.

        Parameters
        ----------
        corrector_id:
            Species id to add.  Either a recognised strong corrector
            (``"NaOH"``, ``"KOH"``), or a species that appears in at least
            one :class:`~PyOMES.reactions.EquilibriumReaction` stoichiometry.
        ph_target:
            Target pH.
        liquid_key:
            Key of the liquid phase in :attr:`phases`.
        max_add_mol_per_L:
            Search upper bound in mol per litre of liquid.
        tol_pH:
            Convergence tolerance on pH.

        Returns
        -------
        float
            Moles added (always ≥ 0).

        Raises
        ------
        ValueError
            If *corrector_id* is not a known strong corrector or in any
            equilibrium stoichiometry, shifts pH in the wrong direction, or
            the target is unreachable within *max_add_mol_per_L*.
        """
        from scipy.optimize import brentq  # lazy — not needed at module load

        liq   = self.phases[liquid_key]
        V     = float(liq.V_L)
        max_n = max_add_mol_per_L * V

        # Resolve the species actually added to n_mol. Strong correctors like
        # NaOH are represented by their charge-carrying ion (Na+) so that no
        # dissolution reaction is needed in the database.
        _dose_id = _STRONG_CORRECTOR_ION.get(corrector_id, corrector_id)

        # Suppress the ConservationMonitor for the duration of the probe and
        # brentq search: every internal advance adds an unbalanced strong ion
        # (e.g. Na+) to n_mol, which would trigger spurious charge/element
        # warnings.  The monitor is restored and re-baselined to the corrected
        # initial state in a finally block so simulation drift is tracked
        # cleanly from there regardless of whether the correction succeeds.
        _saved_monitor = self._conservation_monitor
        self._conservation_monitor = None
        try:
            # 1. Validate: corrector must participate in at least one equilibrium
            #    (skipped for known strong correctors whose ion is always tracked)
            if self.reaction_system is None:
                raise ValueError("CV has no reaction_system; cannot evaluate pH.")
            if corrector_id not in _STRONG_CORRECTOR_ION:
                eq_species = {
                    e.species.id
                    for rxn in self.reaction_system.single_phase_equilibria
                    for e in rxn.stoichiometry
                }
                if corrector_id not in eq_species:
                    raise ValueError(
                        f"{corrector_id!r} does not appear in any EquilibriumReaction "
                        f"stoichiometry and therefore cannot shift pH through speciation.\n"
                        f"Species covered by equilibria: {sorted(eq_species)}"
                    )

            # 2. Baseline pH via zero-timestep advance
            snap0 = self.snapshot_state()
            self.advance(dt_h=0.0, t_h=0.0)
            pH_base = float(liq.pH)
            self.restore_state(snap0)

            if abs(pH_base - ph_target) < tol_pH:
                self.advance(dt_h=0.0, t_h=0.0)
                print(f"pH already at target ({pH_base:.4f}); equilibrated in place.")
                return 0.0

            # 3. Probe direction: tiny addition reveals which way pH moves
            eps = max_n * 1e-6
            snap0 = self.snapshot_state()
            liq.n_mol[_dose_id] = liq.n_mol.get(_dose_id, 0.0) + eps
            self.advance(dt_h=0.0, t_h=0.0)
            pH_probe = float(liq.pH)
            self.restore_state(snap0)

            dpH = pH_probe - pH_base
            if abs(dpH) < 1e-12:
                raise ValueError(
                    f"Adding {corrector_id!r} produced no measurable pH change at "
                    "the current composition. The species may be fully buffered or "
                    "decoupled from the proton balance at this pH."
                )

            needs_rise     = ph_target > pH_base
            species_raises = dpH > 0
            if needs_rise != species_raises:
                direction = "raises" if species_raises else "lowers"
                needed    = "raise"  if needs_rise    else "lower"
                raise ValueError(
                    f"Adding {corrector_id!r} {direction} pH "
                    f"(baseline {pH_base:.3f} → probe {pH_probe:.3f}), "
                    f"but reaching ph_target={ph_target:.3f} requires a {needed}. "
                    "Choose a corrector that moves pH in the correct direction."
                )

            # 4. Bracket and root-find with brentq
            def _residual(n_add: float) -> float:
                snap = self.snapshot_state()
                liq.n_mol[_dose_id] = liq.n_mol.get(_dose_id, 0.0) + n_add
                self.advance(dt_h=0.0, t_h=0.0)
                pH = float(liq.pH)
                self.restore_state(snap)
                return pH - ph_target

            r_hi = _residual(max_n)
            if (r_hi * (pH_base - ph_target)) > 0:
                raise ValueError(
                    f"Target pH {ph_target:.3f} not reached within "
                    f"{max_add_mol_per_L:.2f} mol/L of {corrector_id!r} "
                    f"(pH at max dose: {ph_target + r_hi:.3f}). "
                    "Increase max_add_mol_per_L."
                )

            n_opt = brentq(_residual, 0.0, max_n, xtol=1e-12, rtol=1e-10)

            # 5. Apply optimal amount and leave CV in equilibrated state.
            liq.n_mol[_dose_id] = liq.n_mol.get(_dose_id, 0.0) + n_opt
            self.advance(dt_h=0.0, t_h=0.0)
            print(
                f"pH correction: added {n_opt * 1000:.4f} mmol {corrector_id!r} "
                f"({n_opt / V * 1000:.4f} mmol/L).\n"
                f"pH: {pH_base:.4f} → {float(liq.pH):.4f}  (target {ph_target:.4f})"
            )
            return n_opt
        finally:
            # Restore the monitor and rebias it to the corrected state so the
            # subsequent simulation tracks real drift from the new baseline.
            self._conservation_monitor = _saved_monitor
            if _saved_monitor is not None:
                _saved_monitor.reset()

    def __repr__(self):
        phases_str = ", ".join(f"{k}: {type(v).__name__}" for k, v in self.phases.items())
        n_iface = len(self.internal_interfaces)
        n_bounds = len(self.boundaries)
        has_rxn = self.reaction_system is not None
        return (f"ControlVolume(label={self.label!r}, "
                f"phases={{{phases_str}}}, "
                f"interfaces={n_iface}, "
                f"boundaries={n_bounds}, "
                f"reaction_system={has_rxn})")
