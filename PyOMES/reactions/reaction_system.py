# -*- coding: utf-8 -*-
"""ReactionSystem — single attach point for all reactions on a ControlVolume.

A :class:`ReactionSystem` holds one or more reaction declarations
(:class:`~PyOMES.reactions.kinetic.KineticReaction`,
:class:`~PyOMES.reactions.blackbox.BlackBoxReactionModel`, or anything
satisfying :class:`~PyOMES.reactions.equilibrium.EquilibriumConstraint`
— :class:`~PyOMES.reactions.equilibrium.EquilibriumReaction`,
:class:`~PyOMES.chemistry.partition.HenryEquilibrium`,
:class:`~PyOMES.chemistry.partition.KspEquilibrium`,
:class:`~PyOMES.chemistry.partition.RaoultEquilibrium`) and
**pre-buckets them by classification at construction time**. The
buckets are internal; consumers read either the unified list
(``system.reactions``) or the type-specific projections via public
properties (e.g. ``system.cross_phase_equilibria`` for the gas-liquid
link).

Replaces the deleted ``ReactionSet`` class. Differences from the
previous container:

- **Single attach point.** ``cv.reaction_system`` is the only reaction
  attachment on a CV. ``cv.reaction_model`` is gone.
- **Internal pre-bucketing.** Reactions are sorted into
  ``_kinetic_reactions``, ``_single_phase_equilibria``,
  ``_cross_phase_equilibria``, ``_precipitation_equilibria``, and
  ``_blackbox_models`` at ``__init__`` — ``KineticReaction`` and
  ``BlackBoxReactionModel`` by ``isinstance``, everything else via
  :func:`~PyOMES.reactions.equilibrium.classify_equilibrium_constraint`
  (EQUILIBRIUM_CONSTRAINT_UNIFICATION CP3) rather than a hard-coded
  ``isinstance(rxn, EquilibriumReaction)`` check — so a single
  ``HenryEquilibrium`` instance can be constructed once and attached
  here directly, in the same list as a ``KineticGasLiquidLink``'s
  ``partition_models=`` dict, with no separate declaration to keep in
  sync. Downstream callers read the bucket directly instead of
  partitioning on demand.
- **No public ``partition()`` or ``has_equilibrium`` API.** Partition
  is internal; mixed-set attach is always correct because the
  ``compute_rates`` path only iterates the kinetic + blackbox buckets.
- **Immutable post-attach.** No ``add()`` / ``remove()`` methods.
  Mutation requires constructing a new ``ReactionSystem``.

The speciation engine is constructed lazily on first need and cached
on ``system._engine`` (introduced in checkpoint 3 of the
``state-unification`` phase, when the engine starts writing derived
species back to ``n_mol`` directly). In C2 the existing
``SpeciationPropertySolver`` in ``cv.property_solvers`` still owns the
user-facing engine; ``system._engine`` is reserved for the C3+ flow.
"""

from __future__ import annotations

import warnings
from typing import Dict, List, Optional, Sequence, Union

from .equilibrium import EquilibriumConstraint, classify_equilibrium_constraint
from .kinetic import KineticReaction
from .blackbox import BlackBoxReactionModel
from .environment import ReactionEnvironment
from ._shared import fmt_stoichiometry_string


ReactionLike = Union[KineticReaction, EquilibriumConstraint, BlackBoxReactionModel]

_KNOWN_SOLVERS = frozenset({"charge_balance", "newton_raphson"})


class ReactionSystem:
    """A collection of reactions attached to a single ControlVolume.

    Parameters
    ----------
    reactions : sequence of KineticReaction / EquilibriumConstraint /
                BlackBoxReactionModel
        All reactions in this system — ``EquilibriumConstraint`` covers
        ``EquilibriumReaction``, ``HenryEquilibrium``, ``KspEquilibrium``,
        and ``RaoultEquilibrium`` uniformly. Pre-bucketed at construction
        (``KineticReaction``/``BlackBoxReactionModel`` by ``isinstance``,
        equilibrium constraints via
        :func:`~PyOMES.reactions.equilibrium.classify_equilibrium_constraint`).
        Each reaction must already be stoichiometrically validated (the
        declaration classes do this at their own construction).
    label : str
        Human-readable label (optional).

    Notes
    -----
    A ``ReactionSystem`` can hold kinetic and equilibrium reactions
    in one container — this is the canonical attach point on a
    :class:`~PyOMES.core.control_volume.ControlVolume`. The internal
    bucketing routes each reaction to the right consumer:

    - kinetic reactions and black-box models contribute to
      :meth:`compute_rates`,
    - single-phase equilibria are routed to the speciation engine,
    - cross-phase equilibria are exposed via
      :attr:`cross_phase_equilibria` for consumption by
      :class:`~PyOMES.core.gas_liquid_link.KineticGasLiquidLink`.

    Attempting to mix kinds was an error in the old ``ReactionSet``
    (``compute_rates`` raised on equilibria); the pre-bucketing here
    makes that impossible — equilibria are simply not in the kinetic
    iteration path.
    """

    def __init__(
        self,
        reactions: Sequence[ReactionLike],
        *,
        label: str = "",
        solver: str = "charge_balance",
    ):
        self.reactions: List[ReactionLike] = list(reactions)
        self.label = str(label)

        solver = str(solver)
        if solver not in _KNOWN_SOLVERS:
            raise ValueError(
                f"ReactionSystem: unknown solver {solver!r}. "
                f"Valid values: {sorted(_KNOWN_SOLVERS)}. "
                f"Pass solver= to __init__; it cannot be changed after construction."
            )
        self._solver = solver

        # Pre-bucket — the system maintains these as private read-only
        # projections of the source list. Mutation would require
        # constructing a new ReactionSystem (immutable post-attach by
        # design). KineticReaction/BlackBoxReactionModel by isinstance;
        # everything else via classify_equilibrium_constraint() so any
        # EquilibriumConstraint sibling (HenryEquilibrium, KspEquilibrium,
        # RaoultEquilibrium) is recognised uniformly with
        # EquilibriumReaction — not a hard-coded isinstance check.
        self._kinetic_reactions: List[KineticReaction] = []
        self._single_phase_equilibria: List[EquilibriumConstraint] = []
        self._cross_phase_equilibria: List[EquilibriumConstraint] = []
        self._precipitation_equilibria: List[EquilibriumConstraint] = []
        self._blackbox_models: List[BlackBoxReactionModel] = []

        for rxn in self.reactions:
            if isinstance(rxn, KineticReaction):
                self._kinetic_reactions.append(rxn)
            elif isinstance(rxn, BlackBoxReactionModel):
                self._blackbox_models.append(rxn)
            elif isinstance(rxn, EquilibriumConstraint):
                try:
                    kind = classify_equilibrium_constraint(rxn)
                except ValueError as exc:
                    raise TypeError(
                        f"ReactionSystem: {type(rxn).__name__} entry {rxn!r} "
                        f"could not be classified: {exc}"
                    ) from exc
                if kind == "solid_liquid":
                    self._precipitation_equilibria.append(rxn)
                elif kind == "gas_liquid":
                    self._cross_phase_equilibria.append(rxn)
                else:
                    self._single_phase_equilibria.append(rxn)
            else:
                raise TypeError(
                    f"ReactionSystem: unrecognised reaction kind "
                    f"{type(rxn).__name__} for entry {rxn!r}. "
                    f"Expected KineticReaction, an EquilibriumConstraint-"
                    f"conforming type (EquilibriumReaction, HenryEquilibrium, "
                    f"KspEquilibrium, RaoultEquilibrium, ...), or "
                    f"BlackBoxReactionModel."
                )

        # Raise immediately on hard species conflicts (same id, different
        # atoms/charge/MW across reactions). Soft conflicts (same data,
        # distinct objects — e.g. builder-created vs database instances)
        # are silently accepted here; ControlVolume also ignores them.
        from ..chemistry.species_check import check_species_consistency
        check_species_consistency(self.reactions, soft_conflicts="ignore")

        self._engine = None
        self._engine_config: Dict[str, object] = {
            "use_activity": False,
            "activity_model": "davies",
        }
        # AccuracyMonitor reference, propagated to the engine on
        # first build (or immediately if already built).
        # ``cv.reaction_system.attach_monitor(monitor)`` is the
        # canonical attachment path post-C4.
        self._accuracy_monitor = None
        # ConservationMonitor reference (state-unification C6).
        # ``cv.reaction_system.attach_conservation_monitor(monitor)``
        # is the canonical attachment path; the CV's advance() body
        # invokes ``monitor.check_step(self.phases)`` once per step.
        self._conservation_monitor = None

    # ── Public bucket projections ─────────────────────────────────────

    @property
    def kinetic_reactions(self) -> List[KineticReaction]:
        """Read-only view of the kinetic reactions in this system."""
        return list(self._kinetic_reactions)

    @property
    def single_phase_equilibria(self) -> List[EquilibriumConstraint]:
        """Read-only view of the single-phase (acid-base) equilibrium
        constraints in this system (the ones consumed by
        :meth:`BisectionChemicalEquilibriumEngine.from_reactions`/
        :meth:`NRChemicalEquilibriumEngine.from_reactions`)."""
        return list(self._single_phase_equilibria)

    @property
    def cross_phase_equilibria(self) -> List[EquilibriumConstraint]:
        """Read-only view of the gas-liquid equilibrium constraints in
        this system — partition declarations consumed by
        :meth:`~PyOMES.core.gas_liquid_link.KineticGasLiquidLink.derive_speciation_keys`.
        Includes both cross-phase ``EquilibriumReaction`` declarations
        and ``HenryEquilibrium``/``RaoultEquilibrium`` instances."""
        return list(self._cross_phase_equilibria)

    @property
    def precipitation_equilibria(self) -> List[EquilibriumConstraint]:
        """Read-only view of solid-liquid equilibrium constraints in this system.

        A constraint is placed here when its stoichiometry has any
        :class:`~PyOMES.reactions.stoichiometry.StoichiometryEntry` with
        ``phase="solid"`` — i.e. an ``EquilibriumReaction`` or
        ``KspEquilibrium`` mineral dissolution/precipitation constraint,
        with ``log_K = log10(Ksp)`` (solid activity is 1 by convention).
        Folded into the flat list passed to
        :meth:`NRChemicalEquilibriumEngine.from_reactions` by :attr:`engine`,
        which auto-classifies them back out into its active-set loop.
        """
        return list(self._precipitation_equilibria)

    @property
    def blackbox_models(self) -> List[BlackBoxReactionModel]:
        """Read-only view of the black-box reaction models in this
        system (opaque external simulators contributing source
        terms via :meth:`compute_rates`)."""
        return list(self._blackbox_models)

    # ── Speciation-engine attachment ──────────────────────────────────

    @property
    def engine(self):
        """Lazily build (or return) the speciation engine.

        Triggered on first access. The engine is constructed via
        :meth:`BisectionChemicalEquilibriumEngine.from_reactions`/
        :meth:`NRChemicalEquilibriumEngine.from_reactions` using the system's
        equilibrium constraints, plus whatever defaults are pinned in
        :attr:`_engine_config` (set by :meth:`configure_engine`). The
        ``newton_raphson`` solver path also includes
        :attr:`precipitation_equilibria` — ``NRChemicalEquilibriumEngine`` auto-
        classifies solid-liquid items back out of the flat list
        (EQUILIBRIUM_CONSTRAINT_UNIFICATION CP2) — since ``charge_balance``
        (:class:`~PyOMES.chemical_equilibrium.engine.BisectionChemicalEquilibriumEngine`) has no
        precipitation support, that path omits them.

        Returns ``None`` when the system declares no equilibria —
        callers must check before solving.
        """
        if self._engine is not None:
            return self._engine

        if self._solver == "newton_raphson":
            equilibria = (
                list(self._single_phase_equilibria)
                + list(self._cross_phase_equilibria)
                + list(self._precipitation_equilibria)
            )
            if not equilibria:
                return None
            from ..chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
            self._engine = NRChemicalEquilibriumEngine.from_reactions(
                equilibria,
                use_activity=self._engine_config["use_activity"],
                activity_model=self._engine_config["activity_model"],
            )
        else:
            equilibria = list(self._single_phase_equilibria) + list(
                self._cross_phase_equilibria
            )
            if not equilibria:
                return None
            from ..chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
            self._engine = BisectionChemicalEquilibriumEngine.from_reactions(
                equilibria,
                use_activity=self._engine_config["use_activity"],
                activity_model=self._engine_config["activity_model"],
            )
        # Propagate any pre-attached AccuracyMonitor to the
        # freshly-built engine.
        if self._accuracy_monitor is not None:
            self._engine._accuracy_monitor = self._accuracy_monitor
        return self._engine

    def configure_engine(
        self,
        *,
        use_activity: bool = None,
        activity_model: str = None,
    ) -> None:
        """Override engine-construction defaults before first build.

        Call this before any :meth:`engine` access or
        :meth:`advance` step. Tweaks the lazy-build settings; the
        engine is still constructed on first access.

        Parameters
        ----------
        use_activity : bool, optional
            Apply activity-coefficient corrections (Davies / SIT).
            Default ``False``.
        activity_model : str, optional
            ``"davies"`` (default), ``"sit"``, or ``"ideal"``.

        .. warning::
            The solver backend (``"charge_balance"`` vs
            ``"newton_raphson"``) **cannot** be changed here — it must
            be set at construction time via
            ``ReactionSystem(..., solver="newton_raphson")``.
            Calling :meth:`configure_engine` after the engine has been
            built (i.e. after the first :attr:`engine` access or
            :meth:`advance` call) raises :exc:`RuntimeError`.

        Raises
        ------
        RuntimeError
            If the engine has already been built (lazy access
            already triggered). Configure before the first solve.
        """
        if self._engine is not None:
            raise RuntimeError(
                "ReactionSystem.configure_engine: engine has already "
                "been built. Configure before the first .engine access "
                "or .advance() call, or call .attach_engine(...) "
                "with a fully-constructed replacement."
            )
        if use_activity is not None:
            self._engine_config["use_activity"] = bool(use_activity)
        if activity_model is not None:
            self._engine_config["activity_model"] = str(activity_model)

    def attach_engine(self, engine) -> None:
        """Attach a pre-constructed speciation engine.

        Replaces any lazy-built engine. Use this when the default
        ``BisectionChemicalEquilibriumEngine.from_reactions`` construction doesn't
        match what you need (e.g. for tests injecting a stub, or
        for non-standard chemistry configurations).
        """
        self._engine = engine
        if self._accuracy_monitor is not None:
            self._engine._accuracy_monitor = self._accuracy_monitor

    def attach_monitor(self, monitor) -> None:
        """Attach an :class:`AccuracyMonitor` to this system's engine.

        The monitor is propagated to :attr:`_engine` immediately if
        the engine is already built, or stored for propagation on
        first build otherwise. Replaces the pre-C4 attachment path
        that hung the monitor off ``cv.property_solvers[*]._accuracy_monitor``.
        """
        self._accuracy_monitor = monitor
        if self._engine is not None:
            self._engine._accuracy_monitor = monitor

    def attach_conservation_monitor(self, monitor) -> None:
        """Attach a :class:`ConservationMonitor` for element + charge
        accounting.

        The monitor is invoked once per
        :meth:`ControlVolume.advance` step (called from the CV's
        sequential body, parallel to the speciation-side
        AccuracyMonitor hooks). The CV populates the monitor's
        species registry from the system's reaction stoichiometries
        when this method is called from
        :meth:`ControlVolume.__init__`; direct callers must call
        :meth:`monitor.set_species_registry(...)` themselves if
        the system's reactions don't cover all n_mol species.

        Introduced in state-unification C6.
        """
        self._conservation_monitor = monitor

    # ── ReactionModel protocol ────────────────────────────────────────

    def compute_rates(
        self, env: ReactionEnvironment
    ) -> Dict[str, Dict[str, float]]:
        """Sum source terms across kinetic constituents and black-box
        models.

        Equilibrium constituents (both single-phase and cross-phase)
        are silently skipped — they have no rate by construction and
        are routed elsewhere by the pre-bucketing. Mixed-kind
        systems are always safe to integrate.

        Parameters
        ----------
        env : ReactionEnvironment
            Current local conditions (shared by all reactions).

        Returns
        -------
        dict of dict
            ``{phase_key: {species_id: mol_per_h}}``.
        """
        combined: Dict[str, Dict[str, float]] = {}
        for rxn in self._kinetic_reactions:
            terms = rxn.compute_rates(env)
            for phase, species_rates in terms.items():
                if phase not in combined:
                    combined[phase] = {}
                for sp, rate in species_rates.items():
                    combined[phase][sp] = combined[phase].get(sp, 0.0) + rate
        for rxn in self._blackbox_models:
            terms = rxn.compute_rates(env)
            for phase, species_rates in terms.items():
                if phase not in combined:
                    combined[phase] = {}
                for sp, rate in species_rates.items():
                    combined[phase][sp] = combined[phase].get(sp, 0.0) + rate
        return combined

    # ── Aggregation helpers ──────────────────────────────────────────

    @property
    def species_ids(self) -> list:
        """Sorted list of all species IDs across all reactions."""
        ids = set()
        for rxn in self.reactions:
            species_ids = getattr(rxn, "species_ids", None)
            if species_ids is not None:
                ids.update(species_ids)
        return sorted(ids)

    @property
    def phases(self) -> list:
        """Sorted list of all phase keys across all reactions."""
        phases = set()
        for rxn in self.reactions:
            rxn_phases = getattr(rxn, "phases", None)
            if rxn_phases is not None:
                phases.update(rxn_phases)
        return sorted(phases)

    def plot_speciation(
        self,
        anchor_id: str,
        *,
        pH_range: tuple = (0.0, 14.0),
        T_K: Optional[float] = 298.15,
        n_points: int = 500,
        ax=None,
        title: Optional[str] = None,
    ) -> tuple:
        """Plot molar-fraction vs pH for the acid-base ladder containing *anchor_id*.

        Scans :attr:`single_phase_equilibria` for reactions of the form
        ``acid ⇌ base + H⁺``, builds the connected chain automatically
        from *anchor_id*, and plots α fractions analytically.

        Delegates to :func:`PyOMES.reactions.plots.plot_speciation`.
        Matplotlib is imported lazily.

        Parameters
        ----------
        anchor_id : str
            Species ID of any member of the target ladder (e.g. ``"CO2"``,
            ``"AceticAcid"``).
        pH_range : (float, float)
            pH axis limits.  Default ``(0, 14)``.
        T_K : float or None
            Temperature for Van 't Hoff correction.  ``None`` uses stored
            ``log_K`` values.  Default ``298.15`` K.
        n_points : int
            Number of pH points.  Default 500.
        ax : matplotlib Axes, optional
            Axes to draw on; a new figure is created when ``None``.
        title : str, optional
            Axes title.

        Returns
        -------
        fig, ax
            The matplotlib Figure and primary Axes objects.
        """
        from .plots import plot_speciation
        return plot_speciation(
            self,
            anchor_id,
            pH_range=pH_range,
            T_K=T_K,
            n_points=n_points,
            ax=ax,
            title=title,
        )

    def show_reactions(self, n: Optional[int] = None) -> None:
        """Print reactions in stoichiometry-string format.

        Parameters
        ----------
        n : int, optional
            Number of reactions to show.  ``None`` (default) prints all.
        """
        reactions = self.reactions[:n] if n is not None else self.reactions
        skipped = 0
        for rxn in reactions:
            if isinstance(rxn, BlackBoxReactionModel):
                skipped += 1
                continue
            # KineticReaction is the only "->" family here; every
            # EquilibriumConstraint sibling (EquilibriumReaction,
            # HenryEquilibrium, KspEquilibrium, RaoultEquilibrium) is "<->".
            arrow = "->" if isinstance(rxn, KineticReaction) else "<->"
            print(f"{rxn.label or type(rxn).__name__}:")
            print(f"  {fmt_stoichiometry_string(rxn.stoichiometry, arrow=arrow)}")
        if skipped:
            print(
                f"({skipped} BlackBoxReactionModel(s) skipped"
                " — no explicit stoichiometry)"
            )

    def show_balance(self, *, elements: Optional[Sequence[str]] = None) -> None:
        """Print per-element residuals for every reaction in the system.

        Calls ``.show_balance()`` on each kinetic reaction and each
        equilibrium constraint (single-phase, cross-phase, precipitation)
        that implements it — ``EquilibriumReaction``/``KineticReaction``
        do; ``HenryEquilibrium``/``KspEquilibrium``/``RaoultEquilibrium``
        don't (no elemental balance to validate — they're parameter
        conversions, not user-declared stoichiometry strings) and are
        skipped, like :class:`BlackBoxReactionModel`.

        Parameters
        ----------
        elements : sequence of str, optional
            Elements to check.  Defaults to every element found across
            the participants of each individual reaction.
        """
        buckets = (
            self._kinetic_reactions
            + self._single_phase_equilibria
            + self._cross_phase_equilibria
            + self._precipitation_equilibria
        )
        skipped = 0
        for rxn in buckets:
            if not hasattr(rxn, "show_balance"):
                skipped += 1
                continue
            rxn.show_balance(elements=elements)
            print()
        skipped += len(self._blackbox_models)
        if skipped:
            print(
                f"({skipped} reaction(s) skipped — no show_balance() "
                "(BlackBoxReactionModel and PartitionModel-derived "
                "EquilibriumConstraint types have no elemental balance to validate)"
            )

    def __len__(self) -> int:
        return len(self.reactions)

    def __iter__(self):
        return iter(self.reactions)

    def __repr__(self) -> str:
        lbl = f" [{self.label}]" if self.label else ""
        return (
            f"ReactionSystem("
            f"{len(self._kinetic_reactions)} kinetic, "
            f"{len(self._single_phase_equilibria)} single-phase eq, "
            f"{len(self._cross_phase_equilibria)} cross-phase eq, "
            f"{len(self._precipitation_equilibria)} precipitation eq, "
            f"{len(self._blackbox_models)} black-box, "
            f"solver={self._solver!r}"
            f"{lbl})"
        )
