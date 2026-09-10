# State Unification Refactor — Design Note

> **Status: Shipped 2026-05-22** on the `state-unification`
> branch. Tag: `state-unification-shipped`. Final suite: 793 / 0 in
> 137s. See the companion
> [STATE_UNIFICATION_CHECKLIST.md](STATE_UNIFICATION_CHECKLIST.md)
> for the per-checkpoint implementation log (C1–C6, plus the C4
> sub-commits 4a–4e), the sub-decision deltas, and the BSM2
> sentinel re-baselines (pH 3.340171 → 3.303580 across the phase;
> max drift inside model tolerances).
>
> **Design resolved:** 2026-05-20. Scheduled after
> `INTEGRATOR_REMOVAL`, before `SIMULATION_CLASS`. One feature branch,
> six checkpoints (C1–C6). Scope locked.

## Context

The chemistry-unification track (phases 1–4, shipped through
2026-05-16) unified the **declaration** of kinetic and equilibrium
reactions through a single `Reaction(kind=...)` class. The track
deliberately preserved the runtime split — kinetic reactions on
`cv.reaction_model`, equilibrium reactions routed via
`SpeciationEngine.from_reactions()` and exposed through
`cv.property_solvers[...]`. The split is defended on stiffness
grounds in
[`SPECIATION_REACTIONMODEL_BOUNDARY.md`](../design/SPECIATION_REACTIONMODEL_BOUNDARY.md);
that defence holds.

What didn't hold under scrutiny was the *type machinery* the track
left on top of the split:

- [`PropertySolver`](../../src/core/property_solver.py) protocol with
  one implementer (`SpeciationPropertySolver`).
- [`PropertyResult`](../../src/core/property_solver.py) dataclass
  carrying `pH`, `ionic_strength`, `species`, `alphas`, `logH_guess`,
  `diagnostics`, `raw`.
- [`SpeciationPropertySolver`](../../src/core/speciation_solver.py)
  adapter wrapping the engine.
- `ReactionSet.partition()` public API + `has_equilibrium` warning
  at CV attach time.
- `chem_env` dict carrying `t_h`, `T_K`, `strong_kwargs`,
  `logH_guess`.
- B2 resolution: "equilibrium outputs are properties, not state" —
  storage separation between `phase.n_mol` (conserved totals) and
  `PropertyResult` (derived species).

Walking the codebase against PHREEQC's reference design surfaced that
this entire type surface is **speciation-shaped infrastructure with
one implementer and growing complexity**, not load-bearing on the
stiffness argument. The "future viscosity / density solvers" defence
for `PropertySolver` is six months old with zero concrete demand
materialising. B2 is a particular operator-splitting choice; PHREEQC
ships unified storage with element-total conservation and has done so
for 30+ years.

This refactor collapses the type machinery by **unifying state on
`Phase`**, while preserving operator splitting as the default
integration strategy.

## Resolved design

### Storage — unified `n_mol`, no B2 separation

- `Phase.n_mol` holds **all species** in one dict: conserved totals
  (`TIC`, `CT_NH_T`, strong ions) AND equilibrium-derived (`H⁺`,
  `OH⁻`, `HCO3⁻`, `CO3²⁻`, `NH3`, `NH4⁺`, …).
- No `_conserved` / `_derived` marker sets; no protected-write
  protocol. Apply-flux can write to any species. Engine redistributes
  via element-total conservation after each operation. PHREEQC-style.
- `Phase.pH` is a `@property` computed as
  `-log10(n_mol["H+"] / V_L)`. Not stored. Single source of truth.
- `Phase.properties: dict` (or named fields) for genuine scalar
  derived properties — viscosity, density. Populated by
  `PropertyCalculator`s.

### Reactions — independent classes, no shared base

- `Reaction(kind=...)` deleted. Replaced by three **independent
  concrete classes** (no ABC, no shared base):
  - `KineticReaction` — single reaction with `rate_fn` and validated
    stoichiometry.
  - `EquilibriumReaction` — single reaction with `log_K`,
    `dH_J_per_mol`, `T_ref_K`, `is_cross_phase`, `total_id`.
  - `BlackBoxReactionModel` — opaque wrapper for an external
    network solver (FBA, GSM); produces source terms for many
    species from one solve.
- Shared validation (atom/charge closure, species_ids extraction,
  phase enumeration) lives as free functions in
  `src/reactions/_shared.py`.
- Type hints use `Union[KineticReaction, EquilibriumReaction]` or
  `Union[KineticReaction, EquilibriumReaction, BlackBoxReactionModel]`
  where the broader set is needed.
- Hard break — no `kind="..."` deprecation alias. Builders migrate
  in C1.
- **Future option (deferred, not in scope):** a `ReactionDeclaration`
  Protocol layer could be added later for flexibility — structural
  typing would let third-party reaction implementers satisfy the
  declaration contract without subclassing or modifying VLsim's
  class set. Useful if a plugin / extension API surface emerges. See
  *Open items deferred* below for the trigger conditions to revisit.

### Container — `ReactionSystem` replaces `ReactionSet`

- `cv.reaction_system: ReactionSystem` is the **single attach point**
  for all reactions on a CV.
- `ReactionSystem` pre-buckets its reactions by `isinstance` at
  construction time:
  - `system._kinetic_reactions: list[KineticReaction]`
  - `system._single_phase_equilibria: list[EquilibriumReaction]`
    (engine input)
  - `system._cross_phase_equilibria: list[EquilibriumReaction]`
    (read by `KineticGasLiquidLink.derive_speciation_keys` via
    public property `system.cross_phase_equilibria`)
  - `system._blackbox_models: list[BlackBoxReactionModel]`
- The speciation engine is built **lazily on first need**, cached on
  `system._engine`. Per-solve temperature handled via van 't Hoff
  inside the engine (no cache invalidation needed).
- Reaction list is **immutable post-attach**. No `add()` / `remove()`
  methods. Mutation requires constructing a new `ReactionSystem`.
- `ReactionSet.partition()` and `has_equilibrium` warning gone —
  partition is internal to `ReactionSystem`; mixed-set attach is
  always correct.

### `cv.advance()` flow

```
1. external feeds + boundaries        (apply_flux to any n_mol entry)
2. equilibrium refresh                (engine writes derived to n_mol)
3. PropertyCalculators (if any)       (viscosity, density → phase.properties)
4. kinetic + blackbox integration     (single Euler post-INTEGRATOR_REMOVAL)
5. internal transfer                  (link reads n_mol + properties directly)
```

- `t_h` ownership moves to the **orchestrator** (Simulation in
  `SIMULATION_CLASS`); passed per-call as a direct argument:
  `cv.advance(dt_h, t_h)`. Forwarded into `ReactionEnvironment` and
  boundary flux evaluation. CV does not store time.
- `chem_env` dict **disappears entirely**. Strong-ion totals are
  regular species in `n_mol`. `logH_guess` derived from
  `n_mol["H+"]` of the previous step. `T_K` lives on Phase.
- `ReactionEnvironment` reads pH from `phase.pH` and derived
  concentrations directly from `phase.n_mol`. No `PropertyResult`
  threading.
- `AdvanceResult.properties` field removed (no `PropertyResult`s to
  carry). `AdvanceResult.reaction_sources` and
  `AdvanceResult.transfer_record` survive.

### Property layer — narrow `PropertyCalculator`

- New `PropertyCalculator` protocol — narrow scope, **scalar
  derived properties only** (viscosity, density,
  heat capacity, compressibility, …).
  ```python
  @runtime_checkable
  class PropertyCalculator(Protocol):
      key: str
      def compute(self, phase: Phase, T_K: float, P_atm: float) -> float: ...
  ```
- `cv.property_calculators: List[PropertyCalculator]` — attached to
  CV, empty by default.
- Output stored on `phase.properties: Dict[str, float]` keyed by
  calculator `key`.
- Invoked **once, before kinetic reactions** in `advance()`, so
  kinetic rate laws can read `phase.properties["viscosity"]`. No
  post-step re-evaluation in v1 (add `update_after: bool` opt-in
  only if concrete need materialises).
- Speciation is NOT a `PropertyCalculator` — it's state completion,
  handled by `ReactionSystem`'s internal engine.

### Monitoring

- `AccuracyMonitor` (existing) injection path changes:
  ```python
  cv.reaction_system.attach_monitor(self._accuracy_monitor)
  ```
  Replaces the current attribute-injection through `property_solvers`
  list. `ReactionSystem` propagates to the internal engine on first
  build.
- **New `ConservationMonitor`** introduced. Element + charge
  accounting per `advance()` step.
  - Default policy: **`warn`**, mirroring `AccuracyMonitor`.
    Configurable via `VLsim.config.warnings.conservation`.
  - Thresholds: per-step `1e-8 × max(total_element_mol, 1)`;
    cumulative drift `1e-6 × max(total_element_mol, 1)`. Same shape
    for charge residual. Relative + max-clamp avoids spurious
    triggers in dilute systems.
  - Throttled via existing config throttle pattern (post
    `chemistry-unification-4`).
  - Catches: feeds without matching counterions, kinetic rates
    that violate stoichiometric balance, accumulated float roundoff.

### Solver strategy — operator splitting stays default

Operator splitting remains the default. The `StepSolver` protocol is
unchanged. DAE seam preserved for future opt-in. See
[`DAE_VS_OPERATOR_SPLITTING.md`](../upcoming/DAE_VS_OPERATOR_SPLITTING.md) for
the trigger conditions to revisit.

Crucially: the **unified `n_mol`** design makes the DAE seam clean
when it eventually arrives. Operator splitting and DAE share the
same data model; switching becomes a `StepSolver` choice, not a
data-model surgery.

## Checkpoint slicing

Six checkpoints, single feature branch (`state-unification` per
[[feedback_branching_convention]]), each ships a working state.

### C1 — Reaction subclass split

- Delete `Reaction` (kind-dispatched single class).
- Introduce `KineticReaction`, `EquilibriumReaction` as independent
  classes with their own `__init__` (no shared base).
- `BlackBoxReactionModel` unchanged (already independent).
- Shared validation in `src/reactions/_shared.py`.
- Builder migration: BSM2, ADM1, `models/vlmodels/`, demos, tests.
  Hard break, no shim.
- ~50–100 file touches, mechanical.

### C2 — `ReactionSet` → `ReactionSystem`

- Rename. Pre-bucketing by `isinstance` at construction.
- `partition()` and `has_equilibrium` become internal.
- `cv.reaction_system` (renamed from `cv.reaction_model`) as the
  attach point.
- `KineticGasLiquidLink.derive_speciation_keys` reads
  `system.cross_phase_equilibria` directly.
- Engine still constructed at attach time in this checkpoint
  (lazy build introduced in C3).

### C3 — Phase unification (B2 removal)

- `Phase.n_mol` widens to hold all species (totals + derived).
- `Phase.pH` becomes a `@property` computed from `n_mol["H+"]`.
- Engine's `solve()` writes derived species directly back to
  `phase.n_mol` via a privileged path
  (`phase._refresh_derived(values)`).
- Element-total conservation extracted automatically from the
  reaction graph by the engine — no special "TIC" / "CT_NH_T"
  declarations needed (but legacy declarations still work as
  conservation hints).
- `KineticGasLiquidLink._effective_henry` reads alphas inline
  from `n_mol` ratios — no separate `alphas` dict needed.
- BSM2 / ADM1 builders may need updates for declarations that
  conflated totals with molecular forms.

### C4 — Property layer collapse

- Delete `PropertySolver` protocol, `PropertyResult` dataclass,
  `SpeciationPropertySolver` adapter.
- Delete `cv.property_solvers` list. Replace with direct engine
  access via `cv.reaction_system._engine` (lazy-built).
- `chem_env` dict → `t_h: float` direct argument to
  `cv.advance(dt_h, t_h)`.
- `ReactionEnvironment` reads pH from `phase.pH`, derived
  concentrations from `phase.n_mol`. No `PropertyResult` threading.
- `AccuracyMonitor` injection moves to
  `cv.reaction_system.attach_monitor(monitor)`.
- `AdvanceResult.properties` field removed.
- ~150 tests directly affected (every test reading `result.properties`,
  every test constructing `chem_env`).

### C5 — `PropertyCalculator` protocol

- Introduce protocol + `cv.property_calculators` list.
- `phase.properties: Dict[str, float]` field.
- Invocation point: before kinetic reactions in `advance()`.
- **No concrete implementers in this checkpoint** — the protocol
  exists for future viscosity / density / heat-capacity models.
- Validates that the seam works without committing to a particular
  property model.

### C6 — `ConservationMonitor`

- Element + charge accounting on each `advance()` step.
- Configurable via `VLsim.config.warnings.conservation`.
- `cv.reaction_system.attach_conservation_monitor(monitor)` wiring
  pattern mirrors `AccuracyMonitor`.
- Default `warn` policy. Throttled.

Each checkpoint ships as a single commit (or small commit cluster)
with checklist update + per-checkpoint test re-baseline as needed.

## Resolved decisions

1. **Storage unified on Phase, B2 separation dropped.** Apply-flux
   writes to any species; equilibrium engine redistributes via
   element-total conservation. PHREEQC-style.
2. **Reaction taxonomy: three independent classes**
   (`KineticReaction`, `EquilibriumReaction`, `BlackBoxReactionModel`).
   No ABC, no shared base. Free-function utilities for shared
   validation.
3. **Hard break on `Reaction(kind=...)`** — no deprecation alias.
   Builders migrate in C1.
4. **`ReactionSystem` is the new container name** —
   `cv.reaction_system` attaches all reactions, internal pre-bucket
   by type.
5. **Time ownership moves to orchestrator.** `cv.advance(dt_h, t_h)`
   takes `t_h` as direct argument. `chem_env` dict disappears.
6. **`PropertyCalculator` runs before reactions** (so kinetic rates
   can read viscosity / density). Single invocation per step; no
   `update_after` in v1.
7. **`PropertySolver` / `PropertyResult` / `SpeciationPropertySolver`
   removed.** Speciation engine called directly from
   `cv.advance()` via `ReactionSystem`. Engine retains its
   richer surface for standalone use (`engine.solve(**kwargs)`).
8. **`AccuracyMonitor` injection via `attach_monitor` method**
   on `ReactionSystem`. Replaces attribute injection through
   `property_solvers` list.
9. **`ConservationMonitor` introduced, `warn` default.** Configurable
   thresholds via `VLsim.config.warnings.conservation`.
10. **Operator splitting stays default.** DAE deferred per
    [`DAE_VS_OPERATOR_SPLITTING.md`](../upcoming/DAE_VS_OPERATOR_SPLITTING.md).
    Unified `n_mol` keeps the DAE seam clean for future opt-in.

## What's load-bearing vs collapsing

### Load-bearing (kept)

- **Operator splitting** as default — stiffness argument from
  [`SPECIATION_REACTIONMODEL_BOUNDARY.md`](../design/SPECIATION_REACTIONMODEL_BOUNDARY.md)
  holds.
- **`SpeciationEngine` internals** — Newton + Davies activity model
  + ionic-strength fixed-point + Van 't Hoff temperature correction
  + `AccuracyMonitor` integration.
- **`StepSolver` protocol** — extension seam for future DAE add-on.
- **Reaction declaration shape** — `stoichiometry`, `Species`,
  `StoichiometryEntry`, atom/charge validation, `log_K`,
  `dH_J_per_mol`, `T_ref_K`, `is_cross_phase`.
- **Cross-phase equilibrium reactions** —
  `EquilibriumReaction(is_cross_phase=True)` for partition
  declarations consumed by `KineticGasLiquidLink`.
- **`AccuracyMonitor`** — same accuracy heuristics
  (pH change, Newton iters, scipy step rejections, ionic-strength
  regime), same throttle, same config integration.

### Collapsing (removed)

- **B2 separation** between conserved totals and derived species
  storage.
- **`PropertySolver`** protocol, **`PropertyResult`** dataclass,
  **`SpeciationPropertySolver`** adapter,
  **`cv.property_solvers`** list.
- **`ReactionSet.partition()`** public API + **`has_equilibrium`**
  attach-time warning machinery.
- **`chem_env`** dict (collapses to `t_h: float` argument).
- **`Reaction(kind="...")`** single class (replaced by independent
  subclasses).
- **`ReactionSet`** name (renamed to `ReactionSystem`).
- **`AdvanceResult.properties`** field.
- **`PropertyResult.alphas`** as a separate channel — alphas
  derived inline from `n_mol` ratios.
- **`logH_guess` threading through `chem_env`** — derived from
  `n_mol["H+"]` at previous step.

## Migration scope

Clean breaks, no shims:

- **Delete**: `src/core/property_solver.py`,
  `src/core/speciation_solver.py`, `chem_env` field across
  `ControlVolume.advance` and its callers, `Reaction` class,
  `ReactionSet` class.
- **Rewrite**: `src/core/control_volume.py` (advance flow, attach
  paths), `src/core/phases.py` (`n_mol` semantics, `pH` property,
  `_refresh_derived` privileged path, `properties` dict),
  `src/reactions/__init__.py` (exports), `src/reactions/builder.py`,
  `src/speciation/engine.py` (`from_reactions` consumes new
  `EquilibriumReaction`; `solve` writes back to `n_mol`),
  `src/core/gas_liquid_link.py` (alphas read inline).
- **Migrate (no shim)**: BSM2 / ADM1 builders in
  `models/vlmodels/`, all four fermenter demos, every test
  constructing `Reaction(kind=...)`, `chem_env=...`, or reading
  `result.properties`.

**Test impact:** ~150–200 of 821 tests directly affected, expected
post-phase total similar (deletions and additions roughly balance).

**Net code delta:** estimate −300 to −500 lines (collapse of
property machinery + chem_env dict + Reaction kind-branching
validation exceeds additions for `ReactionSystem` + `ConservationMonitor`
+ `PropertyCalculator`).

## Open items deferred

Not in scope; revisit when triggered:

- **`PropertyCalculator.update_after` opt-in** for post-step
  re-evaluation. Add when first concrete need surfaces (none today).
- **Precipitation reactions** (cross-phase liquid-solid with `Ksp`).
  Current `is_cross_phase` binary survives; richer typing
  (`cross_phase_kind: {"partition", "solubility"}`) deferred until
  first precipitation model lands.
- **Three-phase reactions** (gas-liquid-solid in one
  stoichiometry). Not currently supported, not a regression.
- **DAE solver add-on.** Trigger conditions in
  [`DAE_VS_OPERATOR_SPLITTING.md`](../upcoming/DAE_VS_OPERATOR_SPLITTING.md).
  `StepSolver` seam preserved.
- **`ReactionDeclaration` Protocol layer over the independent
  reaction classes.** Would add structural typing as a unified
  "any single-reaction declaration" interface on top of
  `KineticReaction` and `EquilibriumReaction`, enabling third-party
  reaction implementers (plugins, external model integrations) to
  satisfy the contract without subclassing or modifying VLsim's
  class set. Deferred because the current taxonomy is closed
  (kinetic + equilibrium covers the chemistry consensus) and no
  concrete third-party / plugin demand exists today. Revisit if:
  - A plugin / extension API for custom reaction kinds becomes a
    project goal.
  - External model integration patterns benefit from
    structural-conformance rather than VLsim-class subclassing.
  - Testing patterns benefit from mock-via-Protocol satisfaction.

  Discussed in 2026-05-20 design exploration; the Protocol approach
  (option c) was considered alongside ABC (option a) and the
  chosen independent-classes (option b). Adding the Protocol layer
  later is additive — the independent classes would continue to
  exist; the Protocol just names their shared shape and admits
  outside implementers. No re-architecture required.

## Relationship to other phases

- **Preceded by**
  [[integrator-removal-phase-designed]] —
  `INTEGRATOR_REMOVAL` ships first (clean, independent, ~150 LOC
  subtraction). Its single-Euler kinetic sub-step becomes the
  baseline for C4's kinetic integration. **Sequencing
  changed 2026-05-20** from the prior plan (where
  `INTEGRATOR_REMOVAL` was scheduled after `SIMULATION_CLASS`).
- **Followed by**
  [[simulation-class-phase-design-resolved]] —
  `SIMULATION_CLASS` builds on the new `cv.advance(dt_h, t_h)`
  signature, unified `n_mol`, and no-`chem_env` shape. **Two
  decisions in `SIMULATION_CLASS.md` require revision before
  implementation**:
  - **Decision 11** (`chem_env` permanence) — overridden;
    `chem_env` collapses to `t_h: float` argument.
  - **Decision 4** (records taxonomy) — likely fine; may add
    `conservation_records` channel if `ConservationMonitor`
    feeds into `BatchResult`. Probably defer to a recorder
    extension.
- **Builds on**
  [[project_cv_refactor]]'s chemistry-unification phases 1–4 —
  declaration unification (Phase 1), `alphas` channel (Phase 2),
  cross-phase equilibrium (Phase 3), AccuracyMonitor (Phase 4).
- **Aligned with**
  [[project_speciation_reactionmodel_boundary]] — operator
  splitting stays default; that note's stiffness defence holds.
  The "type split feasible but marginal" verdict in that note
  is upgraded to "type split happens as part of this refactor"
  (C1) because the file is already open for unrelated reasons.

## Cross-references

- [`SPECIATION_REACTIONMODEL_BOUNDARY.md`](../design/SPECIATION_REACTIONMODEL_BOUNDARY.md) —
  prior design note defending operator splitting on stiffness
  grounds.
- [`DAE_VS_OPERATOR_SPLITTING.md`](../upcoming/DAE_VS_OPERATOR_SPLITTING.md) —
  solver strategy note, DAE seam preservation, trigger conditions.
- [`INTEGRATOR_REMOVAL.md`](INTEGRATOR_REMOVAL.md) — preceding
  phase, single-Euler kinetic sub-step baseline.
- [`SIMULATION_CLASS.md`](../upcoming/SIMULATION_CLASS.md) — following phase
  (decisions 11 and 4 reconciled inline post-ship).
- [`CHEMISTRY_UNIFICATION.md`](../upcoming/CHEMISTRY_UNIFICATION.md) and
  [`CHEMISTRY_UNIFICATION_PLAN.md`](../upcoming/CHEMISTRY_UNIFICATION_PLAN.md) —
  declaration-unification track this refactor builds on.

## How to apply (start-of-work checklist)

When starting work on this phase:

- Follow [[feedback_branching_convention]] — create
  `state-unification` branch off `main`, write
  `STATE_UNIFICATION_CHECKLIST.md` modelled on
  `CHEMISTRY_UNIFICATION_4_CHECKLIST.md`, update the README's
  "Currently in flight" section.
- Per-checkpoint commit discipline applies (lesson from
  chemistry-unification-3's lost-edits incident).
- C1 starts with the inventory step: grep all
  `Reaction(kind=`, `chem_env`, `PropertyResult`,
  `PropertySolver`, `property_solvers`, `ReactionSet.partition`,
  `has_equilibrium` callsites and classify by checkpoint.
- BSM2 reference golden test
  ([`test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py))
  needs re-baselining after C3 (numerical: speciation now writes
  to `n_mol` directly) and re-checking after C4 (structural:
  `chem_env` gone). Document drifts inline per the
  chemistry-unification phases' precedent.
