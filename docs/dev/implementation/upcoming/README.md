# Phases — Upcoming

Planning notes for refactor work that is **not yet started**. Each
note describes the scope, design choices to resolve, and triggers
that would justify picking it up. Notes here are stable —
implementation only begins when one is moved through to active work
(i.e. a checklist file is added) and ultimately to
[../shipped/](../shipped/) once shipped.

## Currently in flight

*(none — `simulation-class` shipped 2026-05-27 as the third and
final of the three sequenced phases.)*

## Design discussions (pre-phase, not yet a checklist)

- **[RESERVOIR_TYPE.md](RESERVOIR_TYPE.md)** — 2026-07-10, revised
  2026-07-10. Now concludes with `FlowBoundary`, a real unifying protocol
  for `PhaseInterface`/`CVLink`/`ExternalBoundary`, enabled by a
  `PhaseCarrier` structural type (`ControlVolume` and the sibling
  `Reservoir` type both satisfy it). `Reservoir` (passive audit
  accumulator + `Phase`-borrowed multi-phase accounting, deliberately not
  `ControlVolume`-shaped) is the piece that gives `ExternalBoundary` a
  symmetric "other side," which is what makes the unification possible.
  Extends `MASS_EXCHANGE_ARCHITECTURE.md` §7/§12 Q1. `PartitionModel`'s
  relationship to `FlowBoundary` is explicitly deferred — shared open
  thread with `PHENOMENA_PROTOCOL.md`. Open questions listed in the note;
  no branch, no checklist, no code yet.
- **[PHENOMENA_PROTOCOL.md](PHENOMENA_PROTOCOL.md)** — 2026-07-10, Phase 1
  of two. A three-tier `Phenomena → {EquilibriumPhenomena, KineticPhenomena}`
  taxonomy grouping reactions and transfer models by evaluate-at-a-point
  vs. algebraic-constraint shape (per `MASS_EXCHANGE_ARCHITECTURE.md`
  §14.2's litmus test) rather than by chemistry/transport domain.
  `EquilibriumPhenomena` is largely already real (`EquilibriumConstraint` +
  `LAYER1_GAP_CLOSURE`'s Newton-fold); `KineticPhenomena` is supported by
  `SimultaneousAdaptiveSolver`'s existing combined-RHS evaluation.
  `ReactionSystem` renamed `PhenomenaSystem` in this note (role unchanged —
  bucket, own the engine, execute kinetic rates, provide the
  `ControlVolume` attachment point); `Phenomena` connects to `ControlVolume`
  through it, not directly. No behavior change scoped here — pure taxonomy.
  Followed by
  [KINETIC_TRANSFER_GENERALIZATION.md](KINETIC_TRANSFER_GENERALIZATION.md).
- **[KINETIC_TRANSFER_GENERALIZATION.md](KINETIC_TRANSFER_GENERALIZATION.md)** —
  2026-07-10, Phase 2 of two (depends on `PHENOMENA_PROTOCOL.md`). Gives
  `KineticTransferModel` an optional custom rate function (currently
  constrained to linear `kLa` relaxation) so it can actually satisfy
  `KineticPhenomena` instead of staying inert config deferred to
  `KineticGasLiquidLink`. `PhaseInterface`, which the generalized class
  would satisfy, is now `FlowBoundary`'s same-CV case per
  `RESERVOIR_TYPE.md` §5.1 — framing note only, no design change. Higher-risk half — touches the shipped
  `transfer_models=` API. Five open questions listed in the note
  (signature, `partition_model` optionality, `KineticGasLiquidLink`'s
  role, `EquilibriumTransferModel` parity, regression scope); no branch,
  no checklist, no code yet.

## Recently shipped

- `step-solver-interface-refinement` (2026-07-08) — Stage 1 of the
  solver-interface-refinement plan (see below). Ownership-guard
  `OrchestrationWarning` on `cv.advance()`; `MonolithicODESolver` rejects a
  per-CV `solver=` it can never consult; `SequentialAdvanceSolver` reified
  (`solver=None` is now sugar for a real class); both `Simultaneous*` solvers
  generalized off the gas/liquid assumption (only `"liquid"` required now;
  every registered `internal_interfaces` entry contributes, not just the
  first); shared swappable `clamp_fn` module
  (`proportional_clamp`/`floor_clamp`, `clamp_fn=None` disables clamping);
  `negative_mole`/`clamp_invoked` `AccuracyMonitor` diagnostics; unified
  `StateVector` class replacing three near-duplicate packing
  implementations; a custom-`StepSolver`/`SystemSolver` demo. `_StateVector`
  deleted (no shim). 1950 → 2003 tests. Tag
  `step-solver-interface-refinement-shipped`. See
  [../shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md](../shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md)
  and
  [STEP_SOLVER_REFINEMENT_CHECKLIST.md](STEP_SOLVER_REFINEMENT_CHECKLIST.md).

- `simulation-class` (2026-05-27) — third and final of the three
  sequenced phases. New `Simulation` class as the single
  orchestration pathway for all CV-based models. Subsumed
  `MultiCVSystem` (deleted); absorbed the `run_batch` time-loop
  body; promoted controllers and profiles to first-class concepts
  with structured per-step records; ported the legacy `src/control/`
  controller machinery onto CV-shaped state via the Pattern 1
  collapse (`compute(state, dt_h) -> ControlAction` — no
  intermediate `Commands` or separate `Actuator`); added
  RunContext-based lifecycle gating with Pattern B unchecked
  setters for the orchestrator-mediated mutation path. Final
  suite 977/0; tag `simulation-class-shipped`. See
  [../shipped/SIMULATION_CLASS.md](../shipped/SIMULATION_CLASS.md)
  and
  [../shipped/SIMULATION_CLASS_CHECKLIST.md](../shipped/SIMULATION_CLASS_CHECKLIST.md).

- `state-unification` (2026-05-22) — second of the three
  sequenced phases. Collapsed speciation-shaped type machinery
  (`PropertySolver` / `PropertyResult` / `SpeciationPropertySolver`
  / `chem_env` / `ReactionSet.partition` / `Reaction`
  kind-branching) and unified `Phase` storage on PHREEQC-style
  `n_mol` (totals + derived). `cv.reaction_system` is the single
  reaction attach point; the new narrow `PropertyCalculator`
  protocol replaces `PropertySolver`; the engine writes derived
  species back to `n_mol` via `phase._refresh_derived(values)`. C6
  adds `ConservationMonitor` (element + charge accounting per
  step). Final suite 793 / 0; tag `state-unification-shipped`.
  BSM2 sentinels re-baselined across the phase (pH 3.340171 →
  3.303580). See
  [../shipped/STATE_UNIFICATION.md](../shipped/STATE_UNIFICATION.md)
  and
  [../shipped/STATE_UNIFICATION_CHECKLIST.md](../shipped/STATE_UNIFICATION_CHECKLIST.md).

- `integrator-removal` (2026-05-20) — deleted the `Integrator`
  Protocol + `RK4Integrator` + `EulerIntegrator` (~150 LOC
  subtraction); the sequential `cv.advance()` body's reaction
  sub-step collapses to a single forward Euler evaluation;
  `StepSolver` becomes the sole user-facing extension surface for
  time integration. BSM2 golden re-baselined (largest drift S_h2
  at ~6e-3 relative; pH 1.6e-7 relative). 816 standalone tests
  passing (down from 821 by 5 deleted tests on the removed
  surface). First of three sequenced phases (resequenced
  2026-05-20): `INTEGRATOR_REMOVAL` → `STATE_UNIFICATION` →
  `SIMULATION_CLASS`. See
  [../shipped/INTEGRATOR_REMOVAL.md](../shipped/INTEGRATOR_REMOVAL.md)
  and
  [../shipped/INTEGRATOR_REMOVAL_CHECKLIST.md](../shipped/INTEGRATOR_REMOVAL_CHECKLIST.md).

**Chemistry unification is fully complete.** Six precursor / phase
branches shipped, followed by the trigger-gated `chemistry-unification-3b`
(shipped 2026-06-03). All seven branches:

- `reaction-protocol-cleanup` (2026-05-12) — see
  [../shipped/REACTION_PROTOCOL_CLEANUP_CHECKLIST.md](../shipped/REACTION_PROTOCOL_CLEANUP_CHECKLIST.md).
- `bsm2-reference-test` (2026-05-12) — golden-trajectory regression
  for the BSM2 reference model so chemistry-unification numerics
  drift is visible; see
  [../shipped/BSM2_REFERENCE_TEST.md](../shipped/BSM2_REFERENCE_TEST.md).
- `chemistry-unification-1` (2026-05-13) — declared equilibrium
  reactions consumed by the speciation engine, `chem_env` shrunk to
  its minimal contract, stateless snapshot model. See
  [../shipped/CHEMISTRY_UNIFICATION_1_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_1_CHECKLIST.md).
- `chemistry-unification-2` (2026-05-15) — `PropertyResult.alphas`
  channel; gas-liquid link collapses inline `f_molecular` to a single
  alpha read; `SpeciationCorrection` cascade deleted; dormant
  `_check_f_molecular_consistency` removed. Leak 1 (and Leak 3)
  closed. See
  [../shipped/CHEMISTRY_UNIFICATION_2_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_2_CHECKLIST.md).
- `chemistry-unification-3` (2026-05-15) — canonical-name dedup in
  `_compute_species_eq`; suffix-based `ionic_strength_from_speciation`
  (unmasks a latent NH4+/I bug introduced by Phase 1, BSM2 sentinels
  re-baselined); cross-phase equilibrium reaction infrastructure
  (`Reaction.is_cross_phase`, `log_K=None` for cross-phase, engine
  filters them silently, `KineticGasLiquidLink.derive_speciation_keys`
  wired into `ControlVolume.__init__`). Leak 2 infrastructure shipped;
  BSM2/ADM1 builder migration deferred to Phase 3b alongside the
  unified Species-ID emission convention. See
  [../shipped/CHEMISTRY_UNIFICATION_3_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_3_CHECKLIST.md).
- `chemistry-unification-4` (2026-05-16) — `AccuracyMonitor`
  attached to every `ControlVolume`; `AccuracyWarning(UserWarning)`
  category; `WarningConfig` thresholds + throttle on
  `VLsim.config.warnings`; `VLSIM_WARNINGS` env-var with flat
  presets. Six cheap per-step checks (pH change, Newton iters,
  charge residual, dt vs τ_min stub, scipy step rejections,
  ionic-strength regime). The legacy `_warned_high_I` set on
  `ChemicalEquilibriumEngine` is gone — its ionic-strength threshold
  warning migrates to `AccuracyWarning` via the property solver.
  See
  [../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md).

- `chemistry-unification-3b` (2026-06-03) — `ChemistryDatabase` rollout,
  unified Species-ID emission, BSM2/ADM1 cross-phase migration. Tag
  `chemistry-unification-3b-shipped`. See
  [../shipped/CHEMISTRY_UNIFICATION_3B_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_3B_CHECKLIST.md).

- `partition-model` (2026-06-16) — `PartitionModel` protocol +
  `HenryPartition`; `ChemistryDatabase.partition_models`; H₂S alpha
  correction (2× error at pH=pKa fixed); `_HENRY_PARAMS` + `henry_mol_L_atm()`
  deleted; `vfa_volatility` flag removed from ADM1 builders. 1390→1428 tests.
  Tag `partition-model-shipped`. See
  [../shipped/PARTITION_MODEL.md](../shipped/PARTITION_MODEL.md).

## Priority order

> All phases listed below have now shipped. The priority list is
> preserved as a history of scope decisions and trigger rationale.
> See [OPEN_WORK.md](../OPEN_WORK.md) for currently open items.

1. **[../shipped/CHEMISTRY_UNIFICATION.md](../shipped/CHEMISTRY_UNIFICATION.md)** (design)
   + **[../shipped/CHEMISTRY_UNIFICATION_PLAN.md](../shipped/CHEMISTRY_UNIFICATION_PLAN.md)**
   (plan) — unify the `Reaction` and speciation declaration surfaces
   so equilibrium reactions become a mode of the existing reaction
   framework. *All seven branches shipped (`chemistry-unification-1`
   through `chemistry-unification-4` plus `chemistry-unification-3b`).*
   *Shipped fully 2026-06-03.*

2. **[../shipped/RUN_HISTORY.md](../shipped/RUN_HISTORY.md)** —
   `StreamingFileRecorder` (Parquet chunks + `load_run`), `SparseRecorder`,
   and `SummaryRecorder` + `SummaryResult`. Built on the `Recorder` /
   `BatchRecorder` / `BatchResult` schema shipped in `simulation-class`.
   *Shipped 2026-06-08; tag `run-history-shipped`.  1174→1226 tests.*

3. **[../shipped/SPECIATION_LEVEL_RETIREMENT.md](../shipped/SPECIATION_LEVEL_RETIREMENT.md)** —
   retired the `level` (1, 2, 2.5, 3) parameter on `ChemicalEquilibriumEngine`
   and the three parallel solver implementations. Deleted 6 level/adapter
   files (2,448 lines); unified engine routes directly to
   `solve_from_equilibrium_set` or `solve_acid_base`.
   *Shipped 2026-06-03; tag `speciation-level-retirement-shipped`.*

4. **[../shipped/CUFERMENTER_SUNSET.md](../shipped/CUFERMENTER_SUNSET.md)** —
   retired `CUFermentationSpeciation` and dropped the BioSTEAM bridge
   (Option C). Deleted 9,093 lines: `unit.py`, `loops.py`, `sim/`,
   `solvers/`, `actuators/`, `FermenterResult`, and 7 legacy test files.
   `FeedState` survives.
   *Shipped 2026-06-01; tag `cufermenter-sunset-shipped`.*

## Recently surfaced (2026-05-27 audit)

Trigger-gated planning notes added from the
post-`simulation-class` audit. Not yet ranked into the
numbered priority list above — see each note's
*Trigger conditions* section for the criteria that would pull
it in. All four notes share the audit's framing
("fermenter-shaped framework surfaces dressed up as general")
and naturally bundle if a single "generalisation pass" is
taken.

- **[../shipped/FRAMEWORK_POLISH.md](../shipped/FRAMEWORK_POLISH.md)** —
  P1–P9 all done (~200 LOC). Fixed silent bug: `GasFeed.y` param path
  was discarding writes.
  *Shipped 2026-06-01; tag `framework-polish-shipped`.*
- **[../shipped/PARAM_PATH_DISPATCHER.md](../shipped/PARAM_PATH_DISPATCHER.md)** —
  `MutableScalar` / `MutableDict` descriptors; typed `ParamPath` builder;
  reflective walker. All `@property + _set_X_unchecked` triplets removed;
  unknown paths raise `ParamPathError`. Subsumed FRAMEWORK_POLISH P8.
  *Shipped 2026-06-08; tag `param-path-dispatcher-shipped`.*

## Recently surfaced (2026-05-28 exploration session)

Trigger-gated planning notes added from the 2026-05-28
post-simulation-class exploration session (demos / exports /
checkpointing / logging / FBA / chemistry-database
questions). Cross-cutting interaction: the three "recorder
landscape" notes (RESULT_EXPORT, HPC_CHECKPOINTING,
RUN_HISTORY's sub-deliverables) all touch the same surface;
shipping them together gives a coherent crash-resilient
recorder story.

- **[../shipped/DEMO_RESTRUCTURE.md](../shipped/DEMO_RESTRUCTURE.md)** —
  *Shipped 2026-05-29 with deviations from option β; see the
  reconciliation banner at the top of the planning note for the
  actual layout.* Four fermenter demos moved to `demos/builder/`;
  new `demos/model_api/` umbrella holds `chemistry/` and
  `model_construction/`. Q1 / Q2 demos shipped alongside; Q7
  chemistry-database closed by `chemistry-unification-3b`.
- **[../shipped/RESULT_EXPORT.md](../shipped/RESULT_EXPORT.md)** —
  `BatchResult.to_dataframe / to_csv / to_parquet / to_wide`; 5 NumPy
  channels; pandas + pyarrow optional; demo at
  `demos/model_api/export_results.py`.
  *Shipped 2026-06-01; tag `result-export-shipped`.  32 tests.*
- **[../shipped/HPC_CHECKPOINTING.md](../shipped/HPC_CHECKPOINTING.md)** —
  `sim.save_checkpoint(path, mode="data+env")` API with four
  fidelity levels (`data` / `data+env` / `source` / `full`); `data`/`data+env`
  implemented, `source`/`full` raise `NotImplementedError`.
  *Shipped 2026-06-08; tag `hpc-checkpointing-shipped`. 78 tests.* This
  entry previously read "not yet implemented" — stale; corrected
  2026-07-09 during a scoping-review hygiene pass after cross-checking
  against `git log` found the tag already merged to `main`.

## Recently surfaced (2026-06-11 solver architecture)

Six design docs added from the 2026-06-11 solver architecture session.
Phases A–E all shipped 2026-06-11/12; docs moved to `shipped/`.
See **[SOLVER_ARCHITECTURE.md](SOLVER_ARCHITECTURE.md)** for the anchor doc
and Phase F+G placeholders (SUNDIALS opt-in, no work scheduled).

- **[../shipped/CV_COMPUTE_INTERFACE.md](../shipped/CV_COMPUTE_INTERFACE.md)** — Phase A. *Shipped 2026-06-11.*
- **[../shipped/CONTROLLER_STATE_PROTOCOL.md](../shipped/CONTROLLER_STATE_PROTOCOL.md)** — Phase B. *Shipped 2026-06-12.*
- **[../shipped/SYSTEM_SOLVER_PROTOCOL.md](../shipped/SYSTEM_SOLVER_PROTOCOL.md)** — Phase C. *Shipped 2026-06-12.*
- **[../shipped/IMPLICIT_TRANSPORT.md](../shipped/IMPLICIT_TRANSPORT.md)** — Phase D. *Shipped 2026-06-12.*
- **[../shipped/MONOLITHIC_ODE.md](../shipped/MONOLITHIC_ODE.md)** — Phase E. *Shipped 2026-06-12.*

## Solver interface refinement (multi-stage, Stage 1 shipped)

Multi-round design review of `StepSolver`/`SystemSolver` customizability,
triggered by auditing how flexible `cv.advance()`/`Simulation` solver
selection actually is. Extends the shipped `SOLVER_ARCHITECTURE.md` two-axis
design (Phases A–E) — not the SUNDIALS/DAE Phase F/G track.

- **[../shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md](../shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md)** —
  the design note. 9 findings, 8 shipped as Stage 1 (below); item 5
  (composition/integrator/clamp_fn factoring) explicitly deferred; the
  interleaving/multi-CV-aware solver tier discussion is Stages 2–4 of the
  plan, not yet resolved.
- **[STEP_SOLVER_REFINEMENT_PLAN.md](STEP_SOLVER_REFINEMENT_PLAN.md)** —
  the 4-stage sequencing plan. Stage 0 (rename) and Stage 1 (the bundle)
  shipped; Stage 2 (SIA z-staleness design decision, §12 Q7) and Stage 3
  (reactive D_eff transport) still pending — see the plan doc for current
  status.
- **[STEP_SOLVER_REFINEMENT_CHECKLIST.md](STEP_SOLVER_REFINEMENT_CHECKLIST.md)** —
  Stage 1's full checkpoint-by-checkpoint implementation log (10
  checkpoints, all shipped), including several scope corrections found
  during implementation — kept here rather than archived since it documents
  decisions (e.g. the checkpoint 4/7b split) relevant to Stages 2–4 still
  to come.

## Upcoming phases

### NR Precipitation (two-phase sequence)

Design discussion 2026-06-23. Active-set precipitation equilibrium built on top
of the NR speciation engine shipped in `nr-speciation-engine`.

- **[NR_PRECIPITATION_SPECIATION.md](NR_PRECIPITATION_SPECIATION.md)** —
  Phase 1: speciation layer only. Outer active-set loop in
  `NRChemicalEquilibriumEngine.solve()`; `precipitation_equilibria` bucket on
  `ReactionSystem`; `element_stoichiometry` cross-component mass balance fix;
  `"minerals"` key in output dict; `Ca_plus_plus` / `Mg_plus_plus` in
  `common_species`; `05_precipitation_equilibrium.ipynb` demo. No CV changes.

- **[NR_PRECIPITATION_CV_INTEGRATION.md](NR_PRECIPITATION_CV_INTEGRATION.md)** —
  Phase 2: CV/SolidPhase integration. `_read_from_phases` sums solid
  contribution; `SolidPhase` writeback; consistent phase-type validation across
  all three phase types on `ControlVolume`. Depends on Phase 1 shipping first.

- **[MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md](MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md)** —
  a separate, standalone-by-design track (Fe/Ca/phosphate/citrate networks
  `NRTableau` can't represent). **Partial skeleton, stalled 2026-06-29** —
  see the status banner at the top of that doc; committed directly to `main`
  without a branch or checklist, only one of seven modules tested. Not part
  of the two-phase NR Precipitation sequence above.

### Layer 1 gap closure (two-phase sequence)

Design discussion 2026-07-01 (see `MASS_EXCHANGE_ARCHITECTURE.md` §14 and
`CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md` §3/§18). Folds gas-liquid VLE
(and, per §14, recognizes Ksp as already folded) into the same simultaneous
Newton system as acid-base, closing the "⚠ gap" row in
`MASS_EXCHANGE_ARCHITECTURE.md` §10.4. Reconciles three previously-conflicting
precipitation design threads (`NR_PRECIPITATION_SPECIATION`,
`NR_PRECIPITATION_CV_INTEGRATION`, `MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN`)
without changing any of their shipped/planned numerics.

- **[EQUILIBRIUM_CONSTRAINT_UNIFICATION.md](../shipped/EQUILIBRIUM_CONSTRAINT_UNIFICATION.md)** —
  Phase 1: **shipped 2026-07-02** (tag `equilibrium-constraint-unification-shipped`).
  Declaration-API generalization: `EquilibriumConstraint` protocol +
  `HenryEquilibrium`/`KspEquilibrium`/`RaoultEquilibrium` sibling types (not
  alternate `EquilibriumReaction` constructors — a mid-flight refinement);
  single `classify_equilibrium_constraint()` auto-classification path
  replacing the `is_cross_phase` silent-filter and `precipitation_reactions=`
  kwarg patchwork; fixes the `HenryPartition`/reaction double-declaration bug
  (`ReactionSystem` bucketing reworked to reach it); per-entry
  activity-dispatch helper (groundwork for Phase 2). No solve-time behavior
  change; no rename.

- **[LAYER1_GAP_CLOSURE.md](../shipped/LAYER1_GAP_CLOSURE.md)** —
  Phase 2: **shipped 2026-07-03** (tag `layer1-gap-closure-shipped`). Actual
  Newton-system folding: gas-liquid (Henry) and Raoult/evaporation rows enter
  `g(y,z)=0`; `step_internal_transfer()` scope-filter fix (found no filter
  existed at all beforehand, not a partial one); `WaterVapourBoundary`/
  `VentWaterLoss` retired; precipitation regression verification (numerics
  unchanged, but surfaced a real "CO₂ stripping promotes CaCO₃ scaling"
  effect); closes with the `SpeciationEngine` → `ChemicalEquilibriumEngine`
  rename, plus the `src/speciation/` → `src/chemical_equilibrium/` module
  rename (§18's own deferred decision, executed at this phase's close).

---

## How to start one

1. Read the relevant note in this folder. It captures the design
   thinking and any decisions already pinned.
2. Write a checklist file (e.g. `PHASE6_CHECKLIST.md`) modelled on
   [../shipped/PHASE5_CHECKLIST.md](../shipped/PHASE5_CHECKLIST.md).
   Drop it into this folder while it is being worked.
3. Create a feature branch off `main` named for the phase
   (e.g. `chemistry-unification`).  All implementation lands on that
   branch.
4. Update the "Currently in flight" line above with a pointer to the
   checklist.
5. Push the branch periodically so the work is backed up to GitHub
   while in progress.
6. When the work ships, move both the original note and the
   checklist to [../shipped/](../shipped/), add a
   "Status: Shipped" banner to each, and update the priority list
   above to remove the entry.  Then follow the branching and tagging
   convention below to ship the branch back to `main`.

## Branching and tagging convention

This repo uses a feature-branch-per-phase workflow with explicit
shipping tags.  The convention is:

1. **One branch per phase.**  Create the branch off `main` at the
   start of a phase (e.g. `git checkout -b chemistry-unification`).
   All implementation, doc updates, and intermediate commits land on
   that branch.  Push it to `origin` so the work is backed up.
2. **Merge with `--no-ff` when shipping.**  Once the checklist is
   complete and tests are green, merge the branch into `main` with
   an explicit merge commit:
   ```
   git checkout main
   git merge --no-ff <phase-branch> -m "Merge <phase-branch>: <one-line summary>"
   ```
   The `--no-ff` flag keeps the phase boundary visible in `git log`
   as a fork-and-rejoin shape, even when a fast-forward would
   otherwise be possible.
3. **Tag the shipping commit.**  Create an immovable named pointer
   on the phase's final commit so it can be checked out later by
   name:
   ```
   git tag <phase-name>-shipped <commit-hash>
   ```
   Use names like `phase-7-shipped`, `chemistry-unification-shipped`,
   etc.  These appear on the GitHub "Tags" page and are addressable
   via URLs like `https://github.com/.../tree/<tag-name>`.
4. **Push the branch *and* the tags.**  `git push` ships the merge;
   `git push --tags` ships the tag.  Both are needed — tags are not
   pushed by default.
5. **Delete the feature branch after merge.**  The work now lives on
   `main` and is permanently addressable via the tag, so the branch
   pointer is redundant:
   ```
   git branch -d <phase-branch>
   git push origin --delete <phase-branch>
   ```

This convention preserves every commit as a permanent checkpoint
while keeping `main` the current truth and giving each phase a
human-readable name for browse-back-in-time.
