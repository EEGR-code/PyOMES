# CUFermentationSpeciation Sunset — Design Note

> **Status: SHIPPED 2026-06-01.** Tag: `cufermenter-sunset-shipped`.
> Option C chosen (drop BioSTEAM bridge). All legacy island code deleted.
> Final test count: 913 (from 1021). See "Future BioSTEAM integration"
> section below for the recommended re-integration path.

---

> ~~**Status:** Trigger-gated. Surfaced 2026-05-27 at the end of the
> `simulation-class` phase as the *deferred legacy cleanup* that
> phase explicitly held out of scope. Pick up when: (a) the
> BioSTEAM integration is actively in use and needs a new bridge
> on top of `Simulation`, OR (b) the parallel-codebase maintenance
> cost of the legacy island crosses the threshold.~~

## Context

`CUFermentationSpeciation` is a 3,428-line `biosteam.Unit`-style
fermenter class in
[`models/vlmodels/fermenter/unit.py`](../../models/vlmodels/fermenter/unit.py).
It predates the CV refactor (Phases 1–7) and the
`simulation-class` phase. Its `simulate(feed, tau)` method drives
its own time-loop body (`_calc_ODE_with_headspace`) and consumes
the legacy `FermenterState` / `FermenterCommands` types directly.

A thin BioSTEAM-flowsheet adapter
[`models/vlmodels/fermenter/biosteam.py`](../../models/vlmodels/fermenter/biosteam.py)
inherits from BioSTEAM's `NRELBatchBioreactor` and delegates to
`CUFermentationSpeciation`. This is the **whole reason
`CUFermentationSpeciation` still exists** — it's the
BioSTEAM-compatible entry point for industrial process modelling
flowsheets. Without it, VLsim's BioSTEAM bridge stops working.

After the `simulation-class` phase shipped (2026-05-27),
`CUFermentationSpeciation` is the **only** remaining consumer of:

- `run_batch` in
  [`models/vlmodels/fermenter/config/factory.py`](../../models/vlmodels/fermenter/config/factory.py)
- Legacy controllers in
  [`src/control/loops.py`](../../src/control/loops.py)
  (~1,150 lines) and the legacy
  [`src/control/system.py:ControlSystem`](../../src/control/system.py)
- [`models/vlmodels/fermenter/types.py`](../../models/vlmodels/fermenter/types.py)
  (`FermenterState`, `FermenterCommands`)
- The legacy solver-strategy island in [`src/solvers/`](../../src/solvers/)
  (`CoupledSolver`, `SolverFactory`)
- The legacy types island in [`src/sim/`](../../src/sim/)
  (`Ledger`, `RunResult`, `ControlAction`, parallel `FermenterState`)
- The legacy controller subpackages in
  [`src/control/controllers/`](../../src/control/controllers/)
  (`ph.py`, `pressure.py` — pre-CV controllers consuming
  `FermenterState` / `FermenterCommands`),
  [`src/control/actuators/`](../../src/control/actuators/)
  (`relief_valve.py`), and
  [`src/control/builders/`](../../src/control/builders/)
  (`defaults.py` — imports legacy `loops.PHController`).
  CV-native equivalents shipped in `src/control/cv_loops.py`;
  these subpackages only remain consumed by
  `CUFermentationSpeciation`. The cross-import from
  [`src/control/cv_loops.py:1202`](../../src/control/cv_loops.py#L1202)
  reusing `mass_flow_kg_s` from `actuators/relief_valve.py`
  is a hidden retention hazard — the helper functions must
  move (or inline into `cv_loops.py`) before the actuators
  subpackage can be deleted.
- A stale reference doc at
  [`docs/reference/split_solver_reference.py`](../../docs/reference/split_solver_reference.py)
  references `FermenterState` and predates the CV refactor.
  Delete or annotate as historical at sunset time.
- *(Historical note — no longer applies.)* 3 of 4 demos in
  `demos/` (`cstr_fermenter.py`, `fed_batch_fermenter.py`,
  `microplate_fermenter.py`) were broken at module load
  immediately after `simulation-class` shipped. They were
  migrated onto the new `Simulation` orchestrator in commit
  `5fd8581` (post-ship, on `main`); all four demos now run.
  They no longer consume any legacy-island surface and so do
  not block this sunset. The four migrated demos subsequently
  moved to `demos/builder/` on 2026-05-29 — see
  [DEMO_RESTRUCTURE.md](DEMO_RESTRUCTURE.md).

These were intentionally retained at C14 of `simulation-class`
because deleting `CUFermentationSpeciation` requires the BioSTEAM
bridge question to be answered first. That question is the
centrepiece of this phase.

## The BioSTEAM bridge question

The reason `CUFermentationSpeciation` exists is to plug into
BioSTEAM's `bst.System` flowsheet model — converting BioSTEAM
streams to `FeedState`, delegating simulation to the standalone
core, and writing results back to outlet streams.

Three options for what replaces it:

### Option A — Rebuild the BioSTEAM bridge on top of `Simulation`

Write a new `BioSTEAMFermenter` (or similar) that:
- Inherits from `bst.Unit` (or `NRELBatchBioreactor` if the
  bioreactor-specific surface is wanted).
- Converts BioSTEAM streams to whatever `Simulation` needs (a
  configured `ControlVolume` plus controller / profile / solver
  / recorder wiring, typically via `FermenterBuilder.build_simulation`).
- Calls `sim.run(tau_h, n_steps)`.
- Writes the resulting `BatchResult` back to outlet streams.

This is the *architecturally correct* answer — BioSTEAM users
get the same orchestrator everyone else uses, and the legacy
island's tendrils all retract.

Cost: new BioSTEAM adapter class (~150–300 lines of bridge code),
test suite for the bridge, plus the deletion sweep below.

### Option B — Keep `CUFermentationSpeciation` indefinitely

Acknowledge that the legacy fermenter is *the* BioSTEAM bridge
and let it coexist with `Simulation` permanently. The framework
already supports two paths (legacy `Fermenter.simulate()` for
BioSTEAM users; new `Simulation.run()` for everyone else). The
maintenance cost is real but bounded, and the surface area only
grows if BioSTEAM-specific features are added.

Cost: ongoing maintenance of ~4,300 lines of legacy code; no new
work, but no cleanup either.

### Option C — Drop BioSTEAM integration entirely

`CUFermentationSpeciation` and its BioSTEAM adapter are deleted.
Anyone wanting VLsim-via-BioSTEAM rebuilds it externally.

Cost: feature regression for any external user of the BioSTEAM
bridge. Acceptable only if the actual use of the bridge is
near-zero.

## Migration scope (under Option A)

If Option A is picked, the sunset deletion sweep covers
approximately the following surface:

### Deletions

- [`models/vlmodels/fermenter/unit.py`](../../models/vlmodels/fermenter/unit.py)
  — the 3,428-line `CUFermentationSpeciation` class.
- [`models/vlmodels/fermenter/biosteam.py`](../../models/vlmodels/fermenter/biosteam.py)
  — the legacy adapter (~254 lines), replaced by the new
  `Simulation`-backed adapter.
- [`models/vlmodels/fermenter/types.py`](../../models/vlmodels/fermenter/types.py)
  — `FermenterState`, `FermenterCommands`.
- `run_batch` in
  [`models/vlmodels/fermenter/config/factory.py`](../../models/vlmodels/fermenter/config/factory.py)
  — its body and its `BatchResult` dataclass (the new framework's
  per-CV nested `BatchResult` in
  [`src/core/recorder.py`](../../src/core/recorder.py) replaces
  it).
- [`src/control/loops.py`](../../src/control/loops.py) — legacy
  FermenterState-based concrete controllers (~1,150 lines).
  `src/control/cv_loops.py` already contains the CV-native
  replacements.
- [`src/control/system.py`](../../src/control/system.py) — legacy
  `ControlSystem` (~155 lines). `Simulation._invoke_controllers`
  replaced its orchestration role at C9 of simulation-class.
- [`src/sim/`](../../src/sim/) — entire package (~190 lines:
  `Ledger`, `FermenterState`, `RunResult`, `ControlAction` —
  parallel/legacy types only consumed by
  `CUFermentationSpeciation`).
- [`src/solvers/`](../../src/solvers/) — entire package (~115
  lines: `CoupledSolver`, `SolverFactory` — the legacy strategy
  dispatch only consumed by `CUFermentationSpeciation`).
- The 3 originally-broken demos (`cstr_fermenter.py`,
  `fed_batch_fermenter.py`, `microplate_fermenter.py`) were
  rewritten onto the new `Simulation` orchestrator in commit
  `5fd8581` and now live at
  [`demos/builder/cstr_fermenter.py`](../../demos/builder/cstr_fermenter.py),
  [`demos/builder/fed_batch_fermenter.py`](../../demos/builder/fed_batch_fermenter.py),
  [`demos/builder/microplate_fermenter.py`](../../demos/builder/microplate_fermenter.py)
  post the 2026-05-29 demo restructure. No further action
  required at sunset time.

### Rewrites

- The new `BioSTEAMFermenter` adapter (Option A) — fresh module,
  likely in `models/vlmodels/fermenter/biosteam.py` or a new
  filename.
- `models/vlmodels/__init__.py` and
  `models/vlmodels/fermenter/__init__.py` — drop
  `CUFermentationSpeciation` exports; add the new adapter.
- Any callers of `Fermenter` / `CUFermentationSpeciation` in
  [`src/chemistry/registry.py`](../../src/chemistry/registry.py),
  [`src/chemistry/recipe.py`](../../src/chemistry/recipe.py),
  [`src/chemistry/compounds.py`](../../src/chemistry/compounds.py),
  [`src/speciation/strong_ions.py`](../../src/speciation/strong_ions.py),
  [`src/stream_adapter.py`](../../src/stream_adapter.py) (if any
  are non-trivial) — refactored or deleted.

### Test impact

Active test files exercising the legacy path:

- [`tests/standalone/test_integration.py`](../../tests/standalone/test_integration.py)
  — uses `F.simulate(feed, tau)` style. Significant rewrite or
  deletion needed.
- [`tests/standalone/test_feed_state.py`](../../tests/standalone/test_feed_state.py)
  — likely tests the `FeedState` adapter used by
  `CUFermentationSpeciation`. May survive if `FeedState` is
  retained as a stream-conversion type.
- [`tests/standalone/test_fermenter_construction.py`](../../tests/standalone/test_fermenter_construction.py)
  — tests `F.control_system` construction. Migrates to the new
  framework or deletes.
- Any other test that constructs a `Fermenter` via
  `FermenterFactory.create_volume(...)` + `run_batch(...)` —
  `test_factory.py` notably (dozens of tests). These would either
  migrate to `Simulation.run()` or delete (the new framework
  already has equivalent coverage in `test_simulation.py`).

Estimated ~50–100 tests touched; net test-count change depends
on how much migrate vs delete. Rough ballpark: −60 to −100 tests
deleted (legacy-only behaviour) plus +10 to +30 new (BioSTEAM
bridge coverage), final total similar or slightly down from the
post-`simulation-class` 977.

## Out of scope

- `HPLCColumn.simulate()` — independent concern; its sunset is
  the CONTAINER_LAYERING phase if it gets picked up.
- Renaming `FermenterBuilder.build_simulation()` → `build()`. The
  parallel coexistence pattern (Option B from the C8 / C11
  discussions) is the intended end state for the builder until
  this sunset clears the legacy `build()` path.

## Trigger conditions

Three conditions, any of which makes this work valuable:

1. **Active BioSTEAM use** — a real flowsheet project is using
   `CUFermentationSpeciation`/`BioSTEAMFermenter` and would
   benefit from the new `Simulation` machinery (controllers,
   profiles, lifecycle gating, recorders) being available via
   the BioSTEAM entry point.
2. **Maintenance cost spike** — a chemistry / speciation /
   numerics change that has to be made *twice* (once on
   `Simulation`'s path, once on `CUFermentationSpeciation`'s
   path). The legacy island has had several near-misses across
   `chemistry-unification-3`, `state-unification`, and
   `simulation-class`; each phase carved an explicit out-of-scope
   note for it. The cumulative drag is real but manageable so
   far.
3. **No BioSTEAM in this project's future** — if the BioSTEAM
   integration is unused and unlikely to be picked up, Option C
   (drop the bridge entirely) shrinks this to a pure deletion
   sweep with no new bridge to build.

## Relationship to other phases

- **Mostly independent of** the trigger-gated
  `chemistry-unification-3b`, `CONTAINER_LAYERING`, and
  `SPECIATION_LEVEL_RETIREMENT` phases. None of them depend on
  the legacy fermenter island.
- **Soft interaction with [RUN_HISTORY.md](RUN_HISTORY.md).**
  The sunset deletes `run_batch` in
  `models/vlmodels/fermenter/config/factory.py`. RUN_HISTORY's
  historical body still uses `run_batch` as its baseline
  (already flagged in that note's reframing banner as
  historical). If the sunset ships first, RUN_HISTORY's
  historical body references deleted code; the reframing
  banner already absorbs that fact but a final pruning would
  be in order at that point.
- **Follows** `simulation-class` (shipped 2026-05-27). All the
  machinery the new BioSTEAM bridge would use (`Simulation`,
  `FermenterBuilder.build_simulation`, CV-native controllers,
  RunContext gating) is in place.
- **Would close** the three-broken-demos thread that
  `simulation-class` C12 explicitly deferred — under Option A
  with rewrites, the three demos become canonical examples of
  the new framework.

## Open design questions (resolved 2026-06-01)

1. **Option A vs B vs C** — chose **Option C** (drop the bridge).
   BioSTEAM was not installed, not used in any active test or demo,
   and had no active flowsheet consumers.
2. N/A — Option A was not chosen.
3. N/A — Option A was not chosen.
4. **`FeedState`** — retained in
   [`src/stream_adapter.py`](../../src/stream_adapter.py) as a
   standalone feed-composition type (also used by
   `strong_ions.py` and suite-wide test fixtures).
   `FermenterResult` was deleted (it was BioSTEAM-result-specific).

---

## Future BioSTEAM integration

If BioSTEAM integration becomes valuable again, the recommended
re-integration path is a new `BioSTEAMFermenter` wrapping
`Simulation`. All the necessary infrastructure is already in place:

**Stream conversion** — `FeedState.from_biosteam_stream(stream)`
in [`src/stream_adapter.py`](../../src/stream_adapter.py) converts
a BioSTEAM `Stream` to a `FeedState` (duck-typed, no top-level
`import biosteam`). This is the intended feed-side entry point.

**Fermenter construction** — `FermenterBuilder.build_simulation()`
returns a `Simulation` pre-wired with CV, controllers, profiles,
solver, and recorder. The bridge adapter would call this after
converting the BioSTEAM inlet streams to `FeedState` / initial
liquid-phase concentrations.

**Result write-back** — `BatchResult` (from `VLsim.core.recorder`)
carries the final broth state; the adapter writes molar inventories
back to BioSTEAM outlet streams.

**Suggested implementation sketch:**

```python
class BioSTEAMFermenter(bst.Unit):
    def _run(self):
        feed = FeedState.from_biosteam_stream(self.ins[0])
        sim = (
            FermenterBuilder()
            .vessel(...)
            .transfer_kinetic(...)
            ...
            .build_simulation()
        )
        result = sim.run(tau_h=self.tau_h, n_steps=self.n_steps)
        _write_result_to_stream(result, self.outs[0])
```

**Design choices to revisit:**
- Inherit from `bst.Unit` or `NRELBatchBioreactor`? The latter
  adds titre/productivity parameters that may be useful.
- Expose `F.sim` attribute so the user can call `.run()` themselves?
  More flexible; less BioSTEAM-idiomatic.
- Add BioSTEAM as an optional dev dependency to enable bridge tests.

