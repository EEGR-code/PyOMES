# RUN_HISTORY — Shipped 2026-06-08

> **Shipped:** tag `run-history-shipped`. `StreamingFileRecorder` +
> `load_run` (C1); `SparseRecorder` + `SummaryRecorder` + `SummaryResult`
> (C2). 1174 → 1226 tests.

# Unified Run History — Design Note (historical)

> **Reframed 2026-05-27 (post-simulation-class):** The
> `simulation-class` phase shipped the `Recorder` Protocol +
> default `BatchRecorder` (C3) which already produces a unified
> time-indexed result (`BatchResult`) and is pluggable. This
> note is **no longer about rewriting the result schema** — it
> is about whether a **richer default recorder** is worth
> shipping. The body text below predates that reframing and
> uses the legacy `run_batch` orchestrator as its baseline.
> Treat it as a historical sketch of the unified-history idea;
> the actual implementation in `src/core/recorder.py` is the
> current authoritative shape.
>
> **Concrete sub-deliverables (added 2026-05-28).** The
> "richer default recorder" reframing is made concrete below
> under [§Concrete sub-deliverables](#concrete-sub-deliverables-added-2026-05-28).
> Three new pluggable recorder implementations
> (`StreamingFileRecorder`, `SparseRecorder`, `SummaryRecorder`)
> plus the Parquet-rotated-chunks format decision.
>
> **Cross-link surfaces (added 2026-05-28).**
> - **[CUFERMENTER_SUNSET.md](../shipped/CUFERMENTER_SUNSET.md)**
>   deleted `run_batch` outright on 2026-06-01
>   (tag `cufermenter-sunset-shipped`). The historical body
>   below now references deleted code. The promised body prune
>   is overdue — treat everything below §Concrete
>   sub-deliverables as historical context only.
> - **[PARAM_PATH_DISPATCHER.md](PARAM_PATH_DISPATCHER.md)**
>   defines the typed/string `ParamPath` shape that a richer
>   recorder would naturally use for the `params_changed`
>   channel. If RUN_HISTORY ships before that dispatcher
>   rewrite, the recorder bakes in the current closed-vocabulary
>   strings.
> - **[FRAMEWORK_POLISH.md P9](FRAMEWORK_POLISH.md#p9-zoh-hold-placeholder-disambiguation--severity-low)**
>   answers the small data-model question of "never-fired vs
>   fired-with-no-effect" `ControlAction` records. RUN_HISTORY's
>   redesigned recorder must take the same answer; if P9 ships
>   first, RUN_HISTORY consumes it.
> - **[HPC_CHECKPOINTING.md](HPC_CHECKPOINTING.md)**
>   `Simulation.save_checkpoint` (shipped in `hpc-checkpointing`)
>   always pickles the recorder because `StreamingFileRecorder`
>   does not yet exist. **Once `StreamingFileRecorder` ships here,
>   add an `isinstance` check in `save_checkpoint`:** if the
>   attached recorder is a `StreamingFileRecorder`, skip pickling
>   its mid-array buffer — the flushed Parquet chunks on disk
>   already cover the time-series state. The `include_recorder=True`
>   default stays correct for `BatchRecorder`.

## What this is

A planning note capturing the case for a **unified history object**
that aggregates the per-step `AdvanceResult` outputs of a
simulation into a single time-indexed result the orchestrator
returns to the caller.

The idea originally surfaced during Phase 7 (GLV removal) when
the question of where `pH` / `ionic_strength` should live came
up. Phase 7 chose the minimal answer (drop the convenience
accessors; callers read `result.properties["speciation"].pH` off
each step's `AdvanceResult`); the unified-history idea was kept
as the architecturally cleaner option for a future phase to pick
up. The `simulation-class` phase (shipped 2026-05-27) then
delivered the **structural** half — a `Recorder` Protocol +
`BatchRecorder` producing a per-CV-nested `BatchResult` — which
covers most of what this note originally proposed.

This note now captures the residual question: **what would a
richer default recorder look like?**

## Why it surfaced

After Phase 6 shipped, `GasLiquidVolume` retained two convenience
accessors — `glv.pH` and `glv.ionic_strength` — backed by an
internal `_last_result` cache populated each `advance()` call.
These accessors let callers read post-step properties between
timesteps without holding a reference to the `AdvanceResult`
themselves.

When Phase 7 set out to delete `GasLiquidVolume`, three options
were on the table:

- **A.** Re-home the accessors on `ControlVolume` (`cv.pH`,
  `cv.ionic_strength`) backed by a CV-level `_last_result` cache.
- **B.** Move to a `RunResult` history object the orchestrator
  carries.
- **C.** Drop the convenience accessors entirely; callers read
  `result.properties["speciation"].pH` off the `AdvanceResult`.

C was picked because the live caller surface for `glv.pH` is small
(3 non-test sites + 2 tests) and Chemistry Unification is expected
to install the architecturally proper `phase.pH` (with staleness
guardrails) later anyway — option A would just be temporary scaffolding.

This file captures **why option B is still architecturally
attractive** so a future phase can pick it up cleanly.

## The idea

Today every solver path returns a single `AdvanceResult` per
timestep.  The orchestrator (typically `run_batch` in
[../../models/vlmodels/fermenter/config/factory.py](../../models/vlmodels/fermenter/config/factory.py))
collects those into ad-hoc parallel arrays:

```python
# Today — orchestrator builds parallel arrays manually
t_h: List[float] = []
gas_mol: Dict[str, List[float]] = {}
liquid_mol: Dict[str, List[float]] = {}
pH: List[float] = []
ionic_strength: List[float] = []
P_atm: List[float] = []
transfer_flow: List[Dict[str, float]] = []
boundary_records: List[List[ExternalFluxRecord]] = []
advance_results: List[AdvanceResult] = []

for i in range(n_steps):
    result = cv.advance(dt_h, chem_env=ctx)
    t_h.append(...)
    gas_mol["O2"].append(...)
    pH.append(result.properties["speciation"].pH)
    # … etc.
```

The proposal is a single result object that absorbs this
bookkeeping:

```python
@dataclass
class RunResult:
    """Time-indexed history of a simulation run."""
    t_h: np.ndarray
    advance_results: List[AdvanceResult]
    # Convenience accessors derived from advance_results:
    def pH(self) -> np.ndarray: ...
    def ionic_strength(self) -> np.ndarray: ...
    def gas_mol(self, species: str) -> np.ndarray: ...
    def liquid_mol(self, species: str) -> np.ndarray: ...
    def P_atm(self) -> np.ndarray: ...
    def transfer_flow(self, species: str) -> np.ndarray: ...
    def boundary_record(self, label: str) -> List[ExternalFluxRecord]: ...
```

`run_batch` then becomes:

```python
def run_batch(cv, tau_h, n_steps, ...) -> RunResult:
    advance_results = []
    for i in range(n_steps):
        result = cv.advance(dt_h, chem_env=ctx)
        advance_results.append(result)
    return RunResult(t_h=np.arange(...), advance_results=advance_results)
```

Callers reading `result.pH()` get a vector instead of a scalar; the
underlying data already lives on each `AdvanceResult`, so the
accessors are pure derivations.

## Trade-offs vs current state

**Pros:**

- One clear return type for "what happened during the run",
  replacing today's ~8 parallel arrays + 3 list fields on
  `BatchSimulationResult`.
- Per-step state and cumulative history are cleanly separated:
  `AdvanceResult` is per-step; `RunResult` is the aggregate.
- New time-series quantities (e.g. controller setpoint history,
  derived metabolic rates) hook in by adding accessors to
  `RunResult`, not by adding fields to a result dataclass and
  threading them through the loop.
- The pH-controller pattern (today: `pH = glv.pH`) becomes
  `pH = run_result.pH()[-1]` for the most recent step — explicit,
  no hidden cache.
- Eliminates the per-step `_last_result` cache that today exists on
  `GasLiquidVolume` (and would otherwise have to live on
  `ControlVolume` under option A).

**Cons:**

- Significant API change.  Every caller of `run_batch` and every
  consumer of `BatchSimulationResult` updates.
- Mid-simulation reads (controllers reading pH between steps) need
  a different access pattern.  Today's pattern is `glv.pH` in the
  step callback; the equivalent under `RunResult` is either
  "controller receives the just-returned `AdvanceResult` directly"
  or "controller queries `run_result` from outside the loop" —
  neither is a one-line drop-in.
- `BatchSimulationResult` (in
  [../../models/vlmodels/fermenter/config/factory.py](../../models/vlmodels/fermenter/config/factory.py))
  already does some of this aggregation.  A `RunResult` overlaps
  with that dataclass; merging or replacing would be a non-trivial
  decision.
- Performance: derived accessors (`pH()` returning a numpy array)
  re-walk `advance_results` on each call.  For long runs with many
  reads this adds up — caching becomes necessary, which adds back
  some of the staleness concerns the design tries to eliminate.

## What it would change

**Source surface:**

- New `src/core/run_result.py` (or similar) holding the
  `RunResult` dataclass and accessors.
- [../../models/vlmodels/fermenter/config/factory.py](../../models/vlmodels/fermenter/config/factory.py):
  `run_batch` returns `RunResult`; `BatchSimulationResult` either
  becomes `RunResult` or is deprecated in favour of it.
- Step callbacks: today `step_callback(t_h, glv, advance_result)`
  receives the GLV; equivalent under `RunResult` is
  `step_callback(t_h, cv, advance_result)` plus a way to read
  cumulative history if needed (probably skip this — callbacks
  should be stateless).

**Caller surface:**

- Plotting / analysis code in `systems/` that reads parallel arrays
  off `BatchSimulationResult` switches to `run_result.pH()` etc.
- ~12 systems scripts touched, each a small mechanical change.
- Tests covering `run_batch` migrate to `RunResult` shape.

**Documentation:**

- [../architecture.md](../architecture.md) gains a "Run lifecycle"
  section.
- [../class_diagrams.md](../class_diagrams.md) gains `RunResult`
  in layer 2 alongside `AdvanceResult` and `MultiCVAdvanceResult`.

## Concrete sub-deliverables (added 2026-05-28)

Surfaced from the 2026-05-28 exploration session (intermediate
result logging + HPC checkpointing + result export questions).
The "richer default recorder" reframing maps to three concrete
pluggable `Recorder` implementations plus a format decision.

### `StreamingFileRecorder` — write to disk as the run progresses

A `Recorder` implementation that buffers per-step state in an
Arrow row-group and flushes a Parquet chunk to disk every N
steps. Default `flush_every_n_steps=1000` (configurable). On
crash / interrupt, the most recent flushed chunk is the
recoverable state; the in-flight buffer is lost.

Use cases:
- **Intermediate result logging.** A long simulation can be
  inspected during the run by reading the flushed chunks
  (which any Parquet-aware tool — pandas, Polars, DuckDB —
  can do).
- **Half of HPC checkpointing.** Combined with the
  `HPC_CHECKPOINTING.md` work, the flushed chunks form the
  time-series side of a resumable checkpoint; the `Simulation`
  state pickle is the orchestrator side.
- **Memory-bounded long runs.** Today's `BatchRecorder`
  pre-allocates the full result array — a 10⁶-step run keeps
  everything resident. Streaming caps the in-memory buffer at
  N steps.

### `SparseRecorder` — record every Mth step only

A `Recorder` implementation that records only every M-th step
(plus the initial and final states). Reduces in-memory footprint
and post-run analysis cost for runs where intermediate detail
isn't needed.

### `SummaryRecorder` — final-state only, plus accumulated counters

A `Recorder` implementation that captures only the initial and
final `CVSnapshot`, plus running counters (total moles in/out
per boundary, cumulative integral of selected channels). For
parameter-sweep / Monte-Carlo studies where only the endpoint
matters and storing thousands of full trajectories is wasteful.

### Format decision: Parquet, rotated chunks, JSON manifest

Pinned 2026-05-28 (efficiency analysis under the trigger
session). Rationale and the alternative-format table live in
the planning-note discussion summary; the decision is:

- **Parquet** as the on-disk format. Columnar, compressed
  (snappy default; zstd optional for long runs), interoperable
  with pandas / Polars / DuckDB / Spark.
- **Rotated chunks** — one Parquet file per flush, named with
  a monotonically increasing chunk index
  (`chunk_00000.parquet`, `chunk_00001.parquet`, …). Avoids
  the append-to-single-file ergonomics issues.
- **JSON manifest** (`manifest.json`) alongside the chunks,
  recording: schema version, CV/phase/species inventory, chunk
  boundaries (start/end step index per chunk), Simulation
  metadata at start-of-run (Python version, package version,
  `OMP_NUM_THREADS` for reproducibility cross-reference with
  HPC_CHECKPOINTING).
- **Reader helper** alongside the recorder: `load_run(path) →
  BatchResult` that reassembles chunks into a unified
  in-memory result. Symmetry with the `BatchRecorder` shape.

`pyarrow` becomes an optional dependency; the recorder raises
a clear `ImportError` pointing at `pip install pyarrow` if not
installed.

### Open design questions for the sub-deliverables

1. **Async flush.** Default sync; should the streaming
   recorder offer an async flush option via a background
   thread? Probably not for v1 (complicates crash semantics);
   revisit if I/O latency dominates real runs.
2. **Compression default.** Snappy (faster, larger files) or
   zstd (slower, smaller)? Lean snappy as default; zstd as
   `compression="zstd"` kwarg.
3. **Chunk size policy.** Step-count threshold (current
   recommendation: 1000) or memory-size threshold (e.g. 1 MB
   buffer)? Step-count is simpler; memory-size is more
   adaptive to wide vs narrow CVs.
4. **Per-CV-per-phase chunking.** One Parquet file per CV per
   phase per flush, or one file per flush combining all CVs?
   Per-CV-per-phase scales better with many CVs; combined is
   simpler. Probably combined for v1.

## Triggers — when to actually do this

Don't do this speculatively.  The right time is when one of these
pulls becomes concrete:

- **A controller framework that needs explicit cumulative history.**
  If [src/control/](../../src/control/) grows to support
  controllers that depend on rolling-window pH or ionic-strength
  averages, the parallel-array bookkeeping on
  `BatchSimulationResult` becomes a real friction point and a
  unified `RunResult` earns its keep.
- **Multi-CV history aggregation.**  Today
  `MultiCVSystem.advance_all` returns a `MultiCVAdvanceResult`
  per-step; aggregating those into a multi-CV time series is
  ad-hoc.  A `MultiCVRunResult` parallel to `RunResult` would slot
  in naturally.
- **Performance / memory introspection.**  If users start running
  long simulations and want to query "show me the pH trajectory
  from this run" without re-running, `RunResult` provides the
  right surface for that introspection.
- **Pressure to simplify the public API.**  Today's
  `BatchSimulationResult` is mostly parallel arrays for one CV;
  users running multi-CV systems already get a different result
  shape.  A unified `RunResult` could collapse this.

If none of these is in flight, the per-step `AdvanceResult` plus
ad-hoc orchestrator aggregation pays its way.

## Relationship to other phases

- [GLV_REMOVAL.md](GLV_REMOVAL.md) — Phase 7 explicitly deferred
  this idea (option B in §1) in favour of dropping the GLV
  accessors entirely (option C).  This file is the home for the
  deferred idea.
- [CHEMISTRY_UNIFICATION.md](CHEMISTRY_UNIFICATION.md) — installs
  `phase.pH` and `cv.current_pH()` with staleness guardrails.
  These are *per-step* mid-simulation reads; they are
  complementary to a `RunResult` rather than replaced by it
  (`RunResult.pH()` would still be valuable for plotting / analysis
  even with `phase.pH` available).
- [CONTAINER_LAYERING.md](CONTAINER_LAYERING.md) — the larger
  unbundling.  A unified `RunResult` is one of the layering
  improvements that this file's "concrete starting points" step 3
  flags ("Move to a `RunResult` object the orchestrator carries").
  Picking up `RunResult` would partially satisfy that step.

## Scope estimate

Smaller than [CHEMISTRY_UNIFICATION.md](CHEMISTRY_UNIFICATION.md),
larger than a typical small-phase like Phase 7:

- New dataclass + accessors (~150 lines).
- Migrate `run_batch` and `BatchSimulationResult` (one file).
- Migrate ~12 systems scripts (mechanical, mostly plotting code).
- Tests for `RunResult`.
- Documentation refresh.

Could ship as one PR with a checkpoint structure modelled on
[../shipped/PHASE5_CHECKLIST.md](../shipped/PHASE5_CHECKLIST.md)
or [../shipped/PHASE6_CHECKLIST.md](../shipped/PHASE6_CHECKLIST.md).

## Open design questions to resolve before starting

1. **`RunResult` vs `BatchSimulationResult` — replace or merge?**
   `BatchSimulationResult` already aggregates per-step data; a new
   `RunResult` overlaps with it.  Decision needed: deprecate
   `BatchSimulationResult` in favour of `RunResult`, or extend
   `BatchSimulationResult` to subsume `RunResult`'s accessor
   surface.
2. **Where does `RunResult` live?**  [../../src/core/](../../src/core/)
   if it's a core abstraction, or
   [../../models/vlmodels/fermenter/config/](../../models/vlmodels/fermenter/config/)
   if it's fermenter-specific.  Probably core, since it's solver-
   independent and would be useful for any CV-based simulation.
3. **Multi-CV variant.**  `MultiCVRunResult` aggregating
   per-CV per-step `AdvanceResult` objects — useful or
   over-engineered until there's a multi-zone driver?
4. **Lazy vs eager accessors.**  `pH()` returning a fresh numpy
   array each call vs. computing once and caching.  Caching is
   faster but reintroduces staleness concerns; lazy is simpler
   but slower for long runs.  Probably lazy by default with an
   optional `.materialize()` for callers that want the cached form.
5. **Step-callback signature.**  Today the callback receives the
   GLV; under `RunResult` it could receive the just-built
   `AdvanceResult` plus the CV.  Whether callbacks should also be
   able to read cumulative history is open — leaning toward "no",
   since callbacks should be stateless.
