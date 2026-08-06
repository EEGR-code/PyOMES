# Simulation Class — Design Note (SHIPPED)

## Status: Shipped 2026-05-27

> Implementation landed across 15 checkpoints on the
> `simulation-class` branch (commits `b5f7d5e` C0 through `a29b50a`
> C14, plus C15 doc-shipping). Tag: `simulation-class-shipped`.
> Final standalone suite: 977 / 0.
>
> Major outcomes:
> - New `Simulation` class is the single orchestration pathway for
>   all CV-based models.
> - `MultiCVSystem` deleted (subsumed at C14).
> - Pattern 1 controller protocol (single `compute(state, dt_h) →
>   ControlAction`); legacy `Commands` and `Actuator` deleted.
> - RunContext-based lifecycle gating across CVs, Phases, links,
>   `_LockableList` wrappers; legacy
>   `KineticGasLiquidLink._simulation_running` machinery deleted.
> - Pattern B unchecked-setter dispatch for orchestrator-mediated
>   mutation via the `_resolve_param_path` semantic resolver in
>   `Simulation._apply_param_change`.
> - `FermenterBuilder.build_simulation()` + new `.profile()` /
>   `.recorder()` fluents alongside the legacy `build()`.
>
> Legacy bits retained for the `CUFermentationSpeciation` island
> (out-of-scope per the design doc): `run_batch` in
> `models/vlmodels/fermenter/config/factory.py`, legacy controllers
> in `src/control/loops.py` / `src/control/system.py`,
> `models/vlmodels/fermenter/types.py`. Their deletion folds into
> the CUFermentationSpeciation sunset (separate future phase).

> ## ✅ Reconciled with shipped `STATE_UNIFICATION` (2026-05-22)
>
> Sequencing held: this phase ships **third**, after
> `INTEGRATOR_REMOVAL` (shipped 2026-05-20) and
> [`STATE_UNIFICATION`](../phases-shipped/STATE_UNIFICATION.md)
> (shipped 2026-05-22, tag `state-unification-shipped`). The two
> decisions flagged at design time have been reconciled below; the
> body of this design note has **not** been rewritten — instead
> each affected decision carries an inline reconciliation note
> overriding the original framing.
>
> - **Decision 11 (`chem_env` permanence) — superseded.**
>   `STATE_UNIFICATION` C4 deleted `chem_env` and `chem_env_fn`
>   system-wide. `cv.advance(dt_h, t_h)` is the canonical
>   signature. `Simulation` owns the wall-clock `t_h` accumulator
>   and passes it directly; there is no `chem_env` parameter, no
>   `chem_env_fn` pass-through, no auto-injection. Strong-ion
>   totals are regular `n_mol` species (named like `Na+`, `Cl-`,
>   or unnamed lumps `S_cat`/`S_an` with `atoms={}`); seeding goes
>   through `seed_bsm2_strong_ions` / `seed_adm1_strong_ions`
>   helpers (not via a per-step `chem_env_fn`). See the inline
>   reconciliation in decision 11.
> - **Decision 4 (records taxonomy) — extended.**
>   `STATE_UNIFICATION` C6 added `ConservationMonitor` (per-step
>   element + charge accounting), attached on
>   `cv.reaction_system.attach_conservation_monitor(monitor)` and
>   auto-wired in `ControlVolume.__init__`. When `Simulation` is
>   implemented, the parallel-channel taxonomy gains a
>   `conservation_records` channel alongside `accuracy_records`,
>   mirroring the existing `AccuracyMonitor` plumbing one-for-one.
>   See the inline reconciliation in decision 4.
>
> Other decisions (1, 2, 3, 5, 6, 7, 8, 9, 10, 12) are unaffected
> by `STATE_UNIFICATION` and still apply as written.

## What this is

Introduce a new `Simulation` class as the single orchestration
pathway for all CV-based models in the package. `Simulation`
subsumes today's `MultiCVSystem` (deleted), absorbs the time-loop
body currently in `run_batch` (deleted), promotes controllers and
profiles to first-class concepts with structured per-step records,
and ports the legacy `src/control/` controller machinery onto
CV-shaped state.

The motivation is concrete: today's demo controllers (in
[../../demos/_controllers.py](../../demos/_controllers.py)) bypass
the boundary recording surface — they mutate the CV directly via
`cv.apply_external_flux(...)` and `link.set_kLa(...)` and their
actions do not appear in any structured `BatchResult` field. The
legacy `src/control/` machinery in
[../../src/control/](../../src/control/) (`ControlSystem`,
`Controller`/`Actuator` protocols, `ControlAction`, sampling
periods, cascade DO controllers, matplotlib diagnostics, registry
validation) is feature-rich but orphaned — bound to the legacy
`Fermenter` (`CUFermentationSpeciation`) class and the
`FermenterState` dataclass, invisible to anything built on the CV
refactor (Phases 1–7).

This phase finishes the migration the CV refactor started: bring
the controller machinery onto CV-native state, make controllers
and profiles first-class on a single orchestrator class, and
delete the parallel demo machinery.

## Relationship to other phases / planning notes

- [CONTAINER_LAYERING.md](CONTAINER_LAYERING.md) — recommended
  preserving the three-axis layering (topology / unit-physics /
  integration) and keeping `MultiCVSystem` as a pure topology
  primitive. This phase **inverts that recommendation** for the
  topology/orchestration axis: `Simulation` owns topology
  internally rather than composing with a separate
  `MultiCVSystem`. The decision is deliberate (see *Resolved
  decisions* §1 below) and the CONTAINER_LAYERING note's "MultiCVSystem
  stays" stance is explicitly retired by this work.
- [RUN_HISTORY.md](RUN_HISTORY.md) — the pluggable `Recorder` shape
  in this phase is the architecturally correct home for unified
  run-history. Once `Simulation` ships, the run-history idea
  becomes "ship a richer default `Recorder`" rather than "rewrite
  the result schema."
- [STATE_UNIFICATION.md](../phases-shipped/STATE_UNIFICATION.md)
  (shipped 2026-05-22) — deleted `chem_env` / `chem_env_fn`
  system-wide; `cv.advance(dt_h, t_h)` is the canonical signature.
  `Simulation` owns the `t_h` accumulator and passes it directly,
  no per-step env dict. Strong-ion seeding moved to
  `seed_bsm2_strong_ions` / `seed_adm1_strong_ions` helpers at CV
  construction time.

---

## The shape of `Simulation`

### Construction

```python
class Simulation:
    def __init__(
        self,
        cvs: Dict[str, ControlVolume],
        links: List[CVLink] = None,
        controllers: List[Controller] = None,
        profiles: List[Profile] = None,
        solver: StepSolver | Dict[str, StepSolver] | "MultiCVStepSolver" = None,
        recorder: Recorder = None,
        label: str = "",
    ): ...
```

Post-`STATE_UNIFICATION`: no `chem_env` / `chem_env_fn`
parameters. `Simulation` owns the `t_h` accumulator and passes it
to `cv.advance(dt_h, t_h)` each step. Strong-ion totals are
regular `n_mol` species, seeded at CV construction via
`seed_bsm2_strong_ions(liquid, **kwargs)` /
`seed_adm1_strong_ions(liquid, **kwargs)`.

Single-CV usage:

```python
sim = Simulation(
    cvs={"main": cv},
    controllers=[ph_ctrl, do_ctrl],
    profiles=[temperature_ramp, vvm_schedule],
    solver=EulerSnapshotSolver(),
)
result = sim.run(tau_h=10, n_steps=1000)
```

Multi-CV usage:

```python
sim = Simulation(
    cvs={"sparger": cv_sparger, "bulk": cv_bulk},
    links=[liquid_circ_a_to_b, liquid_circ_b_to_a],
    controllers=[ph_ctrl],
    solver={"sparger": ScipyODESolver(method="Radau"),
            "bulk": EulerSnapshotSolver()},
)
result = sim.run(tau_h=48, n_steps=4800)
```

The link-endpoint validation that today's `MultiCVSystem.__init__`
performs (every link references existing CV keys and existing
phase keys) carries over directly.

### Public API

```python
class Simulation:
    # — Construction —
    def __init__(self, cvs, links, controllers, profiles, solver,
                 recorder, label): ...

    # — Properties / accessors —
    @property
    def cv_keys(self) -> List[str]: ...
    @property
    def total_mol(self) -> Dict[str, float]: ...  # sum across all CVs/phases
    def __getitem__(self, key: str) -> ControlVolume: ...
    def __contains__(self, key: str) -> bool: ...

    # — Run —
    def run(self, tau_h: float, n_steps: int) -> Result: ...

    # — Snapshot —
    def snapshot(self) -> "Simulation": ...

    # — Internal —
    def _step(self, dt_h: float) -> StepRecord: ...  # primitive used by .run()
```

Notes:

- **`.iterate()` is deferred.** The pluggable `Recorder` covers
  streaming/sparse use cases without requiring a generator API.
  `.iterate()` can land later as a thin wrapper around `_step` if
  a real workflow needs it.
- **`.run()` is always destructive.** No `reset` kwarg. State,
  recorder, and controller integrals reset at the top of every
  `.run()` call. Resume capability deferred until a use case
  demands it (additive change when needed).
- **`_step()` is private.** Public users drive simulations via
  `.run()`; tests and orchestrator internals can call `_step`
  directly.

### Result of `.run()`

The return type is whatever the `Recorder` produces from its
`finalize()` method. Default `BatchRecorder` returns a
`BatchResult` with the same schema as today's
[`models/vlmodels/fermenter/config/factory.py:BatchResult`](../../models/vlmodels/fermenter/config/factory.py)
**plus** new parallel channels for controller and profile actions:

```python
@dataclass
class BatchResult:
    # — Existing fields (preserved bit-for-bit) —
    t_h: np.ndarray
    gas_mol: Dict[str, np.ndarray]
    liquid_mol: Dict[str, np.ndarray]
    pH: np.ndarray
    ionic_strength: np.ndarray
    P_atm: np.ndarray
    transfer_records: List[LinkFlowRecord]              # was transfer_flow
    boundary_records: List[List[ExternalFluxRecord]]
    advance_results: List[AdvanceResult]
    runtime_s: float

    # — New fields (from this phase) —
    controller_actions: List[List[ControlAction]]       # per-step list per controller
    profile_actions: List[List[ProfileRecord]]          # per-step list per profile
```

For multi-CV simulations a `MultiCVBatchResult` variant (or
expanded `BatchResult`) records per-CV channels keyed by CV key.
Exact shape deferred to checkpoint 1 of the implementation, since
it depends on how the recorder pre-allocation handles the
single-vs-multi case.

---

## Controller protocol — port from `src/control/`

The CV-native controller protocol consumes a typed `CVSnapshot`
(replacing `FermenterState`) and returns a structured
`ControlAction` directly. There is no intermediate `Commands`
type and no separate `Actuator` protocol — the legacy
`compute(state, cmd) → cmd` plus `rates(state, cmd) → dict`
two-step collapses into a single method, with controller-internal
state (setpoint, integral, last-sample time, held command) carried
on `self`. See *Resolved decisions* §2 below for the design
rationale.

```python
@dataclass(frozen=True)
class CVSnapshot:
    """Read-only view of a CV at one timestep, given to controllers."""
    cv_key: str
    t_h: float
    pH: Optional[float]              # from cv.phases["liquid"].pH (None if no H+ species)
    ionic_strength: Optional[float]  # from liquid.speciation["IonicStrength"] / liquid.properties
    T_K: float
    V_liq_L: float
    V_gas_L: float
    P_gas_atm: float
    n_gas_mol: Dict[str, float]
    n_liq_mol: Dict[str, float]
    y_gas: Dict[str, float]
    species_properties: Dict[str, float]  # forwarded from liquid.properties (PropertyCalculator outputs); per-species derived moles now in n_liq_mol
    sensors: Dict[str, Any]               # extension point (DO_mol_L, etc. derived)

@dataclass(frozen=True)
class SimulationSnapshot:
    """Snapshot across all CVs for multi-CV controllers."""
    t_h: float
    cvs: Dict[str, CVSnapshot]

class Controller(Protocol):
    """Reads snapshot, returns structured action. State carried on self."""
    def compute(self, state: CVSnapshot | SimulationSnapshot,
                dt_h: float) -> ControlAction: ...
    # Optional:
    sample_period_h: Optional[float]
    sample_period_s: Optional[float]
    def reset(self) -> None: ...
```

The orchestrator builds a fresh `CVSnapshot` (or
`SimulationSnapshot` for multi-CV controllers) at every step from
the post-`cv.advance` state. Single-CV controllers receive
`CVSnapshot`; multi-CV controllers receive `SimulationSnapshot`
and pick which CV(s) to act on.

### ControlAction — the record

```python
@dataclass(frozen=True)
class ControlAction:
    """Auditable record of one controller's per-step action."""
    controller_label: str
    target_cv_key: str
    t_h: float
    dt_h: float
    # Flux-producing actions (e.g. pH dosing)
    flux_applied: Dict[str, Dict[str, float]]   # {phase_key: {species: mol/h}}
    # Parameter-adjusting actions (e.g. kLa change)
    params_changed: Dict[str, Any]              # {param_path: new_value}
    # Aggregated totals over this step
    vented_mol: Dict[str, float]
    dosed_mol: Dict[str, float]
```

`ControlAction` is the unified output type for both flux-producing
controllers (pH dosing → `flux_applied`) and parameter-adjusting
controllers (DO via kLa → `params_changed`). The legacy
[../../src/control/system.py:ControlSystem.action()](../../src/control/system.py#L99)
already produces a `ControlAction(vented_mol, dosed_mol)`; this
shape expands it to cover parameter mutations recorded as
structured deltas rather than direct setter calls.

### Sampling-period support

Keep verbatim from
[../../src/control/system.py:ControlSystem.step()](../../src/control/system.py#L56)
— if a controller exposes `sample_period_h` (or
`sample_period_s`), its `compute()` is only called at sampling
instants; between samples, the command is held (zero-order hold).

This matches real digital control hardware and lets different
controllers run at different rates in the same simulation
(e.g. fast pressure loop, slow pH loop). Implementation cost is
modest (~10 lines on the orchestrator side); the legacy code is
the reference.

### Sense / compute / record per step

```
1. Build CVSnapshot for each CV (post-advance state)
2. For each controller (respecting sample_period_h):
     action = controller.compute(snapshot, dt_h)   # returns ControlAction directly
3. Apply ControlAction:
     - flux_applied → cv.apply_external_flux(phase, flux, dt)
     - params_changed → setter calls through the lifecycle-aware path
4. Recorder.record_step(... controller_actions=[ControlAction, ...] ...)
```

Controllers carry their compute-vs-rates internal split themselves
(setpoint, integral, derivative, sampling-hold state are
controller-private). The framework-visible surface is one snapshot
in, one `ControlAction` out.

---

## Profile protocol

New `Profile` protocol modelled on the legacy
[../../models/vlmodels/fermenter/profiles.py:ProfileSet](../../models/vlmodels/fermenter/profiles.py)
but promoted to per-profile granularity and recording per action.

```python
class Profile(Protocol):
    label: str
    def apply(self, t_h: float, sim: Simulation) -> ProfileRecord: ...

@dataclass(frozen=True)
class ProfileRecord:
    profile_label: str
    t_h: float
    targets: Dict[str, Any]    # what was changed and to what value
```

Concrete profiles (`TemperatureRamp`, `VVMSchedule`,
`SetpointTrajectory`) implement `apply` and return a record
describing the mutation. The orchestrator stores
`profiles: List[Profile]` and iterates per step, recording each
action.

`ProfileSet` (the existing class) survives optionally as a helper
for grouping related profiles, but the orchestrator's primary
surface is the flat list.

---

## Recorder protocol

```python
class Recorder(Protocol):
    def record_init(self, sim: Simulation, n_steps: int) -> None: ...
    def record_step(
        self,
        step_index: int,
        t_h: float,
        dt_h: float,
        advance_results: Dict[str, AdvanceResult],   # per CV
        controller_actions: List[ControlAction],
        profile_actions: List[ProfileRecord],
    ) -> None: ...
    def finalize(self) -> Any: ...    # returns Result (BatchResult, MultiCVBatchResult, …)
```

Default `BatchRecorder`:
- Pre-allocates arrays sized `(n_steps + 1,)` at `record_init`
  based on each CV's `n_mol.keys()`.
- Stores arrays per-CV (single-CV uses cv_key `"main"` by
  convention; multi-CV uses the user's keys).
- Returns a `BatchResult` from `finalize()` with the schema
  documented above.

Future recorders (out of scope for this phase, but the shape
allows them):
- `SparseRecorder` — record every Nth step.
- `StreamingRecorder` — write to Parquet/HDF5 incrementally.
- `SummaryRecorder` — record only the final state plus summary
  statistics.
- `CompositeRecorder` — stack multiple recorders.

---

## Lifecycle gating

While `sim._context.is_running is True` (set at the top of
`.run()`, cleared on exit), the following raise `RuntimeError`:

- **Internal interface setters** —
  `KineticGasLiquidLink.set_kLa(...)`,
  `set_kLa_with_co2_ratio(...)`, and any future setter on any
  `PhaseInterface` subclass.
- **`cv.boundaries` mutation** — `append`, `pop`, `clear`,
  reassignment.
- **`cv.reaction_system` reassignment.**
- **`cv.property_calculators` mutation.**
- **Phase parameter changes** — `phase.V_L = ...`, `phase.T_K = ...`.
- **Simulation-level mutation** —
  `sim.controllers.append(...)`, `sim.profiles.append(...)`,
  `sim.solver = ...`, `sim.recorder = ...`.

### Mechanism: `RunContext`

A small `RunContext` dataclass holds the run-bound state shared
across all objects owned by a single `Simulation`:

```python
@dataclass
class RunContext:
    is_running: bool = False
    label: str = ""
    # Future fields (deferred until consumers appear):
    # current_t_h, current_step_index, recorder, ...
```

The Simulation creates exactly one `RunContext` in `__init__`.
Every lockable owned object (`ControlVolume`, `Phase`,
`PhaseInterface` subclasses, `_LockableList` wrappers) holds a
reference to that *one* context. Lockable mutators consult
`self._context.is_running`; the Simulation flips a single flag
once at run-entry and once at run-exit.

This generalises the legacy `_enter_simulation()` /
`_exit_simulation()` pattern on `KineticGasLiquidLink` — those
two methods and the per-object `_simulation_running` bool are
migrated to consume the unified `RunContext`. The existing
`UserWarning` becomes a hard `RuntimeError`.

Design properties:
- **No reference cycles.** Owned objects reference the
  `RunContext` leaf, never the `Simulation`.
- **Multi-Simulation safe by construction.** Two `Simulation`s
  with disjoint CV sets have disjoint `RunContext`s; one's
  `is_running` flip never affects the other's lockables.
- **One source of truth.** No distributed per-object
  `_simulation_running` bools that can drift if a code path
  forgets to flip one.
- **Test-friendly.** Construct `RunContext(is_running=True)`
  manually, assign to a CV/phase/link, verify the mutator
  raises — no mocking of `Simulation` required.
- **Snapshot semantics are clean.** `cv.snapshot()` returns a CV
  with `_context = None`; the snapshot is structurally
  independent and not bound to anyone's run until added to a
  new `Simulation`.

The principle: during a `.run()`, *only* the controller and
profile protocols can change anything. Both produce structured
records. Direct setter access from user code (e.g. a Jupyter cell
inadvertently mutating live state) is fenced off entirely.

---

## Migration map

### Deletions (clean breaks, no shims)

- [../../models/vlmodels/fermenter/config/factory.py:run_batch](../../models/vlmodels/fermenter/config/factory.py)
  — replaced by `Simulation.run()`. `FermenterFactory.create_volume`
  stays (it builds the CV).
- [../../src/core/multi_cv.py](../../src/core/multi_cv.py) entire
  module — `MultiCVSystem` and `MultiCVAdvanceResult` subsumed by
  `Simulation`. The `_apply_links` logic moves into
  `Simulation._step`. *(Alternative considered:
  `MultiCVAdvanceResult` survives as a renamed `MultiCVResult`
  type — see [RESULT_TYPE_RENAME.md](RESULT_TYPE_RENAME.md), an
  unfinished 2026-05-18 exploration. Locked design is the deletion
  above; resolve at implementation time if the alternative becomes
  attractive.)*
- [../../demos/_controllers.py](../../demos/_controllers.py) —
  demo `PHController`, `DOController`, `make_controller_callback`
  all deleted. Demos use the ported `src/control/` controllers
  instead.
- [../../models/vlmodels/fermenter/types.py:FermenterState](../../models/vlmodels/fermenter/types.py)
  and `FermenterCommands` — the live `FermenterState` (consumed
  by `src/control/` and demo `_bootstrap.py`) is replaced by
  `CVSnapshot`. `FermenterCommands` is dropped along with the
  `Commands` collapse (decision 2). Final file disposition (delete
  the whole module, or strip to whatever `unit.py` legacy
  consumers still need) settled at C14.
  *(Note: [../../src/sim/state.py:FermenterState](../../src/sim/state.py)
  is a separate class in the `CUFermentationSpeciation` legacy
  island; not touched here — see Out-of-scope.)*

### Rewrites

- [../../src/control/](../../src/control/) — migrate from
  `FermenterState` to `CVSnapshot` / `SimulationSnapshot`.
  Collapse `Controller.compute(state, cmd, dt_h) → cmd` plus
  `Actuator.rates(state, cmd) → dict` into a single
  `Controller.compute(state, dt_h) → ControlAction` (see
  decision 2). The legacy `Actuator` protocol and `Commands` type
  are deleted; controllers carry their internal state on `self`
  and return `ControlAction` directly. `ControlSystem`
  restructures around the same single-method shape, orchestrated
  by `Simulation`.
- [../../src/control/state_builder.py](../../src/control/state_builder.py)
  — rewrite to build `CVSnapshot` from `(cv, advance_result)`
  instead of `FermenterState` from a `Fermenter`.
- [../../src/control/loops.py](../../src/control/loops.py) — port
  `PHController`, `DOAgitationController`, `DOCascadeController`,
  `PressureReliefController`, `InstantPressureReliefController`,
  `SmoothPressureReliefController` to read from `CVSnapshot`.
  Compound-ID validation against the chemistry registry survives.
  Diagnostics (matplotlib plotting, time-series storage) survive.
- [../../models/vlmodels/fermenter/profiles.py](../../models/vlmodels/fermenter/profiles.py)
  — adapt `ProfileSet` items to implement the new `Profile`
  protocol with structured `ProfileRecord` output.
- [../../models/vlmodels/fermenter/config/builder.py:FermenterBuilder](../../models/vlmodels/fermenter/config/builder.py)
  — evolves to produce a `Simulation` directly. Fluent methods
  for controllers, profiles, solver, recorder. `build()` returns
  `Simulation`, not `ControlVolume`. (Post-`STATE_UNIFICATION`:
  no `chem_env_fn` fluent method; strong-ion seeding is a CV
  construction concern via `seed_*_strong_ions` helpers.)
- [../../demos/batch_fermenter.py](../../demos/batch_fermenter.py),
  [../../demos/cstr_fermenter.py](../../demos/cstr_fermenter.py),
  [../../demos/fed_batch_fermenter.py](../../demos/fed_batch_fermenter.py),
  [../../demos/microplate_fermenter.py](../../demos/microplate_fermenter.py)
  — switch from `run_batch(cv, ..., step_callback=...)` to
  `Simulation(...).run(...)`. Controllers come from
  `VLsim.control` (the ported legacy ones), not
  `demos/_controllers.py`.
- [../../models/vlmodels/adm1/base.py](../../models/vlmodels/adm1/base.py),
  [../../models/vlmodels/adm1/bsm2.py](../../models/vlmodels/adm1/bsm2.py)
  — replace `run_batch(cv, ...)` with
  `Simulation(cvs={"main": cv}).run(...)`. Strong-ion seeding
  already migrated in `STATE_UNIFICATION` C4 to
  `seed_adm1_strong_ions` / `seed_bsm2_strong_ions` at CV
  construction; no per-step `chem_env_fn` plumbing remains.
- [../../tests/standalone/test_multi_cv.py](../../tests/standalone/test_multi_cv.py)
  — rewrite to test `Simulation` multi-CV behaviour.
- [../../tests/standalone/test_controllers.py](../../tests/standalone/test_controllers.py)
  — migrate from `FermenterState` to `CVSnapshot`. Test cascade
  controllers, sampling periods, registry validation.

### Out of scope

- **Legacy `Fermenter` (`CUFermentationSpeciation`) in
  [../../models/vlmodels/fermenter/unit.py](../../models/vlmodels/fermenter/unit.py).**
  This class has a BioSTEAM dependency and its own
  `.simulate(feed, tau)` method. It survives; its sunset is a
  separate, BioSTEAM-aware migration. Documented as the surviving
  legacy orchestration path.
- **`HPLCColumn.simulate()` in
  [../../models/vlmodels/hplc/column.py](../../models/vlmodels/hplc/column.py).**
  HPLC packs its own state vector and has a 1D-field topology
  that doesn't fit the CV graph. CONTAINER_LAYERING.md
  acknowledges this. HPLC remains a documented exception until /
  unless it's re-expressed as a CV-graph (separate, large piece
  of work, no current pull).
- **`src/sim/` and `src/solvers/` packages.**
  Both are consumed exclusively by the legacy
  `CUFermentationSpeciation` orchestrator in
  [../../models/vlmodels/fermenter/unit.py](../../models/vlmodels/fermenter/unit.py)
  (`Ledger`, `FermenterState`, `RunResult`, `CoupledSolver`,
  `SolverFactory`). They live or die with the BioSTEAM-fermenter
  legacy island and are not consolidated by this phase. Note that
  [../../src/sim/state.py:FermenterState](../../src/sim/state.py)
  is a separate class from
  [../../models/vlmodels/fermenter/types.py:FermenterState](../../models/vlmodels/fermenter/types.py)
  — the former is the `CUFermentationSpeciation` snapshot, the
  latter is the one `src/control/` reads and gets replaced by
  `CVSnapshot` (see decision 2 and the deletion entry above).
- **Container layering more broadly.**
  [CONTAINER_LAYERING.md](CONTAINER_LAYERING.md)'s topology /
  unit-physics / integration separation is partially served by
  this phase (integration is already pluggable via `StepSolver`;
  topology is now owned by `Simulation`; unit-physics descriptors
  like `GasLiquidUnit` are not introduced). Full unbundling stays
  deferred.

---

## Resolved decisions

Captured here so a future maintainer doesn't re-litigate
settled questions.

### 1. `Simulation` owns topology (vs. composes with `MultiCVSystem`)

Three options were on the table during design:
- **A** — Rename `MultiCVSystem` → `Simulation`, extend with
  controllers/profiles/recorder.
- **B** — Keep `MultiCVSystem` as topology primitive, introduce
  new `Simulation` that composes with it.
- **C** — New `Simulation` class subsumes `MultiCVSystem` entirely;
  delete `MultiCVSystem`.

**Picked C.** Rationale: single conceptual unit easier to reason
about; no compositional seam to design around; clean break aligned
with the project's no-shims discipline. Inverts
[CONTAINER_LAYERING.md](CONTAINER_LAYERING.md)'s "MultiCVSystem
stays as topology primitive" recommendation — that note is
explicitly retired by this work for the topology/orchestration
axis. Three-axis layering for integration (`StepSolver`) is
preserved.

### 2. Controller protocol: port legacy `src/control/`, collapse `Commands` and `Actuator`

Three options at the protocol layer:
- **Q1.x** — Port legacy `src/control/` machinery onto CV state.
- **Q1.y** — New minimal CV-native protocol with
  `step(snapshot) → ControllerRecord`.
- **Q1.z** — Controllers as state-dependent boundaries.

**Picked Q1.x with one modification.** The legacy machinery is
feature-rich (sampling periods, cascade controllers, diagnostics,
registry validation) and represents serious controller design that
the codebase already paid for — the port brings these features
onto the CV path without re-implementing them. But the legacy
`Controller.compute(state, cmd) → cmd` plus
`Actuator.rates(state, cmd) → dict` two-step **collapses into a
single `Controller.compute(snapshot, dt_h) → ControlAction`**.
The intermediate `Commands` dataclass and the `Actuator` protocol
are dropped; controller-internal state (setpoint, integral,
sampling-hold) lives on `self`, and the framework-visible record
is the structured `ControlAction`.

Reasoning for the collapse (resolved 2026-05-26 during the
SIMULATION_CLASS implementation discussion):

- **Cross-model framework.** The legacy `FermenterCommands` names
  fermenter-shaped fields (`valve_open`, `air_mol_h`,
  `acid_mol_h`) that do not fit ADM1 / BSM2 / HPLC. Surveying
  serious process-simulation tooling, the dominant pattern is
  per-actuator typed objects (Cantera's `MassFlowController` /
  `Valve`, Modelica's typed connectors, IDAES `manipulated_vars`
  lists), not a flat shared `Commands` struct. VLsim's
  `ControlAction` with
  `flux_applied: Dict[phase, Dict[species, mol/h]]` and
  `params_changed: Dict[param_path, value]` is already a path-
  addressed action shape on the output side; symmetrising the
  input side by dropping `Commands` keeps the framework
  cross-model without an adapter layer.
- **Zero-order hold** between sampling instants is preserved by
  the controller storing its last action on `self` — no
  orchestrator-visible `Commands` is required.
- **No legacy leak.** `vlmodels/fermenter/types.py:FermenterCommands`
  was kept only to service the legacy port; deleting it removes
  fermenter-shaped vocabulary from the framework layer.

Sampling-period support kept (`sample_period_h` /
`sample_period_s`) — not legacy baggage, but the ideal approach
for honest digital-controller simulation.

### 3. State snapshot for controllers: typed `CVSnapshot`

Three options:
- Raw `(cv, result)` pair — controllers reach into structure.
- Typed `CVSnapshot` dataclass — named field access.
- Dict-keyed flat accessor.

**Picked typed snapshot.** Rationale: the named-field discipline
is what makes legacy controllers readable; preserving it is most
of the value of porting them. `SimulationSnapshot` wraps
`{cv_key: CVSnapshot}` for multi-CV controllers.

### 4. Records taxonomy: parallel channels

> **✅ Reconciled post-`STATE_UNIFICATION` (2026-05-22).**
> `STATE_UNIFICATION` C6 shipped `ConservationMonitor` (element +
> charge accounting per `cv.advance` step), attached on
> `cv.reaction_system.attach_conservation_monitor(monitor)` and
> auto-wired in `ControlVolume.__init__`. The parallel-channel
> taxonomy below is extended with a `conservation_records` channel
> alongside `accuracy_records`, mirroring the `AccuracyMonitor`
> plumbing one-for-one. The recorder default surfaces both
> summaries (`VLsim.print_accuracy_summary()` /
> `VLsim.print_conservation_summary()`) at the end of a `.run()`.

Three options:
- Parallel channels (`boundary_records`, `controller_actions`,
  `profile_actions`).
- Unified `intervention_records` with `kind` discriminator.
- Tiered (structured for some, flat for others).

**Picked parallel channels.** Each intervention type has different
structured content; flattening loses per-type schema. Easy to scan
chronologically via a derived accessor if needed.

### 5. Profile protocol: list of `Profile`, not `ProfileSet`

`Profile` becomes a Protocol with `apply(t_h, sim) → ProfileRecord`.
`Simulation` stores `profiles: List[Profile]`. `ProfileSet` survives
as optional grouping sugar. Symmetric with `controllers: List[Controller]`.

### 6. Recorder: pluggable, default `BatchRecorder`

Pros (schema flexibility, pre-allocation by topology, future
streaming/sparse recorders, future-proofing run-history
unification) outweigh cons (one more concept, default semantics
required). Default `BatchRecorder` reproduces today's
`BatchResult` schema bit-for-bit, plus new `controller_actions`
and `profile_actions` channels.

### 7. Step API: `.run()` only

Three options:
- `.run()` only.
- `.run()` + `.step()` (per-step primitive).
- `.run()` + `.step()` + `.iterate()` (generator).

**Picked `.run()` only** for v1. `_step()` exists as a private
primitive. `.iterate()` is deferred — most of its use cases are
served by a streaming/sparse recorder. Easy to add later if a
real workflow demands it (purely additive).

### 8. Solver: single OR per-CV dict (with type surface left open for a future third arm)

```python
solver=EulerSnapshotSolver()                        # one solver, all CVs
solver={"sparger": ScipyODESolver(...),             # per-CV
        "bulk": EulerSnapshotSolver()}
```

Carries forward today's `MultiCVSystem.advance_all(solvers=...)`
flexibility. Per-CV is theoretical today (no model uses multi-CV
in production), but cheap to support and real engineering cases
exist (sparger + bulk fermenter with different stiffness; two-stage
AD).

**Type-surface note (added with [MULTICV_STEP_SOLVER.md](../phases-upcoming/MULTICV_STEP_SOLVER.md)):**
The annotation is `StepSolver | Dict[str, StepSolver] | "MultiCVStepSolver"`
rather than `StepSolver | Dict[str, StepSolver]`. The third arm is a
**deliberate placeholder**: `MultiCVStepSolver` is a future protocol
for a system-level integrator that sees all CVs + all links at once
(see MULTICV_STEP_SOLVER.md for the architectural rationale —
compartmental fermenter models and HPLC-as-N-CVs both pull in this
direction, but neither is an active trigger today). No implementation
ships in this phase; the placeholder costs nothing now and prevents a
backwards-compatible-break later for any user doing
`isinstance(solver, (StepSolver, dict))`-style dispatch. Internal
dispatch in `Simulation._step` should branch on `(StepSolver, dict,
else)` so the third arm raises a clean `NotImplementedError` until
someone wires it through.

### 9. Lifecycle gating: hard `RuntimeError`, comprehensive lock, `RunContext` mechanism

During `.run()`, attempts to mutate locked state raise
`RuntimeError`. No `UserWarning` fallback (no legacy users to
appease). The set of locked mutators is documented above.

**Mechanism (resolved 2026-05-26 during the SIMULATION_CLASS
implementation discussion):** the gating is implemented via a
small `RunContext` dataclass that owned objects hold a reference
to (not via owner back-references to the `Simulation`, and not
via per-object `_simulation_running` flags). See *Lifecycle
gating* above for the full design. Three alternatives considered:

- **(a) Per-object flag** — each lockable holds its own
  `_simulation_running: bool` (today's pattern on
  `KineticGasLiquidLink`). Distributed state; the Simulation has
  to walk and flip N flags at every `.run()` boundary; risk of
  flag drift if a code path forgets one.
- **(b) Back-reference to `Simulation`** — each lockable holds
  `_owner_simulation: Simulation`; mutators consult
  `self._owner_simulation._is_running`. Introduces a reference
  cycle (`Sim → CV → Sim`), complicates GC and pickling, and
  forces a strict 1:1 CV-to-Simulation ownership model.
- **(c) `RunContext` (picked)** — each lockable holds
  `_context: Optional[RunContext]`, where `RunContext` is a small
  leaf dataclass created by the Simulation. Single source of
  truth (one flag flipped once); no reference cycles;
  multi-Simulation safe by construction; testable without
  mocking the `Simulation`; natural growth path (the
  `RunContext` can later carry `current_t_h`,
  `current_step_index`, `recorder`, etc. as consumers appear).

### 10. Resume / `reset` kwarg: deferred entirely

`.run()` always destructive. Add `reset: bool = True` later if a
real workflow needs it — purely additive change, existing code
unaffected.

### 11. `t_h` accumulator, no `chem_env`

> **✅ Reconciled post-`STATE_UNIFICATION` (2026-05-22).**
> `STATE_UNIFICATION` C4 deleted `chem_env` and `chem_env_fn`
> system-wide. The original framing of this decision (passing
> `chem_env` through to `cv.advance` as permanent plumbing) is
> superseded; the resolved form is captured here.

`Simulation` owns a wall-clock `t_h` accumulator and passes it
directly to `cv.advance(dt_h, t_h)` each step. There is no
`chem_env` parameter on `Simulation.__init__`, no `chem_env_fn`
pass-through, and no auto-injection machinery. Strong-ion totals
are regular `n_mol` species (named like `Na+`, `Cl-`, or unnamed
lumps `S_cat` / `S_an` with `atoms={}`); seeding goes through
the helpers `seed_bsm2_strong_ions(liquid, **kwargs)` /
`seed_adm1_strong_ions(liquid, **kwargs)` (shipped in
`STATE_UNIFICATION` C4) at CV construction, not via a per-step
`chem_env_fn`. `logH_guess` derives from `liquid.n_mol["H+"]` at
the previous step (the speciation engine handles this internally);
`T_K` reads from `Phase`.

The narrow case the original decision retained `chem_env_fn` for —
state-dependent strong-ion totals across a `ScipyODESolver` macro
step — is the *only* legitimate need still uncovered. If a real
workflow surfaces it post-`SIMULATION_CLASS`, a small dedicated
hook (`Simulation(strong_ion_profile=callable)` or a step-level
profile entry) is the additive-extension path. Until then,
`Simulation` does not carry the surface.

### 12. Naming: `Simulation`

Plain. Discoverable. No connotation collision worth worrying about.

---

## Test impact

Rough estimate:

- **New tests** for `Simulation` construction, single-CV run,
  multi-CV run, controller integration, profile integration,
  lifecycle gating, recorder pluggability. ~80–120 tests.
- **Rewritten tests**:
  [test_multi_cv.py](../../tests/standalone/test_multi_cv.py)
  (~30 tests, migrate from `MultiCVSystem` to `Simulation`),
  [test_controllers.py](../../tests/standalone/test_controllers.py)
  (~50 tests, migrate from `FermenterState` to `CVSnapshot`),
  [test_factory.py](../../tests/standalone/test_factory.py)
  (test_run_batch_* → test_simulation_run_*).
- **Existing demo tests**: BSM2 reference golden test
  ([test_bsm2_reference.py](../../tests/standalone/test_bsm2_reference.py))
  must continue to pass bit-for-bit — pure structural change for
  ADM1/BSM2, no numerical drift expected. Other demo / integration
  tests update mechanically.
- **Approximate total touched**: ~200–250 tests.

Current standalone test count is 793 (post-`STATE_UNIFICATION`,
shipped 2026-05-22). Expected post-phase: ~830–870 (net +40 from
new coverage minus some redundant tests deleted alongside demo
`_controllers.py`).

---

## Scope estimate

Larger than any single chemistry-unification phase. Comparable in
size to Phase 7 (GLV removal) plus a chemistry-unification phase
combined — touches:

- 1 new module (`src/core/simulation.py`, ~600 lines including
  `CVSnapshot`, `SimulationSnapshot`, `Recorder` protocol, default
  `BatchRecorder`)
- 1 new module (`src/control/snapshot.py` for the snapshot
  builder, ~150 lines)
- Rewrite of `src/control/` interfaces + state_builder + loops
  (~800 lines touched across 4 files)
- Profile protocol introduction (~100 lines in profiles.py)
- Deletion of `src/core/multi_cv.py`, `run_batch`,
  `demos/_controllers.py` (~700 lines deleted)
- Migration of 4 demo files, ADM1 + BSM2 builders, FermenterBuilder
- ~200–250 tests touched

Net code delta: ballpark −400 to −600 lines (deletions exceed
additions because the demo controller / `make_controller_callback`
/ `run_batch` body / `MultiCVSystem` adds up to ~1000 lines,
replaced by ~600 lines of `Simulation` + ~150 lines of snapshot
builder + ~150 lines of new `Profile` machinery).

---

## Checklist sketch — implementation checkpoints

Rough sequencing for the eventual `SIMULATION_CLASS_CHECKLIST.md`.
Each checkpoint is a commit (per the per-checkpoint discipline
established in chemistry-unification-3 / -4).

1. **Skeleton.** New `src/core/simulation.py` with `Simulation`
   class, validators, construction, `__repr__`, `snapshot()`.
   Just enough to import. No `.run()` yet.

2. **Snapshot builders.** `CVSnapshot`, `SimulationSnapshot`
   dataclasses; builder function consuming `(cv, advance_result)`.
   Standalone tested.

3. **Recorder protocol + `BatchRecorder` default.** Pre-allocation
   from topology + `n_steps`. Producing `BatchResult` shape
   bit-for-bit equivalent to today (without controller/profile
   channels yet).

4. **`Simulation._step` + `Simulation.run`.** Single-CV path only.
   Absorbs today's `run_batch` body. Hardcoded empty controllers
   / profiles. Passes existing demo simulations bit-for-bit
   (compared against current `BatchResult` outputs).

5. **Multi-CV path.** `_step` extended to apply links and advance
   each CV. Per-CV solver dispatch. `test_multi_cv.py` migrated.

6. **Lifecycle gating.** `_is_running` flag; `RuntimeError` on
   mutators. Wire each locked mutator to check the back-reference.

7. **Controller port — interfaces + snapshot consumer.** Migrate
   `src/control/interfaces.py`, `state_builder.py` from
   `FermenterState` to `CVSnapshot`. `ControlSystem.action` returns
   `ControlAction` (extended shape with `params_changed`).

8. **Controller port — concrete loops.** Migrate `PHController`,
   `DOAgitationController`, `DOCascadeController`,
   pressure-relief variants. Registry validation, diagnostics,
   sampling periods preserved.

9. **Wire controllers into `Simulation`.** Sense/compute/actuate/
   record sequence per step. Controller actions recorded in
   `BatchResult.controller_actions`.

10. **Profile protocol + concrete profiles.** `Profile` protocol;
    port `ProfileSet` items to per-profile granularity.
    `BatchResult.profile_actions` channel.

11. **`FermenterBuilder` migration.** Builder produces `Simulation`,
    not `ControlVolume`. Fluent methods for all new fields.

12. **Demo migration.** All four fermenter demos rewritten.
    `demos/_controllers.py` deleted.

13. **ADM1 / BSM2 migration.** Replace `run_batch(...)` calls with
    `Simulation(...).run(...)`. BSM2 reference golden test must
    pass bit-for-bit.

14. **Delete `run_batch` and `MultiCVSystem`.** Verify nothing
    else imports them; remove `src/core/multi_cv.py` and the
    `run_batch` function from `factory.py`. `BatchResult` survives
    (still produced by `BatchRecorder.finalize()`).

15. **Doc updates.** [../architecture.md](../architecture.md)
    section on the orchestrator rewritten around `Simulation`;
    [../class_diagrams_simple.md](../class_diagrams_simple.md)
    updated. This file (`SIMULATION_CLASS.md`) and its checklist
    move to [../phases-shipped/](../phases-shipped/).

Each checkpoint commits. Per the lesson from
chemistry-unification-3, per-checkpoint commits are
load-bearing, not optional.

---

## Open questions deferred to checkpoint time

These are implementation specifics, not architecture:

- **Exact field set on `CVSnapshot`.** What lives on `sensors`?
  Probably: `DO_mol_L` (derived from O₂ in liquid + V_liq),
  any extension fields. Settle when first controller migrates.
- **Multi-CV `BatchResult` shape.** Per-CV channels keyed by CV
  key, or one big flat result? Likely per-CV; settle at
  checkpoint 5.
- **`Commands` class shape for the ported protocol.** Today's
  `FermenterCommands` has a fixed field set; the ported version
  needs to be CV-agnostic. Could remain a dataclass with all-CV
  fields, or become a dict. Settle at checkpoint 7.
- **Strong-ion profiling post-port.** Post-`STATE_UNIFICATION`,
  `chem_env` is gone; the rare "state-dependent strong-ion totals
  across a `ScipyODESolver` macro step" use case (if any caller
  surfaces it) becomes either a `Profile` that mutates n_mol or
  a small dedicated `strong_ion_profile` hook. Settle at
  checkpoint 10.

None of these unblock the design; all are local decisions made
when the corresponding checkpoint lands.
