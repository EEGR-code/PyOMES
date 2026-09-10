# Simulation Class — Checklist (SHIPPED)

> **Status: Shipped 2026-05-27** — 15 checkpoints landed across
> 19 commits on the `simulation-class` branch (C0 design lock-in
> through C14 deletion sweep, plus C15 doc-shipping). Final
> standalone suite: **977 / 0**. Tag: `simulation-class-shipped`.
>
> Per-checkpoint test deltas (793 baseline + 184 new from this
> phase — net -36 from C14's MultiCVSystem test deletion):
>
> | C# | Commit | Δ tests | What |
> |---|---|---|---|
> | C0 (×2) | b5f7d5e, 24518ef | — | Design lock-in (RunContext, Pattern 1, Pattern B) |
> | C1 | b2c42b5 | +22 | Simulation skeleton + RunContext |
> | C2 | 53442fd | +16 | CVSnapshot / SimulationSnapshot + builders |
> | C3 | 5b1e3b3 | +19 | Recorder protocol + BatchRecorder + BatchResult |
> | C4 | 66d0675 | +16 | Simulation.run / _step single-CV (bit-for-bit vs run_batch) |
> | C5 | 5d34b6c | +9 | Multi-CV + per-CV solver dispatch (bit-for-bit vs MultiCVSystem) |
> | C6 | 924417f | +27 | Lifecycle gating via RunContext |
> | C7 | 246edef | +8 | Controller protocol (Pattern 1 collapse) |
> | C8a | 47c9ae5 | +20 | PHController CV-native |
> | C8b | e7bbb2a | +16 | DOAgitationController CV-native |
> | C8c | 1e6b321 | +16 | DOCascadeController CV-native |
> | C8d | 513a927 | +16 | Pressure-relief × 3 CV-native |
> | C9 | 15bb957 | +11 | Wire controllers into Simulation.run |
> | C10 | aa16fe5 | +11 | Profile protocol + concrete profiles |
> | C11 | 0605555 | +11 | FermenterBuilder.build_simulation() |
> | C12 | 6d46c46 | 0 (demos) | Migrate batch_fermenter demo; delete _controllers.py |
> | C13 | a5b43ca | +2 | ADM1 uses Simulation orchestrator |
> | C14 | a29b50a | -36 | Delete MultiCVSystem + test_multi_cv.py |
> | C15 | this commit | 0 (docs) | Doc shipping; move docs/phases-{upcoming → shipped} |

Working checklist for the `simulation-class` branch. Conceptual
framing and resolved decisions are in
[SIMULATION_CLASS.md](SIMULATION_CLASS.md). This file is the
implementation log. Modelled on
[../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md).

## Goal

Introduce a new `Simulation` class as the **single orchestration
pathway** for all CV-based models in the package. The phase:

1. **Subsumes** `MultiCVSystem` (deleted) and absorbs the time-loop
   body of `run_batch` (deleted).
2. **Ports** the legacy `src/control/` controller machinery onto
   CV-shaped state (`FermenterState` → `CVSnapshot`).
3. **Promotes** controllers and profiles to first-class concepts
   with structured per-step records that show up in `BatchResult`.
4. **Makes the `Recorder` pluggable**, with a default
   `BatchRecorder` that reproduces today's `BatchResult` schema
   bit-for-bit plus new `controller_actions` and `profile_actions`
   channels.
5. **Adds hard `RuntimeError` lifecycle gating** during `.run()` —
   topology and parameters are frozen; controllers and profiles
   change state only through their formal protocols.

After this branch ships:

- New module [`src/core/simulation.py`](../../src/core/simulation.py)
  defines `Simulation`, the public construction API, `.run()`, and
  the private `_step()` primitive.
- New module [`src/core/snapshot.py`](../../src/core/snapshot.py)
  defines `CVSnapshot` and `SimulationSnapshot` dataclasses + the
  `build_snapshot` function consuming `(cv, advance_result)`.
- New module [`src/core/recorder.py`](../../src/core/recorder.py)
  defines the `Recorder` protocol, `BatchRecorder` (default), and
  `BatchResult` (preserving today's fields + new channels).
- [`src/control/`](../../src/control/) controllers consume
  `CVSnapshot`/`SimulationSnapshot` instead of `FermenterState` and
  return `ControlAction` directly from a single `compute(state, dt_h)`
  method. The legacy `Actuator` protocol and the intermediate
  `Commands` dataclass are deleted (see SIMULATION_CLASS.md
  decision 2).
- [`models/vlmodels/fermenter/profiles.py`](../../models/vlmodels/fermenter/profiles.py)
  exposes a `Profile` Protocol; concrete items implement
  `apply(t_h, sim) → ProfileRecord`.
- [`models/vlmodels/fermenter/config/builder.py`](../../models/vlmodels/fermenter/config/builder.py)
  `FermenterBuilder.build()` returns a `Simulation`, not a CV.
- All four fermenter demos use `Simulation(...).run(...)`. The
  parallel demo controllers in `demos/_controllers.py` are
  **deleted**.
- ADM1 and BSM2 builders construct a `Simulation` and call
  `.run()` directly. (Strong-ion seeding was migrated to
  `seed_adm1_strong_ions` / `seed_bsm2_strong_ions` at CV
  construction in STATE_UNIFICATION C4; the
  `make_*_chem_env_fn` factories no longer exist.)
- [`src/core/multi_cv.py`](../../src/core/multi_cv.py) **deleted**.
  [`run_batch`](../../models/vlmodels/fermenter/config/factory.py)
  **deleted**.
  [`make_controller_callback`](../../demos/_controllers.py)
  **deleted** (file deleted entirely).
- BSM2 reference golden test
  ([`test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py))
  passes bit-for-bit — the change is structural, not numerical.
- ~830–870 standalone tests passing (was 793 post-`STATE_UNIFICATION`;
  net +40 from new coverage minus deleted parallel-machinery tests).

## Out of scope

Everything reserved for later phases or trigger-gated deferral:

- **Legacy `Fermenter` (`CUFermentationSpeciation`) sunset.** The
  3388-line class in
  [`models/vlmodels/fermenter/unit.py`](../../models/vlmodels/fermenter/unit.py)
  with its BioSTEAM dependency and own `.simulate(feed, tau)`
  shape survives. Its sunset is its own future migration; the
  BioSTEAM wrapper needs rewriting against the new orchestrator
  first. Documented as the surviving legacy orchestration path.
- **`HPLCColumn.simulate()` migration.** HPLC's 1D-field topology
  doesn't fit the CV graph. CONTAINER_LAYERING.md acknowledges
  this. Documented exception.
- **`src/sim/` and `src/solvers/` packages.** Both are consumed
  exclusively by the legacy `CUFermentationSpeciation` orchestrator
  in [`models/vlmodels/fermenter/unit.py`](../../models/vlmodels/fermenter/unit.py)
  (`Ledger`, `FermenterState`, `RunResult`, `CoupledSolver`,
  `SolverFactory`). They live or die with the BioSTEAM-fermenter
  legacy island and are not consolidated by this phase.
  [`src/sim/state.py:FermenterState`](../../src/sim/state.py) is a
  separate class from
  [`models/vlmodels/fermenter/types.py:FermenterState`](../../models/vlmodels/fermenter/types.py)
  — the former lives in this island; the latter is the one
  `src/control/` reads and gets replaced by `CVSnapshot`.
- **`.iterate()` generator API.** Deferred. The pluggable Recorder
  covers streaming/sparse use cases without a generator. Easy to
  add later as a thin wrapper around `_step`.
- **`.run(..., reset=False)` resume capability.** Deferred. `.run()`
  always destructive in v1. Purely additive when added.
- **`.step(dt_h)` public method.** `_step` stays private. Tests and
  orchestrator internals use it; public users drive via `.run()`.
- **HPLC, legacy Fermenter, or BioSTEAM integration.** All
  external-to-Simulation.
- **Run-history unification.** Reframed by this phase to "ship a
  richer default Recorder"; the design note
  [`RUN_HISTORY.md`](RUN_HISTORY.md) stays trigger-gated.
- **Container-layering unit-physics descriptors (`GasLiquidUnit`,
  etc.).** This phase covers only the topology/orchestration
  axis; unit-physics axis stays in
  [`CONTAINER_LAYERING.md`](CONTAINER_LAYERING.md).
- **`chemistry-unification-3b`.** Independent of this phase; stays
  trigger-gated.

## Resolved decisions

The 12 design decisions are pinned in
[SIMULATION_CLASS.md §Resolved decisions](SIMULATION_CLASS.md).
Do not re-litigate in this branch. The most relevant for
implementation flow:

- `Simulation` subsumes `MultiCVSystem` entirely (Option C).
- Controllers receive a typed `CVSnapshot` / `SimulationSnapshot`,
  not raw `(cv, result)`.
- **Controller protocol collapsed (decision 2, reconciled
  2026-05-26).** Single `Controller.compute(snapshot, dt_h) →
  ControlAction`; the legacy `Actuator` protocol and `Commands`
  dataclass are dropped. Controller-internal state (setpoint,
  integral, sampling-hold) lives on `self`. See SIMULATION_CLASS.md
  decision 2 for the cross-model rationale.
- Records use parallel channels (`boundary_records`,
  `controller_actions`, `profile_actions`, `accuracy_records`, and
  `conservation_records`), not a unified intervention channel.
  `accuracy_records` shipped in chemistry-unification-4;
  `conservation_records` shipped in STATE_UNIFICATION C6
  (`ConservationMonitor`).
- **Lifecycle gating uses a `RunContext` object (resolved
  2026-05-26).** A small `RunContext` dataclass is planted on
  the Simulation at construction; every lockable owned object
  (CVs, Phases, `PhaseInterface` subclasses, `_LockableList`
  wrappers) holds a reference to *that one context*. Mutators
  consult `self._context.is_running`. One flag flipped once per
  `.run()`. Hard `RuntimeError`, not `UserWarning`. No
  back-references to the `Simulation` itself (no reference
  cycles); no distributed per-object `_simulation_running` bools
  (no flag drift). The legacy `KineticGasLiquidLink._simulation_running`
  flag and `_enter_simulation`/`_exit_simulation` methods are
  migrated to consume the unified RunContext at C6.
- `.run()` always destructive; no `reset` kwarg in v1.
- **`Simulation` owns `t_h`, no `chem_env` (decision 11,
  reconciled post-STATE_UNIFICATION).** `Simulation` accumulates
  wall-clock `t_h` and passes it directly to
  `cv.advance(dt_h, t_h)`. No `chem_env` / `chem_env_fn`
  parameters anywhere on the orchestrator surface. Strong-ion
  totals are regular `n_mol` species seeded at CV construction via
  `seed_bsm2_strong_ions` / `seed_adm1_strong_ions`.

## Checkpoints

Ordered so each leaves the test suite runnable. The skeleton
(1) lands first so later checkpoints can import. Snapshot
builders (2) and recorder (3) land next as standalone pieces.
`_step` + `.run()` single-CV (4) is the first user-facing
milestone. Multi-CV (5) and lifecycle gating (6) follow. The
controller port (7–9) is the largest sub-sequence and lands
in three checkpoints. Profiles (10), builder migration (11),
demo migration (12), ADM1/BSM2 migration (13), deletions (14),
and docs (15) close out.

**Commit discipline:** commit immediately after each checkpoint
finishes (suite green for that checkpoint's slice). Mid-flight
edits were silently lost between checkpoints during
chemistry-unification-3; per-checkpoint commits are
load-bearing, not optional. See
[../shipped/CHEMISTRY_UNIFICATION_3_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_3_CHECKLIST.md)
for the precedent.

### 1. `Simulation` skeleton

- [ ] New module
      [`src/core/simulation.py`](../../src/core/simulation.py)
      defining a small `RunContext` dataclass alongside the
      `Simulation` class:

      ```python
      @dataclass
      class RunContext:
          """Run-bound state shared across all objects owned by a
          single Simulation.run(). Owned objects hold a reference;
          mutators consult it to gate lifecycle-locked behaviour."""
          is_running: bool = False
          label: str = ""
          # Future fields (deferred until consumers appear):
          # current_t_h, current_step_index, recorder, ...
      ```

      The Simulation creates exactly one `RunContext` in `__init__`
      and propagates it to every owned mutable object at hookup
      time. Owned objects never reference the `Simulation`
      directly — only its `RunContext`. This avoids reference
      cycles, makes multi-Simulation scenarios trivially safe
      (each Sim has its own context), and gives owned objects
      one explicit dependency that's mockable in tests.

- [ ] `Simulation` class skeleton:
      - `__init__(self, cvs, links=None, controllers=None,
        profiles=None, solver=None, recorder=None, label="")` —
        no `chem_env` / `chem_env_fn` parameters (decision 11,
        post-STATE_UNIFICATION).
      - Field initialisation with `dict(cvs)`, `list(links or [])`,
        etc.
      - `self._context = RunContext(is_running=False,
        label=self.label)` — the run-bound state container.
      - For each CV in `self.cvs.values()`: assert `cv._context is
        None` (raise on double-ownership; see sanity check below),
        then `cv._context = self._context`. Wiring of nested
        objects (Phases, links, interfaces) is C6's problem;
        C1 only wires the CV layer.
      - Link-endpoint validation (carried verbatim from
        [`src/core/multi_cv.py:75-101`](../../src/core/multi_cv.py#L75)):
        every link references existing CV keys and existing phase
        keys on the source/sink CVs.
      - `_t_h: float = 0.0` accumulator — Simulation owns the
        wall-clock time and passes it to `cv.advance(dt_h, t_h)`
        each step (decision 11).
      - `_recorder` resolved at construction: if `recorder is None`,
        defaults to `None` and is filled in at `.run()` time with
        a `BatchRecorder()` instance (forward-compat for
        checkpoint 3).
      - Accessor helpers: `__getitem__`, `__contains__`,
        `cv_keys` property.
      - `total_mol()` — sum across all CVs/phases (carried from
        [`MultiCVSystem.total_mol`](../../src/core/multi_cv.py#L117)).
      - `snapshot()` — deep-copy each CV via `cv.snapshot()` (the
        snapshotted CV must have `_context = None`; the new
        Simulation that consumes the snapshot wires a fresh
        context). Share links/controllers/profiles/recorder/solver.
      - `__repr__` similar in style to `MultiCVSystem.__repr__`.

- [ ] Add `_context: Optional[RunContext] = None` attribute to
      [`ControlVolume`](../../src/core/control_volume.py) (one-line
      change in `__init__`; no behaviour yet — C6 adds the
      consuming guards).
- [ ] Update [`ControlVolume.snapshot()`](../../src/core/control_volume.py)
      so the returned CV has `_context = None`. Independence from
      any Simulation that previously owned the original.
- [ ] No `.run()` method yet. No `_step()` method yet.
- [ ] [`src/core/__init__.py`](../../src/core/__init__.py) — add
      `Simulation` and `RunContext` to exports. **Do not** remove
      `MultiCVSystem` from exports yet (deleted at checkpoint 14).
- [ ] New test module
      [`tests/standalone/test_simulation.py`](../../tests/standalone/test_simulation.py)
      with: construction tests, accessor tests, link validation
      tests (the four `KeyError` cases from `MultiCVSystem`),
      `total_mol` correctness, `snapshot()` independence (and
      verify the snapshotted CV's `_context` is `None`),
      `__repr__` shape, **`RunContext` wiring**: every CV's
      `_context` points at the Simulation's `_context` post-
      construction; double-ownership raises (constructing a
      second `Simulation` with a CV already owned by another
      raises). ~15 tests.

Sanity check: `from VLsim.core import Simulation, RunContext; sim = Simulation(cvs={"main": cv})` succeeds;
`cv._context is sim._context` and `cv._context.is_running is False`. Attempting
`Simulation(cvs={"main": cv})` a second time on the already-owned `cv` raises a clear `RuntimeError`.
`pytest tests/standalone/test_simulation.py -q` passes. Existing 793 tests still pass.

### 2. `CVSnapshot` and `SimulationSnapshot` dataclasses

- [ ] New module
      [`src/core/snapshot.py`](../../src/core/snapshot.py)
      with:

      ```python
      @dataclass(frozen=True)
      class CVSnapshot:
          cv_key: str
          t_h: float
          pH: Optional[float]
          ionic_strength: Optional[float]
          T_K: float
          V_liq_L: float
          V_gas_L: float
          P_gas_atm: float
          n_gas_mol: Dict[str, float]
          n_liq_mol: Dict[str, float]
          y_gas: Dict[str, float]
          species_properties: Dict[str, float]
          sensors: Dict[str, Any]

      @dataclass(frozen=True)
      class SimulationSnapshot:
          t_h: float
          cvs: Dict[str, CVSnapshot]
      ```
- [ ] `build_cv_snapshot(cv, cv_key, t_h, advance_result)` function
      that pulls fields off the CV's phases directly:
      `cv.phases["liquid"].pH` (catches `ValueError` if no
      `n_mol["H+"]` species, surfaces `None`),
      `cv.phases["liquid"].speciation.get("IonicStrength")` or
      `liquid.properties.get("ionic_strength")`,
      `liquid.n_mol` / `gas.n_mol` for species inventories,
      `liquid.properties` forwarded onto `species_properties`.
      Post-STATE_UNIFICATION the `AdvanceResult.properties` field is
      gone (deleted in C4b); pH and derived chemistry-state live on
      the CV's phases directly.
- [ ] `build_simulation_snapshot(sim, t_h, results: Dict[str, AdvanceResult])`
      function that aggregates per-CV snapshots.
- [ ] `sensors` field: derived quantities — `DO_mol_L` (O₂ in liquid
      ÷ V_liq_L), `T_C` (T_K − 273.15). Extension point for future
      derived properties; not coupled to any controller.
- [ ] Tests in `test_simulation.py`: snapshot building from a
      fresh CV, from a CV after one advance step, with and without
      speciation. ~10 tests.

Sanity check: `snap = build_cv_snapshot(cv, "main", t_h=0.0, advance_result=None)` returns a frozen
`CVSnapshot` with reasonable defaults. Round-trip via `dataclasses.asdict` and back. Frozen-ness
enforced: `snap.pH = 7.0` raises.

### 3. `Recorder` protocol + `BatchRecorder` default

- [ ] New module
      [`src/core/recorder.py`](../../src/core/recorder.py) with:

      ```python
      class Recorder(Protocol):
          def record_init(self, sim: "Simulation", n_steps: int) -> None: ...
          def record_step(
              self,
              step_index: int,
              t_h: float,
              dt_h: float,
              advance_results: Dict[str, AdvanceResult],
              controller_actions: List[ControlAction],
              profile_actions: List[ProfileRecord],
          ) -> None: ...
          def finalize(self) -> Any: ...
      ```
- [ ] `BatchRecorder` implementation:
      - Pre-allocates per-CV arrays at `record_init` based on
        `cv.phases[*].n_mol.keys()` and `n_steps`.
      - Records initial state (t=0) into index 0; per-step state
        into indices 1..n_steps.
      - Stores `transfer_records`, `boundary_records`,
        `controller_actions`, `profile_actions`,
        `advance_results` lists.
      - `finalize()` returns a `BatchResult` dataclass.
- [ ] `BatchResult` dataclass (move from
      [`models/vlmodels/fermenter/config/factory.py`](../../models/vlmodels/fermenter/config/factory.py)
      to `src/core/recorder.py`; preserve existing fields **plus**
      new parallel channels: `controller_actions: List[List[ControlAction]]`,
      `profile_actions: List[List[ProfileRecord]]`,
      `accuracy_records: List[List[AccuracyWarning]]` (shipped in
      chemistry-unification-4), and
      `conservation_records: List[List[ConservationRecord]]`
      (shipped in STATE_UNIFICATION C6 via `ConservationMonitor`).
      `ControlAction` and `ProfileRecord` stubs defined in
      `src/control/actions.py` (new module — note `src/sim/` is
      out-of-scope per the legacy-island deferral) and imported here.
- [ ] For multi-CV: `BatchResult.gas_mol` and `liquid_mol` become
      `Dict[cv_key, Dict[species, np.ndarray]]`. Single-CV
      simulations populate under `"main"` (or the sole key)
      — exact ergonomic shape settled here.
- [ ] Tests for `BatchRecorder` standalone (no Simulation
      needed): pre-allocation correctness, record_step
      bookkeeping, finalize shape. ~12 tests.

Sanity check: Construct a `BatchRecorder`, call `record_init` with a fake `Simulation` stub and `n_steps=10`,
verify array sizes. Call `record_step` 10 times with synthetic AdvanceResult inputs, verify time-series
contents. `finalize()` returns a `BatchResult` with the expected shape.

### 4. `Simulation._step` and `Simulation.run` — single-CV path

- [ ] On `Simulation`, implement `_step(dt_h: float, t_h: float) ->
      Dict[str, AdvanceResult]`:
      - For each CV in `self.cvs`, call
        `cv.advance(dt_h, t_h, solver=self.solver, ...)` —
        passing `t_h` directly (decision 11, post-STATE_UNIFICATION;
        no `chem_env=ctx` argument).
      - Single-CV: ignore `self.links` (empty list anyway).
      - Returns `{cv_key: AdvanceResult}`.
- [ ] On `Simulation`, implement `run(tau_h: float, n_steps: int)
      -> Any`:
      - Build default `BatchRecorder()` if `self._recorder is None`.
      - Call `recorder.record_init(self, n_steps)`.
      - Compute `dt_h = tau_h / n_steps`; build `t` grid.
      - Run initial speciation solve at t=0 to populate `pH[0]`
        and `I[0]` in the recorder (carried from
        [`run_batch`](../../models/vlmodels/fermenter/config/factory.py#L546-L557)).
      - `self._t_h = 0.0` (or carry over if a resume kwarg ever
        lands; v1 always destructive).
      - For each step `i in 1..n_steps`:
        - Call `results = self._step(dt_h, self._t_h)`.
        - `self._t_h += dt_h`.
        - Record via `recorder.record_step(i, t[i], dt_h, results,
          controller_actions=[], profile_actions=[])`.
      - Return `recorder.finalize()`.
- [ ] Runtime measurement: `time.perf_counter()` around the loop,
      stored on the result. (Carried from `run_batch`.)
- [ ] Newly-appearing species handling (carried from
      [`run_batch`](../../models/vlmodels/fermenter/config/factory.py#L607-L615)):
      if a species appears during the run, the recorder
      allocates an array on first observation.
- [ ] Tests in `test_simulation.py`: single-CV run bit-for-bit
      equivalent to `run_batch` for representative fermenter
      configurations. Compare against today's `BatchResult` output
      for at least 3 demos (batch, CSTR, fed-batch). ~15 tests.

Sanity check: For each fermenter demo, run today's `run_batch(cv, tau_h, n_steps)` and
`Simulation(cvs={"main": cv}).run(tau_h, n_steps)` on the same CV configuration. Assert all time-series
arrays match bit-for-bit. (At this point controllers are still external — demos that use `step_callback`
will fail this; that's checkpoint 9's job.)

### 5. Multi-CV path

- [ ] On `Simulation._step`, extend to apply inter-CV links **before**
      advancing CVs (carried from
      [`MultiCVSystem._apply_links`](../../src/core/multi_cv.py#L192-L224)).
      Returns `link_records: List[LinkFlowRecord]` alongside CV
      results.
- [ ] Per-CV solver dispatch:
      - If `self.solver` is a `StepSolver` instance, use it for all CVs.
      - If `self.solver` is a `Dict[str, StepSolver]`, use per-CV
        lookup (CVs absent from the dict default to the sequential
        body).
- [ ] On `BatchResult`, add `link_records: List[List[LinkFlowRecord]]`
      (parallel channel — list per step, like boundary_records).
- [ ] On `BatchRecorder`, record link flows per step.
- [ ] Tests: multi-CV simulation reproducing
      [`tests/standalone/test_multi_cv.py`](../../tests/standalone/test_multi_cv.py)'s
      cases on `Simulation` instead of `MultiCVSystem`. ~20 tests.
      Old `test_multi_cv.py` not yet deleted (it tests `MultiCVSystem`
      directly until checkpoint 14).

Sanity check: A two-CV system advancing with `MultiCVSystem.advance_all` for one step produces the same
final state as `Simulation(cvs={...}, links=[...])._step(dt_h, t_h)`. Per-CV solver dict works — one CV
uses Scipy, the other uses Euler; both advance.

### 6. Lifecycle gating

Built on the `RunContext` object planted at C1. The Simulation
flips one flag (`self._context.is_running`); every lockable
mutator across the package consults its own `self._context`
reference.

- [ ] `Simulation.run()` brackets: at entry, set
      `self._context.is_running = True`; in a `try/finally`,
      clear at exit. One flag, flipped once per run.
- [ ] Propagate `_context` to **nested** lockable objects. C1
      wired only the CV layer; C6 extends to:
      - `Phase` objects: when a `ControlVolume` is constructed
        (or when `_context` is set on it), forward to each phase
        via `phase._context = self._context`. Snapshot returns
        phases with `_context = None`.
      - `PhaseInterface` subclasses (notably
        [`KineticGasLiquidLink`](../../src/core/gas_liquid_link.py)):
        same forwarding pattern. The CV propagates `_context` to
        its `internal_interfaces` collection.
      - `_LockableList` wrappers for
        `ControlVolume.boundaries` and
        `ControlVolume.property_calculators`: the wrapper
        carries a `_context` reference and consults it on every
        mutating method (`append`, `pop`, `clear`, `__setitem__`,
        `extend`).
- [ ] Migrate
      [`KineticGasLiquidLink._enter_simulation` /
      `_exit_simulation` / `_simulation_running`](../../src/core/gas_liquid_link.py#L341)
      to the unified `RunContext` mechanism. Delete the
      per-object `_simulation_running` bool and the two enter/exit
      methods. `_warn_if_simulating(...)` becomes
      `_raise_if_running(action: str)` checking `self._context`,
      raising `RuntimeError` (no more `UserWarning` — clean break,
      no legacy callers).
- [ ] Locked surfaces (from
      [SIMULATION_CLASS.md §Lifecycle gating](SIMULATION_CLASS.md)):
      - `KineticGasLiquidLink.set_kLa`,
        `set_kLa_with_co2_ratio`, and any future
        `PhaseInterface` setter — `RuntimeError` via
        `_raise_if_running`.
      - `ControlVolume.boundaries.append/pop/clear` — wrap the
        list in `_LockableList` that consults
        `self._context.is_running`.
      - `ControlVolume.reaction_system` reassignment — property
        with guard (renamed from `reaction_model` in
        STATE_UNIFICATION).
      - `ControlVolume.property_calculators` mutation — same
        `_LockableList` pattern (renamed from `property_solvers`
        in STATE_UNIFICATION C5).
      - `Phase.V_L`, `Phase.T_K` setters — property with guard.
      - `Simulation.controllers`, `Simulation.profiles` mutation,
        `Simulation.solver`/`Simulation.recorder` reassignment —
        the Simulation guards its own surface by checking
        `self._context.is_running` (no `_context` indirection
        needed; the Simulation literally owns the context).
- [ ] Standardise the error shape across all lockable mutators:
      ```
      RuntimeError(
          f"Cannot mutate {qualified_attr} while Simulation"
          f" {self._context.label!r} is running. Use a "
          f"Controller or Profile to modify state mid-run."
      )
      ```
- [ ] For each gated public mutator, add an underscore-prefixed
      `_<name>_unchecked` sibling that bypasses the gate. The
      public mutator becomes `if self._context and
      self._context.is_running: raise ...; self._<name>_unchecked(*args)`.
      The orchestrator's controller-action apply path (C9) calls
      the unchecked variants; user code calls only the public
      versions. This is **Pattern B** of the privileged-mutation
      design space — see C9 below for the rationale (public /
      private split is the dominant idiom in pragmatic Python:
      pandas, SQLAlchemy, Cantera).
- [ ] Tests: every lockable mutator raises `RuntimeError` when
      called while its `_context.is_running` is `True`; the same
      mutator succeeds when `_context is None` or
      `_context.is_running is False`. Multi-Simulation isolation
      test: two `Simulation`s with disjoint CV sets run
      sequentially without bleed-over (one's `is_running` flip
      doesn't affect the other's CVs). ~15 tests.

Sanity check: A controller-free demo runs cleanly. A test that mutates `cv.boundaries.append(...)` mid-run
raises `RuntimeError` with a clear message naming the attribute and the Simulation's label. The same
mutation outside `.run()` works normally. After the run, `sim._context.is_running is False` and a fresh
mutation succeeds.

### 7. Controller port — interfaces + snapshot consumer

Applies decision 2 (Pattern 1 collapse): single
`Controller.compute(snapshot, dt_h) → ControlAction`. No
intermediate `Commands` type; no separate `Actuator` protocol.

- [ ] Rewrite
      [`src/control/interfaces.py`](../../src/control/interfaces.py)
      to consume `CVSnapshot` / `SimulationSnapshot` and return
      `ControlAction` directly:

      ```python
      class Controller(Protocol):
          def compute(
              self,
              state: CVSnapshot | SimulationSnapshot,
              dt_h: float,
          ) -> ControlAction: ...
          # Optional:
          sample_period_h: Optional[float]
          sample_period_s: Optional[float]
          def reset(self) -> None: ...
      ```

      The legacy `Actuator` protocol is **deleted**;
      `compute()` performs both the control logic and the
      rate-to-action translation. Controllers carry their
      internal state (setpoint, integral, derivative,
      last-sample time, held action) on `self`.
- [ ] Rewrite
      [`src/control/state_builder.py`](../../src/control/state_builder.py)
      to build `CVSnapshot` from `(cv, advance_result)` rather
      than `FermenterState` from a `Fermenter`. The old
      `FermenterState`-based path is deleted entirely (it's
      bound to legacy `Fermenter` which is out of scope for
      this branch; if it needs to survive, it survives in
      [`models/vlmodels/fermenter/unit.py`](../../models/vlmodels/fermenter/unit.py)
      as a private helper for the legacy class).
- [ ] **Delete the `Commands` concept entirely.** No replacement
      dataclass; no new `src/control/commands.py` module.
      [`vlmodels/fermenter/types.py:FermenterCommands`](../../models/vlmodels/fermenter/types.py)
      is no longer used by the new framework; its remaining
      legacy consumers in
      [`models/vlmodels/fermenter/unit.py`](../../models/vlmodels/fermenter/unit.py)
      (lines 340, 410, 634, 1079) are part of the
      CUFermentationSpeciation island. Final file disposition
      (delete vs unit-private inlining) settled at C14.
- [ ] Define `ControlAction` in **new** module
      `src/control/actions.py` (not `src/sim/control.py` — that
      module is part of the legacy `src/sim/` island, out of
      scope per the Out-of-scope list). Shape (extended from
      the legacy `ControlAction(vented_mol, dosed_mol)`):

      ```python
      @dataclass(frozen=True)
      class ControlAction:
          controller_label: str
          target_cv_key: str
          t_h: float
          dt_h: float
          flux_applied: Dict[str, Dict[str, float]]  # {phase: {species: mol/h}}
          params_changed: Dict[str, Any]              # {param_path: new_value}
          vented_mol: Dict[str, float]
          dosed_mol: Dict[str, float]
      ```
- [ ] Tests: `state_builder` produces correct `CVSnapshot` from
      a CV in known state. `ControlAction` shape. ~10 tests.

Sanity check: `build_cv_snapshot(cv, "main", t_h, advance_result)` produces a snapshot whose `pH`,
`DO_mol_L`, `V_liq_L`, etc. match direct reads off the CV. Old `FermenterState`-consuming code in
`src/control/` is removed (existing tests that depended on `FermenterState` now fail — they're
rewritten in checkpoint 9).

### 8. Controller port — concrete loops

- [ ] Migrate every concrete controller in
      [`src/control/loops.py`](../../src/control/loops.py) to read
      from `CVSnapshot` and return `ControlAction` (Pattern 1):
      - `PHController` — reads `state.pH`, `state.V_liq_L`,
        `state.sensors["CT_P_mol_L"]` (phosphate cap hint); returns
        a `ControlAction` with `flux_applied={"liquid": {acid_id: -mol/h}}`
        or equivalent base dosing. Internal state on `self`:
        setpoint, integral, last-sample time.
      - `DOAgitationController` — reads `state.sensors["DO_mol_L"]`;
        returns `params_changed={"link.kLa": new_kLa}` (or similar
        path).
      - `DOCascadeController` — reads `state.sensors["DO_mol_L"]`;
        returns `params_changed` mixing kLa, gas-feed composition,
        and VVM; carries tier-state on `self`.
      - `InstantPressureReliefController`,
        `SmoothPressureReliefController`,
        `PressureReliefController` — read `state.P_gas_atm`,
        `state.n_gas_mol`, `state.T_K`, `state.V_gas_L`,
        `state.y_gas`; return `flux_applied={"gas": {...}}` for
        vented moles. `vented_mol` populated for audit.
- [ ] Preserve all feature surfaces:
      - Sampling-period support
        ([`src/control/system.py:79-91`](../../src/control/system.py#L79))
      - Cascade tier-active diagnostics
        ([`src/control/loops.py:DOCascadeController`](../../src/control/loops.py))
      - matplotlib `plot_diagnostics()` on each controller
      - `enable_diagnostics` flag + diagnostic time-series storage
      - Registry validation of acid/base compound IDs
        ([`src/control/loops.py:152-166`](../../src/control/loops.py#L152))
      - `requires_segmented_ode` flag
      - `control_tags` for categorisation
      - `reset()` methods
- [ ] Restructure
      [`src/control/system.py:ControlSystem`](../../src/control/system.py)
      around the collapsed single-method shape. The legacy
      `step(state, cmd, dt_h, t_h)` → `cmd` plus
      `rates(state, cmd)` → `dict` plus
      `action(state, dt_h)` → `ControlAction` three-step
      collapses to one call per controller per step:
      `ctrl.compute(snapshot, dt_h)` → `ControlAction`. The
      `ControlSystem` retains sampling-period gating, cascade
      coordination, and registry validation, but stops shuttling
      a `Commands` dataclass between sub-calls.
- [ ] **Multi-CV-aware controllers**: a controller can declare
      `target_cv_key: Optional[str] = None`. If set,
      `ControlSystem` calls `controller.compute(sim_snap.cvs[key], dt_h)`
      with the per-CV snapshot. If None, controller receives the
      full `SimulationSnapshot` and picks its own CV (the rare
      cross-CV case).
- [ ] Tests in
      [`tests/standalone/test_controllers.py`](../../tests/standalone/test_controllers.py):
      port from `FermenterState`-based fixtures to
      `CVSnapshot`-based fixtures. Each controller now returns
      `ControlAction` directly; assertions check
      `action.flux_applied` / `action.params_changed` content.
      Cover cascade behaviour, sampling periods, anti-windup,
      diagnostics. ~50 tests rewritten.

Sanity check: A `PHController` with a setpoint of 6.5 and `CVSnapshot(pH=7.0)` returns a `ControlAction`
with positive `flux_applied["liquid"]["H3PO4"]` (or equivalent acid). `DOCascadeController` reports
`_active_tier=1` at low setpoints, escalates to tier 2 when RPM saturates; tier is visible on the
returned action's `params_changed`. Sampling-period controllers return an unchanged held action
between samples. `pytest tests/standalone/test_controllers.py -q` passes.

### 9. Wire controllers into `Simulation`

- [ ] In `Simulation._step` (and `.run()`'s per-step loop), after
      `cv.advance(...)` returns for each CV:
      1. Build `sim_snap = build_simulation_snapshot(self, t_h, results)`.
      2. Sample-period filter: for each controller, determine if
         it should fire this step (based on `sample_period_h` and
         last-fire time). Between samples, the controller's held
         `ControlAction` from `self` is reused (zero-order hold).
      3. Invoke each firing controller in one call:
         `action = controller.compute(sim_snap or per_cv_snap, dt_h)`
         — returns `ControlAction` directly (Pattern 1; no
         intermediate `Commands`, no separate `Actuator`).
      4. Apply the `ControlAction`:
         - `flux_applied` → `cv.apply_external_flux(phase, flux, dt_h)`.
         - `params_changed` → setter calls through the lifecycle-aware
           path (which is itself locked, but Simulation has owner
           rights — uses a privileged setter helper).
      5. Pass `controller_actions: List[ControlAction]` into the
         recorder's `record_step`.
- [ ] `Simulation` stores a `ControlSystem` internally (constructed
      from `self.controllers`). Reused per step.
- [ ] `Simulation.run()` calls `controller.reset()` on each
      controller at the top of every run (before
      `recorder.record_init`). `.run()` is always destructive
      (decision 10), so controller internal state — setpoint
      integrals, sampling-hold buffers, last-fire timestamps —
      resets unconditionally. Controllers without a `reset()`
      method are tolerated (skipped via `hasattr` check).
- [ ] Privileged mutation path for `params_changed` application
      (resolved 2026-05-26): use **Pattern B — public/private
      underscore split** (matches pandas, SQLAlchemy, Cantera
      and the dominant idiom in pragmatic scientific Python).
      Each gated mutator gains an `_<name>_unchecked` sibling
      that bypasses the `RunContext` gate; the public mutator
      (`set_kLa`, etc.) keeps the gate and delegates to the
      unchecked variant. The orchestrator's apply path resolves
      `param_path` strings (e.g. `"links.gas_liquid.kLa.O2"`)
      to the target object and calls the unchecked variant
      directly.

      Pattern B was picked over the alternatives:
      - **Pattern A (capability tokens on every mutator
        signature)** — too heavy for non-safety-critical sim.
      - **Pattern C (`with cv._context.allow_overrides():`
        scope-flip)** — the gate becomes porous: any `Profile.apply`
        or other mutation running inside the `with` window slips
        through inadvertently. No bounded blast radius.

      Concretely, every lockable mutator (C6 list) gets an
      underscore-prefixed unchecked sibling. New helper
      `_resolve_param_path(sim, path) -> (obj, attr, setter_name)`
      walks `"cv_key.attr.subattr..."` paths down `self.cvs` to
      the target object and returns the unchecked setter to
      call. Lives in `src/core/simulation.py` (private module
      helper). Tests cover the path-resolver alongside the apply
      loop.
- [ ] Tests in `test_simulation.py`: controllers integrated into a
      Simulation produce expected acid dosing and kLa changes;
      `controller_actions` appears in the BatchResult and contains
      structured ControlAction records (not empty); a sampling-period
      controller fires at the right times. ~25 tests.

Sanity check: A simulation with a `PHController(setpoint=6.5)` reaches and maintains pH ≈ 6.5. The
result has `len(result.controller_actions) == n_steps`, each entry containing a `ControlAction` with
populated `flux_applied`. A `DOCascadeController` test shows tier escalation in the recorded actions.

### 10. Profile protocol + concrete profiles

- [ ] In
      [`models/vlmodels/fermenter/profiles.py`](../../models/vlmodels/fermenter/profiles.py),
      define the `Profile` Protocol:

      ```python
      class Profile(Protocol):
          label: str
          def apply(self, t_h: float, sim: "Simulation") -> ProfileRecord: ...
      ```
- [ ] Migrate `ProfileSet` items (the concrete profiles inside
      today's class) to standalone classes implementing `Profile`.
      Each one's `apply()` mutates state through the privileged
      setter helper (so it doesn't trip lifecycle gating) and
      returns a `ProfileRecord(profile_label, t_h, targets)`.
- [ ] `ProfileSet` survives as an optional helper for grouping
      related profiles, but the orchestrator interfaces with a
      flat list of profiles.
- [ ] In `Simulation._step` (and `.run()`'s per-step loop), before
      controller invocation: for each profile, call
      `profile.apply(t_h, self)`, collect `profile_actions: List[ProfileRecord]`,
      pass to `recorder.record_step`. (No `ChemEnvProfile` —
      `chem_env` was deleted in STATE_UNIFICATION C4. If a future
      caller needs state-dependent strong-ion totals across a
      Scipy macro step — the one residual case noted in
      SIMULATION_CLASS.md decision 11 — a small dedicated
      `strong_ion_profile` hook is the additive-extension path,
      not a `chem_env`-override mechanism.)
- [ ] Tests: a `TemperatureRamp` profile changes
      `cv.phases["liquid"].T_K` between two specified times; the
      action appears in `BatchResult.profile_actions`. ~10 tests.

Sanity check: A demo simulation with a `TemperatureRamp` profile shows T_K rising over the configured
window; `result.profile_actions[i]` contains a record with `targets={"liquid.T_K": <value>}` for each
step in the ramp window.

### 11. `FermenterBuilder` migration

- [ ] Rewrite
      [`models/vlmodels/fermenter/config/builder.py:FermenterBuilder`](../../models/vlmodels/fermenter/config/builder.py)
      so `build()` returns a `Simulation`, not a `ControlVolume`.
- [ ] New fluent methods on `FermenterBuilder`:
      - `.controller(ctrl)` — append to controllers list.
      - `.profile(profile)` — append to profiles list.
      - `.solver(solver)` — set solver (one or per-CV dict).
      - `.recorder(recorder)` — set recorder.
- [ ] Drop any `.chem_env_fn(fn)` fluent method on the builder
      (decision 11; `chem_env_fn` was deleted system-wide in
      STATE_UNIFICATION C4e). Strong-ion seeding is a CV-construction
      concern via `seed_bsm2_strong_ions` / `seed_adm1_strong_ions`,
      not a builder-level fluent.
- [ ] `FermenterFactory.create_volume(...)` stays — it builds the
      CV; the builder wraps it in a Simulation.
- [ ] Tests in
      [`tests/standalone/test_builder.py`](../../tests/standalone/test_builder.py):
      build chain produces a `Simulation` with the expected
      controllers/profiles/solver/recorder attached. ~15 tests
      rewritten.

Sanity check: `sim = (FermenterBuilder().vessel(...).chemistry(...).controller(ph_ctrl).build())`
returns a `Simulation`. `sim.run(tau_h, n_steps)` produces a `BatchResult`.

### 12. Demo migration + delete `demos/_controllers.py`

- [ ] Rewrite each fermenter demo:
      - [`demos/batch_fermenter.py`](../../demos/batch_fermenter.py)
      - [`demos/cstr_fermenter.py`](../../demos/cstr_fermenter.py)
      - [`demos/fed_batch_fermenter.py`](../../demos/fed_batch_fermenter.py)
      - [`demos/microplate_fermenter.py`](../../demos/microplate_fermenter.py)

      Each demo replaces `run_batch(cv, ..., step_callback=...)` with
      `Simulation(cvs={"main": cv}, controllers=[...], solver=...).run(...)`.
      Controllers come from `VLsim.control` (the ported legacy ones),
      not the deleted demo controllers.
- [ ] **Delete**
      [`demos/_controllers.py`](../../demos/_controllers.py)
      entirely. `PHController`, `DOController`,
      `make_controller_callback` all go.
- [ ] Update demo imports to use `from VLsim.control import PHController, DOAgitationController`
      etc. instead of the deleted demo classes.
- [ ] Run each demo end-to-end after rewriting. Compare
      qualitatively against the previous output (numerical drift is
      expected for pH-controlled runs because the ported
      `PHController` has slightly different default tuning; the
      structural correctness is what's being verified). Tests in
      [`tests/standalone/test_factory.py`](../../tests/standalone/test_factory.py)
      need similar updates.

Sanity check: `python demos/batch_fermenter.py` produces a non-NaN time series with controllers acting.
The same for CSTR, fed-batch, microplate demos. No imports of the deleted `demos/_controllers.py`
remain in the repo (`grep -r "_controllers" demos/ tests/` returns empty).

### 13. ADM1 / BSM2 migration

- [ ] Update
      [`models/vlmodels/adm1/base.py`](../../models/vlmodels/adm1/base.py)'s
      docstring example and any test fixture using `run_batch(...)` →
      `Simulation(cvs={"main": cv}).run(...)`. No `chem_env_fn`
      parameter (deleted in STATE_UNIFICATION C4); strong-ion
      seeding already moved to `seed_adm1_strong_ions` at CV
      construction.
- [ ] Update
      [`models/vlmodels/adm1/bsm2.py`](../../models/vlmodels/adm1/bsm2.py)
      and
      [`tests/standalone/test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py).
      Strong-ion seeding via `seed_bsm2_strong_ions`. BSM2 reference
      golden test must pass **bit-for-bit** — the change is
      structural, not numerical.
- [ ] Tests: the BSM2 reference test passes against the new
      `Simulation` path. Re-baseline only if structural
      differences in float ordering produce numerical drift —
      should not happen, but document if it does.

Sanity check: `pytest tests/standalone/test_bsm2_reference.py -q` passes. `pytest tests/standalone -q`
passes overall.

### 14. Delete `run_batch` and `MultiCVSystem`

- [ ] Delete `run_batch` function from
      [`models/vlmodels/fermenter/config/factory.py`](../../models/vlmodels/fermenter/config/factory.py).
      `FermenterFactory.create_volume` stays. `BatchResult`
      lives in `src/core/recorder.py` (moved at checkpoint 3);
      delete the legacy import or alias here.
- [ ] Delete [`src/core/multi_cv.py`](../../src/core/multi_cv.py)
      entire module — `MultiCVSystem`, `MultiCVAdvanceResult`.
- [ ] Update
      [`src/core/__init__.py`](../../src/core/__init__.py) to
      remove `MultiCVSystem` and `MultiCVAdvanceResult` exports.
- [ ] Update `__init__.py` files anywhere `MultiCVSystem` was
      re-exported (e.g. `VLsim/__init__.py`).
- [ ] **Search** the repo for stale references to `run_batch`,
      `MultiCVSystem`, `MultiCVAdvanceResult`,
      `make_controller_callback` — all must be gone except
      historical-record mentions in
      [`docs/shipped/`](../shipped/).
- [ ] Delete
      [`tests/standalone/test_multi_cv.py`](../../tests/standalone/test_multi_cv.py)
      (its coverage migrated into `test_simulation.py` at
      checkpoint 5).
- [ ] Delete any `test_factory.py` tests that specifically tested
      `run_batch` (the equivalent Simulation tests live in
      `test_simulation.py`).
- [ ] Delete
      [`models/vlmodels/fermenter/types.py`](../../models/vlmodels/fermenter/types.py)
      (resolved 2026-05-26: remove if entirely legacy after the
      port). Both `FermenterState` and `FermenterCommands` lose
      their new-framework consumers during C7–C12. To enable
      deletion, inline `FermenterCommands` references in the
      remaining legacy
      [`unit.py`](../../models/vlmodels/fermenter/unit.py)
      call-sites (lines 340, 410, 634, 1079) — either as a
      private class at the top of `unit.py` or inline at each
      `getattr(unit, '_cmd', None) or FermenterCommands()`
      use-site. Small CUFermentationSpeciation edits are in
      scope at C14 specifically to enable the deletion (the
      class's runtime behaviour is unchanged).
- [ ] Update [`demos/_bootstrap.py`](../../demos/_bootstrap.py)
      docstring that imports `FermenterState` from the deleted
      module.

Sanity check: `grep -r "run_batch\|MultiCVSystem\|make_controller_callback\|vlmodels.fermenter.types" src/ models/ demos/ tests/standalone/`
returns empty. Test suite passes. Documentation in `docs/shipped/` still mentions these as
historical references; that's expected.

### 15. Documentation updates

- [ ] Update [`docs/architecture.md`](../architecture.md):
      - Replace the "Orchestrator" section (currently describing
        `_calc_ODE_with_headspace` from the legacy fermenter — a
        stale reference) with a "Simulation" section describing
        the new class.
      - Update the "Multi-CV Support" section to describe
        `Simulation` (not `MultiCVSystem`).
      - Sweep any other references to `run_batch`/`MultiCVSystem`.
- [ ] Update **both** class-diagram files in parallel — they are
      intentionally maintained side-by-side, full vs simplified:
      [`docs/class_diagrams.md`](../class_diagrams.md) and
      [`docs/class_diagrams_simple.md`](../class_diagrams_simple.md).
      Add `Simulation`, `RunContext`, `CVSnapshot`,
      `SimulationSnapshot`, `Recorder`, `BatchRecorder` in their
      natural layer. Remove `MultiCVSystem` /
      `MultiCVAdvanceResult` from both.
- [ ] Update [`docs/solvers.md`](../solvers.md): the existing
      reference to `cv.advance()` is unchanged, but
      "`run_batch()` uses `EulerSnapshotSolver` by default and
      accepts a `solver=` override" → "`Simulation(...).run(...)`
      uses `EulerSnapshotSolver` by default and accepts a
      `solver=` parameter at construction".
- [ ] Move
      [`docs/upcoming/SIMULATION_CLASS.md`](SIMULATION_CLASS.md)
      and this checklist to
      [`docs/shipped/`](../shipped/). Add a status
      banner at the top of each:
      `> **Status: Shipped YYYY-MM-DD** — N checkpoints landed on the simulation-class branch (N commits). M standalone tests passing. Tag: simulation-class-shipped.`
- [ ] Update
      [`docs/upcoming/README.md`](README.md):
      - Remove `SIMULATION_CLASS.md` from the priority list.
      - Update "Currently in flight" section. CONTAINER_LAYERING
        and RUN_HISTORY entries' references to the
        topology/orchestration axis can be tightened now that the
        work has actually shipped.
- [ ] Update
      [`docs/shipped/README.md`](../shipped/README.md)
      with the new entry.

Sanity check: Architecture doc no longer references `run_batch` or `MultiCVSystem` in current-state
descriptions. Class diagram reflects the new shape. The phase docs sit in `shipped/` with a
clear status banner.

## Ship

- [ ] Final test sweep: `pytest tests/standalone -q` green.
- [ ] BSM2 reference golden test green (bit-for-bit).
- [ ] All four fermenter demos run end-to-end.
- [ ] `git checkout main && git merge --no-ff simulation-class -m "Merge simulation-class: ..."`.
- [ ] `git tag simulation-class-shipped <commit-hash>`.
- [ ] `git push && git push --tags`.
- [ ] `git branch -d simulation-class && git push origin --delete simulation-class`.
- [ ] Update [`docs/upcoming/README.md`](README.md)
      "Currently in flight" to remove the simulation-class entry.
- [ ] Update auto-memory: status note in
      [`project_simulation_class.md`](../../C:/Users/k2473520/.claude/projects/c--Users-k2473520-VLcode/memory/project_simulation_class.md)
      moves from "design resolved, ready for branch" to "shipped
      YYYY-MM-DD". Update
      [`project_cv_refactor.md`](../../C:/Users/k2473520/.claude/projects/c--Users-k2473520-VLcode/memory/project_cv_refactor.md)
      and `MEMORY.md` similarly.

## Open questions deferred to checkpoint time

Implementation specifics, not architecture. Settled when the
corresponding checkpoint lands:

- **Exact field set on `CVSnapshot`** (checkpoint 2). Initial set
  matches the legacy `FermenterState` translated to CV-native
  shape; add fields as controllers need them.
- **Multi-CV `BatchResult` shape** (checkpoint 5). Per-CV channels
  keyed by CV key, or one big flat result? Likely per-CV.
- **`RunContext` eventual field set** (no fixed checkpoint). Today
  only `is_running` and `label`. Likely additions over time as
  consumers appear: `current_t_h` (moves from `Simulation._t_h`),
  `current_step_index` (today on `Recorder`), `recorder`
  reference (so monitors can ask "what step am I in?"). Each
  field lands when a real consumer wants it; no anticipation
  needed in C1.

Resolved during the 2026-05-26 design discussion (no longer
open):
- *Privileged mutation path for `params_changed`* — Pattern B
  (public/private underscore split). See C9.
- *`vlmodels/fermenter/types.py` disposition* — delete at C14,
  inline `FermenterCommands` into `unit.py` legacy callers to
  enable. See C14.
- *`controller.reset()` timing* — automatic at top of every
  `.run()` (matches always-destructive `.run()` semantic, decision
  10). See C9.
- *Lifecycle gating mechanism* — `RunContext` object held by
  every lockable. See C1, C6, and SIMULATION_CLASS.md decision 9.

None of these unblock the design; all are local decisions at
their corresponding checkpoints.
