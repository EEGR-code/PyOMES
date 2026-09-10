# Solvers — Summary

PyOMES has two orthogonal solver axes. **Axis 1** (this doc's main
subject) is per-CV physics: three `StepSolver` rungs that advance one
`ControlVolume` by one macro timestep `dt_h`. **Axis 2** (see
[Axis 2 — SystemSolver](#axis-2--systemsolver-whole-system-integration)
below) is whole-system integration across CVs, links, and controllers.
They compose orthogonally: `Simulation(solver=..., system_solver=...)`.

## Axis 1 — StepSolver: at a glance

| Solver | Semantics | Lives in | Where it earns its keep |
|---|---|---|---|
| `SequentialAdvanceSolver` (`solver=None` default) | Sequential operator split: speciation → feed → react → transfer; each sub-step sees post-previous state | [PyOMES/core/solvers.py](../PyOMES/core/solvers.py) | Default path; small `dt_h`; simple models without strong feed/reaction coupling |
| `SimultaneousEulerSolver` | Simultaneous explicit Euler: all sub-systems read the same frozen snapshot, deltas summed and applied together, with swappable clamping | [PyOMES/core/solvers.py](../PyOMES/core/solvers.py) | Larger `dt_h` where operator-splitting bias matters; matches the pre-refactor semantics that ADM1 was validated against |
| `SimultaneousAdaptiveSolver` | `solve_ivp` adaptive; explicit (DOP853) or implicit (Radau, BDF) methods; optional speciation freeze (DAE-style, BSM2/PyADM1 convention) | [PyOMES/core/solvers.py](../PyOMES/core/solvers.py) | Stiff models like BSM2 anaerobic digestion (fast H₂ kinetics); long macro timesteps that explicit Euler cannot survive |

All three implement the same `StepSolver` protocol (one method,
`solve_step(cv, dt_h, t_h, external_source_terms) -> AdvanceResult`)
and plug into `cv.advance(...)` via the `solver=` keyword — `solver=None`
is pure sugar for `SequentialAdvanceSolver()`, not a separate code path.
`Simulation(cvs=..., solver=...).run(tau_h, n_steps)` accepts a
single `StepSolver`, a `Dict[cv_key, StepSolver]` for per-CV
dispatch, or `None` (default — every CV uses
`SequentialAdvanceSolver`).

All three require only a `"liquid"` phase (speciation and reactions
are liquid-scoped) — a gas phase, solid phase, additional phases, or
none of those beyond liquid are all supported; a CV with no gas phase
runs under any of the three exactly as one with gas+liquid does.

`cv.advance()` on a CV already owned by a `Simulation` emits
`OrchestrationWarning` (a promotable `UserWarning`) — see
[Ownership guard](#ownership-guard-on-cvadvance) below.

---

## `SequentialAdvanceSolver` — sequential operator split (the default)

The foundation, and `solver=None`'s target — `cv.advance(dt_h, t_h)`
and `cv.advance(dt_h, t_h, solver=SequentialAdvanceSolver())` are
byte-identical calls. Promoted from ~70 lines of inline
`ControlVolume.advance()` logic to a real `StepSolver`-conforming class
so the default is reachable, wrappable, and placeable in a
config-driven dispatch dict like any other solver. Each call performs
five ordered sub-steps; each sub-step operates on the state left by
the previous one:

```
1. Run property calculators       →  phase.properties[key]
2. Run speciation engine          →  derived species in phase.n_mol (via _refresh_derived)
3. Apply external source terms    (feed, dosing, inter-CV transport)
4. Integrate reactions            (single forward Euler over dt_h; sees post-feed state)
5. Internal transfer              (every registered PhaseInterface, reads canonical n_mol)
```

This is a Lie-Trotter operator splitting with O(`dt`) error. The
ordering follows [ORDERING.md](phases-shipped/ORDERING.md): property
calculators and speciation run first so derived properties and
molecular fractions in `n_mol` are current for the reaction and
transfer sub-steps; feed runs before reactions because in a
continuously fed CV, material that arrives in the timestep is
immediately available for reaction. `cv.advance(dt_h, t_h)` is the
canonical signature post-state-unification — the legacy `chem_env`
dict is gone.

It is **not** equivalent to a single explicit-Euler step — see
[ORDERING_CRITIQUE.md](phases-shipped/ORDERING_CRITIQUE.md) Review 2. Reactions
see post-feed moles, transfer sees post-feed-post-reaction moles. For
substrate-limited Monod kinetics this is a feature; for cases where
the snapshot semantics are needed, use `SimultaneousEulerSolver`.

**Strengths:** simple, transparent, no extra dependencies, returns a
clean `AdvanceResult`. The behaviour is identical to what each
`PropertyCalculator`, `ReactionModel`, and `PhaseInterface` defines
on its own.

**Weaknesses:** the ordering of feed/reaction/transfer is baked in. No
adaptive step control.

**Clamping:** `clamp_fn=proportional_clamp` by default (each sub-step's
delta dict is clamped independently as it's computed — this solver is
sequential, not simultaneous, so there's no single combined dict to
clamp once the way `SimultaneousEulerSolver` does). Pass
`clamp_fn=floor_clamp`, a bespoke composite (see
[PyOMES/core/clamping.py](../PyOMES/core/clamping.py)'s module docstring for
the pattern), or `clamp_fn=None` to disable clamping entirely — the
last genuinely allows a species to go negative, which
`AccuracyMonitor.check_negative_mole` then flags. See
[Clamping](#clamping-non-negativity-handling) below.

### Ownership guard on `cv.advance()`

If a CV is owned by a `Simulation` (wired into `Simulation(cvs=...)`),
calling `cv.advance(...)` directly — rather than through
`sim.run(...)` — emits `OrchestrationWarning`: inter-CV links,
controllers, and profiles will **not** be applied for that step. This
remains a legitimate way to advance an *unowned* CV (unit tests,
notebooks exploring a bare CV's behaviour) and is the trusted internal
building block every `SystemSolver` calls (via
`cv._advance_unchecked(...)`, which bypasses the guard). For an owned
CV, prefer `sim.run(...)` — even for a single CV with zero links —
since that's where `t_h` bookkeeping, the recorder, and lifecycle
safety actually live. `OrchestrationWarning` is a `UserWarning`
subclass; promote or silence it the same way as `AccuracyWarning`
(`warnings.simplefilter("error", PyOMES.core.OrchestrationWarning)`).

---

## `SimultaneousEulerSolver` — simultaneous explicit Euler

A common solver passed to `Simulation.run` for fermenter-shaped
CVs. Phase-generic — covers whatever phases the CV has (only
`"liquid"` required), and every registered `internal_interfaces` entry
contributes, not just a single gas-liquid link. One macro step:

```
1. SNAPSHOT every phase the CV has (frozen pre-step state)
2. COMPUTE all sub-system fluxes from the snapshot, in any order:
     - boundaries (gas feed, vent, membrane)
     - speciation (engine.solve(phases=...) — writes derived n_mol)
     - property calculators (viscosity, density, …)
     - every registered internal interface (gas-liquid Henry, liquid-solid Ksp, …)
     - reaction model
3. COLLECT into per-phase delta dicts
4. CLAMP — clamp_fn scales removal fluxes so no species goes negative
5. APPLY all deltas to live state in one batch
6. Final speciation on post-step live state (for reporting)
```

All sub-systems read the same `t = t₀` state. No sub-system's output
feeds another within the same step. This is a true single-step
explicit Euler that differs from `cv.advance()` by O(`dt²`)
cross-terms.

**Strengths:**
- No operator-splitting bias between feed, reaction, and transfer.
- Swappable clamping prevents negative moles even with large `dt_h`
  and aggressive removal rates.
- Predictable cost per step (single evaluation of every sub-system).
- Matches the pre-refactor gas-liquid semantics that ADM1 and BSM2
  cases have been validated against.

**Weaknesses:**
- Still O(`dt`) overall — large `dt_h` produces visible truncation
  error.
- No adaptive control; the user must pick a stable `dt_h` themselves.
- For very stiff systems (BSM2 H₂, fast pH-coupled gas transfer)
  explicit Euler is unstable at any reasonable `dt_h` — fall back to
  `SimultaneousAdaptiveSolver`.

**Clamping:** `clamp_fn=proportional_clamp` by default, applied once to
the combined multi-source delta dict (step 4 above). Same
`clamp_fn=`/`None` options as `SequentialAdvanceSolver` — see
[Clamping](#clamping-non-negativity-handling) below.

---

## `SimultaneousAdaptiveSolver` — adaptive `solve_ivp`

Wraps `scipy.integrate.solve_ivp` to integrate the full system as a
single ODE over `[0, dt_h]`. Phase-generic, like `SimultaneousEulerSolver`
— only `"liquid"` is required. State is packed via
`StateVector` ([PyOMES/core/state_vector.py](../PyOMES/core/state_vector.py))
into a flat array covering every phase the CV has, alphabetically by
phase key then species key, with `H⁺` excluded as algebraic. One
derivative `f(t, y)` evaluates:

```
1. Unpack y → temp phase objects (one per phase the CV has)
2. Speciation (algebraic — engine.solve(phases=...) at every f evaluation,
   or cached if freeze_speciation; derived species land in phase.n_mol)
3. Boundary fluxes (instantaneous=True path)
4. Every registered internal interface (instantaneous=True path)
5. Reaction rates
6. Sum into dy/dt
```

`solve_ivp` then chooses sub-steps adaptively to satisfy `rtol` /
`atol`.

**Configuration:**
- `method`: `"DOP853"` (explicit, high-order, default), `"Radau"` /
  `"BDF"` (implicit, for stiff systems), `"LSODA"` (auto-switching),
  etc.
- `rtol` / `atol`: solver tolerances. Defaults `1e-6` / `1e-9`.
- `max_step`: caps the internal step (hours). Default `1.0`.
- `freeze_speciation`: if `True`, speciation is solved once at the
  start of the macro step and reused for all internal sub-steps. This
  matches BSM2/PyADM1's DAE convention where ion states are frozen
  during integration and updated algebraically after.

**Strengths:**
- Handles genuinely stiff systems that explicit Euler cannot
  (BSM2 anaerobic digestion with H₂ kinetics, pH-coupled CO₂
  transfer).
- Adaptive step control means much larger macro `dt_h` is acceptable.
- High-order methods (DOP853) give better accuracy per derivative
  evaluation when the system is non-stiff.

**Weaknesses:**
- Per-step cost is variable — many derivative evaluations may be
  required for stiff sections.
- No swappable `clamp_fn`; relies on `floor_nonnegative` (a raw
  state-array floor, conceptually distinct from `clamp_fn` — see
  [Clamping](#clamping-non-negativity-handling)) after each unpack.
  With aggressive removal rates the integrator may waste steps near
  zero.
- Diagnostics are sparser than the snapshot path (no per-boundary
  records, no transfer record) — this is intentional, since `solve_ivp`
  evaluates each sub-system many times per macro step and per-step
  records would be misleading.
- Pulls in `scipy` as a dependency for any model that uses it.

---

## Clamping — non-negativity handling

All three `StepSolver`s restore a physical invariant (non-negative
inventory) that the raw numerics don't guarantee on their own —
[PyOMES/core/clamping.py](../PyOMES/core/clamping.py) holds the shared,
swappable implementations:

| Function | Used by | Behaviour |
|---|---|---|
| `proportional_clamp` (default for both discrete-step solvers) | `SequentialAdvanceSolver`, `SimultaneousEulerSolver` | Scales a removal rate down so the result lands at exactly zero, preserving relative stoichiometry across whatever species share the delta dict |
| `floor_clamp` | Either discrete-step solver, passed explicitly | Per-species floor (optionally at `eps` instead of exact zero) — simpler, does not preserve relative stoichiometry |
| `floor_nonnegative` | `SimultaneousAdaptiveSolver` (always, not swappable) | Floors a raw ODE state array — a different failure mode (adaptive-integrator overshoot on a trial state), not a `(deltas, current_mol, dt_h)` decision |

`clamp_fn` is a plain callable, `clamp_fn(deltas, current_mol, dt_h) ->
deltas`, operating correctly on an **arbitrary subset** of species (not
whole-dict ownership) — this is what lets a bespoke composite match an
external reference model's specific per-species convention using the
shared functions as ingredients:

```python
def bsm2_style_clamp(deltas, current_mol, dt_h):
    result = dict(deltas)
    for sp in NEVER_CLAMP:                # species the reference model lets go negative
        result[sp] = deltas[sp]
    for sp in FLOOR_AT_EPSILON:            # reference model floors at eps, not exact zero
        result[sp] = floor_clamp({sp: deltas[sp]}, current_mol, dt_h, eps=1e-9)[sp]
    remaining = {k: v for k, v in deltas.items() if k not in NEVER_CLAMP | FLOOR_AT_EPSILON}
    result.update(proportional_clamp(remaining, current_mol, dt_h))
    return result

SimultaneousEulerSolver(clamp_fn=bsm2_style_clamp)
```

Passing `clamp_fn=None` to either discrete-step solver disables
clamping entirely — a genuine way to let a species go negative for
diagnostic purposes. `AccuracyMonitor.check_negative_mole` then flags
it post-step; independently, `AccuracyMonitor.check_clamp_invoked`
compares each delta dict before/after `clamp_fn` runs and flags
species whose rate was scaled down, naming the before/after values —
useful beyond visibility, since it signals that the chosen
`dt_h`/kinetics combination is aggressive enough that clamping is
doing real work.

---

## Axis 2 — SystemSolver: whole-system integration

Axis 1 (above) is per-CV physics. **Axis 2** is the orthogonal
question: given a whole `Simulation` (multiple CVs, inter-CV links,
controllers, profiles), how does one macro step get orchestrated?
Selected via `Simulation(system_solver=...)`; `None` (default) uses
the transparent `ExplicitEulerSystemSolver` path.

| Solver | Sequence per macro step | Earns its keep when |
|---|---|---|
| `ExplicitEulerSystemSolver` (default) | `profiles(t) → links(dt) → cv.advance(dt) → controllers` | τ_CFL ≫ dt_h and slow controllers — the common case |
| `StrangSplittingSystemSolver` | `profiles(t) → links(dt/2) → cv.advance(dt) → links(dt/2) → controllers` | Free 2nd-order splitting upgrade over Euler, ~zero extra cost |
| `MultirateSystemSolver` | `profiles(t) → cv.advance(dt) → [links(dt/M)] × M → controllers` | Fast inter-CV circulation (τ_CFL comparable to dt_h) without a global step-size cut |
| `ImplicitTransportSystemSolver` | Sparse implicit link-flux solve, unconditionally stable | Many-CV compartmental transport models (HPLC, packed beds) with no tight controllers |
| `MonolithicODESolver` | Single `solve_ivp` co-integrating CV species, inter-CV transport, and controller differential state, with outer-loop event scheduling at controller period boundaries | Tight feedback loops (T_c < dt_h), multi-rate digital controllers, co-integrated integral states (PI/PID, observers) |

All five implement `SystemSolver.advance_system(sim, dt_h, t_h)` and
call `cv._advance_unchecked(...)` internally (the trusted-orchestrator
path — this is why direct `cv.advance()` calls warn but
`sim.run(...)` never does), except `MonolithicODESolver`, which
integrates every CV directly via `cv.compute_rhs()` inside its shared
`solve_ivp` call and therefore **raises** if a per-CV `Simulation(solver=...)`
is also configured — there is no legitimate way for it to be silently
correct, since `MonolithicODESolver` never consults it.

Per-CV Axis-1 choice (`solver=`) and Axis-2 choice (`system_solver=`)
are otherwise fully orthogonal — pick any `StepSolver` per CV under
any `SystemSolver`. Full design rationale, the DAE/SUNDIALS Phase F/G
placeholders, and the "Identified extensions" (SIA, reactive D_eff) are
in [`docs/design/SOLVER_ARCHITECTURE.md`](design/SOLVER_ARCHITECTURE.md).

**Known duplication (documented, not fixed):** two independent
controller-event-firing mechanisms exist —
`EventScheduler` (priority-queue based, used by
`StrangSplittingSystemSolver`/`MultirateSystemSolver`/
`ImplicitTransportSystemSolver`) and `MonolithicODESolver`'s own inline
`ctrl_schedule` list. Both correctly implement periodic-controller
firing at `T_c` boundaries; they simply aren't the same object. No
fix scheduled — flagged here so a future solver addition doesn't
accidentally introduce a third mechanism instead of reusing one of
these two.

Writing your own `SystemSolver` is the same shape as writing your own
`StepSolver` — one method, composing `sim._apply_links`/
`cv._advance_unchecked`/`sim._invoke_controllers` in whatever order
your model needs. See
[demos/features/SolverProtocols/](../demos/features/SolverProtocols/) —
[`0_README.ipynb`](../demos/features/SolverProtocols/0_README.ipynb)
for the protocol overview, then
[`01_writing_a_custom_solver.ipynb`](../demos/features/SolverProtocols/01_writing_a_custom_solver.ipynb)
for a worked example at both axes.

---

## Choosing a solver

| If your model has… | Start with |
|---|---|
| Mild kinetics, small `dt_h`, no special needs | `cv.advance()` (`SequentialAdvanceSolver`, the default) |
| Aggressive feeds or strong feed/reaction coupling, larger `dt_h` | `SimultaneousEulerSolver` |
| Stiff fast kinetics (BSM2 H₂), pH-coupled transfer, long timesteps | `SimultaneousAdaptiveSolver(method="Radau")` or `SimultaneousAdaptiveSolver(method="BDF")` |
| Reproducing BSM2 / PyADM1 validation results | `SimultaneousAdaptiveSolver(freeze_speciation=True, method="Radau")` |

If a model works with the default and the results look right, there
is no reason to switch. The non-default solvers exist because real
models exist that the default cannot integrate accurately or stably.

---

## Accuracy + conservation monitoring

Every `ControlVolume` ships with two per-step monitors, both
attached on the `cv.reaction_system` surface and auto-wired in
`ControlVolume.__init__`. Each emits a `UserWarning` subclass when
a configured threshold is crossed.

`AccuracyMonitor` (chemistry-unification-4) — hooks the speciation
engine for solver-quality signals:

| Check | Hook point | What it flags |
|---|---|---|
| pH change per step | End-of-step speciation in both solvers | Operator splitting losing pH information; consider smaller `dt_h` or `SimultaneousAdaptiveSolver` |
| Newton iterations | `BisectionChemicalEquilibriumEngine.solve()` | Warm-start far from solution; sharp composition change or stiffness |
| Charge-balance residual | `BisectionChemicalEquilibriumEngine.solve()` | DAE drift after a kinetic step; splitting error larger than expected |
| Scipy step rejection rate | `SimultaneousAdaptiveSolver.solve_step()` | Adaptive solver fighting stiffness; consider `method="Radau"`/`"BDF"` |
| Ionic-strength regime | `BisectionChemicalEquilibriumEngine.solve()` | I exceeds the validity range of Davies (~0.5 mol/L) or the ideal-solution assumption (~0.10 mol/L) |

`ConservationMonitor` (state-unification C6) — element + charge
accounting per `cv.advance` step:

| Check | What it flags |
|---|---|
| Per-step element drift | Stoichiometric imbalance in a kinetic step (rate law producing/consuming atoms inconsistently) |
| Cumulative element drift | Slow leak across many steps that no single step trips on its own |
| Per-step charge drift | Mismatched cation/anion source terms (e.g. dosing one without its counter-ion) |
| Cumulative charge drift | Slow charge imbalance accumulating across the run |

The species registry is auto-populated from the reaction
stoichiometries; species in `phase.n_mol` without a registry entry
are skipped (typical for unnamed strong-ion lumps `S_cat`/`S_an`
with `atoms={}`).

Tune thresholds + throttle on the package singleton before
running:

```python
import PyOMES, warnings

# Fine-grained
PyOMES.config.warnings.pH_change_threshold = 0.5
PyOMES.config.warnings.throttle = "silent"

# Preset
PyOMES.config.warnings = PyOMES.WarningConfig.production()

# Promote to error / silence per-category via Python's standard
# warnings machinery
warnings.simplefilter("error", PyOMES.AccuracyWarning)
```

In CI / batch runs, the `VLSIM_WARNINGS` environment variable
selects from flat presets:

```bash
$ VLSIM_WARNINGS=silent     python run.py    # all AccuracyWarning suppressed
$ VLSIM_WARNINGS=verbose    python run.py    # every check emits every step
$ VLSIM_WARNINGS=production python run.py    # first 3 per category, then summary
```

After a run, `PyOMES.print_accuracy_summary()` and
`PyOMES.print_conservation_summary()` each print a one-line
emission count keyed by category (useful in `throttle="once"`
mode for "did anything fire?").

### Diagnostic mode — planned

Rigorous accuracy verification (step-halving Richardson
extrapolation, cross-solver comparison) is a planned future
entry point. The design note targets
`cv.run_with_diagnostics(t_end_h, dt_h, *, mode="step_halving")`
returning a per-species relative-error report; the
implementation is deferred until a concrete diagnostic workflow
appears. Use the cheap checks above in the interim; manually
running a model at `dt_h` and `dt_h/2` and comparing endpoints
covers most use cases.

---

## Solver promotion (shipped, Phase 6) + interface refinement (shipped)

All three `StepSolver`s live at the same layer.  `SimultaneousEulerSolver`
and `SimultaneousAdaptiveSolver` were promoted in
[Phase 6 of the CV refactor](phases-shipped/SOLVER_PROMOTION.md) so any
`ControlVolume` can opt into any strategy via `cv.advance(solver=...)`
— originally restricted to CVs with `"gas"` and `"liquid"` phases;
generalized to any CV with at least `"liquid"` in the
`step-solver-interface-refinement` phase (see
[STEP_SOLVER_INTERFACE_REFINEMENT.md](phases-shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md)).
The [`simulation-class` phase](phases-shipped/SIMULATION_CLASS.md)
(shipped 2026-05-27) exposes per-CV solver choice through
`Simulation(cvs=..., solver={cv_key: StepSolver})` — a single
`StepSolver` applies to every CV, a dict dispatches per CV, and a
future `MultiCVStepSolver` placeholder is reserved for a
system-level integrator (raises `NotImplementedError` until
implemented). Phase 7 had already deleted `GasLiquidVolume` (see
[phases-shipped/PHASE7_CHECKLIST.md](phases-shipped/PHASE7_CHECKLIST.md));
callers now hold a plain `ControlVolume` and either pass the
solver directly to `cv.advance(...)` or hand it to
`Simulation`.

The `StepSolver` protocol takes a `ControlVolume` and returns the
unified `AdvanceResult`, which carries `transfer_record` (a
`LinkFlowRecord` from one active internal-interface transfer — the
first, if more than one produced nonzero flux) and
`boundary_records` (per-boundary `ExternalFluxRecord`s).  The legacy
`GasLiquidAdvanceResult` was dropped in Phase 6.

Items shipped by `step-solver-interface-refinement` (2026-07):
`SequentialAdvanceSolver` reification, the `OrchestrationWarning`
ownership guard, `MonolithicODESolver` rejecting a per-CV `solver=` it
can never consult, generalizing both `Simultaneous*` solvers off the
gas/liquid assumption, the shared swappable `clamp_fn` module, the
`negative_mole`/`clamp_invoked` diagnostics, and unified state-vector
packing (`PyOMES/core/state_vector.py`). Items still deferred — the
multi-CV-aware solver tier (SIA, see
[SOLVER_ARCHITECTURE.md](design/SOLVER_ARCHITECTURE.md) "Identified
extensions"), orchestrator-level adaptive macro `dt_h`, decoupling
`StepSolver` from the concrete `ControlVolume` class, and the
composition/integrator/clamp_fn three-dimension factoring (item 5) —
are catalogued in
[STEP_SOLVER_INTERFACE_REFINEMENT.md](phases-shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md)
and [SOLVER_PROMOTION.md](phases-shipped/SOLVER_PROMOTION.md) "Future
considerations".
