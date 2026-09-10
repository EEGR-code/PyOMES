# Solver Architecture — Two-Axis Integration Strategy

## Status

Design reference, 2026-06-11. Updated 2026-06-29.

**Phases A–E shipped** (2026-06-11/12). Phases F and G remain opt-in
placeholders; trigger conditions not yet met. Two extensions identified by the
mass-exchange architecture review (§ Identified extensions) are candidates for
the next solver phase.

This document supersedes `DAE_VS_OPERATOR_SPLITTING.md` and
`MULTICV_STEP_SOLVER.md`, which have been replaced. Companion doc:
[`MASS_EXCHANGE_ARCHITECTURE.md`](MASS_EXCHANGE_ARCHITECTURE.md) — covers
solver coupling from the reactive transport perspective (SNIA/SIA/DSA, reactive
D_eff, cv.advance() contract).

## Context

VLsim advances its physical system using **operator splitting**: at each macro
timestep, the speciation engine resolves fast equilibria algebraically, kinetic
reactions integrate forward, then inter-CV links transfer material between CVs.
This is the bioengineering field norm — PHREEQC, BSM2, and ADM1 reference
implementations all use the same approach — and is correct and efficient for
the primary use cases the codebase targets.

Two categories of model have been identified where the current approach is
structurally inadequate:

- **Transport-stiff multi-CV systems** (compartmental fermenter models, HPLC
  column cells): fast inter-CV circulation imposes a CFL stability constraint
  that forces impractically small macro `dt_h`.
- **Tightly controlled systems** (fast pH/DO feedback, multi-rate digital
  controllers): controller integral states need co-integration with the process
  for accuracy; discrete update clocks need correct event handling.

A two-axis architecture addresses both without changing the operator-splitting
default for the common case.

## The two axes

**Axis 1 — per-CV physics strategy**

How each `ControlVolume` advances its own intra-CV state: speciation solve,
kinetic integration, property calculators. Extension point: `StepSolver`
protocol, `cv.advance(solver=...)`.

| Implementation | Notes |
|---|---|
| Sequential Euler (default) | Current; operator splitting, O(dt) splitting error |
| `SimultaneousAdaptiveSolver` | Stiff per-CV ODE; handles fast kinetic sub-processes |
| `DAEStepSolver` *(Phase F, opt-in)* | Speciation as algebraic constraint; SUNDIALS IDA |

**Axis 2 — system integrator**

How the system advances all differential state jointly: inter-CV transport
coupling, controller integral states, and discrete controller update events.
Extension point: `SystemSolver` protocol, `Simulation._step()`.

| Implementation | Notes |
|---|---|
| Explicit Euler (default) | Current; CFL-bounded, controllers operator-split |
| `ExplicitEulerSystemSolver` *(Phase C)* | Same behaviour; validates protocol wiring |
| `StrangSplittingSystemSolver` *(Phase C)* | Second-order splitting; near-zero cost over Euler |
| `MultirateSystemSolver` *(Phase C)* | Adaptive subcycling; removes CFL instability |
| `ImplicitTransportSystemSolver` *(Phase D)* | IMEX; unconditionally stable for transport |
| `MonolithicODESolver` *(Phase E)* | Co-integrates CV + controller state; full multi-rate |
| `DAESystemSolver` *(Phase G, opt-in)* | System-wide DAE; SUNDIALS IDA |

The two axes are **orthogonal and composable**. Any Axis 1 choice can be
paired with any Axis 2 choice. `ImplicitTransportSystemSolver` (Axis 2) with
`SimultaneousAdaptiveSolver` per CV (Axis 1) is valid; so is `StrangSplittingSystemSolver`
with the default sequential Euler.

## Stiffness taxonomy

Three distinct stiffness types appear in VLsim models. They require different
treatments at different levels of the architecture and must not be conflated.

| Type | Source | Location | Addressed by |
|---|---|---|---|
| Speciation stiffness | Acid-base equilibria 10⁺ orders faster than kinetics | Intra-CV | Operator splitting (current); `DAEStepSolver` Phase F (opt-in) |
| Kinetic stiffness | Fast vs slow reactions within one CV | Intra-CV | `SimultaneousAdaptiveSolver` Axis 1 |
| Transport stiffness | Fast inter-CV circulation, τ_circ < `dt_h` | Inter-CV | `ImplicitTransportSystemSolver` Phase D |

## Why operator splitting is the right default for speciation

Speciation timescales are 10+ orders of magnitude faster than kinetic
bioreaction timescales (microseconds vs hours). The charge balance is
effectively instantaneous at kinetic resolution. Treating it algebraically —
explicit index reduction — is what PHREEQC, BSM2, and ADM1 reference codes do;
operator-split outputs are directly comparable to published reference
trajectories.

Alternatives are worse:
- Fake infinite rates → reintroduces the stiffness the algebraic solve avoids.
- Real DAE solver (SUNDIALS IDA) → replaces 1D brentq + Davies fixed-point
  with a harder convergence problem; loses `dt`-independence of equilibrium
  satisfaction; requires a C-compiler on installation.

A friction-free middle ground exists: `SimultaneousAdaptiveSolver` with
`freeze_speciation=False` re-runs the speciation engine at each ODE sub-step
(algebraic substitution rather than operator splitting). Slower per step,
no new dependency, appropriate when splitting error is measurable but SUNDIALS
is not warranted.

## Controller hybrid dynamics

Pattern 1 controllers (`compute(state) -> ControlAction`) currently update
once per macro step. Controllers with **integral state** (PI/PID, Luenberger
observers, MPC internal models) have differential equations governing that
state — `di/dt = e(t)` — and produce actions that depend on the current
integral value.

Two properties interact:
- The **integral state** is continuous: it should be co-integrated alongside
  CV species concentrations for accuracy in tight feedback loops.
- The **control output** is discrete: a digital controller fires at a fixed
  update period T_c and holds output constant (zero-order hold, ZOH) until
  the next update.

This is a hybrid continuous/discrete system. The system integrator (Axis 2)
owns both: it integrates controller state continuously as part of the system
state vector, and fires controller updates as events at T_c boundaries using
`scipy.solve_ivp(events=...)`.

For slow regulatory control (T_c >> `dt_h`), the current macro-step behaviour
is correct and requires no change. For tight loops (T_c < `dt_h`), the Phase E
`MonolithicODESolver` provides both co-integration and event scheduling.

Controllers are Axis 2 objects — they observe multi-CV state and act system-wide.
Their differential state belongs at the system integrator level, not inside any
individual CV's `advance()`. Multi-rate digital update clocks are handled by
event scheduling inside Axis 2; no third axis is required.

## Two tracks

```
Friction-free (pure scipy):   A → B → C → D
                                       ↘
                                         E

Opt-in (SUNDIALS):            A → F
                              E → G
```

Phases A–E deliver the practical target without new dependencies. Phases F and
G are placeholders whose interface seams are preserved by Phase A.

## Phase roadmap

| Phase | Name | Scope | Axis | Status |
|---|---|---|---|---|
| A | CV_COMPUTE_INTERFACE | CV + engine interface additions | 1 + 2 foundation | ✅ shipped 2026-06-11 (`cv-compute-interface-shipped`) |
| B | CONTROLLER_STATE_PROTOCOL | Controller differential-state + update-period protocol | 2 | ✅ shipped 2026-06-12 |
| C | SYSTEM_SOLVER_PROTOCOL | SystemSolver protocol, event scheduling, simple solvers | 2 | ✅ shipped 2026-06-12 |
| D | IMPLICIT_TRANSPORT | IMEX implicit-transport solver | 2 | ✅ shipped 2026-06-12 |
| E | MONOLITHIC_ODE | Co-integrated CV + controller ODE + multi-rate events | 2 | ✅ shipped 2026-06-12 |
| F | PER_CV_DAE *(opt-in)* | Per-CV speciation as DAE; SUNDIALS | 1 | ⏸ placeholder — trigger condition not met |
| G | SYSTEM_DAE *(opt-in)* | System-wide DAE across all CVs + controllers; SUNDIALS | 1 + 2 | ⏸ placeholder — depends on E + F |

## Identified extensions

Two solver extensions surfaced during the mass-exchange architecture review
([`MASS_EXCHANGE_ARCHITECTURE.md`](MASS_EXCHANGE_ARCHITECTURE.md) §10). Neither
was part of the original A–E roadmap; both fit naturally as Axis 2 additions.

### SIA — Sequential Iterative SystemSolver

The current default is **SNIA** (Sequential Non-Iterative): transport applied as
a pre-correction before `cv.advance()`, O(dt) splitting error. A
`SequentialIterativeSystemSolver` would iterate `(transport → cv.advance())`
until `‖Δy‖ < tol`, achieving O(dt²) accuracy at the cost of multiple CV
advances per macro step.

Scope:
- New `SystemSolver` implementation wrapping any inner `SystemSolver`
- Convergence loop over the existing `advance_system` → `cv.advance` path
- Tolerance and max-iteration parameters
- Compatible with any per-CV StepSolver (Axis 1 orthogonality preserved)
- Corrects the z-staleness issue in `cv.advance()` for large source-term steps

### Reactive inter-zone transport (D_eff correction)

Species in local chemical equilibrium within each zone cannot diffuse
independently without violating thermodynamic consistency. The correct treatment
is to transport conserved totals with an effective diffusivity:

```
D_eff,CT = Σ_i  (∂z_i/∂CT) × D_i
```

where `∂z_i/∂CT` is the column for CT in the gray-box Jacobian `∂z/∂y`,
now exposed by `GrayBoxEngineProtocol.jacobian_dz_dy()` (shipped in
`speciation-engine-protocols-shipped`).

Scope:
- Extend `DispersiveFlow.compute_flow()` to accept an optional engine reference
- When the engine satisfies `GrayBoxEngineProtocol`, compute D_eff per conserved
  total before assembling the flux dict
- Requires MonolithicODE (Phase E) or SIA iteration to be consistent
- Activated by a flag or by the presence of the gray-box engine; default
  behaviour (per-species kLa without D_eff) unchanged

---

## Phases F and G — placeholders

Both require SUNDIALS (`scikits.odes`) and are only warranted after Phases A–E
are complete and a concrete trigger has surfaced.

**Phase F — PER_CV_DAE**

Trigger: operator splitting demonstrably accumulates error faster than tolerance
for a specific model — detectable by comparing trajectories at `dt_h` vs
`dt_h / 2` and observing non-negligible divergence in a quantity of interest.

Scope when triggered:
- `DAEStepSolver` implementing `StepSolver` via `scikits.odes` IDA
- Uses `cv.compute_differential_rhs` + `cv.compute_algebraic_residual`
  (both stubbed in Phase A) and `engine.algebraic_species()` (added in Phase A)
  for differential/algebraic state-vector split
- Initial condition consistency: project onto constraint surface at `t=0` via
  existing speciation engine
- Optional install group: `pip install vlsim[dae]`

**Phase G — SYSTEM_DAE**

Trigger: Phase F in use AND tightly coupled multi-CV speciation is the binding
bottleneck.

Scope when triggered:
- `DAESystemSolver` extending `SystemSolver` protocol
- All CVs' `compute_algebraic_residual` forms part of the full system residual;
  no per-CV index reduction
- Controller differential and algebraic states co-integrated in the joint system
- Requires `compute_differential_rhs` activated across all CVs (Phase A stub
  promoted to a real implementation)

## Resolved positions

- Operator splitting is the correct default for speciation stiffness. Per-CV
  DAE (Phase F) is an opt-in for edge cases only, with a concrete trigger
  condition required before scheduling.
- The two axes are orthogonal. Per-CV physics choice (Axis 1) and system
  integration strategy (Axis 2) are independent decisions made per-model.
- Controllers are Axis 2 objects. Their differential state belongs at the system
  integrator level. Multi-rate digital update clocks are handled by event
  scheduling inside Axis 2 without a third axis.
- The friction-free track (Phases A–E) should be complete before any SUNDIALS
  work begins.
- Phase A's additions (`compute_differential_rhs`, `compute_algebraic_residual`
  stub, `engine.algebraic_species`) preserve the Phase F/G paths without
  requiring any revisit of the CV interface when they are eventually scheduled.

## Cross-references

- [`MASS_EXCHANGE_ARCHITECTURE.md`](MASS_EXCHANGE_ARCHITECTURE.md) —
  companion doc covering solver coupling from the reactive transport side:
  SNIA/SIA/DSA, reactive D_eff, cv.advance() contract and engine scope
- [`SPECIATION_REACTIONMODEL_BOUNDARY.md`](SPECIATION_REACTIONMODEL_BOUNDARY.md) —
  design reference on the BisectionChemicalEquilibriumEngine ↔ ReactionModel boundary; runtime
  routing detail and naming asymmetry rationale
- [`CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md`](CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md) —
  engine/solver split; black/gray/white-box engine protocols;
  `GrayBoxEngineProtocol.jacobian_dz_dy()` and z-update strategies for the
  reactive D_eff extension
- [`../shipped/CV_COMPUTE_INTERFACE.md`](../implementation/shipped/CV_COMPUTE_INTERFACE.md) — Phase A
- [`../shipped/CONTROLLER_STATE_PROTOCOL.md`](../implementation/shipped/CONTROLLER_STATE_PROTOCOL.md) — Phase B
- [`../shipped/SYSTEM_SOLVER_PROTOCOL.md`](../implementation/shipped/SYSTEM_SOLVER_PROTOCOL.md) — Phase C
- [`../shipped/IMPLICIT_TRANSPORT.md`](../implementation/shipped/IMPLICIT_TRANSPORT.md) — Phase D
- [`../shipped/MONOLITHIC_ODE.md`](../implementation/shipped/MONOLITHIC_ODE.md) — Phase E
- [`../../src/core/control_volume.py`](../../src/core/control_volume.py) —
  `ControlVolume.advance()` and `compute_rhs()`
- [`../../src/core/solvers.py`](../../src/core/solvers.py) —
  `StepSolver` protocol, `SimultaneousAdaptiveSolver`
- [`../../src/core/system_solver.py`](../../src/core/system_solver.py) —
  `SystemSolver` protocol, shipped Axis 2 implementations
- [`../../src/core/simulation.py`](../../src/core/simulation.py) —
  `Simulation._step()` dispatch
- [`../../src/core/links.py`](../../src/core/links.py) —
  `DispersiveFlow` (target for reactive D_eff extension)
