# Class Diagrams — Simplified

Compact, interface-only class diagrams of the PyOMES core. Mirrors the
structure of [architecture.md](architecture.md): one diagram per
section, each diagram showing **only the public surface and the
relationships** that connect classes. Attributes and methods are
omitted unless they're the single defining signal for that class
(typically the one method on a protocol, or one attribute that
disambiguates two siblings).

Use this file when you need the picture of how the pieces fit. For
the full annotated signatures see [class_diagrams.md](class_diagrams.md).

---

## 1 — ControlVolume at the centre

The CV owns all thermodynamic state and orchestrates a single
timestep via `advance(dt_h, t_h, ...)`. Phases hold mole inventories
and derived scalar properties; an optional `ReactionSystem` carries
the chemistry; an optional list of `PropertyCalculator`s feeds scalar
derived properties (viscosity, density, …); two monitors watch the
solve.

```mermaid
classDiagram
    direction TB

    class ControlVolume {
        +advance(dt_h, t_h, ...) AdvanceResult
    }
    class Phase {
        <<protocol>>
    }
    class PhaseInterface {
        <<protocol>>
    }
    class ExternalBoundary {
        <<protocol>>
    }
    class ReactionSystem
    class PropertyCalculator {
        <<protocol>>
    }
    class AccuracyMonitor
    class ConservationMonitor

    ControlVolume "1" o-- "many" Phase : phases
    ControlVolume "1" o-- "many" PhaseInterface : internal_interfaces
    ControlVolume "1" o-- "many" ExternalBoundary : boundaries
    ControlVolume --> ReactionSystem : reaction_system
    ControlVolume "1" o-- "many" PropertyCalculator : property_calculators
    ReactionSystem ..> AccuracyMonitor : attach_monitor
    ReactionSystem ..> ConservationMonitor : attach_conservation_monitor
```

---

## 2 — Phases (storage)

State-unification put **everything** in one canonical `n_mol` dict:
totals, derived species (`CO2aq`, `H+`, `OH-`, …), strong ions
(`Na+`, `Cl-`, `S_cat`, `S_an`). The speciation engine writes derived
values back via `_refresh_derived(values)`. `phase.properties` is a
separate dict for scalar derived properties (viscosity, density, …)
written by `PropertyCalculator` instances.

```mermaid
classDiagram
    direction TB

    class Phase {
        <<protocol>>
        +n_mol: dict
        +properties: dict
        +T_K: float
        +V_L: float
    }
    class GasPhase
    class LiquidPhase
    class SolidPhase

    Phase <|.. GasPhase
    Phase <|.. LiquidPhase
    Phase <|.. SolidPhase
```

---

## 3 — Reactions & chemistry

Three independent declaration classes (no shared base) carry the
chemistry. `ReactionSystem` is the **single attach point** on a CV;
it pre-buckets its contents by `isinstance` and lazily builds the
`BisectionChemicalEquilibriumEngine` from its equilibrium reactions. The engine writes
derived species directly to `phase.n_mol`. Scalar derived properties
(viscosity, density, …) live on a separate, narrower
`PropertyCalculator` protocol.

```mermaid
classDiagram
    direction TB

    class ReactionModel {
        <<protocol>>
        +compute_rates(env) dict
    }
    class KineticReaction
    class EquilibriumReaction
    class BlackBoxReactionModel
    class ReactionSystem
    class BisectionChemicalEquilibriumEngine
    class PropertyCalculator {
        <<protocol>>
        +key: str
        +compute(phase, T_K, P_atm) float
    }
    class ReactionBuilder

    ReactionModel <|.. KineticReaction
    ReactionModel <|.. BlackBoxReactionModel
    ReactionModel <|.. ReactionSystem
    ReactionSystem "1" *-- "many" KineticReaction
    ReactionSystem "1" *-- "many" EquilibriumReaction
    ReactionSystem "1" *-- "many" BlackBoxReactionModel
    ReactionSystem --> BisectionChemicalEquilibriumEngine : engine (lazy)
    BisectionChemicalEquilibriumEngine ..> EquilibriumReaction : from_reactions
    ReactionBuilder ..> KineticReaction : builds
```

Note: `EquilibriumReaction` is deliberately **not** a `ReactionModel`
— it's algebraic, not rate-producing.

---

## 4 — Time-stepping

`cv.advance(...)` is the default sequential path. For stiffer or
larger-`dt_h` work, pass a `StepSolver` via `solver=` to replace the
body with a snapshot Euler or an adaptive `solve_ivp`-based pass.

```mermaid
classDiagram
    direction LR

    class StepSolver {
        <<protocol>>
        +solve_step(cv, dt_h, t_h, ...) AdvanceResult
    }
    class SimultaneousEulerSolver
    class SimultaneousAdaptiveSolver
    class AdvanceResult

    StepSolver <|.. SimultaneousEulerSolver
    StepSolver <|.. SimultaneousAdaptiveSolver
    StepSolver ..> AdvanceResult : produces
```

---

## 5 — Multi-CV topology and orchestration

`Simulation` is the single orchestration class — it owns the
dict of CVs, the inter-CV link list, controllers, profiles, a
step solver (or per-CV dict), and a recorder. `Simulation.run`
drives the per-step loop. A shared `RunContext` flips
`is_running` on entry and exit so lockable mutators across the
system gate themselves consistently. (The legacy `MultiCVSystem`
was deleted in the `simulation-class` phase, 2026-05-27.)

```mermaid
classDiagram
    direction TB

    class Simulation {
        +run(tau_h, n_steps) BatchResult
    }
    class RunContext {
        +is_running: bool
        +label: str
    }
    class ControlVolume
    class CVLink {
        <<protocol>>
    }
    class AdvectiveLink
    class DiffusiveLink
    class KineticGasLiquidLink

    Simulation "1" o-- "many" ControlVolume : cvs
    Simulation "1" o-- "many" CVLink : links
    Simulation "1" *-- "1" RunContext : _context
    CVLink <|.. AdvectiveLink
    CVLink <|.. DiffusiveLink
    CVLink <|.. KineticGasLiquidLink
```

`KineticGasLiquidLink` is dual-typed — it's both a `CVLink` (carries
material between CVs) and a `PhaseInterface` (acts as an internal
gas-liquid interface inside a single CV).

---

## 6 — External boundaries

State-dependent fluxes across the system edge — feed schedules,
vents, membranes, controller-driven dosing. Boundaries are a runtime
protocol; an `apply_boundary(boundary, cv, dt_h)` helper invokes
one and returns the diagnostic record.

```mermaid
classDiagram
    direction TB

    class ExternalBoundary {
        <<protocol>>
        +compute_flux(cv, dt_h, ...) dict
    }
    class GasFeed
    class PressureReliefVent

    ExternalBoundary <|.. GasFeed
    ExternalBoundary <|.. PressureReliefVent
```

---

## 7 — Monitoring

Two per-CV monitors auto-attached on the `ReactionSystem` surface.
`AccuracyMonitor` watches the speciation engine for solver-quality
signals; `ConservationMonitor` watches element + charge balance
across phases per `cv.advance` step. Both emit a `UserWarning`
subclass and share a throttle plumbing
(`once` / `first_N` / `always` / `silent`).

```mermaid
classDiagram
    direction LR

    class AccuracyMonitor
    class ConservationMonitor
    class AccuracyWarning
    class ConservationWarning

    AccuracyMonitor ..> AccuracyWarning : emits
    ConservationMonitor ..> ConservationWarning : emits
```

---

## 8 — Control system and profiles

Pattern 1 collapse (`simulation-class` 2026-05-27): each
`Controller` exposes a single `compute(state, dt_h) ->
ControlAction` method consuming a `CVSnapshot` (or
`SimulationSnapshot`). The legacy `Commands` dataclass and
`Actuator` protocol are deleted; internal state lives on the
controller. `Profile`s are open-loop time-varying mutators with
the same shape — `apply(t_h, sim) -> ProfileRecord`. The
`Simulation` orchestrator invokes both each step and routes
their outputs through Pattern B unchecked-setter dispatch.

```mermaid
classDiagram
    direction TB

    class Controller {
        <<protocol>>
        +compute(state, dt_h) ControlAction
    }
    class Profile {
        <<protocol>>
        +apply(t_h, sim) ProfileRecord
    }
    class PHController
    class DOAgitationController
    class DOCascadeController
    class PressureReliefController
    class InstantPressureReliefController
    class SmoothPressureReliefController
    class TemperatureRamp
    class VVMSchedule
    class SetpointTrajectory
    class ControlAction
    class ProfileRecord
    class CVSnapshot

    Controller <|.. PHController
    Controller <|.. DOAgitationController
    Controller <|.. DOCascadeController
    Controller <|.. PressureReliefController
    Controller <|.. InstantPressureReliefController
    Controller <|.. SmoothPressureReliefController
    Profile <|.. TemperatureRamp
    Profile <|.. VVMSchedule
    Profile <|.. SetpointTrajectory
    Controller ..> CVSnapshot : reads
    Controller ..> ControlAction : produces
    Profile ..> ProfileRecord : produces
```

`Simulation._step` calls `build_simulation_snapshot(...)` after
the CV advance to assemble a `SimulationSnapshot` (or hands the
single `CVSnapshot` directly for single-CV simulations); each
controller's `compute(state, dt_h)` returns a `ControlAction`,
which the orchestrator applies via `cv.apply_external_flux`
(`flux_applied`) and the `_resolve_param_path` semantic dispatch
to Pattern B unchecked setters (`params_changed`).

---

## Relationship key

- `<|--` solid arrow, hollow head — class inheritance.
- `<|..` dashed arrow, hollow head — protocol implementation.
- `*--` solid line, filled diamond — composition (lifecycle ownership).
- `o--` solid line, hollow diamond — aggregation (uses, doesn't own).
- `-->` solid arrow — association (holds reference).
- `..>` dashed arrow — dependency (uses transiently / produces).

---

## See also

- [architecture.md](architecture.md) — narrative description of the
  same picture, with the `advance()` sub-step ordering, directory
  structure, and orchestrator + control-system layers.
- [class_diagrams.md](class_diagrams.md) — full annotated diagrams
  with every public attribute and method.
- [solvers.md](solvers.md) — semantics of each `StepSolver`.
