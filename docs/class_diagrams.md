# Class Diagrams

Four interface-level UML class diagrams covering the post-Phase-7 core
of PyOMES. Each layer is a separate diagram so the picture stays
readable; classes that span layers (e.g. `KineticGasLiquidLink`,
`ExternalBoundary`) appear in whichever layer they primarily belong
and are referenced by name in the others.

The four layers reflect the architectural framing in
[CONTAINER_LAYERING.md](phases-upcoming/CONTAINER_LAYERING.md):

1. **Topology and containers** — how state is decomposed in space.
2. **Integration / time-stepping** — how time advances.
3. **Properties and reactions** — chemistry-flavoured plug-ins.
4. **External boundaries** — fluxes across the system edge.

Each diagram is **interface-level**: class names, key public
attributes, and key public methods, plus the relationships between
them. Private helpers and one-off utility methods are omitted to
keep the diagrams readable.

The diagrams render natively in GitHub, VS Code's markdown preview,
and most modern markdown viewers via the Mermaid `classDiagram`
syntax. No external tools required.

---

## 1 — Topology and containers

How state is decomposed in space, and what holds it together.
`Phase` is the leaf state; `ControlVolume` aggregates phases,
internal interfaces, external `boundaries`, the
`reaction_system` (chemistry), and the optional list of
`property_calculators` (scalar derived properties). A fermenter
is a plain `ControlVolume` with `{"gas", "liquid"}` phases plus
a `KineticGasLiquidLink` as an internal interface —
`GasLiquidVolume` was deleted in Phase 7. `Simulation` (shipped
in `simulation-class` 2026-05-27) is the single orchestration
class — it owns the dict of CVs, the inter-CV link list,
controllers, profiles, and dispatches a per-CV `StepSolver`
through its `solver=` parameter. The legacy `MultiCVSystem` was
deleted at C14 of that phase.

State-unification (2026-05-21) collapsed the speciation-side
type machinery: `cv.advance(dt_h, t_h)` replaces the legacy
`(dt_h, chem_env)` signature; derived species (`H+`, `OH-`,
`CO2aq`, `HCO3-`, `NH3`, `NH4+`, strong ions) live directly in
`phase.n_mol`; pH derives from `n_mol["H+"]`. The
`KineticGasLiquidLink` reads alphas inline from `phase.n_mol`
ratios — no `PropertyResult.alphas` channel.

```mermaid
classDiagram
    direction TB

    class Phase {
        <<protocol>>
        +T_K: float
        +V_L: float
        +n_mol: dict
        +total_mol() dict
        +apply_flux(flux, dt_h)
        +snapshot() Phase
    }
    class GasPhase {
        +n_total: float
        +P_atm: float
        +y: dict
        +p_atm: dict
        +P_water_sat_atm: float
        +P_total_wet_atm: float
    }
    class LiquidPhase {
        +concentrations_mol_L: dict
        +concentrations_g_L: dict
        +pH: float
        +ionic_strength: float
        +properties: dict
        +_refresh_derived(values)
    }
    class SolidPhase {
        +mass_g: dict
    }
    Phase <|.. GasPhase
    Phase <|.. LiquidPhase
    Phase <|.. SolidPhase

    class PhaseInterface {
        <<protocol>>
        +phase_a_key: str
        +phase_b_key: str
        +compute_flux(state_a, state_b, dt_h) dict
    }
    class HenryEquilibriumInterface {
        +henry: dict
        +equilibrium_species: set
        +compute_flux(...) dict
    }
    class KineticGasLiquidLink {
        +gas_cv_key: str
        +gas_phase_key: str
        +liquid_cv_key: str
        +liquid_phase_key: str
        +henry: dict
        +kLa: dict
        +equilibrium_species: set
        +speciation_keys: dict
        +molecular_driving_force: set
        +henry_params: dict
        +label: str
        +derive_speciation_keys(reaction_system)
        +set_kLa(species, value)
        +set_kLa_with_co2_ratio(kLa_O2, co2_ratio)
        +set_transfer_mode(species, mode, kLa)
        +validate() list
        +compute_flux(state_a, state_b, dt_h, instantaneous) dict
        +compute_flow(cvs, dt_h, instantaneous) dict
    }
    PhaseInterface <|.. HenryEquilibriumInterface
    PhaseInterface <|.. KineticGasLiquidLink

    class CVLink {
        <<protocol>>
        +source_cv_key: str
        +source_phase_key: str
        +sink_cv_key: str
        +sink_phase_key: str
        +label: str
        +compute_flow(cvs, dt_h) dict
    }
    class AdvectiveLink {
        +Q_L_per_h: float
        +species_filter: list
    }
    class DiffusiveLink {
        +kLa: dict
        +V_eff_L: float
    }
    CVLink <|.. AdvectiveLink
    CVLink <|.. DiffusiveLink
    CVLink <|.. KineticGasLiquidLink

    class ControlVolume {
        +phases: dict
        +internal_interfaces: list
        +boundaries: list
        +reaction_system: ReactionSystem
        +property_calculators: list
        +label: str
        +phase_keys: list
        +advance(dt_h, t_h, external_source_terms, solver) AdvanceResult
        +step_internal_transfer(dt_h) TransferDiagnostics
        +compute_reaction_rates(t_h) dict
        +apply_external_flux(phase_key, flux, dt_h)
        +total_mol() dict
        +snapshot() ControlVolume
    }
    ControlVolume "1" o-- "many" Phase : phases
    ControlVolume "1" o-- "many" PhaseInterface : internal_interfaces
    ControlVolume "1" o-- "many" ExternalBoundary : boundaries

    class Simulation {
        +cvs: dict
        +links: list
        +controllers: _LockableList
        +profiles: _LockableList
        +solver: Any
        +recorder: Any
        +cv_keys: list
        +run(tau_h, n_steps) BatchResult
        +total_mol() dict
        +snapshot() Simulation
    }
    Simulation "1" o-- "many" ControlVolume : cvs
    Simulation "1" o-- "many" CVLink : links
    Simulation "1" *-- "1" RunContext : _context

    class RunContext {
        +is_running: bool
        +label: str
    }

    class ExternalBoundary {
        <<protocol>>
    }
    class StepSolver {
        <<protocol>>
    }
```

### Relationship key

- `<|--` solid arrow with hollow head — class inheritance.
- `<|..` dashed arrow with hollow head — protocol implementation.
- `*--` solid line with filled diamond — composition (lifecycle ownership).
- `o--` solid line with hollow diamond — aggregation (uses, doesn't own lifecycle).
- `-->` solid arrow — association (holds reference).

---

## 2 — Integration / time-stepping

How a single timestep advances. `StepSolver` strategies operate on a
`ControlVolume` directly (snapshot or adaptive) and are dispatched
via `cv.advance(solver=...)`; the default sequential body of
`ControlVolume.advance()` integrates the reaction sub-step with a
single forward Euler evaluation (no pluggable inner integrator).
The diagnostic dataclasses on the right are what each path returns
to the orchestrator — Phase 6 unified the old
`GasLiquidAdvanceResult` into the single `AdvanceResult`, which
carries `transfer_record` and `boundary_records`.

See [docs/solvers.md](solvers.md) for the semantics of each
`StepSolver`.

```mermaid
classDiagram
    direction LR

    class StepSolver {
        <<protocol>>
        +solve_step(cv, dt_h, t_h, external_source_terms) AdvanceResult
    }
    class SimultaneousEulerSolver {
        +solve_step(...) AdvanceResult
    }
    class SimultaneousAdaptiveSolver {
        +rtol: float
        +atol: float
        +max_step: float
        +method: str
        +freeze_speciation: bool
        +solve_step(...) AdvanceResult
    }
    StepSolver <|.. SimultaneousEulerSolver
    StepSolver <|.. SimultaneousAdaptiveSolver

    class TransferDiagnostics {
        +species_before: dict
        +species_after: dict
        +residual: dict
        +fluxes: list
        +total_mol_before: float
        +total_mol_after: float
        +total_residual: float
        +max_abs_residual: float
        +is_conserved(tol) bool
        +summary() str
    }
    class AdvanceResult {
        +transfer: TransferDiagnostics
        +reaction_sources: dict
        +transfer_record: LinkFlowRecord
        +boundary_records: list
    }
    class LinkFlowRecord {
        +link_label: str
        +source: str
        +sink: str
        +flow_mol_per_h: dict
        +total_mol_transferred: float
    }
    class ExternalFluxRecord {
        +boundary_label: str
        +phase_key: str
        +flux_mol_per_h: dict
        +dt_h: float
        +mol_applied: dict
        +total_mol_applied: float
        +summary() str
    }

    AdvanceResult *-- TransferDiagnostics
    AdvanceResult *-- LinkFlowRecord
    AdvanceResult *-- ExternalFluxRecord

    SimultaneousEulerSolver ..> AdvanceResult : produces
    SimultaneousAdaptiveSolver ..> AdvanceResult : produces

    class ControlVolume {
        <<from layer 1>>
    }
    class Simulation {
        <<from layer 1>>
    }
    ControlVolume ..> StepSolver : dispatched via solver=
    ControlVolume ..> AdvanceResult : produces
    Simulation ..> BatchResult : produces (via Recorder.finalize)
```

---

## 3 — Properties and reactions

Two narrow chemistry-flavoured plug-ins on the CV. **Property
calculators** compute one scalar derived property each (viscosity,
density, …) and write it to `phase.properties[key]`. **The reaction
system** holds every reaction declaration on the CV and pre-buckets
them by isinstance at construction so each one routes to the right
consumer. Speciation (acid-base equilibria) is *not* a property —
the engine writes derived molecular species straight back into
`phase.n_mol`, sharing one canonical store with the kinetic
sub-step.

`KineticGasLiquidLink` reads the molecular fraction it needs for the
speciation-corrected Henry constant directly from
`liquid.n_mol[mol_key]` (resolved via `_CANONICAL_NAMES` against the
gas key), then divides by the canonical total — the speciation
engine has already populated both for the current step.

> **State as of `state-unification` (2026-05-22, shipped).** This
> phase collapsed the speciation-shaped type machinery left over
> from earlier rounds. `PropertySolver` / `PropertyResult` /
> `SpeciationPropertySolver` are deleted; their slot is taken by
> the narrow [`PropertyCalculator`](../PyOMES/core/property_calculator.py)
> protocol (one `key: str`, one `compute(phase, T_K, P_atm) -> float`,
> no `chem_env`, no result struct). The unified `Reaction` class
> with `kind=` branching is replaced by three independent
> declaration classes:
> [`KineticReaction`](../PyOMES/reactions/kinetic.py),
> [`EquilibriumReaction`](../PyOMES/reactions/equilibrium.py), and
> [`BlackBoxReactionModel`](../PyOMES/reactions/blackbox.py) — no
> shared base, shared validation lives as free functions in
> `_shared.py`. `ReactionSet` is replaced by
> [`ReactionSystem`](../PyOMES/reactions/reaction_system.py), which
> pre-buckets reactions into `_kinetic_reactions`,
> `_single_phase_equilibria`, `_cross_phase_equilibria`, and
> `_blackbox_models` at `__init__`. `cv.reaction_system` is the
> single attach point (the old `cv.reaction_model` is gone);
> `cv.property_solvers` becomes `cv.property_calculators`. The
> system lazily builds its `BisectionChemicalEquilibriumEngine` on first access, and
> the engine writes derived species (`CO2aq`, `NH3`, `HAc`, …)
> straight to `phase.n_mol` via the privileged
> `phase._refresh_derived(values)` hook — eliminating the
> `PropertyResult.alphas` channel and the dict-of-dicts plumbing
> that fed it. `chem_env` is deleted system-wide;
> `cv.advance(dt_h, t_h)` is the new entry signature.
>
> Monitoring is attached on the same surface: an
> [`AccuracyMonitor`](../PyOMES/monitoring/accuracy.py) hooks engine
> internals via `cv.reaction_system.attach_monitor(monitor)`, and a
> [`ConservationMonitor`](../PyOMES/monitoring/conservation.py) (new
> in C6) checks element + charge balance per step via
> `cv.reaction_system.attach_conservation_monitor(monitor)`. The CV
> auto-attaches a default `ConservationMonitor` in `__init__` and
> populates its species registry from the reaction stoichiometries.
>
> **Still deferred:** retirement of the engine `level` parameter
> (currently legacy default `level=2`, hidden behind
> `configure_engine`) — folds into the future
> `ChemistryDatabase` track. See
> [phases-upcoming/SPECIATION_LEVEL_RETIREMENT.md](phases-upcoming/SPECIATION_LEVEL_RETIREMENT.md).

```mermaid
classDiagram
    direction TB

    class PropertyCalculator {
        <<protocol>>
        +key: str
        +compute(phase, T_K, P_atm) float
    }

    class BisectionChemicalEquilibriumEngine {
        +use_warmstart: bool
        +n_solve_calls: int
        +from_reactions(equilibrium_rxns, ...)$ BisectionChemicalEquilibriumEngine
        +solve(phases?, **kwargs) dict
    }

    class ReactionModel {
        <<protocol>>
        +compute_rates(env) dict
    }
    class KineticReaction {
        +stoichiometry: list
        +rate_fn: callable
        +balance_elements: tuple
        +label: str
        +species_ids: list
        +phases: list
        +is_cross_phase: bool
        +compute_rates(env) dict
    }
    class EquilibriumReaction {
        +stoichiometry: list
        +log_K: float?
        +dH_J_per_mol: float?
        +T_ref_K: float
        +total_id: str?
        +balance_elements: tuple
        +label: str
        +species_ids: list
        +phases: list
        +is_cross_phase: bool
    }
    class BlackBoxReactionModel {
        +external_model
        +flux_mapping: dict
        +balance_elements: tuple
        +balance_atol: float
        +on_imbalance: str
        +imbalance_history: list
        +compute_rates(env) dict
    }
    class ReactionSystem {
        +reactions: list
        +label: str
        +kinetic_reactions: list
        +single_phase_equilibria: list
        +cross_phase_equilibria: list
        +blackbox_models: list
        +species_ids: list
        +phases: list
        +engine: BisectionChemicalEquilibriumEngine
        +configure_engine(use_activity, activity_model, level)
        +attach_engine(engine)
        +attach_monitor(monitor)
        +attach_conservation_monitor(monitor)
        +compute_rates(env) dict
    }
    class FluxEntry {
        +species_id: str
        +phase: str
        +atoms: dict
        +sign: float
    }
    class ReactionBuilder {
        +aerobic_growth(...)$ KineticReaction
        +from_coefficients(...)$ KineticReaction
    }
    class ReactionEnvironment {
        +T_K: float
        +V_L: float
        +pH: float?
        +has_pH: bool
        +concentrations: dict
        +properties: dict
        +t_h: float
        +S(species_id, default) float
        +X(organism_id, default) float
        +prop(key, default) float
    }
    class StoichiometryEntry {
        +species: Species
        +phase: str
        +coefficient: float
    }
    class Species {
        <<frozen>>
        +id: str
        +atoms: Mapping
        +charge: int
        +MW: float
    }
    ReactionModel <|.. KineticReaction
    ReactionModel <|.. ReactionSystem
    ReactionModel <|.. BlackBoxReactionModel
    ReactionSystem "1" *-- "many" KineticReaction : kinetic
    ReactionSystem "1" *-- "many" EquilibriumReaction : equilibria
    ReactionSystem "1" *-- "many" BlackBoxReactionModel : blackbox
    ReactionSystem --> BisectionChemicalEquilibriumEngine : engine (lazy)
    BisectionChemicalEquilibriumEngine ..> EquilibriumReaction : from_reactions
    KineticReaction "1" *-- "many" StoichiometryEntry
    EquilibriumReaction "1" *-- "many" StoichiometryEntry
    StoichiometryEntry --> Species : species
    BlackBoxReactionModel "1" *-- "many" FluxEntry : flux_mapping
    ReactionBuilder ..> KineticReaction : builds
    ReactionModel ..> ReactionEnvironment : reads

    class ControlVolume {
        <<from layer 1>>
    }
    ControlVolume ..> PropertyCalculator : property_calculators
    ControlVolume ..> ReactionSystem : reaction_system
    ControlVolume ..> ReactionEnvironment : builds for compute_rates

    class LiquidPhase {
        <<from layer 1>>
    }
    BisectionChemicalEquilibriumEngine ..> LiquidPhase : writes derived species via _refresh_derived

    class KineticGasLiquidLink {
        <<from layer 1>>
    }
    KineticGasLiquidLink ..> LiquidPhase : reads molecular n_mol via _CANONICAL_NAMES
    KineticGasLiquidLink ..> ReactionSystem : derives speciation_keys from cross_phase_equilibria
```

---

## 4 — External boundaries

State-dependent fluxes across the system edge. `ExternalBoundary`
is a runtime-checkable protocol, not a base class — implementations
just need `phase_key`, `label`, and
`compute_flux(cv, dt_h, instantaneous)`.  The boundaries layer is
where feed schedules, vents, membranes, and controller-driven dosing
get plugged in.  `apply_boundary` is a free-function helper that
invokes a boundary against a CV and returns the diagnostic record.

```mermaid
classDiagram
    direction TB

    class ExternalBoundary {
        <<protocol>>
        +phase_key: str
        +label: str
        +compute_flux(cv, dt_h, instantaneous) dict
    }
    class GasFeed {
        +vvm_min: float
        +y: dict
        +P_inlet_atm: float
        +phase_key: str
        +liquid_phase_key: str
        +label: str
        +V_liq_override: float
        +T_override_K: float
    }
    class PressureReliefVent {
        +P_set_atm: float
        +mode: str
        +k_vent_per_h: float
        +smooth_width_atm: float
        +phase_key: str
        +label: str
    }
    class MembraneGasBoundary {
        +permeability: dict
        +area_m2: float
        +external_atmosphere: dict
        +phase_key: str
        +label: str
    }
    class LiquidFeed {
        +Q_L_per_h: float
        +feed_conc_mol_L: dict
        +phase_key: str
        +label: str
    }
    class LiquidDrain {
        +Q_L_per_h: float
        +species_filter: set
        +phase_key: str
        +label: str
    }
    class ProportionalGasOutlet {
        +k_p_L_per_h_per_atm: float
        +P_atm: float
        +include_water_vapour: bool
        +last_water_loss_mol_per_h: float
        +phase_key: str
        +label: str
    }
    ExternalBoundary <|.. GasFeed
    ExternalBoundary <|.. PressureReliefVent
    ExternalBoundary <|.. MembraneGasBoundary
    ExternalBoundary <|.. LiquidFeed
    ExternalBoundary <|.. LiquidDrain
    ExternalBoundary <|.. ProportionalGasOutlet
    %% WaterVapourBoundary/VentWaterLoss retired (LAYER1_GAP_CLOSURE CP4) —
    %% replaced by transfer_models=EquilibriumTransferModel(RaoultEquilibrium()),
    %% which properly conserves mass between the liquid and gas H2O pools.

    class ExternalFluxRecord {
        +boundary_label: str
        +phase_key: str
        +flux_mol_per_h: dict
        +dt_h: float
        +mol_applied: dict
    }
    class apply_boundary {
        <<function>>
        +apply_boundary(boundary, cv, dt_h) ExternalFluxRecord
    }
    apply_boundary ..> ExternalBoundary : invokes
    apply_boundary ..> ExternalFluxRecord : returns

    class ControlVolume {
        <<from layer 1>>
    }
    class SimultaneousEulerSolver {
        <<from layer 2>>
    }
    ControlVolume "1" --> "many" ExternalBoundary : boundaries
    SimultaneousEulerSolver ..> ExternalBoundary : evaluates from snapshot
```

---

## 5 — Monitoring

Two per-CV monitors attached on the `ReactionSystem` surface watch
the speciation + integration path for drift. Both share a throttle
plumbing (`once` / `first_N` / `always` / `silent`) configured via
`PyOMES.config.warnings`, and both expose a module-level summary
counter for end-of-run reporting.

- [`AccuracyMonitor`](../PyOMES/monitoring/accuracy.py) — emitted from
  inside the `BisectionChemicalEquilibriumEngine` solve loop when the Newton residual,
  charge balance, or iteration count exceeds the configured
  thresholds. Attached via
  `cv.reaction_system.attach_monitor(monitor)`; the system
  propagates the monitor to the engine on first build.
- [`ConservationMonitor`](../PyOMES/monitoring/conservation.py)
  (state-unification C6) — invoked once per
  `ControlVolume.advance` step. Tracks element totals and net
  charge across all phases; emits a `ConservationWarning` when
  drift exceeds the per-step or cumulative threshold. The species
  registry (`{species_id: Species}`) is populated automatically by
  `ControlVolume.__init__` from the reaction stoichiometries.

```mermaid
classDiagram
    direction TB

    class AccuracyWarning
    class AccuracyMonitor {
        +step_count: int
        +reset()
        +set_species_registry(registry)
        +check_solve(...)
    }
    class ConservationWarning
    class ConservationMonitor {
        +step_count: int
        +reset()
        +set_species_registry(registry)
        +check_step(phases) None
    }
    AccuracyMonitor ..> AccuracyWarning : emits
    ConservationMonitor ..> ConservationWarning : emits

    class ReactionSystem {
        <<from layer 3>>
    }
    class BisectionChemicalEquilibriumEngine {
        <<from layer 3>>
    }
    class ControlVolume {
        <<from layer 1>>
    }
    ReactionSystem ..> AccuracyMonitor : attach_monitor
    ReactionSystem ..> ConservationMonitor : attach_conservation_monitor
    BisectionChemicalEquilibriumEngine ..> AccuracyMonitor : invokes from solve()
    ControlVolume ..> ConservationMonitor : invokes check_step in advance()
```

---

## How to use these diagrams

- **Reading a new piece of code:** start in layer 1 to see where it
  fits structurally, then drill into 2/3/4 for the specific machinery.
- **Adding a new unit:** layers 1 and 4 are usually what you touch.
  Build a `ControlVolume` with the right phases, internal interfaces,
  and boundaries (the fermenter pattern is `{"gas", "liquid"}` phases
  plus a `KineticGasLiquidLink`).
- **Adding a new integration strategy:** layer 2 — implement the
  `StepSolver` protocol (operating on a `ControlVolume`) and pass it
  via `cv.advance(solver=...)` or
  `Simulation(cvs=..., solver={cv_key: solver}).run(...)`.
- **Adding a new chemistry plug-in:** layer 3 — implement
  `PropertyCalculator` (one scalar value, `phase.properties[key]`)
  or contribute a reaction declaration (`KineticReaction`,
  `EquilibriumReaction`, or `BlackBoxReactionModel`) to a
  `ReactionSystem`.

Diagrams are kept at interface-level deliberately: private helpers
(`_clamp_deltas`, `_build_reaction_environment`,
`_SnapshotCV`, etc.) are implementation detail and would clutter
the picture. Read the source for those; the diagrams are the
**contract** view.

---

## Maintenance

These diagrams are hand-curated, not auto-generated. They will drift
if the code changes and no one updates them. When making a change
that affects a class's public surface (new protocol method, renamed
attribute, new wrapper class, removed method), update the relevant
layer here in the same PR. The diagrams are short — touching them
is cheap.

Out of scope for these diagrams (intentionally):

- The `PyOMES/chemistry/` package (compounds, registry, `EquilibriumSet`,
  `ThermodynamicConfig`).
- The `PyOMES/chemical_equilibrium/` engine internals (`Level1`, `Level2`,
  `Level2_5`, `StrongIonsSolver`).
- The `PyOMES/control/` package (`ControlSystem`, `PIController`,
  actuators).
- The `models/vlmodels/` packages (fermenter unit, ADM1, HPLC).
- The `PyOMES/numerics/` package (advection / dispersion schemes).

If a future phase pulls any of those into the architectural picture,
add a new layer rather than expanding an existing one — the
[CONTAINER_LAYERING.md](phases-upcoming/CONTAINER_LAYERING.md) framing argues for
keeping concerns separable.
