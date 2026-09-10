# Architecture: CV-Centric PyOMES Simulation

This document describes the internal architecture of the PyOMES simulation
module as of v0.12.5.  (The core package, previously distributed as
`fermenter`, was renamed to `PyOMES` in this version; the `vlmodels`
package of concrete units kept its name.)

## Overview

The simulation is structured around a **ControlVolume** (CV) that owns all
thermodynamic state (mole inventories, equilibrium, derived properties). The
fermenter orchestrator is a thin loop that applies external fluxes to the CV,
calls `cv.advance()` each timestep, and handles controller actions (pressure
relief, pH dosing, DO control).

This design was chosen to support multi-zone reactor models, where multiple CVs
represent different zones within a single unit operation (e.g. a sparger zone
with high kLa and a bulk liquid zone with lower kLa), connected by material
transport links.

## Key Abstractions

### ControlVolume (`PyOMES/core/control_volume.py`)

A CV holds:

- **Phases** — `GasPhase` and `LiquidPhase` objects containing mole
  inventories in a single canonical `n_mol` dict (totals + derived
  species both live here, PHREEQC-style).
- **Internal interfaces** — `KineticGasLiquidLink` for kLa-based
  gas-liquid partitioning.
- **Reaction system** (optional) — one `ReactionSystem` holding all
  reaction declarations on the CV, pre-bucketed by kind.
- **Property calculators** (optional) — list of `PropertyCalculator`
  instances each computing one scalar derived property (viscosity,
  density, …) onto `phase.properties[key]`.
- **Monitors** — `AccuracyMonitor` and `ConservationMonitor`
  attached on the reaction system surface (the CV auto-attaches a
  default `ConservationMonitor` in `__init__`).

The `advance(dt_h, t_h, external_source_terms=None, solver=None)`
method sequences:

1. Property calculators (viscosity, density, …) — write to
   `phase.properties`.
2. Speciation solve via `reaction_system.engine.solve(phases=...)` —
   writes derived molecular species (`CO2aq`, `NH3`, `HAc`, `H+`, …)
   back to `phase.n_mol` via the privileged
   `phase._refresh_derived(values)` hook.
3. Apply `external_source_terms` (orchestrator-supplied fluxes).
4. Integrate kinetic + black-box reactions (post-feed state; reads
   pre-step properties via `env.prop(key)`).
5. Internal kinetic transfer (`KineticGasLiquidLink`, kLa, Henry
   constants corrected via molecular fractions resolved from
   `phase.n_mol` against the system's `cross_phase_equilibria`).
6. `ConservationMonitor.check_step(phases)` — element + charge drift
   guard.

The CV returns an `AdvanceResult` containing transfer diagnostics
and reaction source terms. `chem_env` is gone; pH is read from
`liquid.n_mol["H+"]` and ionic strength is recomputable from the
canonical `n_mol` directly.

### Gas-liquid CV pattern

After Phase 7 (`phases-shipped/PHASE7_CHECKLIST.md`),
`GasLiquidVolume` is deleted.  A fermenter is now a plain
`ControlVolume` whose `phases` dict is
`{"gas": GasPhase, "liquid": LiquidPhase}` and whose
`internal_interfaces` list contains a `KineticGasLiquidLink` acting
as a `PhaseInterface`.  Callers construct one via
`FermenterFactory.create_volume(...)` (or `FermenterBuilder().build()`).

Both `SimultaneousEulerSolver` and `SimultaneousAdaptiveSolver` operate directly on
a `ControlVolume` and return the unified `AdvanceResult`, which
carries `transfer_record` and `boundary_records`. Boundaries live
on `cv.boundaries`. Controllers that need the gas-liquid link reach
into `cv.internal_interfaces` directly with an `isinstance` check
against `KineticGasLiquidLink`. pH is read inline from
`liquid.n_mol["H+"]` (the speciation engine refreshes it at the
top of every `advance` step). The `AdvanceResult.properties`
channel is gone — derived species are first-class in `n_mol` and
scalar properties live on `phase.properties`.

See `CV_UPDATE.md`, `SOLVER_PROMOTION.md`,
`SIMULATION_CLASS.md` (shipped 2026-05-27; retired `MultiCVSystem`
in favour of the unified `Simulation` orchestrator), and
`CONTAINER_LAYERING.md` for the broader layering question of
whether the fermenter pattern and `HPLCColumn` should be
unbundled into orthogonal topology / unit / integration axes.

### Phases (`PyOMES/core/phases.py`)

`GasPhase` and `LiquidPhase` hold mole inventories in a single
canonical `n_mol` dict and volume/temperature. State-unification
collapsed totals + derived species into one store (PHREEQC-style):
the speciation engine writes derived molecular species (`H+`,
`OH-`, `CO2aq`, `HCO3-`, `NH3`, `NH4+`, …) back to `n_mol` via the
privileged `_refresh_derived(values)` hook. Strong ions (named
species like `Na+`/`Cl-`, plus unnamed lumps `S_cat`/`S_an` with
`atoms={}`) live alongside them in `n_mol`. Gas-phase derived
quantities (pressure, mole fractions, partial pressures) are
computed from the ideal gas law; liquid-phase concentrations
(mol/L and g/L) from moles and volume.

`phase.properties` is a separate dict for scalar derived properties
(viscosity, density, heat capacity, …), written by
`PropertyCalculator` instances and read by rate laws via
`env.prop(key)`.

Species without a Henry constant in the equilibrium interface
(e.g. biomass, substrates) are stored in the liquid phase but are
not touched by the equilibrium solver. They participate in
inter-zone transport via advective links.

### Property Calculators (`PyOMES/core/property_calculator.py`)

A `PropertyCalculator` is a narrow protocol with a `key: str` and
`compute(phase, T_K, P_atm) -> float`. Each calculator computes
**one scalar derived property** (viscosity, density, heat capacity,
…) and writes it to `phase.properties[key]`. Rate laws read it via
`env.prop("viscosity")`.

Speciation is *not* a `PropertyCalculator`. Acid-base equilibria are
state-completion handled by the `BisectionChemicalEquilibriumEngine` attached to
`cv.reaction_system`; the engine writes derived molecular species
straight to `phase.n_mol`, not to `phase.properties`.

Property calculators are registered via
`cv.property_calculators=[...]` and evaluated automatically at the
top of every `advance()` step (before the speciation solve and
reaction integration).

### Equilibrium pathways

Two parallel equilibrium pathways exist in the codebase:

1. **`BisectionChemicalEquilibriumEngine` + `KineticGasLiquidLink`** — the canonical
   pathway used by every CV-based model (ADM1, BSM2, any
   `Simulation`-orchestrated CV). Acid-base equilibria are declared
   as `EquilibriumReaction` objects on the `ReactionSystem` and
   solved by the engine in step 1 of `cv.advance()`; the engine
   writes derived molecular species back to `phase.n_mol` via the
   privileged `_refresh_derived(values)` hook. Gas-liquid
   partitioning is kLa-driven through `KineticGasLiquidLink`, with
   alpha-correction read inline from `phase.n_mol` against the
   reaction system's `cross_phase_equilibria` bucket.

2. **`ProcessCoupledEquilibrator` + `HenryEquilibriumInterface`**
   (`PyOMES/equilibria/`, `PyOMES/core/gl_equilibrium.py`) — a parallel
   single-shot CO₂/O₂/N₂ partition solver used **only** by the
   legacy `CUFermentationSpeciation` class in
   `models/vlmodels/fermenter/unit.py`. Slated for removal in the
   `CUFERMENTER_SUNSET` phase. Do not use for new work.

## Reaction Framework (`PyOMES/reactions/`)

### Three independent declaration classes

A reaction is declared as one of three independent classes — no
shared base, shared validation lives as free functions in
`_shared.py`:

- `KineticReaction` — stoichiometric template + callable rate law.
  Stoichiometry is validated at construction for elemental balance
  (default `("C", "H", "O")`); `StoichiometryError` carries a
  per-species breakdown.
- `EquilibriumReaction` — stoichiometric template + `log_K` (and
  optional Van 't Hoff parameters `dH_J_per_mol` / `T_ref_K`).
  Routed to the `BisectionChemicalEquilibriumEngine` as an algebraic constraint
  (single-phase) or to the `KineticGasLiquidLink` as a partition
  declaration (cross-phase, `log_K` optional). Has no
  `compute_rates`.
- `BlackBoxReactionModel` — adapter for opaque external simulators
  (FBA, genome-scale, proprietary). Stoichiometry is not
  inspectable, so elemental balance is checked at runtime against
  the configured `flux_mapping` with a policy (warn, raise, log,
  ignore).

### ReactionSystem

`ReactionSystem` is the **single attach point** for every reaction
on a CV (`cv.reaction_system`). At construction it pre-buckets its
input into `_kinetic_reactions`, `_single_phase_equilibria`,
`_cross_phase_equilibria`, and `_blackbox_models` by `isinstance` —
each consumer pulls from its own bucket, and mixed-kind systems
are always safe to integrate (equilibria are simply not in the
kinetic iteration path).

It also owns the lazy speciation engine via the `engine` property:
on first access, `BisectionChemicalEquilibriumEngine.from_reactions(...)` is built
from the bucketed equilibria with the settings pinned by
`configure_engine(use_activity=, activity_model=, level=)`. Tests
needing a pre-built engine inject it via `attach_engine(engine)`.
Monitors attach on the same surface: `attach_monitor(monitor)` for
the `AccuracyMonitor` (propagated to the engine on build) and
`attach_conservation_monitor(monitor)` for the `ConservationMonitor`.

### ReactionBuilder

`ReactionBuilder.aerobic_growth()` derives O₂/CO₂/H₂O (and N-source)
stoichiometric coefficients from substrate formula, biomass formula,
and yield and returns a validated `KineticReaction`.
`ReactionBuilder.from_coefficients()` accepts explicit
stoichiometry for custom reactions.

### ReactionModel protocol

`KineticReaction`, `ReactionSystem`, and `BlackBoxReactionModel`
all implement `compute_rates(env) -> {phase: {species: mol/h}}`.
`EquilibriumReaction` deliberately does not — it's algebraic, not
a rate producer.

## Orchestration (`PyOMES/core/simulation.py`)

### Simulation

`Simulation` is the single orchestration class for all CV-based
models (shipped in the `simulation-class` phase, 2026-05-27). It
holds a named dict of ControlVolumes, an optional list of inter-CV
`CVLink` objects, controllers, profiles, a step solver (or per-CV
dict), and a recorder. The constructor wires every CV (and its
nested phases, links, lockable lists) to a shared `RunContext`
that the lifecycle gate consults during `.run()`.

`Simulation.run(tau_h, n_steps)` is the user-facing entry point.
The per-step body (`_step`) is structured as:

1. **Profiles** (open-loop time-varying mutators) fire *before*
   integration so the advance sees updated state. Each profile's
   `apply(t_h, sim)` returns a `ProfileRecord` for audit.
2. **Inter-CV link flows** are computed and applied symmetrically
   (source loses, sink gains), yielding one `LinkFlowRecord` per
   link.
3. **Each CV advances** via `cv.advance(dt_h, t_h, solver=...)`
   with its dispatched solver. The CV refreshes property
   calculators, runs the speciation solve into `phase.n_mol`,
   applies source terms, integrates kinetic reactions one Euler
   sub-step, then runs internal kinetic transfer.
4. **A `SimulationSnapshot`** is built from the post-advance
   state.
5. **Controllers** (respecting sample-period gating; zero-order
   hold between samples) consume the snapshot and return a
   `ControlAction` describing the mutation they want applied.
6. **The orchestrator applies each action**: `flux_applied` via
   `cv.apply_external_flux`; `params_changed` via the Pattern B
   semantic-path dispatch to unchecked setters
   (`link._set_kLa_unchecked`, `GasFeed` attribute writes, etc.).
7. The recorder receives per-step records and finalises into a
   `BatchResult` (per-CV nested arrays).

The legacy `_calc_ODE_with_headspace()` method on the
`CUFermentationSpeciation` BioSTEAM wrapper still exists for
backward compatibility; it is **out of scope** for the new
orchestrator and follows its own time-loop pattern.

### Lifecycle gating (RunContext)

A small `RunContext` dataclass holds the run-bound state shared
across all objects owned by a single `Simulation`. The Simulation
flips `self._context.is_running = True` at run-entry and clears
it in a `try/finally`. Every lockable owned object — CVs,
Phases, `KineticGasLiquidLink`, `_LockableList` wrappers for
`cv.boundaries` / `cv.property_calculators` /
`sim.controllers` / `sim.profiles` — holds a reference to that
one context and consults `self._context.is_running` from a
`raise_if_running` guard before mutating.

For controllers and profiles that need to mutate state during a
run (the legitimate case), the framework provides
**Pattern B unchecked siblings** to every gated mutator
(e.g. `_set_kLa_unchecked`, `_set_T_K_unchecked`,
`_set_reaction_system_unchecked`). User code calls the public
gated method; the orchestrator's action-apply path uses the
unchecked siblings.

### Links (`links.py`)

`AdvectiveLink` models bulk liquid (or gas) circulation. All species in the
source phase are carried at their local concentration proportional to the
volumetric flow rate: `flux_i = C_i × Q`.

`DiffusiveLink` models concentration-driven transfer of specific species:
`flux_i = kLa_i × V_eff × (C_source - C_sink)`.

Both link types include safety clamping to prevent transferring more moles than
exist in the source within one timestep.

### Example: Two-Zone Fermenter

```python
from PyOMES.core import (
    ControlVolume, GasPhase, LiquidPhase,
    Simulation, AdvectiveLink,
)

# Sparger zone: small, high kLa, has the gas phase
cv_sparger = ControlVolume(
    phases={"gas": GasPhase(...), "liquid": LiquidPhase(...)},
    internal_interfaces=[link_sparger],
    reaction_system=my_reaction_system,
)

# Bulk zone: larger, lower kLa, liquid only
cv_bulk = ControlVolume(
    phases={"liquid": LiquidPhase(...)},
    reaction_system=my_reaction_system,
)

# Bidirectional circulation
sim = Simulation(
    cvs={"sparger": cv_sparger, "bulk": cv_bulk},
    links=[
        AdvectiveLink("sparger", "liquid", "bulk", "liquid", Q_L_per_h=500),
        AdvectiveLink("bulk", "liquid", "sparger", "liquid", Q_L_per_h=500),
    ],
)
result = sim.run(tau_h=24.0, n_steps=2400)
# result.liquid_mol["sparger"]["O2"] etc.
```

## Control System (`PyOMES/control/`)

Controllers in the new framework consume a typed `CVSnapshot`
(single-CV) or `SimulationSnapshot` (multi-CV) and return a
structured `ControlAction` from a single
`compute(state, dt_h) -> ControlAction` method. Internal state
(setpoint, integral, sampling-hold buffer) lives on the controller
instance.

1. `Simulation.run` calls `ctrl.reset()` at run-entry on every
   controller (decision 10: always destructive).
2. Per step (post-advance), the orchestrator builds the snapshot,
   filters by sample-period (zero-order hold between samples),
   and invokes `ctrl.compute(snapshot, dt_h)`.
3. The returned `ControlAction` carries `flux_applied`
   (`{phase: {species: mol/h}}`), `params_changed`
   (`{semantic_path: value}`), `vented_mol`, and `dosed_mol`.
4. The orchestrator applies the action: `flux_applied` via
   `cv.apply_external_flux`; `params_changed` via the C9
   semantic-path resolver to the Pattern B unchecked setters.

CV-native concrete controllers ship in `PyOMES/control/cv_loops.py`:
`PHController`, `DOAgitationController` (+`DOController` alias),
`DOCascadeController`, `InstantPressureReliefController`,
`SmoothPressureReliefController`, `PressureReliefController`.
Profiles (open-loop time-varying mutators) ship in
`PyOMES/control/cv_profiles.py`: `TemperatureRamp`, `VVMSchedule`,
`SetpointTrajectory`.

The legacy `PyOMES/control/system.py:ControlSystem` and
`PyOMES/control/loops.py` (FermenterState-based controllers) remain
alive for the still-out-of-scope `CUFermentationSpeciation`
BioSTEAM wrapper; their deletion is a follow-up phase.

## Directory layout

```
PyOMES/
  __init__.py
  units.py                       # Shared constants and unit conversions
  config.py                      # PyOMES.config — WarningConfig, env-var presets
  stream_adapter.py              # FeedState, FermenterResult (BioSTEAM-flavoured)
  chemistry/                     # Compound registry, recipes, Species declarations
    compounds.py, recipe.py, registry.py, species.py, species_check.py,
    common_species.py, equilibria.py, thermo_params.py, types.py, chem_recipe.py
  core/                          # Framework core: phases, CV, orchestration
    phases.py                    # GasPhase, LiquidPhase, SolidPhase
    interfaces.py                # PhaseInterface, TransferDiagnostics, AdvanceResult
    control_volume.py            # ControlVolume.advance(dt_h, t_h, …)
    gas_liquid_link.py           # KineticGasLiquidLink (kLa + alpha)
    gl_equilibrium.py            # HenryEquilibriumInterface — legacy (CUFermenter island)
    property_calculator.py       # PropertyCalculator protocol
    solvers.py                   # SimultaneousEulerSolver, SimultaneousAdaptiveSolver, StepSolver
    links.py                     # CVLink, AdvectiveLink, DiffusiveLink, LinkFlowRecord
    boundaries.py                # GasFeed, Vent, MembraneGasBoundary, LiquidFeed, …
    simulation.py                # Simulation orchestrator + run loop + RunContext
    snapshot.py                  # CVSnapshot, SimulationSnapshot + builders
    recorder.py                  # Recorder protocol, BatchRecorder, BatchResult
    lifecycle.py                 # RunContext, _LockableList, raise_if_running
  reactions/                     # Reaction framework
    stoichiometry.py             # StoichiometryEntry + elemental balance
    kinetic.py                   # KineticReaction
    equilibrium.py               # EquilibriumReaction (log_K, single- and cross-phase)
    blackbox.py                  # BlackBoxReactionModel
    reaction_system.py           # ReactionSystem (pre-bucketed container)
    builder.py                   # ReactionBuilder (aerobic_growth, from_coefficients)
    environment.py               # ReactionEnvironment
    protocols.py                 # ReactionModel protocol
    _shared.py                   # Shared validation helpers
  speciation/                    # Acid-base speciation engine
    engine.py                    # BisectionChemicalEquilibriumEngine + level dispatch (1, 2, 2.5)
    acid_base.py                 # _CANONICAL_NAMES, ladder definitions
    activity.py, activity_models.py, sit.py, strong_ions.py, …
    chemistry_level1.py, chemistry_level2.py, carbamate_level25.py,
      chemistry_level3.py        # level=3 is a stub
    legacy_adapter.py            # speciate_multi_acids_TIC for level ≤ 2 callers
    api.py, factory.py, associations_level2.py
  equilibria/                    # Legacy CO₂/O₂/N₂ partition (CUFermenter island)
    engine.py                    # ProcessCoupledEquilibrator
    coupled.py, factory.py, interfaces.py, vle.py, peng_robinson.py
  control/                       # Controllers and profiles
    actions.py                   # ControlAction, ProfileRecord (new framework)
    cv_loops.py                  # CV-native PHController, DO controllers, pressure-relief
    cv_profiles.py               # TemperatureRamp, VVMSchedule, SetpointTrajectory
    interfaces.py                # Controller Protocol (compute → ControlAction)
    loops.py                     # Legacy FermenterState controllers (CUFermenter island)
    system.py                    # Legacy ControlSystem (CUFermenter island)
    controllers/, actuators/, builders/   # Legacy support (CUFermenter island)
  monitoring/                    # AccuracyMonitor, ConservationMonitor
  kinetics/                      # Plug-in kinetic model protocol + worked models
    core.py, mapping.py, models/yeast_acetate_v1.py, …
  numerics/                      # Shared numerical methods
    spatial.py                   # Advection (upwind, TVD) + dispersion
  properties/                    # Physical property models (viscosity)
  param/                         # Parameter sweep utilities
  sim/                           # Legacy types (Ledger, FermenterState) — CUFermenter island
  solvers/                       # Legacy solver dispatch (CoupledSolver) — CUFermenter island

models/
  vlmodels/
    __init__.py
    headspace.py
    adm1/                        # ADM1 + BSM2 — built on Simulation
      base.py, bsm2.py, bsm2_direct.py
    fermenter/                   # Legacy CUFermentationSpeciation + builder/factory/profiles
      unit.py                    # CUFermentationSpeciation (~3400 LOC, CUFermenter island)
      biosteam.py                # BioSTEAM adapter (CUFermenter island)
      types.py                   # FermenterState, FermenterCommands (CUFermenter island)
      profiles.py
      config/
        builder.py               # FermenterBuilder (fluent; .build() and .build_simulation())
        factory.py               # FermenterFactory + legacy run_batch (CUFermenter island)
        configs.py, kinetics.py
    hplc/                        # HPLC column model (own .simulate() loop; not on Simulation)
      column.py
```

The "CUFermenter island" annotation marks modules slated for
removal in the trigger-gated `CUFERMENTER_SUNSET` phase. New code
must not depend on island-tagged surfaces; the canonical
orchestration pathway is `Simulation` + a `ControlVolume` built
via `FermenterBuilder.build_simulation()`.
