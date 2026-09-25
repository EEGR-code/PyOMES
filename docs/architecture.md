# Architecture: CV-Centric PyOMES Simulation

This document describes the internal architecture of the PyOMES simulation
module as of v0.12.5.  Concrete models (ADM1/BSM2, the HPLC column)
live in the separate `vlmodels` package under `models/`.

## Overview

The simulation is structured around a **ControlVolume** (CV) that owns all
thermodynamic state (mole inventories, equilibrium, derived properties). The
`Simulation` orchestrator owns the time loop: each step it applies open-loop
profiles and inter-CV links, advances every CV, then applies the actions its
controllers return (pressure relief, pH dosing, DO control).

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
- **Internal interfaces** — a `KineticGasLiquidLink` for gas-liquid
  transfer, kinetic (kLa-limited) or at equilibrium per species. It is
  built from the CV's `transfer_models=` argument (one
  `KineticTransferModel` or `EquilibriumTransferModel` per species).
- **Boundaries** — exchanges with the outside: gas and liquid feeds,
  liquid drains, headspace overpressure outlets, and gas-permeable
  membranes (`cv.boundaries`).
- **Reaction system** (optional) — one `ReactionSystem` holding all
  reaction declarations on the CV, pre-bucketed by kind.
- **Property calculators** (optional) — list of `PropertyCalculator`
  instances each computing one scalar derived property (viscosity,
  density, …) onto `phase.properties[key]`.
- **Monitors** — `AccuracyMonitor` and `ConservationMonitor`
  attached on the reaction system surface (the CV auto-attaches a
  default `ConservationMonitor` in `__init__`).

`advance(dt_h, t_h, external_source_terms=None, solver=None)` hands
the step to a `StepSolver` (see "Step solvers" below). With
`solver=None` it uses `SequentialAdvanceSolver`, which sequences:

1. Speciation solve via `reaction_system.engine.solve(phases=...)` —
   writes derived molecular species (`CO2aq`, `HCO3-`, `NH3`, `H+`, …)
   back to `phase.n_mol` via the privileged
   `phase._refresh_derived(values)` hook.
2. Property calculators (viscosity, density, …) — write to
   `phase.properties`, after speciation so they can read derived
   species.
3. Apply `external_source_terms` (orchestrator-supplied fluxes).
4. Apply boundary fluxes (`cv.boundaries`).
5. Integrate kinetic + black-box reactions (post-feed state; reads
   properties via `env.prop(key)`).
6. Internal transfer (`KineticGasLiquidLink`, with Henry constants
   corrected via molecular fractions resolved from `phase.n_mol`
   against the system's `cross_phase_equilibria`).
7. `ConservationMonitor.check_step(phases)` — element + charge drift
   guard.

The CV returns an `AdvanceResult` containing transfer diagnostics
and reaction source terms. pH is read from `liquid.n_mol["H+"]` and
ionic strength is recomputable from the canonical `n_mol` directly.

### Gas-liquid CV pattern

A stirred tank is a plain `ControlVolume` whose `phases` dict is
`{"gas": GasPhase, "liquid": LiquidPhase}` and whose
`internal_interfaces` list contains a `KineticGasLiquidLink` acting
as a `PhaseInterface`.  Callers construct one via
`StirredTankBuilder().build()` (or `StirredTankFactory.create_volume(...)`).

Controllers never touch the link directly: they read a snapshot and
change link parameters such as kLa through `params_changed` paths
(see "Control System" below). Derived species are first-class in
`n_mol`; scalar properties live on `phase.properties`.

### Step solvers (`PyOMES/core/solvers.py`)

A `StepSolver` advances one CV by one step and returns an
`AdvanceResult`, which carries `transfer_record` and
`boundary_records`. Three are provided:

- `SequentialAdvanceSolver` (the default) — runs the sub-steps one
  after another, in the order listed above.
- `SimultaneousEulerSolver` — computes every sub-system's deltas
  from one frozen snapshot of the state and applies them together,
  so the result does not depend on sub-step order.
- `SimultaneousAdaptiveSolver` — integrates the same combined
  right-hand side with `scipy.integrate.solve_ivp`, solving
  speciation algebraically at every internal step.

Only `SequentialAdvanceSolver` runs property calculators; under the
two simultaneous solvers `phase.properties` is not refreshed.

### Phases (`PyOMES/core/phases.py`)

`GasPhase` and `LiquidPhase` hold mole inventories in a single
canonical `n_mol` dict and volume/temperature. Totals and derived
species share that one store (PHREEQC-style): the speciation engine writes derived molecular species (`H+`,
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

Species with no gas-liquid transfer model and no declared
equilibrium (e.g. biomass, substrates) are stored in the liquid
phase but are not touched by the speciation engine or the
gas-liquid link. They participate in
inter-zone transport via advective links.

### Property Calculators (`PyOMES/core/property_calculator.py`)

A `PropertyCalculator` is a narrow protocol with a `key: str` and
`compute(phase, T_K, P_atm) -> float`. Each calculator computes
**one scalar derived property** (viscosity, density, heat capacity,
…) and writes it to `phase.properties[key]`. Rate laws read it via
`env.prop("viscosity")`.

Speciation is *not* a `PropertyCalculator`. Acid-base equilibria are
state-completion handled by the speciation engine of
`cv.reaction_system` (see "Equilibrium pathway" below); the engine
writes derived molecular species straight to `phase.n_mol`, not to
`phase.properties`.

Property calculators are registered via
`cv.property_calculators=[...]` and evaluated by
`SequentialAdvanceSolver` on every step, after the speciation solve
and before reaction integration.

### Equilibrium pathway

Equilibria are declared on the `ReactionSystem` (as
`EquilibriumReaction`, `HenryEquilibrium`, `KspEquilibrium`, …) and
solved by the system's speciation engine at the start of each step;
the engine writes derived molecular species back to `phase.n_mol`
via the privileged `_refresh_derived(values)` hook. The engine is
chosen by `ReactionSystem(solver=...)`:

- `"charge_balance"` (default) — `BisectionChemicalEquilibriumEngine`,
  which solves the charge balance for pH with a bracketed root-finder.
  It has no precipitation support.
  ADM1, BSM2 and `StirredTankBuilder` models use it.
- `"newton_raphson"` — `NRChemicalEquilibriumEngine`, a tableau
  Newton-Raphson solver that also handles gas-liquid rows folded
  into the tableau and mineral precipitation.

Gas-liquid transfer runs separately, through `KineticGasLiquidLink`:
kLa-limited or instantaneous per species, with the Henry constant
alpha-corrected inline from `phase.n_mol` against the reaction
system's `cross_phase_equilibria` bucket.

## Reaction Framework (`PyOMES/reactions/`)

### Three kinds of declaration

A reaction is declared as one of three kinds — no shared base
class; shared validation lives as free functions in `_shared.py`:

- `KineticReaction` — stoichiometric template + callable rate law.
  Stoichiometry is validated at construction for elemental balance
  (by default over every element that appears in it);
  `StoichiometryError` carries a per-species breakdown.
- Equilibrium constraints — `EquilibriumReaction` (stoichiometric
  template + `log_K`, with optional Van 't Hoff parameters
  `dH_J_per_mol` / `T_ref_K`) and the specialised `HenryEquilibrium`,
  `RaoultEquilibrium` and `KspEquilibrium`, all conforming to the
  `EquilibriumConstraint` protocol. Each is classified by the phases
  it spans: single-phase constraints go to the speciation engine,
  gas-liquid ones to the `KineticGasLiquidLink` as partition
  declarations, and solid-liquid ones to the precipitation loop of
  the NR engine. None has `compute_rates`.
- `BlackBoxReactionModel` — adapter for opaque external simulators
  (FBA, genome-scale, proprietary). Stoichiometry is not
  inspectable, so elemental balance is checked at runtime against
  the configured `flux_mapping` with a policy (warn, raise, log,
  ignore).

### ReactionSystem

`ReactionSystem` is the **single attach point** for every reaction
on a CV (`cv.reaction_system`). At construction it pre-buckets its
input into `_kinetic_reactions`, `_blackbox_models` (by
`isinstance`), and `_single_phase_equilibria`,
`_cross_phase_equilibria`, `_precipitation_equilibria` (by
`classify_equilibrium_constraint`) — each consumer pulls from its
own bucket, and mixed-kind systems are always safe to integrate
(equilibria are simply not in the kinetic iteration path).

It also owns the lazy speciation engine via the `engine` property.
On first access the engine selected by `ReactionSystem(solver=...)`
is built with `from_reactions(...)`: the Bisection engine from the
single-phase and gas-liquid constraints, the NR engine from those
plus the solid-liquid ones. Activity settings are pinned beforehand
by `configure_engine(use_activity=, activity_model=)`; the solver
choice itself is fixed at construction. Models and tests needing a
pre-built engine inject it via `attach_engine(engine)`.
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
models. It holds a named dict of ControlVolumes, an optional list
of inter-CV `CVLink` objects, controllers, profiles, a step solver
(or per-CV dict), and a recorder. The constructor wires every CV
(and its nested phases, links, lockable lists) to a shared
`RunContext` that the lifecycle gate consults during `.run()`.
`save_checkpoint` / `load_checkpoint` write and restore a run's
full state.

`Simulation.run(tau_h, n_steps)` is the user-facing entry point.
If a `system_solver=` is given (`PyOMES/core/system_solver.py`:
explicit Euler, Strang splitting, multirate, implicit transport, or
a monolithic ODE solve), it advances the whole multi-CV system each
step. Otherwise the default per-step body is:

1. **Profiles** (open-loop time-varying mutators) fire *before*
   integration so the advance sees updated state. Each profile's
   `apply(t_h, sim)` returns a `ProfileRecord` for audit.
2. **Inter-CV link flows** are computed and applied symmetrically
   (source loses, sink gains), yielding one `LinkFlowRecord` per
   link.
3. **Each CV advances** via `cv.advance(dt_h, t_h, solver=...)`
   with its dispatched step solver (by default the sequential
   body listed under ControlVolume: speciation, property
   calculators, source terms and boundary fluxes, reactions,
   internal transfer).
4. **A `SimulationSnapshot`** is built from the post-advance
   state.
5. **Controllers** (respecting sample-period gating; zero-order
   hold between samples) consume the snapshot and return a
   `ControlAction` describing the mutation they want applied.
6. **The orchestrator applies each action**: `flux_applied` via
   `cv.apply_external_flux`; each `params_changed` key (a path
   string such as `"internal_interfaces[KineticGasLiquidLink].kLa.O2"`,
   a parsed `ParamPath`, or a class-level descriptor handle such as
   `KineticGasLiquidLink.kLa["O2"]`) is walked from the CV and the
   value written through the target's unchecked setter.
7. The recorder receives per-step records and finalises into a
   `BatchResult` (per-CV nested arrays).

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
run (the legitimate case), every gated mutator has an **unchecked
sibling**. Most gated attributes (phase `T_K`, the link's `kLa`,
a gas feed's composition) are declared with the `MutableScalar` /
`MutableDict` descriptors in `PyOMES/control/descriptors.py`, which
supply the gate and the unchecked write in one declaration; the rest
have explicit `_set_*_unchecked` methods (e.g.
`_set_reaction_system_unchecked`, `Simulation._set_controllers_unchecked`).
User code goes through the gated public attribute; the orchestrator's
action-apply path uses the unchecked write.

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
    Simulation, AdvectiveLink, KineticTransferModel,
)

# Sparger zone: small, high kLa, has the gas phase
cv_sparger = ControlVolume(
    phases={"gas": GasPhase(...), "liquid": LiquidPhase(...)},
    transfer_models={"O2": KineticTransferModel(o2_henry, k_transfer=150.0)},
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

Controllers consume a typed `CVSnapshot`
(single-CV) or `SimulationSnapshot` (multi-CV) and return a
structured `ControlAction` from a single
`compute(state, dt_h) -> ControlAction` method. Internal state
(setpoint, integral, sampling-hold buffer) lives on the controller
instance.

1. `Simulation.run` calls `ctrl.reset()` at run-entry on every
   controller, so a re-run always starts from fresh controller state.
2. Per step (post-advance), the orchestrator builds the snapshot,
   filters by sample-period (zero-order hold between samples),
   and invokes `ctrl.compute(snapshot, dt_h)`.
3. The returned `ControlAction` carries `flux_applied`
   (`{phase: {species: mol/h}}`), `params_changed`
   (`{semantic_path: value}`), `vented_mol`, and `dosed_mol`.
4. The orchestrator applies the action: `flux_applied` via
   `cv.apply_external_flux`; `params_changed` via the path walker
   described under Simulation, step 6.

Concrete controllers ship in `PyOMES/control/cv_loops.py`:
`PHController`, `DOAgitationController` (+`DOController` alias),
`DOCascadeController`, `InstantPressureReliefController`,
`SmoothPressureReliefController`, `PressureReliefController`.
Profiles (open-loop time-varying mutators) ship in
`PyOMES/control/cv_profiles.py`: `TemperatureRamp`, `VVMSchedule`,
`SetpointTrajectory`.

## Directory layout

```
PyOMES/
  __init__.py
  units.py                       # Shared constants and unit conversions
  config.py                      # PyOMES.config — WarningConfig, env-var presets
  compounds.py                   # ChemicalRegistry, Chemical — standalone compound database
  chemistry/                     # species.py (Species), common_species.py (inorganic aqueous species),
                                 # species_check.py (cross-reaction consistency), partition.py
                                 # (phase-partition protocols, ideal-gas VLE model)
  databases/                     # ChemistryDatabase (species + reactions + ThermoFramework bundle) and
                                 # stock databases: aqueous, anaerobic digestion, basic bioprocess
  core/                          # Framework core: phases, CV, orchestration
    phases.py                    # GasPhase, LiquidPhase, SolidPhase
    interfaces.py                # PhaseInterface, TransferDiagnostics, AdvanceResult
    control_volume.py            # ControlVolume.advance(dt_h, t_h, …)
    gas_liquid_link.py           # KineticGasLiquidLink (kinetic or equilibrium transfer + alpha)
    transfer_models.py           # KineticTransferModel, EquilibriumTransferModel (per-species transfer)
    property_calculator.py       # PropertyCalculator protocol
    solvers.py                   # StepSolver protocol; SequentialAdvanceSolver (default),
                                 # SimultaneousEulerSolver, SimultaneousAdaptiveSolver
    clamping.py                  # Non-negativity clamping for discrete-step solvers
    state_vector.py              # State-vector packing for ODE-style solvers
    system_solver.py             # SystemSolver protocol + multi-CV integrators (splitting, multirate, …)
    system_env.py                # SystemEnv — read-only system state view for controllers
    links.py                     # CVLink, AdvectiveLink, DiffusiveLink, LinkFlowRecord
    boundaries.py                # External boundaries: gas/liquid feeds, drains, overpressure outlets,
                                 # membranes
    simulation.py                # Simulation orchestrator: run loop, action dispatch, checkpoints
    snapshot.py                  # CVSnapshot, SimulationSnapshot + builders
    recorder.py                  # Recorder protocol, BatchRecorder → BatchResult, plus streaming,
                                 # sparse and summary recorders
    lifecycle.py                 # RunContext, _LockableList, raise_if_running
  reactions/                     # Reaction framework
    stoichiometry.py             # StoichiometryEntry + elemental balance
    environment.py               # ReactionEnvironment
    protocols.py                 # ReactionModel protocol
    _shared.py                   # Shared validation helpers
    reaction_system.py           # ReactionSystem (pre-bucketed container)
    blackbox.py                  # BlackBoxReactionModel
    kinetic/                     # reaction.py (KineticReaction), rate_laws.py (GrowthKinetics, Monod,
                                 # Contois, Andrews, ...), builder.py (ReactionBuilder: aerobic_growth,
                                 # monod_aerobic_growth, from_coefficients)
    equilibrium/                 # constraint.py (EquilibriumConstraint, vant_hoff_log_K,
                                 # classify_equilibrium_constraint), reaction.py (EquilibriumReaction:
                                 # log_K, single- and cross-phase), interphase.py (HenryEquilibrium,
                                 # RaoultEquilibrium, KspEquilibrium), plots.py (plot_vant_hoff,
                                 # plot_speciation)
  chemical_equilibrium/          # Aqueous equilibrium solvers (acid-base, complexation, folded gas-liquid
                                 # and precipitation rows) built on thermo/'s activity models
    protocols.py                 # Engine protocols, EquilibriumResult
    numerical_gradient.py        # NumericalGradientEquilibriumEngine (wraps any engine)
    engines/
      bisection/                 # engine.py (BisectionChemicalEquilibriumEngine), acid_base.py,
                                 # equilibria.py (EquilibriumSet/EquilibriumDef charge-balance format),
                                 # ionic_strength.py (suffix-based ionic strength)
      nr/                        # engine.py (NRChemicalEquilibriumEngine), tableau.py, solver.py
      phreeqc.py                 # PHREEQCChemicalEquilibriumEngine (optional phreeqpython)
  thermo/                        # Stateless models and conventions (no solvers; depends only on units.py)
    framework.py                 # ThermoFramework (liquid_activity + gas_eos), THERMO_IDEAL, THERMO_DAVIES
    equilibrium_constants.py     # van 't Hoff helpers (vant_hoff_K, vant_hoff_log_K)
    liquid/                      # protocols.py (LiquidPhaseModel, ActivityModel, DifferentiableLiquidModel),
                                 # ideal.py, davies.py, sit.py (activity models), water_properties.py,
                                 # factory.py (make_activity_model)
    gas/                         # protocols.py (GasEOS), ideal.py (IdealGasEOS), peng_robinson.py
                                 # (PengRobinsonEOS, CriticalProperties, BIOGAS_SPECIES, BIOGAS_KIJ)
  control/                       # Controllers and profiles
    actions.py                   # ControlAction, ProfileRecord
    cv_loops.py                  # PHController, DO controllers, pressure-relief controllers
    cv_profiles.py               # TemperatureRamp, VVMSchedule, SetpointTrajectory
    interfaces.py                # Controller protocol (compute → ControlAction), ControllerBase
    descriptors.py               # MutableScalar, MutableDict (lifecycle-gated attributes)
    param_path.py                # ParamPath — parsed params_changed paths
  monitoring/                    # AccuracyMonitor, ConservationMonitor
  numerics/                      # Shared numerical methods
    spatial.py                   # Advection (upwind, TVD) + dispersion
  properties/                    # Physical property models (viscosity)
  templates/
    stirred_tank/                # builder.py (StirredTankBuilder), factory.py (StirredTankFactory),
                                 # configs.py (config dataclasses), profiles.py (time profiles)

models/                          # Installable vlmodels package (pyproject.toml, setup.py)
  vlmodels/
    __init__.py
    headspace.py                 # Ideal-gas headspace helpers, clamp()
    adm1/                        # base.py (ADM1), bsm2.py (BSM2) — ControlVolume builders;
                                 # bsm2_direct.py (BSM2DirectModel, a native-unit ReactionModel)
    hplc/                        # HPLC column model (own .simulate() loop; not on Simulation)
      column.py
```

The canonical orchestration pathway is `Simulation` + a
`ControlVolume` built via
`PyOMES.templates.stirred_tank.StirredTankBuilder.build_simulation()`.
