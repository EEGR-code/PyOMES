# PyOMES core library

Subpackage-by-subpackage breakdown of `PyOMES/`. See the top-level
[README.md](../README.md) for installation, a quick start, and the
high-level architecture (`ControlVolume`, the chemical equilibrium
engine, the reaction framework, control loops). This file is the
detailed map for browsing the source; each subpackage also carries
its own module docstring with more detail than is repeated here.

| Subpackage | Contents |
|---|---|
| [`core/`](core/) | Foundational types: `Phase` (`GasPhase`, `LiquidPhase`, `SolidPhase`), `PhaseInterface`, `ControlVolume`, boundaries (`GasFeed`, `LiquidFeed`, `LiquidDrain`, `PressureReliefVent`, `MembraneGasBoundary`), CV-to-CV links (`AdvectiveLink`, `DiffusiveLink`), ODE step solvers (`SequentialAdvanceSolver`, `SimultaneousEulerSolver`, `SimultaneousAdaptiveSolver`), the `Simulation` orchestrator, and result/recorder types (`BatchResult`, `SummaryResult`). |
| [`chemistry/`](chemistry/) | Species definitions and phase-partition models (`PartitionModel`, `MultispeciesVLEPartition`) used to construct feeds and initial conditions. |
| [`chemical_equilibrium/`](chemical_equilibrium/) | Aqueous acid-base speciation: `BisectionChemicalEquilibriumEngine` (and Newton-Raphson engine variants) solve equilibrium from a set of declared reactions rather than a fixed tier of hardcoded chemistry. Includes activity-coefficient models (ideal, Davies, SIT) and an optional PHREEQC bridge (`phreeqc` extra). |
| [`reactions/`](reactions/) | Reaction declarations, grouped by kind. `kinetic/` holds `KineticReaction`, the growth rate laws (`Monod`, `Contois`, `Andrews`, and five others) and `ReactionBuilder`, which constructs stoichiometrically validated `KineticReaction` objects from an organism formula and balance mode; elemental-balance errors (`StoichiometryError`) raise at construction. `equilibrium/` holds `EquilibriumReaction`, the `EquilibriumConstraint` protocol, and `HenryEquilibrium`, `RaoultEquilibrium` and `KspEquilibrium`, gas-liquid and solid-liquid equilibrium constraints that double as `PartitionModel` instances. At the top level, `ReactionSystem` aggregates reactions for a `ControlVolume`; `BlackBoxReactionModel` wraps external kinetic functions behind the `ReactionModel` protocol. |
| [`databases/`](databases/) | `ChemistryDatabase` — a frozen bundle of `ThermoFramework` + species + `ReactionSystem` + partition models that a `ControlVolume` accepts. Stock databases (`AQUEOUS_DEFAULT`, `BIOPROCESS_BASIC`, `AD_BASIC`) compose by `.extend()`, not mutation. |
| [`control/`](control/) | Feedback control loops: `PHController`, `DOAgitationController`, `DOCascadeController`, and pressure-relief controllers (`PressureReliefController` and its instant/smooth variants). |
| [`properties/`](properties/) | Physical property models — currently viscosity correlations behind the `ViscosityModel` protocol. |
| [`numerics/`](numerics/) | Spatial discretization schemes (upwind, TVD, dispersion) for spatially resolved models. |
| [`monitoring/`](monitoring/) | Cheap per-step accuracy and mass-conservation checks that warn when a simulation runs outside the regime its solver/activity model is reliable in. |
| [`templates/stirred_tank/`](templates/stirred_tank/) | `StirredTankBuilder` — the fluent builder used in the top-level Quick Start — plus its factory, presets, and kinetics helpers. |
| [`thermo/`](thermo/) | `ThermoFramework` and van 't Hoff helpers at the top level; [`liquid/`](thermo/liquid/) holds the liquid-phase activity models (`IdealLiquidModel`, `DaviesLiquidModel`, `SITLiquidModel`) and the water property correlations they use, which back `chemical_equilibrium/`; [`gas/`](thermo/gas/) holds the gas-phase equations of state (`IdealGasEOS`, `PengRobinsonEOS`). |
| [`compounds.py`](compounds.py) | `ChemicalRegistry` / `Chemical` — standalone named-compound database (molecular weights, atom compositions), decoupled from `Species`. Backs the stirred-tank template's default organism/substrate composition lookup. |
| [`units.py`](units.py) | Shared unit conversions and physical constants. |
| [`config.py`](config.py) | Package-level accuracy-warning thresholds and throttling, configurable via the `VLSIM_WARNINGS` environment variable. |

For concrete, runnable models built on top of this library (anaerobic
digestion, HPLC columns), see the separate [`models/`](../models/)
package and its [README](../models/README.md).
