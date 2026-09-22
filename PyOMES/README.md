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
| [`chemistry/`](chemistry/) | `ChemicalRegistry` / compound database and species definitions used to construct feeds and initial conditions. |
| [`chemical_equilibrium/`](chemical_equilibrium/) | Aqueous acid-base speciation: `BisectionChemicalEquilibriumEngine` (and Newton-Raphson engine variants) solve equilibrium from a set of declared reactions rather than a fixed tier of hardcoded chemistry. Includes activity-coefficient models (ideal, Davies, SIT) and an optional PHREEQC bridge (`phreeqc` extra). |
| [`reactions/`](reactions/) | `ReactionBuilder` constructs stoichiometrically validated `KineticReaction` / `EquilibriumReaction` objects from an organism formula and balance mode; elemental-balance errors (`StoichiometryError`) raise at construction. `ReactionSystem` aggregates reactions for a `ControlVolume`; `BlackBoxReactionModel` wraps external kinetic functions behind the `ReactionModel` protocol. |
| [`control/`](control/) | Feedback control loops: `PHController`, `DOAgitationController`, `DOCascadeController`, and pressure-relief controllers (`PressureReliefController` and its instant/smooth variants). |
| [`properties/`](properties/) | Physical property models — currently viscosity correlations behind the `ViscosityModel` protocol. |
| [`numerics/`](numerics/) | Spatial discretization schemes (upwind, TVD, dispersion) for spatially resolved models. |
| [`monitoring/`](monitoring/) | Cheap per-step accuracy and mass-conservation checks that warn when a simulation runs outside the regime its solver/activity model is reliable in. |
| [`templates/stirred_tank/`](templates/stirred_tank/) | `StirredTankBuilder` — the fluent builder used in the top-level Quick Start — plus its factory, presets, and kinetics helpers. |
| [`thermo/`](thermo/) | `ThermoFramework` and liquid-phase activity models (`IdealLiquidModel`, `DaviesLiquidModel`, `SITLiquidModel`), water property correlations that back `chemical_equilibrium/`, and gas-phase equations of state (`IdealGasEOS`, `PengRobinsonEOS`). |
| [`stream_adapter.py`](stream_adapter.py) | `FeedState` — the canonical feed/broth composition type used throughout PyOMES. |
| [`units.py`](units.py) | Shared unit conversions and physical constants. |
| [`config.py`](config.py) | Package-level accuracy-warning thresholds and throttling, configurable via the `VLSIM_WARNINGS` environment variable. |

For concrete, runnable models built on top of this library (anaerobic
digestion, HPLC columns), see the separate [`models/`](../models/)
package and its [README](../models/README.md).
