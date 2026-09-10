# PyOMES

**Package:** PyOMES (installed as the `PyOMES` Python distribution; the
core library was previously distributed as `fermenter`)
**Version:** 0.12.5
**Domain:** Bioprocess simulation (fermentation, anaerobic digestion)
**Language:** Python ≥ 3.10
**Key dependencies:** NumPy, SciPy (no BioSTEAM required)

---

## Overview

A modular framework for simulating fermentation bioreactors and related bioprocesses. Core capabilities:

- Multi-phase, multi-zone reactor modeling (gas headspace, liquid broth, solid material)
- Aqueous speciation chemistry (pH, ionic strength, carbonate/phosphate/sulfate systems)
- Gas-liquid equilibrium (Henry's law + coupled VLE)
- Stoichiometrically validated reaction networks with fail-fast elemental balance checks
- Feedback control systems (pressure relief, pH dosing, dissolved oxygen)
- ODE integration for kinetic models (single Euler in the sequential body; SciPy adaptive via `SimultaneousAdaptiveSolver`)

---

## Installation

For local development, install the package in editable mode from the repository root:

```bash
python -m pip install -e .
```

Install optional dependency groups as needed:

```bash
python -m pip install -e ".[test]"
python -m pip install -e ".[phreeqc]"
python -m pip install -e ".[export]"
python -m pip install -e ".[all]"
```

The concrete model implementations used by several examples live in `models/`.
When running scripts directly from a checkout, include that directory on
`PYTHONPATH`:

```bash
export PYTHONPATH="$PWD/models:$PYTHONPATH"
```

---

## Quick Start

```python
from vlmodels.fermenter import create_standalone_fermenter
from PyOMES.stream_adapter import FeedState

F = create_standalone_fermenter(
    tau=2,           # residence time (h)
    T=305.15,        # temperature (K)
    V_total_L=2.0,
    organism_id="Yeast",
    balance_basis="CHO",
)

feed = FeedState.from_mass_concentrations({"AceticAcid": 1.0, "Yeast": 0.1})
result = F.simulate(feed, tau=2)
# result.time_series → dict of numpy arrays (pH, pressure, composition, ...)
```

---

## Running Tests

The default pytest configuration runs the standalone suite:

```bash
python -m pytest
```

The test extra installs the optional packages used by the broader test and
export paths:

```bash
python -m pip install -e ".[test]"
```

---

## Repository Layout

```
PyOMES/                 # Main package
  core/                 # Foundational abstractions: ControlVolume, GasPhase, LiquidPhase, PhaseInterface
  chemistry/            # Chemical database, compound registry, solution recipe builders
  speciation/           # Acid-base equilibrium engine (3 levels of detail)
  equilibria/           # Coupled gas-liquid VLE (Henry's law + speciation)
  reactions/            # Reaction / ReactionSet with stoichiometry validation
  control/              # ControlSystem, pressure/pH/DO controllers and actuators
  properties/           # Physical property models (viscosity correlations)
  sim/                  # FermenterState, ControlAction, RunResult types
  solvers/              # ODE integration strategies (CoupledSolver, SimultaneousAdaptiveSolver, SimultaneousEulerSolver)
  kinetics/             # KineticModel protocol and plug-in system
  numerics/             # Spatial discretization (upwind, TVD, dispersion)
  stream_adapter/       # Lightweight FeedState / FermenterResult (no BioSTEAM)

models/                 # Concrete implementations (imported as `vlmodels`)
  fermenter/            # CUFermentationSpeciation — main production fermenter unit
  adm1/                 # Anaerobic digestion (ADM1, BSM2 variants)
  hplc/                 # HPLC column chromatography model

tests/
  standalone/           # 25 unit + integration test modules
  legacy/               # Carbonate/pH/buffer validation against NIST standards
```

---

## Core Abstractions

### ControlVolume
The central simulation object. Owns a dict of `Phase` objects (gas, liquid), a list of `PhaseInterface` objects (transfer mechanisms), and optional `PropertySolver` and `ReactionModel` instances.

```
ControlVolume.advance(dt_h, chem_env) →  AdvanceResult
  1. Evaluate property solvers    (pH, viscosity, ionic strength) →  property_results
  2. Apply external fluxes        (feed, control dosing)
  3. Integrate reactions          (post-feed; receives property_results)
  4. Step internal interfaces     (gas-liquid mass transfer, equilibrium;
                                   receives property_results)
```

The new ordering eliminates the legacy `_last_properties` cache.  See
`CV_UPDATE.md` for the full rationale and `ORDERING.md` for the
feed-before-reactions argument.

### Speciation Engine
Three tiers of acid-base chemistry, selected at construction time:

| Level | Systems included | Notes |
|-------|-----------------|-------|
| 1 | Water + core acid-base | Fastest |
| 2 | Level 1 + phosphate, bisulfate, ion pairing | Default |
| 2.5 | Level 2 + carbamate equilibrium | Most detailed |

Activity models: Ideal or Davies (ionic strength correction).

### Reaction Framework
`Reaction` validates elemental balance (C, H, O, N) at construction — simulation objects can never be built with a stoichiometrically inconsistent reaction. `ReactionSet` aggregates multiple reactions; `BlackBoxReactionModel` wraps external kinetic functions behind the `ReactionModel` protocol.

### Control System
`ControlSystem.step()` follows a sense → compute → actuate cycle:

| Controller | Action |
|-----------|--------|
| `PressureReliefController` | Vent gas phase to maintain set-point pressure |
| `PHController` | Dose acid or base; re-runs speciation after dosing |
| `DOController` / `DOCascadeController` | Adjust kLa or agitation rate |

### Multi-CV Systems
`MultiCVSystem` connects multiple `ControlVolume` objects via `AdvectiveLink` (bulk circulation) or `DiffusiveLink` (concentration-driven transfer). Used for multi-zone fermenters (e.g., sparger zone with high kLa vs. bulk zone).

---

## Key Design Patterns

- **Protocol-based extensibility** — `ReactionModel`, `PhaseInterface`, `PropertySolver`, `ViscosityModel` are `typing.Protocol` types; no inheritance required.
- **Immutable snapshots** — `Phase.snapshot()` produces independent copies for logging and diagnostics.
- **Fail-fast validation** — stoichiometric and configuration errors raise at construction, not at runtime.
- **Lightweight I/O types** — `FeedState` / `FermenterResult` decouple the framework from BioSTEAM `Stream` objects, enabling standalone use.

---

## Test Coverage

| Area | Test files |
|------|-----------|
| Core phases and interfaces | `test_core.py`, `test_gas_liquid_link.py`, `test_gas_liquid_volume.py` |
| Reaction stoichiometry | `test_reactions.py` |
| Speciation / chemistry | `test_speciation.py`, `test_compounds.py` |
| Controllers | `test_controllers.py` |
| Multi-CV links | `test_multi_cv.py` |
| Fermenter construction | `test_fermenter_construction.py`, `test_builder.py`, `test_factory.py`, `test_configs.py` |
| Full integration | `test_integration.py` |
| Sub-models | `test_hplc_column.py`, `test_headspace.py`, `test_kinetics.py` |
| Physical properties | `test_viscosity.py`, `test_spatial_schemes.py` |
| Legacy validation | Carbonate, pH, NIST buffer benchmarks |

---

## Optional Features

| Extra | Dependencies | Purpose |
|-------|--------------|---------|
| `phreeqc` | `phreeqpython>=1.6` | PHREEQC-backed chemistry comparisons and validation |
| `export` | `pandas`, `pyarrow` | Tabular export workflows |
| `test` | `pytest`, `pandas`, `pyarrow`, `phreeqpython>=1.6` | Development and validation suite |
| `all` | `phreeqpython>=1.6`, `pandas`, `pyarrow` | All optional runtime features |

---

## License

This project is distributed under the license in [LICENSE](LICENSE).
