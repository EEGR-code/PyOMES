# PyOMES

**Package:** PyOMES (installed as the `PyOMES` Python distribution)
**Version:** 0.12.5
**Domain:** Bioprocess simulation (fermentation, anaerobic digestion)
**Language:** Python ≥ 3.10
**Core dependencies:** NumPy, SciPy

---

## What is PyOMES?

A modular framework for simulating fermentation bioreactors and related bioprocesses. 

Core capabilities:
- Warnings - built-in diagnostics flag common model problems proactively, across every capability below
- Reactions - define (bio)chemical reaction phenomena explicitly or use third-party models as plug-ins
- Phase equilibrium - predict gas-liquid-solid partition behaviour using built-in thermodynamic models (Henry's law + coupled VLE, precipitation equilibrium)
- Solution chemistry - track solution pH and chemical ionization states
- Process dynamics - predict time-resolved process behaviour using ODE and ODE-DAE formulations (SciPy-backed solver)
- Control - implement and assess control frameworks using custom or pre-built controller definitions (e.g. pressure relief, pH dosing, and dissolved oxygen controllers)
- Multi-zone systems - multi-phase, multi-zone equipment modeling


---

## Installation

PyOMES isn't published on PyPI yet, so install it from a local copy of the repository:

1. Download the repository — either clone it, or download and extract a ZIP of it — and move into that directory:
   ```bash
   git clone https://github.com/MGuo-Lab/PyOMES.git
   cd PyOMES
   ```
2. Install the package in editable mode:
   ```bash
   python -m pip install -e .
   ```

Install optional dependency groups as needed (see [Optional Features](#optional-features) for what each one provides):

```bash
python -m pip install -e ".[test]"
python -m pip install -e ".[phreeqc]"
python -m pip install -e ".[export]"
python -m pip install -e ".[all]"
```

Concrete model implementations (anaerobic digestion, HPLC columns, ...) live in a separate `models/` package — see [models/README.md](models/README.md) for installing and using it.

---

## Quick Start

PyOMES models are built from `ControlVolume` / `Simulation` machinery. For new users with standard modelling interests, working with this machinery directly may require unnecessary learning. To enable rapid deployment of standard model formulations, templates are available via the `templates` subpackage.

For example, here's a 2 L batch fermenter with pH control, built with the `templates` subpackage's fluent `StirredTankBuilder`:

```python
from PyOMES.templates.stirred_tank import StirredTankBuilder
from PyOMES.control import PHController

builder = (
    StirredTankBuilder()
    .vessel(V_total_L=2.0, T_K=305.15)
    .gas_feed(vvm_min=1.0, composition={"O2": 0.21, "N2": 0.79})
    .transfer_kinetic(kLa_O2=150.0)
    .chemistry()
    .organism("Yeast")
    .substrate("AceticAcid", mu_max=0.5, Ks=5e-3, yield_gX_gS=0.36)
    .controller(PHController(setpoint=5.0, Kp=0.5, Ki=0.0))
)
sim = builder.build_simulation()
result = sim.run(tau_h=5.0, n_steps=1000)

print(f"Final pH: {result.pH['main'][-1]:.2f}")
# result.liquid_mol["main"][species] -> time series (mol)
```

See [docs/tutorials/templates/](docs/tutorials/templates/) for more worked examples (chemostat, fed-batch, and microplate topologies), and [docs/tutorials/](docs/tutorials/) more broadly for tutorials that build up the underlying `Phase` / `PhaseInterface` / `ControlVolume` plumbing this builder wraps.

---

## Repository Layout

```
PyOMES/                 # Core library — see PyOMES/README.md for a subpackage-by-subpackage breakdown
models/                 # Concrete model implementations, installed separately as `vlmodels` — see models/README.md
tests/                  # Test suite: standalone/, validation/, legacy/, performance/
docs/                   # Tutorials and design/development documentation
```

---

## Core Abstractions

### ControlVolume
The central simulation object. Owns a dict of `Phase` objects (gas, liquid, solid), a list of `PhaseInterface` objects (transfer mechanisms), and optional property and reaction models. Multiple `ControlVolume` objects can be connected through a `Simulation`'s `links` (see Multi-CV Systems below) to model multi-zone reactors.

### Chemical Equilibrium
`BisectionChemicalEquilibriumEngine` solves aqueous acid-base speciation from a set of declared equilibrium reactions, rather than a fixed tier of hardcoded chemistry. Activity correction is optional: ideal by default, or a Davies/SIT ionic-strength correction, selected on the `StirredTankBuilder` via `chemistry(use_activity=..., activity_model=...)`.

### Reaction Framework
`ReactionBuilder` constructs stoichiometrically validated `KineticReaction` / `EquilibriumReaction` objects from an organism formula, substrate, and balance mode (e.g. `"CHO"`, `"CHON"`) — elemental-balance errors (`StoichiometryError`) raise at construction, not simulation time. `ReactionSystem` aggregates reactions for a `ControlVolume`; `BlackBoxReactionModel` wraps external kinetic functions behind the `ReactionModel` protocol.

### Control System
Each controller follows a sense → compute → actuate cycle:

| Controller | Action |
|-----------|--------|
| `PressureReliefController` | Vent gas phase to maintain set-point pressure |
| `PHController` | Dose acid or base; re-runs speciation after dosing |
| `DOAgitationController` / `DOCascadeController` | Adjust kLa or agitation rate |

### Multi-CV Systems
`Simulation` connects multiple `ControlVolume` objects via a `links` list of `AdvectiveLink` (bulk circulation) or `DiffusiveLink` (concentration-driven transfer) objects. Used for multi-zone reactors (e.g., a sparger zone with high kLa vs. a bulk zone).

---

## Key Design Patterns

- **Protocol-based extensibility** — `ReactionModel`, `PhaseInterface`, `PropertyCalculator`, `ViscosityModel` are `typing.Protocol` types; no inheritance required.
- **Immutable snapshots** — `Phase.snapshot()` produces independent copies for logging and diagnostics.
- **Fail-fast validation** — stoichiometric and configuration errors raise at construction, not at runtime.
- **Lightweight I/O types** — `FeedState` decouples feed-composition data from any particular process-simulation tool, enabling standalone use.

---

## Running Tests

```bash
python -m pip install -e ".[test]"
python -m pytest
```

By default (see `pyproject.toml`), `pytest` runs `tests/standalone/` and `tests/validation/`. `tests/legacy/` and `tests/performance/` are not part of the default run — `tests/legacy/` predates the current reaction-driven chemical equilibrium engine, and `tests/performance/` is reserved for future runtime benchmarks.

### Test Coverage

As of this writing, `tests/standalone/` and `tests/validation/` cover:

| Area | Representative test files |
|------|-----------|
| Core control-volume mechanics | `test_core.py`, `test_cv_advance.py`, `test_cv_compute_interface.py`, `test_clamping.py`, `test_state_vector.py` |
| Gas-liquid transfer & boundaries | `test_gas_liquid_link.py`, `test_transfer_models.py`, `test_membrane.py`, `test_boundaries.py`, `test_feed_state.py` |
| Simulation orchestration | `test_simulation.py` |
| ODE solvers | `test_system_solver.py`, `test_implicit_transport_solver.py`, `test_monolithic_ode_solver.py`, `test_simultaneous_adaptive_solver_jac.py` |
| Aqueous chemical equilibrium / speciation | `test_speciation.py`, `test_speciation_protocols.py`, `test_nr_speciation_engine.py`, `test_equilibrium_classification.py`, `test_activity_dispatch.py`, `test_strong_ions.py` |
| PHREEQC cross-validation | `test_phreeqc_engine.py`, `tests/validation/speciation/test_phreeqc_nr_agreement.py` |
| NIST / analytical reference validation | `tests/validation/speciation/test_carbonate_phosphate_benchmarks.py`, `test_iron_oxidation.py`, `test_saturation_index.py` |
| Chemistry database & species | `test_chemistry_database.py`, `test_compounds.py`, `test_species.py`, `test_partition_model.py`, `test_thermo_framework.py`, `test_liquid_phase_model.py` |
| Reaction stoichiometry & kinetics | `test_reactions.py`, `test_stoichiometry.py`, `test_equilibrium_constraint.py`, `test_kinetics.py` |
| Control loops | `test_controller_state_protocol.py`, `test_descriptors.py`, `test_param_path.py` |
| StirredTank templates | `test_builder.py`, `test_configs.py` |
| Physical properties & numerics | `test_viscosity.py`, `test_spatial_schemes.py` |
| Monitoring | `test_accuracy_monitor.py`, `test_conservation_monitor.py` |
| Sub-models (HPLC, ADM1/BSM2) | `test_hplc_column.py`, `test_bsm2_reference.py` |

This table groups ~60 test modules thematically rather than listing all of them — browse [`tests/standalone/`](tests/standalone/) and [`tests/validation/`](tests/validation/) for the current, definitive list.

---

## Optional Features

| Extra | Dependencies | Purpose |
|-------|--------------|---------|
| `phreeqc` | `phreeqpython>=1.6` | PHREEQC-backed chemistry comparisons and validation |
| `export` | `pandas`, `pyarrow` | Tabular export workflows |
| `test` | `pytest`, `pandas`, `pyarrow`, `phreeqpython>=1.6` | Development and validation suite |
| `all` | `phreeqpython>=1.6`, `pandas`, `pyarrow` | All optional runtime features |

---

## How to Cite

If you use PyOMES in your work, please cite:

> Errington, E., Vinestock, T., Lee, J., & Guo, M. (2026). *PyOMES: an open-source framework for biochemical process modelling*. arXiv:2608.06360 [q-bio.QM]. https://arxiv.org/abs/2608.06360

```bibtex
@misc{errington2026pyomes,
  title         = {PyOMES: an open-source framework for biochemical process modelling},
  author        = {Errington, Ethan and Vinestock, Tom and Lee, Jaewook and Guo, Miao},
  year          = {2026},
  eprint        = {2608.06360},
  archivePrefix = {arXiv},
  primaryClass  = {q-bio.QM},
  url           = {https://arxiv.org/abs/2608.06360}
}
```

A machine-readable citation is also available in [CITATION.cff](CITATION.cff).

---

## License

This project is distributed under the license in [LICENSE](LICENSE).
