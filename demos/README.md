# PyOMES demos

Runnable tutorials for the `Simulation` orchestrator shipped in
[`simulation-class`](../docs/dev/implementation/shipped/SIMULATION_CLASS.md)
(2026-05-27).

## Prerequisite

Install the package in editable mode from the repo root:

```bash
pip install -e .
```

The local [_bootstrap.py](_bootstrap.py) adds `models/` to `sys.path`
so the demos can import `vlmodels.fermenter.config`. The
`pip install -e .` exposes the `PyOMES` package directly.

## Layout

```
demos/
  builder/                       quick-start: FermenterBuilder fluent API
    batch_fermenter.py
    cstr_fermenter.py
    fed_batch_fermenter.py
    microplate_fermenter.py
  model_api/                     direct framework usage (no builder)
    chemistry/                   reaction / equilibrium / blackbox declarations
      reaction_system.py
      fba/
        fba_toy.py
        fba_ecoli_core.py
        ecoli_core.json
    D2Cworkshop/          manual Simulation assembly from primitives
      basic_layout/         foundational examples (raw construction)
        Example1_mtp_well.ipynb
        Example2_batch_fermenter.ipynb
        raw_construction.py
      updated_layout/       extended scenarios (in development)
  features/                      single-class/single-feature deep dives
    ChemicalEquilibriumProtocol/
      0_README.ipynb           architecture overview (start here)
      01_bisection_engine_basics.ipynb
      02_nr_engine_basics.ipynb
      03_phreeqc_engine_basics.ipynb
    SolverProtocols/
      0_README.ipynb           architecture overview (start here)
      01_writing_a_custom_solver.ipynb   writing your own StepSolver/SystemSolver
  usecases/                      scenario-first, five-minute worked examples
    0_README.ipynb           notebook index (start here)
    01_predict_ph_simple_liquid.ipynb
```

The principle: top-level subtrees name **which framework layer**
the demo teaches. `features/` is the one exception — it's organized by
**class**, not layer, for focused API walkthroughs that don't fit the
layer taxonomy (e.g. "how do I construct and call this one engine class").

| Subtree | When to read it |
|---|---|
| [`builder/`](builder/) | You want to run a fermenter and don't care how it's wired. Five-line fluent API to a configured `Simulation`. Start here. |
| [`model_api/chemistry/`](model_api/chemistry/) | You want to know how to declare reactions — kinetic with a rate law, single-phase equilibria with `log_K`, cross-phase partitions, or opaque black-box solvers (FBA). |
| [`model_api/D2Cworkshop/`](model_api/D2Cworkshop/) | You want to know how to build a `Simulation` end-to-end without the builder — explicit phases, links, boundaries, controllers. |
| [`features/`](features/) | You already know which class you're using and want its instantiation/call conventions and gotchas — e.g. [`features/ChemicalEquilibriumProtocol/0_README.ipynb`](features/ChemicalEquilibriumProtocol/0_README.ipynb) for the architecture overview, then the three notebooks covering `BisectionChemicalEquilibriumEngine`, `NRChemicalEquilibriumEngine`, and `PHREEQCChemicalEquilibriumEngine`'s distinct `solve()` conventions and capabilities; or [`features/SolverProtocols/0_README.ipynb`](features/SolverProtocols/0_README.ipynb) for the `StepSolver`/`SystemSolver` protocol overview, then [`01_writing_a_custom_solver.ipynb`](features/SolverProtocols/01_writing_a_custom_solver.ipynb) for writing your own. |
| [`usecases/`](usecases/) | You have a concrete scenario in mind ("I have a sample, I want to know X") and want a short, worked answer rather than a full API tour — e.g. [`usecases/01_predict_ph_simple_liquid.ipynb`](usecases/01_predict_ph_simple_liquid.ipynb) for predicting the pH of a liquid-only sample with no gas or solid phase. |

The chemistry and construction demos compose: `raw_construction.py`
imports its reactions from `chemistry/reaction_system.py` rather
than redeclaring them, demonstrating that chemistry definitions
are reusable artifacts.

## Reading order

If you're new to the framework:

1. **[`builder/batch_fermenter.py`](builder/batch_fermenter.py)** —
   the canonical build → wire → run → report pattern under the
   simplest topology. Five lines of fluent API.
2. **[`model_api/chemistry/reaction_system.py`](model_api/chemistry/reaction_system.py)** —
   what the builder hides for chemistry. Declare an aerobic-growth
   kinetic reaction, an acid-base equilibrium, and a cross-phase
   CO₂ partition; inspect the post-bucketing `ReactionSystem`.
3. **[`model_api/D2Cworkshop/basic_layout/raw_construction.py`](model_api/D2Cworkshop/basic_layout/raw_construction.py)** —
   what the builder hides for topology. Same batch fermenter as
   step 1, but built by hand from `GasPhase` + `LiquidPhase` +
   `KineticGasLiquidLink` + `ControlVolume` + `Simulation`, with
   the chemistry imported from step 2.
4. **[`builder/cstr_fermenter.py`](builder/cstr_fermenter.py)** —
   adds inlet/outlet boundaries and a DO controller alongside pH.
5. **[`builder/fed_batch_fermenter.py`](builder/fed_batch_fermenter.py)** —
   volume-grows-over-time wrinkle with a substrate-balance check.
6. **[`builder/microplate_fermenter.py`](builder/microplate_fermenter.py)** —
   chemistry stripped to a minimum, membrane gas exchange instead
   of sparging.
7. **[`model_api/chemistry/fba/fba_toy.py`](model_api/chemistry/fba/fba_toy.py)**
   and
   **[`model_api/chemistry/fba/fba_ecoli_core.py`](model_api/chemistry/fba/fba_ecoli_core.py)** —
   dynamic FBA wired through `BlackBoxReactionModel`. Both use
   original hand-crafted pedagogical networks (attribution notes in
   each file document scope and limitations); treat as protocol
   examples (how to wire an FBA solver), not as biochemistry
   references.

## Running

Run any demo from the repo root:

```bash
python demos/builder/batch_fermenter.py
python demos/model_api/chemistry/reaction_system.py
python demos/model_api/D2Cworkshop/basic_layout/raw_construction.py
# ... and so on
```

Each prints a short report. The builder demos report initial → final
concentrations, pH, controller diagnostics. The `model_api/`
demos vary by topic (the chemistry demo prints a bucket inspector;
the FBA demos print trajectory samples).

## Pattern in the new framework

The canonical shape across the builder demos:

```python
from vlmodels.fermenter.config import FermenterBuilder
from PyOMES.control.cv_loops import PHController

cv = (
    FermenterBuilder()
    .vessel(...)
    .transfer_kinetic(...)
    .chemistry(...)
    .organism(...)
    .substrate(...)
    .build()
)
cv.boundaries.append(PressureReliefVent(...))

sim = Simulation(cvs={"main": cv}, controllers=[PHController(...)])
result = sim.run(tau_h=5.0, n_steps=1000)
result.liquid_mol["main"]["AceticAcid"]   # ndarray
result.pH["main"]                          # ndarray
```

`FermenterBuilder.build_simulation_and_run(tau_h, n_steps)` is the
one-shot equivalent when you don't need to inspect or attach custom
boundaries — see [`batch_fermenter.py`](builder/batch_fermenter.py)
for that form. The `model_api/` demos unfold each of those builder
steps into explicit construction.

## Helpers

- **[_bootstrap.py](_bootstrap.py)** — adds `models/` to `sys.path`
  so `vlmodels` resolves without a separate install. Imported at
  the top of every demo via a one-line sys.path tweak that points
  each demo back at this folder.

## ADM1 / BSM2

The BSM2 reference model ships as a regression test rather than a
demo — see
[tests/standalone/test_bsm2_reference.py](../tests/standalone/test_bsm2_reference.py).
That test compares trajectory sentinels against captured
pyadm1-style reference values, which is the canonical "BSM2 vs
reference" comparison artifact.

## Restructuring history

This three-subtree layout shipped 2026-05-29 alongside the Q1 / Q2
demos under `model_api/`. The previous flat layout (four fermenter
demos at the top level, FBA under `demos/api/fba/`) is described
in [docs/dev/implementation/shipped/DEMO_RESTRUCTURE.md](../docs/dev/implementation/shipped/DEMO_RESTRUCTURE.md)
for context.
