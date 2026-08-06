# DARE2CYCLE Workshop Notebooks

Demonstration notebooks prepared for the DARE2CYCLE (D2C) project meeting.
Notebooks are organised into subfolders by purpose.

## Subfolders

| Folder | Contents |
|---|---|
| [`basic_layout/`](basic_layout/) | Foundational examples — raw construction of PyOMES models from primitives |
| [`updated_layout/`](updated_layout/) | Design-basis versions — explains *why* each parameter was chosen, not just *how* to use the API |

## `basic_layout/` notebooks

| File | Description |
|---|---|
| [`Example1_mtp_well.ipynb`](basic_layout/Example1_mtp_well.ipynb) | Closed-system MTP well — sealed headspace, equilibrium gas transfer, equilibrium chemistry loaded from `AQUEOUS_DEFAULT` database, Monod aerobic growth on acetic acid |
| [`Example2_batch_fermenter.ipynb`](basic_layout/Example2_batch_fermenter.ipynb) | Step-by-step assembly of a 0-D sparged batch fermenter — exposes every layer the builder hides: species, kinetics, acid-base equilibria, gas-liquid transfer models, and `ControlVolume` construction |
| [`Example3_CSTR.ipynb`](basic_layout/Example3_CSTR.ipynb) | Extends Example2 with `LiquidFeed` + `LiquidDrain` boundaries to model a CSTR at a specified dilution rate; verifies convergence to the analytical Monod steady state |

## How to run

Launch Jupyter from the repo root:

```
jupyter lab demos/model_api/D2Cworkshop/
```

Each notebook is self-contained: it locates the repo root automatically and adds
the `models/` directory to `sys.path`, so no installation step is required beyond
the standard PyOMES dependencies.

## `updated_layout/` notebooks

Each notebook mirrors one `basic_layout` example but leads with a pre-code design
basis section and closes with a validation section that checks explicit analytical
predictions against the simulation output. See
[`updated_layout/0_README.ipynb`](updated_layout/0_README.ipynb) for the reading-order
guide.

| File | Design-basis question answered |
|---|---|
| [`0_README.ipynb`](updated_layout/0_README.ipynb) | Overview and reading-order guide |
| [`Example1_mtp_well.ipynb`](updated_layout/Example1_mtp_well.ipynb) | Does the headspace contain enough O₂ to support the intended growth, or will O₂ run out before the substrate? |
| [`Example2_batch_fermenter.ipynb`](updated_layout/Example2_batch_fermenter.ipynb) | Does the chosen kLa supply enough oxygen to sustain growth without DO crashing? Is the pH setpoint achievable? |
| [`Example3_CSTR.ipynb`](updated_layout/Example3_CSTR.ipynb) | Where on the washout curve should the reactor operate? How does dilution rate set steady-state biomass and substrate conversion? |

## Key API points demonstrated

- `ControlVolume(transfer_models={...})` — canonical gas-liquid transfer via `KineticTransferModel` / `EquilibriumTransferModel`
- `ReactionBuilder.monod_aerobic_growth` — Monod kinetics and elemental-balance stoichiometry in one call
- `AQUEOUS_DEFAULT.reactions` — pulling standard equilibrium chemistry from the package database
- `Simulation.run()` → `BatchResult` — time-series output, summary statistics
