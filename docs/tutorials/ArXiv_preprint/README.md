# ArXiv preprint tutorials

Curated, user-facing worked examples. Unlike [`demos/`](../../../demos/),
which is the full example set organised by framework layer (for
contributors exploring the API), this directory holds a smaller,
hand-picked selection intended for newcomers working through the
package end to end.

These three were migrated from [`demos/usecases/`](../../../demos/usecases/)
as the subset most closely covered by the project's ArXiv preprint.
The story continues in `demos/usecases/03_grow_ecoli_on_acetic_acid.ipynb`,
which stayed there.

## Notebooks

| Notebook | Situation | Uses |
|---|---|---|
| [01_predict_ph_simple_liquid.ipynb](01_predict_ph_simple_liquid.ipynb) | I'm making up a defined growth medium from KH₂PO₄ (phosphate buffer) and NH₄Cl (nitrogen source), no gas headspace or solid phase to track — what pH does that land at, across the range of doses used in practice? How does that compare against PHREEQC? | `NRChemicalEquilibriumEngine`, `PHREEQCChemicalEquilibriumEngine` (optional) |
| [02_kinetic_co2_equilibration_microplate_well.ipynb](02_kinetic_co2_equilibration_microplate_well.ipynb) | Pure water, in direct contact with a large atmospheric reservoir (O₂/N₂/CO₂) across a gas-liquid interface with a finite mass-transfer coefficient (kLa) rather than an instantaneous equilibrium — how does pH evolve over time as dissolved CO₂ approaches its Henry's-law equilibrium, and how does kLa itself set the timescale to get there? First notebook with genuinely kinetic (rate-limited) gas transfer. | `ControlVolume`, `Simulation`, `KineticTransferModel` |
| [03_cstr_dilution_rate_sweep.ipynb](03_cstr_dilution_rate_sweep.ipynb) | Same organism/substrate as `demos/usecases/03_grow_ecoli_on_acetic_acid.ipynb`, but now run as a chemostat — a CSTR fed and drained at the same volumetric flow rate, so the working volume holds steady while biomass and substrate settle onto a dilution-rate-dependent steady state. How does the dilution rate affect the reactor's volumetric productivity, and where does washout kick in? | `ControlVolume`, `Simulation`, `LiquidFeed`, `LiquidDrain` |

Launch from the repo root with `jupyter lab docs/tutorials/ArXiv_preprint/`.

## Regenerating

These notebooks are generated from
[`_generate_notebooks.py`](_generate_notebooks.py), split out from
[`demos/usecases/_generate_notebooks.py`](../../../demos/usecases/_generate_notebooks.py)
so this folder regenerates independently of the rest of `demos/usecases/`.
Edit the relevant section in this folder's script and re-run:

```bash
python docs/tutorials/ArXiv_preprint/_generate_notebooks.py
```
