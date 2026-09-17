# ArXiv preprint tutorials

Curated, user-facing worked examples: a smaller, hand-picked selection
intended for newcomers working through the package end to end, as opposed
to the full topic-organized tutorial set one level up.

These three were migrated from `demos/usecases/`, once the full example set
organised by framework layer (retired entirely 2026-09-17, its content
migrated into `docs/tutorials/` or removed), as the subset most closely
covered by the project's ArXiv preprint. A fourth notebook in that original
series, growing the same organism/substrate as a **batch** culture rather
than the chemostat below, was removed the same day as unmigrated,
low-value content.

## Notebooks

| Notebook | Situation | Uses |
|---|---|---|
| [01_predict_ph_simple_liquid.ipynb](01_predict_ph_simple_liquid.ipynb) | I'm making up a defined growth medium from KH₂PO₄ (phosphate buffer) and NH₄Cl (nitrogen source), no gas headspace or solid phase to track — what pH does that land at, across the range of doses used in practice? How does that compare against PHREEQC? | `NRChemicalEquilibriumEngine`, `PHREEQCChemicalEquilibriumEngine` (optional) |
| [02_kinetic_co2_equilibration_microplate_well.ipynb](02_kinetic_co2_equilibration_microplate_well.ipynb) | Pure water, in direct contact with a large atmospheric reservoir (O₂/N₂/CO₂) across a gas-liquid interface with a finite mass-transfer coefficient (kLa) rather than an instantaneous equilibrium — how does pH evolve over time as dissolved CO₂ approaches its Henry's-law equilibrium, and how does kLa itself set the timescale to get there? First notebook with genuinely kinetic (rate-limited) gas transfer. | `ControlVolume`, `Simulation`, `KineticTransferModel` |
| [03_cstr_dilution_rate_sweep.ipynb](03_cstr_dilution_rate_sweep.ipynb) | Inoculate the notebook 02 vessel with *E. coli* growing on acetic acid, run as a chemostat — a CSTR fed and drained at the same volumetric flow rate, so the working volume holds steady while biomass and substrate settle onto a dilution-rate-dependent steady state. How does the dilution rate affect the reactor's volumetric productivity, and where does washout kick in? | `ControlVolume`, `Simulation`, `LiquidFeed`, `LiquidDrain` |

Launch from the repo root with `jupyter lab docs/tutorials/ArXiv_preprint/`.

## Regenerating

These notebooks are generated from
[`_generate_notebooks.py`](_generate_notebooks.py), originally split out from
`demos/usecases/_generate_notebooks.py` (deleted along with the rest of
`demos/` on 2026-09-17) so this folder could regenerate independently of the
rest of that set. Edit the relevant section in this folder's script and
re-run:

```bash
python docs/tutorials/ArXiv_preprint/_generate_notebooks.py
```
