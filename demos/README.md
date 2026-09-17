# PyOMES demos

The examples that used to live here — the `StirredTankBuilder` quick-starts,
the direct-API chemistry/topology demos, the protocol deep-dives, and the
D2C Workshop notebooks — have moved into topic-organized subfolders under
[`docs/tutorials/`](../docs/tutorials/): `templates/`, `reactions/`,
`protocols/`, `D2C_workshop/`, and `results/`. See that folder's README for
the full index and what each one covers.

What's left here is content that hasn't been migrated:

## Prerequisite

Install the package in editable mode from the repo root:

```bash
pip install -e .
```

## Layout

```
demos/
  aerobic_fermentation_stoichiometry.ipynb   CHNOSP elemental-balance stoichiometry demo
  usecases/                                  scenario-first, five-minute worked examples
    0_README.ipynb                             notebook index (start here)
    02_kinetic_co2_equilibration_microplate_well.ipynb
    03_grow_ecoli_on_acetic_acid.ipynb
```

| File/subtree | What it's for |
|---|---|
| [`aerobic_fermentation_stoichiometry.ipynb`](aerobic_fermentation_stoichiometry.ipynb) | Derives all stoichiometric coefficients for aerobic growth on glucose from an elemental (CHNOSP) balance, comparing two published biomass elemental compositions (Roels extended, Upcraft). |
| [`usecases/`](usecases/) | You have a concrete scenario in mind ("I have a sample, I want to know X") and want a short, worked answer rather than a full API tour — e.g. [`usecases/03_grow_ecoli_on_acetic_acid.ipynb`](usecases/03_grow_ecoli_on_acetic_acid.ipynb) for growing *E. coli* on acetic acid in a batch bioreactor. A curated subset of `usecases/` closest to the project's ArXiv preprint lives in [`docs/tutorials/ArXiv_preprint/`](../docs/tutorials/ArXiv_preprint/) instead. |

## Running

Both are notebooks — launch from the repo root:

```bash
jupyter lab demos/
```

## ADM1 / BSM2

The BSM2 reference model ships as a regression test rather than a
demo — see
[tests/standalone/test_bsm2_reference.py](../tests/standalone/test_bsm2_reference.py).
That test compares trajectory sentinels against captured
pyadm1-style reference values, which is the canonical "BSM2 vs
reference" comparison artifact.

## Restructuring history

This folder used to hold `builder/`, `model_api/`, and `features/` subtrees
(organised by framework layer) alongside `usecases/` (organised by scenario).
Those three layer-organised subtrees were migrated into topic-organized
`docs/tutorials/` subfolders in the `tutorials-reorg` phase (see
[`docs/dev/implementation/upcoming/TUTORIALS_REORG_CHECKLIST.md`](../docs/dev/implementation/upcoming/TUTORIALS_REORG_CHECKLIST.md),
moved to `shipped/` once complete) — the eventual goal being to retire this
directory entirely once its remaining content also has a home elsewhere.
The layer-organised layout itself shipped 2026-05-29; the flat layout before
that (four fermenter demos at the top level, FBA under `demos/api/fba/`) is
described in [docs/dev/implementation/shipped/DEMO_RESTRUCTURE.md](../docs/dev/implementation/shipped/DEMO_RESTRUCTURE.md)
for context.
