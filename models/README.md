# models/ (`vlmodels`)

Concrete bioprocess/wastewater model implementations built on top of
the core [`PyOMES`](../PyOMES/) framework. This directory is a
separate installable package (see [`pyproject.toml`](pyproject.toml)),
imported as `vlmodels`, so it can depend on `PyOMES` without `PyOMES`
itself depending on any specific model.

## Installing

From the repository root, with `PyOMES` already installed
(see the top-level [README.md](../README.md#installation)):

```bash
python -m pip install -e ./models
```

If you'd rather run scripts directly from a checkout without
installing `vlmodels`, add this directory to `PYTHONPATH` instead:

```bash
export PYTHONPATH="$PWD/models:$PYTHONPATH"
```

## Contents

| Package | Provides |
|---|---|
| [`adm1/`](vlmodels/adm1/) | Anaerobic Digestion Model No. 1 (ADM1) as a reusable reaction model (`build_adm1_reactions`, `build_adm1_cv`), plus a BSM2-canonical variant structurally consistent with the BSM2 Matlab/Simulink reference (Rosén & Jeppsson, 2006). |
| [`hplc/`](vlmodels/hplc/) | HPLC column model: 1D advection-dispersion with Langmuir adsorption, discretised into finite-volume cells and integrated with `scipy.integrate.solve_ivp`. |
| [`headspace.py`](vlmodels/headspace.py) | Ideal-gas headspace utility functions (total moles, mole fractions, headspace pressure) shared by other models. |

For worked examples that exercise these models end to end, see
[`docs/tutorials/`](../docs/tutorials/).
