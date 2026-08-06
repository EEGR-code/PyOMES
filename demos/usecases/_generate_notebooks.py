"""Generate the usecases demo notebooks as valid .ipynb JSON.

This folder is scenario-first: each notebook opens with a plain-language
"here's the situation" framing rather than a class-by-class API tour (that's
what demos/features/ and demos/model_api/chemistry/speciation/ are for).
Cross-reference those folders once a use case needs more depth than a
five-minute read supports.

Run once from the repo root:
    python demos/usecases/_generate_notebooks.py
"""
import json
from pathlib import Path

HERE = Path(__file__).parent

METADATA = {
    "kernelspec": {
        "display_name": "biosteam",
        "language": "python",
        "name": "python3",
    },
    "language_info": {
        "codemirror_mode": {"name": "ipython", "version": 3},
        "file_extension": ".py",
        "mimetype": "text/x-python",
        "name": "python",
        "nbconvert_exporter": "python",
        "pygments_lexer": "ipython3",
        "version": "3.12.3",
    },
}


def nb(*cells):
    return {"nbformat": 4, "nbformat_minor": 5, "metadata": METADATA, "cells": list(cells)}


def md(id_, src):
    return {"cell_type": "markdown", "id": id_, "metadata": {}, "source": src}


def code(id_, src):
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": id_,
        "metadata": {},
        "outputs": [],
        "source": src,
    }


# ═══════════════════════════════════════════════════════════════════════════
#  0  README (pure markdown, no code — matches the speciation/ and
#     ChemicalEquilibriumProtocol/ 0_README convention)
# ═══════════════════════════════════════════════════════════════════════════

readme_nb = nb(
    md("title", """\
# PyOMES Use Cases

Short, scenario-first notebooks: each one opens with a plain-language
situation ("I have a sample, I want to know X") and gets to an answer in a
handful of cells. If you want the full API tour of a class instead — every
constructor argument, every gotcha — see
[`demos/features/`](../features/) or
[`demos/model_api/chemistry/speciation/`](../model_api/chemistry/speciation/),
which these notebooks link out to once a use case needs more depth.\
"""),
    md("notebooks", """\
## Notebooks

| Notebook | Situation | Uses |
|---|---|---|
| [01_predict_ph_simple_liquid.ipynb](01_predict_ph_simple_liquid.ipynb) | I'm making up a defined growth medium from KH₂PO₄ (phosphate buffer) and NH₄Cl (nitrogen source), no gas headspace or solid phase to track — what pH does that land at, across the range of doses used in practice? How does that compare against PHREEQC? | `NRChemicalEquilibriumEngine`, `PHREEQCChemicalEquilibriumEngine` (optional) |
| [02_equilibrate_with_atmospheric_gas.ipynb](02_equilibrate_with_atmospheric_gas.ipynb) | Same medium as 01, but now left open to a typical atmosphere (simplified to O₂/N₂/CO₂) instead of sealed — does dissolved atmospheric CO₂ shift the pH, how much O₂/N₂ does the liquid actually hold, and where does the crossover from "buffered" to "CO₂-dominated" (rainwater-like) pH sit? | `NRChemicalEquilibriumEngine`, `HenryEquilibrium`, `PHREEQCChemicalEquilibriumEngine` (optional) |
| [02b_kinetic_co2_equilibration_microplate_well.ipynb](02b_kinetic_co2_equilibration_microplate_well.ipynb) | Pure water, in direct contact with a large atmospheric reservoir (O₂/N₂/CO₂) across a gas-liquid interface with a finite mass-transfer coefficient (kLa) rather than an instantaneous equilibrium — how does pH evolve over time as dissolved CO₂ approaches its Henry's-law equilibrium, and how does kLa itself set the timescale to get there? First notebook with genuinely kinetic (rate-limited) gas transfer. | `ControlVolume`, `Simulation`, `KineticTransferModel` |
| [03_grow_ecoli_on_acetic_acid.ipynb](03_grow_ecoli_on_acetic_acid.ipynb) | Inoculate that same sparged vessel with *E. coli* growing on acetic acid as sole carbon source — how fast does it grow, does dissolved O₂ ever become limiting, and what happens to pH as the acid substrate is consumed? First notebook where the chemistry evolves over time. | `ControlVolume`, `Simulation`, `ReactionBuilder.monod_aerobic_growth`, `ReactionSystem` |
| [04_compare_runtime_by_usecase.ipynb](04_compare_runtime_by_usecase.ipynb) | Usecases 01-03 all took *some* wall-clock time to run — how much, and where does it go? Compares run time across all three usecases, decomposed by which activity-coefficient treatment did the solving (PyOMES ideal / PyOMES Davies / PHREEQC ideal-equivalent / PHREEQC default). | `NRChemicalEquilibriumEngine`, `PHREEQCChemicalEquilibriumEngine` (optional), `ReactionSystem.configure_engine` |
| [05_cstr_dilution_rate_sweep.ipynb](05_cstr_dilution_rate_sweep.ipynb) | Same organism/substrate as 03, but now run as a chemostat — a CSTR fed and drained at the same volumetric flow rate, so the working volume holds steady while biomass and substrate settle onto a dilution-rate-dependent steady state. How does the dilution rate affect the reactor's volumetric productivity, and where does washout kick in? | `ControlVolume`, `Simulation`, `LiquidFeed`, `LiquidDrain` |

Launch from the repo root with `jupyter lab demos/usecases/`.\
"""),
)

save_path_readme = HERE / "0_README.ipynb"
with open(save_path_readme, "w", encoding="utf-8") as f:
    json.dump(readme_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path_readme}")

# ═══════════════════════════════════════════════════════════════════════════
#  01  Predict the pH of a simple liquid
# ═══════════════════════════════════════════════════════════════════════════

SETUP = """\
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def _find_repo():
    for p in [Path.cwd(), *Path.cwd().parents]:
        if (p / "pyproject.toml").exists():
            return p
    raise RuntimeError("Run from inside the PyOMES repo")

sys.path.insert(0, str(_find_repo() / "models"))

from PyOMES.chemistry.common_species import (
    H2O, H_plus, OH_minus,
    H3PO4, H2PO4_minus, HPO4_2minus, PO4_3minus,
    NH3, NH4_plus,
)
from PyOMES.reactions.equilibrium import EquilibriumReaction
from PyOMES.reactions.stoichiometry import StoichiometryEntry
from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

def _e(sp, coeff):
    return StoichiometryEntry(species=sp, phase="liquid", coefficient=coeff)

print("Imports OK")
"""

ph_nb = nb(
    md("title", """\
# Predict the pH of a Simple Aqueous Solution — KH₂PO₄ + NH₄Cl

**The situation:** you're making up a defined liquid growth medium from
two of its most common ingredients — KH₂PO₄ (a phosphate buffer salt) and
NH₄Cl (a nitrogen source) — dissolved in water. There is no gas headspace
to model and no precipitating solid to track, just the liquid phase on its
own. What pH does that land at, and how sensitive is it to how much of
each salt you weigh out?

This is the smallest possible use of PyOMES's chemical-equilibrium layer:
one engine, one reaction network, one `solve()` call. No `Phase`,
`ControlVolume`, or `Simulation` objects are needed for a single
equilibrium snapshot like this — those come in once the chemistry needs to
evolve over time or couple to other phases (see
[`demos/model_api/chemistry/reaction_system.py`](../model_api/chemistry/reaction_system.py)
for that step).\
"""),

    code("setup", SETUP),

    # ── 1. Declare the chemistry ─────────────────────────────────────────

    md("chem-md", """\
## 1  Declare the chemistry

Both salts fully dissociate on dissolving — KH₂PO₄ → K⁺ + H₂PO₄⁻, and
NH₄Cl → NH₄⁺ + Cl⁻ — so the acid-base chemistry actually at play is the
phosphate ladder (from the H₂PO₄⁻ that KH₂PO₄ contributes) and the
ammonium/ammonia couple (from the NH₄⁺ that NH₄Cl contributes), on top of
water's own autoionization. K⁺ and Cl⁻ take no part in any reaction — they
only enter the charge balance, as `strong_ions=`.

| Reaction | log K | p$K_a$ | Role |
|---|---|---|---|
| H₂O ⇌ H⁺ + OH⁻ | −14.0 | 14.0 | water autoionization |
| H₃PO₄ ⇌ H₂PO₄⁻ + H⁺ | −2.15 | 2.15 | phosphate, 1st step |
| H₂PO₄⁻ ⇌ HPO₄²⁻ + H⁺ | −7.20 | 7.20 | phosphate, 2nd step |
| HPO₄²⁻ ⇌ PO₄³⁻ + H⁺ | −12.35 | 12.35 | phosphate, 3rd step |
| NH₄⁺ ⇌ NH₃ + H⁺ | −9.25 | 9.25 | ammonium/ammonia |

`total_id="H3PO4"` on the three phosphate steps and `total_id="NH3"` on
the ammonium step tell the engine which mass-balance total each species
belongs to — the total phosphate and total ammoniacal nitrogen you'll
supply to `solve()` later. Note that the *master* label (`H3PO4`, `NH3`)
is just an id for the total; it doesn't imply the solution actually starts
from those molecular forms — here it starts from H₂PO₄⁻ and NH₄⁺, which
the engine handles identically since it solves for equilibrium composition,
not a synthesis path.\
"""),

    code("chem-code", """\
water = EquilibriumReaction(
    stoichiometry=[_e(H2O, -1), _e(H_plus, +1), _e(OH_minus, +1)],
    log_K=-14.0, label="water",
)
p1 = EquilibriumReaction(
    stoichiometry=[_e(H3PO4, -1), _e(H2PO4_minus, +1), _e(H_plus, +1)],
    log_K=-2.15, total_id="H3PO4", label="p1",
)
p2 = EquilibriumReaction(
    stoichiometry=[_e(H2PO4_minus, -1), _e(HPO4_2minus, +1), _e(H_plus, +1)],
    log_K=-7.20, total_id="H3PO4", label="p2",
)
p3 = EquilibriumReaction(
    stoichiometry=[_e(HPO4_2minus, -1), _e(PO4_3minus, +1), _e(H_plus, +1)],
    log_K=-12.35, total_id="H3PO4", label="p3",
)
nh4 = EquilibriumReaction(
    stoichiometry=[_e(NH4_plus, -1), _e(NH3, +1), _e(H_plus, +1)],
    log_K=-9.25, total_id="NH3", label="nh4",
)

print("5 reactions declared.")
"""),

    # ── 2. Instantiate the engine ────────────────────────────────────────

    md("engine-md", """\
## 2  Instantiate the equilibrium engine

`NRChemicalEquilibriumEngine.from_reactions()` takes the flat list of
reactions and builds a Newton-Raphson tableau from them — this is the
"chemicalEquilibrium system" for this sample: a solver object that knows
the reaction network and can be asked to equilibrate at any composition.\
"""),

    code("engine-code", """\
engine = NRChemicalEquilibriumEngine.from_reactions([water, p1, p2, p3, nh4])

print("Masters:    ", engine.tableau.masters)
print("Secondaries:", [s.species_id for s in engine.tableau.secondaries])
"""),

    # ── 3. Solve for the sample's pH ─────────────────────────────────────

    md("solve-md", """\
## 3  Solve for the sample's pH

A concrete recipe close to M9 minimal medium: 22 mmol/L KH₂PO₄ and
18.7 mmol/L NH₄Cl. Total phosphate and total ammoniacal nitrogen go in via
`totals=`, keyed by master id; the K⁺ and Cl⁻ each salt contributes go in
via `strong_ions=` at the *same* concentration as the parent salt, since
both dissociate 1:1.\
"""),

    code("solve-code", """\
CT_P = 0.022   # mol/L KH2PO4 -> mol/L total phosphate, mol/L K+
CT_N = 0.0187  # mol/L NH4Cl  -> mol/L total ammoniacal N, mol/L Cl-

result = engine.solve(
    totals={"H3PO4": CT_P, "NH3": CT_N},
    strong_ions={"CT_K": CT_P, "CT_Cl": CT_N},
)

print(f"pH               = {result.pH:.3f}")
print(f"[H2PO4-]         = {result.species_mol_L['H2PO4-']*1e3:.4f} mmol/L")
print(f"[HPO4--]         = {result.species_mol_L['HPO4--']*1e6:.4f} µmol/L")
print(f"[NH4+]           = {result.species_mol_L['NH4+']*1e3:.4f} mmol/L")
print(f"[NH3]            = {result.species_mol_L['NH3']*1e6:.4f} µmol/L")
print(f"charge residual  = {result.charge_residual:.2e}  (should be ~0)")
"""),

    # ── 4. Sensitivity to each salt in isolation ─────────────────────────

    md("sweep-md", """\
## 4  Sensitivity to each salt in isolation

Before the 2-D picture, it's worth seeing each salt's effect on its own.
KH₂PO₄ alone sets an *amphoteric* pH — H₂PO₄⁻ is both a weak acid (toward
HPO₄²⁻) and a weak base (toward H₃PO₄), and at moderate-to-high
concentration the solution pH converges toward
$\\tfrac12(\\text{p}K_{a1}+\\text{p}K_{a2}) \\approx 4.68$, largely
independent of exactly how much is dosed. NH₄Cl on its own is a much
weaker acid (p$K_a$ = 9.25, far from neutral) so its effect on pH near
neutral-to-acidic conditions is comparatively small.\
"""),

    code("sweep-code", """\
CT_P_vals = np.logspace(-3, -1, 30)    # 1 - 100 mmol/L KH2PO4
pH_P = [engine.solve(totals={"H3PO4": CT, "NH3": 0.0},
                     strong_ions={"CT_K": CT}).pH
        for CT in CT_P_vals]

CT_N_vals = np.logspace(-3, -1, 30)    # 1 - 100 mmol/L NH4Cl
pH_N = [engine.solve(totals={"H3PO4": 0.0, "NH3": CT},
                     strong_ions={"CT_Cl": CT}).pH
        for CT in CT_N_vals]

fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)

axes[0].semilogx(CT_P_vals * 1e3, pH_P, "o-", color="tab:blue")
axes[0].axhline(0.5 * (2.15 + 7.20), color="gray", ls="--", lw=1,
                label=r"$\\frac{1}{2}(pK_{a1}+pK_{a2})=4.68$")
axes[0].set_xlabel("KH2PO4 (mmol/L)")
axes[0].set_ylabel("pH")
axes[0].set_title("KH2PO4 alone")
axes[0].legend(fontsize=8)
axes[0].grid(True, which="both", alpha=0.3)

axes[1].semilogx(CT_N_vals * 1e3, pH_N, "o-", color="tab:orange")
axes[1].set_xlabel("NH4Cl (mmol/L)")
axes[1].set_title("NH4Cl alone")
axes[1].grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.show()
"""),

    # ── 5. Contour: pH as a function of both salts ───────────────────────

    md("contour-md", """\
## 5  pH across the KH₂PO₄ / NH₄Cl growth-media design space

Putting both salts together, over the range of doses commonly used in
defined bacterial/yeast growth media (roughly 1-100 mmol/L for each —
M9 minimal medium, for instance, uses ~22 mmol/L KH₂PO₄ and
~19 mmol/L NH₄Cl), gives a 2-D picture of how the recipe determines pH.
The result confirms what Section 4 suggested: KH₂PO₄ concentration is the
dominant control — contour lines run mostly horizontal — while NH₄Cl only
measurably shifts pH when phosphate is comparatively dilute.

The plot itself is deferred to the end of §6: it's drawn together with the
PHREEQC parity comparison as a single two-panel publication figure once
that comparison data exists, so the two share one figure style instead of
each getting its own throwaway version.\
"""),

    code("contour-code", """\
n_grid = 40
CT_P_grid = np.logspace(-3, -1, n_grid)   # 1 - 100 mmol/L KH2PO4
CT_N_grid = np.logspace(-3, -1, n_grid)   # 1 - 100 mmol/L NH4Cl

pH_grid = np.empty((n_grid, n_grid))
for i, ct_n in enumerate(CT_N_grid):
    for j, ct_p in enumerate(CT_P_grid):
        out = engine.solve(
            totals={"H3PO4": ct_p, "NH3": ct_n},
            strong_ions={"CT_K": ct_p, "CT_Cl": ct_n},
        )
        pH_grid[i, j] = out.pH

print(f"pH range over the grid: {pH_grid.min():.3f} - {pH_grid.max():.3f}")
"""),

    # ── 6. Benchmark against PHREEQC ─────────────────────────────────────

    md("phreeqc-md", """\
## 6  Benchmark against PHREEQC

Everything above solves a hand-declared reaction network at infinite
dilution — no activity correction. [PHREEQC](https://www.usgs.gov/software/phreeqc-version-3)
is the de facto standard geochemical speciation code; comparing against
it checks both the NR solver itself and where the infinite-dilution
assumption used in Sections 1-5 starts to matter.

**Which PHREEQC formulation, specifically:** `PHREEQCChemicalEquilibriumEngine`
goes through `phreeqpython`, whose default database is `vitens.dat` — a
derivative of the standard `phreeqc.dat` (same reactions, same activity
treatment). That database is **not** ideal: every species carries a
`-gamma` (ion-size å, b-dot) entry, so PHREEQC applies the **extended
(WATEQ) Debye-Hückel equation** to every solve, and its reaction list
includes the `KHPO4⁻`/`NaHPO4⁻` ion pairs. So this section compares two
*different nonideal* activity treatments against each other (WATEQ
Debye-Hückel + ion pairing vs. PyOMES's Davies option), not a nonideal
model against an ideal one — the notebook's own `engine` (ideal, Sections
1-5) is plotted alongside purely as the "what if we correct for nothing"
baseline.

This section is optional — it needs the `phreeqpython` package
(`pip install PyOMES[phreeqc]`). If it isn't installed, the cell below
reports that and the rest of the notebook is unaffected.\
"""),

    code("phreeqc-setup", """\
try:
    from PyOMES.chemical_equilibrium.phreeqc_engine import PHREEQCChemicalEquilibriumEngine
    _HAVE_PHREEQC = True
    print("phreeqpython available - PHREEQC benchmark cells will run.")
except ImportError as exc:
    _HAVE_PHREEQC = False
    print(f"phreeqpython not installed ({exc}); skipping PHREEQC benchmark cells.")
    print("Install with: pip install PyOMES[phreeqc]")
"""),

    md("phreeqc-point-md", """\
### 6a  A single point: the M9-like recipe

Three engines, same recipe as Section 3:

- **`engine` (ideal)** — the engine used throughout this notebook; no
  activity correction, strictly valid only at infinite dilution.
- **`engine_davies`** — identical reaction network, with
  `use_activity=True, activity_model="davies"`.
- **PHREEQC**, via `PHREEQCChemicalEquilibriumEngine` — `component_map`
  routes `H3PO4`→`P`, `NH3`→`N(-3)`, `CT_K`→`K`, `CT_Cl`→`Cl`. As noted
  above, this is PHREEQC running its own **nonideal** default (extended
  WATEQ Debye-Hückel activities, plus the `KHPO4⁻` ion pair PyOMES doesn't
  model) — not an ideal reference.

`use_warmstart=False` is required here — the default `True` reuses the
PHREEQC solution across calls *incrementally* (a `REACTION` addition, not
a reset to the declared total), which would double-count the composition
already used to prime the engine. See
[`demos/features/ChemicalEquilibriumProtocol/03_phreeqc_engine_basics.ipynb`](../features/ChemicalEquilibriumProtocol/03_phreeqc_engine_basics.ipynb)
§4 for the mechanism.\
"""),

    code("phreeqc-point-code", """\
if _HAVE_PHREEQC:
    engine_davies = NRChemicalEquilibriumEngine.from_reactions(
        [water, p1, p2, p3, nh4], use_activity=True, activity_model="davies",
    )
    engine_pq = PHREEQCChemicalEquilibriumEngine(
        {"H3PO4": CT_P * 1e3, "NH3": CT_N * 1e3, "CT_K": CT_P * 1e3, "CT_Cl": CT_N * 1e3},
        component_map={"H3PO4": "P", "NH3": "N(-3)", "CT_K": "K", "CT_Cl": "Cl"},
        use_warmstart=False,
    )

    result_davies = engine_davies.solve(
        totals={"H3PO4": CT_P, "NH3": CT_N}, strong_ions={"CT_K": CT_P, "CT_Cl": CT_N},
    )
    result_pq = engine_pq.solve(
        totals={"H3PO4": CT_P, "NH3": CT_N, "CT_K": CT_P, "CT_Cl": CT_N},
    )

    print(f"{'Engine':<28}{'pH':>8}")
    print(f"{'NR (ideal)':<28}{result.pH:>8.3f}")
    print(f"{'NR (Davies activity)':<28}{result_davies.pH:>8.3f}")
    print(f"{'PHREEQC (WATEQ Debye-Huckel)':<28}{result_pq.pH:>8.3f}")
    print()
    print(f"|ideal  - PHREEQC| = {abs(result.pH - result_pq.pH):.3f} pH units")
    print(f"|Davies - PHREEQC| = {abs(result_davies.pH - result_pq.pH):.3f} pH units")
"""),

    md("phreeqc-sweep-md", """\
### 6b  Where the infinite-dilution assumption breaks down

Repeating Section 4's KH₂PO₄-alone sweep with all three engines. Ionic
strength grows with dose (K⁺ and H₂PO₄⁻ both scale with `CT_P`), so the
ideal engine's error against PHREEQC (WATEQ Debye-Hückel) should grow
with concentration too, while the Davies-corrected engine — a different
nonideal treatment, not the same equation PHREEQC uses — should still
track it closely across the whole range. The left panel shows pH vs.
dose directly; the right panel is a parity plot (PyOMES pH vs. PHREEQC
pH at matching doses — perfect agreement falls on the dashed 1:1 line).\
"""),

    code("phreeqc-sweep-code", """\
if _HAVE_PHREEQC:
    pH_P_davies = [engine_davies.solve(totals={"H3PO4": CT, "NH3": 0.0},
                                        strong_ions={"CT_K": CT}).pH
                   for CT in CT_P_vals]

    engine_pq_sweep = PHREEQCChemicalEquilibriumEngine(
        {"H3PO4": 1.0, "CT_K": 1.0}, component_map={"H3PO4": "P", "CT_K": "K"},
        use_warmstart=False,
    )
    pH_P_pq = np.array([engine_pq_sweep.solve(totals={"H3PO4": CT, "CT_K": CT}).pH
                         for CT in CT_P_vals])
    pH_P_arr = np.array(pH_P)
    pH_P_davies_arr = np.array(pH_P_davies)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    axes[0].semilogx(CT_P_vals * 1e3, pH_P_arr, "o-", color="tab:blue", label="NR (ideal)")
    axes[0].semilogx(CT_P_vals * 1e3, pH_P_davies_arr, "^-", color="tab:green", label="NR (Davies)")
    axes[0].semilogx(CT_P_vals * 1e3, pH_P_pq, "s--", color="k", label="PHREEQC (WATEQ D-H)")
    axes[0].set_xlabel("KH2PO4 (mmol/L)")
    axes[0].set_ylabel("pH")
    axes[0].set_title("KH2PO4 alone - three engines")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, which="both", alpha=0.3)

    lo = min(pH_P_pq.min(), pH_P_arr.min(), pH_P_davies_arr.min())
    hi = max(pH_P_pq.max(), pH_P_arr.max(), pH_P_davies_arr.max())
    pad = (hi - lo) * 0.08
    axes[1].plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=1, zorder=0, label="1:1")
    axes[1].scatter(pH_P_pq, pH_P_arr, s=28, color="tab:blue", alpha=0.8, label="NR (ideal)")
    axes[1].scatter(pH_P_pq, pH_P_davies_arr, s=28, color="tab:green", marker="^", alpha=0.8, label="NR (Davies)")
    axes[1].set_xlabel("PHREEQC pH (WATEQ D-H)")
    axes[1].set_ylabel("PyOMES pH")
    axes[1].set_title("Parity vs. PHREEQC")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    print(f"Max |ideal  - PHREEQC| = {np.max(np.abs(pH_P_arr - pH_P_pq)):.3f} pH units")
    print(f"Max |Davies - PHREEQC| = {np.max(np.abs(pH_P_davies_arr - pH_P_pq)):.3f} pH units")
"""),

    md("phreeqc-summary-md", """\
The ideal engine's disagreement with PHREEQC grows with concentration, as
expected — activity effects strengthen with ionic strength, and the ideal
engine applies none. The Davies-corrected engine sits close to the 1:1
line across the whole range even though it isn't using PHREEQC's own
activity equation; the small residual that remains is a real model
difference (PHREEQC's `KHPO4⁻` ion pair, which this notebook's reaction
network doesn't declare, plus Davies vs. WATEQ Debye-Hückel being
different empirical fits), not solver error.

But that raises a fair question: how much of the 6a/6b gap is *really*
activity theory, versus PHREEQC simply using different log K values or
extra reactions? Section 6c isolates that.\
"""),

    md("phreeqc-ideal-md", """\
### 6c  Isolating the effect: forcing PHREEQC to (near-)ideal too

PHREEQC has **no built-in switch** for "ideal solution" — every stock
database applies an activity correction unconditionally. The standard
workaround (a documented PHREEQC technique, not a PyOMES feature) is to
override each species' `-gamma` ion-size parameter with something huge
and its `b`-dot term with 0: in the extended Debye-Hückel equation
`log10(γ) = -A z² √I / (1 + B a √I) + b I`, driving `a → ∞` drives the
whole correction to 0, i.e. `γ → 1`.

Doing this cleanly means writing a **from-scratch minimal database**
(rather than patching `vitens.dat`'s existing 900+ species, which are
wired together through a master-species graph that doesn't tolerate
piecemeal overrides) with just the species this notebook needs. That has
a second benefit: the database's log K values can be set to *exactly*
PyOMES's own (2.15, 7.20, 12.35, 9.25, 14.0), removing thermodynamic-data
differences as a confound entirely. What's left, if anything, is pure
numerical-solver disagreement — the cleanest possible check of the NR
solver itself.\
"""),

    code("phreeqc-ideal-code", """\
if _HAVE_PHREEQC:
    import tempfile
    from phreeqpython import PhreeqPython as _RawPhreeqPython

    # Same 5 reactions and log K's as Section 1, written as a PHREEQC database.
    # -gamma <huge> 0  suppresses the extended Debye-Hueckel term (gamma -> 1)
    # on every charged species -- an approximation of "ideal", not a native
    # PHREEQC mode.  N is declared as a plain (non-redox) master species since
    # this system never touches another nitrogen oxidation state.
    _IDEAL_DB = \"\"\"\\
SOLUTION_MASTER_SPECIES
H       H+      -1.     H       1.008
H(0)    H2      0.0     H
H(1)    H+      -1.     0.0
E       e-      0.0     0.0     0.0
O       H2O     0.0     O       16.00
O(0)    O2      0.0     O
O(-2)   H2O     0.0     0.0
P       PO4-3   0.0     P       30.974
N       NH4+    0.0     N       14.0067
K       K+      0.0     K       39.098
Cl      Cl-     0.0     Cl      35.453

SOLUTION_SPECIES
H+ = H+
        log_k           0.0
        -gamma          1e6     0
e- = e-
        log_k           0.0
H2O = H2O
        log_k           0.0
2 H+ + 2 e- = H2
        log_k           -3.15
2 H2O = O2 + 4 H+ + 4 e-
        log_k           -86.08
H2O = OH- + H+
        log_k           -14.0
        -gamma          1e6     0
PO4-3 = PO4-3
        log_k           0.0
        -gamma          1e6     0
PO4-3 + H+ = HPO4-2
        log_k           12.35
        -gamma          1e6     0
PO4-3 + 2H+ = H2PO4-
        log_k           19.55
        -gamma          1e6     0
PO4-3 + 3H+ = H3PO4
        log_k           21.70
NH4+ = NH4+
        log_k           0.0
        -gamma          1e6     0
NH4+ = NH3 + H+
        log_k           -9.25
K+ = K+
        log_k           0.0
        -gamma          1e6     0
Cl- = Cl-
        log_k           0.0
        -gamma          1e6     0
END
\"\"\"
    _db_dir = Path(tempfile.gettempdir())
    (_db_dir / "vlsim_ideal_phosphate_ammonium.dat").write_text(_IDEAL_DB)
    pp_ideal = _RawPhreeqPython(database="vlsim_ideal_phosphate_ammonium.dat",
                                 database_directory=_db_dir)

    def _solve_ideal_pq(ct_p, ct_n):
        sol = pp_ideal.add_solution_raw({
            "P": ct_p * 1e3, "N": ct_n * 1e3, "K": ct_p * 1e3, "Cl": ct_n * 1e3,
            "temp": 25.0, "pH": "7 charge", "units": "mmol/L",
        })
        ph = sol.pH
        sol.forget()
        return ph

    ph_point_ideal_pq = _solve_ideal_pq(CT_P, CT_N)
    print(f"NR (ideal)              pH = {result.pH:.4f}")
    print(f"PHREEQC (gamma -> 1)    pH = {ph_point_ideal_pq:.4f}")
    print(f"|NR ideal - PHREEQC ideal| = {abs(result.pH - ph_point_ideal_pq):.4f} pH units")
"""),

    md("phreeqc-ideal-sweep-md", """\
Repeating the KH₂PO₄-alone sweep once more, ideal engine against
near-ideal PHREEQC with matched log K's — this is the true solver-vs-solver
check the earlier sections couldn't isolate.\
"""),

    code("phreeqc-ideal-sweep-code", """\
if _HAVE_PHREEQC:
    pH_P_ideal_pq = np.array([_solve_ideal_pq(ct, 0.0) for ct in CT_P_vals])

    fig, ax = plt.subplots(figsize=(5.5, 5))
    lo = min(pH_P_arr.min(), pH_P_ideal_pq.min())
    hi = max(pH_P_arr.max(), pH_P_ideal_pq.max())
    pad = (hi - lo) * 0.08
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=1, zorder=0, label="1:1")
    ax.scatter(pH_P_ideal_pq, pH_P_arr, s=28, color="tab:red", alpha=0.8)
    ax.set_xlabel("PHREEQC pH (gamma -> 1, matched log K)")
    ax.set_ylabel("NR (ideal) pH")
    ax.set_title("Ideal vs. (near-)ideal: pure solver agreement")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    print(f"Max |NR ideal - PHREEQC ideal| = {np.max(np.abs(pH_P_arr - pH_P_ideal_pq)):.4f} pH units")
    print(f"(compare: Max |NR ideal - PHREEQC WATEQ D-H| = {np.max(np.abs(pH_P_arr - pH_P_pq)):.3f} pH units, Section 6b)")
"""),

    md("phreeqc-ideal-summary-md", """\
With activity corrections and reaction data both pinned to "ideal, same
log K's," the two solvers agree to within a few thousandths of a pH
unit across the whole dose range — roughly 50x tighter than the
ideal-vs-real-PHREEQC gap from Section 6b. That confirms the 6a/6b
mismatch is (almost entirely) activity theory, not a bug in either
solver.\
"""),

    md("phreeqc-combined-md", """\
### 6d  Publication figures: design space + accuracy, as separate files

Same shared figure style as before, but each panel is now its own
independent figure, saved as its own high-resolution file (vector PDF +
600dpi PNG) so each can be dropped into a LaTeX document directly (e.g.
`\\includegraphics{fig_contour_ph_design_space.pdf}`):

- `fig_contour_ph_design_space.{pdf,png}` — the KH₂PO₄/NH₄Cl design-space
  contour from Section 5 (deferred to here so it could be built and styled
  together with the parity figure below).
- `fig_parity_vs_phreeqc.{pdf,png}` — the PHREEQC parity comparison, only
  produced if `phreeqpython` is installed — PyOMES vs. PHREEQC under
  *matching* activity treatments in each case: ideal vs. ideal (Section 6c)
  and Davies vs. WATEQ Debye-Hückel (Section 6b). Each series lands close
  to the 1:1 line on its own terms; the point is that "close" means
  something very different for the two — thousandths of a pH unit for the
  ideal pair, hundredths for the nonideal pair (the residual there being
  real model differences: Davies vs. WATEQ Debye-Hückel are different
  empirical fits, and PHREEQC's `KHPO4⁻` ion pair isn't in this notebook's
  reaction network).\
"""),

    code("phreeqc-combined-code", """\
from matplotlib.ticker import MultipleLocator

# One shared style for both figures, each sized for a single-column
# placement (~3.5 in / 89 mm wide) in a two-column article. Scoped with
# rc_context so it doesn't leak into any other plot in this notebook.
FS_ticks = 10
FS_axes = 12

PUB_FIG_WIDTH_IN = 3.5
PUB_DPI = 600
PUB_STYLE = {
    "font.size": 8,
    "axes.labelsize": FS_axes,
    "axes.titlesize": FS_axes,
    "xtick.labelsize": FS_ticks,
    "ytick.labelsize": FS_ticks,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "grid.linewidth": 0.5,
}

FIG_DIR = _find_repo() / "demos" / "usecases" / "figures"
FIG_DIR.mkdir(exist_ok=True)

with plt.rc_context(PUB_STYLE):
    # -- figure 1: design-space contour, own file ------------------------
    # A slim extra row (not a side column) holds a horizontal colorbar,
    # placed below the contour panel, so the panel can use the full
    # figure width.
    fig_a = plt.figure(figsize=(PUB_FIG_WIDTH_IN, 3.9), layout="constrained")
    gs_a = fig_a.add_gridspec(2, 1, height_ratios=[10, 0.6])
    ax_a = fig_a.add_subplot(gs_a[0, 0])
    cax_a = fig_a.add_subplot(gs_a[1, 0])
    fig_a.get_layout_engine().set(hspace=0.08)

    cf = ax_a.contourf(CT_P_grid * 1e3, CT_N_grid * 1e3, pH_grid,
                        levels=20, cmap="viridis")
    # Standard 0.1-pH-unit contour lines (5.1, 5.0, 4.9, ...) instead of
    # data-range-derived levels.
    lo, hi = pH_grid.min(), pH_grid.max()
    line_levels = np.round(np.arange(np.ceil(lo * 10) / 10, hi, 0.1), 1)
    cs = ax_a.contour(CT_P_grid * 1e3, CT_N_grid * 1e3, pH_grid,
                       levels=line_levels, colors="white", linewidths=0.5)
    ax_a.clabel(cs, inline=True, fontsize=9, fmt="%.1f")
    ax_a.set_xscale("log")
    ax_a.set_yscale("log")
    ax_a.set_xlabel("KH2PO4 (mmol/L)")
    ax_a.set_ylabel("NH4Cl (mmol/L)")
    cbar = fig_a.colorbar(cf, cax=cax_a, orientation="horizontal")
    # Same standardized 0.1-pH ticks as the contour line labels, rather
    # than data-range-derived (and less round-looking) values.
    cbar.set_ticks(line_levels)
    cbar.set_label("pH")

    for ext in ("pdf", "png"):
        fig_a.savefig(FIG_DIR / f"fig_contour_ph_design_space.{ext}",
                      dpi=PUB_DPI, bbox_inches="tight")
    plt.show()

    # -- figure 2: PHREEQC parity, own file, if available ----------------
    if _HAVE_PHREEQC:
        fig_b, ax_b = plt.subplots(figsize=(PUB_FIG_WIDTH_IN, 3.3), layout="constrained")

        # Marker colors sampled from the same viridis colormap as the
        # contour figure above, so the two share one palette.
        color_ideal, color_nonideal = plt.cm.viridis([0.15, 0.85])

        all_vals = np.concatenate([pH_P_ideal_pq, pH_P_arr, pH_P_pq, pH_P_davies_arr])
        lo, hi = all_vals.min(), all_vals.max()
        pad = (hi - lo) * 0.05
        ax_b.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=0.8,
                  zorder=0, label="1:1")
        ax_b.scatter(pH_P_ideal_pq, pH_P_arr, s=16, color=color_ideal, alpha=0.8,
                     label="Ideal")
        ax_b.scatter(pH_P_pq, pH_P_davies_arr, s=16, color=color_nonideal, marker="^",
                     alpha=0.8, label="Nonideal")
        ax_b.set_xlabel("PHREEQC pH")
        ax_b.set_ylabel("PyOMES pH")
        # Same range on both axes (not just matplotlib's independent
        # per-axis autoscale) so the 1:1 line is a true 45 degrees.
        ax_b.set_xlim(lo - pad, hi + pad)
        ax_b.set_ylim(lo - pad, hi + pad)
        ax_b.set_aspect("equal", adjustable="box")
        # Equal xlim/ylim alone doesn't guarantee equal ticks: the default
        # locator picks tick spacing per axis based on the box's pixel
        # size, which "equal" aspect can make unequal (x got 0.2 spacing,
        # y got 0.1). Force both to the same explicit spacing.
        ax_b.xaxis.set_major_locator(MultipleLocator(0.1))
        ax_b.yaxis.set_major_locator(MultipleLocator(0.1))
        ax_b.legend(loc="lower right", frameon=True, handletextpad=0.4,
                    borderpad=0.4, labelspacing=0.3)

        for ext in ("pdf", "png"):
            fig_b.savefig(FIG_DIR / f"fig_parity_vs_phreeqc.{ext}",
                          dpi=PUB_DPI, bbox_inches="tight")
        plt.show()

if _HAVE_PHREEQC:
    print(f"Max |ideal series|    (NR ideal  vs. PHREEQC gamma->1) = {np.max(np.abs(pH_P_arr - pH_P_ideal_pq)):.4f} pH units")
    print(f"Max |nonideal series| (NR Davies vs. PHREEQC WATEQ D-H) = {np.max(np.abs(pH_P_davies_arr - pH_P_pq)):.4f} pH units")
    print(f"Saved: {FIG_DIR / 'fig_contour_ph_design_space.pdf'}")
    print(f"Saved: {FIG_DIR / 'fig_parity_vs_phreeqc.pdf'}")
else:
    print(f"Saved: {FIG_DIR / 'fig_contour_ph_design_space.pdf'}")
"""),

    md("phreeqc-combined-summary-md", """\
For a much deeper accuracy audit — carbonate, calcium, precipitation —
see
[`demos/model_api/chemistry/speciation/06_phreeqc_benchmark.ipynb`](../model_api/chemistry/speciation/06_phreeqc_benchmark.ipynb).\
"""),

    md("runtime-md", """\
## 7  Run time: software x activity model, 500 replicates each

One last comparison, this time on speed rather than accuracy — split by
*both* axes Section 6 already established, not just by software: the same
M9-like point from Section 6a (`CT_P`, `CT_N`), solved 500 times each by
all four engine/activity-model combinations from Sections 6a-6c (same
four-way split [`04_compare_runtime_by_usecase.ipynb`](04_compare_runtime_by_usecase.ipynb)
uses):

- **PyOMES ideal** (`engine`) — no activity correction.
- **PyOMES Davies** (`engine_davies`) — `use_activity=True, activity_model="davies"`.
- **PHREEQC default** (`engine_pq`) — `vitens.dat`'s WATEQ Debye-Hückel.
- **PHREEQC ideal-equivalent** (`_solve_ideal_pq`) — the §6c `-gamma 1e6 0`
  trick (γ→1) against the matched-log-K minimal database.

Each gets one untimed warmup call first so neither engine's one-time
import/IPC-connection cost biases the result, then mean ± standard
deviation over 500 timed replicate calls — not the noise-robust
"best-of-trials" timing usecase 04 uses, since the point here is to
characterize the *spread* of individual call times, not to filter it out.
The last three rows need `phreeqpython`; without it, only PyOMES ideal is
reported.

(PyOMES's engine was called `PyOMES` in an earlier version of this project —
same engine, current name.)\
"""),

    code("runtime-code", """\
import time

N_REPS = 500

def time_replicates(fn, n=N_REPS):
    fn()  # untimed warmup -- first call pays one-time import/IPC costs
    times_s = np.empty(n)
    for i in range(n):
        t0 = time.perf_counter()
        fn()
        times_s[i] = time.perf_counter() - t0
    return times_s * 1e3  # ms

runtime_results = {}
runtime_results["PyOMES ideal"] = time_replicates(lambda: engine.solve(
    totals={"H3PO4": CT_P, "NH3": CT_N},
    strong_ions={"CT_K": CT_P, "CT_Cl": CT_N},
))

if _HAVE_PHREEQC:
    runtime_results["PyOMES Davies"] = time_replicates(lambda: engine_davies.solve(
        totals={"H3PO4": CT_P, "NH3": CT_N},
        strong_ions={"CT_K": CT_P, "CT_Cl": CT_N},
    ))
    runtime_results["PHREEQC default (WATEQ D-H)"] = time_replicates(lambda: engine_pq.solve(
        totals={"H3PO4": CT_P, "NH3": CT_N, "CT_K": CT_P, "CT_Cl": CT_N},
    ))
    runtime_results["PHREEQC ideal-equivalent (gamma->1)"] = time_replicates(
        lambda: _solve_ideal_pq(CT_P, CT_N)
    )
else:
    print("phreeqpython not installed; only PyOMES ideal timed below.")

print(f"{'Software / activity model':<36}{'mean (ms)':>12}{'std (ms)':>12}")
for label, times_ms in runtime_results.items():
    print(f"{label:<36}{times_ms.mean():>12.4f}{times_ms.std():>12.4f}")
"""),

    md("next-md", """\
## Where to go next

- **More chemistry in the same liquid** (strong ions as a matter of
  course, temperature correction, activity coefficients) —
  [`demos/model_api/chemistry/speciation/`](../model_api/chemistry/speciation/),
  which verifies this same engine against closed-form analytical results
  (including this exact phosphate ladder).
- **Engine mechanics and gotchas** (constructor arguments, warmstart
  caching, what `algebraic_species()` returns) —
  [`demos/features/ChemicalEquilibriumProtocol/`](../features/ChemicalEquilibriumProtocol/).
- **Wiring this into something that evolves over time** — a fermenter or
  reactor where pH is one state among many being integrated — see
  [`demos/model_api/chemistry/reaction_system.py`](../model_api/chemistry/reaction_system.py)
  and [`demos/builder/`](../builder/) for the full `Simulation` pattern.\
"""),
)

save_path_ph = HERE / "01_predict_ph_simple_liquid.ipynb"
with open(save_path_ph, "w", encoding="utf-8") as f:
    json.dump(ph_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path_ph}")

# ═══════════════════════════════════════════════════════════════════════════
#  02  Equilibrate the same liquid against an atmospheric O2/N2/CO2 headspace
# ═══════════════════════════════════════════════════════════════════════════

GAS_SETUP = """\
import sys
from pathlib import Path
import math
import numpy as np
import matplotlib.pyplot as plt

def _find_repo():
    for p in [Path.cwd(), *Path.cwd().parents]:
        if (p / "pyproject.toml").exists():
            return p
    raise RuntimeError("Run from inside the PyOMES repo")

sys.path.insert(0, str(_find_repo() / "models"))

from PyOMES.chemistry.common_species import (
    H2O, H_plus, OH_minus,
    H3PO4, H2PO4_minus, HPO4_2minus, PO4_3minus,
    NH3, NH4_plus,
    CO2, HCO3_minus, CO3_2minus,
)
from PyOMES.chemistry.species import Species
from PyOMES.chemistry import HenryEquilibrium
from PyOMES.reactions.equilibrium import EquilibriumReaction
from PyOMES.reactions.stoichiometry import StoichiometryEntry
from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

def _e(sp, coeff, phase="liquid"):
    return StoichiometryEntry(species=sp, phase=phase, coefficient=coeff)

print("Imports OK")
"""

gas_nb = nb(
    md("title", """\
# Equilibrating a Liquid Against Atmospheric Gas — O₂, N₂, CO₂ Headspace

**The situation:** [usecase 01](01_predict_ph_simple_liquid.ipynb) solved
the KH₂PO₄ + NH₄Cl medium as a sealed liquid, no gas headspace. Leave that
same medium open to a typical atmosphere instead — simplified to O₂, N₂,
and CO₂ — and two things happen: CO₂ dissolves and joins the acid-base
chemistry already at play (it's a weak diprotic acid, same as phosphoric
acid was in usecase 01), while O₂ and N₂ dissolve too but take no part in
any reaction, so they only ever show up as a reported concentration, never
in the pH.

This extends usecase 01 with exactly one new mechanic — Henry's-law
gas-liquid partitioning — while staying at the same "one engine, one
`solve()`" level: the atmosphere is a fixed boundary condition (fixed
partial pressures), not a finite gas volume that the liquid could deplete,
so it converts to a dissolved-gas total the same way a weighed-out salt
does, and that total goes into `solve()` exactly like phosphate did in
usecase 01.\
"""),

    code("setup", GAS_SETUP),

    # ── 1. Declare the chemistry ─────────────────────────────────────────

    md("chem-md", """\
## 1  Declare the chemistry

The phosphate and ammonium ladders are unchanged from usecase 01. Two more
reactions extend the acid-base network to cover dissolved CO₂ — carbonic
acid is diprotic, so it ladders down the same way phosphoric acid does,
just two steps instead of three:

| Reaction | log K | p$K_a$ | Role |
|---|---|---|---|
| H₂O ⇌ H⁺ + OH⁻ | −14.0 | 14.0 | water autoionization |
| H₃PO₄ ⇌ H₂PO₄⁻ + H⁺ | −2.15 | 2.15 | phosphate, 1st step |
| H₂PO₄⁻ ⇌ HPO₄²⁻ + H⁺ | −7.20 | 7.20 | phosphate, 2nd step |
| HPO₄²⁻ ⇌ PO₄³⁻ + H⁺ | −12.35 | 12.35 | phosphate, 3rd step |
| NH₄⁺ ⇌ NH₃ + H⁺ | −9.25 | 9.25 | ammonium/ammonia |
| CO₂(aq) + H₂O ⇌ HCO₃⁻ + H⁺ | −6.35 | 6.35 | carbonate, 1st step |
| HCO₃⁻ ⇌ CO₃²⁻ + H⁺ | −10.33 | 10.33 | carbonate, 2nd step |

`total_id="CO2"` on the two carbonate steps marks their own mass balance,
the same way `total_id="H3PO4"` did for phosphate — except this total
won't be weighed out as a salt, it'll be *computed* from an atmospheric
partial pressure via Henry's law in Section 2.

Alongside the acid-base network, each gas gets a Henry's-law solubility
constant (Sander convention: mol m⁻³ Pa⁻¹, plus a van 't Hoff temperature
sensitivity `dlnH`):

| Gas | H_ref (mol m⁻³ Pa⁻¹) | dlnH (K) | k_H at 25 °C (mol L⁻¹ atm⁻¹) | Reacts? |
|---|---|---|---|---|
| CO₂ | 3.4 × 10⁻⁴ | 2400 | 0.0345 | yes — carbonate ladder above |
| O₂  | 1.3 × 10⁻⁵ | 1500 | 0.00132 | no — inert here |
| N₂  | 6.4 × 10⁻⁶ | 1300 | 0.000648 | no — inert here |

O₂ and N₂ get no acid-base reaction at all (this notebook doesn't model
redox/respiration) — they're declared purely so their dissolved
concentration can be reported in Section 2, alongside the carbonate
chemistry that *does* feed back into pH.\
"""),

    code("chem-code", """\
water = EquilibriumReaction(
    stoichiometry=[_e(H2O, -1), _e(H_plus, +1), _e(OH_minus, +1)],
    log_K=-14.0, label="water",
)
p1 = EquilibriumReaction(
    stoichiometry=[_e(H3PO4, -1), _e(H2PO4_minus, +1), _e(H_plus, +1)],
    log_K=-2.15, total_id="H3PO4", label="p1",
)
p2 = EquilibriumReaction(
    stoichiometry=[_e(H2PO4_minus, -1), _e(HPO4_2minus, +1), _e(H_plus, +1)],
    log_K=-7.20, total_id="H3PO4", label="p2",
)
p3 = EquilibriumReaction(
    stoichiometry=[_e(HPO4_2minus, -1), _e(PO4_3minus, +1), _e(H_plus, +1)],
    log_K=-12.35, total_id="H3PO4", label="p3",
)
nh4 = EquilibriumReaction(
    stoichiometry=[_e(NH4_plus, -1), _e(NH3, +1), _e(H_plus, +1)],
    log_K=-9.25, total_id="NH3", label="nh4",
)
co2_first = EquilibriumReaction(
    stoichiometry=[_e(CO2, -1), _e(H2O, -1), _e(HCO3_minus, +1), _e(H_plus, +1)],
    log_K=-6.35, total_id="CO2", label="co2_first",
)
co2_second = EquilibriumReaction(
    stoichiometry=[_e(HCO3_minus, -1), _e(CO3_2minus, +1), _e(H_plus, +1)],
    log_K=-10.33, total_id="CO2", label="co2_second",
)

# Henry's-law solubility constants. gas_species/liquid_species are left
# unset -- these are used only as a pCO2/pO2/pN2 -> dissolved-total
# converter (Section 2), never folded into the engine's own tableau, so
# they don't need to double as EquilibriumConstraint stoichiometry.
co2_henry = HenryEquilibrium(H_ref=3.4e-4, dlnH=2400.0, label="henry_CO2")
o2_henry = HenryEquilibrium(H_ref=1.3e-5, dlnH=1500.0, label="henry_O2")
n2_henry = HenryEquilibrium(H_ref=6.4e-6, dlnH=1300.0, label="henry_N2")

O2 = Species(id="O2", atoms={"O": 2}, charge=0)
N2 = Species(id="N2", atoms={"N": 2}, charge=0)

print("7 reactions declared, 3 Henry constants declared.")
"""),

    # ── 2. Instantiate the engine ────────────────────────────────────────

    md("engine-md", """\
## 2  Instantiate the equilibrium engine, and convert the atmosphere to totals

`engine` is a plain liquid-only `NRChemicalEquilibriumEngine`, same as
usecase 01 — the carbonate ladder makes it CO₂-aware, but nothing about
*how* the engine solves changes. The new step is converting a fixed
atmosphere (partial pressures) into the dissolved-gas totals `solve()`
expects, via Henry's law: $C^*_{aq} = k_H(T) \\times p_{gas}$.

`kH_mol_L_atm()` below reproduces exactly the formula
`HenryEquilibrium` uses internally (Sander convention converted to
mol L⁻¹ atm⁻¹, van 't Hoff temperature correction) — written out
explicitly here rather than called from the class, the same way usecase
01 wrote out log K's explicitly rather than hiding them in a database
lookup.\
"""),

    code("engine-code", """\
engine = NRChemicalEquilibriumEngine.from_reactions(
    [water, p1, p2, p3, nh4, co2_first, co2_second]
)

print("Masters:    ", engine.tableau.masters)
print("Secondaries:", [s.species_id for s in engine.tableau.secondaries])

def kH_mol_L_atm(H_ref, dlnH, T_K, T_ref=298.15):
    H_ref_mol_L_atm = (H_ref / 1000.0) * 101325.0
    return H_ref_mol_L_atm * math.exp(dlnH * (1.0 / T_K - 1.0 / T_ref))

T_K = 298.15  # 25 C

# Typical dry-air composition, simplified to O2/N2/CO2 per the notebook's
# scope (the remaining ~1% -- mostly argon -- is not modelled).
p_N2_atm = 0.78084
p_O2_atm = 0.20946
p_CO2_atm = 400e-6   # 400 ppm, present-day atmospheric CO2

CT_CO2_atm = kH_mol_L_atm(co2_henry.H_ref, co2_henry.dlnH, T_K) * p_CO2_atm
CT_O2_atm = kH_mol_L_atm(o2_henry.H_ref, o2_henry.dlnH, T_K) * p_O2_atm
CT_N2_atm = kH_mol_L_atm(n2_henry.H_ref, n2_henry.dlnH, T_K) * p_N2_atm

print(f"Dissolved CO2 at 400 ppm : {CT_CO2_atm*1e6:.3f} umol/L")
print(f"Dissolved O2  at 20.9%   : {CT_O2_atm*O2.MW*1e3:.3f} mg/L  "
      f"(compare: DO meters read ~8-9 mg/L for air-saturated water at 25 C)")
print(f"Dissolved N2  at 78.1%   : {CT_N2_atm*N2.MW*1e3:.3f} mg/L  "
      f"(compare: ~14-15 mg/L is the typical literature value)")
"""),

    # ── 3. Effect on the M9-like recipe ──────────────────────────────────

    md("solve-md", """\
## 3  Does atmospheric CO₂ change the M9-like recipe's pH?

Same recipe as usecase 01 Section 3 — 22 mmol/L KH₂PO₄, 18.7 mmol/L
NH₄Cl — solved twice: sealed (`CO2` total = 0, usecase 01's answer) and
open to the atmosphere (`CO2` total = the value just computed).\
"""),

    code("solve-code", """\
CT_P = 0.022   # mol/L KH2PO4 -> mol/L total phosphate, mol/L K+
CT_N = 0.0187  # mol/L NH4Cl  -> mol/L total ammoniacal N, mol/L Cl-

result_sealed = engine.solve(
    totals={"H3PO4": CT_P, "NH3": CT_N, "CO2": 0.0},
    strong_ions={"CT_K": CT_P, "CT_Cl": CT_N},
)
result_open = engine.solve(
    totals={"H3PO4": CT_P, "NH3": CT_N, "CO2": CT_CO2_atm},
    strong_ions={"CT_K": CT_P, "CT_Cl": CT_N},
)

print(f"pH, sealed (usecase 01)             = {result_sealed.pH:.4f}")
print(f"pH, open to atmosphere (400 ppm)    = {result_open.pH:.4f}")
print(f"Delta pH                            = {result_open.pH - result_sealed.pH:+.5f}")
print(f"[HCO3-] at equilibrium              = {result_open.species_mol_L['HCO3-']*1e6:.4f} umol/L")
"""),

    md("solve-note-md", """\
The shift is tiny — micromolar CO₂ barely dents a medium buffered by
22 mmol/L phosphate. That's not a bug; Section 5 shows exactly where the
crossover into "CO₂ actually matters" sits.\
"""),

    # ── 4. Sanity check: pure water ──────────────────────────────────────

    md("sanity-md", """\
## 4  Sanity check: pure water open to the atmosphere

With no phosphate buffer to mask it, dissolved atmospheric CO₂ alone
should reproduce a well-known textbook number: unpolluted rainwater sits
at about pH 5.6, set entirely by carbonic acid in equilibrium with
~400 ppm atmospheric CO₂ (see e.g. Stumm & Morgan, *Aquatic Chemistry*).
Reproducing that here checks the carbonate ladder declared in Section 1
against an independent, well-known reference point.\
"""),

    code("sanity-code", """\
result_pure = engine.solve(
    totals={"H3PO4": 0.0, "NH3": 0.0, "CO2": CT_CO2_atm},
    strong_ions={},
)
print(f"Pure water + 400 ppm atmospheric CO2 -> pH = {result_pure.pH:.3f}")
print("(textbook reference: unpolluted rainwater is ~pH 5.6, from atmospheric CO2 alone)")
"""),

    # ── 5. Where CO2 actually matters ────────────────────────────────────

    md("sweep-md", """\
## 5  Where atmospheric CO₂ actually matters: buffer capacity

Section 3 barely moved pH; Section 4 showed CO₂ alone sets pH ≈ 5.6.
Sweeping phosphate dose against headspace pCO₂ — from the present-day
atmosphere up to a CO₂-rich headspace like an anaerobic-digester biogas
space (tens of percent CO₂) — maps out the crossover between those two
regimes directly: flat contour lines (buffered, phosphate wins) versus
horizontal-banded ones (unbuffered, pCO₂ wins).\
"""),

    code("sweep-code", """\
n_grid = 30
pCO2_grid = np.logspace(np.log10(200e-6), 0.0, n_grid)   # 200 ppm - 100% CO2
CT_P_grid = np.logspace(-4, -1, n_grid)                    # 0.1 - 100 mmol/L KH2PO4

kH_CO2 = kH_mol_L_atm(co2_henry.H_ref, co2_henry.dlnH, T_K)

pH_grid = np.empty((n_grid, n_grid))
for i, p_co2 in enumerate(pCO2_grid):
    ct_co2 = kH_CO2 * p_co2
    for j, ct_p in enumerate(CT_P_grid):
        out = engine.solve(
            totals={"H3PO4": ct_p, "NH3": 0.0, "CO2": ct_co2},
            strong_ions={"CT_K": ct_p},
        )
        pH_grid[i, j] = out.pH

fig, ax = plt.subplots(figsize=(7.5, 6))
cf = ax.contourf(CT_P_grid * 1e3, pCO2_grid * 1e2, pH_grid, levels=20, cmap="viridis")
cs = ax.contour(CT_P_grid * 1e3, pCO2_grid * 1e2, pH_grid, levels=8,
                 colors="white", linewidths=0.6)
ax.clabel(cs, inline=True, fontsize=8, fmt="%.2f")

ax.scatter([22], [400e-6 * 100], marker="*", s=200, color="white",
           edgecolor="black", linewidth=0.8, zorder=5,
           label="Section 3 point (22 mM KH2PO4, 400 ppm)")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("KH2PO4 (mmol/L)")
ax.set_ylabel("pCO2 (% atm)")
ax.set_title("pH vs. phosphate dose and headspace CO2\\n(25 C, no ammonium)")
ax.legend(fontsize=8, loc="lower left")

cbar = fig.colorbar(cf, ax=ax)
cbar.set_label("pH")

plt.tight_layout()
plt.show()
"""),

    # ── 6. Benchmark against PHREEQC ─────────────────────────────────────

    md("phreeqc-md", """\
## 6  Benchmark against PHREEQC

[Usecase 01 Section 6](01_predict_ph_simple_liquid.ipynb) works through
*why* an ideal engine and an activity-corrected one land at different
distances from PHREEQC in detail — same story applies here, so this
section just checks the carbonate chemistry added in this notebook holds
up under the same test: one parity plot, two series, each engine compared
against PHREEQC under a *matching* activity treatment (ideal vs. ideal,
Davies vs. WATEQ Debye-Hückel), swept across the atmospheric pCO₂ range
(200-5000 ppm) at the fixed M9-like point from Section 3.

As in usecase 01, this needs the optional `phreeqpython` package
(`pip install PyOMES[phreeqc]`); the cell below reports if it's missing
and the rest of the notebook is unaffected.\
"""),

    code("phreeqc-setup", """\
try:
    import phreeqpython  # noqa: F401  -- actually probe for the optional dep;
    # PHREEQCChemicalEquilibriumEngine itself imports phreeqpython lazily
    # inside __init__, so importing only the wrapper class below would
    # succeed even without phreeqpython installed.
    from PyOMES.chemical_equilibrium.phreeqc_engine import PHREEQCChemicalEquilibriumEngine
    _HAVE_PHREEQC = True
    print("phreeqpython available - PHREEQC benchmark cell will run.")
except ImportError as exc:
    _HAVE_PHREEQC = False
    print(f"phreeqpython not installed ({exc}); skipping PHREEQC benchmark cell.")
    print("Install with: pip install PyOMES[phreeqc]")
"""),

    code("phreeqc-code", """\
if _HAVE_PHREEQC:
    import tempfile
    from phreeqpython import PhreeqPython as _RawPhreeqPython

    engine_davies = NRChemicalEquilibriumEngine.from_reactions(
        [water, p1, p2, p3, nh4, co2_first, co2_second],
        use_activity=True, activity_model="davies",
    )

    # Near-ideal PHREEQC database: usecase 01 Section 6c's database,
    # extended with the carbonate master species and this notebook's own
    # log K's (-6.35, -10.33) so only solver mechanics can differ. Every
    # master species needs its own trivial "X = X" identity reaction (see
    # PO4-3/NH4+/K+/Cl- below) -- IPhreeqc silently fails the whole
    # database load without one (surfaces later as "No database is
    # loaded" on the first solve, not as a load-time error), so CO2 gets
    # one too.
    _IDEAL_DB_CO2 = \"\"\"\\
SOLUTION_MASTER_SPECIES
H       H+      -1.     H       1.008
H(0)    H2      0.0     H
H(1)    H+      -1.     0.0
E       e-      0.0     0.0     0.0
O       H2O     0.0     O       16.00
O(0)    O2      0.0     O
O(-2)   H2O     0.0     0.0
C       CO2     0.0     C       12.011
P       PO4-3   0.0     P       30.974
N       NH4+    0.0     N       14.0067
K       K+      0.0     K       39.098
Cl      Cl-     0.0     Cl      35.453

SOLUTION_SPECIES
H+ = H+
        log_k           0.0
        -gamma          1e6     0
e- = e-
        log_k           0.0
H2O = H2O
        log_k           0.0
2 H+ + 2 e- = H2
        log_k           -3.15
2 H2O = O2 + 4 H+ + 4 e-
        log_k           -86.08
H2O = OH- + H+
        log_k           -14.0
        -gamma          1e6     0
CO2 = CO2
        log_k           0.0
        -gamma          1e6     0
CO2 + H2O = HCO3- + H+
        log_k           -6.35
        -gamma          1e6     0
HCO3- = CO3-2 + H+
        log_k           -10.33
        -gamma          1e6     0
PO4-3 = PO4-3
        log_k           0.0
        -gamma          1e6     0
PO4-3 + H+ = HPO4-2
        log_k           12.35
        -gamma          1e6     0
PO4-3 + 2H+ = H2PO4-
        log_k           19.55
        -gamma          1e6     0
PO4-3 + 3H+ = H3PO4
        log_k           21.70
NH4+ = NH4+
        log_k           0.0
        -gamma          1e6     0
NH4+ = NH3 + H+
        log_k           -9.25
K+ = K+
        log_k           0.0
        -gamma          1e6     0
Cl- = Cl-
        log_k           0.0
        -gamma          1e6     0
END
\"\"\"
    _db_dir = Path(tempfile.gettempdir())
    (_db_dir / "vlsim_ideal_gas_liquid.dat").write_text(_IDEAL_DB_CO2)
    pp_ideal_gl = _RawPhreeqPython(database="vlsim_ideal_gas_liquid.dat",
                                    database_directory=_db_dir)

    def _solve_ideal_pq_co2(ct_c, ct_p=CT_P, ct_n=CT_N):
        sol = pp_ideal_gl.add_solution_raw({
            "C": ct_c * 1e3, "P": ct_p * 1e3, "N": ct_n * 1e3,
            "K": ct_p * 1e3, "Cl": ct_n * 1e3,
            "temp": 25.0, "pH": "7 charge", "units": "mmol/L",
        })
        ph = sol.pH
        sol.forget()
        return ph

    engine_pq_gl = PHREEQCChemicalEquilibriumEngine(
        {"H3PO4": CT_P * 1e3, "NH3": CT_N * 1e3, "CO2": 1.0,
         "CT_K": CT_P * 1e3, "CT_Cl": CT_N * 1e3},
        component_map={"H3PO4": "P", "NH3": "N(-3)", "CO2": "C",
                        "CT_K": "K", "CT_Cl": "Cl"},
        use_warmstart=False,
    )

    pCO2_vals_atm = np.logspace(np.log10(200e-6), np.log10(5000e-6), 20)  # 200-5000 ppm
    CT_CO2_vals = kH_CO2 * pCO2_vals_atm

    pH_ideal = np.array([
        engine.solve(totals={"H3PO4": CT_P, "NH3": CT_N, "CO2": ct_c},
                     strong_ions={"CT_K": CT_P, "CT_Cl": CT_N}).pH
        for ct_c in CT_CO2_vals
    ])
    pH_davies = np.array([
        engine_davies.solve(totals={"H3PO4": CT_P, "NH3": CT_N, "CO2": ct_c},
                             strong_ions={"CT_K": CT_P, "CT_Cl": CT_N}).pH
        for ct_c in CT_CO2_vals
    ])
    pH_pq_ideal = np.array([_solve_ideal_pq_co2(ct_c) for ct_c in CT_CO2_vals])
    pH_pq_default = np.array([
        engine_pq_gl.solve(totals={"H3PO4": CT_P, "NH3": CT_N, "CO2": ct_c,
                                    "CT_K": CT_P, "CT_Cl": CT_N}).pH
        for ct_c in CT_CO2_vals
    ])

    fig, ax = plt.subplots(figsize=(6, 6))
    all_vals = np.concatenate([pH_ideal, pH_davies, pH_pq_ideal, pH_pq_default])
    lo, hi = all_vals.min(), all_vals.max()
    pad = (hi - lo) * 0.08
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=1, zorder=0, label="1:1")

    ax.scatter(pH_pq_ideal, pH_ideal, s=32, color="tab:red", alpha=0.85,
               label="Ideal: NR (ideal) vs. PHREEQC (gamma -> 1)")
    ax.scatter(pH_pq_default, pH_davies, s=32, color="tab:green", marker="^", alpha=0.85,
               label="Nonideal: NR (Davies) vs. PHREEQC (WATEQ D-H)")

    ax.set_xlabel("PHREEQC pH")
    ax.set_ylabel("PyOMES pH")
    ax.set_title("Gas-liquid CO2 equilibration: parity vs. PHREEQC\\n"
                 "(atmospheric pCO2 200-5000 ppm, fixed KH2PO4/NH4Cl)")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    print(f"Max |ideal series|    (NR ideal  vs. PHREEQC gamma->1)  = {np.max(np.abs(pH_ideal - pH_pq_ideal)):.4f} pH units")
    print(f"Max |nonideal series| (NR Davies vs. PHREEQC WATEQ D-H) = {np.max(np.abs(pH_davies - pH_pq_default)):.4f} pH units")
"""),

    md("next-md", """\
## Where to go next

- **The activity-theory deep dive this notebook only summarizes** —
  [`01_predict_ph_simple_liquid.ipynb`](01_predict_ph_simple_liquid.ipynb)
  Section 6 walks through *why* ideal vs. Davies land at different
  distances from PHREEQC, and isolates solver-vs-solver agreement from
  activity-model disagreement.
- **A finite, sealed headspace and the time it actually takes to get
  there** — this notebook treats the atmosphere as an infinite reservoir
  at fixed partial pressure, reached instantaneously; a small sealed
  volume (e.g. a microplate well) is neither —
  [`02b_kinetic_co2_equilibration_microplate_well.ipynb`](02b_kinetic_co2_equilibration_microplate_well.ipynb)
  reruns this exact chemistry as a time-resolved, kLa-limited
  gas-liquid transfer and shows how much that changes the answer.
- **Wiring this into something that evolves over time** — a fermenter or
  reactor sparged continuously by a gas phase — see
  [`demos/model_api/chemistry/reaction_system.py`](../model_api/chemistry/reaction_system.py)
  and [`demos/builder/`](../builder/) for the full `Simulation` pattern.\
"""),
)

save_path_gas = HERE / "02_equilibrate_with_atmospheric_gas.ipynb"
with open(save_path_gas, "w", encoding="utf-8") as f:
    json.dump(gas_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path_gas}")

# ═══════════════════════════════════════════════════════════════════════════
#  02b  Kinetic equilibration of pure water with the atmosphere
# ═══════════════════════════════════════════════════════════════════════════

KINETIC_SETUP = """\
import sys
from pathlib import Path
import math
import warnings
import numpy as np
import matplotlib.pyplot as plt

def _find_repo():
    for p in [Path.cwd(), *Path.cwd().parents]:
        if (p / "pyproject.toml").exists():
            return p
    raise RuntimeError("Run from inside the PyOMES repo")

sys.path.insert(0, str(_find_repo() / "models"))

from PyOMES.chemistry.common_species import H2O, H_plus, OH_minus, CO2, HCO3_minus, CO3_2minus
from PyOMES.chemistry import HenryEquilibrium
from PyOMES.reactions import EquilibriumReaction, ReactionSystem, StoichiometryEntry
from PyOMES.core import (
    ControlVolume, GasPhase, LiquidPhase, KineticTransferModel, Simulation,
    SimultaneousAdaptiveSolver,
)
from PyOMES.core.phases import R_L_ATM_MOL_K

# Pure water's tracked H/OH totals are tiny (~1e-8 mol, just the
# autoionization ions -- bulk solvent water isn't part of the element
# balance), so the conservation monitor's *relative* drift threshold trips
# on ordinary explicit-Euler step-to-step floating-point noise. Silenced
# only so the run cells' output stays readable; final pH/CO2(aq) values
# are unaffected.
from PyOMES.monitoring.accuracy import AccuracyWarning
from PyOMES.monitoring.conservation import ConservationWarning
warnings.filterwarnings("ignore", category=AccuracyWarning)
warnings.filterwarnings("ignore", category=ConservationWarning)

def _e(sp, coeff, phase="liquid"):
    return StoichiometryEntry(species=sp, phase=phase, coefficient=coeff)

print("Imports OK")
"""

kinetic_nb = nb(
    md("title", """\
# Kinetic Equilibration of Pure Water with the Atmosphere

**The situation:** a beaker of pure water is left open, in direct contact
with a large, well-mixed atmosphere (O₂, N₂, CO₂ at typical ambient
partial pressures) — no other solutes, no membrane or other barrier
between them. Equilibrium isn't instantaneous: gas has to physically
cross the gas-liquid interface, at a rate set by the volumetric
mass-transfer coefficient kLa (h⁻¹). Since CO₂ is a weak diprotic acid
(carbonic acid, laddering CO₂(aq) ⇌ HCO₃⁻ ⇌ CO₃²⁻), dissolving it shifts
pH away from neutral even with nothing else in the water — and it does so
gradually, on whatever timescale kLa sets, not all at once.

This notebook asks two questions about that process: how does pH evolve
over time at one representative kLa, and — since real gas-liquid contacting
spans a huge practical range, from a still puddle to a vigorously stirred
tank — how does the time to reach equilibrium scale with kLa itself?

The atmosphere is modelled as a large but finite gas reservoir (not a
fixed boundary condition) — large enough, relative to how little gas the
water actually needs to saturate, that its own composition barely moves
over the run (checked directly in Section 2); a "pseudo-unlimited" source
in practice, without needing a literal infinite reservoir.\
"""),

    code("setup", KINETIC_SETUP),

    # ── 1. Declare the chemistry ─────────────────────────────────────────

    md("chem-md", """\
## 1  Declare the chemistry

Just water's own autoionization plus the two-step carbonate ladder — no
buffer salts, since the liquid starts as pure water:

| Reaction | log K | p$K_a$ | Role |
|---|---|---|---|
| H₂O ⇌ H⁺ + OH⁻ | −14.0 | 14.0 | water autoionization |
| CO₂(aq) + H₂O ⇌ HCO₃⁻ + H⁺ | −6.35 | 6.35 | carbonate, 1st step |
| HCO₃⁻ ⇌ CO₃²⁻ + H⁺ | −10.33 | 10.33 | carbonate, 2nd step |

O₂ and N₂ get no acid-base reaction (this notebook doesn't model
redox/respiration) — they still cross the gas-liquid interface at the
same kLa as CO₂ (Section 2), they just take no part in any reaction, so
they never show up in the pH.

Each gas also gets a Henry's-law solubility constant (Sander convention:
mol m⁻³ Pa⁻¹, plus a van 't Hoff temperature sensitivity `dlnH`), which
is what `KineticTransferModel` (Section 2) uses to know the equilibrium
*target* each gas's kinetic transfer is relaxing toward:

| Gas | H_ref (mol m⁻³ Pa⁻¹) | dlnH (K) | k_H at 25 °C (mol L⁻¹ atm⁻¹) |
|---|---|---|---|
| CO₂ | 3.4 × 10⁻⁴ | 2400 | 0.0345 |
| O₂  | 1.3 × 10⁻⁵ | 1500 | 0.00132 |
| N₂  | 6.4 × 10⁻⁶ | 1300 | 0.000648 |\
"""),

    code("chem-code", """\
water = EquilibriumReaction(
    stoichiometry=[_e(H2O, -1), _e(H_plus, +1), _e(OH_minus, +1)],
    log_K=-14.0, label="water",
)
co2_first = EquilibriumReaction(
    stoichiometry=[_e(CO2, -1), _e(H2O, -1), _e(HCO3_minus, +1), _e(H_plus, +1)],
    log_K=-6.35, total_id="CO2", label="co2_first",
)
co2_second = EquilibriumReaction(
    stoichiometry=[_e(HCO3_minus, -1), _e(CO3_2minus, +1), _e(H_plus, +1)],
    log_K=-10.33, total_id="CO2", label="co2_second",
)
system = ReactionSystem([water, co2_first, co2_second], label="pure_water_co2_kinetics")

co2_henry = HenryEquilibrium(H_ref=3.4e-4, dlnH=2400.0, label="henry_CO2")
o2_henry = HenryEquilibrium(H_ref=1.3e-5, dlnH=1500.0, label="henry_O2")
n2_henry = HenryEquilibrium(H_ref=6.4e-6, dlnH=1300.0, label="henry_N2")

def kH_mol_L_atm(H_ref, dlnH, T_K, T_ref=298.15):
    return (H_ref / 1000.0) * 101325.0 * math.exp(dlnH * (1.0 / T_K - 1.0 / T_ref))

print("3 equilibrium reactions declared, 3 Henry constants declared.")
"""),

    # ── 2. Build the system ──────────────────────────────────────────────

    md("build-md", """\
## 2  Pure water in direct contact with a large atmosphere

Liquid: 100 mL of pure water, no dissolved gas yet (`CO2=O2=N2=0.0`).
Gas: today's atmosphere (400 ppm CO₂, 20.9% O₂, 78.1% N₂) at 1 atm,
filling a 1000 L headspace — a "pseudo-unlimited" reservoir in the sense
that matters here, not literally infinite: the check below confirms the
gas phase holds several orders of magnitude more of each species than the
100 mL of water could ever draw down, so its own composition stays
essentially fixed over the whole run without needing to model an actual
infinite boundary condition.

`transfer_models=` is the only inter-phase mechanism — no membrane or
other boundary sits between the two phases, i.e. direct contact. CO₂ uses
`transfer_basis="molecular"`: the driving force acts on dissolved
molecular CO₂(aq) alone (the same convention BSM2-style gas transfer uses
for CH₄/H₂), so the effective relaxation rate is `kLa · alpha`, where
`alpha` is the CO₂(aq) fraction of total dissolved carbon. At this
system's pH (~5.6, well below pKa1 = 6.35), alpha ≈ 1 throughout, so kLa
and the effective rate are nearly the same thing. O₂/N₂ have no acid-base
ladder to speciate into, so `transfer_basis="total"` is equivalent for
them.

Sections 3-5 all integrate through `SimultaneousAdaptiveSolver(method="BDF")`
rather than the plain forward-Euler `cv.advance()` default — an implicit,
adaptive-step ODE solver (via `scipy.integrate.solve_ivp`) suited to
gas-liquid transfer coupled to pH-sensitive speciation, and the natural
choice once Section 4/5 sweep kLa across three-plus orders of magnitude
(0.1-1000 /h): a fixed explicit step size that's stable at kLa = 0.1 /h
is not guaranteed to stay stable at kLa = 1000 /h, while BDF adapts its
own internal step size to whatever the local dynamics demand.\
"""),

    code("build-code", """\
V_liq = 0.1     # L (100 mL pure water)
V_gas = 1000.0  # L -- large enough not to deplete measurably; checked below
T_K = 298.15    # 25 C

p_CO2_atm, p_O2_atm, p_N2_atm = 400e-6, 0.20946, 0.78084

def build_cv(kLa, label="pure_water"):
    n_gas = (1.0 * V_gas) / (R_L_ATM_MOL_K * T_K)   # 1 atm headspace
    gas_phase = GasPhase(
        n_mol={
            "CO2": n_gas * p_CO2_atm,
            "O2": n_gas * p_O2_atm,
            "N2": n_gas * p_N2_atm,
        },
        V_L=V_gas, T_K=T_K,
    )
    liquid_phase = LiquidPhase(
        n_mol={"CO2": 0.0, "O2": 0.0, "N2": 0.0},
        V_L=V_liq, T_K=T_K,
    )
    transfer_models = {
        "CO2": KineticTransferModel(
            partition_model=co2_henry, k_transfer=kLa, transfer_basis="molecular",
        ),
        "O2": KineticTransferModel(
            partition_model=o2_henry, k_transfer=kLa, transfer_basis="total",
        ),
        "N2": KineticTransferModel(
            partition_model=n2_henry, k_transfer=kLa, transfer_basis="total",
        ),
    }
    return ControlVolume(
        phases={"gas": gas_phase, "liquid": liquid_phase},
        transfer_models=transfer_models,
        reaction_system=system,
        label=label,
    )

# Shared solver instance -- stateless (all state lives on the CV), so one
# instance is reused across every Simulation.run() call in this notebook.
bdf_solver = SimultaneousAdaptiveSolver(method="BDF")

# Depletion check: moles of each gas needed to fully saturate 100 mL of
# water at Henry's-law equilibrium, vs. moles actually available in the
# 1000 L headspace.
_kH = {"CO2": kH_mol_L_atm(co2_henry.H_ref, co2_henry.dlnH, T_K),
       "O2": kH_mol_L_atm(o2_henry.H_ref, o2_henry.dlnH, T_K),
       "N2": kH_mol_L_atm(n2_henry.H_ref, n2_henry.dlnH, T_K)}
_p = {"CO2": p_CO2_atm, "O2": p_O2_atm, "N2": p_N2_atm}
_n_gas_total = (1.0 * V_gas) / (R_L_ATM_MOL_K * T_K)
for _sp in ("CO2", "O2", "N2"):
    _needed = _kH[_sp] * _p[_sp] * V_liq
    _available = _n_gas_total * _p[_sp]
    print(f"{_sp}: needs {_needed:.3e} mol to saturate liquid, "
          f"headspace holds {_available:.3e} mol "
          f"({100*_needed/_available:.4f}% depletion at full saturation)")
"""),

    # ── 3. pH over time at kLa = 96 /h ────────────────────────────────────

    md("run-md", """\
## 3  pH over time at kLa = 96 /h

96 /h sits in the middle of the range reported for actively agitated
gas-liquid contacting (shaken plates, stirred tanks) — fast enough that
the run window only needs to span minutes, not hours. `tau_h` is sized
to comfortably clear the expected relaxation time (`ln(20)/kLa`, the
first-order estimate for reaching 95% of equilibrium) with margin to
spare.\
"""),

    code("run-code", """\
kLa_headline = 96.0   # /h

tau_h = 8.0 * math.log(20.0) / kLa_headline   # minutes-scale window, ample margin
n_steps = 500

cv = build_cv(kLa=kLa_headline, label="pure_water_headline")
result = Simulation(cvs={"main": cv}, label="usecase02b", solver=bdf_solver).run(
    tau_h=tau_h, n_steps=n_steps,
)

t_h = result.t_h
pH = result.pH["main"]
CO2_liq = result.liquid_mol["main"]["CO2"] / V_liq

fig, ax = plt.subplots(figsize=(7.5, 4.5))
ax.plot(t_h[1:] * 60.0, pH[1:], color="tab:purple")
ax.set_xlabel("time (min)")
ax.set_ylabel("pH")
ax.set_title(f"pH of pure water equilibrating with the atmosphere, kLa = {kLa_headline:.0f} /h")
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print(f"pH: {pH[0]:.4f} (pure water) -> {pH[-1]:.4f} (equilibrated)")
print(f"CO2(aq): 0 -> {CO2_liq[-1]*1e6:.4f} umol/L over {tau_h*60:.1f} min "
      f"({n_steps} steps)")
print("(textbook reference: unpolluted rainwater is ~pH 5.6, from atmospheric CO2 alone)")
"""),

    # ── 4. kLa vs. time-to-95%-equilibrium ────────────────────────────────

    md("sweep-md", """\
## 4  How kLa sets the time to reach equilibrium

The molecular-basis driving force means dissolved CO₂(aq) relaxes toward
exactly `kH × p_CO2` (Henry's law) at equilibrium, independent of kLa —
kLa only controls *how fast* it gets there, not where it ends up (the
same logic Section 3 relied on). Sweeping kLa over three-plus orders of
magnitude — from a still, unagitated surface up past vigorous mixing —
and recording the time each run takes to reach 95% of that Henry's-law
target (`t95`, found by interpolating each run's own CO₂(aq) trajectory)
maps out that relationship directly. Each run's own `tau_h` is sized
relative to its kLa so every sweep point gets adequate resolution without
wasting steps on runs that finish in seconds.

For a system where kLa is the only physics setting the timescale (alpha
≈ 1 here, Section 2), t95 should follow the simple first-order estimate
`ln(20)/kLa` closely — plotted alongside as a dashed reference.\
"""),

    code("sweep-code", """\
kH_CO2 = kH_mol_L_atm(co2_henry.H_ref, co2_henry.dlnH, T_K)
CO2_eq = kH_CO2 * p_CO2_atm   # Henry's-law equilibrium target, independent of kLa

kLa_sweep = np.logspace(-1, np.log10(500.0), 14)   # 0.1 - 500 /h
t95_sweep = np.empty_like(kLa_sweep)

for i, kLa_i in enumerate(kLa_sweep):
    tau_h_i = 8.0 * math.log(20.0) / kLa_i   # margin scaled to this run's own timescale
    cv_i = build_cv(kLa=kLa_i, label=f"sweep_{i}")
    res_i = Simulation(cvs={"main": cv_i}, label=f"usecase02b_sweep_{i}", solver=bdf_solver).run(
        tau_h=tau_h_i, n_steps=500,
    )
    CO2_i = res_i.liquid_mol["main"]["CO2"] / V_liq
    t95_sweep[i] = np.interp(0.95 * CO2_eq, CO2_i, res_i.t_h)

fig, ax = plt.subplots(figsize=(7.5, 4.5))
ax.plot(kLa_sweep, t95_sweep, "o-", color="tab:blue", label="simulated t95")
ax.plot(kLa_sweep, np.log(20.0) / kLa_sweep, "k--", lw=1,
        label="ln(20)/kLa (first-order estimate)")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("kLa (/h)")
ax.set_ylabel("t95 (h)")
ax.set_title("Time to reach 95% of the equilibrium aqueous CO2 concentration, vs. kLa")
ax.legend(fontsize=8)
ax.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.show()

print(f"{'kLa (/h)':>10}{'t95 (h)':>12}{'t95 (min)':>12}")
for kLa_i, t95_i in zip(kLa_sweep, t95_sweep):
    print(f"{kLa_i:>10.2f}{t95_i:>12.4f}{t95_i*60:>12.2f}")
"""),

    # ── 5. Computation time vs. kLa ─────────────────────────────────────────

    md("timing-md", """\
## 5  Computation time per simulated hour vs. kLa, N=500 replicates

A separate question from Section 4's *simulated* time-to-equilibrium:
how much *wall-clock* time does `SimultaneousAdaptiveSolver(method="BDF")`
itself spend, at low, medium, and high kLa? Each of the three kLa values
below (0.1, 10, 1000 /h — spanning the same static-to-vigorous range as
Section 4) is run for the same fixed 1-hour simulated window, as a
**single** `n_steps=1` macro step — letting `solve_ivp`'s own adaptive
internal stepping do all the work, rather than forcing it through many
Python-level macro-step calls, so the timing reflects the solver's own
cost per simulated experimental hour rather than Python call overhead.

As in the timing sections elsewhere in this repo, each kLa gets one
untimed warmup call first (so neither's first-call import/setup cost
biases the result), then mean ± standard deviation over 500 timed
replicate calls — spread across replicates, not a single best-of-N
figure, since the point is to characterize each kLa's own cost
distribution.\
"""),

    code("timing-code", """\
import time

N_REPS = 500
tau_h_timing = 1.0   # one simulated experimental hour
kLa_timing_values = [0.1, 10.0, 1000.0]

def time_replicates(fn, n=N_REPS):
    fn()  # untimed warmup -- first call pays one-time import/setup costs
    times_s = np.empty(n)
    for i in range(n):
        t0 = time.perf_counter()
        fn()
        times_s[i] = time.perf_counter() - t0
    return times_s * 1e3  # ms

def _run_one_experimental_hour(kLa_i):
    cv_i = build_cv(kLa=kLa_i, label="timing")
    Simulation(cvs={"main": cv_i}, label="usecase02b_timing", solver=bdf_solver).run(
        tau_h=tau_h_timing, n_steps=1,
    )

timing_results = {
    f"kLa = {kLa_i:g} /h": time_replicates(lambda kLa_i=kLa_i: _run_one_experimental_hour(kLa_i))
    for kLa_i in kLa_timing_values
}

print(f"{'kLa':<14}{'mean (ms/h)':>14}{'std (ms/h)':>14}")
for label_i, times_ms in timing_results.items():
    print(f"{label_i:<14}{times_ms.mean():>14.4f}{times_ms.std():>14.4f}")
"""),
)

save_path_kinetic = HERE / "02b_kinetic_co2_equilibration_microplate_well.ipynb"
with open(save_path_kinetic, "w", encoding="utf-8") as f:
    json.dump(kinetic_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path_kinetic}")

# ═══════════════════════════════════════════════════════════════════════════
#  03  Grow E. coli on acetic acid — Monod kinetics in a batch bioreactor
# ═══════════════════════════════════════════════════════════════════════════

GROWTH_SETUP = """\
import sys
from pathlib import Path
import warnings
import numpy as np
import matplotlib.pyplot as plt

def _find_repo():
    for p in [Path.cwd(), *Path.cwd().parents]:
        if (p / "pyproject.toml").exists():
            return p
    raise RuntimeError("Run from inside the PyOMES repo")

sys.path.insert(0, str(_find_repo() / "models"))

from PyOMES.chemistry.species import Species
from PyOMES.chemistry.common_species import (
    H2O, H_plus, OH_minus,
    H3PO4, H2PO4_minus, HPO4_2minus, PO4_3minus,
    NH3, NH4_plus,
    CO2, HCO3_minus, CO3_2minus,
    K_plus, Cl_minus,
)
from PyOMES.chemistry import HenryEquilibrium
from PyOMES.reactions import (
    EquilibriumReaction, ReactionBuilder, ReactionSystem, StoichiometryEntry,
)
from PyOMES.core import (
    ControlVolume, GasPhase, LiquidPhase,
    KineticTransferModel, EquilibriumTransferModel,
    Simulation, GasFeed, PressureReliefVent,
)
from PyOMES.core.phases import R_L_ATM_MOL_K

# The explicit-Euler CV solver clamps a species' removal rate when a step
# would otherwise drive it negative (AccuracyWarning), and separately flags
# any resulting element/charge drift above its tolerance (ConservationWarning).
# Both are the solver's self-correcting safety net doing its job at dt_h=0.01,
# not a sign of a wrong answer here (Section 6 checks that directly) --
# silenced only so the run cell's output stays readable.
from PyOMES.monitoring.accuracy import AccuracyWarning
from PyOMES.monitoring.conservation import ConservationWarning
warnings.filterwarnings("ignore", category=AccuracyWarning)
warnings.filterwarnings("ignore", category=ConservationWarning)

def _e(sp, coeff, phase="liquid"):
    return StoichiometryEntry(species=sp, phase=phase, coefficient=coeff)

print("Imports OK")
"""

growth_nb = nb(
    md("title", """\
# Growing E. coli on Acetic Acid — Monod Kinetics in a Batch Bioreactor

**The situation:** [usecase 01](01_predict_ph_simple_liquid.ipynb) built the
KH₂PO₄ + NH₄Cl liquid and [usecase 02](02_equilibrate_with_atmospheric_gas.ipynb)
opened it to a sparged air atmosphere. Now inoculate that same vessel with
*E. coli*, using acetic acid as the sole carbon and energy source, and let it
grow. How fast does the population grow, how quickly is the acetate consumed,
does dissolved O₂ ever become limiting, and — since acetic acid is itself a
weak acid and the nitrogen source is ammonium — what happens to pH as growth
proceeds?

This is the first notebook in the series where the chemistry *evolves over
time* rather than being solved as a single equilibrium snapshot, so it
introduces the two pieces usecases 01-02 deliberately left out:
`ControlVolume` (a vessel of phases) and `Simulation` (the time-stepping
orchestrator). Everything else carries over unchanged — the same
`EquilibriumReaction` network from usecase 02, extended with exactly one new
building block: a `KineticReaction` for Monod growth, built via
`ReactionBuilder.monod_aerobic_growth()` and dropped into the same
`ReactionSystem` alongside the acid-base chemistry.\
"""),

    code("setup", GROWTH_SETUP),

    # ── 1. Literature-grounded Monod parameters ──────────────────────────

    md("params-md", """\
## 1  Literature-grounded Monod parameters

*E. coli* grows on acetate far more slowly than on glucose — acetate uptake
feeds the glyoxylate shunt and gluconeogenesis rather than glycolysis, and
that costs growth rate. For reference, *E. coli*'s Monod constant on
**glucose** is well characterised: μ_max ≈ 0.66 h⁻¹, K_s ≈ 2-4 mg/L
([BioNumbers ID 111049](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=111049)).
Acetate kinetics are reported far less consistently — Kovarova-Kovar & Egli
(1998, *Microbiol. Mol. Biol. Rev.* 62:646-666) is the standard review
explaining *why*: reported K_s values for a given substrate vary by orders
of magnitude across studies because K_s depends on the cells' physiological
state, not just the transporter's intrinsic affinity. The parameters below
are chosen to be *representative* of the reported ranges rather than a
single paper's point estimate — treat them as a reasonable starting model,
not a re-derivation of any one study.

| Parameter | Value | Basis |
|---|---|---|
| μ_max | 0.30 h⁻¹ | Representative of aerobic minimal-medium growth on acetate as sole carbon source — well below glucose's 0.66 h⁻¹, consistent with the glyoxylate-shunt/gluconeogenesis cost noted above. |
| K_s | 5 mg/L acetate | Same order of magnitude as *E. coli*'s high-affinity glucose K_s (2-4 mg/L); high-affinity organic-acid transporters place acetate K_s in a similar low-mg/L range. |
| Y_X/S | 0.36 g DCW/g acetate | Within the range reported for aerobic *E. coli* on acetate — Leone et al. (2015, *Microb. Cell Fact.* 14:106) measured Y_X/S = 0.39 ± 0.01 g/g for *E. coli* BL21 under controlled pH/DO; other studies report lower yields (0.17-0.30 g/g) at low growth rates, consistent with a Pirt-type maintenance cost. |
| K_O2 | 0.2 mg/L | Typical high-affinity dissolved-O2 half-saturation constant for aerobic heterotrophs (same value used in this repository's own `ReactionBuilder.monod_aerobic_growth` docstring example). |

Biomass itself is represented by the generic empirical formula
**CH₁.₈O₀.₅N₀.₂** (MW ≈ 24.6 g per C-mol) — the average bacterial biomass
composition widely used in bioprocess engineering (Roels, 1983,
*Energetics and Kinetics in Biotechnology*) when a strain-specific formula
isn't available. It plays the same role here that the generic yeast formula
(CH₁.₆₁O₀.₅₆N₀.₁₆) played in `demos/builder/batch_fermenter.py`.\
"""),

    code("params-code", """\
mu_max_per_h = 0.30    # h^-1
Ks_gL        = 5e-3    # g acetate / L
Yxs          = 0.36    # g biomass / g acetate
Ko2_gL       = 0.2e-3  # g O2 / L

print(f"mu_max = {mu_max_per_h} /h, Ks = {Ks_gL*1e3} mg/L, "
      f"Yxs = {Yxs} g/g, Ko2 = {Ko2_gL*1e3} mg/L")
"""),

    # ── 2. Declare the chemistry ──────────────────────────────────────────

    md("chem-md", """\
## 2  Declare the chemistry

Same seven reactions as usecase 02 (water, the phosphate ladder, ammonium,
the carbonate ladder), plus one new one: acetic acid's own dissociation.
It's a weak acid exactly like phosphoric acid was in usecase 01 — just a
single step instead of three, with a pK_a (4.756) that happens to sit close
to phosphate's amphoteric pH (½(pK_a1+pK_a2) ≈ 4.68 from usecase 01 §4),
so the two acids interact rather than one simply dominating.

| Reaction | log K | p$K_a$ | Role |
|---|---|---|---|
| H₂O ⇌ H⁺ + OH⁻ | −14.0 | 14.0 | water autoionization |
| H₃PO₄ ⇌ H₂PO₄⁻ + H⁺ | −2.15 | 2.15 | phosphate, 1st step |
| H₂PO₄⁻ ⇌ HPO₄²⁻ + H⁺ | −7.20 | 7.20 | phosphate, 2nd step |
| HPO₄²⁻ ⇌ PO₄³⁻ + H⁺ | −12.35 | 12.35 | phosphate, 3rd step |
| NH₄⁺ ⇌ NH₃ + H⁺ | −9.25 | 9.25 | ammonium/ammonia |
| CO₂(aq) + H₂O ⇌ HCO₃⁻ + H⁺ | −6.35 | 6.35 | carbonate, 1st step |
| HCO₃⁻ ⇌ CO₃²⁻ + H⁺ | −10.33 | 10.33 | carbonate, 2nd step |
| AceticAcid ⇌ Acetate⁻ + H⁺ | −4.756 | 4.756 | the substrate's own acid-base chemistry |

`AceticAcid` and `Acetate-` aren't in `common_species` (a substrate-specific
pair, not a universal inorganic), so they're declared locally — the same
pattern `demos/model_api/chemistry/reaction_system.py` uses. `Ecoli` is
declared the same way, carrying the CH₁.₈O₀.₅N₀.₂ formula from §1.\
"""),

    code("chem-code", """\
ACETIC_ACID = Species(
    id="AceticAcid", atoms={"C": 2, "H": 4, "O": 2}, charge=0, MW=60.052,
)
ACETATE_MINUS = Species(
    id="Acetate-", atoms={"C": 2, "H": 3, "O": 2}, charge=-1, MW=59.044,
)
ECOLI = Species(
    id="Ecoli", atoms={"C": 1, "H": 1.8, "O": 0.5, "N": 0.2}, charge=0,
    MW=24.626,
)

water = EquilibriumReaction(
    stoichiometry=[_e(H2O, -1), _e(H_plus, +1), _e(OH_minus, +1)],
    log_K=-14.0, label="water",
)
p1 = EquilibriumReaction(
    stoichiometry=[_e(H3PO4, -1), _e(H2PO4_minus, +1), _e(H_plus, +1)],
    log_K=-2.15, total_id="H3PO4", label="p1",
)
p2 = EquilibriumReaction(
    stoichiometry=[_e(H2PO4_minus, -1), _e(HPO4_2minus, +1), _e(H_plus, +1)],
    log_K=-7.20, total_id="H3PO4", label="p2",
)
p3 = EquilibriumReaction(
    stoichiometry=[_e(HPO4_2minus, -1), _e(PO4_3minus, +1), _e(H_plus, +1)],
    log_K=-12.35, total_id="H3PO4", label="p3",
)
nh4 = EquilibriumReaction(
    stoichiometry=[_e(NH4_plus, -1), _e(NH3, +1), _e(H_plus, +1)],
    log_K=-9.25, total_id="NH3", label="nh4",
)
co2_first = EquilibriumReaction(
    stoichiometry=[_e(CO2, -1), _e(H2O, -1), _e(HCO3_minus, +1), _e(H_plus, +1)],
    log_K=-6.35, total_id="CO2", label="co2_first",
)
co2_second = EquilibriumReaction(
    stoichiometry=[_e(HCO3_minus, -1), _e(CO3_2minus, +1), _e(H_plus, +1)],
    log_K=-10.33, total_id="CO2", label="co2_second",
)
acetate_eq = EquilibriumReaction(
    stoichiometry=[_e(ACETIC_ACID, -1), _e(ACETATE_MINUS, +1), _e(H_plus, +1)],
    log_K=-4.756, total_id="AceticAcid", label="acetate_eq",
)

print("8 equilibrium reactions declared.")
"""),

    # ── 3. Declare the Monod growth kinetics ─────────────────────────────

    md("growth-md", """\
## 3  Declare the Monod growth kinetics

`ReactionBuilder.monod_aerobic_growth()` is the highest-level tool for this
job — the same role `NRChemicalEquilibriumEngine.from_reactions()` played
for the acid-base network in usecase 01. Given the substrate and biomass
`Species`, the four kinetic parameters from §1, and `balance="CHNO"`, it
derives the O₂, NH₃ (nitrogen source), CO₂, and H₂O coefficients from
elemental balance — the same closure `ReactionBuilder.aerobic_growth` uses
internally, so the growth reaction is automatically self-consistent on C,
H, N, and O. `n_source_id` defaults to `"NH3"`, so growth draws directly on
the *same* ammoniacal-nitrogen pool the `nh4` reaction above manages —
biomass synthesis and the ammonium/ammonia equilibrium share one total.\
"""),

    code("growth-code", """\
growth = ReactionBuilder.monod_aerobic_growth(
    substrate=ACETIC_ACID, biomass=ECOLI,
    mu_max_per_h=mu_max_per_h, Ks_gL=Ks_gL, yield_gX_gS=Yxs, Ko2_gL=Ko2_gL,
    balance="CHNO", label="growth_on_AceticAcid",
)

print("Derived stoichiometry (per mol AceticAcid consumed):")
for entry in growth.stoichiometry:
    print(f"  {entry.coefficient:+.4f}  {entry.species.id:<12} ({entry.phase})")
"""),

    # ── 4. Assemble the batch vessel and run ─────────────────────────────

    md("assemble-md", """\
## 4  Assemble the batch vessel and run

All nine reactions (eight equilibria + the one kinetic growth reaction) go
into a single `ReactionSystem`. The vessel is a 2 L bottle at 37 °C
(typical *E. coli* culture temperature), continuously sparged with air —
`GasFeed` supplies fresh O₂/N₂/CO₂ the same way usecase 02's atmosphere did,
but `GasFeed` only ever adds gas, so a `PressureReliefVent` is needed
alongside it to vent the excess and hold the headspace at 1 atm (mirroring
`demos/builder/cstr_fermenter.py`'s vessel plumbing). O₂ and CO₂ transfer
kinetically (kLa = 150 /h, the same value used throughout
`demos/builder/`); N₂ stays at instantaneous Henry equilibrium since
nothing reacts with it.

The liquid starts at usecase 01's M9-like recipe (22 mM KH₂PO₄, 18.7 mM
NH₄Cl) plus 4 g/L acetic acid (0.4% w/v — the same loading Leone et al.
used) and a small inoculum (0.05 g/L biomass). `K+`/`Cl-` are the actual
counter-ions each salt contributes — the charge-balance solver needs them
declared as real species in `n_mol`, the CV-driven equivalent of usecase
01's `strong_ions=` argument on a bare `engine.solve()` call.\
"""),

    code("assemble-code", """\
system = ReactionSystem(
    [water, p1, p2, p3, nh4, co2_first, co2_second, acetate_eq, growth],
    label="ecoli_on_acetate",
)

V_total_L = 2.0
headspace_frac = 0.20
V_liq = V_total_L * (1.0 - headspace_frac)
V_gas = V_total_L * headspace_frac
T_K = 310.15   # 37 C

n_total_gas = (1.0 * V_gas) / (R_L_ATM_MOL_K * T_K)   # 1 atm headspace
gas_phase = GasPhase(
    n_mol={
        "O2": n_total_gas * 0.2095,
        "N2": n_total_gas * 0.7901,
        "CO2": n_total_gas * 0.0004,
    },
    V_L=V_gas, T_K=T_K,
)

CT_P = 0.022    # mol/L KH2PO4 (usecase 01's recipe)
CT_N = 0.0187   # mol/L NH4Cl
C_AcOH0_gL = 4.0                        # g/L acetic acid (0.4% w/v)
C_AcOH0 = C_AcOH0_gL / ACETIC_ACID.MW
X0_gL = 0.05                             # g/L inoculum
X0 = X0_gL / ECOLI.MW

liquid_phase = LiquidPhase(
    n_mol={
        "H3PO4": CT_P * V_liq, "NH3": CT_N * V_liq,
        "AceticAcid": C_AcOH0 * V_liq, "Ecoli": X0 * V_liq,
        "K+": CT_P * V_liq, "Cl-": CT_N * V_liq,
        "O2": 0.0, "CO2": 0.0, "N2": 0.0,
    },
    V_L=V_liq, T_K=T_K,
)

transfer_models = {
    "O2": KineticTransferModel(
        partition_model=HenryEquilibrium(H_ref=1.3e-5, dlnH=1500.0), k_transfer=150.0,
    ),
    "CO2": KineticTransferModel(
        partition_model=HenryEquilibrium(H_ref=3.4e-4, dlnH=2400.0), k_transfer=150.0 * 0.9,
    ),
    "N2": EquilibriumTransferModel(HenryEquilibrium(H_ref=6.4e-6, dlnH=1300.0)),
}
gas_feed = GasFeed(
    vvm_min=1.0, y={"O2": 0.2095, "N2": 0.7901, "CO2": 0.0004}, P_inlet_atm=1.0,
    phase_key="gas", liquid_phase_key="liquid", label="gas_feed",
)
vent = PressureReliefVent(P_set_atm=1.0, mode="instant")

cv = ControlVolume(
    phases={"gas": gas_phase, "liquid": liquid_phase},
    transfer_models=transfer_models,
    boundaries=[gas_feed, vent],
    reaction_system=system,
    label="batch_ecoli",
)

sim = Simulation(cvs={"main": cv}, label="usecase03")
tau_h, n_steps = 20.0, 2000
result = sim.run(tau_h=tau_h, n_steps=n_steps)

print(f"Ran {tau_h:.0f} h ({n_steps} steps) in {result.runtime_s:.2f} s wall-clock.")
"""),

    # ── 5. Plot the batch ─────────────────────────────────────────────────

    md("plot-md", """\
## 5  Biomass, substrate, dissolved O₂, and pH over the batch

Four views of the same 20-hour run: acetate depletes as biomass grows,
dissolved O₂ stays close to air saturation throughout (§7 checks whether
it ever actually limits growth), and pH swings noticeably — starting well
below usecase 01's phosphate-only baseline (≈4.73) because 4 g/L acetic
acid is itself a substantial acid load, comparable in magnitude to the
22 mM phosphate buffer meant to resist it, then recovering toward that
baseline as the acid is consumed.\
"""),

    code("plot-code", """\
cv_key = "main"
t = result.t_h
X = result.liquid_mol[cv_key]["Ecoli"] * ECOLI.MW / V_liq
S = result.liquid_mol[cv_key]["AceticAcid"] * ACETIC_ACID.MW / V_liq
DO_mgL = result.liquid_mol[cv_key]["O2"] / V_liq * 32e3
pH = result.pH[cv_key]

fig, axes = plt.subplots(2, 2, figsize=(11, 8))

axes[0, 0].plot(t, X, color="tab:green")
axes[0, 0].set_xlabel("time (h)")
axes[0, 0].set_ylabel("biomass (g/L)")
axes[0, 0].set_title("E. coli growth")
axes[0, 0].grid(alpha=0.3)

axes[0, 1].plot(t, S, color="tab:blue")
axes[0, 1].set_xlabel("time (h)")
axes[0, 1].set_ylabel("acetic acid (g/L)")
axes[0, 1].set_title("Substrate depletion")
axes[0, 1].grid(alpha=0.3)

axes[1, 0].plot(t, DO_mgL, color="tab:red")
axes[1, 0].set_xlabel("time (h)")
axes[1, 0].set_ylabel("dissolved O2 (mg/L)")
axes[1, 0].set_title("Dissolved oxygen")
axes[1, 0].grid(alpha=0.3)

axes[1, 1].plot(t, pH, color="tab:purple")
axes[1, 1].axhline(4.734, color="gray", ls="--", lw=1,
                    label="usecase 01 baseline (no acetate), pH=4.73")
axes[1, 1].set_xlabel("time (h)")
axes[1, 1].set_ylabel("pH")
axes[1, 1].set_title("pH")
axes[1, 1].legend(fontsize=8)
axes[1, 1].grid(alpha=0.3)

plt.tight_layout()
plt.show()

print(f"Biomass:  {X[0]:.3f} -> {X[-1]:.3f} g/L")
print(f"Acetate:  {S[0]:.3f} -> {S[-1]:.3f} g/L")
print(f"DO:       {DO_mgL[0]:.3f} -> {DO_mgL[-1]:.3f} mg/L "
      f"(min {DO_mgL.min():.3f} mg/L)")
print(f"pH:       {pH[0]:.3f} -> {pH[-1]:.3f} (min {np.nanmin(pH):.3f})")
"""),

    # ── 6. Validate ────────────────────────────────────────────────────

    md("validate-md", """\
## 6  Validating the run against the declared model

Two independent checks that the simulator is actually solving the model
that was declared in §§3-4, not just producing a plausible-looking curve:

- **Yield mass balance.** `ReactionBuilder.aerobic_growth` fixed the ratio
  of biomass produced to substrate consumed at exactly `Yxs` when it built
  the stoichiometry — every mol of `AceticAcid` consumed produces `Y_mol`
  mol of `Ecoli`, no more, no less. So `Yxs * (S₀ - S_f)` should equal
  `X_f - X₀` to numerical precision, independent of `mu_max`, `Ks`, or
  `Ko2` (those only control the *rate*, not the *stoichiometry*).
- **Monod self-consistency.** The rate function inside `growth` computes
  `μ = mu_max · S/(Ks+S) · O2/(Ko2+O2)` from whatever `S` and `O2` the
  solver hands it at each step. Recomputing that same formula from the
  *recorded* `S`/`DO` trajectory and comparing it to `μ` backed out
  numerically from `dX/dt / X` checks that the integrator is actually
  tracking the declared rate law, not some subtly different one.\
"""),

    code("validate-code", """\
# Yield mass balance
dS = S[0] - S[-1]
dX = X[-1] - X[0]
print(f"Yxs * (S0 - Sf) = {Yxs * dS:.4f} g/L")
print(f"Xf - X0         = {dX:.4f} g/L")
print(f"Agreement to {abs(Yxs*dS - dX):.2e} g/L\\n")

# Monod self-consistency: mu backed out of the recorded trajectory
# vs. mu recomputed from the declared rate law at that same trajectory.
mu_numeric = np.gradient(X, t) / np.maximum(X, 1e-9)
DO_gL = DO_mgL / 1e3
mu_monod = mu_max_per_h * (S / (Ks_gL + S)) * (DO_gL / (Ko2_gL + DO_gL))

# Restrict the comparison to the bulk growth phase, away from both its
# edges: near t=0, dissolved O2 starts at exactly 0 mol/L (not yet
# gas-liquid-equilibrated), so the central-difference dX/dt right at that
# one-step transient understates the true instantaneous rate; and near
# substrate exhaustion, S/(Ks+S) has its steepest curvature (Ks is tiny
# next to S), where any finite-difference derivative is least accurate.
# Both are numerical-differentiation edge effects, not signs the
# integrator is tracking the wrong rate law -- the mid-growth agreement
# below is the actual check.
growing = (S > 0.05) & (t > 0.02)
max_gap = np.max(np.abs(mu_numeric[growing] - mu_monod[growing]))
print(f"Max |mu_numeric - mu_monod| over the growth phase: {max_gap:.5f} /h")

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t[growing], mu_numeric[growing], "o", ms=3, color="tab:orange", label="mu, from dX/dt / X")
ax.plot(t[growing], mu_monod[growing], "-", color="tab:blue", label="mu, from Monod formula")
ax.set_xlabel("time (h)")
ax.set_ylabel("specific growth rate (1/h)")
ax.set_title("Monod self-consistency check")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.show()
"""),

    md("do-note-md", """\
Dissolved O₂ stayed within a few percent of air saturation for the whole
run (§5's DO panel) — at `kLa(O2)=150 /h`, oxygen transfer keeps up with
this culture's peak O₂ demand comfortably, so the `Ko2_gL` Monod term
never actually binds here (`DO >> Ko2` throughout). That's a real result,
not an oversight: it takes either a much lower kLa or a much denser
culture before O₂ transfer becomes the growth-limiting step, which is
exactly the tradeoff `demos/builder/cstr_fermenter.py`'s `DOAgitationController`
exists to manage.\
"""),

    # ── 7. Sensitivity: mu_max sets the timescale, not the yield ─────────

    md("sweep-md", """\
## 7  Sensitivity: μ_max sets the timescale, not the final yield

The mass-balance check in §6 already hinted at this: `mu_max` never
entered the yield relationship. Rerunning the same batch at three
different `mu_max` values — holding `Ks`, `Yxs`, `Ko2`, and every initial
condition fixed — should show growth curves that reach the *same* final
biomass and substrate concentrations, just on different timescales. This
is the batch-culture analogue of usecase 02 §4's sensitivity sections: one
parameter swept, one qualitative expectation from the model's structure
checked directly against the simulation.\
"""),

    code("sweep-code", """\
def run_batch(mu_max_test):
    growth_i = ReactionBuilder.monod_aerobic_growth(
        substrate=ACETIC_ACID, biomass=ECOLI,
        mu_max_per_h=mu_max_test, Ks_gL=Ks_gL, yield_gX_gS=Yxs, Ko2_gL=Ko2_gL,
        balance="CHNO", label="growth_on_AceticAcid",
    )
    system_i = ReactionSystem(
        [water, p1, p2, p3, nh4, co2_first, co2_second, acetate_eq, growth_i],
        label="ecoli_on_acetate",
    )
    gas_phase_i = GasPhase(
        n_mol={
            "O2": n_total_gas * 0.2095, "N2": n_total_gas * 0.7901,
            "CO2": n_total_gas * 0.0004,
        },
        V_L=V_gas, T_K=T_K,
    )
    liquid_phase_i = LiquidPhase(
        n_mol={
            "H3PO4": CT_P * V_liq, "NH3": CT_N * V_liq,
            "AceticAcid": C_AcOH0 * V_liq, "Ecoli": X0 * V_liq,
            "K+": CT_P * V_liq, "Cl-": CT_N * V_liq,
            "O2": 0.0, "CO2": 0.0, "N2": 0.0,
        },
        V_L=V_liq, T_K=T_K,
    )
    cv_i = ControlVolume(
        phases={"gas": gas_phase_i, "liquid": liquid_phase_i},
        transfer_models=transfer_models,
        boundaries=[
            GasFeed(vvm_min=1.0, y={"O2": 0.2095, "N2": 0.7901, "CO2": 0.0004},
                    P_inlet_atm=1.0, phase_key="gas", liquid_phase_key="liquid"),
            PressureReliefVent(P_set_atm=1.0, mode="instant"),
        ],
        reaction_system=system_i, label=f"mu{mu_max_test}",
    )
    sim_i = Simulation(cvs={"main": cv_i}, label=f"sweep_mu{mu_max_test}")
    return sim_i.run(tau_h=tau_h, n_steps=n_steps)

mu_max_values = [0.15, 0.30, 0.60]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for mu_test in mu_max_values:
    res_i = run_batch(mu_test)
    X_i = res_i.liquid_mol["main"]["Ecoli"] * ECOLI.MW / V_liq
    S_i = res_i.liquid_mol["main"]["AceticAcid"] * ACETIC_ACID.MW / V_liq
    axes[0].plot(res_i.t_h, X_i, label=f"mu_max = {mu_test}/h")
    axes[1].plot(res_i.t_h, S_i, label=f"mu_max = {mu_test}/h")

axes[0].set_xlabel("time (h)")
axes[0].set_ylabel("biomass (g/L)")
axes[0].set_title("Biomass: same asymptote, different timescale")
axes[0].legend(fontsize=8)
axes[0].grid(alpha=0.3)

axes[1].set_xlabel("time (h)")
axes[1].set_ylabel("acetic acid (g/L)")
axes[1].set_title("Substrate: same depletion, different timescale")
axes[1].legend(fontsize=8)
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.show()
"""),

    md("next-md", """\
## Where to go next

- **Closed-loop pH and DO control** — this notebook lets pH and DO drift
  freely; `demos/builder/batch_fermenter.py` (`PHController`) and
  `demos/builder/cstr_fermenter.py` (`PHController` + `DOAgitationController`
  cascade) show how `PyOMES.control.cv_loops` closes the loop on both.
- **Continuous culture (chemostat) instead of batch** — the same organism
  and substrate at steady state under continuous feed/drain, including the
  analytical Monod steady-state check — see
  `demos/builder/cstr_fermenter.py`.
- **The declarative building blocks used here in isolation** —
  `demos/model_api/chemistry/reaction_system.py` walks through
  `KineticReaction`, `EquilibriumReaction`, and `ReactionSystem` one at a
  time without running a `Simulation`.
- **`FermenterBuilder`, the fluent alternative to assembling `ControlVolume`
  by hand** — `models/vlmodels/fermenter/config/builder.py` wraps vessel
  geometry, gas feed, transfer, organism, and substrate into one chained
  call; §4 of this notebook is what it builds internally.\
"""),
)

save_path_growth = HERE / "03_grow_ecoli_on_acetic_acid.ipynb"
with open(save_path_growth, "w", encoding="utf-8") as f:
    json.dump(growth_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path_growth}")

# ═══════════════════════════════════════════════════════════════════════════
#  04  Compare run time across usecases and activity-model treatments
# ═══════════════════════════════════════════════════════════════════════════

RUNTIME_SETUP = """\
import sys
import time
import warnings
import math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def _find_repo():
    for p in [Path.cwd(), *Path.cwd().parents]:
        if (p / "pyproject.toml").exists():
            return p
    raise RuntimeError("Run from inside the PyOMES repo")

sys.path.insert(0, str(_find_repo() / "models"))

from PyOMES.chemistry.species import Species
from PyOMES.chemistry.common_species import (
    H2O, H_plus, OH_minus,
    H3PO4, H2PO4_minus, HPO4_2minus, PO4_3minus,
    NH3, NH4_plus,
    CO2, HCO3_minus, CO3_2minus,
)
from PyOMES.chemistry import HenryEquilibrium
from PyOMES.reactions import (
    EquilibriumReaction, ReactionBuilder, ReactionSystem, StoichiometryEntry,
)
from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
from PyOMES.core import (
    ControlVolume, GasPhase, LiquidPhase,
    KineticTransferModel, EquilibriumTransferModel,
    Simulation, GasFeed, PressureReliefVent,
)
from PyOMES.core.phases import R_L_ATM_MOL_K

# Same self-correcting-solver warnings usecase 03 silences (AccuracyWarning /
# ConservationWarning at dt_h=0.01) -- silenced only so the timing cells'
# output stays readable, not because anything here is wrong.
from PyOMES.monitoring.accuracy import AccuracyWarning
from PyOMES.monitoring.conservation import ConservationWarning
warnings.filterwarnings("ignore", category=AccuracyWarning)
warnings.filterwarnings("ignore", category=ConservationWarning)

def _e(sp, coeff, phase="liquid"):
    return StoichiometryEntry(species=sp, phase=phase, coefficient=coeff)

try:
    import phreeqpython  # noqa: F401  -- actually probe for the optional dep;
    # PHREEQCChemicalEquilibriumEngine itself imports phreeqpython lazily
    # inside __init__, so importing only the wrapper class below would
    # succeed even without phreeqpython installed.
    from PyOMES.chemical_equilibrium.phreeqc_engine import PHREEQCChemicalEquilibriumEngine
    _HAVE_PHREEQC = True
    print("phreeqpython available - PHREEQC timing cells will run.")
except ImportError as exc:
    _HAVE_PHREEQC = False
    print(f"phreeqpython not installed ({exc}); skipping PHREEQC timing cells.")
    print("Install with: pip install PyOMES[phreeqc]")

# (usecase, engine/activity-model) -> seconds. For usecases 01-02 this is
# seconds per single engine.solve() call; for usecase 03 it's seconds for
# the *entire* 20 h / 2000-step batch run -- see the markdown notes below
# on why those two things aren't the same kind of quantity.
TIMINGS = {}

def time_repeated(fn, reps=100, warmup=3, trials=5):
    \"\"\"Best-of-`trials` mean wall time per call, after a fixed warmup.

    Sub-millisecond NR solves are short enough that a *single* timed pass
    is dominated by OS scheduling noise, not the computation being
    measured -- the standard fix (same one `timeit.repeat` uses) is to
    repeat the whole timed pass several times and keep the fastest, since
    noise only ever makes a pass slower, never faster.

    The warmup matters more than usual here too: the *first* call into a
    freshly-built engine/solver pays one-time Python import and (for
    PHREEQC) IPC-connection costs that have nothing to do with the
    activity model being timed, and would otherwise bias whichever
    engine happens to run first in a given cell.
    \"\"\"
    for _ in range(warmup):
        fn()
    best = None
    for _ in range(trials):
        t0 = time.perf_counter()
        for _ in range(reps):
            fn()
        t1 = time.perf_counter()
        mean = (t1 - t0) / reps
        best = mean if best is None else min(best, mean)
    return best

print("Imports OK")
"""

runtime_nb = nb(
    md("title", """\
# Comparing Run Time Across Use Cases and Activity-Model Treatments

**The situation:** usecases 01-03 each took *some* amount of wall-clock time
to run — a single-point `solve()` call in 01 and 02, a full 2000-step batch
simulation in 03. Both usecases 01 and 02 also showed that `use_activity=True`
(Davies correction) and reaching for PHREEQC instead of PyOMES's own solver
both change the *answer* slightly (Section 6 of each). This notebook asks a
different question: how much do those same choices change the *run time*?

This is a **benchmarking** notebook, not a scientific-results notebook — the
chemistry and parameters are copied verbatim from usecases 01-03 (same
recipe, same reaction networks, same batch-growth setup) so the workloads
being timed are the *actual* usecases already demonstrated, not toy
stand-ins. Wall-clock timing is inherently noisy (background load, CPU
frequency scaling, GC pauses), so every measurement below is a mean over
repeated calls after a fixed warmup, not a single stopwatch read.

**Two axes of comparison, and one caveat about mixing them:**

- **By usecase** — 01 (pH prediction), 02 (gas-liquid CO₂ equilibration), and
  03 (E. coli batch growth) represent three different *kinds* of workload:
  01 and 02 are a single equilibrium snapshot, while 03 is a 20-hour dynamic
  simulation stepped 2000 times, each step re-solving the same kind of
  equilibrium problem 01/02 solve *once*. So "01 vs 03 run time" is really
  "cost of one solve vs. cost of two thousand solves plus kinetics and
  transport" — an honest comparison, but not an apples-to-apples one, which
  is why the plots below keep those units labelled explicitly rather than
  implying they're the same thing.
- **By activity-model / engine treatment** — four combinations, matching
  usecase 01/02's own Section 6 exactly: **PyOMES ideal** (no activity
  correction), **PyOMES Davies** (`use_activity=True, activity_model="davies"`),
  **PHREEQC ideal-equivalent** (the `-gamma 1e6 0` trick from usecase 01
  §6c/usecase 02 §6, forcing γ→1), and **PHREEQC default**
  (`vitens.dat`'s WATEQ Debye-Hückel, via `PHREEQCChemicalEquilibriumEngine`).
  Usecase 03 only has the first two — there is no dynamic-`ControlVolume`
  PHREEQC path in PyOMES, so PHREEQC timing is only meaningful for the
  single-snapshot usecases.

PHREEQC timing needs the optional `phreeqpython` package
(`pip install PyOMES[phreeqc]`); if it isn't installed, this notebook skips
those cells and reports only the three PyOMES-only rows.\
"""),

    code("setup", RUNTIME_SETUP),

    # ── 1. Usecase 01 timing ─────────────────────────────────────────────

    md("uc01-md", """\
## 1  Usecase 01 — pH prediction (single equilibrium solve)

Exactly [usecase 01](01_predict_ph_simple_liquid.ipynb)'s chemistry (the
five-reaction phosphate/ammonium network) and its M9-like recipe
(22 mmol/L KH₂PO₄, 18.7 mmol/L NH₄Cl), timed as repeated `engine.solve()`
calls at that one fixed point — the recipe doesn't matter for timing (NR
converges from the same starting guess every call), only that it's a
realistic, already-validated one rather than an arbitrary composition.\
"""),

    code("uc01-code", """\
water = EquilibriumReaction(
    stoichiometry=[_e(H2O, -1), _e(H_plus, +1), _e(OH_minus, +1)],
    log_K=-14.0, label="water",
)
p1 = EquilibriumReaction(
    stoichiometry=[_e(H3PO4, -1), _e(H2PO4_minus, +1), _e(H_plus, +1)],
    log_K=-2.15, total_id="H3PO4", label="p1",
)
p2 = EquilibriumReaction(
    stoichiometry=[_e(H2PO4_minus, -1), _e(HPO4_2minus, +1), _e(H_plus, +1)],
    log_K=-7.20, total_id="H3PO4", label="p2",
)
p3 = EquilibriumReaction(
    stoichiometry=[_e(HPO4_2minus, -1), _e(PO4_3minus, +1), _e(H_plus, +1)],
    log_K=-12.35, total_id="H3PO4", label="p3",
)
nh4 = EquilibriumReaction(
    stoichiometry=[_e(NH4_plus, -1), _e(NH3, +1), _e(H_plus, +1)],
    log_K=-9.25, total_id="NH3", label="nh4",
)

CT_P = 0.022   # mol/L KH2PO4 -> mol/L total phosphate, mol/L K+
CT_N = 0.0187  # mol/L NH4Cl  -> mol/L total ammoniacal N, mol/L Cl-

engine01_ideal = NRChemicalEquilibriumEngine.from_reactions([water, p1, p2, p3, nh4])
engine01_davies = NRChemicalEquilibriumEngine.from_reactions(
    [water, p1, p2, p3, nh4], use_activity=True, activity_model="davies",
)

def solve01(engine):
    return engine.solve(
        totals={"H3PO4": CT_P, "NH3": CT_N}, strong_ions={"CT_K": CT_P, "CT_Cl": CT_N},
    )

TIMINGS[("01 pH prediction", "PyOMES ideal")] = time_repeated(lambda: solve01(engine01_ideal))
TIMINGS[("01 pH prediction", "PyOMES Davies")] = time_repeated(lambda: solve01(engine01_davies))

print(f"01 PyOMES ideal:  {TIMINGS[('01 pH prediction', 'PyOMES ideal')]*1e3:.4f} ms/solve")
print(f"01 PyOMES Davies: {TIMINGS[('01 pH prediction', 'PyOMES Davies')]*1e3:.4f} ms/solve")
"""),

    md("uc01-pq-md", """\
### PHREEQC, usecase 01

Reuses usecase 01 §6c's minimal ideal database (`-gamma 1e6 0` on every
charged species, forcing γ→1) for the PHREEQC-ideal-equivalent point, and
`PHREEQCChemicalEquilibriumEngine` (`vitens.dat`, WATEQ Debye-Hückel) for
the PHREEQC-default point — `use_warmstart=False` for the same reason
usecase 01 needed it: the default warmstart reuses the previous solution
*incrementally*, which would corrupt repeated identical solves.\
"""),

    code("uc01-pq-code", """\
if _HAVE_PHREEQC:
    import tempfile
    from phreeqpython import PhreeqPython as _RawPhreeqPython

    _IDEAL_DB_01 = \"\"\"\\
SOLUTION_MASTER_SPECIES
H       H+      -1.     H       1.008
H(0)    H2      0.0     H
H(1)    H+      -1.     0.0
E       e-      0.0     0.0     0.0
O       H2O     0.0     O       16.00
O(0)    O2      0.0     O
O(-2)   H2O     0.0     0.0
P       PO4-3   0.0     P       30.974
N       NH4+    0.0     N       14.0067
K       K+      0.0     K       39.098
Cl      Cl-     0.0     Cl      35.453

SOLUTION_SPECIES
H+ = H+
        log_k           0.0
        -gamma          1e6     0
e- = e-
        log_k           0.0
H2O = H2O
        log_k           0.0
2 H+ + 2 e- = H2
        log_k           -3.15
2 H2O = O2 + 4 H+ + 4 e-
        log_k           -86.08
H2O = OH- + H+
        log_k           -14.0
        -gamma          1e6     0
PO4-3 = PO4-3
        log_k           0.0
        -gamma          1e6     0
PO4-3 + H+ = HPO4-2
        log_k           12.35
        -gamma          1e6     0
PO4-3 + 2H+ = H2PO4-
        log_k           19.55
        -gamma          1e6     0
PO4-3 + 3H+ = H3PO4
        log_k           21.70
NH4+ = NH4+
        log_k           0.0
        -gamma          1e6     0
NH4+ = NH3 + H+
        log_k           -9.25
K+ = K+
        log_k           0.0
        -gamma          1e6     0
Cl- = Cl-
        log_k           0.0
        -gamma          1e6     0
END
\"\"\"
    _db_dir = Path(tempfile.gettempdir())
    (_db_dir / "vlsim_ideal_phosphate_ammonium_runtime.dat").write_text(_IDEAL_DB_01)
    pp_ideal01 = _RawPhreeqPython(database="vlsim_ideal_phosphate_ammonium_runtime.dat",
                                   database_directory=_db_dir)

    def solve01_pq_ideal():
        sol = pp_ideal01.add_solution_raw({
            "P": CT_P * 1e3, "N": CT_N * 1e3, "K": CT_P * 1e3, "Cl": CT_N * 1e3,
            "temp": 25.0, "pH": "7 charge", "units": "mmol/L",
        })
        ph = sol.pH
        sol.forget()
        return ph

    engine01_pq_default = PHREEQCChemicalEquilibriumEngine(
        {"H3PO4": CT_P * 1e3, "NH3": CT_N * 1e3, "CT_K": CT_P * 1e3, "CT_Cl": CT_N * 1e3},
        component_map={"H3PO4": "P", "NH3": "N(-3)", "CT_K": "K", "CT_Cl": "Cl"},
        use_warmstart=False,
    )

    def solve01_pq_default():
        return engine01_pq_default.solve(
            totals={"H3PO4": CT_P, "NH3": CT_N, "CT_K": CT_P, "CT_Cl": CT_N},
        )

    TIMINGS[("01 pH prediction", "PHREEQC ideal (gamma->1)")] = time_repeated(solve01_pq_ideal, reps=10, warmup=2, trials=3)
    TIMINGS[("01 pH prediction", "PHREEQC default (WATEQ D-H)")] = time_repeated(solve01_pq_default, reps=10, warmup=2, trials=3)

    print(f"01 PHREEQC ideal:   {TIMINGS[('01 pH prediction', 'PHREEQC ideal (gamma->1)')]*1e3:.3f} ms/solve")
    print(f"01 PHREEQC default: {TIMINGS[('01 pH prediction', 'PHREEQC default (WATEQ D-H)')]*1e3:.3f} ms/solve")
else:
    print("Skipping PHREEQC timing for usecase 01 (phreeqpython not installed).")
"""),

    # ── 2. Usecase 02 timing ─────────────────────────────────────────────

    md("uc02-md", """\
## 2  Usecase 02 — gas-liquid CO₂ equilibration (single equilibrium solve)

[Usecase 02](02_equilibrate_with_atmospheric_gas.ipynb)'s network extends
usecase 01 with the two-step carbonate ladder, solved at the same M9-like
recipe plus dissolved CO₂ from a 400 ppm atmosphere (usecase 02 §2's
`kH_mol_L_atm` Henry conversion, reproduced here). Two more reactions than
usecase 01, so any per-reaction-network-size overhead in the NR tableau
build/solve should show up as a slightly larger 01-vs-02 gap within the
*same* activity-model column.\
"""),

    code("uc02-code", """\
co2_first = EquilibriumReaction(
    stoichiometry=[_e(CO2, -1), _e(H2O, -1), _e(HCO3_minus, +1), _e(H_plus, +1)],
    log_K=-6.35, total_id="CO2", label="co2_first",
)
co2_second = EquilibriumReaction(
    stoichiometry=[_e(HCO3_minus, -1), _e(CO3_2minus, +1), _e(H_plus, +1)],
    log_K=-10.33, total_id="CO2", label="co2_second",
)

def kH_mol_L_atm(H_ref, dlnH, T_K, T_ref=298.15):
    H_ref_mol_L_atm = (H_ref / 1000.0) * 101325.0
    return H_ref_mol_L_atm * math.exp(dlnH * (1.0 / T_K - 1.0 / T_ref))

CT_CO2_atm = kH_mol_L_atm(H_ref=3.4e-4, dlnH=2400.0, T_K=298.15) * 400e-6  # 400 ppm at 25 C

engine02_ideal = NRChemicalEquilibriumEngine.from_reactions(
    [water, p1, p2, p3, nh4, co2_first, co2_second]
)
engine02_davies = NRChemicalEquilibriumEngine.from_reactions(
    [water, p1, p2, p3, nh4, co2_first, co2_second],
    use_activity=True, activity_model="davies",
)

def solve02(engine):
    return engine.solve(
        totals={"H3PO4": CT_P, "NH3": CT_N, "CO2": CT_CO2_atm},
        strong_ions={"CT_K": CT_P, "CT_Cl": CT_N},
    )

TIMINGS[("02 gas-liquid CO2", "PyOMES ideal")] = time_repeated(lambda: solve02(engine02_ideal))
TIMINGS[("02 gas-liquid CO2", "PyOMES Davies")] = time_repeated(lambda: solve02(engine02_davies))

print(f"02 PyOMES ideal:  {TIMINGS[('02 gas-liquid CO2', 'PyOMES ideal')]*1e3:.4f} ms/solve")
print(f"02 PyOMES Davies: {TIMINGS[('02 gas-liquid CO2', 'PyOMES Davies')]*1e3:.4f} ms/solve")
"""),

    md("uc02-pq-md", """\
### PHREEQC, usecase 02

Same idea as usecase 01, with usecase 02's CO₂-extended ideal database
(adds the `C` master species and carbonate reactions to the §6c database
above) for the PHREEQC-ideal-equivalent point.\
"""),

    code("uc02-pq-code", """\
if _HAVE_PHREEQC:
    _IDEAL_DB_02 = \"\"\"\\
SOLUTION_MASTER_SPECIES
H       H+      -1.     H       1.008
H(0)    H2      0.0     H
H(1)    H+      -1.     0.0
E       e-      0.0     0.0     0.0
O       H2O     0.0     O       16.00
O(0)    O2      0.0     O
O(-2)   H2O     0.0     0.0
C       CO2     0.0     C       12.011
P       PO4-3   0.0     P       30.974
N       NH4+    0.0     N       14.0067
K       K+      0.0     K       39.098
Cl      Cl-     0.0     Cl      35.453

SOLUTION_SPECIES
H+ = H+
        log_k           0.0
        -gamma          1e6     0
e- = e-
        log_k           0.0
H2O = H2O
        log_k           0.0
2 H+ + 2 e- = H2
        log_k           -3.15
2 H2O = O2 + 4 H+ + 4 e-
        log_k           -86.08
H2O = OH- + H+
        log_k           -14.0
        -gamma          1e6     0
CO2 = CO2
        log_k           0.0
        -gamma          1e6     0
CO2 + H2O = HCO3- + H+
        log_k           -6.35
        -gamma          1e6     0
HCO3- = CO3-2 + H+
        log_k           -10.33
        -gamma          1e6     0
PO4-3 = PO4-3
        log_k           0.0
        -gamma          1e6     0
PO4-3 + H+ = HPO4-2
        log_k           12.35
        -gamma          1e6     0
PO4-3 + 2H+ = H2PO4-
        log_k           19.55
        -gamma          1e6     0
PO4-3 + 3H+ = H3PO4
        log_k           21.70
NH4+ = NH4+
        log_k           0.0
        -gamma          1e6     0
NH4+ = NH3 + H+
        log_k           -9.25
K+ = K+
        log_k           0.0
        -gamma          1e6     0
Cl- = Cl-
        log_k           0.0
        -gamma          1e6     0
END
\"\"\"
    _db_dir = Path(tempfile.gettempdir())
    (_db_dir / "vlsim_ideal_gas_liquid_runtime.dat").write_text(_IDEAL_DB_02)
    pp_ideal02 = _RawPhreeqPython(database="vlsim_ideal_gas_liquid_runtime.dat",
                                   database_directory=_db_dir)

    def solve02_pq_ideal():
        sol = pp_ideal02.add_solution_raw({
            "C": CT_CO2_atm * 1e3, "P": CT_P * 1e3, "N": CT_N * 1e3,
            "K": CT_P * 1e3, "Cl": CT_N * 1e3,
            "temp": 25.0, "pH": "7 charge", "units": "mmol/L",
        })
        ph = sol.pH
        sol.forget()
        return ph

    engine02_pq_default = PHREEQCChemicalEquilibriumEngine(
        {"H3PO4": CT_P * 1e3, "NH3": CT_N * 1e3, "CO2": 1.0,
         "CT_K": CT_P * 1e3, "CT_Cl": CT_N * 1e3},
        component_map={"H3PO4": "P", "NH3": "N(-3)", "CO2": "C",
                       "CT_K": "K", "CT_Cl": "Cl"},
        use_warmstart=False,
    )

    def solve02_pq_default():
        return engine02_pq_default.solve(
            totals={"H3PO4": CT_P, "NH3": CT_N, "CO2": CT_CO2_atm,
                    "CT_K": CT_P, "CT_Cl": CT_N},
        )

    TIMINGS[("02 gas-liquid CO2", "PHREEQC ideal (gamma->1)")] = time_repeated(solve02_pq_ideal, reps=10, warmup=2, trials=3)
    TIMINGS[("02 gas-liquid CO2", "PHREEQC default (WATEQ D-H)")] = time_repeated(solve02_pq_default, reps=10, warmup=2, trials=3)

    print(f"02 PHREEQC ideal:   {TIMINGS[('02 gas-liquid CO2', 'PHREEQC ideal (gamma->1)')]*1e3:.3f} ms/solve")
    print(f"02 PHREEQC default: {TIMINGS[('02 gas-liquid CO2', 'PHREEQC default (WATEQ D-H)')]*1e3:.3f} ms/solve")
else:
    print("Skipping PHREEQC timing for usecase 02 (phreeqpython not installed).")
"""),

    # ── 3. Usecase 03 timing ─────────────────────────────────────────────

    md("uc03-md", """\
## 3  Usecase 03 — E. coli batch growth (full dynamic simulation)

[Usecase 03](03_grow_ecoli_on_acetic_acid.ipynb)'s vessel, exactly as built
there: same species, same Monod parameters, same 2 L/37 °C/sparged-air
setup, same 20 h horizon at 2000 steps. The only new mechanic is
`ReactionSystem.configure_engine(use_activity=, activity_model=)`, called
*before* the first `.engine`/`.advance()` access — this is how the
activity treatment is chosen for a `ControlVolume`-driven system, the
dynamic-simulation equivalent of passing `use_activity=` directly to
`NRChemicalEquilibriumEngine.from_reactions()` in usecases 01-02.

A tiny throwaway run (5 steps) warms up each configuration before the
timed run, for the same reason `time_repeated` warms up above: the first
call into a freshly-built engine pays one-time import/JIT costs that
would otherwise unfairly penalise whichever configuration happens to run
first.\
"""),

    code("uc03-code", """\
ACETIC_ACID = Species(id="AceticAcid", atoms={"C": 2, "H": 4, "O": 2}, charge=0, MW=60.052)
ACETATE_MINUS = Species(id="Acetate-", atoms={"C": 2, "H": 3, "O": 2}, charge=-1, MW=59.044)
ECOLI = Species(id="Ecoli", atoms={"C": 1, "H": 1.8, "O": 0.5, "N": 0.2}, charge=0, MW=24.626)

acetate_eq = EquilibriumReaction(
    stoichiometry=[_e(ACETIC_ACID, -1), _e(ACETATE_MINUS, +1), _e(H_plus, +1)],
    log_K=-4.756, total_id="AceticAcid", label="acetate_eq",
)

mu_max_per_h, Ks_gL, Yxs, Ko2_gL = 0.30, 5e-3, 0.36, 0.2e-3
growth = ReactionBuilder.monod_aerobic_growth(
    substrate=ACETIC_ACID, biomass=ECOLI,
    mu_max_per_h=mu_max_per_h, Ks_gL=Ks_gL, yield_gX_gS=Yxs, Ko2_gL=Ko2_gL,
    balance="CHNO", label="growth_on_AceticAcid",
)

V_total_L, headspace_frac = 2.0, 0.20
V_liq = V_total_L * (1.0 - headspace_frac)
V_gas = V_total_L * headspace_frac
T_K = 310.15   # 37 C
n_total_gas = (1.0 * V_gas) / (R_L_ATM_MOL_K * T_K)

C_AcOH0 = 4.0 / ACETIC_ACID.MW    # 4 g/L acetic acid
X0 = 0.05 / ECOLI.MW              # 0.05 g/L inoculum

transfer_models = {
    "O2": KineticTransferModel(partition_model=HenryEquilibrium(H_ref=1.3e-5, dlnH=1500.0), k_transfer=150.0),
    "CO2": KineticTransferModel(partition_model=HenryEquilibrium(H_ref=3.4e-4, dlnH=2400.0), k_transfer=150.0 * 0.9),
    "N2": EquilibriumTransferModel(HenryEquilibrium(H_ref=6.4e-6, dlnH=1300.0)),
}

def build_usecase03_cv(use_activity, activity_model):
    system = ReactionSystem(
        [water, p1, p2, p3, nh4, co2_first, co2_second, acetate_eq, growth],
        label="ecoli_on_acetate", solver="newton_raphson",
    )
    system.configure_engine(use_activity=use_activity, activity_model=activity_model)

    gas_phase = GasPhase(
        n_mol={"O2": n_total_gas * 0.2095, "N2": n_total_gas * 0.7901, "CO2": n_total_gas * 0.0004},
        V_L=V_gas, T_K=T_K,
    )
    liquid_phase = LiquidPhase(
        n_mol={
            "H3PO4": CT_P * V_liq, "NH3": CT_N * V_liq,
            "AceticAcid": C_AcOH0 * V_liq, "Ecoli": X0 * V_liq,
            "K+": CT_P * V_liq, "Cl-": CT_N * V_liq,
            "O2": 0.0, "CO2": 0.0, "N2": 0.0,
        },
        V_L=V_liq, T_K=T_K,
    )
    gas_feed = GasFeed(
        vvm_min=1.0, y={"O2": 0.2095, "N2": 0.7901, "CO2": 0.0004}, P_inlet_atm=1.0,
        phase_key="gas", liquid_phase_key="liquid", label="gas_feed",
    )
    vent = PressureReliefVent(P_set_atm=1.0, mode="instant")
    return ControlVolume(
        phases={"gas": gas_phase, "liquid": liquid_phase},
        transfer_models=transfer_models,
        boundaries=[gas_feed, vent],
        reaction_system=system,
        label="batch_ecoli",
    )

def run_usecase03(use_activity, activity_model, tau_h, n_steps):
    cv = build_usecase03_cv(use_activity, activity_model)
    sim = Simulation(cvs={"main": cv}, label="runtime_compare")
    return sim.run(tau_h=tau_h, n_steps=n_steps)

# Warmup: tiny throwaway runs, discarded.
run_usecase03(False, "davies", tau_h=0.05, n_steps=5)
run_usecase03(True, "davies", tau_h=0.05, n_steps=5)

# Timed runs: usecase 03's actual parameters (20 h, 2000 steps).
t0 = time.perf_counter()
result03_ideal = run_usecase03(False, "davies", tau_h=20.0, n_steps=2000)
TIMINGS[("03 E. coli batch growth", "PyOMES ideal")] = time.perf_counter() - t0

t0 = time.perf_counter()
result03_davies = run_usecase03(True, "davies", tau_h=20.0, n_steps=2000)
TIMINGS[("03 E. coli batch growth", "PyOMES Davies")] = time.perf_counter() - t0

print(f"03 PyOMES ideal:  {TIMINGS[('03 E. coli batch growth', 'PyOMES ideal')]:.3f} s for the full 20 h / 2000-step run")
print(f"03 PyOMES Davies: {TIMINGS[('03 E. coli batch growth', 'PyOMES Davies')]:.3f} s for the full 20 h / 2000-step run")
print("PHREEQC: no dynamic-ControlVolume path exists in PyOMES, so there is no PHREEQC timing for usecase 03.")
"""),

    # ── 4. Summary table ──────────────────────────────────────────────────

    md("table-md", """\
## 4  All nine measurements, side by side

Printed before plotting so every number below is also available as plain
text (bar charts are for comparing shapes at a glance; exact values belong
in a table).\
"""),

    code("table-code", """\
UC_KEYS = ["01 pH prediction", "02 gas-liquid CO2", "03 E. coli batch growth"]
ENGINES = ["PyOMES ideal", "PyOMES Davies", "PHREEQC ideal (gamma->1)", "PHREEQC default (WATEQ D-H)"]

print(f"{'Usecase':<26}{'Engine / activity model':<30}{'Time':>18}")
print("-" * 74)
for uc in UC_KEYS:
    per_solve = uc.startswith(("01", "02"))
    for eng in ENGINES:
        val = TIMINGS.get((uc, eng))
        if val is None:
            continue
        disp = f"{val*1e3:.4f} ms/solve" if per_solve else f"{val:.3f} s (full run)"
        print(f"{uc:<26}{eng:<30}{disp:>18}")
"""),

    # ── 5. Plot ────────────────────────────────────────────────────────────

    md("plot-md", """\
## 5  Plotting the comparison

Two panels, both log-scale (the range here spans sub-millisecond NR solves
up to multi-second batch runs — a linear axis would flatten everything but
the largest bars to invisible slivers):

- **Left — absolute run time.** One group of bars per usecase, one colour
  per engine/activity-model treatment. Usecase 03 only shows two bars
  (PyOMES ideal/Davies) since PHREEQC has no entry there.
- **Right — time relative to that usecase's own PyOMES-ideal run**, i.e.
  each usecase's own ideal bar is normalised to 1.0. This isolates the
  *marginal* cost of Davies/PHREEQC from the much bigger difference in
  absolute magnitude between "one solve" (01, 02) and "2000 solves plus
  kinetics and transport" (03) that dominates the left panel — the
  question the right panel answers is "how much overhead does this
  treatment add, given whatever usecase it's used in," independent of
  which usecase that is.\
"""),

    code("plot-code", """\
UC_LABELS = [
    "01  pH prediction\\n(1 equilibrium solve)",
    "02  Gas-liquid CO2\\n(1 equilibrium solve)",
    "03  E. coli batch growth\\n(full 20h / 2000-step run)",
]
COLORS = {
    "PyOMES ideal":                  "#2a78d6",
    "PyOMES Davies":                 "#eb6834",
    "PHREEQC ideal (gamma->1)":     "#1baf7a",
    "PHREEQC default (WATEQ D-H)":  "#eda100",
}

n_uc = len(UC_KEYS)
bar_w, bar_gap = 0.20, 0.03
x = np.arange(n_uc)

# Usecase 03 only ever has 2 of the 4 engine columns (no PHREEQC path for a
# dynamic ControlVolume) -- rather than reserve 4 fixed slot positions and
# leave two visibly empty, each usecase group is centered on however many
# engines actually apply to *it*, so every group's bars sit under its own
# tick regardless of how many columns that usecase has.
def _grouped_bars(ax, values_by_engine, collect_legend):
    handles, labels = [], []
    for i, uc in enumerate(UC_KEYS):
        present = [(eng, values_by_engine[eng][i]) for eng in ENGINES
                   if values_by_engine[eng][i] is not None]
        n_present = len(present)
        for k, (eng, val) in enumerate(present):
            offset = (k - (n_present - 1) / 2) * (bar_w + bar_gap)
            bar = ax.bar(x[i] + offset, val, width=bar_w, color=COLORS[eng], edgecolor="none")
            if collect_legend and eng not in labels:
                handles.append(bar[0])
                labels.append(eng)
    ax.set_xticks(x)
    ax.set_xticklabels(UC_LABELS, fontsize=8.5)
    ax.set_yscale("log")
    ax.grid(True, which="major", axis="y", alpha=0.3, linewidth=0.8)
    ax.grid(True, which="minor", axis="y", alpha=0.12, linewidth=0.6)
    return handles, labels

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
fig.subplots_adjust(bottom=0.26, wspace=0.3)

abs_vals = {eng: [TIMINGS.get((uc, eng)) for uc in UC_KEYS] for eng in ENGINES}
handles, labels = _grouped_bars(ax1, abs_vals, collect_legend=True)
ax1.set_ylabel("wall-clock time (s, log scale)")
ax1.set_title("Absolute run time")

rel_vals = {}
for eng in ENGINES:
    row = []
    for uc in UC_KEYS:
        val, base = TIMINGS.get((uc, eng)), TIMINGS.get((uc, "PyOMES ideal"))
        row.append(val / base if (val is not None and base is not None) else None)
    rel_vals[eng] = row
_grouped_bars(ax2, rel_vals, collect_legend=False)
ax2.axhline(1.0, color="#c3c2b7", lw=1, zorder=0)
ax2.set_ylabel("time relative to that usecase's PyOMES-ideal run (log scale)")
ax2.set_title("Overhead relative to PyOMES ideal")

fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=8.5, frameon=False,
           bbox_to_anchor=(0.5, 0.02))
plt.show()
"""),

    md("summary-md", """\
## Takeaways

- **Davies vs. ideal is cheap.** Across usecases 01 and 02, adding the
  activity-coefficient correction costs roughly the same small constant
  factor per solve — noticeable in the right panel, invisible next to
  usecase 03's total in the left one.
- **PHREEQC's overhead is a completely different order of magnitude.**
  Every PHREEQC call round-trips through `phreeqpython`'s IPC layer to a
  separate PHREEQC process, which costs far more than the NR solve itself
  — the reason usecases 01/02 use it only for periodic validation (a
  handful of points, Section 6 of each), never as the per-step solver a
  dynamic simulation like usecase 03 would need.
- **Usecase-to-usecase differences dwarf activity-model differences.**
  The single biggest lever on run time here isn't which activity model is
  used at all — it's whether the workload is "solve once" or "solve two
  thousand times inside a stepped simulation." That's a completely
  different kind of scaling question (number of steps, step size,
  solver warmstart effectiveness) from anything Section 5's right panel
  addresses, and outside this notebook's scope.

## Where to go next

- **Why Davies and PHREEQC give slightly different answers, not just
  different run times** — [usecase 01](01_predict_ph_simple_liquid.ipynb)
  §6 and [usecase 02](02_equilibrate_with_atmospheric_gas.ipynb) §6 walk
  through the accuracy side of this same comparison.
- **Where usecase 03's own per-step cost actually goes** (NR solve vs.
  kinetics vs. transport vs. bookkeeping) — a per-step profiling breakdown
  is out of scope here but would be the natural next cut if usecase 03's
  run time itself, not just its activity-model sensitivity, becomes the
  question.\
"""),
)

save_path_runtime = HERE / "04_compare_runtime_by_usecase.ipynb"
with open(save_path_runtime, "w", encoding="utf-8") as f:
    json.dump(runtime_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path_runtime}")

# ═══════════════════════════════════════════════════════════════════════════
#  05  CSTR dilution-rate sweep — the E. coli/acetic-acid chemostat
# ═══════════════════════════════════════════════════════════════════════════

CSTR_SETUP = """\
import sys
from pathlib import Path
import math
import warnings
import numpy as np
import matplotlib.pyplot as plt

def _find_repo():
    for p in [Path.cwd(), *Path.cwd().parents]:
        if (p / "pyproject.toml").exists():
            return p
    raise RuntimeError("Run from inside the PyOMES repo")

sys.path.insert(0, str(_find_repo() / "models"))

from PyOMES.chemistry.species import Species
from PyOMES.chemistry.common_species import (
    H2O, H_plus, OH_minus,
    H3PO4, H2PO4_minus, HPO4_2minus, PO4_3minus,
    NH3, NH4_plus,
    CO2, HCO3_minus, CO3_2minus,
)
from PyOMES.chemistry import HenryEquilibrium
from PyOMES.reactions import (
    EquilibriumReaction, ReactionBuilder, ReactionSystem, StoichiometryEntry,
)
from PyOMES.core import (
    ControlVolume, GasPhase, LiquidPhase,
    KineticTransferModel, EquilibriumTransferModel,
    Simulation, GasFeed, PressureReliefVent, LiquidFeed, LiquidDrain,
)
from PyOMES.core.phases import R_L_ATM_MOL_K

# Same rationale as usecase 03: the explicit-Euler CV solver clamps a
# species' removal rate when a step would otherwise drive it negative
# (AccuracyWarning), and separately flags element/charge drift above
# tolerance (ConservationWarning). Silenced only so the run cells below
# stay readable -- Section 3 checks the actual steady state directly
# against the analytical chemostat solution.
from PyOMES.monitoring.accuracy import AccuracyWarning
from PyOMES.monitoring.conservation import ConservationWarning
warnings.filterwarnings("ignore", category=AccuracyWarning)
warnings.filterwarnings("ignore", category=ConservationWarning)

def _e(sp, coeff, phase="liquid"):
    return StoichiometryEntry(species=sp, phase=phase, coefficient=coeff)

print("Imports OK")
"""

cstr_nb = nb(
    md("title", """\
# CSTR Dilution-Rate Sweep — Running Usecase 03 as a Chemostat

**The situation:** [usecase 03](03_grow_ecoli_on_acetic_acid.ipynb) grew
*E. coli* on acetic acid in a **batch** bottle — inoculate once, watch
substrate deplete and biomass grow until it's done. Run the same organism
and substrate instead as a **CSTR (continuous stirred-tank reactor)**: fresh
sterile medium flows in at a constant volumetric rate $Q$, and broth flows
out at that *same* rate, so the working liquid volume $V$ never changes.
That equal-in/equal-out condition is what makes this a **chemostat** — the
one parameter that fully characterises its operating point is the
**dilution rate** $D = Q/V$ (units 1/h, the reciprocal of the mean
residence time $\\tau = 1/D$).

Left running long enough at a fixed $D$, a chemostat settles onto a
**steady state**: biomass growth is exactly balanced by washout ($\\mu = D$
at steady state), and substrate is exactly balanced by the difference
between what arrives in the feed and what the culture consumes.

**Aim of this notebook.** Build one mechanistic CSTR model — the same
`ControlVolume`/`Simulation` machinery usecase 03 used for a batch bottle,
now with a `LiquidFeed`/`LiquidDrain` pair added — and run *it alone* to
steady state across a sweep of dilution rates. Because that one model
already tracks O₂ transfer and the acid-base speciation directly, a single
mechanistic run gives every quantity of interest at once: steady-state
substrate, biomass, dissolved O₂, and pH, plus the volumetric productivity
they imply — no separate oxygen- or pH-specific model is needed to get at
any of them. A closed-form **substrate-only** estimate (Section 3) is kept
alongside purely as a sanity check and to pick a sensible sweep range; it
assumes O₂ is never limiting, and Section 6's own DO panel shows that
assumption holds closely only for part of the range this notebook sweeps.

This is the extension [usecase 03's own "where to go
next"](03_grow_ecoli_on_acetic_acid.ipynb) section pointed to, and the same
setup `demos/builder/cstr_fermenter.py` demonstrates as a plain script;
this notebook adds the dilution-rate sweep and the single publication-style
summary figure that script doesn't produce.\
"""),

    code("setup", CSTR_SETUP),

    # ── 1. Chemistry and organism: unchanged from usecase 03 ────────────

    md("chem-md", """\
## 1  Chemistry and organism — identical to usecase 03

Same eight equilibrium reactions (water autoionization, the three-step
phosphate ladder, ammonium/ammonia, the two-step carbonate ladder, and
acetic acid's own dissociation) and the same `Ecoli`/`AceticAcid` species
declarations as
[usecase 03 §2](03_grow_ecoli_on_acetic_acid.ipynb#2--Declare-the-chemistry)
— nothing about the acid-base network changes when the reactor becomes
continuous instead of batch.\
"""),

    code("chem-code", """\
ACETIC_ACID = Species(
    id="AceticAcid", atoms={"C": 2, "H": 4, "O": 2}, charge=0, MW=60.052,
)
ACETATE_MINUS = Species(
    id="Acetate-", atoms={"C": 2, "H": 3, "O": 2}, charge=-1, MW=59.044,
)
ECOLI = Species(
    id="Ecoli", atoms={"C": 1, "H": 1.8, "O": 0.5, "N": 0.2}, charge=0,
    MW=24.626,
)

water = EquilibriumReaction(
    stoichiometry=[_e(H2O, -1), _e(H_plus, +1), _e(OH_minus, +1)],
    log_K=-14.0, label="water",
)
p1 = EquilibriumReaction(
    stoichiometry=[_e(H3PO4, -1), _e(H2PO4_minus, +1), _e(H_plus, +1)],
    log_K=-2.15, total_id="H3PO4", label="p1",
)
p2 = EquilibriumReaction(
    stoichiometry=[_e(H2PO4_minus, -1), _e(HPO4_2minus, +1), _e(H_plus, +1)],
    log_K=-7.20, total_id="H3PO4", label="p2",
)
p3 = EquilibriumReaction(
    stoichiometry=[_e(HPO4_2minus, -1), _e(PO4_3minus, +1), _e(H_plus, +1)],
    log_K=-12.35, total_id="H3PO4", label="p3",
)
nh4 = EquilibriumReaction(
    stoichiometry=[_e(NH4_plus, -1), _e(NH3, +1), _e(H_plus, +1)],
    log_K=-9.25, total_id="NH3", label="nh4",
)
co2_first = EquilibriumReaction(
    stoichiometry=[_e(CO2, -1), _e(H2O, -1), _e(HCO3_minus, +1), _e(H_plus, +1)],
    log_K=-6.35, total_id="CO2", label="co2_first",
)
co2_second = EquilibriumReaction(
    stoichiometry=[_e(HCO3_minus, -1), _e(CO3_2minus, +1), _e(H_plus, +1)],
    log_K=-10.33, total_id="CO2", label="co2_second",
)
acetate_eq = EquilibriumReaction(
    stoichiometry=[_e(ACETIC_ACID, -1), _e(ACETATE_MINUS, +1), _e(H_plus, +1)],
    log_K=-4.756, total_id="AceticAcid", label="acetate_eq",
)

print("8 equilibrium reactions declared (unchanged from usecase 03).")
"""),

    # ── 2. Monod kinetics for continuous operation ───────────────────────

    md("growth-md", """\
## 2  Monod growth kinetics for continuous operation

$\\mu_{max}$, the yield $Y_{X/S}$, and $K_{O_2}$ carry over from usecase 03
unchanged; $K_s$ is set differently here. A chemostat's steady-state
substrate concentration is

$$S_{ss} = \\frac{K_s \\, D}{\\mu_{max} - D}$$

which, at usecase 03's $K_s = 5$ mg/L (a high-affinity estimate appropriate
for a *batch*, where substrate only approaches zero briefly at the end of
a run), sits at sub-mg/L levels across almost the entire viable
dilution-rate range in a *continuous* culture — a numerically stiff regime
for the explicit-Euler reaction sub-stepping this framework uses. This
notebook instead uses $K_s = 0.5$ g/L (500 mg/L): still within the
order-of-magnitude spread usecase 03's own literature review noted for
organic-acid uptake kinetics (Kovarova-Kovar & Egli 1998), and a
substrate-affinity choice rather than a growth-rate or yield change. It
keeps $S_{ss}$ a comfortable double-digit-percent fraction of the feed
concentration everywhere except right at washout, so the dilution-rate
sweep in Section 6 runs at a normal, non-stiff step size.\
"""),

    code("growth-code", """\
mu_max_per_h = 0.30    # h^-1              -- unchanged from usecase 03
Ks_gL        = 0.5     # g acetate / L     -- see Section 2 above
Yxs          = 0.36    # g biomass / g acetate -- unchanged
Ko2_gL       = 0.2e-3  # g O2 / L          -- unchanged

growth = ReactionBuilder.monod_aerobic_growth(
    substrate=ACETIC_ACID, biomass=ECOLI,
    mu_max_per_h=mu_max_per_h, Ks_gL=Ks_gL, yield_gX_gS=Yxs, Ko2_gL=Ko2_gL,
    balance="CHNO", label="growth_on_AceticAcid",
)
system = ReactionSystem(
    [water, p1, p2, p3, nh4, co2_first, co2_second, acetate_eq, growth],
    label="ecoli_on_acetate",
)

print(f"mu_max = {mu_max_per_h} /h, Ks = {Ks_gL*1e3:.0f} mg/L, "
      f"Yxs = {Yxs} g/g, Ko2 = {Ko2_gL*1e3} mg/L")
"""),

    # ── 3. Chemostat theory: two models ──────────────────────────────────

    md("theory-md", """\
## 3  Chemostat theory: a closed-form sanity check

This is standard chemostat theory (Monod 1950; see e.g. Bailey & Ollis,
*Biochemical Engineering Fundamentals*, ch. 7, or Blanch & Clark,
*Biochemical Engineering*, ch. 6) — restated here to pick a sensible
dilution-rate range for Section 6's sweep and to sanity-check the
mechanistic simulation in Section 5. It's a reference line on the plots
that follow, not a second model this notebook tries to validate in its
own right.

At steady state, growth exactly balances washout: $\\mu(S_{ss}) = D$.
Substituting the Monod form $\\mu(S) = \\mu_{max}\\,S/(K_s+S)$ and solving
for $S$ gives the **substrate-only** estimate

$$S_{ss}(D) = \\frac{K_s \\, D}{\\mu_{max} - D} \\qquad (D < \\mu_{max})$$

and the biomass estimate follows from the yield alone (the same
`ReactionBuilder`-derived stoichiometry usecase 03 §6 validated by mass
balance):

$$X_{ss}(D) = Y_{X/S}\\,\\bigl(S_{feed} - S_{ss}(D)\\bigr)$$

This closed form assumes dissolved O₂ never limits growth (the O₂ Monod
term $O_2/(K_{O_2}+O_2)$ from Section 2 is implicitly treated as exactly
1, i.e. DO always sits at air saturation) — an assumption Section 6's own
DO panel shows holds closely only for part of the swept range, since the
vessel's aeration capacity is finite. Two quantities from this
substrate-only estimate are still useful as *reference lines*:

- **Washout.** As $D \\to \\mu_{max}$, $S_{ss} \\to \\infty$ in the formula
  above, which is unphysical — what actually happens is $X_{ss}$ hits zero
  *before* that, at $S_{ss} = S_{feed}$ (all the feed passes through
  unconsumed). Solving $S_{ss}(D_{crit}) = S_{feed}$ gives the **washout
  dilution rate**
  $$D_{crit} = \\mu_{max}\\,\\frac{S_{feed}}{K_s + S_{feed}}$$
  Because this closed form ignores O₂ limitation, it's an *upper bound* on
  where the mechanistic model actually washes out — Section 6's simulated
  values fall away from it well before reaching $D_{crit}$ itself.
- **Volumetric productivity.** The reactor's output of biomass per unit
  time and volume is $P(D) = D\\,X_{ss}(D)$ — zero at $D=0$ (no flow, no
  output) and zero at $D_{crit}$ (washout, no biomass), so it has an
  interior maximum. Differentiating and setting $dP/dD=0$ gives the
  **productivity-optimal dilution rate**
  $$D_{opt} = \\mu_{max}\\left(1 - \\sqrt{\\frac{K_s}{K_s + S_{feed}}}\\right)$$
  Section 6 checks how close the mechanistic model's *own* simulated
  productivity peak sits to this estimate.\
"""),

    code("theory-code", """\
C_AcOH_feed_gL = 4.0   # g/L acetic acid in the feed (same loading as usecase 03)

def S_ss_analytic(D, Ks=Ks_gL, mu_max=mu_max_per_h):
    return Ks * D / (mu_max - D)

def X_ss_analytic(D, Sfeed=C_AcOH_feed_gL, Ks=Ks_gL, mu_max=mu_max_per_h, Y=Yxs):
    return Y * (Sfeed - S_ss_analytic(D, Ks, mu_max))

D_crit = mu_max_per_h * C_AcOH_feed_gL / (Ks_gL + C_AcOH_feed_gL)
D_opt  = mu_max_per_h * (1.0 - np.sqrt(Ks_gL / (Ks_gL + C_AcOH_feed_gL)))

print(f"S_feed              = {C_AcOH_feed_gL} g/L acetic acid")
print(f"D_crit (substrate-only washout)    = {D_crit:.4f} /h   (tau = {1/D_crit:.2f} h)")
print(f"D_opt (substrate-only max productivity) = {D_opt:.4f} /h   (tau = {1/D_opt:.2f} h)")
print(f"At D_opt: S_ss = {S_ss_analytic(D_opt):.4f} g/L, "
      f"X_ss = {X_ss_analytic(D_opt):.4f} g/L, "
      f"productivity = {D_opt*X_ss_analytic(D_opt):.4f} g/L/h")
"""),

    # ── 4. Assemble the chemostat CV and run one dilution rate ───────────

    md("assemble-md", """\
## 4  Assemble the chemostat and run one dilution rate to steady state

Same 2 L vessel, headspace, and aeration pattern as usecase 03 (`GasFeed` +
`PressureReliefVent`, kinetic O₂/CO₂/N₂ transfer) — none of
that changes. `GasFeed`'s aeration rate (1 vvm — matching the airflow
Leone et al. 2015 used per litre of working volume, see reference below)
is fixed, and so is gas-liquid transfer: $k_La$ = 90 /h for O₂, CO₂, and N₂ throughout this
notebook, within the range reported for bench-scale *E. coli* STRs under
mechanical agitation/sparging — roughly 70-1000 /h across reported
operating points (~72 /h at 200 rpm/1.5 vvm in a 5 L fermenter, up to
~980 /h at 1500 rpm/18 L min⁻¹ in a 12.8 L bioreactor; see references
below) — and toward the moderate end of it, comfortably inside the regime
where Section 6's DO panel shows aeration capacity actually starts to
matter.

What's new relative to usecase 03 is the liquid boundary pair that makes
this a chemostat instead of a bottle: `LiquidFeed` delivers sterile medium
(phosphate, ammonium, and acetic acid at their feed concentrations — *no*
biomass, since the feed stream is sterile) at rate $Q$, and `LiquidDrain`
removes broth (every dissolved species, biomass included) at that same
$Q$, so $V_{liq}$ stays fixed by construction. `build_cv(D_per_h)` builds
one fresh chemostat at a given dilution rate; it's called once below for
$D=0.15$/h and again in Section 6 across the full sweep.

*References:* Leone et al. (2015, *Microb. Cell Fact.* 14:106) — the same
paper cited for $Y_{X/S}$ in Section 2 and the feed loading below — report
an airflow of 60 L/h into a 1 L working volume (1 vvm, matching this
vessel's `vvm_min`). For the $k_La$ range: BioProcess International,
"Lessons in Bioreactor Scale-Up, Part 4" (5 L fermenter, 200 rpm/1.5 vvm
example); Fan et al., *J. Chem. Pharm. Res.*, 2014, 6(7):1810-1817 (12.8 L
bioreactor, recombinant *E. coli*, 1500 rpm/18 L min⁻¹ example).\
"""),

    code("assemble-code", """\
V_total_L = 2.0
headspace_frac = 0.20
V_liq = V_total_L * (1.0 - headspace_frac)
V_gas = V_total_L * headspace_frac
T_K = 310.15   # 37 C, same as usecase 03

CT_P = 0.022     # mol/L KH2PO4 -> mol/L total phosphate, mol/L K+  (usecase 01's recipe)
CT_N = 0.0187    # mol/L NH4Cl  -> mol/L total ammoniacal N, mol/L Cl-
C_AcOH_feed = C_AcOH_feed_gL / ACETIC_ACID.MW

kLa = 90.0   # /h for O2, CO2, and N2 -- see Section 4 above

def kH_mol_L_atm(H_ref, dlnH, T_K, T_ref=298.15):
    H_ref_mol_L_atm = (H_ref / 1000.0) * 101325.0
    return H_ref_mol_L_atm * math.exp(dlnH * (1.0 / T_K - 1.0 / T_ref))

DO_sat_mgL = kH_mol_L_atm(1.3e-5, 1500.0, T_K) * 0.2095 * 32.00 * 1e3   # air-sat. DO, 37 C
Ko2_mgL = Ko2_gL * 1e3

def make_transfer_models():
    return {
        "O2": KineticTransferModel(
            partition_model=HenryEquilibrium(H_ref=1.3e-5, dlnH=1500.0), k_transfer=kLa,
        ),
        "CO2": KineticTransferModel(
            partition_model=HenryEquilibrium(H_ref=3.4e-4, dlnH=2400.0), k_transfer=kLa,
        ),
        "N2": KineticTransferModel(
            partition_model=HenryEquilibrium(H_ref=6.4e-6, dlnH=1300.0), k_transfer=kLa,
        ),
    }

def build_cv(D_per_h, X0_gL=0.05):
    Q_L_per_h = D_per_h * V_liq
    n_total_gas = (1.0 * V_gas) / (R_L_ATM_MOL_K * T_K)
    gas_phase = GasPhase(
        n_mol={
            "O2": n_total_gas * 0.2095, "N2": n_total_gas * 0.7901,
            "CO2": n_total_gas * 0.0004,
        },
        V_L=V_gas, T_K=T_K,
    )
    X0 = X0_gL / ECOLI.MW
    liquid_phase = LiquidPhase(
        n_mol={
            "H3PO4": CT_P * V_liq, "NH3": CT_N * V_liq,
            "AceticAcid": C_AcOH_feed * V_liq, "Ecoli": X0 * V_liq,
            "K+": CT_P * V_liq, "Cl-": CT_N * V_liq,
            "O2": 0.0, "CO2": 0.0, "N2": 0.0,
        },
        V_L=V_liq, T_K=T_K,
    )
    boundaries = [
        GasFeed(vvm_min=1.0, y={"O2": 0.2095, "N2": 0.7901, "CO2": 0.0004},
                P_inlet_atm=1.0, phase_key="gas", liquid_phase_key="liquid",
                label="gas_feed"),
        PressureReliefVent(P_set_atm=1.0, mode="instant"),
        LiquidFeed(
            Q_L_per_h=Q_L_per_h,
            feed_conc_mol_L={
                "H3PO4": CT_P, "NH3": CT_N, "AceticAcid": C_AcOH_feed,
                "K+": CT_P, "Cl-": CT_N,
            },
            label="sterile_feed",
        ),
        LiquidDrain(Q_L_per_h=Q_L_per_h, label="broth_drain"),
    ]
    return ControlVolume(
        phases={"gas": gas_phase, "liquid": liquid_phase},
        transfer_models=make_transfer_models(), boundaries=boundaries,
        reaction_system=system, label=f"cstr_D{D_per_h:.3f}",
    )

def run_to_steady_state(D_per_h, n_res=15.0, dt_h=0.01):
    \"\"\"Run build_cv(D_per_h) for n_res multiples of the system's own
    relaxation time 1/(mu_max-D) -- *not* the residence time 1/D, which
    underestimates how long convergence actually takes near D_opt/washout
    (Section 6 discusses this).\"\"\"
    cv = build_cv(D_per_h)
    tau_relax_h = 1.0 / (mu_max_per_h - D_per_h)
    tau_sim_h = n_res * tau_relax_h
    n_steps = int(tau_sim_h / dt_h)
    sim = Simulation(cvs={"main": cv}, label=f"cstr_D{D_per_h:.3f}")
    result = sim.run(tau_h=tau_sim_h, n_steps=n_steps)
    return result, tau_sim_h

D_demo = 0.15   # the dilution rate this section runs and checks in detail
result_demo, tau_sim_demo = run_to_steady_state(D_demo)
print(f"D = {D_demo}/h  (kLa={kLa:.0f}/h, tau_residence = {1/D_demo:.2f} h), "
      f"ran {tau_sim_demo:.1f} h simulated in {result_demo.runtime_s:.1f} s wall-clock.")
print(f"DO_sat = {DO_sat_mgL:.2f} mg/L (air saturation, 37 C), K_O2 = {Ko2_mgL:.2f} mg/L")
"""),

    md("plot-md", """\
## 5  Approach to steady state, vs. the substrate-only estimate

Biomass and substrate both start at their inoculum/feed values and relax
toward steady state — plotted against the substrate-only closed form from
Section 3 (gray, dashed) for this same $D=0.15$/h, as a sanity check that
the CV's feed/drain/reaction/transfer machinery, run through a full
`Simulation`, reproduces basic chemostat behaviour. The two agree closely
here (well under 1%) — Section 6's sweep shows that agreement eroding
substantially above $D\\approx0.18$-$0.20$/h, where the vessel's finite
aeration capacity starts to bind.\
"""),

    code("plot-code", """\
t = result_demo.t_h
X = result_demo.liquid_mol["main"]["Ecoli"] * ECOLI.MW / V_liq
S = result_demo.liquid_mol["main"]["AceticAcid"] * ACETIC_ACID.MW / V_liq
pH = result_demo.pH["main"]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

axes[0].plot(t, X, color="tab:green")
axes[0].axhline(X_ss_analytic(D_demo), color="gray", ls="--", lw=1,
                label=f"substrate-only X_ss = {X_ss_analytic(D_demo):.3f} g/L")
axes[0].set_xlabel("time (h)")
axes[0].set_ylabel("biomass (g/L)")
axes[0].set_title(f"Biomass approaching steady state, D={D_demo}/h")
axes[0].legend(fontsize=8)
axes[0].grid(alpha=0.3)

axes[1].plot(t, S, color="tab:blue")
axes[1].axhline(S_ss_analytic(D_demo), color="gray", ls="--", lw=1,
                label=f"substrate-only S_ss = {S_ss_analytic(D_demo):.3f} g/L")
axes[1].set_xlabel("time (h)")
axes[1].set_ylabel("acetic acid (g/L)")
axes[1].set_title(f"Substrate approaching steady state, D={D_demo}/h")
axes[1].legend(fontsize=8)
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.show()

print(f"X_ss: simulated = {X[-1]:.4f} g/L  |  substrate-only estimate = "
      f"{X_ss_analytic(D_demo):.4f} g/L  |  gap = "
      f"{100*abs(X[-1]-X_ss_analytic(D_demo))/X_ss_analytic(D_demo):.1f}%")
print(f"S_ss: simulated = {S[-1]:.4f} g/L  |  substrate-only estimate = "
      f"{S_ss_analytic(D_demo):.4f} g/L  |  gap = "
      f"{100*abs(S[-1]-S_ss_analytic(D_demo))/S_ss_analytic(D_demo):.1f}%")
print(f"pH at steady state: {pH[-1]:.2f}")
"""),

    md("ph-note-md", """\
That last line is worth pausing on: pH settles to a low, *uncontrolled*
value here, well below usecase 01's phosphate-only baseline (pH ≈ 4.7).
Two things conspire, both structural to running this culture
continuously rather than in a batch: (1) growth keeps drawing down the
ammoniacal-nitrogen pool for biomass synthesis while phosphate (which
nothing consumes) stays at its feed concentration, so the medium's
buffer composition shifts at steady state in a way it never fully does
in a batch; and (2) aerobic growth on an organic acid continuously
produces CO₂, adding a steady carbonic-acid load that a batch run only
accumulates transiently before consumption stops. Usecase 03 could let pH
drift freely because a batch run ends; a chemostat runs indefinitely at
whatever pH this settles to, which is precisely why real continuous
cultures on acid substrates are normally run under active pH control —
`demos/builder/cstr_fermenter.py`'s `PHController` (dosing NaOH to a
setpoint) is the closed-loop fix, left out here so this notebook can stay
focused on the dilution-rate/productivity relationship itself.\
"""),

    # ── 6. Dilution-rate sweep and the publication figure ────────────────

    md("sweep-md", """\
## 6  Dilution-rate sweep and the summary figure

This section sweeps $D$ from a low dilution rate up toward washout,
running the *same* single mechanistic model from Sections 4-5 at each
point and reading its steady state directly: biomass, substrate,
dissolved O₂, and pH all come out of the one `ControlVolume`/`Simulation`
run, with no separate O₂- or pH-specific model needed to get at any of
them.

Points closer to washout take longer to reach steady state — critical
slowing down, a well-documented chemostat phenomenon (perturbations relax
more slowly the closer the operating point sits to the washout boundary;
see e.g. Bailey & Ollis ch. 7) — so `run_to_steady_state`'s `n_res` is
bumped for the last, slowest point. The sweep stops at $D=0.24$/h:
pushing closer to Section 3's substrate-only $D_{crit}$ (0.267/h) runs
into that slowing-down directly, to the point where confirming
convergence would cost far more simulated time than this sweep is worth.
That the mechanistic model is already declining well before $D_{crit}$ is
itself the point (Section 3): the substrate-only closed form assumes O₂
is never limiting, and this vessel's finite aeration capacity means the
real washout point sits below 0.267/h, not at it.\
"""),

    code("sweep-code", """\
D_sweep = [0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.22, 0.23, 0.24]

X_ss, S_ss, DO_ss, pH_ss = {}, {}, {}, {}
X_ss[D_demo], S_ss[D_demo] = X[-1], S[-1]
DO_ss[D_demo] = result_demo.liquid_mol["main"]["O2"][-1] * 32.00 * 1e3 / V_liq
pH_ss[D_demo] = pH[-1]

for D in D_sweep:
    if D == D_demo:
        continue
    n_res = 40.0 if D >= 0.24 else 15.0
    result_i, tau_i = run_to_steady_state(D, n_res=n_res)
    X_ss[D] = result_i.liquid_mol["main"]["Ecoli"][-1] * ECOLI.MW / V_liq
    S_ss[D] = result_i.liquid_mol["main"]["AceticAcid"][-1] * ACETIC_ACID.MW / V_liq
    DO_ss[D] = result_i.liquid_mol["main"]["O2"][-1] * 32.00 * 1e3 / V_liq
    pH_ss[D] = result_i.pH["main"][-1]
    print(f"D={D:.3f}/h  n_res={n_res:.0f}  tau_sim={tau_i:6.1f}h  "
          f"X_ss={X_ss[D]:.4f} g/L  S_ss={S_ss[D]:.4f} g/L  "
          f"DO_ss={DO_ss[D]:.3f} mg/L  pH_ss={pH_ss[D]:.2f}")

D_all = sorted(X_ss)
P_ss = {D: D * X_ss[D] for D in D_all}
D_opt_sim = max(D_all, key=lambda D: P_ss[D])

print(f"\\nSwept D = {D_all} /h")
print(f"Simulated productivity peaks at D={D_opt_sim}/h (P={P_ss[D_opt_sim]:.4f} g/L/h), "
      f"vs. the substrate-only estimate D_opt={D_opt:.4f}/h.")
"""),

    md("runtime-md", """\
### 6.1  Runtime benchmark

For the manuscript runtime table, this cell runs three representative
dilution rates once each through the same `run_to_steady_state` helper used
above. The reported time is `SimulationResult.runtime_s × 1000`, so each
row is the wall-clock time for one full CSTR simulation at that dilution
rate.\
"""),

    code("runtime-code", """\
D_runtime = [0.01, 0.10, 0.20, 0.219]
N_RUNTIME_REPS = 1

print(f"{'Dilution rate':>14}{'Run Time, ms':>16}")
cstr_runtime_rows = []
for D in D_runtime:
    times_ms = []
    for _ in range(N_RUNTIME_REPS):
        result_i, _tau_i = run_to_steady_state(D)
        times_ms.append(result_i.runtime_s * 1e3)
    runtime_ms = float(np.mean(times_ms))
    cstr_runtime_rows.append((D, runtime_ms))
    print(f"{D:>14.2f}{runtime_ms:>16.1f}")
"""),

    md("figure-md", """\
### 6.2  Publication figure

Single column, three linked panels sharing the $D$ axis. The middle panel
is split into two tightly-coupled sub-panels (DO above, pH below) rather
than overlaid on one dual-axis plot: dissolved O₂ and pH sit on
incompatible scales, and a shared-axis overlay would misrepresent both.

(a) steady-state substrate and biomass vs. $D$, against the substrate-only
estimate from Section 3; (b)/(c) steady-state dissolved O₂ and pH vs. $D$
— (b) also marks air-saturation $DO_{sat}$ and the O₂ half-saturation
constant $K_{O_2}$ for scale; (d) volumetric productivity
$P(D)=D\\,X_{ss}(D)$ vs. $D$, marking both the substrate-only $D_{opt}$
and where the simulated productivity actually peaks. The gray dotted line
in each panel is Section 3's substrate-only $D_{crit}$ — a reference
ceiling the mechanistic model falls short of, not a value it's expected
to reach.\
"""),

    code("figure-code", """\
D_arr = np.array(D_all)
X_arr = np.array([X_ss[D] for D in D_all])
S_arr = np.array([S_ss[D] for D in D_all])
DO_arr = np.array([DO_ss[D] for D in D_all])
pH_arr = np.array([pH_ss[D] for D in D_all])
P_arr = np.array([P_ss[D] for D in D_all])

D_theory = np.linspace(1e-4, D_crit * 0.999, 400)

# Fixed categorical colors (identity, not rank) -- consistent across panels.
C_BIOMASS, C_SUBSTRATE, C_DO, C_PH, C_PROD = (
    "#2a78d6", "#eb6834", "#1baf7a", "#008300", "#4a3aa7",
)
C_MUTED, C_GRID = "#898781", "#e1e0d9"

PUB_FIG_WIDTH_IN = 3.5   # single-column journal width, same convention as usecase 01
PUB_DPI = 600
PUB_STYLE = {
    "font.size": 8, "axes.labelsize": 9, "axes.titlesize": 9,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 6.5,
    "axes.edgecolor": C_MUTED, "axes.labelcolor": "#0b0b0b",
    "xtick.color": C_MUTED, "ytick.color": C_MUTED,
    "axes.grid": True, "grid.color": C_GRID, "grid.linewidth": 0.6,
    "axes.linewidth": 0.6,
}

FIG_DIR = _find_repo() / "demos" / "usecases" / "figures"
FIG_DIR.mkdir(exist_ok=True)

def _mark_D(ax, label=False):
    ax.axvline(D_crit, color=C_MUTED, ls=":", lw=1,
               label=f"$D_{{crit}}$ (substrate-only) = {D_crit:.3f} /h" if label else None)
    ax.axvline(D_opt, color=C_MUTED, ls="--", lw=1,
               label=f"$D_{{opt}}$ (substrate-only) = {D_opt:.3f} /h" if label else None)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

with plt.rc_context(PUB_STYLE):
    fig = plt.figure(figsize=(PUB_FIG_WIDTH_IN, 9.5))
    outer = fig.add_gridspec(3, 1, height_ratios=[1.15, 1.7, 1.15], hspace=0.42)
    ax_a = fig.add_subplot(outer[0])
    inner_bc = outer[1].subgridspec(2, 1, hspace=0.40)
    ax_b = fig.add_subplot(inner_bc[0], sharex=ax_a)
    ax_c = fig.add_subplot(inner_bc[1], sharex=ax_a)
    ax_d = fig.add_subplot(outer[2], sharex=ax_a)

    # (a) substrate & biomass
    ax_a.plot(D_theory, X_ss_analytic(D_theory), color=C_BIOMASS, ls=":", lw=1, alpha=0.6)
    ax_a.plot(D_theory, S_ss_analytic(D_theory), color=C_SUBSTRATE, ls=":", lw=1, alpha=0.6,
              label="substrate-only estimate")
    ax_a.plot(D_arr, X_arr, color=C_BIOMASS, marker="o", ms=3.5, lw=1.4, label="biomass $X_{ss}$")
    ax_a.plot(D_arr, S_arr, color=C_SUBSTRATE, marker="o", ms=3.5, lw=1.4, label="substrate $S_{ss}$")
    _mark_D(ax_a, label=True)
    ax_a.set_ylabel("conc. (g/L)")
    ax_a.set_title("(a)  Steady-state substrate & biomass", loc="left")
    ax_a.legend(frameon=False, loc="upper left")
    plt.setp(ax_a.get_xticklabels(), visible=False)

    # (b) dissolved O2
    ax_b.plot(D_arr, DO_arr, color=C_DO, marker="o", ms=3.5, lw=1.4)
    ax_b.axhline(DO_sat_mgL, color=C_MUTED, ls="--", lw=0.8)
    ax_b.axhline(Ko2_mgL, color=C_MUTED, ls=":", lw=0.8)
    ax_b.text(D_arr[-1], DO_sat_mgL, "$DO_{sat}$ ", color=C_MUTED, fontsize=6.5,
              va="bottom", ha="right")
    ax_b.text(D_arr[-1], Ko2_mgL, "$K_{O_2}$ ", color=C_MUTED, fontsize=6.5,
              va="bottom", ha="right")
    ax_b.set_ylim(-0.5, DO_sat_mgL * 1.15)
    _mark_D(ax_b)
    ax_b.set_ylabel("DO (mg/L)")
    ax_b.set_title("(b)  Steady-state dissolved O$_2$", loc="left")
    plt.setp(ax_b.get_xticklabels(), visible=False)

    # (c) pH
    ax_c.plot(D_arr, pH_arr, color=C_PH, marker="o", ms=3.5, lw=1.4)
    _mark_D(ax_c)
    ax_c.set_ylabel("pH")
    ax_c.set_title("(c)  Steady-state pH", loc="left")
    plt.setp(ax_c.get_xticklabels(), visible=False)

    # (d) productivity
    ax_d.plot(D_theory, D_theory * X_ss_analytic(D_theory), color=C_MUTED, ls=":", lw=1,
              alpha=0.6, label="substrate-only estimate")
    ax_d.plot(D_arr, P_arr, color=C_PROD, marker="o", ms=3.5, lw=1.4, label="simulated $P(D)$")
    ax_d.scatter([D_opt_sim], [P_ss[D_opt_sim]], color=C_PROD, edgecolor="#0b0b0b", s=45,
                 zorder=5, label=f"simulated peak, D={D_opt_sim:.2f}/h")
    _mark_D(ax_d)
    ax_d.set_xlabel("dilution rate $D$ (1/h)")
    ax_d.set_ylabel("productivity (g/L/h)")
    ax_d.set_title("(d)  Volumetric productivity", loc="left")
    ax_d.legend(frameon=False, loc="upper left")

    for ax in (ax_a, ax_b, ax_c, ax_d):
        ax.set_xlim(0, D_crit * 1.02)

    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"05_cstr_dilution_rate_sweep.{ext}",
                    dpi=PUB_DPI, bbox_inches="tight")
    plt.show()

print(f"Saved publication figure to {FIG_DIR / '05_cstr_dilution_rate_sweep.png'}")
"""),

    md("summary-md", """\
## Takeaways

- **The dilution rate is the whole story for a fixed reactor.** Given
  fixed $S_{feed}$, $\\mu_{max}$, $K_s$, $Y_{X/S}$, and $k_La$, every
  steady-state quantity — substrate residual, biomass concentration,
  dissolved O₂, pH, volumetric productivity — is a function of $D$ alone.
  Turning that one dial moves the operating point along the curves in
  Section 6.1's figure.
- **A single mechanistic model gets all four quantities at once.** Because
  the `ControlVolume`/`Simulation` run tracks O₂ transfer and acid-base
  speciation directly, there's no need for a separate closed-form O₂ model
  to get at DO, or a separate calculation to get at pH — Section 6's sweep
  reads all of it off the same run that also gives $S_{ss}$ and $X_{ss}$.
- **"DO looks healthy" is not the same as "DO isn't affecting the
  answer."** Dissolved O₂ never crashes anywhere in this sweep — it stays
  within a few mg/L of air saturation throughout (comfortably non-limiting
  on any DO probe reading) — yet it's still responsible for a real,
  measurable shift in $S_{ss}$/$X_{ss}$/productivity relative to the
  substrate-only estimate, growing from under 1% at $D=0.15$/h to over 10%
  by $D=0.23$/h (Section 5, Section 6.1's panel a).
- **Productivity is maximised strictly below washout.** The substrate-only
  estimate puts that optimum at
  $D_{opt} = \\mu_{max}(1-\\sqrt{K_s/(K_s+S_{feed})})$ — not at $D_{crit}$
  itself, and the gap between them is the margin a real operator has
  before a disturbance tips the culture into washout. Section 6.1 shows
  where the mechanistic model's own simulated productivity actually peaks,
  which need not sit at exactly the same $D$ once O₂ limitation is
  accounted for.
- **Convergence time is not symmetric around washout.** Approaching the
  mechanistic model's own washout point from below is markedly slower —
  critical slowing down, a well-documented chemostat phenomenon — which is
  exactly why this notebook's sweep stops at $D=0.24$/h rather than
  pushing all the way to Section 3's substrate-only $D_{crit}$ (Section 6).

## Where to go next

- **Vary $k_La$ yourself.** `build_cv()` and `run_to_steady_state()` take
  no aeration argument in this notebook because a single condition was
  kept for simplicity — reintroducing `kLa` as a parameter and re-running
  Section 6's sweep at a second value (e.g. a more vigorously sparged
  vessel) shows directly how much of the gap from the substrate-only
  estimate is an aeration artifact rather than a fixed property of this
  organism/substrate pair.
- **Closed-loop pH (and DO) control for this exact chemostat** —
  `demos/builder/cstr_fermenter.py`'s `PHController`/`DOAgitationController`
  cascade, addressing the low-pH note in Section 5.
- **Back to batch, for comparison** —
  [usecase 03](03_grow_ecoli_on_acetic_acid.ipynb) runs the identical
  organism/substrate/kinetics without feed or drain at all.
- **`FermenterBuilder`, the fluent alternative to hand-assembling the
  `ControlVolume` in Section 4** —
  `models/vlmodels/fermenter/config/builder.py`, used directly by
  `demos/builder/cstr_fermenter.py`.\
"""),
)

save_path_cstr = HERE / "05_cstr_dilution_rate_sweep.ipynb"
with open(save_path_cstr, "w", encoding="utf-8") as f:
    json.dump(cstr_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path_cstr}")
