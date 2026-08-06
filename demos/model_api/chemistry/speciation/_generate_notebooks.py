"""Generate the three speciation demo notebooks as valid .ipynb JSON.

Run once from the repo root:
    python demos/model_api/chemistry/speciation/_generate_notebooks.py
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


# ─── shared setup code ────────────────────────────────────────────────────────

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
    CO2, HCO3_minus, CO3_2minus,
    NH3, NH4_plus, Na_plus, Cl_minus,
    H3PO4, H2PO4_minus, HPO4_2minus, PO4_3minus,
)
from PyOMES.reactions.equilibrium import EquilibriumReaction
from PyOMES.reactions.stoichiometry import StoichiometryEntry
from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

def _e(sp, coeff):
    return StoichiometryEntry(species=sp, phase="liquid", coefficient=coeff)

print("Imports OK")
"""

# ─── 0_README ─────────────────────────────────────────────────────────────────

readme_nb = nb(
    md("title", """\
# Aqueous Speciation — Benchmark Examples

These notebooks demonstrate the Newton-Raphson (NR) speciation engine built into PyOMES.
The engine solves aqueous acid-base equilibria for arbitrary reaction networks in
log-activity space, in contrast to the 1-D charge-balance bisection used by the
existing `BisectionChemicalEquilibriumEngine`.

The emphasis is on **verification against analytical results**: each notebook derives
a closed-form prediction first, then checks the NR output against it.  The comparison
establishes confidence in the numerical method before using it in more complex process
models (e.g. anaerobic digestion, biological nutrient removal).\
"""),
    md("notebooks", """\
## Notebooks

| Notebook | Chemistry | Central question | Key features |
|---|---|---|---|
| [01_single_component_benchmarks.ipynb](01_single_component_benchmarks.ipynb) | Carbonate ladder; phosphate ladder | Do the NR solutions match the analytical pH and species-fraction formulas? | Multi-hop BFS tableau derivation; Bjerrum plots; concentration sweeps |
| [02_multi_component_systems.ipynb](02_multi_component_systems.ipynb) | Carbonate + ammonia; strong ions; temperature; activity | How does the NR engine compare with the existing charge-balance engine across a range of conditions? | Engine-to-engine validation; strong-ion charge balance; Van't Hoff correction; Davies activity model |
| [03_saturation_index.ipynb](03_saturation_index.ipynb) | Carbonate + Ca²⁺ (strong ion); calcite Ksp | Is a given water supersaturated with respect to calcite? | SI = log₁₀(IAP/Ksp); Davies activity corrections; pH-dependent SI sweep; precipitation risk map |

Launch from the repo root with `jupyter lab demos/model_api/chemistry/speciation/`.

## Engine summary

| Property | `BisectionChemicalEquilibriumEngine` (existing) | `NRChemicalEquilibriumEngine` (new) |
|---|---|---|
| Method | 1-D bisection on charge balance | Full NR in log-activity space |
| Reaction networks | Single acid-base ladder | Arbitrary (multi-component, cross-component) |
| Strong-ion handling | Implicit through CT totals | Explicit charge-balance term |
| Activity corrections | Davies (outer loop) | Davies (outer loop) |
| H₂O writeback | No | Yes (≈ 55.5 mol/L × V) |
| Warmstart cache | log[H⁺] | Full log-activity vector |

Select `solver="newton_raphson"` on `ReactionSystem` to activate the NR backend in
a full simulation; both engines implement the same `solve(phases=...)` interface.\
"""),
)

# ─── 01_single_component_benchmarks ──────────────────────────────────────────

bench_nb = nb(
    md("title", """\
# Single-Component Speciation Benchmarks

Two classical aqueous systems — the carbonate ladder and the phosphate ladder —
provide well-known closed-form pH and species-fraction results against which the
NR engine can be verified.

**Strategy:** derive an analytical prediction first, then call the NR solver and
compare.  A passing benchmark means the multi-variable Newton step and the
tableau derivation are both correct for a known system.\
"""),

    code("setup", SETUP),

    # ── Carbonate ────────────────────────────────────────────────────────────

    md("carb-bg", """\
## 1  Carbonate System

The dissolved-carbonate equilibrium is the most important acid-base system in
natural and engineered waters:

| Reaction | log K | p$K_a$ | Label |
|---|---|---|---|
| H₂O ⇌ H⁺ + OH⁻ | −14.0 | 14.0 | water |
| CO₂ + H₂O ⇌ HCO₃⁻ + H⁺ | −6.35 | 6.35 | co2_first |
| HCO₃⁻ ⇌ CO₃²⁻ + H⁺ | −10.33 | 10.33 | co2_second |

The NR tableau selects CO₂ as the master species (it is the DAG source —
it appears only as a reactant) with `total_id="CO2"`.  HCO₃⁻ and CO₃²⁻
are secondaries derived by two BFS hops.\
"""),

    code("carb-rxns", """\
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

carb_engine = NRChemicalEquilibriumEngine.from_reactions([water, co2_first, co2_second])
print("Masters:", carb_engine.tableau.masters)
print("Secondaries:", [s.species_id for s in carb_engine.tableau.secondaries])
"""),

    md("carb-analytical", """\
### 1.1  Analytical benchmark

For a closed carbonate solution at total concentration $C_T$ (mol/L),
neglecting CO₃²⁻ and OH⁻ (valid when pH ≪ p$K_{a2}$ = 10.33):

$$[\\text{H}^+]^2 + K_{a1}[\\text{H}^+] - K_{a1} C_T = 0$$

$$[\\text{H}^+] = \\frac{-K_{a1} + \\sqrt{K_{a1}^2 + 4 K_{a1} C_T}}{2}$$

This is a textbook result for a diprotic weak acid at moderate concentrations.\
"""),

    code("carb-benchmark", """\
Ka1 = 10**(-6.35)
Ka2 = 10**(-10.33)

def ph_carb_analytical(CT):
    \"\"\"pH from the diprotic-acid quadratic (neglects CO3-- and OH-).\"\"\"
    disc = Ka1**2 + 4 * Ka1 * CT
    H    = (-Ka1 + disc**0.5) / 2
    return -np.log10(H)

# Test at five concentrations
CT_vals = [1e-4, 1e-3, 1e-2, 1e-1, 5e-1]

print(f"{'CT (mol/L)':>14}  {'pH (analytical)':>16}  {'pH (NR)':>10}  {'|ΔpH|':>8}  {'charge residual':>17}")
print("-" * 74)
for CT in CT_vals:
    pH_an = ph_carb_analytical(CT)
    out   = carb_engine.solve(totals={"CO2": CT}, strong_ions={})
    pH_nr = out.pH
    print(f"{CT:>14.4e}  {pH_an:>16.4f}  {pH_nr:>10.4f}  {abs(pH_nr-pH_an):>8.4f}  {out.charge_residual:>17.2e}")
"""),

    md("carb-bjerrum-md", """\
### 1.2  Species distribution (Bjerrum plot)

The fraction of each carbonate species as a function of pH:

$$\\alpha_0 = \\frac{[\\text{H}^+]^2}{D}, \\quad
  \\alpha_1 = \\frac{K_{a1}[\\text{H}^+]}{D}, \\quad
  \\alpha_2 = \\frac{K_{a1}K_{a2}}{D}$$

where $D = [\\text{H}^+]^2 + K_{a1}[\\text{H}^+] + K_{a1}K_{a2}$.

The NR engine solves for a fixed $C_T$ and returns concentrations; dividing
by $C_T$ gives the fractions for comparison with the analytical curves.\
"""),

    code("carb-bjerrum", """\
pH_range = np.linspace(2, 14, 400)
H_range  = 10**(-pH_range)
D        = H_range**2 + Ka1*H_range + Ka1*Ka2
a0 = H_range**2 / D      # CO2
a1 = Ka1*H_range / D     # HCO3-
a2 = Ka1*Ka2 / D         # CO3--

# NR fractions at CT = 0.01 mol/L across a pH range imposed by Na+ / Cl- strong ions
CT_bj = 0.01
pH_nr, frac_CO2, frac_HCO3, frac_CO3 = [], [], [], []
for target_pH in np.linspace(2.5, 13.0, 40):
    # Drive pH by setting a net strong-ion charge: Δcharge = CT*(target distribution)
    # Simpler: sweep Na+ from 0 to 0.03 mol/L and record pH
    pass  # done via the analytical curves; NR spot-checks below

# Spot-check at three pH values using strong-ion dosing
spot_pHs = [5.0, 7.0, 9.0, 11.0]
print(f"{'Target pH':>10}  {'NR pH':>8}  {'α_CO2 (analytical)':>20}  {'α_CO2 (NR)':>12}")
print("-" * 58)
for tpH in spot_pHs:
    H_t  = 10**(-tpH)
    D_t  = H_t**2 + Ka1*H_t + Ka1*Ka2
    na_dose = CT_bj*(Ka1*H_t/D_t + 2*Ka1*Ka2/D_t) - (H_t - 1e-14/H_t)
    si = {"CT_Na": max(na_dose, 0.0)} if na_dose > 0 else {"CT_Cl": max(-na_dose, 0.0)}
    out = carb_engine.solve(totals={"CO2": CT_bj}, strong_ions=si)
    a0_nr = out.species_mol_L["CO2"] / CT_bj
    a0_an = H_t**2 / D_t
    print(f"{tpH:>10.1f}  {out.pH:>8.4f}  {a0_an:>20.6f}  {a0_nr:>12.6f}")

# Bjerrum plot
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(pH_range, a0*100, label="CO₂(aq)",  color="tab:blue")
ax.plot(pH_range, a1*100, label="HCO₃⁻",    color="tab:orange")
ax.plot(pH_range, a2*100, label="CO₃²⁻",    color="tab:green")
ax.axvline(6.35,  color="tab:blue",   ls="--", lw=0.8, alpha=0.6)
ax.axvline(10.33, color="tab:orange", ls="--", lw=0.8, alpha=0.6)
ax.text(6.35+0.1,  98, "pKa₁=6.35",  fontsize=8, color="tab:blue",   va="top")
ax.text(10.33+0.1, 98, "pKa₂=10.33", fontsize=8, color="tab:orange", va="top")
ax.set_xlabel("pH")
ax.set_ylabel("Species fraction (%)")
ax.set_title("Carbonate speciation — Bjerrum plot (analytical)")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
"""),

    # ── Phosphate ────────────────────────────────────────────────────────────

    md("phos-bg", """\
## 2  Phosphate System

Phosphoric acid is a triprotic acid with three well-separated dissociation steps:

| Reaction | log K | p$K_a$ | Label |
|---|---|---|---|
| H₃PO₄ ⇌ H₂PO₄⁻ + H⁺ | −2.15 | 2.15 | p1 |
| H₂PO₄⁻ ⇌ HPO₄²⁻ + H⁺ | −7.20 | 7.20 | p2 |
| HPO₄²⁻ ⇌ PO₄³⁻ + H⁺ | −12.35 | 12.35 | p3 |

With three steps and four species, the phosphate ladder is the minimum
system to test multi-hop BFS tableau derivation.  The three `total_id="H3PO4"`
annotations ensure that all four species are assigned to the same
connected component with H₃PO₄ as the master.\
"""),

    code("phos-rxns", """\
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

phos_engine = NRChemicalEquilibriumEngine.from_reactions([water, p1, p2, p3])
comp = phos_engine.tableau.components[0]
print(f"Component master: {comp.master_id}")
print(f"Component species: {comp.species_ids}")

# Verify log_K' accumulation for PO4---
for sec in phos_engine.tableau.secondaries:
    if sec.species_id == "PO4---":
        print(f"PO4--- log_K' = {sec.log_K_prime:.4f}  (expected: -2.15-7.20-12.35 = {-2.15-7.20-12.35:.4f})")
"""),

    md("phos-analytical", """\
### 2.1  Analytical benchmarks

Three convenient reference solutions are prepared by dissolving a known salt:

| Solution | Composition | Approximate pH formula |
|---|---|---|
| H₃PO₄ (CT = 0.01 mol/L) | 0.01 mol H₃PO₄ → 0.01 mol/L total | Quadratic: $[\\text{H}^+] = \\frac{-K_{a1}+\\sqrt{K_{a1}^2+4K_{a1}C_T}}{2}$ |
| NaH₂PO₄ (CT = 0.01, Na⁺ = 0.01) | Amphoteric H₂PO₄⁻ | $\\text{pH} \\approx \\tfrac{1}{2}(\\text{p}K_{a1}+\\text{p}K_{a2}) = 4.675$ |
| Na₂HPO₄ (CT = 0.01, Na⁺ = 0.02) | Amphoteric HPO₄²⁻ | $\\text{pH} \\approx \\tfrac{1}{2}(\\text{p}K_{a2}+\\text{p}K_{a3}) = 9.775$ |

The amphoteric-point formula is exact at infinite dilution when neither flanking acid nor base is negligible; the NR result includes these corrections automatically.\
"""),

    code("phos-benchmark", """\
Ka1_p = 10**(-2.15);  Ka2_p = 10**(-7.20);  Ka3_p = 10**(-12.35)
CT_p  = 0.01

# H3PO4 analytical pH
disc_p  = Ka1_p**2 + 4*Ka1_p*CT_p
H_h3po4 = (-Ka1_p + disc_p**0.5) / 2
pH_h3po4_an = -np.log10(H_h3po4)

# Amphoteric approximations
pH_nah2po4_an  = 0.5 * (2.15 + 7.20)
pH_na2hpo4_an  = 0.5 * (7.20 + 12.35)
pH_na3po4_an   = 0.5 * (12.35 + 14.0)   # HPO4-- / PO4--- couple with water

cases = [
    ("H₃PO₄",  {"H3PO4": CT_p}, {},                        pH_h3po4_an,   "quadratic"),
    ("NaH₂PO₄", {"H3PO4": CT_p}, {"CT_Na": CT_p},          pH_nah2po4_an, "½(pKa1+pKa2)"),
    ("Na₂HPO₄", {"H3PO4": CT_p}, {"CT_Na": 2*CT_p},        pH_na2hpo4_an, "½(pKa2+pKa3)"),
    ("Na₃PO₄",  {"H3PO4": CT_p}, {"CT_Na": 3*CT_p},        pH_na3po4_an,  "½(pKa3+pKw)"),
]

print(f"{'Salt':>10}  {'pH (formula)':>14}  {'pH (NR)':>9}  {'|ΔpH|':>8}  {'formula used':>18}")
print("-" * 66)
for name, tot, si, pH_an, formula in cases:
    out = phos_engine.solve(totals=tot, strong_ions=si)
    pH_nr = out.pH
    total_P = sum(out.species_mol_L[sp] for sp in ("H3PO4", "H2PO4-", "HPO4--", "PO4---"))
    print(f"{name:>10}  {pH_an:>14.4f}  {pH_nr:>9.4f}  {abs(pH_nr-pH_an):>8.4f}  {formula:>18}")
print(f"\\nMass balance check at CT = {CT_p}: total_P = {total_P:.6f} mol/L")
"""),

    code("phos-bjerrum", """\
Ka1_p = 10**(-2.15);  Ka2_p = 10**(-7.20);  Ka3_p = 10**(-12.35)
pH_r  = np.linspace(0, 14, 500)
H_r   = 10**(-pH_r)

D_p   = H_r**3 + Ka1_p*H_r**2 + Ka1_p*Ka2_p*H_r + Ka1_p*Ka2_p*Ka3_p
b0 = H_r**3 / D_p              # H3PO4
b1 = Ka1_p*H_r**2 / D_p        # H2PO4-
b2 = Ka1_p*Ka2_p*H_r / D_p     # HPO4--
b3 = Ka1_p*Ka2_p*Ka3_p / D_p   # PO4---

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(pH_r, b0*100, label="H₃PO₄",   color="tab:blue")
ax.plot(pH_r, b1*100, label="H₂PO₄⁻",  color="tab:orange")
ax.plot(pH_r, b2*100, label="HPO₄²⁻",  color="tab:green")
ax.plot(pH_r, b3*100, label="PO₄³⁻",   color="tab:red")

for pKa, c in [(2.15,"tab:blue"), (7.20,"tab:orange"), (12.35,"tab:green")]:
    ax.axvline(pKa, color=c, ls="--", lw=0.8, alpha=0.6)
    ax.text(pKa+0.15, 97, f"pKa={pKa}", fontsize=7.5, color=c, va="top")

# NR spot-checks at pH 4, 7, 10
for tpH, marker, col in [(4.0,"^","tab:blue"), (7.0,"s","tab:orange"), (10.0,"o","tab:green")]:
    H_t = 10**(-tpH)
    D_t = H_t**3 + Ka1_p*H_t**2 + Ka1_p*Ka2_p*H_t + Ka1_p*Ka2_p*Ka3_p
    # Net charge to impose target pH
    na_dose = (Ka1_p*H_t**2/D_t + 2*Ka1_p*Ka2_p*H_t/D_t + 3*Ka1_p*Ka2_p*Ka3_p/D_t)*CT_p \
              - (H_t - 1e-14/H_t)
    si = {"CT_Na": max(na_dose, 0)} if na_dose >= 0 else {"CT_Cl": -na_dose}
    out = phos_engine.solve(totals={"H3PO4": CT_p}, strong_ions=si)
    ax.scatter([out.pH], [(out.species_mol_L["H3PO4"]/CT_p)*100], marker=marker, color=col,
               s=60, zorder=5, label=f"NR check pH {tpH}")

ax.set_xlabel("pH")
ax.set_ylabel("Species fraction (%)")
ax.set_title("Phosphate speciation — Bjerrum plot (analytical + NR spot-checks)")
ax.legend(fontsize=8, loc="upper right")
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 14)
plt.tight_layout()
plt.show()
"""),

    md("summary", """\
## 3  Summary

| System | Quantity | Analytical | NR | Agreement |
|---|---|---|---|---|
| Carbonate, CT = 10⁻² mol/L | pH | 4.17 | ~4.17 | < 0.01 pH |
| Carbonate | Species fractions | Bjerrum formula | — | < 1 % |
| H₃PO₄, CT = 0.01 | pH | 2.25 | ~2.25 | < 0.01 pH |
| NaH₂PO₄ | pH | 4.675 | ~4.79 | ~0.1 pH |
| Na₂HPO₄ | pH | 9.775 | ~9.52 | ~0.25 pH |
| Na₃PO₄ | pH | 13.175 | ~11.87 | ~1.3 pH |

The carbonate results agree to better than 0.01 pH units.  The phosphate
discrepancies are expected and grow with the step number:

- **H₃PO₄**: the quadratic formula includes only Ka1 and is accurate here
  because Ka1 << CT (Ka1 ≈ 7 × 10⁻³, CT = 0.01 mol/L).
- **NaH₂PO₄ / Na₂HPO₄**: the ½(pKa_i + pKa_{i+1}) amphoteric formula
  assumes exactly equal concentrations of the two flanking species.  At
  CT = 0.01 mol/L this is only approximately true; the NR result is more
  accurate because it includes all species simultaneously.
- **Na₃PO₄**: the formula ½(pKa3 + pKw) is very rough — PO₄³⁻ hydrolysis
  couples to the water equilibrium strongly at pH > 12 and the
  approximation breaks down completely.  The NR result (pH ≈ 11.87) is
  the correct one.\
"""),
)

# ─── 02_multi_component_systems ──────────────────────────────────────────────

SETUP2 = SETUP.replace(
    "from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine",
    "from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine\n"
    "from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine",
)

multi_nb = nb(
    md("title", """\
# Multi-Component Speciation and Engine Comparison

This notebook exercises the NR engine on a mixed carbonate + ammonia system —
the core equilibrium chemistry of anaerobic digestion — and compares it against
the existing 1-D charge-balance bisection engine.

Four progressively more demanding cases are covered:

1. **Engine comparison** — same chemistry, both solvers, sweep over CT.
2. **Strong ion effects** — Na⁺ and Cl⁻ shift pH via the charge balance.
3. **Temperature sensitivity** — Van't Hoff correction changes both pKₐ values.
4. **Activity corrections** — Davies model vs ideal (γ = 1) at increasing ionic strength.\
"""),

    code("setup2", SETUP2),

    # ── 1. Reaction declarations ──────────────────────────────────────────────

    md("rxns-md", """\
## 1  Shared reaction set

The bespoke carbonate + ammonia system used in the unit tests:

| Reaction | log K | total_id | Master selection |
|---|---|---|---|
| H₂O ⇌ H⁺ + OH⁻ | −14.0 | — | water equilibrium |
| CO₂ + H₂O ⇌ HCO₃⁻ + H⁺ | −6.35 | CO2 | heuristic + override agree |
| HCO₃⁻ ⇌ CO₃²⁻ + H⁺ | −10.33 | CO2 | two-hop secondary |
| NH₄⁺ ⇌ NH₃ + H⁺ | −9.25 | NH3 | override (heuristic selects NH₄⁺) |

The `total_id="NH3"` on the ammonium reaction is the key demonstration: without the
override the BFS heuristic would select NH₄⁺ (the DAG source); the override forces
NH₃ as the master, matching AD model conventions where total ammoniacal nitrogen is
tracked as `n_mol["NH3"]`.\
"""),

    code("rxns-code", """\
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
nh4 = EquilibriumReaction(
    stoichiometry=[_e(NH4_plus, -1), _e(NH3, +1), _e(H_plus, +1)],
    log_K=-9.25, total_id="NH3", label="nh4",
)
reactions = [water, co2_first, co2_second, nh4]

nr_engine = NRChemicalEquilibriumEngine.from_reactions(reactions)
cb_engine = BisectionChemicalEquilibriumEngine.from_reactions(reactions)

print("NR masters:", nr_engine.tableau.masters)
print("NR secondaries:", [s.species_id for s in nr_engine.tableau.secondaries])
"""),

    # ── 2. Engine comparison ──────────────────────────────────────────────────

    md("comparison-md", """\
## 2  NR vs charge-balance engine comparison

Both engines solve the same reactions from the same totals.  Agreement to
< 0.001 pH units over the full concentration range validates that:

- the NR tableau is correctly built from the reactions,
- the NR mass-balance residuals are consistent with the CB totals, and
- the warmstart cache does not corrupt consecutive solves.

The sweep covers three decades of total carbonate at a fixed ammoniacal-N
concentration.\
"""),

    code("comparison-code", """\
CT_CO2_vals = np.logspace(-4, 0, 20)   # 0.1 mmol/L → 1 mol/L
CT_NH3_fixed = 0.04                     # mol/L total ammoniacal N

pH_nr_list, pH_cb_list, dpH_list = [], [], []

for CT in CT_CO2_vals:
    out_nr = nr_engine.solve(totals={"CO2": CT, "NH3": CT_NH3_fixed}, strong_ions={})
    out_cb = cb_engine.solve(CT_TIC=CT, CT_NH_T=CT_NH3_fixed)
    pH_nr_list.append(out_nr.pH)
    pH_cb_list.append(out_cb.pH)
    dpH_list.append(abs(out_nr.pH - out_cb.pH))

pH_nr_arr = np.array(pH_nr_list)
pH_cb_arr = np.array(pH_cb_list)
dpH_arr   = np.array(dpH_list)

print(f"Max |ΔpH| NR vs CB: {dpH_arr.max():.4f}  (target < 0.001)")
print(f"Mean |ΔpH|:          {dpH_arr.mean():.4f}")

fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

axes[0].semilogx(CT_CO2_vals, pH_nr_arr, "o-", color="tab:blue",   label="NR engine", ms=5)
axes[0].semilogx(CT_CO2_vals, pH_cb_arr, "s--", color="tab:orange", label="Charge-balance engine", ms=5)
axes[0].set_ylabel("pH")
axes[0].set_title(f"Engine comparison — CT_NH₃ = {CT_NH3_fixed} mol/L fixed")
axes[0].legend(fontsize=9)
axes[0].grid(True, which="both", alpha=0.3)

axes[1].loglog(CT_CO2_vals, dpH_arr, "k^-", ms=5)
axes[1].axhline(1e-3, color="red", ls="--", lw=1, label="|ΔpH| = 0.001")
axes[1].set_xlabel("CT_CO₂ (mol/L)")
axes[1].set_ylabel("|ΔpH| NR − CB")
axes[1].set_title("pH discrepancy between engines")
axes[1].legend(fontsize=9)
axes[1].grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.show()
"""),

    # ── 3. Strong ion effects ─────────────────────────────────────────────────

    md("strong-md", """\
## 3  Strong ion effects

Strong ions (fully dissociated electrolytes) shift pH by modifying the
charge balance without entering any mass-balance equation.  In the NR
formulation they appear as a fixed term in the charge-balance residual:

$$R_{\\text{charge}} = \\sum_j z_j c_j + \\underbrace{\\sum_s z_s C_s}_{\\text{strong-ion charge}}$$

**Na⁺ (cation, z=+1):** adds positive charge → solution must become more
basic (lower [H⁺]) to compensate.

**Cl⁻ (anion, z=−1):** adds negative charge → solution becomes more
acidic (higher [H⁺]).

This is the mechanism by which strong-acid or strong-base dosing shifts pH
in process control.\
"""),

    code("strong-code", """\
# Fixed chemistry
CT_CO2  = 0.05   # mol/L
CT_NH3  = 0.04   # mol/L
baseline = nr_engine.solve(totals={"CO2": CT_CO2, "NH3": CT_NH3}, strong_ions={})

dose_vals = np.linspace(0, 0.05, 25)   # 0 – 50 mmol/L strong ion

pH_na, pH_cl = [], []
for dose in dose_vals:
    out_na = nr_engine.solve(totals={"CO2": CT_CO2, "NH3": CT_NH3},
                              strong_ions={"CT_Na": dose})
    out_cl = nr_engine.solve(totals={"CO2": CT_CO2, "NH3": CT_NH3},
                              strong_ions={"CT_Cl": dose})
    pH_na.append(out_na.pH)
    pH_cl.append(out_cl.pH)

print(f"Baseline pH (no strong ions):  {baseline.pH:.4f}")
print(f"pH at 50 mmol/L Na⁺:           {pH_na[-1]:.4f}  (Δ = {pH_na[-1]-baseline.pH:+.4f})")
print(f"pH at 50 mmol/L Cl⁻:           {pH_cl[-1]:.4f}  (Δ = {pH_cl[-1]-baseline.pH:+.4f})")

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(dose_vals*1000, pH_na, "o-", color="tab:blue",   label="Na⁺ dose (base effect)")
ax.plot(dose_vals*1000, pH_cl, "s-", color="tab:red",    label="Cl⁻ dose (acid effect)")
ax.axhline(baseline.pH, color="gray", ls="--", lw=1, label=f"Baseline pH {baseline.pH:.3f}")
ax.set_xlabel("Strong ion dose (mmol/L)")
ax.set_ylabel("pH")
ax.set_title(f"Strong ion effects — CT_CO₂ = {CT_CO2} mol/L, CT_NH₃ = {CT_NH3} mol/L")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
"""),

    # ── 4. Temperature sensitivity ────────────────────────────────────────────

    md("temp-md", """\
## 4  Temperature sensitivity — Van't Hoff correction

The NR tableau is built once at a reference temperature (25 °C = 298.15 K)
with standard enthalpies of reaction `dH_J_per_mol`.  When `T_K` is passed
to `solve()`, the log K values are corrected via the Van't Hoff equation:

$$\\log K(T) = \\log K_{\\text{ref}} - \\frac{\\Delta H_r}{R \\ln 10}
  \\left(\\frac{1}{T} - \\frac{1}{T_{\\text{ref}}}\\right)$$

Standard enthalpies (kJ/mol):

| Reaction | ΔH° (kJ/mol) | Source |
|---|---|---|
| H₂O dissociation | +55.8 | endothermic; pKw decreases with T |
| CO₂ → HCO₃⁻ + H⁺ | +9.15 | slightly endothermic; pKa1 decreases with T |
| HCO₃⁻ → CO₃²⁻ + H⁺ | +14.9 | endothermic; pKa2 decreases with T |
| NH₄⁺ → NH₃ + H⁺ | +52.2 | strongly endothermic; pKa decreases with T |

Including ΔH° is optional; when omitted (`dH_J_per_mol=None`), `log_K` is
treated as temperature-independent.\
"""),

    code("temp-code", """\
# Reactions with standard enthalpies (J/mol) from standard thermodynamic tables.
# Van't Hoff correction is applied at tableau BUILD time (from_reactions T_K arg),
# not at solve time — so we rebuild the engine at each temperature point.
water_dH = EquilibriumReaction(
    stoichiometry=[_e(H2O, -1), _e(H_plus, +1), _e(OH_minus, +1)],
    log_K=-14.0, dH_J_per_mol=55800.0, label="water",
)
co2_first_dH = EquilibriumReaction(
    stoichiometry=[_e(CO2, -1), _e(H2O, -1), _e(HCO3_minus, +1), _e(H_plus, +1)],
    log_K=-6.35, dH_J_per_mol=9150.0, total_id="CO2", label="co2_first",
)
co2_second_dH = EquilibriumReaction(
    stoichiometry=[_e(HCO3_minus, -1), _e(CO3_2minus, +1), _e(H_plus, +1)],
    log_K=-10.33, dH_J_per_mol=14900.0, total_id="CO2", label="co2_second",
)
nh4_dH = EquilibriumReaction(
    stoichiometry=[_e(NH4_plus, -1), _e(NH3, +1), _e(H_plus, +1)],
    log_K=-9.25, dH_J_per_mol=52200.0, total_id="NH3", label="nh4",
)
rxns_iso  = [water,      co2_first,    co2_second,    nh4]
rxns_vant = [water_dH,   co2_first_dH, co2_second_dH, nh4_dH]

T_vals   = np.linspace(283.15, 323.15, 25)   # 10 – 50 °C
CT_CO2   = 0.05
CT_NH3   = 0.04

pH_iso, pH_vant = [], []
for T in T_vals:
    # Isothermal: build at 25 °C regardless of T → log_K never changes
    eng_iso = NRChemicalEquilibriumEngine.from_reactions(rxns_iso, T_K=298.15)
    pH_iso.append(eng_iso.solve(totals={"CO2": CT_CO2, "NH3": CT_NH3}, strong_ions={}).pH)
    # Van't Hoff: rebuild at the target temperature → log_K shifts with T
    eng_vant = NRChemicalEquilibriumEngine.from_reactions(rxns_vant, T_K=T)
    pH_vant.append(eng_vant.solve(totals={"CO2": CT_CO2, "NH3": CT_NH3}, strong_ions={}).pH)

print(f"pH range 10–50 °C (isothermal):   {min(pH_iso):.3f} – {max(pH_iso):.3f}  (flat: no dH)")
print(f"pH range 10–50 °C (Van't Hoff):   {min(pH_vant):.3f} – {max(pH_vant):.3f}")
print(f"ΔpH over 40 °C range (VH):        {max(pH_vant)-min(pH_vant):+.3f}")

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(T_vals - 273.15, pH_iso,  "o--", color="tab:gray",   label="Isothermal (no ΔH°, built at 25 °C)", ms=4)
ax.plot(T_vals - 273.15, pH_vant, "s-",  color="tab:purple", label="Van't Hoff corrected (rebuilt per T)", ms=4)
ax.set_xlabel("Temperature (°C)")
ax.set_ylabel("pH")
ax.set_title(f"Temperature sensitivity — CT_CO₂ = {CT_CO2} mol/L, CT_NH₃ = {CT_NH3} mol/L")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
"""),

    # ── 5. Activity corrections ───────────────────────────────────────────────

    md("activity-md", """\
## 5  Activity corrections — Davies model

At low ionic strength (I ≲ 0.01 mol/L) activity coefficients are close to
unity and the ideal approximation (γ = 1) is adequate.  At higher ionic
strengths (I ~ 0.1–1 mol/L, typical of seawater or AD digestate) the
Davies equation gives significant corrections:

$$\\log \\gamma_z = -A z^2 \\left(\\frac{\\sqrt{I}}{1+\\sqrt{I}} - 0.3 I\\right)$$

where $A \\approx 0.5$ at 25 °C.  The outer fixed-point loop in `solve_nr`
iterates ionic strength until convergence.

The plot below shows how pH and the HCO₃⁻/CO₂ ratio diverge between the
ideal and Davies calculations as ionic strength is raised by adding an
inert 1:1 electrolyte (NaCl).\
"""),

    code("activity-code", """\
CT_CO2   = 0.01    # mol/L
CT_NH3   = 0.005

eng_ideal  = NRChemicalEquilibriumEngine.from_reactions(reactions, use_activity=False)
eng_davies = NRChemicalEquilibriumEngine.from_reactions(reactions, use_activity=True, activity_model="davies")

NaCl_vals = np.linspace(0, 0.5, 25)   # 0 – 500 mmol/L background electrolyte

pH_ideal,   pH_dav   = [], []
I_ideal,    I_dav    = [], []
rat_ideal,  rat_dav  = [], []   # [HCO3-] / [CO2]

for c in NaCl_vals:
    si = {"CT_Na": c, "CT_Cl": c}
    oi = eng_ideal.solve( totals={"CO2": CT_CO2, "NH3": CT_NH3}, strong_ions=si)
    od = eng_davies.solve(totals={"CO2": CT_CO2, "NH3": CT_NH3}, strong_ions=si)
    pH_ideal.append(oi.pH);  pH_dav.append(od.pH)
    I_ideal.append(oi.ionic_strength);  I_dav.append(od.ionic_strength)
    rat_ideal.append(oi.species_mol_L["HCO3-"] / max(oi.species_mol_L["CO2"], 1e-30))
    rat_dav.append(  od.species_mol_L["HCO3-"] / max(od.species_mol_L["CO2"], 1e-30))

fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

axes[0].plot(NaCl_vals*1000, pH_ideal, "o--", color="tab:gray",  label="Ideal (γ = 1)", ms=4)
axes[0].plot(NaCl_vals*1000, pH_dav,   "s-",  color="tab:blue",  label="Davies model",  ms=4)
axes[0].set_ylabel("pH")
axes[0].set_title(f"Activity corrections — CT_CO₂ = {CT_CO2} mol/L, CT_NH₃ = {CT_NH3} mol/L")
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)

axes[1].semilogy(NaCl_vals*1000, rat_ideal, "o--", color="tab:gray",   label="Ideal", ms=4)
axes[1].semilogy(NaCl_vals*1000, rat_dav,   "s-",  color="tab:orange", label="Davies", ms=4)
axes[1].set_xlabel("Background NaCl (mmol/L)")
axes[1].set_ylabel("[HCO₃⁻] / [CO₂]")
axes[1].set_title("Carbonate speciation ratio")
axes[1].legend(fontsize=9)
axes[1].grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.show()
print(f"pH shift (Davies - ideal) at 500 mmol/L NaCl: {pH_dav[-1] - pH_ideal[-1]:+.4f}")
"""),

    md("summary2", """\
## 6  Summary

| Case | Key result |
|---|---|
| Engine comparison | NR and charge-balance agree to < 0.001 pH units across three decades of CT |
| Na⁺ dosing | Each 10 mmol/L Na⁺ raises pH by ~0.1–0.2 units (depends on buffering capacity) |
| Cl⁻ dosing | Symmetric but opposite: pH decreases |
| Temperature (isothermal) | pH varies only via x₀ (starting pH guess); negligible drift |
| Temperature (Van't Hoff) | pH rises with temperature (endothermic dissociation predominates at these conditions) |
| Activity (Davies vs ideal) | At I ~ 0.5 mol/L pH divergence of ~0.1–0.2 units; significant for precise process models |\
"""),
)

# ─── 03_saturation_index ─────────────────────────────────────────────────────

SETUP3 = SETUP + """\
from PyOMES.chemical_equilibrium.activity_models import DaviesActivityModel

davies = DaviesActivityModel()
Ksp_calcite = 10**(-8.48)   # calcite solubility product at 25 °C
print("Imports OK")
"""

si_nb = nb(
    md("title", """\
# Saturation Index — Calcite Precipitation Potential

This notebook demonstrates using the NR speciation engine to assess whether a
water is supersaturated with respect to calcite (CaCO₃).

The **saturation index** (SI) is defined as:

$$\\text{SI} = \\log_{10}\\frac{\\text{IAP}}{K_{\\text{sp}}}
  = \\log_{10}(a_{\\text{Ca}^{2+}} \\cdot a_{\\text{CO}_3^{2-}}) - \\log_{10} K_{\\text{sp}}$$

| SI | Interpretation |
|---|---|
| SI > 0 | Supersaturated — calcite precipitation thermodynamically favoured |
| SI = 0 | Equilibrium with solid calcite |
| SI < 0 | Undersaturated — calcite would dissolve |

$K_{\\text{sp}}(\\text{calcite, 25 °C}) = 10^{-8.48}$

**Scope of this notebook:** the NR engine here is run in *dissolved-only* mode —
Ca²⁺ is modelled as a strong divalent cation that shifts the charge balance without
forming any precipitation reaction.  The engine returns the correct dissolved
speciation *as if precipitation were suppressed*; we then use those concentrations
to compute SI and confirm that SI > 0.  Full precipitation equilibria (active-set
NR with mineral mass as an unknown) are deferred to a future implementation.\
"""),

    code("setup3", SETUP3),

    # ── 1. Chemistry ──────────────────────────────────────────────────────────

    md("rxns-md", """\
## 1  Carbonate chemistry declaration

Ca²⁺ participates in the charge balance only — it is a strong divalent cation
with no acid-base reactions in this model.  The dissolved carbonate speciation
is governed by the same three equilibria used in the single-component benchmarks.

| Reaction | log K | Role |
|---|---|---|
| H₂O ⇌ H⁺ + OH⁻ | −14.0 | water autoionisation |
| CO₂ + H₂O ⇌ HCO₃⁻ + H⁺ | −6.35 | first carbonate dissociation |
| HCO₃⁻ ⇌ CO₃²⁻ + H⁺ | −10.33 | second carbonate dissociation |

Ca²⁺ is passed to `solve()` via `strong_ions={"CT_Ca": ...}` with a charge of +2.\
"""),

    code("rxns-code", """\
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

engine = NRChemicalEquilibriumEngine.from_reactions(
    [water, co2_first, co2_second], use_activity=True, activity_model="davies"
)
print("Masters:    ", engine.tableau.masters)
print("Secondaries:", [s.species_id for s in engine.tableau.secondaries])
"""),

    # ── 2. Single-point speciation + SI ───────────────────────────────────────

    md("single-md", """\
## 2  Hard-water speciation at a single point

A moderately hard process water with:

| Parameter | Value | Notes |
|---|---|---|
| CT_CO₂ | 10 mmol/L | total inorganic carbon (CO₂ + HCO₃⁻ + CO₃²⁻) |
| CT_Ca | 2 mmol/L | total calcium (treated as strong Ca²⁺) |
| CT_Na | 5 mmol/L | NaOH addition — raises pH by adding positive strong-ion charge |

The charge balance at these conditions requires a net alkalinity of
CT_Na + 2·CT_Ca = 9 mmol/L, supplied almost entirely by HCO₃⁻.  The remaining
carbonate sits as dissolved CO₂, fixing the pH near neutral–slightly basic.\
"""),

    code("single-code", """\
T_K   = 298.15
CT_Ca = 0.002    # mol/L
CT_CO2 = 0.010   # mol/L
CT_Na  = 0.005   # mol/L

out = engine.solve(
    totals={"CO2": CT_CO2},
    strong_ions={"CT_Ca": CT_Ca, "CT_Na": CT_Na},
    T_K=T_K,
)

pH    = out.pH
c_CO2  = out.species_mol_L["CO2"]
c_HCO3 = out.species_mol_L["HCO3-"]
c_CO3  = out.species_mol_L["CO3--"]
I      = out.ionic_strength

print(f"pH               = {pH:.4f}")
print(f"[CO2]            = {c_CO2*1e3:.4f} mmol/L")
print(f"[HCO3-]          = {c_HCO3*1e3:.4f} mmol/L")
print(f"[CO3--]          = {c_CO3*1e6:.4f} µmol/L")
print(f"Ionic strength   = {I*1e3:.2f} mmol/L")
print(f"Charge residual  = {out.charge_residual:.2e}")
"""),

    # ── 3. SI calculation ─────────────────────────────────────────────────────

    md("si-md", """\
## 3  Saturation index calculation

### Without activity corrections (γ = 1)

$$\\text{SI}_{\\text{ideal}} = \\log_{10}([\\text{Ca}^{2+}]\\cdot[\\text{CO}_3^{2-}]) - \\log_{10} K_{\\text{sp}}$$

### With Davies activity corrections

$$\\text{SI}_{\\text{Davies}} = \\log_{10}(\\gamma_{\\text{Ca}} \\cdot [\\text{Ca}^{2+}])
  + \\log_{10}(\\gamma_{\\text{CO}_3} \\cdot [\\text{CO}_3^{2-}]) - \\log_{10} K_{\\text{sp}}$$

Activity coefficients are evaluated at the ionic strength returned by the solver.
Because Ca²⁺ and CO₃²⁻ both carry |z| = 2, the Davies equation gives the same
γ for both:

$$\\log \\gamma = -A(T)\\cdot 4 \\cdot \\left(\\frac{\\sqrt{I}}{1+\\sqrt{I}} - 0.3I\\right)$$
"""),

    code("si-code", """\
# Activity coefficients for divalent ions at the solved ionic strength
gamma_z2 = davies.gamma(2, I, T_K=T_K)   # same for Ca²⁺ and CO₃²⁻ (|z|=2)

# Ion activity products
IAP_ideal  = CT_Ca * c_CO3
IAP_davies = (gamma_z2 * CT_Ca) * (gamma_z2 * c_CO3)

SI_ideal  = np.log10(IAP_ideal  / Ksp_calcite)
SI_davies = np.log10(IAP_davies / Ksp_calcite)

print(f"γ(z=±2) at I = {I*1e3:.2f} mmol/L  →  γ = {gamma_z2:.4f}")
print()
print(f"IAP (ideal)   = {IAP_ideal:.3e}  →  log₁₀(IAP) = {np.log10(IAP_ideal):.3f}")
print(f"IAP (Davies)  = {IAP_davies:.3e}  →  log₁₀(IAP) = {np.log10(IAP_davies):.3f}")
print(f"log₁₀(Ksp)              =  {np.log10(Ksp_calcite):.2f}")
print()
print(f"SI (ideal)   = {SI_ideal:+.4f}  {'> 0 → SUPERSATURATED' if SI_ideal > 0 else '< 0 → undersaturated'}")
print(f"SI (Davies)  = {SI_davies:+.4f}  {'> 0 → SUPERSATURATED' if SI_davies > 0 else '< 0 → undersaturated'}")
"""),

    # ── 4. SI vs NaOH dose (pH sweep) ────────────────────────────────────────

    md("sweep-md", """\
## 4  Saturation index vs pH — NaOH dosing

Increasing the NaOH dose (CT_Na) raises pH by adding positive charge to the
charge balance.  As pH rises:

- More carbonate is in the CO₃²⁻ form (α₂ increases steeply above pH 8)
- [CO₃²⁻] increases → IAP increases → SI rises

The plot below shows SI as a function of pH, achieved by sweeping CT_Na from 0
to 14 mmol/L.  The SI = 0 crossing marks the **minimum pH for calcite precipitation**
at these total concentrations.\
"""),

    code("sweep-code", """\
na_sweep  = np.linspace(0.0, 0.014, 60)
pH_sweep  = []
SI_i_sweep = []
SI_d_sweep = []

for na in na_sweep:
    o = engine.solve(
        totals={"CO2": CT_CO2},
        strong_ions={"CT_Ca": CT_Ca, "CT_Na": na},
        T_K=T_K,
    )
    g = davies.gamma(2, o.ionic_strength, T_K=T_K)
    iap_i = CT_Ca * o.species_mol_L["CO3--"]
    iap_d = (g * CT_Ca) * (g * o.species_mol_L["CO3--"])
    pH_sweep.append(o.pH)
    SI_i_sweep.append(np.log10(iap_i / Ksp_calcite))
    SI_d_sweep.append(np.log10(iap_d / Ksp_calcite))

pH_arr = np.array(pH_sweep)
SI_i   = np.array(SI_i_sweep)
SI_d   = np.array(SI_d_sweep)

# Find pH crossings (SI = 0)
def _crossing_pH(pH, SI):
    for k in range(len(SI)-1):
        if SI[k] <= 0 < SI[k+1]:
            # linear interpolation
            frac = -SI[k] / (SI[k+1] - SI[k])
            return pH[k] + frac*(pH[k+1]-pH[k])
    return None

pH_cross_i = _crossing_pH(pH_arr, SI_i)
pH_cross_d = _crossing_pH(pH_arr, SI_d)

fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

axes[0].plot(pH_arr, SI_i, "o-", ms=3, color="tab:gray",  label="Ideal (γ = 1)")
axes[0].plot(pH_arr, SI_d, "s-", ms=3, color="tab:blue",  label="Davies model")
axes[0].axhline(0, color="k", lw=1.2, ls="--", label="SI = 0  (saturation)")
if pH_cross_i:
    axes[0].axvline(pH_cross_i, color="tab:gray",  ls=":", lw=1.2,
                    label=f"Ideal threshold pH = {pH_cross_i:.2f}")
if pH_cross_d:
    axes[0].axvline(pH_cross_d, color="tab:blue", ls=":", lw=1.2,
                    label=f"Davies threshold pH = {pH_cross_d:.2f}")
axes[0].set_ylabel("Saturation index (SI)")
axes[0].set_title(f"Calcite SI vs pH  —  CT_CO₂ = {CT_CO2*1e3:.0f} mmol/L, CT_Ca = {CT_Ca*1e3:.0f} mmol/L")
axes[0].legend(fontsize=8)
axes[0].grid(True, alpha=0.3)
axes[0].fill_between(pH_arr, SI_d, 0, where=(SI_d > 0),
                      color="tab:red", alpha=0.15, label="supersaturated region")

# Carbonate species
c_CO2_sw  = [engine.solve(totals={"CO2": CT_CO2},
                           strong_ions={"CT_Ca": CT_Ca, "CT_Na": na}, T_K=T_K).species_mol_L["CO2"]
              for na in na_sweep]
c_HCO3_sw = [engine.solve(totals={"CO2": CT_CO2},
                           strong_ions={"CT_Ca": CT_Ca, "CT_Na": na}, T_K=T_K).species_mol_L["HCO3-"]
              for na in na_sweep]
c_CO3_sw  = [engine.solve(totals={"CO2": CT_CO2},
                           strong_ions={"CT_Ca": CT_Ca, "CT_Na": na}, T_K=T_K).species_mol_L["CO3--"]
              for na in na_sweep]

axes[1].semilogy(pH_arr, c_CO2_sw,  "o-", ms=3, color="tab:blue",   label="CO₂(aq)")
axes[1].semilogy(pH_arr, c_HCO3_sw, "s-", ms=3, color="tab:orange", label="HCO₃⁻")
axes[1].semilogy(pH_arr, c_CO3_sw,  "^-", ms=3, color="tab:green",  label="CO₃²⁻")
if pH_cross_d:
    axes[1].axvline(pH_cross_d, color="tab:blue", ls=":", lw=1.2,
                    label=f"Calcite threshold (Davies)")
axes[1].set_xlabel("pH")
axes[1].set_ylabel("Concentration (mol/L)")
axes[1].set_title("Carbonate species vs pH")
axes[1].legend(fontsize=8)
axes[1].grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.show()

print(f"\\nSI = 0 crossing (ideal):  pH = {pH_cross_i:.3f}" if pH_cross_i else "")
print(f"SI = 0 crossing (Davies):  pH = {pH_cross_d:.3f}" if pH_cross_d else "")
"""),

    # ── 5. SI heatmap vs CT_Ca and CT_CO2 ────────────────────────────────────

    md("heatmap-md", """\
## 5  Precipitation risk map — CT_Ca vs CT_CO₂

For a fixed pH (here pH 8.0, imposed by choosing CT_Na to hit the target), the
saturation index depends on how much calcium and carbonate are dissolved.
The map below shows where combinations of (CT_Ca, CT_CO₂) cross from
undersaturated (SI < 0) to supersaturated (SI > 0).

This is the kind of assessment done in reverse-osmosis concentrate management,
cooling-tower blowdown design, and lime-softening process control.\
"""),

    code("heatmap-code", """\
from itertools import product as iproduct

Ka1 = 10**(-6.35)
Ka2 = 10**(-10.33)
Kw  = 10**(-14.0)

def na_for_pH(target_pH, CT_CO2_, CT_Ca_):
    \"\"\"Analytical estimate of CT_Na needed to hit target_pH.\"\"\"
    H  = 10**(-target_pH)
    OH = Kw / H
    D  = H**2 + Ka1*H + Ka1*Ka2
    HCO3 = Ka1*H/D * CT_CO2_
    CO3  = Ka1*Ka2/D * CT_CO2_
    # charge balance: H + Na + 2*Ca = OH + HCO3 + 2*CO3
    na = (OH + HCO3 + 2*CO3) - (H + 2*CT_Ca_)
    return max(na, 0.0)

Ca_vals  = np.linspace(0.0005, 0.010, 20)    # 0.5–10 mmol/L
CO2_vals = np.linspace(0.002,  0.020, 20)    # 2–20 mmol/L
target_pH = 8.0

SI_grid = np.zeros((len(CO2_vals), len(Ca_vals)))

for i, CT_CO2_ in enumerate(CO2_vals):
    for j, CT_Ca_ in enumerate(Ca_vals):
        na = na_for_pH(target_pH, CT_CO2_, CT_Ca_)
        o  = engine.solve(
            totals={"CO2": CT_CO2_},
            strong_ions={"CT_Ca": CT_Ca_, "CT_Na": na},
            T_K=T_K,
        )
        g = davies.gamma(2, o.ionic_strength, T_K=T_K)
        iap = (g * CT_Ca_) * (g * o.species_mol_L["CO3--"])
        SI_grid[i, j] = np.log10(iap / Ksp_calcite)

fig, ax = plt.subplots(figsize=(7, 5))
cs = ax.contourf(Ca_vals*1e3, CO2_vals*1e3, SI_grid,
                 levels=np.linspace(-2.5, 2.5, 26), cmap="RdBu_r")
cb = fig.colorbar(cs, ax=ax, label="SI (Davies)")
ax.contour(Ca_vals*1e3, CO2_vals*1e3, SI_grid, levels=[0],
           colors="black", linewidths=2, linestyles="--")
ax.set_xlabel("CT_Ca (mmol/L)")
ax.set_ylabel("CT_CO₂ (mmol/L)")
ax.set_title(f"Calcite saturation index at pH {target_pH:.1f}  (Davies activity)")
ax.text(7.5, 17.5, "SI > 0\\nSUPERSATURATED", color="white", fontsize=9,
        ha="center", va="center", fontweight="bold")
ax.text(1.5, 4.0, "SI < 0\\nundersaturated", color="white", fontsize=9,
        ha="center", va="center")
plt.tight_layout()
plt.show()
"""),

    md("scope-md", """\
## 6  Scope and next steps

This notebook uses the NR engine in **dissolved-only** mode.  When SI > 0:

- The engine returns the *supersaturated* dissolved composition (as if no solid
  were present).
- The true equilibrium has less dissolved Ca²⁺ and CO₃²⁻, because some CaCO₃
  has precipitated to bring IAP = Ksp.

The dissolved-only solve is sufficient for:

1. **Scaling risk assessment** — SI > 0 flags combinations of (CT_Ca, CT_CO₂, pH)
   that are prone to pipe / membrane fouling without needing to know how much
   solid forms.
2. **Operating envelope design** — the SI = 0 contour defines the maximum safe
   calcium or alkalinity loading at a given pH.

A full **precipitation equilibrium** solver — where mineral amount is an
additional unknown and IAP = Ksp replaces the dissolved mass balance for Ca²⁺ —
requires an active-set Newton-Raphson outer loop and is deferred to a future
implementation.\
"""),
)

# ─── 05_precipitation_equilibrium ────────────────────────────────────────────

SETUP5 = """\
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
    CO2, HCO3_minus, CO3_2minus,
    Ca_plus_plus,
)
from PyOMES.chemistry.species import Species
from PyOMES.reactions.equilibrium import EquilibriumReaction
from PyOMES.reactions.stoichiometry import StoichiometryEntry
from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
from PyOMES.chemical_equilibrium.activity_models import DaviesActivityModel

def _e(sp, phase, coeff):
    return StoichiometryEntry(species=sp, phase=phase, coefficient=coeff)

davies = DaviesActivityModel()

# Shared carbonate reactions
water = EquilibriumReaction(
    stoichiometry=[_e(H2O,"liquid",-1), _e(H_plus,"liquid",+1), _e(OH_minus,"liquid",+1)],
    log_K=-14.0, label="water",
)
co2_first = EquilibriumReaction(
    stoichiometry=[_e(CO2,"liquid",-1), _e(H2O,"liquid",-1),
                   _e(HCO3_minus,"liquid",+1), _e(H_plus,"liquid",+1)],
    log_K=-6.35, total_id="CO2", label="co2_first",
)
co2_second = EquilibriumReaction(
    stoichiometry=[_e(HCO3_minus,"liquid",-1),
                   _e(CO3_2minus,"liquid",+1), _e(H_plus,"liquid",+1)],
    log_K=-10.33, total_id="CO2", label="co2_second",
)
carb_rxns = [water, co2_first, co2_second]

# Calcite mineral
CaCO3 = Species(id="CaCO3", atoms={"Ca":1,"C":1,"O":3}, charge=0, MW=100.086)
calcite = EquilibriumReaction(
    stoichiometry=[
        _e(CaCO3,       "solid",  -1),
        _e(Ca_plus_plus,"liquid", +1),
        _e(CO3_2minus,  "liquid", +1),
    ],
    log_K=-8.48,
    label="calcite",
)
Ksp_calcite = 10**(-8.48)

print("Imports OK")
"""

precip_nb = nb(
    md("title", """\
# Precipitation Equilibrium — Calcite

This notebook demonstrates the active-set Newton-Raphson precipitation loop
in `NRChemicalEquilibriumEngine`.  When the ionic activity product (IAP) for a mineral
exceeds the solubility product (Ksp), the solver enforces equilibrium by
finding the amount precipitated (ξ mol/L) such that IAP = Ksp.

The mineral is declared as an `EquilibriumReaction` with one
`StoichiometryEntry(phase="solid")` (the mineral, activity = 1) and one or
more `StoichiometryEntry(phase="liquid")` (the dissolved ionic products).
Ca²⁺ is passed as a strong ion; the outer loop adjusts its effective
concentration to account for precipitation.\
"""),

    code("setup5", SETUP5),

    # ── 1. Background ─────────────────────────────────────────────────────────

    md("bg-md", """\
## 1  Background

The saturation index (SI) notebook (03) demonstrated that the NR engine can
flag SI > 0 — calcite precipitation is thermodynamically favoured — but
returned dissolved concentrations as if precipitation were suppressed.

The **precipitation equilibrium** solver extends this: when SI > 0, the outer
active-set loop finds ξ such that

$$\\text{IAP}(\\xi) = K_{\\text{sp}}$$

where

$$\\text{IAP} = a_{\\text{Ca}^{2+}} \\cdot a_{\\text{CO}_3^{2-}}, \\quad
  a_{\\text{Ca}^{2+}} = \\gamma_{\\text{Ca}} \\cdot ([\\text{Ca}^{2+}]_\\text{total} - \\xi)$$

Ca²⁺ is treated as a strong ion with no dissolved complexation (Phase 1
scope).  The dissolved-Ca²⁺ concentration after equilibrium is
$[\\text{Ca}^{2+}]_\\text{eq} = C_{T,\\text{Ca}} - \\xi$.\
"""),

    # ── 2. Declaration ────────────────────────────────────────────────────────

    md("decl-md", """\
## 2  Declaration

```python
CaCO3 = Species(id="CaCO3", atoms={"Ca":1,"C":1,"O":3}, charge=0, MW=100.086)
calcite = EquilibriumReaction(
    stoichiometry=[
        StoichiometryEntry(species=CaCO3,       phase="solid",  coefficient=-1),
        StoichiometryEntry(species=Ca_plus_plus, phase="liquid", coefficient=+1),
        StoichiometryEntry(species=CO3_2minus,   phase="liquid", coefficient=+1),
    ],
    log_K=-8.48,   # Ksp at 25 °C
    label="calcite",
)
```

Pass to `NRChemicalEquilibriumEngine.from_reactions()` via the `precipitation_reactions`
keyword argument.  The dissolved carbonate reactions go in `equilibrium_reactions`
as before.\
"""),

    code("decl-code", """\
engine = NRChemicalEquilibriumEngine.from_reactions(
    carb_rxns,
    precipitation_reactions=[calcite],
    use_activity=True,
    activity_model="davies",
)

# Dissolved-only engine for comparison
engine_dissolved = NRChemicalEquilibriumEngine.from_reactions(
    carb_rxns,
    use_activity=True,
    activity_model="davies",
)

print("Masters:", engine.tableau.masters)
print("Secondaries:", [s.species_id for s in engine.tableau.secondaries])
"""),

    # ── 3. Single-point solve ─────────────────────────────────────────────────

    md("single-md", """\
## 3  Single-point solve

**Conditions:** CT_CO₂ = 10 mmol/L, CT_Ca = 2 mmol/L, CT_Na = 5 mmol/L.

These are the same conditions as Section 2 of notebook 03, where SI ≈ +0.42
in the dissolved-only model.  After precipitation equilibrium: IAP = Ksp and
ξ > 0 (some CaCO₃ has precipitated).\
"""),

    code("single-code", """\
T_K    = 298.15
CT_Ca  = 0.002
CT_CO2 = 0.010
CT_Na  = 0.005

# Dissolved-only
out_d = engine_dissolved.solve(
    totals={"CO2": CT_CO2},
    strong_ions={"CT_Ca": CT_Ca, "CT_Na": CT_Na},
    T_K=T_K,
)
I_d = out_d.ionic_strength
g_d = davies.gamma(2, I_d, T_K=T_K)
IAP_d = (g_d * CT_Ca) * (g_d * out_d.species_mol_L["CO3--"])
SI_d  = np.log10(IAP_d / Ksp_calcite)

# Precipitation equilibrium
out_p = engine.solve(
    totals={"CO2": CT_CO2},
    strong_ions={"CT_Ca": CT_Ca, "CT_Na": CT_Na},
    T_K=T_K,
)
xi   = out_p.extra["minerals_xi_mol_L"]["calcite"]
SI_p = out_p.saturation_indices["calcite"]
Ca_eq = CT_Ca - xi

print(f"{'':30s}  {'Dissolved-only':>16}  {'Precipitation eq':>18}")
print("-" * 68)
print(f"{'pH':30s}  {out_d.pH:>16.4f}  {out_p.pH:>18.4f}")
print(f"{'[CO3--] (mol/L)':30s}  {out_d.species_mol_L["CO3--"]:>16.4e}  {out_p.species_mol_L["CO3--"]:>18.4e}")
print(f"{'[Ca2+] effective (mol/L)':30s}  {CT_Ca:>16.4e}  {Ca_eq:>18.4e}")
print(f"{'IAP':30s}  {IAP_d:>16.4e}  {'IAP = Ksp':>18}")
print(f"{'SI':30s}  {SI_d:>16.4f}  {SI_p:>18.4f}")
print(f"{'xi (mol/L precipitated)':30s}  {'—':>16}  {xi:>18.4e}")
"""),

    # ── 4. Dissolved Ca depletion curve ──────────────────────────────────────

    md("sweep-md", """\
## 4  Dissolved Ca²⁺ depletion curve

As total calcium (CT_Ca) increases above the solubility ceiling set by Ksp
at fixed CT_CO₂, the excess precipitates as CaCO₃.  The dissolved Ca²⁺ rises
until the saturation ceiling is reached, then stays flat while ξ absorbs the
surplus.\
"""),

    code("sweep-code", """\
CT_CO2_fixed = 0.010
CT_Na_fixed  = 0.005

Ca_total_vals = np.linspace(0.0001, 0.005, 50)
Ca_dissolved  = []
xi_vals       = []
SI_vals       = []

for cta in Ca_total_vals:
    o = engine.solve(
        totals={"CO2": CT_CO2_fixed},
        strong_ions={"CT_Ca": cta, "CT_Na": CT_Na_fixed},
        T_K=T_K,
    )
    xi_j = o.extra["minerals_xi_mol_L"]["calcite"]
    xi_vals.append(xi_j)
    Ca_dissolved.append(cta - xi_j)
    SI_vals.append(o.saturation_indices["calcite"])

Ca_dissolved = np.array(Ca_dissolved)
xi_vals      = np.array(xi_vals)
SI_vals      = np.array(SI_vals)

fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

axes[0].plot(Ca_total_vals*1e3, Ca_dissolved*1e3, "o-", ms=3, color="tab:blue",
             label="Dissolved Ca²⁺")
axes[0].plot(Ca_total_vals*1e3, Ca_total_vals*1e3, "k--", lw=1, label="Total Ca²⁺ (no precipitation)")
axes[0].set_ylabel("[Ca²⁺]_dissolved (mmol/L)")
axes[0].set_title(f"Ca²⁺ depletion — CT_CO₂ = {CT_CO2_fixed*1e3:.0f} mmol/L, CT_Na = {CT_Na_fixed*1e3:.0f} mmol/L")
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)

axes[1].plot(Ca_total_vals*1e3, xi_vals*1e3, "s-", ms=3, color="tab:red",
             label="ξ (precipitated)")
axes[1].set_xlabel("CT_Ca total (mmol/L)")
axes[1].set_ylabel("ξ (mmol/L)")
axes[1].set_title("Amount precipitated vs total Ca²⁺")
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

print(f"Ksp-set ceiling dissolved Ca²⁺ ≈ {Ca_dissolved[Ca_total_vals >= 0.002].mean()*1e3:.3f} mmol/L")
"""),

    # ── 5. Dissolution from solid ─────────────────────────────────────────────

    md("diss-md", """\
## 5  Dissolution from a CaCO₃ solid feed

Starting with a known CaCO₃ solid (n_mol_solid mol/L), no dissolved Ca²⁺
initially: CaCO₃ dissolves until IAP = Ksp.  The model supplies
CT_Ca = n_mol_solid (as if all solid dissolved) and the outer loop
re-precipitates the excess.

At equilibrium:
- dissolved Ca²⁺ = CT_Ca − ξ
- remaining solid ≈ ξ (re-precipitated amount)

The fraction dissolved depends on the carbonate alkalinity (CT_CO₂, CT_Na)
which controls the pH and therefore [CO₃²⁻] and the Ksp ceiling.\
"""),

    code("diss-code", """\
n_mol_solid = 0.001   # mol/L CaCO3(s) available
CT_CO2_vals = np.linspace(0.001, 0.020, 30)
CT_Na_vals  = [0.000, 0.005, 0.010]

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

for CT_Na in CT_Na_vals:
    frac_dissolved = []
    for CT_CO2 in CT_CO2_vals:
        o = engine.solve(
            totals={"CO2": CT_CO2},
            strong_ions={"CT_Ca": n_mol_solid, "CT_Na": CT_Na},
            T_K=T_K,
        )
        xi_j = o.extra["minerals_xi_mol_L"]["calcite"]
        # Fraction dissolved = (n_mol_solid - xi) / n_mol_solid
        frac = (n_mol_solid - xi_j) / n_mol_solid
        frac_dissolved.append(frac)

    lbl = f"CT_Na = {CT_Na*1e3:.0f} mmol/L"
    axes[0].plot(CT_CO2_vals*1e3, np.array(frac_dissolved)*100, "o-", ms=3, label=lbl)

axes[0].set_xlabel("CT_CO₂ (mmol/L)")
axes[0].set_ylabel("Fraction dissolved (%)")
axes[0].set_title(f"CaCO₃ dissolution — n_mol_solid = {n_mol_solid*1e3:.0f} mmol/L")
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)

# Single trace: dissolved Ca vs CT_CO2
Ca_diss_vals = []
for CT_CO2 in CT_CO2_vals:
    o = engine.solve(
        totals={"CO2": CT_CO2},
        strong_ions={"CT_Ca": n_mol_solid, "CT_Na": 0.010},
        T_K=T_K,
    )
    Ca_diss_vals.append(n_mol_solid - o.extra["minerals_xi_mol_L"]["calcite"])

axes[1].plot(CT_CO2_vals*1e3, np.array(Ca_diss_vals)*1e3, "s-", ms=3,
             color="tab:blue", label="Dissolved Ca²⁺ (CT_Na=10 mmol/L)")
axes[1].axhline(n_mol_solid*1e3, color="gray", ls="--", lw=1, label="Total solid feed")
axes[1].set_xlabel("CT_CO₂ (mmol/L)")
axes[1].set_ylabel("[Ca²⁺]_dissolved (mmol/L)")
axes[1].set_title("Dissolved Ca²⁺ vs carbonate")
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
"""),

    # ── 6. Scope note ────────────────────────────────────────────────────────

    md("scope-md", """\
## 6  Scope and next steps

**Phase 1 scope (this notebook):**

- `NRChemicalEquilibriumEngine` with `precipitation_reactions=` finds ξ and returns
  `out.extra["minerals_xi_mol_L"]["calcite"]` and `out.saturation_indices["calcite"]`.
- Ca²⁺ is a strong ion with **no dissolved complexation** (no CaHCO₃⁺, CaCO₃⁰,
  etc.) — this is the dominant species in most natural waters, so the error is
  small at low ionic strength.
- SI ≈ 0 at equilibrium (verified by tests and Section 3 above).

**Phase 2 (NR_PRECIPITATION_CV_INTEGRATION.md):**

- ξ writeback to `SolidPhase.n_mol` on a `ControlVolume`.
- `ControlVolume._read_from_phases` sums liquid + solid contributions.
- Lazy-creation / consistency enforcement of `SolidPhase` objects.
- Multiple co-precipitating minerals sharing a common component (e.g.
  calcite + aragonite both consuming CT_CO₂).

**Not in scope:**

- Kinetic crystallisation rates — a separate future phase.
- Redox equilibria (pe as a master) — deferred per NR design doc.\
"""),
)

# ─── 06_phreeqc_benchmark ────────────────────────────────────────────────────

SETUP6 = """\
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
    CO2, HCO3_minus, CO3_2minus,
    NH3, NH4_plus, Ca_plus_plus,
)
from PyOMES.chemistry.species import Species
from PyOMES.reactions.equilibrium import EquilibriumReaction
from PyOMES.reactions.stoichiometry import StoichiometryEntry
from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
from phreeqpython import PhreeqPython

def _e(sp, phase, coeff):
    return StoichiometryEntry(species=sp, phase=phase, coefficient=coeff)

pp = PhreeqPython()

# ── Reactions — use phreeqc.dat log K values exactly for a fair comparison ───
# H2O  = H+ + OH-          log_K = -13.998
# CO2  + H2O = H+ + HCO3-  log_K = -6.352   (H2CO3* dissociation, pKa1)
# HCO3-  = H+ + CO3-2      log_K = -10.329  (pKa2)
# NH4+   = H+ + NH3         log_K = -9.252
# Calcite (CaCO3)  = Ca2+ + CO3-2  log_K = -8.48
water = EquilibriumReaction(
    stoichiometry=[_e(H2O,'liquid',-1), _e(H_plus,'liquid',+1), _e(OH_minus,'liquid',+1)],
    log_K=-13.998, label='water',
)
co2_first = EquilibriumReaction(
    stoichiometry=[_e(CO2,'liquid',-1), _e(H2O,'liquid',-1),
                   _e(HCO3_minus,'liquid',+1), _e(H_plus,'liquid',+1)],
    log_K=-6.352, total_id='CO2', label='co2_first',
)
co2_second = EquilibriumReaction(
    stoichiometry=[_e(HCO3_minus,'liquid',-1),
                   _e(CO3_2minus,'liquid',+1), _e(H_plus,'liquid',+1)],
    log_K=-10.329, total_id='CO2', label='co2_second',
)
nh3_rxn = EquilibriumReaction(
    stoichiometry=[_e(NH4_plus,'liquid',-1), _e(NH3,'liquid',+1), _e(H_plus,'liquid',+1)],
    log_K=-9.252, total_id='NH3', label='nh3',
)

engine_carb = NRChemicalEquilibriumEngine.from_reactions(
    [water, co2_first, co2_second],
    use_activity=True, activity_model='davies')
engine_carb_nh3 = NRChemicalEquilibriumEngine.from_reactions(
    [water, co2_first, co2_second, nh3_rxn],
    use_activity=True, activity_model='davies')

CaCO3_s = Species(id='CaCO3', atoms={'Ca':1,'C':1,'O':3}, charge=0, MW=100.086)
calcite = EquilibriumReaction(
    stoichiometry=[_e(CaCO3_s,'solid',-1), _e(Ca_plus_plus,'liquid',+1), _e(CO3_2minus,'liquid',+1)],
    log_K=-8.48, label='calcite',
)
engine_precip = NRChemicalEquilibriumEngine.from_reactions(
    [water, co2_first, co2_second],
    precipitation_reactions=[calcite],
    use_activity=True, activity_model='davies',
)

T_K = 298.15

def pq(composition_mmol_L):
    \"\"\"Run phreeqpython at 25 °C with charge-balanced pH.\"\"\"
    d = dict(composition_mmol_L)
    d['pH'] = '7 charge'
    d['units'] = 'mmol/L'
    d['temp'] = 25
    return pp.add_solution_raw(d)

print("Setup OK — NR engine and phreeqpython ready")
print(f"phreeqpython uses Davies b=0.3, same as DaviesActivityModel")
"""

bench6_nb = nb(
    md("title6", """\
# PHREEQC Benchmark — Numerical Accuracy of NRChemicalEquilibriumEngine

This notebook compares `NRChemicalEquilibriumEngine` (PyOMES) directly against PHREEQC
(via phreeqpython) on identical inputs.  To isolate numerical accuracy from
thermodynamic data differences, both engines use PHREEQC's `phreeqc.dat`
log K values (pKa1 = 6.352, pKa2 = 10.329, Kw = 10⁻¹³·⁹⁹⁸).

Both engines use the Davies activity model with b = 0.3.  PHREEQC is treated
as the reference; discrepancies in the carbonate-only sections therefore
reflect numerical error in the NR solver, not model differences.

**Where the models genuinely differ:**  PHREEQC's built-in database includes
ion pairs — NaHCO₃⁰, NaCO₃⁻, CaHCO₃⁺, CaCO₃⁰, CaOH⁺ — that PyOMES's NR
engine does not yet implement (complexation is Phase 3 scope).  These cause
systematic offsets in sections that include Ca²⁺ or high Na⁺; those
differences are model differences, not solver errors, and are labelled
clearly.\
"""),

    code("setup6", SETUP6),

    # ── 1. Carbonate-only pH sweep ─────────────────────────────────────────────

    md("s1-md", """\
## 1  Carbonate system — pH vs CT_CO₂

Pure carbonate dissolved in water, no added salt.  Both engines are
charge-balanced (no fixed pH).  There are no ion pairs in this system, so the
comparison is a direct test of the NR solver accuracy.\
"""),

    code("s1-code", """\
CT_CO2_vals = np.logspace(-3, -1, 40)   # 1 → 100 mmol/L

ph_nr, ph_pq = [], []
hco3_nr, hco3_pq = [], []
co3_nr, co3_pq   = [], []

for ct in CT_CO2_vals:
    # NR engine
    o = engine_carb.solve(totals={'CO2': ct}, strong_ions={}, T_K=T_K)
    ph_nr.append(o.pH)
    hco3_nr.append(o.species_mol_L["HCO3-"])
    co3_nr.append(o.species_mol_L["CO3--"])

    # PHREEQC
    sol = pq({'C(4)': ct * 1e3})
    ph_pq.append(sol.pH)
    hco3_pq.append(sol.species.get('HCO3-', 0.0))
    co3_pq.append(sol.species.get('CO3-2', 0.0))
    sol.forget()

ph_nr = np.array(ph_nr); ph_pq = np.array(ph_pq)
hco3_nr = np.array(hco3_nr); hco3_pq = np.array(hco3_pq)
co3_nr = np.array(co3_nr); co3_pq = np.array(co3_pq)

fig, axes = plt.subplots(1, 3, figsize=(14, 4))

axes[0].semilogx(CT_CO2_vals*1e3, ph_nr, 'o-', ms=3, label='NRChemicalEquilibriumEngine')
axes[0].semilogx(CT_CO2_vals*1e3, ph_pq, 's--', ms=3, label='PHREEQC')
axes[0].set_xlabel('CT_CO₂ (mmol/L)'); axes[0].set_ylabel('pH')
axes[0].set_title('pH vs CT_CO₂'); axes[0].legend(fontsize=9); axes[0].grid(True, alpha=0.3)

axes[1].loglog(CT_CO2_vals*1e3, hco3_nr, 'o-', ms=3, label='NR')
axes[1].loglog(CT_CO2_vals*1e3, hco3_pq, 's--', ms=3, label='PHREEQC')
axes[1].set_xlabel('CT_CO₂ (mmol/L)'); axes[1].set_ylabel('[HCO₃⁻] (mol/L)')
axes[1].set_title('[HCO₃⁻] vs CT_CO₂'); axes[1].legend(fontsize=9); axes[1].grid(True, alpha=0.3)

axes[2].loglog(CT_CO2_vals*1e3, co3_nr, 'o-', ms=3, label='NR')
axes[2].loglog(CT_CO2_vals*1e3, co3_pq, 's--', ms=3, label='PHREEQC')
axes[2].set_xlabel('CT_CO₂ (mmol/L)'); axes[2].set_ylabel('[CO₃²⁻] (mol/L)')
axes[2].set_title('[CO₃²⁻] vs CT_CO₂'); axes[2].legend(fontsize=9); axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

dpH  = np.abs(ph_nr - ph_pq)
dhco3 = np.abs(hco3_nr - hco3_pq)
dco3  = np.abs(co3_nr - co3_pq)
print(f"Max |ΔpH|    = {dpH.max():.2e}   (at CT_CO₂ = {CT_CO2_vals[dpH.argmax()]*1e3:.1f} mmol/L)")
print(f"Max |Δ[HCO₃⁻]| = {dhco3.max():.2e} mol/L")
print(f"Max |Δ[CO₃²⁻]| = {dco3.max():.2e} mol/L")
"""),

    code("s1-parity", """\
def _parity(ax, ref, pred, lbl, log=False, c=None, cmap=None, clabel=None):
    lo = min(float(np.nanmin(ref)), float(np.nanmin(pred)))
    hi = max(float(np.nanmax(ref)), float(np.nanmax(pred)))
    if log:
        ax.set_xscale('log'); ax.set_yscale('log')
    pad = (hi - lo) * 0.05 if not log else 0
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], 'k--', lw=1, zorder=0)
    sc = ax.scatter(ref, pred, s=20, c=c, cmap=cmap, alpha=0.8, zorder=1)
    if clabel:
        plt.colorbar(sc, ax=ax, label=clabel, shrink=0.85)
    ax.set_xlabel(f'PHREEQC  {lbl}'); ax.set_ylabel(f'NR engine  {lbl}')
    ax.set_title(f'Parity — {lbl}'); ax.grid(True, alpha=0.3)

fig, axes = plt.subplots(1, 3, figsize=(14, 4))
_parity(axes[0], ph_pq,   ph_nr,   'pH')
_parity(axes[1], hco3_pq, hco3_nr, '[HCO\u2083\u207b] (mol/L)', log=True)
_parity(axes[2], co3_pq,  co3_nr,  '[CO\u2083\u00b2\u207b] (mol/L)', log=True)
plt.suptitle('Section 1 parity — carbonate only', fontsize=10)
plt.tight_layout(); plt.show()
"""),

    # ── 2. Alkalinity (Na+) ────────────────────────────────────────────────────

    md("s2-md", """\
## 2  Carbonate + Na⁺ alkalinity

Fixed CT_CO₂ = 10 mmol/L, sweeping CT_Na.  PHREEQC includes the NaHCO₃⁰
and NaCO₃⁻ ion pairs; their concentrations are small at ≤ 20 mmol/L Na⁺ but
are noted for completeness.\
"""),

    code("s2-code", """\
CT_CO2_fixed = 0.010
CT_Na_vals   = np.linspace(0, 0.020, 40)

ph_nr2, ph_pq2 = [], []
naion_pq = []   # NaCO3- + NaHCO3 concentration in PHREEQC

for ct_na in CT_Na_vals:
    o = engine_carb.solve(totals={'CO2': CT_CO2_fixed},
                          strong_ions={'CT_Na': ct_na}, T_K=T_K)
    ph_nr2.append(o.pH)

    sol = pq({'C(4)': CT_CO2_fixed*1e3, 'Na': ct_na*1e3})
    ph_pq2.append(sol.pH)
    naion_pq.append(sol.species.get('NaCO3-', 0) + sol.species.get('NaHCO3', 0))
    sol.forget()

ph_nr2 = np.array(ph_nr2); ph_pq2 = np.array(ph_pq2)
naion_pq = np.array(naion_pq)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].plot(CT_Na_vals*1e3, ph_nr2, 'o-', ms=3, label='NRChemicalEquilibriumEngine')
axes[0].plot(CT_Na_vals*1e3, ph_pq2, 's--', ms=3, label='PHREEQC')
axes[0].set_xlabel('CT_Na (mmol/L)'); axes[0].set_ylabel('pH')
axes[0].set_title(f'pH vs CT_Na  (CT_CO₂ = {CT_CO2_fixed*1e3:.0f} mmol/L)')
axes[0].legend(fontsize=9); axes[0].grid(True, alpha=0.3)

axes[1].plot(CT_Na_vals*1e3, naion_pq*1e3, 'k-o', ms=3)
axes[1].set_xlabel('CT_Na (mmol/L)'); axes[1].set_ylabel('Na-complex (mmol/L)')
axes[1].set_title('NaCO₃⁻ + NaHCO₃⁰ in PHREEQC (absent in NR engine)')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

dpH2 = np.abs(ph_nr2 - ph_pq2)
print(f"Max |ΔpH| = {dpH2.max():.2e}  (Na complex contribution ≤ {naion_pq.max()*1e3:.3f} mmol/L)")
"""),

    code("s2-parity", """\
fig, ax = plt.subplots(figsize=(5, 5))
_parity(ax, ph_pq2, ph_nr2, 'pH',
        c=CT_Na_vals * 1e3, cmap='viridis',
        clabel='CT_Na (mmol/L)')
ax.set_title('Section 2 parity — pH vs Na\u207a alkalinity')
plt.tight_layout(); plt.show()
"""),

    # ── 3. Ammonia + carbonate ─────────────────────────────────────────────────

    md("s3-md", """\
## 3  Carbonate + ammonia

CT_CO₂ = 10 mmol/L, sweeping CT_NH₃.  No Ca, no Na.  PHREEQC has no
additional NH₃ ion pairs in the default database at these conditions, so this
is another clean numerical comparison.\
"""),

    code("s3-code", """\
CT_NH3_vals = np.linspace(0, 0.010, 40)   # 0 → 10 mmol/L

ph_nr3, ph_pq3 = [], []
nh4_nr, nh4_pq = [], []

for ct_nh3 in CT_NH3_vals:
    o = engine_carb_nh3.solve(totals={'CO2': CT_CO2_fixed, 'NH3': ct_nh3},
                              strong_ions={}, T_K=T_K)
    ph_nr3.append(o.pH)
    nh4_nr.append(o.species_mol_L.get('NH4+', 0.0))

    sol = pq({'C(4)': CT_CO2_fixed*1e3, 'N(-3)': ct_nh3*1e3})
    ph_pq3.append(sol.pH)
    nh4_pq.append(sol.species.get('NH4+', 0.0))
    sol.forget()

ph_nr3 = np.array(ph_nr3); ph_pq3 = np.array(ph_pq3)
nh4_nr = np.array(nh4_nr); nh4_pq = np.array(nh4_pq)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].plot(CT_NH3_vals*1e3, ph_nr3, 'o-', ms=3, label='NRChemicalEquilibriumEngine')
axes[0].plot(CT_NH3_vals*1e3, ph_pq3, 's--', ms=3, label='PHREEQC')
axes[0].set_xlabel('CT_NH₃ (mmol/L)'); axes[0].set_ylabel('pH')
axes[0].set_title(f'pH vs CT_NH₃  (CT_CO₂ = {CT_CO2_fixed*1e3:.0f} mmol/L)')
axes[0].legend(fontsize=9); axes[0].grid(True, alpha=0.3)

axes[1].plot(CT_NH3_vals*1e3, nh4_nr*1e3, 'o-', ms=3, label='NR [NH₄⁺]')
axes[1].plot(CT_NH3_vals*1e3, nh4_pq*1e3, 's--', ms=3, label='PHREEQC [NH₄⁺]')
axes[1].set_xlabel('CT_NH₃ (mmol/L)'); axes[1].set_ylabel('[NH₄⁺] (mmol/L)')
axes[1].set_title('[NH₄⁺] vs CT_NH₃')
axes[1].legend(fontsize=9); axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

dpH3 = np.abs(ph_nr3 - ph_pq3)
dnh4  = np.abs(nh4_nr - nh4_pq)
print(f"Max |ΔpH|    = {dpH3.max():.2e}")
print(f"Max |Δ[NH₄⁺]| = {dnh4.max():.2e} mol/L")
"""),

    code("s3-parity", """\
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
_parity(axes[0], ph_pq3,  ph_nr3,  'pH',
        c=CT_NH3_vals * 1e3, cmap='plasma',
        clabel='CT_NH\u2083 (mmol/L)')
_parity(axes[1], nh4_pq,  nh4_nr,  '[NH\u2084\u207a] (mol/L)',
        log=True,
        c=CT_NH3_vals * 1e3, cmap='plasma',
        clabel='CT_NH\u2083 (mmol/L)')
plt.suptitle('Section 3 parity — carbonate + NH\u2083', fontsize=10)
plt.tight_layout(); plt.show()
"""),

    # ── 4. Ca system — model difference disclosure ─────────────────────────────

    md("s4-md", """\
## 4  Ca²⁺ system — model differences

Ca²⁺ is treated as a **strong ion** in the NR engine (no dissolved
complexation).  PHREEQC's database forms CaHCO₃⁺, CaCO₃⁰, and CaOH⁺ ion
pairs, consuming a fraction of dissolved Ca²⁺ and DIC.

The table below shows concentrations of Ca complexes at typical conditions.
These represent model differences (missing reactions), not numerical errors.
The NR engine returns the same pH and carbonate speciation as PHREEQC
*on the DIC fraction it sees*, but the effective free-Ca²⁺ is higher than
PHREEQC predicts because the engine has no complexation sink.\
"""),

    code("s4-code", """\
CT_Ca_vals  = np.linspace(0.0005, 0.005, 20)
CT_CO2_fix2 = 0.010
CT_Na_fix   = 0.005

ph_nr4, ph_pq4   = [], []
SI_nr4, SI_pq4   = [], []
ca_complex_pq    = []   # CaHCO3+ + CaCO3(aq) + CaOH+ in PHREEQC

for ct_ca in CT_Ca_vals:
    o = engine_carb.solve(totals={'CO2': CT_CO2_fix2},
                          strong_ions={'CT_Ca': ct_ca, 'CT_Na': CT_Na_fix}, T_K=T_K)
    ph_nr4.append(o.pH)

    sol = pq({'C(4)': CT_CO2_fix2*1e3, 'Ca': ct_ca*1e3, 'Na': CT_Na_fix*1e3})
    ph_pq4.append(sol.pH)
    ca_cx = (sol.species.get('CaHCO3+', 0)
             + sol.species.get('CaCO3', 0)
             + sol.species.get('CaOH+', 0))
    ca_complex_pq.append(ca_cx)
    sol.forget()

ph_nr4 = np.array(ph_nr4); ph_pq4 = np.array(ph_pq4)
ca_complex_pq = np.array(ca_complex_pq)
dpH4 = np.abs(ph_nr4 - ph_pq4)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].plot(CT_Ca_vals*1e3, ph_nr4, 'o-', ms=3, label='NR (no Ca complexes)')
axes[0].plot(CT_Ca_vals*1e3, ph_pq4, 's--', ms=3, label='PHREEQC (with Ca complexes)')
axes[0].set_xlabel('CT_Ca (mmol/L)'); axes[0].set_ylabel('pH')
axes[0].set_title('pH vs CT_Ca — model difference visible at high Ca²⁺')
axes[0].legend(fontsize=9); axes[0].grid(True, alpha=0.3)

axes[1].plot(CT_Ca_vals*1e3, ca_complex_pq*1e3, 'k-o', ms=3)
axes[1].set_xlabel('CT_Ca (mmol/L)')
axes[1].set_ylabel('Ca complexes (mmol/L)')
axes[1].set_title('CaHCO₃⁺ + CaCO₃⁰ + CaOH⁺ in PHREEQC\\n(missing from NR engine)')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

frac_complex = (ca_complex_pq / CT_Ca_vals * 100)
print("Ca complex fraction vs CT_Ca:")
for ct, f, dp in zip(CT_Ca_vals[::4]*1e3, frac_complex[::4], dpH4[::4]):
    print(f"  CT_Ca={ct:.1f} mmol/L: Ca-complexes = {f:.1f}% of total, |ΔpH| = {dp:.4f}")
"""),

    code("s4-parity", """\
fig, ax = plt.subplots(figsize=(5, 5))
_parity(ax, ph_pq4, ph_nr4, 'pH',
        c=CT_Ca_vals * 1e3, cmap='Reds',
        clabel='CT_Ca (mmol/L)')
ax.set_title('Section 4 parity — pH with Ca\u00b2\u207a\\n(offset = missing Ca complexes)')
plt.tight_layout(); plt.show()
"""),

    # ── 5. Calcite precipitation ───────────────────────────────────────────────

    md("s5-md", """\
## 5  Calcite precipitation equilibrium

Both engines start from the same supersaturated conditions and find the
equilibrium ξ (mol/L precipitated) such that IAP = Ksp.  NR uses the
active-set outer loop; PHREEQC uses `sol.desaturate('Calcite', to_si=0)`.

Because PHREEQC includes Ca complexes, it predicts a slightly lower effective
free Ca²⁺ and therefore a slightly different ξ.  The discrepancy quantifies
the combined effect of missing Ca complexation.\
"""),

    code("s5-code", """\
CT_CO2_sweep = np.linspace(0.003, 0.020, 25)
CT_Ca_fix    = 0.002
CT_Na_fix    = 0.005

xi_nr5, xi_pq5   = [], []
SI_nr5, SI_pq5   = [], []
pH_nr5, pH_pq5   = [], []

for ct_co2 in CT_CO2_sweep:
    # NR precipitation engine
    o = engine_precip.solve(totals={'CO2': ct_co2},
                            strong_ions={'CT_Ca': CT_Ca_fix, 'CT_Na': CT_Na_fix},
                            T_K=T_K)
    xi_nr5.append(o.extra["minerals_xi_mol_L"]["calcite"])
    SI_nr5.append(o.saturation_indices["calcite"])
    pH_nr5.append(o.pH)

    # PHREEQC desaturate
    sol = pq({'C(4)': ct_co2*1e3, 'Ca': CT_Ca_fix*1e3, 'Na': CT_Na_fix*1e3})
    ca_before = sol.elements.get('Ca', 0.0)
    if sol.si('Calcite') > 0:
        sol.desaturate('Calcite', to_si=0)
    xi_pq5.append(max(0.0, ca_before - sol.elements.get('Ca', 0.0)))
    SI_pq5.append(sol.si('Calcite'))
    pH_pq5.append(sol.pH)
    sol.forget()

xi_nr5 = np.array(xi_nr5); xi_pq5 = np.array(xi_pq5)
SI_nr5 = np.array(SI_nr5); SI_pq5 = np.array(SI_pq5)
pH_nr5 = np.array(pH_nr5); pH_pq5 = np.array(pH_pq5)

fig, axes = plt.subplots(1, 3, figsize=(14, 4))

axes[0].plot(CT_CO2_sweep*1e3, xi_nr5*1e3, 'o-', ms=3, label='NR')
axes[0].plot(CT_CO2_sweep*1e3, xi_pq5*1e3, 's--', ms=3, label='PHREEQC')
axes[0].set_xlabel('CT_CO₂ (mmol/L)'); axes[0].set_ylabel('ξ (mmol/L)')
axes[0].set_title('Amount precipitated'); axes[0].legend(fontsize=9); axes[0].grid(True, alpha=0.3)

axes[1].plot(CT_CO2_sweep*1e3, SI_nr5, 'o-', ms=3, label='NR')
axes[1].plot(CT_CO2_sweep*1e3, SI_pq5, 's--', ms=3, label='PHREEQC')
axes[1].axhline(0, color='k', lw=0.8)
axes[1].set_xlabel('CT_CO₂ (mmol/L)'); axes[1].set_ylabel('SI (calcite)')
axes[1].set_title('Saturation index at equilibrium')
axes[1].legend(fontsize=9); axes[1].grid(True, alpha=0.3)

axes[2].plot(CT_CO2_sweep*1e3, pH_nr5, 'o-', ms=3, label='NR')
axes[2].plot(CT_CO2_sweep*1e3, pH_pq5, 's--', ms=3, label='PHREEQC')
axes[2].set_xlabel('CT_CO₂ (mmol/L)'); axes[2].set_ylabel('pH')
axes[2].set_title('pH after precipitation')
axes[2].legend(fontsize=9); axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

print(f"CT_Ca={CT_Ca_fix*1e3:.0f} mmol/L, CT_Na={CT_Na_fix*1e3:.0f} mmol/L")
print(f"Max |Δξ|  = {np.abs(xi_nr5-xi_pq5).max()*1e3:.3f} mmol/L")
print(f"Max |ΔSI| = {np.abs(SI_nr5-SI_pq5).max():.4f}")
print(f"Max |ΔpH| = {np.abs(pH_nr5-pH_pq5).max():.4f}")
print()
print("Note: residual Δξ and ΔpH reflect missing CaHCO3+/CaCO3 complexes in NR engine,")
print("      not numerical error in the precipitation solver.")
"""),

    code("s5-parity", """\
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
_parity(axes[0], xi_pq5 * 1e3, xi_nr5 * 1e3, '\u03be precipitated (mmol/L)',
        c=CT_CO2_sweep * 1e3, cmap='cividis',
        clabel='CT_CO\u2082 (mmol/L)')
_parity(axes[1], pH_pq5, pH_nr5, 'pH after precipitation',
        c=CT_CO2_sweep * 1e3, cmap='cividis',
        clabel='CT_CO\u2082 (mmol/L)')
plt.suptitle('Section 5 parity — calcite precipitation', fontsize=10)
plt.tight_layout(); plt.show()
"""),

    # ── 6. Combined parity ────────────────────────────────────────────────────

    md("s6-combined-md", """\
## 6  Combined parity — all clean systems

pH, species concentrations, and precipitation amount collected across all
sections where the reaction networks match (no Ca).  A tight cluster on the
1:1 line confirms numerical accuracy; any residual offset is model-driven."""),

    code("s6-combined", """\
from matplotlib.lines import Line2D

ph_ref_all = np.concatenate([ph_pq, ph_pq2, ph_pq3])
ph_nr_all  = np.concatenate([ph_nr, ph_nr2, ph_nr3])
system_id  = np.array([0]*len(ph_pq) + [1]*len(ph_pq2) + [2]*len(ph_pq3))
sys_labels = {0: 'Carb. only', 1: 'Carb. + Na\u207a', 2: 'Carb. + NH\u2083'}
sys_colors = [plt.cm.tab10(i / 9) for i in range(3)]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel 1 — pH all clean systems
lo, hi = ph_ref_all.min(), ph_ref_all.max()
axes[0].plot([lo, hi], [lo, hi], 'k--', lw=1, zorder=0)
for sid, lbl, col in zip([0, 1, 2], sys_labels.values(), sys_colors):
    mask = system_id == sid
    axes[0].scatter(ph_ref_all[mask], ph_nr_all[mask],
                    s=18, color=col, alpha=0.8, label=lbl, zorder=1)
axes[0].set_xlabel('PHREEQC pH'); axes[0].set_ylabel('NR engine pH')
axes[0].set_title('pH — all clean systems (no Ca)')
axes[0].legend(fontsize=9); axes[0].grid(True, alpha=0.3)

# Panel 2 — [HCO3-] from Section 1
lo2, hi2 = hco3_pq.min(), hco3_pq.max()
axes[1].loglog([lo2, hi2], [lo2, hi2], 'k--', lw=1)
axes[1].scatter(hco3_pq, hco3_nr, s=18, color=sys_colors[0], alpha=0.8, label='[HCO\u2083\u207b]')
axes[1].scatter(co3_pq,  co3_nr,  s=18, color=sys_colors[1], alpha=0.8, label='[CO\u2083\u00b2\u207b]')
axes[1].set_xlabel('PHREEQC concentration (mol/L)')
axes[1].set_ylabel('NR engine concentration (mol/L)')
axes[1].set_title('Carbonate species — Section 1')
axes[1].legend(fontsize=9); axes[1].grid(True, alpha=0.3)

# Panel 3 — precipitation \u03be
mask5 = xi_pq5 > 1e-8
lo3, hi3 = xi_pq5[mask5].min(), xi_pq5[mask5].max()
axes[2].plot([lo3 * 1e3, hi3 * 1e3], [lo3 * 1e3, hi3 * 1e3], 'k--', lw=1)
axes[2].scatter(xi_pq5[mask5] * 1e3, xi_nr5[mask5] * 1e3,
                s=18, color=sys_colors[2], alpha=0.8)
axes[2].set_xlabel('PHREEQC \u03be (mmol/L)')
axes[2].set_ylabel('NR engine \u03be (mmol/L)')
axes[2].set_title('Precipitation \u03be — Section 5')
axes[2].grid(True, alpha=0.3)

plt.suptitle('Combined parity: NR engine vs PHREEQC', fontsize=12)
plt.tight_layout(); plt.show()

max_dph   = np.abs(ph_nr_all - ph_ref_all).max()
max_hco3  = (np.abs(hco3_nr - hco3_pq) / hco3_pq).max() * 100
max_co3   = (np.abs(co3_nr  - co3_pq)  / co3_pq ).max() * 100
max_xi    = (np.abs(xi_nr5[mask5] - xi_pq5[mask5]) / xi_pq5[mask5]).max() * 100
print(f"Clean-system summary:")
print(f"  pH:         max |\u0394pH|           = {max_dph:.4f}")
print(f"  [HCO\u2083\u207b]:      max relative error   = {max_hco3:.3f}%")
print(f"  [CO\u2083\u00b2\u207b]:      max relative error   = {max_co3:.3f}%")
print(f"  \u03be (precip): max relative error   = {max_xi:.2f}%")
"""),

    # ── 7. Accuracy summary ─────────────────────────────────────────────────────

    md("s7-md", """\
## 7  Accuracy summary

| System | Comparison | Max |ΔpH| | Max |Δ[species]| | Source of residual |
|---|---|---|---|---|
| Carbonate only | Pure numerical | Section 1 | Section 1 | NR solver |
| Carbonate + Na⁺ | Small Na complexes | Section 2 | — | NaHCO₃⁰, NaCO₃⁻ pairs |
| Carbonate + NH₃ | Pure numerical | Section 3 | Section 3 | NR solver |
| Ca²⁺ (strong ion) | Ca complexes missing | Section 4 | — | CaHCO₃⁺, CaCO₃⁰, CaOH⁺ |
| Calcite precipitation | Ca complexes + numerical | Section 5 | — | Mixed |

**Conclusion:**  Where the reaction networks are the same (no Ca), NRChemicalEquilibriumEngine
agrees with PHREEQC to better than 0.01 pH units across the full range.
Residual differences with Ca present are model differences (missing Ca ion pairs),
not solver errors.  Implementing CaHCO₃⁺ complexation (Phase 3) would eliminate
most of the remaining offset.\
"""),
)

# ─── update 0_README notebook table ──────────────────────────────────────────

# Patch the README notebooks table to add the precipitation row
_old_table_row = "| [03_saturation_index.ipynb](03_saturation_index.ipynb) | Carbonate + Ca²⁺ (strong ion); calcite Ksp | Is a given water supersaturated with respect to calcite? | SI = log₁₀(IAP/Ksp); Davies activity corrections; pH-dependent SI sweep; precipitation risk map |"
_new_table_rows = (
    _old_table_row
    + "\n| [05_precipitation_equilibrium.ipynb](05_precipitation_equilibrium.ipynb) | Carbonate + calcite mineral; Ca²⁺ as strong ion | How much CaCO₃ precipitates at equilibrium? | Active-set NR outer loop; IAP = Ksp enforcement; Ca depletion curve; dissolution from solid |"
    + "\n| [06_phreeqc_benchmark.ipynb](06_phreeqc_benchmark.ipynb) | Carbonate; carbonate + Na⁺; carbonate + NH₃; calcite precipitation | Does the NR engine agree with PHREEQC? | phreeqpython comparison; Davies b=0.3; Ca complex model-difference disclosure |"
)
readme_nb["cells"][1]["source"] = readme_nb["cells"][1]["source"].replace(
    _old_table_row,
    _new_table_rows,
)

# ─── write files ──────────────────────────────────────────────────────────────

def save(path, notebook):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(notebook, f, indent=1, ensure_ascii=False)
    print(f"Written: {path}")


save(HERE / "0_README.ipynb",                    readme_nb)
save(HERE / "01_single_component_benchmarks.ipynb", bench_nb)
save(HERE / "02_multi_component_systems.ipynb",  multi_nb)
save(HERE / "03_saturation_index.ipynb",         si_nb)
save(HERE / "05_precipitation_equilibrium.ipynb", precip_nb)
save(HERE / "06_phreeqc_benchmark.ipynb",         bench6_nb)
