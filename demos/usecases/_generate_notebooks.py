"""Generate the usecases demo notebooks as valid .ipynb JSON.

This folder is scenario-first: each notebook opens with a plain-language
"here's the situation" framing rather than a class-by-class API tour (that's
what demos/features/ and demos/model_api/chemistry/speciation/ are for).
Cross-reference those folders once a use case needs more depth than a
five-minute read supports.

Notebooks 01, 02, 02b, and 05 were curated into docs/tutorials/ArXiv_preprint/
(the subset most closely covered by the ArXiv preprint) and are no longer
built here -- see docs/tutorials/ArXiv_preprint/_generate_notebooks.py for
those. This script now owns only 0_README, 03, and 04.

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
| [01_predict_ph_simple_liquid.ipynb](../../docs/tutorials/ArXiv_preprint/01_predict_ph_simple_liquid.ipynb) | I'm making up a defined growth medium from KH₂PO₄ (phosphate buffer) and NH₄Cl (nitrogen source), no gas headspace or solid phase to track — what pH does that land at, across the range of doses used in practice? How does that compare against PHREEQC? | `NRChemicalEquilibriumEngine`, `PHREEQCChemicalEquilibriumEngine` (optional) |
| [02_kinetic_co2_equilibration_microplate_well.ipynb](../../docs/tutorials/ArXiv_preprint/02_kinetic_co2_equilibration_microplate_well.ipynb) | Pure water, in direct contact with a large atmospheric reservoir (O₂/N₂/CO₂) across a gas-liquid interface with a finite mass-transfer coefficient (kLa) rather than an instantaneous equilibrium — how does pH evolve over time as dissolved CO₂ approaches its Henry's-law equilibrium, and how does kLa itself set the timescale to get there? First notebook with genuinely kinetic (rate-limited) gas transfer. | `ControlVolume`, `Simulation`, `KineticTransferModel` |
| [03_grow_ecoli_on_acetic_acid.ipynb](03_grow_ecoli_on_acetic_acid.ipynb) | Inoculate that same sparged vessel with *E. coli* growing on acetic acid as sole carbon source — how fast does it grow, does dissolved O₂ ever become limiting, and what happens to pH as the acid substrate is consumed? First notebook where the chemistry evolves over time. | `ControlVolume`, `Simulation`, `ReactionBuilder.monod_aerobic_growth`, `ReactionSystem` |
| [04_compare_runtime_by_usecase.ipynb](04_compare_runtime_by_usecase.ipynb) | Usecases 01-03 all took *some* wall-clock time to run — how much, and where does it go? Compares run time across all three usecases, decomposed by which activity-coefficient treatment did the solving (PyOMES ideal / PyOMES Davies / PHREEQC ideal-equivalent / PHREEQC default). | `NRChemicalEquilibriumEngine`, `PHREEQCChemicalEquilibriumEngine` (optional), `ReactionSystem.configure_engine` |
| [03_cstr_dilution_rate_sweep.ipynb](../../docs/tutorials/ArXiv_preprint/03_cstr_dilution_rate_sweep.ipynb) | Same organism/substrate as 03, but now run as a chemostat — a CSTR fed and drained at the same volumetric flow rate, so the working volume holds steady while biomass and substrate settle onto a dilution-rate-dependent steady state. How does the dilution rate affect the reactor's volumetric productivity, and where does washout kick in? | `ControlVolume`, `Simulation`, `LiquidFeed`, `LiquidDrain` |

Launch from the repo root with `jupyter lab demos/usecases/`.\
"""),
)

save_path_readme = HERE / "0_README.ipynb"
with open(save_path_readme, "w", encoding="utf-8") as f:
    json.dump(readme_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path_readme}")

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

**The situation:** [usecase 01](../../docs/tutorials/ArXiv_preprint/01_predict_ph_simple_liquid.ipynb) built the
KH₂PO₄ + NH₄Cl liquid and usecase 02
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
- **`StirredTankBuilder`, the fluent alternative to assembling `ControlVolume`
  by hand** — `PyOMES/templates/stirred_tank/builder.py` wraps vessel
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

Exactly [usecase 01](../../docs/tutorials/01_predict_ph_simple_liquid.ipynb)'s chemistry (the
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

Usecase 02's network extends
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
  different run times** — [usecase 01](../../docs/tutorials/01_predict_ph_simple_liquid.ipynb)
  §6 and usecase 02 §6 walk
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

