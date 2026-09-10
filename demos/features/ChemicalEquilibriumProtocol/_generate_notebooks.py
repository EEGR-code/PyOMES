"""Generate the ChemicalEquilibriumProtocol feature notebooks as valid .ipynb JSON.

This folder covers the `ChemicalEquilibriumEngineProtocol` family — the
three peer engine implementations (Bisection, NR, PHREEQC) satisfy the same
solve() contract but differ in scope and calling convention. Each notebook
covers one engine specifically; see 10_engine_protocol_hierarchy.ipynb in
demos/model_api/chemistry/speciation/ for a cross-engine comparison.

Run once from the repo root:
    python demos/features/ChemicalEquilibriumProtocol/_generate_notebooks.py
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
#  0  Architecture overview (pure markdown, no code — matches the
#     model_api/chemistry/speciation/0_README.ipynb convention)
# ═══════════════════════════════════════════════════════════════════════════

overview_nb = nb(
    md("title", """\
# `ChemicalEquilibriumProtocol` — Architecture Overview

`ChemicalEquilibriumEngineProtocol` is PyOMES's structural contract for
"solve aqueous equilibrium chemistry, return an `EquilibriumResult`."
Three peer engines satisfy it — `BisectionChemicalEquilibriumEngine`,
`NRChemicalEquilibriumEngine`, `PHREEQCChemicalEquilibriumEngine` — and
**none of them is "the" default**. Each trades off scope, calling
convention, and dependency footprint differently; picking the right one
is a decision about your chemistry's shape, not a question of which is
"newest" or "most capable" in every respect.

This folder's notebooks each cover one engine's instantiation/call
mechanics in isolation. This notebook is the map: what the protocol
hierarchy actually is, how the three engines compare, and which one to
reach for.\
"""),

    md("notebooks", """\
## Notebooks in this folder

| Notebook | Engine | Central question |
|---|---|---|
| [01_bisection_engine_basics.ipynb](01_bisection_engine_basics.ipynb) | `BisectionChemicalEquilibriumEngine` | How do I construct it via `from_reactions()` and call `solve(CT_TIC=..., ...)`? |
| [02_nr_engine_basics.ipynb](02_nr_engine_basics.ipynb) | `NRChemicalEquilibriumEngine` | How do I use `solve(totals={...})`, and what does gas-liquid folding / precipitation actually look like? |
| [03_phreeqc_engine_basics.ipynb](03_phreeqc_engine_basics.ipynb) | `PHREEQCChemicalEquilibriumEngine` | How do I construct it via `component_map` instead of `from_reactions()`, and what does PHREEQC's species-name translation look like? |

Launch from the repo root with
`jupyter lab demos/features/ChemicalEquilibriumProtocol/`.\
"""),

    md("hierarchy", """\
## The protocol hierarchy: three capability tiers

`ChemicalEquilibriumEngineProtocol` is the **black-box** tier — every
engine satisfies at least this much. Two further tiers add capability
without breaking the tier below (`GrayBoxEngineProtocol` and
`WhiteBoxEngineProtocol` both structurally extend the black-box
interface):

| Tier | Protocol | Adds beyond the tier below | Method(s) added |
|---|---|---|---|
| Black box | `ChemicalEquilibriumEngineProtocol` | — (base tier) | `solve()`, `algebraic_species()`, `reset_cache()`, `reset_counters()`, `n_solve_calls` |
| Gray box | `GrayBoxEngineProtocol` | Total sensitivity at the converged point | `jacobian_dz_dy()` |
| White box | `WhiteBoxEngineProtocol` | Residual + split partial Jacobians, for DAE coupling | `residual()`, `jacobian_dg_dz()`, `jacobian_dg_dy()` |

**Which engines satisfy which tier, natively:**

| Engine | Black box | Gray box | White box |
|---|---|---|---|
| `BisectionChemicalEquilibriumEngine` | ✅ | ❌ | ❌ |
| `NRChemicalEquilibriumEngine` | ✅ | ✅ (via `retain_jacobian=True`) | ✅ (via `retain_jacobian=True`, and only when the tableau has no folded gas-liquid secondaries) |
| `PHREEQCChemicalEquilibriumEngine` | ✅ | ❌ | ❌ |

Any black-box engine can be *promoted* to gray-box without modifying it,
by wrapping it in `NumericalGradientEquilibriumEngine` — it computes
`jacobian_dz_dy()` via central finite differences around whatever
`solve(totals=...)` call the wrapped engine already supports. This is
how `BisectionChemicalEquilibriumEngine` and
`PHREEQCChemicalEquilibriumEngine` gain gray-box capability in practice;
there is no equivalent finite-difference promotion to white-box (that
tier needs the analytic residual/Jacobian structure only `NRChemicalEquilibriumEngine`
actually builds).

See
[`../../model_api/chemistry/speciation/10_engine_protocol_hierarchy.ipynb`](../../model_api/chemistry/speciation/10_engine_protocol_hierarchy.ipynb)
for a fully executed, code-level walkthrough of all three tiers,
including `NumericalGradientEquilibriumEngine` wrapping and a
`jacobian_dz_dy()` cost comparison.\
"""),

    md("declaration-side", """\
## The other half: how chemistry gets *into* an engine

The protocol hierarchy above is about the **solve side** — what an
engine returns and what capability tiers it exposes. A separate,
related unification covers the **declaration side**: `EquilibriumReaction`,
`HenryEquilibrium`, `KspEquilibrium`, and `RaoultEquilibrium` all satisfy
one shared `EquilibriumConstraint` protocol
(`EQUILIBRIUM_CONSTRAINT_UNIFICATION`), so a single flat list of declared
constraints can feed `from_reactions()`. This is exactly how
`BisectionChemicalEquilibriumEngine`/`NRChemicalEquilibriumEngine` accept
chemistry — but note that `PHREEQCChemicalEquilibriumEngine` opts out of
this entirely: it has no `from_reactions()` at all, and takes chemistry
via `component_map` into PHREEQC's own database instead (see
[03_phreeqc_engine_basics.ipynb](03_phreeqc_engine_basics.ipynb), Section 1).\
"""),

    md("comparison", """\
## Engine comparison

| Property | `BisectionChemicalEquilibriumEngine` | `NRChemicalEquilibriumEngine` | `PHREEQCChemicalEquilibriumEngine` |
|---|---|---|---|
| Construction | `from_reactions(equilibrium_reactions, ...)` | `from_reactions(equilibrium_reactions, ...)` | `__init__(components, component_map=...)` — no `from_reactions()` |
| Chemistry source | Declared `EquilibriumReaction`/`EquilibriumConstraint` list | Declared `EquilibriumReaction`/`EquilibriumConstraint` list | PHREEQC's own thermodynamic database |
| `solve()` totals | Named kwargs (`CT_TIC=`, `CT_NH_T=`, ...) | `totals={master_id: mol_L}` + `strong_ions={...}` | `totals={component_id: mol_L}` (no `strong_ions=` split) |
| Method | 1-D bisection on the charge balance | Full Newton-Raphson in log-activity space | External PHREEQC solve (`phreeqpython`) |
| Reaction networks | Single/independent acid-base ladders only | Arbitrary cross-component networks | Whatever PHREEQC's database supports |
| Gas-liquid folding | ❌ (diverted to `cross_phase_constraints`, never solved) | ✅ (simultaneous, via folded tableau secondaries) | N/A (PHREEQC's own gas-phase handling, outside this comparison) |
| Precipitation (Ksp) | ❌ (diverted, never solved) | ✅ (auto-detected, nested active-set loop) | Whatever PHREEQC's database supports |
| `charge_residual` in result | Always `None` | Populated (last Newton residual row) | Not populated |
| External dependency | None | None | `phreeqpython` (optional: `pip install PyOMES[phreeqc]`) |
| Warmstart | `logH_warmstart`/`I_warmstart` cache | Full log-activity vector cache | PHREEQC `Solution` object reused via `REACTION` — see 03's documented drift gotcha |

Every engine returns the same `EquilibriumResult` type and satisfies
`ChemicalEquilibriumEngineProtocol`, so code written against the
black-box interface (`solve()`, `algebraic_species()`, `reset_cache()`,
`reset_counters()`, `n_solve_calls`) is portable across all three —
the differences above are what you hit the moment you need
engine-specific *construction* or a *capability* only some of them have.\
"""),

    md("which-one", """\
## Which engine should I use?

- **Your chemistry is one acid-base ladder, or a few independent ones**
  (no species whose mass balance couples to more than one "total"), and
  you don't need gas-liquid/precipitation folding: `BisectionChemicalEquilibriumEngine`.
  Simplest, no external dependency, historically validated.
- **Your chemistry has cross-component coupling, gas-liquid partition
  terms that must be solved simultaneously with the acid-base system,
  or mineral precipitation**: `NRChemicalEquilibriumEngine`. This is
  also the engine to reach for if you need gray/white-box Jacobian
  access natively (`retain_jacobian=True`).
- **You want PHREEQC's own thermodynamic database and ion-pair library**
  (more extensive than PyOMES's declared-reaction chemistry, at the cost
  of an external dependency and no `from_reactions()`), or you're
  cross-validating PyOMES's own engines against an independent reference:
  `PHREEQCChemicalEquilibriumEngine`.

All three can be promoted to gray-box via `NumericalGradientEquilibriumEngine`
if you need `jacobian_dz_dy()` and the engine doesn't have it natively.\
"""),

    md("cross-references", """\
## Cross-references

- [`docs/dev/ideas/CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md`](../../../docs/dev/ideas/CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md) —
  full protocol-hierarchy design record; engine/solver split; §18 naming
  history (`SpeciationEngine` → `ChemicalEquilibriumEngine` →
  `BisectionChemicalEquilibriumEngine`).
- [`docs/dev/implementation/shipped/EQUILIBRIUM_CONSTRAINT_UNIFICATION.md`](../../../docs/dev/implementation/shipped/EQUILIBRIUM_CONSTRAINT_UNIFICATION.md) —
  the declaration-side unification (`EquilibriumConstraint`,
  `classify_equilibrium_constraint()`).
- [`docs/dev/implementation/shipped/LAYER1_GAP_CLOSURE.md`](../../../docs/dev/implementation/shipped/LAYER1_GAP_CLOSURE.md) —
  gas-liquid/precipitation folding into the NR tableau; the CP6 engine
  rename.
- [`../../model_api/chemistry/speciation/10_engine_protocol_hierarchy.ipynb`](../../model_api/chemistry/speciation/10_engine_protocol_hierarchy.ipynb) —
  executed, code-level black/gray/white-box walkthrough.
- [`../../model_api/chemistry/speciation/02_multi_component_systems.ipynb`](../../model_api/chemistry/speciation/02_multi_component_systems.ipynb) —
  Bisection-vs-NR accuracy comparison across conditions.
- [`../../model_api/chemistry/speciation/06_phreeqc_benchmark.ipynb`](../../model_api/chemistry/speciation/06_phreeqc_benchmark.ipynb) —
  NR-vs-PHREEQC accuracy comparison.\
"""),
)

save_path_overview = HERE / "0_README.ipynb"
with open(save_path_overview, "w", encoding="utf-8") as f:
    json.dump(overview_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path_overview}")


SETUP = """\
import sys
from pathlib import Path

def _find_repo():
    for p in [Path.cwd(), *Path.cwd().parents]:
        if (p / "pyproject.toml").exists():
            return p
    raise RuntimeError("Run from inside the PyOMES repo")

sys.path.insert(0, str(_find_repo() / "models"))

from PyOMES.chemistry.common_species import (
    H2O, H_plus, OH_minus,
    CO2, HCO3_minus, CO3_2minus,
    NH3, NH4_plus,
)
from PyOMES.reactions.equilibrium import EquilibriumReaction
from PyOMES.reactions.stoichiometry import StoichiometryEntry
from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine

def _e(sp, coeff, phase="liquid"):
    return StoichiometryEntry(species=sp, phase=phase, coefficient=coeff)

print("Imports OK")
"""

basics_nb = nb(
    md("title", """\
# `BisectionChemicalEquilibriumEngine` — Instantiation and Use

Three engines satisfy `ChemicalEquilibriumEngineProtocol`:
`BisectionChemicalEquilibriumEngine`, `NRChemicalEquilibriumEngine`, and
`PHREEQCChemicalEquilibriumEngine`. None is "the" default — pick the one
that matches your chemistry's shape. This notebook covers
`BisectionChemicalEquilibriumEngine`: PyOMES's original speciation engine, a
1-D bisection solver over the charge balance, scanning `pH_min`–`pH_max`
for the root. It solves single acid-base ladders (and independent ladders
that only interact through the charge balance, e.g. carbonate + ammonia)
declared as `EquilibriumReaction` instances.

**Contrast with `NRChemicalEquilibriumEngine`:** the NR engine builds a
log-linear tableau and solves the full multi-variable Newton system, so it
handles arbitrary cross-component networks (species whose mass balance
couples to more than one "total") and gas-liquid/solid-liquid folding.
`BisectionChemicalEquilibriumEngine` cannot do either of those — it is the right
choice when your chemistry is a small number of decoupled acid-base ladders
and you want the simpler, historically-validated bisection path. See
[`02_nr_engine_basics.ipynb`](02_nr_engine_basics.ipynb) for the NR engine's
own instantiation/call conventions, including the two capabilities this
engine structurally cannot have, or
[`03_phreeqc_engine_basics.ipynb`](03_phreeqc_engine_basics.ipynb) for the
PHREEQC-backed engine, which forgoes declared `EquilibriumReaction`
networks entirely in favour of PHREEQC's own database.

This notebook focuses on the mechanics of the engine itself: constructing
it, calling `solve()`, and reading back an `EquilibriumResult`. For
engine-to-engine accuracy comparisons against `NRChemicalEquilibriumEngine`,
see
[`../../model_api/chemistry/speciation/02_multi_component_systems.ipynb`](../../model_api/chemistry/speciation/02_multi_component_systems.ipynb)
and
[`../../model_api/chemistry/speciation/10_engine_protocol_hierarchy.ipynb`](../../model_api/chemistry/speciation/10_engine_protocol_hierarchy.ipynb).\
"""),

    code("setup", SETUP),

    # ── 1. Instantiation ──────────────────────────────────────────────────────

    md("s1-md", """\
## 1  Instantiation via `from_reactions()`

`BisectionChemicalEquilibriumEngine.from_reactions()` classifies each declared
`EquilibriumReaction` into `"water"`, `"acid"`, or `"cation_acid"` from its
stoichiometry, and folds them into an internal `EquilibriumSet`
(`engine._equilibrium_set`). Reactions that share a `total_id` (e.g. a
polyprotic ladder) are merged into a single `EquilibriumDef` automatically.

Any item that is *not* a single-phase acid-base equilibrium (a gas-liquid
`HenryEquilibrium`, a solid-liquid `KspEquilibrium`, etc.) is **not silently
dropped** — it is diverted to `engine.cross_phase_constraints` instead of
being folded into the tableau. Section 5 below shows why that matters.\
"""),

    code("s1-code", """\
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

engine = BisectionChemicalEquilibriumEngine.from_reactions([water, co2_first, co2_second])

print("Engine constructed:  ", type(engine).__name__)
print("EquilibriumSet names:", [d.name for d in engine._equilibrium_set])
print("Cross-phase constraints (empty — no gas/solid items declared):",
      engine.cross_phase_constraints)
"""),

    # ── 2. solve() and total-key mapping ──────────────────────────────────────

    md("s2-md", """\
## 2  `solve()` and the total-key mapping

Unlike `NRChemicalEquilibriumEngine.solve(totals={...})`,
`BisectionChemicalEquilibriumEngine.solve()` takes totals as **named keyword
arguments**, using a fixed mapping from species id to a conventional total
key (`_TOTAL_KEY_OVERRIDES` in `engine.py`):

| Species id | Total kwarg |
|---|---|
| `CO2` | `CT_TIC` |
| `NH3` / `NH4+` | `CT_NH_T` |
| `H3PO4` | `CT_P` |
| `SO4--` | `CT_SO4` |

`solve()` returns an immutable `EquilibriumResult` — the same return type
`NRChemicalEquilibriumEngine` and `PHREEQCChemicalEquilibriumEngine` use.
`species_mol_L` is restricted to a fixed canonical tuple this engine has
always written back (`_CANONICAL_WRITEBACK_SPECIES`); anything else lands in
`.extra`.

**Field coverage differs by engine.** `charge_residual` is populated by the
NR solver (it's the last row of the Newton residual) but `solve_acid_base`
and `solve_from_equilibrium_set` — the two paths behind
`BisectionChemicalEquilibriumEngine.solve()` — never set it, so it's always `None`
here. `pH`, `species_mol_L`, and `ionic_strength` are populated on every
engine; don't assume every `EquilibriumResult` field is.\
"""),

    code("s2-code", """\
result = engine.solve(CT_TIC=0.01)   # 10 mmol/L total inorganic carbon

print(f"pH               = {result.pH:.4f}")
print(f"species_mol_L    = {result.species_mol_L}")
print(f"charge_residual  = {result.charge_residual!r}  (always None for this engine)")
print(f"ionic_strength   = {result.ionic_strength:.4e} mol/L")
print(f"alphas           = {result.alphas}")
"""),

    md("s2-gotcha-md", """\
### 2.1  Gotcha: NR-style `totals=` is silently ignored

`solve()` accepts `**kwargs` and only ever looks for the specific total
keywords above — it never reads a `totals=` dict. Passing NR-style
`totals={"CO2": 0.01}` compiles and runs, but every total defaults to `0.0`
because `CT_TIC` was never set. There is no error; the result is just wrong.
This is exactly the kind of engine-convention mismatch this notebook exists
to make obvious.\
"""),

    code("s2-gotcha-code", """\
result_wrong = engine.solve(totals={"CO2": 0.01})   # WRONG for this engine
result_right = engine.solve(CT_TIC=0.01)             # correct

print(f"pH with totals={{'CO2': 0.01}} (ignored): {result_wrong.pH:.4f}  ← defaults to CT_TIC=0")
print(f"pH with CT_TIC=0.01 (correct):           {result_right.pH:.4f}")
"""),

    # ── 3. Strong ions and activity model ─────────────────────────────────────

    md("s3-md", """\
## 3  Strong ions and the Davies activity model

Strong ions are also passed as direct keyword arguments (`CT_Na`, `CT_Cl`,
`CT_Ca`, ... — see the list in `engine.py`), not a nested `strong_ions=`
dict. Activity corrections are configured once at construction time via
`use_activity=` / `activity_model=`, exactly as for `NRChemicalEquilibriumEngine`.\
"""),

    code("s3-code", """\
baseline = engine.solve(CT_TIC=0.01)
dosed    = engine.solve(CT_TIC=0.01, CT_Na=0.01)

print(f"pH baseline (no strong ions):  {baseline.pH:.4f}")
print(f"pH with 10 mmol/L Na+ dosing:  {dosed.pH:.4f}")

# Like-for-like ideal-vs-Davies comparison: same totals AND same strong-ion
# composition (10 mmol/L NaCl — net-neutral charge), only the activity
# model differs between the two engines.
engine_davies = BisectionChemicalEquilibriumEngine.from_reactions(
    [water, co2_first, co2_second], use_activity=True, activity_model="davies",
)
ideal_nacl  = engine.solve(CT_TIC=0.01, CT_Na=0.01, CT_Cl=0.01)
davies_nacl = engine_davies.solve(CT_TIC=0.01, CT_Na=0.01, CT_Cl=0.01)

print(f"\\npH ideal (γ=1),   +10 mmol/L NaCl: {ideal_nacl.pH:.4f}")
print(f"pH Davies-corrected, +10 mmol/L NaCl: {davies_nacl.pH:.4f}")
print(f"Ionic strength (Davies run):          {davies_nacl.ionic_strength*1e3:.2f} mmol/L")
"""),

    # ── 4. Multi-component (independent ladders) ──────────────────────────────

    md("s4-md", """\
## 4  Multi-component: independent acid-base ladders

`from_reactions()` can bind more than one acid family at once — each
`total_id` group becomes its own `EquilibriumDef`, and the bisection solves
the *combined* charge balance across all of them simultaneously. This is
different from cross-component *coupling* (e.g. a species whose mass balance
spans two totals, which this engine cannot represent): here, carbonate and
ammonia only interact through the shared pH.\
"""),

    code("s4-code", """\
nh4 = EquilibriumReaction(
    stoichiometry=[_e(NH4_plus, -1), _e(NH3, +1), _e(H_plus, +1)],
    log_K=-9.25, total_id="NH3", label="nh4",
)

engine_multi = BisectionChemicalEquilibriumEngine.from_reactions([water, co2_first, co2_second, nh4])
result_multi = engine_multi.solve(CT_TIC=0.01, CT_NH_T=0.04)

print("EquilibriumSet names:", [d.name for d in engine_multi._equilibrium_set])
print(f"pH                  = {result_multi.pH:.4f}")
print(f"species_mol_L       = {result_multi.species_mol_L}")
print()
print("For a full NR-vs-charge-balance accuracy sweep on this exact carbonate")
print("+ ammonia system, see 02_multi_component_systems.ipynb in the")
print("model_api/chemistry/speciation/ demo family.")
"""),

    # ── 5. Cross-phase constraints gotcha ─────────────────────────────────────

    md("s5-md", """\
## 5  Gotcha: gas-liquid/solid-liquid items are diverted, not solved

If a `HenryEquilibrium` (or `KspEquilibrium`/`RaoultEquilibrium`) is included
in the reaction list, `from_reactions()` classifies it as `"gas_liquid"` (or
`"solid_liquid"`) and routes it to `engine.cross_phase_constraints` instead
of folding it into the acid-base tableau. `solve()` never reads
`cross_phase_constraints` — the item is exposed for a caller such as
`KineticGasLiquidLink` to pick up, not solved by this engine at all.

Declaring a Henry term and expecting it to affect the liquid-phase pH here
is a real mistake this API shape invites; the item is preserved (not
discarded) precisely so a caller can detect it and route it elsewhere —
but `BisectionChemicalEquilibriumEngine.solve()` itself will not fold it in. Use
`NRChemicalEquilibriumEngine` (see `LAYER1_GAP_CLOSURE`'s gas-liquid folding)
when you need the Henry/Raoult term to enter the same solve.\
"""),

    code("s5-code", """\
from PyOMES.chemistry import HenryEquilibrium

co2_henry = HenryEquilibrium(H_ref=3.4e-4, dlnH=2400.0, gas_species="CO2", liquid_species="CO2",
                              label="co2_henry")

engine_mixed = BisectionChemicalEquilibriumEngine.from_reactions([water, co2_first, co2_second, co2_henry])

print("Cross-phase constraints found:", [c.label for c in engine_mixed.cross_phase_constraints])

result_mixed = engine_mixed.solve(CT_TIC=0.01)
print(f"pH with Henry term declared:  {result_mixed.pH:.4f}")
print(f"pH without it (Section 2):    {result.pH:.4f}")
print("(identical — the Henry term never entered the charge balance)")
"""),

    # ── 6. Convenience helpers ─────────────────────────────────────────────────

    md("s6-md", """\
## 6  Convenience helpers

- `get_CO2aq_from_totals(**kwargs)` — one-shot CO2(aq) lookup, checking both
  `species_mol_L` and `extra` since "CO2aq" is not in the canonical
  writeback tuple.
- `logH_warmstart` / `I_warmstart` — the current warm-start cache (used to
  seed the next bisection scan when `use_warmstart=True`, the default).
- `reset_cache()` / `reset_counters()` — clear the warm-start cache / the
  `n_solve_calls` counter.
- `algebraic_species()` — species ids written to `phase.n_mol` by
  `_refresh_derived`; the algebraic state `z` in a DAE formulation.\
"""),

    code("s6-code", """\
print("get_CO2aq_from_totals:", engine.get_CO2aq_from_totals(CT_TIC=0.01))
print("logH_warmstart:        ", engine.logH_warmstart)
print("I_warmstart:           ", engine.I_warmstart)
print("n_solve_calls so far:  ", engine.n_solve_calls)

engine.reset_cache()
engine.reset_counters()
print("\\nafter reset_cache()/reset_counters():")
print("logH_warmstart:        ", engine.logH_warmstart)
print("n_solve_calls:         ", engine.n_solve_calls)

print("\\nalgebraic_species():  ", sorted(engine.algebraic_species()))
"""),

    md("summary", """\
## Summary

| Topic | Key takeaway |
|---|---|
| Construction | `from_reactions()` classifies each item; non-acid-base items go to `cross_phase_constraints`, not the tableau |
| `solve()` totals | Fixed keyword mapping (`CT_TIC`, `CT_NH_T`, `CT_P`, `CT_SO4`) — **not** `totals={...}` |
| Strong ions | Direct kwargs (`CT_Na`, `CT_Ca`, ...), not a nested dict |
| Result type | `EquilibriumResult`, same as the NR and PHREEQC engines; `species_mol_L` restricted to a fixed canonical tuple |
| Multi-component | Independent acid ladders solve together via the shared charge balance; true cross-component coupling needs `NRChemicalEquilibriumEngine` |
| Gas-liquid items | Diverted to `cross_phase_constraints`, never folded into `solve()` — a common source of silent no-ops |
| Helpers | `get_CO2aq_from_totals`, warm-start cache properties, `reset_cache`/`reset_counters`, `algebraic_species()` |\
"""),
)

save_path = HERE / "01_bisection_engine_basics.ipynb"
with open(save_path, "w", encoding="utf-8") as f:
    json.dump(basics_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path}")


# ═══════════════════════════════════════════════════════════════════════════
#  02  NRChemicalEquilibriumEngine basics
# ═══════════════════════════════════════════════════════════════════════════

SETUP_NR = """\
import sys
from pathlib import Path

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
from PyOMES.chemistry import HenryEquilibrium
from PyOMES.reactions.equilibrium import EquilibriumReaction
from PyOMES.reactions.stoichiometry import StoichiometryEntry
from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

def _e(sp, coeff, phase="liquid"):
    return StoichiometryEntry(species=sp, phase=phase, coefficient=coeff)

print("Imports OK")
"""

nr_nb = nb(
    md("title", """\
# `NRChemicalEquilibriumEngine` — Instantiation and Use

Three engines satisfy `ChemicalEquilibriumEngineProtocol`:
`BisectionChemicalEquilibriumEngine`, `NRChemicalEquilibriumEngine`, and
`PHREEQCChemicalEquilibriumEngine`. None is "the" default — pick the one
that matches your chemistry's shape. This notebook covers
`NRChemicalEquilibriumEngine`: a Newton-Raphson solver over a log-linear
tableau in log-activity space. Unlike
`BisectionChemicalEquilibriumEngine` (see
[`01_bisection_engine_basics.ipynb`](01_bisection_engine_basics.ipynb)),
it handles **arbitrary cross-component networks** (species whose mass
balance couples to more than one "total"), **gas-liquid folding**
(Henry/Raoult rows solved simultaneously with the acid-base system,
not diverted), and **solid-liquid precipitation** (Ksp, via a nested
active-set loop).

This notebook focuses on the mechanics of the engine itself:
constructing it, calling `solve()`, and reading back an
`EquilibriumResult` — including the two capabilities
`BisectionChemicalEquilibriumEngine` structurally cannot have. See
[`03_phreeqc_engine_basics.ipynb`](03_phreeqc_engine_basics.ipynb) for
the third engine, which forgoes declared reaction networks entirely.
For a deep accuracy comparison against the Bisection engine and the full
protocol hierarchy (black/gray/white-box), see
[`../../model_api/chemistry/speciation/02_multi_component_systems.ipynb`](../../model_api/chemistry/speciation/02_multi_component_systems.ipynb)
and
[`../../model_api/chemistry/speciation/10_engine_protocol_hierarchy.ipynb`](../../model_api/chemistry/speciation/10_engine_protocol_hierarchy.ipynb).\
"""),

    code("setup", SETUP_NR),

    # ── 1. Instantiation ──────────────────────────────────────────────────────

    md("s1-md", """\
## 1  Instantiation via `from_reactions()`

`NRChemicalEquilibriumEngine.from_reactions()` takes **one flat list**
of any `EquilibriumConstraint`-conforming item — `EquilibriumReaction`,
`HenryEquilibrium`, `KspEquilibrium`, `RaoultEquilibrium` — and
auto-classifies each via `classify_equilibrium_constraint()`:
acid-base items build the Newton tableau graph; solid-liquid items are
auto-detected as precipitation reactions (no separate
`precipitation_reactions=` kwarg required); gas-liquid items with a
real `log_K` are **folded directly into the tableau** as gas-phase
secondaries. This is the opposite of
`BisectionChemicalEquilibriumEngine`'s behaviour, which diverts
everything except single-phase acid-base items to
`cross_phase_constraints` — see Section 5.\
"""),

    code("s1-code", """\
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

engine = NRChemicalEquilibriumEngine.from_reactions([water, co2_first, co2_second], T_K=298.15)

print("Engine constructed:  ", type(engine).__name__)
print("Tableau masters:     ", engine._tableau.masters)
print("gas_liquid_species(): ", engine.gas_liquid_species())
"""),

    # ── 2. solve() and the totals= dict ────────────────────────────────────────

    md("s2-md", """\
## 2  `solve()` and the `totals=` dict

Unlike `BisectionChemicalEquilibriumEngine.solve(CT_TIC=..., CT_NH_T=...)`,
`NRChemicalEquilibriumEngine.solve()` takes totals as a single **`totals=`
dict, keyed by master id** (the species each connected component is
tracked against — see `engine._tableau.masters`), plus a separate
`strong_ions=` dict for non-reactive charge carriers (`CT_Na`, `CT_Cl`, ...).
Both engines return the same `EquilibriumResult` type, but
`charge_residual` — `None` for the Bisection engine — is populated here:
it's literally the last row of the Newton residual, so it always exists
for this solver.\
"""),

    code("s2-code", """\
result = engine.solve(totals={"CO2": 0.01}, strong_ions={})   # 10 mmol/L total inorganic carbon

print(f"pH               = {result.pH:.4f}")
print(f"species_mol_L    = {result.species_mol_L}")
print(f"charge_residual  = {result.charge_residual:.2e}  (populated — last Newton residual row)")
print(f"ionic_strength   = {result.ionic_strength:.4e} mol/L")
"""),

    md("s2-gotcha-md", """\
### 2.1  Gotcha: Bisection-style named kwargs are silently ignored

`solve()` only ever reads `totals=`/`strong_ions=`/`phases=`/`T_K=`/
`V_liq_L=`/`V_gas_L=` from its `**kwargs` — a Bisection-style
`CT_TIC=0.01` is accepted syntactically (absorbed into `**kwargs` and
never read) but has no effect; `totals` defaults to `{}` and every
total is `0.0`. Exactly the mirror image of Section 2.1 in the
Bisection notebook.\
"""),

    code("s2-gotcha-code", """\
result_wrong = engine.solve(CT_TIC=0.01)                        # WRONG for this engine
result_right = engine.solve(totals={"CO2": 0.01}, strong_ions={})  # correct

print(f"pH with CT_TIC=0.01 (ignored):        {result_wrong.pH:.4f}  ← defaults to totals={{}}")
print(f"pH with totals={{'CO2': 0.01}} (correct): {result_right.pH:.4f}")
"""),

    # ── 3. Strong ions and activity model ─────────────────────────────────────

    md("s3-md", """\
## 3  Strong ions and the Davies activity model

Strong ions go in the `strong_ions=` dict (`{"CT_Na": 0.01, ...}`), not
direct kwargs. Activity corrections are configured once at construction
time via `use_activity=`/`activity_model=`, exactly as for
`BisectionChemicalEquilibriumEngine`.\
"""),

    code("s3-code", """\
baseline = engine.solve(totals={"CO2": 0.01}, strong_ions={})
dosed    = engine.solve(totals={"CO2": 0.01}, strong_ions={"CT_Na": 0.01})

print(f"pH baseline (no strong ions):  {baseline.pH:.4f}")
print(f"pH with 10 mmol/L Na+ dosing:  {dosed.pH:.4f}")

engine_davies = NRChemicalEquilibriumEngine.from_reactions(
    [water, co2_first, co2_second], use_activity=True, activity_model="davies", T_K=298.15,
)
strong_nacl = {"CT_Na": 0.01, "CT_Cl": 0.01}
ideal_nacl  = engine.solve(totals={"CO2": 0.01}, strong_ions=strong_nacl)
davies_nacl = engine_davies.solve(totals={"CO2": 0.01}, strong_ions=strong_nacl)

print(f"\\npH ideal (γ=1),   +10 mmol/L NaCl: {ideal_nacl.pH:.4f}")
print(f"pH Davies-corrected, +10 mmol/L NaCl: {davies_nacl.pH:.4f}")
print(f"Ionic strength (Davies run):          {davies_nacl.ionic_strength*1e3:.2f} mmol/L")
"""),

    # ── 4. Precipitation (Ksp) ─────────────────────────────────────────────────

    md("s4-md", """\
## 4  Solid-liquid: precipitation (Ksp) is auto-classified

A solid-liquid `EquilibriumReaction` — one `StoichiometryEntry(phase="solid")`
(the mineral) plus one or more `phase="liquid"` dissolved products, with
`log_K` set to `log10(Ksp)` — is auto-detected from the same flat list
passed to `from_reactions()` and routed into `engine._precipitation_reactions`,
which drives an outer active-set loop around the inner Newton solve. No
separate `precipitation_reactions=` kwarg is needed (that kwarg still
exists but is deprecated — see `EQUILIBRIUM_CONSTRAINT_UNIFICATION` CP2).\
"""),

    code("s4-code", """\
CaCO3_solid = Species(id="CaCO3(s)", atoms={"Ca": 1, "C": 1, "O": 3}, charge=0, MW=100.086)
calcite = EquilibriumReaction(
    stoichiometry=[
        _e(CaCO3_solid, -1, phase="solid"),
        _e(Ca_plus_plus, +1),
        _e(CO3_2minus, +1),
    ],
    log_K=-8.48, label="calcite",   # log10(Ksp) at 25 C
)

engine_precip = NRChemicalEquilibriumEngine.from_reactions(
    [water, co2_first, co2_second, calcite], use_activity=True, activity_model="davies", T_K=298.15,
)
print("Auto-detected precipitation reactions:",
      [r.label for r in engine_precip._precipitation_reactions])

# Supersaturated w.r.t. calcite: high CT_CO2 + Ca2+ dosing
out_precip = engine_precip.solve(
    totals={"CO2": 0.010}, strong_ions={"CT_Ca": 0.002, "CT_Na": 0.005},
)
xi = out_precip.extra["minerals_xi_mol_L"]["calcite"]
si = out_precip.saturation_indices["calcite"]
print(f"xi (mol/L precipitated) = {xi:.6f}")
print(f"saturation index (~0 at equilibrium) = {si:.4f}")
"""),

    # ── 5. Gas-liquid folding ──────────────────────────────────────────────────

    md("s5-md", """\
## 5  Gas-liquid folding: the capability Bisection cannot have

A fully-parameterized `HenryEquilibrium` (`gas_species`/`liquid_species`
set, so it carries a real `log_K`) is folded directly into the Newton
tableau as a gas-phase secondary attached to the acid-base component its
liquid form belongs to — solved **simultaneously** with the acid-base
system, not as a sequential post-speciation correction
(`KineticGasLiquidLink`'s SNIA path). This is the opposite of
`BisectionChemicalEquilibriumEngine.from_reactions()`, which diverts the
same `HenryEquilibrium` to `cross_phase_constraints` and never solves it
(Section 5 of the Bisection notebook).

Folding a gas-liquid component changes `solve()`'s calling contract:
`totals` for that component must be the **total across both phases**
(mol/L liquid-volume-normalized), and `V_liq_L`/`V_gas_L` (or a
`phases={"liquid":..., "gas":...}` dict) become required. The result
gains a `partial_pressures_atm` field (atm, not mol/L — kept separate
from `species_mol_L` to avoid mislabeling a pressure as a
concentration).\
"""),

    code("s5-code", """\
co2_henry = HenryEquilibrium(H_ref=3.4e-4, dlnH=2400.0, gas_species="CO2", liquid_species="CO2",
                              label="co2_henry")

engine_folded = NRChemicalEquilibriumEngine.from_reactions(
    [water, co2_first, co2_second, co2_henry], T_K=298.15,
)
print("gas_liquid_species(): ", engine_folded.gas_liquid_species())

V_liq, V_gas, T_K = 1.0, 0.2, 298.15
n_total_C = 0.05   # mol, across BOTH phases
out_folded = engine_folded.solve(
    totals={"CO2": n_total_C / V_liq}, strong_ions={},
    V_liq_L=V_liq, V_gas_L=V_gas,
)

C_liq_total = (
    out_folded.species_mol_L["CO2"] + out_folded.species_mol_L["HCO3-"]
    + out_folded.species_mol_L["CO3--"]
)
n_liq = C_liq_total * V_liq
n_gas = out_folded.partial_pressures_atm["CO2"] * V_gas / (0.0820574 * T_K)

print(f"pH                        = {out_folded.pH:.4f}")
print(f"partial_pressures_atm     = {out_folded.partial_pressures_atm}")
print(f"n_liq (dissolved C, mol)  = {n_liq:.6f}")
print(f"n_gas (gas-phase C, mol)  = {n_gas:.6f}")
print(f"n_liq + n_gas             = {n_liq + n_gas:.6f}  (must equal n_total_C = {n_total_C})")
"""),

    md("s5-gotcha-md", """\
### 5.1  Gotcha: `retain_jacobian=True` refuses a folded gas tableau

The white-box Jacobian machinery (`residual()`, `jacobian_dg_dz()`,
`jacobian_dz_dy()`) keys its internal caches by bare `species_id`, which
collides whenever a gas secondary shares an id with its liquid parent
(the common Henry case — both literally `"CO2"`). Constructing with
`retain_jacobian=True` against a tableau with folded gas secondaries
raises `NotImplementedError` rather than silently producing a wrong
Jacobian; the main `solve()` path above is unaffected. See
`10_engine_protocol_hierarchy.ipynb` for the white-box path against a
non-folded (acid-base-only) tableau.\
"""),

    code("s5-gotcha-code", """\
try:
    NRChemicalEquilibriumEngine.from_reactions(
        [water, co2_first, co2_second, co2_henry], T_K=298.15, retain_jacobian=True,
    )
except NotImplementedError as exc:
    print(f"NotImplementedError: {str(exc)[:120]}...")
"""),

    md("summary", """\
## Summary

| Topic | Key takeaway |
|---|---|
| Construction | One flat list; `classify_equilibrium_constraint()` auto-routes acid-base / gas-liquid / solid-liquid |
| `solve()` totals | `totals={master_id: mol/L}` dict + `strong_ions={...}` dict — **not** named kwargs |
| Precipitation | Solid-liquid items auto-detected from the flat list; drives a nested active-set loop; `saturation_indices` / `extra["minerals_xi_mol_L"]` on the result |
| Gas-liquid folding | Fully-parameterized `HenryEquilibrium`/`RaoultEquilibrium` solved *simultaneously*, not diverted; needs `V_liq_L`/`V_gas_L` (or `phases=`); result gains `partial_pressures_atm` |
| `retain_jacobian` | Raises `NotImplementedError` if the tableau has folded gas secondaries (id-collision guard) |
| Result type | `EquilibriumResult`, same as Bisection/PHREEQC; `charge_residual` populated (Bisection always returns `None`) |\
"""),
)

save_path_nr = HERE / "02_nr_engine_basics.ipynb"
with open(save_path_nr, "w", encoding="utf-8") as f:
    json.dump(nr_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path_nr}")


# ═══════════════════════════════════════════════════════════════════════════
#  03  PHREEQCChemicalEquilibriumEngine basics
# ═══════════════════════════════════════════════════════════════════════════

SETUP_PHREEQC = """\
import sys
from pathlib import Path

def _find_repo():
    for p in [Path.cwd(), *Path.cwd().parents]:
        if (p / "pyproject.toml").exists():
            return p
    raise RuntimeError("Run from inside the PyOMES repo")

sys.path.insert(0, str(_find_repo() / "models"))

from PyOMES.chemical_equilibrium.phreeqc_engine import (
    PHREEQCChemicalEquilibriumEngine, phreeqc_to_vlsim,
)

print("Imports OK")
"""

phreeqc_nb = nb(
    md("title", """\
# `PHREEQCChemicalEquilibriumEngine` — Instantiation and Use

Three engines satisfy `ChemicalEquilibriumEngineProtocol`:
`BisectionChemicalEquilibriumEngine`, `NRChemicalEquilibriumEngine`, and
`PHREEQCChemicalEquilibriumEngine`. None is "the" default — pick the one
that matches your chemistry's shape. This notebook covers
`PHREEQCChemicalEquilibriumEngine`: a black-box wrapper around
[`phreeqpython`](https://github.com/Vitens/phreeqpython) (an **optional**
dependency — `pip install PyOMES[phreeqc]`), delegating to PHREEQC's own
thermodynamic database and ion-pair library rather than PyOMES-declared
`EquilibriumReaction` networks.

**Contrast with `BisectionChemicalEquilibriumEngine`/`NRChemicalEquilibriumEngine`:**
this is the most structurally different of the three. It has no
`from_reactions()` classmethod at all — chemistry comes from PHREEQC's
own database via a `component_map` (PyOMES id → PHREEQC element key), not
from declared `EquilibriumReaction`/`HenryEquilibrium`/`KspEquilibrium`
instances. It satisfies only the black-box tier of the protocol hierarchy
(`solve()`, `algebraic_species()`, `reset_cache()`, `reset_counters()`) —
wrap it with `NumericalGradientEquilibriumEngine` for `jacobian_dz_dy()`
rather than expecting it directly, unlike `NRChemicalEquilibriumEngine`
which supports `retain_jacobian=True` natively.

If `phreeqpython` is not installed, constructing this engine raises a
clear `ImportError` with installation instructions rather than failing
obscurely later — this notebook assumes it's already installed (it's
installed in this kernel's environment).

For a deep accuracy comparison against the other two engines, see
[`../../model_api/chemistry/speciation/06_phreeqc_benchmark.ipynb`](../../model_api/chemistry/speciation/06_phreeqc_benchmark.ipynb)
and
[`../../model_api/chemistry/speciation/10_engine_protocol_hierarchy.ipynb`](../../model_api/chemistry/speciation/10_engine_protocol_hierarchy.ipynb).
See also [`01_bisection_engine_basics.ipynb`](01_bisection_engine_basics.ipynb)
and [`02_nr_engine_basics.ipynb`](02_nr_engine_basics.ipynb) for the other
two engines.\
"""),

    code("setup", SETUP_PHREEQC),

    # ── 1. Instantiation ──────────────────────────────────────────────────────

    md("s1-md", """\
## 1  Instantiation: `component_map`, not `from_reactions()`

Construction takes a priming composition (`{vlsim_id: mmol/L}`) plus a
`component_map` translating each PyOMES component id to a PHREEQC element
key (oxidation-state notation where the element has more than one, e.g.
`"C(4)"` for carbon(IV) as CO2/carbonate, `"N(-3)"` for reduced nitrogen
as NH3/NH4+). The priming solve both validates the mapping immediately
(a bad key raises during PHREEQC's own parse, not silently) and caches a
conservative superset of species for `algebraic_species()`.\
"""),

    code("s1-code", """\
engine = PHREEQCChemicalEquilibriumEngine(
    {"CO2": 1.0, "Ca": 0.5},                       # priming composition, mmol/L
    component_map={"CO2": "C(4)", "Ca": "Ca"},
)

print("Engine constructed:  ", type(engine).__name__)
print("component_map:       ", engine.component_map)
print("n_solve_calls:       ", engine.n_solve_calls)
"""),

    # ── 2. solve() and the totals= dict ────────────────────────────────────────

    md("s2-md", """\
## 2  `solve()` and the `totals=` dict

`solve()` uses the **same `totals=` dict convention as
`NRChemicalEquilibriumEngine`** — `{vlsim_component_id: mol_L}` — not
`BisectionChemicalEquilibriumEngine`'s named kwargs. There is no
`strong_ions=` split: every component (including strong ions like Ca) goes
through the same `component_map`/`totals` path, since PHREEQC's own
database handles the full ion-pairing network.\
"""),

    code("s2-code", """\
result = engine.solve(totals={"CO2": 0.010, "Ca": 0.002})

print(f"pH               = {result.pH:.4f}")
print(f"logH             = {result.logH:.4f}")
print(f"ionic_strength   = {result.ionic_strength:.4e} mol/L")
print(f"n_solve_calls    = {engine.n_solve_calls}")
"""),

    # ── 3. Species name translation ────────────────────────────────────────────

    md("s3-md", """\
## 3  Species name translation: `phreeqc_to_vlsim`

PHREEQC's own species names use charge-number notation (`Ca+2`, `CO3-2`);
PyOMES convention repeats the sign (`Ca++`, `CO3--`). The default
`species_map=phreeqc_to_vlsim` translates every key in `sol.species`
before it lands in `result.species_mol_L` — this is what lets a calcium
carbonate system (declared only via `component_map`, no hand-built
`EquilibriumReaction` network) come back with names consistent with the
rest of PyOMES.\
"""),

    code("s3-code", """\
print("Ca+2  ->", phreeqc_to_vlsim("Ca+2"))
print("CO3-2 ->", phreeqc_to_vlsim("CO3-2"))
print("HCO3- ->", phreeqc_to_vlsim("HCO3-"), " (pass-through, already singly-charged)")

print()
print("From the solve() above:")
for sp in sorted(result.species_mol_L):
    if sp in ("Ca++", "CO3--", "CaCO3", "CaHCO3+", "HCO3-", "CO2"):
        print(f"  {sp:<10} {result.species_mol_L[sp]:.6e} mol/L")
"""),

    # ── 4. Warmstart cache and a real gotcha ───────────────────────────────────

    md("s4-md", """\
## 4  Gotcha: warmstart drifts from a fresh solve

`use_warmstart=True` (the default) reuses the PHREEQC `Solution` object
across calls via `sol.change(...)`, which issues a PHREEQC `REACTION`
block — an **incremental** addition relative to whatever composition the
solution already holds — rather than resetting to an absolute total. A
fresh engine (`use_warmstart=False`) creates a brand-new `Solution` on
every call instead, always starting from the declared total.

For a single one-shot `solve()` call these agree closely. Across a
*sequence* of calls with changing totals, they can drift apart by more
than the "should agree" tolerance used in this repo's own test suite
(`abs=0.02` pH) — this is a known, currently-undocumented-elsewhere
limitation of the warmstart path, not a numerical-precision artifact.
Reach for `use_warmstart=False` (or `reset_cache()` between calls) if you
need exact agreement with a fresh solve at every step.\
"""),

    code("s4-code", """\
warm = PHREEQCChemicalEquilibriumEngine({"CO2": 1.0}, component_map={"CO2": "C(4)"}, use_warmstart=True)
cold = PHREEQCChemicalEquilibriumEngine({"CO2": 1.0}, component_map={"CO2": "C(4)"}, use_warmstart=False)

print(f"{'CT_CO2':>8} {'warm pH':>10} {'cold pH':>10} {'|diff|':>8}")
for ct in (0.001, 0.005, 0.010, 0.050):
    out_w = warm.solve(totals={"CO2": ct})
    out_c = cold.solve(totals={"CO2": ct})
    print(f"{ct:>8.3f} {out_w.pH:>10.4f} {out_c.pH:>10.4f} {abs(out_w.pH - out_c.pH):>8.4f}")
"""),

    # ── 5. algebraic_species(), reset_cache(), reset_counters() ────────────────

    md("s5-md", """\
## 5  `algebraic_species()`, `reset_cache()`, `reset_counters()`

- `algebraic_species()` — a **conservative superset** discovered at the
  priming solve, translated through `species_map`. Species with
  negligible concentration at the priming point are still included;
  they'll just have near-zero values in practice. Unlike the Bisection/NR
  engines (whose algebraic species come from the declared reaction
  graph), this set comes from whatever PHREEQC's database considers
  possible for the priming composition.
- `reset_cache()` — discards the warmstart `Solution`; the next
  `solve()` builds a fresh one (a one-shot escape from Section 4's
  drift, at the cost of losing the warmstart speed benefit for that call).
- `reset_counters()` — resets `n_solve_calls` to zero.\
"""),

    code("s5-code", """\
alg = engine.algebraic_species()
print("algebraic_species() sample:", sorted(alg)[:10])
print("'Ca++' in algebraic_species():", "Ca++" in alg)
print("'CO3--' in algebraic_species():", "CO3--" in alg)

print("\\nn_solve_calls before reset:", engine.n_solve_calls)
engine.reset_counters()
print("n_solve_calls after reset_counters():", engine.n_solve_calls)

engine.reset_cache()
print("_sol after reset_cache():", engine._sol)
"""),

    md("summary", """\
## Summary

| Topic | Key takeaway |
|---|---|
| Construction | Direct `__init__(components, component_map=...)` — **no** `from_reactions()`; chemistry comes from PHREEQC's own database, not declared `EquilibriumReaction` networks |
| Dependency | Optional (`phreeqpython`); raises a clear `ImportError` at construction if missing |
| `solve()` totals | `totals={component_id: mol_L}` dict — same shape as `NRChemicalEquilibriumEngine`, no `strong_ions=` split |
| Name translation | `species_map=phreeqc_to_vlsim` (default) converts `Ca+2`→`Ca++`, `CO3-2`→`CO3--`, etc.; pass `None` for raw PHREEQC names |
| Protocol tier | Black-box only — wrap with `NumericalGradientEquilibriumEngine` for `jacobian_dz_dy()`, unlike NR's native `retain_jacobian=True` |
| Warmstart gotcha | `use_warmstart=True` reacts incrementally (PHREEQC `REACTION` block) rather than resetting to an absolute total — drifts from a fresh solve across a sequence of calls; use `use_warmstart=False` or `reset_cache()` for exact agreement |
| `algebraic_species()` | Conservative superset from the priming solve, not derived from a declared reaction graph |\
"""),
)

save_path_phreeqc = HERE / "03_phreeqc_engine_basics.ipynb"
with open(save_path_phreeqc, "w", encoding="utf-8") as f:
    json.dump(phreeqc_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path_phreeqc}")
