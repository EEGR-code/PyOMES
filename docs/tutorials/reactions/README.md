# Reactions

Direct usage of the `PyOMES.reactions` API surface — how to declare
`KineticReaction`, `EquilibriumReaction`, and `BlackBoxReactionModel`
instances, and how a `ReactionSystem` pre-buckets them. No builder, no
topology — this is the layer the [`templates/`](../templates/) builder
demos hide, and what [`D2C_workshop/`](../D2C_workshop/) wires into a full
`Simulation`.

## Files

| File | What it shows |
|---|---|
| [reaction_system.ipynb](reaction_system.ipynb) | Aerobic-growth `KineticReaction` (via `ReactionBuilder.aerobic_growth`), acid-base `EquilibriumReaction` (`log_K=-pKa`), cross-phase CO₂ partition. Prints the `ReactionSystem` bucket inventory. |
| [chemistry_database.py](chemistry_database.py) | `ChemistryDatabase` lifecycle: import a stock database (`AD_BASIC`), extend it with a custom species, override its `ThermoFramework`. **Currently broken** — `dataclasses.replace(..., activity_model=...)` no longer matches `ThermoFramework`'s constructor; pre-existing, not caused by this move (see the phase checklist's follow-ups). Stays `.py` until fixed — a notebook needs genuine output, not a crash. |
| [partition_model.py](partition_model.py) | `PartitionModel`/`HenryPartition` inspection, temperature dependence, H₂S alpha correction. **Currently broken** — calls a `.beta()` method that doesn't exist on `HenryEquilibrium`, the class that replaced the now-deprecated `HenryPartition`; pre-existing, not caused by this move. Stays `.py` until fixed. |
| [fba/fba_toy.ipynb](fba/fba_toy.ipynb) | Dynamic FBA on a 7-reaction toy network, wired through `BlackBoxReactionModel`. The dFBA coupling pattern is faithful to Mahadevan et al 2002; the network shape is hand-crafted (attribution warning in the notebook). |
| [fba/fba_ecoli_core.ipynb](fba/fba_ecoli_core.ipynb) | Same FBA wiring against an 18-reaction E. coli subset loaded from JSON. Network is an original hand-crafted pedagogical construction (note in the notebook documents scope and limitations); covers aerobic growth plus mixed-acid fermentation (PFL + ADH anaerobic route). |
| [fba/ecoli_core.json](fba/ecoli_core.json) | Stoichiometry data for `fba_ecoli_core.ipynb`. |

## Reading order

1. **`reaction_system.ipynb`** — start here. Three reactions declared, each
   demonstrating a different kind.
2. **`fba/fba_toy.ipynb`** — pick this up if you want to wire an external
   solver (dynamic FBA, FBA, CBM, ...) through `BlackBoxReactionModel`.
   The protocol is generic; FBA is just one example.
3. **`fba/fba_ecoli_core.ipynb`** — adds JSON-loaded networks and shows how
   to write a small loader against the `JSONFBASolver` schema.

## Running

```bash
jupyter lab docs/tutorials/reactions/
```

The two `.py` files (broken, see above) run from the repo root:

```bash
python docs/tutorials/reactions/chemistry_database.py
python docs/tutorials/reactions/partition_model.py
```

## Pattern: chemistry as a reusable artifact

`reaction_system.ipynb` builds three factory functions:

```python
make_aerobic_growth_on_acetate(...)  # KineticReaction
make_acetate_dissociation(...)       # EquilibriumReaction (single-phase)
make_co2_partition()                 # EquilibriumReaction (cross-phase)
```

A caller in the same process could populate a `ReactionSystem` from them:

```python
rxn_system = ReactionSystem([
    make_aerobic_growth_on_acetate(),
    make_acetate_dissociation(),
    make_co2_partition(),
])
```

[`D2C_workshop/raw_construction.py`](../D2C_workshop/raw_construction.py)
declares this exact same pattern independently (a by-value copy, not an
import) — each tutorial folder stays self-sufficient rather than depending
on a sibling one.

The same reactions could be reused across different topologies (batch,
CSTR, microplate, multi-CV networks) without redeclaring stoichiometry,
balance closure, or rate laws.
