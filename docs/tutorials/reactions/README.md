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
| [reaction_system.py](reaction_system.py) | Aerobic-growth `KineticReaction` (via `ReactionBuilder.aerobic_growth`), acid-base `EquilibriumReaction` (`log_K=-pKa`), cross-phase CO₂ partition. Runnable: prints the `ReactionSystem` bucket inventory. Importable: each reaction is a factory function callable from other demos. |
| [chemistry_database.py](chemistry_database.py) | `ChemistryDatabase` lifecycle: import a stock database (`AD_BASIC`), extend it with a custom species, override its `ThermoFramework`. **Currently broken** — `dataclasses.replace(..., activity_model=...)` no longer matches `ThermoFramework`'s constructor; pre-existing, not caused by this move (see the phase checklist's follow-ups). |
| [partition_model.py](partition_model.py) | `PartitionModel`/`HenryPartition` inspection, temperature dependence, H₂S alpha correction. **Currently broken** — calls a `.beta()` method that doesn't exist on `HenryEquilibrium`, the class that replaced the now-deprecated `HenryPartition`; pre-existing, not caused by this move. |
| [fba/fba_toy.py](fba/fba_toy.py) | Dynamic FBA on a 7-reaction toy network, wired through `BlackBoxReactionModel`. The dFBA coupling pattern is faithful to Mahadevan et al 2002; the network shape is hand-crafted (attribution warning in the file). |
| [fba/fba_ecoli_core.py](fba/fba_ecoli_core.py) | Same FBA wiring against an 18-reaction E. coli subset loaded from JSON. Network is an original hand-crafted pedagogical construction (note in the file documents scope and limitations); covers aerobic growth plus mixed-acid fermentation (PFL + ADH anaerobic route). |
| [fba/ecoli_core.json](fba/ecoli_core.json) | Stoichiometry data for `fba_ecoli_core.py`. |

## Reading order

1. **`reaction_system.py`** — start here. Three reactions declared, each
   demonstrating a different kind. Run it to see the bucket inspector.
2. **`fba/fba_toy.py`** — pick this up if you want to wire an external
   solver (dynamic FBA, FBA, CBM, ...) through `BlackBoxReactionModel`.
   The protocol is generic; FBA is just one example.
3. **`fba/fba_ecoli_core.py`** — adds JSON-loaded networks and shows how
   to write a small loader against the `JSONFBASolver` schema.

## Running

From the repo root:

```bash
python docs/tutorials/reactions/reaction_system.py
python docs/tutorials/reactions/fba/fba_toy.py
python docs/tutorials/reactions/fba/fba_ecoli_core.py
```

## Pattern: chemistry as a reusable artifact

`reaction_system.py` exposes three factory functions:

```python
make_aerobic_growth_on_acetate(...)  # KineticReaction
make_acetate_dissociation(...)       # EquilibriumReaction (single-phase)
make_co2_partition()                 # EquilibriumReaction (cross-phase)
```

[`D2C_workshop/`](../D2C_workshop/)'s `raw_construction.py` imports this
module and calls those factories to populate a `ReactionSystem`:

```python
import reaction_system as chem

rxn_system = ReactionSystem([
    chem.make_aerobic_growth_on_acetate(),
    chem.make_acetate_dissociation(),
    chem.make_co2_partition(),
])
```

The same reactions could be reused across different topologies (batch,
CSTR, microplate, multi-CV networks) without redeclaring stoichiometry,
balance closure, or rate laws.
