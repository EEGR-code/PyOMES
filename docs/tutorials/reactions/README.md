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
| [aerobic_fermentation_stoichiometry.ipynb](aerobic_fermentation_stoichiometry.ipynb) | Derives aerobic-growth-on-glucose stoichiometric coefficients from an elemental (CHNOSP) balance by hand, comparing two published biomass elemental compositions (Roels extended, Upcraft) — the derivation behind what `ReactionBuilder.aerobic_growth`'s `balance=` argument does internally. Section 6 runs the resulting stoichiometry as a dynamic batch simulation via `StirredTankBuilder`. |
| [reaction_system.ipynb](reaction_system.ipynb) | Aerobic-growth `KineticReaction` (via `ReactionBuilder.aerobic_growth`), acid-base `EquilibriumReaction` (`log_K=-pKa`), cross-phase CO₂ partition. Prints the `ReactionSystem` bucket inventory. |
| [chemistry_database.ipynb](chemistry_database.ipynb) | `ChemistryDatabase` lifecycle: import a stock database (`AD_BASIC`), extend it with a custom species, override its `ThermoFramework` (`liquid_activity=DaviesLiquidModel()`). |
| [partition_model.ipynb](partition_model.ipynb) | `PartitionModel`/`HenryEquilibrium` inspection, temperature dependence, H₂S alpha correction (`partition_ratio(alpha=...)`), extending a database with a custom `HenryEquilibrium`. |
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
4. **`aerobic_fermentation_stoichiometry.ipynb`** — a self-contained detour
   into where `ReactionBuilder.aerobic_growth`'s coefficients actually come
   from; read independently of 1-3, not a prerequisite for any of them.

## Running

```bash
jupyter lab docs/tutorials/reactions/
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
