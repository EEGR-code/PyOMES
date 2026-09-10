# Model-API demos

Direct usage of the `PyOMES.core` and `PyOMES.reactions` API surfaces —
no `FermenterBuilder`. Two subtrees, organised by **what you're
defining**:

- [`chemistry/`](chemistry/) — reaction declarations
  (`KineticReaction`, `EquilibriumReaction`, `BlackBoxReactionModel`)
  and how a `ReactionSystem` pre-buckets them.
- [`D2Cworkshop/`](D2Cworkshop/) — manual assembly
  of a `Simulation` from `Phase` + `Link` + `ControlVolume`
  primitives, with chemistry pulled in from the `chemistry/` demos.

The two compose: `D2Cworkshop/basic_layout/raw_construction.py` imports
its reactions from `chemistry/reaction_system.py` rather than
redeclaring them, making the reuse pattern visible in code.

## Demos

### [`chemistry/`](chemistry/) — define reactions

| File | What it shows |
|---|---|
| [reaction_system.py](chemistry/reaction_system.py) | Aerobic-growth `KineticReaction` (via `ReactionBuilder.aerobic_growth`), acid-base `EquilibriumReaction` (`log_K=-pKa`), cross-phase CO₂ partition. Runnable: prints the `ReactionSystem` bucket inventory. Importable: each reaction is a factory function callable from other demos. |
| [fba/fba_toy.py](chemistry/fba/fba_toy.py) | Dynamic FBA on a 7-reaction toy network, wired through `BlackBoxReactionModel`. The dFBA coupling pattern is faithful to Mahadevan et al 2002; the network shape is hand-crafted (attribution warning in the file). |
| [fba/fba_ecoli_core.py](chemistry/fba/fba_ecoli_core.py) | Same FBA wiring against an 18-reaction E. coli subset loaded from JSON. Network is an original hand-crafted pedagogical construction (note in the file documents scope and limitations); covers aerobic growth plus mixed-acid fermentation (PFL + ADH anaerobic route). |
| [fba/ecoli_core.json](chemistry/fba/ecoli_core.json) | Stoichiometry data for `fba_ecoli_core.py`. |

### [`D2Cworkshop/`](D2Cworkshop/) — assemble a Simulation

| File | What it shows |
|---|---|
| [raw_construction.py](D2Cworkshop/basic_layout/raw_construction.py) | Mirrors [`builder/batch_fermenter.py`](../builder/batch_fermenter.py) end-to-end — same batch sparged topology, same aerobic growth on acetic acid, same PI pH control — but builds every object explicitly. Imports chemistry from [reaction_system.py](chemistry/reaction_system.py). |

## Reading order

1. **`chemistry/reaction_system.py`** — start here. Three reactions
   declared, each demonstrating a different kind. Run it to see
   the bucket inspector.
2. **`D2Cworkshop/basic_layout/raw_construction.py`** — see those
   reactions wired into a complete `Simulation`. Side-by-side with
   [`builder/batch_fermenter.py`](../builder/batch_fermenter.py)
   shows what the builder hides.
3. **`chemistry/fba/fba_toy.py`** — pick this up if you want to
   wire an external solver (dynamic FBA, FBA, CBM, ...) through
   `BlackBoxReactionModel`. The protocol is generic; FBA is just
   one example.
4. **`chemistry/fba/fba_ecoli_core.py`** — adds JSON-loaded
   networks and shows how to write a small loader against the
   `JSONFBASolver` schema.

## Running

From the repo root:

```bash
python demos/model_api/chemistry/reaction_system.py
python demos/model_api/D2Cworkshop/basic_layout/raw_construction.py
python demos/model_api/chemistry/fba/fba_toy.py
python demos/model_api/chemistry/fba/fba_ecoli_core.py
```

## Pattern: chemistry as a reusable artifact

`reaction_system.py` exposes three factory functions:

```python
make_aerobic_growth_on_acetate(...)  # KineticReaction
make_acetate_dissociation(...)       # EquilibriumReaction (single-phase)
make_co2_partition()                 # EquilibriumReaction (cross-phase)
```

`raw_construction.py` imports the module and calls those factories
to populate a `ReactionSystem`:

```python
import reaction_system as chem

rxn_system = ReactionSystem([
    chem.make_aerobic_growth_on_acetate(),
    chem.make_acetate_dissociation(),
    chem.make_co2_partition(),
])
```

The same reactions could be reused across different topologies
(batch, CSTR, microplate, multi-CV networks) without redeclaring
stoichiometry, balance closure, or rate laws. Adding a future
`chemistry_database.py` demo (per
[CHEMISTRY_UNIFICATION_PLAN.md](../../docs/dev/implementation/shipped/CHEMISTRY_UNIFICATION_PLAN.md))
extends the same pattern to library-loaded chemistry.
