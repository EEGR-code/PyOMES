# Partition Constraint Relocation — Design Note

> **Status: Shipped 2026-09-24** — merged into `main` via `git merge --no-ff`
> as commit `afbc546`, tagged `partition-constraint-relocation-shipped`. Six
> checkpoints: the move of `HenryEquilibrium`/`RaoultEquilibrium`/`KspEquilibrium`
> (with `_resolve_species`) from `chemistry/partition.py` to the new
> `reactions/phase_equilibria.py`, exported from `reactions/__init__.py`; docstring
> paths; notebooks and scripts; docs and live notes; a layering test; and tests for
> `RaoultEquilibrium`'s custom parameters. Open questions 1-2 were settled as a new
> file (`phase_equilibria.py`) and exporting all three classes; question 5 (renaming
> `test_partition_model.py`) stayed no. The note below is the audit as originally
> written; where it disagrees with the checklist's "Re-verification" section, the
> checklist is correct. The main corrections: the package graph does not become a
> DAG (two other cycles remain, now logged in `OPEN_WORK.md`); the old file also had
> an unused `TYPE_CHECKING` import of `thermo`; checkpoints pickle these classes, so
> ones saved earlier will not load; the "31 external files" break down differently
> (26 needed edits); and live notes and docstrings the note did not list carried
> stale paths. Full suite green post-merge: 2098 passed, 0 failed (2086 plus two
> layering tests and ten `RaoultEquilibrium` tests). See
> [`PARTITION_CONSTRAINT_RELOCATION_CHECKLIST.md`](PARTITION_CONSTRAINT_RELOCATION_CHECKLIST.md).

> Written before implementation, from a conversational audit. Surfaced 2026-09-23
> in a conversational review of `chemistry/`
> (the same one that produced
> [`EQUILIBRIUM_SET_RELOCATION.md`](EQUILIBRIUM_SET_RELOCATION.md)). It
> develops the second, still-open cause of the package-level
> `chemistry` <-> `reactions` cycle logged in `OPEN_WORK.md`'s
> "Package-level layering" entry (logged 2026-09-22,
> `chemistry-reactions-kinetics-cleanup` checkpoint 13, decision D7) into
> a concrete plan. That entry named two possible fixes and left the
> choice open; this note picks one and records the evidence.

## Motivation

`chemistry/partition.py` holds two different kinds of thing:

- **Phase-partition protocols and helpers** — `PartitionModel`,
  `MultispeciesPartitionModel`, `MultispeciesVLEPartition`. Consumed by
  `core/` (`gas_liquid_link.py`, `transfer_models.py`) and
  `databases/database.py`.
- **Equilibrium constraints** — `HenryEquilibrium`, `RaoultEquilibrium`,
  `KspEquilibrium`. The module's own docstring (`partition.py:16-22`)
  says each "also satisfies
  `PyOMES.reactions.equilibrium.EquilibriumConstraint`", so the same
  instance can go into `KineticGasLiquidLink`'s `partition_models=` dict
  *and* into the reaction list fed to the speciation engine. To build
  that second interface they import `reactions.stoichiometry.
  StoichiometryEntry`/`_parse_stoichiometry` and
  `reactions.equilibrium.vant_hoff_log_K` at module level
  (`partition.py:43-44`).

Those two imports are the `chemistry -> reactions` edge. `reactions ->
chemistry` is the opposite, legitimate edge (`reactions/equilibrium.py:49`
and `reactions/stoichiometry.py:25-26` import `Species`/`common_species`).
Both directions existing at once is the package-level cycle.
`tests/standalone/test_import_graph_acyclic.py` passes because it checks
only the module-level graph and ignores function-level and
`TYPE_CHECKING` imports.

## Which fix

`OPEN_WORK.md` listed two directions: (a) move `StoichiometryEntry` and
its helpers into `chemistry/`, or (b) move the three `*Equilibrium`
classes into `reactions/` and leave the `PartitionModel` protocol in
`chemistry/`. This note takes **(b)**. Evidence:

- `PartitionModel` is a `runtime_checkable` structural `Protocol`. The
  three classes satisfy it by shape, not inheritance, so nothing forces
  them to share a package with it.
- `core/gas_liquid_link.py:97` and `core/transfer_models.py:57`
  (`TYPE_CHECKING`-only) import only `PartitionModel`. Every concrete
  `HenryEquilibrium(...)` in those files is a docstring example or comment.
  `core/` depends on the shape, not the implementations.
- `reactions/reaction_system.py`, `reactions/equilibrium.py` and the
  engines under `chemical_equilibrium/` mention the three classes only in
  docstrings, comments and error-message strings. They handle them through
  the `EquilibriumConstraint` protocol. No logic there changes.
- (a) was not evaluated in depth. `StoichiometryEntry` appears in roughly
  80 files by name (a count that includes docs, so not like-for-like with
  the ~31 below), and moving it would drag `_parse_stoichiometry` and
  `validate_balance` along. (b) moves only the pieces that are already
  `reactions/`-shaped.

## What moves and what stays

Line numbers are as of 2026-09-23 and will drift.

| Piece | Lines | Needs `reactions/`? | Destination |
|---|---|---|---|
| `PartitionModel` (protocol) | 75-102 | No | stays |
| `MultispeciesPartitionModel` (protocol) | 515-557 | No | stays |
| `MultispeciesVLEPartition` | 566-631 | No (only `units`) | stays |
| `_kH_mol_L_atm_from_ref` | 560-563 | No | stays (shared helper) |
| `HenryEquilibrium` | 105-250 | Yes (`StoichiometryEntry`) | moves |
| `RaoultEquilibrium` | 259-385 | Yes (same) | moves |
| `KspEquilibrium` | 390-510 | Yes (`_parse_stoichiometry`, `vant_hoff_log_K`) | moves |
| `_resolve_species` | 54-71 | No, but only Henry/Raoult call it | moves with them |

After the split, the remaining `chemistry/partition.py` imports nothing
from `reactions/`. With the `reactions -> chemistry` edge unchanged, the
package graph becomes a DAG (based on the imports read during this
review; a package-level test at the end of the phase turns that from a
reading into a checked fact).

The moved code's imports become:
- `StoichiometryEntry`, `_parse_stoichiometry` from `.stoichiometry`
- `vant_hoff_log_K` from `.equilibrium`
- `Species`, `common_species` from `..chemistry` (the direction
  `stoichiometry.py` already uses)
- `_kH_mol_L_atm_from_ref` from `..chemistry.partition` — `HenryEquilibrium`
  and `MultispeciesVLEPartition` both use it, so importing it from its one
  home avoids duplicating it.

## Destination in `reactions/`

Open question 1. Two options:

- **New file** (placeholder name `reactions/phase_equilibrium.py`).
  Recommended: merging into `equilibrium.py` takes it from 352 to roughly
  700 lines, and these are a distinct family (gas-liquid and solid-liquid
  constraints that double as partition models) from `EquilibriumReaction`.
- **Merge into `equilibrium.py`**, next to `EquilibriumConstraint` and
  `EquilibriumReaction`, which the three classes conform to. Fits the repo's
  "type plus siblings in one file" precedent, at the cost of file size.

## Consumer updates

Counted directly on 2026-09-23.

**Real code imports of the concrete classes (production):**
- `databases/anaerobic_digestion.py:33` — module-level
  `from ..chemistry.partition import HenryEquilibrium`; constructs five
  instances at import time.
- `databases/bioprocess_basic.py:31` — same import; two instances.
  (`databases -> reactions` already exists via `database.py`'s import of
  `reactions.reaction_system`, so neither `databases/` change adds a new
  package edge.)
- `templates/stirred_tank/factory.py:47` —
  `from PyOMES.chemistry.partition import HenryEquilibrium, PartitionModel`.
  Split into two lines; `PartitionModel` stays.
- `models/vlmodels/adm1/base.py:1055` — lazy
  `from PyOMES.chemistry import RaoultEquilibrium`.
- `chemistry/__init__.py` — drop the three imports and their `__all__`
  entries; update the module docstring if it names them.
- `reactions/__init__.py` — add the three, or leave them deep-only (open
  question 2).

**Docstring / prose references (repoint, no logic change):**
- `reactions/reaction_system.py` (lines 9-11 and prose),
  `reactions/equilibrium.py` (lines 77-79, 118-120 and prose),
  `chemical_equilibrium/engines/bisection/engine.py:147-149`, prose in
  `chemical_equilibrium/engines/nr/{tableau,solver,engine}.py`,
  `core/gas_liquid_link.py` (a `>>> from PyOMES.chemistry import
  HenryEquilibrium` example at line 68 plus prose), `databases/
  anaerobic_digestion.py` docstring, and `partition.py`'s own module
  docstring (which needs rewriting to describe only what stays).

**Outside `PyOMES/` and `docs/dev/` (31 files by class-name match):**
- 17 test files under `tests/standalone/` and `tests/validation/`
  (`test_walker`, `test_simulation`, `test_partition_model`,
  `test_nr_tableau_gas_liquid`, `test_membrane`, `test_gas_liquid_link`,
  `test_equilibrium_constraint_dual_role`, `test_equilibrium_constraint`,
  `test_cv_advance`, `test_chemistry_database`, `test_raoult_h2o_fold_cp4`,
  `test_precipitation_gas_liquid_cp5`, `test_nr_gas_liquid_cp2`,
  `test_equilibrium_classification`, `test_transfer_models`,
  `test_step_internal_transfer_scope_filter`, and
  `tests/validation/speciation/test_iron_oxidation.py`).
- Notebooks: `tests/validation/speciation/{07,08}_*.ipynb`,
  `docs/tutorials/reactions/partition_model.ipynb`,
  `docs/tutorials/protocols/ChemicalEquilibriumProtocol/{0_README,01,02,03}*.ipynb`,
  `docs/tutorials/ArXiv_preprint/{02,03}_*.ipynb`.
- Scripts: `docs/tutorials/D2C_workshop/raw_construction.py`,
  `docs/tutorials/ArXiv_preprint/_generate_notebooks.py` (which is itself
  scheduled for removal by `NOTEBOOK_GENERATOR_REMOVAL.md`; if that ships
  first, one fewer file here).
- `docs/architecture.md` and `docs/tutorials/reactions/README.md`.
- `pyproject.toml` configures no doctests (`testpaths` only), so stale
  `>>>` examples go out of date without failing tests.

Shipped-phase docs under `docs/dev/implementation/shipped/` that mention
`chemistry.partition` are historical logs and are left alone, per this
repo's convention.

## No re-export shim

Re-exporting the three classes from `chemistry/__init__.py` would put the
`chemistry -> reactions` edge straight back, defeating the point. It is
also a real import-order risk: `chemistry/__init__.py` runs whenever
anything imports `chemistry.species`, including `reactions/stoichiometry.py`
while it is itself mid-import, so a package `__init__` importing from
`reactions/` can hit a partially initialised module depending on import
order. The existing acyclic-graph test would not catch that. This repo also
avoids compatibility shims by convention, so every call site is updated in
place and the old import path simply stops working.

## Checkpoint order

Rough shape; the real checklist decides granularity.

1. Create the new `reactions/` file, remove the three classes and
   `_resolve_species` from `partition.py`, update production imports in
   `PyOMES/`, `models/` and `templates/`.
2. Update the 17 test files. Steps 1 and 2 need to land together for a
   green suite, since there is no shim.
3. Update notebooks, tutorials, the generator script and docs; repoint the
   docstring cross-refs.
4. Add a package-level layering test asserting `chemistry` never imports
   `reactions` (function-level and `TYPE_CHECKING` imports included), so
   the fix cannot silently regress.
5. Close out the `OPEN_WORK.md` "Package-level layering" entry.

Whether notebooks are executed by CI (and so must be re-run, not just
edited) was not checked; verify when writing the checklist.

## Open questions

1. **New file or merge into `equilibrium.py`?** See "Destination in
   `reactions/`". Leaning new file, name to be decided.
2. **Export from `reactions/__init__.py`?** Today it exports
   `EquilibriumReaction` but not `EquilibriumConstraint`, `vant_hoff_log_K`
   or `classify_equilibrium_constraint`, so deep-only is the consistent
   default. These three are user-facing constructors (unlike those
   helpers), which argues for exporting them. Decide when the file is
   named.
3. **Sequencing against `EQUILIBRIUM_SET_RELOCATION.md`.** Independent.
   After both ship, `chemistry/` has no outbound dependency except `units`.
   Either order works; doing `EquilibriumSet` first is smaller and lets
   this phase's layering test assert a cleaner `chemistry/`.
4. **Sequencing against `EXPLICIT_SPECIES_RESOLUTION.md` Phase 2.** That
   phase rewrites `_resolve_species`, which moves here. If Phase 2 lands
   first, it's one function fewer to move; if this lands first, Phase 2's
   edits happen in the new file. No conflict either way.
5. **Rename `test_partition_model.py`?** It exercises the three moved
   classes as much as the protocol. Left as-is unless the checklist finds a
   reason.

## Trigger conditions

No hard trigger; nothing is broken today (the module-level graph is
acyclic). Reasonable pickup points: before adding any further
`reactions/`-shaped type to `chemistry/`, or when a layering test at
package granularity is wanted. It is a phase of its own, larger than
`EQUILIBRIUM_SET_RELOCATION.md` because of the ~31 external files, though
nearly all of that is mechanical.

## How to start one

Per this folder's usual convention
([README.md](../upcoming/README.md#how-to-start-one)): resolve open questions 1-2,
write a checklist file (`PARTITION_CONSTRAINT_RELOCATION_CHECKLIST.md`),
cut a branch off `main` (suggested name: `partition-constraint-relocation`).
Logic change is small (two imports per moved class, one file split);
risk sits in missing a call site, which the full suite plus the new
layering test should catch.
