# Reactions Subfolder Structure — Design Note

> Design discussion, 2026-09-24. No branch, no checklist, no code yet. Written
> from a planning conversation that audited `PyOMES/reactions/` (13 flat `.py`
> files, 3,840 lines) and asked whether grouping files by the kind of reaction
> they implement would make the package easier for a new reader to navigate.
> Direction approved by the repo owner the same day: **`kinetic/` and
> `equilibrium/` subfolders, `blackbox.py` stays flat**. Modelled on the
> `chemical_equilibrium/engines/` reorganisation
> ([`../shipped/CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER.md`](../shipped/CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER.md)).

## Re-verification of the kickoff brief (2026-09-24, on `main` at `61d481a`)

The brief that started this note gave a file list, an intra-package import
graph and rough consumer counts, and asked for all of it to be re-checked.

**Method.** Line counts with `wc -l`. Imports read from every file in
`reactions/`. Consumers counted with a script that walks the whole working tree
(every file type, including gitignored `notes/` and `scratch/` and hidden
`.github/` and `.claude/`; excluding `.git/`, `__pycache__/`, `.pytest_cache/`
and `PyOMES.egg-info/`), matching `reactions.<module>`, `reactions/<module>` and
`reactions\<module>`. Real import statements were then separated from prose
with a second pattern (`from|import` + `PyOMES.` or relative dots +
`reactions.<module>`). Bare `<module>.py` mentions were counted separately
because they are ambiguous (`builder.py` also names
`templates/stirred_tank/builder.py`).

### Confirmed as the brief says

- 13 files, 3,840 lines, every per-file count exact: `__init__` 102,
  `stoichiometry` 334, `environment` 86, `protocols` 60, `_shared` 175,
  `kinetic` 180, `rate_laws` 525, `builder` 368, `equilibrium` 352,
  `phase_equilibria` 477, `blackbox` 228, `reaction_system` 568, `plots` 385.
- `builder.py` imports `kinetic` and `rate_laws` (plus `stoichiometry`,
  `environment`, `chemistry.species`, `chemistry.common_species`).
- `phase_equilibria.py` imports `equilibrium.vant_hoff_log_K` and
  `chemistry.partition._kH_mol_L_atm_from_ref`.
- `reactions` <-> `chemical_equilibrium` is a function-level package cycle:
  `reaction_system.py:261,273` import the engines;
  `engines/bisection/engine.py:203`, `engines/nr/engine.py:239` and
  `engines/nr/tableau.py:401` import `reactions.equilibrium`.
- Neither guard test looks at the layout of `reactions/`.
  `test_import_graph_acyclic.py` checks module-level imports only.
  `test_package_layering.py` checks `chemistry/` only; its two
  `reactions.equilibrium` strings (lines 71, 73) are synthetic detector inputs,
  not imports, and stay valid.

### Corrected

1. **`plots.py` is not lazy-imported both ways.** `equilibrium.py:61` and
   `reaction_system.py:54` import it at module level. `plots.py` imports
   `equilibrium` and `reaction_system` only under `TYPE_CHECKING` (lines 16-18).
   Only matplotlib (and numpy in `plot_speciation`) is imported lazily.
2. **`_shared.py` has three users, not two:** `reaction_system.py:53` imports
   `fmt_stoichiometry_string`, besides `kinetic.py` and `equilibrium.py`.
3. **`rate_laws.py` imports nothing from `reactions/`.** `builder.py` is its only
   in-package user; `kinetic.py` does not import it (a rate law is any callable).
4. **`phase_equilibria.py` has more dependencies than listed:** also
   `..units` (twice), `..chemistry.species`, `..chemistry.common_species` and
   `.stoichiometry`.
5. **`plots.py` serves one group, not two.** `plot_speciation` plots acid-base
   α fractions from a `ReactionSystem`'s equilibrium ladder. Both functions are
   equilibrium plots; `ReactionSystem` is just one caller.
6. **Consumer counts** differ, most of all for `reaction_system` (Sphinx-style
   docstring references). Files outside `reactions/`:

   | Module | Brief | All | Live (excl. `shipped/`, `ideas/`) | With a real import |
   |---|---|---|---|---|
   | `equilibrium` | ~51 | 54 | 46 | 40 (1 in `scratch/`) |
   | `stoichiometry` | ~42 | 46 | 39 | 37 (1 in `scratch/`) |
   | `reaction_system` | ~16 | 33 | 22 | 8 |
   | `phase_equilibria` | ~10 | 13 | 11 | 3 |
   | `kinetic` | ~9 | 7 | 6 | 3 |
   | `builder` | ~1 | 8 | 3 | 1 |
   | `environment` | ~4 | 8 | 5 | 4 |
   | `rate_laws` | ~5 | 6 | 4 | 3 |
   | `protocols` | ~2 | 6 | 2 | 0 |
   | `blackbox` | ~5 | 4 | 1 | 0 |
   | `plots` | ~2 | 3 | 2 | 0 |
   | `_shared` | — | 3 | 0 | 0 |

   No hits in `notes/`, `.github/` or `.claude/`. `pyproject.toml` and
   `setup.py` name no module paths (`find_packages(include=["PyOMES",
   "PyOMES.*"])` picks up new subpackages automatically).

### Found along the way

- **Almost every deep import pulls a name the package root already exports.**
  Across the 40 files importing `reactions.equilibrium`, the names are
  `EquilibriumReaction` (45 uses), `classify_equilibrium_constraint` (15),
  `EquilibriumConstraint` (9) and `vant_hoff_log_K` (8). Only the last three
  are not exported, and they appear in 8 files: the three engine files above,
  `test_chemistry_database.py`, `test_equilibrium_classification.py`,
  `test_equilibrium_constraint.py`, `test_equilibrium_constraint_dual_role.py`
  and `test_partition_model.py`. Every import of `kinetic`, `rate_laws`,
  `builder` and `phase_equilibria` pulls exported names only.
- **Doc errors.** `docs/architecture.md:372-382` lists the package without
  `rate_laws.py` or `plots.py`. `README.md:106` and `PyOMES/README.md:15` say
  `ReactionBuilder` builds `KineticReaction` / `EquilibriumReaction` objects;
  it builds only `KineticReaction` (`builder.py` imports nothing from
  `equilibrium`).
- **Pickled checkpoints.** `Simulation.save_checkpoint` with `"data"` fidelity
  pickles the full object graph (`core/simulation.py:374-399`), so moved classes
  change their pickled module path. The repo holds no `.pkl`/`.pickle` files.

## Motivation

A new reader listing `reactions/` sees 13 files side by side. Three are named
after reaction types (`kinetic`, `equilibrium`, `blackbox`), but others say
nothing about which type they serve: `builder.py` builds only kinetic
reactions, `rate_laws.py` holds only growth kinetics for them, `plots.py` holds
only equilibrium plots, and `phase_equilibria.py` is a second family of
equilibrium constraints. `protocols.py` holds one protocol (`ReactionModel`)
while the other two live elsewhere (`EquilibriumConstraint` in
`equilibrium.py`, `GrowthKinetics` in `rate_laws.py`).

## Groups and seams

| Group | Files | Imports from `reactions/` |
|---|---|---|
| Kinetic | `kinetic`, `rate_laws`, `builder` | `stoichiometry`, `environment`, `_shared`, each other |
| Equilibrium | `equilibrium`, `phase_equilibria`, `plots` | `stoichiometry`, `_shared`, each other (`plots` -> `reaction_system` under `TYPE_CHECKING` only) |
| Black-box | `blackbox` | `environment` |
| Shared | `stoichiometry`, `environment`, `protocols`, `_shared`, `reaction_system`, `__init__` | — |

**No group imports another group.** Only the shared files are imported across
groups, and only `reaction_system.py` and `__init__.py` import every group. The
split follows an existing seam, as the engines split did, and creates no import
cycle inside `reactions/` or between packages.

## Target layout

```
reactions/
    __init__.py            public exports (names unchanged)
    stoichiometry.py
    environment.py
    protocols.py           ReactionModel (the cross-type protocol)
    _shared.py
    reaction_system.py
    blackbox.py
    kinetic/
        __init__.py        docstring only
        reaction.py        KineticReaction                     (was kinetic.py)
        rate_laws.py       GrowthKinetics, Monod, ...
        builder.py         ReactionBuilder
    equilibrium/
        __init__.py        docstring only
        constraint.py      EquilibriumConstraint, vant_hoff_log_K,
                           classify_equilibrium_constraint     (split from equilibrium.py)
        reaction.py        EquilibriumReaction                 (split from equilibrium.py)
        interphase.py      HenryEquilibrium, RaoultEquilibrium,
                           KspEquilibrium                      (was phase_equilibria.py)
        plots.py           plot_vant_hoff, plot_speciation
```

Top level: 7 files and 2 folders. Kinetic group about 1,070 lines, equilibrium
group about 1,210, shared files about 1,550.

**Splitting `equilibrium.py`.** Lines 67-152 (the protocol, `vant_hoff_log_K`,
`classify_equilibrium_constraint`, plus `_LOG10_E` and the `R` import) go to
`constraint.py`; `EquilibriumReaction` (155-352) goes to `reaction.py`. The
halves do not depend on each other: `classify_equilibrium_constraint` never
refers to `EquilibriumReaction` (it classifies by phase tags) and
`EquilibriumReaction` never calls `vant_hoff_log_K`. `constraint.py` needs only
`_shared.phases_from_entries` and `units`; `reaction.py` keeps the other
`_shared` helpers, `_parse_stoichiometry`, `Species` and `plots`. The engines,
`ReactionSystem` and `interphase.py` depend on `constraint.py` only.

## Decisions (2026-09-24)

1. **Layout B above:** `kinetic/` and `equilibrium/` subfolders. Chosen over
   keeping the package flat (with a file map in the `__init__.py` docstring,
   about 3 files touched), and over also giving `blackbox` a folder.
2. **`blackbox.py` stays flat.** One file, nothing supporting it; same reasoning
   as `engines/phreeqc.py`. It can become a folder if it grows.
3. **Folder names reuse the old module names** (`kinetic`, `equilibrium`),
   matching `KineticReaction` / `EquilibriumReaction`. Consequence: the old
   dotted paths now name packages, so a stale
   `from PyOMES.reactions.equilibrium import EquilibriumReaction` fails with
   `ImportError: cannot import name`, not `ModuleNotFoundError`, and searches
   for stale paths must tell `reactions.equilibrium` from
   `reactions.equilibrium.constraint` (match the old path followed by a
   non-`.` character).
4. **`equilibrium.py` is split into `constraint.py` and `reaction.py`,** as its
   own checkpoint straight after the move, so the move stays a pure relocation.
   Cost accepted: the files that import the constraint names (8 listed above,
   plus `reaction_system.py` and `interphase.py`) are edited twice.
5. **`phase_equilibria.py` is renamed `interphase.py`.** Each of its classes
   relates the same species across two phases, which is also why they double as
   `PartitionModel`s. Rejected: `cross_phase.py` (a cross-phase
   `EquilibriumReaction` is also cross-phase, and `ReactionSystem`'s
   `_cross_phase_equilibria` bucket holds both), `partition.py` (confusable with
   `chemistry/partition.py`), keeping the name (stutters under `equilibrium/`).
   The module docstring says it holds the named physical-law constraints that
   double as `PartitionModel`s. No extra cost: every file naming the old module
   is rewritten by the move anyway.
6. **No shims.** Subpackage `__init__.py` files are docstring-only, as in
   `chemical_equilibrium/engines/`. Old paths simply stop working.
7. **Package root exports unchanged.** `EquilibriumConstraint`,
   `vant_hoff_log_K` and `classify_equilibrium_constraint` stay deep-only, now
   at `PyOMES.reactions.equilibrium.constraint`.
8. **Consumers outside `PyOMES/` use the short form.** Tests, notebooks,
   generator scripts and `models/` switch every deep import of an exported name
   to `from PyOMES.reactions import ...`, so they no longer depend on file
   layout. Code inside `PyOMES/` keeps deep imports (importing the package root
   from inside the package runs all of `reactions/__init__.py`, a larger cycle
   risk). Same rule as the partition relocation's Decision 2. *Recommended
   2026-09-24; the owner asked for a recommendation and may overturn it.*
9. **Absolute imports in the moved files.** Anything a moved file imports from
   outside its own folder becomes absolute (`from PyOMES.reactions.stoichiometry
   import ...`, `from PyOMES.chemistry.species import ...`,
   `from PyOMES.units import ...`); imports between files in the same folder
   stay relative (`from .rate_laws import Monod`). Staying files that import
   into the folders keep their relative style (`from .equilibrium.constraint
   import ...` in `reaction_system.py`). The three engine lazy imports, which
   must be rewritten anyway, become absolute
   (`from PyOMES.reactions.equilibrium.constraint import ...`). The other deep
   relative imports in `chemical_equilibrium/engines/` stay as they are, and
   `OPEN_WORK.md`'s "Relative imports three or more dots deep" entry is updated
   for the three lines this changes.
10. **`docs/architecture.md`, `README.md` and `PyOMES/README.md` are corrected
    in this phase**: the architecture tree gains the missing files and the new
    layout, and both READMEs stop saying `ReactionBuilder` builds
    `EquilibriumReaction` objects. Those lines are edited anyway.
11. **Defaults (not put to a vote):** the 37 files that deep-import
    `stoichiometry` are left alone (that path does not change); the one
    `scratch/` notebook (gitignored) is not updated; `plot_vant_hoff` moves
    unchanged although it has no caller (already in `OPEN_WORK.md`, whose entry
    is repointed); old pickled checkpoints that hold moved classes will not load
    (accepted, as in the partition relocation); the dated comment at
    `test_bsm2_reference.py:351` (`src/reactions/equilibrium.py`) is a record of
    a past re-baseline and is left as written.

## What breaks and what has to change

**Deep paths that stop working** (files with a real import statement):

| Old path | New path | Import files | Of which need a non-exported name |
|---|---|---|---|
| `reactions.equilibrium` | `.equilibrium.reaction` / `.equilibrium.constraint` | 40 | 8 |
| `reactions.kinetic` | `.kinetic.reaction` | 3 (`adm1/base.py`, `adm1/bsm2.py`, `aerobic_fermentation_stoichiometry.ipynb`) | 0 |
| `reactions.rate_laws` | `.kinetic.rate_laws` | 3 (`templates/stirred_tank/{__init__,builder}.py`, `test_rate_laws.py`) + own doctest | 0 |
| `reactions.phase_equilibria` | `.equilibrium.interphase` | 3 (`databases/{anaerobic_digestion,bioprocess_basic}.py`, `templates/stirred_tank/factory.py`) | 0 |
| `reactions.builder` | `.kinetic.builder` | 1 (`test_rate_laws.py:325`) | 0 |
| `reactions.plots` | `.equilibrium.plots` | 0 (prose only: `README.md:170` extras table) | 0 |

The 40 `equilibrium` import files: 16 under `tests/standalone/`, 11 under
`tests/validation/` (tests, 7 notebooks, the generator's string templates), 4
under `docs/tutorials/` (including the ArXiv generator's string templates), 3
engine files, 3 `databases/` files, 2 `models/` files, 1 in `scratch/`.

**Unchanged paths:** `stoichiometry`, `reaction_system`, `environment`,
`protocols`, `_shared`, `blackbox`.

**Scale.** The union of live files naming any moved module is 65: 8 in
`reactions/`, 13 elsewhere in `PyOMES/`, 29 tests, 11 docs (`OPEN_WORK.md`,
`architecture.md`, upcoming notes, tutorials), 2 in `models/`, `README.md`, and
the 1 gitignored `scratch/` file that is left alone.

**Import rewrites inside `reactions/`:**

- `kinetic/reaction.py`: `..chemistry.species`, `.stoichiometry`,
  `.environment`, `._shared` -> absolute.
- `kinetic/builder.py`: `.kinetic` -> `.reaction`; `.rate_laws` stays;
  `.stoichiometry`, `.environment`, `..chemistry.species`,
  `..chemistry.common_species` -> absolute.
- `kinetic/rate_laws.py`: no imports to change; doctest path at line 19.
- `equilibrium/constraint.py` and `reaction.py`: `..chemistry.species`,
  `..units`, `.stoichiometry`, `._shared` -> absolute; `.plots` stays.
- `equilibrium/interphase.py`: `.equilibrium` -> `.constraint`;
  `..chemistry.*`, `..units`, `.stoichiometry` -> absolute.
- `equilibrium/plots.py`: `..units` -> absolute; `TYPE_CHECKING` imports:
  `.equilibrium` -> `.reaction`, `.reaction_system` -> absolute.
- Staying files: `__init__.py` (five import lines) and `reaction_system.py`
  (`.equilibrium`, `.kinetic`, `.plots`). `_shared.py` and `protocols.py` change
  in docstrings only.

**Prose and docstrings to repoint:** docstring cross-references inside
`reactions/`, `chemical_equilibrium/engines/bisection/engine.py`,
`chemistry/partition.py`, `core/gas_liquid_link.py`, `core/transfer_models.py`,
`thermo/equilibrium_constants.py`; live docs `OPEN_WORK.md`,
`docs/architecture.md`, `README.md`, `PyOMES/README.md`, and the upcoming notes
that cite the moved paths (`EXPLICIT_SPECIES_RESOLUTION.md`,
`KINETIC_TRANSFER_GENERALIZATION.md`, `PHENOMENA_PROTOCOL.md`,
`RESERVOIR_TYPE.md`, `README.md`). Bare `*.py` mentions in those files need
reading one by one, since some name other files. `docs/dev/implementation/shipped/`
and `docs/dev/ideas/` are left alone.

## Risks

- **Notebooks and lazy imports fail only when run.** The engine imports are
  function-level; notebook cells are not run by pytest. Re-run every changed
  notebook (from the scratch directory, so no stray output folders land in the
  repo — the engines phase hit this).
- **Generator string templates** (`tests/validation/speciation/_generate_notebooks.py`,
  `docs/tutorials/ArXiv_preprint/_generate_notebooks.py`) hold import lines as
  text. Edit each template and its notebooks together so they do not drift.
- **Line endings.** Every file under `reactions/`, the docs, most tests and the
  notebooks are CRLF. Bulk edits must preserve each file's line endings and
  trailing bytes; check for bare LF afterwards.
- **Silent-skip guards.** Check whether any notebook wraps a moved import in
  `try/except ImportError` (the engines phase found two that would have
  silently skipped cells).
- **Pickled checkpoints** holding `KineticReaction`, `EquilibriumReaction`,
  `HenryEquilibrium`/`RaoultEquilibrium`/`KspEquilibrium` or rate-law
  instances will not load (Decision 11).

## Checkpoint shape

Rough; the checklist decides granularity. Full suite before and after each
(baseline 2098 passed on `main`).

1. Re-verify this note's inventory with a fresh search; record baseline.
2. Move the kinetic group into `kinetic/`; rewrite its imports (Decision 9),
   `__init__.py`, `reaction_system.py` and every consumer (Decision 8).
3. Move `equilibrium.py` whole to `equilibrium/reaction.py`,
   `phase_equilibria.py` to `equilibrium/interphase.py`, `plots.py` to
   `equilibrium/plots.py`; rewrite imports and consumers.
4. Split `equilibrium/reaction.py` into `constraint.py` and `reaction.py`;
   repoint the constraint-name importers.
5. Docstrings, prose and live docs (Decision 10), including subpackage
   docstrings that say what each folder holds.
6. Sweep every old path in every form (dotted, slash, relative, bare filename)
   across `.py`, `.ipynb`, `.md`, `.toml`, `.yml`, `.json`; confirm each old
   path fails to import; re-run changed notebooks.

## Open questions

1. **Guard the seam with a test?** A small test in the style of
   `test_package_layering.py` asserting that `kinetic/` and `equilibrium/` never
   import each other, `blackbox.py` imports neither, and the folders import
   `reaction_system` only under `TYPE_CHECKING`. Recommended: yes, as the last
   checkpoint, so a later edit cannot quietly blur the groups.
2. **Deep `stoichiometry` imports next to rewritten lines.** Many notebooks and
   tests have `from PyOMES.reactions.stoichiometry import StoichiometryEntry`
   on the line after the `equilibrium` import that Decision 8 rewrites.
   Recommended: where the two are adjacent, fold both into one short-form line
   (`from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry`);
   leave other `stoichiometry` importers alone.
3. **Sequencing with other notes.** `EXPLICIT_SPECIES_RESOLUTION.md` Phase 2
   edits `phase_equilibria.py` (to be `interphase.py`) and Phase 0 edits
   `stoichiometry.py` (does not move); `PHENOMENA_PROTOCOL.md` renames
   `ReactionSystem`; `NOTEBOOK_GENERATOR_REMOVAL.md` retires both generators.
   Recommended: independent; whichever lands second works on the new paths. If
   the generator removal lands first, the template edits disappear.

## How to start one

Per this folder's convention ([README.md](README.md#how-to-start-one)): resolve
the open questions, write `REACTIONS_SUBFOLDER_STRUCTURE_CHECKLIST.md`, and cut
a branch off `main` named `reactions-subfolder-structure`. The change is file
moves, one file split and import rewrites, with no behaviour change; the risk is
a missed call site, which the full suite, a notebook re-run and the old-path
sweep should catch.
