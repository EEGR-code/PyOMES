# Phase Kickoff Checklist — reactions-subfolder-structure

> Checklist for [`REACTIONS_SUBFOLDER_STRUCTURE.md`](REACTIONS_SUBFOLDER_STRUCTURE.md),
> the source of truth for motivation, groups, target layout and decisions; do not
> restate it here. Where this checklist and the note disagree, this checklist
> wins: the "Re-verification" section below corrects the note from a fresh search.
> See [`README.md`](README.md)'s "Branching and tagging convention". Modelled on
> [`../shipped/PARTITION_CONSTRAINT_RELOCATION_CHECKLIST.md`](../shipped/PARTITION_CONSTRAINT_RELOCATION_CHECKLIST.md)
> and
> [`../shipped/CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER_CHECKLIST.md`](../shipped/CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER_CHECKLIST.md).

**Working rules**

- Checkpoints 1-3 are pure refactors: file moves, one file split and import
  rewrites, with no behaviour change. Checkpoint 4 changes only docstrings,
  comments and prose. Checkpoint 5 changes nothing unless the sweep finds a
  missed site. Checkpoint 6 adds a test.
- The repo owner runs every `git` command (branch, add, commit, push, merge,
  tag). Read-only git (`status`, `log`, `diff`, `show`) is fine. Files are moved
  with a plain filesystem move, not `git mv`. Git commands are handed over one
  per line (PowerShell 5.1: never chained with `&&`); commit messages use two
  `-m` flags (a strapline and one body paragraph) and no attribution lines.
- One checkpoint at a time: run the full suite first, edit, run the sanity
  checks, report what changed and the results as they are, stop, and hand over
  the commit commands. Do not start the next checkpoint until told to.
- Full suite: `python -m pytest -p no:cacheprovider -q`.
- No re-export shims, aliases or compatibility layers. Old import paths stop
  working; every call site is updated in place. Subpackage `__init__.py` files
  are docstring-only.
- Every notebook changed in a checkpoint is re-run in that checkpoint, from the
  scratchpad directory (so no output folders land in the repo). Notebook edits
  touch only `source` lines, never saved outputs.
- Nearly every file involved is CRLF. Bulk edits use a script that preserves
  each file's line endings and exact trailing bytes; check for bare LF after
  every edit. Known oddities, to be preserved as they are:
  `tests/validation/speciation/10_engine_protocol_hierarchy.ipynb` already has 3
  bare LFs among 1,111 CRLFs; `docs/tutorials/ArXiv_preprint/01_predict_ph_simple_liquid.ipynb`
  and `docs/tutorials/reactions/aerobic_fermentation_stoichiometry.ipynb` end
  in `\r\n}` with no final newline.
- Leave `docs/dev/implementation/shipped/` and `docs/dev/ideas/` untouched.
- Docs and docstrings describe current behaviour only: no phase or checkpoint
  labels, no pointers to this checklist or the design note.
- Anything that looks like a bug or dead code beyond scope is logged in
  `OPEN_WORK.md`, not fixed.

**Decisions** (settled 2026-09-24, before starting; 1-11 are the note's, in brief)

1. Layout B: `kinetic/` and `equilibrium/` subfolders under `reactions/`.
2. `blackbox.py` stays flat.
3. Folder names reuse the old module names. `from PyOMES.reactions.equilibrium
   import EquilibriumReaction` and the same for `kinetic` therefore fail with
   `ImportError: cannot import name`; the other four old paths (`rate_laws`,
   `builder`, `phase_equilibria`, `plots`) fail with `ModuleNotFoundError`.
   Searches for stale paths must match the old path followed by a non-`.`
   character.
4. `equilibrium.py` is split into `constraint.py` and `reaction.py` in its own
   checkpoint, straight after the move. The constraint-name importers are
   edited twice (accepted).
5. `phase_equilibria.py` is renamed `equilibrium/interphase.py`.
6. No shims; subpackage `__init__.py` files are docstring-only.
7. Package root exports unchanged. `EquilibriumConstraint`, `vant_hoff_log_K`
   and `classify_equilibrium_constraint` stay deep-only, at
   `PyOMES.reactions.equilibrium.constraint`.
8. Consumers outside `PyOMES/` (tests, notebooks, generator scripts, `models/`)
   switch every deep import of an exported name to
   `from PyOMES.reactions import ...`. Code inside `PyOMES/` keeps deep imports.
   The note records this as a recommendation; the owner confirmed it on
   2026-09-24 over deep imports everywhere and the short form everywhere.
9. Moved files import anything outside their own folder absolutely
   (`from PyOMES.reactions.stoichiometry import ...`,
   `from PyOMES.chemistry.species import ...`, `from PyOMES.units import ...`);
   imports within a folder stay relative. Staying files keep their relative
   style (`from .equilibrium.constraint import ...`). The three engine lazy
   imports become absolute; the other deep relative imports in
   `chemical_equilibrium/engines/` stay, and the `OPEN_WORK.md` "Relative
   imports three or more dots deep" entry is updated (14 lines across 6 files
   become 11 across 6).
10. `docs/architecture.md`, `README.md` and `PyOMES/README.md` are corrected in
    this phase (missing files in the tree; `ReactionBuilder` builds only
    `KineticReaction`).
11. Defaults: the gitignored `scratch/` notebook is not updated; `plot_vant_hoff`
    moves unchanged and its `OPEN_WORK.md` entry is repointed; old pickled
    checkpoints holding moved classes will not load; the dated comment at
    `test_bsm2_reference.py:351` is left as written. The note's default of
    leaving every `stoichiometry` importer alone is narrowed by decision 13.
12. **Seam guard** (open question 1): yes, as the last checkpoint, as a second
    test in `tests/standalone/test_package_layering.py`. The file's import
    resolver is generalised once to return the full dotted target plus a
    "`TYPE_CHECKING` only" flag, and both tests use it; the `chemistry/` test's
    behaviour is unchanged. Rules: `kinetic/` and `equilibrium/` never import
    each other; `blackbox.py` imports neither; neither folder imports
    `reaction_system` except under `TYPE_CHECKING`; neither folder imports the
    `PyOMES.reactions` package root. A synthetic self-check makes an empty result
    meaningful.
13. **`stoichiometry` imports** (open question 2): in every file this phase
    edits anyway, every deep `StoichiometryEntry` import switches to the short
    form, folded into the rewritten line where the two sit in the same block
    and adjacent, on its own `from PyOMES.reactions import StoichiometryEntry`
    line otherwise. Aliases are kept (`StoichiometryEntry as E`). Code inside
    `PyOMES/` keeps its deep imports (decision 8). That is 31 live files, 45 lines
    (see Re-verification). `tests/standalone/test_stoichiometry.py` is left
    alone: it imports the private `_parse_stoichiometry` and tests that module
    directly.
14. **Sequencing** (open question 3): independent of
    `EXPLICIT_SPECIES_RESOLUTION.md`, `PHENOMENA_PROTOCOL.md` and
    `NOTEBOOK_GENERATOR_REMOVAL.md`; whichever lands second works on the new
    paths. No other phase branch exists (only `main`, checked 2026-09-24).
15. The note gets an entry in `upcoming/README.md`'s "Design discussions",
    committed with this checklist, and replaced by a "Recently shipped" entry at
    shipping.

## Pre-flight

- [x] `git status -sb` clean, `main` level with `origin/main` at `f88f372`
      (checked 2026-09-24)
- [x] Full-suite baseline on `main`: **2098 passed**, 0 failed, 166 warnings,
      1m50s (`python -m pytest -p no:cacheprovider -q`)
- [ ] Branch created off `main`: `reactions-subfolder-structure`
- [ ] This checklist and the `upcoming/README.md` entry committed on that branch
      as the first commit

## During

- [x] Plan doc exists in `docs/dev/implementation/upcoming/`
      ([`REACTIONS_SUBFOLDER_STRUCTURE.md`](REACTIONS_SUBFOLDER_STRUCTURE.md))
- [ ] Checkpoints tracked below, one commit each
- [ ] **If work stalls:** add a status banner to the top of the plan doc at once

### Re-verification of the design note (2026-09-24, on `main` at `f88f372`)

**Method.** A script walked the whole working tree, including gitignored
`notes/` and `scratch/`, hidden `.github/` and `.claude/`, and
`PyOMES.egg-info/`, skipping only `.git/`, `__pycache__/` and `.pytest_cache/`.
It read `.py`, `.ipynb`, `.md`, `.toml`, `.yml`/`.yaml`, `.json`, `.txt`,
`.cfg`, `.ini` and `.rst` files and matched `reactions.<module>`,
`reactions/<module>` and `reactions\<module>` (not followed by a word
character), relative `from .<module>` inside `reactions/`, and bare
`<module>.py`. Each hit was classified: in `.py` files by AST and tokenizer
(module-level import, lazy import, `TYPE_CHECKING` import, string or docstring,
string template holding an import line, comment); in notebooks by cell type
(code import, code comment, markdown, saved output). A second AST script listed
the names each file imports from each moved module, including generator string
templates and notebook code cells. Cross-checked with plain `grep -r` (which also
sees ignored and hidden paths) and the Grep tool. Line counts with `wc -l`, line
endings by byte count.

**Confirmed as the note says**

- 13 files, 3,840 lines; every per-file count exact. Every file in `reactions/`
  is CRLF with a trailing `\r\n`.
- Intra-package imports: `builder.py` imports `kinetic`, `rate_laws`,
  `stoichiometry`, `environment`, `chemistry.species`, `chemistry.common_species`
  (lines 37-42); `rate_laws.py` imports nothing from `reactions/`;
  `phase_equilibria.py` imports `.equilibrium.vant_hoff_log_K` (:49),
  `.stoichiometry` (:50), `..chemistry` (:44-46), `..units` twice (:47-48);
  `equilibrium.py:61` and `reaction_system.py:54` import `plots` at module level;
  `plots.py:17-18` imports `equilibrium` and `reaction_system` under
  `TYPE_CHECKING` only; `_shared.py` has three users (`kinetic.py:38`,
  `equilibrium.py:52`, `reaction_system.py:53`); `blackbox.py` imports only
  `environment`. No group imports another group.
- The `reactions` <-> `chemical_equilibrium` cycle is function-level only:
  `reaction_system.py:261,273`; `engines/bisection/engine.py:203`,
  `engines/nr/engine.py:239`, `engines/nr/tableau.py:401`.
- `test_package_layering.py:71,73` are synthetic detector inputs, not imports.
- Per-module file counts outside `reactions/` (all / live / live with a real
  import), with `PyOMES.egg-info/` excluded as in the note: `kinetic` 7/6/3,
  `builder` 8/3/1, `phase_equilibria` 13/11/3, `plots` 3/2/0,
  `reaction_system` 33/22/8, `environment` 8/5/4, `protocols` 6/2/0,
  `_shared` 3/0/0, `blackbox` 4/1/0.
- Exactly 8 files import a non-exported name from `reactions.equilibrium`: the
  three engine files, `test_chemistry_database.py`,
  `test_equilibrium_classification.py`, `test_equilibrium_constraint.py`,
  `test_equilibrium_constraint_dual_role.py`, `test_partition_model.py`. Inside
  `reactions/`, `phase_equilibria.py` (`vant_hoff_log_K`) and
  `reaction_system.py` (`EquilibriumConstraint`,
  `classify_equilibrium_constraint`) do too, and `plots` names
  (`plot_vant_hoff`, `plot_speciation`) are imported only inside `reactions/`.
  Every import of `kinetic`, `rate_laws`, `builder` and `phase_equilibria`
  outside `reactions/` pulls exported names only.
- The doc errors: `docs/architecture.md:372-382` lacks `rate_laws.py` and
  `plots.py`; `README.md:106` and `PyOMES/README.md:15` say `ReactionBuilder`
  builds `EquilibriumReaction` objects (`builder.py` imports nothing from
  `equilibrium`). `PyOMES/README.md:15` also says the rate laws "live in
  `rate_laws.py`", which the move makes stale. `README.md:170` names
  `PyOMES.reactions.plots`.
- No `.pkl`/`.pickle`/`.joblib` files anywhere. No `mock.patch`, `importlib`,
  `sys.modules` or `__import__` string names a `reactions` module. No hits in
  `notes/`, `.github/`, `.claude/`, or any `.toml`/`.yml`/`.json` file. No
  notebook markdown cell, code comment or saved output names a moved path.
- **No silent-skip guard wraps a moved import.** Every `except ImportError` in
  `.py`/`.ipynb` guards `phreeqpython`/the PHREEQC engine, matplotlib, pytest or
  optional recorder back ends.
- Packaging: `setup.py`'s `find_packages(include=["PyOMES", "PyOMES.*"])` will
  find the new subpackages, and the editable install's finder
  (`__editable___pyomes_0_12_5_finder.py`) maps the whole `PyOMES` directory, so
  no reinstall is needed. `PyOMES.egg-info/SOURCES.txt` (gitignored, generated)
  lists the old files and is already stale (no `phase_equilibria.py`); left
  alone.
- The `OPEN_WORK.md` "Relative imports three or more dots deep" entry is
  accurate today: 14 lines across 6 files.

**Discrepancies** (none changes a decision; all are folded into the checkpoints)

1. **`reactions.equilibrium` has 39 files with a real import, not 40.**
   `tests/standalone/` has 15, not 16: the note's second pattern matched the two
   synthetic strings in `test_package_layering.py`. The other groups match the
   note (11 under `tests/validation/`, 4 under `docs/tutorials/`, 3 engine files,
   3 `databases/` files, 2 `models/` files, 1 in `scratch/`). All/live counts
   (54/46) match.
2. **`reactions.rate_laws` has 2 files with a real import, not 3**:
   `templates/stirred_tank/__init__.py:30` and `test_rate_laws.py:27`.
   `templates/stirred_tank/builder.py:342` is a `::` code example inside a
   docstring (checkpoint 4), like `rate_laws.py:19`.
3. **Both halves of the split need `phases_from_entries`.**
   `EquilibriumReaction.phases` (:271) calls it, not only
   `classify_equilibrium_constraint` (:147). `constraint.py` also needs
   `StoichiometryEntry` (the protocol's annotation at :88), `math`, and the
   `typing` names `Literal`, `Optional`, `Protocol`, `Sequence`,
   `runtime_checkable`. The protocol begins at line 66 (the decorator), not 67.
4. **The Bisection engine depends on both halves, not only `constraint.py`.**
   Its lazy import (`engines/bisection/engine.py:203-206`) pulls
   `EquilibriumReaction` as well, which it uses in
   `isinstance(rxn, EquilibriumReaction)` at :233. After the split it needs one
   import from `constraint` and one from `reaction`. The NR engine, the tableau,
   `reaction_system.py` and `interphase.py` depend on `constraint.py` only, as the
   note says.
5. **The module docstring of `equilibrium.py` (:2-40) describes only
   `EquilibriumReaction`**, and `classify_equilibrium_constraint`'s docstring
   (:117, :121) names `EquilibriumReaction` by short `:class:` reference. The
   split needs a new docstring for `constraint.py` and a full-path reference in
   the classifier's docstring.
6. **Docstring and comment paths the note's prose list misses:**
   `engines/nr/engine.py:199`, `engines/nr/tableau.py:15,366`,
   `databases/anaerobic_digestion.py:11`,
   `templates/stirred_tank/builder.py:342` (docstring example), and
   `tests/standalone/test_rate_laws.py:2,128,290` (module docstring path, and
   bare `rate_laws.py`/`builder.py` in two docstrings). The note's
   "inside `reactions/`" covers `builder.py:7`, `equilibrium.py:11,77-79,118-120,300`,
   `kinetic.py:155`, `phase_equilibria.py:9,97,247,363,365,385`,
   `protocols.py:10,21`, `rate_laws.py:6,19`,
   `reaction_system.py:5,7-11,27,74,452` and `_shared.py:6-7`.
7. **The live union is 67 files, not 65**: 8 in `reactions/`, 14 elsewhere in
   `PyOMES/` (including `PyOMES/README.md` and `databases/aqueous.py`), 29 tests,
   12 docs (the 7 dev docs below plus 4 tutorial notebooks and the ArXiv
   generator), 2 in `models/`, `README.md`, and the gitignored `scratch/` file.
   Decision 13 adds no file to that set: every file it touches already imports
   `equilibrium` or `kinetic`.
8. **`stoichiometry` import files are 36, not 37**, with 50 import lines: 3 in
   `PyOMES/databases/` (relative, kept deep by decision 8), the `scratch/` copy,
   `test_stoichiometry.py`, and 31 live files outside `PyOMES/` with 45 lines, all
   importing only `StoichiometryEntry`. Those 31 are decision 13's scope: 30 also
   import `equilibrium` (every live `equilibrium` importer outside `PyOMES/`
   except `test_equilibrium_constraint_dual_role.py` and
   `test_partition_model.py`), and `aerobic_fermentation_stoichiometry.ipynb`
   imports `kinetic`. 40 of the 50 lines sit next to an `equilibrium` import (39
   live, in 27 files).
9. **Some cited line numbers are already stale**, so live docs are re-derived,
   not renamed: `OPEN_WORK.md:141` cites `equilibrium.py:92` (the function is at
   :94); `RESERVOIR_TYPE.md:155` cites `equilibrium.py:115-119`;
   `EXPLICIT_SPECIES_RESOLUTION.md:53-54` cite `equilibrium.py:221,231` and
   `kinetic.py:90,96`; `KINETIC_TRANSFER_GENERALIZATION.md:32` cites
   `kinetic.py:48-80`; `PHENOMENA_PROTOCOL.md:53` cites
   `phase_equilibria.py:189-218`.
10. **Doctests are not run by the suite** (`pyproject.toml` sets only
    `testpaths`; no `--doctest-modules`). The `>>>` examples touched by the phase
    (`rate_laws.py:19`, and the module examples in the moved files) are run by
    hand at checkpoint 4.
11. **Live docs to repoint** (the note lists files, not lines):
    `OPEN_WORK.md:136,141,150,171,179,208,461,468,470,594,609,669,676,696,720-724`
    and the "three or more dots" entry at :748-757;
    `EXPLICIT_SPECIES_RESOLUTION.md:53-54,61,63,66,70-71,73,75,183`;
    `KINETIC_TRANSFER_GENERALIZATION.md:32,210,213`;
    `PHENOMENA_PROTOCOL.md:53,98,252,254,259`; `RESERVOIR_TYPE.md:155`;
    `upcoming/README.md:27,37` (`:153` is a shipped-phase entry and stays);
    `docs/architecture.md:372-382`; `README.md:106,170`; `PyOMES/README.md:15`.
    Each needs reading: some describe past work (for example `OPEN_WORK.md:461,
    470` narrate a shipped move) and may be history, not a live path.
12. Consumer import styles inside `PyOMES/` (kept, decision 8): absolute in
    `templates/stirred_tank/__init__.py:30` and `factory.py:48`; relative
    `..reactions.<module>` in the three `databases/` files (`anaerobic_digestion.py:33-34`,
    `aqueous.py:21`, `bioprocess_basic.py:31-32`).

### Checkpoints

- [x] 1. **Kinetic group into `kinetic/`.** Move `kinetic.py` to
      `kinetic/reaction.py`, `rate_laws.py` to `kinetic/rate_laws.py`,
      `builder.py` to `kinetic/builder.py`; add a docstring-only
      `kinetic/__init__.py`. Imports (decision 9): `reaction.py` and `builder.py`
      make their `..chemistry.*`, `.stoichiometry`, `.environment`, `._shared`
      imports absolute; `builder.py`'s `.kinetic` becomes `.reaction` and
      `.rate_laws` stays. Staying files: `reactions/__init__.py` (three import
      lines) and `reaction_system.py:50` (`.kinetic.reaction`). Consumers:
      `templates/stirred_tank/__init__.py:30` (deep,
      `PyOMES.reactions.kinetic.rate_laws`); outside `PyOMES/`, short form
      (decision 8): `models/vlmodels/adm1/{base,bsm2}.py` (kinetic line only; their
      `equilibrium` and `stoichiometry` lines wait for checkpoint 2),
      `test_rate_laws.py:27,325`, and `aerobic_fermentation_stoichiometry.ipynb`
      (with its `StoichiometryEntry` import, decision 13). Docstrings wait for
      checkpoint 4.
      Sanity: `kinetic.py`, `rate_laws.py`, `builder.py` gone from the top level;
      AST of each moved file equal to the original once import statements are
      removed; `from PyOMES.reactions.kinetic import KineticReaction` raises
      `ImportError` and `import PyOMES.reactions.rate_laws` /
      `PyOMES.reactions.builder` raise `ModuleNotFoundError`; fresh-interpreter
      import of `PyOMES`, `PyOMES.reactions`, `PyOMES.reactions.kinetic.builder`,
      `PyOMES.templates` and `models.vlmodels.adm1.bsm2` in several orders;
      the aerobic notebook re-run from the scratchpad; no bare LF in changed files;
      full suite **2098 passed**.
      _Notes: done 2026-09-24. Suite before: **2098 passed** (2m07s); after:
      **2098 passed**, 0 failed, 166 warnings, 2m07s. The three files were moved
      with a plain filesystem move; `kinetic/rate_laws.py` is byte-identical to the
      old `rate_laws.py`, and `kinetic/reaction.py` and `kinetic/builder.py` differ
      from the originals only in their import lines (4 and 5 lines; the ASTs with
      imports removed are identical). `kinetic/__init__.py` is docstring-only
      (checkpoint 4 reviews its wording). Edited lines elsewhere, all imports:
      `reactions/__init__.py` (3), `reaction_system.py:50`,
      `templates/stirred_tank/__init__.py:30` (deep), `adm1/base.py:43` and
      `adm1/bsm2.py:46` (short form; their `equilibrium` and `stoichiometry` lines
      wait for checkpoint 2), `test_rate_laws.py:27,325` (short form), and cell 17
      of `aerobic_fermentation_stoichiometry.ipynb`, where the adjacent
      `stoichiometry` and `kinetic` imports became one
      `from PyOMES.reactions import KineticReaction, StoichiometryEntry` line
      (decision 13); only that cell's `source` changed. Edits were scripted with
      exact-count byte replacements: no bare LF in any changed file, trailing
      bytes unchanged. Old paths: `from PyOMES.reactions.kinetic import
      KineticReaction` raises `ImportError`; `PyOMES.reactions.rate_laws` and
      `PyOMES.reactions.builder` raise `ModuleNotFoundError`. Nine
      fresh-interpreter import orders pass, each asserting that the root,
      subpackage and `templates.stirred_tank` names are the same objects.
      **Deviation (notebook re-run):** the tool sandbox caps every process at
      about 15 % of one CPU core (a Windows job object; a bare busy loop gets the
      same share), and the full run had not finished after about 30 minutes, so
      the notebook was run without its simulation: code cells before cell 18,
      then one `build_sim(...)` call for the first composition, which constructs
      the `KineticReaction`, the `ReactionSystem` and the stirred-tank
      `Simulation` without calling `sim.run` (10 cells plus the call, 2.4 s, from
      the scratchpad; no files left in the repo). Cells 18-20 (the 60 h
      simulation and its two report/plot cells) were not run. CI does not
      execute notebooks, so it does not cover this either.
      **Found after committing:** with `core.autocrlf=true`, git stored the moved
      files under their new paths as LF (the old paths were stored as CRLF), so
      `git diff`/`git blame` show `kinetic/reaction.py` and `kinetic/builder.py` as
      whole-file rewrites; `git blame -w` and `git log -M --follow` see through it.
      Working-tree files are unchanged (CRLF). Left as is: re-storing them as CRLF
      would be a second whole-file rewrite. From checkpoint 2 on, each moved file is
      staged so git stores it with the same line endings as its old path
      (`git -c core.autocrlf=false add` where the old path was stored as CRLF).
      In checkpoint 2 that is `equilibrium/reaction.py` only; `phase_equilibria.py`
      and `plots.py` were already stored as LF._
- [x] 2. **Equilibrium group into `equilibrium/`, no split yet.** Move
      `equilibrium.py` whole to `equilibrium/reaction.py`, `phase_equilibria.py` to
      `equilibrium/interphase.py`, `plots.py` to `equilibrium/plots.py`; add a
      docstring-only `equilibrium/__init__.py`. Imports (decision 9):
      `..chemistry.*`, `..units`, `.stoichiometry`, `._shared` become absolute;
      `reaction.py`'s `.plots` stays; `interphase.py`'s `.equilibrium` becomes
      `.reaction`; `plots.py`'s `TYPE_CHECKING` imports become `.reaction` and
      absolute `PyOMES.reactions.reaction_system`. Staying files:
      `reactions/__init__.py` (two import lines), `reaction_system.py:49,54`
      (`.equilibrium.reaction`, `.equilibrium.plots`). Package consumers (deep):
      the three engine lazy imports become absolute
      `from PyOMES.reactions.equilibrium.reaction import ...`; the three
      `databases/` files (`..reactions.equilibrium.reaction`,
      `..reactions.equilibrium.interphase`); `templates/stirred_tank/factory.py:48`.
      Outside `PyOMES/`: 32 live files (short form for `EquilibriumReaction`;
      the 5 tests needing a constraint name import it from
      `PyOMES.reactions.equilibrium.reaction` for now), 30 of them with
      `StoichiometryEntry` imports (decision 13): 15 standalone tests, 3
      validation tests, 7 validation notebooks, 3 tutorial notebooks, 2 `models/`
      files, and the two generators' string templates
      (`tests/validation/speciation/_generate_notebooks.py:60,1133,1479` and
      `docs/tutorials/ArXiv_preprint/_generate_notebooks.py:69`), each kept
      identical to its notebooks' cells, checked by script.
      Sanity: AST of each moved file equal to the original minus imports; old
      paths fail (`reactions.equilibrium` with `ImportError`,
      `reactions.phase_equilibria` and `reactions.plots` with
      `ModuleNotFoundError`); the three engine lazy paths exercised (NR engine,
      tableau and Bisection engine solves in the suite; also called directly);
      fresh-interpreter import of every entry point in several orders, including
      `PyOMES.chemical_equilibrium` before `PyOMES.reactions`; all 10 changed
      notebooks re-run from the scratchpad; generator templates match their
      notebook cells; both generators compile; no bare LF beyond the known 3;
      full suite **2098 passed**.
      _Notes: done 2026-09-24. Suite before: **2098 passed** (2m01s); after:
      **2098 passed**, 0 failed, 166 warnings, 1m52s. Moves by plain filesystem
      move; each moved file's AST with imports removed is identical to the
      original, and line endings and trailing bytes are unchanged (all CRLF). New
      docstring-only `equilibrium/__init__.py` (wording reviewed in checkpoint 4).
      Package edits, all import lines: `reaction.py` (4), `interphase.py` (7),
      `plots.py` (3, two of them under `TYPE_CHECKING`), `reactions/__init__.py`
      (2), `reaction_system.py:49,54`, the three engine lazy imports (now
      absolute `PyOMES.reactions.equilibrium.reaction`), `databases/`
      (`anaerobic_digestion.py:33-34`, `aqueous.py:21`, `bioprocess_basic.py:31-32`)
      and `stirred_tank/factory.py:48`. Outside `PyOMES/`, 32 files by script: 38
      adjacent `equilibrium`/`stoichiometry` pairs folded into
      `from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry`
      (aliases kept), 26 constraint-name imports pointed at
      `PyOMES.reactions.equilibrium.reaction` for now, and 4 lone
      `StoichiometryEntry` imports given their own short-form line. With
      checkpoint 1's notebook line that is all 45 decision-13 lines. Existing
      short-form lines next to a rewritten one (for example
      `from PyOMES.reactions import HenryEquilibrium` in `02_nr_engine_basics`)
      were left separate. **Added by the owner's decision (2026-09-24):** in
      `adm1/base.py` and `adm1/bsm2.py` only, the four adjacent
      `PyOMES.reactions` lines, including the unchanged
      `reaction_system` path, became one parenthesised short-form import; the
      other six files that deep-import `ReactionSystem` are left alone. Checks:
      old paths fail (`reactions.equilibrium` with `ImportError`,
      `reactions.phase_equilibria` and `reactions.plots` with
      `ModuleNotFoundError`); 12 fresh-interpreter import orders pass, including
      `PyOMES.chemical_equilibrium` and `PyOMES.chemistry.partition` first, each
      asserting that root and deep names are the same objects; the three engine
      lazy imports called directly (`build_tableau`, and `from_reactions` on both
      engines). Notebooks: only `source` changed, one import line per notebook.
      Generators: each rewritten line appears in its template; the ArXiv cell
      matches its template verbatim; the validation generator composes cells from
      pieces, so no cell matched verbatim before either, and every line of each
      cell is in the generator, as before (one docstring line in `06` already
      differed, unchanged). All 10 changed notebooks re-run from the scratchpad,
      every code cell, OK in 0.0-5.0 s each; `phreeqpython` is installed, so the
      PHREEQC-guarded cells ran. No files left in the repo._
- [x] 3. **Split `equilibrium/reaction.py`.** `constraint.py` gets
      `EquilibriumConstraint`, `vant_hoff_log_K`, `classify_equilibrium_constraint`,
      `_LOG10_E`, and the imports they need (`math`, `typing` names,
      `R_J_PER_MOL_K`, `StoichiometryEntry`, `phases_from_entries`; discrepancy 3);
      `reaction.py` keeps `EquilibriumReaction` and its imports (including
      `phases_from_entries`). New module docstring for `constraint.py`; the
      classifier's `EquilibriumReaction` reference gets its full path
      (discrepancy 5). Repoint to `constraint`: `reaction_system.py:49`,
      `interphase.py`, the NR engine and tableau lazy imports, the Bisection
      engine (two lines: `constraint` and `reaction`, discrepancy 4), and the 5
      tests.
      Sanity: the definitions in the two files, taken together, have the same AST
      as the pre-split file's (imports and module docstring removed); no unused
      import in either file; `vant_hoff_log_K` results bit-identical before and
      after on a grid of constraints and temperatures; full suite **2098 passed**.
      _Notes: done 2026-09-24. Suite before: **2098 passed** (1m50s); after:
      **2098 passed**, 0 failed, 166 warnings, 1m48s. `constraint.py` holds lines 66-152 of the pre-split file
      (the protocol, `vant_hoff_log_K`, `classify_equilibrium_constraint`) plus
      `_LOG10_E`, with its own module docstring and only the imports it uses
      (`math`, five `typing` names, `R_J_PER_MOL_K`, `StoichiometryEntry`,
      `phases_from_entries`); `reaction.py` keeps its module docstring and
      `EquilibriumReaction`, and drops `math`, `R_J_PER_MOL_K` and the `typing`
      names it no longer uses. Every top-level definition of the pre-split file is
      in exactly one of the two files with an identical AST; the one text change is
      the classifier docstring's reference to `EquilibriumReaction`, now a full
      path (discrepancy 5). No unused import in either file (AST check; `pyflakes`
      is not installed). Repointed to `constraint`: `interphase.py`,
      `reaction_system.py:49`, the NR engine and tableau lazy imports, and 26
      imports in the 5 tests (`test_chemistry_database.py` 1,
      `test_equilibrium_classification.py` 9, `test_equilibrium_constraint.py` 11,
      `test_equilibrium_constraint_dual_role.py` 3, `test_partition_model.py` 2).
      The Bisection engine's lazy import is now two lines, `constraint` and
      `reaction` (discrepancy 4). `equilibrium/__init__.py`'s file list gained
      `constraint.py`. While writing `constraint.py`'s docstring, a first draft
      said the engines accept any conforming type; the Bisection engine accepts
      only `EquilibriumReaction` for single-phase items
      (`bisection/engine.py:233`), so the docstring says so. Checks: a fingerprint
      of `vant_hoff_log_K` (48 values over 7 constraints and 8 temperatures,
      including one 1e-11 K off the reference) and `classify_equilibrium_constraint`
      (7 labels) is bit-identical before and after (SHA-256 `91bb33a9...`); the old
      paths `from PyOMES.reactions.equilibrium.reaction import
      EquilibriumConstraint` (and the other two names) raise `ImportError`; 8
      fresh-interpreter import orders pass, each asserting that root, deep and
      re-imported names are the same objects; the three engine lazy imports called
      directly, and the Bisection engine still rejects a non-constraint with
      `ValueError`; no bare LF. No notebook or generator imports a constraint
      name, so none was re-run._
- [ ] 4. **Docstrings, comments and live docs** (decision 10). Every dotted,
      slash and bare path in discrepancies 6 and 11, re-derived against the new
      files (line citations re-read, not shifted). Subpackage `__init__.py`
      docstrings say what each folder holds; `interphase.py`'s module docstring
      says it holds the named physical-law constraints that double as
      `PartitionModel`s. `docs/architecture.md` tree shows the new layout with
      `rate_laws.py` and `plots.py`; `README.md:106` and `PyOMES/README.md:15` stop
      saying `ReactionBuilder` builds `EquilibriumReaction` objects;
      `README.md:170` names `PyOMES.reactions.equilibrium.plots`. `OPEN_WORK.md`:
      the entries listed, the `plot_vant_hoff` entry repointed, and the
      three-dots entry updated to 11 lines across 6 files.
      Sanity: every changed relative link resolves; the touched `>>>` examples
      run by hand; changed `.py` files have identical ASTs with docstrings
      removed (no code change); full suite **2098 passed**.
- [ ] 5. **Sweep.** Search every old path in every form (dotted, slash,
      backslash, relative, bare filename; old path followed by a non-`.`
      character, decision 3) across `.py`, `.ipynb`, `.md`, `.toml`, `.yml`,
      `.json`, with plain `grep` and the Grep tool. Remaining hits only in
      `shipped/`, `docs/dev/ideas/`, the design note, this checklist, gitignored
      paths, `test_bsm2_reference.py:351`, `upcoming/README.md:153`, and bare
      `builder.py` mentions that name `templates/stirred_tank/builder.py`. AST
      import audit of every `.py` and notebook code cell: every `PyOMES`/`models`
      import resolves, including each imported name. Each old path fails to
      import with the expected error. Every notebook changed in the phase re-run
      once more from the scratchpad.
      Sanity: full suite **2098 passed**.
- [ ] 6. **Seam guard** (decision 12). Generalise `test_package_layering.py`'s
      resolver; add the `reactions/` layout test and its synthetic self-check
      cases. Also run the new check against a copy with an injected
      `kinetic -> equilibrium` import and a module-level `reaction_system` import
      in `equilibrium/plots.py`, to see it fail.
      Sanity: both tests pass; full suite passes with the count up by the number
      of test functions added.

## Shipping

- [ ] Full test suite green on the branch
- [ ] `git checkout main`
- [ ] `git merge --no-ff reactions-subfolder-structure -m "Merge reactions-subfolder-structure: <summary>"`
- [ ] `git tag reactions-subfolder-structure-shipped` on the merge commit
- [ ] `git push` and `git push --tags` (separate commands)
- [ ] `git branch -d reactions-subfolder-structure` and
      `git push origin --delete reactions-subfolder-structure`
- [ ] Full suite on `main` after the merge
- [ ] Move the design note and this checklist to `docs/dev/implementation/shipped/`
      (plain filesystem move); add "Shipped" banners; repoint the moved docs' own
      relative links
- [ ] Update `upcoming/README.md`: remove the "Design discussions" entry, add a
      "Recently shipped" entry
- [ ] Repoint any live link to the moved docs (search for
      `REACTIONS_SUBFOLDER_STRUCTURE` across `docs/` and `OPEN_WORK.md`)
