# Phase Kickoff Checklist — partition-constraint-relocation

> Checklist for [`PARTITION_CONSTRAINT_RELOCATION.md`](PARTITION_CONSTRAINT_RELOCATION.md),
> the source of truth for motivation and design; do not restate it here. Where
> this checklist and the note disagree, this checklist wins: the note was written
> from a conversational audit and the "Re-verification" section below corrects
> it. See [`README.md`](README.md)'s "Branching and tagging convention". Modelled
> on
> [`EQUILIBRIUM_SET_RELOCATION_CHECKLIST.md`](../shipped/EQUILIBRIUM_SET_RELOCATION_CHECKLIST.md).

**Working rules**

- Checkpoint 1 is a pure move: no behaviour change beyond the import path.
  Checkpoints 2-4 change only docstrings, comments, notebooks, scripts and
  prose. Checkpoint 5 adds one test and edits `OPEN_WORK.md`. Checkpoint 6 adds
  tests for `RaoultEquilibrium`'s custom parameters and one more `OPEN_WORK.md`
  entry; it changes no code.
- The repo owner runs every `git` command (branch, add, commit, push). At each
  checkpoint: run the full suite first, edit, run the sanity check, report,
  stop, and hand over the commands (two `-m` flags, no attribution lines).
- No shims, aliases or re-exports in `chemistry/` (decision 3 below).
- Leave `docs/dev/implementation/shipped/` and `docs/dev/ideas/` untouched; they
  describe past states. Live notes under `upcoming/` and `OPEN_WORK.md` are
  repointed (checkpoint 4).
- Docs and docstrings describe current behaviour only: no phase or checkpoint
  labels, no pointers to this checklist or the design note.
- Anything that looks like a bug or dead code beyond scope is logged in
  `OPEN_WORK.md`, not fixed.

**Decisions** (settled 2026-09-23, before starting)

1. The three classes go to a new file, `PyOMES/reactions/phase_equilibria.py`
   (not merged into `equilibrium.py`: about 420 lines would move, taking it from
   352 to about 770, and the classes are a distinct family). The plural name says
   the file holds more than one kind of phase equilibrium (gas-liquid Henry and
   Raoult, solid-liquid Ksp).
2. `HenryEquilibrium`, `RaoultEquilibrium` and `KspEquilibrium` are exported from
   `reactions/__init__.py`, so the short user-facing form becomes
   `from PyOMES.reactions import HenryEquilibrium`. `EquilibriumConstraint`,
   `vant_hoff_log_K` and `classify_equilibrium_constraint` stay deep-only.
   Tests, notebooks and scripts use the short form; code inside `PyOMES/` follows
   the import style of its neighbours.
3. No re-export from `chemistry/`. Consequence, accepted: a checkpoint written by
   `Simulation.save_checkpoint` before this phase pickles the classes under
   `PyOMES.chemistry.partition` and will fail to load afterwards. `load_checkpoint`
   only compares the package version, which does not change, so it gives no
   warning. The repo has no such pickle fixtures.
4. `_kH_mol_L_atm_from_ref` stays in `chemistry/partition.py` (its other user,
   `MultispeciesVLEPartition`, stays), and `phase_equilibria.py` imports it from
   there. Importing an underscore name across packages has precedent:
   `partition.py` imports `_parse_stoichiometry` from `reactions/` today.
5. The moved code's design-note pointers are reworded, not copied: module
   docstring "`EQUILIBRIUM_CONSTRAINT_UNIFICATION` CP1", `§14.1` (twice),
   "`NR_PRECIPITATION_SPECIATION`" and "wrapping litmus test". Docstring wording
   only; no code changes.
6. The layering test is narrow: no module under `chemistry/` imports anything in
   `PyOMES` except `units` and `chemistry` itself, at any depth (function-level
   and `TYPE_CHECKING` included). It does not assert an acyclic package graph;
   two other package cycles exist (see Re-verification) and are logged, not fixed.
7. `test_partition_model.py` keeps its name (`README.md:152` lists it by name).
8. Sequencing: `EquilibriumSet` shipped first, so after this phase `chemistry/`
   has no outbound dependency except `units`. `EXPLICIT_SPECIES_RESOLUTION.md`
   Phase 2 is independent; if it lands later, its edits happen in
   `phase_equilibria.py`, where `_resolve_species` moves unchanged.
9. `RaoultEquilibrium`'s water values (`P_sat_ref`, `dH_vap`, `T_ref`,
   `C_water_mol_L`) are already constructor arguments, so no API change is made.
   Its defaults stay as they are (a private constant for two, inline literals for
   two). This phase adds tests that custom values work (checkpoint 6) and logs the
   lack of one source for water properties in `OPEN_WORK.md`. Making the defaults
   public named constants, or consolidating the copies of the water values, are
   separate changes and are out of scope here.

## Pre-flight

- [x] `git status -sb` clean, `main` level with `origin/main` (checked
      2026-09-23)
- [x] Full-suite baseline on `main`: **2086 passed**, 0 failed, 166 warnings,
      2m41s (`python -m pytest -p no:cacheprovider -q`)
- [x] Branch created off `main`: `partition-constraint-relocation`
- [x] This checklist committed on that branch as the first commit (`5af2cc0`)

## During

- [x] Plan doc exists in `docs/dev/implementation/upcoming/`
      ([`PARTITION_CONSTRAINT_RELOCATION.md`](PARTITION_CONSTRAINT_RELOCATION.md))
- [ ] Checkpoints tracked below, one commit each
- [ ] **If work stalls:** add a status banner to the top of the plan doc at once

### Re-verification of the design note (2026-09-23, on `main`)

Searched every file type (`.py`, `.ipynb`, `.md`, `.toml`, `.yml`, `.json`) with
plain `grep` (which also covers gitignored and hidden paths: `notes/`, `scratch/`,
`PyOMES.egg-info/`, `.github/`, `.claude/`) and the Grep tool, for
`HenryEquilibrium`, `RaoultEquilibrium`, `KspEquilibrium`, `_resolve_species`,
`_kH_mol_L_atm_from_ref`, `_P_SAT_REF`, `_T_REF_WATER`, and every spelling of the
`chemistry.partition` path. 64 files name at least one moved class.

**Confirmed as the note says**

- Every `partition.py` line number in the note's table (54-71, 105-250, 259-385,
  390-510, 515-557, 560-563, 566-631; module docstring 16-22; imports 43-44).
- The `chemistry -> reactions` edge is exactly `partition.py:43-44` (plus a
  `TYPE_CHECKING` import at :52). `reactions -> chemistry` is the legitimate
  opposite edge (`builder.py:41-42`, `equilibrium.py:49`, `kinetic.py:35`,
  `reaction_system.py:48`, `stoichiometry.py:25-26`; the note names only two).
- Production import sites: `databases/anaerobic_digestion.py:33`,
  `databases/bioprocess_basic.py:31`, `templates/stirred_tank/factory.py:47`,
  `models/vlmodels/adm1/base.py:1055` (lazy). `core/` imports only
  `PartitionModel` (`gas_liquid_link.py:97`, `transfer_models.py:57`,
  `databases/database.py:30`, all of which stay).
- `reactions/reaction_system.py`, `reactions/equilibrium.py` and
  `chemical_equilibrium/` handle the classes only through the
  `EquilibriumConstraint` protocol; they name them in docstrings, comments and
  error-message strings only. No `isinstance` or import anywhere.
- `tests/standalone/test_import_graph_acyclic.py` passes because it counts only
  module- and class-level imports. Every entry-point import (`PyOMES`,
  `.chemistry`, `.reactions`, `.reactions.stoichiometry`, `.reactions.equilibrium`,
  `.chemistry.partition`, `.databases`, `.templates`, `.core`) works in a fresh
  interpreter today, so the note's import-order risk is latent, not live.
- No mock-patch or `importlib` string references to the module path; the install
  is editable, so no reinstall is needed. `PyOMES.egg-info/SOURCES.txt` lists
  `partition.py` (still exists).

**Discrepancies (all folded into the checkpoints below)**

- **"The package graph becomes a DAG" is false.** Package-level graph with every
  import counted (including function-level and `TYPE_CHECKING`): today the
  strongly connected groups are {`chemical_equilibrium`, `chemistry`, `reactions`}
  and {`control`, `core`}. Simulating the move removes `chemistry` from the first,
  leaving `chemical_equilibrium` <-> `reactions` (function-level only:
  `reaction_system.py:261,273` one way; `bisection/engine.py:203`,
  `nr/engine.py:239`, `nr/tableau.py:401` the other) and `control` <-> `core`
  (module-level both ways, e.g. `control/__init__.py:11`, `core/phases.py:25`).
  The move adds no new package edge (`databases -> reactions`,
  `templates -> reactions` and `reactions -> chemistry` all exist). So checkpoint 5
  narrows the `OPEN_WORK.md` entry and records the two remaining cycles instead of
  closing it, and the layering test is scoped to `chemistry` (decision 6).
- **Checkpoints pickle these classes.** `Simulation.save_checkpoint` pickles the
  full object graph, which includes `HenryEquilibrium` via `partition_models`
  (decision 3).
- **The "31 external files" count matches by class name, but the composition
  differs.**
  - `docs/architecture.md:363` is a false positive (`HenryEquilibriumInterface`,
    an unrelated legacy class). It needs no class edit, only the new file added to
    its `reactions/` tree.
  - Real edits outside `PyOMES/`: 17 test files (16 in `tests/standalone/`, plus
    `tests/validation/speciation/test_iron_oxidation.py`); 7 notebooks
    (`ArXiv_preprint/02,03`, `ChemicalEquilibriumProtocol/01,02`,
    `reactions/partition_model`, `tests/validation/speciation/07,08`); 2 scripts
    (`ArXiv_preprint/_generate_notebooks.py` at lines 828 and 1209, and
    `D2C_workshop/raw_construction.py:56`); `models/vlmodels/adm1/base.py:1055`.
  - No edit needed: `ChemicalEquilibriumProtocol/03_phreeqc_engine_basics.ipynb`
    and `0_README.ipynb`, and the markdown cells in the other notebooks, which use
    bare class names with no path; `docs/tutorials/reactions/README.md:17` (bare
    names).
  - Notebooks are not run by CI (`.github/workflows/tests.yml` runs bare `pytest`
    over `tests/standalone` and `tests/validation`) and no test executes one.
    Saved outputs contain no module path. Only the ArXiv notebooks have a
    generator, so `_generate_notebooks.py` and notebooks 02/03 must be edited
    together. `07/08` are hand-authored.
  - Gitignored and left alone: `scratch/ArXiv_preprint/02,03` (untracked copies
    with the same import, which will break) and `notes/chemistry_notes.md`
    (personal notes with `partition.py` line references).
- **Stale docstring paths the note missed** (path forms that go stale; bare-name
  prose and error strings stay accurate and are not edited):
  `core/gas_liquid_link.py:118,573` and `core/transfer_models.py:70,112`
  (`~PyOMES.chemistry.HenryEquilibrium`), and `databases/anaerobic_digestion.py:11`,
  besides the ones the note lists (`gas_liquid_link.py:68`,
  `bisection/engine.py:147-149`, `reactions/equilibrium.py:77-79,118-120`,
  `reactions/reaction_system.py:9-11`). `chemistry/__init__.py`'s docstring does
  not name the classes, so it needs no edit.
- **Live notes with `partition.py` references the note does not mention:**
  `upcoming/EXPLICIT_SPECIES_RESOLUTION.md` (`partition.py:54-71`, `63-65`,
  `145-146`, `245-246`, `310-311`, `380-381`, and the Phase 2 heading),
  `upcoming/PHENOMENA_PROTOCOL.md:53,259` (a link labelled `partition.py:221-253`
  and a link to the file), `upcoming/README.md:27,37,46`, and `OPEN_WORK.md:170,
  206,465,590`, which list `chemistry/partition.py` among files.
  `upcoming/README.md:212` sits in a shipped-phase entry and stays.
- **The moved region also holds** `_P_SAT_REF` and `_T_REF_WATER` (255-256), used
  only by `RaoultEquilibrium`'s defaults; the note's table omits them. The
  `TYPE_CHECKING` block (`ThermoFramework`, `StoichiometryEntry`) is unused or
  redundant and is dropped rather than carried over.
- `PyOMES/README.md:13` still says `chemistry/` holds "acid-base equilibrium
  sets", a leftover from `equilibrium-set-relocation`; fixed in checkpoint 4
  because this phase changes what `chemistry/` holds again.
- `equilibrium.py` after a merge would be about 770 lines, not the note's ~700.
- **Found while working (checkpoint 1):** `RaoultEquilibrium` can already be
  built with custom `P_sat_ref`, `dH_vap`, `T_ref` and `C_water_mol_L`
  (`RaoultEquilibrium(P_sat_ref=0.0313, dH_vap=43990.0).P_sat(310.15)` differs
  from the default), but no test does so. The water values also exist in three
  places with no shared source: `phase_equilibria.py` (`_P_SAT_REF`, and inline
  44011.0 and 55.51), `chemical_equilibrium/engines/nr/engine.py:55`
  (`_C_WATER_MOL_L`, derived from density and molar mass, about 55.51) and
  `models/vlmodels/adm1/base.py:1050` (a literal 55.51).
- Unrelated, left alone: `upcoming/README.md:9` has a stray `+-*` line.

### Checkpoints

- [x] 1. Move the classes. Steps 1 and 2 of the note land together, since there is
      no shim.
      Create `PyOMES/reactions/phase_equilibria.py` holding `_resolve_species`,
      `_P_SAT_REF`, `_T_REF_WATER`, `HenryEquilibrium`, `RaoultEquilibrium` and
      `KspEquilibrium`, with a module docstring describing the file (dual
      `PartitionModel`/`EquilibriumConstraint` role; the activity-correction
      section moves here from `partition.py`) and the design-note pointers
      reworded (decision 5). Imports, in the `..` style of `reactions/`:
      `from ..chemistry import common_species`,
      `from ..chemistry.species import Species`,
      `from ..chemistry.partition import _kH_mol_L_atm_from_ref`,
      `from ..units import R_J_PER_MOL_K as _R_J_MOL` and `R_L_ATM_PER_MOL_K`,
      `from .equilibrium import vant_hoff_log_K`,
      `from .stoichiometry import StoichiometryEntry, _parse_stoichiometry`.
      `chemistry/partition.py`: delete the moved code and the `reactions` and
      `TYPE_CHECKING` imports; rewrite the module docstring to describe only the
      two protocols and `MultispeciesVLEPartition`; the remaining imports are
      `math`, `dataclass`/`field`, `typing` and `R_L_ATM_PER_MOL_K`.
      `chemistry/__init__.py`: drop the three names and their `__all__` entries.
      `reactions/__init__.py`: import and export the three (module docstring and
      `__all__`). Production imports: `databases/anaerobic_digestion.py:33`,
      `databases/bioprocess_basic.py:31`, `templates/stirred_tank/factory.py:47`
      (split: `HenryEquilibrium` from `reactions`, `PartitionModel` stays),
      `models/vlmodels/adm1/base.py:1055`. Tests: update the import in the 17 files
      (split mixed lines such as `from PyOMES.chemistry import HenryEquilibrium,
      PartitionModel` and `Species, HenryEquilibrium`).
      Sanity: fresh-interpreter `import` of each entry point listed under
      Re-verification, in several orders; import-graph guard green;
      `test_partition_model.py`, `test_equilibrium_constraint*.py`,
      `test_hpc_checkpointing.py` and `test_bsm2_reference.py` pass; full suite
      **2086 passed**; `chemistry/partition.py` contains no `reactions`; no
      `.py` outside `docs/` still imports the classes from `PyOMES.chemistry`.
      _Notes: done 2026-09-23. Suite before: **2086 passed** (3m33s); after:
      **2086 passed**, 0 failed, 166 warnings, 4m02s. The six targeted files give
      147 passed. Every entry point imports in a fresh interpreter, including
      `phase_equilibria` first and `partition` before `reactions`;
      `chemistry/` code has no `reactions` import (only a `>>>` example in
      `common_species.py:23`, which is a docstring, not an import). Neither
      `partition.py` nor `phase_equilibria.py` has an unused import. The moved code
      is verbatim except the module docstring, the `_resolve_species` docstring
      (it claimed to skip "the extra hop through `reactions`", no longer true in
      `reactions/`), the reworded design-note pointers (decision 5), and the dropped
      `TYPE_CHECKING` block. One correction while rewriting: the old module
      docstring said only "single-ion" `KspEquilibrium` satisfies
      `EquilibriumConstraint`; the class satisfies it for any stoichiometry, and
      only its `PartitionModel` role is single-ion, so the new docstring says that.
      Production imports follow their neighbours: deep
      `..reactions.phase_equilibria` in `databases/` and `factory.py`, package-level
      `PyOMES.reactions` in `adm1/base.py`. The 17 test files were rewritten by a
      script (77 import lines, 5 of them split into two lines) because they use CRLF
      line endings and the Edit tool would have risked bare LFs; the diff shows no
      other change and no LF-only lines. Remaining old-path references, all for
      later checkpoints: the docstring paths listed in checkpoint 2, 7 notebooks
      plus the gitignored `scratch/` copies, and the two `docs/tutorials` scripts._
- [x] 2. Repoint docstring and comment paths in `PyOMES/`: `~PyOMES.chemistry.
      partition.X` and `~PyOMES.chemistry.HenryEquilibrium` become
      `~PyOMES.reactions.phase_equilibria.X` in `bisection/engine.py:147-149`,
      `databases/anaerobic_digestion.py:11`, `reactions/equilibrium.py:77-79,
      118-120`, `reactions/reaction_system.py:9-11`,
      `core/gas_liquid_link.py:118,573` and `core/transfer_models.py:70,112`; the
      `>>>` example at `gas_liquid_link.py:68` imports from `PyOMES.reactions`.
      Nothing else in these files changes.
      Sanity: no `PyOMES.chemistry.(partition.)?(Henry|Raoult|Ksp)` left under
      `PyOMES/`; the `gas_liquid_link.py` example runs; full suite **2086 passed**.
      _Notes: done 2026-09-24. Suite before: **2086 passed** (8m01s); after:
      **2086 passed**, 0 failed, 166 warnings, 5m44s. The
      `gas_liquid_link.py:68-84` example (construct `HenryEquilibrium` from
      `PyOMES.reactions`, build a `KineticGasLiquidLink`) runs unchanged. No other
      text in the nine edited lines changed._
- [x] 3. Notebooks and scripts: change the import in the 7 notebooks' code cells,
      `_generate_notebooks.py:828,1209` (kept identical to notebooks 02/03) and
      `raw_construction.py:56`. Source cells only; saved outputs are not re-run
      into the repo.
      Sanity: run each edited notebook's cells that construct the classes (or the
      whole notebook where it runs quickly) from a scratch copy, and
      `raw_construction.py` if it is short; confirm generator and notebook cells
      match by script; full suite **2086 passed**.
      _Notes: done 2026-09-24. Suite before: **2086 passed** (3m07s); after:
      **2086 passed**, 0 failed, 166 warnings, 2m27s (docs/scripts/notebooks are
      outside `testpaths`, so no test-count change was expected). Each of the 7
      notebooks got exactly one changed line (verified by diff against a backup
      taken before editing), splitting the moved name onto its own
      `from PyOMES.reactions import ...` line beside the existing
      `from PyOMES.chemistry import ...`/`from PyOMES.reactions import ...` line,
      same pattern as checkpoint 1's test-file rewrite. Editing was scripted for
      the same CRLF reason as checkpoint 1; the script additionally had to
      preserve each file's exact trailing bytes after the final `}` byte-for-byte
      (three different endings across the 7: none, `\r\n`, and a lone trailing
      `\r` with no `\n` — all three now match their originals), and CRLF counts
      match old-to-new for every file. `raw_construction.py`'s single import line
      was merged into its existing `from PyOMES.reactions import (...)` block, and
      the file was run end-to-end (`python docs/tutorials/D2C_workshop/
      raw_construction.py`): exit 0, 0.474 s wall-clock, prints unchanged from a
      run before this checkpoint. `_generate_notebooks.py`'s two sites were first
      merged into their existing `from PyOMES.reactions import (...)` blocks (to
      match `raw_construction.py`'s style), which did not match the corresponding
      notebook cells; switched to a separate `from PyOMES.reactions import
      HenryEquilibrium` line instead (matching the scripted notebook edit's shape),
      confirmed identical to the `KINETIC_SETUP`/`CSTR_SETUP` string constants by
      direct string comparison against notebooks 02 and 03's first code cell — not
      run (would regenerate every notebook the script produces, not just 02/03,
      risking unrelated diffs from output/timestamp drift). The generator was not
      executed for that reason; its two other sanity checks (syntax parse, string
      match) covered the intended check. The three false-positive grep hits this
      produced (`from PyOMES.chemistry import X\nfrom PyOMES.reactions import
      Henry...` spanning the split) were confirmed benign — the two matched
      substrings sit either side of the new line break, not a leftover old-style
      import. Gitignored `scratch/ArXiv_preprint/02,03` untouched, as decided._
- [ ] 4. Docs and live notes: `docs/architecture.md` (add `phase_equilibria.py` to
      the `reactions/` tree); `PyOMES/README.md:13` (drop "and acid-base
      equilibrium sets", and mention the constraints under `reactions/` if that
      row lists contents); repoint `partition.py` references in
      `EXPLICIT_SPECIES_RESOLUTION.md` (paths and line numbers re-read against the
      new file), `PHENOMENA_PROTOCOL.md:53,259` and `upcoming/README.md:27,37,46`;
      add the new file beside `chemistry/partition.py` in `OPEN_WORK.md:170,206,
      465,590` where the sentence is about the moved code.
      Sanity: repo-wide search for `chemistry.partition` and `chemistry/partition`
      finds only shipped docs, `docs/dev/ideas/`, the design note, this checklist,
      gitignored paths, and references to what stays (`PartitionModel`,
      `MultispeciesVLEPartition`); every changed link resolves; full suite
      unchanged.
- [ ] 5. Layering test and `OPEN_WORK.md`. New
      `tests/standalone/test_package_layering.py`: parse every `.py` under
      `PyOMES/chemistry/` with `ast.walk` (all imports at any depth, including
      `TYPE_CHECKING`), resolve relative imports, and assert the set of `PyOMES`
      subpackages imported is a subset of `{chemistry, units}`. It includes a
      synthetic self-check (a function-level import, a `TYPE_CHECKING` import and
      a `from PyOMES import reactions` form must each be caught) so an empty result
      is meaningful. `OPEN_WORK.md`: narrow the "Package-level layering" entry
      (the `chemistry` <-> `reactions` cause is gone) and record the two remaining
      cycles with their edges (`control` <-> `core`, `chemical_equilibrium` <->
      `reactions`), unfixed.
      Sanity: the new test passes; full suite **2087 passed** (or one more per
      test function added).
- [ ] 6. Tests for `RaoultEquilibrium`'s custom parameters (decision 9), added to
      the Raoult section of `tests/standalone/test_partition_model.py`. No source
      change. Cases: `P_sat(T_ref)` equals a custom `P_sat_ref`; `P_sat` at another
      temperature follows Clausius-Clapeyron with a custom `dH_vap` and `T_ref`;
      `partition_ratio` scales with `C_water_mol_L` and inversely with `P_sat`;
      `log_K`, `dH_J_per_mol` and `T_ref_K` report the custom values; default
      construction is unchanged. If a case exposes a defect, it is logged in
      `OPEN_WORK.md`, not fixed here.
      `OPEN_WORK.md`: add an entry for the missing single source of truth for
      water properties (the three locations under Re-verification), including the
      option of public named constants for the defaults, unfixed.
      Sanity: the new tests pass; full suite passes with the count up by the number
      of test functions added.

## Shipping

- [ ] Full test suite green on the branch
- [ ] `git checkout main`
- [ ] `git merge --no-ff partition-constraint-relocation -m "Merge partition-constraint-relocation: <summary>"`
- [ ] `git tag partition-constraint-relocation-shipped` on the merge commit
- [ ] `git push && git push --tags`
- [ ] `git branch -d partition-constraint-relocation` and
      `git push origin --delete partition-constraint-relocation`
- [ ] Full suite on `main` after the merge
- [ ] Move the design note and this checklist to `docs/dev/implementation/shipped/`
      (plain filesystem move); add "Shipped" banners
- [ ] Update `upcoming/README.md`: remove the "Design discussions" entry, add a
      "Recently shipped" entry
- [ ] Repoint any live link to the moved docs (search for
      `PARTITION_CONSTRAINT_RELOCATION` across `docs/` and `OPEN_WORK.md`)
