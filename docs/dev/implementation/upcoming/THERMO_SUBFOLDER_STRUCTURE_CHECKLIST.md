
# Phase Kickoff Checklist — thermo-subfolder-structure

> Checklist for [`THERMO_SUBFOLDER_STRUCTURE.md`](THERMO_SUBFOLDER_STRUCTURE.md),
> the source of truth for motivation, groups, target layout and decisions; do not
> restate it here. Where this checklist and the note disagree, this checklist
> wins: the "Re-verification" section below corrects the note from a fresh search.
> See [`README.md`](README.md)'s "Branching and tagging convention". Modelled on
> [`../shipped/REACTIONS_SUBFOLDER_STRUCTURE_CHECKLIST.md`](../shipped/REACTIONS_SUBFOLDER_STRUCTURE_CHECKLIST.md)
> and
> [`../shipped/CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER_CHECKLIST.md`](../shipped/CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER_CHECKLIST.md).

**Working rules**

- Checkpoints 1-2 are pure refactors: file moves, two file splits and import
  rewrites, with no behaviour change. Checkpoint 3 changes `GasEOS` from a base
  class to a `Protocol` (decision 11) and adds a test. Checkpoint 4 replaces two
  private helpers with one public one, bit-identically (decision 12), and adds a
  test. Checkpoint 5 changes only docstrings, comments and prose. Checkpoint 6
  changes nothing unless the sweep finds a missed site. Checkpoint 7 adds a test.
- The repo owner runs every `git` command (branch, add, commit, push, merge,
  tag). Read-only local git (`status`, `log`, `diff`, `show`, `ls-files`) is
  fine; nothing that touches the network. Files are moved with a plain
  filesystem move, not `git mv`. Git commands are handed over one per line
  (PowerShell 5.1: never chained with `&&`); commit messages use two `-m` flags
  (a strapline and one body paragraph) and no attribution lines. The repo is in
  OneDrive: a commit may ask "Rename ... index.lock ... failed. Should I try
  again? (y/n)"; answer `y`. No git command is run while a commit may be waiting.
- Line endings when staging (`core.autocrlf=true`): a file added under a new
  path is stored as LF unless staged with autocrlf off. Before handing over
  staging commands, check how the old path is stored (`git ls-files --eol`). A
  new file that comes from an `i/crlf` file is staged with
  `git -c core.autocrlf=false add <path>`; one from an `i/lf` file with a plain
  `git add`. Check `git diff --cached -M --stat` before committing. Edited files
  keep their stored line endings (git does not convert a file whose index copy
  already has CRLF).
- One checkpoint at a time: run the full suite first, edit, run the sanity
  checks, report what changed and the results as they are, record notes here,
  stop, and hand over the commit commands. Do not start the next checkpoint until
  told to.
- Full suite: `python -m pytest -p no:cacheprovider -q`.
- No re-export shims, aliases or compatibility layers. Old import paths stop
  working; every call site is updated in place. Subpackage `__init__.py` files
  are docstring-only.
- Proof of no behaviour change: the AST of each moved or split definition is
  identical to the original once imports are removed; numeric output of
  everything split is bit-identical (fingerprint by SHA-256, before and after);
  doc-only edits leave the AST identical with docstrings removed. Each old path
  is shown to fail to import.
- Every notebook changed in a checkpoint is re-run in that checkpoint, from the
  scratchpad directory (so no output folders land in the repo). Notebook edits
  touch only `source` lines, never saved outputs. No notebook is expected to
  change (see Re-verification).
- The tool sandbox caps each process at about 15 % of one CPU core. If a
  notebook's long cells do not finish, run it up to the slow cells plus a
  construction check and record that here as a deviation. Notebooks and the
  suite are never run at the same time.
- Nearly every file involved is CRLF in the working tree. Bulk edits use a script
  that preserves each file's line endings and exact trailing bytes; check for
  bare LF after every edit.
- Leave `docs/dev/implementation/shipped/` and `docs/dev/ideas/` untouched.
- Docs and docstrings describe current behaviour only: no phase or checkpoint
  labels, no pointers to this checklist or the design note.
- Anything that looks like a bug or dead code beyond scope is logged in
  `OPEN_WORK.md`, not fixed.

**Decisions** (settled 2026-09-24, before starting; 1-10 are the note's, in brief)

1. Layout: `liquid/` and `gas/` subfolders under `thermo/`; `framework.py` and
   `equilibrium_constants.py` stay at the top level.
2. `gas/ideal.py` is its own file, matching `liquid/ideal.py`.
3. `water_properties.py` and `factory.py` go into `liquid/`.
4. File names drop the suffix the folder supplies: `davies.py`, `sit.py`,
   `peng_robinson.py`. Class names are unchanged. None of the old module names
   is reused, so every old path fails with `ModuleNotFoundError`.
5. No shims; subpackage `__init__.py` files are docstring-only.
6. Package-root exports unchanged.
7. Consumers outside `PyOMES/` (here, three test files) switch every deep import
   of an exported name to `from PyOMES.thermo import ...`; only
   `DifferentiableLiquidModel` keeps a deep import, at
   `PyOMES.thermo.liquid.protocols`. Code inside `PyOMES/` keeps deep imports.
8. Moved files import within their folder relatively (`from .water_properties
   import ...`) and anything outside it absolutely (`from PyOMES.units import
   ...`). Staying files keep their relative style
   (`from .liquid.davies import DaviesLiquidModel`).
9. The moves are pure relocations. Behaviour-touching changes are separate
   checkpoints (decisions 11 and 12). Wiring `ThermoFramework.gas_eos` into
   `core/` is out of scope and stays in `OPEN_WORK.md`.
10. Docstrings of moved code describe current behaviour: phase labels and
    design-doc references are rewritten in the docstring checkpoint, along with
    `water_properties.py`'s stale `speciation/` path and `framework.py`'s
    `KineticGasLiquidLink` claim; the `OPEN_WORK.md` count is updated.
11. **`GasEOS` becomes a `@runtime_checkable` `Protocol`** (open question 1), in
    its own checkpoint after the moves. `IdealGasEOS` stops subclassing it and
    `PengRobinsonEOS` satisfies it structurally. The protocol's docstring states
    plainly that the two implementations' `partial_pressures_atm` return
    different quantities (partial pressures vs fugacities); the planned split
    into `partial_pressures_atm` and `fugacities_atm` stays in `OPEN_WORK.md`.
12. **One public `water_kg_per_L(T_K)` in `liquid/water_properties.py`** (open
    question 2), in its own checkpoint. It replaces the two private `_kg_per_L`
    copies (`davies.py`, `sit.py`), and `ionic_strength_molal_from_molar` uses
    it; all three are bit-identical today. Not exported from the package root
    (decision 6); every caller is in `liquid/`. SIT's two inline
    `water_density_kg_per_m3(T_K) / 1000.0` computations (no fallback, so a
    switch would change results for a bad density) are left and recorded in the
    `OPEN_WORK.md` entry, which is otherwise closed. The near-identical Davies and
    SIT Jacobian loops are left alone.
13. **`DifferentiableLiquidModel` is not exported from the package root** (open
    question 3).
14. **Seam guard** (open question 4): yes, as the last checkpoint, as a third
    test in `tests/standalone/test_package_layering.py`, reusing its `_imports()`
    resolver, with a synthetic self-check. Rule: any import from a module in
    `liquid/` or `gas/` that resolves into `PyOMES.thermo` must stay inside that
    module's own folder, in every form (module-level, lazy, `TYPE_CHECKING`,
    `from .. import x`, plain `import`). That forbids the other side,
    `framework`, `equilibrium_constants` and the package root.
15. **Sequencing** (open question 5): resolved. `reactions-subfolder-structure`
    shipped on 2026-09-24 (tag `reactions-subfolder-structure-shipped`); its
    effect on the note is recorded under Re-verification. No other phase branch
    exists (only `main`, checked 2026-09-24).
16. The note is already listed in `upcoming/README.md`'s "Design discussions";
    that entry is replaced by a "Recently shipped" entry at shipping.
17. New files from a split keep the stored line endings of the file they come
    from: the `liquid/` pieces of `liquid_phase_model.py` (`i/crlf`) are staged
    with autocrlf off; the `gas/` pieces of `gas_eos.py` (`i/lf`) with a plain
    `git add`.

## Pre-flight

- [x] `git status -sb` clean, `main` level with `origin/main` at `8c89048` as of
      the last fetch (checked 2026-09-24, no fetch run)
- [x] Full-suite baseline on `main`: **2100 passed**, 0 failed, 166 warnings,
      1m48s (`python -m pytest -p no:cacheprovider -q`)
- [x] Branch created off `main`: `thermo-subfolder-structure`
- [x] This checklist committed on that branch as the first commit (`ce2500a`)

## During

- [x] Plan doc exists in `docs/dev/implementation/upcoming/`
      ([`THERMO_SUBFOLDER_STRUCTURE.md`](THERMO_SUBFOLDER_STRUCTURE.md))
- [ ] Checkpoints tracked below, one commit each
- [ ] **If work stalls:** add a status banner to the top of the plan doc at once

### Re-verification of the design note (2026-09-24, on `main` at `8c89048`)

**Method.** A script walked the whole working tree, including gitignored
`notes/`, `scratch/` and `PyOMES.egg-info/`, and hidden `.github/` and
`.claude/`, skipping only `.git/`, `__pycache__/` and `.pytest_cache/`. It read
`.py`, `.ipynb`, `.md`, `.toml`, `.yml`/`.yaml`, `.json`, `.txt`, `.cfg`, `.ini`
and `.rst` files and matched `thermo.<module>`, `thermo/<module>` and
`thermo\<module>` (not followed by a word character), relative
`from .<module> import`, and bare `<module>` / `<module>.py` for the four
distinctive names (bare `factory.py` only on lines that also say `thermo`).
Each hit was classified: in `.py` files by tokenizer and AST (module-level
import, lazy import, `TYPE_CHECKING` import, string or docstring, comment,
code); in notebooks by cell type (code, markdown, saved output). The same script
listed every import of `PyOMES.thermo` or a relative `thermo` in every `.py`
file and notebook code cell, including generator string templates. Cross-checked
with plain `grep -r` (63 old-path occurrences in 17 files) and the Grep tool (the
same, minus the 3 in gitignored `PyOMES.egg-info/SOURCES.txt`; it counts matching
lines, so the 9 occurrences in `THERMODYNAMIC_MODEL_ARCHITECTURE.md` show as 7).
Separate searches for dynamic imports (`importlib`, `__import__`, `sys.modules`,
`mock.patch`) naming `thermo`, for pickles, and for the names each split piece
uses. Line counts with `wc -l`; line endings by byte count and
`git ls-files --eol`.

**Confirmed as the note says**

- 8 files, 1,459 lines; every per-file count exact.
- Intra-package imports: `liquid_phase_model.py:34` and `sit_liquid_model.py:27`
  import `.water_properties`; `factory.py:6-7` imports `.liquid_phase_model`
  (`ActivityModel`, `DaviesLiquidModel`, `IdealLiquidModel`) and
  `.sit_liquid_model`; `framework.py:20-21` imports `.liquid_phase_model` and
  `.gas_eos`; `water_properties.py`, `gas_eos.py` and `equilibrium_constants.py`
  import nothing from `thermo/`. Outside `thermo/`, only `..units`
  (`framework.py:19`, `equilibrium_constants.py:24`, `gas_eos.py:37`). Neither
  group imports the other.
- **Consumers of moved paths:** exactly three files, all tests, 11 import lines,
  all lazy (inside test functions): `test_nr_gas_liquid_cp2.py` (7 lines),
  `test_gas_eos.py:28,79,93`, `test_liquid_phase_model.py:244`. Every name they
  pull is exported from the root except `DifferentiableLiquidModel`.
  `test_gas_eos.py:2` names the old path in its module docstring.
- Unmoved paths: `thermo.framework` is imported by `databases/aqueous.py:24`
  (module level), `databases/database.py:31` (`TYPE_CHECKING`) and
  `test_thermo_framework.py:138,145` (`THERMO_IDEAL`, `THERMO_DAVIES`, not
  exported, so they stay deep); `thermo.equilibrium_constants` by
  `engines/bisection/acid_base.py:43`, `engines/nr/engine.py:47`,
  `engines/nr/tableau.py:43` and `test_equilibrium_constants.py:15`.
- Everything else imports from the package root: `PyOMES/__init__.py:20`, the
  engines (`acid_base.py:42,887`, `bisection/engine.py:44`, `nr/engine.py:46`),
  `docs/tutorials/reactions/chemistry_database.ipynb` (cell 6), validation
  notebooks `03`, `04`, `05`, the validation generator's templates (`:767`,
  `:1134`), `test_saturation_index.py:25`, and the tests. **No notebook or
  generator names a moved path**, in code, markdown or saved output, so no
  notebook changes in this phase.
- No `mock.patch`, `importlib`, `__import__` or `sys.modules` string names a
  `thermo` module; no `.pkl`/`.pickle`/`.joblib` file anywhere. No hit in
  `notes/`, `scratch/`, `.github/`, `.claude/`, or any `.toml`/`.yml`/`.json`
  file. `test_package_layering.py:148` is synthetic detector input.
- Nothing in `PyOMES/` outside `thermo/` constructs `IdealGasEOS` or
  `PengRobinsonEOS` or reads `ThermoFramework.gas_eos`; no `isinstance(...,
  GasEOS)` and no `GasEOS()` call anywhere. `PengRobinsonEOS` does not subclass
  `GasEOS` (its docstring at `gas_eos.py:239` says it implements it).
- **Split boundaries.** `liquid_phase_model.py`: `DifferentiableLiquidModel`
  41-84, `LiquidPhaseModel` 87-135, `ActivityModel` 138-143 (the protocols,
  41-143), `_kg_per_L` 146-151, `_LN10` 154, `IdealLiquidModel` 157-182,
  `DaviesLiquidModel` 185-304 (uses `_kg_per_L` at :288 and `_LN10` at :297).
  `gas_eos.py`: `GasEOS` 40-47, `IdealGasEOS` 50-63, and 66-506 for
  Peng-Robinson (`OMEGA_A`/`OMEGA_B`, `CriticalProperties`, `BIOGAS_SPECIES`,
  `BIOGAS_KIJ`, `_get_kij`, `_kappa`, `_alpha`, `_a_pure`, `_b_pure`,
  `_solve_cubic_Z`, `PengRobinsonEOS`), which references neither `GasEOS` nor
  `IdealGasEOS` except in that docstring. Needed imports: `gas/ideal.py`:
  `dataclass`, `Dict`, `R_L_ATM_PER_MOL_K as R`, `.protocols.GasEOS`;
  `gas/peng_robinson.py`: `math`, `dataclass`, `field`, `Dict`, `Tuple`, `R`;
  `gas/protocols.py`: `Dict`.
- The two `_LN10` values (`float(np.log(10.0))` in Davies, `math.log(10.0)` in
  SIT) are the same float; the two `_kg_per_L` bodies are identical.
- Doctests are not run by the suite (`pyproject.toml` sets only `testpaths`).
- Packaging: `setup.py`'s `find_packages(include=["PyOMES", "PyOMES.*"])` will
  find `thermo.liquid` and `thermo.gas`; `PyOMES` imports from the working tree
  (editable install), so no reinstall is needed. `PyOMES.egg-info/SOURCES.txt`
  (gitignored, generated, already stale: it lacks `gas_eos.py` and `factory.py`)
  is left alone.
- `docs/dev/ideas/THERMODYNAMIC_MODEL_ARCHITECTURE.md` still uses `src/thermo/`
  paths (9 occurrences); left alone. `OPEN_WORK.md:20,35` narrate past moves
  into `thermo/gas_eos.py`.

**Discrepancies** (none changes a decision; all are folded into the checkpoints)

1. **Line endings.** The note says "All eight files are LF". In the working tree
   all eight are CRLF with a trailing `\r\n`. In the index, `__init__.py`,
   `framework.py`, `liquid_phase_model.py` and `sit_liquid_model.py` are
   `i/crlf`; `equilibrium_constants.py`, `factory.py`, `gas_eos.py` and
   `water_properties.py` are `i/lf`. Staging per the working rule and
   decision 17: `liquid/sit.py`, `liquid/protocols.py`, `liquid/ideal.py`,
   `liquid/davies.py` with autocrlf off; `liquid/water_properties.py`,
   `liquid/factory.py` and the three `gas/` pieces with a plain `git add`.
2. **Line numbers moved by the reactions phase.** It merged two import lines in
   `test_nr_gas_liquid_cp2.py`, so every thermo import there sits 2 lines
   earlier than the note says: `liquid_phase_model` at 345-347, 352, 356, 378,
   389; `sit_liquid_model` at 351, 369; `DifferentiableLiquidModel` at 346 and
   352. Lines 345-347 import `DaviesLiquidModel` and `DifferentiableLiquidModel`
   together, so under decision 7 they become two imports (short form and deep).
   The `docs/architecture.md` `thermo/` block is at 397-400, not 393-396.
   `engines/nr/tableau.py` imports `equilibrium_constants` at :43, not :42.
3. **Baseline:** 2100 passed, not 2098 (the reactions phase added two tests).
4. **`liquid_phase_model.py:32` imports `Optional` and never uses it.** No piece
   of the split needs it, so it is dropped (an import-only difference; the
   definitions' ASTs are unaffected). `liquid/protocols.py` needs `numpy` for
   `DifferentiableLiquidModel`'s `np.ndarray` return annotation, plus `Dict`,
   `Protocol`, `runtime_checkable`; `liquid/ideal.py` needs `dataclass` and
   `Dict`; `liquid/davies.py` needs `numpy`, `dataclass`, `Dict` and the three
   `water_properties` names (via `_kg_per_L` for `water_density_kg_per_m3`).
5. **The mol/L → mol/kg conversion exists five times, not three.** Besides
   `ionic_strength_molal_from_molar` and the two `_kg_per_L` copies,
   `sit_liquid_model.py:167` (`gamma_all`) and `:288` (`compute_gammas`)
   compute `water_density_kg_per_m3(T_K) / 1000.0` inline, with no fallback.
   Decision 12 leaves those two and records them.
6. **Docstring rewrites the note does not list** (decision 10), in addition to
   `liquid_phase_model.py:4`, `sit_liquid_model.py:227`, `gas_eos.py:19`,
   `water_properties.py:5-6` and `framework.py:35-36`:
   `liquid_phase_model.py:45-48` ("Added by CP2 of `LAYER1_GAP_CLOSURE` per §8.4
   of `THERMODYNAMIC_MODEL_ARCHITECTURE.md`", "(this phase)"), `:58` ("§10.4 of
   `MASS_EXCHANGE_ARCHITECTURE.md`"), `:101-104` (two "(in CP2)" labels; the
   `NRTLLiquidModel` "(future)" line is reviewed), `:262` (§8.4 pointer);
   `sit_liquid_model.py:218` (§8.4 pointer), `:226` ("per the design doc").
   Reviewed, not necessarily changed (rule: text about the code stays):
   "backward-compat" wording at `sit_liquid_model.py:92,124-125,273` (check
   whether `compute_gammas` is still called as described) and `framework.py:42-47,
   54` ("Backward-compatible properties"); `gas_eos.py:239` ("Implements the
   `GasEOS` protocol"), true after checkpoint 3. The `OPEN_WORK.md`
   development-history entry counts 8 such lines in `thermo/` (by its own rule:
   lines naming a design doc under `docs/dev/` or a `CP<n>`/phase label); 7 are
   the lines above and the 8th is `equilibrium_constants.py:16`, a pointer to
   `OPEN_WORK.md` that stays. Expected after checkpoint 5: 1.
7. **Live `OPEN_WORK.md` lines the note's list misses:** the note names five
   entries; the constants entries also cite `thermo/gas_eos.py` at `:208` (the
   "Let a simulation set the value of `R`" consumer list) and `:462`, `:468`
   (the "Sweep the package for each fundamental constant" update; `:462`
   narrates a past move and gets "(now ...)", `:468` is a current file list).
   Full list: `:120,124,132` (mol/L entry, rewritten in checkpoint 4), `:208`,
   `:256` (development-history entry's example), `:370` (charge entry), `:462`,
   `:468`, `:584` (gas EOS entry), `:684` (water constants entry). `:20` and
   `:35` stay as history.
8. **Doctests to run by hand:** the module example at `gas_eos.py:19-22` and the
   class example at `:253-256` (checkpoint 5). No `>>>` example in the liquid
   files.
9. `upcoming/README.md:26,29` are this phase's own "Design discussions" entry,
   replaced at shipping (decision 16).
10. Observation, left alone: `test_partition_model.py:483`'s docstring says
    "IdealGasEOS path" about `MultispeciesVLEPartition`, which takes no EOS
    (`chemistry/partition.py:118-120`); it imports nothing from `thermo`. Loose
    wording, not a path.

### Checkpoints

- [x] 1. **Liquid group into `liquid/`.** Split `liquid_phase_model.py` into
      `liquid/protocols.py` (module docstring and lines 41-143),
      `liquid/ideal.py` (157-182) and `liquid/davies.py` (`_kg_per_L`, `_LN10`,
      185-304); move `sit_liquid_model.py` to `liquid/sit.py`, and
      `water_properties.py` and `factory.py` into `liquid/`; add a
      docstring-only `liquid/__init__.py`. New short module docstrings for
      `ideal.py` and `davies.py` (wording reviewed in checkpoint 5). Imports
      (decisions 8 and discrepancy 4): `factory.py` imports `.protocols`,
      `.ideal`, `.davies`, `.sit`. Staying files: `thermo/__init__.py` (liquid,
      SIT, factory and water blocks) and `framework.py:20`. Tests (decision 7):
      `test_nr_gas_liquid_cp2.py` (short form, with `DifferentiableLiquidModel`
      deep at `PyOMES.thermo.liquid.protocols`; 345-347 become two imports) and
      `test_liquid_phase_model.py:244` (short form). Docstrings wait for
      checkpoint 5.
      Sanity: old files gone; each definition's AST identical to the original,
      and the three split files together hold every definition of
      `liquid_phase_model.py` exactly once; fingerprint of Ideal, Davies and SIT
      (`gamma_all`, `gamma`, `jacobian_dgamma_dx`, SIT `compute_gammas`) on a
      grid of compositions and temperatures, `make_activity_model` for every
      accepted string, the water functions, and `THERMO_DAVIES`, bit-identical
      before and after; `import PyOMES.thermo.liquid_phase_model`,
      `.sit_liquid_model`, `.water_properties`, `.factory` raise
      `ModuleNotFoundError`; fresh-interpreter imports of `PyOMES`,
      `PyOMES.thermo`, `PyOMES.thermo.liquid.davies`, `PyOMES.databases`,
      `PyOMES.chemical_equilibrium` in several orders, each asserting root and
      deep names are the same objects; no bare LF; full suite **2100 passed**.
      _Notes: done 2026-09-24. Suite before: **2100 passed** (2m09s); after:
      **2100 passed**, 0 failed, 166 warnings, 2m11s. By script, with byte-level
      line handling: `liquid/protocols.py` holds the old module docstring and
      lines 41-143, `liquid/ideal.py` lines 157-182, `liquid/davies.py`
      `_kg_per_L`, `_LN10` and lines 185-304; `ideal.py` and `davies.py` have new
      one-line module docstrings and `liquid/__init__.py` is docstring-only (all
      three reviewed in checkpoint 5). Each new file imports only what it uses
      (AST check: no unused or undefined names), which drops the unused
      `Optional` (discrepancy 4). `sit_liquid_model.py` → `liquid/sit.py` and
      `water_properties.py` → `liquid/water_properties.py` by plain filesystem
      move, byte-identical; `factory.py` → `liquid/factory.py` with its two
      import lines rewritten to four (`.protocols`, `.ideal`, `.davies`, `.sit`).
      Staying files, import lines only: `thermo/__init__.py` (the liquid, SIT,
      factory and water blocks, now `.liquid.<module>`) and `framework.py:20`
      (now three lines). Tests: `test_nr_gas_liquid_cp2.py` (7 import lines;
      345-347 became a short-form `DaviesLiquidModel` line plus a deep
      `PyOMES.thermo.liquid.protocols` line for `DifferentiableLiquidModel`, and
      352 is deep; the other five are short form) and
      `test_liquid_phase_model.py:244` (short form; it now sits next to an
      existing short-form line at :243, left as two lines). Checks: the top-level
      statements of `liquid_phase_model.py`, minus imports and module docstring,
      are the same 7 ASTs as those of the three new files taken together, each
      once; `factory.py`, `__init__.py` and `framework.py` have identical ASTs
      apart from imports; a fingerprint of 13,668 values (Ideal, Davies and SIT
      `gamma`, `gamma_all`, `jacobian_dgamma_dx`, SIT `compute_gammas` with
      default and custom `epsilon`, over 11 temperatures including 1e-3 K, −5 K,
      1500 K and NaN and 5 compositions including zero and negative ones; the
      four water functions; `make_activity_model` for 16 inputs including the
      error cases; `THERMO_DAVIES`/`THERMO_IDEAL`; the SIT tables; 7 of the values
      are recorded exceptions at extreme temperatures) is bit-identical before and
      after (SHA-256 `88bed7f6...`, two runs before). `import
      PyOMES.thermo.liquid_phase_model`, `.sit_liquid_model`, `.water_properties`
      and `.factory` raise `ModuleNotFoundError`. Nine fresh-interpreter import
      orders pass (`PyOMES`, `PyOMES.thermo`, `.liquid`, `.liquid.davies`,
      `.liquid.factory`, `.liquid.sit`, `.framework`, `PyOMES.databases`,
      `PyOMES.chemical_equilibrium` first), each asserting that the 19 root names,
      the engines' `make_activity_model`/`ActivityModel` and
      `THERMO_DAVIES.liquid_activity` are the new modules' objects and that
      Davies and SIT still satisfy `DifferentiableLiquidModel`. No bare LF in any
      new or changed file; every one ends in `\r\n`. No notebook changed. The
      first run of the edit script stopped on a wrong boundary check (line 157 is
      the `@dataclass` decorator, not the class line) after creating only the
      empty `liquid/` folder, which was removed before the corrected run.
      Staging: `liquid/sit.py`, `protocols.py`, `ideal.py`, `davies.py` with
      autocrlf off (decision 17); the rest with a plain `git add`._
- [x] 2. **Gas group into `gas/`.** Split `gas_eos.py` into `gas/protocols.py`
      (`GasEOS`, 40-47, with the interface part of the module docstring,
      including the partial-pressure vs fugacity distinction), `gas/ideal.py`
      (`IdealGasEOS`, 50-63, still subclassing `GasEOS` via `.protocols`) and
      `gas/peng_robinson.py` (66-506, with the usage example and references);
      add a docstring-only `gas/__init__.py`. `..units` becomes
      `from PyOMES.units import R_L_ATM_PER_MOL_K as R` in `ideal.py` and
      `peng_robinson.py`. Staying files: `thermo/__init__.py` (gas block) and
      `framework.py:21` (`.gas.protocols`). Tests: `test_gas_eos.py:28,79,93`
      (short form); its docstring at :2 waits for checkpoint 5.
      Sanity: AST of each definition identical; fingerprint of `IdealGasEOS` and
      `PengRobinsonEOS` (`pressure_atm`, `partial_pressures_atm`, custom `kij`,
      species outside `BIOGAS_SPECIES`) on a grid of amounts, temperatures and
      volumes bit-identical; `import PyOMES.thermo.gas_eos` raises
      `ModuleNotFoundError`; fresh-interpreter import orders; no bare LF; full
      suite **2100 passed**.
      _Notes: done 2026-09-24, after checkpoint 1 was committed as `13e6139`.
      Suite before: **2100 passed** (2m08s); after: **2100 passed**, 0 failed,
      166 warnings, 2m10s. By script, with byte-level line handling:
      `gas/protocols.py` holds `GasEOS` (lines 40-47) under a new title line plus
      the old docstring's interface paragraph (lines 4-11, with the
      partial-pressure vs fugacity distinction); `gas/ideal.py` holds
      `IdealGasEOS` (50-63), still subclassing `GasEOS`, imported with
      `from .protocols import GasEOS`, under a new one-line docstring;
      `gas/peng_robinson.py` holds lines 66-506 under a new title line plus the
      old docstring's pressure-range note, usage example and references (lines
      13-28). All three keep the source file's blank line between docstring and
      `from __future__`. `..units` became `from PyOMES.units import
      R_L_ATM_PER_MOL_K as R` in `ideal.py` and `peng_robinson.py`; `gas/protocols.py`
      imports only `Dict`. `gas/__init__.py` is docstring-only. Docstring
      wording, including the usage example's old import path, waits for
      checkpoint 5. Staying files, import lines only: `thermo/__init__.py` (the
      gas block became `.gas.protocols`, `.gas.ideal` and a parenthesised
      `.gas.peng_robinson` import) and `framework.py:23`
      (`.gas.protocols`). `test_gas_eos.py:28,79,93` use the short form; its
      docstring at :2 waits for checkpoint 5. Checks: the 14 top-level statements
      of `gas_eos.py`, minus imports and module docstring, are the same ASTs as
      those of the three new files taken together, each once; `__init__.py` and
      `framework.py` differ from `HEAD` only in imports, and `test_gas_eos.py` is
      identical with imports removed at every depth (its imports are inside the
      test functions); no unused or undefined names in the new files. A 2,768-value
      fingerprint (`IdealGasEOS` and `PengRobinsonEOS` with the default, empty and
      extended `kij`, a species set extended with argon; `pressure_atm` and
      `partial_pressures_atm` over 6 temperatures, 5 volumes including 0 and
      1e-6 L, 5 amounts and 8 mixtures including unknown species, zero amounts
      and an empty dict; the `BIOGAS_*` tables; `ThermoFramework(gas_eos=...)`;
      3 recorded `AttributeError`s for `partial_pressures_atm(None)` on
      Peng-Robinson, which the ideal EOS accepts) is bit-identical before and after
      (SHA-256 `1fba3358...`, two runs before), and the checkpoint-1 liquid
      fingerprint is unchanged. `import PyOMES.thermo.gas_eos` raises
      `ModuleNotFoundError`. Nine fresh-interpreter import orders pass (`PyOMES`,
      `PyOMES.thermo`, `.gas`, `.gas.ideal`, `.gas.peng_robinson`, `.framework`,
      `.liquid`, `PyOMES.databases`, `PyOMES.chemical_equilibrium` first), each
      asserting that the root gas names are the new modules' objects,
      `IdealGasEOS` still subclasses `GasEOS` and `PengRobinsonEOS` does not.
      No bare LF in any new or changed file. No notebook changed. Staging: the
      `gas/` files with a plain `git add` (`gas_eos.py` was stored LF,
      decision 17)._
- [x] 3. **`GasEOS` as a `Protocol`** (decision 11). `@runtime_checkable class
      GasEOS(Protocol)` with the same two method signatures and `...` bodies;
      `IdealGasEOS` drops the base class. Known behaviour changes, all
      intended: `GasEOS()` now raises `TypeError` (nothing calls it);
      `isinstance`/`issubclass` of `PengRobinsonEOS` against `GasEOS` become
      true. Add one test to `test_gas_eos.py` asserting both implementations
      satisfy `GasEOS` and neither subclasses it.
      Sanity: the checkpoint-2 fingerprint bit-identical; `ThermoFramework(
      gas_eos=PengRobinsonEOS(...))` still constructs; full suite **2101
      passed**.
      _Notes: done 2026-09-24, after checkpoint 2 was committed as `97346d2`.
      Suite before: **2100 passed** (2m20s); after: **2101 passed**, 0 failed,
      166 warnings, 2m11s (the one new test). `gas/protocols.py`: `GasEOS` is now
      `@runtime_checkable class GasEOS(Protocol)` with the same two method
      signatures, `...` bodies and one-line method docstrings; its class
      docstring says implementations satisfy it structurally and that
      `partial_pressures_atm` returns partial pressures for `IdealGasEOS` and
      fugacities for `PengRobinsonEOS`. `gas/ideal.py`: `IdealGasEOS` no longer
      subclasses `GasEOS`, and the now-unused `from .protocols import GasEOS` is
      gone. `test_gas_eos.py`: new `TestGasEOSProtocol` with
      `test_both_satisfy_gas_eos_without_subclassing` (it would fail on the
      previous code on both counts). The signatures of both implementations match
      the protocol's exactly (`inspect.signature`). Checks: the 2,768-value gas
      fingerprint is bit-identical to checkpoint 2's (it includes
      `isinstance(IdealGasEOS(), GasEOS)`, still true); `GasEOS()` now raises
      `TypeError: Protocols cannot be instantiated` (nothing in the repo calls
      it); `isinstance`/`issubclass` against `GasEOS` are true for both
      implementations and false for `object()`; `IdealGasEOS.__mro__` is
      `(IdealGasEOS, object)`; `ThermoFramework(gas_eos=PengRobinsonEOS(...))`
      constructs; per-definition AST comparison against `HEAD`: `IdealGasEOS`
      identical apart from its bases, module docstrings unchanged, only `GasEOS`
      changed. No bare LF. The module docstring's "abstract interface" wording is
      reviewed in checkpoint 5._
- [x] 4. **One `water_kg_per_L`** (decision 12). Add
      `water_kg_per_L(T_K)` to `liquid/water_properties.py` (body of the old
      `_kg_per_L`); `davies.py` and `sit.py` drop `_kg_per_L` and import it;
      `ionic_strength_molal_from_molar` uses it. Add one test to
      `test_liquid_phase_model.py` (value at 298.15 K; fallback of 1.0 for a
      non-physical temperature). `OPEN_WORK.md`'s mol/L entry is closed for
      these three and records SIT's two inline copies (discrepancy 5).
      Sanity: the checkpoint-1 liquid fingerprint bit-identical, plus
      `ionic_strength_molal_from_molar` over a grid including non-finite and
      non-positive inputs and a non-physical temperature; no `_kg_per_L` left;
      full suite **2102 passed**.
      _Notes: done 2026-09-24, after checkpoint 3 was committed as `614dc53`.
      Suite before: **2101 passed** (2m15s); after: **2102 passed**, 0 failed,
      166 warnings, 4m09s (the one new test; the slower run is sandbox load, not
      the change). `liquid/water_properties.py`: new `water_kg_per_L(T_K)`, whose
      body is the old `_kg_per_L` verbatim, placed after
      `water_density_kg_per_m3`; `ionic_strength_molal_from_molar` now returns
      `float(I / water_kg_per_L(T_K))` after its unchanged `I` guard.
      `davies.py` and `sit.py` drop `_kg_per_L` and import the helper; each
      Jacobian's `dIm_dImolL` line calls it; `davies.py` no longer imports
      `water_density_kg_per_m3`, and `sit.py` still does, for its two inline
      copies. `test_liquid_phase_model.py`: new `test_water_kg_per_L` (value at
      298.15 K equals the density / 1000, about 0.997; 1.0 at 2273.15 K, where
      the density is negative, and at NaN), with a deep import, since the helper is
      not exported. The two private copies were identical to the helper, so the
      Jacobians are unchanged by construction. `ionic_strength_molal_from_molar`
      used to test the kg/L value rather than the density: the two can differ
      only when the density is positive but below about 2.5e-321 kg/m³ (dividing
      it by 1000 underflows to 0), which no representable temperature produces;
      elsewhere `I / (rho / 1000)` and `I / 1.0` are the same floats as before.
      Checks: a new 800-value fingerprint of the conversion and its users
      (density, both old private copies vs the helper, `ionic_strength_molal_from_molar`
      for 9 ionic strengths including NaN, ±inf, negative, 0 and 1e-300, and the
      Davies and SIT `gamma_all` and Jacobians, over 16 temperatures of which 8
      give a bad density: around the correlation's pole at T_C = −68.12963, 1500 K,
      2273.15 K, NaN and ±inf) is bit-identical before and after (SHA-256
      `5e9b452d...`, two runs before), and the 13,668-value checkpoint-1 liquid
      fingerprint is unchanged. Per-definition AST against `HEAD`: `_kg_per_L`
      removed from both model files, `water_kg_per_L` added,
      `ionic_strength_molal_from_molar` changed, and in `DaviesLiquidModel` and
      `SITLiquidModel` only the call name on the `dIm_dImolL` line. No unused or
      undefined names; no `_kg_per_L` left in any `.py` or notebook; no bare LF.
      `OPEN_WORK.md`: the "Three copies" entry is now "SIT converts mol/L to
      mol/kg-water inline, without the bad-density fallback". It records that the
      three copies share the helper, and that SIT's `gamma_all` and
      `compute_gammas` still divide by `water_density_kg_per_m3(T_K) / 1000.0`
      inline: at a negative density (for example 1500 K) their ε sum uses negative
      molalities while the Debye–Hückel term uses 1 kg/L (measured: at 2273.15 K
      SIT γ for 0.05 M NaCl is 0.8214 against 0.8239 at 298.15 K; NaN at NaN).
      Nothing links to the old heading. Two of this checkpoint's edits were
      rejected once in the edit prompt and then re-applied unchanged on the
      owner's go-ahead._
- [x] 5. **Docstrings, comments and live docs** (decision 10, discrepancies 6-8).
      Every path and label listed, re-derived against the new files (line
      citations re-read, not shifted). Subpackage `__init__.py` docstrings say
      what each folder holds; module docstrings of the new files reviewed;
      `framework.py`'s `gas_eos` parameter says what the field does today
      (nothing reads it). Live docs: `docs/architecture.md:397-400` (the
      `thermo/` block shows both folders), `PyOMES/README.md:22` (names the
      subfolders), `OPEN_WORK.md` at the lines in discrepancy 7 and the
      development-history count, `test_gas_eos.py:2`.
      Sanity: changed `.py` files have identical ASTs with docstrings removed;
      the two gas `>>>` examples run by hand; every changed relative link
      resolves; no bare LF; full suite **2102 passed**.
      _Notes: done 2026-09-24, after checkpoint 4 was committed as `ca98358`.
      Suite before: **2102 passed** (2m40s); after: **2102 passed**, 0 failed,
      166 warnings, 2m33s. Every claim in a rewritten docstring was checked
      against the code first. **`liquid/protocols.py`:** module title names the
      three protocols; the `gas_eos.py` path is now a cross-reference to
      `GasEOS`; the dual-protocol note says all three models satisfy both
      protocols and names the real users (`HenryEquilibrium` calls `gamma_all`,
      `interphase.py:184`; the NR solver calls the per-ion `gamma` in its outer
      ionic-strength loop; the Bisection engine calls `gamma` or SIT's
      `compute_gammas`); `DifferentiableLiquidModel` loses "CP2 of
      `LAYER1_GAP_CLOSURE`", both § pointers and "(this phase)" and says in the
      present tense that nothing in the package calls it (only tests do) and why
      the NR solver does not need it; `LiquidPhaseModel` loses both "(in CP2)"
      labels, and the `NRTLLiquidModel` "(future)" bullet became "none is
      implemented". **`davies.py`, `sit.py`:** the Jacobian docstrings lose the
      § pointer and "(per the design doc)", and point at
      `PyOMES.thermo.liquid.protocols.DifferentiableLiquidModel` (the SIT one
      named the deleted `liquid_phase_model`). **SIT `compute_gammas` was
      documented as called by the NR engine; its only caller is the Bisection
      engine's acid-base solver** (`acid_base.py:243-267`), so its docstring,
      the class docstring's bullet and the `ION_CHARGES` comment now say so and
      drop "backward-compat". **`water_properties.py`:** the stale
      `PyOMES/speciation/` sentence became "used by the Davies and SIT activity
      models in this folder" (nothing outside `thermo/` imports these functions).
      **Subpackage docstrings:** `liquid/__init__.py` and `gas/__init__.py` list
      each file in the style of `reactions/kinetic/__init__.py`; the liquid one
      says `DifferentiableLiquidModel` and `water_kg_per_L` are not exported.
      **`gas/protocols.py`:** "abstract interface" became "protocol".
      **`gas/peng_robinson.py`:** the usage example imports from `PyOMES.thermo`;
      the class says it satisfies `GasEOS` (full cross-reference) and returns
      fugacities. **`framework.py`:** `gas_eos` now says nothing in the package
      reads it (the `KineticGasLiquidLink` claim was false); "Backward-compatible
      properties" became "Derived properties", which says the Bisection and NR
      engines copy `use_activity`/`activity_model` into their own attributes
      (`bisection/engine.py:106-107`, `nr/engine.py:152-153`), and the section
      comment matches. **`test_gas_eos.py`:** module docstring names
      `PyOMES.thermo.gas`, drops "checkpoint 11, decision D4", "this one file"
      and "the checkpoint-1 audit", and lists the protocol test. **Live docs:**
      `docs/architecture.md`'s `thermo/` block now lists `framework.py`,
      `equilibrium_constants.py`, `liquid/` and `gas/` in the style of the
      `reactions/` block; `PyOMES/README.md:22` names both subfolders, with links;
      `OPEN_WORK.md` at the lines in discrepancy 7: `:209` (the `R` consumers list
      now names `gas/ideal.py` and `gas/peng_robinson.py`), `:374` (`liquid/davies.py`,
      `liquid/sit.py`), `:466` and `:472` (dated narration, "now split into
      `thermo/gas/`"), `:588` (`gas/ideal.py`), `:688`
      (`liquid/water_properties.py`); the development-history entry keeps its
      dated 2026-09-20 figures, swaps its dead `liquid_phase_model.py` example
      for `core/control_volume.py:386`, and gains an update: `thermo/` is at 1
      line (the entry's rule, recounted: only `equilibrium_constants.py:16`, a
      deliberate pointer to `OPEN_WORK.md`). `:20`, `:35` stay as history.
      **Checks:** every changed `.py` file has an identical AST to `HEAD` with
      docstrings removed (no code change; the `sit.py` comment is not in the
      AST); the two `>>>` blocks in `gas/peng_robinson.py` run with
      `doctest.testmod` (7 examples, 0 failed); the three new README link targets
      exist; no bare LF, every file ends in `\r\n`; no new line over 100
      characters (the two over it are code lines that already were). A sweep for
      the old module names, `PyOMES/speciation`, the `KineticGasLiquidLink` claim
      and "backward-compat" over live files leaves only intended hits and, outside
      `thermo/`, other packages' own "backward compatible" wording (out of scope).
      **Observation, not changed:** a markdown cell in
      `docs/tutorials/reactions/chemistry_database.ipynb` (cell 5) calls
      `use_activity`/`activity_model` "derived read-only properties for
      backward-compatible inspection"; not false, outside this checkpoint's list,
      and editing it would mean a notebook re-run, so it is left._
- [x] 6. **Sweep.** Search every old path in every form (dotted, slash,
      backslash, relative, bare filename, Sphinx cross-reference) across every
      file type, with plain `grep` and the Grep tool, matching full paths for the
      ambiguous new names (`protocols.py`, `ideal.py`, `factory.py`, `sit.py`).
      Remaining hits only in `shipped/`, `docs/dev/ideas/`, the design note, this
      checklist, gitignored paths, and `OPEN_WORK.md:20,35` (history). AST import
      audit of every `.py` and notebook code cell: every `PyOMES`/`models` import
      resolves, including each imported name. Each old path fails to import.
      Sanity: full suite **2102 passed**.
      _Notes: done 2026-09-24, after checkpoint 5 was committed as `56a1e43`. The
      sweep found no missed site, so the only file changed besides this checklist
      is `OPEN_WORK.md` (one new entry, below). Full suite **2102 passed**, 0
      failed, 166 warnings, 2m07s. **Old-path search:** plain `grep -rnE` over the
      whole tree (gitignored and hidden paths included) and every file type, for
      `thermo.`/`thermo/`/`thermo\` followed by any old module or by a new module
      name without its folder (`protocols`, `ideal`, `davies`, `sit`,
      `peng_robinson`, `water_properties`, `factory`), relative
      `from .<old module>`, and bare `liquid_phase_model`, `sit_liquid_model`,
      `gas_eos` (with or without `.py`, not inside a longer name): 138 lines in
      16 files; 31 in `shipped/`, 15 in `docs/dev/ideas/`, 81 in this phase's note
      and checklist, 3 in gitignored `PyOMES.egg-info/SOURCES.txt`, and 8 live,
      all intended (the `gas_eos` field name in `framework.py:36,54` and
      `architecture.md:398`; `OPEN_WORK.md:20,35` (history) and `:466,472`
      (dated narration with "now `thermo/gas/`"); this phase's own entry at
      `upcoming/README.md:29`, replaced at shipping). No relative import of a
      moved module is left in `thermo/`'s top level. The Grep tool finds 135 lines
      in 15 files, exactly the same minus the 3 gitignored ones. (A first run
      used `\w` inside a bracket expression, which extended `grep` does not
      support, so it also matched `test_liquid_phase_model.py`; rerun with
      `[:alnum:]_`.) **AST import audit:** every `.py` file and notebook code
      cell, including lazy, `TYPE_CHECKING` and docstring/template import lines:
      3,480 `PyOMES`/`models` imports in 180 files across 93 modules, each module
      imported and each name resolved. 23 problems, **none involving `thermo`**:
      16 in gitignored `scratch/` notebooks and synthetic strings in
      `test_package_layering.py`, as in the reactions phase; the known
      `tests/run_tests.py:105`, `test_simulation.py:2186` and
      `batch_fermenter.ipynb` cell 10 (does not parse; logged); and four
      docstring examples with pre-rename paths (`adm1/bsm2.py:24`,
      `adm1/bsm2_direct.py:20` (logged), `hplc/column.py:48`, twice). This audit
      reads docstring examples, which the reactions phase's did not, so the last
      group is new: logged in `OPEN_WORK.md` as "Docstrings in `models/` and
      `numerics/` still use pre-rename module paths", together with
      `numerics/spatial.py:55` (a `fermenter.models...` cross-reference found by
      the follow-up search). The first audit run missed `bsm2.py:24` (a bug in
      the script's handling of multi-line template imports); fixed and rerun
      before the numbers above. **Old paths:** `import PyOMES.thermo.<m>` raises
      `ModuleNotFoundError` for all five old modules and for the five new module
      names used without their folder (`davies`, `ideal`, `protocols`, `sit`,
      `peng_robinson`); `from PyOMES.thermo import` `liquid_phase_model`,
      `gas_eos`, `DifferentiableLiquidModel` and `water_kg_per_L` raise
      `ImportError` (the last two are not exported, decisions 12-13).
      **Notebooks:** the branch changes no `.ipynb`, `.toml`, `.yml` or `.json`
      file (`git diff --name-only main..HEAD`: 17 `.py`, 4 `.md`), so none was
      re-run._
- [ ] 7. **Seam guard** (decision 14). Add `_thermo_layout_violations()` and two
      tests to `test_package_layering.py`: a synthetic self-check (crossings in
      every import form, allowed same-folder imports, files outside the rule's
      scope) and the real-files test (after asserting the new modules are found).
      Update the module docstring. Also inject, in memory, a
      `liquid -> gas` import, a lazy `framework` import in `gas/`, and a
      `TYPE_CHECKING` package-root import in `liquid/`, to see each fail.
      Sanity: all layering tests pass; full suite **2104 passed**.

## Shipping

- [ ] Full test suite green on the branch
- [ ] `git checkout main`
- [ ] `git merge --no-ff thermo-subfolder-structure -m "Merge thermo-subfolder-structure: group thermo/ into liquid/ and gas/ subfolders"`
- [ ] `git tag thermo-subfolder-structure-shipped` on the merge commit
- [ ] `git push` and `git push --tags` (separate commands)
- [ ] `git branch -d thermo-subfolder-structure` and
      `git push origin --delete thermo-subfolder-structure`
- [ ] Full suite on `main` after the merge
- [ ] Move the design note and this checklist to `docs/dev/implementation/shipped/`
      (plain filesystem move); add "Shipped" banners; repoint the moved docs' own
      relative links
- [ ] Update `upcoming/README.md`: remove the "Design discussions" entry, add a
      "Recently shipped" entry
- [ ] Repoint any live link to the moved docs (search for
      `THERMO_SUBFOLDER_STRUCTURE` across `docs/` and `OPEN_WORK.md`)
