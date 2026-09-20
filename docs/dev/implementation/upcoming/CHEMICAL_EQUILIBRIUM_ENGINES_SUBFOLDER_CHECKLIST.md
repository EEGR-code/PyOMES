
# Phase Kickoff Checklist — chemical-equilibrium-engines-subfolder

> Checklist for [`CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER.md`](CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER.md),
> which is the source of truth for goals, audit, target layout and Decisions
> 1-10. All open questions there are resolved; do not re-litigate them. See
> `README.md`'s "Branching and tagging convention" for the branch/merge/tag
> rationale. Modelled on
> [`../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md`](../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md).

**Working rules for this phase**

- Pure refactor: no behaviour or numerical change. Anything that looks like a
  bug or improvement is noted, not fixed. **Exception: Part D (checkpoint 14)
  deliberately changes numbers** (unifying the gas constant); it runs last and
  is verified by a recorded before/after shift, not by bit-identity.
- The repo owner runs every `git add`, `git commit` and `git push`. At each
  checkpoint: edit, run the checkpoint's tests, report, then stop for review
  and commit.
- No shims at old paths, no deprecation periods (Decisions 2, 9, 10).
- Before deleting or moving anything, re-verify the audit claim with a fresh
  search over `.py`, `.ipynb` and `.md` (excluding the gitignored `scratch/`).
  If the plan doc is wrong, stop and report.
- Leave `docs/dev/implementation/shipped/` and `docs/dev/ideas/` untouched,
  including their stale `src/speciation/` paths.
- Use `git mv` for every move in Part C.

## Pre-flight

- [x] `git status -sb` clean (no stray uncommitted work left over from a
      previous task) — only the untracked plan doc
      `CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER.md`, which is committed together
      with this checklist as the branch's first commit
- [x] `git log origin/main..main --oneline` empty (nothing unpushed sitting
      around from earlier work) — verified 2026-09-20
- [ ] Branch created off current `main`:
      `git checkout -b chemical-equilibrium-engines-subfolder`
- [ ] This checklist file (and the plan doc) committed on that branch as the
      first commit — so partial work is never silent

## During

- [x] Plan/design doc exists in `docs/dev/implementation/upcoming/`
      (`CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER.md`)
- [ ] Checkpoints tracked below as they land, one commit per checkpoint
- [ ] **If work stalls or is paused before shipping:** add a status banner
      to the top of the plan doc *immediately* — what's built, what's
      tested, why it stopped, which branch/commit it's on. Don't leave this
      for a future audit to discover.

### Checkpoints

- [x] 1. Inventory every import site (`.py` and `.ipynb`, plus the two
      notebook generators' string templates and lazy in-function imports) of
      the modules being moved or deleted. Record counts in "Checkpoint 1
      inventory" below. Baseline: full test suite green.

**Part A — thermo consolidation** (retire `activity_models.py` / `sit.py`)

- [x] 2. Add `ActivityModel` and `make_activity_model` to `thermo/` (in
      `liquid_phase_model.py` or a small `thermo/factory.py` — decide here);
      export from `PyOMES.thermo`. Fix the stale `PyOMES/speciation/`
      docstring paths in `liquid_phase_model.py`. `activity_models.py` still
      exists after this checkpoint.
      _Decided: `ActivityModel` protocol in `thermo/liquid_phase_model.py`
      (next to its sibling `LiquidPhaseModel`); `make_activity_model` in new
      `thermo/factory.py` (top-level imports, no import cycle). Old-vs-new
      factory checked equal on 10 valid + 3 error inputs;
      `test_liquid_phase_model`, `test_thermo_framework`, `test_speciation`,
      `test_nr_tableau_gas_liquid`: 116 passed. Left for you:
      `thermo/water_properties.py:6` also has a stale `PyOMES/speciation/`
      mention (outside the plan's scope)._
- [x] 3. Repoint importers at `PyOMES.thermo`: `engine.py`, `nr_engine.py`,
      `acid_base.py`, `activity.py` (`thermo.water_properties`),
      `test_nr_tableau_gas_liquid.py`, and the `DaviesActivityModel` users in
      `tests/validation/speciation/` (tests, generator, 3 notebooks →
      `DaviesLiquidModel`). Remove `IdealActivityModel` /
      `DaviesActivityModel` from `__init__.py` and `__all__` (Decision 9).
      _Notes: package files use the relative form `..thermo` (matches the
      neighbouring `..units` / `..core` imports). Notebooks 03/04/05 were
      edited as raw JSON (two lines each; the Edit tool refuses `.ipynb` and
      `nbformat` is not installed), then their code cells were executed in a
      scratch dir: all three run and reproduce the saved values. Two comments
      in `acid_base.py` that said "DaviesActivityModel" now say
      "DaviesLiquidModel". Still mentioning the old name in prose only, left
      for C12 (decide then): `06_phreeqc_benchmark.ipynb` (2 places) and
      generator line 1543._
- [x] 4. Delete `sit.py` and `activity_models.py`; delete the two alias tests
      in `test_liquid_phase_model.py`; repoint the factory test there. Delete
      `davies_log10_gamma` / `davies_gamma` and fix `activity.py`'s docstring
      (Decision 6). Run `test_liquid_phase_model.py`,
      `test_thermo_framework.py`, `test_speciation.py`,
      `tests/validation/speciation/`.
      _Notes: of the four shim-only tests in `test_liquid_phase_model.py`, the
      two alias tests were deleted; the water-helpers test was repointed to
      `PyOMES.thermo` (not deleted, so the `debye_huckel_A(298.15) ≈ 0.509`
      assertion keeps a home), as was the SIT-epsilon test; both factory tests
      repointed to `PyOMES.thermo`. `activity.py` also lost its now-unused
      `debye_huckel_A` / `ionic_strength_molal_from_molar` import. Suite goes
      from 2066 to 2064 tests (two alias tests deleted)._

**Part B — remove `strong_ions.py`**

- [x] 5. Search `.ipynb` files (earlier audit covered `.py` and `.md` only)
      for `strong_ions_from_feed_molL` / `PyOMES.chemical_equilibrium.strong_ions`.
      If clear: delete `strong_ions.py` and
      `tests/standalone/test_strong_ions.py`; remove the `__init__.py` import
      and `__all__` entry (Decision 10); drop the stale comment in
      `chemistry/registry.py`. `FeedState` and its fixtures stay.
      _Notes: `.ipynb` search clean (zero hits). Also removed
      `test_strong_ions.py` from the test-coverage table in the root
      `README.md` (line 149), which the plan did not list. Suite: 2064 →
      2057 (the deleted file held 7 tests). The `rich_feed` fixture in
      `tests/standalone/conftest.py` was used only by the deleted test and is
      now unused (whole-repo search, all file types, no dynamic lookups); the
      same goes for the hand-copied `_rich_feed` / `"rich_feed"` entry in
      `tests/run_tests.py` (a non-pytest runner that is already broken, see
      `OPEN_WORK.md`). Both left in place per the plan ("fixtures stay") —
      candidate for a later cleanup. `docs/architecture.md:394` still lists `strong_ions.py`
      (C12)._
- [x] 6. Update the "Adjacent, out of scope" section of
      `STRONG_ION_INFERENCE_GENERALIZATION.md` (`strong_ions.py` removed;
      `SALT_DISSOCIATION_MAP` now has no consumer). Run `test_feed_state.py`
      and the speciation tests.
      _Notes: docs-only; added a dated update paragraph at the end of that
      section. `test_feed_state.py`, `test_speciation.py`,
      `test_speciation_protocols.py`, `tests/validation/speciation/`: 112
      passed. The same doc's `nr_engine.py` / `nr_solver.py` / `engine.py`
      path citations are left for C12._

**Part C — `engines/` subfolder**

- [x] 7. De-duplicate `acid_base._vant_hoff_K` / `nr_tableau._vant_hoff_log_K`
      into one helper in `PyOMES/thermo/` (Decision 3; filename decided here,
      e.g. `thermo/equilibrium_constants.py`). Own commit, own tests. Repoint
      `acid_base.py` and `nr_tableau.py`.
      _Notes: new `thermo/equilibrium_constants.py` holds the correction once
      (`vant_hoff_delta_ln_K`) plus two thin wrappers, `vant_hoff_K` and
      `vant_hoff_log_K`, each keeping the edge-case rules of the original it
      replaces (they differ: the K-space one also guards K_ref <= 0 and
      T_K <= 0). Not exported from `PyOMES.thermo.__init__`. `nr_engine.py`
      now imports the helper directly (it used to import the private name from
      `nr_tableau`). Removed as unused: `_R_J_PER_MOLK` (acid_base),
      `_R_J_MOL_K`, `_LOG10_E` and `import numpy as np` (nr_tableau).
      **Bit-identical, verified:** an 814-value fingerprint of
      `solve_acid_base` (all 8 enthalpies, 7 temperatures) and the NR engine
      (water/carbonate, plus calcite Ksp and a Henry row; 5 temperatures) was
      byte-identical before and after (same SHA-256). New
      `tests/standalone/test_equilibrium_constants.py` (23 tests) pins the
      helpers to verbatim copies of the two originals with exact equality.
      **Not folded in:** `reactions/equilibrium.py:vant_hoff_log_K` (the third
      copy) uses a log10(e) constant 1 ulp different from the one `nr_tableau`
      used, so merging it would change ~20% of corrected values by up to ~4e-15
      in log10 K; logged in `OPEN_WORK.md` instead._
- [x] 8. `git mv` the NR files into `engines/nr/` (`engine.py`, `tableau.py`,
      `solver.py`, `__init__.py`); rewrite relative imports (`..units` →
      `...units`, `..core.phases` → `...core.phases`, ...). Run `test_nr_*`
      and `test_equilibrium_classification.py`.
      _Notes: three `git mv` renames (`nr_engine.py` → `engines/nr/engine.py`,
      `nr_tableau.py` → `tableau.py`, `nr_solver.py` → `solver.py`); new
      `engines/__init__.py` and `engines/nr/__init__.py` are docstring-only
      (no re-exports, so no import side effects; revisit at C11).
      **Plan correction:** the plan said `..units` → `...units`, but `engines/nr/`
      is two levels deeper, so it is `....units` (four dots); `.protocols` →
      `...protocols`. The flat `engines/phreeqc.py` (C10) will be three dots.
      Plan doc fixed. 10 relative imports rewritten inside the moved files
      (incl. lazy `engine.py:239`, `tableau.py:408`). Old dotted paths
      rewritten in 31 tracked `.py`/`.ipynb` files (89 replacements: 69
      `nr_engine`, 16 `nr_tableau`, 4 `nr_solver`), covering tests, both
      generators' string templates, 8 validation notebooks, 2 tutorial
      notebooks, `reactions/reaction_system.py`, and docstring cross-references
      in the package. Bare filename mentions in prose (`nr_solver.py` in
      comments, e.g. `engines/nr/engine.py:86`, `tableau.py:113`,
      `protocols.py:73`, `thermo/liquid_phase_model.py:53`,
      `adm1/base.py:1048`) and current `.md` docs are left for C12. Verified:
      814-value fingerprint byte-identical before/after (same SHA-256); old
      path now raises `ModuleNotFoundError`; 10 changed notebooks run;
      generators compile; full suite 2080 passed (unchanged)._
- [x] 9. `git mv` the Bisection files into `engines/bisection/` (`engine.py`,
      `acid_base.py`, `__init__.py`). Run `test_speciation*.py` and
      `test_bisection_chemical_equilibrium_engine_alias.py`.
      _Notes: two `git mv` renames plus a docstring-only
      `engines/bisection/__init__.py`. Inside the moved files: `.acid_base`
      stays (sibling); `..thermo` → `....thermo`; `.protocols` → `...protocols`;
      `.activity` → `...activity` (acid_base, incl. its lazy import); lazy
      `..chemistry.equilibria` / `..reactions.equilibrium` → four dots
      (engine). Staying files repointed: `api.py`, `factory.py`, package
      `__init__.py` (import + the `solve_acid_base` re-export). Old dotted paths
      rewritten in 23 tracked `.py`/`.ipynb` files (38 replacements: 34
      `engine`, 4 `acid_base`): tests, `reactions/reaction_system.py`,
      `adm1/{base,bsm2}.py`, generator template, 3 notebooks, docstrings. The
      `unittest.mock.patch.object(_eng_mod, "solve_acid_base")` in
      `test_speciation.py` still works (patches the name in the moved module).
      Verified: new 271-value fingerprint of the Bisection path (engine class at
      4 temperatures ± activity, `get_CO2aq_from_totals`, the legacy
      `ChemicalEquilibriumEngine` alias, `SpeciationFactory` + adapter in both
      policies, direct `solve_acid_base`, and "package re-export is the same
      callable") byte-identical before/after; the checkpoint-8 NR fingerprint
      still byte-identical; old paths raise `ModuleNotFoundError`; 3 changed
      notebooks run (from the scratch dir); generator compiles; full suite 2080
      passed (unchanged). Housekeeping: a stray untracked `figures/` folder
      created by my checkpoint-8 notebook run (cwd was the repo root) was
      found and removed; later notebook runs use the scratch dir._
- [x] 9b. _Added during Part C (plan Decision 16)._ Delete the orphaned,
      Bisection-only `api.py` (`SpeciationEngineAdapter`) and `factory.py`
      (`SpeciationFactory`) and their two package-level exports; reword the two
      `protocols.py` docstrings that cited them (`n_iter` note and
      `EquilibriumResult.to_dict`). Log knock-on effects outside this phase.
      _Notes: prompted by the question of whether `factory.py` belonged in
      `engines/bisection/`; a fresh search (`.py`, `.ipynb`, `.md`, config)
      found zero callers, and the package has no outside users yet, so
      deleting was chosen over moving or generalising (engine selection
      already exists in `ReactionSystem.engine`; the adapter's typed API is
      Bisection-shaped). **Restore point: commit `2e5554a`**
      (`git show 2e5554a:PyOMES/chemical_equilibrium/api.py`, likewise
      `factory.py`). Top level of `chemical_equilibrium/` is now `__init__`,
      `protocols`, `activity`, `activity_dispatch`, `numerical_gradient` and
      `phreeqc_engine` (leaves in C10): everything left is engine-agnostic.
      Verified: Bisection fingerprint minus the removed adapter lines
      (265 values) byte-identical to the checkpoint-9 baseline; NR fingerprint
      byte-identical; `SpeciationFactory`, `SpeciationEngineAdapter` and both
      module paths now fail to import while every other package-level name
      still imports; full suite 2080 passed (unchanged: neither file had
      tests).
      **Effects outside this phase's scope (logged, not touched):**
      (1) `chemistry/types.py`: `AqueousEquilibrium` is now unused;
      `AqueousTotalsUser.to_engine()` has no caller; `AqueousTotals` and
      `AqueousTotalsUser` are reachable only via `SolutionRecipe.to_totals_user()`
      / `to_totals()` (`chemistry/recipe.py`), which nothing outside `recipe.py`
      calls; all three still exported from `PyOMES.chemistry`.
      (2) Stale wording: `chemistry/recipe.py:5` ("standalone speciation
      interface") and `chemistry/registry.py:18` (lists `AqueousTotalsUser` as
      a consumer). (3) `EquilibriumResult.to_dict()` was added for the
      adapter's `raw=` field; it stays (harmless, tested once in
      `test_speciation_protocols.py`) but now has no in-package caller.
      (4) `docs/architecture.md:398` still lists `api.py, factory.py` (C12).
      (1)–(2) are recorded in the recipe-layer section of
      `STRONG_ION_INFERENCE_GENERALIZATION.md`, next to the
      `SALT_DISSOCIATION_MAP` question they now join._
- [x] 10. `git mv phreeqc_engine.py engines/phreeqc.py`. Run
      `tests/validation/speciation/`.
      _Notes: one `git mv` rename. The file's only relative import was
      `.protocols` → `..protocols` (three dots would apply to a parent-package
      import such as `..units`; the file has none — `phreeqpython` is an
      absolute lazy import inside `__init__`, unchanged). Old dotted path
      rewritten in 8 files (8 replacements, exactly the checkpoint-1 inventory):
      `bisection/engine.py` and `phreeqc.py` docstrings, `test_phreeqc_engine.py`,
      `test_phreeqc_nr_agreement.py`, the ArXiv generator template,
      notebooks `ArXiv/01`, `validation/10` and `protocols/03`. **Silent-skip
      trap checked:** ArXiv `01` and validation `10` wrap the import in
      `try/except ImportError` and would have *silently skipped* their PHREEQC
      cells on a wrong path; confirmed `_HAVE_PHREEQC = True` and
      `HAS_PHREEQC = True` with the engine's `__module__` now
      `...engines.phreeqc`. Verified: new 182-value PHREEQC fingerprint
      (real `phreeqpython 1.6.2` solves, 3 temperatures × 2 systems × 2 calls
      for the cache path, plus name-translation functions) byte-identical
      before/after; Bisection and NR fingerprints still byte-identical; old
      path raises `ModuleNotFoundError`; full suite 2080 passed (unchanged)._
- [x] 10b. _Added during Part C; unrelated to the engine layout._ Make the
      PHREEQC availability flag consistent and make the ArXiv guard actually
      work. _Notes: found while checking checkpoint 10's silent-skip trap.
      Two notebooks guarded the same import with different names —
      `_HAVE_PHREEQC` (ArXiv `01_predict_ph_simple_liquid.ipynb` and its
      `_generate_notebooks.py`, 9 uses each, from 2026-09-14) and `HAS_PHREEQC`
      (validation `10_engine_protocol_hierarchy.ipynb`, 3 uses, hand-written,
      no generator) — with no repo convention behind either (no other
      availability flags exist; the package raises `ImportError` at the point
      of use and pytest uses `importorskip`). The names also hid a real
      difference: `engines/phreeqc.py` imports `phreeqpython` lazily inside
      `__init__`, so the ArXiv guard (which only imported the engine module)
      reported `True` even when `phreeqpython` was missing, and the PHREEQC
      cells then crashed instead of being skipped as the notebook's own text
      promises. Pre-existing; not caused by this phase.
      **Fix:** standardised on `HAS_PHREEQC` (PEP 8 constant style; matches the
      correct guard) by renaming the 9 uses in the ArXiv generator template and
      the 9 in the notebook, and added notebook 10's `import phreeqpython`
      presence check to the ArXiv guard. Raw byte edits with exact-count
      assertions (generator is CRLF, notebook JSON `\n` escapes handled);
      validation notebook 10 untouched. **Verified:** before the fix, with
      `phreeqpython` simulated absent, the notebook failed at code cell 7 with
      the flag `True`; after, it completes all 13 cells with the flag `False`,
      and with `phreeqpython` present the flag is `True` and it completes as
      before; the generator template and the notebook's guard cell are
      identical; all 14 changed notebook lines are `source` lines (no saved
      output touched); generator compiles; no `_HAVE_PHREEQC` left in
      `.py`/`.ipynb`; full suite 2080 passed (unchanged). The historical
      mention of `_HAVE_PHREEQC` in checkpoint 10's notes above is left as the
      record of what was true then._
- [x] 11. Update the package `__init__.py` re-exports, then external callers:
      `PyOMES/reactions/reaction_system.py`, `models/vlmodels/adm1/{base,bsm2}.py`,
      both notebook generators, notebooks, tests. Finish with a repo-wide
      search (`.py`, `.ipynb`, `.md`) for each old module path to confirm
      nothing still points at it.
      _Notes: the external callers were repointed checkpoint by checkpoint
      (C8–C10), so this was mainly verification. **Re-exports:** the package
      `__init__.py` is unchanged from C9b (`BisectionChemicalEquilibriumEngine`,
      `ChemicalEquilibriumEngine`, `EquilibriumResult`,
      `ionic_strength_from_speciation`, `warn_if_high_ionic_strength`,
      `solve_acid_base`, all resolving to the new locations). **Decision:**
      `engines/__init__.py`, `engines/nr/__init__.py` and
      `engines/bisection/__init__.py` stay docstring-only — the NR and PHREEQC
      engines were not exported from the package root before the move either,
      so re-exporting them would be new API, and it is trivial to add later.
      (Observation, not acted on: the root exports only the Bisection engine, so
      NR — the newer engine — needs a deep import; whether the root should
      export all three is a separate API question.) Fixed a stale docstring I
      had written in `engines/__init__.py` (it still listed `api` and `factory`,
      deleted in 9b, and omitted `activity_dispatch`).
      **Verification:** (1) packaging — `setup.py`'s
      `find_packages(include=["PyOMES", "PyOMES.*"])` discovers
      `engines`, `engines.nr` and `engines.bisection`, and their `__init__.py`
      files are tracked, so CI's `pip install -e .` will see them; (2) sweep of
      every old module path in every form (dotted, slash, relative, bare
      filename) over tracked `.py`/`.ipynb`/`.md`: **zero old dotted, slash or
      relative-import paths in code, notebooks or generators**; the only code
      hits are 5 prose comments naming `nr_solver.py` / `nr_engine.py`
      (`engines/nr/engine.py:86`, `engines/nr/tableau.py:113`,
      `protocols.py:73`, `thermo/liquid_phase_model.py:53`, `adm1/base.py:1048`)
      — C12; current `.md` docs still to fix in C12: `docs/architecture.md:394`,
      `OPEN_WORK.md:57,59,164`, `upcoming/README.md:60`,
      `STRONG_ION_INFERENCE_GENERALIZATION.md` (~20 refs),
      `NR_PRECIPITATION_CV_INTEGRATION.md:144`; historical docs untouched; (3)
      **AST import audit** (script, not grep): parsed every tracked `.py` and
      every notebook code cell, resolved 2,065 `PyOMES`/`models` imports across
      229 files including lazy in-function and relative ones, and checked that
      each module exists and each imported name resolves — **no problems
      attributable to this phase**. It flagged only two pre-existing items:
      `tests/run_tests.py:107` (`create_standalone_fermenter`, already in
      `OPEN_WORK.md`) and `test_simulation.py:2186`, a test that asserts a
      deleted module stays deleted; (4) no dynamic/string imports, `patch()` or
      `sys.modules` targets mention the old paths; (5) full suite 2080 passed
      (unchanged)._
- [x] 12. Update docstring cross-references (~15 files using
      `PyOMES.chemical_equilibrium.engine...`-style paths) and current docs,
      including `OPEN_WORK.md` (cites `activity_models.make_activity_model`) and
      `docs/architecture.md` (lists `activity_models.py`, `sit.py`). Leave
      `shipped/` and `docs/dev/ideas/` alone.
      _Notes: the dotted-path docstring references (~89 + 38 + 8 replacements)
      were already rewritten with the moves (C8–C10), so this checkpoint handled
      the remaining prose and the current docs. **Code comments (5):**
      `engines/nr/engine.py:86`, `engines/nr/tableau.py:113`,
      `protocols.py:72-73`, `thermo/liquid_phase_model.py:53`,
      `models/vlmodels/adm1/base.py:1048` now name the real files. **Tutorial
      notebook:** two markdown cells in `01_bisection_engine_basics.ipynb` said
      bare `engine.py` (now ambiguous: there is one under `bisection/` and one
      under `nr/`); edited as exact-count raw bytes, markdown only, still valid
      JSON. **Docs:** `docs/architecture.md` — replaced the stale `speciation/`
      block (which described files such as `chemistry_level1.py` and
      `legacy_adapter.py` that no longer exist, plus the deleted
      `activity_models.py`/`sit.py`/`strong_ions.py`/`api.py`/`factory.py`)
      with the real `chemical_equilibrium/` tree including `engines/`, and
      added a `thermo/` line so a reader can find where the activity models
      went (a small addition: `thermo/` was not in that tree at all); the
      separate stale `equilibria/`/`sim/`/`solvers/` "CUFermenter island"
      entries are untouched (an `OPEN_WORK.md` item, whose entry now records
      that this block is fixed). `OPEN_WORK.md` — the `use_activity` entry now
      cites `PyOMES.thermo.make_activity_model` and
      `engines/{bisection,nr}/engine.py` (and notes the deleted `factory.py`);
      the constants entry cites `engines/nr/solver.py`.
      `upcoming/README.md` and `NR_PRECIPITATION_CV_INTEGRATION.md` — paths.
      **`STRONG_ION_INFERENCE_GENERALIZATION.md`:** every file-and-line
      citation (about 20 references) was re-derived against the moved code, not
      just renamed: the NR-engine ones moved down one line (checkpoint 7 added
      an import), the solver and Bisection-engine ones did not; a script then
      verified all 9 cited ranges contain what the doc says (including that the
      two solver `_STRONG_CHARGES` copies really sit inside `_ionic_strength`
      and `solve_nr`); a short dated note at the top records the re-derivation.
      **Deliberately left:** `shipped/` and `docs/dev/ideas/` (historical);
      `test_accuracy_monitor.py:560` ("The legacy engine.py RuntimeWarning was
      replaced…", a narrative of a past phase); `OPEN_WORK.md:21` and
      `architecture.md` `equilibria/` lines (long-deleted island); my own
      dated notes that name removed files. Verified: changed code files
      compile; full suite 2080 passed (unchanged)._
- [x] 12b. _Added during Part C (requested after checkpoint 12)._ **Documentation
      sweep and rewrite of `PyOMES/chemical_equilibrium/`**: remove
      development-history references (phase and checkpoint labels, pointers to
      old design documents, "this phase" wording) from docstrings and comments,
      unless there is a good reason for one to remain. Docs should describe what
      the code does now and why, in the present tense.
      **Scope:** the 9 `.py` files under `PyOMES/chemical_equilibrium/`. A
      survey (2026-09-20, by script, classifying each hit by token type)
      found **83 distinct lines** across 9 files, in 7 categories:
      (A) 33 lines naming a design doc — `LAYER1_GAP_CLOSURE` ×20,
      `MASS_EXCHANGE_ARCHITECTURE` ×3, `MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN` ×3,
      `PARTITION_MODEL` ×2, `NR_PRECIPITATION_CV_INTEGRATION` ×2,
      `EQUILIBRIUM_CONSTRAINT_UNIFICATION` ×2, `CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE` ×1;
      (B) 33 checkpoint labels (`CP1`/`CP2`…, `C4d`, "Stage 16");
      (C) 30 phase/era names (`chemistry-unification-3b`, "Phase 2", "Layer 1");
      (D) 1 retired "Level 1 foundation" (`acid_base.py` header);
      (E) 7 "design doc" / `§` pointers; (F) 17 changelog-style phrases
      ("`--- NEW:`", "was deleted in", "Fixed since", "this phase's folding");
      (G) 7 old-project-name mentions (`vlsim`, all in `engines/phreeqc.py`).
      By kind: 91 docstring hits, 23 comments, **14 runtime strings** (7
      distinct error/warning messages). Heaviest files: `engines/nr/tableau.py`
      (29 lines), `engines/nr/engine.py` (16), `engines/nr/solver.py` (10),
      `engines/bisection/acid_base.py` (8), `phreeqc.py` (7).
      **Rules.** (1) Keep the *reasoning*, drop the *label*: where a comment
      records a constraint that is still true and non-obvious (for example why
      `retain_jacobian=True` is unsupported with folded gas-liquid rows, or why
      the tableau has one master per component), rewrite it as a present-tense
      statement of that constraint. (2) A pointer to a design document may stay
      **only** if that document is the canonical, still-current explanation, is
      too long to summarise in two or three lines, and lives at a stable path
      (`docs/dev/ideas/` or `shipped/`, never `upcoming/`, which moves on ship);
      each retained pointer is listed in this entry's notes with its reason.
      (3) Text that is *about* the code (a class named after a protocol, the word
      "legacy" for a code path that genuinely exists today) is not history and is
      reviewed, not removed. (4) No code changes. **Runtime strings** (the 7
      messages, e.g. the `NotImplementedError` in `engines/nr/engine.py` and the
      `ConfigurationError` texts in `engines/nr/tableau.py`) may be reworded to
      state the limitation directly, but must keep the fragments tests match on:
      **`V_liq_L`** (3 tests in `test_nr_tableau_gas_liquid.py`) and
      **`bridge`** (`ConfigurationError` test) — and each rewritten message is
      listed in the notes.
      **Decision needed at kickoff:** `engines/phreeqc.py` calls PyOMES's
      former name `vlsim` in its docs and in a public function,
      `phreeqc_to_vlsim` (also the default `species_map`, imported by
      `test_phreeqc_engine.py`). Rewording the docs is in scope; *renaming the
      function* (e.g. to `phreeqc_to_pyomes`) is a code change. The package has
      no outside users, so a rename without an alias is cheap (one module, one
      test file, the protocol tutorial notebook); default if not decided:
      docs only, function name kept, rename logged in `OPEN_WORK.md`.
      **Verification (planned):** an AST comparison per file with docstrings
      removed — the code must be identical except the explicitly listed
      runtime-string edits; the full suite; a re-run of the survey script
      (target: 0 lines, or only the listed retained pointers); the
      `V_liq_L` / `bridge` tests by name.
      **Out of scope, logged separately:** the same style of reference exists
      elsewhere — 108 lines in 32 files under the rest of `PyOMES/` (`core/` 67,
      `chemistry/` 14, `thermo/` 8, `templates/` 6, `reactions/` 5,
      `control/` 4, other 4) — and in tests (module docstrings such as
      "Tests for CP2 of LAYER1_GAP_CLOSURE", and files named after checkpoints:
      `test_nr_gas_liquid_cp2.py`, `test_raoult_h2o_fold_cp4.py`,
      `test_precipitation_gas_liquid_cp5.py`). See `OPEN_WORK.md`.
      **RESULT (done).** The survey went from **83 lines to 4**, and all 4 are
      justified: two are ordinary English ("is used to convert", in
      `solver.py` and `tableau.py`) and two name the real function
      `phreeqc_to_vlsim`. **No design-document pointer was retained**: every one
      failed the "canonical, current, too long to summarise, stable path" test
      (`LAYER1_GAP_CLOSURE`, `MASS_EXCHANGE_ARCHITECTURE`,
      `MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN`,
      `NR_PRECIPITATION_CV_INTEGRATION`, `EQUILIBRIUM_CONSTRAINT_UNIFICATION`,
      `PARTITION_MODEL` and `CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE` are either
      shipped/historical or live in `upcoming/`), and each surrounding statement
      was already self-explanatory once the label was dropped. Beyond the
      planned categories, a second, looser scan found more history-flavoured
      wording my patterns had missed: three **"CHANGE (Update N)"** changelog
      headers in `acid_base.py`, "behave exactly as before", "(legacy)/(new)"
      labels, "separated for future extensions", "it has always written back",
      "deferred" (now "is not performed"), and `protocols.py` calling the
      Bisection solver "legacy" although it is `ReactionSystem`'s default solver.
      **Statements that were false, not just stale, and are now corrected:**
      `engines/nr/engine.py`'s module "Scope" said precipitation was "deferred
      (see design doc)" although the engine has a full precipitation loop and
      gas-liquid folding; `activity_dispatch.py` promised that a later phase
      "wires this dispatch into `build_tableau()`" (that phase shipped without
      doing so — the function is called only by its own tests); `acid_base.py`
      promised the `_HA`/`_A-` fallback "is removed in the PARTITION_MODEL phase"
      (it shipped and the fallback is still live for the BSM2 VFA rows) and
      described a `_CANONICAL_NAMES` path that no longer exists;
      `tableau.py` said a component's total would read the gas phase "once CP2
      wires" it (it already does); `phreeqc.py` had a comment describing the
      opposite of what the next line does ("convert PHREEQC keys back").
      **Reasoning kept, label dropped**, for example: why `retain_jacobian=True`
      is unsupported with folded gas-liquid rows; why `SecondaryEntry.c_key`
      needs a `":gas"` suffix; why totals span both phases; why gas-phase write
      back is not implemented in `apply_to_phases`.
      **Runtime messages reworded: 7** (1 Bisection engine, 2 NR engine, 1
      solver, 3 tableau), all listed by the AST tool; the fragments tests match
      on survive (`V_liq_L`, `bridge`); the 22 tests in
      `test_nr_tableau_gas_liquid.py` pass.
      **Verification.** (1) AST comparison of all 14 files against `HEAD`, with
      docstrings removed: **no code changes**; only string literals differ, and
      only those 7. (2) Bit-level fingerprints for all three engines (NR,
      Bisection, PHREEQC — 814 + 265 + 182 values) byte-identical. (3) Full
      suite 2080 passed (unchanged). (4) Both surveys re-run. (5) The line-number
      citations in `STRONG_ION_INFERENCE_GENERALIZATION.md` shifted again with
      these edits (the `engines/nr/engine.py` and `engines/bisection/engine.py`
      ranges; the solver ones did not), so they were re-derived and re-verified
      with a stricter checker that requires each range to start on the defining
      line and end on its closing line — the checkpoint-12 checker was too
      lenient (it passed ranges that were one line off).
      **Decision applied:** `phreeqc_to_vlsim` (and the local `vlsim_name`) kept;
      docs reworded; rename logged in `OPEN_WORK.md`. **Logged, not touched
      (`OPEN_WORK.md`):** `activity_dispatch.py` as an orphan candidate (same
      situation as the deleted `api.py`/`factory.py`); the same false or stale
      statements found outside this directory
      (`chemistry/equilibria.py:510-511`, `core/gas_liquid_link.py:928,958,966`,
      `reactions/reaction_system.py:250`); two saved notebook outputs that
      quote the old warning text; the wider 108-line cleanup.

**Part D — gas-constant unification** (after Part C; see the plan doc's
"Gas-constant definitions (Part D)" audit and Decisions 11–15)

- [ ] 13. _Value-preserving._ Re-verify the plan's audit table with a fresh
      search. Replace `cv_loops._R_UNIV` (bit-identical copy) with an import
      from `units`; make `test_equilibrium_constants.py` import `R` instead of
      holding a literal. Add the guard test with an allowlist of every copy that
      remains (Decision 15). Bit-identical (fingerprint).
- [ ] 14. _Numerics-changing, own commit._ Derive `R_L_ATM_PER_MOL_K` from
      `R_J_PER_MOL_K`; repoint `core/phases.py` (rename all 25 files, Decision
      13), `partition.py`, `peng_robinson.py`, `cv_loops.py`, `plots.py`, the
      three ADM1/BSM2 files (Decision 14 — **confirm keep-or-unify first**),
      test literals and `Example1_mtp_well.ipynb` `R_LA`; empty the guard
      allowlist. Record a gas-liquid + BSM2 fingerprint before, the shift after
      (expect ~4.1e-7 relative on `nRT/V`, ~3.2e-7 on ADM1 `Ka(T)`), and
      re-baseline sentinels with a dated before/after comment.
- [ ] 15. Full suite green, then ship (see below).

### Checkpoint 1 inventory

Recorded 2026-09-20 from fresh searches over tracked `.py`, `.ipynb` and `.md`
files (`scratch/` is gitignored and excluded; the plan doc and this checklist
are excluded from the counts).

**Baseline:** full suite (`python -m pytest -x`, configured `testpaths`,
including `tests/validation/`) = **2066 passed, 0 failed, 0 skipped, 335
warnings, 1m43s**. This is the number to compare against at checkpoint 13.

**Absolute-path references** (`chemical_equilibrium.<module>` or `/<module>`),
by module, excluding historical docs (occurrences / files):

| Module | Occ. | Files | Fate |
|---|---|---|---|
| `nr_engine` | 73 | 32 | move (C8) |
| `engine` | 35 | 22 | move (C9) |
| `nr_tableau` | 16 | 4 | move (C8) |
| `activity_models` | 12 | 8 | delete (A4) |
| `phreeqc_engine` | 8 | 8 | move (C10) |
| `nr_solver` | 6 | 4 | move (C8) |
| `acid_base` | 4 | 2 | move (C9) |
| `sit` | 2 | 1 | delete (A4) |
| `strong_ions` | 1 | 1 | delete (B5) |
| **Total** | **157** | | |

By location (all matches, including the 11 in `docs/dev/ideas/`, which are
left alone; `docs/dev/implementation/shipped/` has none):

| Location | Occ. | Files |
|---|---|---|
| `PyOMES/chemical_equilibrium/` (docstrings and lazy imports, absolute form) | 16 | 6 |
| Other `PyOMES/` packages (`reactions/`, `core/`, `chemistry/types.py`) | 12 | 7 |
| `models/vlmodels/adm1/{base,bsm2}.py` | 2 | 2 |
| `tests/standalone/` | 85 | 16 |
| `tests/validation/**/test_*.py` | 5 | 3 |
| `tests/validation/speciation/_generate_notebooks.py` (inside string templates) | 8 | 1 |
| `tests/validation/speciation/*.ipynb` | 14 | 8 |
| `docs/tutorials/ArXiv_preprint/_generate_notebooks.py` | 2 (string templates) | 1 |
| `docs/tutorials/**/*.ipynb` | 5 | 4 |
| Current docs (`OPEN_WORK.md`, `STRONG_ION_INFERENCE_GENERALIZATION.md`, `NR_PRECIPITATION_CV_INTEGRATION.md`) | 8 | 3 |
| `docs/dev/ideas/` (historical, leave alone) | 11 | 4 |

Location rows sum to 168 (157 non-historical + 11 in `docs/dev/ideas/`).

**Relative imports inside `chemical_equilibrium/` that change:**

- Moving files import staying files: `engine.py`, `nr_engine.py`,
  `phreeqc_engine.py` → `.protocols`; `acid_base.py` → `.activity`. These gain
  a dot after the move.
- Moving files import each other: `engine.py` → `.acid_base`;
  `nr_engine.py` → `.nr_tableau`, `.nr_solver`; `nr_solver.py` → `.nr_tableau`.
  The engines do not import each other (only docstrings mention the other
  engines' paths), confirming the plan's audit.
- Staying files import moving files: `__init__.py` (`.engine`, `.acid_base`),
  `api.py` (`.engine`), `factory.py` (`.engine`).
- Parent-relative imports inside moving files gain a dot (`..` → `...`):
  `nr_tableau.py:44` (`..units`) and lazy `nr_tableau.py:435`
  (`..reactions.equilibrium`); `nr_solver.py:49` (`..core.phases`);
  lazy `engine.py:202-203` (`..chemistry.equilibria`, `..reactions.equilibrium`);
  lazy `nr_engine.py:238` (`..reactions.equilibrium`).
- Part A files: `activity.py:13`, `acid_base.py:42,921`, `engine.py:44`,
  `nr_engine.py:45` import `.activity_models`; `__init__.py:7` re-exports its
  aliases. Nothing imports `.sit`.

**External lazy / relative imports** (only fail when executed):
`reactions/reaction_system.py:268` (`..chemical_equilibrium.nr_engine`) and
`:280` (`..chemical_equilibrium.engine`); `reactions/equilibrium.py:18,281,283`,
`reactions/_shared.py:173`, `reactions/__init__.py:12`,
`core/gas_liquid_link.py:928,958` (`acid_base`), `core/property_calculator.py:17`,
`chemistry/types.py:33` (mostly docstring/comment mentions — check each at C11/C12).
`core/solvers.py:653` imports `..chemical_equilibrium.protocols`, which stays.

**Package-level imports** (`from PyOMES.chemical_equilibrium import ...`):
only `ChemicalEquilibriumEngine` and `BisectionChemicalEquilibriumEngine`
(2 test files). Both stay exported. Nothing imports `IdealActivityModel`,
`DaviesActivityModel` or `strong_ions_from_feed_molL` from the package root.

**Part A claims re-verified:**

- `davies_log10_gamma` / `davies_gamma`: only the two definitions and one
  self-call in `activity.py`; zero hits anywhere else in `.py`, `.ipynb`, `.md`.
- `SITActivityModel`: only `sit.py` and one test. Nothing in the package imports
  `sit.py`.
- `DaviesActivityModel` importers: `tests/validation/speciation/_generate_notebooks.py`
  (lines 768, 1136), `test_saturation_index.py:26`, notebooks `03`, `04`, `05`
  (import line + usage each), `test_liquid_phase_model.py:182`. Prose-only
  mentions (no import): notebook `06_phreeqc_benchmark.ipynb` (lines 45 and 130),
  generator line 1543, `acid_base.py` comments (lines 249, 255).
- `make_activity_model` users: `engine.py:44,306`, `nr_engine.py:45,366,965`,
  `acid_base.py:921`, `test_nr_tableau_gas_liquid.py:379`,
  `test_liquid_phase_model.py:186,265`, `OPEN_WORK.md`. `ActivityModel`
  protocol: `acid_base.py:42,202,870` (type hints), docstring in `nr_solver.py:374`.
- Water-property re-exports from `activity_models`: only `activity.py:13` in
  the package, plus one test (see below).

**Part B claims re-verified:** `strong_ions_from_feed_molL` /
`chemical_equilibrium.strong_ions` appear in `strong_ions.py`, `__init__.py`
(import + `__all__`), `test_strong_ions.py` (9 uses), and the comment at
`chemistry/registry.py:148`. **Zero hits in any `.ipynb`.**
`SALT_DISSOCIATION_MAP` is used only by `strong_ions.py:16,78`, defined at
`registry.py:117` and exported from `chemistry/__init__.py` (stays, Decision 8).

**Discrepancies from the plan doc (none block the phase; decisions needed at
the noted checkpoint):**

1. **`test_liquid_phase_model.py` has four shim-only tests, not two**, plus two
   factory tests: `test_davies_activity_model_alias` (l.181),
   `test_water_helpers_still_importable_from_speciation` (l.190),
   `test_sit_activity_model_alias` (l.260),
   `test_sit_epsilon_importable_from_speciation` (l.270); factory tests
   `test_make_activity_model_returns_davies_liquid_model` (l.185) and
   `test_make_activity_model_returns_sit` (l.264). The plan names only "two
   alias tests" and "the factory test". Proposal for C4: delete the two alias
   tests and the water-helpers test (pure re-export checks); repoint the SIT
   epsilon test to `PyOMES.thermo.sit_liquid_model`; repoint both factory tests.
2. **`sit.py` also re-exports `SIT_EPSILON`, `ION_CHARGES`, `_get_epsilon`**
   (the plan mentions only `SITLiquidModel` and the alias). Only the test above
   uses them; canonical definitions are in `thermo/sit_liquid_model.py`.
3. **A third van 't Hoff function exists:** `reactions/equilibrium.py:92`
   `vant_hoff_log_K(constraint, T_K)` (tested in `test_equilibrium_constraint.py`).
   Decision 3 names only `_vant_hoff_K` (`acid_base.py`, 9 call sites) and
   `_vant_hoff_log_K` (`nr_tableau.py`, 5 call sites, plus `nr_engine.py:651`).
   Check at C7 whether the third is the same maths; whether to fold it in is
   a decision for you.
4. **Current docs citing moved paths beyond the plan's C12 list:**
   `STRONG_ION_INFERENCE_GENERALIZATION.md` (~17 `nr_engine.py` / `nr_solver.py` /
   `engine.py` references with line numbers), `NR_PRECIPITATION_CV_INTEGRATION.md:144`,
   `upcoming/README.md:60,62`, `adm1/base.py:1048` comment,
   `protocols.py:72-73`, `thermo/liquid_phase_model.py:53`, `nr_tableau.py:117`,
   `nr_engine.py:85` comments, `01_bisection_engine_basics.ipynb:130` ("see the
   list in `engine.py`"), and `docs/architecture.md:391-394` (also lists the stale
   `speciation/` folder name). Plan cites `OPEN_WORK.md` and `architecture.md` only.
5. Plan says the lazy import is at `nr_tableau.py:435`; that line is
   `from ..reactions.equilibrium import`, correct as a lazy relative import
   that needs `...` after the move (line 410 is a docstring path).

## Shipping

- [ ] Full test suite green on the branch
- [ ] `git checkout main`
- [ ] `git merge --no-ff chemical-equilibrium-engines-subfolder -m "Merge chemical-equilibrium-engines-subfolder: move engines into engines/ subfolder, retire thermo shims, remove strong_ions"`
- [ ] `git tag chemical-equilibrium-engines-subfolder-shipped <commit-hash>`
- [ ] `git push && git push --tags` (both — tags are not pushed by default)
- [ ] `git branch -d chemical-equilibrium-engines-subfolder` and
      `git push origin --delete chemical-equilibrium-engines-subfolder`
- [ ] Move the plan doc + this checklist to `docs/dev/implementation/shipped/`, add a
      "Shipped" banner to both
- [ ] Update `docs/dev/implementation/upcoming/README.md`'s "Recently shipped" list
- [ ] Note the deliberate breaking changes (old engine module paths,
      `IdealActivityModel` / `DaviesActivityModel`, `strong_ions_from_feed_molL`)
      in release notes / changelog if one is kept
