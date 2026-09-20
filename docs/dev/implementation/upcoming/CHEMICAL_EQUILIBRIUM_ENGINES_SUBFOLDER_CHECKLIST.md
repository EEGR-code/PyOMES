# Phase Kickoff Checklist — chemical-equilibrium-engines-subfolder

> Checklist for [`CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER.md`](CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER.md),
> which is the source of truth for goals, audit, target layout and Decisions
> 1-10. All open questions there are resolved; do not re-litigate them. See
> `README.md`'s "Branching and tagging convention" for the branch/merge/tag
> rationale. Modelled on
> [`../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md`](../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md).

**Working rules for this phase**

- Pure refactor: no behaviour or numerical change. Anything that looks like a
  bug or improvement is noted, not fixed.
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

- [ ] 5. Search `.ipynb` files (earlier audit covered `.py` and `.md` only)
      for `strong_ions_from_feed_molL` / `PyOMES.chemical_equilibrium.strong_ions`.
      If clear: delete `strong_ions.py` and
      `tests/standalone/test_strong_ions.py`; remove the `__init__.py` import
      and `__all__` entry (Decision 10); drop the stale comment in
      `chemistry/registry.py`. `FeedState` and its fixtures stay.
- [ ] 6. Update the "Adjacent, out of scope" section of
      `STRONG_ION_INFERENCE_GENERALIZATION.md` (`strong_ions.py` removed;
      `SALT_DISSOCIATION_MAP` now has no consumer). Run `test_feed_state.py`
      and the speciation tests.

**Part C — `engines/` subfolder**

- [ ] 7. De-duplicate `acid_base._vant_hoff_K` / `nr_tableau._vant_hoff_log_K`
      into one helper in `PyOMES/thermo/` (Decision 3; filename decided here,
      e.g. `thermo/equilibrium_constants.py`). Own commit, own tests. Repoint
      `acid_base.py` and `nr_tableau.py`.
- [ ] 8. `git mv` the NR files into `engines/nr/` (`engine.py`, `tableau.py`,
      `solver.py`, `__init__.py`); rewrite relative imports (`..units` →
      `...units`, `..core.phases` → `...core.phases`, ...). Run `test_nr_*`
      and `test_equilibrium_classification.py`.
- [ ] 9. `git mv` the Bisection files into `engines/bisection/` (`engine.py`,
      `acid_base.py`, `__init__.py`). Run `test_speciation*.py` and
      `test_bisection_chemical_equilibrium_engine_alias.py`.
- [ ] 10. `git mv phreeqc_engine.py engines/phreeqc.py`. Run
      `tests/validation/speciation/`.
- [ ] 11. Update the package `__init__.py` re-exports, then external callers:
      `PyOMES/reactions/reaction_system.py`, `models/vlmodels/adm1/{base,bsm2}.py`,
      both notebook generators, notebooks, tests. Finish with a repo-wide
      search (`.py`, `.ipynb`, `.md`) for each old module path to confirm
      nothing still points at it.
- [ ] 12. Update docstring cross-references (~15 files using
      `PyOMES.chemical_equilibrium.engine...`-style paths) and current docs,
      including `OPEN_WORK.md` (cites `activity_models.make_activity_model`) and
      `docs/architecture.md` (lists `activity_models.py`, `sit.py`). Leave
      `shipped/` and `docs/dev/ideas/` alone.
- [ ] 13. Full suite green, then ship (see below).

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
