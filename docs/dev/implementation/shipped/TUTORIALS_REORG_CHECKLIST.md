# Phase Kickoff Checklist — tutorials-reorg

> **Status: Shipped 2026-09-17** — merged into `main` via `git merge --no-ff`
> as commit `1970c8e`, tagged `tutorials-reorg-shipped`. All 10 checkpoints
> plus the pre-ship `pyproject.toml` testpaths fix landed. Full suite green
> post-merge: 2029 passed, 36 skipped. Five deferred bugs/cleanup items split
> into their own phase — see
> [`../upcoming/TUTORIALS_FOLLOWUPS_CHECKLIST.md`](../upcoming/TUTORIALS_FOLLOWUPS_CHECKLIST.md).

> Reorganize demos/ examples into topic-based docs/tutorials/ subdirectories, and split
> validation/performance content into tests/validation/ and tests/performance/. See the
> design discussion this was worked out from (conversation, not a doc) for full rationale
> on naming and scope decisions — summarized in this checklist's Checkpoints.

## Pre-flight

- [ ] `git status -sb` clean (no stray uncommitted work left over from a previous task)
- [ ] `git log origin/main..main --oneline` empty (nothing unpushed sitting around from earlier work)
- [ ] Branch created off current `main`: `git checkout -b tutorials-reorg`
- [ ] This checklist file committed on that branch as the first commit

## During

- [ ] Checkpoints tracked below as they land, one commit per checkpoint
- [ ] **If work stalls or is paused before shipping:** add a status banner to the top of this
      file immediately — what's built, what's tested, why it stopped, which commit it's on.

### Checkpoints

- [x] 1. Delete `docs/publications/ArXiv_preprint/` (figures + README); fix its one
      cross-reference in `docs/tutorials/README.md`. Sanity check: `grep -r "publications/ArXiv_preprint"`
      returns nothing.
- [x] 2. Remove `04_compare_runtime_by_usecase.ipynb` — delete the file; edit
      `demos/usecases/_generate_notebooks.py` (drop the 04 section, its README row, fix the
      module docstring) and regenerate `0_README.ipynb`; hand-edit the `"runtime-md"` cell in
      the already-committed `docs/tutorials/01_predict_ph_simple_liquid.ipynb` directly (do not
      regenerate — `02_kinetic_co2_equilibration_microplate_well.ipynb` has baked outputs a
      full regenerate would wipe) and fix the same cell's source in
      `docs/tutorials/_generate_notebooks.py`; trim stale mentions in `docs/tutorials/README.md`,
      `OPEN_WORK.md`, `NOTEBOOK_GENERATOR_REMOVAL.md`; create `tests/performance/` (empty, short
      README). Sanity check: `grep -r "04_compare_runtime_by_usecase"` returns nothing.
- [x] 3. Move `demos/model_api/chemistry/speciation/` → `tests/validation/speciation/`
      (including its self-contained generator, updating its two cosmetic path strings). Fix
      relative links in the ~10 files that reference it (`demos/usecases/*`,
      `docs/tutorials/{01_predict_ph_simple_liquid.ipynb,_generate_notebooks.py}`,
      `demos/features/ChemicalEquilibriumProtocol/_generate_notebooks.py` + its 3 rendered
      notebooks, `NOTEBOOK_GENERATOR_REMOVAL.md`, `NR_PRECIPITATION_SPECIATION.md`). Write
      `test_*.py` files under `tests/validation/speciation/` asserting each notebook's core
      numerical claim — new coverage for 07/08 (iron oxidation, currently zero), tightened
      coverage for 01/03/06 (currently loose/partial); 02/05/10 already well covered, no new
      test needed. Notebooks are kept as the plotted/narrative companion (same pairing as
      `demos/builder/batch_fermenter.py`+`.ipynb`). Sanity check: `pytest tests/validation/`
      and `pytest tests/standalone/` both green; `grep -r "model_api/chemistry/speciation"`
      returns nothing.
- [x] 4. Copy `demos/model_api/chemistry/{reaction_system.py,fba/,chemistry_database.py,
      partition_model.py}` → `docs/tutorials/reactions/`; add a short README. Sanity check:
      `python docs/tutorials/reactions/reaction_system.py` runs clean.
- [x] 5. Copy `demos/features/{ChemicalEquilibriumProtocol,SolverProtocols}/` →
      `docs/tutorials/protocols/` (notebooks + READMEs, not generators); recompute relative
      links for the new depth and the `speciation/` move; write a top-level
      `docs/tutorials/protocols/README.md`. Sanity check: every notebook's relative links
      resolve (no 404s when opened).
- [x] 6. Copy `demos/builder/` → `docs/tutorials/templates/`; fix `README.md`'s link depths and
      its `D2Cworkshop` cross-link. Sanity check: all 4 `.py` files run clean from the new path.
- [x] 7. Delete `demos/model_api/D2Cworkshop/basic_layout/` entirely; copy
      `demos/model_api/D2Cworkshop/updated_layout/` → `docs/tutorials/D2C_workshop/`; replace
      the stale `README.md` stub with `0_README.ipynb`-derived content; fix
      `raw_construction.py`'s pre-existing `_chem_dir` path bug (`parents[1]` currently resolves
      to a directory that doesn't exist). **Do not point the fix at `docs/tutorials/reactions/`**
      — inline `raw_construction.py`'s own copy of the three reaction-factory functions
      (`make_aerobic_growth_on_acetate`, `make_acetate_dissociation`, `make_co2_partition`)
      instead of `import reaction_system as chem`, so `D2C_workshop/` is self-sufficient and
      doesn't depend on a sibling tutorial folder (tutorials are standalone by default). This
      also removes the only reason `reactions/reaction_system.py` needed to stay a `.py` module
      — see checkpoint 10. Sanity check: `raw_construction.py` runs standalone with no imports
      outside its own folder (besides `PyOMES` itself).
- [x] 8. Verify: run/open every copied file from checkpoints 4–7 in its new location, confirm
      no broken links/imports. **Gate — do not proceed to 9 until this passes.** PASSED — see
      commit for details (1 real bug found and fixed: a checkpoint-6 typo, `../../templates/`
      should have been `../templates/`; the two already-tracked follow-up bugs are the only
      remaining script failures, both expected).
- [x] 9. Delete migrated originals from `demos/` (`builder/`, `features/`, `model_api/`
      entirely); rewrite `demos/README.md` for the much smaller remaining tree
      (`usecases/` + `aerobic_fermentation_stoichiometry.ipynb`); remove
      `demos/model_api/README.md`. **Gap found mid-checkpoint:** `demos/model_api/
      export_results.py` and `solver_comparison.py` were never in `demos/model_api/README.md`'s
      documented file table, so no earlier investigation surfaced them — caught only because
      `git rm -r demos/model_api` swept them up unexpectedly. Neither fit an existing
      `docs/tutorials/` folder's topic. Resolved: `solver_comparison.py` → converted to a
      notebook, `docs/tutorials/protocols/SolverProtocols/02_comparing_system_solvers.ipynb`
      (fills a real gap — that folder had no engine/solver *comparison* table, unlike its
      `ChemicalEquilibriumProtocol` sibling). `export_results.py` → converted to a notebook in a
      new `docs/tutorials/results/` folder (post-processing a `BatchResult` is orthogonal to
      every existing folder's topic). Both conversions executed for real (not fabricated
      output) via a shared exec-and-capture-stdout harness, matching the "every notebook carries
      genuine outputs" convention. After the repo-wide re-run of checkpoint 8's link-check
      (worth doing again after any deletion, not just after copies), two more stale forward-
      references surfaced and were fixed: `demos/usecases/0_README.ipynb`'s `../features/` link
      (target now gone) and `docs/tutorials/ArXiv_preprint/01_predict_ph_simple_liquid.ipynb`'s
      link to the now-deleted `demos/model_api/chemistry/reaction_system.py`.
- [x] 10. Convert `docs/tutorials/reactions/*.py` to notebook form. Converted the 3 unblocked
      files — `reaction_system.py` → `reaction_system.ipynb`, `fba/fba_toy.py` →
      `fba/fba_toy.ipynb`, `fba/fba_ecoli_core.py` → `fba/fba_ecoli_core.ipynb` — each executed
      for real via a shared exec-and-capture-stdout harness (not fabricated output), verified
      against the original `.py` scripts' actual printed output before deleting them.
      `chemistry_database.py`/`partition_model.py` stay `.py`, still blocked on their
      pre-existing bugs. Note: `fba_ecoli_core.py` used `Path(__file__)` to locate
      `ecoli_core.json`, which doesn't exist in a notebook — changed to a bare relative path
      (Jupyter's default cwd is the notebook's own directory). Updated all cross-references
      (`D2C_workshop/{raw_construction.py,README.md}`, `ArXiv_preprint/`, `reactions/README.md`)
      and re-ran the full link-checker + every remaining `.py` script — clean, only the two
      already-tracked bugs fail. Revisit whether `templates/`'s three plain-script files
      (`cstr_fermenter.py`, `fed_batch_fermenter.py`, `microplate_fermenter.py`, from checkpoint
      6) should get the same treatment as a future follow-up.

- [x] Pre-ship fix: `pyproject.toml` `testpaths`. Found while doing a final pre-ship
      verification (running `pytest tests/` broadly rather than the specific paths used at each
      checkpoint): CI runs bare `pytest` (`.github/workflows/tests.yml`), and
      `[tool.pytest.ini_options] testpaths` was `["tests/standalone"]` only — meaning the new
      `tests/validation/speciation/test_*.py` coverage from checkpoint 3b would **never have run
      in CI**, silently defeating the whole point of writing real pytest assertions instead of
      leaving that content as notebooks-only. Fixed: `testpaths = ["tests/standalone",
      "tests/validation"]`. Verified bare `pytest` now collects both — 2029 passed, 36 skipped,
      matching the explicit `pytest tests/validation/ tests/standalone/` runs from earlier
      checkpoints exactly. Also fixed `docs/tutorials/results/README.md` and its notebook to
      install via the project's existing `pip install -e ".[export]"` extras group instead of an
      ad-hoc `pandas pyarrow` install.
## Follow-up phase

Five bugs/cleanup items were found during this phase but deliberately not fixed here — they're
debugging/investigation work, not file moves, and didn't block this phase's own goal. Split into
their own phase, [`../upcoming/TUTORIALS_FOLLOWUPS_CHECKLIST.md`](../upcoming/TUTORIALS_FOLLOWUPS_CHECKLIST.md),
so this phase could ship as one clean, complete unit: the two broken `reactions/` scripts
(`chemistry_database.py`, `partition_model.py`), `D2C_workshop/raw_construction.py`'s
control-loop issue, iron-oxidation notebooks' missing test coverage, and the orphaned
`demos/usecases/03_cstr_dilution_rate_sweep.ipynb` duplicate.

## Shipping

- [x] Full test suite green on the branch (`pytest`)
- [x] `git checkout main`
- [x] `git merge --no-ff tutorials-reorg -m "Merge tutorials-reorg: ..."`
- [x] `git tag tutorials-reorg-shipped <commit-hash>`
- [x] `git push && git push --tags`
- [x] `git branch -d tutorials-reorg` and `git push origin --delete tutorials-reorg`
- [x] Move this checklist to `docs/dev/implementation/shipped/`, add a "Shipped" banner
