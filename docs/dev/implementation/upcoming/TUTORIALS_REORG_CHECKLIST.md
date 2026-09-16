# Phase Kickoff Checklist — tutorials-reorg

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
- [ ] 10. Convert `docs/tutorials/reactions/*.py` to notebook form, for consistency with the
      rest of `docs/tutorials/` (all notebooks except this folder). Gated on: (a) checkpoint 7's
      de-coupling landing, so `reaction_system.py` has no importers left; (b) the two
      pre-existing bugs below being fixed first — a notebook's outputs should be genuine baked
      results, not a crash, so `chemistry_database.py`/`partition_model.py` can't be converted
      until they actually run. `fba/fba_toy.py` and `fba/fba_ecoli_core.py` have no such
      blocker and could convert independently. Revisit whether `templates/`'s three plain-script
      files (`cstr_fermenter.py`, `fed_batch_fermenter.py`, `microplate_fermenter.py`, from
      checkpoint 6) should get the same treatment once this lands.

## Follow-ups found during this phase (not fixed here, tracked so they aren't lost)

- [ ] `tests/validation/speciation/{07_iron_oxidation,08_iron_oxidation_and_precipitation}.ipynb`
      have zero pytest coverage. Their own P1 prediction (>70% Fe2+ conversion) is already
      marked `[CHECK]` (failed) in the notebook's own output — actual conversion is 0.2%,
      likely because the demo equilibrates to pH 6.5 while the Singer-Stumm rate
      constant/narrative assumes pH ~4.7. Needs investigation before writing tests against it.
- [ ] `docs/tutorials/reactions/chemistry_database.py` crashes: `dataclasses.replace(...,
      activity_model=...)` — `ThermoFramework.__init__()` no longer accepts `activity_model`.
      Pre-existing, confirmed broken in `demos/model_api/chemistry/chemistry_database.py` too
      (not caused by this phase's move) — likely drift from the StirredTankBuilder/ThermoFramework
      refactor.
- [ ] `demos/usecases/03_cstr_dilution_rate_sweep.ipynb` appears to be an orphaned duplicate —
      not produced by `demos/usecases/_generate_notebooks.py`, and not referenced anywhere in
      the repo (confirmed by grep). Likely a leftover from before this notebook was originally
      curated into `docs/tutorials/` (pre-dating this phase). Probably safe to delete, but
      left untouched since it's outside this phase's scope.
- [ ] `docs/tutorials/reactions/partition_model.py` crashes: `h2s.beta(...)` —
      `HenryEquilibrium` (replacement for the deprecated `HenryPartition`) has no `.beta()`
      method. Pre-existing, confirmed broken in `demos/model_api/chemistry/partition_model.py`
      too.
- [ ] `docs/tutorials/D2C_workshop/raw_construction.py` runs (its import bug is now fixed —
      see checkpoint 7) but produces pH 12.089 against a `PHController(setpoint=5.0)`, plus
      `ConservationWarning`s for O/C/H and charge balance exceeding their stated thresholds.
      Unlike the two bugs above, this was never previously observable at all — the script has
      never run successfully before (confirmed: the original `demos/` copy raises
      `ModuleNotFoundError` before reaching the simulation). Likely a PHController
      tuning issue or a stoichiometry mismatch; not investigated. Does not block this phase —
      confirmed none of `Example1_mtp_well.ipynb`/`Example2_batch_fermenter.ipynb`/
      `Example3_CSTR.ipynb` (the actual tutorial content in this folder) import or depend on
      `raw_construction.py` — it's an independent, standalone comparison script.

## Shipping

- [ ] Full test suite green on the branch (`pytest`)
- [ ] `git checkout main`
- [ ] `git merge --no-ff tutorials-reorg -m "Merge tutorials-reorg: ..."`
- [ ] `git tag tutorials-reorg-shipped <commit-hash>`
- [ ] `git push && git push --tags`
- [ ] `git branch -d tutorials-reorg` and `git push origin --delete tutorials-reorg`
- [ ] Move this checklist to `docs/dev/implementation/shipped/`, add a "Shipped" banner
