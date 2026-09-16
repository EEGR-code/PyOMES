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
- [ ] 3. Move `demos/model_api/chemistry/speciation/` → `tests/validation/speciation/`
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
- [ ] 4. Copy `demos/model_api/chemistry/{reaction_system.py,fba/,chemistry_database.py,
      partition_model.py}` → `docs/tutorials/reactions/`; add a short README. Sanity check:
      `python docs/tutorials/reactions/reaction_system.py` runs clean.
- [ ] 5. Copy `demos/features/{ChemicalEquilibriumProtocol,SolverProtocols}/` →
      `docs/tutorials/protocols/` (notebooks + READMEs, not generators); recompute relative
      links for the new depth and the `speciation/` move; write a top-level
      `docs/tutorials/protocols/README.md`. Sanity check: every notebook's relative links
      resolve (no 404s when opened).
- [ ] 6. Copy `demos/builder/` → `docs/tutorials/templates/`; fix `README.md`'s link depths and
      its `D2Cworkshop` cross-link. Sanity check: all 4 `.py` files run clean from the new path.
- [ ] 7. Delete `demos/model_api/D2Cworkshop/basic_layout/` entirely; copy
      `demos/model_api/D2Cworkshop/updated_layout/` → `docs/tutorials/D2C_workshop/`; replace
      the stale `README.md` stub with `0_README.ipynb`-derived content; fix
      `raw_construction.py`'s two docstring links **and** the pre-existing `_chem_dir` path bug
      (`parents[1]` → should resolve to the new `reactions/` sibling). Sanity check:
      `raw_construction.py` actually imports successfully (it doesn't today).
- [ ] 8. Verify: run/open every copied file from checkpoints 4–7 in its new location, confirm
      no broken links/imports. **Gate — do not proceed to 9 until this passes.**
- [ ] 9. Delete migrated originals from `demos/` (`builder/`, `features/`, `model_api/`
      entirely); rewrite `demos/README.md` for the much smaller remaining tree
      (`usecases/` + `aerobic_fermentation_stoichiometry.ipynb`); remove
      `demos/model_api/README.md`.

## Shipping

- [ ] Full test suite green on the branch (`pytest`)
- [ ] `git checkout main`
- [ ] `git merge --no-ff tutorials-reorg -m "Merge tutorials-reorg: ..."`
- [ ] `git tag tutorials-reorg-shipped <commit-hash>`
- [ ] `git push && git push --tags`
- [ ] `git branch -d tutorials-reorg` and `git push origin --delete tutorials-reorg`
- [ ] Move this checklist to `docs/dev/implementation/shipped/`, add a "Shipped" banner
