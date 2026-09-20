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

- [ ] 1. Inventory every import site (`.py` and `.ipynb`, plus the two
      notebook generators' string templates and lazy in-function imports) of
      the modules being moved or deleted. Record counts in "Checkpoint 1
      inventory" below. Baseline: full test suite green.

**Part A — thermo consolidation** (retire `activity_models.py` / `sit.py`)

- [ ] 2. Add `ActivityModel` and `make_activity_model` to `thermo/` (in
      `liquid_phase_model.py` or a small `thermo/factory.py` — decide here);
      export from `PyOMES.thermo`. Fix the stale `PyOMES/speciation/`
      docstring paths in `liquid_phase_model.py`. `activity_models.py` still
      exists after this checkpoint.
- [ ] 3. Repoint importers at `PyOMES.thermo`: `engine.py`, `nr_engine.py`,
      `acid_base.py`, `activity.py` (`thermo.water_properties`),
      `test_nr_tableau_gas_liquid.py`, and the `DaviesActivityModel` users in
      `tests/validation/speciation/` (tests, generator, 3 notebooks →
      `DaviesLiquidModel`). Remove `IdealActivityModel` /
      `DaviesActivityModel` from `__init__.py` and `__all__` (Decision 9).
- [ ] 4. Delete `sit.py` and `activity_models.py`; delete the two alias tests
      in `test_liquid_phase_model.py`; repoint the factory test there. Delete
      `davies_log10_gamma` / `davies_gamma` and fix `activity.py`'s docstring
      (Decision 6). Run `test_liquid_phase_model.py`,
      `test_thermo_framework.py`, `test_speciation.py`,
      `tests/validation/speciation/`.

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

_Filled in when checkpoint 1 lands._

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
