# Phase Kickoff Checklist — activity-model-parameter

> Checklist for [`ACTIVITY_MODEL_PARAMETER.md`](ACTIVITY_MODEL_PARAMETER.md), the
> source of truth for motivation, design, inventory and decisions; do not restate
> it here. Where this checklist and the note disagree, this checklist wins. See
> [`README.md`](README.md)'s "Branching and tagging convention". Modelled on
> [`../shipped/THERMO_SUBFOLDER_STRUCTURE_CHECKLIST.md`](../shipped/THERMO_SUBFOLDER_STRUCTURE_CHECKLIST.md).

**Working rules**

- Checkpoints 1-2 change signatures and call sites with no change in engine
  output, except the one intended change: an invalid activity-model name fails
  when the engine is built instead of at its first solve. Checkpoint 1 also adds
  tests. Checkpoint 3 changes notebook `source` lines and generator code only,
  plus re-running one notebook. Checkpoint 4 changes only prose and
  `OPEN_WORK.md`.
- The repo owner runs every `git` command (branch, add, commit, push, merge,
  tag). Read-only local git (`status`, `log`, `diff`, `show`, `ls-files`) is
  fine; nothing that touches the network. Git commands are handed over one per
  line (PowerShell 5.1: never chained with `&&`); commit messages use two `-m`
  flags (a strapline and one body paragraph) and no attribution lines. The repo
  is in OneDrive: a commit may ask "Rename ... index.lock ... failed. Should I
  try again? (y/n)"; answer `y`. No git command is run while a commit may be
  waiting.
- One checkpoint at a time: run the full suite first, edit, run the sanity
  checks, report what changed and the results as they are, record notes here,
  stop, and hand over the commit commands. Do not start the next checkpoint until
  told to. Edits are shown one per message, each preceded by a one-line summary.
- Full suite: `python -m pytest -p no:cacheprovider -q`.
- No compatibility aliases or shims: `use_activity=` and `thermo=` stop being
  accepted by the engines, and every call site is updated in place.
- Proof of no behaviour change: the engine-output fingerprint (below) is
  identical after checkpoints 1, 2 and 3; the BSM2 sentinel tests pass unchanged.
  Prose-only edits leave the AST identical with docstrings removed.
- Notebook edits touch only `source`, never saved outputs, except
  `docs/tutorials/reactions/chemistry_database.ipynb`, which is re-run because
  one saved output prints the removed name. Notebooks are re-run from the
  scratchpad directory so no output folders land in the repo. The tool sandbox
  caps each process at about 15 % of one CPU core: if a notebook does not finish,
  record how far it ran as a deviation. Notebooks and the suite are never run at
  the same time. Notebooks are never regenerated from the generator scripts.
- Most files involved are CRLF in the working tree. Bulk edits use a script that
  preserves each file's line endings and exact trailing bytes; check for bare LF
  after every edit.
- Leave `docs/dev/implementation/shipped/` and `docs/dev/ideas/` untouched.
- Docs and docstrings describe current behaviour only: no phase or checkpoint
  labels, no pointers to this checklist or the design note.
- Anything that looks like a bug or dead code beyond scope is logged in
  `OPEN_WORK.md`, not fixed.

## Pre-flight

- [ ] `git status -sb` clean apart from this phase's docs (the design note, this
      checklist, the `upcoming/README.md` entry and the new `OPEN_WORK.md` entry
      "A CV's `chemistry_db` activity model never reaches its speciation engine")
- [ ] `git log origin/main..main --oneline` empty
- [ ] Branch created off current `main`: `git switch -c activity-model-parameter`
- [ ] Those docs committed on the branch as the first commit
- [ ] Baseline full suite recorded here (pass / fail / skip counts)
- [ ] Fingerprint baseline captured. A scratchpad script (not committed) builds
      each engine over a fixed set of cases — pure water, carbonate, acetate,
      ammonium, phosphate, a mixed BSM2-like liquid, and one high-ionic-strength
      case — for ideal, Davies and SIT, through the **old** parameters, solves at
      two temperatures, and records the SHA-256 of every result value (pH, ionic
      strength, every species). Hash recorded here. The script takes the
      parameter style as an argument so the same cases run through the new
      parameter after each checkpoint.

## Checkpoints

- [ ] **1. Resolver, framework and engines.**
      `thermo/liquid/factory.py` (`make_activity_model(activity_model)`: a name
      in `"ideal"`/`"davies"`/`"sit"`, case-insensitive, or a model object passed
      through; `ValueError` listing the accepted names otherwise),
      `thermo/liquid/__init__.py` docstring, `thermo/framework.py` (remove
      `use_activity`), both engines (`activity_model=` in `__init__` and
      `from_reactions`, resolved at construction and stored as the object on
      `engine.activity_model`; `thermo=`, `use_activity` and `_liquid_activity`
      removed), `acid_base.py` and `nr/solver.py` (ideal checks by type; default
      model from the resolver), `reactions/reaction_system.py`
      (`_engine_config`, `configure_engine`), `monitoring/accuracy.py`
      (`check_ionic_strength(I_molL, activity_model)`), and their tests.
      New tests: each accepted name and a model object resolve to the right
      class; an unknown name raises at engine construction; a model object that
      is not `IdealLiquidModel` and has no `name` takes the ionic-strength loop.
      Sanity: fingerprint identical; suite green; `use_activity` absent from
      `PyOMES/` and the edited tests.
- [ ] **2. Stirred tank and models.** `templates/stirred_tank/builder.py`
      (`chemistry(activity_model="ideal")`), `configs.py`
      (`ChemistryConfig.activity_model`), `factory.py` (pass-through),
      `models/vlmodels/adm1/base.py` and `bsm2.py`, and their tests
      (`test_builder.py`, `test_configs.py`). Sanity: fingerprint identical; BSM2
      sentinels pass; suite green; `use_activity` absent from `PyOMES/`,
      `models/` and `tests/standalone/`.
- [ ] **3. Notebooks and generators.** The 9 notebooks listed in the design note
      (edit `source` only) and both `_generate_notebooks.py` scripts. Re-run
      `chemistry_database.ipynb` after changing its cell 6 to print and assert on
      `activity_model`. Sanity: every code cell of every edited notebook parses;
      each notebook's diff touches only the intended `source` lines (plus the
      re-run notebook's outputs); fingerprint identical; `use_activity` absent
      from `docs/tutorials/` and `tests/`.
- [ ] **4. Docs and close-out.** `README.md` and `docs/architecture.md`
      (parameter description, checked against the code); delete the `OPEN_WORK.md`
      entry "`chemical_equilibrium`'s `use_activity`/`activity_model` split could
      be one parameter". Sanity: repo-wide sweep finds no `use_activity` outside
      `docs/dev/implementation/shipped/` and `docs/dev/ideas/`; relative links in
      edited files resolve; `git diff --stat` shows only `.md` files.

## Notes

<!-- Baseline counts, fingerprint hash, deviations and findings go here as each
checkpoint lands. -->

## Shipping

- [ ] Full test suite green on the branch
- [ ] `git switch main`
- [ ] `git merge --no-ff activity-model-parameter -m "Merge activity-model-parameter: <summary>"`
- [ ] `git tag activity-model-parameter-shipped`
- [ ] `git push origin main` and `git push origin activity-model-parameter-shipped`
- [ ] `git branch -d activity-model-parameter`
- [ ] Move the design note and this checklist to `docs/dev/implementation/shipped/`,
      add a "Shipped" banner to both
- [ ] Update `docs/dev/implementation/upcoming/README.md`: remove the "Design
      discussions" entry and add one to "Recently shipped"
