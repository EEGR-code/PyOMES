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

- [x] `git status -sb` clean apart from this phase's docs (the design note, this
      checklist, the `upcoming/README.md` entry and the new `OPEN_WORK.md` entry
      "A CV's `chemistry_db` activity model never reaches its speciation engine")
- [x] `git log origin/main..main --oneline` empty
- [x] Branch created off current `main`: `git switch -c activity-model-parameter`
- [x] Those docs committed on the branch as the first commit (`bd8b372`)
- [x] Baseline full suite recorded here (pass / fail / skip counts)
- [x] Fingerprint baseline captured. A scratchpad script (not committed) builds
      each engine over a fixed set of cases — pure water, carbonate, acetate,
      ammonium, phosphate, a mixed BSM2-like liquid, and one high-ionic-strength
      case — for ideal, Davies and SIT, through the **old** parameters, solves at
      two temperatures, and records the SHA-256 of every result value (pH, ionic
      strength, every species). Hash recorded here. The script takes the
      parameter style as an argument so the same cases run through the new
      parameter after each checkpoint.

## Checkpoints

- [x] **1. Resolver, framework and engines** (with checkpoint 2's package and
      model call sites folded in; see Notes).
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
- [x] **2. Stirred tank and models.** Done inside checkpoint 1 (see Notes). `templates/stirred_tank/builder.py`
      (`chemistry(activity_model="ideal")`), `configs.py`
      (`ChemistryConfig.activity_model`), `factory.py` (pass-through),
      `models/vlmodels/adm1/base.py` and `bsm2.py`, and their tests
      (`test_builder.py`, `test_configs.py`). Sanity: fingerprint identical; BSM2
      sentinels pass; suite green; `use_activity` absent from `PyOMES/`,
      `models/` and `tests/standalone/`.
- [x] **3. Notebooks and generators.** The 9 notebooks listed in the design note
      (edit `source` only) and both `_generate_notebooks.py` scripts. Re-run
      `chemistry_database.ipynb` after changing its cell 6 to print and assert on
      `activity_model`. Sanity: every code cell of every edited notebook parses;
      each notebook's diff touches only the intended `source` lines (plus the
      re-run notebook's outputs); fingerprint identical; `use_activity` absent
      from `docs/tutorials/` and `tests/`.
- [x] **4. Docs and close-out.** `README.md` and `docs/architecture.md`
      (parameter description, checked against the code); the three usage
      examples that pass a `speciation_level` argument `chemistry()` and
      `ChemistryConfig` no longer have (`templates/stirred_tank/builder.py:17`,
      `templates/stirred_tank/factory.py:20` docstrings and
      `docs/tutorials/templates/README.md:61`; found in checkpoint 1); delete the `OPEN_WORK.md`
      entry "`chemical_equilibrium`'s `use_activity`/`activity_model` split could
      be one parameter". Sanity: repo-wide sweep finds no `use_activity` outside
      `docs/dev/implementation/shipped/` and `docs/dev/ideas/`; relative links in
      edited files resolve; `git diff --stat` shows only `.md` files.

## Notes

**Pre-flight (2026-09-25, on `bd8b372`).** Baseline suite: 2104 passed, 0
failed, 166 warnings (129 s). Fingerprint baseline (old parameters):
`c649ad00ad8fc1501f458122aeedd932d9def639ad168b65bbe78a38876ac763` over 1,452
values, identical on repeated runs. It covers four routes, each for ideal /
Davies / SIT: the Bisection engine built directly (7 cases × 25 and 37 °C); the
NR engine from `from_reactions` (3 cases × 2 temperatures, plus a calcite
precipitation case); `ReactionSystem.configure_engine` with both solvers,
solving from a phase; and `StirredTankBuilder.chemistry()` with two substrates
(the only way the factory calls `configure_engine`), recording which model class
the reaction system resolves plus one advance step. The first three routes give
different values for each model, so the hash would catch a model mix-up. A
default single-substrate tank declares no equilibria and has no engine, so its
activity setting changes no number; hence the class check for route 4.

**Checkpoint 1 (2026-09-25).** Suite: 2130 passed, 0 failed, 166 warnings
(2104, minus 4 removed tests — `test_default_use_activity_false` and the three
`thermo=` engine tests — plus 29 in the new
`tests/standalone/test_activity_model_argument.py` and 1 new monitor test).
Fingerprint through the new parameter: `c649ad00…a38876ac763`, identical to the
baseline. No `use_activity` left in `PyOMES/`, `models/` or the tests outside the
generator and notebooks (checkpoint 3). Every edited file keeps its stored line
endings.

- **Deviation: checkpoint 2 folded into 1.** The stirred-tank factory and the
  ADM1/BSM2 builders pass `use_activity=` to `configure_engine` and the engines,
  so they break as soon as the engines change; a separate checkpoint 1 could not
  have a green suite.
- **Settled while implementing** (owner's decisions, 2026-09-25):
  - `make_activity_model` raises `TypeError` for an object with no callable
    `gamma`, so a non-model fails when the engine is built. `name` is not
    required: it is only a label.
  - `solve_from_equilibrium_set`'s `activity_model` is required; its silent
    `None` → ideal branch is gone (its one caller always passed a model).
  - The engines' default stays `"ideal"`, with no warning: it is a documented
    default, not a fallback. A condition-based warning (ideal at high ionic
    strength) is the existing `OPEN_WORK.md` entry "No engine emits a
    high-ionic-strength warning".
  - `ReactionSystem.configure_engine` resolves the model immediately (a bad
    name fails at that call) and stores the resolved object.
  - `ChemistryConfig.activity_model` takes a name or a model object (option A):
    it is checked in `__post_init__` and stored as given; `to_dict()` keeps a
    model object as the object instead of letting `asdict` flatten it into a
    dict. A named model round-trips through JSON; a model object round-trips in
    Python only, as its docstring says.
- The resolver's accepted names are now exactly the models' `name` labels; the
  old unused aliases (`"daviesliquidmodel"` and so on) are gone.
- `configure_engine`'s docstring and error message referred to a
  `ReactionSystem.advance()` that does not exist; they now say the engine is
  built on first `engine` access (a CV's first step).
- The three models keep no state between calls (Ideal and Davies are frozen
  dataclasses; nothing in `sit.py` assigns to `self`), so building the model
  once per engine instead of once per solve is equivalent, as the fingerprint
  confirms.

**Checkpoint 3 (2026-09-25).** 9 notebooks and both generators edited (+26/−30
lines). Notebooks were edited at the JSON level after checking that all 9
re-serialise byte-for-byte, so only the targeted `source` lines changed; a script
verified no saved output changed. Code cells: the pair became `activity_model=`
(`02_multi_component_systems` got `activity_model="ideal"` for its
`use_activity=False` engine). Prose: `01_bisection_engine_basics` and
`02_nr_engine_basics` cell 8 describe the single argument;
`chemistry_database` cell 5 no longer calls the properties "backward-compatible"
(that was never accurate), and cell 6 prints and asserts on
`thermo.activity_model`. Every code cell of every edited notebook parses.

- **How the notebooks were run.** `nbclient` is not installed, so notebooks were
  run by a scratchpad script that executes code cells in order in one namespace
  (from the scratchpad directory, `MPLBACKEND=Agg`) and captures each cell's
  stdout. `chemistry_database.ipynb` was run in full: cells 2, 4, 8 and 10
  reproduced their saved output exactly, and only cell 6's saved text was
  replaced (now `ideal` / `davies`). The other 8 notebooks were run from the
  start up to their last edited cell, read-only: every cell with a saved output
  (24 cells) printed exactly that output; `01_predict_ph_simple_liquid` has no
  saved outputs and ran without error to cell 15. No notebook was regenerated.
- Fingerprint identical (`c649ad00…a38876ac763`); suite 2130 passed; no
  `use_activity` left in `PyOMES/`, `models/`, `tests/` or `docs/tutorials/`.

**Checkpoint 4 (2026-09-25).** `README.md` and `docs/architecture.md` describe the
single `activity_model` argument. `OPEN_WORK.md`: the entry
"`chemical_equilibrium`'s `use_activity`/`activity_model` split could be one
parameter" is deleted, and the "package root exports only the Bisection engine"
entry no longer points at it. Repo-wide sweep: no `use_activity` or
`speciation_level` outside `shipped/`, `ideas/` and this phase's own docs (the
design note, this checklist and its `upcoming/README.md` entry, which move to
`shipped/` or are replaced at shipping). Relative links resolve; line endings
unchanged.

- **The three `speciation_level` examples** were wrong in more ways than the
  argument: this tank declares no equilibria, so it has no speciation engine,
  its pH is `nan`, and `chemistry(...)` has no effect on it. The `.chemistry(...)`
  / `chemistry=` line is removed from all three rather than rewritten.
  `builder.py`'s example also ended with `result.pH[-1]`, which raises `KeyError`
  (`pH` is keyed by CV); it now ends with `result.liquid_mol["main"]["Yeast"][-1]`.
  The tutorials README example used `PressureReliefVent` without importing it;
  the import is added. All three examples run (with `n_steps` reduced for speed).
- **Deviation:** the checklist expected checkpoint 4's diff to be `.md` only;
  `builder.py` and `factory.py` also change, docstrings only (AST identical with
  docstrings removed).
- **Not fixed, noted:** the builder example adds no inoculum, so its biomass stays
  0; the tutorials README example attaches a `PHController` to a tank that
  computes no pH. Both predate this phase and are example-quality issues, not
  wrong API.

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
