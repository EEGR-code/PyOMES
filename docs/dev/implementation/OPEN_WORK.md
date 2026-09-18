# Open work

Standalone follow-up items surfaced during other phases — not yet
scoped as their own phase, no branch, no checklist. Referenced from
[`upcoming/README.md`](upcoming/README.md)'s "Priority order" section.

## README.md is broadly stale, well beyond any single phase

Surfaced 2026-09-15 while fixing `stirred-tank-template`'s checkpoint 9
doc cross-references. Root `README.md` describes an architecture from
well before several shipped phases:

- **`speciation/` package** (Overview/Repository Layout) — renamed to
  `chemical_equilibrium/` at the close of `LAYER1_GAP_CLOSURE`
  (shipped 2026-07-03).
- **`PropertySolver`/`ReactionModel`/`MultiCVSystem`** (Core
  Abstractions, Key Design Patterns) — retired by `state-unification`
  (2026-05-22) and `simulation-class` (2026-05-27); `Simulation` is now
  the single orchestration pathway.
- **`create_standalone_fermenter`/`CUFermentationSpeciation`** (Quick
  Start, Repository Layout) — a different, already-deleted BioSTEAM-era
  API; retired entirely by `cufermenter-sunset` (shipped 2026-06-01,
  9,093 lines deleted). Unrelated to the `StirredTankBuilder`/
  `StirredTankFactory` template this migration moved.
- **Test Coverage table** references test files that may no longer
  exist (`test_fermenter_construction.py`, `test_factory.py`,
  `test_gas_liquid_volume.py`) — not verified, flagging alongside the
  rest rather than auditing file-by-file here.

Needs a proper documentation-refresh pass against current `main` —
Repository Layout, Core Abstractions, Quick Start, and the Test
Coverage table all need re-deriving from the actual codebase, not a
line-by-line patch. Deliberately **not** done as part of
`stirred-tank-template` — that phase's checkpoint 9 only touched the
one Repository Layout line (`models/fermenter/` — now removed, since
that directory no longer exists) that was directly about its own move.

## `tests/run_tests.py` imports a deleted `create_standalone_fermenter`

Surfaced 2026-09-18 during a README.md audit (see the "README.md is
broadly stale" entry above). `tests/run_tests.py` still does
`from PyOMES import PressureReliefController, PHController,
create_standalone_fermenter` and calls it in `_fermenter_minimal()` /
`_fermenter_full()`. That name was never restored after
`cufermenter-sunset` (2026-06-01) — it isn't exported from
`PyOMES/__init__.py` and doesn't exist under `models/vlmodels/`
either, so this script currently fails on import. It isn't part of
the pytest suite (`pyproject.toml`'s `testpaths` only covers
`tests/standalone` and `tests/validation`), so it doesn't show up as
a CI failure — likely why it's gone unnoticed. Needs its own pass:
either migrate it to `StirredTankBuilder` (mirroring the
`stirred-tank-template` migration already done for the demo
notebooks) or delete it if it's fully superseded by
`tests/standalone`.

## `chemical_equilibrium`'s `use_activity`/`activity_model` split could be one parameter

Surfaced 2026-09-18 during the README.md audit, while checking
`StirredTankBuilder.chemistry()`'s `use_activity: bool` +
`activity_model: str = "davies"` signature for the README's Quick
Start rewrite. The split is threaded from
`PyOMES.chemical_equilibrium.activity_models.make_activity_model(use_activity,
activity_model)` through `engine.py`, `factory.py`, and
`nr_engine.py` — `activity_model` is only meaningful when
`use_activity=True`, and `"ideal"` is not itself a valid value for
`activity_model` (it's only reachable via `use_activity=False`).
Consider collapsing this into a single `activity_model: str`
parameter that accepts `"ideal"` alongside `"davies"`/`"sit"`,
removing the separate boolean gate. Touches the builder, the three
engines above, and their callers — worth scoping as its own small
phase rather than a drive-by fix.

## `demos/_bootstrap.py` — resolved 2026-09-16, deleted

Surfaced 2026-09-15 (see above): seven `model_api/` files still
imported `_bootstrap` with zero actual `vlmodels` use. Verified
2026-09-16 the same way — grep for `vlmodels` in each of the seven
found no real usage — and removed the dead `import _bootstrap` (plus
the now-pointless `sys.path.insert` scaffolding) from all seven. That
left two genuine holdouts that hadn't surfaced in the checkpoint-8
sweep: `demos/aerobic_fermentation_stoichiometry.ipynb` and
`demos/builder/batch_fermenter.ipynb`, both still importing
`vlmodels.fermenter.config.FermenterBuilder` — a module that no longer
exists (`models/vlmodels/` has no `fermenter/` subpackage any more;
`stirred-tank-template` moved it to
`PyOMES.templates.stirred_tank.StirredTankBuilder`, an identical
fluent API). Both notebooks were broken imports, not merely stale
comments. Migrated both to `StirredTankBuilder`, at which point
`demos/_bootstrap.py` had zero remaining importers and was deleted,
along with the stale doc pointers to it in `demos/README.md`.

## `docs/architecture.md` still describes deleted CUFermenter-era code

Surfaced 2026-09-15, same session as above. While fixing
`stirred-tank-template`'s checkpoint 9 (the doc directly referenced the
old `builder.py`/`factory.py` paths and `FermenterBuilder`/
`FermenterFactory` — fixed), found the doc still tags several *other*
`PyOMES/` subpackages "CUFermenter island" as if awaiting the
trigger-gated `CUFERMENTER_SUNSET` phase — `PyOMES/equilibria/
gl_equilibrium.py`, `PyOMES/control/{loops.py,system.py,controllers/,
actuators/,builders/}`, `PyOMES/sim/` — but `cufermenter-sunset`
already shipped 2026-06-01 (9,093 lines deleted, see
`../shipped/CUFERMENTER_SUNSET.md`). Worth checking whether these
paths still exist at all; if not, the whole "CUFermenter island"
framing throughout this doc is dead and should come out, not just be
re-pathed. Not touched — same reasoning as the README.md entry above,
needs its own pass rather than a line-by-line patch mid-migration.

## Four files mention FermenterBuilder/FermenterFactory only in prose, not imports

Surfaced 2026-09-15, same session. Not in `stirred-tank-template`'s
checkpoint 8 (demo imports) or 9 (doc cross-references) scope — these
reference the old names as comparison points in comments/docstrings,
not broken imports, so the code still runs correctly:
`demos/model_api/README.md`, `demos/model_api/chemistry/
reaction_system.py`, `demos/model_api/D2Cworkshop/{basic_layout,
updated_layout}/raw_construction.py`. Left untouched; sweep together
in a follow-up pass if one happens.

The other two originally on this list, `demos/usecases/
_generate_notebooks.py` and `docs/tutorials/_generate_notebooks.py`,
were fixed 2026-09-16 as a side effect of the `demos/_bootstrap.py`
cleanup above — both emit a "Where to go next" markdown cell into a
committed notebook (`demos/usecases/03_grow_ecoli_on_acetic_acid.ipynb`
and `docs/tutorials/03_cstr_dilution_rate_sweep.ipynb`) that pointed at
`models/vlmodels/fermenter/config/builder.py`; repointed to
`PyOMES.templates.stirred_tank.StirredTankBuilder` /
`PyOMES/templates/stirred_tank/builder.py` in both the generator source
and (by hand-editing the specific cell, not by rerunning the generator —
see note below) the committed notebooks. `docs/tutorials/
03_cstr_dilution_rate_sweep.ipynb` had the same stale pointer
independently (it isn't generated by either script) and was fixed the
same way.

**Trap found in the process:** rerunning `docs/tutorials/
_generate_notebooks.py` after a source fix is not safe to commit
as-is — `02_kinetic_co2_equilibration_microplate_well.ipynb` ships
with real executed outputs (plots, timing numbers) baked in, and a
full regenerate silently replaces them with the blank
`execution_count: null` / no-output template, plus at least one
unrelated string got mangled by existing replace logic in that script
("VLsim" in a §7 note became "PyOMES", breaking the sentence). Prefer
hand-editing the specific stale cell in the committed `.ipynb` (e.g.
via a notebook-aware editor that touches only that cell) over
regenerating from source when the target notebook may carry baked
outputs. `demos/usecases/03_grow_ecoli_on_acetic_acid.ipynb` and
`0_README.ipynb` were safe to regenerate (no baked outputs, confirmed
by a 1-line/no-op diff) — but that's a property of each notebook, not
the method, and needs re-checking per file. (`04_compare_runtime_by_usecase.ipynb`,
also confirmed safe at the time, was deleted entirely in the
`tutorials-reorg` phase — see `TUTORIALS_REORG_CHECKLIST.md`.)
