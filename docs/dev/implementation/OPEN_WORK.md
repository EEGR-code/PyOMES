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

## `demos/_bootstrap.py` may be dead weight outside the builder demos

Surfaced 2026-09-15, same session as above. `stirred-tank-template`'s
checkpoint 8 confirmed and removed `import _bootstrap` from the five
demos that used it for `vlmodels.fermenter.config` (now
`PyOMES.templates.stirred_tank`, resolved via the editable install —
no `models/` path needed). Seven other files still `import _bootstrap`
(`model_api/chemistry/{reaction_system,partition_model}.py`,
`model_api/chemistry/fba/{fba_toy,fba_ecoli_core}.py`,
`model_api/D2Cworkshop/{basic_layout,updated_layout}/
raw_construction.py`, `model_api/solver_comparison.py`) — a grep for
`vlmodels` in each of the seven found **zero** actual `vlmodels`
imports. Worth checking whether the import is genuinely load-bearing
in any of them (dynamic/indirect import? a real need not caught by a
plain grep?) or whether it's leftover from an earlier layout and can
be dropped file-by-file. Not investigated further or touched — outside
this migration's scope, and each file needs its own verification
before assuming the import is safe to remove.

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

## Six files mention FermenterBuilder/FermenterFactory only in prose, not imports

Surfaced 2026-09-15, same session. Not in `stirred-tank-template`'s
checkpoint 8 (demo imports) or 9 (doc cross-references) scope — these
reference the old names as comparison points in comments/docstrings,
not broken imports, so the code still runs correctly:
`demos/model_api/README.md`, `demos/model_api/chemistry/
reaction_system.py`, `demos/model_api/D2Cworkshop/{basic_layout,
updated_layout}/raw_construction.py`, `demos/usecases/
_generate_notebooks.py`, `docs/tutorials/_generate_notebooks.py`. The
two `_generate_notebooks.py` scripts also emit this text into their
generated, committed notebooks — a cosmetic staleness, not a breakage.
Left untouched; sweep together in a follow-up pass if one happens.
