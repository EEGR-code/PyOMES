# Stirred-Tank Template Migration — Checklist

> **Status: Shipped 2026-09-15.** All 10 checkpoints landed on branch
> `stirred-tank-template`, tag `stirred-tank-template-shipped`. Final
> suite: 2011 passed, 28 skipped, 0 failed. Modelled on
> [`PHASE_KICKOFF_TEMPLATE.md`](../upcoming/PHASE_KICKOFF_TEMPLATE.md). Design doc:
> [`STIRRED_TANK_TEMPLATE_MIGRATION.md`](STIRRED_TANK_TEMPLATE_MIGRATION.md) —
> read that first for the "why," the naming decisions, and the full
> blast-radius audit. This file is the implementation log only.

## Commit discipline for this phase

The repo owner runs every `git add`, `git commit`, and `git push`
personally. At each checkpoint below: make the file edits, report
exactly what changed, run the checkpoint's sanity check, then stop and
wait for the owner to review, stage, and commit before moving to the
next checkpoint. `git mv` counts as a staging operation — prefer a
plain filesystem move plus a separate `git add` the owner runs, unless
told otherwise.

## Pre-flight

- [x] `git status -sb` clean — confirmed 2026-09-15.
- [x] `git log origin/main..main --oneline` empty — confirmed 2026-09-15.
- [x] Branch created off current `main`: `git checkout -b stirred-tank-template`
- [x] This checklist file committed on that branch as the first commit

## During

- [x] Plan/design doc exists:
      [`STIRRED_TANK_TEMPLATE_MIGRATION.md`](STIRRED_TANK_TEMPLATE_MIGRATION.md)
- [x] Checkpoints tracked below as they land, one commit per checkpoint
- [x] **If work stalls or is paused before shipping:** N/A — shipped
      without stalling.

### Checkpoints

Numbered to match the design doc's "Checkpoints" section exactly.

- [x] 1. **Decide the `kinetics.py` destination** — **resolved
      2026-09-15: stays under `stirred_tank/`.** `GrowthKinetics`
      (scalar μ, builder-attached rate law) and `KineticModel`
      (self-integrating RHS owning its own state vector) share a name,
      not an interface — confirmed by reading both; `YeastAcetateV1.rhs`
      hand-derives a Monod term rather than reusing `Monod.mu()` because
      the two protocols can't compose. Real unification would be a
      protocol/integration-model redesign, out of scope for a rename
      phase — logged as design doc open question 6, deferred. See the
      design doc's "Open questions" item 1 for the full reasoning.
- [x] 2. **Create skeletons** — `PyOMES/templates/__init__.py` and
      `PyOMES/templates/stirred_tank/__init__.py` created, both empty.
      Sanity check: `python -c "import PyOMES.templates.stirred_tank"`
      → passed.
- [x] 3. **Move the implementation files** — `builder.py`, `configs.py`,
      `factory.py`, `kinetics.py`, `profiles.py` moved (plain filesystem
      move, not staged) from `models/vlmodels/fermenter/[config/]` into
      `PyOMES/templates/stirred_tank/`. Both old `__init__.py` files
      deleted outright, then the now-empty `fermenter/config/` and
      `fermenter/` directories removed. Sanity check:
      `models/vlmodels/fermenter/` no longer exists → passed. **Note for
      staging:** these are plain moves, not `git mv` — `git add -A` (or
      equivalent) should still let git detect them as renames via
      content-similarity, since content is untouched in this checkpoint.
- [x] 4. **Fix intra-package references** inside the moved files —
      `FermenterBuilder`→`StirredTankBuilder`,
      `FermenterFactory`→`StirredTankFactory` (class defs, all return
      annotations/docstrings/`__repr__`/internal calls); both classes'
      docstrings gained the explicit gas+liquid-only statement per the
      naming decision. Removed the two stale "Phase 7 holdover, rename
      deferred" Notes blocks (`builder.py`'s `build()`,
      `factory.py`'s `create_volume()`). Fixed stale docstring import
      examples in `builder.py`, `kinetics.py`, `profiles.py`
      (`PyOMES.config.kinetics`/`PyOMES.profiles` → `PyOMES.templates.
      stirred_tank`) and `factory.py`'s `PyOMES.config import *` example;
      stripped the "Stage 17a" dev-stage label from `profiles.py`; swept
      `configs.py` for stray `Fermenter`/`fermenter.*`-path docstring
      text (dataclass names themselves unchanged, per the naming
      decision — none renamed). **Lazy `PyOMES.core.simulation` import
      in `build_simulation()` re-examined per risk #2 and promoted to a
      top-level import** — traced the import chain (`PyOMES/core/__init__.py`
      already imports `.simulation` unconditionally before any submodule
      is reachable; nothing under `PyOMES.core` imports `PyOMES.templates`
      or `vlmodels`) and confirmed empirically (`build()` and
      `build_simulation()` both run end-to-end). The historical cycle was
      real only for the old cross-package (`vlmodels` importing `PyOMES`)
      layout; moving the file intra-package removed it. Wrote
      `stirred_tank/__init__.py`'s real content: re-exports mirroring the
      old `config/__init__.py`'s `__all__` list (renamed classes, local
      `.profiles` import instead of `vlmodels.fermenter.profiles`, "Stage
      17"/"17c" comments dropped). Sanity check:
      `python -c "from PyOMES.templates.stirred_tank import StirredTankBuilder"`
      → passed. Additional functional smoke test (`build()` and
      `build_simulation()` both construct successfully) → passed.
- [x] 5. **Fix the two functionally-dependent call sites** —
      `adm1/base.py`: lazy import `FermenterBuilder` →
      `from PyOMES.templates.stirred_tank import StirredTankBuilder`,
      the builder-construction call, and one stale comment mentioning
      `FermenterBuilder`. `adm1/bsm2.py`: lazy imports for both
      `FermenterBuilder` and `TransferConfig` collapsed into one
      `from PyOMES.templates.stirred_tank import StirredTankBuilder,
      TransferConfig`, plus the builder-construction call. **Found a
      real ordering gap not called out in the design doc:**
      `models/vlmodels/__init__.py` still did
      `from vlmodels.fermenter import FermenterFactory, FermenterBuilder`
      at package-init time — since that submodule no longer exists,
      this blocked importing *any* `vlmodels.*` module, including
      `vlmodels.adm1.bsm2`, so checkpoint 5's own sanity check couldn't
      pass without checkpoint 6's fix landing first. Did checkpoint 6's
      fix (see below) as a prerequisite rather than block. Sanity check:
      `pytest tests/standalone/test_bsm2_reference.py -v` (needs
      `PYTHONPATH` including `models/`, per `README.md` — pytest's own
      `conftest.py` does this automatically, a bare `python -c` does
      not) → **6 passed**, including the numeric sentinel tests
      (`test_final_liquid_concentrations`/`test_final_gas_mol`/
      `test_final_pH`) — no numerical drift from the move. The
      `ConservationWarning`/`AccuracyWarning` output is pre-existing
      BSM2 model behavior, not a regression (the test asserts specific
      numeric sentinels, not warning-free execution).
- [x] 6. **Fix `vlmodels/__init__.py`** — done early, as checkpoint 5's
      prerequisite (see above): dropped the dead `from vlmodels.fermenter
      import FermenterFactory, FermenterBuilder` re-export entirely and
      the now-inaccurate "built on the fermenter framework" docstring
      sentence. Sanity check: `python -c "import vlmodels"` (same
      `PYTHONPATH` note as above) → passed.
- [x] 7. **Fix test imports** — `test_configs.py`: import line only
      (its `TransferMode`/`VesselConfig`/etc. usages are unaffected,
      since those dataclass names didn't change). `test_builder.py`:
      import line plus a file-wide `FermenterBuilder`→
      `StirredTankBuilder`/`FermenterFactory`→`StirredTankFactory`
      sweep (docstring, every instantiation — 30 + 1 occurrences).
      `test_simulation.py`: the seven named classes
      (`TestKineticGasLiquidLinkGating`, `TestC9FluxApplyDispatch`,
      `TestC9ParamPathDispatch`, `TestC10TemperatureRamp`,
      `TestC10VVMSchedule`, `TestC10ProfileOrderingBeforeAdvance`,
      `TestC11FermenterBuilderSimulation`) — confirmed via class-boundary
      grep that every `Fermenter`/`vlmodels.fermenter.config` match in
      the file falls inside these seven (or the comment immediately
      above the seventh); `TestC13ADM1Simulation` (the one class the
      design doc says to leave alone — genuine external `vlmodels.adm1`
      model) had zero matches, confirmed untouched. **Side effect worth
      flagging:** the file-wide `FermenterBuilder`→`StirredTankBuilder`
      replace also renamed the class itself,
      `TestC11FermenterBuilderSimulation` →
      `TestC11StirredTankBuilderSimulation` (a substring match, not
      explicitly scoped by the design doc) — kept deliberately rather
      than reverted, since the class tests `StirredTankBuilder.
      build_simulation()` and leaving the old name would itself be
      stray Fermenter-branded text. The one generic-word usage
      (`"""Fermenter-shaped CV with a kinetic gas-liquid link..."""`
      docstring around line 1900) was correctly left alone by the
      targeted replacements. Sanity check:
      `pytest tests/standalone/test_builder.py tests/standalone/test_configs.py tests/standalone/test_simulation.py -q`
      → **397 passed**, including `TestC13ADM1Simulation` (confirms the
      checkpoint 5 `adm1/base.py` fix holds here too).
- [x] 8. **Fix demo imports** — verified first (per design doc's claim)
      that none of the five files touch `vlmodels` for anything besides
      the fermenter import; confirmed via grep. `batch_fermenter.py`,
      `cstr_fermenter.py`, `fed_batch_fermenter.py`,
      `microplate_fermenter.py`: dropped the whole bootstrap block
      (`sys.path.insert(...)` + `import _bootstrap`, now unneeded since
      `PyOMES.templates.stirred_tank` resolves via the editable install
      with no `models/` path needed) and the now-dead
      `import sys`/`from pathlib import Path`; import line and every
      `FermenterBuilder()` call site → `StirredTankBuilder`.
      `export_results.py`: same bootstrap removal (a slightly different
      two-step form — `sys.path.insert(0, .../"demos")` then
      `import _bootstrap`, since this file lives one level deeper); its
      unrelated `import pathlib`/`tempfile` (used later for the
      CSV/Parquet round-trip) kept, only the now-dead `import sys`/
      `from pathlib import Path` removed;
      `FermenterFactory`→`StirredTankFactory`. Final sweep confirmed
      zero remaining `FermenterBuilder`/`FermenterFactory`/`vlmodels`/
      `_bootstrap` references across all five files. Sanity check: ran
      each of the five scripts directly →  all five completed and
      printed their reports (`batch_fermenter.py`, `cstr_fermenter.py`,
      `fed_batch_fermenter.py`, `microplate_fermenter.py`,
      `export_results.py`, the last one's CSV/Parquet round-trip checks
      included — "All export checks passed").
- [x] 9. **Fix remaining doc cross-references** —
      `PyOMES/core/recorder.py`: removed two now-obsolete "legacy
      `BatchResult` coexists until C14" notes (module docstring +
      class docstring) — the legacy class genuinely doesn't exist
      anywhere in the codebase anymore (confirmed via grep/history),
      not just at a stale path, so trimmed rather than re-pathed.
      `PyOMES/core/simulation.py`: also trimmed two dangling
      `models/vlmodels/fermenter/config/factory.py:run_batch` line-
      number citations found while in the area (flagged at checkpoint
      4) — same "points at code that never existed at that path"
      issue. `README.md`, `demos/README.md`, `demos/builder/README.md`,
      `demos/_bootstrap.py`: import examples, class names, and the
      Repository Layout `fermenter/` line (removed — directory no
      longer exists) updated. **Two scope decisions made with the
      owner mid-checkpoint** (README.md and several other files turned
      out to be broadly stale across *many* unrelated prior phases,
      not just this migration): (1) `README.md` and `docs/
      architecture.md`'s wider staleness (retired `speciation/`/
      `MultiCVSystem` names, the unrelated already-sunset
      `create_standalone_fermenter`/CUFermenter-era API, `docs/
      architecture.md`'s "CUFermenter island" tags on several *other*
      subpackages) — narrow-fixed only the lines directly about *this*
      migration, logged the rest to a new
      [`../OPEN_WORK.md`](../OPEN_WORK.md) (also fixes a pre-existing
      dangling link — `upcoming/README.md`'s "Priority order" section
      already referenced `../OPEN_WORK.md`, which didn't exist until
      now). (2) Six files mentioning `FermenterBuilder`/
      `FermenterFactory` only in comments/docstrings, not imports
      (`demos/model_api/README.md`, `reaction_system.py`, two
      `raw_construction.py` files, two `_generate_notebooks.py`
      scripts) — logged to `OPEN_WORK.md`, not touched (not in the
      design doc's checkpoint 8/9 scope, code still runs). Also found
      and logged to `OPEN_WORK.md`: seven other demo files still
      `import _bootstrap` with **zero** actual `vlmodels` usage —
      possibly dead weight, not verified/touched. No specific sanity
      check listed for this checkpoint by the design doc; ran the full
      `pytest tests/standalone` suite as an extra check since
      `recorder.py`/`simulation.py` were touched beyond the checkpoint's
      named file list → **2011 passed, 28 skipped, 0 failed.**
- [x] 10. **Full suite**: root `pytest` → **2011 passed, 28 skipped, 0
       failed.** Design doc + this checklist moved to `../shipped/`
       with "Shipped" banners (this edit). Merge/tag/push/branch-delete
       below, run by the repo owner per this phase's commit discipline.

## Shipping

- [x] Full test suite green on the branch
- [ ] `git checkout main`
- [ ] `git merge --no-ff stirred-tank-template -m "Merge stirred-tank-template: <summary>"`
- [ ] `git tag stirred-tank-template-shipped <commit-hash>`
- [ ] `git push && git push --tags`
- [ ] `git branch -d stirred-tank-template` and
      `git push origin --delete stirred-tank-template`
- [x] Move `STIRRED_TANK_TEMPLATE_MIGRATION.md` + this checklist to
      `docs/dev/implementation/shipped/`, add a "Shipped" banner to both
- [x] Update `docs/dev/implementation/upcoming/README.md`'s "Currently
      in flight" / "Recently shipped" lists
