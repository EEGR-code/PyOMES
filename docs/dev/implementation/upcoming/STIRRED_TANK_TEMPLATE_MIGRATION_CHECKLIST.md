# Stirred-Tank Template Migration — Checklist

> Working checklist for the `stirred-tank-template` branch. Modelled on
> [`PHASE_KICKOFF_TEMPLATE.md`](PHASE_KICKOFF_TEMPLATE.md). Design doc:
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
- [ ] Branch created off current `main`: `git checkout -b stirred-tank-template`
- [ ] This checklist file committed on that branch as the first commit

## During

- [x] Plan/design doc exists:
      [`STIRRED_TANK_TEMPLATE_MIGRATION.md`](STIRRED_TANK_TEMPLATE_MIGRATION.md)
- [ ] Checkpoints tracked below as they land, one commit per checkpoint
- [ ] **If work stalls or is paused before shipping:** add a status
      banner to the top of the design doc immediately — what's built,
      what's tested, why it stopped, which branch/commit it's on.

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
- [ ] 3. **Move the implementation files** (four or five, per
      checkpoint 1) into `PyOMES/templates/stirred_tank/`, preserving
      history. Delete the two old `__init__.py` files outright. Sanity
      check: `models/vlmodels/fermenter/` no longer exists.
- [ ] 4. **Fix intra-package references** inside the moved files
      (`builder.py`, `factory.py`, `configs.py`, `kinetics.py`,
      `profiles.py`) — renames, deferred-rename docstring cleanup, lazy
      `PyOMES.core.simulation` import re-examined, real
      `stirred_tank/__init__.py` content written. Sanity check:
      `python -c "from PyOMES.templates.stirred_tank import StirredTankBuilder"`
      succeeds on its own.
- [ ] 5. **Fix the two functionally-dependent call sites**:
      `adm1/base.py:976` and `adm1/bsm2.py:804-805` (plus its
      `TransferConfig` import). Sanity check:
      `pytest tests/standalone/test_bsm2_reference.py -v` in isolation.
- [ ] 6. **Fix `vlmodels/__init__.py`**. Sanity check:
      `python -c "import vlmodels"`.
- [ ] 7. **Fix test imports** — `test_builder.py`, `test_configs.py`
      (import lines only), and the seven named classes in
      `test_simulation.py`. Sanity check:
      `pytest tests/standalone/test_builder.py tests/standalone/test_configs.py tests/standalone/test_simulation.py -v`.
- [ ] 8. **Fix demo imports** — the four `demos/builder/*.py` files
      plus `export_results.py`; update import line, drop
      `import _bootstrap` from all five. Sanity check: run each script
      directly, confirm it still prints its report.
- [ ] 9. **Fix remaining doc cross-references** —
      `PyOMES/core/recorder.py`'s stale docstring mention, `README.md`,
      `demos/README.md`, `demos/builder/README.md`.
- [ ] 10. **Full suite**: root `pytest`. Then ship per the
       branching/tagging convention (`--no-ff` merge, tag
       `stirred-tank-template-shipped`), moving the design doc + this
       checklist to `../shipped/` with a "Shipped" banner.

## Shipping

- [ ] Full test suite green on the branch
- [ ] `git checkout main`
- [ ] `git merge --no-ff stirred-tank-template -m "Merge stirred-tank-template: <summary>"`
- [ ] `git tag stirred-tank-template-shipped <commit-hash>`
- [ ] `git push && git push --tags`
- [ ] `git branch -d stirred-tank-template` and
      `git push origin --delete stirred-tank-template`
- [ ] Move `STIRRED_TANK_TEMPLATE_MIGRATION.md` + this checklist to
      `docs/dev/implementation/shipped/`, add a "Shipped" banner to both
- [ ] Update `docs/dev/implementation/upcoming/README.md`'s "Currently
      in flight" / "Recently shipped" lists
