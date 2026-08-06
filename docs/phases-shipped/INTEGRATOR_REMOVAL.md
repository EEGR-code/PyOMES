# Integrator Removal — Design Note

> **Status: Shipped 2026-05-20** on the `integrator-removal`
> branch. Tag: `integrator-removal-shipped`. See the companion
> [INTEGRATOR_REMOVAL_CHECKLIST.md](INTEGRATOR_REMOVAL_CHECKLIST.md)
> for the per-checkpoint implementation log, the §Inventory
> classification (52 callsites bucketed at checkpoint 1), and the
> BSM2 sentinel re-baseline (max drift S_h2 ~6e-3 relative,
> others 1e-7 to 5e-6).

## Status

Design resolved 2026-05-18 through exploration. Scope locked. Small
phase — single concern, narrow code surface (~150 LOC subtraction
across `src/`, `tests/`, `docs/`).

Implementation shipped through `integrator-removal` per the
[branching and tagging convention](../phases-upcoming/README.md#branching-and-tagging-convention) —
off `main`, merged `--no-ff`, tagged `integrator-removal-shipped`.

## What this is

Remove the two-protocol architecture for ODE integration. Today
the codebase has:

- `Integrator` Protocol (`step(rhs, t0, y0, dt) -> y_final`) in
  [src/reactions/integrators.py](../../src/reactions/integrators.py),
  with `RK4Integrator` and `EulerIntegrator` implementations. Used
  in exactly one place: `ControlVolume._integrate_reactions` inside
  the sequential `cv.advance()` body.
- `StepSolver` Protocol (`solve_step(cv, dt_h, chem_env, ...) ->
  AdvanceResult`) in
  [src/core/solvers.py](../../src/core/solvers.py), with
  `EulerSnapshotSolver` and `ScipyODESolver`. Used via
  `cv.advance(solver=...)`; replaces the entire sequential body
  when supplied. Both implementations bypass `Integrator` entirely.

This phase deletes the `Integrator` Protocol, both implementations,
the `integrator` kwarg on `ControlVolume.__init__`, and the
`cv.integrator` attribute. The sequential `cv.advance()` body's
reaction sub-step collapses from a 4-stage RK4 to a single Euler
evaluation. `StepSolver` becomes the sole user-facing extension
surface for time integration.

## Why

1. **The two-protocol architecture has no concrete user.**
   `Integrator` has one production caller. No third-party
   `Integrator` implementation exists in the repo or its
   dependants. No production model passes `integrator=` to
   `ControlVolume`. The Protocol was extracted during Phases 1–5
   to preserve numerical parity with the pre-refactor
   `_calc_ODE_with_headspace` RK4 step (the `RK4Integrator`
   docstring still says so verbatim). That parity goal is long
   since met.

2. **RK4 inside Lie-Trotter is wasted accuracy.** The sequential
   body is a first-order operator split (feed → reactions →
   transfer, each evaluated from previous-sub-step state). RK4
   inside the reactions sub-step gives accurate substrate
   depletion within `dt_h`, but the outer splitting is O(dt)
   regardless — global error is dominated by splitting
   cross-terms, not by the reaction sub-step's internal accuracy.
   Callers who care about accuracy reach for `ScipyODESolver`
   (adaptive, with stiff methods available), not for the RK4
   inside the sequential body. See
   [docs/phases-shipped/ORDERING_CRITIQUE.md](../phases-shipped/ORDERING_CRITIQUE.md)
   Review 2 for the operator-splitting discussion.

3. **The split is non-standard versus the wider literature.**
   Canonical bioprocess simulators (Aspen, gPROMS, WEST, BSM2
   MATLAB reference, PyADM1) all integrate monolithically with a
   stiff DAE solver — one integrator, one RHS, no inner/outer
   protocol split. Reactive-transport codes (PHREEQC, CrunchFlow)
   use operator splitting but with **same-order** inner solvers
   per sub-step. The mixed-order character (RK4 inside, Euler
   outside) is bespoke to VLsim and corresponds to no published
   scheme. After removal, `ScipyODESolver` is the standard
   monolithic-DAE path and `EulerSnapshotSolver` is textbook
   forward Euler — both are recognisable to anyone reading the
   code.

4. **Future readers pay the dead-code cost.** Any phase that
   touches `cv.advance()` or the time-loop has to consider the
   `Integrator` indirection. `SIMULATION_CLASS` is the next phase
   to do exactly this. Removing the indirection first (or as a
   companion phase) keeps the future surface honest.

## Scope — concrete deletions

- [src/reactions/integrators.py](../../src/reactions/integrators.py)
  — delete entirely (~100 LOC).
- [src/reactions/__init__.py](../../src/reactions/__init__.py) —
  remove `Integrator`, `RK4Integrator`, `EulerIntegrator` exports
  and module-level docstring bullet.
- [src/core/control_volume.py](../../src/core/control_volume.py):
  - `integrator` kwarg on `__init__` — delete.
  - `self.integrator` attribute and default `RK4Integrator()`
    initialisation — delete.
  - `_integrate_reactions` — rewrite the stepping call from
    `self.integrator.step(rhs, 0.0, y0, dt_h)` to a single Euler
    evaluation: `y_final = y0 + dt_h * rhs(0.0, y0)`. Species
    vector map and `rhs` construction unchanged.
  - `snapshot()` — remove the `integrator=self.integrator`
    argument.
- [tests/standalone/test_reactions.py:540-556](../../tests/standalone/test_reactions.py#L540-L556)
  — delete the RK4Integrator/EulerIntegrator unit-test block.
  Other tests in this file remain.
- [docs/architecture.md:31](../../docs/architecture.md#L31) —
  remove the "Integrator (optional)" bullet from the
  ControlVolume description.
- [docs/architecture.md:283](../../docs/architecture.md#L283) —
  remove `integrators.py` from the directory tree.
- [docs/class_diagrams.md:219-226](../../docs/class_diagrams.md#L219-L226)
  — remove the Integrator / RK4Integrator / EulerIntegrator class
  block.
- [docs/solvers.md](../../docs/solvers.md) — update the "RK4 by
  default" mention in the `ControlVolume.advance()` description
  and any other references to inner reaction integration.

## Sentinel risk

The change shifts numerical behaviour for any caller that:

1. Runs `cv.advance()` **without** an explicit `solver=`
   argument, **and**
2. Has a non-trivial `reaction_model`.

Pre-change: reactions advance with RK4 over `[0, dt_h]`
(4 RHS evaluations). Post-change: reactions advance with one Euler
step (1 RHS evaluation × `dt_h`). Drift magnitude scales with how
much reaction rates change with species concentrations within one
`dt_h` — Monod kinetics with substantial substrate depletion per
step will see the largest drift.

**Known sentinel callsite — re-baseline expected:**

- [tests/standalone/test_bsm2_reference.py:119](../../tests/standalone/test_bsm2_reference.py#L119)
  — BSM2 golden trajectory uses bare `cv.advance(dt_h, chem_env)`.
  Despite [docs/solvers.md](../../docs/solvers.md)'s recommendation
  to use `ScipyODESolver(freeze_speciation=True, method="Radau")`
  for "reproducing BSM2 / PyADM1 validation results" (which
  describes a *user-facing* recommendation), the in-repo test
  exercises the default sequential body. Re-baseline the golden
  values, document the drift inline in the test file (as previous
  chemistry-unification phases did when sentinels shifted), and
  keep `RTOL_SENTINEL` at its current `1e-9`.

**Known unaffected callsites:**

- [models/vlmodels/fermenter/config/factory.py:589](../../models/vlmodels/fermenter/config/factory.py#L589)
  — passes `solver=solver` (`EulerSnapshotSolver` or
  `ScipyODESolver` per config). Bypasses `_integrate_reactions`.
  All four fermenter demos route through this factory.
- BSM2 production model itself — only the *test* uses bare
  `cv.advance()`; the model definition in
  [models/vlmodels/adm1/bsm2.py](../../models/vlmodels/adm1/bsm2.py)
  is solver-agnostic and the production caller would pass
  `ScipyODESolver`.

**Investigate during implementation:**

- [models/vlmodels/fermenter/unit.py:2217](../../models/vlmodels/fermenter/unit.py#L2217)
  — legacy `CUFermentationSpeciation` path uses bare
  `_cv.advance(dt_h, chem_env)`. This class is marked out-of-scope
  for `SIMULATION_CLASS` migration; any sentinel drift here is
  bounded to the legacy unit and is not blocking.
- [tests/standalone/test_cv_advance.py](../../tests/standalone/test_cv_advance.py)
  — many tests call bare `cv.advance()` but most are structural
  (presence of fields, basic mass-balance). Classify each test
  during checkpoint 1: "structural — no shift expected" vs
  "numerical — shift expected, may need re-baseline."
- [tests/standalone/test_multi_cv.py](../../tests/standalone/test_multi_cv.py)
  and [tests/standalone/test_gas_liquid_link.py](../../tests/standalone/test_gas_liquid_link.py)
  use `sys.advance_all(dt_h)` without `solvers=`, which routes to
  bare `cv.advance()` per CV. Same classification needed.

## Sketched checkpoints

1. **Inventory** — grep all `cv.advance(` callsites, classify
   each as "structural / no shift" vs "numerical / re-baseline."
   Commit the inventory as a note in the checklist file.
2. **Code deletion** — delete `integrators.py`, the
   `__init__.py` exports, the `integrator` kwarg + attribute, and
   rewrite `_integrate_reactions` to one-step Euler. Commit.
3. **Test deletion** — remove the RK4Integrator/EulerIntegrator
   unit-test block. Commit.
4. **Run the standalone suite** — collect all numerical
   sentinel failures. Each failure should be attributable to a
   callsite from step 1's "numerical / re-baseline" list. If a
   failure shows up in the "structural / no shift" list, stop and
   investigate — it indicates the change has a wider blast radius
   than expected.
5. **Re-baseline BSM2 golden** — update
   `test_bsm2_reference.py` golden values; document drift inline.
6. **Re-baseline other numerical tests** — case-by-case,
   inline-documented.
7. **Docs sweep** — `architecture.md`, `class_diagrams.md`,
   `solvers.md` updates.
8. **Final test run + tag.**

## Relationship to other phases

- **[SIMULATION_CLASS.md](SIMULATION_CLASS.md)** preserves
  bit-for-bit sentinel parity ("the change is structural, not
  numerical" per
  [[project_simulation_class]]). Bundling
  integrator removal into `SIMULATION_CLASS` would conflate two
  different sentinel-shift sources (structural reorganisation +
  numerical scheme change), making post-merge drift hard to
  attribute. Keep them separate.
- **Ordering recommendation:** Ship `SIMULATION_CLASS` first
  (structural, sentinels unchanged), then `integrator-removal`
  (numerical, sentinels re-baselined). The reverse order also
  works but inverts the cleaner causal story for any future
  bisect against the BSM2 golden test.
- **[CONTAINER_LAYERING.md](CONTAINER_LAYERING.md)** and
  **[CHEMISTRY_UNIFICATION_PLAN.md](CHEMISTRY_UNIFICATION_PLAN.md)
  Phase 3b** are independent of this work.

## Out of scope

- Replacing `EulerSnapshotSolver`'s explicit-Euler scheme with
  anything else. Different concern.
- Changing the operator-split ordering (feed → reactions →
  transfer). That's
  [docs/phases-shipped/ORDERING.md](../phases-shipped/ORDERING.md)
  territory and not on the table.
- Removing the `_integrate_reactions` method itself. The method
  remains; only its inner stepping collapses to a single Euler
  evaluation. (A future phase could fold the reaction sub-step
  directly into the sequential body without a dedicated method,
  but that's a separate readability win, not scope here.)
- Touching `ScipyODESolver` — it's the standard monolithic-DAE
  pattern, no changes needed.
- Touching the `StepSolver` Protocol itself. After this phase it
  becomes the sole time-integration extension surface, but its
  shape is already correct.

## Trigger

Scheduled. Pick up after `SIMULATION_CLASS` ships. The phase is
small (1–2 days) and self-contained; running it as a follow-on
keeps the SIMULATION_CLASS golden-test discipline intact and
sequences the two kinds of sentinel shift cleanly.
