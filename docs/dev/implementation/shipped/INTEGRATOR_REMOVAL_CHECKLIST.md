# Integrator Removal — Checklist

> **Status: Shipped 2026-05-20** — all 8 checkpoints landed on
> the `integrator-removal` branch (7 commits, one per checkpoint
> after C0; C6 was a documented no-op).
> **816 standalone tests passing**, down from 821 baseline by 5
> deleted tests (4 RK4/Euler unit tests in C3 + 1
> snapshot-integrator test in C2). BSM2 sentinels re-baselined at
> C5: largest drift S_h2 at ~6e-3 relative, others 1e-7 to 5e-6,
> pH 1.6e-7 relative. Tag: `integrator-removal-shipped`.

Implementation log for `integrator-removal`. The design framing —
why this phase existed, what got deleted, what numerical drift to
expect, the relationship to `STATE_UNIFICATION` and
`SIMULATION_CLASS` — lives in
[INTEGRATOR_REMOVAL.md](INTEGRATOR_REMOVAL.md). This file is the
file-level edit log, the §Inventory classification (52 callsites
bucketed at C1), and the §Test run record (3 failures at C4, all in
`test_bsm2_reference.py`; re-baselined at C5). Modelled on
[CHEMISTRY_UNIFICATION_4_CHECKLIST.md](CHEMISTRY_UNIFICATION_4_CHECKLIST.md).

## Goal

Delete the two-protocol architecture for ODE integration:

- The `Integrator` Protocol, `RK4Integrator`, and `EulerIntegrator`
  in [`src/reactions/integrators.py`](../../src/reactions/integrators.py)
  disappear entirely (~100 LOC).
- `ControlVolume.__init__` loses its `integrator` kwarg and the
  `self.integrator` attribute.
- `ControlVolume._integrate_reactions` keeps its method shape but
  collapses its inner stepping call
  `self.integrator.step(rhs, 0.0, y0, dt_h)` to a single forward
  Euler evaluation `y0 + dt_h * rhs(0.0, y0)`.
- `StepSolver` (`EulerSnapshotSolver`, `ScipyODESolver`) is left
  untouched and becomes the sole user-facing extension surface for
  time integration.

After this branch ships:

- [`src/reactions/integrators.py`](../../src/reactions/integrators.py)
  is gone. `src/reactions/__init__.py` no longer exports `Integrator`,
  `RK4Integrator`, or `EulerIntegrator`; the module-level docstring
  bullet is removed.
- [`src/core/control_volume.py`](../../src/core/control_volume.py)
  has no `integrator` kwarg, no `self.integrator`, and no
  `integrator=self.integrator` plumbing in `snapshot()`. The
  reaction sub-step inside the sequential `cv.advance()` body is one
  Euler evaluation.
- The RK4Integrator/EulerIntegrator unit-test block in
  [`tests/standalone/test_reactions.py`](../../tests/standalone/test_reactions.py)
  is deleted. Other tests in that file remain.
- [`tests/standalone/test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py)
  golden trajectory values are re-baselined; drift is documented
  inline at the test's location per the chemistry-unification
  precedent. `RTOL_SENTINEL` stays at `1e-9`.
- [`docs/architecture.md`](../../docs/architecture.md),
  [`docs/class_diagrams.md`](../../docs/class_diagrams.md), and
  [`docs/solvers.md`](../../docs/solvers.md) drop the
  Integrator/RK4 references; the directory tree in `architecture.md`
  loses `integrators.py`.

## Out of scope

Reaffirmed from [INTEGRATOR_REMOVAL.md](INTEGRATOR_REMOVAL.md)
"Out of scope":

- **Touching `StepSolver` or its implementations.** After this
  phase it is the sole time-integration extension surface, but its
  shape is already correct.
- **Replacing `EulerSnapshotSolver`'s explicit-Euler scheme with
  anything else.** Different concern.
- **Changing the operator-split ordering** (feed → reactions →
  transfer). That's `docs/shipped/ORDERING.md` territory.
- **Removing `_integrate_reactions` itself.** The method stays; only
  its inner stepping collapses. Folding the method directly into
  the sequential body is a separate readability win.
- **Migrating the legacy `CUFermentationSpeciation` path** in
  [`models/vlmodels/fermenter/unit.py`](../../models/vlmodels/fermenter/unit.py)
  beyond observing whatever drift it experiences. That class is
  marked out-of-scope for `SIMULATION_CLASS`; any drift here is
  bounded to the legacy unit and is not blocking.

## Resolved decisions

Pinned in [INTEGRATOR_REMOVAL.md](INTEGRATOR_REMOVAL.md) (2026-05-18,
sequencing revised 2026-05-20):

- **Sequencing: this phase ships first.** Originally scheduled
  after `SIMULATION_CLASS`. Resequenced 2026-05-20 to ship ahead
  of `STATE_UNIFICATION` and `SIMULATION_CLASS`. The ~150 LOC
  subtraction lands cleanly on the current code surface; deferring
  it past `STATE_UNIFICATION` would mean re-doing the
  `_integrate_reactions` rewrite against whatever shape
  STATE_UNIFICATION leaves behind. Cheaper to ship now.

- **Sentinel re-baseline expected for `test_bsm2_reference.py`.**
  Bare `cv.advance(dt_h, chem_env)` (no explicit `solver=`) is the
  only known production-shaped callsite that will shift
  numerically. Drift magnitude scales with how much reaction rates
  change with species concentrations across one `dt_h` — Monod
  kinetics with substantial substrate depletion per step show the
  largest drift. Re-baseline values, document drift inline,
  preserve `RTOL_SENTINEL = 1e-9`.

- **Other `cv.advance()` callsites: classify, then run.** Most
  tests calling bare `cv.advance()` are structural (presence of
  fields, basic mass balance) and unaffected. A minority are
  numerical and will need re-baseline. The classification (output
  of checkpoint 1, embedded below in §Inventory) gates the
  interpretation of checkpoint 4's test-run failures: any failure
  in the "structural / no shift" bucket means the change has a
  wider blast radius than expected, stop and investigate.

- **`StepSolver` left untouched.** No behavioural change to
  `EulerSnapshotSolver` or `ScipyODESolver`. The fermenter factory
  ([`models/vlmodels/fermenter/config/factory.py:589`](../../models/vlmodels/fermenter/config/factory.py#L589))
  passes `solver=solver` explicitly, so all four fermenter demos
  route around `_integrate_reactions` and see no drift.

- **`_integrate_reactions` method shape preserved.** The method
  signature, species vector map, and `rhs` construction stay
  identical; only the line
  `y_final = self.integrator.step(rhs, 0.0, y0, dt_h)` collapses to
  `y_final = y0 + dt_h * rhs(0.0, y0)`. Folding the method directly
  into `advance()` is left as a future readability pass.

## Checkpoints

The checkpoints are ordered so each leaves the test suite in a
runnable state. The inventory lands first (1) so checkpoint 4's
interpretation has a key. Code deletion + RHS collapse land next
(2). The targeted test-block deletion follows (3). The full suite
run (4) drives any re-baselines (5, 6). Docs (7) and ship (8) close
out.

**Commit discipline:** commit immediately after each checkpoint
finishes (suite green for that checkpoint's slice). Mid-flight edits
were silently lost between checkpoints during
`chemistry-unification-3`; per-checkpoint commits are load-bearing,
not optional.

### 1. Inventory — classify every bare `cv.advance(` callsite

- [x] Grep all `cv.advance(`, `_cv.advance(`, `.advance(` callsites
      across `src/`, `tests/`, `models/`, and `demos/`.
- [x] For each callsite, classify as one of:
      - **Bypasses `_integrate_reactions` (no shift):** explicit
        `solver=` argument, regardless of `solver` type.
      - **Structural / no shift expected:** uses bare
        `cv.advance(dt_h, chem_env)` but asserts only on
        non-numerical properties (presence of fields, exception
        behaviour, basic mass conservation under simple kinetics).
      - **Numerical / re-baseline expected:** uses bare
        `cv.advance(dt_h, chem_env)` and asserts on golden values
        or tight tolerances against pre-change RK4 trajectories.
      - **Multi-CV routing:** `sys.advance_all(dt_h)` without
        `solvers=` — classify per the CV's reaction_model load.
- [x] Record the inventory inline in this checklist under
      §Inventory below. (See §Inventory section for the
      per-file/per-line breakdown.)
- [x] Document any callsite where the classification is genuinely
      ambiguous; resolve at checkpoint 4 by reading the test
      assertions rather than guessing now. (None ambiguous —
      every callsite has a clear bucket. Bordering case is
      `test_membrane.py:145` where steady-state shape over 2000
      reaction-active steps could in principle shift; assertion
      is on tail variance not absolute values, classified
      structural.)

Sanity check: the §Inventory section is populated. Every callsite
known from [INTEGRATOR_REMOVAL.md](INTEGRATOR_REMOVAL.md)'s
"Known sentinel callsite", "Known unaffected callsites", and
"Investigate during implementation" lists is accounted for.

### 2. Code deletion + RHS collapse

- [ ] Delete
      [`src/reactions/integrators.py`](../../src/reactions/integrators.py)
      entirely.
- [ ] [`src/reactions/__init__.py`](../../src/reactions/__init__.py):
      remove the
      `from .integrators import RK4Integrator, EulerIntegrator, Integrator`
      line, the three entries from `__all__`, and the
      "Integrators" bullet from the module-level docstring.
- [ ] [`src/core/control_volume.py`](../../src/core/control_volume.py):
      - Remove the `integrator: Optional[Any] = None` kwarg from
        `__init__`.
      - Remove the `self.integrator = ...` initialisation block
        (incl. the lazy `from ..reactions.integrators import
        RK4Integrator` import).
      - Remove the `integrator` entry from the `__init__` docstring's
        Parameters block.
      - In `_integrate_reactions`, replace
        `y_final = self.integrator.step(rhs, 0.0, y0, dt_h)` with
        `y_final = y0 + dt_h * np.asarray(rhs(0.0, y0), dtype=float)`.
      - Remove the `integrator=self.integrator` argument from
        `snapshot()`.
      - Update the `_integrate_reactions` docstring's "Uses the CV's
        integrator (default RK4) ..." paragraph to describe the new
        single-Euler step honestly.
- [ ] Check that no other src/ file imports `Integrator`,
      `RK4Integrator`, or `EulerIntegrator`. Grep before committing.

Sanity check: `python -c "from VLsim.reactions import RK4Integrator"`
raises `ImportError`. `python -c "import VLsim"` succeeds. A
trivial `ControlVolume()` construction with no `reaction_model`
runs without errors.

### 3. Delete the RK4/Euler unit-test block

- [ ] [`tests/standalone/test_reactions.py`](../../tests/standalone/test_reactions.py)
      lines around 540–556: delete the RK4Integrator/EulerIntegrator
      unit-test block. Verify by greping for `RK4Integrator` and
      `EulerIntegrator` in the file — should return nothing after
      the deletion. Other tests in this file remain.
- [ ] Grep `tests/` more broadly for any other import of
      `RK4Integrator`, `EulerIntegrator`, or `Integrator` from
      `VLsim.reactions` — fix or delete each one. Likely only the
      block above, but verify.

Sanity check:
`python -m pytest tests/standalone/test_reactions.py -q` passes.
The file still has its non-integrator tests; the test count for
this one file drops by however many tests were in the deleted block.

### 4. Run the standalone suite — categorise failures

- [x] Run `python -m pytest tests/standalone -q` and collect the full
      list of failures.
- [x] Cross-reference each failure with the §Inventory classification
      from checkpoint 1.
- [x] **Expected:** every failure attributes to a callsite in the
      "numerical / re-baseline expected" or "multi-CV routing"
      bucket. — Confirmed: all 3 failures attribute to
      `test_bsm2_reference.py:119`.
- [x] **Unexpected:** any failure in the "structural / no shift" or
      "bypasses `_integrate_reactions`" bucket. Stop and
      investigate. — Not triggered: zero failures in those
      buckets.
- [x] Record the failure-vs-inventory cross-reference inline in
      this checklist under §Test run below.

Sanity check: a written cross-reference exists in §Test run. Every
failure is either expected (re-baseline in checkpoint 5/6) or
flagged for investigation.

### 5. Re-baseline BSM2 golden

- [x] [`tests/standalone/test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py):
      update golden trajectory values to match the post-change
      single-Euler output. Keep `RTOL_SENTINEL = 1e-9`.
- [x] Document the drift inline at the test's location, mirroring
      the format used in previous chemistry-unification phases when
      sentinels shifted. Inline note records: (i) which species
      shifted most (S_h2 at ~6e-3 relative, in both liquid and gas
      phases), (ii) magnitude per species (1e-7 to 6e-3), (iii)
      attribution to RK4 → single Euler with pointer to
      `INTEGRATOR_REMOVAL.md`'s "Monod kinetics with substantial
      substrate depletion" framing.
- [x] Re-run `python -m pytest
      tests/standalone/test_bsm2_reference.py -q` — 6/6 pass.

Sanity check: `test_bsm2_reference.py` passes. The inline drift
note is present and attributes the shift to this phase.

### 6. Re-baseline other numerical tests

- [x] No-op for this phase. Checkpoint 4 confirmed all 3
      remaining failures (`TestBSM2Sentinels::test_final_*`) were
      in the same trajectory captured at
      `test_bsm2_reference.py:119`, all covered by the
      checkpoint-5 re-baseline. No other test in the standalone
      suite asserts on a sentinel that shifts under RK4 → single
      Euler.

Sanity check: deferred to checkpoint 8's final test run.

### 7. Documentation sweep

- [x] [`docs/architecture.md`](../../docs/architecture.md): removed
      the "Integrator (optional)" bullet from the ControlVolume
      description; removed `integrators.py` from the
      `src/reactions/` directory tree; updated the
      `_calc_ODE_with_headspace` orchestrator description to say
      "single forward Euler over `dt_h`" instead of "RK4".
- [x] [`docs/class_diagrams.md`](../../docs/class_diagrams.md):
      removed the Integrator / RK4Integrator / EulerIntegrator
      mermaid class block; removed `+integrator: Integrator` from
      the ControlVolume class block; removed the
      `ControlVolume ..> Integrator` relation arrow; rewrote the
      "two time-stepping families" paragraph to drop the
      `Integrator` half.
- [x] [`docs/solvers.md`](../../docs/solvers.md): rewrote the
      `ControlVolume.advance()` sub-step list — "(RK4 by default)"
      becomes "(single forward Euler over `dt_h`)".
- [x] [`docs/upcoming/STATE_UNIFICATION.md`](STATE_UNIFICATION.md):
      no live references (only `[[integrator-removal-phase-designed]]`
      memory link, which is intentional and stays).
- [x] [`docs/upcoming/SIMULATION_CLASS.md`](SIMULATION_CLASS.md):
      no live references (only "system-level integrator"
      placeholder for `MultiCVStepSolver`, intentional).
- [x] [`docs/design/SPECIATION_REACTIONMODEL_BOUNDARY.md`](../design/SPECIATION_REACTIONMODEL_BOUNDARY.md):
      L25 solver table previously claimed "RK4 / scipy ODE" for
      the kinetic side; updated to "single Euler / scipy ODE".
- [x] [`README.md`](../../README.md): "ODE integration for
      kinetic models (RK4, Euler, SciPy)" became "single Euler in
      the sequential body; SciPy adaptive via `ScipyODESolver`".
- [x] Stale `.pyc` removed:
      `src/reactions/__pycache__/integrators.cpython-312.pyc`.

Sanity check: `grep -rn 'RK4Integrator\|EulerIntegrator' docs/
README.md tests/ src/ models/` returns only intentional matches —
the inline re-baseline note in `test_bsm2_reference.py:189` and
this checklist + design doc + upcoming/README.md "Currently
in flight" entry (all of which will move to shipped/ at
checkpoint 8).

### 8. Final test run + ship

- [ ] Move [`INTEGRATOR_REMOVAL.md`](INTEGRATOR_REMOVAL.md) and this
      checklist to
      [`../shipped/`](../shipped/) with "Shipped"
      banners at the top of each, mirroring the chemistry-unification
      precedent.
- [ ] [`docs/upcoming/README.md`](README.md): remove
      `INTEGRATOR_REMOVAL.md` from the priority list. Move the
      "Currently in flight" entry to a shipped pointer under
      "Recently shipped" (or whatever banner section the README has
      at the time).
- [ ] Final test sweep: `python -m pytest tests/standalone -q`.
      Confirm pass count: baseline post chemistry-unification-4 is
      821; the deletion in checkpoint 3 drops the count by the
      RK4/Euler unit-test block (~2-6 tests, exact number recorded
      under §Test run). No new tests added in this phase.
- [ ] Ship via the convention in [README.md](README.md):
      ```
      git checkout main
      git merge --no-ff integrator-removal \
          -m "Merge integrator-removal: delete Integrator Protocol; single-Euler reaction sub-step"
      ```
- [ ] Tag: `git tag integrator-removal-shipped <commit-hash>`.
- [ ] Push: `git push && git push --tags`.
- [ ] Delete branch:
      `git branch -d integrator-removal` and
      `git push origin --delete integrator-removal`.

## Final test count expectation

Baseline 821 → ~821 − N, where N is the number of unit tests in the
deleted RK4/Euler block at
[`tests/standalone/test_reactions.py:540-556`](../../tests/standalone/test_reactions.py#L540).
Confirmed empirically at checkpoint 8.

## Risk notes

- **Drift attribution depends on the inventory.** If checkpoint 1
  miscategorises a test as structural and that test fails at
  checkpoint 4, the workflow says "investigate" — but the test
  may simply be numerical-with-a-loose-tolerance and a legitimate
  re-baseline target. Resolve by reading the test's actual
  assertions, not by re-categorising on the fly.
- **`_integrate_reactions` returning `y_final < 0`.** RK4 sometimes
  hides overshoot via averaging across sub-stages; single Euler is
  more prone to undershoot below zero on a fast-depleting species.
  The existing `max(0.0, ...)` clamp on the writeback survives,
  but downstream tests that assert on tight mass-balance at the
  per-step level may notice the clamp activating more often. Watch
  for this at checkpoint 4.
- **`CUFermentationSpeciation` legacy path.** Likely drifts
  numerically but its tests (if any) are out-of-scope per
  `SIMULATION_CLASS`. Any drift here is recorded inline but not
  treated as blocking.
- **Multi-CV `sys.advance_all(dt_h)` without `solvers=`.**
  [`tests/standalone/test_multi_cv.py`](../../tests/standalone/test_multi_cv.py)
  and
  [`tests/standalone/test_gas_liquid_link.py`](../../tests/standalone/test_gas_liquid_link.py)
  route to bare `cv.advance()` per CV. Most such tests are
  structural (presence-of-link tests) but a minority check
  trajectories; the inventory should disambiguate.

## Inventory

Populated at checkpoint 1 (2026-05-20). Source grep:
`\.advance\(|\.advance_all\(` across `*.py` files in `src/`,
`tests/`, `models/`. Docstring-only mentions (e.g.
`src/core/solvers.py:15`, `src/core/multi_cv.py:10`, etc.) excluded.

### Bypasses `_integrate_reactions` (no shift — explicit `solver=`)

- `models/vlmodels/fermenter/config/factory.py:589` — production
  fermenter factory; passes `solver=solver` (`EulerSnapshotSolver`
  or `ScipyODESolver` per config). All four fermenter demos route
  through this factory.
- `tests/standalone/test_cv_advance.py:709, 730` —
  `TestCVAdvanceSolverDispatch` cases, explicit
  `solver=EulerSnapshotSolver()`.
- `tests/standalone/test_multi_cv.py:613` —
  `TestAdvanceAllSolvers.test_per_cv_solver_dispatched`, explicit
  `solvers={"zone": EulerSnapshotSolver()}`.

### Production / src callsites — pass-through

- `src/core/multi_cv.py:181` — `MultiCVSystem.advance_all` forwards
  `solver=solvers.get(key)` (or `None`) to per-CV `cv.advance()`.
  Not a test assertion site; the behavioural shift propagates to
  callers that pass no `solvers=`.

### Structural / no shift expected

All callsites in this bucket either have no `reaction_model`
attached (so `_integrate_reactions` short-circuits on the
`n_species == 0` guard) or have an A→B-style equal-stoichiometry
reaction whose conservation laws hold for both RK4 and Euler, with
assertions only on direction / steady-state shape / total-moles
conservation — none of which a single-Euler reaction sub-step
would break.

- `tests/standalone/test_cv_advance.py` — 26 callsites total in
  this file. All structural:
  - L162, L167, L174 — `TestAdvanceNoReaction`: returns
    `AdvanceResult`, transfer diagnostics present, no reaction
    sources. No reaction_model.
  - L184, L196 — speciation property solver attached, no
    reaction_model. Asserts on pH presence + range.
  - L219, L226 — pre-step properties; no reaction_model. Asserts
    `pH_with == pH_without` at `abs=1e-9`.
  - L243, L252 — external sources only. No reaction_model.
    Asserts on mol-delta at `abs=1e-12`/`1e-15`.
  - L269, L276, L284, L290, L302, L312 —
    `TestAdvanceWithReaction`: A→B equal-C reaction. Direction
    asserts (`A_after < A_before`), C-conservation
    (`A1+B1 == A0+B0 abs=1e-12`), reaction_sources sign, O2
    untouched (`abs=1e-15`). Equal-stoichiometry C conservation
    holds under both RK4 and single Euler; per-species
    trajectories shift but these assertions don't see it.
  - L393, L405, L415 — `TestAdvanceMassBalance`: no-reaction
    conservation, A→B C-conservation, external accounting.
  - L514, L526, L542, L552, L566 —
    `TestKineticGasLiquidLinkInCV`: `KineticGasLiquidLink` only,
    no reaction_model. Asserts AdvanceResult presence, O2
    transfer direction, O2 total conservation (`abs=1e-12`),
    transfer diagnostics conservation, pH in range.
  - L601 — `test_kinetic_approaches_equilibrium`: 500 steps of
    gas-liquid transfer, no reaction_model. Asserts Henry
    equilibrium at `rel=0.02`.
  - L634, L643 — `TestCVBoundaries`: `LiquidFeed` boundary, no
    reaction_model. Asserts feed-applied amount + no-boundary
    conservation.

- `tests/standalone/test_membrane.py` — 3 callsites:
  - L128 — `MembraneGasBoundary` only, no reaction_model. Asserts
    O2 inventory increased.
  - L145 — well plate with O2-consuming reaction. Asserts
    steady-state shape (`max(tail)-min(tail) / mean < 0.05`) over
    2000 steps. Steady-state location is set by membrane
    permeability balancing reaction rate; the Euler/RK4 transient
    shape may shift but the steady-state assertion is robust.
  - L177 — well plate with A→CO2 reaction. Asserts
    `gas.CO2 > 0` after 500 steps. Direction assertion, robust.

- `tests/standalone/test_multi_cv.py` — 11 callsites, all
  structural:
  - L275, L285, L295 — `TestAdvanceAll`: result-type checks,
    link_records populated, conservation under no-reaction.
  - L324 — `TestMixingConvergence`: 500 steps of bidirectional
    flow, no reaction_model. Asserts concentration convergence
    at `abs=0.005`.
  - L346 — same family, mass conservation `abs=1e-10`. No
    reaction.
  - L391 — `TestReactiveAndNonreactive`: A→B equal-C reaction in
    one CV. Asserts B accumulates in non-reactive zone
    (direction) + A+B conservation `abs=1e-10`.
  - L427, L453 — `TestThreeZoneRing`: no reaction. Convergence +
    conservation asserts.
  - L485, L493 — `TestSingleCVSystem.test_single_cv_no_links_matches_direct_advance`:
    asserts that direct `cv1.advance()` and
    `sys.advance_all()` give identical results `abs=1e-14`.
    Critical observation — both routes use the same
    `_integrate_reactions` path; both shift identically, the
    identity holds.
  - L521 — `TestDiffusiveLinkInSystem`: no reaction. Convergence
    assert `abs=0.01`.

- `tests/standalone/test_gas_liquid_link.py` — 5 callsites:
  - L429 — no reaction_model, asserts mass conservation across
    gas-liquid transfer `abs=1e-12`.
  - L449 — same, equilibrium species case `abs=1e-12`.
  - L515 — `test_exponential_form_prevents_overshoot`: no
    reaction_model. Asserts monotone increase of n_liq.
  - L665 — `test_kinetic_o2_dissolution_approaches_equilibrium`:
    no reaction_model. Henry equilibrium `rel=0.01`,
    conservation `abs=1e-10`.
  - L699 — `test_kinetic_co2_stripping_approaches_equilibrium`:
    no reaction_model. Direction + conservation.

### Numerical / re-baseline expected

- `tests/standalone/test_bsm2_reference.py:119` — known sentinel
  per `INTEGRATOR_REMOVAL.md`. Bare `cv.advance()` over a full
  ADM1 reaction set with Monod kinetics and substantial substrate
  depletion per `DT_H`. Golden trajectories at the file's tail
  will shift; re-baseline at checkpoint 5.

### Production / non-test bare `cv.advance()` (drift bounded)

- `models/vlmodels/fermenter/unit.py:2217` — legacy
  `CUFermentationSpeciation` path, bare `_cv.advance(dt_h,
  chem_env)`. Out-of-scope for `SIMULATION_CLASS` migration; any
  drift here is bounded to the legacy unit and is not blocking
  per [INTEGRATOR_REMOVAL.md](INTEGRATOR_REMOVAL.md) "Sentinel
  risk → Investigate during implementation". Watch for failures
  in any test that exercises this class; the inventory above
  found no such test using `_cv.advance` directly.

### Summary

- 4 callsites bypass `_integrate_reactions` (explicit `solver=`).
- 1 production pass-through (`MultiCVSystem.advance_all`).
- 45 callsites in the structural bucket (26 in `test_cv_advance.py`,
  3 in `test_membrane.py`, 11 in `test_multi_cv.py`, 5 in
  `test_gas_liquid_link.py`).
- 1 callsite in the numerical / re-baseline bucket
  (`test_bsm2_reference.py`).
- 1 production legacy callsite (`fermenter/unit.py:2217`,
  CUFermentationSpeciation), no test exercising it found.

**Checkpoint 4 expectation:** the BSM2 reference test
(`test_bsm2_reference.py`) fails on sentinel comparisons and gets
re-baselined at checkpoint 5. Any other failure (esp. in any of the
26 `test_cv_advance.py` callsites or the 11 `test_multi_cv.py`
callsites) indicates the change has a wider blast radius than
expected and triggers the "stop and investigate" path.

## Test run

### Checkpoint 4 — initial suite run

`python -m pytest tests/standalone -q --tb=no` after checkpoints
2 + 3 landed:

- **Passed: 813**
- **Failed: 3** — all in `tests/standalone/test_bsm2_reference.py`:
  - `TestBSM2Sentinels::test_final_liquid_concentrations`
  - `TestBSM2Sentinels::test_final_gas_mol`
  - `TestBSM2Sentinels::test_final_pH`
- **Runtime: 283s** (BSM2 reference is by far the slowest test;
  the rest finish in seconds).

**Test count arithmetic:** baseline 821 (post chemistry-unification-4)
− 4 (deleted `TestIntegrators` block in checkpoint 3) − 1 (deleted
`test_snapshot_preserves_integrator` in checkpoint 2) = **816
expected**. Actual 813 passed + 3 failed = 816. ✓

**Cross-reference against §Inventory:**

- All 3 failures attribute to the single numerical/re-baseline
  callsite predicted in §Inventory:
  `tests/standalone/test_bsm2_reference.py:119`. The three
  `TestBSM2Sentinels::test_*` cases all read from the same
  trajectory captured by the L119 bare `cv.advance()` loop.
- **Zero failures in the 45-callsite structural bucket.** The
  "stop and investigate" path is not triggered — the inventory
  classification held under the actual test run.
- **Zero failures from the production legacy callsite
  (`fermenter/unit.py:2217`).** No test in the standalone suite
  exercises the CUFermentationSpeciation path through that line.
  Drift bounded as design predicted.

Re-baseline at checkpoint 5 covers all 3 sentinel failures
(same trajectory, same golden file).
