# BSM2 Reference Test — Checklist

> **Status: Shipped 2026-05-12** — all 5 checkpoints landed. 714
> standalone tests passing (708 baseline + 6 new BSM2 tests). Tag:
> `bsm2-reference-test-shipped`.

Working checklist for the `bsm2-reference-test` branch.

## Why this branch exists

The BSM2 reference implementation
([`models/vlmodels/adm1/bsm2.py`](../../models/vlmodels/adm1/bsm2.py))
has zero direct test coverage. The
[REACTION_PROTOCOL_CLEANUP](../shipped/REACTION_PROTOCOL_CLEANUP_CHECKLIST.md)
branch explicitly carved BSM2 regression coverage *out* of scope as a
separate concern (line 35–36). It is now picked up as that separate
concern, immediately before
[CHEMISTRY_UNIFICATION_PLAN.md](CHEMISTRY_UNIFICATION_PLAN.md)
begins — because `chemistry-unification-1` will perturb the
declaration surface that BSM2 sits on top of (`Reaction.kind`,
equilibrium reactions consumed by `SpeciationEngine`, charge-balance
validation), and we want any numerical drift to be **visible** rather
than silent.

## Goal

Lock the current end-of-trajectory state of a short BSM2 simulation
behind hardcoded golden values so chemistry-unification phases must
either:

- match the golden values bit-for-bit (refactor-only changes), or
- update the goldens with an explicit, reviewable diff (numerics-
  changing changes), in the same commit that introduces the change.

## Out of scope

- **`bsm2_direct.py`.** The COD-direct variant does not go through
  the speciation solver and is not on the chemistry-unification
  critical path. A separate test branch can pick it up later if a
  driver appears.
- **Literature validation.** This test detects drift from *current
  numerics*, not correctness against Rosén & Jeppsson 2006. Use a
  separate validation suite if/when literature comparison is needed.
- **Long-horizon trajectories.** A short (~1 hour sim) trajectory is
  sufficient to exercise hydrolysis + uptake + speciation + gas
  transfer paths. Days-long runs are out of scope until a driver
  appears.
- **Touching `class_diagrams.md`.** The deferred-refactor note at
  [class_diagrams.md:305-323](../class_diagrams.md#L305-L323) is
  already strong; nothing to add.

## Resolved decisions

- **Initial state: builder defaults + minimal substrate seed.** Use
  `build_bsm2_reactions()` and `build_bsm2_cv()` with default args;
  seed the liquid phase with a small documented "starter culture"
  because the builder leaves the liquid phase essentially empty. The
  seed values are arbitrary but realistic — the test validates
  *self-consistency over time*, not the seed itself.
- **Sentinels + invariants, not full snapshot.** Hardcode ~6
  sentinel species concentrations + pH + gas-phase partial pressures
  at end-of-trajectory. Add invariants (pH in AD range, non-negative
  concentrations) checked over the full trajectory. Avoids the
  ~25-number update burden of a full snapshot every time chemistry-
  unification touches numerics, while still catching drift.
- **Hardcoded goldens, not a fixture file.** Easier to review the
  diff when a future change intentionally updates them. Fixture
  files lose the "this was deliberate" signal in code review.
- **In-code docstring note on `SpeciationPropertySolver`** about its
  scheduled absorption into the unified `Reaction` framework. The
  existing note in `class_diagrams.md` is strong but invisible to
  anyone reading the class itself.

## Checkpoints

### 1. Add forward-note docstring to `SpeciationPropertySolver`

- [x] [`src/core/speciation_solver.py:37`](../../src/core/speciation_solver.py#L37) —
      append a "Forward note" paragraph to the class docstring,
      cross-referencing
      [`class_diagrams.md`](../class_diagrams.md#L305) and
      [`CHEMISTRY_UNIFICATION.md`](CHEMISTRY_UNIFICATION.md). Note
      that the class is scheduled for absorption into the unified
      `Reaction` framework once `chemistry-unification-1` lands, and
      that `PropertySolver` itself survives for genuinely property-
      flavoured concerns (viscosity, density).

### 2. Add golden-trajectory test for BSM2 reference

- [x] Create
      [`tests/standalone/test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py)
      with:
      - module-scoped fixture that builds the BSM2 CV, seeds liquid
        substrates, runs ~100 steps at `dt=0.01h`, and returns the
        final state plus history snapshots,
      - `TestBSM2Sentinels` class asserting hardcoded end-of-
        trajectory values for `S_ac`, `S_pro`, `S_h2`, `S_ch4`,
        `CO2`, `NH3`, `pH`, and gas-phase `S_ch4`/`S_h2`/`CO2` mole
        counts at `rel_tol=1e-9`,
      - `TestBSM2Invariants` class asserting pH is physically valid
        (1 < pH < 13, loose because BSM2 default `CT_cation=0`
        produces an unbuffered acidic trajectory with these seed
        VFAs — drift detection, not operational realism, is the
        purpose), all concentrations remain non-negative, and the
        speciation result is populated at the endpoint.

### 3. Capture current numerics and lock in goldens

- [x] Run the test once with placeholder goldens to print the actual
      values from the current numerics. Replace placeholders with
      the printed values. Re-run; test passes.

### 4. Final test sweep

- [x] Run the full standalone test suite. Actual: 708 → 714 (six
      new tests: three sentinels, three invariants).
- [x] `python -m pytest tests/standalone -q` is green.

### 5. Ship

- [x] Move this doc + the planning entry to `shipped/`,
      replace "In progress" with "Shipped <date>" + test count.
- [x] Update [`README.md`](README.md) priority list (remove "in
      flight" if applicable; this branch sits between
      `reaction-protocol-cleanup` and `chemistry-unification-1`).
- [x] Merge `--no-ff` into `main`, tag `bsm2-reference-test-shipped`,
      delete the branch.

## Final test count

708 standalone tests pre-branch → 714 post-branch (`TestBSM2Sentinels`
with 3 tests for liquid concentrations, gas mole counts, and pH;
`TestBSM2Invariants` with 3 tests for pH validity, non-negativity, and
endpoint speciation population).
