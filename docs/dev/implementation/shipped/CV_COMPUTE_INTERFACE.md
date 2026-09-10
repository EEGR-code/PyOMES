# CV_COMPUTE_INTERFACE — Phase A

## Status

**Shipped 2026-06-11.** Tag: `cv-compute-interface-shipped`. 1226→1256 tests.
See [`../shipped/SOLVER_ARCHITECTURE.md`](../shipped/SOLVER_ARCHITECTURE.md)
for the full two-axis context.

## Goal

Add a clean computational interface to `ControlVolume` and `SpeciationEngine`
that system-level and DAE solvers can program against. Refactor `cv.advance()`
to call `compute_rhs()` internally so its behaviour is unchanged while the
building blocks become independently reachable. No new user-visible behaviour.
All existing tests pass without modification.

## New `ControlVolume` methods

### `compute_rhs(t_h) → Dict[str, Dict[str, float]]`

ODE interface. Runs speciation by index reduction (engine writes derived
species to `phase.n_mol`), then evaluates kinetic reaction rates. Returns
`{phase_key: {species_id: rate_mol_per_h}}`. Does not mutate CV state beyond
the speciation write-back (which is idempotent for a given totals state).

`cv.advance()` is refactored to call this for its Euler step.

### `compute_differential_rhs(t_h, y_state, z_state) → Dict[str, Dict[str, float]]`

```python
def compute_differential_rhs(
    self,
    t_h: float,
    y_state: Dict[str, Dict[str, float]],  # {phase_key: {species: n_mol}} — differential totals
    z_state: Dict[str, Dict[str, float]],  # {phase_key: {species: n_mol}} — algebraic derived
) -> Dict[str, Dict[str, float]]:
```

DAE interface. Evaluates kinetic reaction rates at explicit (`y_state`,
`z_state`) without calling the speciation engine. Required by `DAEStepSolver`
(Phase F) and `DAESystemSolver` (Phase G). Ships as a real implementation in
Phase A: merge y and z into a working concentration dict, build
`ReactionEnvironment`, call `reaction_system.compute_rates`. Does not touch
`phase.n_mol`.

### `compute_algebraic_residual(t_h, y_state, z_state) → Dict[str, Dict[str, float]]`

DAE interface. Returns the speciation charge-balance residual at explicit
(`y_state`, `z_state`). Ships as a `NotImplementedError` stub. Phase F
activates the real implementation when `DAEStepSolver` is written. Adding the
stub now means Phase F does not require a breaking change to the CV interface.

### `compute_jacobian(t_h) → Optional[np.ndarray]`

Returns the Jacobian of `compute_rhs` with respect to the flattened per-CV
species vector, shape `(S_total, S_total)`, or `None` if not implemented.
Ships returning `None`. Individual model implementations may override for
performance. Used by `MonolithicODESolver` (Phase E) to avoid finite
differences when available.

### `snapshot_state() → Dict[str, Dict[str, float]]`

Returns `{phase_key: dict(phase.n_mol)}` — a deep copy of the current mole
inventory across all phases. Used by multi-sub-step solvers for trial-state
evaluation:

```python
snap = cv.snapshot_state()
cv._set_trial_state(trial_vec)
rhs = cv.compute_rhs(t_h)
cv.restore_state(snap)
```

### `restore_state(snapshot) → None`

Writes back a snapshot taken by `snapshot_state()`. Overwrites `phase.n_mol`
entries directly. Pair with `snapshot_state()` for reversible trial-state
evaluation.

## New `SpeciationEngine` method

### `algebraic_species() → frozenset[str]`

Returns the set of species IDs the engine writes to `phase.n_mol` via
`_refresh_derived`. This is the algebraic state `z` for any DAE formulation;
the complement within `phase.n_mol` is the differential state `y`.

Implementation: inspect the output species of each entry in `_equilibrium_set`.
Returns `frozenset()` when no `EquilibriumSet` is bound (legacy kwargs path;
all species are differential in that case).

## `advance()` refactor

The sequential body of `cv.advance()` currently computes the reaction Euler
step inline via `_integrate_reactions()`. After Phase A, step 3 delegates to
`compute_rhs()`:

```python
# Step 3 (reactions) — after Phase A:
if self.reaction_system is not None:
    rhs = self.compute_rhs(t_h)
    for pk, sp_rates in rhs.items():
        for sp, rate in sp_rates.items():
            phase = self.phases[pk]
            phase.n_mol[sp] = max(0.0, phase.n_mol.get(sp, 0.0) + rate * dt_h)
```

`_integrate_reactions` may be retained as a thin wrapper or removed; its
logic moves into `compute_rhs`. User-visible behaviour of `cv.advance()` is
unchanged.

## Implementation checklist

- [ ] `cv.compute_rhs(t_h)` — real implementation
- [ ] `cv.compute_differential_rhs(t_h, y_state, z_state)` — real implementation
- [ ] `cv.compute_algebraic_residual(t_h, y_state, z_state)` — stub (`NotImplementedError`)
- [ ] `cv.compute_jacobian(t_h)` — stub (returns `None`)
- [ ] `cv.snapshot_state()` — real implementation
- [ ] `cv.restore_state(snapshot)` — real implementation
- [ ] `engine.algebraic_species()` — real implementation on `SpeciationEngine`
- [ ] `cv.advance()` refactored to call `compute_rhs()` for the Euler step body
- [ ] `_integrate_reactions` refactored or removed
- [ ] All existing tests pass without modification
- [ ] New tests: `compute_rhs` values match deltas produced by `cv.advance()`
      at the same state
- [ ] New tests: `snapshot_state` / `restore_state` round-trip fidelity
- [ ] New tests: `algebraic_species` returns correct set for BSM2/ADM1 engines
- [ ] New tests: `compute_differential_rhs` gives same kinetic rates as
      `compute_rhs` when called with consistent y + z state

## Resolved positions

### `advance()` calls `compute_reaction_rates`, not `compute_rhs`

The implementation checklist and advance() refactor snippet in this document
stated that `cv.advance()` would call `compute_rhs()` for its Euler step.
**This was not implemented that way**, and for good reason.

`compute_rhs` re-runs the speciation engine before evaluating rates.  In
`advance()`, feeds are applied in step 2 *after* speciation has already run
in step 1.  Calling `compute_rhs` in step 3 would re-solve speciation on the
*post-feed* totals, shifting H⁺ and pH within the same macro step.  This
breaks the operator-splitting contract that VLsim models depend on:

> **Speciation is pinned at the pre-step totals for the duration of one
> `advance()` call.**  External sources injected in step 2 do not affect
> the same-step pH.

`test_advance_speciation_reflects_pre_step_state` in `test_cv_advance.py`
encodes this contract explicitly and caught the regression.

**What was implemented instead:** `advance()` step 3 calls
`compute_reaction_rates(t_h)` directly (no speciation re-run).  Speciation is
run once, in step 1, and its output is pinned for the entire step.
`compute_rhs` is the right interface for *external* callers (system-level
solvers in Phase C+) that need a self-contained RHS evaluator and are
responsible for managing their own speciation timing.

**Does this need a follow-up refinement?** No.  The boundary is clean:

- `advance()` (sequential, operator-split) → uses `compute_reaction_rates`
- Phase C+ system solvers (own the full step) → use `compute_rhs`

The only change needed to this document is this note.  The checklist item
"`cv.advance()` refactored to call `compute_rhs()` for the Euler step body"
should be read as: *the Euler step body is now expressed in terms of
`compute_reaction_rates` + inline dt multiply, with `compute_rhs` providing
the independently reachable equivalent for external callers.*  The building
blocks are independently reachable; `advance()` just uses the narrower one
that preserves the operator-splitting invariant.

## Cross-references

- [`SOLVER_ARCHITECTURE.md`](SOLVER_ARCHITECTURE.md) — two-axis context
- [`SYSTEM_SOLVER_PROTOCOL.md`](SYSTEM_SOLVER_PROTOCOL.md) — Phase C; first consumer
- [`../../src/core/control_volume.py`](../../src/core/control_volume.py)
- [`../../src/speciation/engine.py`](../../src/speciation/engine.py)
