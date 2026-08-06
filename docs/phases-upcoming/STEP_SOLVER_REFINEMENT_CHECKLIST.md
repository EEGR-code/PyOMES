# Step-Solver Interface Refinement — Checklist

> **Status: Shipped 2026-07-08.** All 10 checkpoints (0 through 9) landed
> on branch `step-solver-interface-refinement`, tag
> `step-solver-interface-refinement-shipped`. Working checklist for
> that branch (Stage 1 of
> [STEP_SOLVER_REFINEMENT_PLAN.md](STEP_SOLVER_REFINEMENT_PLAN.md)).
> Conceptual framing is in
> [`../phases-shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md`](../phases-shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md)
> (items 1, 2, 3, 4, 6, 7, 8, 9 — item 5 explicitly deferred, item 10
> shipped as checkpoint 0 above). This file is the implementation log:
> file-level edits, ordered checkpoints, sanity checks. Modelled on
> [`../phases-shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md`](../phases-shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md).
> 1950 → 2003 tests (27 skipped throughout), 0 failed.

## Goal

Close two live correctness gaps in the shipped two-axis solver
architecture (`cv.advance()`'s missing ownership guard;
`MonolithicODESolver` silently discarding per-CV `solver=`), and ship
five independent, low-risk cleanups the same design review surfaced
(reifying the sequential default, unifying state-vector packing,
generalizing off the gas/liquid assumption, a swappable clamping
module with its diagnostics, and a "write your own solver" demo).
None of this changes numerics for any existing caller — every item is
additive or a structural rename with byte-identical default behavior.

## Out of scope

- **Item 5 (composition/integrator/clamp_fn as three independently
  swappable dimensions).** Explicitly deferred by the source design
  note — a bigger lift than anything else in this list, to be built
  only when a concrete model needs a (composition, method) pair the
  two shipped `StepSolver`s can't express.
- **Adaptive macro-stepping** (`AdvanceResult` "actual dt consumed"
  field) and **decoupling `StepSolver` from the concrete
  `ControlVolume` class** — both listed as open questions with "no
  trigger yet" in the source note. Not touched here.
- **The interleaving/multi-CV solver tier** (SIA, reactive D_eff) —
  Stages 2–4 of the parent plan, not this branch.
- **`MASS_EXCHANGE_ARCHITECTURE.md` §11 priority-2 item** (engine-scope
  construction-time validation via `algebraic_species()`) — unrelated
  track, not part of this bundle.

## Resolved decisions

Pinned during the 2026-07-06 scoping review (see
`STEP_SOLVER_REFINEMENT_PLAN.md`'s *Verification notes* for the
code-level evidence behind each):

- **`OrchestrationWarning` lives in `src/core/control_volume.py`**,
  next to `ControlVolume`, mirroring the existing precedent of
  `AccuracyWarning` in `monitoring/accuracy.py` and
  `ConservationWarning` in `monitoring/conservation.py` (warning
  categories are defined next to the code that emits them, not in a
  shared warnings module). **Correction found at checkpoint 1:**
  top-level `src/__init__.py` (the `VLsim` package) does not export
  anything from `src/core/` at all — `ControlVolume`, `Simulation`,
  `StepSolver`, `AccuracyWarning`'s core-adjacent siblings, none of it.
  `AccuracyWarning`/`ConservationWarning` are exported at top level
  because they come from `src/monitoring/`, a package top-level
  already re-exports from; `core/` primitives are accessed via
  `VLsim.core.X` only. So `OrchestrationWarning` is exported from
  `src/core/__init__.py` only, not top-level `src/__init__.py`.
- **`_advance_unchecked()` is the trusted internal entry point** —
  every `SystemSolver.advance_system()` that currently calls the
  public `cv.advance(dt_h, t_h, solver=solver)` must switch to it, or
  the ownership guard will fire on every orchestrated step. Four call
  sites confirmed (checkpoint 2 lists them by file:line).
- **`SequentialAdvanceSolver` goes in `src/core/solvers.py`**,
  alongside `SimultaneousEulerSolver`/`SimultaneousAdaptiveSolver` —
  same file, same `StepSolver`-conforming shape, consistent with
  where the other two rungs of Axis 1 already live.
- **Clamp-module scope for this phase**: wire `clamp_fn` into
  `SequentialAdvanceSolver` and `SimultaneousEulerSolver` only, per
  the source note. `SimultaneousAdaptiveSolver`'s unpack-time floor
  stays conceptually separate (different failure mode — adaptive
  overshoot, not discrete-dt overdraw) but should be named/swappable
  too rather than a bare `max(0.0, …)`, per the note's own guidance.

## Checkpoints

Ordered so state-vector/solver-object changes land before the
clamp/diagnostics work that builds on them, and the demo notebook
comes last so it can show the final API shape rather than an
intermediate one.

**Commit discipline:** commit after each checkpoint leaves the suite
green — per-checkpoint commits are the repo's established convention
(see `CHEMISTRY_UNIFICATION_4_CHECKLIST.md`'s note on why this is
load-bearing, not optional).

### 0. Close out Stage 0 and open the branch — done 2026-07-06

- [x] Confirm `rename-simultaneous-step-solvers` is merged to `main`
      and tagged `rename-simultaneous-step-solvers-shipped` (see
      `STEP_SOLVER_REFINEMENT_PLAN.md` → *How to start Stage 1*).
      Correction versus the plan: the branch had never actually been
      committed (branch and `main` were at the same commit; the rename
      was sitting entirely as uncommitted working-tree changes). Split
      into two commits before merging: the rename itself (`92a1ee7`)
      and an unrelated stale-doc-status fix in
      `MASS_EXCHANGE_ARCHITECTURE.md` (`5abed00`) that was sitting in
      the same uncommitted diff. Merged `--no-ff` (`d6c680f`), tagged,
      pushed, local+remote branch deleted (remote never existed).
- [x] `git checkout -b step-solver-interface-refinement` off the
      updated `main`.

Sanity check: `git log --oneline -3` on the new branch shows the
rename merge commit as an ancestor; `pytest tests/standalone -q`
passes at 1950/0 before any new edits.

### 1. Reify `SequentialAdvanceSolver`

- [x] [`src/core/solvers.py`](../../src/core/solvers.py) — new class
      `SequentialAdvanceSolver` implementing `StepSolver`. Moved the
      inline body of `ControlVolume.advance()` (former lines 601–670)
      into `solve_step(self, cv, dt_h, t_h=0.0, external_source_terms=None)`
      verbatim, renaming `self.<attr>` → `cv.<attr>` throughout.
- [x] [`src/core/control_volume.py`](../../src/core/control_volume.py)
      `advance()` — replaced the dispatch/inline-body fallback with
      `solver = solver if solver is not None else SequentialAdvanceSolver()`
      followed by one unconditional `solver.solve_step(...)` call.
      Docstring updated to stop describing a separate inline path.
- [x] [`src/core/__init__.py`](../../src/core/__init__.py) — export
      `SequentialAdvanceSolver` alongside `SimultaneousEulerSolver` /
      `SimultaneousAdaptiveSolver`. **Correction versus the plan:**
      top-level `src/__init__.py` (the `VLsim` package) does not
      export `ControlVolume`, `Simulation`, or any solver at all —
      those all live under `VLsim.core` only. The *Resolved decisions*
      instruction to also touch `src/__init__.py` was wrong; not done.
      Also fixed a pre-existing gap found while touching this file:
      `SimultaneousAdaptiveSolver` was imported in `core/__init__.py`
      but missing from `__all__` — added.

Sanity check: added
`TestCVAdvanceSolverDispatch::test_default_is_sugar_for_sequential_advance_solver`
in `test_cv_advance.py`, confirming `cv.advance(dt_h, t_h)` and
`cv.advance(dt_h, t_h, solver=SequentialAdvanceSolver())` produce
identical `n_mol` state and `reaction_sources`. Full suite:
1950 → 1951 passed, 27 skipped, 0 failed.

### 2. Ownership guard on `cv.advance()` + fix internal call sites — done

- [x] [`src/core/control_volume.py`](../../src/core/control_volume.py)
      — added `class OrchestrationWarning(UserWarning)` (module level,
      next to `ControlVolume`). Split `advance()` into a thin public
      wrapper (checks `self._context is not None`, warns with
      `stacklevel=2` if so, then calls `_advance_unchecked`) and
      `_advance_unchecked(self, dt_h, t_h=0.0, *,
      external_source_terms=None, solver=None)` holding the dispatch
      logic from checkpoint 1.
- [x] Updated all four internal callers to bypass the guard:
      - [`src/core/simulation.py:859`](../../src/core/simulation.py#L859)
        (`Simulation._step_default`)
      - [`src/core/system_solver.py:245`](../../src/core/system_solver.py#L245)
        (`StrangSplittingSystemSolver.advance_system`)
      - [`src/core/system_solver.py:350`](../../src/core/system_solver.py#L350)
        (`MultirateSystemSolver.advance_system`)
      - [`src/core/system_solver.py:563`](../../src/core/system_solver.py#L563)
        (`ImplicitTransportSystemSolver.advance_system`)

      All four now call `cv._advance_unchecked(dt_h, t_h, solver=solver)`.
      `MonolithicODESolver` needed no change — confirmed it never calls
      `cv.advance()` at all.
- [x] [`src/core/__init__.py`](../../src/core/__init__.py) — exported
      `OrchestrationWarning`. **Correction versus the plan** (same as
      checkpoint 1): not exported from top-level `src/__init__.py` —
      that package doesn't re-export anything from `core/`.

Sanity check: added `TestOrchestrationWarningGuard` (4 tests) to
`test_system_solver.py` — bare unowned CV emits nothing;
`Simulation`-owned CV calling `cv.advance()` directly emits exactly
one `OrchestrationWarning`; `cv._advance_unchecked()` on the same
owned CV emits nothing; a real `sim.run()` under all five shipped
`SystemSolver`s (`None`/Explicit/Strang/Multirate/ImplicitTransport/
Monolithic) on a 2-CV linked simulation emits **zero**
`OrchestrationWarning`s. Full suite: 1951 → 1955 passed, 27 skipped,
0 failed.

### 3. `MonolithicODESolver` validates per-CV `solver=` — done

- [x] [`src/core/system_solver.py`](../../src/core/system_solver.py)
      `MonolithicODESolver.advance_system()` — raises `ValueError` at
      entry if `sim.solver is not None`. **Decision made:** checked at
      `advance_system()` entry, not at construction/setter time —
      `Simulation.solver`/`Simulation.system_solver` are independent
      settable properties with no cross-validation between them
      (`simulation.py` lines 270–292), and either can be (re)assigned
      in any order before `sim.run()`. Entry-time checking is correct
      regardless of assignment order; construction-time checking would
      only catch the case where both happen to be passed to `__init__`
      together and would miss later reassignment via the property
      setters.

Sanity check: added `TestMonolithicRejectsPerCVSolver` (4 tests) to
`test_monolithic_ode_solver.py` — raises for a single `StepSolver`
instance, raises for a per-CV dict, does not raise when `sim.solver`
is `None` (the common case), and raises via a real `sim.run()` call
regardless of whether `system_solver` or `solver` was assigned first.
Full suite: 1955 → 1959 passed, 27 skipped, 0 failed.

### 4. Generalize off the gas/liquid assumption — `SimultaneousEulerSolver` done, `SimultaneousAdaptiveSolver` paused (see below)

- [x] [`src/core/solvers.py`](../../src/core/solvers.py)
      `SimultaneousEulerSolver` — replaced the `raise ValueError` guard
      requiring exactly `{"gas", "liquid"}` with a narrower check that
      only `"liquid"` must be present (speciation/reactions are
      liquid-scoped — matches `SequentialAdvanceSolver`'s existing
      requirement). Snapshot construction generalized to
      `{k: p.snapshot() for k, p in cv.phases.items()}`; per-phase
      delta dict `deltas: Dict[str, Dict[str, float]]` replaces the
      hardcoded `gas_deltas`/`liq_deltas` pair; apply/clamp loop over
      `cv.phases` generically.
- [x] Replaced the `next(iface for iface in cv.internal_interfaces if
      isinstance(iface, KineticGasLiquidLink))` search with a loop over
      **every** `cv.internal_interfaces` entry, routing each one's flux
      via its own `phase_a_key`/`phase_b_key` rather than assuming
      `"gas"`/`"liquid"`. **Bonus fix found while doing this:** the old
      code only ever processed the *first* matching interface — a CV
      with more than one registered `PhaseInterface` silently dropped
      every interface after the first. Now all of them contribute to
      `deltas`; `AdvanceResult.transfer_record` (a single optional
      slot, not a list) still only captures one for backward
      compatibility — documented in a code comment,
      `step_internal_transfer`'s `TransferDiagnostics` remains the
      complete multi-interface audit trail.
- [x] Verified the boundaries-handling loop needed no *further* change
      beyond generalizing which `deltas[phase_key]` sub-dict each
      boundary's flux is routed into (was hardcoded to two named
      dicts).

Sanity check: added `test_snapshot_solver_generalizes_to_liquid_only_cv`
(liquid-only CV, no gas phase, now succeeds under
`SimultaneousEulerSolver` — previously raised) and updated
`test_snapshot_solver_fails_fast_on_missing_phase` →
`test_snapshot_solver_fails_fast_on_missing_liquid_phase` (a gas-only
CV still correctly fails fast, since `"liquid"` is still required).
All existing gas-liquid trajectories unchanged (snapshot construction
and delta routing are behaviorally identical for the `{"gas",
"liquid"}` case, just generalized). Full suite: 1959 → 1960 passed,
27 skipped, 0 failed.

**`SimultaneousAdaptiveSolver` — found to be a substantially bigger,
riskier task than "mechanical," paused for a decision rather than
guessed at:**

Reading the full `solve_step()` (not just skimming) surfaced that,
unlike `SimultaneousEulerSolver`, this solver's entire ODE hot loop is
hand-unrolled around exactly two named phases for performance:

- `_StateVector` (the class this checkpoint's own item 4 — state-vector
  unification — is meant to eventually replace) is hardcoded to
  `cv.phases["gas"]`/`cv.phases["liquid"]` throughout: species
  ordering, `pack`/`unpack_to_dicts`/`write_to_cv`. Generalizing
  `SimultaneousAdaptiveSolver` off gas/liquid properly means
  generalizing `_StateVector` first — which is checkpoint 7's job, not
  this one. Doing it twice (once ad hoc here, once properly in
  checkpoint 7's unified `state_vector.py`) is wasted, riskier work.
- The `f(t, y)` RHS closure captures `gas_phase`/`liq_phase`/`gas_V_L`/
  `liq_V_L`/etc. as closure locals for hot-loop performance, and calls
  `iface.compute_flux(..., instantaneous=True)` — an `instantaneous=`
  kwarg that is **not** part of the `PhaseInterface` protocol
  (confirmed by reading `interfaces.py`); it's specific to
  `KineticGasLiquidLink.compute_flux`'s actual signature. Since
  `KineticGasLiquidLink` is currently the *only* concrete
  `PhaseInterface` implementation in the codebase (confirmed —
  `KineticTransferModel`/`EquilibriumTransferModel` are per-species
  config objects consumed by `_build_transfer_link`, which still
  constructs a `KineticGasLiquidLink` under the hood; there is no
  second interface class), generalizing the *iteration* here without
  also generalizing what `instantaneous=` means for an arbitrary
  interface would be a change that compiles but doesn't actually widen
  what's usable.
- The Jacobian machinery (`use_engine_jacobian`, `_alg_override`,
  `GrayBoxEngineProtocol` detection) is also two-phase-shaped.

**Decided (user, 2026-07-06): option (b) — sequence after checkpoint 7.**
Checkpoint order for the rest of this phase is now: 5 (clamp module) →
6 (diagnostics) → 7 (state-vector unification) → **7b**
(`SimultaneousAdaptiveSolver` gas/liquid generalization, built on the
new unified `state_vector.py` instead of hand-rolling a throwaway
generalization of `_StateVector` now) → 8 (demo) → 9 (ship). See
checkpoint 7b below, inserted after checkpoint 7 in this file.

### 5. Shared clamp module

- [ ] **Before writing code:** confirm scope against
      [`src/core/phases.py`](../../src/core/phases.py)'s `apply_flux`
      (lines ~262/419/517 on `GasPhase`/`LiquidPhase`/`SolidPhase`) —
      each already has its own unconditional `max(0.0, current + rate
      * dt)` floor baked into the Phase method itself, used by *every*
      caller (feeds, boundaries, links), not just
      `SequentialAdvanceSolver`.

      **Decided:** `apply_flux` gained a `clamp: bool = True` keyword
      (default preserves every existing caller's behavior exactly).
      `clamp_fn` applies **per sub-step** for `SequentialAdvanceSolver`
      (feeds, each boundary, reactions independently — this solver is
      sequential/Lie-Trotter, there's no single combined dict the way
      `SimultaneousEulerSolver` has one), then that sub-step's
      `apply_flux` call passes `clamp=False` so the built-in floor
      doesn't silently re-clamp (or un-defeat a deliberate
      `clamp_fn=None`) on top. This *is* "bypassing `apply_flux`'s
      built-in floor for the feed/boundary sub-steps" — proven
      numerically equivalent to the old behavior first: for a single
      species clamped independently, `proportional_clamp`'s
      rate-scaling and `apply_flux`'s old result-flooring produce the
      *identical* final `n_mol` (both reduce to "land exactly on
      zero"), so the default configuration is byte-identical to
      pre-checkpoint-5, confirmed by the full suite still passing
      unchanged (see below).
- [x] New module [`src/core/clamping.py`](../../src/core/clamping.py) —
      `proportional_clamp` (moved verbatim from `solvers.py`'s
      `_clamp_deltas`, now deleted), `floor_clamp` (new, with an `eps`
      parameter for reference models that floor at a small epsilon
      rather than exact zero). Both operate on an arbitrary subset of
      species — confirmed by `test_operates_on_arbitrary_species_subset`.
- [x] `SimultaneousEulerSolver.__init__` — added
      `clamp_fn=proportional_clamp` constructor parameter; step 4 now
      calls `self.clamp_fn(...)` (skipped entirely when
      `clamp_fn=None`) and step 5's `apply_flux` calls pass
      `clamp=False`.
- [x] `SequentialAdvanceSolver.__init__` — added
      `clamp_fn=proportional_clamp` constructor parameter; `clamp_fn=None`
      disables clamping (each sub-step passes its raw, possibly
      negative-inducing deltas straight through with `clamp=False`).
- [ ] `SimultaneousAdaptiveSolver`'s hardcoded unpack-time floor —
      **not done, folded into checkpoint 7b instead.** On inspection
      this floor operates on the raw ODE state array/dict directly
      (`max(0.0, float(y[i]))`), not on a `(deltas, current_mol, dt_h)`
      triple — a genuinely different shape than `clamp_fn`, not a
      drop-in `floor_clamp` call. Naming it properly means touching
      the same hot-loop code checkpoint 7b is already rewriting for
      the gas/liquid generalization; doing it twice (once here, once
      there) would be wasted work on fragile, performance-sensitive
      code. Added as an explicit bullet to checkpoint 7b.

Sanity check: added `tests/standalone/test_clamping.py` (13 tests) —
`proportional_clamp`/`floor_clamp` as pure functions (including the
arbitrary-subset requirement); `clamp_fn=None` on both
`SequentialAdvanceSolver` and `SimultaneousEulerSolver` produces
negative `n_mol` on a deliberately overdrawing rate law (proves
clamping was doing real work before); default `clamp_fn` behaves
identically to pre-checkpoint-5 on the same scenario; `floor_clamp`
swappable in as an alternative; the design note's own
`bsm2_style_clamp` composite example built out of `proportional_clamp`
as an ingredient, confirmed to clamp `FLOOR_AT_EPSILON` species to
`eps` while leaving `NEVER_CLAMP` species free to go negative and
`Normal` species landing at exactly zero via the shared default. Full
suite: 1960 → 1973 passed, 27 skipped, 0 failed — confirms the
default-`clamp_fn` behavioral-equivalence claim above across every
existing gas-liquid/BSM2/ADM1 trajectory, not just the new scenario.

### 6. Negative-mole check + clamp-invoked diagnostic — done

- [x] `AccuracyMonitor.check_negative_mole(n_mol)` — scans a post-step
      `n_mol` dict for values below `-negative_mole_tolerance` (new
      `WarningConfig` field, default `1e-9` mol); emits `AccuracyWarning`
      category `negative_mole` naming the worst offender and how many
      species were affected. Wired into both `SequentialAdvanceSolver`
      and `SimultaneousEulerSolver`, only when `self.clamp_fn is None`
      (guarded via a `_monitor_negative_mole(cv, n_mol)` module helper,
      mirroring the existing `_monitor_pH_post_step` pattern).
- [x] `AccuracyMonitor.check_clamp_invoked(deltas_before, deltas_after, dt_h)` —
      compares a delta dict before/after `clamp_fn`; emits
      `AccuracyWarning` category `clamp_invoked` naming every changed
      species and its before→after rate (beyond
      `clamp_invoked_tolerance`, new `WarningConfig` field, default
      `1e-9` mol/h), in the source note's message format. Wired into
      every `clamp_fn` call site in both solvers (three per-substep
      sites in `SequentialAdvanceSolver`; the one combined-dict site in
      `SimultaneousEulerSolver`) via a `_monitor_clamp_invoked(cv,
      before, after, dt_h)` module helper.

Sanity check: added 8 unit tests to `test_accuracy_monitor.py`
(`TestCheckNegativeMole`, `TestCheckClampInvoked` — threshold
crossing, no-op below threshold, worst-offender naming, unchanged-rate
no-op) and 4 integration tests to `test_clamping.py` — a deliberately
aggressive `dt_h`/kinetics combination with default clamping emits
`clamp_invoked` (mild kinetics emits nothing); the same scenario with
`clamp_fn=None` instead emits `negative_mole` and no `clamp_invoked`
(confirmed mutually exclusive — one only fires when clamping ran, the
other only when it didn't, exactly per the original sanity-check
design). Full suite: 1973 → 1985 passed, 27 skipped, 0 failed.

### 7. Unify state-vector packing — `system_solver.py` pair done, `_StateVector` folded into checkpoint 7b

- [x] New module [`src/core/state_vector.py`](../../src/core/state_vector.py)
      — `StateVector` class with `pack()`/`unpack(y, floor=False)`/
      `index_map`. Parameters: `cvs` (dict, scope), `exclude_species`
      (frozenset, default empty), `ctrl_list` (optional controller
      differential-state extension).
- [x] [`src/core/system_solver.py`](../../src/core/system_solver.py) —
      `_pack_state`/`_unpack_state`/`_state_index_map`/
      `_pack_extended_state`/`_unpack_extended_state` are now thin
      wrapper functions delegating to `StateVector` internally, **kept
      at their original signatures** so every existing caller
      (`ImplicitTransportSystemSolver`'s matrix assembly,
      `MonolithicODESolver`, `_build_rhs`) needed zero changes.
      **Scope correction found while implementing:** the checklist's
      "replace all three call sites" undersold what "call site" meant
      for `_StateVector` — its `gas_species`/`liq_species`/`n_gas`/
      `n_liq` attributes are read directly (not through a clean
      function call) at ~15 separate points throughout
      `SimultaneousAdaptiveSolver`'s ODE hot loop and Jacobian
      machinery. Migrating it here would mean rewriting all of those
      call sites anyway — indistinguishable from checkpoint 7b's own
      scope. Left `_StateVector` untouched in this checkpoint; folded
      its migration into checkpoint 7b (see that checkpoint's updated
      task list) so the hot loop is rewritten once, not twice — exactly
      the reason 7b was deferred in the first place.
- [x] `_unpack_extended_state`'s `n_cv` parameter is now vestigial
      (`StateVector` recomputes its own offset from `sim.cvs`) but kept
      in the signature for compatibility with existing callers;
      documented in its docstring. A regression test
      (`test_unpack_extended_state_ignores_stale_n_cv_arg`) confirms a
      wrong/stale `n_cv` can't corrupt the result.

Sanity check: added `tests/standalone/test_state_vector.py` (12 tests)
— `StateVector` pack/unpack roundtrip, CV→phase→species alphabetical
ordering, `exclude_species`, controller-state extension, and a
regression class confirming all five `system_solver.py` thin wrappers
still produce identical output to their pre-refactor bodies. Every
`ImplicitTransportSystemSolver`/`MonolithicODESolver` trajectory test
passes unchanged (pure refactor — call sites untouched). Full suite:
1985 → 1997 passed, 27 skipped, 0 failed.

Sanity check: `SimultaneousAdaptiveSolver`, every `SystemSolver`, and
`MonolithicODESolver`'s controller co-integration all produce
identical trajectories to pre-change on their respective existing
tests (this is a pure refactor — any numeric diff is a bug in the
unification, not an intended change).

### 7b. Generalize `SimultaneousAdaptiveSolver` off the gas/liquid assumption — done

Deferred from checkpoint 4 (2026-07-06 decision) to run after
checkpoint 7 so this uses the unified `state_vector.py` instead of
hand-rolling a second, throwaway generalization of `_StateVector`.

- [x] [`src/core/solvers.py`](../../src/core/solvers.py)
      `SimultaneousAdaptiveSolver.solve_step()` — replaced the
      `{"gas", "liquid"}` `ValueError` guard with the same `"liquid"`-only
      requirement `SimultaneousEulerSolver` has.
- [x] Replaced `_StateVector(cv)` with
      `StateVector({"_cv": cv}, exclude_species=frozenset({"H+"}))`
      (single-CV scope, `"_cv"` a placeholder key never surfaced to the
      caller) plus two new `StateVector` methods needed for this hot
      loop's per-phase index arithmetic: `phase_species(cv_key)` →
      `{phase_key: [species, ...]}` and `phase_offset(cv_key)` →
      `{phase_key: starting_index}`, both added to `state_vector.py` in
      this checkpoint (not checkpoint 7 — the need wasn't apparent
      until this solver's actual access pattern was rewritten).
      **`_StateVector` itself is deleted** (no backwards-compat shim,
      per house convention) — it had exactly one caller, now migrated.
- [x] Rewrote the `f(t, y)` RHS closure's phase-specific locals
      (`gas_phase`/`liq_phase`/`gas_V_L`/`liq_V_L`/`ext_gas_rates`/
      `ext_liq_rates`/`boundary_rates_gas`/`boundary_rates_liq`) to loop
      over `phase_layout`/`phase_offset` generically — same
      per-phase-dict-then-one-indexed-pass pattern already used in
      `SimultaneousEulerSolver`'s generalized `solve_step` (checkpoint
      4), for consistency between the two solvers.
- [x] Replaced the single `link = next(iface for iface in
      cv.internal_interfaces if isinstance(iface, KineticGasLiquidLink))`
      search with a loop over every `cv.internal_interfaces` entry,
      accumulating into a per-phase `transfer_rates` dict (same
      bonus-fix pattern as checkpoint 4 — every registered interface
      now contributes, not just the first).
- [x] Resolved the `instantaneous=True` kwarg question as **option
      (a)**: added `instantaneous: bool = False` to the `PhaseInterface`
      protocol declaration in `interfaces.py`, mirroring
      `ExternalBoundary.compute_flux`'s already-existing identically-named
      parameter (confirmed present on all six `ExternalBoundary`
      implementations). Zero behavior change —
      `KineticGasLiquidLink.compute_flux` already accepted this kwarg;
      the protocol declaration was just formalizing an existing
      undocumented contract.
- [x] Verified the Jacobian machinery: `sv.n_gas`/`sv.n_liq`/
      `sv.liq_species` became `liq_offset`/`n_liq`/`liq_species` derived
      once from `phase_offset.get("liquid", 0)` /
      `phase_layout.get("liquid", [])` before `jac_callable` is defined
      (not phase-generic *within* the Jacobian itself — the analytical
      dz_dy path is inherently liquid-component-scoped by the
      speciation engine's own contract — but no longer hardcoded to
      assume gas occupies indices `[0, n_gas)`).
- [x] **Carried over from checkpoint 5:** added
      `floor_nonnegative(y) -> y` to `clamping.py` (vectorized
      `np.maximum(y, 0.0)`), replacing three remaining bare inline
      `max(0.0, float(y[...]))` floors — one in `f(t,y)`'s state
      unpack, one in the Jacobian's component-totals computation, one
      in `StateVector.unpack`'s `floor=True` branch (checkpoint 7's
      code, updated here since this checkpoint is what actually needed
      it). Kept conceptually separate from `clamp_fn` per the source
      note (documented in `floor_nonnegative`'s own docstring).

Sanity check: existing gas-liquid CV trajectories under
`SimultaneousAdaptiveSolver` (including the Jacobian-enabled path,
`use_engine_jacobian=True`, `test_simultaneous_adaptive_solver_jac.py`'s
full 24-test suite) are unchanged. Added
`test_adaptive_solver_fails_fast_on_missing_liquid_phase` and a
3-way-parametrized (`DOP853`/`BDF`/`Radau`)
`test_adaptive_solver_generalizes_to_liquid_only_cv` to
`test_cv_advance.py` — a liquid-only CV (no gas phase) now runs under
every method, previously raised `ValueError` for all of them. Plus 2
new `floor_nonnegative` unit tests in `test_clamping.py`. Full suite:
1997 → 2003 passed, 27 skipped, 0 failed.

### 8. Custom-solver demo notebook — done

- [x] Originally shipped as `demos/model_api/custom_solver_demo.py` (a
      `.py` script, matching the sibling `solver_comparison.py`
      convention). **Corrected post-ship** (2026-07-08, at user request)
      to an actual notebook — the checklist item's own title always
      said "notebook" — at
      [`demos/features/SolverProtocols/01_writing_a_custom_solver.ipynb`](../../demos/features/SolverProtocols/01_writing_a_custom_solver.ipynb),
      generated from
      [`demos/features/SolverProtocols/_generate_notebooks.py`](../../demos/features/SolverProtocols/_generate_notebooks.py)
      following the `demos/features/ChemicalEquilibriumProtocol/`
      generation-script convention (`nb()`/`md()`/`code()` helpers,
      executed via `jupyter nbconvert --execute` under the `biosteam`
      kernel, real captured outputs committed). The `.py` script was
      deleted (clean break, no shim) rather than kept alongside a
      duplicate notebook; `docs/solvers.md`'s cross-reference updated
      to point at the notebook. **Expanded further** (same day, second
      user request) with a companion
      [`0_README.ipynb`](../../demos/features/SolverProtocols/0_README.ipynb)
      architecture-overview notebook (pure markdown, condensed from
      `docs/solvers.md`) — matching `ChemicalEquilibriumProtocol`'s
      `0_README.ipynb` + numbered-notebook shape — and the walkthrough
      renamed from `step_and_system_solver_demo.ipynb` to
      `01_writing_a_custom_solver.ipynb` to fit the new two-notebook
      numbering.
- [x] `SpeciationAfterFeedStepSolver`, a custom `StepSolver` that
      deliberately reverses `SequentialAdvanceSolver`'s documented
      pre-feed-speciation ordering (feed → speciation → reactions →
      transfer, instead of speciation → feed → reactions → transfer).
      Chose a real documented ordering decision to invert (see
      `docs/phases-shipped/ORDERING.md`) rather than an arbitrary
      toy reordering, so the demo teaches an actual design trade-off.
- [x] `AsymmetricThreeStageSystemSolver`, a custom `SystemSolver` —
      asymmetric 3-stage link/CV-advance interleaving (default weights
      `(0.2, 0.5, 0.3)`) following the identical extension pattern
      `StrangSplittingSystemSolver` uses (`sim._apply_links`/
      `cv._advance_unchecked`/`sim._invoke_controllers`), demonstrating
      that a bespoke N-stage split is no more than composing those same
      three calls differently.
- [x] Both run against a real `Simulation(...)`, compared against a
      shipped solver on the same model (`SequentialAdvanceSolver` for
      Axis 1; `ExplicitEulerSystemSolver` for Axis 2).

Sanity check: notebook executed end-to-end via `jupyter nbconvert
--execute` (biosteam kernel, matching the `ChemicalEquilibriumProtocol`
sibling notebooks), 0 error outputs, JSON validity confirmed by
`json.load`. Captured outputs identical to the original script's
manual run: Axis 1 shows pH staying at 7.000 under the shipped default
(speciation pinned pre-feed) vs. shifting to 12.301 under the custom
solver (speciation sees the post-feed state) for the same strong-base
feed and CV — concretely demonstrating that ordering changes the
physical answer, not just internal bookkeeping. Axis 2 shows both
solvers conserving mass and agreeing to within ~6% at `dt_h=0.1`
(finite-dt splitting-error difference, same trade-off
`StrangSplittingSystemSolver` demonstrates against Euler). No dedicated
pytest wrapper — matches the `ChemicalEquilibriumProtocol` notebooks'
own precedent (none of them have one either).

### 9. Documentation + ship

- [x] [`docs/solvers.md`](../../docs/solvers.md) — documented
      `SequentialAdvanceSolver`, `clamp_fn=`, `OrchestrationWarning`,
      the gas/liquid generalization, added the missing Axis-2 section
      this doc had never had, and added a new Clamping section
      consolidating `clamp_fn`/`floor_nonnegative` in one place.
- [x] Fixed the three (not four — one of the four listed items,
      `StepSolver`/`SystemSolver` decoration, and the `_solver_for` fix
      it enables, count as one finding, not two) "minor
      documentation/consistency findings": added `@runtime_checkable`
      to `StepSolver` (`src/core/solvers.py`) and changed
      `Simulation._solver_for` to use `isinstance(solver, SystemSolver)`
      instead of `hasattr(solver, "advance_system")`; corrected
      `compute_jacobian`'s docstring (it isn't called by
      `MonolithicODESolver` *or* `SimultaneousAdaptiveSolver` — grepped
      and confirmed nothing calls this stub at all; the actual
      analytical-Jacobian path uses `engine.jacobian_dz_dy()`, a
      different object); documented the `EventScheduler` vs.
      `MonolithicODESolver`-inline-`ctrl_schedule` duplication in
      `docs/solvers.md`'s new Axis-2 section (no fix, as scoped).
- [x] [`docs/phases-shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md`](../phases-shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md)
      — added a "Shipped" banner; moved from `phases-upcoming/` via
      `git mv`. Also fixed a pre-existing broken same-dir link to
      `SOLVER_ARCHITECTURE.md` (which lives in `docs/design/`, not
      this file's directory, in either its old or new location) while
      touching the file's cross-references for the move.
- [x] [`docs/phases-upcoming/STEP_SOLVER_REFINEMENT_PLAN.md`](STEP_SOLVER_REFINEMENT_PLAN.md) —
      Stage 1 marked shipped; Stages 2–4 left in place (still pending).
- [x] [`docs/phases-upcoming/README.md`](README.md) — replaced the
      "Recently surfaced" entry with a "Solver interface refinement"
      section reflecting Stage 1 shipped / Stages 2–4 pending; added a
      `step-solver-interface-refinement` entry to "Recently shipped".
- [x] Final test sweep: `python -m pytest tests/standalone -q` →
      **2003 passed, 27 skipped, 0 failed** (baseline 1950 passed
      post-rename-merge; +53 net new across checkpoints 1, 2, 3, 4, 5,
      6, 7, 7b).
- [x] Ship: merged `--no-ff` to `main` (commit `de0f898`). Full suite
      re-run on `main` post-merge: 2003 passed, 27 skipped, 0 failed.
- [x] Tag: `step-solver-interface-refinement-shipped` created on `de0f898`.
- [ ] Push: **deliberately deferred at user's request** — merge and tag
      are local only; `main` is 14 commits ahead of `origin/main`
      unpushed. Push when ready with `git push && git push --tags`.
- [ ] Delete branch: **deferred** — `step-solver-interface-refinement`
      still exists locally (no remote copy was ever pushed, same as
      Stage 0). Delete after pushing, with
      `git branch -d step-solver-interface-refinement`.

## Final test count expectation

1950 → 2003 passed, 27 skipped, 0 failed. Confirmed empirically at
checkpoint 9. Purely additive — no existing behavior changes, so no
test needed rewriting (two tests were *renamed*/re-scoped to reflect
intentionally *widened* behavior: `test_snapshot_solver_fails_fast_on_missing_phase`
→ `..._missing_liquid_phase` in checkpoint 4, and the equivalent for
`SimultaneousAdaptiveSolver` in checkpoint 7b — both because the
gas/liquid generalization deliberately made a previously-failing case
succeed).

## Risk notes

- **Checkpoint 2 is the one with real blast radius.** Missing any of
  the four internal call sites means `Simulation.run()` spams
  `OrchestrationWarning` on every step for real users, which is worse
  than the silent bug this checkpoint fixes. The sanity check's
  "zero warnings across all five `SystemSolver`s in a real
  integration test" is load-bearing — don't skip it.
- **Checkpoint 5's `apply_flux` question needs an answer before code
  is written, not during.** The three `Phase` subclasses' baked-in
  floor is easy to miss if checkpoint 5 starts from
  `SequentialAdvanceSolver`'s reified body without re-reading
  `phases.py`. Get the scope decision in writing (in this file) before
  starting the checkpoint.
- **Checkpoint 7 (state-vector unification) touches the widest
  surface for the lowest immediate payoff.** If time-boxed, it's the
  safest checkpoint to defer to a follow-up commit within the same
  branch without blocking the others — nothing downstream in this
  checklist depends on it.
- **Checkpoint 8 (demo)** depends on checkpoint 1's reified default
  existing so the demo can show `SequentialAdvanceSolver` as a real
  class rather than "the thing that runs when `solver=None`."
