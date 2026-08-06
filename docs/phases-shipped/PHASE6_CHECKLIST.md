# Phase 6 Checklist — Promote Solvers to `ControlVolume`

> **Status: Shipped 2026-05-01.** All five checkpoints landed cleanly;
> 747 standalone tests pass (`python -m pytest tests/standalone -q`)
> — 740 baseline + 7 new tests added across Checkpoints 1 and 3.
> All checkboxes below are ticked as the implementation log.  The
> remaining GLV simplification work is captured in
> [GLV_REMOVAL.md](../phases-upcoming/GLV_REMOVAL.md); the larger
> container-layering story stays in
> [CONTAINER_LAYERING.md](../phases-upcoming/CONTAINER_LAYERING.md).

Working checklist for Phase 6 of the
[CV refactor](CV_UPDATE.md). The five-phase CV
refactor shipped 2026-04-29; this phase is the smallest natural
follow-on and the first slice of the larger
[CONTAINER_LAYERING.md](../phases-upcoming/CONTAINER_LAYERING.md) factoring.

Conceptual framing is in
[SOLVER_PROMOTION.md](SOLVER_PROMOTION.md). This file is the
implementation plan.

## Goal

Promote `EulerSnapshotSolver` and `ScipyODESolver` from being
`GasLiquidVolume`-attached strategies to being usable directly via
`cv.advance(solver=...)` on any `ControlVolume`. As a side effect:

- `boundaries` becomes a `ControlVolume`-level concept.
- `GasLiquidAdvanceResult` is dropped; the unified `AdvanceResult`
  carries `transfer_record` and `boundary_records`.
- `GasLiquidVolume` shrinks to a near-transparent wrapper holding
  only the `pH`/`ionic_strength` accessors and the `.create(...)`
  convenience factory.
- `MultiCVSystem.advance_all` exposes per-CV solver choice via a
  `solvers={cv_key: solver}` parameter.
- The remaining `context` → `chem_env` naming holdovers
  (`context_fn`, `add_context`, `apply_context`) are renamed.

The post-Phase-6 GLV is essentially a 30-line shim that a future
phase ([GLV_REMOVAL.md](../phases-upcoming/GLV_REMOVAL.md)) can delete cleanly.

## Out of scope

- **System-wide simultaneous integration across multiple CVs.** The
  multi-CV-aware solver tier (snapshot Euler over a multi-zone
  graph) is not delivered. Per-CV solver choice via
  `MultiCVSystem.advance_all` is. The simultaneous-integration
  capability is the deferred [CONTAINER_LAYERING.md](../phases-upcoming/CONTAINER_LAYERING.md)
  work; **kept noted for future development**.
- **Orchestrator-level adaptive macro `dt_h`.** Solver-internal
  sub-stepping (e.g. `ScipyODESolver`'s adaptive `solve_ivp`) is
  unchanged. Multi-rate integration (different `dt_h` per zone) and
  reactive `dt_h` shrinking on solver failure are **kept noted for
  future development**; no target doc yet.
- **Removal of `GasLiquidVolume`.** The class is reduced to a thin
  shim in this phase but not deleted. Caller migration is captured
  in [GLV_REMOVAL.md](../phases-upcoming/GLV_REMOVAL.md).
- **Generalising `EulerSnapshotSolver` to arbitrary phase
  combos.** It will continue to expect `"gas"` and `"liquid"` keys;
  a fail-fast validation will be added so non-conforming CVs error
  loudly.
- **Refactoring the sequential body of `cv.advance()` into a
  `SequentialAdvanceSolver` strategy.** The sequential body grows a
  boundaries pass but is not packaged as a strategy.
- **`cv.compute_reaction_rates` (read-only sibling of `advance`)**
  stays untouched. It has no integration aspect.
- **`KineticGasLiquidLink` rename.** Named for its physics, not its
  container. Stays.

## Resolved decisions

- **`GasLiquidAdvanceResult` is dropped (Path A).** Its fields move
  onto `AdvanceResult`: `transfer_record: Optional[LinkFlowRecord]`
  and `boundary_records: List[ExternalFluxRecord]`. Rationale:
  matches the no-shims precedent and aligns with the long-term
  direction of removing GLV entirely.
- **Boundaries are a `ControlVolume`-level concept.** Stored as
  `cv.boundaries` (mutable list, like
  `cv.internal_interfaces`). Applied in the sequential body of
  `cv.advance` and in both `StepSolver` implementations.
- **`StepSolver.solve_step` takes a `ControlVolume` directly** (not
  a `GasLiquidVolume`). Returns `AdvanceResult`.
- **`MultiCVSystem.advance_all` extension is included** in this
  phase: gains a `solvers={cv_key: solver}` parameter, per-CV
  forwarding to `cv.advance(solver=...)`. ~5 lines + 1 test.
- **Naming pass is folded in** (context_fn → chem_env_fn,
  add_context → add_chem_env, apply_context → apply_chem_env). All
  caller sites updated.
- **`EulerSnapshotSolver` keeps its `"gas"` / `"liquid"` phase-key
  expectation.** A fail-fast validation will be added at the entry
  to `solve_step` so misconfigured CVs error with a useful message.
- **The legacy `_apply_boundaries` method on `GasLiquidVolume`** is
  dropped (no callers; was already deprecated).
- **No backwards-compat shims, aliases, or deprecation wrappers.**
  Clean breaks; callers updated together.

## Pre-flight

- [x] Working tree on a Phase 6 branch (or `main` if you'd rather
      stage and review before branching).
- [x] `python -m pytest tests/standalone -q` passes (740 tests
      expected). **Stop and investigate if not — a red baseline
      invalidates every checkpoint below.**
- [x] Backup confirmed at the existing path
      `C:\Users\k2473520\VLcode_backup_2026-04-27_153216`. (Take a
      fresh backup before starting if you want a Phase-6-specific
      restore point.)
- [x] Audit grep for direct `solve_step(glv, ...)` calls outside
      `src/core/solvers.py` and `src/core/gas_liquid_volume.py`:
      ```
      Grep "solve_step\(" in tests/, models/, systems/, src/
      ```
      Expected: no live callers reach into the solver protocol
      directly. Any that do will need updating in Checkpoint 2.

---

# Checkpoints

Each checkpoint is **independently green** —
`python -m pytest tests/standalone -q` after each, stop and report
on red.

## Checkpoint 1 — Add `boundaries` to `ControlVolume`

**Why first:** the boundary-handling responsibility is the largest
behavioural change. Doing it before any solver-protocol work means
later checkpoints don't have to reason about boundaries living on
two different classes during the transition.

**Edits in
[../../src/core/control_volume.py](../../src/core/control_volume.py):**

- [x] `__init__` gains a `boundaries: Optional[List] = None`
      parameter, stored as `self.boundaries = list(boundaries or [])`.
      Position after `internal_interfaces` for symmetry.
- [x] Update the `__init__` docstring to describe the new parameter.
- [x] Inside `advance()`, between the existing "external source
      terms" pass and the "reactions" pass, add a boundaries pass:
      ```python
      # 2b. Apply state-dependent boundary fluxes.
      for boundary in self.boundaries:
          flux = boundary.compute_flux(self, dt_h)
          phase_key = boundary.phase_key
          if phase_key in self.phases:
              self.phases[phase_key].apply_flux(flux, dt_h)
      ```
      (Boundaries are computed from current state and applied via
      `apply_flux`. No diagnostic record is collected here — that's
      the snapshot solver's job and is added in Checkpoint 4.)
- [x] Update `advance()` docstring to document the new sequence
      and reference [ORDERING.md](ORDERING.md).
- [x] Update `snapshot()` to copy the boundaries list:
      `boundaries=list(self.boundaries)`.
- [x] Add `self.boundaries` to `__repr__`.

**Edits in
[../../src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py):**

- [x] In `GasLiquidVolume.__init__`, pass `boundaries=boundaries`
      through to `ControlVolume(...)` instead of storing
      `self._boundaries` separately. Drop `self._boundaries`.
- [x] Update the `boundaries` property to read
      `self._cv.boundaries` instead of `self._boundaries`.

**Tests:**

- [x] In a suitable test file (likely
      [../../tests/standalone/test_cv_advance.py](../../tests/standalone/test_cv_advance.py)
      or a new `test_cv_boundaries.py`), add a test:
      `cv.advance(...)` on a single-phase liquid CV with a
      `LiquidFeed` boundary applies the feed flux without a
      `GasLiquidVolume` wrapper. Confirms boundaries work
      CV-natively.
- [x] Add a test: a CV with no boundaries advances exactly as before
      (regression guard).

**Verify:**

- [x] `python -m pytest tests/standalone -q` → 740+ passing (count
      may rise with new tests).

---

## Checkpoint 2 — Update `StepSolver` protocol to take `ControlVolume`

**Why second:** boundaries are now CV-resident from C1, so the
solver can read them from `cv.boundaries` rather than reaching into
GLV.

**Edits in
[../../src/core/solvers.py](../../src/core/solvers.py):**

- [x] Update the `StepSolver` protocol signature:
      ```python
      class StepSolver(Protocol):
          def solve_step(
              self,
              cv: "ControlVolume",
              dt_h: float,
              chem_env: Dict[str, Any],
              external_source_terms: Optional[Dict[str, Dict[str, float]]] = None,
          ) -> "AdvanceResult": ...
      ```
      Note return type is now `AdvanceResult` (will gain
      `transfer_record` / `boundary_records` fields in Checkpoint
      4 — until then, return a `GasLiquidAdvanceResult` from inside
      this checkpoint, but its caller is now CV-typed).
- [x] Update `EulerSnapshotSolver.solve_step` to take `cv`:
  - Replace `glv.gas_phase` / `glv.liquid_phase` with
    `cv.phases["gas"]` / `cv.phases["liquid"]`.
  - Replace `glv._boundaries` with `cv.boundaries`.
  - Replace `glv.property_solvers` with `cv.property_solvers`.
  - Replace `glv.reaction_model` with `cv.reaction_model`.
  - Replace `glv._transfer_link` with the appropriate
    `cv.internal_interfaces[0]` (or, more robustly, find the
    first `KineticGasLiquidLink` in `cv.internal_interfaces`).
  - Add a fail-fast validation at the top of `solve_step`:
    ```python
    if "gas" not in cv.phases or "liquid" not in cv.phases:
        raise ValueError(
            f"EulerSnapshotSolver expects a CV with 'gas' and "
            f"'liquid' phases; got {list(cv.phases.keys())}"
        )
    ```
- [x] Update `ScipyODESolver.solve_step` similarly. The internal
      `_context_fn` capture may need updating; rename also touches
      this in Checkpoint 5.

**Edits in
[../../src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py):**

- [x] Update `GasLiquidVolume.advance` to call
      `self._solver.solve_step(self._cv, dt_h, chem_env, external_source_terms=...)`
      instead of passing `self`. Boundaries no longer pass through
      as a parameter — they're on `cv.boundaries` now.

**Tests:**

- [x] No new tests required for this checkpoint specifically; the
      existing test suite exercises the solver-via-GLV path and
      will catch any regression.

**Verify:**

- [x] `python -m pytest tests/standalone -q` → 740+ passing.

---

## Checkpoint 3 — Add `solver=` dispatch to `cv.advance()` and `MultiCVSystem.advance_all`

**Why third:** with boundaries on CV (C1) and the protocol
CV-shaped (C2), the dispatch is now mechanical.

**Edits in
[../../src/core/control_volume.py](../../src/core/control_volume.py):**

- [x] `advance()` gains a `solver=None` keyword argument.
- [x] When `solver is None`, run the existing sequential body
      unchanged (now boundary-aware from C1).
- [x] When `solver is not None`, dispatch:
      ```python
      if solver is not None:
          return solver.solve_step(
              self, dt_h, chem_env,
              external_source_terms=external_source_terms,
          )
      ```
      Place the dispatch at the top of `advance()` so it short-
      circuits the sequential body.
- [x] Update `advance()` docstring to describe the dispatch
      semantics.

**Edits in
[../../src/core/multi_cv.py](../../src/core/multi_cv.py):**

- [x] `MultiCVSystem.advance_all` gains a
      `solvers: Optional[Dict[str, "StepSolver"]] = None`
      parameter.
- [x] Inside the per-CV advance loop, look up the solver:
      ```python
      cv_results: Dict[str, Any] = {}
      for key, cv in self.cvs.items():
          ctx = chem_envs.get(key, {})
          ext = external_source_terms.get(key, None)
          solver = solvers.get(key) if solvers else None
          cv_results[key] = cv.advance(
              dt_h, chem_env=ctx,
              external_source_terms=ext,
              solver=solver,
          )
      ```
- [x] Update docstring.

**Edits in
[../../src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py):**

- [x] Replace `GasLiquidVolume.advance` body with a forward through
      `cv.advance(solver=...)`:
      ```python
      def advance(self, dt_h, chem_env=None, external_source_terms=None):
          chem_env = dict(chem_env or {})
          result = self._cv.advance(
              dt_h, chem_env=chem_env,
              external_source_terms=external_source_terms,
              solver=self._solver,
          )
          self._last_result = result
          return result
      ```

**Tests:**

- [x] Add a test: bare `ControlVolume` with two phases (`"gas"`,
      `"liquid"`) and a `KineticGasLiquidLink` running under
      `EulerSnapshotSolver` via `cv.advance(solver=...)` produces
      results equivalent to today's `GasLiquidVolume.advance` (same
      phases, same link, same chem_env). Proves the abstraction
      holds without GLV.
- [x] Add a test: `MultiCVSystem.advance_all(solvers={"zone": EulerSnapshotSolver()})`
      dispatches the solver per CV. Easiest with one CV in the
      system; one assertion on the result type.
- [x] Add a test: `cv.advance(solver=EulerSnapshotSolver())` on a
      CV missing `"gas"` raises a `ValueError` with a useful
      message (the C2 fail-fast).

**Verify:**

- [x] `python -m pytest tests/standalone -q` → 740+ passing.

---

## Checkpoint 4 — Drop `GasLiquidAdvanceResult`; unify on `AdvanceResult`

**Why fourth:** with the dispatch in place (C3), result-type
unification is the cleanup. After this, GLV's own result type is
gone and `AdvanceResult` carries everything.

**Edits in
[../../src/core/interfaces.py](../../src/core/interfaces.py):**

- [x] Add to `AdvanceResult`:
      ```python
      transfer_record: Optional["LinkFlowRecord"] = None
      boundary_records: List["ExternalFluxRecord"] = field(default_factory=list)
      ```
- [x] Update the `AdvanceResult` docstring to describe the new
      fields and clarify the difference between `transfer`
      (mass-balance `TransferDiagnostics`) and `transfer_record`
      (gas-liquid link `LinkFlowRecord`) — they sound similar but
      mean different things.
- [x] Add the necessary imports at the top of the file
      (`LinkFlowRecord`, `ExternalFluxRecord`).

**Edits in
[../../src/core/solvers.py](../../src/core/solvers.py):**

- [x] Both solvers stop creating `GasLiquidAdvanceResult` and
      instead return `AdvanceResult` directly with all fields
      populated. The `gas_cv_result` and `liquid_cv_result` content
      simply collapses into the unified result.
- [x] Drop the
      `from .gas_liquid_volume import GasLiquidAdvanceResult`
      import (and any equivalent in the scipy path).

**Edits in
[../../src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py):**

- [x] Delete the `GasLiquidAdvanceResult` class.
- [x] Update `advance()` return type annotation to `AdvanceResult`.
- [x] Drop the `_apply_boundaries` legacy method (no live callers).
- [x] Update `_last_result` type annotation to
      `Optional[AdvanceResult]`. The `pH` and `ionic_strength`
      accessors keep working as-is (they read
      `self._last_result.properties` which both types share).

**Edits in
[../../src/core/__init__.py](../../src/core/__init__.py):**

- [x] Drop `GasLiquidAdvanceResult` from imports and `__all__`.

**Caller updates:**

- [x] Grep for `GasLiquidAdvanceResult` and `gas_cv_result` /
      `liquid_cv_result` references across the live tree:
      ```
      Grep "GasLiquidAdvanceResult|gas_cv_result|liquid_cv_result"
        in src/, models/, tests/, systems/
      ```
- [x] Update each caller to read from the unified `AdvanceResult`.
      Likely sites:
  - [../../models/vlmodels/fermenter/config/factory.py](../../models/vlmodels/fermenter/config/factory.py)
  - [../../tests/standalone/test_gas_liquid_volume.py](../../tests/standalone/test_gas_liquid_volume.py)
  - [../../tests/standalone/test_factory.py](../../tests/standalone/test_factory.py)
  - Possibly [../../src/chemistry/thermo_params.py](../../src/chemistry/thermo_params.py)
- [x] Tests reading the old fields update to read the new unified
      fields.

**Verify:**

- [x] `python -m pytest tests/standalone -q` → 740+ passing.
- [x] Final grep:
      ```
      Grep "GasLiquidAdvanceResult|gas_cv_result|liquid_cv_result"
      ```
      Expected: zero hits in the live tree.

---

## Checkpoint 5 — Naming pass: `context_fn` → `chem_env_fn`, `add_context` → `add_chem_env`, `apply_context` → `apply_chem_env`

**Why last:** mechanical rename, kept separate from structural work
to avoid muddling diffs. Now that the structural work is shipped,
the rename runs as a self-contained sweep.

**Renames:**

- [x] **`context_fn` → `chem_env_fn`:**
  - [x] [../../src/core/solvers.py](../../src/core/solvers.py):
        `ScipyODESolver._context_fn` attribute and
        `context_fn = self._context_fn` capture inside
        `solve_step`.
  - [x] [../../models/vlmodels/fermenter/config/factory.py](../../models/vlmodels/fermenter/config/factory.py):
        `context_fn` parameter on `run_batch` and any helper
        functions; `context_fn = ...` assignments inside.
  - [x] Any `def context_fn(t, i):` factory definitions in
        `systems/` example scripts.
  - [x] Tests that pass `context_fn=...` as a keyword argument.

- [x] **`add_context` → `add_chem_env`** in
      [../../models/vlmodels/fermenter/profiles.py](../../models/vlmodels/fermenter/profiles.py)
      and all callers.

- [x] **`apply_context` → `apply_chem_env`** in
      [../../models/vlmodels/fermenter/profiles.py](../../models/vlmodels/fermenter/profiles.py)
      and all callers.

**Approach:**

- [x] Use `Grep` to enumerate every occurrence of each old name
      across the live tree:
      ```
      Grep "\\bcontext_fn\\b"
      Grep "\\badd_context\\b"
      Grep "\\bapply_context\\b"
      ```
      Word boundary anchors avoid touching unrelated identifiers
      like `context` (which is now `chem_env` everywhere relevant).
- [x] Use `Edit` with `replace_all=true` for the bulk rename
      file-by-file.
- [x] Update any docstrings, comments, and parameter names
      referring to `context_fn` etc.
- [x] Update planning notes /
      [../phases-shipped/](../phases-shipped/) docs only if they
      describe the *current* code; historical plans are frozen.

**Verify:**

- [x] `python -m pytest tests/standalone -q` → 740+ passing.
- [x] Final grep:
      ```
      Grep "\\bcontext_fn\\b|\\badd_context\\b|\\bapply_context\\b"
      ```
      Expected: zero hits in the live tree (excluding
      `phases-shipped/` historical docs and any `legacy/` content).

---

# Final verification

- [x] `python -m pytest tests/standalone -q` — 740+ passing.
- [x] No remaining `_boundaries` or `_apply_boundaries` references
      on `GasLiquidVolume`:
      ```
      Grep "_boundaries|_apply_boundaries" in src/core/gas_liquid_volume.py
      ```
      Expected: zero hits.
- [x] No remaining `GasLiquidAdvanceResult` references in live tree.
- [x] No remaining `context_fn`, `add_context`, `apply_context`
      identifiers in live tree.
- [x] [../solvers.md](../solvers.md) "Future direction" section
      updated to past-tense (the promotion has shipped).
- [x] [../architecture.md](../architecture.md) `GasLiquidVolume`
      section updated: GLV is now a near-transparent shim; solvers
      live on CV.
- [x] Move [SOLVER_PROMOTION.md](SOLVER_PROMOTION.md) and this
      checklist to [../phases-shipped/](../phases-shipped/) with
      "Status: Shipped" banners.
- [x] Update [README.md](../phases-upcoming/README.md) priority list (drop
      SOLVER_PROMOTION; CHEMISTRY_UNIFICATION moves to #1;
      CONTAINER_LAYERING stays #2 with GLV_REMOVAL added as #3).
- [x] Update memory files (`MEMORY.md` index hook,
      `project_cv_refactor.md`).

---

# Future considerations (notes for later phases)

Items surfaced during Phase 6 that are deliberately deferred:

1. **System-wide simultaneous integration across multiple CVs.**
   `MultiCVSystem.advance_all(solvers={...})` lets each CV pick its
   own solver, but the multi-CV-aware solver tier (snapshot Euler
   over the whole graph) is genuinely new capability. Trigger:
   a multi-zone model where the operator-splitting bias between
   inter-zone flow and within-zone kinetics measurably matters.
   See [CONTAINER_LAYERING.md](../phases-upcoming/CONTAINER_LAYERING.md).

2. **Orchestrator-level adaptive macro `dt_h`.** Two flavours:
   reactive shrinking (solver fails → orchestrator retries with
   smaller `dt_h`) and multi-rate (different `dt_h` per zone).
   Neither has infrastructure today. The `StepSolver` protocol
   could grow a richer return type for the reactive case; the
   multi-rate case is much harder. No target doc yet.

3. **Boundary protocol evolution.** `ExternalBoundary` is currently
   a runtime-checkable Protocol with `phase_key`, `label`, and
   `compute_flux(cv, dt_h, instantaneous)`. If a future model needs
   boundaries that span multiple phases, or boundaries that
   participate in the multi-CV-aware solver tier, the protocol
   may need extension.

4. **Generalising `EulerSnapshotSolver` to arbitrary phase
   combos.** Today it expects `"gas"` and `"liquid"`. A solver
   that snapshots and clamps an arbitrary `Dict[str, Phase]`
   would be more reusable. Bigger structural change; gated on a
   real use case.

5. **`GasLiquidVolume` removal.** Phase 6 reduces GLV to a
   ~30-line shim. A future phase deletes it entirely and migrates
   every caller. See [GLV_REMOVAL.md](../phases-upcoming/GLV_REMOVAL.md).

---

# Audit appendix — files touched

Reference list for the verification step.

**Source — touched in Phase 6:**

- [../../src/core/control_volume.py](../../src/core/control_volume.py) —
  add `boundaries` parameter; add boundaries pass to sequential
  `advance()`; add `solver=` dispatch.
- [../../src/core/solvers.py](../../src/core/solvers.py) —
  protocol takes CV; both solvers updated; return `AdvanceResult`;
  add fail-fast for snapshot solver.
- [../../src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py) —
  drop `_boundaries` (moves to CV); drop `GasLiquidAdvanceResult`;
  drop `_apply_boundaries`; `advance()` becomes a forward.
- [../../src/core/interfaces.py](../../src/core/interfaces.py) —
  `AdvanceResult` gains `transfer_record` and `boundary_records`.
- [../../src/core/multi_cv.py](../../src/core/multi_cv.py) —
  `advance_all` gains `solvers={cv_key: solver}` parameter.
- [../../src/core/__init__.py](../../src/core/__init__.py) —
  drop `GasLiquidAdvanceResult` export.

**Models / systems — touched in Phase 6:**

- [../../models/vlmodels/fermenter/config/factory.py](../../models/vlmodels/fermenter/config/factory.py) —
  `context_fn` rename; possibly `gas_cv_result` / `liquid_cv_result`
  reads if any.
- [../../models/vlmodels/fermenter/profiles.py](../../models/vlmodels/fermenter/profiles.py) —
  `add_context` and `apply_context` method renames.
- [../../models/vlmodels/fermenter/config/builder.py](../../models/vlmodels/fermenter/config/builder.py) —
  possibly `context_fn` references.
- Any `systems/` script using `context_fn` or the renamed profile
  methods (audit during C5).

**Tests — touched in Phase 6:**

- [../../tests/standalone/test_cv_advance.py](../../tests/standalone/test_cv_advance.py) —
  new tests for `cv.advance(boundaries=...)` and
  `cv.advance(solver=...)`.
- [../../tests/standalone/test_multi_cv.py](../../tests/standalone/test_multi_cv.py) —
  new test for `advance_all(solvers={...})`.
- [../../tests/standalone/test_gas_liquid_volume.py](../../tests/standalone/test_gas_liquid_volume.py) —
  drop `gas_cv_result` / `liquid_cv_result` reads; possibly
  `_boundaries` references.
- [../../tests/standalone/test_factory.py](../../tests/standalone/test_factory.py) —
  `context_fn` rename in any test that uses it.

**Source — *not* touched:**

- `src/core/control_volume.py::compute_reaction_rates` — read-only;
  no integration aspect; out of scope per resolved decisions.
- `src/core/gas_liquid_link.py` — already dual-protocol from
  Phase 3; no changes needed.
- `src/core/multi_cv.py` apart from the small `advance_all`
  extension.
- Anywhere the `chem_env` rename has already landed.
