# Phase 5 Checklist — Simplify `GasLiquidVolume`

> **Status: Shipped 2026-04-29.** All five checkpoints landed cleanly;
> 740 standalone tests pass (`python -m pytest tests/standalone -q`).
> All checkboxes below are ticked as the implementation log. The
> open follow-on layering work is captured in
> [CONTAINER_LAYERING.md](../upcoming/CONTAINER_LAYERING.md).

Working checklist for Phase 5 of the [CV refactor](CV_UPDATE.md).
Phases 1–4 were complete before this work began; the tree was green
(`python -m pytest tests/standalone -q` → 740 passing).

## Goal

Reduce `GasLiquidVolume` to a thin wrapper that holds **one**
`ControlVolume` containing both `GasPhase` and `LiquidPhase` with
`KineticGasLiquidLink` as an internal `PhaseInterface` (the dual
protocol added in Phase 3). Remove the two-CV split, the internal
`MultiCVSystem`, and the GLV-shaped tangles inside the snapshot and
scipy solvers.

The public surface (`glv.advance`, `glv.boundaries`, `glv.gas_phase`,
`glv.liquid_phase`, `glv.pH`, `glv.ionic_strength`, `glv.total_mol`,
`glv.apply_external_flux`, `glv.snapshot`) is preserved so callers in
`models/vlmodels/**` and `systems/**` keep working.

## Out of scope

- Renaming `GasLiquidVolume` to `Vessel`.
- Promoting `EulerSnapshotSolver` / `ScipyODESolver` onto
  `CV.advance(solver=...)`.
- Deleting `MultiCVSystem` (it has genuine multi-zone users —
  see [CONTAINER_LAYERING.md](../upcoming/CONTAINER_LAYERING.md)).
- Renaming `context_fn` or `profiles.py::add_context` /
  `apply_context` (deferred to a separate naming pass).
- Any reaction-model, speciation, or chemistry-engine refactor.
- Unifying `cv.advance()` and `EulerSnapshotSolver` semantics.

## Resolved decisions

- **No backwards-compat shims, aliases, or deprecation wrappers.**
  Clean breaks at each step, callers updated together.
- **`MultiCVSystem` stays.** Phase 5 only stops using it inside GLV.
  The class itself remains exported and tested; multi-zone users
  (`tests/standalone/test_multi_cv.py`,
  `src/chemistry/thermo_params.py`) are untouched.
- **`glv.gas_cv` and `glv.liquid_cv` properties are removed.** No
  aliases. The grep covered the live tree; only test files reach in.
- **`GasLiquidVolume` remains a class** holding internal CV +
  boundaries + configurable solver — not reduced to a bare factory
  function. Boundaries and solver are real responsibilities that need
  a home.
- **`KineticGasLiquidLink.compute_flux` gains an `instantaneous=False`
  keyword** so `ScipyODESolver` doesn't need to keep one foot in
  `compute_flow`. This is a small additive extension to the new
  PhaseInterface signature, contained to the link class.

## Pre-flight

- [x] Working tree is on `main` (or a Phase 5 branch you intend to
      use). Check `git status` shows the same set of staged/modified
      files reported at the start of the session — no surprises.
- [x] `python -m pytest tests/standalone -q` passes (740 tests).
      **Stop and investigate if not — a red baseline invalidates every
      checkpoint below.**
- [x] Backup confirmed at
      `C:\Users\k2473520\VLcode_backup_2026-04-27_153216`.

---

# Checkpoints

Each checkpoint is **independently green**. Run
`python -m pytest tests/standalone -q` after each one; if red, stop
and report rather than rolling forward.

## Checkpoint 1 — Add `property_solvers` / `reaction_model` accessors and migrate consumers

**Why first:** Lets the rest of the refactor read these from a stable
public surface, so when GLV's internals change in Checkpoint 5 the
external callers don't need a second edit.

**Edits:**

- [x] [src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py)
      — add two read-only properties on `GasLiquidVolume`:
      ```python
      @property
      def property_solvers(self) -> List:
          return self._cv_liquid.property_solvers

      @property
      def reaction_model(self) -> Any:
          return self._cv_liquid.reaction_model
      ```
      (After Checkpoint 5 these forward to `self._cv.property_solvers`
      / `self._cv.reaction_model`.)
- [x] [src/chemistry/thermo_params.py:374](../../src/chemistry/thermo_params.py#L374)
      `_apply_to_speciation_solver`: replace
      `cv_liq = getattr(glv, "_cv_liquid", None); ... for ps in getattr(cv_liq, "property_solvers", [])`
      with `for ps in getattr(glv, "property_solvers", [])`.
- [x] [src/chemistry/thermo_params.py:658](../../src/chemistry/thermo_params.py#L658)
      `collect_thermo_params` GLV branch: same substitution.
- [x] [src/chemistry/thermo_params.py:836](../../src/chemistry/thermo_params.py#L836)
      — already uses `glv._last_result`; **no change needed**.
- [x] [models/vlmodels/fermenter/profiles.py:512](../../models/vlmodels/fermenter/profiles.py#L512)
      `TemperatureSetpointTarget.apply`: replace
      `for ps in getattr(self._glv._cv_liquid, "property_solvers", [])`
      with `for ps in getattr(self._glv, "property_solvers", [])`.
- [x] [models/vlmodels/fermenter/config/factory.py:548-555](../../models/vlmodels/fermenter/config/factory.py#L548-L555)
      initial speciation seed: replace
      `if hasattr(glv, "_cv_liquid") and glv._cv_liquid.property_solvers:`
      and the inner `glv._cv_liquid.phases` read with
      `property_solvers = getattr(glv, "property_solvers", []); if property_solvers: ... _solver.solve(glv.phases, initial_ctx, dt_h)`.
- [x] [tests/standalone/test_builder.py](../../tests/standalone/test_builder.py)
      — replace every `glv._cv_liquid.property_solvers` →
      `glv.property_solvers`, every `glv._cv_liquid.reaction_model` →
      `glv.reaction_model`. Lines 172, 180, 196, 207, 208, 217, 221,
      232, 241, 264.
- [x] [tests/standalone/test_factory.py](../../tests/standalone/test_factory.py)
      — same substitutions. Lines 124, 125, 204, 217, 218, 234, 241.

**Verify:**

- [x] `python -m pytest tests/standalone -q` → 740 passing.

---

## Checkpoint 2 — Extend `KineticGasLiquidLink.compute_flux` with `instantaneous` kwarg

**Why:** `ScipyODESolver` needs the instantaneous-rate path that
`compute_flow(..., instantaneous=True)` currently provides. Without
this extension, the scipy solver would have to keep calling
`compute_flow` directly even after Phase 5, defeating the
single-protocol cleanup.

**Edits:**

- [x] [src/core/gas_liquid_link.py](../../src/core/gas_liquid_link.py)
      `KineticGasLiquidLink.compute_flux`: change signature to
      ```python
      def compute_flux(
          self, state_a, state_b, dt_h,
          property_results=None, *, instantaneous: bool = False,
      ) -> Dict[str, float]:
      ```
      and forward `instantaneous=instantaneous` into the
      `compute_flow(...)` call at the end of the method.
- [x] Update the docstring to mention the new kwarg and that it
      mirrors `compute_flow`'s behaviour.
- [x] No other call sites need updating yet — `compute_flux` is only
      called from `ControlVolume.step_internal_transfer`, which passes
      no `instantaneous` and so picks up the default `False`.

**Verify:**

- [x] `python -m pytest tests/standalone -q` → 740 passing.

---

## Checkpoint 3 — Rewrite `EulerSnapshotSolver` against single-CV phases

**Why:** Removes the GLV-shaped `_cv_gas` / `_cv_liquid` /
`_gas_cv_key` / `_liq_cv_key` reaches and the synthetic two-CV
`_snap_cvs` dict. After this, the snapshot solver only needs the
phases dict and the link as a `PhaseInterface`.

**Note:** This checkpoint runs *before* Checkpoint 5 deletes the
two-CV split. So during Checkpoint 3, GLV still has `_cv_gas` and
`_cv_liquid` — but we read from `glv.gas_phase` / `glv.liquid_phase`
(which work in both worlds) instead of `glv._cv_gas.phases["gas"]` /
`glv._cv_liquid.phases["liquid"]`.

**Edits in [src/core/solvers.py](../../src/core/solvers.py)
`EulerSnapshotSolver.solve_step`:**

- [x] Replace `snap_gas = glv._cv_gas.phases["gas"].snapshot()` and
      `snap_liq = glv._cv_liquid.phases["liquid"].snapshot()` with
      `snap_gas = glv.gas_phase.snapshot()` and
      `snap_liq = glv.liquid_phase.snapshot()`.
- [x] Drop the `snap_cvs = {glv._gas_cv_key: ..., glv._liq_cv_key: ...}`
      construction. Replace
      `transfer_flow = glv._transfer_link.compute_flow(snap_cvs, dt_h, speciation_override=spec_result)`
      with
      `transfer_flow = glv._transfer_link.compute_flux(snap_gas, snap_liq, dt_h, property_results={"speciation": spec_result} if spec_result is not None else None)`.
- [x] Replace the speciation call's phase argument: change
      `solver.solve({"liquid": snap_liq}, chem_env, dt_h)` to read
      from `glv.property_solvers` and pass `{"liquid": snap_liq}`
      unchanged (same content; just stop reaching through
      `_cv_liquid`).
- [x] Replace the reaction snapshot:
      `snap_liq_cv = glv._cv_liquid.snapshot()` →
      `snap_liq_cv = ControlVolume(phases={"liquid": snap_liq}, property_solvers=list(glv.property_solvers), reaction_model=glv.reaction_model)`.
      Or, simpler: build a minimal CV inline that satisfies
      `compute_reaction_rates`. **Pick whichever is shorter to read.**
- [x] Replace the apply step:
      `glv._cv_gas.phases["gas"].apply_flux(gas_deltas, dt_h)` →
      `glv.gas_phase.apply_flux(gas_deltas, dt_h)`. Same for liquid.
- [x] Replace the final speciation loop: iterate
      `for s in glv.property_solvers` and call
      `s.solve({"gas": glv.gas_phase, "liquid": glv.liquid_phase}, chem_env, dt_h)` —
      i.e. pass the full phases dict, not a fabricated subset.
- [x] In the `LinkFlowRecord`, replace
      `f"{glv._gas_cv_key}.gas"` / `f"{glv._liq_cv_key}.liquid"` with
      `"gas"` / `"liquid"` (the phase keys are sufficient now that
      there is only one CV).

**Verify:**

- [x] `python -m pytest tests/standalone -q` → 740 passing.
- [x] Spot-check: run
      `python -m pytest tests/standalone/test_gas_liquid_volume.py tests/standalone/test_factory.py tests/standalone/test_builder.py -q`
      to ensure the GLV-heavy tests are happy.

---

## Checkpoint 4 — Rewrite `ScipyODESolver` against single-CV phases

**Why:** Same reasoning as Checkpoint 3, for the adaptive path.

**Edits in [src/core/solvers.py](../../src/core/solvers.py)
`ScipyODESolver.solve_step`:**

- [x] Drop the captures `gas_cv_key = glv._gas_cv_key`,
      `liq_cv_key = glv._liq_cv_key`, `cv_liquid = glv._cv_liquid`.
      Replace any uses with `glv.gas_phase` / `glv.liquid_phase` /
      `glv.property_solvers` / `glv.reaction_model`.
- [x] Drop the `_snap_cvs = {gas_cv_key: ..., liq_cv_key: ...}`
      construction inside the closure.
- [x] Replace the in-derivative transfer call:
      `transfer_flow = link.compute_flow(_snap_cvs, dt_h, speciation_override=spec_result, instantaneous=True)`
      → `transfer_flow = link.compute_flux(_temp_gas, _temp_liq, dt_h, property_results={"speciation": spec_result} if spec_result is not None else None, instantaneous=True)`.
- [x] Replace the in-derivative speciation call's phases dict:
      `prop_solver.solve({"liquid": _temp_liq}, ctx, dt_h)` is fine
      to leave as-is (same content; no `_cv_liquid` reach involved).
      Actually verify there's no `cv_liquid` reach inside the closure
      — the speciation solver is captured as `prop_solver` already.
- [x] Replace the final speciation loop:
      `for s in cv_liquid.property_solvers: res = s.solve(cv_liquid.phases, chem_env, dt_h)` →
      `for s in glv.property_solvers: res = s.solve({"gas": glv.gas_phase, "liquid": glv.liquid_phase}, chem_env, dt_h)`.
- [x] Sanity-check the `_StateVector` class — it already operates on
      `glv.gas_phase` / `glv.liquid_phase` and doesn't reach into
      `_cv_*`. **Should need no changes.**

**Verify:**

- [x] `python -m pytest tests/standalone -q` → 740 passing.

---

## Checkpoint 5 — Rewrite `GasLiquidVolume.__init__` as single-CV

**Why:** This is the structural change. Once Checkpoints 3 and 4 are
green, the solvers no longer depend on the two-CV split, so we can
collapse it.

**Edits in [src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py):**

- [x] Rewrite `__init__` to build a single `ControlVolume`:
      ```python
      self._cv = ControlVolume(
          phases={"gas": gas_phase, "liquid": liquid_phase},
          internal_interfaces=[transfer_link],
          property_solvers=list(property_solvers or []),
          reaction_model=reaction_model,
          label=label,
      )
      self._transfer_link = transfer_link
      self._boundaries = list(boundaries or [])
      if solver is not None:
          self._solver = solver
      else:
          from .solvers import EulerSnapshotSolver
          self._solver = EulerSnapshotSolver()
      self._last_result: Optional[GasLiquidAdvanceResult] = None
      self.label = str(label)
      ```
- [x] Delete `_gas_cv_key`, `_liq_cv_key`, `_cv_gas`, `_cv_liquid`,
      `_phase_to_cv`, and `_system` (the internal `MultiCVSystem`).
- [x] Update `phases` property:
      ```python
      @property
      def phases(self) -> Dict[str, Phase]:
          return self._cv.phases
      ```
- [x] Update `gas_phase` / `liquid_phase` properties:
      ```python
      @property
      def gas_phase(self) -> GasPhase:
          return self._cv.phases["gas"]

      @property
      def liquid_phase(self) -> LiquidPhase:
          return self._cv.phases["liquid"]
      ```
- [x] Update `total_mol`:
      ```python
      def total_mol(self) -> Dict[str, float]:
          return self._cv.total_mol()
      ```
- [x] Update `apply_external_flux`:
      ```python
      def apply_external_flux(self, phase_key, flux_mol_per_h, dt_h):
          if phase_key in self._cv.phases:
              self._cv.apply_external_flux(phase_key, flux_mol_per_h, dt_h)
      ```
- [x] Update `property_solvers` and `reaction_model` properties added
      in Checkpoint 1 to forward to `self._cv` instead of
      `self._cv_liquid`.
- [x] Update `snapshot()` to construct via the single-CV path:
      ```python
      def snapshot(self) -> "GasLiquidVolume":
          return GasLiquidVolume(
              gas_phase=self.gas_phase.snapshot(),
              liquid_phase=self.liquid_phase.snapshot(),
              transfer_link=self._transfer_link,
              property_solvers=list(self._cv.property_solvers),
              reaction_model=self._cv.reaction_model,
              boundaries=list(self._boundaries),
              label=self.label,
              solver=self._solver,
          )
      ```
      Drop the `gas_cv_key` / `liquid_cv_key` arguments — they no
      longer exist.
- [x] Drop the `gas_cv_key` / `liquid_cv_key` constructor arguments
      and the `gas_cv` / `liquid_cv` properties. **Clean break, no
      aliases.**
- [x] Update [src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py)'s
      module docstring and the `GasLiquidVolume` class docstring to
      describe the new single-CV internal shape.
- [x] [src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py)
      `GasLiquidVolume.create` classmethod: drop the
      `gas_cv_key="_gas"` / `liquid_cv_key="_liquid"` locals — replace
      with `gas_cv_key="gas"`, `liquid_cv_key="liquid"` *as
      arguments to `KineticGasLiquidLink`* so the link's CVLink
      protocol fields still resolve. Verify the link's
      `compute_flow` path (which Checkpoint 3 stopped using) still
      works for any external multi-CV caller of the link — but since
      no one outside GLV uses this link, this is purely a hygiene
      check.
- [x] Drop the `from .multi_cv import MultiCVSystem,
      MultiCVAdvanceResult` import. `MultiCVAdvanceResult` is unused;
      `MultiCVSystem` is only referenced from the removed
      `self._system`.

**Test updates:**

- [x] [tests/standalone/test_gas_liquid_volume.py:104-110](../../tests/standalone/test_gas_liquid_volume.py#L104-L110):
      replace
      ```python
      assert glv.gas_cv is not None
      assert "gas" in glv.gas_cv.phases
      ...
      assert glv.liquid_cv is not None
      assert "liquid" in glv.liquid_cv.phases
      ```
      with
      ```python
      assert "gas" in glv.phases
      assert isinstance(glv.phases["gas"], GasPhase)
      assert "liquid" in glv.phases
      assert isinstance(glv.phases["liquid"], LiquidPhase)
      ```
      (Adjust imports if `GasPhase` / `LiquidPhase` aren't already
      imported in that test file.)

**Verify:**

- [x] `python -m pytest tests/standalone -q` → 740 passing.
- [x] Run a representative `systems/` script end-to-end:
      `python systems/batch_fermenter.py` should complete and print
      the same final pH / pressure as on `main`.

---

## Final verification

- [x] `python -m pytest tests/standalone -q` → 740 passing.
- [x] `git status` shows only the files in the audit, plus
      [PHASE5_CHECKLIST.md](PHASE5_CHECKLIST.md) itself.
- [x] No remaining references to `_cv_gas`, `_cv_liquid`,
      `_gas_cv_key`, `_liq_cv_key`, `_phase_to_cv`, or `_system`
      *as attributes* in the live tree (worktree files under
      `.claude/worktrees/**` are out of scope):
      ```
      Grep _cv_gas|_cv_liquid|_gas_cv_key|_liq_cv_key|_phase_to_cv|_system\b
        in src/, models/, tests/, systems/
      ```
      Expected hits: zero in `src/core/gas_liquid_volume.py` and
      `src/core/solvers.py`; zero in `models/vlmodels/**`; zero in
      `tests/standalone/test_*.py`; zero in `systems/**`.
- [x] `gas_cv` / `liquid_cv` removed from
      [src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py).
- [x] [docs/architecture.md](../architecture.md) updated note: GLV
      now holds a single CV; the two-CV pattern is historical.

---

# Audit appendix — files touched

This is the blast radius identified before Checkpoint 1 starts. Used
as the reference list for the verification step above.

**Source — touched in Phase 5:**

- [src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py) —
  rewrite (Checkpoint 1 + 5).
- [src/core/solvers.py](../../src/core/solvers.py) — rewrite both solvers
  (Checkpoints 3 + 4).
- [src/core/gas_liquid_link.py](../../src/core/gas_liquid_link.py) — extend
  `compute_flux` (Checkpoint 2).
- [src/chemistry/thermo_params.py](../../src/chemistry/thermo_params.py) —
  drop `_cv_liquid` reads (Checkpoint 1).

**Models / systems — touched in Phase 5:**

- [models/vlmodels/fermenter/profiles.py](../../models/vlmodels/fermenter/profiles.py) —
  drop `_cv_liquid` read (Checkpoint 1).
- [models/vlmodels/fermenter/config/factory.py](../../models/vlmodels/fermenter/config/factory.py) —
  drop `_cv_liquid` read (Checkpoint 1).

**Tests — touched in Phase 5:**

- [tests/standalone/test_gas_liquid_volume.py](../../tests/standalone/test_gas_liquid_volume.py) —
  drop `gas_cv` / `liquid_cv` assertions (Checkpoint 5).
- [tests/standalone/test_builder.py](../../tests/standalone/test_builder.py) —
  swap `_cv_liquid.x` → `glv.x` (Checkpoint 1).
- [tests/standalone/test_factory.py](../../tests/standalone/test_factory.py) —
  swap `_cv_liquid.x` → `glv.x` (Checkpoint 1).

**Source — *not* touched:**

- [src/core/control_volume.py](../../src/core/control_volume.py) — already
  Phase-3-clean.
- [src/core/multi_cv.py](../../src/core/multi_cv.py) — `MultiCVSystem`
  remains for genuine multi-zone work.
- [src/core/boundaries.py](../../src/core/boundaries.py),
  [systems/_controllers.py](systems/_controllers.py) — only mention
  `GasLiquidVolume` in docstrings.
- [models/vlmodels/adm1/base.py](../../models/vlmodels/adm1/base.py),
  [models/vlmodels/adm1/bsm2.py](../../models/vlmodels/adm1/bsm2.py),
  [models/vlmodels/adm1/bsm2_direct.py](../../models/vlmodels/adm1/bsm2_direct.py)
  — construct GLV via the public surface; no `_cv_*` reach.
- [models/vlmodels/fermenter/config/builder.py](../../models/vlmodels/fermenter/config/builder.py)
  — high-level builder; uses `GasLiquidVolume` as a return type only.
- [tests/standalone/test_multi_cv.py](../../tests/standalone/test_multi_cv.py),
  [tests/standalone/test_gas_liquid_link.py](../../tests/standalone/test_gas_liquid_link.py),
  [tests/standalone/test_membrane.py](../../tests/standalone/test_membrane.py) —
  do not depend on GLV internals.
