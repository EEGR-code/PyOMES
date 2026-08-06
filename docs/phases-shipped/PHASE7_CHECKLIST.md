# Phase 7 Checklist — Remove `GasLiquidVolume`

> **Status: Shipped 2026-05-12** — all five checkpoints landed.  712
> standalone tests passing post-deletion (`python -m pytest
> tests/standalone -q`).

Working checklist for Phase 7 of the
[CV refactor](CV_UPDATE.md).  After Phase 6
(see [PHASE6_CHECKLIST.md](PHASE6_CHECKLIST.md)),
`GasLiquidVolume` was reduced to a near-transparent shim — solvers
and boundaries live on `ControlVolume`, the unified `AdvanceResult`
carries `transfer_record` and `boundary_records`, and the
`context` → `chem_env` rename had shipped.  This phase deleted the
class entirely with mostly-mechanical caller migration.

Conceptual framing is in
[GLV_REMOVAL.md](GLV_REMOVAL.md).  The deferred unified-history-
object idea (option B for the pH question) is captured in
[../phases-upcoming/RUN_HISTORY.md](../phases-upcoming/RUN_HISTORY.md).
This file is the implementation plan.

## Goal

Delete `GasLiquidVolume` as a class.  After this phase, every
caller (`models/vlmodels/`, `systems/`, `tests/`) holds a plain
`ControlVolume` configured with the right phases, internal
interfaces, and boundaries.  The link-wiring logic that today lives
in `GasLiquidVolume.create()` moves into
`FermenterFactory.create_volume()` — no parallel free-function
factory is introduced (since `.create()` has only one live caller).

As side effects:

- The `pH` / `ionic_strength` accessors are dropped, not re-homed.
  Callers read `result.properties["speciation"].pH` off the
  `AdvanceResult` returned by `cv.advance(...)`.  See
  [RUN_HISTORY.md](RUN_HISTORY.md) for why the unified-history
  alternative was deferred.
- The `_last_result` cache is dropped (lived only on GLV; not
  re-homed on CV).
- The `transfer_link` accessor is dropped.  Controllers and other
  callers that need the gas-liquid link reach into
  `cv.internal_interfaces` directly with an `isinstance(...)`
  lookup (option 6b).
- Variables holding the post-removal `ControlVolume` are swept
  from `glv` → `cv`.

After Phase 7, `src/core/gas_liquid_volume.py` is deleted; the
public surface in `src/core/__init__.py` no longer exports
`GasLiquidVolume`.

## Out of scope

- **Chemistry Unification.**  [CHEMISTRY_UNIFICATION.md](CHEMISTRY_UNIFICATION.md)
  installs `phase.pH` / `cv.current_pH()` with staleness
  guardrails.  Phase 7 does not pre-install these; it just drops
  the GLV accessors and lets Chemistry Unification install the
  proper accessors from scratch when it lands.
- **Unified `RunResult` history object.**  Captured in
  [RUN_HISTORY.md](RUN_HISTORY.md) as a deferred follow-on.  Phase
  7 takes the minimal path; callers that need cumulative history
  build their own arrays as today.
- **`KineticGasLiquidLink` rename.**  Kept as-is — name describes
  physics, not container.
- **Resurrecting `GasLiquidAdvanceResult`.**  Dropped in Phase 6.
  Stays dropped.  The unified `AdvanceResult` carries everything.
- **Generalising `EulerSnapshotSolver` to arbitrary phase combos.**
  Still expects `"gas"` and `"liquid"` keys.  Out of scope.
- **`HPLCColumn` reframing.**  Its 1D-cell topology and bespoke
  ODE state vector are unchanged.  See
  [CONTAINER_LAYERING.md](CONTAINER_LAYERING.md).

## Resolved decisions

- **`pH` / `ionic_strength` are dropped (option C).**  Callers
  read `result.properties["speciation"].pH` /
  `result.properties["speciation"].ionic_strength` off each
  `AdvanceResult`.  Live caller surface is small (3 non-test
  sites + 2 tests).  The architecturally cleaner option B
  (unified `RunResult`) is preserved as a deferred follow-on in
  [RUN_HISTORY.md](RUN_HISTORY.md).
- **No new free-function factory.**  `GasLiquidVolume.create()`
  has exactly one live caller
  ([../../models/vlmodels/fermenter/config/factory.py:317](../../models/vlmodels/fermenter/config/factory.py#L317)
  inside `FermenterFactory.create_volume()`).  The link-wiring
  logic inlines into `create_volume`'s body.  No
  `make_gas_liquid_cv()` is introduced.
- **`transfer_link` access pattern is option 6b.**  No accessor
  on `ControlVolume`.  Callers that need the link do the lookup:
  ```python
  link = next(
      i for i in cv.internal_interfaces
      if isinstance(i, KineticGasLiquidLink)
  )
  link.set_kLa("O2", new_value)
  ```
- **Variable naming sweep `glv` → `cv`.**  All caller-side
  variables holding the post-removal `ControlVolume` rename to
  `cv` (or `<unit>_cv` like `fermenter_cv` if a multi-CV
  collision is foreseeable, but the live tree has no such
  collisions).
- **`_last_result` is dropped, not re-homed.**  Tied to the
  pH/ionic_strength decision — without those accessors,
  `_last_result` has no consumers.  Callers retain their own
  reference to the most recent `AdvanceResult` if they need it.
- **`test_gas_liquid_volume.py` is deleted, not repurposed.**
  Most of its coverage is duplicated in
  [../../tests/standalone/test_cv_advance.py](../../tests/standalone/test_cv_advance.py)
  (the post-Phase-6 CV-with-`KineticGasLiquidLink` tests) and
  [../../tests/standalone/test_membrane.py](../../tests/standalone/test_membrane.py)
  (which migrates its `GasLiquidVolume(...)` calls to
  `ControlVolume(...)`).  Any unique coverage moves into
  `test_cv_advance.py`.
- **No backwards-compat shims, aliases, or deprecation wrappers.**
  Clean breaks; callers updated together.

## Pre-flight

- [ ] Working tree on a Phase 7 branch (or `main` if you'd rather
      stage and review before branching).
- [ ] `python -m pytest tests/standalone -q` passes (747 tests
      expected).  **Stop and investigate if not — a red baseline
      invalidates every checkpoint below.**
- [ ] Backup confirmed at the existing path
      `C:\Users\k2473520\VLcode_backup_2026-04-27_153216`.  (Take
      a fresh backup before starting if you want a Phase-7-specific
      restore point.)
- [ ] Audit greps for caller surfaces touched in this phase:
      ```
      Grep "GasLiquidVolume" in src/, models/, systems/, tests/
      Grep "glv\.pH|glv\.ionic_strength" in src/, models/, systems/, tests/
      Grep "glv\.gas_phase|glv\.liquid_phase" in src/, models/, systems/, tests/
      Grep "glv\.transfer_link" in src/, models/, systems/, tests/
      Grep "GasLiquidVolume\.create" in src/, models/, systems/, tests/
      Grep "isinstance.*GasLiquidVolume" in src/, models/, systems/, tests/
      ```
      Expected counts (live tree, after Phase 6):
      - `GasLiquidVolume` references: ~25–35 sites
      - `glv.pH`: 3 non-test + 2 test
      - `glv.ionic_strength`: 0 non-test
      - `glv.gas_phase` / `glv.liquid_phase`: ~15 sites in tests + a
        few in systems
      - `glv.transfer_link`: a small number of controller-flavoured
        sites in `models/vlmodels/` and `src/control/`
      - `GasLiquidVolume.create`: 1 site
        ([../../models/vlmodels/fermenter/config/factory.py:317](../../models/vlmodels/fermenter/config/factory.py#L317))
      - `isinstance(..., GasLiquidVolume)`: ~3–5 sites in tests

---

# Checkpoints

Each checkpoint is **independently green** —
`python -m pytest tests/standalone -q` after each, stop and report
on red.

## Checkpoint 1 — Mark Phase 7 in-flight

**Why first:** doc-only ceremony that signals work has begun.
Updating the in-flight pointer before any code change keeps
[../phases-upcoming/README.md](../phases-upcoming/README.md) accurate as soon as
implementation starts, so anyone reading the planning docs sees
the current state.  No source change; pytest run after is just
a no-op verification that the baseline is still green before any
real edits.

**Edits in
[../phases-upcoming/README.md](../phases-upcoming/README.md):**

- [ ] Update the "Currently in flight" section.  Replace the
      "**None.**" line with a pointer to this checklist:
      ```
      **[PHASE7_CHECKLIST.md](PHASE7_CHECKLIST.md)** — *In progress —
      Checkpoint 2 of 5* (GLV removal).  Design note:
      [GLV_REMOVAL.md](GLV_REMOVAL.md).
      ```
      Refresh the status line after each subsequent checkpoint
      completes (e.g. "Checkpoint 3 of 5"), and clear back to
      "**None.**" once the docs move to `phases-shipped/` in C5.

**Verify:**

- [ ] `python -m pytest tests/standalone -q` → 747 passing
      (unchanged from pre-flight baseline; no source touched).

---

## Checkpoint 2 — Drop `pH` / `ionic_strength` / `_last_result` from `GasLiquidVolume`

**Why second:** the pH-question is the only design-flavoured
change in this phase.  Doing it as an isolated, early checkpoint
means every later checkpoint is a pure mechanical sweep.  The
change is local (one source file plus a handful of caller sites)
and the test surface is small.

**Edits in
[../../src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py):**

- [ ] Remove the `pH` property (currently
      `Optional[float]` reading `_last_result.properties.get("speciation").pH`).
- [ ] Remove the `ionic_strength` property (mirrors `pH`).
- [ ] Remove the `_last_result: Optional[AdvanceResult] = None`
      field assignment in `__init__`.
- [ ] Remove the `self._last_result = result` assignment in
      `advance()`.  `advance()` returns `result` directly (already
      does post-Phase-6).
- [ ] Update class docstring: drop the "convenience accessor"
      paragraph mentioning `glv.pH` / `glv.ionic_strength`.

**Caller migration:**

- [ ] [../../systems/_controllers.py:84](../../systems/_controllers.py#L84):
      `pH = glv.pH` → controller's step callback receives the
      `AdvanceResult` and reads
      `result.properties["speciation"].pH`.  This is a small
      controller-API tweak; the controller's `step` method gains a
      `result` parameter or reads from a passed-in state.
- [ ] [../../systems/twelve_rxn_AD_stage2_3.py:1085](../../systems/twelve_rxn_AD_stage2_3.py#L1085):
      `pH = glv.pH` → same migration; the script's logging code
      reads from the just-returned `AdvanceResult`.
- [ ] [../../tests/standalone/test_gas_liquid_volume.py:217](../../tests/standalone/test_gas_liquid_volume.py#L217)
      and `:221`: drop the two `assert glv.pH is …` tests.  The
      "pH is None before first advance" semantics no longer have a
      meaningful test surface (no cache to be empty).

**Verify:**

- [ ] No remaining `glv.pH` or `glv.ionic_strength` references in
      the live tree:
      ```
      Grep "glv\.pH|glv\.ionic_strength" in src/, models/, systems/, tests/
      ```
      Expected: zero hits.
- [ ] `python -m pytest tests/standalone -q` → 745 passing
      (747 baseline − 2 dropped pH-accessor tests).

---

## Checkpoint 3 — Migrate `glv.gas_phase` / `glv.liquid_phase` callers to `glv.phases[…]`

**Why third:** mechanical caller-side cleanup that doesn't change
any types.  After this checkpoint, callers no longer use the
GLV-specific `gas_phase` / `liquid_phase` named accessors —
they use the CV-compatible `phases["gas"]` / `phases["liquid"]`
dict access that already works on both GLV and CV.  Setting this up
before the type flip in C4 means C4's variable rename
(`glv` → `cv`) is fully mechanical with no accessor changes.

**Caller migration (mechanical replacement):**

- [ ] All `glv.gas_phase` → `glv.phases["gas"]` (or
      `glv["gas"]` via the existing `__getitem__`).  Pick one
      style consistently per file.
- [ ] All `glv.liquid_phase` → `glv.phases["liquid"]` (or
      `glv["liquid"]`).
- [ ] Touched files (audit during this checkpoint):
  - [../../tests/standalone/test_gas_liquid_volume.py](../../tests/standalone/test_gas_liquid_volume.py)
    — many sites (the file is being deleted in C5 but for now its
    asserts run, so update for green).
  - [../../tests/standalone/test_cv_advance.py](../../tests/standalone/test_cv_advance.py)
    — a few sites in the GLV-vs-CV equivalence test.
  - [../../tests/standalone/test_membrane.py](../../tests/standalone/test_membrane.py)
    — assertion code reading membrane-driven state.
  - Any [../../systems/](../../systems/) script that reads
    `glv.gas_phase.X` or `glv.liquid_phase.X`.

**Edits in
[../../src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py):**

- [ ] **No source change in this checkpoint.**  The `gas_phase` /
      `liquid_phase` properties stay on `GasLiquidVolume` until
      C5 (when the whole class deletes).  Removing them now would
      force every caller migration to happen at the same time as
      the source change; splitting them — caller-side migration
      here in C3, source-side property removal in C5 — keeps each
      diff smaller and reviewable independently.

**Verify:**

- [ ] Final grep for `glv\.gas_phase\b|glv\.liquid_phase\b`
      across the live tree.  Expected: zero hits (excluding
      `gas_liquid_volume.py` and any test that tests the GLV
      accessors directly — those tests are dropped in C2 or get
      deleted with the test file in C5).
- [ ] `python -m pytest tests/standalone -q` → 745 passing.

---

## Checkpoint 4 — Switch factories to return `ControlVolume`; sweep `glv` → `cv`

**Why fourth:** with the design-flavoured change (C2) and the
accessor migration (C3) already done, this checkpoint is a single
atomic flip: factories change return types, callers' variable
names sweep to `cv`, and the few remaining GLV-specific surfaces
(`transfer_link.set_kLa`, `isinstance(…, GasLiquidVolume)`)
migrate.

This checkpoint touches the most files but each individual change
is mechanical.

**Edits in
[../../models/vlmodels/fermenter/config/factory.py](../../models/vlmodels/fermenter/config/factory.py):**

- [ ] Replace the `GasLiquidVolume.create(…)` call at line 317
      with explicit `KineticGasLiquidLink` + `ControlVolume`
      construction.  The new body:
      ```python
      from VLsim.core.gas_liquid_link import KineticGasLiquidLink
      from VLsim.core.control_volume import ControlVolume

      link = KineticGasLiquidLink(
          gas_cv_key="gas", gas_phase_key="gas",
          liquid_cv_key="liquid", liquid_phase_key="liquid",
          henry=henry_dict,
          kLa=kLa_dict,
          equilibrium_species=eq_species,
          speciation_keys=speciation_keys or {"CO2": "CO2aq"},
          speciation_corrections=speciation_corrections or {},
          henry_params=henry_raw_params if henry_raw_params else None,
          _label=f"{label}_gl_transfer" if label else "gl_transfer",
      )
      return ControlVolume(
          phases={"gas": gas_phase, "liquid": liquid_phase},
          internal_interfaces=[link],
          boundaries=boundaries or [],
          property_solvers=[spec_solver],
          reaction_model=rxn_model,
          label=label,
      )
      ```
- [ ] Drop `from VLsim.core.gas_liquid_volume import GasLiquidVolume`.
- [ ] Update `FermenterFactory.create_volume`'s return-type
      annotation: `GasLiquidVolume` → `ControlVolume`.
- [ ] Update `BatchSimulationResult` field type annotation from
      `List[AdvanceResult]` (already correct post-Phase-6) — no
      change needed there, but verify `glv` parameter type
      annotations on `run_batch` and helpers update from
      `GasLiquidVolume` to `ControlVolume`, and rename `glv` →
      `cv` in the function signatures and bodies.
- [ ] Update docstrings: every `GasLiquidVolume` reference in the
      file becomes `ControlVolume`.

**Edits in
[../../models/vlmodels/fermenter/config/builder.py](../../models/vlmodels/fermenter/config/builder.py):**

- [ ] Update `FermenterBuilder.build()` return-type annotation
      from `GasLiquidVolume` to `ControlVolume`.  No body change
      required — `build()` calls into `FermenterFactory.create_volume()`
      whose return type just flipped.
- [ ] Drop the `GasLiquidVolume` import if no longer used.

**Edits in
[../../models/vlmodels/adm1/base.py](../../models/vlmodels/adm1/base.py),
[../../models/vlmodels/adm1/bsm2.py](../../models/vlmodels/adm1/bsm2.py),
and any other ADM1 builder:**

- [ ] Rename `build_adm1_glv` → `build_adm1_cv` (and any sibling
      ADM1 builders that follow the same pattern).
- [ ] Update return-type annotations.
- [ ] Update docstrings.
- [ ] Sweep `glv` → `cv` in function bodies.
- [ ] Update callers of these builders in [../../systems/](../../systems/).

**Edits in [../../systems/](../../systems/) scripts:**

- [ ] Sweep `glv` → `cv` variable name across all 12+ scripts.
- [ ] Migrate any remaining `glv.transfer_link.X` to the 6b
      lookup pattern:
      ```python
      from VLsim.core.gas_liquid_link import KineticGasLiquidLink
      link = next(
          i for i in cv.internal_interfaces
          if isinstance(i, KineticGasLiquidLink)
      )
      link.set_kLa(...)
      ```
- [ ] [../../systems/_controllers.py](../../systems/_controllers.py):
      sweep variable names; verify the controller still works
      (the pH read was already migrated in C2 to read from
      `result.properties`).

**Edits in [../../tests/standalone/](../../tests/standalone/):**

- [ ] [../../tests/standalone/test_factory.py](../../tests/standalone/test_factory.py):
      `isinstance(glv, GasLiquidVolume)` → `isinstance(cv, ControlVolume)`.
      Sweep `glv` → `cv`.
- [ ] [../../tests/standalone/test_builder.py](../../tests/standalone/test_builder.py):
      same.
- [ ] [../../tests/standalone/test_membrane.py](../../tests/standalone/test_membrane.py):
      replace `GasLiquidVolume(…)` constructor calls with
      `ControlVolume(…)` plus inline link construction.  Sweep
      `glv` → `cv`.
- [ ] [../../tests/standalone/test_cv_advance.py](../../tests/standalone/test_cv_advance.py):
      the GLV-vs-CV equivalence test in
      `TestCVAdvanceSolverDispatch` no longer has a meaningful
      "GLV path" to compare against.  Either delete that test
      (since the equivalence is now trivial — both paths are the
      same code) or convert it to a self-comparison (snapshot
      solver vs sequential body).

**Edits in [../../src/control/](../../src/control/):**

- [ ] Audit any controllers that take a `GasLiquidVolume` parameter
      and migrate their type annotations and bodies to
      `ControlVolume`.  Use the 6b transfer-link lookup pattern
      where needed.

**Verify:**

- [ ] Final grep for `glv\.transfer_link` across the live tree.
      Expected: zero hits (all migrated to 6b lookup).
- [ ] Final grep for `isinstance.*GasLiquidVolume` across the
      live tree.  Expected: zero hits in the live tree (only the
      class definition itself remains, deleted in C5).
- [ ] Final grep for `\bglv\b` across the live tree.  Expected:
      zero hits in `models/`, `systems/`, `tests/`, `src/`
      excluding `src/core/gas_liquid_volume.py` (the file itself
      is deleted in C5).
- [ ] `python -m pytest tests/standalone -q` → 745 passing.

---

## Checkpoint 5 — Delete the `GasLiquidVolume` class

**Why last:** all callers are migrated by C4.  The class is
unreachable; deleting it is purely cleanup.

**File deletions / edits:**

- [ ] Delete [../../src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py).
- [ ] [../../src/core/__init__.py](../../src/core/__init__.py):
      remove the `GasLiquidVolume` import and the `"GasLiquidVolume"`
      entry from `__all__`.
- [ ] Delete [../../tests/standalone/test_gas_liquid_volume.py](../../tests/standalone/test_gas_liquid_volume.py).
      Most of its coverage duplicates
      [../../tests/standalone/test_cv_advance.py](../../tests/standalone/test_cv_advance.py)
      §`TestCVWithKineticGasLiquidLink` (added in Phase 6's
      Checkpoint 3).  If any unique coverage exists (e.g.
      multi-step kinetic-approach-to-equilibrium), port that one
      test into `test_cv_advance.py` before deleting.

**Documentation refresh:**

- [ ] [../architecture.md](../architecture.md): delete the
      `GasLiquidVolume` section entirely.  Replace with a brief
      "Gas-liquid CV pattern" section explaining that callers
      construct a `ControlVolume` with `{"gas", "liquid"}` phases
      and a `KineticGasLiquidLink` as an internal interface
      (typically via `FermenterFactory.create_volume()`).
- [ ] [../class_diagrams.md](../class_diagrams.md): delete the
      `GasLiquidVolume` box from layer 1.  Update arrows and
      relationships:
      - Layer 1: drop the `_cv`, `_transfer_link`, `_solver`
        composition arrows.
      - Layer 2: drop the `GasLiquidVolume ..> ControlVolume :
        forwards advance() through cv.advance(solver=_solver)`
        arrow.
      - Layer 4: drop the `GasLiquidVolume ..> ControlVolume :
        forwards .boundaries to cv.boundaries` arrow.
      - Update intro text in each layer to drop GLV mentions.
- [ ] [../solvers.md](../solvers.md): final scan for any
      remaining `GasLiquidVolume` mentions; replace with
      `ControlVolume` (the post-Phase-6 doc already says
      "`ControlVolume` (or `GasLiquidVolume`, which today wraps a
      CV)" — drop the parenthetical).
- [ ] Move
      [GLV_REMOVAL.md](GLV_REMOVAL.md) and
      [PHASE7_CHECKLIST.md](PHASE7_CHECKLIST.md) to
      [./](./) with "Status:
      Shipped" banners.  Fix internal cross-links (relative paths
      change).
- [ ] Update [README.md](README.md) priority list (drop
      GLV_REMOVAL; CHEMISTRY_UNIFICATION stays #1;
      CONTAINER_LAYERING stays #2; RUN_HISTORY.md becomes #3).
- [ ] Update [./README.md](./README.md)
      with a Phase 7 entry.
- [ ] Update memory files (`MEMORY.md` index hook,
      `project_cv_refactor.md`).

**Verify:**

- [ ] `python -m pytest tests/standalone -q` — 745+ passing
      (count may rise slightly if any unique
      `test_gas_liquid_volume.py` coverage was ported).
- [ ] No remaining `GasLiquidVolume` references in the live tree:
      ```
      Grep "GasLiquidVolume" in src/, models/, systems/, tests/
      ```
      Expected: zero hits in the live tree (only `phases-shipped/`
      historical docs and `phases-upcoming/RUN_HISTORY.md`'s
      narrative references remain).

---

# Final verification

- [ ] `python -m pytest tests/standalone -q` — 745+ passing.
- [ ] No remaining `GasLiquidVolume`, `glv` variable name,
      `glv.gas_phase` / `glv.liquid_phase` / `glv.transfer_link` /
      `glv.pH` / `glv.ionic_strength` / `glv._last_result` /
      `GasLiquidVolume.create` references in the live tree.
- [ ] [../solvers.md](../solvers.md) has no remaining
      `GasLiquidVolume` mentions other than historical context.
- [ ] [../architecture.md](../architecture.md) `GasLiquidVolume`
      section is gone or rewritten as a "Gas-liquid CV pattern"
      paragraph.
- [ ] [../class_diagrams.md](../class_diagrams.md) no longer
      shows `GasLiquidVolume` in any layer.
- [ ] Move
      [SOLVER_PROMOTION.md](SOLVER_PROMOTION.md) — already in
      `phases-shipped/` — no action.
- [ ] Move
      [GLV_REMOVAL.md](GLV_REMOVAL.md) and this checklist to
      [./](./) with "Status:
      Shipped" banners.
- [ ] Update [README.md](README.md) priority list.
- [ ] Update memory files.

---

# Future considerations (notes for later phases)

Items surfaced during Phase 7 that are deliberately deferred:

1. **Unified `RunResult` history object.**  The architecturally
   cleaner answer to the pH-question (option B).  Captured in
   [RUN_HISTORY.md](RUN_HISTORY.md).  Trigger: a controller
   framework or analysis path that benefits from explicit
   cumulative history rather than ad-hoc parallel arrays in
   `BatchSimulationResult`.
2. **Chemistry Unification's `phase.pH` / `cv.current_pH()`.**
   After [CHEMISTRY_UNIFICATION.md](CHEMISTRY_UNIFICATION.md)
   Phase B lands, callers reading
   `result.properties["speciation"].pH` can switch to the
   guardrailed `phase.pH` accessor.  Phase 7 leaves callers in the
   `result.properties[…]` form; Chemistry Unification migrates
   them.
3. **`fermenter_cv` naming convention for multi-zone scripts.**
   Phase 7 sweeps `glv` → `cv` since the live tree has no
   collisions.  If a future caller mixes single-vessel
   fermentation with `MultiCVSystem` topology, a more specific
   name like `fermenter_cv` or `main_cv` avoids variable
   shadowing in iteration loops.  Convention to apply locally
   when the situation arises; not enforced globally.
4. **`FermenterFactory.create_volume` and `FermenterBuilder.build`
   naming.**  Both still imply a "volume" / a builder that
   "builds" something — appropriate when they returned a
   `GasLiquidVolume`, slightly less so when they return a
   `ControlVolume`.  Renaming is a cosmetic follow-on; not
   required by Phase 7.

---

# Audit appendix — files touched

Reference list for the verification step.

**Source — touched in Phase 7:**

- [../../src/core/gas_liquid_volume.py](../../src/core/gas_liquid_volume.py)
  — drop `pH` / `ionic_strength` / `_last_result` (C2); deleted
  entirely (C5).
- [../../src/core/__init__.py](../../src/core/__init__.py) — drop
  `GasLiquidVolume` export (C5).
- [../../src/control/](../../src/control/) — controller
  type-annotation migrations and 6b transfer-link lookup
  migrations (C4).

**Models / systems — touched in Phase 7:**

- [../../models/vlmodels/fermenter/config/factory.py](../../models/vlmodels/fermenter/config/factory.py)
  — inline `GasLiquidVolume.create()` body, switch return type to
  `ControlVolume`, sweep `glv` → `cv` (C4).
- [../../models/vlmodels/fermenter/config/builder.py](../../models/vlmodels/fermenter/config/builder.py)
  — switch return type to `ControlVolume` (C4).
- [../../models/vlmodels/adm1/base.py](../../models/vlmodels/adm1/base.py),
  [../../models/vlmodels/adm1/bsm2.py](../../models/vlmodels/adm1/bsm2.py)
  — rename `build_adm1_glv` → `build_adm1_cv`, switch return
  types, sweep `glv` → `cv` (C4).
- [../../systems/](../../systems/) — sweep `glv` → `cv` across
  all 12+ scripts; migrate `glv.pH` (C2) and `glv.transfer_link`
  (C4) sites.
- [../../systems/_controllers.py](../../systems/_controllers.py)
  — `glv.pH` migration (C2); `glv` → `cv` sweep (C4).

**Tests — touched in Phase 7:**

- [../../tests/standalone/test_gas_liquid_volume.py](../../tests/standalone/test_gas_liquid_volume.py)
  — drop two `glv.pH` tests (C2); delete the file entirely (C5),
  porting any unique coverage to `test_cv_advance.py`.
- [../../tests/standalone/test_factory.py](../../tests/standalone/test_factory.py)
  — `isinstance(glv, GasLiquidVolume)` → `isinstance(cv,
  ControlVolume)`, `glv` → `cv` sweep (C4).
- [../../tests/standalone/test_builder.py](../../tests/standalone/test_builder.py)
  — same.
- [../../tests/standalone/test_membrane.py](../../tests/standalone/test_membrane.py)
  — replace `GasLiquidVolume(…)` constructors with
  `ControlVolume(…)` + inline link construction; `glv` → `cv`
  sweep (C4).
- [../../tests/standalone/test_cv_advance.py](../../tests/standalone/test_cv_advance.py)
  — drop or rework the GLV-vs-CV equivalence test (C4); port any
  unique `test_gas_liquid_volume.py` coverage (C5).

**Documentation — touched in Phase 7:**

- [../phases-upcoming/README.md](../phases-upcoming/README.md) — update in-flight
  pointer (C1); update priority list (C5).
- [../architecture.md](../architecture.md) — rewrite GLV section
  (C5).
- [../class_diagrams.md](../class_diagrams.md) — drop
  `GasLiquidVolume` from all layers (C5).
- [../solvers.md](../solvers.md) — final GLV-mention scrub (C5).
- [./README.md](./README.md) —
  add Phase 7 entry (C5).
- Move [GLV_REMOVAL.md](GLV_REMOVAL.md) and this checklist to
  [./](./) with status banners
  (C5).

**Source — *not* touched:**

- `src/core/control_volume.py` — Phase 6 already added
  `boundaries` and `solver=` dispatch; nothing more needed here.
- `src/core/solvers.py` — Phase 6 already CV-shaped both solvers.
- `src/core/gas_liquid_link.py` — `KineticGasLiquidLink` stays
  as-is; name kept for physics meaning.
- `src/core/multi_cv.py` — orthogonal to GLV removal.
- `src/speciation/`, `src/equilibria/` — out of scope (Chemistry
  Unification territory).
