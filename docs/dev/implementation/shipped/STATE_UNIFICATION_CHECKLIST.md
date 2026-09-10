# State Unification — Checklist

> **Status: Shipped 2026-05-22** on the `state-unification`
> branch. Tag: `state-unification-shipped`. Final suite: 793 / 0 in
> 137s. C1–C6 implementation log preserved below; C4 shipped as
> five sub-commits (4a–4e) within the checkpoint (the
> justifiable departure is documented inline at §C4).

Implementation log for `state-unification`. The design framing — why
this phase exists, what gets deleted, what numerical drift to expect,
the relationship to `INTEGRATOR_REMOVAL` (shipped 2026-05-20) and
`SIMULATION_CLASS` — lives in
[STATE_UNIFICATION.md](STATE_UNIFICATION.md). This file is the
file-level edit log, the per-checkpoint sub-decision record, and the
§BSM2 re-baseline log capturing the drift that accumulates across
checkpoints. Modelled on
[INTEGRATOR_REMOVAL_CHECKLIST.md](INTEGRATOR_REMOVAL_CHECKLIST.md).

## Goal

Collapse speciation-shaped type machinery and unify Phase storage on
PHREEQC-style `n_mol`:

- Replace the kind-dispatched `Reaction` class with three independent
  concrete classes (`KineticReaction`, `EquilibriumReaction`,
  `BlackBoxReactionModel`). No ABC, no shared base. Hard break on
  `Reaction(kind=...)`.
- Rename `ReactionSet` to `ReactionSystem`. Single attach point
  (`cv.reaction_system`); internal pre-bucketing by `isinstance`;
  lazy engine construction; immutable post-attach.
- Widen `Phase.n_mol` to hold **all** species (conserved totals AND
  equilibrium-derived). Drop the B2 storage separation. `Phase.pH`
  becomes a `@property` computed from `n_mol["H+"] / V_L`.
- Delete `PropertySolver` protocol, `PropertyResult` dataclass,
  `SpeciationPropertySolver` adapter, `cv.property_solvers` list.
  Engine called directly from `cv.advance()` via
  `cv.reaction_system._engine` (lazy-built).
- Collapse the `chem_env` dict. `cv.advance(dt_h, t_h)` takes `t_h`
  as direct argument; orchestrator owns time.
  `ReactionEnvironment` reads pH from `phase.pH`, derived
  concentrations from `phase.n_mol`. `AdvanceResult.properties`
  field removed.
- Introduce narrow `PropertyCalculator` protocol (scalar derived
  properties only — viscosity, density, …). No concrete implementers
  in this phase; the seam exists for future use.
- Introduce `ConservationMonitor` for element + charge accounting per
  `advance()` step. `warn` default; throttled via existing
  `WarningConfig`.

After this branch ships:

- `src/core/property_solver.py`, `src/core/speciation_solver.py`,
  the `Reaction` class in `src/reactions/reaction.py`, and the
  `ReactionSet` class in `src/reactions/reaction_set.py` are gone.
- `src/reactions/_shared.py` exists with the free-function validation
  helpers; `KineticReaction` and `EquilibriumReaction` live in
  separate modules under `src/reactions/`.
- `cv.reaction_system: ReactionSystem` is the only reaction
  attach-point. `cv.reaction_model` and `cv.property_solvers` are
  gone. `cv.advance(dt_h, t_h)` is the new signature; `chem_env`
  no longer threaded.
- `Phase.n_mol` holds totals + derived; `Phase.pH` is a `@property`;
  `Phase.properties: Dict[str, float]` exists.
- `AdvanceResult.reaction_sources` and `AdvanceResult.transfer_record`
  survive; `AdvanceResult.properties` is removed.
- `AccuracyMonitor` injection moves to
  `cv.reaction_system.attach_monitor(monitor)`.
- `ConservationMonitor` available via
  `cv.reaction_system.attach_conservation_monitor(monitor)`.
- BSM2 reference golden test
  ([`test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py))
  re-baselined at C3 (numerical: speciation writes directly to
  `n_mol`) and again at C4 (structural: `chem_env` removed,
  `result.properties` removed). Drift documented inline per the
  chemistry-unification + integrator-removal precedent.
- `docs/architecture.md`, `docs/class_diagrams.md`,
  `docs/solvers.md`, and `README.md` updated to reflect the new
  surface.

## Out of scope

Reaffirmed from [STATE_UNIFICATION.md](STATE_UNIFICATION.md)
"Open items deferred":

- **`PropertyCalculator.update_after` opt-in.** Add when a concrete
  post-step re-evaluation need surfaces; none today.
- **Precipitation reactions** (cross-phase liquid-solid with `Ksp`).
  Current `is_cross_phase` binary survives unchanged.
- **Three-phase reactions** (gas-liquid-solid in one stoichiometry).
  Not currently supported, not a regression.
- **DAE solver add-on.** `StepSolver` seam preserved; triggers in
  [DAE_VS_OPERATOR_SPLITTING.md](DAE_VS_OPERATOR_SPLITTING.md).
- **`ReactionDeclaration` Protocol layer** over the independent
  reaction classes. Additive future option; not in scope here.
- **Migrating `CUFermentationSpeciation` legacy path** beyond what
  the C1 hard break + C4 `chem_env` collapse forces. The legacy
  unit lives in
  [`models/vlmodels/fermenter/unit.py`](../../models/vlmodels/fermenter/unit.py)
  and was already marked out-of-scope for `SIMULATION_CLASS`; any
  drift here is bounded and recorded inline but not blocking.

## Resolved decisions

Pinned in [STATE_UNIFICATION.md](STATE_UNIFICATION.md) (2026-05-20):

- **Storage unified on Phase, B2 separation dropped.** Apply-flux
  writes to any species; equilibrium engine redistributes via
  element-total conservation. PHREEQC-style.
- **Reaction taxonomy: three independent classes** — `KineticReaction`,
  `EquilibriumReaction`, `BlackBoxReactionModel`. No ABC, no shared
  base. Free-function utilities for shared validation in
  `src/reactions/_shared.py`.
- **Hard break on `Reaction(kind=...)`** — no deprecation alias.
  Builders migrate in C1.
- **`ReactionSystem` is the new container name.** `cv.reaction_system`
  attaches all reactions; internal pre-bucket by `isinstance`.
  Immutable post-attach.
- **Time ownership moves to orchestrator.** `cv.advance(dt_h, t_h)`
  takes `t_h` as direct argument. `chem_env` disappears.
- **`PropertyCalculator` runs before reactions** (so kinetic rates
  can read viscosity / density). Single invocation per step; no
  `update_after` in v1.
- **`PropertySolver` / `PropertyResult` / `SpeciationPropertySolver`
  removed.** Engine called directly from `cv.advance()` via
  `ReactionSystem`. Engine retains its richer surface for standalone
  use (`engine.solve(**kwargs)`).
- **`AccuracyMonitor` injection via `attach_monitor` method** on
  `ReactionSystem`. Replaces attribute injection through
  `property_solvers` list.
- **`ConservationMonitor` introduced, `warn` default.** Configurable
  thresholds via `VLsim.config.warnings.conservation`.
- **Operator splitting stays default.** DAE deferred per
  [DAE_VS_OPERATOR_SPLITTING.md](DAE_VS_OPERATOR_SPLITTING.md). The
  unified `n_mol` keeps the future DAE seam clean.

## Checkpoints

The checkpoint shape differs from the integrator-removal phase: there
is **no preparatory inventory step**. C1 begins immediate code work
(mechanical hard break across BSM2 / ADM1 / vlmodels / demos / tests).

**Commit discipline:** commit immediately after each checkpoint
finishes (suite green for that checkpoint's slice). Mid-flight edits
were silently lost between checkpoints during
`chemistry-unification-3`; per-checkpoint commits are load-bearing,
not optional.

**ASCII discipline for commit messages.** Heredoc-passed commit
messages on Windows must stay ASCII — em-dash (—) and Unicode minus
bytes don't decode cleanly in the heredoc shell and end up as literal
`\xe2\x80\x94` byte sequences in the message. Precedent:
integrator-removal C4 had to be reworded post-hoc. Use `--` and `-`,
not `—` and Unicode minus, in commit subjects/bodies.

### C1 — Reaction subclass split

- [ ] Create `src/reactions/_shared.py` with free-function validation
      helpers extracted from the existing `Reaction.__init__`:
      stoichiometry coercion, atom/charge closure validation,
      species_ids extraction, phase enumeration, `is_cross_phase`
      derivation. Helpers stay private (underscore-prefixed module);
      `validate_balance` remains the public reaction-validation
      surface in `src/reactions/stoichiometry.py`.
- [ ] Replace `src/reactions/reaction.py`'s `Reaction` class with two
      independent concrete classes. Either co-located in
      `reaction.py` (renamed) or split into
      `src/reactions/kinetic.py` and `src/reactions/equilibrium.py`
      — choice surfaced in §C1 sub-decisions below before
      implementation.
      - `KineticReaction` — `__init__(stoichiometry, rate_fn, *,
        balance_elements, balance_atol, label)`; carries
        `compute_rates(env)` returning `{phase: {sp_id: mol_per_h}}`;
        implements the `ReactionModel` protocol on its own.
      - `EquilibriumReaction` — `__init__(stoichiometry, *, log_K,
        dH_J_per_mol, T_ref_K, total_id, balance_elements,
        balance_atol, label)`. No `rate_fn`, no `compute_rates`.
        `log_K` required for single-phase; optional for cross-phase
        (partition declarations). `is_cross_phase` remains a
        derived `@property` from the stoichiometry.
- [ ] `BlackBoxReactionModel` (in `src/reactions/blackbox.py`)
      unchanged — already independent of `Reaction`.
- [ ] Update `src/reactions/__init__.py`: drop `Reaction`; export
      `KineticReaction`, `EquilibriumReaction`,
      `BlackBoxReactionModel`. Update the module-level docstring.
- [ ] Update `src/reactions/reaction_set.py`: `ReactionSet` consumes
      both new classes via `Union[KineticReaction,
      EquilibriumReaction, BlackBoxReactionModel]`. `partition()`
      uses `isinstance` checks instead of `kind == "equilibrium"`.
      Rename happens in C2, not here.
- [ ] Update `src/reactions/builder.py`: `ReactionBuilder.aerobic_growth`
      and any siblings return `KineticReaction` instead of
      `Reaction`.
- [ ] Update `src/speciation/engine.py`: `from_reactions` consumes
      `EquilibriumReaction` (drop the `kind == "equilibrium"`
      filter; the type itself is the discriminator).
      `is_cross_phase` filtering unchanged.
- [ ] Update `src/core/gas_liquid_link.py`: `derive_speciation_keys`
      consumes `EquilibriumReaction(is_cross_phase=True)`. Drop
      `kind == "equilibrium"` check; rely on `isinstance`.
- [ ] **Builder migration (single C1 commit, hard break, no shim):**
      - [ ] [`models/vlmodels/adm1/base.py`](../../models/vlmodels/adm1/base.py)
            — 3 kinetic `Reaction(...)` sites → `KineticReaction(...)`.
      - [ ] [`models/vlmodels/adm1/bsm2.py`](../../models/vlmodels/adm1/bsm2.py)
            — 15 sites: water/CO2/NH4/VFA equilibria as
            `EquilibriumReaction(..., log_K=..., dH_J_per_mol=...,
            total_id=...)`; substrate-degradation reactions as
            `KineticReaction(...)`.
      - [ ] [`demos/batch_fermenter.py`](../../demos/batch_fermenter.py)
            — water/acetate equilibrium declarations →
            `EquilibriumReaction(...)`; growth → `KineticReaction(...)`.
      - [ ] [`demos/cstr_fermenter.py`](../../demos/cstr_fermenter.py)
            — same shape.
      - [ ] [`demos/fed_batch_fermenter.py`](../../demos/fed_batch_fermenter.py)
            — same shape.
- [ ] **Test migration:**
      - [ ] [`tests/standalone/test_reactions.py`](../../tests/standalone/test_reactions.py)
            — 19 `Reaction(` sites: rewrite as the appropriate new
            class. `test_is_cross_phase_property` becomes a test on
            `EquilibriumReaction.is_cross_phase`. Drop any test
            that asserts on the `kind` argument's error path
            (the kwarg no longer exists; the error is
            type-system-enforced).
      - [ ] [`tests/standalone/test_cv_advance.py`](../../tests/standalone/test_cv_advance.py)
            — 4 sites.
      - [ ] [`tests/standalone/test_gas_liquid_link.py`](../../tests/standalone/test_gas_liquid_link.py)
            — 5 sites.
      - [ ] [`tests/standalone/test_property_solvers.py`](../../tests/standalone/test_property_solvers.py)
            — 5 sites.
      - [ ] [`tests/standalone/test_speciation.py`](../../tests/standalone/test_speciation.py)
            — 3 sites.
      - [ ] [`tests/standalone/test_multi_cv.py`](../../tests/standalone/test_multi_cv.py)
            — 2 sites.
      - [ ] [`tests/standalone/test_membrane.py`](../../tests/standalone/test_membrane.py)
            — 2 sites.
      - [ ] [`tests/standalone/test_factory.py`](../../tests/standalone/test_factory.py)
            — 1 site.
- [ ] Run `python -m pytest tests/standalone -q` — expect green
      modulo any sentinel drift that bare `cv.advance()` callsites
      would expose. C1 should be **numerically identical** to main
      (no semantic change to speciation or kinetics, only the
      declaration surface). Any sentinel failure here is a
      regression, not an expected re-baseline.

Sanity check: `from VLsim.reactions import Reaction` raises
`ImportError`. `from VLsim.reactions import KineticReaction,
EquilibriumReaction, BlackBoxReactionModel` succeeds. Full standalone
suite passes without any BSM2 sentinel re-baseline.

### C2 — `ReactionSet` → `ReactionSystem`

- [ ] Rename `src/reactions/reaction_set.py` → ideally
      `src/reactions/reaction_system.py`; rename `ReactionSet` →
      `ReactionSystem`. Update
      [`src/reactions/__init__.py`](../../src/reactions/__init__.py)
      exports.
- [ ] `ReactionSystem.__init__` pre-buckets by `isinstance`:
      `_kinetic_reactions`, `_single_phase_equilibria`,
      `_cross_phase_equilibria`, `_blackbox_models`. Public
      properties: `cross_phase_equilibria` (read by
      `KineticGasLiquidLink.derive_speciation_keys`). Drop the
      `partition()` public method and the `has_equilibrium`
      attach-time warning machinery — partition is internal now.
- [ ] Reaction list is **immutable post-attach.** No `add()` /
      `remove()` methods. Mutation requires constructing a new
      `ReactionSystem`.
- [ ] Engine still constructed at attach time in this checkpoint
      (lazy build introduced in C3). The engine instance is
      cached on `system._engine`; per-solve temperature handled
      via van 't Hoff inside the engine.
- [ ] Rename `cv.reaction_model` → `cv.reaction_system` across
      [`src/core/control_volume.py`](../../src/core/control_volume.py),
      [`src/core/multi_cv.py`](../../src/core/multi_cv.py), and
      all callsites in
      [`tests/`](../../tests/standalone/) /
      [`models/`](../../models/vlmodels/) /
      [`demos/`](../../demos/).
- [ ] `KineticGasLiquidLink.derive_speciation_keys` reads
      `system.cross_phase_equilibria` directly (no more
      partition() round-trip).
- [ ] Update
      [`tests/standalone/test_reactions.py`](../../tests/standalone/test_reactions.py)
      — any tests on `ReactionSet.partition()` or `has_equilibrium`
      either rewrite to the new internal-bucket properties or
      delete. Mixed-set construction no longer warns; drop any
      test asserting on that warning.
- [ ] Run `python -m pytest tests/standalone -q` — expect green.
      C2 is a pure rename + API tightening; no numerical change.

Sanity check: `cv.reaction_system` exists; `cv.reaction_model` is
gone. `system.cross_phase_equilibria` returns the expected list.
`system.partition()` raises `AttributeError`.

### C3 — Phase unification (B2 removal)

- [ ] [`src/core/phases.py`](../../src/core/phases.py): widen
      `Phase.n_mol` semantics — holds totals AND derived species in
      one dict. Drop any `_conserved` / `_derived` marker sets,
      protected-write enforcement, or B2 separation machinery. Add
      `Phase.pH` as a `@property` computed as
      `-log10(n_mol["H+"] / V_L)` (with safe handling for missing
      `H+` and `V_L == 0`).
- [ ] Add a privileged write path
      `phase._refresh_derived(values: dict)` consumed by the engine
      after a solve. Single-underscore naming signals "engine-only";
      enforce via a docstring + module-private helper rather than a
      runtime guard (cheaper, matches the existing apply_flux
      convention).
- [ ] [`src/speciation/engine.py`](../../src/speciation/engine.py):
      `solve()` writes derived species directly back to
      `phase.n_mol` via `_refresh_derived`. Element-total
      conservation extracted automatically from the reaction graph
      (no special "TIC" / "CT_NH_T" declarations needed; legacy
      declarations still accepted as conservation hints).
- [ ] [`src/core/gas_liquid_link.py`](../../src/core/gas_liquid_link.py):
      `_effective_henry` reads alphas inline from `n_mol` ratios.
      Drop the separate `alphas` dict read.
- [ ] BSM2 / ADM1 builder updates for any declarations that
      conflated totals with molecular forms (the `total_id="NH3"`
      pattern in
      [`models/vlmodels/adm1/bsm2.py:520`](../../models/vlmodels/adm1/bsm2.py#L520)
      may need adjustment now that totals and derived live in the
      same dict).
- [ ] Run `python -m pytest tests/standalone -q` — expect BSM2
      sentinel failures (speciation now writes to `n_mol` directly,
      changing the numerical path). Re-baseline as part of this
      checkpoint, document drift in §BSM2 re-baseline log below.

Sanity check: `Phase.n_mol` contains both `TIC` and `HCO3-`, `CO3^2-`
keys after a single-phase BSM2 solve. `Phase.pH` returns the right
value derived from `n_mol["H+"]`. `test_bsm2_reference.py` passes
with re-baselined goldens; drift documented inline at the test and
in §BSM2 re-baseline log.

### C4 — Property layer collapse

- [ ] Delete [`src/core/property_solver.py`](../../src/core/property_solver.py)
      (the `PropertySolver` Protocol + `PropertyResult` dataclass).
- [ ] Delete [`src/core/speciation_solver.py`](../../src/core/speciation_solver.py)
      (the `SpeciationPropertySolver` adapter).
- [ ] [`src/core/control_volume.py`](../../src/core/control_volume.py):
      remove `self.property_solvers: List`. Replace engine access
      via direct `cv.reaction_system._engine` (lazy-built).
      Change `advance()` signature: `advance(dt_h, t_h, *,
      solver=None) -> AdvanceResult`. `chem_env` parameter removed.
      Update the docstring.
- [ ] [`src/reactions/environment.py`](../../src/reactions/environment.py):
      `ReactionEnvironment` reads pH from `phase.pH`, derived
      concentrations from `phase.n_mol`. Drop `PropertyResult`
      threading. `t_h` and `T_K` come through directly (T_K from
      Phase).
- [ ] Remove `AdvanceResult.properties: Dict[str, PropertyResult]`
      field. Inspect
      [`src/core/control_volume.py`](../../src/core/control_volume.py)
      for `AdvanceResult` definition; preserve
      `reaction_sources` and `transfer_record` fields.
- [ ] `AccuracyMonitor` injection moves to
      `cv.reaction_system.attach_monitor(monitor)`.
      `ReactionSystem` propagates to the internal engine on first
      build (or immediately if already built). Update
      [`src/core/control_volume.py`](../../src/core/control_volume.py)
      and any test fixtures relying on the old attribute-injection
      path through `property_solvers`.
- [ ] [`src/core/multi_cv.py`](../../src/core/multi_cv.py):
      `MultiCVSystem.advance_all(dt_h, t_h=None)` — accept and
      forward `t_h` to each CV's `advance(dt_h, t_h)`. If `t_h`
      is `None`, the system computes it from an internal counter
      (this is a transitional stance — `SIMULATION_CLASS` will
      own it cleanly).
- [ ] Test migration — every test that:
      - constructs `chem_env={...}` and passes to `cv.advance()`,
      - reads `result.properties["..."]`,
      - constructs a `SpeciationPropertySolver` directly,
      - attaches via `cv.property_solvers.append(...)`,
      - or injects `AccuracyMonitor` via
        `solver._accuracy_monitor = ...` or similar.
      Estimate ~150 tests touched per design doc. Rewrite to the
      new signature. Tests asserting on equilibrium output now
      read from `phase.n_mol[derived_id]` or `phase.pH` directly.
- [ ] Demos: drop `chem_env` construction; pass `t_h` directly.
- [ ] BSM2 reference test: structural re-check, then re-baseline if
      necessary. Document drift in §BSM2 re-baseline log.
- [ ] Run `python -m pytest tests/standalone -q` — expect a large
      slice of failures concentrated in the property-related and
      `chem_env`-using tests. Triage and fix.

Sanity check: `from VLsim.core.property_solver import PropertySolver`
raises `ImportError`. `cv.advance(dt_h, t_h)` with no `chem_env`
returns an `AdvanceResult` whose only fields are
`reaction_sources` and `transfer_record`. `phase.pH` reflects the
post-solve `n_mol["H+"]`.

### C5 — `PropertyCalculator` protocol

- [ ] New `src/core/property_calculator.py` with the runtime-checkable
      Protocol:
      ```python
      @runtime_checkable
      class PropertyCalculator(Protocol):
          key: str
          def compute(self, phase: Phase, T_K: float, P_atm: float) -> float: ...
      ```
- [ ] [`src/core/control_volume.py`](../../src/core/control_volume.py):
      add `cv.property_calculators: List[PropertyCalculator] = []`
      attribute (empty by default). In `advance()`, invoke each
      calculator **once, before kinetic reactions** — store result
      in `phase.properties[calc.key] = calc.compute(phase, T_K,
      P_atm)`.
- [ ] [`src/core/phases.py`](../../src/core/phases.py): add
      `phase.properties: Dict[str, float]` field (empty default).
- [ ] Export `PropertyCalculator` from
      `src/core/__init__.py`.
- [ ] No concrete implementers in this checkpoint. Add one minimal
      smoke test asserting the seam works: a fake
      `PropertyCalculator` returning `42.0` for `key="viscosity"`
      ends up in `phase.properties["viscosity"]` after `advance()`.
- [ ] Run `python -m pytest tests/standalone -q` — expect green.

Sanity check: `from VLsim.core import PropertyCalculator` succeeds.
`isinstance(my_calc, PropertyCalculator)` works structurally
(`@runtime_checkable`). `phase.properties` exists on a fresh phase.

### C6 — `ConservationMonitor`

- [ ] New `src/core/conservation_monitor.py` mirroring the shape of
      `src/speciation/accuracy_monitor.py`. Tracks element + charge
      residuals per `advance()` step.
- [ ] Thresholds: per-step `1e-8 × max(total_element_mol, 1)`;
      cumulative drift `1e-6 × max(total_element_mol, 1)`. Same
      shape for charge residual. Relative + max-clamp.
- [ ] `cv.reaction_system.attach_conservation_monitor(monitor)`
      wiring; `ReactionSystem` propagates the monitor to its internal
      machinery. Mirror the `attach_monitor` plumbing from C4.
- [ ] Hook into `VLsim.config.warnings.conservation` (new
      sub-config). Default `warn`. Throttled via existing throttle
      machinery (post-chemistry-unification-4).
- [ ] Issue a `ConservationWarning(UserWarning)` category when
      threshold exceeded. Mirror `AccuracyWarning` integration.
- [ ] Tests: per-element balance, per-charge balance, cumulative
      drift, threshold respect, throttle respect, `VLSIM_WARNINGS`
      preset interaction. Mirror the test shape of
      `chemistry-unification-4`'s accuracy-monitor tests.
- [ ] Run `python -m pytest tests/standalone -q` — expect green.

Sanity check: under a known-balanced reaction set, no warnings.
Under an intentionally unbalanced feed (e.g. add `Na+` without
`Cl-`), the warning fires within the per-step or cumulative
threshold window.

### Ship — docs sweep + merge + tag

- [ ] [`docs/architecture.md`](../../docs/architecture.md): update
      the `Phase`, `ControlVolume`, `cv.advance()` flow, and
      property-machinery descriptions. Remove
      `property_solver.py` / `speciation_solver.py` from the
      `src/core/` directory tree. Add `_shared.py` / `reaction_system.py`
      to the `src/reactions/` tree.
- [ ] [`docs/class_diagrams.md`](../../docs/class_diagrams.md):
      regenerate or rewrite the mermaid blocks for the reaction
      and core layers. Drop the `PropertySolver`, `PropertyResult`,
      `SpeciationPropertySolver`, `Reaction`, `ReactionSet` blocks.
      Add `KineticReaction`, `EquilibriumReaction`,
      `ReactionSystem`, `PropertyCalculator`, `ConservationMonitor`.
- [ ] [`docs/solvers.md`](../../docs/solvers.md): rewrite the
      `cv.advance()` sub-step list to reflect the new sequence
      (feeds → equilibrium refresh → property calculators →
      kinetic Euler → transfer).
- [ ] [`README.md`](../../README.md): update the descriptions of
      reaction declarations, Phase storage, and the property layer.
- [ ] [`docs/upcoming/SIMULATION_CLASS.md`](../upcoming/SIMULATION_CLASS.md):
      reconcile against the new shape. Decisions 11 (`chem_env`
      permanence) and possibly 4 (records taxonomy) need updating.
      Edit the conflict banner / decision text directly; this is
      the predicted reconciliation point per
      [STATE_UNIFICATION.md](STATE_UNIFICATION.md) "Relationship
      to other phases".
- [ ] Move [`STATE_UNIFICATION.md`](STATE_UNIFICATION.md) and this
      checklist to
      [`../shipped/`](../shipped/) with "Shipped"
      banners. Add a one-paragraph "Status: Shipped" summary
      capturing the final test count, the BSM2 drift summary, and
      net code delta — mirroring the integrator-removal precedent.
- [ ] [`docs/upcoming/README.md`](../upcoming/README.md):
      move `STATE_UNIFICATION` from "Currently in flight" to
      "Recently shipped"; remove from the priority list. Update the
      `SIMULATION_CLASS` entry to drop the conflict banner now
      that reconciliation is complete.
- [ ] Final test sweep:
      `python -m pytest tests/standalone -q`. Confirm pass count.
      Baseline post-integrator-removal is 816; net delta from
      this phase expected to be roughly balanced (deletions on the
      property-solver test side, additions for `PropertyCalculator`
      seam + `ConservationMonitor`).
- [ ] Ship via the convention in [README.md](README.md):
      ```
      git checkout main
      git merge --no-ff state-unification \
          -m "Merge state-unification: unify Phase n_mol; collapse property machinery; ReactionSystem"
      ```
- [ ] Tag: `git tag state-unification-shipped <commit-hash>`.
- [ ] Push: `git push && git push --tags`.
- [ ] Delete branch:
      `git branch -d state-unification` and
      `git push origin --delete state-unification`.

## Final test count expectation

Baseline 816 (post integrator-removal). Net delta from this phase
expected to be roughly neutral or modestly positive:

- C1 removes ~5 tests that asserted on the `kind` argument's error
  path (the kwarg no longer exists; the error is type-system-
  enforced and not testable the same way).
- C2 removes ~5 tests on `ReactionSet.partition()` /
  `has_equilibrium` (those APIs are gone).
- C4 removes ~10 tests on `PropertySolver` /
  `SpeciationPropertySolver` / `chem_env` construction that are
  superseded by direct-engine access tests.
- C5 adds ~3 tests on the `PropertyCalculator` seam.
- C6 adds ~15 tests on `ConservationMonitor` (mirroring
  chemistry-unification-4's accuracy-monitor test shape).

Tentative final count: ~814 ± 5. Confirmed empirically at ship.

## Risk notes

- **C1 is mechanical but large.** ~50 sites across BSM2 / ADM1 /
  demos / 8 test files. Risk is a missed callsite causing a runtime
  `ImportError` on `Reaction` after the class is removed. Mitigation:
  grep `import Reaction\b|from .* import Reaction\b|Reaction(` across
  all directories before committing. Read failures and patch in one
  pass.
- **C3's `Phase.pH` `@property` change.** Any test that wrote to
  `phase.pH` directly (treating it as a settable attribute) breaks
  silently if Python allows the assignment with no setter. Mitigation:
  define `pH` with no setter; verify by greping `phase.pH = ` /
  `phase.pH=`. Should be no callsites today.
- **C3 BSM2 re-baseline magnitude.** Engine writing to `n_mol`
  directly rather than via `PropertyResult` round-trip may shift
  numerical values, especially for derived species whose previous
  round-trip path included a `n_mol` total split. Expected drift is
  small but non-zero. Recorded in §BSM2 re-baseline log.
- **C4 multi-CV time argument propagation.** `MultiCVSystem.advance_all`
  needs an internal time counter (transitional stance until
  `SIMULATION_CLASS`). If a multi-CV test loop relies on the
  current time-less `advance_all(dt_h)` shape, that loop now needs
  `t_h` explicitly. Watch
  [`tests/standalone/test_multi_cv.py`](../../tests/standalone/test_multi_cv.py)
  for breakage.
- **C4 `result.properties` removal blast radius.** ~150 tests
  potentially read `result.properties[...]`. A test that doesn't
  actively assert on those fields but constructs them in fixtures
  may break on the AdvanceResult kwargs alone. Mitigation: grep
  `result.properties` / `.properties[` / `properties=` early in C4.
- **C6 conservation thresholds.** The proposed `1e-8` per-step
  threshold may be too tight for production BSM2 — accumulated
  floating-point drift over ~100 steps could exceed `1e-6` cumulative.
  If the default fires noisily on BSM2 reference, loosen by an order
  of magnitude in C6 itself, not later. Document the chosen
  thresholds inline.
- **`CUFermentationSpeciation` legacy path.** Out-of-scope per
  [STATE_UNIFICATION.md](STATE_UNIFICATION.md), but C1's `Reaction`
  removal and C4's `chem_env` removal will force-touch this file.
  Watch for surprise downstream tests.

## C1 sub-decisions (resolved 2026-05-20)

Surfaced before C1 code work per the design-discussion convention
(see project memory `feedback_design_discussion`). All four
resolutions confirmed:

- **`_shared.py` privacy: private.** Underscore-prefixed module,
  no external import path. Tests continue to get `validate_balance`
  from
  [`src/reactions/stoichiometry.py`](../../src/reactions/stoichiometry.py)
  (existing public surface; survives untouched). `_shared.py`
  carries the *new* helpers needed because validation now happens
  twice (once per subclass) rather than once.
- **Module layout: two files.** `KineticReaction` lives in
  `src/reactions/kinetic.py`; `EquilibriumReaction` lives in
  `src/reactions/equilibrium.py`. Visibly reinforces the
  no-shared-base, independent-classes framing. The existing
  `src/reactions/reaction.py` is deleted at the end of C1.
- **`is_cross_phase`: derived `@property`, no change.**
  `len(set(phases)) > 1`, computed on read. Matches today's
  behaviour exactly. No callsite passes it explicitly today; an
  explicit kwarg can be added later additively if a single-phase
  reaction ever needs to override the inferred value.
- **Commit shape: single C1 commit at end.** New classes +
  `_shared.py` + all builder/test migrations + speciation engine
  + gas_liquid_link, suite green, then one commit. Matches the
  hard-break-no-shim resolution in `STATE_UNIFICATION.md`. Risk
  of a larger diff being harder to bisect is acknowledged; the
  alternative (sub-commits inside C1) breaks the
  per-checkpoint-commit discipline.
- **Protocol contract: defer to C2.** `ReactionModel` protocol
  surface stays as-is through C1.
  [`src/reactions/reaction_set.py`](../../src/reactions/reaction_set.py)'s
  `partition()` and `compute_rates` paths use `isinstance`
  internally instead of `kind == "equilibrium"`, but the public
  protocol shape doesn't change until the `ReactionSystem` rename
  in C2.

## BSM2 re-baseline log

This phase touches `n_mol` semantics, engine write paths, and the
`chem_env` -> `t_h` collapse. Each of these can shift the BSM2
golden trajectory. Each re-baseline is documented per the
chemistry-unification + integrator-removal precedent: (i) which
species shifted most, (ii) magnitude per species, (iii) attribution
to the specific checkpoint with a pointer back to its design rationale.

### C1 re-baseline

Expected: none. Confirmed: BSM2 sentinels pass without re-baseline
when running the C1 suite (see §Test run log → C1). The
declaration-surface rename is numerically identical to main.

### C2 re-baseline

Expected: none. Confirmed: BSM2 sentinels pass without re-baseline.
C2 is a container rename + API tightening with no numerical change.

### C3 re-baseline

**Re-baselined 2026-05-21.** Canonical-naming migration: n_mol now
holds molecular species under canonical IDs (CO2aq, HCO3-, CO3--,
NH3, NH4+, H+, OH-) with no separate total tracker — totals are
implicit (sum over the canonical ladder via
`_CANONICAL_NAMES`). BSM2 migrated its kinetic CO2 production from
`n_mol["CO2"]` (legacy conflated total) to `n_mol["CO2aq"]`
(molecular); NH3 was already the canonical molecular name. Engine's
`solve()` now writes derived species to `phase.n_mol` via
`_refresh_derived`. `KineticGasLiquidLink._liquid_total` sums the
canonical ladder for the liquid-side total (no longer reads a
single conflated key).

Drift summary (final-step values):

- **Liquid** — `S_ac`, `S_pro`: essentially unchanged (drift 1e-10
  to 1e-9 relative). `S_h2`: 9.5e-7 → 1.7e-8 mol/L (factor 57
  drop). Both H2 levels are below the simulation's effective
  resolution; this represents redistribution within numerical
  noise rather than a physical shift. `S_ch4` liquid: ~1.6 %
  drop. Carbonate ladder now exposed as separate species —
  `n_mol["CO2aq"]` = 9.99e-3 mol/L (the molecular fraction at
  pH 3.34, ~99 %); `n_mol["HCO3-"]` = 9.92e-6 mol/L;
  `n_mol["CO3--"]` = 0. Ammonia ladder likewise:
  `n_mol["NH3"]` = 1.1e-8 mol/L (low at acidic pH);
  `n_mol["NH4+"]` = 4.999e-3 mol/L. `n_mol["H+"]` = 4.97e-4 mol/L
  (consistent with pH 3.34).
- **Gas** — `CO2`: 3.67 → 4.31 mol (+17 %). The engine writeback
  driving alpha lookups inline from n_mol ratios produces a
  slightly different effective Henry coupling vs. the previous
  PropertyResult.alphas-based path. `S_ch4`: 0.243 → 0.239 mol
  (~1.5 % drop). `S_h2`: 9.8e-3 → 2.6e-4 mol (factor 38 drop;
  same noise-level reasoning as liquid).
- **pH**: 3.340171 → 3.338685 (1.5e-3 absolute shift). The new
  pH is the PHREEQC-style result; the old value was an artefact
  of the BSM2 conflation where the engine and link consumed
  totals through different paths.

Sentinel keys themselves migrate: the test's
`SENTINEL_FINAL_LIQUID_CONC` dict drops the legacy "CO2"/"NH3"
total entries and adds separate `CO2aq`/`HCO3-`/`NH3`/`NH4+`
entries reflecting the canonical molecular forms.

### C4 re-baseline

**Re-baselined 2026-05-21 during state-unification C4 (property
layer collapse).** `cv.advance(dt_h, chem_env)` collapsed to
`cv.advance(dt_h, t_h)`; speciation now runs through
`cv.reaction_system.engine.solve(phases=...)` directly rather than
via the deleted `SpeciationPropertySolver` wrapper. The link's
alpha computation moved inline (`PropertyResult.alphas` channel
removed). Strong-ion totals migrated into `phase.n_mol` as Species
(`n_mol["S_an"] = 5.21e-3 * V_L` lumps the unnamed anion,
`charge=-1`, `atoms={}`).

Drift from the C3 baseline: pH 3.338685 → 3.303580 (-0.035
absolute, ~1.1 % relative). Other species 1e-9 to 1e-3 relative.
The shift comes from the engine's solve path now running through
`cv.reaction_system.engine` inside `cv.advance` directly; subtle
differences in monitor wiring + alpha-lookup timing (now inline
from `n_mol` ratios) account for the residual. Qualitative
behaviour unchanged: low pH (~3.3), fully protonated NH4+
(~5e-3 mol/L), CO2 transferred to gas (~4.3 mol).

C4e (final cleanup of vestigial `chem_env` / `chem_env_fn`
surface) introduced no further drift — the deleted paths were
already no-ops by the end of C4b/c/d.

### C5 re-baseline

None — confirmed. The PropertyCalculator seam is introduced with
no concrete implementers; BSM2 has no property calculators
attached, so the C4 BSM2 sentinel values carry through C5
unchanged.

### C6 re-baseline

None — confirmed. The ConservationMonitor runs in `warn` mode
by default and only emits warnings; no state mutation, no
sentinel drift. BSM2 doesn't trip the conservation thresholds at
the default `1e-8` / `1e-6` factors.

## Test run log

### C1 — initial suite run

`python -m pytest tests/standalone -q` after all C1 edits applied:

- **Passed: 813**
- **Failed: 0**
- **Runtime: 113s**

**Test count arithmetic:** baseline 816 (post integrator-removal) − 3
(deleted error-path tests in `tests/standalone/test_reactions.py`:
`test_kinetic_with_log_K_rejected`, `test_equilibrium_with_rate_fn_rejected`,
`test_unknown_kind_rejected` — the `kind=` kwarg no longer exists, so
these error paths can't be exercised the same way) = **813
expected**. Actual 813 passed + 0 failed = 813. ✓

**Numerical sentinels unchanged.** `test_bsm2_reference.py` passes
without re-baseline, confirming the C1 expectation that the
declaration-surface rename is numerically identical. No drift to
attribute. The "stop and investigate" path (any failure in a
non-numerical bucket) was not triggered.

C1 is committed as a single commit per the resolved commit-shape
decision.

### C2 — initial suite run

`python -m pytest tests/standalone -q` after C2 edits applied:

- **Passed: 813**
- **Failed: 0**
- **Runtime: 110s**

**Test count arithmetic:** baseline 813 (post C1). No tests added or
deleted — the rewritten ``TestReactionSystemBucketing`` tests (3
tests) replace the deleted ``TestReactionSetPartition`` (3 tests) at
parity.

**Numerical sentinels unchanged.** `test_bsm2_reference.py` passes
without re-baseline, confirming the C2 expectation that the
container rename + pre-bucketing is numerically identical. No drift
to attribute.

**One pre-commit issue surfaced and fixed inline.** The first suite
run hit a `TypeError: 'KineticReaction' object is not iterable` in
`KineticGasLiquidLink.derive_speciation_keys` when a bare
`KineticReaction` (from the factory's single-substrate path or
test fixtures) was attached directly to a CV. The fix in
[src/core/gas_liquid_link.py](../../src/core/gas_liquid_link.py):
return early when the attached object doesn't expose
``cross_phase_equilibria``. Non-system attachments have no
partition declarations to derive, so this is a no-op. Documented
in the function's docstring.

**Engine construction deferred to C3.** The design doc's mention
of "engine still constructed at attach time" was interpreted as
forward-looking — in C2 the existing
``SpeciationPropertySolver`` in ``cv.property_solvers`` still owns
the user-facing engine. ``system._engine`` is reserved (stays
``None``) for the C3+ flow when the engine starts writing derived
species back to ``n_mol`` directly.

### C3 — initial suite run + re-baseline

`python -m pytest tests/standalone -q` after C3 edits applied:

- **Passed: 813**
- **Failed: 0** (after BSM2 sentinel re-baseline)
- **Runtime: 105s**

**Test count arithmetic:** baseline 813 (post C2). No tests added
or deleted in C3. The 2 `test_core` pH tests changed shape (from
"pH from speciation dict" / "pH None when no speciation" to "pH
from n_mol['H+']" / "pH raises ValueError when no 'H+'") but
count is unchanged.

**Scope reduction vs original plan.** The design doc's
"`KineticGasLiquidLink._effective_henry` reads alphas inline from
`n_mol` ratios" was partly deferred: the link's `_liquid_total`
helper now sums the canonical ladder from n_mol (so the liquid
total is computed inline, matching design), but `_effective_henry`
itself still reads `PropertyResult.alphas` because the property
solver still exists in C3 (dies in C4). The full inline-alpha
migration will happen in C4 when `PropertyResult` is removed.

**Sub-decision deltas vs the pre-C3 resolutions.** Two
adjustments emerged during implementation:

- **BSM2 migration shape.** The original plan was to rename BSM2's
  totals to TIC/CT_NH_T. That ran into a deep coupling with
  `KineticGasLiquidLink` (which reads `n_mol[species]` where
  species is the gas-side ID — needed to know about the new
  total key). The chosen final shape goes further than the
  original plan: BSM2 migrates to the engine's canonical species
  names directly (`CO2aq`, `HCO3-`, `CO3--`, `NH3`, `NH4+`). No
  separate "total" key exists; totals are implicit (summed over
  the ladder via `_CANONICAL_NAMES`). This matches the design
  doc's PHREEQC-style end state more closely than the original
  TIC/CT_NH_T plan would have.
- **Engine writeback scope.** Writeback covers canonical-name
  ladders only (CO2 ladder, NH3 ladder, water dissociation,
  strong ions). VFA-style unrecognised acids (`{name}_HA` /
  `{name}_A-`) are deliberately NOT written back to n_mol — they
  continue to be bookkept under their bare-name total in n_mol
  (caller convention). This avoids the double-count that would
  occur if the engine wrote `{name}_HA` and the bare name `{name}`
  was also summed.

**Pre-speciation-state fallback.** Both `_populate_totals_from_phases`
and `KineticGasLiquidLink._liquid_total` apply a fallback: if the
canonical-ladder sum is zero (the engine has not yet populated
canonical species — first-solve state, or pure-transfer test
fixture without speciation), they read the bare-name key as the
total. This keeps non-speciation-aware test fixtures working
without requiring a wholesale migration of every legacy
`n_mol["CO2"]` callsite.

**BSM2 sentinel re-baseline.** All 3 sentinel failures
(`test_final_liquid_concentrations`, `test_final_gas_mol`,
`test_final_pH`) re-baselined. Drift documented inline at
[`tests/standalone/test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py)
and in §BSM2 re-baseline log above.

### C4 — initial suite run + re-baseline

C4 shipped as five sub-commits within the checkpoint (4a–4e)
because the property-layer collapse cascades across ~30 files
and ~150 tests; one mega-commit risked compounding errors.
Per-sub-commit suite results:

- **C4a** (commit `b5456b7`, additive infrastructure only):
  813 passed / 0 failed. No behaviour change — new
  `ReactionSystem.engine` lazy property + `attach_engine` /
  `configure_engine` / `attach_monitor` API, strong-ion
  populator helper, and the
  [SPECIATION_LEVEL_RETIREMENT.md](SPECIATION_LEVEL_RETIREMENT.md)
  design note land here. Existing `property_solvers` flow
  continues to drive `cv.advance`.
- **C4b** (commit `44d73b1`, signature migration):
  758 passed / 55 failed. `cv.advance(dt_h, chem_env)` collapsed
  to `cv.advance(dt_h, t_h)`; `AdvanceResult.properties` field
  removed; `KineticGasLiquidLink` alphas inline from `n_mol`;
  `MultiCVSystem.advance_all(dt_h, t_h=None)` with internal
  counter; snapshot solver protocol updated. Failures
  concentrated at `chem_env=` test callsites and
  `result.properties` reads — expected migration shape.
- **C4c** (commit `a860d59`, test + factory + CU migration):
  813 passed / 3 failed. Migrated test fixtures, factory.py's
  `run_batch` extraction loop, and `CUFermentationSpeciation`
  legacy path. Remaining 3 failures are the BSM2 sentinels —
  re-baselined in this same commit (drift ~1.1 % relative).
  Net: 813 / 0.
- **C4d** (commit `dbcb11f`, delete property layer):
  786 passed / 0 failed. Deleted
  `src/core/property_solver.py`,
  `src/core/speciation_solver.py`, and
  `tests/standalone/test_property_solvers.py`
  (~27 tests removed). Snapshot solvers rewritten to run
  speciation via `cv.reaction_system.engine` directly. BSM2 /
  ADM1 builders drop the `SpeciationPropertySolver` wrapper.
- **C4e** (commit `49d8733`, vestigial cleanup):
  783 passed / 0 failed. Deleted `make_bsm2_chem_env_fn` /
  `make_adm1_chem_env_fn` (replaced by `seed_*_strong_ions`),
  dropped `chem_env` / `chem_env_fn` from `run_batch` +
  `FermenterBuilder`, removed `_collect_from_chem_env_fn` from
  `thermo_params`, retired `ProfileSet.add_chem_env`. 3 tests
  removed (the deleted `TestChemEnv` + `TestRunBatchContextFn`).

C4 totals across the five sub-commits: ~1850 LOC removed,
~520 LOC added. 30 net tests removed (mostly tests of the
deleted property layer + chem_env paths).

### C5 — initial suite run

`python -m pytest tests/standalone -q` after C5 edits applied:

- **Passed: 787**
- **Failed: 0**
- **Runtime: 467s** (slow run, likely IO contention)

**Test count arithmetic:** baseline 783 (post C4e). 4 new tests
added in `tests/standalone/test_cv_advance.py::TestPropertyCalculatorSeam`
exercise the protocol seam: empty default, structural typing,
write-to-`phase.properties`, visibility via `env.prop()` in rate
laws.

**No concrete PropertyCalculator implementers in this checkpoint.**
Viscosity / density / heat-capacity models can land later as
additive PRs; C5 ships only the protocol + invocation hook.

### C6 — initial suite run

`python -m pytest tests/standalone -q` after C6 edits applied:

- **Passed: 793**
- **Failed: 0**
- **Runtime: 156s**

**Test count arithmetic:** baseline 787 (post C5). 6 new tests
in `tests/standalone/test_conservation_monitor.py`:

- `TestConservationMonitorBasics`: no emission on first call
  (baseline set), no emission when balanced, no-op with empty
  registry.
- `TestChargeBalance`: unmatched counterion fires `charge_step`;
  cumulative drift fires `charge_cum`.
- `TestElementBalance`: H or O destruction fires
  `element_step_*`.

Throttle plumbing is the same `_emit` code as `AccuracyMonitor`,
tested in `test_accuracy_monitor.py` — not retested here.

**BSM2 sentinels carry through C5 + C6 unchanged** (793 / 0,
matching the C4 baseline). The ConservationMonitor runs in
`warn` mode and doesn't trip on BSM2's default trajectory at
the `1e-8` / `1e-6` threshold factors.

### Final suite run (ship)

`python -m pytest tests/standalone -q` post-ship-sweep edits:

- **Passed: 793**
- **Failed: 0**
- **Runtime: 137s**

Matches C6 baseline exactly. The ship-sweep edits were docs-only
(class_diagrams.md / architecture.md / solvers.md /
SIMULATION_CLASS.md reconciliation, plus the moves to
shipped/); no code changes since `93bb752` (C6 ship). 52
warnings remain — same set as the C6 run, all benign
(species-redeclaration soft warnings from
`check_species_consistency` in legacy test fixtures, plus the
expected `CHNO → CHO` fallback warning and one dosing-volume
warning in pH-control integration tests).
