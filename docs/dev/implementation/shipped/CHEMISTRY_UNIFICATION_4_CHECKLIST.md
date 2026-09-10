# Chemistry Unification 4 — Checklist

> **Status: Shipped 2026-05-16** — all 10 checkpoints landed
> across the `chemistry-unification-4` branch (10 commits, one per
> checkpoint). 821 standalone tests passing (778 baseline + 43 net
> new in `test_accuracy_monitor.py`). BSM2 reference sentinels
> unchanged — Phase 4 is purely additive on numerics. Five of six
> cheap checks wired; `dt` vs τ_min has a no-op hook with a
> deferred Jacobian-free estimate; splitting-error proxy and
> `cv.run_with_diagnostics(...)` deferred entirely (forward notes
> in `docs/solvers.md`). Tag: `chemistry-unification-4-shipped`.

Working checklist for `chemistry-unification-4`. Conceptual framing is
in [../upcoming/CHEMISTRY_UNIFICATION.md](../upcoming/CHEMISTRY_UNIFICATION.md)
under "Phase F" and the "Equilibrium-Kinetic Coupling and
Accuracy Monitoring" section. The branch slicing is in
[../upcoming/CHEMISTRY_UNIFICATION_PLAN.md](../upcoming/CHEMISTRY_UNIFICATION_PLAN.md)
under "Phase 4 — `chemistry-unification-4`". This file is the
implementation log: file-level edits, ordered checkpoints, and
sanity checks. Modelled on
[CHEMISTRY_UNIFICATION_3_CHECKLIST.md](CHEMISTRY_UNIFICATION_3_CHECKLIST.md).

## Goal

Add automatic numerical-accuracy monitoring to the chemistry stack so
users learn when they're stepping outside the regime where the chosen
solver / activity model is reliable. Two structural pieces share one
package-level configuration surface:

1. **`AccuracyMonitor` + `AccuracyWarning`.** A per-CV monitor that
   runs cheap per-step checks (pH change, Newton iterations, charge
   residual, splitting-error proxy, `dt` vs τ_min, scipy step
   rejections) and emits `AccuracyWarning` when a threshold is
   crossed. One warning category, one meaning: "your numerical setup
   is being used outside the regime where it's reliable."
2. **`VLsim.config` package-level configuration.** A `WarningConfig`
   dataclass holding thresholds + throttle mode; a module-level
   `config` singleton modifiable before runtime; a `VLSIM_WARNINGS`
   environment variable picking from flat presets
   (`silent`/`verbose`/`production`).

After this branch ships:

- [`src/monitoring/accuracy.py`](../../src/monitoring/accuracy.py)
  defines `AccuracyMonitor` (per-CV check runner) and
  `AccuracyWarning(UserWarning)`. The monitor holds per-CV step
  history (`last_pH`, `step_count`, `_warned_categories`) and reads
  global thresholds + throttle off `VLsim.config.warnings` at
  emission time.
- [`src/config.py`](../../src/config.py) defines `WarningConfig`,
  `Config`, the module-level `config` singleton, and `VLSIM_WARNINGS`
  env-var initialisation. Flat presets only: `silent` / `verbose` /
  `production`. Threshold tuning is a Python concern
  (`VLsim.config.warnings.pH_change_threshold = 0.5`); the env-var
  is for CI/batch ergonomics.
- [`src/__init__.py`](../../src/__init__.py) exports
  `AccuracyWarning`, `WarningConfig`, and `config` for users to
  import as `VLsim.config.warnings` / `VLsim.AccuracyWarning` /
  `VLsim.WarningConfig`.
- [`src/core/control_volume.py`](../../src/core/control_volume.py)
  constructs an `AccuracyMonitor` in `__init__` (attached as
  `cv._accuracy_monitor`). The monitor is lazy: thresholds + throttle
  read off `VLsim.config.warnings` per check, so configuration
  changes between simulations take effect with no re-instantiation.
- [`src/core/solvers.py`](../../src/core/solvers.py) `EulerSnapshotSolver`
  and `ScipyODESolver` call into the monitor at the right hook
  points: pH-jump check on the step-closing speciation re-solve;
  scipy step-rejection rate after `solve_ivp` returns; `dt` vs τ_min
  at construction.
- [`src/speciation/engine.py`](../../src/speciation/engine.py) calls
  into the monitor for Newton-iteration count and charge-balance
  residual on each `solve()`. The legacy `_warned_high_I` per-instance
  cache is **deleted**; the existing ionic-strength threshold warning
  migrates to `AccuracyWarning` and uses the new throttle machinery.
- A module-level summary aggregator in `src/monitoring/accuracy.py`
  counts emissions across CVs so the `"once"` throttle's
  end-of-simulation summary (`"AccuracyWarning emitted N times total
  across [categories...]"`) is meaningful in `MultiCVSystem` runs.

## Out of scope

Everything reserved for later phases or trigger-gated deferral:

- **`cv.run_with_diagnostics(...)`.** The plan's "rigorous diagnostics"
  entry point (step-halving, cross-solver comparison) is **deferred
  entirely**. Phase 4 ships only the cheap per-step monitor. The
  doc reference in `docs/solvers.md` becomes a forward note pointing
  at a future "diagnostics phase". Rationale: step-halving + cross-
  solver compare are substantive numerics work (Richardson
  extrapolation, `DiagnosticsResult` schema, per-mode dispatch);
  doing them properly justifies its own phase. Half-implementations
  rot.
- **Parsed `VLSIM_WARNINGS=key=value` env-var.** Flat presets only
  for this phase. If a CI need surfaces for fine-grained env-var
  threshold control later, the parser is additive on top.
- **Severity tiers (`critical`/`warn`/`all`).** No per-check
  severity field; all `AccuracyWarning` emissions are peers.
- **Absorbing configuration / runtime-fallback warnings.** Only the
  one accuracy-shaped warning (`engine.py`'s ionic-strength
  threshold) migrates. `gas_liquid_link.py`'s configuration warnings
  (orphan kLa, overlap, missing speciation key) stay `UserWarning` —
  they're setup-time concerns, not numerical regime breaches. The
  `equilibria/` runtime-fallback warnings, the
  `acid_base.py`/`chemistry_level2.py` convergence warnings, and
  `core/solvers.py:401`'s `freeze_speciation` architectural notice
  all stay as-is.
- **Phase 3b (deferred, trigger-gated).** `ChemistryDatabase` +
  `ThermoFramework` + unified Species-ID emission still wait on a
  concrete run-time chemistry composition pull. Phase 4 does not
  unblock it.
- **`ThermodynamicConfig` deletion / `validate_thermodynamics`
  rewrite.** Still tied to Phase 3b's `ChemistryDatabase` rollout.

## Resolved decisions

Pinned during the planning walkthrough on 2026-05-16:

- **`VLSIM_WARNINGS` env-var: flat presets only.** Values
  `silent` / `verbose` / `production` map to the three classmethod
  constructors on `WarningConfig`. The env-var is for "make this
  run quieter / louder without editing code"; a preset name does
  that cleanly. Per-key threshold overrides
  (`VLSIM_WARNINGS=pH_change=0.5,newton=20`) and severity tiers
  (`VLSIM_WARNINGS=critical`) were considered and rejected — both
  introduce a second config language alongside Python attribute
  assignment, and the in-Python surface already covers fine-grained
  tuning with IDE/typecheck help. Both remain additive future
  extensions.

- **`AccuracyMonitor`: per-CV instance + module-level summary
  aggregator.** Each `ControlVolume` constructs its own
  `AccuracyMonitor` in `__init__` (stored as `cv._accuracy_monitor`).
  The monitor holds per-CV step history (`last_pH`, `step_count`,
  `_warned_categories: set[str]` for "once" throttle). Thresholds
  and throttle mode are **not copied** onto the monitor — they're
  read off `VLsim.config.warnings` at every emission so
  configuration changes propagate without re-instantiation.
  A module-level `_summary_counter: Counter[str]` in
  `monitoring/accuracy.py` aggregates emission counts across CVs
  so `MultiCVSystem` runs get one meaningful end-of-simulation
  summary, not N per-CV summaries. The module-level singleton and
  engine-/solver-attached alternatives were considered and
  rejected: per-CV ownership matches the post-Phase-7 architecture
  (CV already owns `property_solvers`, `reaction_model`,
  `boundaries`); the "once" throttle reads naturally as
  once-per-CV-per-category which is what users expect in
  multi-zone runs.

- **`cv.run_with_diagnostics(...)` deferred entirely.** No stub
  method on `ControlVolume`, no `DiagnosticsResult` dataclass. The
  documentation references in `docs/solvers.md` become forward
  notes ("planned: step-halving and cross-solver comparison
  entry-points for opt-in expensive accuracy verification") with
  no committed signature. Rationale: the cheap per-step monitor is
  Phase 4's load-bearing scope; the rigorous diagnostics are
  substantive numerics work that justifies its own design pass
  when a real use case appears. A stub raising
  `NotImplementedError` was considered and rejected — reserving
  the name invites half-implementations later.

- **Existing warnings: migrate only `engine.py`'s ionic-strength
  threshold to `AccuracyWarning`.** That warning is the
  unambiguous accuracy case in the current codebase ("Davies
  beyond ~0.5 mol/L is unreliable"; "ideal assumption beyond ~0.1
  mol/L is biased"). It migrates from `RuntimeWarning` to
  `AccuracyWarning`; the per-instance `_warned_high_I` cache is
  deleted in favour of the throttle machinery. Two thresholds
  (`ionic_strength_ideal_threshold = 0.10`,
  `ionic_strength_davies_threshold = 0.50`) move onto
  `WarningConfig`. Configuration warnings (`gas_liquid_link.py`)
  and runtime-fallback warnings (`equilibria/`, `acid_base.py`,
  `chemistry_level2.py`) stay where they are — they're not
  numerical regime breaches. Aggressive migration of every
  accuracy-shaped warning was considered and rejected: it dilutes
  the new category's meaning and risks classification disputes
  for borderline cases (`acid_base.py:1057`'s convergence-warning
  is reachable through both legitimate sharp-pH-change and bug
  conditions; migrating it would commit to a meaning the codebase
  isn't sure of yet).

- **Hook points for per-step checks.**
  - **pH-jump:** in `EulerSnapshotSolver.solve_step` after the
    step-closing speciation re-solve; in `ScipyODESolver.solve_step`
    after the post-step `final_props` speciation
    ([`core/solvers.py:628-637`](../../src/core/solvers.py#L628)).
  - **Newton iterations / charge residual:** inside
    `SpeciationEngine.solve()` after the underlying model returns.
    The engine exposes `iters` and `charge_residual` from the
    solver output dict where available; brentq-based solvers
    don't expose iteration counts directly, so this check is
    best-effort and gated on "field present in `out`".
  - **Ionic-strength threshold:** existing site in
    `SpeciationEngine.solve` migrated in place.
  - **`dt` vs τ_min:** computed once at
    `ControlVolume.__init__` from the reaction model's Jacobian
    estimate; emits at that point. The post-Phase-1 reaction model
    exposes max rate magnitude via `reaction_set.partition()`'s
    kinetic bucket; τ_min is a coarse `1 / max_rate` estimate.
    A first-cut implementation may degrade to "no check
    available" when the Jacobian estimate can't be computed; the
    monitor logs no warning in that case.
  - **Scipy step rejections:** after `solve_ivp` returns in
    `ScipyODESolver.solve_step`. The rejection count is
    approximated as `nfev - len(t)` where `t` is the array of
    accepted steps; warning fires if the ratio exceeds a configured
    fraction.
  - **Splitting-error proxy (pH shift across end-of-step re-solve):**
    deferred to a future phase. It needs both the pre-re-solve pH
    and the post-re-solve pH within the same step; currently only
    the final pH is available. Adding the second snapshot requires
    a small solver-internal refactor that isn't justified for
    Phase 4's scope. The plan's list of cheap checks loses one
    item; the other five cover the intended regime.

- **`WarningConfig` schema.**

  ```python
  @dataclass
  class WarningConfig:
      # Thresholds — split by check
      pH_change_threshold:       float = 0.3   # pH units per step
      newton_iters_threshold:    int   = 15    # iters per speciation solve
      charge_residual_threshold: float = 1e-8  # |sum z_i C_i|
      dt_over_tau_min_threshold: float = 0.1   # dt / tau_min
      scipy_rejection_threshold: float = 0.3   # rejected / accepted

      # Ionic-strength regime (migrated from engine.py:281)
      ionic_strength_ideal_threshold:  float = 0.10  # mol/L
      ionic_strength_davies_threshold: float = 0.50  # mol/L

      # Throttle
      throttle: Literal["once", "first_N", "always", "silent"] = "once"
      first_N:  int = 10

      @classmethod
      def silent(cls)     -> "WarningConfig": return cls(throttle="silent")
      @classmethod
      def verbose(cls)    -> "WarningConfig": return cls(throttle="always")
      @classmethod
      def production(cls) -> "WarningConfig": return cls(throttle="first_N", first_N=3)
  ```

  Each threshold corresponds to one named check in the monitor.
  Adding a future check is one new field on `WarningConfig` and one
  new method on `AccuracyMonitor` — additive, no schema migration.

## Checkpoints

The checkpoints are ordered so each leaves the test suite in a
runnable state. The config + warning category land first (1–2) so
later checkpoints can import them. The monitor itself (3) lands next.
Hook integrations (4–6) follow. The ionic-strength migration (7) is
isolated and self-contained. Tests and ship (8–10) close out.

**Commit discipline:** commit immediately after each checkpoint
finishes (suite green for that checkpoint's slice). Mid-flight edits
were silently lost between checkpoints during `chemistry-unification-3`;
per-checkpoint commits are load-bearing, not optional.

### 1. `WarningConfig`, `Config`, `VLSIM_WARNINGS` env-var

- [ ] New module
      [`src/config.py`](../../src/config.py) implementing the
      schema from *Resolved decisions* above. Module-level
      `config = Config()` singleton. `os.environ` check at import
      time selects `WarningConfig.silent()` / `.verbose()` /
      `.production()` if `VLSIM_WARNINGS` is set to one of those
      values (case-insensitive). Unknown values: log a single
      `UserWarning` and fall back to defaults.
- [ ] No imports from `src/monitoring/` here — this module is the
      foundation, monitoring depends on it.
- [ ] Edge cases: empty `VLSIM_WARNINGS=""` falls through silently
      (treated as unset); whitespace stripped before preset lookup.

Sanity check: `python -c "import VLsim; print(VLsim.config.warnings.throttle)"` prints `once`.
`VLSIM_WARNINGS=silent python -c "import VLsim; print(VLsim.config.warnings.throttle)"` prints `silent`.
`VLSIM_WARNINGS=bogus python -c "import VLsim"` warns once and continues.

### 2. `AccuracyWarning` category + skeleton `AccuracyMonitor`

- [ ] New module
      [`src/monitoring/accuracy.py`](../../src/monitoring/accuracy.py)
      with `AccuracyWarning(UserWarning)` and an
      `AccuracyMonitor` class skeleton:

      ```python
      class AccuracyMonitor:
          def __init__(self) -> None:
              self.last_pH: Optional[float] = None
              self.step_count: int = 0
              self._warned_categories: Set[str] = set()

          def _emit(self, category: str, message: str) -> None:
              """Apply throttle + dispatch warnings.warn."""
              cfg = VLsim.config.warnings
              if cfg.throttle == "silent":
                  return
              if cfg.throttle == "once" and category in self._warned_categories:
                  return
              if cfg.throttle == "first_N":
                  count = _summary_counter[category]
                  if count >= cfg.first_N:
                      return
              self._warned_categories.add(category)
              _summary_counter[category] += 1
              warnings.warn(message, AccuracyWarning, stacklevel=3)
      ```

      Plus a module-level `_summary_counter: Counter[str]`.
- [ ] New helper `print_accuracy_summary()` that emits the
      end-of-simulation `"AccuracyWarning emitted N times total
      across [categories...]"` line when called. Wired in for
      the `"once"` throttle mode; left as a public helper for
      users who want to inspect / clear the counter mid-run.
- [ ] `src/monitoring/__init__.py` exporting
      `AccuracyMonitor` and `AccuracyWarning`.
- [ ] [`src/__init__.py`](../../src/__init__.py) — add exports:
      `AccuracyWarning`, `WarningConfig`, `config`. Keep the
      existing `__all__` ordering; add the new entries at the
      end of the relevant group.

Sanity check: `from VLsim import AccuracyWarning, WarningConfig, config` succeeds.
A skeleton test can instantiate `AccuracyMonitor()`, call a no-op `_emit("test", "test")`, and verify
`warnings.warn` is called with `AccuracyWarning`. `_summary_counter["test"]` increments to 1; a second
`_emit("test", "test")` is silently suppressed under throttle="once".

### 3. `AccuracyMonitor` per-check methods

Implement the per-check methods on `AccuracyMonitor`. Each one is a
small function: read inputs, compute the metric, compare against
the configured threshold, call `_emit(category, message)` on
crossing. All thresholds are read from
`VLsim.config.warnings` at call time (not cached on the monitor),
so post-construction config tweaks propagate.

- [ ] `check_pH_jump(self, current_pH: float) -> None`: emits
      `pH_change` category when `abs(current_pH - last_pH) >
      pH_change_threshold`. Sets `last_pH = current_pH` on the
      way out. No-ops on first call (when `last_pH is None`).
      Message includes the magnitude and threshold for diagnostic
      value.
- [ ] `check_newton_iters(self, iters: int) -> None`: emits
      `newton_iters` category when `iters >
      newton_iters_threshold`. No-op if `iters` is `None` or
      negative (indicates the underlying solver didn't expose
      the count — best-effort by design).
- [ ] `check_charge_residual(self, residual: float) -> None`:
      emits `charge_residual` category when `abs(residual) >
      charge_residual_threshold`. No-op on `None` or `NaN`.
- [ ] `check_dt_vs_tau_min(self, dt_h: float, tau_min_h: float)
      -> None`: emits `dt_over_tau_min` when `dt_h / tau_min_h >
      dt_over_tau_min_threshold`. Called once at CV construction
      typically; `tau_min_h=None` (no Jacobian estimate available)
      is a no-op.
- [ ] `check_scipy_rejections(self, n_accepted: int, n_attempted:
      int) -> None`: emits `scipy_rejections` when the rejection
      fraction `(n_attempted - n_accepted) / max(n_accepted, 1)`
      exceeds `scipy_rejection_threshold`. No-op if `n_attempted
      <= 1` (not enough data).
- [ ] `check_ionic_strength(self, I: float, activity_model: str,
      use_activity: bool) -> None`: emits `ionic_strength_high`
      (with sub-category in the message: `ideal` or `davies`)
      when I exceeds the regime threshold. Replaces the
      hand-rolled per-instance dedup in
      [`engine.py:281`](../../src/speciation/engine.py#L281) —
      throttle is the new dedup. Message format mirrors the
      existing wording so log parsers don't regress.
- [ ] Reset method: `reset(self) -> None` clears
      `_warned_categories` and `last_pH`, leaving `step_count`
      tracked. Useful for solver tests and re-runs in the same
      process; not currently called by production code paths.

Sanity check: each method tested in isolation — set
`VLsim.config.warnings.pH_change_threshold = 0.1`,
call `monitor.check_pH_jump(7.0)` (no warn — first call),
call `monitor.check_pH_jump(7.3)` (warns: |delta|=0.3 > 0.1),
call `monitor.check_pH_jump(7.31)` (no warn — once throttle held).
Switch throttle to `"always"`, call `monitor.check_pH_jump(7.6)`
(warns again).

### 4. Construct monitor in `ControlVolume.__init__`

- [ ] [`src/core/control_volume.py`](../../src/core/control_volume.py)
      — in `__init__`, after `reaction_model` is attached and
      `derive_speciation_keys` has run on internal links, construct
      `self._accuracy_monitor = AccuracyMonitor()`. Import is local
      to the constructor to avoid module-load ordering issues.
- [ ] Add a `dt` vs τ_min hook: try to compute a Jacobian-free
      `tau_min_h` estimate from the kinetic reactions in the
      attached `reaction_model`. If the reaction set's kinetic
      bucket is empty or the rate estimate is unavailable,
      `tau_min_h=None` and the check no-ops. Wire the check at
      `__init__` time only — the answer doesn't change per step.
- [ ] Keep the existing `_warnings.warn(...)` species-consistency
      call at line 120 as-is — it's a `UserWarning` for setup
      issues, not an accuracy concern.

Sanity check: constructing a CV with no kinetic reactions does
not emit any `AccuracyWarning`. Constructing a CV with a fast
kinetic reaction and a large user-supplied `dt_h` (passed via
some downstream test runner) emits one
`AccuracyWarning(category="dt_over_tau_min")`.

### 5. Wire pH-jump + scipy-rejection hooks in solvers

- [ ] [`src/core/solvers.py`](../../src/core/solvers.py)
      `EulerSnapshotSolver.solve_step`: after the step-closing
      speciation re-solve writes the final pH back, call
      `cv._accuracy_monitor.check_pH_jump(final_pH)`. Wrap in a
      `hasattr(cv, '_accuracy_monitor')` guard so the test
      fixtures that bypass `ControlVolume.__init__` (if any)
      don't crash.
- [ ] [`src/core/solvers.py`](../../src/core/solvers.py)
      `ScipyODESolver.solve_step` (lines 605–637): after
      `solve_ivp` returns, compute
      `(n_attempted, n_accepted) = (result.nfev, len(result.t))`
      and call
      `cv._accuracy_monitor.check_scipy_rejections(...)`. The
      `nfev` field is the function-evaluation count, which is a
      proxy for attempted steps (each accepted step ≈ N fevals
      depending on the integrator; the ratio is meaningful
      mainly when it spikes). Wire the pH-jump check at the
      same site that writes `phase.n_mol["H+"]` back from the
      final speciation.
- [ ] Don't migrate the `freeze_speciation=True` warning at
      line 401 to `AccuracyWarning` — it's a one-shot
      architectural notice fired at construction, not a
      per-step accuracy check. Stays `UserWarning`.
- [ ] Don't migrate the `result.success=False` warning at
      line 619 to `AccuracyWarning` — it's a solver failure
      (`solve_ivp` returned an error), not a regime warning.
      Stays `RuntimeWarning`.

Sanity check: a 10-step BSM2-style test run with default
`AccuracyWarning` filters silent (no thresholds crossed) under
the default `throttle="once"`. Manually setting
`VLsim.config.warnings.pH_change_threshold = 0.0` and running
the same test causes exactly one `AccuracyWarning` (subsequent
steps suppressed by `"once"`).

### 6. Wire Newton-iters + charge-residual hooks in engine

- [ ] [`src/speciation/engine.py`](../../src/speciation/engine.py)
      `solve()`: after the model returns `out`, look for fields
      `out.get("iters")` and `out.get("charge_residual")` and
      pass them to `cv._accuracy_monitor.check_newton_iters(...)`
      and `check_charge_residual(...)` respectively. The engine
      doesn't hold a reference to its parent CV today; either:
      (a) thread the monitor through `solve()` as a keyword
      argument (`monitor: Optional[AccuracyMonitor] = None`)
      with a no-op default, or
      (b) accept the monitor reference at
      `SpeciationPropertySolver.__init__` time (the solver does
      receive a `cv` handle in `apply_to_cv` already).
      Decide at implementation time; option (b) avoids
      threading and is probably cleaner. Mark with a TODO if
      either solver path can be reached without a monitor.
- [ ] If the underlying solver (`brentq` / `fsolve` / model-
      specific Newton) doesn't surface `iters` or
      `charge_residual` in its output dict, the corresponding
      check is silently a no-op (`None` -> early return on the
      monitor side). This phase does not add machinery to extract
      iteration counts from scipy primitives; that's a separate
      change.

Sanity check: tests that construct a `SpeciationEngine` and
solve a known-difficult charge-balance problem (high I, sharp
pH boundary) confirm that an `AccuracyWarning(category="newton_iters")`
fires when the underlying solver exposes its iteration count.
For solvers that don't expose it, the test asserts no warning is
emitted (silent best-effort, not a regression).

### 7. Migrate ionic-strength threshold warning to `AccuracyWarning`

- [ ] [`src/speciation/engine.py:281-372`](../../src/speciation/engine.py#L281)
      — replace the existing block (lines 321–372, with the
      `_warned_high_I` per-instance cache and the
      hand-rolled threshold dispatch) with a single call:

      ```python
      I = out.get("IonicStrength", None)
      if I is not None and (I_val := _safe_float(I)) is not None:
          cv._accuracy_monitor.check_ionic_strength(
              I_val,
              activity_model=str(self.activity_model),
              use_activity=bool(self.use_activity),
          )
      ```

      where `_safe_float` returns `None` for `NaN` / non-numeric.
- [ ] Delete the `self._warned_high_I: Set[str]` attribute and
      its initialisation in `__init__`. Throttle handles dedup
      now.
- [ ] Move the threshold table (`{"ideal": 0.10, "davies": 0.50}`)
      onto `WarningConfig` (already specified above as
      `ionic_strength_ideal_threshold` and
      `ionic_strength_davies_threshold`). The
      `check_ionic_strength` method on the monitor reads these.
- [ ] Update the message text in `check_ionic_strength` to
      mirror the legacy wording sufficiently — preserve the
      `"Ionic strength X mol/L exceeds..."` phrasing so log
      parsers (if any) don't regress.

Sanity check: the existing
[`tests/standalone/test_speciation.py`](../../tests/standalone/test_speciation.py)
high-I tests (if they exist; otherwise add one in checkpoint 8) that previously
asserted on `RuntimeWarning` now assert on `AccuracyWarning`. Behaviour
is otherwise unchanged: same threshold, same message format, same
dedup (now via throttle).

### 8. Tests — new

New tests under
[`tests/standalone/test_accuracy_monitor.py`](../../tests/standalone/test_accuracy_monitor.py)
(~20 tests):

- [ ] **Config + env-var (~4 tests):**
      - `VLsim.config.warnings.pH_change_threshold` defaults to
        `0.3`.
      - `WarningConfig.silent()`, `.verbose()`, `.production()`
        return configs with the expected `throttle` values.
      - Setting `os.environ["VLSIM_WARNINGS"] = "silent"`,
        reloading `VLsim.config` (or constructing a fresh
        `Config` via the env-var path), yields
        `throttle="silent"`.
      - Unknown env value emits one `UserWarning` and falls
        back to defaults.
- [ ] **Monitor per-check (~8 tests):**
      - `check_pH_jump` no-ops on first call; warns once when
        threshold crossed.
      - `check_pH_jump` resets `last_pH` correctly across calls.
      - `check_newton_iters` warns when threshold crossed;
        no-ops on `None`.
      - `check_charge_residual` warns on `abs(r) >
        threshold`; no-ops on `NaN`.
      - `check_dt_vs_tau_min` warns when ratio crossed;
        no-ops on `tau_min=None`.
      - `check_scipy_rejections` warns on ratio crossed;
        no-ops on `n_attempted <= 1`.
      - `check_ionic_strength` warns at I > 0.5 for davies;
        warns at I > 0.10 for ideal/use_activity=False; no-ops
        for unknown activity models.
      - Each `_emit` increments `_summary_counter` for the
        emitted category.
- [ ] **Throttle modes (~4 tests):**
      - `"once"` emits exactly one warning per category per
        monitor.
      - `"first_N"` emits N then stops (with `_summary_counter`
        as the global gate, so multi-CV runs collectively cap
        at N).
      - `"always"` emits unconditionally.
      - `"silent"` suppresses entirely.
- [ ] **Integration (~4 tests):**
      - Constructing a `ControlVolume` attaches a non-`None`
        `_accuracy_monitor`.
      - Running an `EulerSnapshotSolver` step with a deliberately
        large pH swing emits a `AccuracyWarning(pH_change)`.
      - Running a `SpeciationEngine.solve()` at high I emits a
        `AccuracyWarning(ionic_strength_high)`.
      - `print_accuracy_summary()` outputs the expected
        "N emissions across [categories]" line after a multi-CV
        run.

### 9. Tests — rewrites

- [ ] [`tests/standalone/test_speciation.py`](../../tests/standalone/test_speciation.py)
      — locate any test asserting on the legacy ionic-strength
      `RuntimeWarning` (search for `"Ionic strength"` /
      `RuntimeWarning` adjacency). Rewrite to assert on
      `AccuracyWarning` instead. If no such test exists today,
      this checkpoint becomes a no-op — the integration test in
      checkpoint 8 covers it.
- [ ] Search for any test using `pytest.warns(RuntimeWarning)`
      that incidentally catches the migrated ionic-strength
      warning. If found, update the expected category.
- [ ] BSM2 reference test
      ([`tests/standalone/test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py))
      should pass unchanged — Phase 4 doesn't touch numerics.
      Verify; if any sentinel drifts, investigate before
      re-baselining (Phase 4 has no mechanism for numerical
      drift).

Sanity check: `grep -rn 'Ionic strength.*RuntimeWarning' tests/`
returns empty after this checkpoint.

### 10. Documentation + ship

- [ ] [`docs/solvers.md`](../../docs/solvers.md) — add a brief
      "Accuracy monitoring" section: what the cheap checks are,
      how to configure (Python attribute / env-var), how to
      silence per-category via `warnings.simplefilter`. Forward
      note for `cv.run_with_diagnostics(...)` — list it as a
      "planned future entry point for opt-in expensive
      step-halving and cross-solver comparison".
- [ ] [`docs/upcoming/CHEMISTRY_UNIFICATION.md`](../upcoming/CHEMISTRY_UNIFICATION.md)
      — Phase F section: update to "Resolved 2026-05-XX in
      `chemistry-unification-4`". Note what deferred (the
      diagnostic mode, the splitting-error proxy).
- [ ] [`docs/upcoming/CHEMISTRY_UNIFICATION_PLAN.md`](../upcoming/CHEMISTRY_UNIFICATION_PLAN.md)
      — Phase 4 section: replace with a "Shipped" note pointing
      at the moved checklist. Update the "when Phase 4 ships"
      doc-update bullet — move both planning docs to
      `shipped/` only if Phase 3b is no longer pending
      (it is trigger-gated, so likely keep them in
      `upcoming/` with the chemistry-unification-1/2/3/4
      shipping banners).
- [ ] [`docs/upcoming/README.md`](../upcoming/README.md) — mark
      `chemistry-unification-4` shipped in the priority list.
      The chemistry-unification track is now complete except for
      3b (trigger-gated).
- [ ] [`README.md`](../../README.md) (if it references the
      phase landscape) — update.
- [ ] Move this checklist to
      [`../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md`](../shipped/)
      with a "Shipped" status banner mirroring Phase 3's.
- [ ] Final test sweep:
      `python -m pytest tests/standalone -q`. Baseline post
      Phase 3: 778 tests. Phase 4 estimate: 778 → ~798 (+20
      from `test_accuracy_monitor.py`; 0 net from rewrites).
- [ ] Ship via the convention in [../upcoming/README.md](../upcoming/README.md):
      `git checkout main && git merge --no-ff chemistry-unification-4
      -m "Merge chemistry-unification-4: AccuracyMonitor + AccuracyWarning + WarningConfig"`.
- [ ] Tag: `git tag chemistry-unification-4-shipped <commit-hash>`.
- [ ] Push: `git push && git push --tags`.
- [ ] Delete branch:
      `git branch -d chemistry-unification-4` and
      `git push origin --delete chemistry-unification-4`.

## Final test count expectation

778 → ~798 (estimate: +20 in `test_accuracy_monitor.py`;
±0 from rewrites — Phase 4 is purely additive on the test
side: no existing test asserts on the legacy
`_warned_high_I`-style behaviour because it lived inside the
engine and was dedup'd to one warning per process, which
typical tests don't observe). Final number confirmed
empirically at checkpoint 10.

## Risk notes

- **Monitor → engine threading.** The cleanest wiring (option (b)
  in checkpoint 6) is for `SpeciationPropertySolver` to capture
  the CV's monitor reference at `apply_to_cv` time. Verify that
  `apply_to_cv` is the lifecycle hook today; if it's renamed or
  removed in some Phase 1/2/3 fallout that I haven't traced, the
  checkpoint flags this and proposes (a) — kwarg threading —
  as the fallback.
- **τ_min from kinetic reactions.** A coarse `1 / max_rate`
  estimate at construction is fine for a regime warning, but if
  the kinetic rate envelope can't be cheaply queried from
  `ReactionSet.partition()`'s kinetic bucket, this check becomes
  a no-op for the first version. Document as a known gap;
  Phase 5 (if it exists) can improve the estimate.
- **Throttle interaction with `pytest.warns(...)`.** Python's
  warnings machinery resets filters per test via pytest's
  `_pytest.recwarn` plumbing, but the monitor's
  `_warned_categories` and `_summary_counter` do NOT reset
  automatically. Tests that emit the same `AccuracyWarning`
  multiple times in a session must either use the monitor's
  `reset()` method or instantiate a fresh `AccuracyMonitor` per
  test. Document on the public methods.
- **Cache lifecycle.** The module-level `_summary_counter`
  persists across test boundaries in the same pytest process.
  Add a `conftest.py` fixture or use `@pytest.fixture(autouse=True)`
  to clear it between tests, OR live with cross-test counting
  (acceptable if no test asserts on absolute counts beyond
  "non-zero"). Decide at checkpoint 8 when the tests are
  written.
