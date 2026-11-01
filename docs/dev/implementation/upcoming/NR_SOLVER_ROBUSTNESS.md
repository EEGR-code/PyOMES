# Newton solver robustness: direction and open issues

Status on 2026-11-01: paused by the owner in favour of a run-speed review. Nothing
here is in progress. Measurements are from one set of runs on `AD_BASIC`,
carbonate-ammonia, calcite and the iron oxidation test, so other chemistry is
untested.

**The goal:** a Newton solve either returns a correct answer or says plainly that it
did not, and a failed solve does not quietly distort a simulation. Detection is
done. Correction and the simulation's response are not.

## What already exists

| Piece | Where | Behaviour |
|---|---|---|
| Limit check | `NRChemicalEquilibriumEngine.__init__` | `ValueError` unless `0 < max_log_activity <= 150` |
| Limit-effects report | `_limit_effects`, `extra` | `clamped_log_activities` and `absent_component_species`, with one `RuntimeWarning` per species per engine |
| Step cap | `_capped_step`, `max_newton_step=4.0` on `solve_nr` | Scales the whole Newton step; `float("inf")` turns it off |
| Failure flag | `EquilibriumResult.converged`, `max_residual` | Final residual across every balance below 1e-7 mol/L; other engines leave them `None` |
| Failure warning | `_warn_not_converged` | One `RuntimeWarning` per engine |
| Failed states in the cache | `_is_stranded` | Kept unless a value is non-finite or a species is above the upper limit |

Tests are in `tests/standalone/test_nr_log_activity_clamp.py`.

## Open issues, in suggested order

### 1. Failed warm starts return wrong answers (do this first)

With ideal activity and warm starts across large jumps, 10 of 600 solves failed,
with pH wrong by up to 3.4. Newton is still progressing at about 0.43 log units per
step when it reaches `max_inner` (50), and ideal activity has no outer loop to give
it another pass. Davies had none of these failures. Small time-step changes did not
fail. These answers now return `converged=False` but are still wrong.

**Direction:** when a warm-started solve fails, retry once from a cold start inside
`engine.solve()`, and report and warn only if the retry fails too. Solving again from
the failed state converged in all 10, and the failures inspected also converged from
a cold start, so compare three options: a cold retry, continuing from the failed
state, and more Newton steps for ideal activity. The precipitation loop has no cache,
so the retry does not apply there.

**Check:** failures on several seeds of warm jumps, extra cost on the failing solves,
no change at all on solves that converge, and the forced iron runs below.

### 2. Nothing acts on the flag

All 7 simulation call sites in `core/` call `apply_to_phases` without checking
`converged`. The single warning per engine is the only signal during a run, and it
can come from a rejected ODE trial state.

- **Report through `AccuracyMonitor`** with a count per run. `check_charge_residual`
  and `check_newton_iters` already exist and have no callers, and
  `EquilibriumResult.n_iter` is never filled in.
- **The step-closing solve** at `core/solvers.py` writes to the live state. See the
  owner decision below.
- **ODE right-hand side:** treat a failed solve as a failed step, so the existing
  restore-to-pre-step path runs. How `solve_ivp` reacts is untested.

### 3. A failed Davies solve can still raise a bare `OverflowError`

It comes from `DaviesLiquidModel.gamma` with no chemistry context, before any result
exists. It was seen with the cap on only when Newton steps were forced low. From
reading the code, a large residual state makes the ionic strength enormous; this was
not traced step by step.

### 4. Absent components' species are not held to their zero total

The pin replaces the component's mass balance, so nothing limits its secondaries, and
they still enter the charge balance and ionic strength. At a limit of 20 in 1 mol/L
Na+, `PO4---` reaches 0.108 mol/L from a zero `H3PO4` total and pH falls by 0.167. At
the default limit, the largest amount measured was 2.8e-31 mol/L, but a species with
more protons removed would leave less margin. It is now reported but not corrected.
Dropping these species from the charge balance and ionic strength would fix the
cause, and would change results for any solve that has them.

### 5. Very small totals are not resolved

The tolerance is absolute, so a carbon total of 1e-19 mol/L was reported as about
1e-14, and 1e-12 was correct. Totals near 1e-6 mol/L are resolved to about 1 part in
10,000. pH was unaffected. A tolerance relative to each row's total would fix it, but
would change results and may cost iterations.

### 6. Smaller items

- The line search takes its last trial step even when no halving helped.
- `converged` describes the final inner solve, not whether the mineral loop finished.
- `retain_jacobian=True` keeps a failed state's Jacobian, so `jacobian_dz_dy()` would
  differentiate from it.
- `max_log_activity` still sets both the clamp and the pin, and cannot be set through
  `ReactionSystem`.
- `equilibrate_to_pH` raises `UnicodeEncodeError` on its `→` message when output goes
  to a cp1252 pipe or file.

## Owner decisions needed

1. **When a step-closing solve fails, what does a simulation do?** Warn and keep the
   last good speciation (the first recommendation), stop with an error, or a setting.
2. **Should enforced checks become optional?** Not decided.
3. **Should the ceiling of 150 rise?** With the cap on, a limit of 300 raised no
   overflow warnings in 200 cold solves. Without the cap, overflow warnings start at
   about 154.
4. **Should `reset_cache()` re-enable warnings?** It does not today.

## Lessons to keep

- **Do not discard failed states from the warm-start cache.** Keeping the last good
  state left repeated solves stuck (converged 0 of 17), and clearing the cache every
  time left a forced iron run with pH wrong by 4.8. Keeping the failed state gives
  2e-14. Any simulation response to a failed solve should be tested against this.
- **Judge failure on the final residual, not the inner loop.** Some solves exhaust one
  activity iteration and still finish correct.
- **Check that a test can fail.** Deliberate breakage found several empty tests:
  `pytest.approx(1e-50, rel=...)` accepts 0.0 unless `abs=0` is passed, and inputs can
  stop reaching the clamp once the cap is on. Tests that need a condition assert that
  it happened.
- **Scale the whole Newton step, not each component.** Limiting each component left
  up to 186 of 300 solves unconverged.
- **Separate the workloads.** A sweep of only Davies warm jumps hid the ideal
  failures, and a claim of "no failures" did not survive the next workload.
- **The threshold is close at high strength.** Good Davies solves reached 5.6e-8 at
  several mol/L against the 1e-7 limit.

## Reproducing the measurements

Throwaway scripts are in the gitignored `/scratch/nr_solver_experiments/`. They exist
only on the owner's machine.

| Script | Purpose |
|---|---|
| `fingerprint.py OUT.json` | Exact values for direct solves, calcite and the iron run; compare files byte for byte. Set `PYTHONIOENCODING=utf-8` |
| `step_cap_lab.py` | In-memory `solve_nr` with counters for Newton steps, trials and failed line searches |
| `residual_survey.py` | Final residuals on normal solves, and residual against true error on forced failures |
| `warm_failure.py` | The warm-start failures in issue 1 |
| `failure_visibility2.py` | Forced failures, and what a user sees |
| `cache_options2.py` | Cache rules compared, including the forced iron runs |
| `mut_plugin5.py` | Pytest plugin applying deliberate breakage in memory, chosen with `CLAMP_MUTATION` |

To force a failure, use `functools.partial(solve_nr, max_inner=k)` on the engine
module's `solve_nr`.

The history of how this was found is in the owner's `devlog.md`, outside the repo.
Detailed entries are also in [`../OPEN_WORK.md`](../OPEN_WORK.md).
