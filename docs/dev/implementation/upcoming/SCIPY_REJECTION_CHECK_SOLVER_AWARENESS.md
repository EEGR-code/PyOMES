# scipy Step-Rejection Check: Solver-Family Awareness — Design Note

> Status: design note, not yet started. No branch, no checklist, no code
> yet. Written on `main` (independently of the `tutorials-followups` phase
> that surfaced the motivating gap) — see "How to start one" below when
> picked up.

## Motivation

While fixing `tests/validation/speciation/07_iron_oxidation.ipynb` /
`08_iron_oxidation_and_precipitation.ipynb` (tutorials-followups checkpoint
4), switching to `SimultaneousAdaptiveSolver(method='BDF')` to properly
resolve a stiff Fe²⁺-oxidation transient (t½ ≈ 11 min inside a 2 h run)
eliminated the `ConservationWarning`s the default explicit solver produced,
but introduced a persistent `AccuracyWarning` from
`AccuracyMonitor.check_scipy_rejections`
(`PyOMES/monitoring/accuracy.py:277-309`):

```
SimultaneousAdaptiveSolver attempted 148 function evaluations for 74
accepted steps (ratio 1.00 above threshold 0.30). Suggests stiffness or a
tight tolerance boundary; consider switching to Radau/BDF or relaxing
rtol/atol.
```

Investigation (varying `rtol`/`atol` across three orders of magnitude,
varying `max_step`, enabling `use_engine_jacobian=True`, and trying both
`BDF` and `Radau`) found the ratio stuck at ~1.0 regardless — and Radau's
ratio (6.54) was *higher* than BDF's, the opposite of what "switch to
Radau/BDF" as a remedy would predict if the problem were genuine step
rejection.

## Root cause

`check_scipy_rejections`' `n_attempted` parameter is populated from
`result.nfev` (total RHS evaluation count) at its call site
(`PyOMES/core/solvers.py:913-926`), used as a proxy for "attempted steps."
The docstring is candid that this is approximate ("a rough indicator"). For
an **explicit** method (RK45, DOP853 — this class's default), `nfev` tracks
accepted steps closely because each accepted step costs a fixed number of
stage evaluations, and a truly rejected step costs evaluations without
producing an accepted one — so a high `(nfev - accepted)/accepted` ratio is
a genuine signal that the step-size controller is struggling (stiffness, or
too-tight tolerances), and "switch to an implicit method" is the right fix.

For an **implicit** method (BDF, Radau), every *accepted* step runs an
internal Newton iteration, which costs several extra RHS evaluations by
design — normal, expected overhead, not evidence of rejection. That
structural overhead is what the ~1.0 ratio was actually measuring, which is
also why no tolerance/step-size adjustment budged it, and why Radau (whose
implicit Runge-Kutta stages need even more evaluations per step than BDF's
linear multistep approach) scored worse rather than better. The check has
no notion of solver family, so it applies the explicit-method assumption
(`nfev ≈ accepted`) uniformly, producing a structural false positive for
every implicit-method run of any genuinely fast/stiff reaction — exactly
the class of problem `SimultaneousAdaptiveSolver(method='BDF'/'Radau')`
exists to solve well.

**scipy itself is not at fault here** — `nfev` is exactly what scipy
documents it as, and BDF's internal Newton iteration is correct, intended
behavior for an implicit solver. The gap is in how `check_scipy_rejections`
interprets that number.

## Proposed fix (Option 1 of 3 considered — see Alternatives)

Make the check aware of which solver family produced the statistics it's
judging, and hold each family to a baseline appropriate to it, rather than
one flat ratio for every method.

- `SimultaneousAdaptiveSolver` already knows `self.method` at the
  `check_scipy_rejections` call site (`solvers.py:923-926`, inside the same
  method that already references `self.method` in its `RuntimeWarning`
  message a few lines above). Pass it through:
  `monitor.check_scipy_rejections(n_accepted=..., n_attempted=..., method=self.method)`.
- Classify `method` into `explicit = {"RK45", "RK23", "DOP853"}` /
  `implicit = {"Radau", "BDF"}` inside `check_scipy_rejections` (or a small
  module-level constant next to it), and compare against a per-family
  threshold instead of one flat value.
- `WarningConfig.scipy_rejection_threshold` (`PyOMES/config.py:45`,
  currently a single `float = 0.3`) becomes two fields (or a
  `Dict[str, float]` keyed by family) —
  `scipy_rejection_threshold_explicit` (keep the existing `0.3` default,
  unaffected) and `scipy_rejection_threshold_implicit` (a new, separately
  calibrated default — needs empirical data, see Open Questions).
- `LSODA` (available via `scipy.integrate.solve_ivp` but not currently
  documented as a `SimultaneousAdaptiveSolver` option) is a hybrid
  explicit/implicit method that switches internally based on detected
  stiffness — it doesn't cleanly belong in either bucket. Not a blocker
  for this fix (it isn't a supported `method` value today), but worth
  flagging so whoever adds `LSODA` support later doesn't have to
  rediscover this.

## Alternatives considered

1. **(Chosen)** Solver-family-aware threshold, above. Cheapest fix that
   directly resolves the concrete false positive; doesn't touch scipy
   internals.
2. **Data-driven per-step baseline.** Instead of a flat threshold, estimate
   an expected `nfev`-per-accepted-step baseline for the method in use
   (e.g., implicit methods' Newton overhead roughly scales with state-vector
   size when no analytical Jacobian is supplied — consistent with
   `use_engine_jacobian=True` measurably reducing, but not eliminating, the
   `nfev` count in the motivating investigation) and flag only when the
   observed ratio exceeds *that* baseline by a margin. More accurate, more
   design work, and needs a defensible formula for the baseline — deferred
   unless Option 1's flat implicit threshold proves too coarse in practice.
3. **True rejection count instead of the `nfev` proxy.** Vanilla
   `scipy.integrate.solve_ivp` doesn't expose actual rejected-step counts in
   its public `OdeResult` — getting a real count means subclassing/hooking
   the internal `OdeSolver` step loop. Most accurate, most invasive; not
   proposed here.

## Open questions (need a decision before implementing)

1. **What should the implicit-method default threshold be?** The
   motivating case saw ratio ≈ 1.0–1.04 across a wide sweep of
   tolerances/`max_step`/Jacobian settings for a genuinely well-behaved,
   physically-accurate BDF run (Fe mass balance drift 0.00005%, essentially
   solver-noise-free). A reasonable starting point is picking a default
   comfortably above the range seen in several genuinely-fine implicit
   runs (e.g. `1.5`–`2.0`) rather than guessing from one case — needs a
   small empirical sweep across a few existing implicit-solver test/demo
   cases in this repo (the ArXiv_preprint kinetic CO2 notebook already uses
   `SimultaneousAdaptiveSolver(method="BDF")` and would be a second data
   point) before picking a number.
2. **Should `use_engine_jacobian=True` shift the implicit baseline down?**
   It measurably reduced `nfev` in the motivating investigation (132→~96-102
   across sweeps) without eliminating the ratio — worth checking whether
   Jacobian-enabled implicit runs cluster at a lower, tighter baseline than
   Jacobian-disabled ones, which would argue for a third bucket rather than
   two.
3. Config shape: two named fields vs. a `Dict[str, float]` keyed by method
   family — the former is more discoverable in `WarningConfig`'s existing
   flat-fields style; the latter generalizes better if a third family
   (e.g. `LSODA`) needs its own bucket later.

## Implementation steps (once the questions above are resolved)

1. Resolve Open Questions 1–3.
2. Add the explicit/implicit method classification (module-level constant
   or small helper) next to `check_scipy_rejections` in
   `PyOMES/monitoring/accuracy.py`.
3. Extend `WarningConfig` (`PyOMES/config.py:45`) per the resolved shape;
   keep the existing `0.3` explicit default unchanged (no behavior change
   for any existing explicit-method caller).
4. Update `check_scipy_rejections`'s signature to accept `method: str`, and
   `solvers.py:923-926`'s call site to pass `self.method`.
5. Regression check: existing explicit-method (`DOP853`, the class default)
   test/demo runs must see zero behavior change — this is purely additive
   for explicit callers and only changes the threshold implicit callers are
   judged against.
6. Unit tests: an implicit-method run with a `nfev`/`accepted` ratio between
   the old flat `0.3` and the new implicit threshold no longer warns; an
   explicit-method run at the same ratio still does (regression guard that
   the fix didn't accidentally widen the explicit-method check too).
7. Revisit `07_iron_oxidation.ipynb`/`08_iron_oxidation_and_precipitation.ipynb`
   (tutorials-followups checkpoint 4) once shipped — the notebooks'
   documented-and-suppressed `AccuracyWarning` workaround can likely be
   removed at that point, the cleanest end-to-end confirmation the fix
   closes the gap it was built for.

## How to start one

Per this folder's usual convention
([README.md](README.md#how-to-start-one)): write a checklist file
(`SCIPY_REJECTION_CHECK_SOLVER_AWARENESS_CHECKLIST.md`), cut a branch off
`main` (suggested name: `scipy-rejection-check-solver-awareness`), work the
checkpoints there, and add a "Currently in flight" pointer to this folder's
`README.md`. Small and self-contained — no dependency on any other in-flight
phase. Natural to bundle with
[REACTION_ENVIRONMENT_PHASE_EXPOSURE.md](REACTION_ENVIRONMENT_PHASE_EXPOSURE.md)
and
[PHCONTROLLER_CORRECTOR_VALIDATION.md](PHCONTROLLER_CORRECTOR_VALIDATION.md)
if a single "monitoring/diagnostics polish" pass is ever taken, since all
three surfaced from the same tutorials-followups checkpoint 4 investigation,
but none depends on the others.
