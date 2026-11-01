# Persistent integrator for `SimultaneousAdaptiveSolver`

Status: design agreed with the owner, 2026-11-01. Default-on decided. No code yet.
Checklist: [PERSISTENT_INTEGRATOR_CHECKLIST.md](PERSISTENT_INTEGRATOR_CHECKLIST.md).

## Problem

`SimultaneousAdaptiveSolver.solve_step` calls `solve_ivp` once per macro step. Each
call discards the step size, order and Jacobian, so `n_steps` sets a restart rate, not
just an output grid. A finer grid makes a run slower and less accurate.

Measured with an in-memory prototype, one BDF object across macro steps, no events:

| Model | Shipped | Persistent | Speedup | Error, shipped → persistent |
|---|---|---|---|---|
| Iron, BDF, ideal, n=200 | 3.3 s, 5,770 solves | 0.31 s, 363 solves | ×10 | 6.6e-5 → 3.2e-6 |
| Iron, BDF, Davies, n=200 | 6.8 s | 0.86 s | ×8 | 1.6e-5 → 3.2e-6 |
| BSM2, BDF, n=100 | 4.7 s, 5,768 solves | 0.61 s, 713 solves | ×8 | 5.5e-5 → 1.3e-6 |

Error is measured against a tight reference. A forced restart every 20 steps gave ×6,
every 5 steps ×3, every step ×1. No real controller, profile or feed has been tested.

## Scope

In scope: `SimultaneousAdaptiveSolver`, restart detection, `env.t_h` under this solver,
and a benchmark in `tests/performance/`.

Out of scope: `MonolithicODESolver`, linked CVs on split system solvers (links change
the state each macro step, so these restart every step and gain nothing), Newton
arrays, and the mineral bugs.

## Design

**Stepping.** Each CV keeps one scipy `OdeSolver` object. `solve_step` advances it until
it reaches or passes the macro-step end, then interpolates there. The object may step
past that end; those values are never written back. `n_steps` becomes an output grid
when nothing forces a restart.

**State is per CV.** One solver instance can serve several CVs, so stepper state is
keyed by CV.

**Totals, not species, are floored.** The state holds every speciated species, but the
engine reads only their sums and overwrites them. A flux in proportion to a speciated
amount (a drain, a rate law on a resolved species) can take one entry negative while
the total is fine. Flooring that entry creates material: without re-speciation at each
restart this reached an error of 0.89 on the feed-and-drain case. So a persistent run
does not floor liquid species in `engine.algebraic_species()`; everything else is
floored as before.

**Restart triggers.** The stepper is dropped and rebuilt from the current state when
any of these holds:

1. Any `n_mol` entry, `T_K` or `V_L` differs **exactly** from the snapshot taken at the
   end of the previous `solve_step`, after the closing solve. This catches controller
   doses, links and direct writes without listing them.
2. A `MutableScalar` or `MutableDict` write happened anywhere since then, detected by a
   generation counter in `control/descriptors.py`. State comparison cannot see these:
   profiles change parameters, never amounts.
3. `Simulation` reports a profile record, or a controller action with a parameter
   change, in the last step. This covers `_write_leaf`'s direct write into a plain dict
   (`link.kLa["O2"]`), which bypasses the descriptors. A controller that runs and does
   nothing must not force a restart; doses are caught by trigger 1.
4. `dt_h` or `external_source_terms` differ from the previous call, because `f` closes
   over them.
5. `freeze_speciation=True`, which keeps its per-macro-step meaning.
6. The caller calls the new public `solver.restart()`.
7. The species in a phase differ from those the state vector was built with, or the
   reaction system or engine has been replaced. A species outside the state vector
   would sit in the scratch phases in stale form and be counted in the totals.

**Known gap.** Changes the framework cannot see are not detected: a plain attribute on a
boundary, or state hidden in a `rate_fn` closure. They need `solver.restart()`, and
`docs/solvers.md` will say so. A sudden change inside a rate law is still handled by
the integrator's error control but may be stepped over if very short, unless `max_step`
is set.

**Failure.** If a continued stepper fails, retry once from a fresh stepper at the
macro-step start. If that fails too, keep today's behaviour: restore the state, emit
`RuntimeWarning`, drop the stepper.

**Checkpoints.** Pickling drops stepper state, because it holds closures. A resumed run
restarts once, so it matches a continuous run within tolerance, not bit for bit.

**Accuracy monitor.** `check_scipy_rejections` receives per-macro-step deltas, and is
skipped when the integrator took no step in that macro step.

**Opt-out.** `SimultaneousAdaptiveSolver(restart_each_step=True)` keeps today's
behaviour exactly, through `solve_ivp`.

## Decisions

- **Default-on** (owner, 2026-11-01). Adaptive results change, so earlier adaptive runs
  will not reproduce exactly. Stored adaptive numbers and notebook outputs need review.
- **`env.t_h` becomes continuous**: `t_h` passed to `solve_step` plus time since then.
  Today it runs over `[0, dt_h]` each macro step, and that cannot be kept once a step
  spans macro boundaries. This changes results for time-dependent rate laws under this
  solver, as a bug fix.
- **`Simulation` still passes the end-of-step time as `t_h`**, so `env.t_h` under this
  solver leads by one `dt_h`, as sequential and Euler already do at the start of a step.
  Correcting that affects every solver and is logged separately, not done here.
- **All `solve_ivp` methods are covered**, explicit ones included. The prototype tested
  BDF and Radau only; the benchmark must cover the others before the default flips.

## Acceptance

Before the default flips, each model in the benchmark must show persistent error no
worse than `restart_each_step=True` against a tight reference, and no more engine
solves. This includes the two untested cases: a real controller and a feed-and-drain
model. A model that fails this stops the work until it is understood.

**Open, for the owner:** two kinds of case fail this as written, and the gate needs a
decision before the default flips.

- *Explicit methods.* Restarting holds them well below the requested tolerance by
  forcing short steps: RK45 reaches 3e-8 at `rtol=1e-6` at about 3× the cost, where
  persistent reaches 1.2e-6. At `rtol=1e-7` persistent passes for RK45, RK23 and DOP853,
  with fewer solves than restarting (253 vs 973, 709 vs 997, 487 vs 1,501).
- *A controller dosing every step.* Both restart every step, at the same cost, and
  persistent's error is 16% larger (3.6e-5 vs 3.1e-5). Stepping past the macro-step end
  takes `Fe2+` negative in 97 of 6,020 trial states, and the floor acts on it. Either
  stopping at the macro-step end (3.13e-5) or removing the floor (3.16e-5) closes the
  gap. A run that has just restarted could stop at the macro-step end until the state
  stops changing.

Every restart trigger has a test that fails if the trigger is removed.

## What changes elsewhere

- **Tests:** three tests in `test_simultaneous_adaptive_solver_jac.py` patch
  `scipy.integrate.solve_ivp` to capture `jac`. They move to a small
  `_make_stepper` seam in `core/solvers.py`. Their assertions do not change.
- **Unaffected:** sequential and Euler solvers, and the BSM2 golden test.
- **Docs:** `docs/solvers.md`, `docs/architecture.md`, and the `Simulation.run`
  docstring, which already calls `n_steps` an output grid.
- **Notebooks to re-run by the owner:** `07_iron_oxidation`,
  `08_iron_oxidation_and_precipitation`, `10_engine_protocol_hierarchy`,
  `ArXiv_preprint/02_kinetic_co2_equilibration_microplate_well`, and
  `SolverProtocols/0_README`.
