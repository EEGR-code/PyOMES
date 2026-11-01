# Phase Checklist — persistent-integrator

> **Status: Shipped 2026-11-01** — merged into `main` by pull request, tagged
> `persistent-integrator-shipped`. Full suite on the final tree: 2292 passed, 87 skipped,
> 2 xfailed.

Design: [PERSISTENT_INTEGRATOR.md](PERSISTENT_INTEGRATOR.md).

## Pre-flight

- [x] `git status -sb` clean, `main` level with `origin/main` at `c0613d0`
- [x] Baseline suite at `c0613d0`, as reported by the owner: 2226 passed, 87 skipped,
      2 xfailed
- [x] Branch created off current `main`: `persistent-integrator`
- [x] This checklist and the design note committed as the first commit

## Checkpoints

One commit per checkpoint. The owner runs the full suite after each checkpoint that
changes `PyOMES/` or `tests/`.

- [x] 1. **Benchmark** in `tests/performance/`: iron (ideal, Davies), BSM2, a controller
      model, a feed-and-drain model. Records counters, CPU time and error against a
      stored tight reference. Baselines recorded before any package change. No change
      to `PyOMES/`.
      Recorded on Python 3.14.3, numpy 2.4.4, scipy 1.17.1. Iron BDF: 5,770 solves,
      3.19 s, worst error 6.65e-5. BSM2 BDF: 5,768 solves, 4.59 s, 1.93e-5.
- [x] 2. **Refactor only:** build `f` and `jac` outside the integrator call. Results
      must be **bit-identical** on every benchmark model.
- [x] 3. **Persistent stepping behind an option, default off:** `_make_stepper`, per-CV stepper,
      interpolation at macro-step ends, state-snapshot restart (trigger 1), `dt_h` and
      `external_source_terms` (trigger 4), `freeze_speciation` (trigger 5),
      `solver.restart()` (trigger 6), species and engine changes (trigger 7), totals
      not floored. With the option off, results stay bit-identical. Benchmark reference
      rebuilt on persistent Radau, because restarting does not conserve totals as well.
- [x] 4. **Parameter triggers:** the keyword renamed to `persistent`; generation counter in `control/descriptors.py`
      (trigger 2) and the `Simulation` signal (trigger 3). A test for each trigger that
      fails when it is removed, including `link.kLa["O2"]`.
- [x] 5. **Robustness:** after a restart, stop at the macro-step end until the state stops changing; retry on failure, pickling, one solver across several CVs,
      monitor deltas. Check `tests/standalone/test_package_layering.py` still passes.
- [x] 6. **`env.t_h` continuous under this solver**, with a test. Log `Simulation`'s
      end-of-step `t_h` in `OPEN_WORK.md`.
- [x] 7. **Acceptance run:** every benchmark model meets the gate in the design note,
      all `solve_ivp` methods included. Stop here if one does not.
- [x] 8. **Flip the default** to `persistent=True`. Move the
      three `solve_ivp`-patching tests to the `_make_stepper` seam. Update docs, and state
      plainly in `docs/solvers.md` and the solver docstring that explicit methods,
      including the default `DOP853`, become faster and less accurate at the same `rtol`,
      with both remedies: a tighter `rtol`, or `persistent=False`.
- [x] 9. **Owner re-runs the 4 notebooks that build this solver**; review any changed outputs.

## Implementation log

Results are from Python 3.14.3, numpy 2.4.4, scipy 1.17.1. A checkpoint is ticked once
the owner's full suite run has passed; each entry gives that result.

### Checkpoint 1: benchmark (`5d45471`)

- Added `tests/performance/bench_adaptive_solver.py`: 13 cases, gating on engine solves
  and error, or on a bit-level fingerprint with `--exact`.
- Baseline, restarting: iron BDF 5,770 solves, 3.19 s, error 6.65e-5; BSM2 BDF 5,768
  solves, 4.59 s, 1.93e-5.
- Found: `iron_ph_sampled_bdf` never dosed inside its 0.05 deadband. The band was
  tightened, and the idle version kept as its own case. So an idle controller must not
  force a restart, and the design changed.
- Full suite: not needed, no package or collected test changed.

### Checkpoint 2: refactor (`1c16b9e`)

- `solve_step` gets its right-hand side and Jacobian from `_build_integrand`. Logic only
  moved, plus a `stacklevel` bump so a warning points at the same caller.
- `_make_stepper` moved to checkpoint 3, so nothing unused was added.
- All 13 cases bit-identical.
- Full suite: 2226 passed, 87 skipped, 2 xfailed.

### Checkpoint 3: persistent stepping, opt-in (`2b8a824`)

- Off: all 13 cases bit-identical to the committed baseline.
- On: iron BDF ×12.9, error 6.6e-5 → 3.2e-6; BSM2 ×7.9, 1.9e-5 → 1.3e-6; iron Radau
  ×17.6; iron LSODA ×35; feed and drain ×12.9, 1.2e-4 → 6.6e-6. A dosing controller and
  a temperature ramp restart every step, ×1.0.
- Found: the feed-and-drain error reached 0.89, because flooring individual speciated
  species creates material. Resolved species are now exempt from the floor under
  persistent. The iron-only prototype could not show this.
- Found: restarting drifts total Fe by 9.8e-5 on feed and drain, and 1.6e-6 on iron
  with the floor fixed, against about 1e-7 for persistent. The first Radau reference was
  restarted and carried that bias, so it was rebuilt on persistent Radau. Persistent
  Radau, BDF and LSODA agreed to within 3e-7 on feed and drain and weak acid; on plain
  iron only total Fe supports the choice, since restarted Radau differs by 1.7e-6.
- Found: the extra dosing-controller error (3.1e-5 → 3.6e-5) comes from overshoot past
  the macro-step end taking `Fe2+` negative in 97 of 6,020 evaluations, where the floor
  acts. Stopping at the end gives 3.13e-5, and removing the floor gives 3.16e-5.
- Found: explicit methods get less accurate at the same `rtol`, because restarting
  capped their steps (RK45 median step 0.02 h vs 0.14 h).
- Logged the clamp findings in `OPEN_WORK.md`.
- Full suite: 2256 passed, 87 skipped, 2 xfailed.

### Checkpoint 4: parameter triggers (`c5e3dc9`)

- Keyword renamed from `restart_each_step` to `persistent`, still default `False`.
- Generation counter in `control/descriptors.py`, counting only changes to an existing
  value, so scratch phases and unchanged rewrites do not restart. `Simulation` notes
  changed profile reports, plain-dict writes and the unchecked-setter fallback.
- Owner decisions: explicit methods pass within 10× `rtol`, documented plainly; the
  dosing-controller gap is fixed rather than excused.
- Removed each restart check in memory to test the tests: the species and
  reaction-system checks had no real protection, and now have their own tests. The
  `freeze_speciation` test was not put through this check.
- Default path: all 13 cases bit-identical. Persistent numbers unchanged from
  checkpoint 3.
- Full suite: 2270 passed, 87 skipped, 2 xfailed.

### Checkpoint 5: robustness (`f3b00a7`)

- A new run stops exactly at the macro-step end; after a macro step with no change, the
  stepper is replaced by one with no end point. Failure retries once from a new run.
- Benchmark gate updated to the owner's decision: explicit methods may reach 10× `rtol`.
- Persistent gate passes on all 13 cases. The dosing-controller and temperature-ramp
  cases are now bit-identical to restarting (3.13e-5, 7.13e-5).
- Cost of the extra restart: iron BDF 363 → 396 solves, ×12.5 → ×11.5, error 3.19e-6 →
  4.52e-6; iron Radau ×17.9 → ×13.8; iron LSODA ×35 → ×26. BSM2 stayed at 713 solves.
- Pickling and one solver across several CVs were already covered at checkpoint 3.
- Full suite: 2278 passed, 87 skipped, 2 xfailed.

### Checkpoint 6: continuous `env.t_h` (`e146295`)

- Under `persistent=True`, rate laws read the `t_h` passed in plus the time since, across
  macro steps. `persistent=False` is unchanged, including its reset each step.
- New restart trigger: a `t_h` that does not follow on from the last step, with a
  tolerance for rounding in `Simulation`'s `linspace` grid.
- `env.t_h` still leads true time by one step under `Simulation`. Logged in
  `OPEN_WORK.md`, with the differences between solvers.
- Found by the new tests: scipy evaluates the right-hand side while building a stepper,
  so the time origin has to be set before that, or the first evaluation reads t = 0.
- No benchmark case has a time-dependent rate law, so this rests on the new tests.
- Full suite: no result reported at this commit; the run at checkpoint 8 covers it.

### Checkpoint 7: acceptance (`bbd9dcd`)

- Added to the benchmark: `iron_pulse_bdf`, a controller with quiet spells between doses,
  and a `clamps` command. Added tests: a time-dependent rate law against its closed form
  on 5 methods, linked CVs against restarting, and resuming with a pickled solver.
- Persistent gate passes on all 14 cases; default path bit-identical on all 14.

| Case | Solves, restarting → persistent | Error, restarting → persistent |
|---|---|---|
| Iron BDF | 5,770 → 396 | 6.6e-5 → 4.5e-6 |
| Iron Davies | 4,968 → 537 | 1.6e-5 → 4.0e-6 |
| Iron Radau | 14,544 → 765 | 9.1e-7 → 3.0e-7 |
| Iron LSODA | 25,418 → 617 | 1.3e-5 → 2.2e-6 |
| Feed and drain | 5,950 → 434 | 1.2e-4 → 6.0e-6 |
| BSM2 BDF | 5,768 → 713 | 1.9e-5 → 1.4e-6 |
| Dose every 20 steps | 5,796 → 998 | 6.4e-5 → 9.4e-6 |
| Controller idle | 5,770 → 396 | 6.6e-5 → 4.5e-6 |
| Controller, sampled | 6,843 → 6,224 | 4.6e-5 → 3.2e-5 |
| Controller, every step | 5,880 → 5,880 | identical, 3.1e-5 |
| Temperature ramp | 12,186 → 12,186 | identical, 7.1e-5 |
| Weak acid RK45 | 973 → 213 | 2.6e-8 → 1.1e-6 |
| Weak acid DOP853 | 1,501 → 474 | 7.2e-8 → 1.5e-6 |
| Weak acid RK23 | 997 → 393 | 2.7e-6 → 7.7e-6 |

- BSM2 over 48 h, 480 steps, outside the benchmark: 25,894 → 2,101 solves, 20.5 s → 1.6 s,
  error 2.7e-6 → 1.6e-6 against persistent Radau at `rtol=1e-9`.
- Clamps measured on 10 of the 14 cases, both paths: the Newton log-activity clamp, the
  `apply_flux` floor, the `PHController` dose limit, the vent cap and BSM2's rate-law
  guard changed no value. The state floor acted on `Fe2+` in the iron Davies, sampled
  and pulse cases, and on `Fe(OH)2+` and `FeOH2+` only in restarted feed and drain. Two
  bounds changed values only before or at the start of a run: `dt_h` in
  `equilibrate_to_pH`, and the Newton initial guess.
- Found: `save_checkpoint` fails for any model with `Species`, on any solver. Logged in
  `OPEN_WORK.md`.
- Not covered: a real integrator failure; the stirred-tank template; a controller that
  changes a parameter (`DOAgitationController`) in a full run.
- Full suite: no result reported at this commit; the run at checkpoint 8 covers it.

### Checkpoint 8: default flipped (`c27920e`)

- `persistent=True` is the default. `persistent=False` stays as the way back.
- The three `solve_ivp`-patching tests capture through `_make_stepper`, with the same
  assertions.
- `docs/solvers.md` and the solver docstring state that explicit methods, including the
  default `DOP853`, become faster and less accurate at the same `rtol`, with both remedies.
- Benchmark: `persistent=False` is pinned as the yardstick, and persistent fingerprints are
  stored, so a later refactor can be checked bit for bit on the default path.
- `freeze_speciation` logged in `OPEN_WORK.md` at the owner's request.
- Full suite: 2292 passed, 87 skipped, 2 xfailed.

### Checkpoint 9: notebooks (`d08f1fc`, `f562270`)

- A scan of every notebook's code found 4 that build this solver: `07_iron_oxidation`,
  `08_iron_oxidation_and_precipitation`, `10_engine_protocol_hierarchy` and
  `ArXiv_preprint/02_kinetic_co2_equilibration_microplate_well`. `SolverProtocols/0_README`,
  listed in the design note, only mentions it.
- Ran their code cells in a scratch script, without changing the files, and compared
  printed output with what is saved.
- `07` and `08`: results unchanged, run time 14-15 s → 0.8 s, largest Fe drift 0.0003% →
  0.0000%. The saved `'KOH'` wording in their pH line was already out of date.
- `10`: cell 14 raised `KeyError: 'nfev'`, because it read counts by patching `solve_ivp`.
  Fixed by passing `persistent=False` there. The suite does not run notebooks, so it could
  not catch this. Its PHREEQC cells depend on `phreeqpython` being installed.
- `ArXiv 02`: the saved sweep table was already out of date (it ends at kLa 500, and the
  code sweeps to 1000), and cells 0 and 5 have no saved output. `t95` moves by up to 0.4%
  against `persistent=False`. Against tight runs, which agree across both paths and BDF and
  Radau, the error at kLa 96, 3.46 and 0.1 is 1.1e-4, 1.2e-4 and 1.7e-6 for persistent, and
  4.2e-3, 4.3e-3 and 9.3e-4 for restarting. Final CO2(aq) goes the other way: 13.7790
  against the tight 13.7788 at kLa 96 and 3.46, an error of 1e-5 where restarting matches.
  This model holds about 1e-6 mol of CO2 against `atol=1e-9`.
- Its comment calling the solver stateless was corrected, in the notebook and its generator.
- The owner regenerated all four with `nbconvert`: no errors, and printed results as above.
  All P1-P5 checks in `07` and `08` still pass.
- Plots reviewed by eye for `07`, `10` and `ArXiv 02`, and look right. `07` and `08` had no
  saved plots before, so there is nothing to compare them with, and each file grew from
  about 45 KB to about 300 KB. `08`'s plots were extracted but not viewed; its printed
  results match `07`'s.
- `10`'s PHREEQC cells now say "skipped", because `phreeqpython` is not installed. The old
  saved output was an error message.
- Notebook text that said these run through `solve_ivp` corrected, in the `ArXiv 02`
  generator too.
- Left as they are: the `~70` vs `~15` accepted BDF steps and the `~0.0005%` to
  `~0.00005%` drift in `07` and `08`, which describe earlier measurements on restarting
  runs; and "analytical Jacobian" in `07`, `08` and `10`, which `noteworthy.md` already
  records as not what this route does.

## Shipping

Merged by pull request to `EEGR-code/PyOMES`, as recent work has been, in place of a
direct `--no-ff` merge.

- [x] Full test suite green on the branch
- [x] Pull request merged with a merge commit
- [x] `persistent-integrator-shipped` tagged on the merge commit and pushed
- [x] Feature branch deleted
- [x] The design note and this checklist moved to `docs/dev/implementation/shipped/`,
      with a "Shipped" banner on both
- [x] `docs/dev/implementation/upcoming/README.md`'s "Recently shipped" list updated
