# Phase Checklist — persistent-integrator

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
- [ ] 4. **Parameter triggers:** generation counter in `control/descriptors.py`
      (trigger 2) and the `Simulation` signal (trigger 3). A test for each trigger that
      fails when it is removed, including `link.kLa["O2"]`.
- [ ] 5. **Robustness:** retry on failure, pickling, one solver across several CVs,
      monitor deltas. Check `tests/standalone/test_package_layering.py` still passes.
- [ ] 6. **`env.t_h` continuous under this solver**, with a test. Log `Simulation`'s
      end-of-step `t_h` in `OPEN_WORK.md`.
- [ ] 7. **Acceptance run:** every benchmark model meets the gate in the design note,
      all `solve_ivp` methods included. Stop here if one does not.
- [ ] 8. **Flip the default**; add `restart_each_step=True` as the opt-out. Move the
      three `solve_ivp`-patching tests to the `_make_stepper` seam. Update docs.
- [ ] 9. **Owner re-runs the listed notebooks**; review any changed outputs.

## Shipping

- [ ] Full test suite green on the branch
- [ ] `git checkout main`
- [ ] `git merge --no-ff persistent-integrator -m "Merge persistent-integrator: <summary>"`
- [ ] `git tag persistent-integrator-shipped <commit-hash>`
- [ ] `git push` and `git push --tags`
- [ ] `git branch -d persistent-integrator` and
      `git push origin --delete persistent-integrator`
- [ ] Move the design note and this checklist to `docs/dev/implementation/shipped/`,
      add a "Shipped" banner to both
- [ ] Update `docs/dev/implementation/upcoming/README.md`'s "Recently shipped" list
