# Phase Kickoff Checklist — conservation-ledger

> **Status: approved by the repo owner 2026-11-07; not started.** Revised
> 2026-10-07 with the decisions agreed
> that day (see "Agreed"), **and 2026-11-07** against `main` at `546c776`:
> the persistent integrator (#7) and the NR solver fixes (#1, #3 to #5) since
> the audit, what phase 2
> ([`LIQUID_VOLUME_MODELS.md`](LIQUID_VOLUME_MODELS.md)) needs from the
> record, and shipping through pull requests. No branch, no code. Phase 1 of
> [`MASS_BALANCE_CLOSURE.md`](MASS_BALANCE_CLOSURE.md), followed by
> `variable-liquid-volume` and then `engine-closure`.

> Checklist for the `conservation-ledger` phase of
> [`MASS_BALANCE_CLOSURE.md`](MASS_BALANCE_CLOSURE.md), the source of truth for
> the standard, the decisions and the measurements; do not restate them here.
> Where this checklist and the note disagree, this checklist wins: in
> particular the checkpoints are in a different order from the note's
> provisional list. See [`README.md`](README.md)'s "Branching and tagging
> convention". Modelled on
> [`EXPLICIT_SPECIES_RESOLUTION_CHECKLIST.md`](../shipped/EXPLICIT_SPECIES_RESOLUTION_CHECKLIST.md).

**What the phase leaves.** No amount is kept non-negative unless the model
asks for it. Every change to a control volume's contents is recorded with its
source. A conservation monitor checks every step of every control volume,
under every step and system solver, reports only genuine imbalances and names
their source. The engines still fail to book water and still force the liquid
neutral; after this phase the monitor says so, by name and amount.

**Working rules**

- Standing rules are in `CLAUDE.md` (tests, git, line endings, docs). Only the
  phase-specific ones are below.
- **Results change in checkpoint 2** (no clamp, no floors) and are measured
  before and after. **In every other checkpoint they are unchanged**: every
  `n_mol` trajectory and pH in the measurement is bit-identical, with one named
  exception, the adaptive solver in checkpoint 3c. What changes elsewhere is
  which warnings are raised; each checkpoint records the warnings before and
  after, per case. An unexplained shift stops the checkpoint.
- **"Clamp" means the mass clamp.** In this checklist a clamp is a step
  solver's `clamp_fn` (and `proportional_clamp`), which changes amounts. The
  NR solver's limits on its log-activity iterates and Newton step (#1, #3,
  #4) are a different thing, inside one solve, and are not touched.
- **Measurement scripts are committed** in `tests/performance/` (run by hand,
  not collected by pytest, like `bench_adaptive_solver.py`), so they survive
  between sessions and the owner can re-run them. Their output is not
  committed; the numbers that matter are recorded in Notes. Run with `py`
  from the repo root.
- **Measurement.** `tests/performance/ledger_fingerprint.py` builds and runs
  fixed cases and records every phase's `n_mol` trajectory, pH,
  the warnings raised by category and, from checkpoint 3, the record's totals
  by source. Cases: the user-defined model of
  `tests/standalone/test_user_defined_model.py` under the Bisection and NR
  engines, plain and with a feed, a vent and a dose (1 h in 20 steps); the same
  model under `SimultaneousEulerSolver` and `SimultaneousAdaptiveSolver`;
  `docs/tutorials/D2C_workshop/raw_construction.py`'s `build()` (1 h in 200
  steps); ADM1 (0.2 h in 8 steps); BSM2 as in `test_bsm2_reference.py`
  (100 × 0.01 h); the stirred-tank template as a batch and a fed-batch tank
  (2 h in 100 steps); and two control volumes joined by an advective link, under
  the default system solver and each of the others. The six BSM2 sentinel tests
  are the regression check for BSM2.
- **Clamp recorder.** `tests/performance/clamp_audit.py` (a pytest plugin and
  a script runner) records every clamp and floor that fires. It is used
  throughout checkpoint 2. It reuses the guard counting of
  `bench_adaptive_solver.py clamps` where that fits.
- **Balance audit.** `tests/performance/balance_audit.py` checks every
  reaction of the repo's models over every element its species carry and over
  charge, for the Notes' balance table and checkpoint 6.
- One checkpoint at a time: full suite first (owner), edit, sanity checks and
  measurement, report as it is, record notes here, stop, hand over the commit
  commands. The next checkpoint starts when told to. Edits are shown one per
  message, each preceded by a one-line summary.
- A decision listed under "To settle" is agreed with the owner before its
  checkpoint starts, and recorded in the notes.
- Anything found beyond scope goes in `OPEN_WORK.md`, not fixed.

## Pre-flight

- [x] Branch created off current `main`: `conservation-ledger` (at `6bbe1e6`,
      2026-11-07)
- [x] This checklist committed on the branch (it came with `main`, approved
      in #11)
- [x] `MASS_BALANCE_CLOSURE.md`: status line updated (phase 1 started), its
      phase 1 list follows this checklist's order, and its `upcoming/README.md`
      entry follows
- [x] `OPEN_WORK.md` gains the three items logged by the audit and the
      decisions: `ProportionalGasOutlet.last_water_loss_mol_per_h` is read by
      nothing (fixed in `variable-liquid-volume` checkpoint 5, so the entry
      says so); the flux guards that conserve; whether `clamp_fn` should exist
      at all (next to the open question on what an explicit clamp should do)
- [x] `git status -sb` clean (after the pre-flight's last commit, "Finish
      conservation-ledger pre-flight")
- [x] Baseline full suite recorded here (pass / fail / xfail counts):
      2292 passed, 87 skipped, 2 xfailed, 420 warnings, 0 failed (owner's
      terminal, 101 s, 2026-11-07). Run with `py -m pytest`: in the owner's
      PowerShell, `python` is a Python 3.12 without pytest
- [x] The measurement scripts written again and committed in
      `tests/performance/` (decided 2026-11-08): the 2026-10-07 copies of
      `clamp_audit.py`, `balance_audit.py` and `ledger_fingerprint.py` were in
      an earlier session's scratchpad and no longer exist. Each is described
      in the Working rules and Notes below. Checked against the 2026-10-07
      numbers where they overlap, all reproduced: BSM2 clamps `H2O` in 100 of
      100 steps (2.0e-6 mol); the fed-batch configuration clamps O2 in 39
      steps (0.870 mol); `raw_construction.py` clamps `AceticAcid` in 145
      steps (1.1e-4 mol); the balance table, every row. On the iron
      feed-and-drain case of `bench_adaptive_solver.py`, `persistent=False`
      floors trial states 29 times (`Fe(OH)2+`, `FeOH2+`) and the final state
      once, and `persistent=True` floors nothing
- [x] The audit re-run on current `main` (2026-11-08; see Notes, "Re-measured
      2026-11-08"): everything reproduced, with two corrections to the
      preview table. The clamp table, the preview with
      the clamp off and the balance table in Notes, and the closure note's
      "Measurements" that checkpoint 1's markers carry. Every number that moved
      since `99d7f7a` is recorded in Notes with its cause (the persistent
      integrator, the NR fixes). Checkpoint 1 uses these numbers, not the
      older ones
- [x] Baseline measurement captured, and a second capture compares identical
      (2026-11-08, code at `e561a2a`): all 27 cases of `ledger_fingerprint.py`
      run without error, and the two captures are bit-identical. `compare`
      was checked to catch a one-ulp change and one extra warning. The capture
      files are not committed; the code they measure is unchanged until
      checkpoint 1, so the baseline is re-captured from that commit when a
      later session needs it. Two cases differ from the checklist's
      description, by necessity: the ADM1 case seeds its own starting state
      (`ADM1_INITIAL_MOL_L`; the repo has none, and the smoke test it came from
      seeds only strong ions), and the batch tank is the fed-batch tank
      without its feed (`batch_fermenter.py` seeds nothing, so nothing grows).
      Both tanks report pH as NaN: the template's liquid declares no pH
      equilibria (`OPEN_WORK.md`, "The stirred-tank template's pH controller
      never doses"). ADM1's pH stays at 12.635 because it sees only the strong
      ions (`OPEN_WORK.md`, "ADM1's pH sees only the strong ions", logged
      today)

## Checkpoints

### Tests first

- [ ] **1. Closure tests as strict expected failures.** New
      `tests/standalone/test_mass_balance_closure.py`, on the user-defined
      model, for each of the Bisection and NR engines and each of four cases:
      plain, a `LiquidFeed`, a headspace with a `PressureReliefVent`, and a
      base dose (Na+ + OH- through `PHController`). Each test checks, per
      element and for charge, that the change in the totals equals what
      crossed the boundary, at 1e-12 relative. What crossed is taken from a
      small recording wrapper in the test file (it delegates `compute_flux`
      and keeps what it returned), so the tests do not depend on the record
      and later cross-check it. Tests that fail today are strict expected
      failures whose reason carries today's number and cause (the note's
      "Measurements"); those that pass are plain tests. The two H and O
      markers in `test_user_defined_model.py` stay as they are. Sanity: every
      marker's number reproduced by the measurement; suite green; measurement
      unchanged.

      *Done 2026-11-08.* 35 tests: 11 pass, 24 strict expected failures.
      Pass: C and charge in every case and Na with the dose under Bisection;
      Na with the dose under NR; one check that the feed, vent and dose cases
      move material. Expected failures: H and O under Bisection (water not
      booked: +4.7e-6 mol H, +4.3e-8 relative; with the dose, the forced
      neutrality: 5.0e-3 mol OH- dosed, −9.9e-3 mol H and −4.95e-3 mol O);
      C, H and O under NR (water reset to 55.509 mol/L and the residual
      tolerance: H and O +1.5e-4 relative, +6.4e-5 with the dose; C +1.5e-8
      to +2.4e-8); charge under NR (−4e-14 mol, −4e-13 with the dose).
      Decisions taken in writing them:
      - The baseline is taken before the first solve, as the monitor's will
        be (checkpoint 4), so the first solve is checked too. That is why the
        numbers differ from the closure note's "Measurements", which start
        after it: Bisection's H is +4.3e-8 here against +2.5e-8 there (the
        first solve also forms HCO3-), and NR's is +1.5e-4 against −2.1e-5
        (its first solve resets the 55.5 mol of water given to 55.509 mol).
      - The tolerance is 1e-12 of the larger of the total before and after;
        for charge, which starts at zero, of the ionic inventory
        (sum of |z|·n).
      - All four cases run through `Simulation` (the dose needs a
        controller); the dose is recorded where it reaches
        `ControlVolume.apply_external_flux`.
      - The user-defined model is loaded from `test_user_defined_model.py`
        by path: the suite does not make test modules importable by name.
      Sanity: every marker's number from a scratch measurement of the same
      cases; `ledger_fingerprint.py` bit-identical to the baseline in all 27
      cases; this file and `test_user_defined_model.py` 45 passed, 26 xfailed.
      Full suite: the owner's plain run after the commit is not recorded;
      the suite run with the clamp recorder (2a, trial floor off) gives
      2301 passed and 2 failed (both floor unit tests, see 2a) of 2303, 87
      skipped, 26 xfailed, consistent with 11 new passes and 24 new xfails

### Non-negativity

- [ ] **2. No clamp by default; no silent floor.** Results change here. A
      model's definer chooses a step and a solver that keep its amounts
      physical; the package warns when an amount goes negative and runs what
      it was given. Five commits:
  - [x] **2a. What relies on it.** The owner runs the clamp recorder over the
        full suite and the notebooks (commands handed over); the notes list
        every test, tutorial and notebook where a clamp or floor fires, every
        test that asserts a non-negative amount, what each rate law in the
        tutorials and notebooks does with a negative concentration, and what
        the two engines do with a negative total. With the adaptive solver's
        trial floor patched off in a scratch run (as the 2026-10-07 preview
        did for the clamp), the notes also list every run that then produces a
        `NaN`, a failed solve, a step failure or a negative accepted amount,
        and the rate law or engine path responsible. No code.

        *Done 2026-11-09* (code at `4481154`). The owner ran `clamp_audit.py`
        over the full suite and the 32 notebooks, each with and without
        `--no-trial-floor`; the tutorial scripts were run at pre-flight.

        *Where a clamp or floor fires, default settings.*

        | Run | Guard | Fires | Amount |
        |---|---|---|---|
        | `test_bsm2_reference.py` (module fixture) | `clamp_fn` | 100 | `H2O` 2.0e-6 mol |
        | `test_clamping.py`, 7 tests | `clamp_fn`, `floor_nonnegative` | 1-2 each | they test the clamp and floor |
        | `test_core.py`: `test_apply_flux_floors_at_zero`, `test_floor_at_zero_breaks_conservation` | phase floor | 1 each | they test the floor |
        | `test_state_vector.py`: `test_unpack_with_floor_clamps_at_zero` | unpack floor | 1 | tests the floor |
        | `test_adaptive_solver_persistent.py`: `test_a_species_the_engine_does_not_resolve_is_still_floored` | trial floor | 1 | tests the floor |
        | `test_adaptive_solver_persistent.py`: `test_species_added` | trial floor | 6 | `Na+x` 2e-25 mol |
        | `ArXiv_preprint/03_cstr_dilution_rate_sweep.ipynb` | `clamp_fn` | 3,635 | `NH3` 0.033 mol |
        | `D2C_workshop/Example1_mtp_well.ipynb` | `clamp_fn` | 12,822 | O2, `NH3`, `Glucose`; 1.1e-4 mol |
        | `D2C_workshop/Example2_batch_fermenter.ipynb` | `clamp_fn`; phase floor | 4,053; 2,625 | O2, `NH3`, `Glucose`, 14.7 mol; `NH3` 7.4e-14 mol |
        | `D2C_workshop/Example3_CSTR.ipynb` | `clamp_fn` | 263 | O2, `Glucose`; 0.029 mol |
        | `reactions/fba/fba_ecoli_core.ipynb`, `fba_toy.ipynb` | `clamp_fn` | 70; 2 | O2 2.7e-3 mol each |
        | Tutorial scripts (pre-flight) | `clamp_fn` | | `cstr_fermenter.py`, `fed_batch_fermenter.py`, `microplate_fermenter.py`, `raw_construction.py` |

        Every other test, and 22 of the 28 notebooks that ran, fire nothing.
        Four notebooks did not run, for want of optional packages, not
        because of the code: `03_phreeqc_engine_basics.ipynb` and
        `06_phreeqc_benchmark.ipynb` (phreeqpython),
        `aerobic_fermentation_stoichiometry.ipynb` and
        `01_exporting_results.ipynb` (pandas). The only phase floor outside
        the floor's own tests is in `Example2_batch_fermenter.ipynb`, and its
        total is 7.4e-14 mol.

        *With the trial floor off.* Nothing depends on it in the repo. No
        notebook and no model test changes outcome; no run records a
        non-finite state, and no final or unpack floor fires that did not
        before. The suite gives 2 failures, both tests of the floor itself
        (`test_clamping.py::TestFloorNonnegative::test_floors_negative_elements`
        calls `floor_nonnegative` directly, which the recorder treats as a
        trial floor; `test_a_species_the_engine_does_not_resolve_is_still_floored`
        tests the floor's exemption), and both go with the floor in 2d. The
        one run known to depend on the floor is outside the suite: the iron
        feed-and-drain case of `bench_adaptive_solver.py` under
        `persistent=False` (29 trial floors, one final floor; Pre-flight).

        *Tests that assert a non-negative amount.* Of the floor or clamp:
        `test_clamping.py` (lines 90, 135, 164, 177). Of a model's state:
        `test_bsm2_reference.py::test_concentrations_non_negative` (≥ −1e-12,
        decided: exempts `H2O` until phase 2); `test_adaptive_solver_persistent.py`
        line 177; `test_implicit_transport_solver.py` lines 241, 406
        (≥ −1e-12, −1e-10); `test_system_solver.py` lines 452-453, 475,
        605-606; `test_nr_speciation_engine.py` lines 657, 699 (dissolved Ca
        from precipitation). Each is re-run in 2b to 2e.

        *Rate laws in the tutorials and notebooks with a negative
        concentration.* The template, the builders and `rate_laws.py` return
        0 for a substrate or biomass at or below zero (audit below). The D2C
        growth law (`raw_construction.py`) and the one in
        `reaction_system.ipynb` do the same; D2C's sums `AceticAcid` and
        `Acetate-`, so one form may be negative. The Singer–Stumm law of
        `07_iron_oxidation.ipynb` and `08_iron_oxidation_and_precipitation.ipynb`
        (k·[Fe2+]·[O2]·[OH-]²·V) and the user-defined model's 50·[AcOH]·[O2]·V
        are products of two concentrations that can both be negative, which
        gives a positive rate.

        *The engines with a negative total* (the user-defined model, one solve,
        a total of −1e-18 and −1e-6 mol for `AceticAcid` or `CO2`). The
        Bisection engine accepts it without a word and distributes it over
        the forms, writing negative amounts back (−1e-6 total: `AceticAcid`
        −5.5e-7, `Acetate-` −4.5e-7); the total is kept. The NR engine
        replaces it, without a word and with `converged=True`, by about zero:
        every form is written back as 1e-50 mol, which creates the missing
        amount (+1e-6 mol). The two differ, and NR's creates material; 2d
        replaces both with the defined response.
  - [ ] **2b. The step solvers' default.** `SequentialAdvanceSolver` and
        `SimultaneousEulerSolver` take `clamp_fn=None` by default; a clamp is
        something the user passes. A negative amount warns, naming the
        species, the amount and the step; the warning no longer speaks of
        clamping being disabled. `proportional_clamp`'s docstring stops
        claiming it preserves stoichiometry.
  - [ ] **2c. The floors in the phases.** `Phase.apply_flux` loses its `clamp`
        parameter and no longer floors (internal transfer,
        `apply_external_flux`, inter-CV links); every caller is updated.
        `ControlVolume._integrate_reactions` is deleted.
  - [ ] **2d. The adaptive solver's floors.** Both go (see "Agreed",
        2026-11-08). The floor on the state written back after a step goes;
        a negative end amount warns as in 2b. The floor on each trial state
        (`floor_nonnegative` in `_build_integrand`, with its exemption of
        engine-resolved species under `persistent=True`) goes too: rates,
        boundary fluxes and the engine see the integrator's own trial state.
        Each engine gets a defined response to a negative total in place of
        whatever it does today (2a records that): a total negative only at
        roundoff is read as zero, and a larger one raises or warns naming the
        component, never silently; the threshold is settled before this
        commit. The solver's docstring says that keeping amounts in range is
        the rate laws' and the tolerances' job, and that a rate law guards its
        own domain explicitly where it must (as the built-in ones do for a
        substrate at or below zero).

        Re-audited first against the persistent integrator (#7). The success
        path now writes back `integrand.floor(y_end)` (`solve_step` in
        `core/solvers.py`); `unpack(y0_copy, floor=True)` remains on the
        retry and restore paths. A running integrator carries on from its own
        unfloored state while the phases hold the floored one, so the floor's
        material appears at write-back and is gone again in the next step
        without passing through anything the record would see. The audit
        lists every floor on these paths (all go; a restored state that was
        never negative needs none), and checks that
        `_PersistentRun.can_continue` still matches the written-back state.
        Results move for adaptive runs; measured before and after, both with
        the integrator persisting and restarting.
  - [ ] **2e. The models that rested on the clamp.** The tutorials, as settled
        below with the 2a list in hand. BSM2's
        `test_concentrations_non_negative` exempts `H2O` until BSM2 carries its
        solvent water (the volume phase), unless a sentinel check shows that
        seeding water now is free.

      Sanity, each commit: measurement before and after, per case, with what
      moved, by how much and why; the clamp recorder shows no clamp or floor
      firing anywhere with default settings; the tutorials and tests of 2a
      re-run and compared; BSM2 sentinels; suite green, with each re-baselined
      test listed with its reason.

### The record

- [ ] **3. The control volume records every change to its contents.** Each
      control volume keeps, per step, a list of entries: a source kind, a
      label, the phase, and the moles of each species. Element and charge flux
      are derived from the entries with `cv.species`. Nothing is declared by
      the user. Every entry holds the change as computed; when a clamp was
      passed, one clamp entry holds the difference between what was applied
      and what was computed.

      Two properties `variable-liquid-volume` relies on, so they are built in
      here: **entries are kept one per source and never netted** before they
      are stored, because that phase's monitor threshold is scaled on each
      step's gross turnover (the sum of the entries' magnitudes per element),
      which netting would lose; and **the source kind is open**, not a closed
      set, because that phase adds a make-up source (`make_up_to`) and an
      evaporation debit on `ProportionalGasOutlet`. Four commits:
  - [ ] **3a. The record and the sequential solver.** The entry type and the
        per-CV record; `SequentialAdvanceSolver` writes an entry for
        `external_source_terms`, for each boundary, for the engine write-back
        (the amounts either side of `apply_to_phases`), for reactions, for
        internal transfer and for a passed clamp; `cv.apply_external_flux`
        writes one (controller doses, vent controllers, `apply_boundary`);
        `Simulation._initial_solve` writes the first engine entry.
        `ExternalFluxRecord` is the boundary entry.
  - [ ] **3b. The snapshot solver, links and stepping system solvers.**
        `SimultaneousEulerSolver`; `Simulation._apply_links` and so the
        default, Strang and multirate system solvers, one entry on each side
        of a link, with `LinkFlowRecord` as the link entry;
        `ImplicitTransportSystemSolver`, whose transport is the difference
        across its linear solve.
  - [ ] **3c. The `solve_ivp` solvers.** `SimultaneousAdaptiveSolver`
        integrates each state-dependent boundary's running total as an extra
        state (a vent, drain, membrane or gas outlet; constant-rate sources
        are rate × time), so each boundary's entry is exact. Under the
        persistent integrator (#7) these states live in the running
        integrator across macro steps: each step's entry is the change in the
        running total over that step, and a restart (when `can_continue`
        fails) starts the totals from the record, not from zero. The extra
        states also enter the integrator's error control, so step sizes can
        change. The shift this causes in adaptive runs that have such
        boundaries is measured and shown to the owner before the commit, both
        when the integrator persists and when it restarts every step; a
        control volume without one is bit-identical. `MonolithicODESolver` gets no extra state: it has no
        boundaries, and its links cancel in the simulation-wide total that
        checkpoint 4 checks.
  - [ ] **3d. Cumulative totals.** Each boundary, link and dosing controller
        keeps a running total of what it moved, per species, as its own small
        record that a later reservoir type can adopt (the adaptive solver's
        extra states are these totals). The recorder stores each step's
        entries.

      Sanity, each commit: for every measurement case and every step, the
      entries sum to the change in every species' amount to 1e-12 relative
      (new tests, one per step solver and system solver, including a run with
      a clamp passed); the recording wrappers of checkpoint 1 agree with the
      boundary entries; measurement bit-identical except as named in 3c; suite
      green.

### The monitor

- [ ] **4. The monitor compares the change in totals with the record.** Every
      step, for each element and for charge, per step and cumulatively: the
      change, less what crossed the boundary (boundaries, links, doses,
      external source terms), is the imbalance. A warning names the sources
      that make it up and the amount of each: engine write-back, a passed
      clamp, the reactions' own stoichiometric residual, and whatever the
      record does not account for (a change made between steps, or anything
      else), which is reported as unexplained. A simulation-wide total, in
      which links cancel, is checked the same way; under
      `MonolithicODESolver` that is the only check. The baseline is taken
      before the first solve: in `Simulation.run` before `_initial_solve`, and
      at the first `advance` of an unowned control volume.
      `equilibrate_to_pH` keeps resetting the baseline. The check on the
      liquid's absolute net charge stays as it is. Thresholds and the
      `throttle` setting unchanged (`variable-liquid-volume` trials a
      threshold on turnover once liquids hold their water; this phase does not
      change it). The recorder stores each step's imbalance
      by source, so every step can be read after a run whatever was warned.
      `ConservationWarning` messages are rewritten; tests that match their
      text are found first. Sanity: `test_conservation_monitor.py` rewritten
      for the new inputs; ADM1's vent warnings gone; the warnings each case
      raises before and after, with the reason for each difference (BSM2's
      first solve, 17.7 mol of H, is reported for the first time); measurement
      values unchanged; suite green.
- [ ] **5. The monitor runs for every control volume and every solver.** The
      species registry is `cv.species` whatever the reaction system is (a bare
      `KineticReaction`, a black-box model, none). The check is made where
      every path ends: by `Simulation._step` for each of its control volumes
      and for the simulation as a whole, after links, advance and controller
      actions, whatever the system solver; and by `advance()` for an unowned
      control volume, whatever the step solver. The call in
      `SequentialAdvanceSolver` and
      `ReactionSystem.attach_conservation_monitor` go. Sanity: one test per
      step solver and per system solver that a deliberate leak warns and a
      closed model does not; the fed-batch configuration warns for the yeast
      nitrogen and for nothing else; every test and tutorial that now warns is
      listed with the reason (notebook survey handed to the owner); suite
      green; measurement values unchanged.

### Definitions

- [ ] **6. A reaction that does not balance warns.** `balance_elements` keeps
      its meaning: the elements that must close, or construction raises
      (`None` is all of them). When a control volume first advances, next to
      the check for unresolved species, every reaction is checked over every
      element its species carry and over charge, whatever `balance_elements`
      says, and one warning lists each reaction that does not close with its
      residuals. The warning is `UnbalancedReactionWarning`, a new
      `UserWarning` subclass next to `ConservationWarning` (see "Agreed").
      Sanity: the warning names BSM2's three reactions, ADM1's
      thirteen and the template's CHO growth, with the elements and residuals
      of the audit; every other reaction in the repo is silent, or listed with
      the reason; BSM2's and ADM1's docstrings say they implement the
      published benchmarks, whose stoichiometry is on a COD basis; suite
      green; measurement values unchanged.

### Close-out

- [ ] **7. Close-out.** `docs/solvers.md`, `docs/architecture.md`, the
      `monitoring/` and `core/clamping.py` docstrings and the `SolverProtocols`
      notebooks describe the monitor, the record and the clamp as they now
      are, keeping the mass clamp and the NR solver's iterate limits apart in
      their wording; `OPEN_WORK.md` entries this phase resolves are removed or cut down
      ("The conservation monitor does not net out boundary flows", "Some
      inventory floors are silent", "Amounts are kept non-negative by
      default", keeping its open question on what an explicit clamp should
      do); `MASS_BALANCE_CLOSURE.md` gains this phase's measurements and what
      they change for phases 2 and 3; `RESERVOIR_TYPE.md` §3.1 notes that the
      per-boundary total exists. The three measurement scripts in
      `tests/performance/` are reviewed: kept, with docstrings that describe
      them as they now are, or removed if nothing later needs them. Sanity:
      relative links resolve; suite green.

## Agreed

With the repo owner, 2026-10-07, while reviewing the first draft.

- **Non-negativity comes first.** The clamp is removed from the defaults
  before the record is built, so the record and the monitor are developed and
  measured on runs with no clamp in them. Keeping amounts physical is the
  model definer's job (step, solver, rate laws that limit themselves); the
  package warns when an amount goes negative.
- **The clamp stays as an option** a user passes to a step solver. Whether
  `clamp_fn` should exist at all is logged in `OPEN_WORK.md` with the open
  question on what a correct clamp should do: today's clamp scales each
  overdrawn species on its own and leaves the other species of the same
  reaction in full, so it creates mass.
- **Booking.** Entries hold the change as computed. A passed clamp gets one
  entry, applied minus computed, which is what the clamp created. Under the
  snapshot solver the clamp acts on the combined change, so with a clamp
  passed a boundary's total there is what it computed, not what the clamp
  left.
- **`Phase.apply_flux(clamp=...)`** goes; clamping exists only as a step
  solver's `clamp_fn`.
- ~~**The adaptive solver's two floors are different in kind.** The floor on
  the final state adds material and goes. Evaluating rates at the non-negative
  part of a trial state conserves and stays: without it the hydrogen
  inhibition terms of BSM2 and ADM1 would be evaluated at slightly negative
  concentrations, near their pole.~~ Superseded 2026-11-08 (below): the trial
  floor does not conserve in general, and goes too.
- **Flux guards that conserve are left in place** and logged: the links' caps
  at what the source holds, `KineticGasLiquidLink._clamp_flux` and the floor
  in `_equilibrium_flux`, and the drain, links and vents skipping or flooring
  a negative amount when computing a flux. Each limits a flux that is applied
  equally to both sides.
- **The `solve_ivp` solvers.** The adaptive step solver integrates each
  state-dependent boundary's running total as a state. A reservoir's contents
  would be the same state under another name, so a later reservoir type adopts
  these. The monolithic solver is checked at simulation level, where links
  cancel, and needs nothing extra.
- **Two checks, kept apart.** The static check of declared stoichiometry
  (checkpoint 6) fires once, at the first advance; it matters only where a
  model narrowed or switched off the construction-time check. The check for
  accidental inconsistency is the monitor's (checkpoints 4 and 5): every step,
  every control volume and the simulation as a whole, with each step's
  imbalance kept by the recorder.
- **BSM2's non-negativity test** exempts `H2O` until phase 2, unless seeding
  water now is shown to be free.

With the repo owner, 2026-11-07.

- **The warning category for checkpoint 6** is a new `UserWarning` subclass,
  `UnbalancedReactionWarning`, next to `ConservationWarning`. The
  definition-time warning is about how a model is written; the monitor's is a
  run-time check of the books. Kept apart, a user can silence BSM2's and
  ADM1's known benchmark imbalances without silencing the monitor.
- **What phase 2 needs from the record** (checkpoint 3): entries kept one per
  source and never netted, and an open source kind.
- **Shipping goes through a pull request** to `EEGR-code/PyOMES`, and the
  branch is deleted only once the merge is confirmed (see "Shipping").

With the repo owner, 2026-11-08.

- **Decisions are general, not shaped to one model.** BSM2 and ADM1 are test
  cases for a decision, not its reason.
- **The adaptive solver's trial floor goes** (checkpoint 2d), replacing the
  2026-10-07 entry above. Runge–Kutta and BDF methods preserve every linear
  conservation law (element totals, charge) to roundoff, but only when every
  rate is evaluated at the method's own trial state; a floor, clip or
  transform applied first removes that guarantee. It also creates material
  wherever the engine sums several of the state's species into one total:
  flooring a negative form raises the total the engine sees (measured in
  `OPEN_WORK.md`, "Flooring a resolved species creates material": total Fe
  drifts 9.8e-5 with the floor, 1.2e-7 without). Keeping amounts in range is
  done where it is visible: rate laws guard their own domain explicitly, each
  engine has a defined response to a negative total, tolerances keep trial
  negatives small, and a negative accepted amount warns. Checkpoint 2a
  measures which runs depend on the floor before anything is committed.
- **Enforced positivity, if a model needs it, is a later opt-in** that
  conserves: rejecting a step whose accepted state breaks a sign constraint
  and retrying with a smaller step, as SUNDIALS' CVODE and IDA do, not
  editing the state. Logged in `OPEN_WORK.md`; 2a's evidence decides whether
  it is needed and in what form.

## To settle

**Before 2e: the tutorials whose results rest on the clamp** (numbers in the
audit below). For each tutorial on the 2a list: run as given and warn, give
the model the missing limitation, or pass a clamp explicitly.

**Before 2d: what the engines do with a negative total, and the runs that
need the trial floor.** The roundoff threshold below which a negative total
reads as zero, and whether a larger one raises or warns. For each run 2a finds
failing without the floor: the rate law guards its own domain, the run's
tolerances change, or it is listed as needing the opt-in positivity mechanism.

**Before 3c is committed: the adaptive solver's shift.** The owner sees the
measured shift in adaptive runs that have a state-dependent boundary, with
the integrator persisting and restarting.

## Notes

**Re-checked 2026-11-07, `main` at `546c776`.** Read from the code, not
measured: the numbers below are re-measured at pre-flight. 33 commits since
the audit; the code changes are in `core/solvers.py` (the persistent
integrator, #7), the NR engine and solver (#1, #3 to #5),
`control/descriptors.py` (a parameter-change counter the persistent
integrator reads) and `core/simulation.py`. In the table below, "Adaptive
solver's final floor" now reads `integrand.floor(y_end)` on success and
`unpack(floor=True)` on retry and restore; the rest of the table and of the
findings still hold. Phase 2 (`variable-liquid-volume`) fixes the gas
outlet's water loss.

**2026-11-08.** The table row "Adaptive solver's trial state" says "rates
only; conserves". That is wrong in general: the integrator's state holds
species the engine sums into totals, and flooring one raises the total
(`OPEN_WORK.md`, "Flooring a resolved species creates material", found
2026-11-01). Since #7, `persistent=True` exempts the engine-resolved species;
`persistent=False` does not. Checkpoint 2d now removes the floor.

**Re-measured 2026-11-08, code at `396176b`** (pre-flight). The closure
note's "Measurements" for the user-defined model, which checkpoint 1's
markers carry, are unchanged by the NR fixes and the persistent integrator
(both cases run the sequential solver). Drift over 1 h in 20 steps, relative:
Bisection C −1.65e-16, H +2.46e-8, O +2.46e-8; NR C +1.63e-8, H −2.14e-5,
O −2.14e-5. A dose of 1e-3 mol followed by one solve, under either engine:
Na+ with OH- moves pH from 3.3885 to 3.8664 and loses 2.00e-3 mol H and
1.00e-3 mol O; OH- alone leaves pH at 3.3885 and loses 1.00e-3 mol H and
1.00e-3 mol O (the note gives only the H); H+ alone leaves pH at 3.3885 and
loses 1.00e-3 mol H. The clamp and balance numbers reproduced are listed
under Pre-flight.

The owner re-ran the tutorials (`clamp_audit.py scripts`, with and without
`--no-clamp`, outputs compared line by line). The clamp table below
reproduces in every row: `cstr_fermenter.py` 14,782 clamps of `AceticAcid`
(0.303 mol); `fed_batch_fermenter.py` 611 of `AceticAcid` (28.3 mol) and 108
of O2 (3.95 mol); `microplate_fermenter.py` O2 (1.39e-5 mol) and `Glucose`
(6.7e-9 mol); `raw_construction.py` 145 (1.1e-4 mol); `batch_fermenter.py`
none. With the clamp off, no floor fires in any of them. The preview
reproduces, with two corrections to the 2026-10-07 table: the biomass
1.8004 → 1.7998 g/L is `cstr_fermenter.py`'s, whose output therefore is not
identical (`microplate_fermenter.py`'s biomass is 5.0130 g/L in both runs;
its DO, p(O2) and glucose lines are as the table gives them); and
`raw_construction.py` prints pH 5.100 in both runs, not 5.099 without the
clamp (its yeast moves 3.7782e-2 → 3.7675e-2 mol, as the table gives).
Whether the table misattributed the biomass line or the CSTR changed since
cannot be told from the old run.

**Audit before drafting (2026-10-07, `main` at `99d7f7a`).** Read from the
code, then measured with two scratchpad scripts (`clamp_audit.py`,
`balance_audit.py`). Searched: `PyOMES/`, `models/`, `tests/`,
`docs/tutorials/` (`.py`, `.ipynb`) and the planning docs.

*For checkpoint 3: every path by which an amount in a phase changes.*

| Path | Where | Applied through | Recorded today |
|---|---|---|---|
| `external_source_terms` | each step solver; passed only by direct `cv.advance` callers (six tests, one notebook), never by `Simulation` | sequential: clamp, then `apply_flux(clamp=False)`; snapshot: summed into the combined change; adaptive: a constant rate in the right-hand side | no |
| `cv.boundaries` (`GasFeed`, `PressureReliefVent`, `MembraneGasBoundary`, `LiquidFeed`, `LiquidDrain`, `ProportionalGasOutlet`) | each step solver | as above; adaptive integrates the instantaneous flux | snapshot solver only (`AdvanceResult.boundary_records`), before its clamp |
| Engine write-back | `apply_to_phases` → `_refresh_derived`: sequential once at the start; snapshot and adaptive on a copy, then on the live phases at the end; `Simulation._initial_solve`; `compute_rhs` (monolithic); `equilibrate_to_pH` through `advance(dt_h=0)` | overwrite | no |
| Clamp | `clamp_fn` in the sequential solver (per sub-step) and the snapshot solver (combined change) | scales removals | `AccuracyWarning` only |
| Reactions | each step solver; `compute_rhs` | sequential writes `n_mol` directly | `AdvanceResult.reaction_sources` (before the clamp) |
| Internal transfer | `step_internal_transfer` (sequential); in the combined change or right-hand side otherwise | `apply_flux` with its floor | `TransferDiagnostics` (sequential) |
| Inter-CV links | `Simulation._apply_links`: default (before advance), Strang (two half steps), multirate (M sub-steps) | `Phase.apply_flux` on each side, with its floor | `LinkFlowRecord`, per link, not per CV |
| Inter-CV links, implicit | `ImplicitTransportSystemSolver`: linear solve, `_unpack_state` | overwrite of every amount | no (returns no link records) |
| Inter-CV links, monolithic | in the right-hand side; `_unpack_extended_state` | overwrite | no |
| Controller actions | `Simulation._apply_controller_action` → `cv.apply_external_flux`: `PHController` doses (liquid), three vent controllers (gas), after the advance | `Phase.apply_flux` with its floor | `ControlAction.flux_applied` (the request, not what the floor left) |
| Adaptive solver's final floor | `unpack(floor=True)` on the final or restored state | overwrite; adds material | no |
| Adaptive solver's trial state | `floor_nonnegative` on each trial state before rates are computed; the floored values are also written to the live phases during the integration and replaced at its end | rates only; conserves | no |
| Direct writes | `equilibrate_to_pH`; ADM1's and BSM2's seeding after the CV is built; `restore_state`; user code. Between steps, not through the CV | `n_mol[...] =` | no |

Found on the way: `MonolithicODESolver` applies no boundary and no internal
transfer (its right-hand side has reactions and links only), and never calls
the monitor; `ProportionalGasOutlet.last_water_loss_mol_per_h` is set and read
by nothing, so the water its docstring says the liquid loses is not removed;
profiles change parameters only, never an amount. The first is already in
`HPLC_CV_GRAPH_TEMPLATE.md`; the second goes to `OPEN_WORK.md` at pre-flight.

*For checkpoint 6: how reactions pass the balance check today.* A reaction
checks only the elements in `balance_elements` (`None` means every element its
species carry) and raises if one does not close. Passing `()` checks nothing,
and charge is never checked (`validate_balance(check_charge=...)` has no
caller). BSM2 passes `()` on all 19 kinetic reactions; ADM1 passes `()` on
decay and disintegration and a per-reaction list on the rest. Residuals over
every element and charge:

| Model | Reactions | Do not close |
|---|---|---|
| BSM2 | 27 | 3: `R1_disintegration` (C −0.584, H +0.337, N −0.518, O −0.085), `R2_hyd_ch` (H +2, O +1: one water), `R4_hyd_li` (C +3.03, H +10.1, O +4.08) |
| ADM1, default | 32 | 13: twelve decay reactions (N +0.2 each), `5_disintegration` (C +0.134, H +0.265, N −0.087, O −0.042) |
| ADM1, ethanol | 36 | 15: fourteen decay reactions, `5_disintegration` |
| Stirred tank, yeast on a CHO basis | 1 | 1: N +0.139 (checked C, H, O) |
| Stirred tank, CHNO; the three databases; `raw_construction.py` | 1; 3, 7, 9; 7 | 0 |

None of these reactions has a charge residual. BSM2's other sixteen kinetic
reactions (eight uptake, seven decay, `R3_hyd_pr`) close although unchecked:
its drift in every step comes from three reactions. `BlackBoxReactionModel`
has its own run-time check (`on_imbalance`) and is left alone. Not audited:
the reactions defined in tests and notebooks.

*For checkpoint 2: which runs clamp today.* Clamp recorder, default settings:

| Run | Clamped | Steps | Removal cut by |
|---|---|---|---|
| BSM2, 100 × 0.01 h | `H2O` | 100 of 100 | 2.0e-6 mol |
| `raw_construction.py` | `AceticAcid` (the protonated form) | 145 | 1.1e-4 mol |
| `fed_batch_fermenter.py`, 10 h in 1000 steps | `AceticAcid` | 611 | 28.3 mol |
| | `O2` | 108 | 3.95 mol |
| `cstr_fermenter.py` | `AceticAcid` | 14,782 | 0.303 mol |
| `microplate_fermenter.py` | `O2` | 568 | 1.4e-5 mol |
| | `Glucose` | 2 | 6.7e-9 mol |
| `batch_fermenter.py` | none | | |
| 545 tests: `tests/validation`, `test_simulation.py`, `test_builder.py`, `test_cv_advance.py`, `test_system_solver.py`, `test_user_defined_model.py`, `test_model_species.py`, `test_bsm2_reference.py` | none, except BSM2 above | | |

BSM2 clamps water in every step, not one as `OPEN_WORK.md` says. The
fed-batch tutorial clamps its substrate in most steps, which the note's 2 h
case did not show. No silent floor fired in any of these runs
(`Phase.apply_flux`, the adaptive solver's trial and final floors), and
`ControlVolume._integrate_reactions` has no caller in the repo. Not measured:
the rest of the suite and every notebook (checkpoint 2a).

*A preview with the step solvers' clamp off* (their default patched in the
scratch run; phase floors left on):

| Run | Today | No clamp |
|---|---|---|
| `fed_batch_fermenter.py` | biomass 1.404 g/L; substrate 0.000 g/L | biomass 1.023 g/L; substrate −0.001 g/L |
| `microplate_fermenter.py` | DO 1.74e-4 mol/L; p(O2) 0.1608 atm; biomass 1.8004 g/L | DO −7.45e-2 mol/L; p(O2) 0.1896 atm; biomass 1.7998 g/L; glucose −0.006 g/L |
| `raw_construction.py` | pH 5.100; yeast 3.7782e-2 mol | pH 5.099; yeast 3.7675e-2 mol |
| `cstr_fermenter.py`, `batch_fermenter.py` | | printed output identical |

Every script ran to the end. The microplate result is the template's growth
law at work: it has no oxygen term (`Ko2_gL` exists on
`ReactionBuilder.monod_aerobic_growth`, not on the template), so the clamp is
what limits growth by oxygen today.

*For checkpoint 2: what the rate laws do with a negative concentration.*

| Rate law | Substrate or biomass ≤ 0 | Other negative input |
|---|---|---|
| `rate_laws.py`: Monod, Contois, Andrews, ContoisAndrews, Tessier, Moser, Blackman (shared `make_rate_fn`) | rate 0 | none read |
| `rate_laws.py`: `DualSubstrateMonod` | rate 0 | secondary substrate negative: rate 0 |
| BSM2 `rate_first_order`, `rate_uptake`, `rate_c4_competitive` | rate 0, tested on the sum of an acid's forms, so one form may be negative | `_I_IN` floors NH3 at 1e-12; `_I_h2` and `_I_nh3` are 1/(1 + S/K_I): above 1 for a negative inhibitor, singular at S = −K_I and negative beyond it (K_I for H2 is 3.5e-6 to 1e-5 kg COD/m³) |
| ADM1 `rate_monod`, `rate_inhibited`, `rate_syntroph`, `rate_ace_meth`, `rate_h2_meth`, `rate_ace_srb`, `rate_h2_srb`, `rate_decay`, `rate_first_order` | rate 0 | free NH3 and H2S are floored at 0; `Ki/(Ki + P)` and `Ki_H2/(Ki_H2 + H2)` have the same pole |
| `raw_construction.py` growth law, the user-defined model's `50·[AcOH]·[O2]·V` | read in checkpoint 2a with the notebooks' laws | a product of two negative concentrations is a positive rate |

No built-in law returns a negative rate for a negative substrate. The hazard
is an inhibition term with a negative inhibitor, in BSM2 and ADM1; neither
clamps an inhibitor in the runs measured.

## Shipping

- [ ] Full test suite green on the branch
- [ ] Branch pushed: `git push`
- [ ] Pull request opened against `EEGR-code/PyOMES` (`gh pr create --repo
      EEGR-code/PyOMES --base main --head conservation-ledger ...`) and merged
      with a merge commit (`gh pr merge <n> --repo EEGR-code/PyOMES --merge`)
- [ ] Merge confirmed: `gh pr view <n> --repo EEGR-code/PyOMES --json state`
      prints `MERGED`. Nothing below runs before this: deleting the branch
      while its pull request is open closes the request unmerged
- [ ] `git switch main`, then `git pull`; `git log --oneline -1` shows the
      merge commit
- [ ] `git tag conservation-ledger-shipped`, then
      `git push origin conservation-ledger-shipped`
- [ ] `git branch -d conservation-ledger`, then
      `git push origin --delete conservation-ledger`
- [ ] Move this checklist to `docs/dev/implementation/shipped/` with a
      "Shipped" banner; the design note stays in `upcoming/` until its third
      phase ships
- [ ] Update `docs/dev/implementation/upcoming/README.md`: the note's entry
      says phase 1 has shipped, and "Recently shipped" gains an entry
