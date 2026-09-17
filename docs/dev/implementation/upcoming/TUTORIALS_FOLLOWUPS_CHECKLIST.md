# Phase Kickoff Checklist — tutorials-followups

> Status: branch `tutorials-followups` cut from `main` (with `tutorials-reorg` merged) on
> 2026-09-17. Five bugs/cleanup items found during that phase (see
> [`../shipped/TUTORIALS_REORG_CHECKLIST.md`](../shipped/TUTORIALS_REORG_CHECKLIST.md) once it
> ships) but deliberately not fixed there — they're debugging/investigation work, not file
> moves, and didn't block that phase's own goal. Split into its own phase so `tutorials-reorg`
> could ship as one clean, complete unit. When this is picked up: this file already exists on
> `main` (committed alongside the `tutorials-reorg` merge), so just cut the branch — no need to
> commit the checklist itself as a first commit, it's already there.

## Pre-flight

- [x] `git status -sb` clean (no stray uncommitted work left over from a previous task)
- [x] `git log origin/main..main --oneline` empty (nothing unpushed sitting around from earlier work)
- [x] On current `main`, with `tutorials-reorg` already merged (confirm this file's "Status" line
      above no longer says "no branch cut yet" — update it once the branch exists)
- [x] Branch created off current `main`: `git checkout -b tutorials-followups`

## During

- [ ] Checkpoints tracked below as they land, one commit per checkpoint
- [ ] **If work stalls or is paused before shipping:** add a status banner to the top of this
      file immediately — what's built, what's tested, why it stopped, which commit it's on.

### Checkpoints

- [x] 1. Fix `docs/tutorials/reactions/chemistry_database.py`: crashes on `dataclasses.replace(
      ..., activity_model=...)` — `ThermoFramework.__init__()` no longer accepts `activity_model`
      as a kwarg (pre-existing, confirmed broken in the pre-`tutorials-reorg` `demos/` copy too;
      likely drift from the `stirred-tank-template` `ThermoFramework` refactor). Find the current
      constructor API and update the call. Once it runs clean, convert to a notebook (exec each
      cell for real and capture genuine stdout — see `tutorials-reorg`'s checkpoint 10 for the
      harness pattern) and update `docs/tutorials/reactions/README.md`. Sanity check:
      `python docs/tutorials/reactions/chemistry_database.py` exits 0 before conversion; the
      resulting notebook's cells show genuine (non-error) output. **Root cause:**
      `ThermoFramework` was refactored to hold a `liquid_activity: LiquidPhaseModel` field
      instead of string kwargs — `use_activity`/`activity_model` are now derived read-only
      properties, not settable. Fix: `dataclasses.replace(AD_BASIC.thermo,
      liquid_activity=DaviesLiquidModel())`. Script confirmed exit 0, then converted to
      `chemistry_database.ipynb` via a one-off exec-and-capture-stdout harness (each cell run
      for real against a shared namespace, genuine stdout captured into `outputs`) — output
      verified to match the script's actual run (21 species, pKa 6.3500 → 6.3065, etc.). Old
      `.py` deleted; `README.md` row and "Running" section updated.
- [x] 2. Fix `docs/tutorials/reactions/partition_model.py`: crashes calling `.beta()` on
      `HenryEquilibrium` — that method doesn't exist on the class that replaced the now-deprecated
      `HenryPartition` (pre-existing, confirmed broken in the pre-`tutorials-reorg` `demos/` copy
      too). Find the current equivalent and update the call. Once it runs clean, convert to a
      notebook and update `docs/tutorials/reactions/README.md`, same as checkpoint 1. Sanity
      check: same as checkpoint 1. **Root cause:** `HenryPartition`/`.beta()` are deprecated —
      `HenryPartition(...)` is now a function that emits `DeprecationWarning` and returns a
      `HenryEquilibrium`, whose equivalent method is `.partition_ratio()` (same signature, same
      "ratio >> 1 -> mostly liquid" semantics, just renamed). Updated the script to import and
      construct `HenryEquilibrium` directly (not the deprecated alias) and renamed all `.beta(`
      calls to `.partition_ratio(`. Confirmed exit 0, then converted to `partition_model.ipynb`
      via the same exec-and-capture-stdout harness as checkpoint 1 — output verified to match the
      script's actual run. Old `.py` deleted; `README.md` row and "Running" section updated (no
      `.py` files remain in this folder, so the section collapsed to just `jupyter lab`).
- [x] 3. Investigate `docs/tutorials/D2C_workshop/raw_construction.py`'s control-loop issue:
      runs (an unrelated import bug was fixed in `tutorials-reorg`) but produces pH 12.089
      against a `PHController(setpoint=5.0)`, plus `ConservationWarning`s for O/C/H and charge
      balance exceeding their stated thresholds. Never previously observable — the script never
      ran successfully before `tutorials-reorg` (confirmed: the pre-move `demos/` copy raises
      `ModuleNotFoundError` before reaching the simulation), so this isn't a regression, it's
      untested code's first real run. Likely a `PHController` tuning issue or a stoichiometry
      mismatch; not yet diagnosed. Doesn't block anything else in this folder — confirmed none
      of `Example1_mtp_well.ipynb`/`Example2_batch_fermenter.ipynb`/`Example3_CSTR.ipynb` (the
      actual tutorial content there) import or depend on `raw_construction.py`. Sanity check:
      re-run and confirm pH settles near the 5.0 setpoint with no `ConservationWarning`s.
      **Root cause (two independent issues, confirmed by isolating each):**
      1. **pH runaway.** The `PHController` doses `chemical_id="H3PO4"` (acid) and
         `base_chemical_id="NaOH"` (base). `NaOH` is a recognised strong-corrector alias
         (`ControlVolume.apply_external_flux` resolves it to `Na+`, which shifts pH purely via
         charge balance — no equilibrium reaction needed) but `H3PO4` is not; it only shifts pH
         by actually dissociating, which requires the phosphate ladder
         (`H3PO4 ⇌ H2PO4- ⇌ HPO4-- ⇌ PO4---`) to be declared in the CV's `reaction_system`. The
         manually-built `ReactionSystem` here never declared it (unlike `AD_BASIC`, which
         `templates/batch_fermenter.py` uses via `StirredTankBuilder`), so dosed `H3PO4`
         accumulated as inert neutral acid (confirmed: 0.356 mol by t=5h, zero effect on pH) —
         the acid half of the PI loop was a silent no-op. Base dosing kept working, so pH only
         ever climbed, unopposed, to 12.089. Fix: added `make_phosphate_ladder()` (by-value copy
         of `bioprocess_basic.py`'s `eq_phosphate_1/2/3`, same log_K values: -2.15/-7.20/-12.35),
         wired into the `ReactionSystem`. Confirmed fixed: pH rises from 3.234 to ~5.0 within
         0.5 h and holds 4.72-5.01 for the rest of the 5 h run (sampled at 8 points across the
         run, not just the endpoint).
      2. **ConservationWarnings.** Unrelated to the above — their magnitudes were byte-identical
         before and after the phosphate-ladder fix. Root cause: `ConservationMonitor`'s per-step
         threshold clamps to an absolute `1e-8 mol` for any element pool under 1 mol total
         (`scale = max(abs(baseline), 1.0)` in `PyOMES/monitoring/conservation.py`), and this is
         a small-scale system (`V_liq`=1.6 L, `V_gas`=0.4 L, sub-1-mol O/C/H pools) with a real
         `GasFeed` air sparge at 1 vvm — genuine, expected open-system O/C/N throughflow (~1 min
         headspace turnover), not a stoichiometry bug. Same situation `ArXiv_preprint/
         _generate_notebooks.py` already documents and silences for its own small-pool kinetic
         runs; applied the same fix here (`warnings.filterwarnings("ignore",
         category=ConservationWarning)`, with a comment explaining why). Drive-by cleanup: also
         replaced the file's remaining `HenryPartition` (deprecated alias) with `HenryEquilibrium`,
         same issue as checkpoint 2. Confirmed: script now exits 0 with zero warnings of any kind.
- [x] 4. Investigate and add pytest coverage for
      `tests/validation/speciation/{07_iron_oxidation,08_iron_oxidation_and_precipitation}.ipynb`
      (currently zero coverage). Their own P1 prediction (>70% Fe2+ conversion) is already marked
      `[CHECK]` (failed) in the notebook's own output — actual conversion is 0.2%, likely because
      the demo equilibrates to pH 6.5 while the Singer-Stumm rate constant/narrative assumes pH
      ~4.7. Diagnose that inconsistency first (fix the notebook's setup, or its narrative claim,
      whichever is actually wrong), then write `test_*.py` coverage for the corrected behavior,
      mirroring `tutorials-reorg` checkpoint 3's approach for notebooks 01/03/06 (see
      `tests/validation/speciation/test_carbonate_phosphate_benchmarks.py` and
      `test_saturation_index.py` for the established style/conventions to follow). Sanity check:
      `pytest tests/validation/` green, including new iron-oxidation tests; the notebook's own P1
      cell shows `[PASS]` or its claim is corrected to match reality. **Root cause (this notebook
      lives in `tests/validation/` — the point is a real, checkable prediction, not a demo-pacing
      choice — so this was treated as a genuine physics bug and verified against literature, not
      curve-fit to hit a target):**
      1. **Wrong pH assumed throughout the narrative.** Section 1 described the reaction at the
         medium's *natural* pH 4.68, but `cv.equilibrate_to_pH('KOH', 6.5)` has always deliberately
         run the kinetics at pH 6.5 — every downstream prediction (rate, dominant Fe(III) species)
         described a scenario the code never actually simulated.
      2. **The real bug: literature `k` applied to the wrong quantity.** Singer & Stumm (1970)'s
         rate law is `r = k [Fe2+] p(O2) [OH-]^2` — calibrated against **p(O2) in atm**. But
         `ControlVolume._build_reaction_environment` only ever builds a kinetic `rate_fn`'s
         `env.concentrations` from the **liquid** phase (`control_volume.py:700-748`, `:924-1003`
         — confirmed by reading both call sites; gas-phase state is never merged in). The old
         `K_SS = 5e11` was that literature `k` applied directly to aqueous `[O2]` in mol/L with no
         unit conversion — understating the rate by the Henry's-law factor (`kH_O2 ≈ 1.317e-3
         mol/(L·atm)`, ~760x) regardless of pH, since `[OH-]^2` alone can't close a gap that large.
         Verified the literature constant via web search rather than trusting memory:
         `k = 1.33e12 M^-2 atm^-1 s^-1` (Singer & Stumm 1970, Eq. 22), cited via USGS's PHREEQC
         docs [Example 9](https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-71.htm),
         which implements the identical rate law. Corrected `K_SS = k / kH_O2 ≈ 3.635e18`, computed
         explicitly in the notebook, not applied by fiat.
      3. **O2 was a fixed, depletable 0.75 mM pool** (capping conversion regardless of kinetics)
         **and Fe2+'s stated loading (2.00 mM) never matched what `phase-cv`'s stock-dilution code
         actually produces (0.88 mM).** Fixed both: replaced the fixed O2 pool with a real 1000 L
         air headspace (`GasPhase` + `EquilibriumTransferModel`), same "pseudo-unlimited reservoir"
         idiom as `ArXiv_preprint/_generate_notebooks.py`'s kinetic CO2 demo (verified <0.01%
         reservoir drawdown over the run); corrected the Fe2+ table entry to 0.88 mM.
      4. **O2's initial condition was itself a bug** — caught by noticing an O2(aq) spike right
         after t=0 in the notebook's plot: O2 started at 0 and let the `EquilibriumTransferModel`
         jump it to its ~0.276 mM equilibrium value on the solver's very first internal step, an
         artificial discontinuity, not a physical initial condition (a medium open to air has
         always been in contact with it). Confirmed this alone roughly quadrupled solver work
         (~70 vs ~15 accepted BDF steps for an otherwise-identical run). Fixed by initializing O2
         at its Henry's-law equilibrium value directly.
      5. **Solver.** This reaction is genuinely stiff (t½ ≈ 11 min inside a 2h run) — switched from
         the default explicit solver to `SimultaneousAdaptiveSolver(method='BDF',
         use_engine_jacobian=True)`, a real fix (Fe mass-balance drift dropped ~10x) that also
         eliminates every `ConservationWarning` the explicit solver accumulated from many large
         per-step NR re-solves at this reaction speed — not a suppression. It does trip a
         confirmed-false-positive `AccuracyWarning` (`check_scipy_rejections` judges
         `nfev`/accepted-steps against a flat threshold calibrated for explicit methods; implicit
         methods' Newton-iteration overhead inflates `nfev` by design, unrelated to genuine step
         rejection — confirmed by sweeping rtol/atol, max_step, and both BDF/Radau: ratio stuck at
         ~1.0 regardless, and Radau scored *worse*, the opposite of what a real rejection problem
         would show) — documented and suppressed in both notebooks and the test file.
      Both notebooks now show all 5 predictions `[PASS]` with zero warnings. Added
      `tests/validation/speciation/test_iron_oxidation.py` (9 tests: rate-law pH-direction and
      4th-order-structure regression guards, K_SS-derivation check, O2-initial-condition regression
      guard, and 4 batch-integration checks against P1/P3/P4/P5) — `pytest tests/validation/` green
      (27 passed, 8 skipped pre-existing).
      **Three core-package findings surfaced and logged as separate design notes** (not fixed here,
      out of `tutorials-followups` scope, each independently pick-up-able):
      [`PHCONTROLLER_CORRECTOR_VALIDATION.md`](PHCONTROLLER_CORRECTOR_VALIDATION.md) (from
      checkpoint 3 — `PHController` should warn when its configured corrector can't actually shift
      pH), [`REACTION_ENVIRONMENT_PHASE_EXPOSURE.md`](REACTION_ENVIRONMENT_PHASE_EXPOSURE.md)
      (kinetic `rate_fn`s can't see gas-phase partial pressure, forcing the manual Henry's-law
      conversion above), and
      [`SCIPY_REJECTION_CHECK_SOLVER_AWARENESS.md`](SCIPY_REJECTION_CHECK_SOLVER_AWARENESS.md)
      (the `AccuracyWarning` false positive above — `check_scipy_rejections` needs a
      solver-family-aware threshold).
- [x] 5. Clean up `demos/usecases/03_cstr_dilution_rate_sweep.ipynb` — an orphaned duplicate, not
      produced by `demos/usecases/_generate_notebooks.py` and not referenced anywhere in the repo
      (confirmed by grep during `tutorials-reorg`), likely a leftover from before this notebook
      was originally curated into `docs/tutorials/`. Re-confirm it's still unreferenced, then
      delete. Sanity check: `grep -r "usecases/03_cstr_dilution_rate_sweep"` (outside this
      checklist) returns nothing before deleting. Re-confirmed: zero references outside this
      checklist — even `demos/usecases/0_README.ipynb`'s own notebook table already links to the
      `docs/tutorials/ArXiv_preprint/03_cstr_dilution_rate_sweep.ipynb` copy, not the local file.
      Deleted.

## Shipping

- [x] Full test suite green on the branch (`pytest`) — 2038 passed, 36 skipped, 0 failed.
- [ ] `git checkout main`
- [ ] `git merge --no-ff tutorials-followups -m "Merge tutorials-followups: ..."`
- [ ] `git tag tutorials-followups-shipped <commit-hash>`
- [ ] `git push && git push --tags`
- [ ] `git branch -d tutorials-followups` and `git push origin --delete tutorials-followups`
- [ ] Move this checklist to `docs/dev/implementation/shipped/`, add a "Shipped" banner
- [ ] Update `docs/dev/implementation/upcoming/README.md`'s "Recently shipped" list
