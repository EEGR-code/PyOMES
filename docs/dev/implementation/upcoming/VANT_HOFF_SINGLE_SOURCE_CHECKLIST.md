# Phase Kickoff Checklist — vant-hoff-single-source

> Checklist for [`VANT_HOFF_SINGLE_SOURCE.md`](VANT_HOFF_SINGLE_SOURCE.md), the
> source of truth for motivation, inventory (copies #1–#15), design and decisions;
> do not restate it here. Where this checklist and the note disagree, this
> checklist wins. See [`README.md`](README.md)'s "Branching and tagging
> convention". Modelled on
> [`../shipped/ACTIVITY_MODEL_PARAMETER_CHECKLIST.md`](../shipped/ACTIVITY_MODEL_PARAMETER_CHECKLIST.md).

**Working rules**

- This phase is **not** bit-identical by rule. Numbers may move where the
  formulation stays correct and each shift is attributable either to
  floating-point arithmetic or to the one named change of formulation (dropping
  the 0.01 K skip in copies #5–8). Every checkpoint measures its shifts and
  records them here; an unexplained shift stops the checkpoint.
- **Measurement.** A scratchpad script (not committed) evaluates each touched
  function, old and new, over a fixed grid: T from 273.15 to 373.15 K in 1 K steps
  plus T_ref ± {1e-12, 1e-9, 0.005, 0.0099, 0.0101, 0.02} K; enthalpies (or E/R)
  covering both signs and 0; the reference values each caller actually uses. It
  reports, per function, the number of differing values and the maximum absolute
  and relative shift, split into "within 0.01 K of T_ref" and "elsewhere".
  Arithmetic shifts are expected at the ~1e-15 relative level; anything larger
  outside the 0.01 K band is unexplained. The engine-output fingerprint from the
  previous phase (`activity_fingerprint.py`, which solves at 25 and 37 °C) is
  re-run as a regression check; its shifts are recorded, not required to be zero.
- The repo owner runs every `git` command (branch, add, commit, push, merge,
  tag). Read-only local git (`status`, `log`, `diff`, `show`, `ls-files`) is
  fine; nothing that touches the network. Git commands are handed over one per
  line (PowerShell 5.1: never chained with `&&`); commit messages use two `-m`
  flags (a strapline and one body paragraph) and no attribution lines. The repo
  is in OneDrive: a commit may ask "Rename ... index.lock ... failed. Should I
  try again? (y/n)"; answer `y`. No git command is run while a commit may be
  waiting.
- Files are moved with a plain filesystem move, not `git mv`. Before staging a
  moved or renamed file, check how the old path is stored (`git ls-files --eol`):
  from an `i/crlf` file, stage with `git -c core.autocrlf=false add <path>`; from
  an `i/lf` file, a plain `git add`. Check `git diff --cached -M --stat` before
  committing.
- One checkpoint at a time: run the full suite first, edit, run the sanity checks
  and the measurement, report what changed and the results as they are, record
  notes here, stop, and hand over the commit commands. Do not start the next
  checkpoint until told to. Edits are shown one per message, each preceded by a
  one-line summary; mechanical multi-file rewrites use a script with exact
  replacements and required counts, and the diff is shown.
- Full suite: `python -m pytest -p no:cacheprovider -q`.
- No compatibility aliases or shims: `thermo.equilibrium_constants` stops
  existing; every importer is updated in place.
- Most files involved are CRLF in the working tree. Bulk edits preserve each
  file's line endings and exact trailing bytes; check for bare LF after every edit.
- Leave `docs/dev/implementation/shipped/` and `docs/dev/ideas/` untouched.
- Docs and docstrings describe current behaviour only: no phase or checkpoint
  labels, no pointers to this checklist or the design note.
- Anything that looks like a bug or dead code beyond scope is logged in
  `OPEN_WORK.md`, not fixed.

## Pre-flight

- [x] `git status -sb` clean apart from this phase's docs (the design note, this
      checklist and the `upcoming/README.md` entry)
- [x] `git log origin/main..main --oneline` empty
- [x] Branch created off current `main`: `git switch -c vant-hoff-single-source`
- [x] Those docs committed on the branch as the first commit (`d9bfa63`)
- [x] Baseline full suite recorded here (pass / fail / skip counts)
- [x] Baseline measurement: the script's "old" values captured for every copy
      #1–#15, and the engine fingerprint recorded

## Checkpoints

- [x] **1. Kernel, module and the reactions copies.** New
      `thermo/temperature_correction.py`: the kernel (`ln_correction`, E/R in K)
      and `vant_hoff_K`, `vant_hoff_log_K` moved from `equilibrium_constants.py`
      unchanged in behaviour, plus the Henry-constant function and
      `clausius_clapeyron` (added now, used from checkpoint 4).
      `equilibrium_constants.py` deleted; importers updated (`nr/tableau.py`,
      `nr/engine.py`, `bisection/acid_base.py`); `thermo/__init__.py` docstring and
      exports if needed. `tests/standalone/test_equilibrium_constants.py` renamed
      `test_temperature_correction.py` and extended to the new functions.
      `test_package_layering.py`'s `thermo` rule updated to name
      `temperature_correction` where it names `equilibrium_constants`.
      `arrhenius_factor` in `reactions/kinetic/rate_laws.py`, calling the kernel.
      #2: `constraint.vant_hoff_log_K(constraint, T)` renamed (name settled at the
      checkpoint) to a one-line wrapper over `thermo`; its callers
      (`interphase.py:472`, tests) updated. #3–4 (`plots.py`) call the kernel.
      Sanity: `vant_hoff_K`/`vant_hoff_log_K` bit-identical to the old ones;
      measured shifts for #2–4 are arithmetic only; suite green; the old module
      path fails to import.
- [x] **2. Bisection engine** (#5–6: `EquilibriumDef.pKas_at_T`,
      `WaterDef.Kw_at_T`). Call the kernel through `vant_hoff_K`; keep their
      `correction == "none"` switch and the `max(Ka, 1e-30)` clamp; drop the
      0.01 K skip. Sanity: outside the 0.01 K band, shifts arithmetic only; inside
      it, the shift equals the previously skipped correction; engine fingerprint
      and BSM2 sentinels recorded; suite green.
- [x] **3. `ThermoFramework`** (#7–8): `pKa_at_T`, `Kw_at_T` as thin callers of
      the kernel, 0.01 K skip dropped. Sanity as checkpoint 2;
      `chemistry_database.ipynb` re-run if its saved output shows a changed value.
- [x] **4. Henry and vapour pressure** (#9–11). `test_package_layering.py`: the
      `chemistry/` rule becomes `{chemistry, thermo, units}`, with its docstring.
      `chemistry/partition.py` `_kH_mol_L_atm_from_ref` keeps its unit conversion
      and calls the Henry-constant function; `RaoultEquilibrium.P_sat` and
      `core/phases.water_vapour_P_sat_atm` call `clausius_clapeyron` (their
      parameter values unchanged). Sanity: bit-identical by construction, confirmed
      by measurement; suite green; import graph still acyclic.
- [x] **5. Models** (#12–15). ADM1 `_arrhenius` calls `arrhenius_factor`;
      `_pKa_NH4` (×3 files) and `_pKa_H2S` call the kernel with their own
      `ΔH / _R_J` (rounded R unchanged). Sanity: bit-identical by construction,
      confirmed; BSM2 sentinels pass unchanged; suite green.
- [x] **6. Guard test and close-out.** A guard test (style of
      `test_gas_constant_single_source.py`) fails on the pattern `1/T − 1/T_ref`
      in `PyOMES/` or `models/` outside `thermo/temperature_correction.py`, with a
      self-check that it catches every form found in the inventory. Delete the
      `OPEN_WORK.md` entry "A third van 't Hoff copy in
      `reactions/equilibrium/constraint.py` differs from `thermo` in the last bit";
      update the van 't Hoff parts of "Sweep the package for each fundamental
      constant" and the "Development-history references…" entry's pointer to
      `equilibrium_constants.py:16`; `docs/architecture.md`'s `thermo/` tree line.
      Sanity: repo-wide sweep for `equilibrium_constants` and for the pattern;
      relative links resolve; suite green.
- [x] **7. Rename the test references.** In
      `tests/standalone/test_temperature_correction.py`, `_legacy_vant_hoff_K` →
      `_reference_vant_hoff_K`, `_legacy_vant_hoff_log_K` →
      `_reference_vant_hoff_log_K`, and `TestBitIdenticalToLegacy` →
      `TestBitIdenticalToReference`. They are independent written-out formulas that
      pin the kernel's arithmetic bit for bit, not leftovers; "legacy" dates from the
      earlier refactor they once guarded. Test-only. Sanity: no `_legacy_` left in
      the file; that test file and the full suite pass.

## Notes

**Pre-flight (2026-09-28, on `d9bfa63`).** Baseline suite: 2130 passed, 0
failed, 166 warnings. Engine fingerprint (`activity_fingerprint.py --style new`):
`c649ad00…a38876ac763` over 1,452 values, unchanged from the previous phase.

Baseline measurement (`vh_measure.py capture`, scratchpad, not committed): 26,338
cases, 27,244 values, covering every copy #1–#15 through its real code path; #3
and #4, which compute inline, are read back from the curves `plot_vant_hoff` and
`plot_speciation` draw. A second run compares identical, so the measurement is
repeatable. It also resolves the differences the design note predicts, from the
baseline alone:

- #1 vs #2 (the log10(e) constant): 899 of 4,520 values differ, by at most
  3.55e-15 in log10 K.
- #1 vs #5 (the 0.01 K skip, as pKa): inside the band, up to 4.65e-4 (at the
  grid's largest enthalpy, 80 kJ/mol; 2.9e-4 at 50 kJ/mol); elsewhere at most
  3.55e-15, i.e. arithmetic only.

**Checkpoint 1 (2026-09-28).** Suite: 2139 passed, 0 failed (+9, all new tests in
`test_temperature_correction.py`). Engine fingerprint identical
(`c649ad00…a38876ac763`). `PyOMES.thermo.equilibrium_constants` no longer
imports. Measurement against the baseline:

- #1 `vant_hoff_K`, `vant_hoff_log_K`: bit-identical after the move.
- #2 (now `constraint_log_K_at`): 899 of 4,520 values differ, at most 3.55e-15
  in log10 K — exactly the 899 values where #1 and #2 differed at baseline, so #2
  now equals #1. Arithmetic (the log10(e) constant).
- #3 `plot_vant_hoff`: 42 of 303 curve values differ, at most 8.9e-16; #4
  `plot_speciation`: 62 of 612 fraction values, at most 5.6e-16. Arithmetic
  (÷R then ×log10(e) instead of ÷(R·ln 10)).
- #5–#15: untouched, 0 differ.

Implementation notes:

- `vant_hoff_delta_ln_K` is replaced by the kernel `ln_correction(E_over_R_K, T_K,
  T_ref_K)` (it had no caller outside its own tests). `vant_hoff_K` and
  `vant_hoff_log_K` keep their exact arithmetic, including `vant_hoff_K`'s
  `np.exp`; `henry_constant`, `clausius_clapeyron` and `arrhenius_factor` use
  `math.exp`, as the copies they replace do, so checkpoints 4–5 can be
  bit-identical. Their tests compare with each copy's written-out expression
  exactly.
- `plot_vant_hoff` calls `vant_hoff_log_K` once per grid point: the curve is over
  an array, and `vant_hoff_log_K`'s edge checks take one temperature.
- `arrhenius_factor` is exported from `PyOMES.reactions` with the rate laws.
- **Deviation:** `docs/architecture.md`'s `thermo/` tree line (planned for
  checkpoint 6) is updated here, because it named a file this checkpoint deletes.
  The layering test's synthetic self-check string for `interphase` was also
  updated to the new name.

**Checkpoint 2 (2026-09-28).** `EquilibriumDef.pKas_at_T` and `WaterDef.Kw_at_T`
take their exponent from `ln_correction`; the 0.01 K skip is gone. Suite: 2139
passed, 0 failed; BSM2 sentinels (6) pass unchanged; engine fingerprint identical
(the engines solve at exactly T_ref and at 37 °C, both outside the band).
Measurement:

- #5 `pKas_at_T`: elsewhere 0 of 4,200 differ (bit-identical); inside the band
  210 of 320 differ, at most 4.65e-4 in pKa — the same 210 values and size the
  pre-flight #1-vs-#5 comparison predicted. Checked against the exact formula:
  every band value now equals `-vant_hoff_log_K(...)` to within 1.8e-15 (the
  Ka ↔ pKa round trip).
- #6 `Kw_at_T`: elsewhere 0 of 840 differ; inside the band 42 of 64, at most
  1.07e-3 relative (the same shift in K space: 4.65e-4 × ln 10).
- Every other copy unchanged.

Deviations from this checkpoint's wording, both to keep the attribution clean:

- The copies call the kernel `ln_correction` directly with `math.exp`, not
  `vant_hoff_K`, which uses `np.exp`: the two are not guaranteed equal in the last
  bit, and calling the kernel keeps everything outside the band bit-identical.
  The formula still lives only in the kernel.
- The 0.01 K skip is replaced by the kernel's own 1e-10 K rule (as in
  `vant_hoff_log_K`), not removed outright: at T exactly equal to T_ref the
  correction is exactly zero, but `pKas_at_T` converts pKa → Ka → pKa, and that
  round trip is not exact in floating point (4.76 can come back as
  4.760000000000001), which would shift every solve at 25 °C in the last bit for
  no physical reason.

**Checkpoint 3 (2026-09-28).** `ThermoFramework.pKa_at_T` and `Kw_at_T` take
their exponent from `ln_correction`, with the same two choices as checkpoint 2
(`math.exp`; 0.01 K rule → 1e-10 K). Suite: 2139 passed, 0 failed; engine
fingerprint identical. Measurement:

- #7 `pKa_at_T`: elsewhere 0 of 4,200 differ; band 210 of 320, at most 4.65e-4 in
  pKa; every band value equals the exact formula to within 1.78e-15.
- #8 `Kw_at_T` (fixed ΔH = 55 900 J/mol): elsewhere 0 of 105 differ; band 6 of 8,
  at most 7.49e-4 relative (4.65e-4 × ln 10 × 55.9/80, as expected).
- `chemistry_database.ipynb` calls `pKa_at_T` only at 308.15 K, outside the band,
  so its saved output is unchanged and it was not re-run.

**Checkpoint 4 (2026-09-28).** `chemistry/partition.py`'s
`_kH_mol_L_atm_from_ref` keeps its unit conversion and calls `henry_constant`;
`RaoultEquilibrium.P_sat` and `core/phases.water_vapour_P_sat_atm` call
`clausius_clapeyron` with their own parameters (unchanged). The layering test's
`chemistry/` rule now allows `thermo` (docstring updated; the test is renamed
`test_chemistry_imports_only_units_and_thermo`). Suite: 2139 passed, 0 failed;
layering and import-graph-acyclic tests pass; engine fingerprint identical.
Measurement: #9, #10 and #11 bit-identical (0 of 1,356, 113 and 113 values
differ), as the exact-expression tests from checkpoint 1 predicted.

**Checkpoint 5 (2026-09-28).** ADM1's `_arrhenius` keeps its guard and calls
`arrhenius_factor` (imported from `PyOMES.reactions`); the four pKa functions
(`adm1/base.py` `_pKa_NH4`, `_pKa_H2S`; `bsm2.py` and `bsm2_direct.py`
`_pKa_NH4`) call `ln_correction` with their own `ΔH / _R_J`, so the rounded R
stays theirs. Suite: 2139 passed, 0 failed; BSM2 sentinels (6) pass unchanged;
engine fingerprint identical. Measurement: #12–#15 bit-identical. A sweep of
`PyOMES/` and `models/` for the `1/T − 1/T_ref` pattern now finds only the kernel
(`thermo/temperature_correction.py:41`): all 15 copies are gone.

**Checkpoint 6 (2026-09-28).** New guard
`tests/standalone/test_temperature_correction_single_source.py`: it walks each
module's syntax tree in `PyOMES/` and `models/` for a difference of two
reciprocals `1/a − 1/b` (literal 1 numerators, either order, any denominator) and
requires exactly one, in the kernel. Self-checks cover every form the inventory
found and the near-misses it must ignore (docstrings, comments, strings, sums,
other numerators). Run against the pre-phase code (`d9bfa63`), the same scanner
finds 16 occurrences in 11 files — the 15 inventory rows, #13 holding two — so it
would have caught every copy. Tests are not scanned, on purpose (the references in
`test_temperature_correction.py`); the guard's docstring says so. Suite: 2142
passed, 0 failed (+3).

`OPEN_WORK.md`: deleted the entry "A third van 't Hoff copy in
`reactions/equilibrium/constraint.py` differs from `thermo` in the last bit"; the
R-consumers list and the development-history note now name
`temperature_correction.py`; "Sweep the package for each fundamental constant"
records that the formula is single-source and guarded. Found while checking that
last cross-reference: "No single source for water's physical constants" listed
only `RaoultEquilibrium`'s vapour-pressure values, not the second set in
`core/phases.py` (BSM2: 0.0313 bar, ΔH_vap/R = 5290 K; 2.6 % apart at 25 °C), which
the design note said it covered; the entry now records both. Remaining
`equilibrium_constants` mentions are history (that note, and the thermo phase's
"Recently shipped" entry in `upcoming/README.md`).

**Checkpoint 7 (2026-09-28).** In `test_temperature_correction.py`:
`_legacy_vant_hoff_K` → `_reference_vant_hoff_K`, `_legacy_vant_hoff_log_K` →
`_reference_vant_hoff_log_K`, `TestBitIdenticalToLegacy` →
`TestBitIdenticalToReference`, plus three mentions this checkpoint's list did not
name: the module docstring's ``_legacy_*`` and the two methods
`test_vant_hoff_K_matches_legacy_exactly` / `test_vant_hoff_log_K_matches_legacy_exactly`
(now `…_matches_reference_exactly`). No "legacy" left in the file and no reference
to the old names anywhere else. That file: 32 passed; suite: 2142 passed, 0
failed.

## Shipping

- [ ] Full test suite green on the branch
- [ ] `git switch main`
- [ ] `git merge --no-ff vant-hoff-single-source -m "Merge vant-hoff-single-source: <summary>"`
- [ ] `git tag vant-hoff-single-source-shipped`
- [ ] `git push origin main` and `git push origin vant-hoff-single-source-shipped`
- [ ] `git branch -d vant-hoff-single-source`
- [ ] Move the design note and this checklist to `docs/dev/implementation/shipped/`,
      add a "Shipped" banner to both
- [ ] Update `docs/dev/implementation/upcoming/README.md`: remove the "Design
      discussions" entry and add one to "Recently shipped"
