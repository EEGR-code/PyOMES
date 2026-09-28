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

- [ ] `git status -sb` clean apart from this phase's docs (the design note, this
      checklist and the `upcoming/README.md` entry)
- [ ] `git log origin/main..main --oneline` empty
- [ ] Branch created off current `main`: `git switch -c vant-hoff-single-source`
- [ ] Those docs committed on the branch as the first commit
- [ ] Baseline full suite recorded here (pass / fail / skip counts)
- [ ] Baseline measurement: the script's "old" values captured for every copy
      #1–#15, and the engine fingerprint recorded

## Checkpoints

- [ ] **1. Kernel, module and the reactions copies.** New
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
- [ ] **2. Bisection engine** (#5–6: `EquilibriumDef.pKas_at_T`,
      `WaterDef.Kw_at_T`). Call the kernel through `vant_hoff_K`; keep their
      `correction == "none"` switch and the `max(Ka, 1e-30)` clamp; drop the
      0.01 K skip. Sanity: outside the 0.01 K band, shifts arithmetic only; inside
      it, the shift equals the previously skipped correction; engine fingerprint
      and BSM2 sentinels recorded; suite green.
- [ ] **3. `ThermoFramework`** (#7–8): `pKa_at_T`, `Kw_at_T` as thin callers of
      the kernel, 0.01 K skip dropped. Sanity as checkpoint 2;
      `chemistry_database.ipynb` re-run if its saved output shows a changed value.
- [ ] **4. Henry and vapour pressure** (#9–11). `test_package_layering.py`: the
      `chemistry/` rule becomes `{chemistry, thermo, units}`, with its docstring.
      `chemistry/partition.py` `_kH_mol_L_atm_from_ref` keeps its unit conversion
      and calls the Henry-constant function; `RaoultEquilibrium.P_sat` and
      `core/phases.water_vapour_P_sat_atm` call `clausius_clapeyron` (their
      parameter values unchanged). Sanity: bit-identical by construction, confirmed
      by measurement; suite green; import graph still acyclic.
- [ ] **5. Models** (#12–15). ADM1 `_arrhenius` calls `arrhenius_factor`;
      `_pKa_NH4` (×3 files) and `_pKa_H2S` call the kernel with their own
      `ΔH / _R_J` (rounded R unchanged). Sanity: bit-identical by construction,
      confirmed; BSM2 sentinels pass unchanged; suite green.
- [ ] **6. Guard test and close-out.** A guard test (style of
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

## Notes

<!-- Baseline counts, measurements, deviations and findings go here as each
checkpoint lands. -->

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
