# Phase Kickoff Checklist — explicit-species-resolution

> Checklist for [`EXPLICIT_SPECIES_RESOLUTION.md`](EXPLICIT_SPECIES_RESOLUTION.md),
> the source of truth for the rule, decisions, settled questions and audit
> results; do not restate them here. Where this checklist and the note disagree,
> this checklist wins. See [`README.md`](README.md)'s "Branching and tagging
> convention". Modelled on
> [`../shipped/VANT_HOFF_SINGLE_SOURCE_CHECKLIST.md`](../shipped/VANT_HOFF_SINGLE_SOURCE_CHECKLIST.md).

**Working rules**

- **The user-defined route comes first.** A model can be defined with every
  species and reaction built by the user, with some taken from a database,
  or with a whole database; the first is the baseline and is verified before
  the others. `tests/standalone/test_user_defined_model.py` builds a model
  with no database or `common_species` import, is run first at every
  checkpoint, and grows as entry points gain `species=` / `reactions=` /
  `chemistry_db=`. Tests and examples that exercise something other than
  databases define their own species; a database is used where the point is
  the database route.
- Results are unchanged by rule, with two named exceptions, each in its own
  checkpoint: the Bisection engine's write-back (checkpoint 2) and the yeast
  molar masses (checkpoint 7). Every checkpoint measures and records; an
  unexplained shift stops the checkpoint.
- **Measurement.** A scratchpad script (`species_fingerprint.py`, not
  committed) builds and runs five fixed cases: the stirred-tank template as a
  batch and a fed-batch tank (Monod growth on acetic acid, pH and DO control
  configured, 2 h in 100 steps); `docs/tutorials/D2C_workshop/raw_construction.py`'s
  `build()` (acetate and phosphate equilibria, CO2 partition, `PHController`
  dosing H3PO4 and NaOH; 1 h in 200 steps); ADM1 (0.2 h in 8 steps); and BSM2
  as in `test_bsm2_reference.py` (100 × 0.01 h). It records every phase's
  `n_mol` trajectory, pH and the warnings raised by category, and reports per
  case the number of differing values and the maximum absolute and relative
  shift. Initial amounts are in mol, so a molar-mass change shows only where
  the model uses the molar mass. The script's build functions are the only
  part that follows API changes; each change is recorded in the notes. The
  six BSM2 sentinel tests are the regression check for BSM2.
- **Fallback recorder.** The pytest plugin used for the audits (scratchpad,
  `-p audit_plugin`) is re-run after each fallback removal to confirm the
  fallback is no longer reached.
- The repo owner runs every `git` command (branch, add, commit, push, merge,
  tag). Read-only local git (`status`, `log`, `diff`, `show`, `ls-files`) is
  fine; nothing that touches the network. Git commands are handed over one per
  line (PowerShell 5.1: never chained with `&&`); commit messages use two `-m`
  flags (a strapline and one body paragraph) and no attribution lines. The repo
  is in OneDrive: a commit may ask "Rename ... index.lock ... failed. Should I
  try again? (y/n)"; answer `y`. No git command is run while a commit may be
  waiting.
- One checkpoint at a time: run the full suite first, edit, run the sanity
  checks and the measurement, report what changed and the results as they are,
  record notes here, stop, and hand over the commit commands. Do not start the
  next checkpoint until told to. Edits are shown one per message, each preceded
  by a one-line summary; mechanical multi-file rewrites use a script with exact
  replacements and required counts, and the diff is shown.
- Full suite: `python -m pytest -p no:cacheprovider -q`.
- No compatibility aliases or shims: `compounds.py`, `common_species.py` and
  `_STRONG_CORRECTOR_ION` stop existing; every importer is updated in place,
  including notebooks and `models/`.
- Notebooks are edited as source and re-run from the scratchpad where affected;
  runs that cannot finish under the CPU cap are reported as not run.
- Most files are CRLF in the working tree. Edits preserve each file's line
  endings and trailing bytes; check for mixed endings after every edit.
- Leave `docs/dev/implementation/shipped/` and `docs/dev/ideas/` untouched.
- Docs, docstrings and error messages describe current behaviour only: no
  phase or checkpoint labels, no pointers to this checklist or the design note.
  Tests match some message fragments; search before rewording.
- Anything that looks like a bug or dead code beyond scope is logged in
  `OPEN_WORK.md`, not fixed.

## Pre-flight

- [x] Branch created off current `main`: `explicit-species-resolution`
- [x] Design note updated with the audits and settled questions and committed
      on the branch (`8a0b7cb`)
- [x] This checklist committed on the branch (`8014a50`)
- [x] `git status -sb` clean
- [x] Baseline full suite recorded here (pass / fail / skip counts)
- [x] Baseline measurement captured (stirred tank, D2C script, ADM1, BSM2)

## Checkpoints

### String stoichiometry

- [x] **1. Drop the `common_species` seed.** `_parse_stoichiometry` looks up
      only the caller's `species=`; `_get_common_species` and the `_cs_mod`
      import go. The unknown-id message lists the caller's ids (or says none
      were given) and names the fix. The `species` parameter docs on
      `EquilibriumReaction`, `KineticReaction`, `KspEquilibrium` and
      `_parse_stoichiometry` stop mentioning `common_species`. Callers define
      their own species: `tests/standalone/test_stoichiometry.py` (16 sites;
      `test_unknown_species_message_lists_common` becomes a test of the new
      message); `docs/tutorials/D2C_workshop/Example3_CSTR.ipynb` cell 14
      (re-run). New `tests/standalone/test_user_defined_model.py` (see working
      rules). Sanity: suite green; recorder shows no stoichiometry fallback
      hits; measurement unchanged.

### Baseline route

- [ ] **2. Bisection write-back.** The engine writes back only species it was
      given (its declared equilibria and the strong ions present in its
      input), including the conjugate bases of user-declared acids;
      `_CANONICAL_WRITEBACK_SPECIES` is deleted. Sanity: the four strict
      expected failures in `test_user_defined_model.py` now pass and their
      markers are removed (H and O conserved to 1e-9); recorder shows no
      zero-filled ions created in `n_mol`; measurement recorded for every
      case, shifts explained (e.g. a change in which keys the state vector
      packs); BSM2 sentinels; suite green.

### Stirred-tank template and `Chemical`

- [ ] **3. Species for the template's compounds.** In the database modules
      (not `common_species`): Yeast, Yeast_CHO, AceticAcid, PropionicAcid,
      ButyricAcid, CitricAcid, O2 and N2 in the bioprocess database; CH4 and
      H2 in the anaerobic-digestion database. Yeast and Yeast_CHO keep the
      registry's explicit MWs (26.868, 24.626) for now so this checkpoint
      changes no result; the acids and gases compute MW from atoms. Registry
      entries not carried over (solvent alias `Water`, the salts) are listed
      in the notes. Sanity: each new `Species` equals the registry entry's
      atoms, and its MW equals the registry's except CitricAcid (−0.001);
      suite green; measurement unchanged.
- [ ] **4. `aerobic_growth` takes `Species`.** All six participants (substrate,
      biomass, N source, O2, CO2, H2O) are `Species` arguments; the separate
      id / atoms / MW arguments and `species_overrides` go, as does the
      `common_species` import and the local `O2`. Every caller updated
      (the factory, tests, notebooks; listed in the notes). Sanity: the built
      stoichiometry is identical entry by entry for the template's default
      configurations; suite green; measurement unchanged.
- [ ] **5. Template takes the model's chemistry; no default database.**
      `StirredTankBuilder` and `StirredTankFactory.create_volume` take
      `species=` / `reactions=` / `chemistry_db=`; the `AD_BASIC` fallback
      goes. Partition models come from the database passed; with none and no
      Henry constant, the existing "No partition model" error is raised.
      ADM1 passes `AD_BASIC`; BSM2 passes nothing (it gives every Henry
      constant). Every other builder and `create_volume` caller (audit 3)
      passes the database it gets today, so results are unchanged; a switch to
      `BIOPROCESS_BASIC` in a tutorial is made only after measuring it.
      `transfer_species`'s docstring stops citing `_HENRY_PARAMS` (and the
      OPEN_WORK entry for it is deleted). Sanity: suite green; measurement
      unchanged; ADM1 and BSM2 sentinels unchanged.
- [ ] **6. Names resolve against the species passed; `compounds.py` goes.**
      `OrganismConfig` / `SubstrateConfig` resolve ids (organism, substrate,
      N source) against the template's species and accept a `Species` in
      place of `atoms=` / `MW=`; a miss raises the error decision 4 describes.
      Default ids go from the configs and from `.organism()` / `.substrate()`.
      Deleted: `compounds.py`, the `Chemical` / `ChemicalRegistry` exports in
      `PyOMES/__init__.py`, `tests/standalone/test_compounds.py`, the
      `conftest.py` fixture, `ChemistryConfig.acid_pKas` and its two
      assertions in `test_configs.py`. Callers updated: the template
      tutorials (`.py` and `batch_fermenter.ipynb`), `test_builder.py`,
      `test_configs.py`. Docs: `PyOMES/README.md`'s `compounds.py` row,
      `docs/architecture.md`, `README.md`'s test table, the comments citing
      the registry's Yeast_CHO. Sanity: repo-wide sweep for
      `ChemicalRegistry`, `compounds`, `acid_pKas` (template); suite green;
      recorder shows no registry use; measurement unchanged.
- [ ] **7. Yeast molar masses from atoms.** Yeast and Yeast_CHO drop their
      explicit MW (24.834, 22.593). The tutorials' `MW_yeast = 26.868`
      (`cstr_fermenter.py`, `fed_batch_fermenter.py`) and any other copy of
      26.868 / 24.626 (sweep, including tests) follow. Sanity: the stirred-tank
      measurement moves; the notes record by how much and confirm the shift is
      only in quantities that pass through the yeast MW; ADM1 and BSM2
      unchanged; suite green (tests that pin gram-based yeast values are
      re-baselined with the reason recorded).

### `ControlVolume` species set and warnings

- [ ] **8. `cv.species`.** `ControlVolume` takes `species=` / `reactions=` /
      `chemistry_db=` and builds the read-only `cv.species` in `__init__` from
      those plus the reaction stoichiometries. A new species-level conflict
      check raises `SpeciesConflictError` when the same id arrives with
      different data. `_collect_species_registry` builds from `cv.species`;
      `_common_species_catalog` is deleted. The factory passes the template's
      species through. `Simulation._warn_thermo_mismatch` unchanged. Sanity:
      suite green; the three tests that reached the catalog (audit 3) now
      get those ids from species passed to them; measurement unchanged.
- [ ] **9. Model species.** BSM2 declares `S_cat` / `S_an` as `Species`
      (`atoms={}`, charge ±1) in its species table; ADM1's `_get_species`
      raises on an unknown id instead of building biomass. Sanity: BSM2 and
      ADM1 unchanged; suite green.
- [ ] **10. Correctors from the model's species.** `equilibrate_to_pH` and
      `apply_external_flux` resolve the corrector from the model's species;
      `_STRONG_CORRECTOR_ION` is deleted. The resolution rule is agreed with
      the repo owner before this checkpoint starts (see notes). Callers
      updated: the D2C workshop notebooks and script, the template tutorials'
      `PHController` configuration, `test_cv_advance.py`,
      `test_iron_oxidation.py`, `test_simulation.py`, `conftest.py`. Sanity:
      pH-control results unchanged; suite green.
- [ ] **11. `UnresolvedSpeciesWarning`.** New `UserWarning` subclass in
      `monitoring/`, exported next to `ConservationWarning`. Warns where an
      id enters (CV construction, feeds, dosing, template set-up), naming
      the id and the call, and in the monitor once per id whenever it first
      appears. The notes list every test and notebook that now warns and
      why, including the NR `H2O` and PHREEQC write-backs left in place.
      Sanity: suite green (tests that turn warnings into errors are fixed by
      declaring the species, not by filtering); measurement unchanged.

### Henry and Raoult

- [ ] **12. `Species` objects only.** `HenryEquilibrium` /
      `RaoultEquilibrium` species fields take `Species` or `None`;
      `_resolve_species` and the `"H2O"` defaults go. In the same checkpoint:
      `databases/anaerobic_digestion.py`, ADM1's water link, the two
      engine-basics notebooks and the test sites in audit 2. Sanity: suite
      green; recorder shows no string resolution; measurement unchanged.

### Retire `common_species`

- [ ] **13. Definitions into the databases.** The 26 definitions move into
      the database modules (`aqueous.py`: water, carbonate, ammonia; which
      module owns phosphate, sulfate, sulfide, the spectator ions and the
      metal ions is settled in the notes); each database lists every species
      it offers. Every importer (62 files) is rewritten by script; the module
      and its `chemistry/__init__.py` export are deleted, with the
      `species.py` docstring that cites it. Sanity: repo-wide sweep for
      `common_species` (only history left); every `Species` object identical
      to before by value; suite green; measurement unchanged.
- [ ] **14. Close-out.** Sweep for `_get_common_species`,
      `_common_species_catalog`, `ChemicalRegistry`, `_STRONG_CORRECTOR_ION`,
      `_CANONICAL_WRITEBACK_SPECIES`, `species_overrides`; `README.md`,
      `PyOMES/README.md`, `docs/architecture.md`; the strong-ion and
      `PHController` notes' cross-references; `OPEN_WORK.md` entries this
      phase resolves. Sanity: relative links resolve; suite green.

## Notes

**Pre-flight (2026-09-30, on `8014a50`).** Baseline suite: 2144 passed, 0
failed, 166 warnings. Baseline measurement (`baseline.json`, scratchpad); a
second capture compares identical, so the measurement is repeatable:

| Case | Series | Values | Warnings |
|---|---|---|---|
| `d2c_raw` | 26 | 5,226 | none |
| `st_batch` | 10 | 1,010 | 1 `AccuracyWarning` |
| `st_fedbatch` | 10 | 1,010 | 1 `AccuracyWarning` |
| `adm1` | 56 | 504 | 2 `ConservationWarning` |
| `bsm2` | 45 | 4,500 | 1 `AccuracyWarning`, 10 `ConservationWarning` |

The two stirred-tank cases have pH NaN at every step and no Na+: the template
tank carries no acid-base equilibria (the database is used only for partition
models), so its configured pH controller never doses. They still cover growth,
gas transfer and the template's name resolution; `d2c_raw` covers equilibria
and NaOH dosing (pH 3.23 → 5.00, Na+ 0 → 0.0196 mol). Logged in `OPEN_WORK.md`.

**Checkpoint 1 (2026-09-30).** `_parse_stoichiometry` builds its lookup from
`species=` only; `_get_common_species` and the `common_species` import are
gone. A miss now reads, e.g., `'Na+' is not among the species passed
(available: CO2, H+, H2O, ...). Add a Species for it to species=, e.g.
species={..., 'Na+': Species(id='Na+', ...)}.`, or "no species were passed"
when none were. The `species` parameter docs on `EquilibriumReaction`,
`KineticReaction`, `KspEquilibrium` and `_parse_stoichiometry` say ids are
looked up there only; `_resolve_species`'s docstring no longer cites the
deleted function.

A first version passed `AQUEOUS_DEFAULT.species` in the tests and
`AD_BASIC.species` in the notebook (committed as `1962192`). Reworked in a
follow-up commit on the repo owner's direction: which database a model uses is the user's choice when
defining it, and the baseline to prove first is a model whose species are all
user-defined.

- `test_stoichiometry.py`: the species it names (H2O, H+, OH-, CO2, HCO3-,
  NH3, acetic acid, acetate) are defined in the file; each test passes exactly
  the species its string names; no database import. Tests and section headers
  named for "common species" renamed; the message test checks the caller's ids
  are listed. Two new tests: no `species=` raises even for water, and an id
  missing from `species=` raises.
- `Example3_CSTR.ipynb` cell 14 declares its six water and carbonate species
  next to `GLUCOSE` and `PEKILO` and passes them to its three reactions
  (exact replacements; no saved outputs). Run from the scratchpad against the
  old and new code: every code cell runs in both, printed output identical
  apart from wall-clock time, same warnings (eight `ConservationWarning`s and
  one `AccuracyWarning`, already present).
- New `test_user_defined_model.py`: water, carbonate and acetate equilibria
  and a kinetic acetate oxidation, all string stoichiometry over species
  built in the file, in a liquid-only `ControlVolume`. Checks that the
  reactions and the conservation monitor hold only those `Species` objects,
  that the model advances to a finite pH, oxidises acetate and conserves
  each element. Four checks fail today and are marked strict expected
  failures, with the current behaviour as the reason:
  - the Bisection engine writes nine ids the model never declared into
    `n_mol`, at zero (Ca++, Cl-, Co++, K+, Mg++, Mn++, Mo7O24------, Na+,
    Zn++);
  - it never writes back the declared `Acetate-` (not on its fixed list), so
    `n_mol` keeps all 0.01 mol as `AceticAcid` while `H+` (4.1e-4 mol)
    reflects the dissociation, and the monitor warns about charge;
  - H and O drift by -9.7e-8 and +2.5e-8 relative over 1 h (C is conserved
    to 3e-16).

  That makes the Bisection write-back a prerequisite for the baseline route,
  so it moves from the `ControlVolume` checkpoints to checkpoint 2, and the
  checkpoints after it are renumbered.

Suite: 2151 passed, 4 xfailed, 0 failed (+2 in `test_stoichiometry.py`, +5
and the 4 expected failures in `test_user_defined_model.py`); 178 warnings
(+12, all the new test's charge `ConservationWarning`s). Fallback recorder: no
stoichiometry hits. Measurement: 0 values differ in all five cases; warnings
unchanged. Sweep:
`_get_common_species`, `_cs_mod` and the old message text remain only in the
planning docs.

**Before checkpoint 8: what `reactions=` means on a `ControlVolume`.** The CV
already takes `reaction_system=`. To agree before checkpoint 8 starts:
whether `reactions=` only contributes species to `cv.species`, or is an
alternative to `reaction_system=`.

**Before checkpoint 10: how a corrector resolves.** Today `"NaOH"` doses
`Na+`. Under the rule the corrector must come from the model's species. To
agree before checkpoint 10 starts: dose the ion `Species` directly (the model
names `Na+`), or keep a salt `Species` (`NaOH` exists in the bioprocess
database) that states which ion it adds.

## Shipping

- [ ] Full test suite green on the branch
- [ ] `git switch main`
- [ ] `git merge --no-ff explicit-species-resolution -m "Merge explicit-species-resolution: <summary>"`
- [ ] `git tag explicit-species-resolution-shipped`
- [ ] `git push origin main` and `git push origin explicit-species-resolution-shipped`
- [ ] `git branch -d explicit-species-resolution`
- [ ] Move the design note and this checklist to `docs/dev/implementation/shipped/`,
      add a "Shipped" banner to both
- [ ] Update `docs/dev/implementation/upcoming/README.md`: remove the "Design
      discussions" entry and add one to "Recently shipped"
