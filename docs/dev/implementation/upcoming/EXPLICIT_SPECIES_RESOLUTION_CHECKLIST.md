# Phase Kickoff Checklist — explicit-species-resolution

> Checklist for [`EXPLICIT_SPECIES_RESOLUTION.md`](EXPLICIT_SPECIES_RESOLUTION.md),
> the source of truth for the rule, decisions, settled questions and audit
> results; do not restate them here. Where this checklist and the note disagree,
> this checklist wins. See [`README.md`](README.md)'s "Branching and tagging
> convention". Modelled on
> [`../shipped/VANT_HOFF_SINGLE_SOURCE_CHECKLIST.md`](../shipped/VANT_HOFF_SINGLE_SOURCE_CHECKLIST.md).

**Working rules**

- Results are unchanged by rule, with two named exceptions, each in its own
  checkpoint: the yeast molar masses (checkpoint 6) and the Bisection engine's
  write-back (checkpoint 8). Every checkpoint measures and records; an
  unexplained shift stops the checkpoint.
- **Measurement.** A scratchpad script (not committed) builds and runs, from a
  fixed configuration: the stirred-tank template (a batch and a fed-batch tank
  with Monod growth on acetic acid, pH control on), ADM1 and BSM2 (short runs).
  It records every phase's `n_mol` trajectory and pH, and reports per model the
  number of differing values and the maximum absolute and relative shift. The
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
- [ ] This checklist committed on the branch
- [ ] `git status -sb` clean
- [ ] Baseline full suite recorded here (pass / fail / skip counts)
- [ ] Baseline measurement captured (stirred tank, ADM1, BSM2)

## Checkpoints

### String stoichiometry

- [ ] **1. Drop the `common_species` seed.** `_parse_stoichiometry` looks up
      only the caller's `species=`; `_get_common_species` and the `_cs_mod`
      import go. The unknown-id message lists the caller's ids (or says none
      were given) and names the fix. The `species` parameter docs on
      `EquilibriumReaction`, `KineticReaction`, `KspEquilibrium` and
      `_parse_stoichiometry` stop mentioning `common_species`. Callers:
      `tests/standalone/test_stoichiometry.py` (16 sites pass `species=` from
      `AQUEOUS_DEFAULT.species`; `test_unknown_species_message_lists_common`
      becomes a test of the new message);
      `docs/tutorials/D2C_workshop/Example3_CSTR.ipynb` cell 14 (re-run).
      Sanity: suite green; recorder shows no stoichiometry fallback hits;
      measurement unchanged.

### Stirred-tank template and `Chemical`

- [ ] **2. Species for the template's compounds.** In the database modules
      (not `common_species`): Yeast, Yeast_CHO, AceticAcid, PropionicAcid,
      ButyricAcid, CitricAcid, O2 and N2 in the bioprocess database; CH4 and
      H2 in the anaerobic-digestion database. Yeast and Yeast_CHO keep the
      registry's explicit MWs (26.868, 24.626) for now so this checkpoint
      changes no result; the acids and gases compute MW from atoms. Registry
      entries not carried over (solvent alias `Water`, the salts) are listed
      in the notes. Sanity: each new `Species` equals the registry entry's
      atoms, and its MW equals the registry's except CitricAcid (−0.001);
      suite green; measurement unchanged.
- [ ] **3. `aerobic_growth` takes `Species`.** All six participants (substrate,
      biomass, N source, O2, CO2, H2O) are `Species` arguments; the separate
      id / atoms / MW arguments and `species_overrides` go, as does the
      `common_species` import and the local `O2`. Every caller updated
      (the factory, tests, notebooks; listed in the notes). Sanity: the built
      stoichiometry is identical entry by entry for the template's default
      configurations; suite green; measurement unchanged.
- [ ] **4. Template takes the model's chemistry; no default database.**
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
- [ ] **5. Names resolve against the species passed; `compounds.py` goes.**
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
- [ ] **6. Yeast molar masses from atoms.** Yeast and Yeast_CHO drop their
      explicit MW (24.834, 22.593). The tutorials' `MW_yeast = 26.868`
      (`cstr_fermenter.py`, `fed_batch_fermenter.py`) and any other copy of
      26.868 / 24.626 (sweep, including tests) follow. Sanity: the stirred-tank
      measurement moves; the notes record by how much and confirm the shift is
      only in quantities that pass through the yeast MW; ADM1 and BSM2
      unchanged; suite green (tests that pin gram-based yeast values are
      re-baselined with the reason recorded).

### `ControlVolume` species set and warnings

- [ ] **7. `cv.species`.** `ControlVolume` takes `species=` / `reactions=` /
      `chemistry_db=` and builds the read-only `cv.species` in `__init__` from
      those plus the reaction stoichiometries. A new species-level conflict
      check raises `SpeciesConflictError` when the same id arrives with
      different data. `_collect_species_registry` builds from `cv.species`;
      `_common_species_catalog` is deleted. The factory passes the template's
      species through. `Simulation._warn_thermo_mismatch` unchanged. Sanity:
      suite green; the three tests that reached the catalog (audit 3) now
      get those ids from species passed to them; measurement unchanged.
- [ ] **8. Bisection write-back.** The engine writes back only species it was
      given (its declared equilibria and the strong ions present in its
      input); `_CANONICAL_WRITEBACK_SPECIES` is deleted. Sanity: recorder
      shows no zero-filled ions created in `n_mol`; BSM2 sentinels and ADM1
      measured, shifts recorded and explained (e.g. a change in which keys
      the state vector packs); suite green.
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

**Before checkpoint 7: what `reactions=` means on a `ControlVolume`.** The CV
already takes `reaction_system=`. To agree before checkpoint 7 starts:
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
