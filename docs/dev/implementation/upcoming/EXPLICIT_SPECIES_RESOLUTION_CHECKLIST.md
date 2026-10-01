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
  checkpoint: the Bisection engine's write-back (checkpoint 3), the O2 molar
  mass in `monod_aerobic_growth` (checkpoint 6) and the yeast molar masses
  (checkpoint 10). Every checkpoint measures and records; an
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

- [x] **2. No built-in water.** An `EquilibriumSet` has no water until
      `set_water(pKw=...)` is called (`pKw` required); the Bisection solver
      uses Kw = 0 without it, so water plays no part in the charge balance
      unless a water reaction is declared. Models that relied on the built-in
      water (pKw 14, no temperature correction) declare it with the same
      values: ADM1 (`build_adm1_reactions`) and
      `docs/tutorials/D2C_workshop/raw_construction.py`. Sanity: every model
      and notebook that solves without a water reaction found first (suite,
      measurement cases and notebooks under a recorder); measurement
      bit-identical; suite green.
- [x] **3. Bisection write-back.** The engine writes back the species of its
      declared equilibria, including every acid form (HA and A-), plus H+,
      and OH- when water is declared; strong ions are inputs and are not
      written; `_CANONICAL_WRITEBACK_SPECIES` is deleted. Without declared
      equilibria the result reports everything the solve computed. Models
      whose rate laws mean an acid's total sum its exact forms themselves
      (BSM2's VFA uptake, the D2C growth law). Sanity: the two id checks in
      `test_user_defined_model.py` pass and lose their markers (H and O stay
      strict expected failures, reason corrected to solvent water); no
      zero-filled ions in `n_mol`; measurement with acid totals matches
      checkpoint 2 to roundoff; BSM2 sentinels; suite green.

### Stirred-tank template and `Chemical`

- [x] **4. Species for the template's compounds.** In the database modules
      (not `common_species`): Yeast, Yeast_CHO, AceticAcid, PropionicAcid,
      ButyricAcid, CitricAcid, O2 and N2 in the bioprocess database; CH4 and
      H2 in the anaerobic-digestion database. Yeast and Yeast_CHO keep the
      registry's explicit MWs (26.868, 24.626) for now so this checkpoint
      changes no result; the acids and gases compute MW from atoms. Registry
      entries not carried over (solvent alias `Water`, the salts) are listed
      in the notes. Sanity: each new `Species` equals the registry entry's
      atoms, and its MW equals the registry's except CitricAcid (−0.001);
      suite green; measurement unchanged.
- [x] **5. `aerobic_growth` takes `Species`.** `aerobic_growth(substrate,
      biomass, *, o2, co2, h2o, yield_gX_gS, rate_fn, balance, n_source, label)`
      and `monod_aerobic_growth(..., *, o2, co2, h2o, ..., n_source, ...)` take
      every participant as a `Species`; the id / atoms / MW arguments,
      `species_overrides`, the built-in `O2`, the `common_species` import and
      the silent NH3 default for a CHNO balance go (CHNO without `n_source`
      raises). `monod_aerobic_growth`'s O2 term reads `o2.id` but keeps its
      32.0 g/mol until checkpoint 6. Every caller updated (listed in the
      notes). Sanity: stoichiometry unchanged; suite green; measurement and
      the edited notebooks unchanged.
- [x] **6. O2 molar mass in `monod_aerobic_growth`.** The O2 Monod term
      converts O2 to g/L with `o2.MW` (31.998) instead of 32.0.
      `DualSubstrateMonod` loses its O2 defaults: `secondary_id` is required
      (keyword-only) and `secondary_MW` is required when
      `secondary_in_mol_L=False`. Sanity: the shift in models using `Ko2_gL`
      is measured and recorded (the D2C Example1/2 notebooks, ArXiv 03, and
      any test using `Ko2_gL`); suite green.
- [x] **7. Template takes the model's chemistry; no default database.**
      `StirredTankFactory.create_volume` takes `chemistry_db=` and
      `species=`, and `StirredTankBuilder.chemistry()` takes them as
      keywords; the `AD_BASIC` fallback goes. The model's species are the
      database's plus its own, merged by the new shared
      `chemistry.species_check.merge_species` (same id, different data
      raises `SpeciesConflictError`). Partition models come only from a
      database passed; with none and no Henry constant, the existing "No
      partition model" error is raised, now saying no database was passed.
      Growth takes O2, CO2 and H2O from the model's species by id, the
      template's stated convention (see checkpoint 8). `reactions=` is left
      to checkpoint 11. ADM1 passes `AD_BASIC`; BSM2 passes nothing. Every
      other builder and `create_volume` caller passes `AD_BASIC`, what it got
      before. `transfer_species`'s docstring stops citing `_HENRY_PARAMS`
      (OPEN_WORK entry deleted). Sanity: the user-defined route through the
      template works with no database; suite green; measurement and the
      edited tutorials unchanged.
- [ ] **8. The template's gas phase comes from the model's declarations.**
      The template invents no gases: the headspace holds the species the
      model declares (initial composition, transfer, gas feed, reactions).
      Three commits. 8a: `VesselConfig.gas_composition` (`{id: fraction}`,
      empty by default) replaces `yO2_init` / `yCO2_init` / `yN2_init`; the
      builder sets it with `.initial_gas(...)`; `AIR` is an explicit
      composition in the bioprocess database. 8b: the transfer presets name
      their species; nothing transfers unless declared. 8c: the gas feed's
      composition has no default; growth's O2 / CO2 / H2O can be given where
      growth is declared, falling back to those ids in the model's species;
      the core sites that assume these ids (DO sensor, `GasFeed` default,
      gas-liquid link defaults, controller gas tables) get their own design
      note and OPEN_WORK entries. Sanity: callers state what they relied on,
      so the measurement is unchanged; a template with differently named
      gases builds; suite green.
- [ ] **9. Names resolve against the species passed; `compounds.py` goes.**
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
- [ ] **10. Yeast molar masses from atoms.** Yeast and Yeast_CHO drop their
      explicit MW (24.834, 22.593). The tutorials' `MW_yeast = 26.868`
      (`cstr_fermenter.py`, `fed_batch_fermenter.py`) and any other copy of
      26.868 / 24.626 (sweep, including tests) follow. Sanity: the stirred-tank
      measurement moves; the notes record by how much and confirm the shift is
      only in quantities that pass through the yeast MW; ADM1 and BSM2
      unchanged; suite green (tests that pin gram-based yeast values are
      re-baselined with the reason recorded).

### `ControlVolume` species set and warnings

- [ ] **11. `cv.species`.** `ControlVolume` takes `species=` / `reactions=` /
      `chemistry_db=` and builds the read-only `cv.species` in `__init__` from
      those plus the reaction stoichiometries. A new species-level conflict
      check raises `SpeciesConflictError` when the same id arrives with
      different data. `_collect_species_registry` builds from `cv.species`;
      `_common_species_catalog` is deleted. The factory passes the template's
      species through. `Simulation._warn_thermo_mismatch` unchanged. Sanity:
      suite green; the three tests that reached the catalog (audit 3) now
      get those ids from species passed to them; measurement unchanged.
- [ ] **12. Model species.** BSM2 declares `S_cat` / `S_an` as `Species`
      (`atoms={}`, charge ±1) in its species table; ADM1's `_get_species`
      raises on an unknown id instead of building biomass. Sanity: BSM2 and
      ADM1 unchanged; suite green.
- [ ] **13. Correctors from the model's species.** `equilibrate_to_pH` and
      `apply_external_flux` resolve the corrector from the model's species;
      `_STRONG_CORRECTOR_ION` is deleted. The resolution rule is agreed with
      the repo owner before this checkpoint starts (see notes). Callers
      updated: the D2C workshop notebooks and script, the template tutorials'
      `PHController` configuration, `test_cv_advance.py`,
      `test_iron_oxidation.py`, `test_simulation.py`, `conftest.py`. Sanity:
      pH-control results unchanged; suite green.
- [ ] **14. `UnresolvedSpeciesWarning`.** New `UserWarning` subclass in
      `monitoring/`, exported next to `ConservationWarning`. Warns where an
      id enters (CV construction, feeds, dosing, template set-up), naming
      the id and the call, and in the monitor once per id whenever it first
      appears. The notes list every test and notebook that now warns and
      why, including the NR `H2O` and PHREEQC write-backs left in place.
      Sanity: suite green (tests that turn warnings into errors are fixed by
      declaring the species, not by filtering); measurement unchanged.

### Henry and Raoult

- [ ] **15. `Species` objects only.** `HenryEquilibrium` /
      `RaoultEquilibrium` species fields take `Species` or `None`;
      `_resolve_species` and the `"H2O"` defaults go. In the same checkpoint:
      `databases/anaerobic_digestion.py`, ADM1's water link, the two
      engine-basics notebooks and the test sites in audit 2. Sanity: suite
      green; recorder shows no string resolution; measurement unchanged.

### Retire `common_species`

- [ ] **16. Definitions into the databases.** The 26 definitions move into
      the database modules (`aqueous.py`: water, carbonate, ammonia; which
      module owns phosphate, sulfate, sulfide, the spectator ions and the
      metal ions is settled in the notes); each database lists every species
      it offers. Every importer (62 files) is rewritten by script; the module
      and its `chemistry/__init__.py` export are deleted, with the
      `species.py` docstring that cites it. Sanity: repo-wide sweep for
      `common_species` (only history left); every `Species` object identical
      to before by value; suite green; measurement unchanged.
- [ ] **17. Close-out.** Sweep for `_get_common_species`,
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
  so it moves from the `ControlVolume` checkpoints to checkpoint 3, and the
  checkpoints after it are renumbered.

Suite: 2151 passed, 4 xfailed, 0 failed (+2 in `test_stoichiometry.py`, +5
and the 4 expected failures in `test_user_defined_model.py`); 178 warnings
(+12, all the new test's charge `ConservationWarning`s). Fallback recorder: no
stoichiometry hits. Measurement: 0 values differ in all five cases; warnings
unchanged. Sweep:
`_get_common_species`, `_cs_mod` and the old message text remain only in the
planning docs.

**Checkpoint 2 (2026-09-30).** Found while planning the write-back fix: with
no water reaction declared, the Bisection engine still modelled water (a
built-in pKw 14 in `EquilibriumSet.__init__`) and wrote OH- into `n_mol`. On
the repo owner's decision the built-in water goes, and this checkpoint comes
before the write-back fix so nothing moves and then moves back.

Who relied on it, found before editing: a recorder (scratchpad
`water_plugin.py`) logged every engine built by `from_reactions` with no water
reaction and later solved, over the full suite, the five measurement cases and
every tutorial and validation notebook or script that builds chemistry (each
run from the scratchpad, 5 min cap). Result: ADM1 (`base.py`, via the three
ADM1 tests in `test_simulation.py`), `raw_construction.py`, and four tests in
`test_speciation.py` (`TestCanonicalEmission` ×3, `TestApplyToPhasesWriteback`).
Every notebook ran to completion except `aerobic_fermentation_stoichiometry.ipynb`
(60 h cells, over the cap); its source takes its equilibria from
`BIOPROCESS_BASIC.reactions`, which include water. `batch_fermenter.ipynb`
stops at its known unparseable cell 10. The notebook generators were not run
(they rewrite repo files). ADM1's engine has no acid-base equilibria at all, so
without water its charge balance would have nothing to balance the strong ions
it seeds.

- `equilibria.py`: `_water` starts as `None`; `water` is `Optional[WaterDef]`;
  `set_water`'s `pKw` is required (it defaulted to 14 when omitted); `__repr__`
  and `summary()` say when no water is set. The class docstring's "or load a
  preset" went (the presets were deleted earlier).
- `acid_base.py` `solve_from_equilibrium_set`: Kw = 0 when the set has no water.
- `engine.py` `from_reactions` docstring: without a water reaction, water plays
  no part in the charge balance and OH- is zero.
- ADM1 `build_adm1_reactions` and `raw_construction.py` declare
  `H2O <-> H+ + OH-` with log_K -14 and no enthalpy, i.e. the built-in values.
  `raw_construction.py` gains a `make_water_dissociation()` factory.
- The four `test_speciation.py` tests stay as they are: they check which ids
  are emitted, not values, and now exercise the no-water path. New
  `TestWaterDeclaration` (4 tests, species defined in the test): a new set has
  no water, `set_water` needs `pKw`, no water reaction gives OH- = 0 with
  H+ = A-, and a declared one gives [H+][OH-] = 1e-14.

Measurement: 0 values differ in all five cases. ADM1 raises one more
`ConservationWarning`, a charge residual: declaring water puts OH- in the
conservation monitor's registry, while the Na+ and Cl- ADM1 seeds after the CV
is built are still invisible to it (fixed by `cv.species`, checkpoint 11). The
monitor now sees part of the charge picture where it saw none; the physics is
unchanged. Suite: 2155 passed, 4 xfailed, 0 failed (+4); 181 warnings (+3,
that charge residual in the three ADM1 tests). `test_user_defined_model.py`
unchanged (5 passed, 4 xfailed): it declares water.

Found, not changed: the engine recognises the solvent by a hard-coded id
(`_SOLVENT_IDS = ("H2O",)` in `engine.py`), another name-keyed backend table.
To settle with the write-back checkpoint or log.

**Checkpoint 3 (2026-09-30).** Decided with the repo owner along the way:

- An engine built without declared equilibria (only tests call it that way,
  with totals as keyword arguments) reports every species it computed.
- An acid's forms are separate species in `n_mol`; a model that needs a total
  computes it where it needs it. Found by measuring: BSM2 and the D2C script
  declare `HA <-> A- + H+` but their rate laws read the HA id as the whole
  acid, which only worked because the fixed list never wrote A- back (BSM2's
  comment said "total tracked under the protonated species id"). Two
  alternatives were considered and set aside: rate laws summing HA + A-
  through a helper that rebuilt the conjugate id from a naming pattern, and a
  per-reaction "keep my total under HA" option on `EquilibriumReaction`. A
  kinetic reaction consuming one form is made up from the others at the next
  solve; if a step removes more of that form than is present, the step solver
  clamps and warns. That is general to fast equilibria with stepped kinetics
  and is documented, not special-cased.

Changes:

- `engine.py`: `_CANONICAL_WRITEBACK_SPECIES` deleted. New `_written_species`
  (the declared equilibria's `species_refs`, H+ when anything is declared,
  OH- when water is); `species_mol_L` holds those, everything else goes to
  `extra`; without an `EquilibriumSet`, everything computed. `algebraic_species()`
  returns the same set, so OH- is in it only when water is declared.
  `from_reactions`' docstring gains a "What is written back" paragraph
  (totals are the model's to compute; clamping warns; a smaller step or
  `SimultaneousAdaptiveSolver`, which re-solves the equilibria at every
  right-hand-side evaluation by default, avoids it).
- `EquilibriumReaction` docstring: each declared form is its own amount.
  `protocols.py` `EquilibriumResult` docstring and the `get_CO2aq_from_totals`
  comment updated.
- BSM2: `rate_uptake` / `rate_c4_competitive` take the exact ids summed as the
  substrate, listed at each uptake reaction (`("S_ac", "S_ac-")`,
  `("S_su",)`, ...); the VFA table names each conjugate base explicitly
  (`_BSM2_VFA`, replacing `_BSM2_PKA_VFA`). D2C `raw_construction.py`'s growth
  law sums `"AceticAcid"` and `"Acetate-"`.
- Tests: `test_user_defined_model.py` loses the two markers that now pass; the
  H/O reason is corrected (solvent water: +2.5e-8 relative over 1 h, exactly
  one O and two H per HCO3- formed). `test_bsm2_reference.py` sums
  `("S_ac", "S_ac-")` and `("S_pro", "S_pro-")` for those sentinels
  (`SENTINEL_LIQUID_IDS`); the sentinel values are unchanged and all six pass.
  Stale comments in `test_speciation.py` and `test_nr_speciation_engine.py`
  updated.
- `01_bisection_engine_basics.ipynb`: three markdown cells stop describing the
  fixed list (exact replacements). It runs; its saved outputs in cells 5 and 11
  still print the old `species_mol_L` with zero-filled ions (not re-saved).

Measurement against checkpoint 2, with each acid's total (HA + A-) in place of
HA: BSM2 821 of 3600 values differ, max 7.8e-11 relative (summation order);
D2C 1133 of 3417, max 4.4e-13; ADM1 and both stirred tanks 0. The zero-filled
ions (Ca++, Co++, K+, Mg++, Mn++, Mo7O24------, Zn++, and Na+/Cl- where not
seeded) no longer appear; BSM2 and D2C gain their conjugate-base series.
BSM2 raises one fewer `ConservationWarning` (9). Clamp warnings: no new ones,
and none on an acid form; the two already raised are unchanged (O2 in the
stirred tanks, H2O in BSM2). Coverage checked: the sequential and simultaneous-Euler step
solvers warn on every `clamp_fn` use; the monolithic solver does not floor;
`Phase.apply_flux`'s default floor and `SimultaneousAdaptiveSolver`'s trial-state
floor are silent (logged).

Suite: 2157 passed, 2 xfailed, 0 failed; 166 warnings (-15: the baseline
test's charge warnings and one BSM2 `ConservationWarning`). Logged in
`OPEN_WORK.md`: solvent water never debited or credited, with the engine's
hard-coded `"H2O"` / `"H+"` / `"OH-"` ids (settles the checkpoint 2 finding);
BSM2's nitrogen inhibition reading molecular NH3 as total nitrogen; the silent
floors.

**Checkpoint 4 (2026-09-30).** `databases/bioprocess_basic.py` defines and
lists `AceticAcid`, `PropionicAcid`, `ButyricAcid`, `CitricAcid`, `Yeast`,
`Yeast_CHO`, `O2` and `N2`; `databases/anaerobic_digestion.py` defines and lists
`CH4` and `H2` (and now imports `Species`). Module docstrings updated. MW from
atoms except `Yeast` (26.868) and `Yeast_CHO` (24.626), whose comment says
they are the template's long-standing values, above the formula weights, with
no recorded source.

Against `ChemicalRegistry.default()`: every atom composition equal; MW exactly
equal except `PropionicAcid` (+1.4e-14, last bit: the registry's 74.079 vs the
computed 74.07900000000001) and `CitricAcid` (192.123 vs 192.124, as settled).
Both will reach the template at checkpoint 9; the propionic one shows up only
at roundoff in `test_builder.py`'s two-substrate case.

Registry ids with no `Species` in any database: the `Water` alias (the
database uses `H2O`) and 21 salts, acids and bases (`(NH4)2SO4`,
`AmmoniumMolybdate`, `AmmoniumSulfate`, `CaCl2`, `CaSO4`, `CoCl2`, `H2SO4`,
`HCl`, `KCl`, `KH2PO4`, `KOH`, `MgSO4`, `MnCl2`, `NH4Cl`, `Na2HPO4`, `Na2SO4`,
`NaCl`, `NaH2PO4`, `NaHCO3`, `NaOH`, `ZnSO4`). No repo caller reaches any of
them through the template: its registry lookups are `Yeast`, `AceticAcid`,
`PropionicAcid` and the default N source `NH3` (swept across `.py`, `.ipynb`,
`.md`); `E_coli` and `Glucose` always come with explicit atoms and MW.
`NH4Cl`, `KH2PO4` and `NaOH` are defined in `bioprocess_basic.py` but not listed
in its species; left for checkpoint 16's "every species it offers".

`docs/tutorials/reactions/chemistry_database.ipynb` prints the extended
database's species count: its saved output says 21, a run now gives 30 (not
re-saved). Suite: 2157 passed, 2 xfailed, 0 failed. Measurement: 0 values
differ in all five cases.

**Checkpoint 5 (2026-09-30).** Decided with the repo owner: the O2 term's
hard-coded 32.0 g/mol becomes `o2.MW` in its own checkpoint (6), so this one
stays a pure refactor.

- `builder.py`: new signatures as in the item; ids, atoms and MW come from the
  `Species` passed, which appear in the reaction as the same objects; the
  element-balance check still rejects an `o2` / `co2` / `h2o` with other atoms.
  Module docstring examples and `monod_aerobic_growth`'s docstring updated.
- `factory.py`: `_build_reaction_system` takes the `chemistry_db` that
  `create_volume` holds and takes O2, CO2 and H2O from its species by id
  (raising, with the available ids, if one is missing); in `AD_BASIC` those
  are the same CO2 / H2O objects the builder imported before and an O2 equal
  by value. Organism and substrate still come from the registry (until
  checkpoint 9), wrapped as `Species` with the registry's atoms and MW. An N
  source the registry does not know now raises `KeyError` instead of warning
  and silently using NH3's atoms; no repo caller hits it. The unused
  `warnings` import went.
- `raw_construction.py`: declares `O2` next to its `CO2` (whose comment about
  matching the builder's internal CO2 was no longer true) and passes both.
- Tests: `test_reactions.py`'s builder tests pass `Species` defined in the
  test (`_gases()`, `_nh3()`, `_acid()`); two new tests: the reaction holds the
  very objects passed, and CHNO without `n_source` raises.
  `test_rate_laws.py`'s Monod comparison passes O2 / CO2 / H2O.
- Notebooks (source only; parsed, edited and written back in each file's own
  format, with a byte-level fallback for ArXiv 03, whose formatting json
  cannot reproduce): ArXiv `03_cstr_dilution_rate_sweep` and its generator
  (same text; generator not run) declare `O2` and pass `CO2`/`H2O`/`NH3` from
  their `common_species` import; D2C Example1/2 add `H2O, CO2, NH3` to their
  `common_species` import and `O2` to their `bioprocess_basic` import; D2C
  Example3 declares `O2` next to its water and carbonate species;
  `reactions/reaction_system.ipynb` declares `O2` and `H2O`. A first attempt
  split Example2's single-string cell sources into one element per line; it
  was restored from `HEAD` (and the four D2C / reactions notebooks' CRLF
  working-tree endings restored) before the edit was redone keeping each
  source's structure.

Run from the scratchpad against the `HEAD` export and the working tree: D2C
Example1 (159 lines), Example2 (130), Example3 (35) and `reaction_system`
(26) print identical output; ArXiv 03 differs only in its "Run Time, ms"
column (wall clock). ArXiv 03 saves a figure under
`docs/tutorials/**/figures/`, which is gitignored.

Suite: 2159 passed, 2 xfailed, 0 failed (+2). Measurement: 0 values differ in
all five cases.

**Checkpoint 6 (2026-09-30).** Decided with the repo owner: `DualSubstrateMonod`'s
own O2 defaults (`secondary_id="O2"`, `secondary_MW=32.0`) go in the same
checkpoint.

- `builder.py`: `monod_aerobic_growth` passes `secondary_MW=float(o2.MW)`;
  docstring says `o2.MW`.
- `rate_laws.py`: `secondary_id` is `field(kw_only=True)` with no default;
  `secondary_MW` defaults to `None`, and `__post_init__` raises when
  `secondary_in_mol_L=False` without it; `make_rate_fn` only converts when it
  is used. Docstring updated. Its only non-test caller (`monod_aerobic_growth`)
  already passed both.
- `test_rate_laws.py`: the frozen reference `_reference_monod_rate_fn` takes
  `MW_O2` explicitly (no 32.0 default); the O2-path fingerprint passes the O2
  `Species`' MW and stays bit-identical on its 500-point grid; the zero-`Ko2`
  edge test passes `MW_O2` too. Two new tests: `secondary_id` is required, and
  g/L conversion without `secondary_MW` raises.
- `OPEN_WORK.md`: the "Mapping-based parameters" entry no longer says the O2
  id is fixed.

Measured shift, notebooks run against the `HEAD` export (no measurement case
uses `Ko2_gL`; all five are unchanged): ArXiv 03 steady-state substrate
0.5385 -> 0.5386 g/L and, at D = 0.24/h, 2.4343 -> 2.4344 g/L (its other
differing lines are the wall-clock "Run Time, ms" column); D2C Example1
productivity 0.0734 -> 0.0733 g/L/h; D2C Example2 final biomass 11.091 ->
11.092 g/L. Every other printed line identical. That is the size expected
from moving the O2 term's g/L by 6e-5 relative. Saved outputs in those
notebooks are not re-saved.

Suite: 2161 passed, 2 xfailed, 0 failed (+2).

**Checkpoint 7 (2026-09-30).** Decided with the repo owner: the builder takes
the chemistry through `.chemistry(activity_model, *, chemistry_db=None,
species=None)`; the template's `reactions=` waits for checkpoint 11; the
species-level conflict check is a new shared function now. Asked why growth
looks up only O2, CO2 and H2O, and by fixed ids: those are the growth
reaction's only participants besides substrate, biomass and N source (N2 is
the inert headspace gas), and the fixed ids run through the whole template
(gas phase, vessel fractions, default transfer, gas feed), so they stay the
template's stated convention for now and checkpoint 8 makes them explicit
end to end.

- `species_check.py`: `merge_species(*sources)` merges `{id: Species}`
  mappings and iterables of `Species`, keeps the first of equal ones, raises
  `SpeciesConflictError` on the same id with different atoms, charge or MW,
  and `ValueError` when a mapping keys a `Species` under another id. Exported
  from `PyOMES.chemistry`; 6 tests in `test_species_check.py`.
- `factory.py`: `create_volume(..., chemistry_db=None, species=None)`; no
  `AD_BASIC` import; module docstring states the chemistry comes from the
  model and the fixed gas ids, and its example passes `AD_BASIC`.
  `_build_reaction_system` takes the merged species; a missing O2 / CO2 / H2O
  raises naming the id, what is available and the fix.
- `builder.py`: `.chemistry()` stores `chemistry_db` / `species` only when
  given, so a later `.chemistry(activity_model=...)` keeps them; `build()`
  passes them on; module example passes `AD_BASIC`.
- ADM1 `build_adm1_cv` passes `AD_BASIC` (it relies on its CH4 / H2 / H2S
  Henry models). BSM2 unchanged.
- Tests: `test_builder.py` and `test_simulation.py` pass `AD_BASIC` to every
  builder (39) and `create_volume` (4) they use (scripted, exact counts).
  New: `TestModelChemistry` in `test_builder.py` (no database and no Henry
  constant raises; a species conflicting with the database raises; growth
  without O2 names the fix; a later `.chemistry()` keeps the database) and,
  in `test_user_defined_model.py`, a stirred tank built from the test's own
  O2 / CO2 / H2O / N2 and explicit Henry constants with no database, whose
  growth reaction holds those objects and which runs and grows.
- Tutorials and docs pass `AD_BASIC`: the four template scripts,
  `batch_fermenter.ipynb` (whose builder table also said `.chemistry()` was a
  "pH-active equilibria preset", which it never was),
  `aerobic_fermentation_stoichiometry.ipynb`, `results/01_exporting_results.ipynb`,
  `README.md` and `docs/tutorials/templates/README.md`; `docs/architecture.md`
  says the builder is given the model's chemistry. `rate_laws.py`'s docstring
  fragment builds nothing and is unchanged.

The measurement script's stirred-tank builder passes `AD_BASIC` (when the
installed builder accepts it). Run against the `HEAD` export: the four
template scripts and `01_exporting_results` print identical output;
`batch_fermenter.ipynb` stops at its known unparseable cell 10 in both, with
identical output before it; `aerobic_fermentation_stoichiometry.ipynb` not run
(60 h cells, over the CPU cap).

Suite: 2173 passed, 2 xfailed, 0 failed (+12). Measurement: 0 values differ
in all five cases.

**Checkpoint 8 (2026-10-01).** Decided with the repo owner: no role keywords
for the template's gases (brittle for a general template); the gas phase is
inferred from what the model already declares; the fixed fractions become a
composition keyed by species id; the core sites go to a new design note.

8a:
- `configs.py`: `VesselConfig.gas_composition: Dict[str, float]` (empty by
  default, normalised by its sum, negative fractions raise) replaces the three
  fractions and the N2-balance rule.
- `factory.py`: the initial gas phase is the composition's species, plus
  transfer species at zero as before; module docstring and example updated.
  The sum is taken in the composition's order, so a composition listing O2,
  CO2, N2 reproduces the old arithmetic bit for bit.
- `builder.py`: `.vessel()` loses the fractions; new `.initial_gas(composition)`
  (separate because `.vessel()` replaces all its settings on each call);
  module example uses `AIR`.
- `bioprocess_basic.py`: `AIR = {"O2": 0.2095, "CO2": 0.0004,
  "N2": 1.0 - 0.2095 - 0.0004}`, the old default written out (N2 computed as
  the old balance was, so bit-identical).
- Callers (scripted, exact counts) state what they relied on: ADM1 and BSM2
  `{"N2": 1.0}` (their O2 and CO2 fractions were zero); `microplate_fermenter.py`
  its own fractions; the other template scripts, three notebooks, both
  READMEs and the tests `AIR`; the user-defined template test its own
  `{"O2": 0.21, "N2": 0.79}`. `test_configs.py`'s two `yN2` tests became three
  composition tests. `batch_fermenter.ipynb`'s builder table gains the
  `.initial_gas` row. The notebook-editing helpers moved to a scratchpad
  module (`nbedit.py`) shared by later scripts.

Measurement against checkpoint 7: 0 values differ in all five cases; BSM2 no
longer carries a zero `gas:O2` entry it never declared (ADM1 keeps its O2
entry: its transfer config still declares O2). Tutorials against `HEAD`: the
four template scripts and `01_exporting_results` identical;
`batch_fermenter.ipynb` stops at its known cell 10 in both. Suite: 2174
passed, 2 xfailed, 0 failed (+1).

**Before checkpoint 11: what `reactions=` means on a `ControlVolume` and on the
stirred-tank template.** Both already take `reaction_system=` (the builder
through `.reaction_system()`); the template's `reactions=` was left to this
decision at checkpoint 7 so the two mean the same. To agree before
checkpoint 11 starts:
whether `reactions=` only contributes species to `cv.species`, or is an
alternative to `reaction_system=`.

**Before checkpoint 13: how a corrector resolves.** Today `"NaOH"` doses
`Na+`. Under the rule the corrector must come from the model's species. To
agree before checkpoint 13 starts: dose the ion `Species` directly (the model
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
