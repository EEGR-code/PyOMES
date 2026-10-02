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
- [x] **8. The template's gas phase comes from the model's declarations.**
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
- [x] **9. Names resolve against the species passed; `compounds.py` goes.**
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
- [x] **10. Yeast molar masses from atoms.** Yeast and Yeast_CHO drop their
      explicit MW (24.834, 22.593). The tutorials' `MW_yeast = 26.868`
      (`cstr_fermenter.py`, `fed_batch_fermenter.py`) and any other copy of
      26.868 / 24.626 (sweep, including tests) follow. Sanity: the stirred-tank
      measurement moves; the notes record by how much and confirm the shift is
      only in quantities that pass through the yeast MW; ADM1 and BSM2
      unchanged; suite green (tests that pin gram-based yeast values are
      re-baselined with the reason recorded).

### `ControlVolume` species set and warnings

- [x] **11. `cv.species`.** `ControlVolume` takes `species=` / `reactions=` /
      `chemistry_db=` and builds the read-only `cv.species` in `__init__` from
      those plus the reaction stoichiometries. A new species-level conflict
      check raises `SpeciesConflictError` when the same id arrives with
      different data. `_collect_species_registry` builds from `cv.species`;
      `_common_species_catalog` is deleted. The factory passes the template's
      species through. `Simulation._warn_thermo_mismatch` unchanged. Sanity:
      suite green; the three tests that reached the catalog (audit 3) now
      get those ids from species passed to them; measurement unchanged.
- [x] **12. Model species.** BSM2 declares `S_cat` / `S_an` as `Species`
      (`atoms={}`, charge ±1) in its species table; ADM1's `_get_species`
      raises on an unknown id instead of building biomass. Sanity: BSM2 and
      ADM1 unchanged; suite green.
- [x] **13. Correctors from the model's species.** `equilibrate_to_pH` and
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

8b:
- `configs.py`: `TransferConfig.kinetic(kLa={id: kLa_per_h}, equilibrium=(ids))`
  and `TransferConfig.equilibrium(ids)` replace `default_kinetic(kLa_O2,
  kLa_CO2_ratio)` / `default_equilibrium()`, which built O2 / CO2 / N2 entries;
  `TransferConfig()` (no transfer) is the default. `kLa_CO2_ratio` went from
  `TransferConfig` and its `to_dict` / `from_dict`: nothing read it. The
  `SpeciesTransferConfig` docstring said a missing Henry constant came from "the
  temperature-dependent correlation"; it comes from the `chemistry_db`'s
  partition model.
- `builder.py`: `.transfer_kinetic(kLa, equilibrium=())` and
  `.transfer_equilibrium(species_ids)` take the species; with no transfer
  declared, `_build_configs` and `.transfer_species` start from no transfer
  (they started from O2 / CO2 / N2 equilibrium). Module example updated.
- `factory.py`: docstring says only the species in ``transfer`` transfer;
  example uses `TransferConfig.kinetic`.
- Callers (scripted; old preset calls found by pattern and rewritten one for
  one, CO2's kLa written as the product the preset computed): tests, ADM1, the
  template scripts, four notebooks (including `aerobic_fermentation_stoichiometry`,
  which passes variables), both READMEs. Builder chains in the tests that relied
  on the silent O2 / CO2 / N2 equilibrium default got
  `.transfer_equilibrium(["O2", "CO2", "N2"])` after their `.initial_gas(AIR)`
  prefix (a later transfer call replaces it, `.transfer_species` adds to it, as
  before). In the tutorials and READMEs the products are shown as their exact
  literals (150.0 * 0.9 == 135.0, 90.0 * 1.0 == 90.0, 120.0 * 0.9 == 108.0).
  Test rewrites: "default transfer is equilibrium" became "no transfer unless
  declared"; `TransferConfig`'s preset tests became tests of `kinetic` /
  `equilibrium` / the empty default; the no-database test declares O2
  transfer; the empty-builder repr test uses a bare builder again (the
  checkpoint 7 script had prefixed it); `batch_fermenter.ipynb`'s builder table
  row updated.

The core kLa ratio (`DOAgitationController.kLa_CO2_ratio`,
`KineticGasLiquidLink.set_kLa_with_co2_ratio`) is for the core-gas-ids note.
Measurement against 8a: 0 values differ in all five cases. Tutorials against
`HEAD`: template scripts and `01_exporting_results` identical; `batch_fermenter.ipynb`
stops at its known cell 10 in both. Suite: 2174 passed, 2 xfailed, 0 failed.

8c:
- `configs.py`: `GasFeedConfig.composition` has no default (it was
  `{"O2": 0.21, "N2": 0.79}`); a feed with `vvm_min > 0` and no composition
  raises. `OrganismConfig` gains `o2_id` / `co2_id` / `h2o_id` (default
  `"O2"` / `"CO2"` / `"H2O"`, resolved against the model's species), and
  `resolve()` now returns `dataclasses.replace(self, atoms=..., MW=...)` so
  every field carries over.
- `builder.py`: `.gas_feed()` docstring says the composition has no default;
  `.organism()` takes the three ids.
- `factory.py`: growth looks up the organism's ids; the miss message names the
  role, the id and the fix (`'O2' is not among the species passed` is kept,
  a test matches it); docstrings and the example's gas feed updated.
- Tests: `GasFeedConfig` default test passes a composition, plus "no
  composition raises" and "zero vvm needs none"; `resolve()` keeps the gas ids;
  one builder test's feed states its composition; and in
  `test_user_defined_model.py` a stirred tank whose gases are `O2_aq` /
  `CO2_aq` / `N2_g` / `H2O_l` builds, its gas phase holds exactly the declared
  gases, and its growth reaction uses those objects with no `O2` / `CO2` /
  `H2O` / `N2` anywhere.
- New design note `upcoming/GAS_SPECIES_IN_CORE.md` (with an `upcoming/README.md`
  entry and an `OPEN_WORK.md` pointer): the core sites that still assume the
  ids (DO sensor in `snapshot.py`, `GasFeed`'s and `MembraneGasBoundary`'s
  default air, the link's `{"CO2": "CO2"}` alpha default and
  `set_kLa_with_co2_ratio`, the DO controllers' kLa paths and `kLa_CO2_ratio`,
  the vent physics' gas tables, which give another id gamma 1.35 and air's
  molar mass unless `chemicals` is passed). Four open questions.

The pre-checkpoint suite for 8c ran while the first edits were being made, so
it is not a clean baseline; `HEAD` was unchanged since 8b's final run (2174
passed), which stands as the baseline. Measurement against 8b: 0 values differ
in all five cases. Suite: 2179 passed, 2 xfailed, 0 failed (+5).

9:

Agreed before starting (2026-10-01), refining the item above: `atoms=` stays.
An organism or substrate is
- an id alone, which must be among the model's species;
- an id with `atoms=` (MW computed from `units.ATOMIC_WEIGHTS` unless `MW=`
  is given), or a `Species`, which joins the model's species.

A definition equal to a model species of the same id (atoms, charge, MW)
uses that species; a different one raises `SpeciesConflictError` unless
`overwrite=True`, which replaces it in the model's species. `MW=` without
`atoms=` and a `Species` with `atoms=` / `MW=` raise. The fields are renamed
to `organism` / `substrate` (positional calls unchanged). `n_source_id` keeps
its name and shape (option b): it resolves against the model's species, has
no default, and is required for CHNO; generalising it to other elements
(sulfur etc.) is the new design note `upcoming/GROWTH_STOICHIOMETRY.md`.

- `configs.py`: `OrganismConfig.organism` / `SubstrateConfig.substrate`
  (str or `Species`, required), `overwrite`, `n_source_id=None` (raises for
  CHNO without one); `resolve()` gone from both; `_check_definition` validates
  the combinations; `to_dict()` keeps a `Species` as the object;
  `ChemistryConfig.acid_pKas` gone.
- `factory.py`: `_resolve_definition` applies the rules above against the
  merged model species (new definitions are added to it);
  `_not_among_species` gives the shared "is not among the species passed"
  wording (organism, substrate, growth gases, N source); the N source is
  the model's `Species`, resolved once rather than per substrate.
- `builder.py`: `.organism(organism, ...)` / `.substrate(substrate, ...)` with
  `overwrite=`, no default ids; `__repr__` shows a `Species`' id.
- Deleted: `PyOMES/compounds.py`, the exports in `PyOMES/__init__.py`,
  `tests/standalone/test_compounds.py` (19 tests), the `registry` fixture in
  `conftest.py` (used only there), the `chemistry/__init__.py` docstring
  sentence about it.
- Tests: `test_configs.py` organism / substrate tests rewritten for the new
  fields and checks (62 -> 65); `test_builder.py` CHNO test names
  `n_source_id="NH3"`, two calls pass an id; `test_user_defined_model.py`
  gains `TestStirredTankOrganismAndSubstrate` (9 tests, user species only):
  id resolves to the model's object, unknown id names the fix, id + atoms
  with computed MW, `Species` used as given, equal definition reuses the
  model's, different one raises, `overwrite=True` replaces, CHNO takes the
  model's NH3, unknown N source names the fix.
- Docs: `PyOMES/README.md` row, `docs/architecture.md` layout line,
  `README.md` test table, the Yeast_CHO comments in
  `D2C_workshop/raw_construction.py` and `reactions/reaction_system.ipynb`
  (now cite `PyOMES.databases.bioprocess_basic`), the OPEN_WORK weighed-salt
  entry. OPEN_WORK gains pointers for the growth-stoichiometry note and for
  `batch_fermenter.ipynb` cell 10, which does not compile (also on `main`).
- Sweep: no `ChemicalRegistry`, `PyOMES.compounds`, `compounds.py` or
  template `acid_pKas` outside history docs (`shipped/`, `obsolete/`, the
  strong-ion note's record of a past decision); `organism_id=` /
  `substrate_id=` remain only as rate-law `make_rate_fn` keywords.

Sanity: growth stoichiometry printed from a `git archive HEAD` export and
the working tree for CHO, CHNO, two substrates and an explicit E_coli /
glucose definition is identical except PropionicAcid, whose MW is now the
database's computed 74.07900000000001 (registry literal 74.079), which moves
its growth coefficients in the last digit. Measurement against 8c: 0 values
differ in all five cases. Tutorial scripts (four template scripts,
`raw_construction.py`): output identical to HEAD. `batch_fermenter.ipynb`:
cells 0-9 print the same on both trees; cell 10 fails on both (above).
Registry use: none possible, the module is gone and the suite imports
cleanly. Suite: 2179 passed, 2 xfailed before; 2172 passed, 2 xfailed, 0
failed after (-19 test_compounds, +3 test_configs, +9
test_user_defined_model).

10:

Agreed before starting (2026-10-01): hand-written yeast `Species` in tutorials
and tests that copied 24.626 / 26.868 for CH1.61O0.56(N0.16) drop `MW=` too.
The ArXiv `03_cstr` notebook and its generator keep `Ecoli` CH1.8O0.5N0.2 at
24.626: that is its formula weight (24.6263), not a copy.

- `bioprocess_basic.py`: `Yeast` and `Yeast_CHO` lose `MW=` (now 24.834 and
  22.593); the comment says MW is the formula weight. `AD_BASIC` takes them
  from there.
- Tutorials: `cstr_fermenter.py` / `fed_batch_fermenter.py` read
  `MW_yeast = AD_BASIC.species["Yeast"].MW` instead of a literal;
  `raw_construction.py` and `reaction_system.ipynb` drop `MW=` from their
  Yeast (and the inoculum comment quotes 22.593); D2C `Example3_CSTR.ipynb`
  drops it from PEKILO and its markdown quotes 22.593.
- Tests: `test_reactions.py` (7 places; the two hand-worked parity tests take
  `MW_X` from the biomass Species they build), `test_builder.py` (one
  `.organism(..., atoms=)`), `test_user_defined_model.py` (`YEAST_MW` is the
  computed value). No test pinned a gram-based yeast value; none needed
  re-baselining.
- Sweep: no 26.868 / 24.626 outside `docs/dev/` except the ArXiv `Ecoli`.

Sanity: suite 2172 passed, 2 xfailed before and after. Measurement against 9:
ADM1 and BSM2 0 values differ. st_batch / st_fedbatch: 897 of 1010 differ
(max rel 18 % mid-run, gas CO2); at the end AceticAcid +6.7 % / +4.2 %, CO2
and H2O about -12 %, Yeast (mol) +0.01 %. The fingerprint fixes the inoculum
in mol, so it holds 7.6 % fewer grams and Monod (g/L) grows more slowly, while
each mol of substrate now makes 8.2 % more mol biomass. With Yeast patched
back to 26.868 on the new code, both cases match 9 exactly (0 of 1010), so the
shift is only the yeast MW. d2c_raw: Yeast +9.0 % at the end (inoculum given
in g/L: 24.626 / 22.593), CO2 -7.0 %. Tutorials against a HEAD export: the
four template scripts and `Example3_CSTR.ipynb` print identical output (they
report in grams, where the MW cancels); `raw_construction.py` moves as d2c_raw;
`reaction_system.ipynb`'s growth reaction goes from
`-1.01 O2, +0.878 Yeast, +1.12 CO2, +1.29 H2O` to
`-0.926 O2, +0.957 Yeast, +1.04 CO2, +1.23 H2O` (its saved outputs still show
the old line; outputs are not re-saved here).

11:

Agreed before starting (2026-10-01), settling the question below: `reactions=`
means the model's reactions, which run. On a `ControlVolume` it is shorthand
for `reaction_system=ReactionSystem([...])` and passing both raises; on the
template (when it gains the argument) they run alongside the growth reactions
it builds, and `.reaction_system()` still replaces everything. A database's
reactions are not run: it is a library to pick from.

- `control_volume.py`: `reactions=` and `species=` (keyword, last in the
  signature); `cv.species`, a read-only view of a dict built once in
  `__init__` with `merge_species(chemistry_db.species, species=,
  stoichiometry species)` (a dict, because a `mappingproxy` does not pickle
  and the HPC checkpointing tests pickle CVs); `_collect_species_registry`
  returns `dict(cv.species)`; `_stoichiometry_species` walks the reactions;
  `_common_species_catalog` deleted; `snapshot()` passes `species=`;
  `merge_species` / `check_species_consistency` / `Species` imported at module
  level (`chemistry/` imports nothing from `core/`).
- `factory.py`: passes the template's merged species (`species=model_species`);
  not `chemistry_db`, so `_warn_thermo_mismatch` behaves as before.
- ADM1 (`models/vlmodels/adm1/base.py`): the new conflict check found its
  `H2S` (MW 34.08) against AD_BASIC's (34.080999..., computed). Agreed: ADM1
  uses the shared `H2S` object as it does CO2 / NH3 / H2O, and its table's H2S
  MW (read by `_mw()` in the H2S inhibition term) is that Species' MW. No
  other ADM1 species differs from AD_BASIC's.
- Tests: `test_cv_advance.py` `TestEquilibrateToPH` and
  `tests/validation/speciation/test_iron_oxidation.py` pass their spectator
  ions (Cl- / K+ / Na+, and K+) as `species=`; `test_accuracy_monitor.py:562`
  never advances, so it needs nothing. `test_user_defined_model.py` gains
  `TestModelSpecies` (8): reaction and passed species, read-only, they feed
  the monitor, conflict raises, `reactions=` runs like a `ReactionSystem`,
  both together raise, `snapshot()` keeps them, the template's CV holds the
  model's species.
- Docs: the registry descriptions in `monitoring/conservation.py`,
  `reactions/reaction_system.py` and `docs/solvers.md`. OPEN_WORK:
  `snapshot()` drops `chemistry_db`; the monitor does not net out boundary
  flows.

Sanity: suite 2172 passed, 2 xfailed before; 2180 passed, 2 xfailed, 0 failed
after (+8). Measurement against 10: 0 values differ in all five cases. Warning
counts are the same except ADM1, 3 -> 4 ConservationWarnings: its registry now
holds Na+, Cl-, N2 and O2 (from AD_BASIC). Na+ / Cl- are seeded after the CV
is built, so the old catalog scan never saw them and the monitor reported a
0.112 mol charge residual, which is (0.100 - 0.030) mol/L x 1.6 L of uncounted
Na+ / Cl-: that false warning is gone. N2 is new to the accounting, and the
N2 the vent releases in one step (1.911e-4 mol N) now shows as N drift
(per-step and cumulative), the same way vented CO2 and water already show as
O and H drift.

12:
- `bsm2.py`: `SPECIES` gains `S_cat` / `S_an` (`atoms={}`, charge +1 / -1,
  so MW 0); `build_bsm2_cv` passes `species=SPECIES` to `.chemistry()`, so
  `cv.species` holds the whole table (BSM2 passed no species before, and the
  CV knew only the species in its reactions). Nothing iterates `SPECIES`, so
  the two lumps reach no ThOD or MW sum.
- ADM1 `base.py`: `_get_species` builds CHO biomass only for organisms in
  `ORG` (CHON for `CHON_ORGS`, as before) and raises `KeyError` for any other
  id not in `SPECIES`; docstring says so. Its `_mw()`, and BSM2's `_mw()` /
  `_thod()` / `_atoms()`, keep their biomass fallbacks (every caller passes a
  table id or an organism) and ADM1's `_at()` has no caller: both logged in
  OPEN_WORK, not changed here.
- Tests: new `tests/standalone/test_model_species.py` (6): ADM1 declared
  species and every organism resolve, a misspelt id raises; BSM2 declares the
  lumps, its CV holds every table species, the monitor counts the lumps.
- For checkpoint 14: BSM2's headspace starts as `N2` (checkpoint 8a), which
  BSM2 does not declare, so `UnresolvedSpeciesWarning` will name it.

Sanity: suite 2180 passed, 2 xfailed before; 2186 passed, 2 xfailed, 0 failed
after (+6). Measurement against 11: 0 values differ in all five cases. BSM2
ConservationWarnings 9 -> 8: the 17.71 mol charge residual is gone, which is
the seeded `S_an` (0.00521 mol/L x 3400 L) the monitor could not count; the C,
H, N and O drift warnings are unchanged.

13:

Agreed before starting (2026-10-01/02), settling "how a corrector resolves":
a corrector is a dose, `{species_id: mol per mol of reagent}`, that states
what the real reagent adds (NaOH is `{"Na+": 1, "OH-": 1}`; a plain id is
one mole of it). Every id must be among the model's species; a dose that is
not charge-neutral warns. No `reagent=` keyword: it would only cross-check
two things the user wrote, and hydrates or dissolved-salt solutions would
not match a formula. Moved to design notes: solution dosing, solid
dissolution and engine-fit warnings (`upcoming/DOSING_AGENTS.md`), and strict
mass balance closure (`upcoming/MASS_BALANCE_CLOSURE.md`), which all three
engines lack for dosed H+ / OH- and water. OPEN_WORK gains the fixed liquid
volume (high priority; `fed_batch_fermenter.py` stays at 1600 L).

- New `PyOMES/chemistry/dose.py`: `as_dose`, `check_dose` (unknown id raises
  with the species available; net charge warns); exported from
  `PyOMES.chemistry`.
- `control_volume.py`: `_STRONG_CORRECTOR_ION` deleted;
  `equilibrate_to_pH(dose, ph_target, ...)` checks the dose against
  `cv.species`, adds each species in proportion, and requires at least one
  species in an equilibrium or charged (else raises); its message says
  "added N mmol of {...}". `apply_external_flux` applies ids as given.
- `cv_loops.py` `PHController`: `chemical_id` / `base_chemical_id` ->
  `acid_dose` / `base_dose` (id or composition, no default, at least one
  required); fluxes and `dosed_mol` are per species of the dose;
  `check_species(cv)`. `simulation.py`: `Simulation.__init__` calls
  `check_species` on each controller that has it, for its `target_cv_key` CV
  or the only CV.
- Engines (read, then tested): Bisection and NR compute H+ / OH- from charge
  balance and overwrite them, and recognise cations through fixed strong-ion
  tables; PHREEQC sets pH by charge from its `component_map`. So dosed OH-
  acts only through its counter-ion, and Na+ + OH- reaches a pH with exactly
  the Na+ that Na+ alone needs (new test, Bisection) - the old behaviour.
- The species check found undeclared ions in five tutorials/notebooks, now
  declared and passed as `species=`: `raw_construction.py` (Na+), D2C
  Example1 and Example2 (K+, Cl-, Na+; both CVs each), the two iron-oxidation
  validation notebooks (the medium's ions in no reaction), and
  `test_iron_oxidation.py` / `test_cv_advance.py` already passed theirs. D2C
  Example3 dosed H3PO4 with no phosphate chemistry, so its acid dose could
  never act (the silent-inert-dose bug); agreed: it now declares the
  phosphate ladder (log K -2.15 / -7.20 / -12.35, as the bioprocess
  database) and Na+.
- Callers: the template tutorials, `batch_fermenter.ipynb`, the D2C notebooks
  and script, `README.md`, `docs/tutorials/templates/README.md` (calls that
  relied on the old `"H3PO4"` default now say `acid_dose="H3PO4"`),
  `conftest.py`, `test_simulation.py` (17 one-line constructions, the default
  test becomes "a dose is required" + "a plain id is one mole", base dose
  asserts per species), `test_controller_state_protocol.py`.
- Tests: `test_cv_advance.py` (unknown dose raises, charged dose warns, ions =
  cation alone), `test_simulation.py` `TestPHControllerDoseCheck` (missing
  species raises, charged dose warns, the target CV is the one checked).
- Docs: `bioprocess_basic.py` and `raw_construction.py` docstrings; the
  `PHCONTROLLER_CORRECTOR_VALIDATION.md` note gains an update saying what is
  now done and what remains (can the dose move pH), and its README entry
  follows.
- Measurement script: its stirred-tank builder now passes
  `acid_dose="H3PO4", base_dose={"Na+": 1, "OH-": 1}`.

Sanity: suite 2186 passed, 2 xfailed before; 2192 passed, 2 xfailed, 0
failed after (+6). Measurement against 12: ADM1, BSM2, st_batch, st_fedbatch
0 values differ; d2c_raw 58 of 3618 differ, all in the recorded liquid OH-
on the steps where the base dose fires: the recorder samples n_mol after the
controller doses and before the next speciation overwrites OH-, so the
dosed OH- (4e-4 mol) shows there; pH, Na+ and every other series are
bit-identical. Tutorial scripts (four template scripts, `raw_construction.py`):
output identical to HEAD. Notebooks against HEAD (pH-correction message
normalised): D2C Example1 159/159 and Example2 130/130 lines identical,
Example3 identical except "3 equilibria" -> "6 equilibria" (the acid dose
never fires in that run), iron-oxidation 07 (60 lines) and 08 (62 lines)
identical (NR engine). `batch_fermenter.ipynb` fails at cell 10 on both
(logged).

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
