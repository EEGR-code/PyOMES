# Explicit Species Resolution — Design Note

> Status: branch `explicit-species-resolution` cut; checklist not yet
> written; no code yet. First written 2026-09-22 while investigating
> whether `chemistry/common_species.py` should move to `PyOMES/databases/`.
> Rewritten 2026-09-30 after working through the OPEN_WORK item on
> molar-mass unification between `compounds.py`'s `Chemical` and
> `Species`: that item turned out to be one case of a wider problem, and
> the design decisions below were agreed with the repo owner in the same
> session. Later on 2026-09-30 the note was re-checked against the code,
> the five audits were run (see "Audit results") and the open questions
> were settled (see "Settled questions").

## The rule

**A model only knows the species and reactions it was given.**

Every species a model uses reaches it through one of three front-end
routes, all explicit in the model's own code:

1. **Build it yourself:** individual `Species` instances and reactions.
2. **Pick from a database:** import a shipped database module and pass the
   species and reactions you want from it.
3. **Pass a whole database** (`ChemistryDatabase`).

Any species name used anywhere in the model (`.organism("Yeast")`, a
string stoichiometry, a feed composition, a pH corrector) is looked up
only in that set. The backend never falls back to a module the model was
not given, and there is no default database.

## Decisions (2026-09-30)

1. **No ambient catalog.** Nothing in the backend scans
   `common_species` (or any other module) for names. `common_species.py`
   itself is retired: its definitions move into the database modules,
   which become the shipped definition lists. This answers the original
   relocation question — the module does not move, it goes away.
2. **No default database.** `StirredTankFactory.create_volume`'s
   `chemistry_db or AD_BASIC` fallback goes, as do the default ids
   `"Yeast"` / `"AceticAcid"` on `OrganismConfig` / `SubstrateConfig`.
3. **Flexible entry points.** Model entry points accept `species=` and
   `reactions=` as well as `chemistry_db=`, and they combine (a database
   plus a few extra species of the user's own). When the same id arrives
   with different data, `check_species_consistency` already raises
   `SpeciesConflictError`, so combining cannot silently pick one version.
4. **Names stay where there is a set to resolve them against.** The
   template call `.organism("Yeast")` keeps working when `"Yeast"` is in
   the species passed. A miss raises an error naming the id and listing
   what is available, e.g. `"Yeast" is not among the species passed to
   this model (available: AceticAcid, CO2, ...). Pass a Species(id="Yeast",
   ...) or a database that defines it.` Where no set is in reach (a
   standalone `HenryEquilibrium` / `RaoultEquilibrium`), the argument is a
   `Species` object, imported from a database module in one line. The API
   must stay easy for new users; error messages carry the fix.
5. **No species is ever silently skipped.** An id in a phase's `n_mol`
   that cannot be resolved to a `Species` raises a warning naming the id.
   Warn, not raise: the simulation itself does not need the `Species`
   (see below), only the conservation checks do.

## Where names are resolved today

| Site | What it does now | Under the rule |
|---|---|---|
| `reactions/stoichiometry.py` `_get_common_species` / `_parse_stoichiometry` | String stoichiometry seeds its lookup with every `Species` in `common_species`, then overlays the caller's `species=`. Called by `EquilibriumReaction`, `KineticReaction` and `KspEquilibrium` | Looks up only the caller's `species=` |
| `reactions/equilibrium/interphase.py` `_resolve_species` | `HenryEquilibrium` / `RaoultEquilibrium` species given as strings are resolved from `common_species`; `RaoultEquilibrium` defaults both fields to `"H2O"`. `databases/anaerobic_digestion.py` itself passes `"CO2"` strings, resolved at import | Fields take `Species` objects; no `"H2O"` default |
| `core/control_volume.py` `_collect_species_registry` / `_common_species_catalog` | Conservation-monitor map built from reaction species, plus any `n_mol` id that matches a `common_species` name; everything else skipped. Built once, in `__init__` | Built from the species the CV was given; unresolved ids warn; ids that appear later are checked too |
| `monitoring/conservation.py` `_element_totals` / `_charge_residual` | Skip any id missing from the map, silently | Warn once per unresolved id |
| `reactions/kinetic/builder.py` `aerobic_growth` | Imports `CO2` / `H2O` from `common_species`, builds its own `O2`, and builds substrate / biomass / N-source `Species` from separate id, atoms and MW arguments | Receives all six `Species` from the caller; `species_overrides` becomes redundant |
| `compounds.py` `ChemicalRegistry`, read by `templates/stirred_tank/configs.py` (`OrganismConfig.resolve`, `SubstrateConfig.resolve`) and `factory.py` (`_build_reaction_system`) | A second name → (MW, atoms) table, 34 entries, only reachable through the template; always `ChemicalRegistry.default()`, so `register()` has no effect there | Deleted; the template resolves names against the species passed |
| `templates/stirred_tank/builder.py` `.organism()` / `.substrate()` | Their own defaults `"Yeast"` / `"AceticAcid"` / `n_source_id="NH3"`, separate from the config dataclasses' defaults. The builder has no `chemistry_db` argument, so every built tank gets `AD_BASIC` | No default ids; the builder takes `species=` / `reactions=` / `chemistry_db=` |
| `templates/stirred_tank/factory.py` `create_volume` | Uses `chemistry_db` only to look up partition models for transfer species without a Henry constant; never passes it to the `ControlVolume` | Passes the model's species set through to the CV |
| `templates/stirred_tank/configs.py` `ChemistryConfig.acid_pKas` | Name-keyed pKa table; no runtime code reads it (two assertions in `tests/standalone/test_configs.py` do) | Deleted; acid equilibria come from the reactions passed |
| `core/control_volume.py` `_STRONG_CORRECTOR_ION` | `equilibrate_to_pH` maps `"NaOH"` → `"Na+"`, `"KOH"` → `"K+"` and writes that id into `n_mol`; `apply_external_flux` applies the same map to every external flux, which is how `PHController` doses NaOH | Corrector resolved from the model's species, on both paths |
| `models/vlmodels/adm1/base.py` `_get_species` | Any id not in ADM1's tables becomes a CHO biomass `Species` | Unknown ids raise |
| `chemical_equilibrium/engines/bisection/engine.py` `_CANONICAL_WRITEBACK_SPECIES` | A fixed tuple of 24 ids written into `n_mol` on every solve, zeros included, whatever the model declared; one of them (`Mo7O24------`) has no `Species` anywhere | Writes back only the species the engine was given |

Related name-keyed tables with their own notes, not duplicated here:
`_STRONG_ION_SPECIES_TO_KEY` and the `_STRONG_CHARGES` copies
([STRONG_ION_INFERENCE_GENERALIZATION.md](STRONG_ION_INFERENCE_GENERALIZATION.md)),
and `PHController`'s corrector ids
([PHCONTROLLER_CORRECTOR_VALIDATION.md](PHCONTROLLER_CORRECTOR_VALIDATION.md)).
Both should follow the rule above when they are picked up.

`thermo/gas/peng_robinson.py`'s critical-property table is keyed by name
too, but it holds property data for a species, not its identity, so it is
out of scope.

## How an id reaches `n_mol` without a `Species`

`n_mol` is a plain `{id: moles}` dict; nothing ties its keys to `Species`
objects. Any code that writes a new key creates an amount with no identity
behind it:

1. **Initial conditions typed by hand**, e.g.
   `LiquidPhase(n_mol={"Na+": 0.01, "Cl-": 0.01})`.
2. **Feeds**, whose compositions are keyed by id (e.g. `feed_conc_mol_L`).
3. **pH dosing**, through `_STRONG_CORRECTOR_ION` (above), in both
   `equilibrate_to_pH` and `apply_external_flux`.
4. **Template-created entries**: the stirred-tank factory writes O2, CO2
   and N2 into every tank's gas phase (BSM2's included, at zero), and a
   zero entry for every gas with a partition model, so CH4 or H2 can be in
   `n_mol` without being in any reaction. No shipped module defines a
   `Species` for O2, N2, CH4 or H2.
5. **Lumped model variables**, such as BSM2's `S_cat` / `S_an`, which have
   no `Species` object at all.
6. **Typos** (`"Na +"`, `"acetic_acid"`).
7. **Engine write-back**: the Bisection engine's fixed 24-id tuple
   (above); the NR engine always writes `H2O`; the PHREEQC engine writes
   one entry per species in its PHREEQC database.
8. **Model helpers run after the CV is built**:
   `seed_adm1_strong_ions` (Na+, Cl-), `seed_bsm2_strong_ions` (`S_cat`,
   `S_an`) and ADM1's water seeding in `build_adm1_cv`. Nothing at
   construction sees these; only the monitor's catch-all does.

The simulation carries these amounts correctly either way. What needs the
`Species` is the conservation monitor: element and charge totals need each
id's atoms and charge. So an unresolved id makes the drift checks blind to
that species, which is why it warns rather than raises.

Warnings are raised in two places:

- **Where the id enters:** phase construction, feeds, dosing and template
  set-up check new ids against the model's species and warn, naming the id
  and the call that introduced it.
- **Catch-all in the monitor:** each step, `n_mol` keys are compared with
  the map and each unresolved id warns once, whenever it first appears.
  This also closes an existing gap: today the map is built once in
  `ControlVolume.__init__`, so an id that first appears later (Na+ from a
  mid-run dose) is never looked up at all.

## Phases

Each phase ships on its own and leaves the suite green.

### Phase 0: string stoichiometry

Remove the `common_species` seed in `_parse_stoichiometry`; delete
`_get_common_species` and the `_cs_mod` import. Every string-stoichiometry
call site without `species=` has to pass its species (audit 1 below).
`KspEquilibrium` goes through the same parser and is covered by the same
change. The unknown-id error message stops listing common species and
lists the caller's instead. The `species` parameter docs on the three
reaction classes and `_parse_stoichiometry` stop mentioning
`common_species`.

### Phase 1: stirred-tank template and `Chemical`

- `OrganismConfig` / `SubstrateConfig` resolve their ids against the
  species passed to the template; a `Species` object is accepted directly
  in place of the loose `atoms=` / `MW=` arguments. The default ids go
  from both the config dataclasses and the builder's `.organism()` /
  `.substrate()` methods, including `n_source_id="NH3"`.
- `aerobic_growth` takes `Species` objects for all six participants.
- `create_volume` loses its `AD_BASIC` default; the builder (which has no
  database argument today) and the factory take `species=` / `reactions=`
  / `chemistry_db=`. Partition models for transfer species without a
  Henry constant come from the database passed; with none passed and no
  Henry constant, the existing "No partition model" error is raised.
- ADM1 and BSM2 are built through `StirredTankBuilder`, so this phase
  reaches them. ADM1 depends on the `AD_BASIC` default (its CH4, H2 and H2S
  transfer species have no Henry constant) and has to pass `AD_BASIC`
  explicitly. BSM2 passes a Henry constant for every transfer species and
  never reads the database.
- Delete `compounds.py`, the `Chemical` / `ChemicalRegistry` exports in
  `PyOMES/__init__.py`, `tests/standalone/test_compounds.py`, the
  `conftest.py` fixture, and `ChemistryConfig.acid_pKas` (with its two
  assertions in `test_configs.py`).
- Yeast, Yeast_CHO and the organic acids become `Species` in the
  bioprocess database, all with MW computed from atoms. For Yeast this is
  a deliberate change: the registry's 26.868 (Yeast) and 24.626
  (Yeast_CHO) are treated as a copy error, not a measured value (see
  settled question 4), and become 24.834 and 22.593. Acetic, propionic and
  butyric acid are unchanged; citric acid moves from 192.124 to 192.123.
- The tutorials' hard-coded `MW_yeast = 26.868`
  (`docs/tutorials/templates/cstr_fermenter.py`,
  `fed_batch_fermenter.py`) and the comments that cite the registry's
  Yeast_CHO (`D2C_workshop/raw_construction.py`,
  `reactions/reaction_system.ipynb`) follow.
- Before/after fingerprint of the stirred-tank template, ADM1 and BSM2.
  The yeast MW change moves stirred-tank results by about 8–9 % on
  gram-based quantities; it goes in its own checkpoint so that movement is
  measured and recorded separately from the refactor, which must leave
  results unchanged. BSM2 sentinels must not move. If tutorials change
  database from `AD_BASIC` (the anaerobic-digestion database) to
  `BIOPROCESS_BASIC`, check what that changes before accepting it.

This phase completes the molar-mass unification that was moved here from
OPEN_WORK: there is then one definition per compound.

### Phase 2: `ControlVolume` species set and warnings

- `ControlVolume` takes `species=` / `reactions=` / `chemistry_db=` and
  merges them, together with the species in its reaction stoichiometries,
  into one read-only mapping, `cv.species`, built once in `__init__`.
  `chemistry_db` is still stored exactly as given; no database is built
  from loose arguments. `Simulation._warn_thermo_mismatch`, which reads
  `cv.chemistry_db.thermo`, keeps skipping CVs whose `chemistry_db` is
  `None`.
- The merge needs a species-level conflict check.
  `check_species_consistency` walks only reaction stoichiometries, and
  `ChemistryDatabase.extend` overrides on a key collision without a word,
  so neither catches the same id arriving with different data from
  `species=` and `chemistry_db=`. The new check raises
  `SpeciesConflictError`, as decision 3 requires.
- `_collect_species_registry` builds the monitor's map from `cv.species`;
  `_common_species_catalog` is deleted.
- Unresolved ids warn with `UnresolvedSpeciesWarning`, a `UserWarning`
  subclass in `monitoring/` next to `ConservationWarning`, at entry
  points and in the monitor, as above.
- `equilibrate_to_pH` and `apply_external_flux` resolve their corrector
  from the model's species instead of `_STRONG_CORRECTOR_ION`.
- The Bisection engine writes back only the species it was given: those
  in its declared equilibria and the strong ions present in its input.
  `_CANONICAL_WRITEBACK_SPECIES` is deleted. Checked with the BSM2
  fingerprint.
- BSM2 declares `S_cat` / `S_an` as `Species` (`atoms={}`, charge ±1) in
  its species table, so BSM2 runs do not warn. Only the conservation
  monitor reads them, so results do not change. The strong-ion note then
  only has to derive their charge from the `Species`.
- ADM1's `_get_species` stops turning unknown ids into biomass and raises
  instead.
- Known warnings this phase leaves in place: the NR engine's unconditional
  `H2O` write and the PHREEQC engine's full species list warn for any id
  the model does not declare. The checklist records which tests and
  notebooks show them.

### Phase 3: `HenryEquilibrium` / `RaoultEquilibrium`

Species fields take `Species` objects only; `_resolve_species` and the
`"H2O"` defaults go. Callers import `H2O` (or whichever species) from a
database module. `databases/anaerobic_digestion.py` passes `"CO2"` strings
and is resolved at import, so it changes in the same checkpoint that
removes `_resolve_species`; otherwise importing the template fails and
takes the whole suite with it. `models/vlmodels/adm1/base.py`'s water link
(`RaoultEquilibrium()` on the default) changes in the same checkpoint.

### Phase 4: retire `common_species`

- Move its 26 definitions into the database modules (`aqueous.py` owns
  water, carbonate, ammonia and so on; the others build on it as they do
  today).
- Each shipped database lists every species it offers, including
  spectator ions (Na+, Cl-, K+), not only those in its reactions, and a
  `Species` for every id it has a partition model for (O2 and N2 in the
  bioprocess database; CH4 and H2 in the anaerobic-digestion database).
  Today `aqueous.py` lists only the 8 species in its reactions, and no
  shipped module defines O2, N2, CH4 or H2. `aerobic_growth`'s own `O2`
  goes.
- Update every importer (62 files import `common_species` on 2026-09-30,
  mostly tests and notebooks; 8 more only mention it) and delete the
  module and its `chemistry/__init__.py` export.

Phase 0 and Phase 1 are independent. Phase 2 needs the entry-point
arguments from Phase 1. Phase 3 is independent of Phases 0–2. Phase 4
goes last, once nothing scans the module.

## Verification

- Full suite before and after each checkpoint.
- Before/after result fingerprints for the stirred-tank template, ADM1
  and BSM2. Every checkpoint leaves them unchanged except the Phase 1
  yeast molar-mass checkpoint, whose movement is recorded.
- Notebooks edited as source only and re-run from the scratchpad where
  affected; long runs that cannot finish under the CPU cap are reported as
  not run.

## Audit results (2026-09-30)

Searched: every tracked `.py` file and notebook code cell by syntax tree,
and `.md`, config and notebook text by grep, excluding
`docs/dev/implementation/shipped/` and `obsolete/`; the three
`vars(common_species)` scans and any `getattr` / `import_module` lookups.
The full suite was also run with a recording pytest plugin (kept outside
the repo, no behaviour change) that logged every id reaching a fallback:
2144 passed, 0 failed. Line numbers are as of that date.

1. **String stoichiometry without `species=` (Phase 0).**
   `tests/standalone/test_stoichiometry.py` has 16 sites that take at
   least one id from `common_species` (lines 25, 34, 47, 58, 71, 84, 95,
   103, 108, 113, 120, 128, 135, 144, 221, 244; 47 and 58 pass `_LOCAL`
   but still take `H+` from `common_species`), and
   `test_unknown_species_message_lists_common` (line 186) asserts on the
   old message. `docs/tutorials/D2C_workshop/Example3_CSTR.ipynb` cell 14
   has three `EquilibriumReaction`s with no `species=`. Nothing in
   `PyOMES/`, `models/` or any other test; no `KspEquilibrium` string
   callers.
2. **`HenryEquilibrium` / `RaoultEquilibrium` (Phase 3).** 104
   constructions. 29 `RaoultEquilibrium` on the `"H2O"` default: one in
   `models/vlmodels/adm1/base.py`, the rest in `test_partition_model`,
   `test_raoult_h2o_fold_cp4`, `test_equilibrium_constraint`,
   `test_equilibrium_classification` and `test_simulation`. 11
   `HenryEquilibrium` with bare strings: `databases/anaerobic_digestion.py`,
   the two `ChemicalEquilibriumProtocol` engine-basics notebooks (cell 13),
   and `test_nr_tableau_gas_liquid` (3), `test_equilibrium_classification`,
   `test_equilibrium_constraint`, `test_equilibrium_constraint_dual_role`,
   `test_nr_gas_liquid_cp2` and `test_precipitation_gas_liquid_cp5`. The
   other 64 pass `Species` objects (4) or no species (60) and are
   unaffected. At run time, 30 call sites across 9 test files resolve a
   string.
3. **Implicit lookup, by fallback.**
   - Stoichiometry seed and `_resolve_species`: items 1 and 2.
   - CV catalog: only three tests reach it
     (`test_accuracy_monitor.py:562`, `test_cv_advance.py:1019`,
     `tests/validation/speciation/test_iron_oxidation.py:160`). None
     asserts on conservation, so removing it fails nothing; they warn
     after Phase 2.
   - `ChemicalRegistry`: the four template tutorials (`batch_fermenter.py`,
     `cstr_fermenter.py`, `fed_batch_fermenter.py`,
     `batch_fermenter.ipynb`), most `.organism()` / `.substrate()` calls in
     `test_builder.py`, the `resolve()` tests in `test_configs.py`,
     `test_compounds.py` and the `conftest.py` fixture.
   - `AD_BASIC` default: every builder user except BSM2 (ADM1, the
     template tutorials, `reactions/aerobic_fermentation_stoichiometry.ipynb`,
     most of `test_builder.py`, `test_simulation.py`'s builder tests), plus
     direct `create_volume` calls in `test_simulation.py` (1905, 3298, 3504),
     `test_builder.py:404` and `results/01_exporting_results.ipynb`. The
     default O2/CO2/N2 transfer needs partition models, so each needs a
     database passed.
   - Monitor warnings: no global `filterwarnings`, but seven test files
     turn warnings into errors locally; `test_cv_advance.py` and
     `test_simultaneous_adaptive_solver_jac.py` meet unresolved ids today.
4. **ADM1 / BSM2.** Both declare species in module tables that reuse the
   `common_species` CO2, NH3 and H2O objects (BSM2's `SPECIES`, ADM1's
   `_get_species` with its biomass fallback). `S_cat` / `S_an` have no
   `Species` and are written after the CV is built; the monitor skips
   `S_an` today. They become `Species` in Phase 2 (above). Other ids
   these models carry without a `Species`: O2 and N2 from the factory,
   Na+ and Cl- from `seed_adm1_strong_ions`, and the Bisection engine's
   zero-filled ions.
5. **Engine write-back.** At run time, engine write-back created 16 ids
   in `n_mol`: CO2, CO3--, H+, H2O, HCO3-, NH4+, OH-, Ca++, Cl-, Co++,
   K+, Mg++, Mn++, Mo7O24------, Na+ and Zn++, most from the Bisection
   engine's fixed tuple. The deprecated `{name}_HA` keys go to the
   result's `extra` dict and never reach `n_mol`.

## Settled questions (2026-09-30)

1. **Subset helper: none for now.** Picking items from `db.species` /
   `db.reactions` expresses route 2, and reactions carry their own
   `Species`. A subset method would also have to decide what happens to
   `partition_models`. Revisit if tutorials repeat the same boilerplate.
2. **Warning category: `UnresolvedSpeciesWarning`**, a `UserWarning`
   subclass in `monitoring/`, following `ConservationWarning` and
   `AccuracyWarning`. It can be filtered or promoted to an error on its
   own; seven test files already turn warnings into errors.
3. **Citric acid MW: computed, 192.123.** It matches every other acid and
   the `H3Cit` the iron-oxidation notebooks already define. Nothing run by
   tests or fingerprints uses the registry's `CitricAcid`.
4. **Yeast MW: computed from atoms.** No source is recorded: both values
   first appear, uncommented, in the initial public release. 24.626 is
   exactly the formula weight of the generic biomass formula
   CH1.8O0.5N0.2 (24.6263), and the initial release gave an `Ecoli`
   species with that formula the same MW; 26.868 is 24.626 plus
   0.16 × 14.007 (26.867). So both are taken as copy errors, not measured
   values, and Yeast becomes 24.834, Yeast_CHO 22.593. The stirred-tank
   results this moves are recorded in their own checkpoint (Phase 1).
5. **Merged species set: `cv.species`**, a read-only mapping built once
   in `ControlVolume.__init__` (Phase 2). `chemistry_db` is stored as
   given, and `Simulation`'s thermo-mismatch check keeps skipping CVs
   without one.
6. **Bisection write-back: fixed in Phase 2.** The engine's fixed 24-id
   tuple is backend knowledge the model never gave it, the same kind of
   hidden name list the rule removes. It writes back only the species it
   was given.

## Relationship to other notes

- [STRONG_ION_INFERENCE_GENERALIZATION.md](STRONG_ION_INFERENCE_GENERALIZATION.md):
  another closed name allowlist, in the speciation engines. Same rule;
  it can reuse this note's species set once Phase 2 lands. Phase 2 here
  removes the Bisection engine's fixed write-back tuple and gives
  `S_cat` / `S_an` `Species`; the strong-ion note keeps the
  `_STRONG_ION_SPECIES_TO_KEY` and `_STRONG_CHARGES` tables.
- [PHCONTROLLER_CORRECTOR_VALIDATION.md](PHCONTROLLER_CORRECTOR_VALIDATION.md):
  the controller's corrector ids should resolve the same way as
  `equilibrate_to_pH`'s after Phase 2, which also changes
  `apply_external_flux`, the path `PHController` doses through.
- The molar-mass unification item was moved here from OPEN_WORK on
  2026-09-30 and is completed by Phase 1.

## How to start one

Per this folder's convention ([README.md](README.md#how-to-start-one)).
The branch `explicit-species-resolution` is cut, the audits are run and
the questions settled; what remains is
`EXPLICIT_SPECIES_RESOLUTION_CHECKLIST.md`, with each fallback removal
in its own checkpoint and the suite run after each.
