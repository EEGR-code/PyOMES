# Explicit Species Resolution — Design Note

> Status: design note, not yet started. No branch, no checklist, no code
> yet. First written 2026-09-22 while investigating whether
> `chemistry/common_species.py` should move to `PyOMES/databases/`.
> Rewritten 2026-09-30 after working through the OPEN_WORK item on
> molar-mass unification between `compounds.py`'s `Chemical` and
> `Species`: that item turned out to be one case of a wider problem, and
> the design decisions below were agreed with the repo owner in the same
> session.

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
| `reactions/stoichiometry.py` `_get_common_species` / `_parse_stoichiometry` | String stoichiometry seeds its lookup with every `Species` in `common_species`, then overlays the caller's `species=` | Looks up only the caller's `species=` |
| `reactions/equilibrium/interphase.py` `_resolve_species` | `HenryEquilibrium` / `RaoultEquilibrium` species given as strings are resolved from `common_species`; `RaoultEquilibrium` defaults both fields to `"H2O"` | Fields take `Species` objects; no `"H2O"` default |
| `core/control_volume.py` `_collect_species_registry` / `_common_species_catalog` | Conservation-monitor map built from reaction species, plus any `n_mol` id that matches a `common_species` name; everything else skipped. Built once, in `__init__` | Built from the species the CV was given; unresolved ids warn; ids that appear later are checked too |
| `monitoring/conservation.py` `_element_totals` / `_charge_residual` | Skip any id missing from the map, silently | Warn once per unresolved id |
| `reactions/kinetic/builder.py` `aerobic_growth` | Imports `CO2` / `H2O` from `common_species`, builds its own `O2`, and builds substrate / biomass / N-source `Species` from separate id, atoms and MW arguments | Receives all six `Species` from the caller; `species_overrides` becomes redundant |
| `compounds.py` `ChemicalRegistry`, read by `templates/stirred_tank/configs.py` (`OrganismConfig.resolve`, `SubstrateConfig.resolve`) and `factory.py` (`_build_reaction_system`) | A second name → (MW, atoms) table, 34 entries, only reachable through the template; always `ChemicalRegistry.default()`, so `register()` has no effect there | Deleted; the template resolves names against the species passed |
| `templates/stirred_tank/configs.py` `ChemistryConfig.acid_pKas` | Name-keyed pKa table; nothing reads it | Deleted; acid equilibria come from the reactions passed |
| `core/control_volume.py` `_STRONG_CORRECTOR_ION` | `equilibrate_to_pH` maps `"NaOH"` → `"Na+"`, `"KOH"` → `"K+"` and writes that id into `n_mol` | Corrector resolved from the model's species |

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
3. **pH dosing**, through `_STRONG_CORRECTOR_ION` (above).
4. **Template-created entries**: the stirred-tank factory adds a zero entry
   for every gas with a partition model, so CH4 or H2 can be in `n_mol`
   without being in any reaction.
5. **Lumped model variables**, such as BSM2's `S_cat` / `S_an`.
6. **Typos** (`"Na +"`, `"acetic_acid"`).

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
call site without `species=` has to pass its species (audit below).

### Phase 1: stirred-tank template and `Chemical`

- `OrganismConfig` / `SubstrateConfig` resolve their ids against the
  species passed to the template; a `Species` object is accepted directly
  in place of the loose `atoms=` / `MW=` arguments.
- `aerobic_growth` takes `Species` objects for all six participants.
- `create_volume` loses its `AD_BASIC` default; the builder and factory
  take `species=` / `reactions=` / `chemistry_db=`.
- Delete `compounds.py`, the `Chemical` / `ChemicalRegistry` exports in
  `PyOMES/__init__.py`, `tests/standalone/test_compounds.py`, the
  `conftest.py` fixture, and `ChemistryConfig.acid_pKas`.
- Yeast, Yeast_CHO and the organic acids become `Species` in the
  bioprocess database. Yeast keeps an explicit MW (it is about 8 % above
  its formula weight, 26.868 vs 24.834; Yeast_CHO 24.626 vs 22.593) with
  a comment giving the reason. The acids compute MW from atoms: identical
  for acetic, propionic and butyric acid; citric acid moves by
  0.001 g/mol unless it keeps an explicit MW.
- Stirred-tank results must be unchanged, checked with a before/after
  fingerprint. If the default database changes from `AD_BASIC` (the
  anaerobic-digestion database) to `BIOPROCESS_BASIC` in tutorials, check
  what that changes before accepting it.

This phase resolves OPEN_WORK's molar-mass item: there is then one
definition per compound.

### Phase 2: `ControlVolume` species set and warnings

- `ControlVolume` takes `species=` / `reactions=` / `chemistry_db=` and
  merges them into one species set (`chemistry_db` is stored today but
  nothing reads it).
- `_collect_species_registry` builds from that set; `_common_species_catalog`
  is deleted.
- Unresolved-id warnings at entry points and in the monitor, as above.
- `equilibrate_to_pH` resolves its corrector from the model's species
  instead of `_STRONG_CORRECTOR_ION`.

### Phase 3: `HenryEquilibrium` / `RaoultEquilibrium`

Species fields take `Species` objects only; `_resolve_species` and the
`"H2O"` defaults go. Callers import `H2O` (or whichever species) from a
database module.

### Phase 4: retire `common_species`

- Move its 26 definitions into the database modules (`aqueous.py` owns
  water, carbonate, ammonia and so on; the others build on it as they do
  today).
- Each shipped database lists every species it offers, including
  spectator ions (Na+, Cl-, K+), not only those in its reactions. Today
  `aqueous.py` lists only the 8 species in its reactions.
- Update every importer (66 files on 2026-09-30, mostly tests and
  notebooks) and delete the module and its `chemistry/__init__.py` export.

Phase 0 and Phase 1 are independent. Phase 2 needs the entry-point
arguments from Phase 1. Phase 4 goes last, once nothing scans the module.

## Verification

- Full suite before and after each checkpoint.
- Before/after result fingerprints for the stirred-tank template and BSM2.
- Notebooks edited as source only and re-run from the scratchpad where
  affected; long runs that cannot finish under the CPU cap are reported as
  not run.

## Audits needed before the checklist

1. `EquilibriumReaction(stoichiometry="...")` / `KineticReaction(stoichiometry="...")`
   call sites without `species=` (Phase 0).
2. `HenryEquilibrium` / `RaoultEquilibrium` construction with bare
   strings or defaults (Phase 3). On 2026-09-30, `RaoultEquilibrium(` or a
   string-valued species argument appears 35 times in 7 files, mostly
   tests; the list still needs sorting into bare-string and object calls.
3. Tests and notebooks relying on implicit lookup without importing
   anything. These only show up when each fallback is removed, so each
   removal gets its own checkpoint with the suite run.
4. How the ADM1 / BSM2 models in `models/` declare their species, and
   whether `S_cat` / `S_an` become `Species` here or in the strong-ion
   note.
5. Ids the speciation engines write into `n_mol` (for example the
   deprecated `{name}_HA` keys in OPEN_WORK), which must be declared or
   will warn.

## Open questions

1. **Subset helper.** Route 2 works with plain `species=` / `reactions=`
   arguments. Is a `ChemistryDatabase` method that builds a smaller
   database from chosen ids worth adding, or is picking items from
   `db.species` / `db.reactions` enough?
2. **Warning category.** A dedicated `UnresolvedSpeciesWarning` (easy to
   filter or promote to an error in tests) or plain `UserWarning`.
3. **Citric acid MW.** Keep 192.124 explicitly or compute 192.123.
4. **Yeast MW source.** The comment needs the real reason for the MW
   above formula weight (ash fraction is the likely one); to be confirmed.
5. **Where the merged species set lives** on a `ControlVolume` when both
   `species=` and `chemistry_db=` are given, and how `Simulation`'s
   thermo-mismatch check (which reads `chemistry_db.thermo`) fits with a
   CV built from loose arguments.

## Relationship to other notes

- [STRONG_ION_INFERENCE_GENERALIZATION.md](STRONG_ION_INFERENCE_GENERALIZATION.md):
  another closed name allowlist, in the speciation engines. Same rule;
  it can reuse this note's species set once Phase 2 lands.
- [PHCONTROLLER_CORRECTOR_VALIDATION.md](PHCONTROLLER_CORRECTOR_VALIDATION.md):
  the controller's corrector ids should resolve the same way as
  `equilibrate_to_pH`'s after Phase 2.
- OPEN_WORK's molar-mass item is resolved by Phase 1.

## How to start one

Per this folder's convention ([README.md](README.md#how-to-start-one)):
settle the open questions, run the audits, write
`EXPLICIT_SPECIES_RESOLUTION_CHECKLIST.md`, and cut a branch off `main`
(suggested name: `explicit-species-resolution`). Phases 0 and 1 can be
picked up without settling the Phase 2–4 questions first.
