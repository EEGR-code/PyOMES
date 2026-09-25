# Phases — Upcoming

Planning notes for refactor work that is **not yet started**. Each
note describes the scope, design choices to resolve, and triggers
that would justify picking it up. Notes here are stable —
implementation only begins when one is moved through to active work
(i.e. a checklist file is added) and ultimately to
[../shipped/](../shipped/) once shipped.
+-*
**Finding work in progress.** This folder does not track which phases are
currently being worked on. Each phase lives on its own feature branch (see
[How to start one](#how-to-start-one)), so list the active ones with
`git branch -a`.

**History below the design notes.** Most sections after "Design discussions"
("Recently shipped", the "Recently surfaced" sections, "Priority order")
record work that has already shipped and are kept for context. The sections
that still describe open work are "Open phases" and the pending stages in
"Solver interface refinement".

## Design discussions (pre-phase, not yet a checklist)

- **[ACTIVITY_MODEL_PARAMETER.md](ACTIVITY_MODEL_PARAMETER.md)** —
  2026-09-25. Replaces the `use_activity: bool` + `activity_model: str` pair
  with one `activity_model` parameter taking a name (`"ideal"`, `"davies"`,
  `"sit"`) or a model object, across `make_activity_model`, both engines,
  `ReactionSystem.configure_engine`, `StirredTankBuilder.chemistry()` and
  `ChemistryConfig`. Removes `use_activity` everywhere (including
  `ThermoFramework`) and the engines' `thermo=` argument; engines resolve the
  model once at construction, and pick the ideal fast path by model type.
  Pure refactor, verified by an engine-output fingerprint. Direction approved;
  checklist in [ACTIVITY_MODEL_PARAMETER_CHECKLIST.md](ACTIVITY_MODEL_PARAMETER_CHECKLIST.md);
  no branch or code yet.
- **[EXPLICIT_SPECIES_RESOLUTION.md](EXPLICIT_SPECIES_RESOLUTION.md)** —
  2026-09-22. Surfaced while investigating whether `chemistry/
  common_species.py` should move to `PyOMES/databases/`: three internal
  call sites (`reactions/stoichiometry.py`'s string-stoichiometry
  parser, `reactions/equilibrium/interphase.py`'s `HenryEquilibrium`/
  `RaoultEquilibrium` species fields, `core/control_volume.py`'s
  charge-conservation registry) resolve unrecognized species ids by
  scanning `common_species.py`'s entire module namespace via `vars()`,
  not from anything the model itself declared — so a model can silently
  pick up (or silently drop, depending on name collision) species the
  user never wrote. Proposes removing all three ambient fallbacks in
  favor of explicit resolution only. Phase 0 (stoichiometry.py) is
  cheap and decided; Phase 1 (control_volume.py) found `ControlVolume`
  already accepts `chemistry_db=` but doesn't consult
  `chemistry_db.species`, so it's mostly wiring; Phase 2 (interphase.py)
  needs an API decision (`RaoultEquilibrium.liquid_species` currently
  defaults to the bare string `"H2O"`, resolved ambiently). The
  original relocation question is downstream of this note, not
  parallel — see its own "Relationship to the relocation question"
  section. No branch, no checklist, no code yet.
- **[PHCONTROLLER_CORRECTOR_VALIDATION.md](PHCONTROLLER_CORRECTOR_VALIDATION.md)** —
  2026-09-17. Surfaced while fixing `tutorials-followups` checkpoint 3
  (`raw_construction.py`'s pH runaway): `PHController` should warn when its
  configured `chemical_id`/`base_chemical_id` can't actually shift pH (not a
  recognised strong-corrector alias and not declared in any equilibrium
  reaction) — the same check `ControlVolume.equilibrate_to_pH` already has,
  just not reused on `PHController`'s actual dosing path. Scoped strictly to
  `PHController`'s own two fields; no other controller/boundary is touched.
  Open questions: warn vs. raise, new warning category vs. plain
  `UserWarning`. No branch, no checklist, no code yet.
- **[REACTION_ENVIRONMENT_PHASE_EXPOSURE.md](REACTION_ENVIRONMENT_PHASE_EXPOSURE.md)** —
  2026-09-17. Surfaced while fixing `tutorials-followups` checkpoint 4
  (iron-oxidation Singer-Stumm rate law): a `KineticReaction`'s `rate_fn`
  only ever sees liquid-phase concentrations (`ControlVolume.
  _build_reaction_environment` never merges in gas-phase state), forcing a
  manual Henry's-law conversion for any literature rate law defined in terms
  of a gas partial pressure. Proposes exposing gas (and, by the same
  reasoning, solid) phase state to `ReactionEnvironment` as an additive,
  read-only extension. Open questions: field shape, default for CVs with no
  gas phase. No branch, no checklist, no code yet.
- **[SCIPY_REJECTION_CHECK_SOLVER_AWARENESS.md](SCIPY_REJECTION_CHECK_SOLVER_AWARENESS.md)** —
  2026-09-17. Surfaced in the same checkpoint 4 investigation:
  `AccuracyMonitor.check_scipy_rejections` judges `solve_ivp`'s `nfev`/
  accepted-steps ratio against one flat threshold calibrated for explicit
  methods, producing a confirmed structural false positive for implicit
  solvers (BDF/Radau), whose per-step Newton-iteration overhead inflates
  `nfev` by design. Proposes a solver-family-aware threshold. No branch, no
  checklist, no code yet.
- **[STRONG_ION_INFERENCE_GENERALIZATION.md](STRONG_ION_INFERENCE_GENERALIZATION.md)** —
  2026-09-18. Surfaced while explaining `NRChemicalEquilibriumEngine.solve()`'s
  `strong_ions=` kwarg for `docs/tutorials/ArXiv_preprint/
  01_predict_ph_simple_liquid.ipynb`: the phase-based `solve(phases=...)`
  path derives `strong_ions` from a hardcoded, closed
  `_STRONG_ION_SPECIES_TO_KEY` allowlist (~15 species) rather than
  structurally, from tableau-component non-membership + the species'
  own `.charge` — so an off-allowlist charged species silently drops out
  of the charge balance. Also finds `_STRONG_CHARGES` duplicated
  verbatim across three locations in `engines/nr/engine.py`/`engines/nr/solver.py`, plus
  a fully independent copy of the same allowlist pattern in the older
  `BisectionChemicalEquilibriumEngine` (`engines/bisection/engine.py`/`acid_base.py`) — left
  out of scope, since that engine is still `ReactionSystem`'s *default*
  solver today, not legacy. Scoped as Phase 0 (decided: de-duplicate the
  three copies, no behavior change) + Phase 1 (derive structurally,
  liquid-phase-only by explicit rule); resolves the `S_cat`/`S_an` generic
  charge-lump question (they should be concrete `Species` declarations,
  not a special case) along the way. No branch, no checklist, no code yet.
- **[NOTEBOOK_GENERATOR_REMOVAL.md](NOTEBOOK_GENERATOR_REMOVAL.md)** —
  2026-09-15, re-audited 2026-09-20. Retires the 2 remaining
  `_generate_notebooks.py` scripts (13 notebooks: `ArXiv_preprint` and
  `tests/validation/speciation`; 5 scripts at the original audit, the other
  three have since been deleted by other phases) in favor of every notebook
  being hand-edited and committed with its outputs embedded — the
  biosteam-style convention, chosen over keeping the generators (with a
  CI discipline bolted on) or extracting their duplicated boilerplate
  into a shared module. Confirmed real duplication exists (the
  water/phosphate/ammonium reaction network retyped in 3+ places) and is
  a known, accepted cost of this choice. Adds a CI check for
  never-executed code cells (drift between source and saved outputs; the
  exact rule is still an open question). No Sphinx/mkdocs/jupyter-book
  — explicitly out of scope. No branch, no checklist, no code yet.
- **[RESERVOIR_TYPE.md](RESERVOIR_TYPE.md)** — 2026-07-10, revised
  2026-07-10. Now concludes with `FlowBoundary`, a real unifying protocol
  for `PhaseInterface`/`CVLink`/`ExternalBoundary`, enabled by a
  `PhaseCarrier` structural type (`ControlVolume` and the sibling
  `Reservoir` type both satisfy it). `Reservoir` (passive audit
  accumulator + `Phase`-borrowed multi-phase accounting, deliberately not
  `ControlVolume`-shaped) is the piece that gives `ExternalBoundary` a
  symmetric "other side," which is what makes the unification possible.
  Extends `MASS_EXCHANGE_ARCHITECTURE.md` §7/§12 Q1. `PartitionModel`'s
  relationship to `FlowBoundary` is explicitly deferred — shared open
  thread with `PHENOMENA_PROTOCOL.md`. Open questions listed in the note;
  no branch, no checklist, no code yet.
- **[PHENOMENA_PROTOCOL.md](PHENOMENA_PROTOCOL.md)** — 2026-07-10, Phase 1
  of two. A three-tier `Phenomena → {EquilibriumPhenomena, KineticPhenomena}`
  taxonomy grouping reactions and transfer models by evaluate-at-a-point
  vs. algebraic-constraint shape (per `MASS_EXCHANGE_ARCHITECTURE.md`
  §14.2's litmus test) rather than by chemistry/transport domain.
  `EquilibriumPhenomena` is largely already real (`EquilibriumConstraint` +
  `LAYER1_GAP_CLOSURE`'s Newton-fold); `KineticPhenomena` is supported by
  `SimultaneousAdaptiveSolver`'s existing combined-RHS evaluation.
  `ReactionSystem` renamed `PhenomenaSystem` in this note (role unchanged —
  bucket, own the engine, execute kinetic rates, provide the
  `ControlVolume` attachment point); `Phenomena` connects to `ControlVolume`
  through it, not directly. No behavior change scoped here — pure taxonomy.
  Followed by
  [KINETIC_TRANSFER_GENERALIZATION.md](KINETIC_TRANSFER_GENERALIZATION.md).
- **[KINETIC_TRANSFER_GENERALIZATION.md](KINETIC_TRANSFER_GENERALIZATION.md)** —
  2026-07-10, Phase 2 of two (depends on `PHENOMENA_PROTOCOL.md`). Gives
  `KineticTransferModel` an optional custom rate function (currently
  constrained to linear `kLa` relaxation) so it can actually satisfy
  `KineticPhenomena` instead of staying inert config deferred to
  `KineticGasLiquidLink`. `PhaseInterface`, which the generalized class
  would satisfy, is now `FlowBoundary`'s same-CV case per
  `RESERVOIR_TYPE.md` §5.1 — framing note only, no design change. Higher-risk half — touches the shipped
  `transfer_models=` API. Five open questions listed in the note
  (signature, `partition_model` optionality, `KineticGasLiquidLink`'s
  role, `EquilibriumTransferModel` parity, regression scope); no branch,
  no checklist, no code yet.
- **[DEMO_RECORDERS.md](DEMO_RECORDERS.md)** — planned after `run-history`
  shipped (2026-06-08). A short demo script comparing the four recorder
  variants (`BatchRecorder`, `StreamingFileRecorder`, `SparseRecorder`,
  `SummaryRecorder`) on the same simulation. Small and non-blocking — the
  trigger is a real model exercising the recorders, or a user asking for an
  example. Its proposed location, `demos/model_api/recorder_comparison.py`,
  no longer exists (`demos/` was retired 2026-09-17), so it needs a new home,
  likely under `docs/tutorials/`. No branch, no checklist, no code yet.
## Recently shipped

- `thermo-subfolder-structure` (2026-09-24) — grouped the 8 flat files in
  `PyOMES/thermo/` by phase: `liquid/` (`liquid_phase_model.py` split into
  `protocols.py`, `ideal.py` and `davies.py`; `sit_liquid_model.py` renamed
  `sit.py`; `water_properties.py`, `factory.py`) and `gas/` (`gas_eos.py` split
  into `protocols.py`, `ideal.py` and `peng_robinson.py`), with `framework.py` and
  `equilibrium_constants.py` staying at the top level. Seven checkpoints. The two
  moves were pure refactors (definitions identical by AST, numeric fingerprints
  bit-identical); there are no shims, so the old deep import paths stop working,
  and package-root exports are unchanged. Two small changes rode along, each in
  its own checkpoint: `GasEOS` is now a runtime-checkable `Protocol`, which
  `PengRobinsonEOS` satisfies as its docstring always claimed; and the two private
  `_kg_per_L` copies became one `water_kg_per_L` helper (bit-identical). A third
  test in `tests/standalone/test_package_layering.py` keeps `liquid/` and `gas/`
  from importing each other, `framework`, `equilibrium_constants` or either
  package root. The docstring pass found statements that were false, not just
  stale: SIT's `compute_gammas` is called by the Bisection engine, not NR, and
  nothing reads `ThermoFramework.gas_eos` (it cited a `KineticGasLiquidLink`
  default). The design note's audit needed correcting (line endings, line numbers
  moved by the reactions phase, five copies of the mol/L → mol/kg conversion
  rather than three). Logged in `OPEN_WORK.md`: SIT's two inline conversions lack
  the bad-density fallback, and docstring examples in `models/` and `numerics/`
  still use pre-rename module paths. Full suite green post-merge: 2104 passed, 0
  failed. Tag `thermo-subfolder-structure-shipped`. See
  [`../shipped/THERMO_SUBFOLDER_STRUCTURE_CHECKLIST.md`](../shipped/THERMO_SUBFOLDER_STRUCTURE_CHECKLIST.md).

- `reactions-subfolder-structure` (2026-09-24) — grouped the 13 flat files in
  `PyOMES/reactions/` by the kind of reaction they serve: `kinetic/` (`reaction.py`
  with `KineticReaction`, `rate_laws.py`, `builder.py`) and `equilibrium/`
  (`equilibrium.py` split into `constraint.py`, holding the `EquilibriumConstraint`
  protocol, `vant_hoff_log_K` and `classify_equilibrium_constraint`, and
  `reaction.py`; `phase_equilibria.py` renamed `interphase.py`; `plots.py`), with
  `stoichiometry`, `environment`, `protocols`, `_shared`, `reaction_system` and
  `blackbox` staying at the top level. Six checkpoints, no behaviour change. There
  are no shims, so the old deep import paths stop working; package-root exports
  are unchanged. Tests, notebooks, both notebook generators and `models/` now
  import exported names, including `StoichiometryEntry`, from `PyOMES.reactions`,
  while code inside `PyOMES/` keeps deep imports. A second test in
  `tests/standalone/test_package_layering.py` keeps `kinetic/` and `equilibrium/`
  from importing each other. The design note's audit needed correcting: several
  counts were off by one or two, both halves of the split need
  `phases_from_entries`, and the Bisection engine depends on both halves, not only
  `constraint.py`. `README.md`, `PyOMES/README.md` (which said `ReactionBuilder`
  builds `EquilibriumReaction` objects) and the `docs/architecture.md` tree (which
  lacked `rate_laws.py` and `plots.py`) were corrected on the way. Two working
  notes: the aerobic-fermentation tutorial notebook's 60 h simulation cells were
  not re-run (the tool sandbox caps CPU; every other cell was), and with
  `core.autocrlf=true` a moved file is stored under its new path as LF unless it is
  staged with autocrlf off, which makes git show it as rewritten. Full suite green
  post-merge: 2100 passed, 0 failed. Tag `reactions-subfolder-structure-shipped`.
  See
  [`../shipped/REACTIONS_SUBFOLDER_STRUCTURE_CHECKLIST.md`](../shipped/REACTIONS_SUBFOLDER_STRUCTURE_CHECKLIST.md).

- `partition-constraint-relocation` (2026-09-24) — moved `HenryEquilibrium`/
  `RaoultEquilibrium`/`KspEquilibrium` (plus `_resolve_species`) from
  `chemistry/partition.py` to a new `reactions/phase_equilibria.py`, exported from
  `reactions/__init__.py`. That removed the `chemistry -> reactions` package edge:
  `chemistry/` now imports only `units`, enforced by a new
  `tests/standalone/test_package_layering.py` that also counts function-level and
  `TYPE_CHECKING` imports. Six checkpoints, no behaviour change beyond the import
  path. There is no re-export from `chemistry/`, so the old import path stops
  working, and checkpoints pickled before the phase will not load. The design
  note's audit needed correcting: the package graph does not become acyclic (two
  other package cycles, `control` <-> `core` and `chemical_equilibrium` <->
  `reactions`, are now logged in `OPEN_WORK.md`, unfixed), the "31 external files"
  were 26 needing edits, and live notes and docstrings the note did not list
  carried stale paths. Also added ten tests showing `RaoultEquilibrium`'s water
  values are already configurable (including another solvent), and logged that
  those values exist in three places with no shared source. Full suite green
  post-merge: 2098 passed, 0 failed. Tag `partition-constraint-relocation-shipped`.
  See
  [`../shipped/PARTITION_CONSTRAINT_RELOCATION_CHECKLIST.md`](../shipped/PARTITION_CONSTRAINT_RELOCATION_CHECKLIST.md).

- `equilibrium-set-relocation` (2026-09-23) — moved `EquilibriumSet`/
  `EquilibriumDef`/`WaterDef` (the acid-base charge-balance declaration
  format) from `chemistry/equilibria.py` to
  `chemical_equilibrium/engines/bisection/equilibria.py`, next to the one
  engine that consumes it, and deleted `EquilibriumSet.bsm2_default()` (a
  vestigial preset that duplicated constants `models/vlmodels/adm1/bsm2.py`
  already declares; its only two callers were tests, which now build a small
  synthetic set). Two checkpoints, no behaviour change. The design note's
  consumer audit undercounted: re-verifying it found comment and docstring
  references in `acid_base.py`, `bsm2.py`, `test_speciation.py` and
  `OPEN_WORK.md` that would have gone stale, all updated. Also found and
  logged in `OPEN_WORK.md`: the deprecated `{name}_HA` fallback in
  `_compute_species_eq` is unreachable from repo code (BSM2 does not use it,
  contrary to an earlier entry), and the package mixes three-dot relative
  imports (`engines/bisection/`, `engines/nr/`) with absolute ones elsewhere.
  Full suite green post-merge: 2086 passed, 0 failed. Tag
  `equilibrium-set-relocation-shipped`. See
  [`../shipped/EQUILIBRIUM_SET_RELOCATION_CHECKLIST.md`](../shipped/EQUILIBRIUM_SET_RELOCATION_CHECKLIST.md).

- `chemistry-reactions-kinetics-cleanup` (2026-09-23) — review of `chemistry/`,
  `kinetics/` and `reactions/` in the style of
  `chemical-equilibrium-engines-subfolder`. 20 checkpoints: dead code deleted
  (`kinetics/`, `thermo_params.py`, the recipe layer — `chem_recipe.py`,
  `recipe.py`, `types.py`, `registry.py` — the top-level `equilibria/`
  package, deprecated `HenryPartition`/`RaoultPartition` aliases, the three
  unused `EquilibriumSet` presets, `tests/legacy/`); import deferrals hoisted
  that no longer worked around anything; three Monod implementations
  unified into one; two defects fixed (`ChemistryDatabase.extend()` was
  dropping the solver/label/engine config; `PHController`'s constructor-time
  id validator checked the wrong, since-deleted table — dropped rather than
  fixed, no replacement built); `chemistry/database.py`/`chemistry/databases/`
  moved to `PyOMES/databases/` (D7); `chemistry.__all__` shrunk to names
  still used, plus a new `plots` extra. Three checkpoints added after the
  original D1-D8 audit closed, from a conversational review that also
  produced [`EXPLICIT_SPECIES_RESOLUTION.md`](EXPLICIT_SPECIES_RESOLUTION.md)
  (logged separately, not part of this phase): `_ATOMIC_WEIGHTS` moved from
  `chemistry/species.py` to `PyOMES/units.py` (a physical-constants table,
  not domain data); `chemistry/compounds.py` moved to `PyOMES/compounds.py`
  (no dependency on `Species` or anything else in `chemistry/`); and
  `ChemicalRegistry.IDs`/`FeedState`/`stream_adapter.py` deleted outright
  (confirmed zero consumers anywhere in the repo outside their own tests).
  Findings logged rather than fixed along the way — see `OPEN_WORK.md` for
  van 't Hoff/constants duplication, gas EOS consistency, multi-species
  rate-law parameters, the remaining `chemistry`<->`reactions` package-level
  cycle via `partition.py`, and molar-mass unification between `Chemical`
  and `Species`. Full suite green post-merge: 2086 passed, 0 failed. Tag
  `chemistry-reactions-kinetics-cleanup-shipped`. See
  [`../shipped/CHEMISTRY_REACTIONS_KINETICS_CLEANUP_CHECKLIST.md`](../shipped/CHEMISTRY_REACTIONS_KINETICS_CLEANUP_CHECKLIST.md).

- `chemical-equilibrium-engines-subfolder` (2026-09-20) — moved the three
  chemical-equilibrium engines (Bisection, NR, PHREEQC) out of a flat 16-file
  `chemical_equilibrium/` package into `engines/{bisection,nr,phreeqc}/`
  subpackages, with `protocols.py` and `numerical_gradient.py` staying shared
  at the top level (Part C). Along the way: retired the
  `activity_models.py`/`sit.py` compatibility shims in favour of `PyOMES.thermo`
  (Part A); deleted the orphaned `strong_ions.py` (Part B), `api.py`/`factory.py`
  (checkpoint 9b), and `activity_dispatch.py` (checkpoint 12c); dissolved
  `activity.py` into `engines/bisection/ionic_strength.py`, deleting the
  never-called `warn_if_high_ionic_strength` (checkpoint 12d); rewrote
  `chemical_equilibrium/`'s docstrings to describe current behaviour rather
  than development history (checkpoint 12b). Part D unified every copy of the
  L·atm gas constant onto one value derived in `PyOMES/units.py`, renaming
  `core.phases.R_L_ATM_MOL_K` to `PyOMES.units.R_L_ATM_PER_MOL_K` with no
  alias — the one deliberate numerical change in the phase (~4e-7 relative on
  gas-liquid quantities; BSM2 sentinels re-baselined with a dated note; ADM1's
  own rounded copy of `R` left alone on purpose, pending word on whether it
  matches a published spec). Findings logged rather than fixed along the way
  — the unwired `AccuracyMonitor.check_ionic_strength`, duplicated
  ionic-strength/charge tables across the three engines, the package root
  exporting only the Bisection engine, and others — see `OPEN_WORK.md`. Full
  suite green post-merge: 2072 passed, 0 failed. Tag
  `chemical-equilibrium-engines-subfolder-shipped`. See
  [`../shipped/CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER_CHECKLIST.md`](../shipped/CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER_CHECKLIST.md).

- `tutorials-followups` (2026-09-17) — five debugging/investigation follow-ups
  split out of `tutorials-reorg` so that phase could ship as one clean unit.
  Fixed `reactions/chemistry_database.py` (`ThermoFramework` now takes
  `liquid_activity=` not string kwargs) and `partition_model.py`
  (`HenryPartition`/`.beta()` deprecated in favor of
  `HenryEquilibrium`/`.partition_ratio()`), converting both to notebooks;
  fixed `D2C_workshop/raw_construction.py`'s pH runaway (missing phosphate
  ladder meant the `PHController`'s acid corrector was a silent no-op) and
  its `ConservationWarning`s (small-scale sparged system hitting the
  monitor's absolute per-step clamp); fixed `tests/validation/speciation/
  {07,08}_iron_oxidation*.ipynb`'s Singer-Stumm rate law (the literature
  constant is calibrated against `p(O2)` in atm, but `rate_fn` only ever
  sees aqueous `[O2]` — a genuine, literature-verified physics bug, not a
  demo-pacing choice, since these notebooks live in `tests/validation/`)
  and added `test_iron_oxidation.py` (previously zero coverage); removed
  the orphaned `demos/usecases/03_cstr_dilution_rate_sweep.ipynb`
  duplicate. Surfaced three core-package findings along the way, split into
  their own design notes rather than fixed here (see "Design discussions"
  above): `PHCONTROLLER_CORRECTOR_VALIDATION.md`,
  `REACTION_ENVIRONMENT_PHASE_EXPOSURE.md`,
  `SCIPY_REJECTION_CHECK_SOLVER_AWARENESS.md`. Full suite green post-merge:
  2038 passed, 36 skipped. Tag `tutorials-followups-shipped`. See
  [`../shipped/TUTORIALS_FOLLOWUPS_CHECKLIST.md`](../shipped/TUTORIALS_FOLLOWUPS_CHECKLIST.md).

- `tutorials-reorg` (2026-09-17) — reorganized `demos/` examples into
  topic-based `docs/tutorials/` subdirectories, and split
  validation/performance content into `tests/validation/` and
  `tests/performance/`. 10 checkpoints plus a pre-ship `pyproject.toml`
  `testpaths` fix (CI runs bare `pytest`, and the new validation tests
  weren't on any configured testpath). Full suite green post-merge: 2029
  passed, 36 skipped. Tag `tutorials-reorg-shipped`. See
  [`../shipped/TUTORIALS_REORG_CHECKLIST.md`](../shipped/TUTORIALS_REORG_CHECKLIST.md).

- `stirred-tank-template` (2026-09-15) — moved
  `models/vlmodels/fermenter/` into `PyOMES` core as
  `PyOMES/templates/stirred_tank/`, renaming `FermenterBuilder`/
  `FermenterFactory` to `StirredTankBuilder`/`StirredTankFactory` (no
  phase qualifier — kept unqualified since `PyOMES.core.phases.SolidPhase`
  already exists and a third phase could plausibly be added to this same
  class later; the gas+liquid-only constraint is now stated explicitly
  in both class docstrings instead). `kinetics.py` stays under
  `stirred_tank/` rather than merging into `PyOMES/kinetics/` — the two
  "kinetics" concepts share a name, not an interface (`GrowthKinetics`
  scalar rate laws vs. `KineticModel`'s self-integrating RHS); real
  unification logged as a deferred open question. Both `adm1/base.py`
  and `adm1/bsm2.py` (functionally dependent, not just cosmetic
  references) updated with no numerical drift (BSM2 sentinel tests
  confirmed). Surfaced and logged (not fixed — separately scoped)
  broader pre-existing documentation staleness in `README.md` and
  `docs/architecture.md` predating this phase by several shipped
  phases — see the new [../OPEN_WORK.md](../OPEN_WORK.md). Final suite:
  2011 passed, 28 skipped, 0 failed. Tag
  `stirred-tank-template-shipped`. See
  [../shipped/STIRRED_TANK_TEMPLATE_MIGRATION.md](../shipped/STIRRED_TANK_TEMPLATE_MIGRATION.md)
  and
  [../shipped/STIRRED_TANK_TEMPLATE_MIGRATION_CHECKLIST.md](../shipped/STIRRED_TANK_TEMPLATE_MIGRATION_CHECKLIST.md).

- `step-solver-interface-refinement` (2026-07-08) — Stage 1 of the
  solver-interface-refinement plan (see below). Ownership-guard
  `OrchestrationWarning` on `cv.advance()`; `MonolithicODESolver` rejects a
  per-CV `solver=` it can never consult; `SequentialAdvanceSolver` reified
  (`solver=None` is now sugar for a real class); both `Simultaneous*` solvers
  generalized off the gas/liquid assumption (only `"liquid"` required now;
  every registered `internal_interfaces` entry contributes, not just the
  first); shared swappable `clamp_fn` module
  (`proportional_clamp`/`floor_clamp`, `clamp_fn=None` disables clamping);
  `negative_mole`/`clamp_invoked` `AccuracyMonitor` diagnostics; unified
  `StateVector` class replacing three near-duplicate packing
  implementations; a custom-`StepSolver`/`SystemSolver` demo. `_StateVector`
  deleted (no shim). 1950 → 2003 tests. Tag
  `step-solver-interface-refinement-shipped`. See
  [../shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md](../shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md)
  and
  [../shipped/STEP_SOLVER_REFINEMENT_CHECKLIST.md](../shipped/STEP_SOLVER_REFINEMENT_CHECKLIST.md).

- `simulation-class` (2026-05-27) — third and final of the three
  sequenced phases. New `Simulation` class as the single
  orchestration pathway for all CV-based models. Subsumed
  `MultiCVSystem` (deleted); absorbed the `run_batch` time-loop
  body; promoted controllers and profiles to first-class concepts
  with structured per-step records; ported the legacy `src/control/`
  controller machinery onto CV-shaped state via the Pattern 1
  collapse (`compute(state, dt_h) -> ControlAction` — no
  intermediate `Commands` or separate `Actuator`); added
  RunContext-based lifecycle gating with Pattern B unchecked
  setters for the orchestrator-mediated mutation path. Final
  suite 977/0; tag `simulation-class-shipped`. See
  [../shipped/SIMULATION_CLASS.md](../shipped/SIMULATION_CLASS.md)
  and
  [../shipped/SIMULATION_CLASS_CHECKLIST.md](../shipped/SIMULATION_CLASS_CHECKLIST.md).

- `state-unification` (2026-05-22) — second of the three
  sequenced phases. Collapsed speciation-shaped type machinery
  (`PropertySolver` / `PropertyResult` / `SpeciationPropertySolver`
  / `chem_env` / `ReactionSet.partition` / `Reaction`
  kind-branching) and unified `Phase` storage on PHREEQC-style
  `n_mol` (totals + derived). `cv.reaction_system` is the single
  reaction attach point; the new narrow `PropertyCalculator`
  protocol replaces `PropertySolver`; the engine writes derived
  species back to `n_mol` via `phase._refresh_derived(values)`. C6
  adds `ConservationMonitor` (element + charge accounting per
  step). Final suite 793 / 0; tag `state-unification-shipped`.
  BSM2 sentinels re-baselined across the phase (pH 3.340171 →
  3.303580). See
  [../shipped/STATE_UNIFICATION.md](../shipped/STATE_UNIFICATION.md)
  and
  [../shipped/STATE_UNIFICATION_CHECKLIST.md](../shipped/STATE_UNIFICATION_CHECKLIST.md).

- `integrator-removal` (2026-05-20) — deleted the `Integrator`
  Protocol + `RK4Integrator` + `EulerIntegrator` (~150 LOC
  subtraction); the sequential `cv.advance()` body's reaction
  sub-step collapses to a single forward Euler evaluation;
  `StepSolver` becomes the sole user-facing extension surface for
  time integration. BSM2 golden re-baselined (largest drift S_h2
  at ~6e-3 relative; pH 1.6e-7 relative). 816 standalone tests
  passing (down from 821 by 5 deleted tests on the removed
  surface). First of three sequenced phases (resequenced
  2026-05-20): `INTEGRATOR_REMOVAL` → `STATE_UNIFICATION` →
  `SIMULATION_CLASS`. See
  [../shipped/INTEGRATOR_REMOVAL.md](../shipped/INTEGRATOR_REMOVAL.md)
  and
  [../shipped/INTEGRATOR_REMOVAL_CHECKLIST.md](../shipped/INTEGRATOR_REMOVAL_CHECKLIST.md).

**Chemistry unification is fully complete.** Six precursor / phase
branches shipped, followed by the trigger-gated `chemistry-unification-3b`
(shipped 2026-06-03). All seven branches:

- `reaction-protocol-cleanup` (2026-05-12) — see
  [../shipped/REACTION_PROTOCOL_CLEANUP_CHECKLIST.md](../shipped/REACTION_PROTOCOL_CLEANUP_CHECKLIST.md).
- `bsm2-reference-test` (2026-05-12) — golden-trajectory regression
  for the BSM2 reference model so chemistry-unification numerics
  drift is visible; see
  [../shipped/BSM2_REFERENCE_TEST.md](../shipped/BSM2_REFERENCE_TEST.md).
- `chemistry-unification-1` (2026-05-13) — declared equilibrium
  reactions consumed by the speciation engine, `chem_env` shrunk to
  its minimal contract, stateless snapshot model. See
  [../shipped/CHEMISTRY_UNIFICATION_1_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_1_CHECKLIST.md).
- `chemistry-unification-2` (2026-05-15) — `PropertyResult.alphas`
  channel; gas-liquid link collapses inline `f_molecular` to a single
  alpha read; `SpeciationCorrection` cascade deleted; dormant
  `_check_f_molecular_consistency` removed. Leak 1 (and Leak 3)
  closed. See
  [../shipped/CHEMISTRY_UNIFICATION_2_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_2_CHECKLIST.md).
- `chemistry-unification-3` (2026-05-15) — canonical-name dedup in
  `_compute_species_eq`; suffix-based `ionic_strength_from_speciation`
  (unmasks a latent NH4+/I bug introduced by Phase 1, BSM2 sentinels
  re-baselined); cross-phase equilibrium reaction infrastructure
  (`Reaction.is_cross_phase`, `log_K=None` for cross-phase, engine
  filters them silently, `KineticGasLiquidLink.derive_speciation_keys`
  wired into `ControlVolume.__init__`). Leak 2 infrastructure shipped;
  BSM2/ADM1 builder migration deferred to Phase 3b alongside the
  unified Species-ID emission convention. See
  [../shipped/CHEMISTRY_UNIFICATION_3_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_3_CHECKLIST.md).
- `chemistry-unification-4` (2026-05-16) — `AccuracyMonitor`
  attached to every `ControlVolume`; `AccuracyWarning(UserWarning)`
  category; `WarningConfig` thresholds + throttle on
  `VLsim.config.warnings`; `VLSIM_WARNINGS` env-var with flat
  presets. Six cheap per-step checks (pH change, Newton iters,
  charge residual, dt vs τ_min stub, scipy step rejections,
  ionic-strength regime). The legacy `_warned_high_I` set on
  `ChemicalEquilibriumEngine` is gone — its ionic-strength threshold
  warning migrates to `AccuracyWarning` via the property solver.
  See
  [../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md).

- `chemistry-unification-3b` (2026-06-03) — `ChemistryDatabase` rollout,
  unified Species-ID emission, BSM2/ADM1 cross-phase migration. Tag
  `chemistry-unification-3b-shipped`. See
  [../shipped/CHEMISTRY_UNIFICATION_3B_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_3B_CHECKLIST.md).

- `partition-model` (2026-06-16) — `PartitionModel` protocol +
  `HenryPartition`; `ChemistryDatabase.partition_models`; H₂S alpha
  correction (2× error at pH=pKa fixed); `_HENRY_PARAMS` + `henry_mol_L_atm()`
  deleted; `vfa_volatility` flag removed from ADM1 builders. 1390→1428 tests.
  Tag `partition-model-shipped`. See
  [../shipped/PARTITION_MODEL.md](../shipped/PARTITION_MODEL.md).

## Priority order (historical — all shipped)

> All phases listed below have now shipped. The priority list is
> preserved as a history of scope decisions and trigger rationale.
> For work that is still open, see "Design discussions" and "Open phases"
> in this file, and [OPEN_WORK.md](../OPEN_WORK.md) for smaller follow-up
> items not yet scoped as phases.

1. **[../shipped/CHEMISTRY_UNIFICATION.md](../shipped/CHEMISTRY_UNIFICATION.md)** (design)
   + **[../shipped/CHEMISTRY_UNIFICATION_PLAN.md](../shipped/CHEMISTRY_UNIFICATION_PLAN.md)**
   (plan) — unify the `Reaction` and speciation declaration surfaces
   so equilibrium reactions become a mode of the existing reaction
   framework. *All seven branches shipped (`chemistry-unification-1`
   through `chemistry-unification-4` plus `chemistry-unification-3b`).*
   *Shipped fully 2026-06-03.*

2. **[../shipped/RUN_HISTORY.md](../shipped/RUN_HISTORY.md)** —
   `StreamingFileRecorder` (Parquet chunks + `load_run`), `SparseRecorder`,
   and `SummaryRecorder` + `SummaryResult`. Built on the `Recorder` /
   `BatchRecorder` / `BatchResult` schema shipped in `simulation-class`.
   *Shipped 2026-06-08; tag `run-history-shipped`.  1174→1226 tests.*

3. **[../shipped/SPECIATION_LEVEL_RETIREMENT.md](../shipped/SPECIATION_LEVEL_RETIREMENT.md)** —
   retired the `level` (1, 2, 2.5, 3) parameter on `ChemicalEquilibriumEngine`
   and the three parallel solver implementations. Deleted 6 level/adapter
   files (2,448 lines); unified engine routes directly to
   `solve_from_equilibrium_set` or `solve_acid_base`.
   *Shipped 2026-06-03; tag `speciation-level-retirement-shipped`.*

4. **[../shipped/CUFERMENTER_SUNSET.md](../shipped/CUFERMENTER_SUNSET.md)** —
   retired `CUFermentationSpeciation` and dropped the BioSTEAM bridge
   (Option C). Deleted 9,093 lines: `unit.py`, `loops.py`, `sim/`,
   `solvers/`, `actuators/`, `FermenterResult`, and 7 legacy test files.
   `FeedState` survives.
   *Shipped 2026-06-01; tag `cufermenter-sunset-shipped`.*

## Recently surfaced (2026-05-27 audit)

Trigger-gated planning notes added from the
post-`simulation-class` audit. Not yet ranked into the
numbered priority list above — see each note's
*Trigger conditions* section for the criteria that would pull
it in. All four notes share the audit's framing
("fermenter-shaped framework surfaces dressed up as general")
and naturally bundle if a single "generalisation pass" is
taken.

- **[../shipped/FRAMEWORK_POLISH.md](../shipped/FRAMEWORK_POLISH.md)** —
  P1–P9 all done (~200 LOC). Fixed silent bug: `GasFeed.y` param path
  was discarding writes.
  *Shipped 2026-06-01; tag `framework-polish-shipped`.*
- **[../shipped/PARAM_PATH_DISPATCHER.md](../shipped/PARAM_PATH_DISPATCHER.md)** —
  `MutableScalar` / `MutableDict` descriptors; typed `ParamPath` builder;
  reflective walker. All `@property + _set_X_unchecked` triplets removed;
  unknown paths raise `ParamPathError`. Subsumed FRAMEWORK_POLISH P8.
  *Shipped 2026-06-08; tag `param-path-dispatcher-shipped`.*

## Recently surfaced (2026-05-28 exploration session)

Trigger-gated planning notes added from the 2026-05-28
post-simulation-class exploration session (demos / exports /
checkpointing / logging / FBA / chemistry-database
questions). Cross-cutting interaction: the three "recorder
landscape" notes (RESULT_EXPORT, HPC_CHECKPOINTING,
RUN_HISTORY's sub-deliverables) all touch the same surface;
shipping them together gives a coherent crash-resilient
recorder story.

- **[../shipped/DEMO_RESTRUCTURE.md](../shipped/DEMO_RESTRUCTURE.md)** —
  *Shipped 2026-05-29 with deviations from option β; see the
  reconciliation banner at the top of the planning note for the
  actual layout.* Four fermenter demos moved to `demos/builder/`;
  new `demos/model_api/` umbrella holds `chemistry/` and
  `model_construction/`. Q1 / Q2 demos shipped alongside; Q7
  chemistry-database closed by `chemistry-unification-3b`.
- **[../shipped/RESULT_EXPORT.md](../shipped/RESULT_EXPORT.md)** —
  `BatchResult.to_dataframe / to_csv / to_parquet / to_wide`; 5 NumPy
  channels; pandas + pyarrow optional; demo at
  `demos/model_api/export_results.py`.
  *Shipped 2026-06-01; tag `result-export-shipped`.  32 tests.*
- **[../shipped/HPC_CHECKPOINTING.md](../shipped/HPC_CHECKPOINTING.md)** —
  `sim.save_checkpoint(path, mode="data+env")` API with four
  fidelity levels (`data` / `data+env` / `source` / `full`); `data`/`data+env`
  implemented, `source`/`full` raise `NotImplementedError`.
  *Shipped 2026-06-08; tag `hpc-checkpointing-shipped`. 78 tests.* This
  entry previously read "not yet implemented" — stale; corrected
  2026-07-09 during a scoping-review hygiene pass after cross-checking
  against `git log` found the tag already merged to `main`.

## Recently surfaced (2026-06-11 solver architecture)

Six design docs added from the 2026-06-11 solver architecture session.
Phases A–E all shipped 2026-06-11/12; docs moved to `shipped/`.
See **[SOLVER_ARCHITECTURE.md](../../ideas/SOLVER_ARCHITECTURE.md)** for the anchor doc
and Phase F+G placeholders (SUNDIALS opt-in, no work scheduled).

- **[../shipped/CV_COMPUTE_INTERFACE.md](../shipped/CV_COMPUTE_INTERFACE.md)** — Phase A. *Shipped 2026-06-11.*
- **[../shipped/CONTROLLER_STATE_PROTOCOL.md](../shipped/CONTROLLER_STATE_PROTOCOL.md)** — Phase B. *Shipped 2026-06-12.*
- **[../shipped/SYSTEM_SOLVER_PROTOCOL.md](../shipped/SYSTEM_SOLVER_PROTOCOL.md)** — Phase C. *Shipped 2026-06-12.*
- **[../shipped/IMPLICIT_TRANSPORT.md](../shipped/IMPLICIT_TRANSPORT.md)** — Phase D. *Shipped 2026-06-12.*
- **[../shipped/MONOLITHIC_ODE.md](../shipped/MONOLITHIC_ODE.md)** — Phase E. *Shipped 2026-06-12.*

## Solver interface refinement (multi-stage, Stage 1 shipped)

Multi-round design review of `StepSolver`/`SystemSolver` customizability,
triggered by auditing how flexible `cv.advance()`/`Simulation` solver
selection actually is. Extends the shipped `SOLVER_ARCHITECTURE.md` two-axis
design (Phases A–E) — not the SUNDIALS/DAE Phase F/G track.

- **[../shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md](../shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md)** —
  the design note. 9 findings, 8 shipped as Stage 1 (below); item 5
  (composition/integrator/clamp_fn factoring) explicitly deferred; the
  interleaving/multi-CV-aware solver tier discussion is Stages 2–4 of the
  plan, not yet resolved.
- **[STEP_SOLVER_REFINEMENT_PLAN.md](STEP_SOLVER_REFINEMENT_PLAN.md)** —
  the sequencing plan (Stages 0–4). Stage 0 (rename) and Stage 1 (the
  bundle) shipped; Stage 2 (z-staleness design decision, §12 Q7), Stage 3
  (reactive D_eff transport) and Stage 4 (SIA, blocked on Stage 2) still
  pending — see the plan doc for current status.
- **[../shipped/STEP_SOLVER_REFINEMENT_CHECKLIST.md](../shipped/STEP_SOLVER_REFINEMENT_CHECKLIST.md)** —
  Stage 1's full checkpoint-by-checkpoint implementation log (10
  checkpoints, all shipped), including several scope corrections found
  during implementation. Archived in `shipped/` now that Stage 1 has
  shipped; still worth reading for decisions (e.g. the checkpoint 4/7b
  split) relevant to Stages 2–4 still to come.

## Recently surfaced (2026-07-01 Layer 1 gap closure)

Both phases shipped 2026-07-02/03; docs moved to `shipped/`.

Design discussion 2026-07-01 (see `MASS_EXCHANGE_ARCHITECTURE.md` §14 and
`CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md` §3/§18). Folds gas-liquid VLE
(and, per §14, recognizes Ksp as already folded) into the same simultaneous
Newton system as acid-base, closing the "⚠ gap" row in
`MASS_EXCHANGE_ARCHITECTURE.md` §10.4. Reconciles three previously-conflicting
precipitation design threads (`NR_PRECIPITATION_SPECIATION`,
`NR_PRECIPITATION_CV_INTEGRATION`, `MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN`)
without changing any of their shipped/planned numerics.

- **[EQUILIBRIUM_CONSTRAINT_UNIFICATION.md](../shipped/EQUILIBRIUM_CONSTRAINT_UNIFICATION.md)** —
  Phase 1: **shipped 2026-07-02** (tag `equilibrium-constraint-unification-shipped`).
  Declaration-API generalization: `EquilibriumConstraint` protocol +
  `HenryEquilibrium`/`KspEquilibrium`/`RaoultEquilibrium` sibling types (not
  alternate `EquilibriumReaction` constructors — a mid-flight refinement);
  single `classify_equilibrium_constraint()` auto-classification path
  replacing the `is_cross_phase` silent-filter and `precipitation_reactions=`
  kwarg patchwork; fixes the `HenryPartition`/reaction double-declaration bug
  (`ReactionSystem` bucketing reworked to reach it); per-entry
  activity-dispatch helper (groundwork for Phase 2). No solve-time behavior
  change; no rename.

- **[LAYER1_GAP_CLOSURE.md](../shipped/LAYER1_GAP_CLOSURE.md)** —
  Phase 2: **shipped 2026-07-03** (tag `layer1-gap-closure-shipped`). Actual
  Newton-system folding: gas-liquid (Henry) and Raoult/evaporation rows enter
  `g(y,z)=0`; `step_internal_transfer()` scope-filter fix (found no filter
  existed at all beforehand, not a partial one); `WaterVapourBoundary`/
  `VentWaterLoss` retired; precipitation regression verification (numerics
  unchanged, but surfaced a real "CO₂ stripping promotes CaCO₃ scaling"
  effect); closes with the `SpeciationEngine` → `ChemicalEquilibriumEngine`
  rename, plus the `src/speciation/` → `src/chemical_equilibrium/` module
  rename (§18's own deferred decision, executed at this phase's close).

## Open phases

### NR Precipitation (two-phase sequence, Phase 1 shipped)

Design discussion 2026-06-23. Active-set precipitation equilibrium built on top
of the NR speciation engine shipped in `nr-speciation-engine`.

- **[NR_PRECIPITATION_SPECIATION.md](../shipped/NR_PRECIPITATION_SPECIATION.md)** —
  Phase 1: **shipped 2026-06-23**. Speciation layer only. Outer active-set loop in
  `NRChemicalEquilibriumEngine.solve()`; `precipitation_equilibria` bucket on
  `ReactionSystem`; `element_stoichiometry` cross-component mass balance fix;
  `"minerals"` key in output dict; `Ca_plus_plus` / `Mg_plus_plus` in
  `common_species`; `05_precipitation_equilibrium.ipynb` demo. No CV changes.

- **[NR_PRECIPITATION_CV_INTEGRATION.md](NR_PRECIPITATION_CV_INTEGRATION.md)** —
  Phase 2: CV/SolidPhase integration. `_read_from_phases` sums solid
  contribution; `SolidPhase` writeback; consistent phase-type validation across
  all three phase types on `ControlVolume`. Phase 1 has shipped, so this is
  unblocked; not yet started. **Needs re-derivation before pickup:** the note
  predates `EQUILIBRIUM_CONSTRAINT_UNIFICATION` and `LAYER1_GAP_CLOSURE`
  (which changed how precipitation is classified and folded into the
  solve), and its file paths use the old `src/` layout — re-check its design
  against the current code first.

## Shelved

Work that was started and then deliberately paused. These are not open
phases — re-decide before picking any of them up.

- **[MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md](MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md)** —
  **shelved 2026-09-20.** A separate, standalone-by-design engine for
  coupled multi-metal / citrate / phosphate chemistry that `NRTableau` can't
  represent. A partial skeleton (7 modules, 353 lines, only `components.py`
  tested) was committed directly to `main` on 2026-06-29 without a branch or
  checklist. Nothing else in the package used it, so it was **removed from
  the package on 2026-09-20**; it can be restored from commit `69d517a` (the
  plan's banner has the exact command). Resume when a concrete model needs
  this chemistry. Its Phase 3 (a persistent solid-phase adapter) overlaps
  `NR_PRECIPITATION_CV_INTEGRATION` above, so design the two together.

---

## How to start one

1. Read the relevant note in this folder. It captures the design
   thinking and any decisions already pinned.
2. Write a checklist file (e.g. `PHASE6_CHECKLIST.md`) modelled on
   [../shipped/PHASE5_CHECKLIST.md](../shipped/PHASE5_CHECKLIST.md).
   Drop it into this folder while it is being worked.
3. Create a feature branch off `main` named for the phase
   (e.g. `chemistry-unification`).  All implementation lands on that
   branch.
4. Push the branch periodically so the work is backed up to GitHub
   while in progress.
5. When the work ships, move both the original note and the
   checklist to [../shipped/](../shipped/), add a
   "Status: Shipped" banner to each, and update the priority list
   above to remove the entry.  Then follow the branching and tagging
   convention below to ship the branch back to `main`.

## Branching and tagging convention

This repo uses a feature-branch-per-phase workflow with explicit
shipping tags.  The convention is:

1. **One branch per phase.**  Create the branch off `main` at the
   start of a phase (e.g. `git checkout -b chemistry-unification`).
   All implementation, doc updates, and intermediate commits land on
   that branch.  Push it to `origin` so the work is backed up.
2. **Merge with `--no-ff` when shipping.**  Once the checklist is
   complete and tests are green, merge the branch into `main` with
   an explicit merge commit:
   ```
   git checkout main
   git merge --no-ff <phase-branch> -m "Merge <phase-branch>: <one-line summary>"
   ```
   The `--no-ff` flag keeps the phase boundary visible in `git log`
   as a fork-and-rejoin shape, even when a fast-forward would
   otherwise be possible.
3. **Tag the shipping commit.**  Create an immovable named pointer
   on the phase's final commit so it can be checked out later by
   name:
   ```
   git tag <phase-name>-shipped <commit-hash>
   ```
   Use names like `phase-7-shipped`, `chemistry-unification-shipped`,
   etc.  These appear on the GitHub "Tags" page and are addressable
   via URLs like `https://github.com/.../tree/<tag-name>`.
4. **Push the branch *and* the tags.**  `git push` ships the merge;
   `git push --tags` ships the tag.  Both are needed — tags are not
   pushed by default.
5. **Delete the feature branch after merge.**  The work now lives on
   `main` and is permanently addressable via the tag, so the branch
   pointer is redundant:
   ```
   git branch -d <phase-branch>
   git push origin --delete <phase-branch>
   ```

This convention preserves every commit as a permanent checkpoint
while keeping `main` the current truth and giving each phase a
human-readable name for browse-back-in-time.
