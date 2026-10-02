# Open work

Standalone follow-up items surfaced during other phases — not yet
scoped as their own phase, no branch, no checklist. Referenced from
[`upcoming/README.md`](upcoming/README.md).

## Findings from the explicit-species-resolution work

Found 2026-09-30 while running the audits in
[`upcoming/EXPLICIT_SPECIES_RESOLUTION.md`](shipped/EXPLICIT_SPECIES_RESOLUTION.md)
and working its checkpoints, not fixed.

- **Core code assumes the gas ids `"O2"`, `"CO2"` and `"N2"`.** The DO sensor,
  `GasFeed`'s and the membrane boundary's default air, the gas-liquid link's
  CO2 defaults, the DO controllers' kLa paths and the vent physics' gas
  tables. Scoped as its own design note:
  [`upcoming/GAS_SPECIES_IN_CORE.md`](upcoming/GAS_SPECIES_IN_CORE.md).
- **Aerobic growth ignores elements it does not balance.** The stirred-tank
  tutorials grow `Yeast` (N 0.16) with `balance_basis="CHO"`, so the biomass
  nitrogen is unaccounted for; a nitrogen source carrying S is not
  balanced for S; a charged source leaves the reaction charged. Scoped as
  its own design note:
  [`upcoming/GROWTH_STOICHIOMETRY.md`](upcoming/GROWTH_STOICHIOMETRY.md).
- **`ControlVolume.snapshot()` drops `chemistry_db`.** It rebuilds the CV
  without it (the species set is carried across as `species=`), so a
  snapshot's `chemistry_db` is `None` and
  `Simulation._warn_thermo_mismatch` skips it. Passing it through would let
  that check run on snapshotted simulations, which may warn where it does
  not today.
- **The conservation monitor does not net out boundary flows.** It compares
  element and charge totals between steps, so gas leaving through a vent or
  material entering with a feed reads as drift ("a kinetic reaction or
  boundary that does not close"). In ADM1's fingerprint case the N2 vented
  in one step (1.911e-4 mol N) is reported as N drift, alongside the O and H
  drift from vented CO2 and water. Subtracting each boundary's flux would
  leave only genuine imbalances. Part of [`upcoming/MASS_BALANCE_CLOSURE.md`](upcoming/MASS_BALANCE_CLOSURE.md).
- **The liquid volume never changes (high priority).** `LiquidPhase.apply_flux`
  changes moles only, `LiquidFeed` adds solute without volume, and nothing in
  `PyOMES/` or `models/` sets a phase's `V_L` after construction. So a
  fed-batch run never dilutes: `docs/tutorials/templates/fed_batch_fermenter.py`
  feeds 5 L/h and its liquid is 1600 L at the start and at the end, and every
  concentration it reports is computed on the starting volume. A CSTR with
  matched feed and drain is unaffected. Fixing it means feeds, drains and
  doses carry a volume and the phase's `V_L` follows; solution dosing
  ([`upcoming/DOSING_AGENTS.md`](upcoming/DOSING_AGENTS.md)) depends on it.
  Part of [`upcoming/MASS_BALANCE_CLOSURE.md`](upcoming/MASS_BALANCE_CLOSURE.md).
- **The FBA notebooks declare no species.** `docs/tutorials/reactions/fba/fba_toy.ipynb`
  and `fba_ecoli_core.ipynb` run a `BlackBoxReactionModel`, which exposes no
  stoichiometry, and pass no `species=`, so every metabolite (Glucose, Biomass,
  O2, CO2, Acetate, ...) raises `UnresolvedSpeciesWarning` and is invisible to
  the conservation checks. Declaring them needs formulas (biomass's in
  particular), and the monitor would then check the FBA rates, which are not
  guaranteed to balance elements, so new drift warnings would follow. A
  decision about those demos, not a mechanical fix.
- **ADM1 / BSM2 helpers fall back to biomass for unknown ids.** ADM1's
  `_mw()` returns the biomass MW, and BSM2's `_mw()` / `_thod()` / `_atoms()`
  the biomass MW, ThOD and formula, for any id not in their tables. Every
  current caller passes a table id or an organism, so nothing is wrong
  today, but a misspelt id would be treated as biomass without a word (ADM1's
  `_get_species` now raises instead). ADM1's `_at()` has no caller anywhere
  in the repo (searched `.py`, `.ipynb`, `.md`).

- **The engines never debit or credit solvent water.** The Bisection engine treats
  `H2O` as a fixed solvent (`_SOLVENT_IDS = ("H2O",)` in
  `engines/bisection/engine.py`): water consumed or produced by a declared
  equilibrium such as `CO2 + H2O <-> HCO3- + H+` is not taken from or added
  to `n_mol`, so each HCO3- formed adds 2 H and 1 O. In
  `tests/standalone/test_user_defined_model.py` that is +2.5e-8 relative in H
  and O over 1 h (strict expected failures there). Fixing it moves every
  model's water amounts. Related name-keyed ids in the same engine: the
  solvent is recognised by the id `"H2O"`, the proton by `"H+"`
  (`_H_PLUS_ID`), and the solver reports H+ and OH- under those fixed ids
  whatever ids the model declared. The NR engine writes its own H2O, and all
  three engines solve pH from a charge balance and overwrite (or ignore)
  H+ and OH- in `n_mol`, so H+ or OH- a feed or dose adds leaves the H and
  O books too. Part of [`upcoming/MASS_BALANCE_CLOSURE.md`](upcoming/MASS_BALANCE_CLOSURE.md).
- **BSM2's nitrogen inhibition reads molecular NH3 as total nitrogen.**
  `_I_IN` and `_I_nh3` in `models/vlmodels/adm1/bsm2.py` read
  `concentrations["NH3"]` as S_IN (total inorganic nitrogen), and `_I_nh3`
  then applies the free-NH3 fraction to it again. The speciation engine has
  always written the NH3 / NH4+ split back to `n_mol`, so both see molecular
  NH3 only.
- **Some inventory floors are silent.** The step solvers warn
  (`AccuracyWarning`) when `clamp_fn` scales a reaction, feed or boundary
  flux, but `Phase.apply_flux`'s default `clamp=True` floors at zero without
  a warning, and it is what internal gas-liquid transfer
  (`ControlVolume.step_internal_transfer`), `apply_external_flux` and
  inter-CV links (`Simulation`) use. `SimultaneousAdaptiveSolver` floors
  trial states with `floor_nonnegative`, also without a warning.

- **`build_adm1_cv(..., ethanol=True)` cannot build.** It calls
  `.transfer_species("Ethanol")` without a Henry constant, and the
  anaerobic-digestion database has no Ethanol partition model, so the
  factory raises `ValueError: No partition model for 'Ethanol'`. No test
  or notebook passes `ethanol=True`.
- **The stirred-tank template's pH controller never doses.** A tank built by
  `StirredTankBuilder` carries no acid-base equilibria: the database is used
  only for partition models, and the reaction system is the growth reaction
  alone. So `pH` is NaN at every step and a `PHController` configured on it
  adds nothing. Measured on the fed-batch tutorial's configuration
  (`docs/tutorials/templates/fed_batch_fermenter.py`, 2 h): no Na+ ever
  appears. The batch, CSTR and microplate tutorials configure the same
  controller on the same kind of tank.

## Four small findings from the molar-mass and `plot_vant_hoff` audit

Found 2026-09-30 while checking the molar-mass and `plot_vant_hoff` items
(now resolved), not fixed. Grouped because each is small; they are
unrelated to each other.

- **`species_check.py` understates what a hard conflict is.** The module
  docstring and `check_species_consistency`'s `Raises` section say a hard
  conflict is the same id with different atoms or charge, but the code also
  compares `MW` (as `SpeciesConflictError`'s own docstring in
  `chemistry/species.py` says). The docstrings need `MW` added.
- **`docs/architecture.md` places `vant_hoff_log_K` in `constraint.py`.** The
  `reactions/equilibrium/` line lists it under `constraint.py`, which only
  imports it; it is defined in `thermo/temperature_correction.py` (which the
  same file's `thermo/` entry already lists). `constraint.py`'s own function is
  `constraint_log_K_at`.
- **`plot_speciation` has no pytest coverage.** Its only callers are
  `ReactionSystem.plot_speciation` and four notebooks (D2C Example1 and
  Example2, validation 07 and 08), none collected by pytest. CI installs the
  `[test]` extra, which has no matplotlib, so a test would need
  `pytest.importorskip("matplotlib")` or matplotlib added to that extra.
  `_build_ladder_from_system` needs no matplotlib and could be tested on its own.
- **`ATOMIC_WEIGHTS` labels itself IUPAC 2021 but uses S = 32.065.** That is the
  older (2007) value; the 2021 conventional value is 32.06. Changing it shifts
  the auto-computed MW of every sulfur species by 0.005 g/mol per S atom, so it
  needs a check of what that moves (sulfate and sulfide species, BSM2) before
  either the value or the label is changed.

## `03_phreeqc_engine_basics.ipynb` describes a warmstart drift the engine no longer has

Found 2026-09-28 when re-running the notebook after renaming `phreeqc_to_pyomes`,
not fixed. `PHREEQCChemicalEquilibriumEngine` builds a fresh PHREEQC solution
for every `solve()`; `use_warmstart` is kept only so existing callers can still
pass it. The notebook's Section 4 ("Gotcha: warmstart drifts from a fresh
solve"), the `reset_cache()` note in Section 5 and the "Warmstart gotcha" row in
the summary still describe the older incremental `sol.change(...)` path, and
three stored outputs come from it. A re-run prints a warm/cold `|diff|` of
0.0000 at every total instead of the stored 0.02–0.09, and the Section 2 solve
gives pH 6.1331 instead of the stored 7.2865, with different species
concentrations in Section 3. Fixing it means rewriting those markdown cells to
match the engine and re-saving the outputs.

## `vlmodels` still carries the former project name

Found 2026-09-28 while removing the former project name `vlsim` from the
package, tests and tutorials; the models package was left out of that change.
The "vl" in `vlmodels` (`models/`) is the same prefix. Renaming it means moving
`models/vlmodels/` (a plain filesystem move, staged as `CLAUDE.md` describes),
changing the package name in `models/pyproject.toml` and `models/setup.py`, the
package's own references (one import in `adm1/bsm2_direct.py`, docstring
examples in `adm1/` and `hplc/`), the imports in
`tests/standalone/test_bsm2_reference.py`, `test_hplc_column.py` and
`test_simulation.py`, the path keys in `test_gas_constant_single_source.py`, a
docstring in `PyOMES/numerics/spatial.py`, and the mentions in `README.md`,
`models/README.md`, `docs/architecture.md`, `CLAUDE.md` and the open planning
docs; then reinstalling with `python -m pip install -e ./models`. The HPLC
template plan ([`upcoming/HPLC_CV_GRAPH_TEMPLATE.md`](upcoming/HPLC_CV_GRAPH_TEMPLATE.md))
retires `models/vlmodels/hplc/`, which would shrink the rename if it ships first.

## Property calculators run only under `SequentialAdvanceSolver`

Found 2026-09-25 while re-checking `docs/architecture.md` against the code, not
fixed. `ControlVolume._run_property_calculators()` has one caller in the repo,
`SequentialAdvanceSolver.solve_step` (`core/solvers.py:151`).
`SimultaneousEulerSolver` and `SimultaneousAdaptiveSolver` never call it, and
neither does `MonolithicODESolver`, which integrates every CV through
`cv.compute_rhs()` instead of `advance()`. Under any of those three,
`phase.properties` keeps whatever was last written (or stays empty), so a rate
law reading `env.prop("viscosity")` sees a stale or missing value. The other
system solvers are affected only through the step solver they hand each CV to.
A fix would call the property calculators from the simultaneous solvers' state
snapshot (and from `compute_rhs`), which changes results for any model that
registers a calculator and uses one of those solvers.

## A CV's `chemistry_db` activity model never reaches its speciation engine

Found 2026-09-25 while scoping the `activity_model` parameter change, not fixed.
Every `ChemistryDatabase` carries a `ThermoFramework` whose `liquid_activity` is
the liquid activity model, and `ControlVolume(chemistry_db=...)` stores the
database on `cv.chemistry_db`. Nothing passes that model on:
`ReactionSystem.engine` builds its engine only from the reaction system's own
engine settings (`configure_engine`, default ideal), and the stirred-tank factory
reads only the database's partition models. So a CV given a database whose
framework uses Davies or SIT still solves ideal chemistry, with no warning. The
only reader of the database's framework is `Simulation`, which warns when two
linked CVs' frameworks differ, as if the difference mattered to the chemistry.

No result is wrong today: all three stock databases (`AQUEOUS_DEFAULT`,
`AD_BASIC`, `BIOPROCESS_BASIC`) are ideal, and the only CV in the repo built with
a non-ideal database is in `TestThermoMismatchWarning`
(`tests/standalone/test_chemistry_database.py`), which has no reactions and only
checks the warning. The likely fix is to make the CV's database model the
default for its engine, with an explicit `ReactionSystem` setting taking
precedence, and to test both. That changes results for any CV whose database is
non-ideal and whose reaction system sets nothing, so it is its own change.

## Regenerating tutorial notebooks can wipe baked outputs

Found 2026-09-16, still applicable while the generator scripts exist
(`docs/tutorials/ArXiv_preprint/_generate_notebooks.py`,
`tests/validation/speciation/_generate_notebooks.py`). A full regenerate is not
safe to commit as-is: `02_kinetic_co2_equilibration_microplate_well.ipynb`
ships with real executed outputs (plots, timing numbers), and regenerating
silently replaces them with the blank `execution_count: null` / no-output
template. The script's existing replace logic also mangled at least one
unrelated string ("VLsim" in a §7 note became "PyOMES"). Prefer
hand-editing the specific stale cell in the committed `.ipynb` over
regenerating, and check per notebook whether it carries baked outputs.
Moot once [`upcoming/NOTEBOOK_GENERATOR_REMOVAL.md`](upcoming/NOTEBOOK_GENERATOR_REMOVAL.md)
ships.

## SIT converts mol/L to mol/kg-water inline, without the bad-density fallback

Surfaced 2026-09-20 during `chemical-equilibrium-engines-subfolder`
checkpoint 4 as three copies of the conversion "divide by water density in
kg/L" in `PyOMES/thermo/`. Those three now share one helper:
`thermo/liquid/water_properties.water_kg_per_L(T_K)`, which returns the
density in kg/L and falls back to `1.0` when the density correlation gives a
non-finite or non-positive value (far outside 0–100 °C).
`ionic_strength_molal_from_molar` and the Davies and SIT
`jacobian_dgamma_dx` methods use it (since 2026-09-24,
`thermo-subfolder-structure`; bit-identical to the three copies it replaced).

Two more copies remain, both in `thermo/liquid/sit.py`:
`SITLiquidModel.gamma_all` and `SITLiquidModel.compute_gammas` each compute
`water_density_kg_per_m3(T_K) / 1000.0` inline to turn concentrations into
molalities for the ion-pair (ε) sum, with no fallback. At a temperature where
the density is negative (for example 1500 K or 2273 K) the ε sum therefore
uses negative molalities, while the Debye–Hückel term, which gets its ionic
strength from `ionic_strength_molal_from_molar`, uses 1 kg/L; the two terms of
the same γ disagree. Inside 0–100 °C the result is the same either way.
Switching both to `water_kg_per_L` would change results only at such
temperatures, so it was left out of the pure-refactor phase; it is a one-line
change in each method if wanted.

## Let a simulation set the value of physical constants such as `R`

Raised 2026-09-20 while planning Part D of
`chemical-equilibrium-engines-subfolder` (one source of truth for the gas
constant in `PyOMES/units.py`). Part D makes the constant *single*; this entry
is the follow-up that would make it *choosable per simulation*, from the
simulation file itself, without editing the package or monkey-patching a
module global.

**Why someone would want it.** Matching a published benchmark that uses a
rounded value (the ADM1/BSM2 reference implementations in
`models/vlmodels/adm1/` carry `_R_J = 8.31446` today, and the ADM1 spec quotes
`R` in bar·m³/kmol/K); sensitivity studies on the constant; reproducing an older
run made before Part D shifted `R` by about 4e-7 relative; and comparing against
another tool (PHREEQC, BioSTEAM) that uses its own value.

**Why it is not a small change.** Every module currently reads `R` by importing
the float (`from PyOMES.units import R_J_PER_MOL_K`), which binds the value at
import time. Overriding per simulation means consumers must read the value from
something they are given, not from a module global. The consumers include hot
paths: `Phase.pressure` / partial-pressure maths in `core/phases.py`,
`core/boundaries.py`, `core/solvers.py`, `chemical_equilibrium/engines/nr/solver.py`,
`thermo/gas/ideal.py`, `thermo/gas/peng_robinson.py`, `chemistry/partition.py`,
`reactions/equilibrium/interphase.py`, `thermo/framework.py` and
`thermo/temperature_correction.py`.

**Design sketch (not decided).**

- A small frozen `PhysicalConstants` dataclass in `units.py`, with the current
  CODATA values as a module-level default. Overriding `R_J_PER_MOL_K` alone must
  update the L·atm value too, which is trivial once Part D derives one from the
  other (Decision 12). Exact conversions (`PA_PER_ATM`, `L_PER_M3`) are not
  overridable.
- Carry it on the object that already scopes a simulation's thermodynamics
  (`ThermoFramework` is the obvious candidate) or on `Simulation` /
  `ControlVolume`, and thread it to the consumers above. Decide whether it is
  per-simulation or per-control-volume.
- Record the constants used in the run metadata and in `save_checkpoint`, so a
  saved run says which `R` produced it.
- Keep the default path fast: resolve the value once per step or at
  construction, not per call inside an ODE right-hand side.

**Open questions.** Is the scope `R` only, or a bundle of constants
(reference temperature `298.15`, `273.15` K offset, `g`)? Is a per-simulation
override actually needed, or is "one correct value, documented" enough? The
trigger should be a concrete model that needs a non-default value, which is
also the point at which the ADM1 rounding question (Part D, Decision 14) gets
answered for real (see the ADM1 / BSM2 entry below). Part D has since landed
(checkpoints 13–14, 2026-09-20): `R` now has one literal in `units.py`, plus the
three deliberate ADM1/BSM2 copies. The wider question of which other constants
would be overridable overlaps with "Sweep the package for each fundamental
constant" below. Do that sweep first, since overriding a constant that still has
scattered copies would only change some of them.

## Development-history references in docstrings and tests outside `chemical_equilibrium/`

Raised 2026-09-20 while scoping checkpoint 12b of
`chemical-equilibrium-engines-subfolder`, which rewrites this kind of
reference inside `PyOMES/chemical_equilibrium/` only. Docstrings and comments
across the package cite the phase or checkpoint that introduced a piece of
code (`CP2 of LAYER1_GAP_CLOSURE`, `chemistry-unification-3b`, "Phase 2"), or
point at a design document by name, instead of describing what the code does
now. They go stale as phases ship and their docs move (`upcoming/` to
`shipped/`), and they tell a new reader about the development process rather
than the code.

**Measured 2026-09-20** (script: lines that name any design doc under
`docs/dev/` or a `CP<n>` / `chemistry-unification-*` / `Phase <n>` / "Layer 1"
label), excluding `chemical_equilibrium/`: **108 lines in 32 files** under
`PyOMES/` — `core/` 67 (14 files), `chemistry/` 14, `thermo/` 8, `templates/` 6,
`reactions/` 5, `control/` 4, other 4. For example `core/control_volume.py:386`
("folded gas-liquid row (CP1/CP2 of `LAYER1_GAP_CLOSURE`)").
**Update 2026-09-24:** `thermo/` is down to 1 line (`thermo-subfolder-structure`
rewrote the other 7 docstring lines when it moved them); the one left,
`equilibrium_constants.py:16`, pointed at this file on purpose. **Update
2026-09-28:** that module was replaced by `thermo/temperature_correction.py`,
whose docstring has no such pointer.

**Tests have the same problem, and worse in the file names.** Module
docstrings read "Tests for CP2 of LAYER1_GAP_CLOSURE…", and several files are
named after the checkpoint that added them (`test_nr_gas_liquid_cp2.py`,
`test_raoult_h2o_fold_cp4.py`, `test_precipitation_gas_liquid_cp5.py`,
`test_step_internal_transfer_scope_filter.py` cites CP3) rather than after the
behaviour they test. Renaming test files is a separate, mechanical change.

Suggested approach: reuse the rules and survey script from checkpoint 12b (keep
the still-true reasoning, drop the label; keep a design-doc pointer only when it
is the canonical explanation and lives at a stable path; AST-compare each file to
prove no code change), package by package, `core/` first.

**Found while doing checkpoint 12b (2026-09-20), outside `chemical_equilibrium/`
and therefore left alone.** Some of these are not just stale labels but
statements that are now false:

- `EquilibriumSet.bsm2_default()` carried a comment saying the deprecated
  `_HA`/`_A-` fallback "is removed in the PARTITION_MODEL phase", on the
  grounds that its four VFA rows had no `species_refs`. The method has since
  been deleted, and the premise was wrong anyway: the BSM2 model builds its
  `EquilibriumSet` through `from_reactions()`, which gives every entry
  `species_refs`. See "Deprecated `{name}_HA` fallback is unreachable from
  repo code" below.
- `core/gas_liquid_link.py:928,958,966` cite
  `...bisection.acid_base._CANONICAL_NAMES`, which no longer exists.
- `reactions/reaction_system.py:250` cites `EQUILIBRIUM_CONSTRAINT_UNIFICATION CP2`.
- Two saved validation notebooks (`05_precipitation_equilibrium.ipynb`,
  `06_phreeqc_benchmark.ipynb`) contain *captured stderr* quoting the old
  `DeprecationWarning` text with its `EQUILIBRIUM_CONSTRAINT_UNIFICATION CP2`
  label. It is a recorded past run, not a live message; it refreshes the next time
  those notebooks are re-run.

**Found 2026-09-25 (`activity-model-parameter` phase): `PropertyResult` is cited
as if it still existed.** No class of that name exists anywhere in the repo;
derived species and their fractions live in `phase.n_mol`, written by the
speciation engine. It is still named in:

- **Two runtime warnings users can see**, in `core/gas_liquid_link.py`: the
  orphan `molecular_driving_force` check (line 459) and `set_transfer_mode`
  (line 551) both say the link looks alphas up in `PropertyResult.alphas`. The
  second also tells the user to add a `speciation_keys` entry by hand, although
  the link derives those keys from the declared cross-phase reactions; that advice
  was not re-checked.
- **`models/vlmodels/adm1/bsm2.py`'s module docstring** ("Note" section, lines
  28–36): says `make_bsm2_callback` "was removed in chemistry-unification-1" and
  that equilibrium species "live in `PropertyResult` only". Both halves are wrong
  now: the engine writes those species back into `phase.n_mol`.
- **Comments:** `bsm2.py:735`, `adm1/base.py:987-996` (which also calls wiring the
  VFA/H₂S equilibria "Phase 3 work"), and docstrings or comments in
  `core/gas_liquid_link.py` (425, 599, 658), `core/control_volume.py:563` and
  `core/solvers.py` (341, 599, 643). Several of these only say `PropertyResult`
  is gone, which is history rather than a false claim.

The two warnings are code strings, so rewording them is a small code change;
the rest is docstring and comment text.

## `chemical_equilibrium/activity_dispatch.py` was deleted (checkpoint 12c)

Found 2026-09-20 in checkpoint 12b and deleted the same day in checkpoint 12c of
`chemical-equilibrium-engines-subfolder`. `activity_for_entry()` was called only
by its own 12 tests (`tests/standalone/test_activity_dispatch.py`). Its old
docstring said a later phase "wires this dispatch into `build_tableau()`", but
that phase (`LAYER1_GAP_CLOSURE`) shipped without doing so: the NR solver folds
gas-liquid rows into the tableau as gas-phase secondaries and handles
solid-liquid equilibria in the engine's precipitation loop. Same grounds as
`api.py` and `factory.py` (checkpoint 9b): no callers, and the package has no
outside users.

Restorable from the last commit before the deletion, `076f6b4`:
`git show 076f6b4:PyOMES/chemical_equilibrium/activity_dispatch.py`, and likewise
the test file. If entry-level activity diagnostics are wanted later, a batch
function designed for that job (γ for the whole liquid composition computed once,
not once per entry) is a better starting point than reviving this one.

## No engine emits a high-ionic-strength warning

Found 2026-09-20 while reviewing `chemical_equilibrium/activity.py` (checkpoint
12d). `AccuracyMonitor.check_ionic_strength` (`monitoring/accuracy.py`) is meant to
warn when ionic strength leaves the range an activity model is good for, with
config thresholds (`ionic_strength_ideal_threshold` 0.10 mol/L,
`ionic_strength_davies_threshold` 0.50 mol/L) and throttling. Nothing in the
package calls it: `ReactionSystem` sets `engine._accuracy_monitor`, but no engine
reads it, and the only caller is a stub engine in
`tests/standalone/test_accuracy_monitor.py`. Its docstring ("Replaces the
per-instance `_warned_high_I` dedup…") describes a hook that no longer exists.
The older `warn_if_high_ionic_strength` in `activity.py` never had a caller either
and was deleted in 12d. Wiring the monitor into the engines would add new
warnings, so it is a behaviour change and was not done there. Both the Bisection
and NR engines already return `ionic_strength` in `EquilibriumResult`, so the call
could sit in `ReactionSystem` rather than in each engine.

**Update 2026-09-22 (`chemistry-reactions-kinetics-cleanup` audit, 2026-09-21):**
a companion write-only attribute, checked directly against the current code.
`ReactionSystem._conservation_monitor` (`reactions/reaction_system.py`) is set
by `attach_conservation_monitor()`, called from `ControlVolume.__init__`
(`core/control_volume.py:304-312`) with the CV's own `ConservationMonitor`
instance — but `reaction_system.py` never calls anything on it afterward. The
CV keeps a second reference to the *same* monitor object on `cv._conservation_monitor`
and invokes `.check_step(cv.phases)` on that one directly from its own
`advance()` body — confirmed live: it is the source of the
`ConservationWarning`s the test suite already prints. So
`reaction_system._conservation_monitor` is not dead in the sense of doing
nothing (the object it points to is genuinely used), but the attribute on
`ReactionSystem` specifically is: nothing reads it through that name. Not
changed here (Part A/B pure-refactor scope) — logged alongside the
`_accuracy_monitor` case above since both are the same shape of problem
(`ReactionSystem` provides an attachment point an engine or the CV's own
code never reads back through).

## Ionic strength and ion charges are defined in several places

Found 2026-09-20 in checkpoint 12d, not changed there because unifying them
changes numbers. Three implementations of `I = ½ Σ z²c`:

- Bisection: `engines/bisection/ionic_strength.py` infers charge from the trailing
  `+`/`-` tokens of each output key. It reads an id that writes its charge as a
  number, such as `Fe2+`, as +1 (declared: +2; see
  `tests/validation/speciation/test_iron_oxidation.py:78-80`, which uses the NR
  path, so nothing wrong is computed today).
- NR: `engines/nr/solver.py::_ionic_strength` uses declared `Species.charge` plus
  a `_STRONG_CHARGES` table for strong ions. That table is copied three times
  (`engines/nr/engine.py`, and twice in `engines/nr/solver.py`), and
  Bisection's charge-balance residual hard-codes strong-ion charges again
  (`engines/bisection/acid_base.py`, `solve_from_equilibrium_set`).
- `thermo/`: `gamma_all` in `liquid/davies.py` and `liquid/sit.py`
  computes I from a caller-supplied `charge` dict and clamps negative
  concentrations to 0, which the other two do not. SIT also keeps a fixed
  `ION_CHARGES` table for `compute_gammas`, which Bisection calls with
  approximate compositions.

Declared charges cannot simply replace the suffix rule in Bisection: the
deprecated `{name}_HA` / `{name}_A-` keys (emitted only for entries added to
an `EquilibriumSet` by hand without `species_refs`) and `Cation(inert)` / `Anion(inert)` have no
`Species` objects. A replacement would have the emitter return a `{key: charge}`
map alongside the concentrations. The comment at `tests/standalone/test_bsm2_reference.py:185-203`
shows the risk: a species emitted under a key the rule does not recognise
silently drops out of I (a missed `NH4+` under-counted BSM2's I by about 25%).
Bisection also ignores `CT_Cu`, `CT_Fe2` and `CT_MoO4`, which NR handles.
Also open: whether the package root should keep re-exporting
`ionic_strength_from_speciation`, which has no users in the repo besides tests.

## ADM1 / BSM2 use a rounded gas constant (`_R_J = 8.31446`)

Decided 2026-09-20 during the gas-constant unification in
`chemical-equilibrium-engines-subfolder` (Decision 14): leave these copies alone
and flag them here. `models/vlmodels/adm1/base.py:224`, `bsm2.py:328` and
`bsm2_direct.py:272` each define `_R_J = 8.31446`, 3.2e-7 below
`PyOMES.units.R_J_PER_MOL_K` (CODATA 8.31446261815324). It is used only for the
van 't Hoff `Ka(T)` corrections of NH4+ and the acid–base pKa functions.

Nothing in the files or the repo says whether the rounding is deliberate (for
instance, to reproduce a published ADM1 or BSM2 specification), so it was not
changed without knowing. What it costs: the ADM1/BSM2 models compute `Ka(T)` with
a slightly different R from the engines, which take R from `units` via `thermo/`.
The 2026-07-02 consolidation unified the same rounding in
`chemistry/thermo_params.py` and `thermo/framework.py`, re-baselining the BSM2
sentinels by about 1e-9 relative, which is the argument for doing the same here.

To resolve: find out whether the benchmark specification fixes R. If it does, keep
the value but replace the three bare literals with one commented, named constant.
If it does not, import `R_J_PER_MOL_K`, re-baseline the BSM2 reference sentinels
with a dated before/after comment, and remove the three matching entries from
`_KNOWN_COPIES` in `tests/standalone/test_gas_constant_single_source.py`.

## Sweep the package for each fundamental constant

The gas-constant work (checkpoints 13 and 14 of
`chemical-equilibrium-engines-subfolder`) established the pattern: one authoritative
definition in `PyOMES/units.py`, everything else imports it, and a guard test fails
on a new copy. It covers only R. The same sweep is needed for every other
fundamental or reference constant, across the whole package and not just
`core/`; the principle (owner, 2026-09-20) is that all of it should take such
values from one place. Nothing under `core/` imports `units` today, and `units.py`
defines only R (two units), `PA_PER_ATM` and the time factors.

Survey of `PyOMES/` and `models/` (2026-09-20, number tokens matched by value, so
comments and strings are excluded; every hit still needs reading, because a value
that matches may be a fitted-formula coefficient and not the constant):

| Constant | Literals | Files | In `core/` | Notes |
|---|---|---|---|---|
| 273.15 (0 °C in K) | 19 | 12 | 2 | `units.py` has no name for it |
| 298.15 (25 °C in K, the standard-state reference) | 39 | 26 | 5 | Also appears as a default argument in many signatures; whether to centralise it is a design choice |
| 101325 (Pa per atm) | 5 | 3 | 0 | `units.PA_PER_ATM` already exists and is not used |
| 1.01325 (bar per atm) | 2 | 2 | 1 | no name in `units.py` |
| 1e-14 (Kw at 25 °C) | 3 | 3 | 0 | may be reference data, not a constant |
| 18.015 (molar mass of water) | 6 | 5 | 0 | overlaps with atomic-weight data in `chemistry/species.py` |
| 55.5 (mol/L of water) | 2 | 2 | 0 | |
| 9.81, 96485, k_B, N_A | 0 | 0 | 0 | not used as literals |

Not yet surveyed: water density and other water properties, conversion factors such
as 1000 L/m3, `ln(10)`, and the temperature-dependence coefficients of vapour
pressure and Henry constants (these are data, not constants, and belong in
`thermo/`, but should be checked for a second copy). The temperature-correction
*formula* they feed is already single-source: it is written only in
`thermo/temperature_correction.py`, and
`tests/standalone/test_temperature_correction_single_source.py` fails on a new copy.
The two water vapour-pressure parameter sets still differ (see "No single source for
water's physical constants").

Suggested approach, reusing the R guard: for each constant, define it once in
`units.py` (derive dependent values from one literal, as R does), find copies by
value with the number-token scan, read each hit, replace the true copies with
imports, and extend `test_gas_constant_single_source.py` (renamed to cover
constants generally) with a per-constant allowlist that shrinks to empty. Most of
these are value-preserving, since every copy holds the same number; anything that
is not (a copy with a rounded value, as R had) should be its own numerics-changing
commit with a recorded before/after shift.

Scope of the guard: `test_gas_constant_single_source.py` scans only `PyOMES/` and
`models/`. Tests and tutorial notebooks are outside it, and the R work found three
copies of R in saved notebook code cells only by a manual re-search
(`Example1_mtp_well`, `Example2_batch_fermenter`, `02_nr_engine_basics`). A sweep
should extend the scan to `tests/` and to notebook code cells, or accept a
periodic manual re-search.

**Update 2026-09-22 (`chemistry-reactions-kinetics-cleanup`):** this phase moved
or deleted several files the 2026-09-20 survey table counted literals in —
`chemistry/thermo_params.py` (deleted, checkpoint 7), `chemistry/recipe.py` /
`chem_recipe.py` (deleted, checkpoint 10), `templates/stirred_tank/kinetics.py`
→ `reactions/rate_laws.py` (checkpoint 12; now `reactions/kinetic/rate_laws.py`), `PyOMES/equilibria/` →
`thermo/gas_eos.py` (checkpoint 11; now split into `thermo/gas/`), `chemistry/database.py` /
`chemistry/databases/` → `PyOMES/databases/` (checkpoint 13) — so the table's
per-file counts are stale, though the phase was pure-refactor for these
literals (moved, not edited), so the aggregate totals per constant should be
close to unchanged. Spot check: `101325`/`298.15` alone still appear 18 times
across just `chemistry/partition.py` (since split into it and
`reactions/equilibrium/interphase.py`), `thermo/gas_eos.py` (now `thermo/gas/`),
`chemical_equilibrium/engines/bisection/equilibria.py`,
`reactions/kinetic/rate_laws.py` and `PyOMES/databases/*.py` post-move. Re-running the
survey against the new layout is part of picking this sweep up, not done here.

## One notebook does not parse on Python 3.10, and the notebooks edited for the gas constant were not re-run

Found 2026-09-20 while checking `chemical-equilibrium-engines-subfolder` before
merging (a parse of every notebook code cell). The problem predates that phase.

- `tests/validation/speciation/05_precipitation_equilibrium.ipynb`, cell 6, uses an
  f-string form that only Python 3.12 accepts (`f-string: unmatched '['` on 3.10),
  while `pyproject.toml` declares `requires-python = ">=3.10"`.

No test executes a notebook, so CI does not notice. Separately, the
gas-constant change edited eight notebooks and one generator by text (the import,
the identifier, and three literal values) and validated that each parses, but did
not re-execute them. Their saved outputs still show numbers from before `R`
moved by 4e-7 relative, which is below what most cells print but not verified. They
refresh the next time each notebook is re-run; see also "Regenerating tutorial
notebooks can wipe baked outputs" above.

## Tutorial and docstring examples that run but show less than they say

Found 2026-09-25 (`activity-model-parameter` phase), not fixed. Each of these runs
without error but is misleading:

- **The `StirredTankBuilder` module docstring example** (`templates/stirred_tank/builder.py`)
  builds and runs a tank in one call (`build_simulation_and_run`), so it never
  sets a starting biomass, and its final line, `result.liquid_mol["main"]["Yeast"][-1]`,
  is `0.0`: nothing grows. The template scripts seed biomass by setting
  `cv.phases["liquid"].n_mol["Yeast"]` after `build()` (e.g.
  `docs/tutorials/templates/cstr_fermenter.py:145`), which the one-call form cannot
  do. The builder has no method for initial composition.
- **The "Canonical shape" example in `docs/tutorials/templates/README.md`**
  attaches `PHController(setpoint=5.0)` to a tank that declares no acid-base
  equilibria, so it has no speciation engine and its pH is `nan`: the controller
  has nothing to act on. It also seeds no biomass, as above.
- **`docs/tutorials/reactions/reaction_system.ipynb`, markdown cell 10,** says
  "four properties expose the type-specific projections" and lists
  `kinetic_reactions`, `single_phase_equilibria`, `cross_phase_equilibria` and
  `blackbox_models`. There is a fifth, `precipitation_equilibria`, and the saved
  output of the next cell shows it (`0 precipitation eq`).
- **The batch tutorials sparge a vessel with no gas outlet** (found 2026-10-02).
  `docs/tutorials/templates/batch_fermenter.ipynb` and `batch_fermenter.py` call
  `.gas_feed(vvm_min=1.0, ...)` and run with `build_simulation_and_run`, so the
  tank's only boundary is the `GasFeed`. Measured on the notebook's configuration
  over 5 h: headspace O2 goes from 3.35 to 3606 mol and N2 from 12.6 to 14320 mol
  in 400 L, and dissolved O2 and N2 follow (0.39 to 423.5 mol and 0.74 to
  841.3 mol in 1600 L). Without the gas feed all four stay at their starting
  values. The CSTR and fed-batch scripts append a `PressureReliefVent` after
  `build()`, which the one-call form cannot do, and no builder method adds one.
  The notebook also seeds no biomass and its pH is `nan`, as in the first two
  items.

The first two want a decision on what the canonical example should demonstrate
(a seeded batch with declared acid-base chemistry would make both the growth and
the pH controller real); the third is a one-cell prose fix; the fourth belongs
with the first two, since a canonical batch example also needs a vent.

## Keep the engine fingerprint scripts, or freeze golden values for each engine

Found 2026-09-20 in checkpoints 12c and 12d. Every pure-refactor checkpoint of
`chemical-equilibrium-engines-subfolder` was verified partly by a "fingerprint": a
fixed set of engine outputs (NR 814 values, Bisection 265, PHREEQC 182) compared
before and after. The scripts that produced them were throwaway and are not in the
repo, so later checkpoints had no baseline to compare against and were verified by
other means (an import-graph proof, AST identity of the moved code, and a fresh
purpose-built fingerprint for the gas-constant change).

Only BSM2 has committed golden values (`test_bsm2_reference.py`, which covers the
Bisection engine and the gas-liquid link end to end). A small module of frozen
sentinel values per engine, in the same style, would let any future refactor prove
it is bit-identical with one test run and no scratch scripts.

## Line endings are mixed across the repo

Found 2026-09-20 during `chemical-equilibrium-engines-subfolder`. Among tracked
`.py`, `.md` and `.ipynb` files, 180 are stored with CRLF line endings in the index,
137 with LF, and 2 with both. There is no `.gitattributes`, and `core.autocrlf` is
`true` on the Windows machine used, so git prints "LF will be replaced by CRLF"
warnings on nearly every edit, and `git diff --check` reports trailing whitespace
on every CRLF line it touches. The two files that mix both were already mixed
before that phase, and the phase flipped no file between the two.

A `.gitattributes` (for example `* text=auto eol=lf`) plus one commit that
renormalises the tree would remove the noise. The renormalisation commit rewrites
most files' line endings, so it should be its own commit with nothing else in it,
and `git blame` will need `--ignore-rev` for it.

## CI runs only for `main`, on Linux

`.github/workflows/tests.yml` triggers on pushes to `main` and pull requests to
`main`, runs `pytest` on Ubuntu with Python 3.10, 3.11 and 3.12, and has no feature-branch trigger. A
phase branch is therefore only tested locally (here, Windows with Python 3.12)
until it is merged, and the first Linux and 3.10/3.11 run happens after the merge.
Opening a draft pull request against `main` before merging, or adding
`push: branches: ['**']` to the workflow, would surface platform problems before
they reach `main`.

## `models/vlmodels/adm1/bsm2_direct.py` has no importer and no test

Found 2026-09-20 when a whole-repo import audit could not import it. Nothing in the
repo imports the module (only `docs/architecture.md` lists it), no test covers it,
and its docstring example (`from PyOMES.models.adm1_bsm2_direct import
build_bsm2_direct_model`) points to a module path that does not exist. It also
imports its neighbour as `from vlmodels.adm1.bsm2 import ...`, which works only
when `models/` is on `sys.path` (the tests' `conftest.py` puts it there);
`models/` itself has no `__init__.py`, so `import models.vlmodels.adm1.bsm2_direct`
fails with "No module named 'vlmodels'". It is one of the three files that keep the
rounded `_R_J = 8.31446` (see the ADM1 / BSM2 entry above). Decide whether it is a
supported second BSM2 implementation, in which case it needs a test and a correct
docstring, or dead code to delete.

## The package root exports only the Bisection engine

Observed 2026-09-20 in checkpoint 11 of `chemical-equilibrium-engines-subfolder`
and not acted on. `PyOMES/chemical_equilibrium/__init__.py` exports
`BisectionChemicalEquilibriumEngine` (and the `ChemicalEquilibriumEngine` alias),
`EquilibriumResult`, `ionic_strength_from_speciation` and `solve_acid_base`, but not
the newer NR engine or the PHREEQC engine, which need deep imports
(`PyOMES.chemical_equilibrium.engines.nr.engine`, `...engines.phreeqc`). The engines
were not exported from the root before the move either, so re-exporting them would
be new API. Whether the root should export all three, none, or a small engine
selector is an API decision.

## A species-based replacement for the deleted recipe layer

Logged 2026-09-22, `chemistry-reactions-kinetics-cleanup` checkpoint 10
(decision D1). `chemistry/recipe.py` (`SolutionRecipe`), `chem_recipe.py`
(`ChemSpec`/`CHEM_DB`), `types.py` (`AqueousTotals`/`AqueousTotalsUser`), and
the ion/salt maps in `registry.py` (`SALT_DISSOCIATION_MAP`,
`ION_TO_ENGINE_KEY`, `normalize_ion_label`, `map_user_ions_to_engine`) were
deleted outright: a fresh repo-wide search found no consumer of any of it.
Two things made that possible rather than just tempting: no outside user
depended on the API, and its actual output shape — pooled `CT_*` totals — was
never what the current framework consumes (species amounts in `n_mol`), so
nothing downstream would have worked with it even if it had been kept.

If a "weighed salt → initial condition" convenience is wanted again, it
should be species-based rather than pooled-total-based: given a compound name
(resolved against the `Species` passed to the model, which carry molar
masses) and a mass or stock-solution dose, emit species amounts
directly (`{"Na+": n_mol, "Cl-": n_mol, ...}`), not `CT_Na`/`CT_Cl`-style
pooled totals. This also supersedes the recipe-layer open question in
`docs/dev/implementation/upcoming/STRONG_ION_INFERENCE_GENERALIZATION.md`.

## Gas EOS: `partial_pressures_atm` returns two different quantities, and `ThermoFramework.gas_eos` is read by nothing

Logged 2026-09-22, `chemistry-reactions-kinetics-cleanup` checkpoint 11
(decision D4), moved but not fixed. `IdealGasEOS.partial_pressures_atm`
(`thermo/gas/ideal.py`) returns partial pressures (`y_i × P`);
`PengRobinsonEOS.partial_pressures_atm` returns fugacities
(`f_i = y_i × φ_i × P`) — same method name on the shared `GasEOS` interface,
different physical quantity. `test_gas_eos.py` (added the same checkpoint)
pins this distinction, so it won't drift silently, but it is still a real
API smell.

Separately: the ideal-gas law is hard-coded independently in at least eight
places (`core/phases.py`, `core/boundaries.py`,
`chemical_equilibrium/engines/nr/solver.py`, `chemistry/partition.py`,
`reactions/equilibrium/interphase.py`,
`control/cv_loops.py`, `templates/stirred_tank/factory.py`,
`models/vlmodels/headspace.py`) instead of going through
`ThermoFramework.gas_eos`, which is read by nothing in production despite
being a real typed `Optional[GasEOS]` field since that checkpoint. A later
change should split the method into `partial_pressures_atm` (ideal) and
`fugacities_atm` (non-ideal), and route both gas pressure and fugacity
through `ThermoFramework.gas_eos` so the seven hard-coded call sites have
one thing to switch on. Overlaps "Sweep the package for each fundamental
constant" above, since `R` is one of the things duplicated across those
same seven files.

## Mapping-based parameters for several limiting/inhibiting species; composable inhibition for the growth rate laws

Logged 2026-09-22, `chemistry-reactions-kinetics-cleanup` checkpoint 12
(decisions D5/D6). `reactions/kinetic/rate_laws.py`'s `DualSubstrateMonod` takes
exactly one secondary species through four scalar arguments (`secondary_id`,
`Ko`, `secondary_in_mol_L`, `secondary_MW`), and
`ReactionBuilder.monod_aerobic_growth` exposes only an O2 term (`Ko2_gL`,
on the O2 `Species` it is given). Neither extends to a second or third limiting
species (e.g. NH3 or a phosphate source) without a new class or a new
argument.

It would be good to take a mapping instead — e.g.
`{species_id: <half-saturation, units/molar mass, limits-or-inhibits>}` —
applied by the rate law to any number of species, with inhibition as a
feature attached to a growth model rather than a separate class per
growth-law × inhibition-law pair. The maths constrains the design: `Andrews`
and `ContoisAndrews` use the Haldane form,
`μ = μmax · r / (Ks + r + r²/Ki)`, where the inhibition term sits inside the
denominator — not a plain multiplier on the base law (as a factor it would be
`(Ks + r) / (Ks + r + r²/Ki)`, which needs the base law's own `Ks`), so a
composable design must treat Haldane specially (e.g. an optional `Ki` on the
base law), while multiplicative forms (non-competitive `1/(1 + I/Ki)`,
exponential `exp(-I/Ki)`) compose freely with any base law. Inhibition is
also often by a *different* species (product inhibition), which a
species-keyed mapping covers naturally. Open questions: the shape of the
per-species parameters; whether species come from `Species` objects so molar
masses are looked up instead of re-entered; which of the eight laws have a
meaningful inhibition companion (`Tessier`, `Moser` and `Blackman` have no
standard one).

Not blocking: `tests/standalone/test_rate_laws.py` pins `Andrews` and
`ContoisAndrews` against a frozen Haldane-form reference and the other six
laws at fixed points, specifically so a later composable-inhibition redesign
can be checked against these values rather than against vibes.

## Package-level layering: two package cycles remain

Originally logged 2026-09-22 (`chemistry-reactions-kinetics-cleanup`
checkpoint 13, decision D7) as `chemistry` cycling with `reactions`. That
cycle is gone: `chemistry/` now imports only `units`, at any depth
(function-level and `TYPE_CHECKING` imports included), and
`tests/standalone/test_package_layering.py` enforces it.

`tests/standalone/test_import_graph_acyclic.py` checks only the *module*-level
graph, so it does not see the package-level cycles that remain. Counting every
import statement (checked 2026-09-24):

- `control` <-> `core`, at module level in both directions. `control/__init__.py:11`,
  `cv_loops.py:17` and `interfaces.py:34-35` import `core`; `core/boundaries.py:34`,
  `gas_liquid_link.py:96`, `phases.py:25`, `recorder.py:27` and `simulation.py:32`
  import `control`. Each side also has function-level imports of the other.
- `chemical_equilibrium` <-> `reactions`, at function level only.
  `reactions/reaction_system.py:261,273` import the engines;
  `chemical_equilibrium/engines/bisection/engine.py:203`, `nr/engine.py:239` and
  `nr/tableau.py:401` import `reactions`.

Nothing is broken today: the module-level graph is acyclic. A test asserting an
acyclic *package* graph cannot be added until both cycles are resolved. Which
package should sit lower in each pair is not decided; fixing either is its own
phase.

## No single source for water's physical constants

Found 2026-09-24 while moving `RaoultEquilibrium` to `reactions/phase_equilibria.py` (now
`reactions/equilibrium/interphase.py`),
not fixed. `RaoultEquilibrium`'s four water values (`P_sat_ref`, `dH_vap`, `T_ref`,
`C_water_mol_L`) are ordinary constructor arguments, so a caller can already supply
their own, or another solvent via `gas_species`/`liquid_species`;
`tests/standalone/test_partition_model.py::TestRaoultCustomParameters` pins that.
The *defaults*, though, come from three unconnected places:

- `reactions/equilibrium/interphase.py`: `_P_SAT_REF = 0.03169` atm and `_T_REF_WATER = 298.15`
  K (private module constants), plus the inline literals `dH_vap = 44011.0` J/mol and
  `C_water_mol_L = 55.51` mol/L.
- `chemical_equilibrium/engines/nr/engine.py:55`: `_C_WATER_MOL_L`, derived from
  `_RHO_WATER_G_L = 1000.0` and `_M_WATER_G_MOL = 18.015` (about 55.51).
- `models/vlmodels/adm1/base.py:1050`: a literal `55.51`.

`thermo/liquid/water_properties.py` already has a temperature-dependent
`water_density_kg_per_m3(T_K)` and no vapour-pressure function. Two smaller changes
would help, and neither has been decided: publish the `RaoultEquilibrium` defaults as
named, documented public constants (for example a saturation pressure and an
enthalpy of vaporisation at 25 °C) so a caller can refer to them, and route the
three `C_water` copies through one value. The second is not a pure refactor: the
engine's derived value is 1000 / 18.015 = 55.5093, which differs from the 55.51
literals by about 1.3e-5 relative, so results move slightly wherever they are unified.

**Also (noted 2026-09-28): two water vapour-pressure parameter sets.** Besides
`RaoultEquilibrium`'s defaults above (P_sat = 0.03169 atm, ΔH_vap = 44 011 J/mol,
at 298.15 K), `core/phases.py` has its own for the headspace water vapour pressure
(`water_vapour_P_sat_atm`): the BSM2 values P_sat = 0.0313 bar (= 0.03089 atm) and
ΔH_vap/R = 5290 K (≈ 43 983 J/mol), at 298.15 K, chosen for BSM2 compatibility. The two
differ by about 2.6 % in P_sat at 25 °C. Both now go through the same
`clausius_clapeyron` function, so only the data differs; which set is canonical, and
whether BSM2 compatibility needs its own, is undecided.

## `ReactionBuilder.aerobic_growth` doesn't check the derived yield is achievable

Logged 2026-09-21 during the `chemistry-reactions-kinetics-cleanup` audit
(checkpoint 1), not fixed (Part A/B pure-refactor scope).
`reactions/kinetic/builder.py`'s `aerobic_growth()` derives O2/CO2/H2O stoichiometry
from elemental balance for whatever `yield_gX_gS` it is given, with no check
that the result is physically sensible. Example found during the audit:
acetate at a 0.9 g/g yield returns without error, with CO2 *consumed* and O2
*produced* — the elements still balance (balance validation only checks mass
closure per element, not reaction direction), but the stoichiometry is
non-physical. A fix would check the sign of the derived `nO2`/`nCO2` (or a
thermodynamic yield bound) and raise or warn when the requested yield isn't
achievable by aerobic respiration.

## Deprecated `{name}_HA` fallback is unreachable from repo code

Found 2026-09-23 while deleting `EquilibriumSet.bsm2_default()`, not fixed.
`_compute_species_eq` in `chemical_equilibrium/engines/bisection/acid_base.py`
keeps a deprecated path that emits `{name}_HA` / `{name}_A-` / `{name}_BH+` /
`{name}_B` keys (via `_species_key`) for `EquilibriumDef` entries with no
`species_refs`. `BisectionChemicalEquilibriumEngine.from_reactions()` always
sets `species_refs` (`_add_acid` and `_add_polyprotic_acid` in `engine.py`),
including for BSM2's four VFA rows (`S_ac` gets `['S_ac', 'S_ac-']`, checked by
building the engine from `_build_bsm2_equilibrium_reactions()`). So nothing in
the repo reaches the fallback; only an `EquilibriumSet` built by hand with
`add()` and no `species_refs`, then passed to `solve(equilibrium_set=...)`, does.
An earlier entry in this file claimed BSM2 still used it; that was wrong.

A candidate for deletion, together with the `else` branch of
`_populate_totals_from_phases` (`engine.py`, reads `n_mol` by entry name) and,
if nothing else needs them, the generic-key half of the suffix rule described
under the ionic-strength duplication entry above. Check the tests that build
generic keys directly (`test_speciation.py`) before removing anything.
