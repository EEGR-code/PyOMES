# Open work

Standalone follow-up items surfaced during other phases — not yet
scoped as their own phase, no branch, no checklist. Referenced from
[`upcoming/README.md`](upcoming/README.md).

## `docs/architecture.md` still describes deleted CUFermenter-era code

Surfaced 2026-09-15; re-audited against `main` 2026-09-19. The doc still
frames several subpackages as a "CUFermenter island" awaiting the
trigger-gated `CUFERMENTER_SUNSET` phase, but `cufermenter-sunset` already
shipped 2026-06-01 (9,093 lines deleted, see
[`shipped/CUFERMENTER_SUNSET.md`](shipped/CUFERMENTER_SUNSET.md)). The
paths it names are largely gone, so the whole "CUFermenter island" framing
is dead and should come out, not just be re-pathed:

- **No longer exist:** `PyOMES/sim/`, `PyOMES/solvers/`,
  `PyOMES/core/gl_equilibrium.py`, `PyOMES/control/{loops.py, system.py,
  controllers/, actuators/, builders/}`, and `PyOMES/equilibria/` (deleted
  entirely 2026-09-22; its two real files moved into `thermo/gas_eos.py`, see
  below). Also `models/vlmodels/fermenter/` (only `adm1/`, `hplc/`,
  `headspace.py` remain).
- **Stale references to those:** the "Equilibrium pathways" section
  (~lines 134–154, which presents `ProcessCoupledEquilibrator` +
  `HenryEquilibriumInterface` as a live second pathway); mentions of `CUFermentationSpeciation`
  (~lines 152, 253, 353); the Repository Layout tree (~lines 372–418).
- **Also stale, not mentioned in the original entry:** the layout tree
  still showed `PyOMES/speciation/` (~line 391), which was renamed to
  `chemical_equilibrium/` at the close of `LAYER1_GAP_CLOSURE` (shipped
  2026-07-03). **Partly fixed 2026-09-20:** `chemical-equilibrium-engines-subfolder`
  (checkpoint 12) replaced that one block with the real
  `chemical_equilibrium/` tree (including the new `engines/` layout) and added
  a `thermo/` line. **Partly fixed 2026-09-22:** `chemistry-reactions-kinetics-cleanup`
  (checkpoint 11, decision D4) moved `PyOMES/equilibria/`'s two real files
  (`vle.py`, `peng_robinson.py`) into `thermo/gas_eos.py` and deleted the
  package, so the "Equilibrium pathways" section's item 2
  (`ProcessCoupledEquilibrator` + `HenryEquilibriumInterface`, presented as a
  live second pathway) is gone rather than re-pathed, and the tree's
  `equilibria/` line is removed rather than corrected. The rest of the tree
  (the `sim/` and `solvers/` "CUFermenter island" entries, the
  `CUFermentationSpeciation` mentions at ~244/344, and any other package that
  is missing) still needs the re-derivation described below.

Needs a pass that re-derives the Repository Layout and the equilibrium
section from the actual tree rather than a line-by-line patch.

## `tests/run_tests.py` imports a deleted `create_standalone_fermenter`

Surfaced 2026-09-18 during a README.md audit. `tests/run_tests.py` still does
`from PyOMES import PressureReliefController, PHController,
create_standalone_fermenter` and calls it in `_fermenter_minimal()` /
`_fermenter_full()`. That name was never restored after
`cufermenter-sunset` (2026-06-01) — it isn't exported from
`PyOMES/__init__.py` and doesn't exist under `models/vlmodels/`
either, so this script currently fails on import. It isn't part of
the pytest suite (`pyproject.toml`'s `testpaths` only covers
`tests/standalone` and `tests/validation`), so it doesn't show up as
a CI failure — likely why it's gone unnoticed. Needs its own pass:
either migrate it to `StirredTankBuilder` (mirroring the
`stirred-tank-template` migration already done for the tutorial
notebooks) or delete it if it's fully superseded by
`tests/standalone`.

## `chemical_equilibrium`'s `use_activity`/`activity_model` split could be one parameter

Surfaced 2026-09-18 while checking
`StirredTankBuilder.chemistry()`'s `use_activity: bool` +
`activity_model: str = "davies"` signature. The split is threaded from
`PyOMES.thermo.make_activity_model(use_activity, activity_model)` (it lived in
`chemical_equilibrium/activity_models.py` until
`chemical-equilibrium-engines-subfolder` moved it) through
`engines/bisection/engine.py` and `engines/nr/engine.py` (the
`chemical_equilibrium/factory.py` that used to be in this list was deleted in
that phase) — `activity_model` is only meaningful when
`use_activity=True`, and `"ideal"` is not itself a valid value for
`activity_model` (it's only reachable via `use_activity=False`).
Consider collapsing this into a single `activity_model: str`
parameter that accepts `"ideal"` alongside `"davies"`/`"sit"`,
removing the separate boolean gate. Touches the builder, the three
engines above, and their callers — worth scoping as its own small
phase rather than a drive-by fix.

## Leftover `FermenterBuilder` mentions in the D2C workshop and reactions tutorials

Originally surfaced 2026-09-15 as part of the `stirred-tank-template`
migration; re-audited 2026-09-19. The old `FermenterBuilder`/
`FermenterFactory` names survive only as prose comparison points (the code
runs correctly) in:

- `docs/tutorials/D2C_workshop/raw_construction.py` — lines 8, 100, 381
  (module docstring, a comment, and a printed demo banner).
- `docs/tutorials/reactions/reaction_system.ipynb` — one comment ("…
  FermenterBuilder uses by default").

Both should say `StirredTankBuilder`. Sweep together; for the notebook,
hand-edit the one cell rather than regenerating (see next entry).

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

## A third van 't Hoff copy in `reactions/equilibrium/constraint.py` differs from `thermo` in the last bit

Surfaced 2026-09-20 during `chemical-equilibrium-engines-subfolder`
checkpoint 7, which merged `acid_base._vant_hoff_K` and
`nr_tableau._vant_hoff_log_K` into `PyOMES/thermo/equilibrium_constants.py`.
`reactions/equilibrium/constraint.py:63` `vant_hoff_log_K(constraint, T_K)` has the same
maths and the same edge-case rules as `thermo.equilibrium_constants.vant_hoff_log_K`,
but converts ln K to log10 K with `_LOG10_E = 1.0 / math.log(10.0)`
(`0.43429448190325176`), while the `nr_tableau` version — now the `thermo` one —
uses `np.log10(np.e)` (`0.4342944819032518`). They differ by 1 ulp, so folding
one into the other changes about a fifth of corrected values by up to ~4e-15 in
log10 K (measured on 20,000 random cases). That is negligible physically but is
a numerical change, so it was left out of a pure-refactor phase.

To finish the de-duplication: make `reactions.equilibrium.constraint.vant_hoff_log_K` a
thin wrapper that calls `thermo.equilibrium_constants.vant_hoff_log_K(
constraint.log_K, constraint.dH_J_per_mol, T_K, constraint.T_ref_K)`, accept
the 1-ulp shift, and re-run the suite (including `test_equilibrium_constraint.py`
and the BSM2 sentinels) to confirm nothing depends on the last bit. Also rename
one of the two same-named functions if the wrapper is not kept, to avoid two
`vant_hoff_log_K` with different signatures.

**Update 2026-09-22 (`chemistry-reactions-kinetics-cleanup` audit, 2026-09-21):**
that phase's audit counted at least nine mass-action/van 't Hoff correction
copies across `chemistry/`, `reactions/`, `thermo/` (this entry and the
ADM1/BSM2 one below already covered some of them) plus four more in `models/`,
estimating seven left in the package after the phase — an estimate, not a
re-checked invariant. What did happen, checked directly: deleting
`chemistry/thermo_params.py` at checkpoint 7 (2026-09-22, decision D3 — its
`AcidDefinition.pKas_at_T`/`WaterDefinition.Kw_at_T` had no consumer) removed
two copies outright; `chemical_equilibrium/engines/bisection/equilibria.py`'s own
`EquilibriumDef.pKas_at_T`/`WaterDef.Kw_at_T` (same formula, real consumer)
survived unmerged, just repointed to import `R_J_PER_MOL_K` from `units.py`
directly; checkpoint 4's D8 fix separately collapsed a duplicate inside
`chemistry/partition.py` (`HenryEquilibrium._kH_mol_L_atm`, since moved to
`reactions/equilibrium/interphase.py`, delegates to `chemistry/partition.py`'s
module-level `_kH_mol_L_atm_from_ref` instead of repeating it, so that
pair counts as one implementation, not two — though it is a Henry-constant
correction, not a pKa/Kw one, so whether the original count included it isn't
clear from the audit text). A quick recount by function (not by line) after
the phase, restricted to distinct `math.exp(-dH/R * (1/T - 1/T_ref))`-shaped
implementations, finds eight in `chemistry/`+`reactions/`+`thermo/`
(`equilibria.py` ×2, `partition.py` ×2 (now one each in `chemistry/partition.py`
and `reactions/equilibrium/interphase.py`), `reactions/equilibrium/constraint.py` ×1,
`thermo/equilibrium_constants.py` ×1 canonical, `thermo/framework.py` ×2) —
close to, not exactly, the audit's "seven" estimate; not reconciled further
here. Not fixed in that phase either way (Part A/B were pure-refactor,
numerics-changing consolidation was explicitly out of scope) — the
de-duplication described above is still the fix.

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
`thermo/gas_eos.py`, `chemistry/partition.py`, `reactions/equilibrium/interphase.py`,
`thermo/framework.py` and `thermo/equilibrium_constants.py`.

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
`reactions/` 5, `control/` 4, other 4. Some are in this phase's new code and
its neighbours, for example `thermo/liquid_phase_model.py` ("Added by CP2 of
`LAYER1_GAP_CLOSURE`…").

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

## Rename `phreeqc_to_vlsim` (and drop the old project name `vlsim`)

`chemical_equilibrium/engines/phreeqc.py` still names PyOMES's former project
name in a public function, `phreeqc_to_vlsim`, which is also the default
`species_map`, and in a local variable `vlsim_name`. Checkpoint 12b reworded the
docs but kept the names, because renaming is a code change. It is also imported
by `tests/standalone/test_phreeqc_engine.py` and used by the
`03_phreeqc_engine_basics.ipynb` tutorial. With no outside users a rename
without an alias is cheap (for example `phreeqc_to_pyomes`): one module, one
test file, one notebook.

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
- `thermo/`: `gamma_all` in `liquid_phase_model.py` and `sit_liquid_model.py`
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
`thermo/`, but should be checked for a second copy).

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
`thermo/gas_eos.py` (checkpoint 11), `chemistry/database.py` /
`chemistry/databases/` → `PyOMES/databases/` (checkpoint 13) — so the table's
per-file counts are stale, though the phase was pure-refactor for these
literals (moved, not edited), so the aggregate totals per constant should be
close to unchanged. Spot check: `101325`/`298.15` alone still appear 18 times
across just `chemistry/partition.py` (since split into it and
`reactions/equilibrium/interphase.py`), `thermo/gas_eos.py`,
`chemical_equilibrium/engines/bisection/equilibria.py`,
`reactions/kinetic/rate_laws.py` and `PyOMES/databases/*.py` post-move. Re-running the
survey against the new layout is part of picking this sweep up, not done here.

## Two notebooks do not parse, and the notebooks edited for the gas constant were not re-run

Found 2026-09-20 while checking `chemical-equilibrium-engines-subfolder` before
merging (a parse of every notebook code cell). Both problems predate that phase.

- `docs/tutorials/templates/batch_fermenter.ipynb`, cell 10, fails with
  "unterminated string literal" on every Python version tried (3.10 and 3.12).
- `tests/validation/speciation/05_precipitation_equilibrium.ipynb`, cell 6, uses an
  f-string form that only Python 3.12 accepts (`f-string: unmatched '['` on 3.10),
  while `pyproject.toml` declares `requires-python = ">=3.10"`.

No test executes a notebook, so CI does not notice either. Separately, the
gas-constant change edited eight notebooks and one generator by text (the import,
the identifier, and three literal values) and validated that each parses, but did
not re-execute them. Their saved outputs still show numbers from before `R`
moved by 4e-7 relative, which is below what most cells print but not verified. They
refresh the next time each notebook is re-run; see also "Regenerating tutorial
notebooks can wipe baked outputs" above.

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
selector is an API decision, and it interacts with "`use_activity` /
`activity_model` split could be one parameter" above.

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
(resolved via `PyOMES.compounds.ChemicalRegistry`, which already carries
molar masses) and a mass or stock-solution dose, emit species amounts
directly (`{"Na+": n_mol, "Cl-": n_mol, ...}`), not `CT_Na`/`CT_Cl`-style
pooled totals. This also supersedes the recipe-layer open question in
`docs/dev/implementation/upcoming/STRONG_ION_INFERENCE_GENERALIZATION.md`.

## Gas EOS: `partial_pressures_atm` returns two different quantities, and `ThermoFramework.gas_eos` is read by nothing

Logged 2026-09-22, `chemistry-reactions-kinetics-cleanup` checkpoint 11
(decision D4), moved but not fixed. `IdealGasEOS.partial_pressures_atm`
(`thermo/gas_eos.py`) returns partial pressures (`y_i × P`);
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
`ReactionBuilder.monod_aerobic_growth` exposes only an O2 term (`Ko2_gL`)
with the id `"O2"` fixed. Neither extends to a second or third limiting
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

`thermo/water_properties.py` already has a temperature-dependent
`water_density_kg_per_m3(T_K)` and no vapour-pressure function. Two smaller changes
would help, and neither has been decided: publish the `RaoultEquilibrium` defaults as
named, documented public constants (for example a saturation pressure and an
enthalpy of vaporisation at 25 °C) so a caller can refer to them, and route the
three `C_water` copies through one value. The second is not a pure refactor: the
engine's derived value is 1000 / 18.015 = 55.5093, which differs from the 55.51
literals by about 1.3e-5 relative, so results move slightly wherever they are unified.

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

## Two small loose ends from the `chemistry-reactions-kinetics-cleanup` audit

Logged 2026-09-21/22, not built — grouped here because each is small and
independent, not because they are related to each other.

- **Molar-mass unification.** `PyOMES/compounds.py`'s `Chemical` (BioSTEAM-shaped,
  moved out of `chemistry/` during this phase since it has no dependency on
  `Species` or anything else in `chemistry/`)
  and `chemistry/species.py`'s `Species` (the framework's own) both carry
  molecular weights for overlapping compound sets, maintained independently.
  A candidate for unification, or at least a documented invariant that they
  agree, if the divergence ever causes a real bug (audit found up to sixteen
  compounds with disagreeing molar masses between `Chemical` and other tables
  before this phase's deletions removed most of those tables).
- **`plot_vant_hoff` has no caller.** `reactions/equilibrium/plots.py`'s `plot_vant_hoff`
  (and `reactions/equilibrium/reaction.py`'s `EquilibriumReaction.plot_vant_hoff`
  wrapper) has no caller anywhere in the repo — no pytest coverage, no
  notebook use. A candidate for deletion in a future pass; left alone here
  since trimming plotting helpers wasn't this phase's scope. `plots.py` now
  has its own `plots` extra (`pyproject.toml`, checkpoint 16 of the same
  phase) if it is kept.

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

## Relative imports three or more dots deep

Found 2026-09-23 while moving `EquilibriumSet`, not fixed. 11 lines across 6
files use `from ....thermo import ...`-style imports, all in
`chemical_equilibrium/engines/bisection/` and `chemical_equilibrium/engines/nr/`,
where the extra nesting level made them long. The rest of the package
(`control/`, `templates/stirred_tank/`) uses absolute
`from PyOMES.… import …`, as do the engines' lazy imports of `PyOMES.reactions`.
Both work with the editable install. Converting the 11 lines is mechanical; pick one convention if the package ever gets a style
pass.
