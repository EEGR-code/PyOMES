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
  controllers/, actuators/, builders/}`, and everything in
  `PyOMES/equilibria/` except `vle.py` and `peng_robinson.py` (no
  `engine.py`, `coupled.py`, `factory.py`, `interfaces.py`). Also
  `models/vlmodels/fermenter/` (only `adm1/`, `hplc/`, `headspace.py` remain).
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
  a `thermo/` line. The rest of the tree (e.g. the `equilibria/`, `sim/` and
  `solvers/` "CUFermenter island" entries, and any other package that is
  missing) still needs the re-derivation described below.

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

- `docs/tutorials/D2C_workshop/raw_construction.py` — lines 8, 99, 380
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

## Three copies of the mol/L → mol/kg-water conversion in `PyOMES/thermo/`

Surfaced 2026-09-20 during `chemical-equilibrium-engines-subfolder`
checkpoint 4, while checking where `debye_huckel_A` lives. Left alone there
because that phase is a pure move/delete refactor. The conversion "divide
ionic strength by water density in kg/L" exists three times:

- `water_properties.ionic_strength_molal_from_molar(I_molL, *, T_K)` — the
  public one, exported from `PyOMES.thermo`. Returns `0.0` for a non-finite or
  non-positive `I_molL`, and falls back to returning `I` unchanged when the
  water density is non-finite or non-positive.
- `liquid_phase_model._kg_per_L(T_K)` and `sit_liquid_model._kg_per_L(T_K)` —
  two private copies with identical bodies. Both return the density in kg/L,
  falling back to `1.0` when the density is bad. They are used only to build
  `dIm_dImolL = 1.0 / _kg_per_L(T_K)` in each model's Jacobian method.

The fallbacks agree (a bad density acts as 1 kg/L in all three), so this looks
like harmless duplication rather than a numerical inconsistency, but that is
unverified. A small cleanup would make one public helper in
`water_properties.py` (for example `water_kg_per_L(T_K)`), and have all three
call sites use it. Check the Jacobian tests in
`tests/standalone/test_liquid_phase_model.py` still pass afterwards.

## A third van 't Hoff copy in `reactions/equilibrium.py` differs from `thermo` in the last bit

Surfaced 2026-09-20 during `chemical-equilibrium-engines-subfolder`
checkpoint 7, which merged `acid_base._vant_hoff_K` and
`nr_tableau._vant_hoff_log_K` into `PyOMES/thermo/equilibrium_constants.py`.
`reactions/equilibrium.py:92` `vant_hoff_log_K(constraint, T_K)` has the same
maths and the same edge-case rules as `thermo.equilibrium_constants.vant_hoff_log_K`,
but converts ln K to log10 K with `_LOG10_E = 1.0 / math.log(10.0)`
(`0.43429448190325176`), while the `nr_tableau` version — now the `thermo` one —
uses `np.log10(np.e)` (`0.4342944819032518`). They differ by 1 ulp, so folding
one into the other changes about a fifth of corrected values by up to ~4e-15 in
log10 K (measured on 20,000 random cases). That is negligible physically but is
a numerical change, so it was left out of a pure-refactor phase.

To finish the de-duplication: make `reactions.equilibrium.vant_hoff_log_K` a
thin wrapper that calls `thermo.equilibrium_constants.vant_hoff_log_K(
constraint.log_K, constraint.dH_J_per_mol, T_K, constraint.T_ref_K)`, accept
the 1-ulp shift, and re-run the suite (including `test_equilibrium_constraint.py`
and the BSM2 sentinels) to confirm nothing depends on the last bit. Also rename
one of the two same-named functions if the wrapper is not kept, to avoid two
`vant_hoff_log_K` with different signatures.

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
`equilibria/vle.py`, `chemistry/partition.py`, `thermo/framework.py` and
`thermo/equilibrium_constants.py`.

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
answered for real. **Depends on Part D landing first**, since overriding a
constant that still has ten copies would only change some of them.

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

- `chemistry/equilibria.py:510-511` (in `EquilibriumSet.bsm2_default()`) says
  the deprecated `_HA`/`_A-` fallback "is removed in the PARTITION_MODEL phase".
  That phase shipped and the fallback is still live: the four VFA rows
  (`S_ac`, `S_pro`, `S_bu`, `S_va`) have no `species_refs`, so the BSM2 model
  still uses it. The same false promise was in `acid_base.py` and was rewritten.
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
deprecated `{name}_HA` / `{name}_A-` keys (the VFA rows of
`EquilibriumSet.bsm2_default()`) and `Cation(inert)` / `Anion(inert)` have no
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
