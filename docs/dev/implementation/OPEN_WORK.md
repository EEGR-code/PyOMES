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
  still shows `PyOMES/speciation/` (~line 391), which was renamed to
  `chemical_equilibrium/` at the close of `LAYER1_GAP_CLOSURE` (shipped
  2026-07-03).

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
`PyOMES.chemical_equilibrium.activity_models.make_activity_model(use_activity,
activity_model)` through `engine.py`, `factory.py`, and
`nr_engine.py` — `activity_model` is only meaningful when
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
