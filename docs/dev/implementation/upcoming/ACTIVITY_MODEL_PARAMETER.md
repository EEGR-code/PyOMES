# One `activity_model` parameter — Design Note

> Design discussion, 2026-09-25. Checklist:
> [`ACTIVITY_MODEL_PARAMETER_CHECKLIST.md`](ACTIVITY_MODEL_PARAMETER_CHECKLIST.md).
> No branch or code yet. Written
> from a planning conversation that picked up the `OPEN_WORK.md` entry
> "`chemical_equilibrium`'s `use_activity`/`activity_model` split could be one
> parameter". Direction approved by the repo owner the same day: **replace the
> `use_activity` / `activity_model` pair with one `activity_model` parameter that
> takes a model name or a model object; remove `use_activity` everywhere,
> including `ThermoFramework`; remove the engines' `thermo=` argument; no
> compatibility aliases** (the package has no outside users, so the design is
> chosen as for a new package).

## The problem

The liquid activity model is chosen with two parameters, `use_activity: bool`
and `activity_model: str = "davies"`, which meet in
`thermo/liquid/factory.py`:

```python
def make_activity_model(use_activity: bool, activity_model: str) -> ActivityModel:
    if not use_activity:
        return IdealLiquidModel()
    ...  # "davies" -> DaviesLiquidModel(), "sit" -> SITLiquidModel(), else ValueError
```

- `use_activity=False` means ideal whatever `activity_model` says, so
  `activity_model="sit"` is silently ignored unless the flag is also set.
- `"ideal"` is not an accepted name: `make_activity_model(True, "ideal")` raises.
  Yet the models' own names are `"ideal"`, `"davies"` and `"sit"`, and
  `ThermoFramework.activity_model` returns `"ideal"` for the ideal model.
- Docstrings already assume the single-string design and are wrong today:
  `ReactionSystem.configure_engine` lists `"davies"`, `"sit"` or `"ideal"`; the NR
  engine lists `"davies"`, `"ideal"`.
- The engines build the model inside every `solve()` call, so an invalid pair is
  accepted at construction and fails only at the first solve.
  `test_thermo_framework.py::test_no_thermo_uses_explicit` builds an engine with
  `use_activity=True, activity_model="ideal"` and passes because it never solves.
- A second route exists for model objects: both engines take
  `thermo=ThermoFramework(...)` and use it for nothing but its `liquid_activity`.
  It is not reachable from `ReactionSystem.configure_engine` or
  `StirredTankBuilder.chemistry()`, so a parametrised model (for example
  `SITLiquidModel(epsilon=...)`, whose ε table a name cannot carry) or a
  user-written model can only be used by building an engine by hand and calling
  `attach_engine`.

## Design

- **One parameter:** `activity_model: str | LiquidPhaseModel = "ideal"` wherever
  the pair is taken today. A string is one of `"ideal"`, `"davies"`, `"sit"`
  (case-insensitive); anything else raises `ValueError` naming the accepted
  values. A model object is used as given.
- **One resolver:** `make_activity_model(activity_model)` takes that single
  argument, returns a model object for a name, and returns a model object
  unchanged. It is the only place a name becomes a model.
- **Engines resolve once, at construction.** `BisectionChemicalEquilibriumEngine`
  and `NRChemicalEquilibriumEngine` take `activity_model=` (in `__init__` and
  `from_reactions`), resolve it immediately and store the model object on
  `engine.activity_model`. A bad name fails at construction. Every `solve()`
  uses the stored object. The `thermo=` argument, `engine.use_activity`,
  `engine.thermo` and the private `_liquid_activity` go.
- **`ReactionSystem.configure_engine(activity_model=...)`** resolves the value
  (name or object) immediately, so a bad name fails at that call, and passes the
  resolved model to whichever engine it builds.
- **Stirred tank:** `StirredTankBuilder.chemistry(activity_model="ideal")` and
  `ChemistryConfig.activity_model` take the same type; the factory passes it
  through.
- **`ThermoFramework`** keeps `liquid_activity` as its field and its read-only
  `activity_model` property (the model's name). `use_activity` is removed.
- **`AccuracyMonitor.check_ionic_strength(I_molL, activity_model)`** drops its
  `use_activity` argument and reads the model name from a string or a model
  object.
- **The ideal fast path is chosen by type.** The Bisection and NR solvers skip
  the ionic-strength loop when the model is ideal, and decide that today with
  `getattr(activity_model, "name", "ideal") == "ideal"` (`acid_base.py:535, 616,
  660, 999`, `nr/solver.py:495`). With model objects accepted, that would treat
  a user-written model with no `name` attribute, or one named `"ideal"`, as ideal
  and solve it at I = 0 with no iteration, silently. The five checks become
  `isinstance(activity_model, IdealLiquidModel)`. For the stock models the path
  taken is the same, since each one's name matches its class.
- **`name` stays, as a label only**: a short string for display, messages and
  `ThermoFramework.activity_model`, and what `AccuracyMonitor` looks thresholds
  up by (an unknown name gets no warning). No solver path depends on it.

Behaviour is unchanged for every translated call:
`(use_activity=False, anything)` becomes `"ideal"`, and `(True, "davies")` /
`(True, "sit")` become `"davies"` / `"sit"`. Each builds the same model class the
old pair did, so engine output should be bit-identical. That is proved with a
fingerprint, not assumed. The one behaviour change is intended: an invalid name
now fails when the engine is built rather than at its first solve.

## Out of scope

- **A CV's `chemistry_db` activity model does not reach its engine.** Found while
  scoping this phase and logged in `OPEN_WORK.md` ("A CV's `chemistry_db`
  activity model never reaches its speciation engine"). Fixing it changes results
  for non-ideal databases, so it is not part of this pure refactor.
- `HenryEquilibrium(thermo=...)`, which uses a `ThermoFramework`'s
  `liquid_activity` for its own γ correction, is unrelated to the engine
  parameter and is unchanged.
- The other activity-related `OPEN_WORK.md` entries (ionic-strength duplication,
  the unwired `AccuracyMonitor`, SIT's inline density conversion) are unchanged.

## Inventory (2026-09-25, on `main` at `7cbf115`)

**Method.** `use_activity` searched across every tracked file outside
`docs/dev/` (130 occurrences in 34 files, plus 2 untracked `scratch/` notebooks
that are ignored), then every use read. Engine attribute reads and `thermo=`
callers searched separately. Notebooks inspected cell by cell for saved outputs.

| Area | Files | What changes |
|---|---|---|
| Resolver | `thermo/liquid/factory.py`, `thermo/liquid/__init__.py` (docstring) | single-argument `make_activity_model` |
| Framework | `thermo/framework.py` | remove `use_activity` |
| Engines | `engines/bisection/engine.py`, `engines/nr/engine.py` | parameter, eager resolution, remove `thermo=` |
| Solver internals | `engines/bisection/acid_base.py` (default model at 887-890; ideal checks at 535, 616, 660, 999), `engines/nr/solver.py:495` | default model via the new resolver; ideal checks by type |
| Reaction system | `reactions/reaction_system.py` | `_engine_config`, `configure_engine` |
| Stirred tank | `templates/stirred_tank/{builder,configs,factory}.py` | `chemistry()`, `ChemistryConfig`, pass-through |
| Monitor | `monitoring/accuracy.py` | `check_ionic_strength` signature |
| Models | `models/vlmodels/adm1/{base,bsm2}.py` | call sites |
| Tests | 11 files, about 55 uses (`test_thermo_framework`, `test_accuracy_monitor`, `test_speciation`, `test_nr_speciation_engine`, `test_equilibrium_classification`, `test_chemistry_database`, `test_configs`, `test_precipitation_gas_liquid_cp5`, `test_builder`, validation `test_saturation_index`, `test_phreeqc_nr_agreement`) | call sites; the `thermo=` engine tests become object-`activity_model` tests |
| Notebooks | 9 tracked; 8 carry saved outputs | hand-edited `source` only |
| Generators | `docs/tutorials/ArXiv_preprint/_generate_notebooks.py`, `tests/validation/speciation/_generate_notebooks.py` | call sites, so a regenerate does not restore the old API |
| Docs | `README.md`, `docs/architecture.md` | parameter description |

One notebook prints the removed name in a saved output:
`docs/tutorials/reactions/chemistry_database.ipynb` cell 6 prints
`AD_BASIC.thermo.use_activity` and asserts on it. It changes to
`activity_model`, and that notebook is re-run so the saved output matches.

## Checkpoints (proposed)

1. **Resolver, framework and engines.** `make_activity_model`,
   `ThermoFramework`, both engines, the ideal checks in `acid_base.py` and
   `nr/solver.py`, `ReactionSystem`, `AccuracyMonitor`, and their tests (including
   a new one: a model object that is not `IdealLiquidModel` takes the
   ionic-strength loop). Fix the wrong docstrings.
2. **Stirred tank and models.** Builder, config, factory, ADM1/BSM2, and their
   tests.
3. **Notebooks and generators.** Hand-edit the 9 notebooks and both generators;
   re-run `chemistry_database.ipynb`; parse every edited cell.
4. **Docs and close-out.** `README.md`, `docs/architecture.md`, the old-name
   sweep, and deleting the `OPEN_WORK.md` entry.

**Verification.** Before checkpoint 1, a throwaway script fingerprints engine
output (SHA-256 of every result value) over a fixed set of cases × ideal /
Davies / SIT × Bisection / NR, built through the old parameters. After each
checkpoint the same cases, built through the new parameter, must give the same
hash. Full suite green at every checkpoint, including the BSM2 sentinels. Final
sweep: no `use_activity` left outside `docs/dev/implementation/shipped/` and
`docs/dev/ideas/`.

## Decisions (settled 2026-09-25)

1. **The resolver keeps its name, `make_activity_model`.** It still makes a
   model for a name, and a rename would add churn for no gain.
2. **`engine.activity_model` holds the resolved model object** (its label is
   `engine.activity_model.name`), so the attribute matches what the parameter
   accepts.
3. **The ideal fast path is chosen by type** (`isinstance(..., IdealLiquidModel)`),
   not by `name`; `name` is kept as a label.
