# Thermo Subfolder Structure — Design Note

> Design discussion, 2026-09-24. No branch, no checklist, no code yet. Written
> from a planning conversation that audited `PyOMES/thermo/` (8 flat `.py`
> files, 1,459 lines) and asked whether the liquid-phase and gas-phase models
> should each get a subfolder so the package is easier for a new reader to
> navigate. Direction approved by the repo owner the same day: **`liquid/` and
> `gas/` subfolders, each with a `protocols.py`, an `ideal.py` and one file per
> non-ideal model; `framework.py` and `equilibrium_constants.py` stay flat**.
> Modelled on [`REACTIONS_SUBFOLDER_STRUCTURE.md`](REACTIONS_SUBFOLDER_STRUCTURE.md)
> and the `chemical_equilibrium/engines/` reorganisation
> ([`../shipped/CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER.md`](../shipped/CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER.md)).

## Inventory (2026-09-24, on `main` at `f88f372`)

**Method.** Line counts with `wc -l`. Imports read from every file in
`thermo/`. Consumers counted with a script that walks the whole working tree
(every file type; excluding `.git/`, `__pycache__/`, `.pytest_cache/` and
`PyOMES.egg-info/`), matching `thermo.<module>`, `thermo/<module>` and
`thermo\<module>`, then separating real import statements (`from|import` +
`PyOMES.` or relative dots + `thermo.<module>`) from prose. Bare `<module>.py`
mentions were checked by hand in live files (`factory.py` also names
`templates/stirred_tank/factory.py`).

### Files

| File | Lines | Contents | Imports from `thermo/` |
|---|---|---|---|
| `__init__.py` | 45 | public exports | all of the below except `equilibrium_constants` |
| `framework.py` | 117 | `ThermoFramework`, `THERMO_IDEAL`, `THERMO_DAVIES` | `liquid_phase_model`, `gas_eos` |
| `equilibrium_constants.py` | 76 | `vant_hoff_delta_ln_K`, `vant_hoff_K`, `vant_hoff_log_K` | none |
| `liquid_phase_model.py` | 304 | protocols `LiquidPhaseModel`, `ActivityModel`, `DifferentiableLiquidModel`; `IdealLiquidModel`; `DaviesLiquidModel`; private `_kg_per_L`, `_LN10` | `water_properties` |
| `sit_liquid_model.py` | 313 | `SITLiquidModel`, `SIT_EPSILON`, `ION_CHARGES`; private `_kg_per_L`, `_LN10`, `_get_epsilon` | `water_properties` |
| `water_properties.py` | 73 | `water_dielectric_constant`, `water_density_kg_per_m3`, `debye_huckel_A`, `ionic_strength_molal_from_molar` | none |
| `factory.py` | 25 | `make_activity_model` | `liquid_phase_model`, `sit_liquid_model` |
| `gas_eos.py` | 506 | `GasEOS` (lines 40-47), `IdealGasEOS` (50-63), then `CriticalProperties`, `BIOGAS_SPECIES`, `BIOGAS_KIJ`, PR helpers and `PengRobinsonEOS` (66-506) | none |

Every file imports from outside `thermo/` only `PyOMES.units` (`framework`,
`equilibrium_constants`, `gas_eos`) or nothing. All eight files are LF.

### Consumers (files outside `thermo/` with a real import statement)

| Module | Importers | Names pulled |
|---|---|---|
| `liquid_phase_model` | `tests/standalone/test_nr_gas_liquid_cp2.py` (lines 347, 354, 358, 380, 391) | `DaviesLiquidModel` (exported), `DifferentiableLiquidModel` (**not** exported) |
| `sit_liquid_model` | `test_nr_gas_liquid_cp2.py` (353, 371), `test_liquid_phase_model.py:244` | `SITLiquidModel`, `ION_CHARGES` (both exported) |
| `gas_eos` | `test_gas_eos.py` (28, 79, 93) | `PengRobinsonEOS`, `IdealGasEOS`, `BIOGAS_SPECIES` (all exported) |
| `water_properties` | none | — |
| `factory` | none | — |
| `framework` | `databases/aqueous.py:24`, `databases/database.py:31`, `test_thermo_framework.py` (138, 145) | does not move |
| `equilibrium_constants` | `engines/bisection/acid_base.py:43`, `engines/nr/engine.py:47`, `engines/nr/tableau.py:42`, `test_equilibrium_constants.py:15` | does not move |

Everything else imports from the package root (`from PyOMES.thermo import
...`), including the engines (`make_activity_model`, `ActivityModel`), every
notebook and `_generate_notebooks.py`. Those imports are unaffected. No notebook
imports a deep `thermo` path, and the two notebooks with `try/except
ImportError` guard PHREEQC only. Neither guard test names a `thermo` module:
`test_package_layering.py:148` uses `thermo` as synthetic detector input only.

### Found along the way

- **The gas side is not used by the simulation.** Nothing in `PyOMES/` outside
  `thermo/` constructs `IdealGasEOS` or `PengRobinsonEOS` or reads
  `ThermoFramework.gas_eos`; the only users are tests. The ideal-gas law is
  written out inline in `core/` and elsewhere. Already logged in `OPEN_WORK.md`
  ("Gas EOS: `partial_pressures_atm` returns two different quantities, and
  `ThermoFramework.gas_eos` is read by nothing"). This phase gives the classes
  a clearer home; it does not wire them in.
- **`framework.py`'s docstring overstates the gas default.** It says
  `gas_eos=None` means "uses ideal-gas via the KineticGasLiquidLink default";
  no code path consults `gas_eos` at all.
- **`PengRobinsonEOS` does not subclass `GasEOS`**, although its docstring says
  it "Implements the `GasEOS` protocol"; `IdealGasEOS` does subclass it.
  `GasEOS` is a plain base class whose methods raise `NotImplementedError`, not
  a `Protocol`. Nothing calls `isinstance(..., GasEOS)`, so nothing breaks
  today. See Open question 1.
- **Uneven liquid layout.** Ideal and Davies share a file with the three
  protocols; SIT has its own. The `LiquidPhaseModel` docstring lists
  `NRTLLiquidModel` as a future implementation with no obvious home.
- **Duplication already logged:** `_kg_per_L` is defined identically in
  `liquid_phase_model.py:146` and `sit_liquid_model.py:36` (`OPEN_WORK.md`,
  "Three copies of the mol/L → mol/kg-water conversion"); `gamma_all`'s
  ionic-strength rules (`OPEN_WORK.md`, charge-inference entry). Not logged:
  the Davies and SIT `jacobian_dgamma_dx` bodies differ only in the
  `d_log10_dIm` line.
- **Stale prose.** `water_properties.py`'s docstring still names
  `PyOMES/speciation/` (renamed `chemical_equilibrium/`).
  `liquid_phase_model.py` and `sit_liquid_model.py` docstrings cite phase
  labels and design-doc sections (`CP2 of LAYER1_GAP_CLOSURE`, "§8.4 of
  THERMODYNAMIC_MODEL_ARCHITECTURE.md"), which `OPEN_WORK.md` already counts
  (8 such lines in `thermo/`).
- **`docs/dev/ideas/THERMODYNAMIC_MODEL_ARCHITECTURE.md`** still uses `src/`
  paths and places `GasEOS` in `src/equilibria/vle.py`. It is an ideas doc and
  is left alone, as `ideas/` is in other phases.
- **Pickled checkpoints.** `ThermoFramework` holds a `LiquidPhaseModel`
  instance, so any `"data"`-fidelity checkpoint pickles `DaviesLiquidModel` /
  `SITLiquidModel` / `IdealLiquidModel` by module path. The repo holds no
  `.pkl`/`.pickle` files.

## Motivation

A new reader listing `thermo/` sees seven files with no grouping. The code has
two clear areas: liquid-phase non-ideality (γ) and gas-phase non-ideality (φ),
which `LiquidPhaseModel`'s docstring already presents as matching sides of the
fugacity equality:

    φ_i(T,P,y) × y_i × P  =  γ_i(T,x) × x_i × f_i^ref
          ↑ GasEOS                  ↑ LiquidPhaseModel

`framework.py` is where the two sides meet, and `equilibrium_constants.py`
belongs to neither. File names do not show this today: "the Davies model"
lives in `liquid_phase_model.py`, "the ideal gas" in `gas_eos.py`, and
`water_properties.py` and `factory.py` look like general utilities although
only the liquid models use them.

## Groups and seams

| Group | Files | Imports from `thermo/` |
|---|---|---|
| Liquid | `liquid_phase_model`, `sit_liquid_model`, `water_properties`, `factory` | each other only |
| Gas | `gas_eos` | nothing |
| Joins both | `framework` | liquid and gas |
| Neither | `equilibrium_constants` | nothing |

**Neither group imports the other.** Only `framework.py` (and `__init__.py`)
uses both. The split follows an existing seam and creates no import cycle.

## Target layout

```
thermo/
    __init__.py              public exports (names unchanged)
    framework.py             ThermoFramework, THERMO_IDEAL, THERMO_DAVIES
    equilibrium_constants.py van 't Hoff helpers
    liquid/
        __init__.py          docstring only
        protocols.py         LiquidPhaseModel, ActivityModel,
                             DifferentiableLiquidModel     (split from liquid_phase_model.py)
        ideal.py             IdealLiquidModel              (split from liquid_phase_model.py)
        davies.py            DaviesLiquidModel             (split from liquid_phase_model.py)
        sit.py               SITLiquidModel, SIT_EPSILON,
                             ION_CHARGES                   (was sit_liquid_model.py)
        water_properties.py
        factory.py           make_activity_model
    gas/
        __init__.py          docstring only
        protocols.py         GasEOS                        (split from gas_eos.py)
        ideal.py             IdealGasEOS                   (split from gas_eos.py)
        peng_robinson.py     PengRobinsonEOS, CriticalProperties,
                             BIOGAS_SPECIES, BIOGAS_KIJ    (split from gas_eos.py)
```

Top level: 3 files and 2 folders. Each side has the same shape: the interface,
the ideal model, then one file per non-ideal model. A new model (NRTL, SRK)
gets its own file next to the others.

**Splitting `liquid_phase_model.py`.** The module docstring (unit convention,
absent-means-γ=1 rule, dual-protocol note) and lines 41-143 (the three
protocols) go to `protocols.py`. `IdealLiquidModel` (157-182) goes to
`ideal.py` and needs no imports beyond `dataclass` and typing.
`DaviesLiquidModel` (185-304) goes to `davies.py` with `_kg_per_L`, `_LN10`
and its `water_properties` imports. None of the models import the protocols
(they satisfy them structurally), so the three new files do not depend on each
other.

**Splitting `gas_eos.py`.** `GasEOS` (40-47) and the part of the module
docstring about the interface (including the partial-pressure vs fugacity
distinction) go to `protocols.py`. `IdealGasEOS` (50-63) goes to `ideal.py`,
importing `R_L_ATM_PER_MOL_K` and, while it still subclasses `GasEOS`,
`.protocols`. Lines 66-506 go to `peng_robinson.py` with the references
section; `CriticalProperties` and the `BIOGAS_*` tables stay with the only
model that uses them.

## Decisions (2026-09-24)

1. **Layout above:** `liquid/` and `gas/` subfolders; `framework.py` and
   `equilibrium_constants.py` stay at the top level. Chosen over keeping the
   package flat and over a `liquid/` folder with `gas_eos.py` left flat (the
   `blackbox.py` precedent from the reactions note; rejected because the gas
   side splits naturally into interface, ideal and Peng-Robinson, and the
   matching structure is the point).
2. **`gas/ideal.py` is its own file**, matching `liquid/ideal.py`, although it
   is about 20 lines. Being the obvious place to find the ideal gas matters more
   than file size.

## Proposed defaults (not yet confirmed by the owner)

3. **`water_properties.py` and `factory.py` go into `liquid/`.** Only the
   liquid models use them. `OPEN_WORK.md` floats adding a water vapour-pressure
   function to `water_properties.py` for `RaoultEquilibrium`; that is still a
   property of liquid water, so the location holds.
4. **File names drop the suffix the folder now supplies:** `davies.py`,
   `sit.py`, `peng_robinson.py`, not `davies_liquid_model.py` etc. Class names
   are unchanged.
5. **No shims.** Subpackage `__init__.py` files are docstring-only, as in
   `chemical_equilibrium/engines/` and the reactions note. Old deep paths stop
   working.
6. **Package root exports unchanged.** `DifferentiableLiquidModel` stays
   deep-only, now at `PyOMES.thermo.liquid.protocols` (see Open question 3).
7. **Consumers outside `PyOMES/` use the short form.** The three test files
   switch every deep import of an exported name to
   `from PyOMES.thermo import ...`; only `DifferentiableLiquidModel` keeps a
   deep import. Same rule as the reactions note's Decision 8.
8. **Import style.** Moved files import from their own folder with relative
   imports (`from .water_properties import ...`) and from outside it with
   absolute ones (`from PyOMES.units import ...`). Staying files import into
   the folders relatively (`from .liquid.davies import DaviesLiquidModel` in
   `framework.py` and `__init__.py`). Same rule as the reactions note's
   Decision 9.
9. **The moves are pure relocations.** Behaviour-touching changes
   (`GasEOS` as a `Protocol`, `_kg_per_L` de-duplication) are separate
   checkpoints if taken at all (Open questions 1 and 2). Wiring
   `ThermoFramework.gas_eos` into `core/` is out of scope; it stays in
   `OPEN_WORK.md`.
10. **Docstrings of moved code describe current behaviour.** The phase labels
    and design-doc section references in the Davies and SIT docstrings are
    rewritten in the docstring checkpoint, since that text is being moved
    anyway; `water_properties.py`'s stale `speciation/` path and
    `framework.py`'s `KineticGasLiquidLink` claim are corrected at the same
    time. The `OPEN_WORK.md` count is updated.

## What breaks and what has to change

**Deep paths that stop working** (files with a real import statement):

| Old path | New path | Import files | Of which need a non-exported name |
|---|---|---|---|
| `thermo.liquid_phase_model` | `.liquid.protocols` / `.liquid.ideal` / `.liquid.davies` | 1 (`test_nr_gas_liquid_cp2.py`, 5 lines) | 1 (`DifferentiableLiquidModel`, lines 348 and 354) |
| `thermo.sit_liquid_model` | `.liquid.sit` | 2 (`test_nr_gas_liquid_cp2.py`, `test_liquid_phase_model.py`) | 0 |
| `thermo.gas_eos` | `.gas.protocols` / `.gas.ideal` / `.gas.peng_robinson` | 1 (`test_gas_eos.py`, 3 lines) + own doctest (line 19) | 0 |
| `thermo.water_properties` | `.liquid.water_properties` | 0 | 0 |
| `thermo.factory` | `.liquid.factory` | 0 | 0 |

**Unchanged paths:** `thermo.framework`, `thermo.equilibrium_constants`, and
the package root.

**Import rewrites inside `thermo/`:** `__init__.py` (the liquid, SIT, factory,
water and gas import blocks) and `framework.py` (lines 20-21). Moved files:
`liquid/davies.py` and `liquid/sit.py` keep `from .water_properties import`;
`liquid/factory.py` imports `.ideal`, `.davies`, `.sit` and `.protocols`;
`gas/ideal.py` and `gas/peng_robinson.py` switch `..units` to
`PyOMES.units`.

**Prose and docstrings to repoint:** `liquid_phase_model.py:4` (the `GasEOS`
path), `sit_liquid_model.py:227` (`:class:` cross-reference to
`DifferentiableLiquidModel`), `gas_eos.py:19` (doctest), module docstrings of
the new files, `test_gas_eos.py:2`; live docs
`docs/architecture.md:393-396` (the `thermo/` tree entry), `PyOMES/README.md:22`
(class list; add the subfolders), and the live `OPEN_WORK.md` entries that cite
the moved files (the mol/L → mol/kg entry, the phase-label entry, the
charge-inference entry, the gas EOS entry, the Raoult water-values entry).
`OPEN_WORK.md`'s historical sentences about files moving *into*
`thermo/gas_eos.py` (lines 20, 35) record past events and stay as written.
`docs/dev/implementation/shipped/` and `docs/dev/ideas/` are left alone.

**Scale.** About 15 live files: 8 in `thermo/` (plus 8 new ones), 3 tests,
`docs/architecture.md`, `PyOMES/README.md`, `OPEN_WORK.md` and
`upcoming/README.md`. Nothing outside `thermo/` in `PyOMES/` changes.

## Risks

- **Ambiguous names in searches.** `protocols.py`, `ideal.py` and `factory.py`
  also exist elsewhere (`reactions/protocols.py`,
  `chemical_equilibrium/protocols.py`, `templates/stirred_tank/factory.py`),
  and `ideal.py` twice within `thermo/`. Sweeps must match full paths
  (`thermo/liquid/ideal`), not bare file names.
- **Sphinx-style cross-references** (`:class:~PyOMES.thermo...`) fail
  silently; include them in the old-path sweep.
- **Pickled checkpoints** holding a `ThermoFramework` (and so a liquid model)
  will not load after the move. Accepted, as in the partition relocation and the
  reactions note.
- **Low notebook risk.** No notebook imports a deep `thermo` path, so no notebook
  needs re-running for this phase. Re-check this at checkpoint 1.

## Checkpoint shape

Rough; the checklist decides granularity. Full suite before and after each
(baseline 2098 passed, 0 failed on `main` at `f88f372`).

1. Re-verify this note's inventory with a fresh search; record baseline.
2. Create `liquid/`: split `liquid_phase_model.py` into `protocols.py`,
   `ideal.py`, `davies.py`; move `sit_liquid_model.py` to `sit.py`, plus
   `water_properties.py` and `factory.py`; rewrite `__init__.py`,
   `framework.py` and the two liquid test files.
3. Create `gas/`: split `gas_eos.py` into `protocols.py`, `ideal.py`,
   `peng_robinson.py`; rewrite `__init__.py`, `framework.py` and
   `test_gas_eos.py`.
4. Docstrings, prose and live docs (Decision 10), including subpackage
   docstrings that say what each folder holds.
5. Sweep every old path in every form (dotted, slash, relative, bare filename,
   Sphinx cross-reference) across `.py`, `.ipynb`, `.md`, `.toml`, `.yml`,
   `.json`; confirm each old path fails to import.
6. Whichever of Open questions 1, 2 and 4 are accepted, each as its own
   checkpoint.

## Open questions

1. **Make `GasEOS` a `Protocol`?** `@runtime_checkable`, like
   `LiquidPhaseModel`; `IdealGasEOS` stops subclassing it and
   `PengRobinsonEOS` then satisfies it for real, so both sides of `thermo/`
   work the same way. Nothing checks `isinstance(..., GasEOS)` today, so no
   caller changes. Recommended: yes, as its own checkpoint after the move. The
   alternative is having `PengRobinsonEOS` subclass `GasEOS` as it stands.
2. **De-duplicate `_kg_per_L` while here?** The `OPEN_WORK.md` entry proposes
   one public helper in `water_properties.py`; the two copies have identical
   bodies, so results should not change. Both files are being moved anyway.
   Recommended: yes, as its own checkpoint, closing the `OPEN_WORK.md` entry.
   Leave the near-identical Jacobian loops alone (a shared helper would need to
   take the DH derivative as a parameter; not worth it for two models).
3. **Export `DifferentiableLiquidModel` from the package root?** It is a public
   protocol used only by tests today. Recommended: no, keep root exports
   unchanged (Decision 6); revisit if a production consumer appears.
4. **Guard the seam with a test?** A small test in the style of
   `test_package_layering.py` asserting that `liquid/` and `gas/` never import
   each other and that neither imports `framework`. Recommended: yes, as the
   last checkpoint, so a later edit cannot blur the groups.
5. **Sequencing with `REACTIONS_SUBFOLDER_STRUCTURE.md`.** The two phases touch
   different packages. They overlap only in prose: both edit the
   `docs/architecture.md` package tree, `OPEN_WORK.md` and this folder's
   `README.md`. `equilibrium_constants.py`'s docstring names
   `reactions.equilibrium` and is repointed by the reactions phase, not this
   one. Recommended: independent; whichever lands second rebases the shared
   prose.

## How to start one

Per this folder's convention ([README.md](README.md#how-to-start-one)): resolve
the proposed defaults and open questions, write
`THERMO_SUBFOLDER_STRUCTURE_CHECKLIST.md`, and cut a branch off `main` named
`thermo-subfolder-structure`. The change is file moves, two file splits and
import rewrites, with no behaviour change; the risk is a missed call site,
which the full suite and the old-path sweep should catch.
