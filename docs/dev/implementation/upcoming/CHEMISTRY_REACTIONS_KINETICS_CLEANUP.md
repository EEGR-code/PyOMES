# chemistry / kinetics / reactions cleanup — design discussion

> **Status: decisions D1-D8 settled 2026-09-21.** Branch
> `chemistry-reactions-kinetics-cleanup` exists with this doc and its checklist;
> no code yet. Written from a review of `PyOMES/chemistry/` (3,893 lines),
> `PyOMES/kinetics/` (308) and `PyOMES/reactions/` (2,806), plus the packages they
> touch (`thermo/`, the top-level `equilibria/`, `chemical_equilibrium/`, `core/`,
> `templates/`, `models/`). Checklist:
> [`CHEMISTRY_REACTIONS_KINETICS_CLEANUP_CHECKLIST.md`](CHEMISTRY_REACTIONS_KINETICS_CLEANUP_CHECKLIST.md).

## Commit discipline for this phase

**The repo owner runs every `git add`, `git commit` and `git push`.** At each
checkpoint: edit, run the checkpoint's sanity check, report, then stop for
review and commit.

## Goal

Apply the scrutiny of `chemical-equilibrium-engines-subfolder` to three more
packages: remove dead and duplicated code, remove import workarounds that no
longer work around anything, give shared concerns one home, fix two defects the
review found, and stop the docs describing deleted code. Three parts, in this
order:

- **Part A — hygiene.** Hoist deferred imports that are no longer load-bearing,
  remove unused imports/constants, fix false docstrings. Bit-identical.
- **Part B — deletions, moves and consolidation.** Delete code with no consumer;
  move code to the package that owns its concern; give the rate laws one home.
  Bit-identical (proved by fingerprint, not assumed).
- **Part C — behaviour fixes.** `ChemistryDatabase.extend()` losing the solver;
  the `PHController` id validator. Behaviour changes, no numerics change to
  existing tests.

Parts A and B are pure refactor: anything that looks like a bug is noted, not
fixed. Part C is isolated in its own commits.

## Settled by the owner (2026-09-21)

- **S1.** There are no outside users of the package. No shims, aliases or
  deprecation periods; repo-internal callers are migrated.
- **S2.** `tests/legacy/` is deleted. GitHub history is the archive.
- **S3.** The three Monod implementations become one source.
- Git history before 2026-08-06 is not available and is not needed.

## Decisions (settled 2026-09-21)

| # | Decision | Outcome |
|---|---|---|
| D1 | The recipe layer (`recipe.py`, `chem_recipe.py`, `types.py`, salt/ion maps in `registry.py`) | **Delete.** It is a lookup table plus arithmetic (name → molar mass, ions, pool contributions) with no working consumer, and its output (pooled `CT_*` totals) is not the shape the framework consumes (species amounts in `n_mol`). A species-based replacement is logged in `OPEN_WORK.md`, not built here. Supersedes the open question in `STRONG_ION_INFERENCE_GENERALIZATION.md`. |
| D2 | `PHController` id validation | **Drop** the constructor-time check (it validates against the recipe table, so it warns on real species such as `NH3` and passes non-species such as `HCl`); delete `registry.py`. **Interim:** a documentation warning on `PHController`'s `chemical_id` / `base_chemical_id` says the ids are not validated and must be a strong-corrector alias or a species in a declared equilibrium, or the dose accumulates as inert. The proper check is the existing design note `PHCONTROLLER_CORRECTOR_VALIDATION.md`, which gets one line saying the old check was removed. |
| D3 | `EquilibriumSet`'s unit-conversion arguments (`Ka`, `lnKa`, `dH`+`dH_unit`, `T_ref`+`T_unit`, `Kw`) | **Delete**, and with them all of `thermo_params.py`. `pKas`, `dH_J_per_mol`, `T_ref_K`, `total_key`, `species_refs` stay. Unit flexibility, if wanted, belongs on `EquilibriumReaction`, the public declaration. |
| D4 | The top-level `equilibria/` package | **Move, do not delete.** `GasEOS`, `IdealGasEOS`, `PengRobinsonEOS`, `CriticalProperties`, `BIOGAS_SPECIES`, `BIOGAS_KIJ` go to `thermo/gas_eos.py`; `HenryIdealVLE` is dropped (no consumer; its partition formula ignores its own `eos`). `ThermoFramework.gas_eos` stays as a real typed field. Add tests pinning current behaviour (about +4). Fix stale docstrings. |
| D5 | Where the single Monod lives; where `ReactionBuilder` lives | **`reactions/rate_laws.py`**, holding `GrowthKinetics`, the shared μ → mol/h conversion, `Monod` and `DualSubstrateMonod`. `ReactionBuilder` stays in `reactions/` (public, reactor-agnostic, a stoichiometry calculator). The builder's `Ko2_gL` term becomes `DualSubstrateMonod(secondary_in_mol_L=False, secondary_MW=32.0)`. One edge case changes: zero `Ko2_gL` with zero O2 returns 0.0 instead of raising `ZeroDivisionError`. |
| D6 | The six other rate laws (Contois, Andrews, ContoisAndrews, Tessier, Moser, Blackman) | **Keep all six**, move them into `reactions/rate_laws.py` unchanged, and add behaviour-level tests (about +8) that pin `Andrews` and `ContoisAndrews` as the reference for a later composable-inhibition redesign. |
| D7 | `chemistry/database.py` + `databases/` | **Move** to a new top-level package `PyOMES/databases/` (name provisional) above `reactions/`. This does not remove the package-level `chemistry <-> reactions` cycle, because `partition.py` still needs `reactions/` (see the layering note below). |
| D8 | `MultispeciesPartitionModel` / `MultispeciesVLEPartition` | **Keep; fix defects only.** `PHENOMENA_PROTOCOL.md` uses them as its deliberate "excluded from the taxonomy" example. Fix: remove the unused `IdealGasEOS` import and its false comment; make the docstring say it is ideal-only (no EOS parameter); replace the duplicate `_kH_mol_L_atm_from_ref` with one shared helper. The 11 tests stay. |

### Notes on decisions

**D5 / D6 — passing multi-species information, and composable inhibition
(deferred design note).** Extra limiting substrates are handled by special cases
today. `DualSubstrateMonod` takes exactly one secondary species through four
scalar arguments (`secondary_id`, `Ko`, `secondary_in_mol_L`, `secondary_MW`), and
`ReactionBuilder.monod_aerobic_growth` exposes only an O2 term (`Ko2_gL`) with the
id `"O2"` fixed. Neither extends to a second or third limiting species (for
example NH3 or a phosphate source) without a new class or a new argument.

It would be good to revise this so that parameters tied to several known species
can be passed as one standard structure: for example a mapping
`{species_id: <parameters for that species>}` (half-saturation constant, units or
molar mass, and whether the species limits or inhibits growth), applied by the
rate law to any number of species. Inhibition should be a feature attached to a
growth model, not a separate class per growth-law × inhibition-law pair. The
maths constrains this: `Andrews` and `ContoisAndrews` use the Haldane form,
μ = μmax · r / (Ks + r + r²/Ki), where the inhibition term sits inside the
denominator. That is not a plain multiplier on the base law (as a factor it is
(Ks + r) / (Ks + r + r²/Ki), which needs the base law's own `Ks`), so a composable
design must treat Haldane specially, for example as an optional `Ki` on the base
law, while multiplicative forms (non-competitive 1/(1 + I/Ki), exponential
exp(−I/Ki)) compose freely with any base law. Inhibition is also often by a
*different* species (product inhibition), which a species-keyed mapping covers.
Open questions: the shape of the per-species parameters; whether species come
from `Species` objects so molar masses are looked up instead of re-entered; which
of the six laws have a meaningful inhibition companion (Tessier, Moser and
Blackman have no standard one).

Not part of this phase. Checkpoint 12 moves the laws as they are, with tests that
pin their current values so the later redesign can be proved to reproduce them
exactly.

**D7 — layering (deferred phase).** Moving `databases/` above `reactions/` fixes
one cause of the package-level `chemistry <-> reactions` cycle. A second cause
remains after Part A hoists the deferred imports in `partition.py`: `Henry`,
`Raoult` and `KspEquilibrium` are equilibrium constraints, a `reactions/` concept
(they need `StoichiometryEntry` and `vant_hoff_log_K`), yet live in `chemistry/`.
The module-level graph is acyclic and the import guard test (checkpoint 2)
enforces that, so nothing breaks. A later layering phase would make the package
graph acyclic too, by either moving `StoichiometryEntry` and its helpers into
`chemistry/`, or moving the three `*Equilibrium` constraint classes into
`reactions/` and leaving the `PartitionModel` protocol (which `core/` imports)
in `chemistry/`. It touches many imports (tests, notebooks) and is its own phase.
Logged in `OPEN_WORK.md` at checkpoint 17.

**D4 — gas EOS consistency (deferred).** The ideal-gas law is hard-coded in at
least seven places (`core/phases.py`, `core/boundaries.py`,
`chemical_equilibrium/engines/nr/solver.py`, `chemistry/partition.py`,
`control/cv_loops.py`, `templates/stirred_tank/factory.py`,
`models/vlmodels/headspace.py`); `ThermoFramework.gas_eos` is read by nothing.
Inside the interface itself, `IdealGasEOS.partial_pressures_atm` returns partial
pressures while `PengRobinsonEOS.partial_pressures_atm` returns fugacities
(y·φ·P): same name, different quantity. This phase moves the code and pins
current behaviour, including that difference, without changing it. A later change
should separate the two methods (`partial_pressures_atm` vs `fugacities_atm`) and
route gas pressure and fugacity through `ThermoFramework.gas_eos`; it overlaps the
constants sweep already in `OPEN_WORK.md`.

## Target layout

```
PyOMES/chemistry/        species.py, common_species.py, species_check.py,
                         compounds.py, partition.py, equilibria.py
PyOMES/reactions/        unchanged, plus rate_laws.py (NEW)
PyOMES/databases/        NEW (name provisional): database.py (ChemistryDatabase),
                         aqueous.py, bioprocess_basic.py, anaerobic_digestion.py
PyOMES/thermo/           unchanged, plus gas_eos.py (NEW)
PyOMES/templates/        stirred_tank/ without kinetics.py
Deleted                  PyOMES/kinetics/, PyOMES/equilibria/, chemistry/thermo_params.py,
                         recipe.py, chem_recipe.py, types.py, registry.py, tests/legacy/
```

## Audit

Method: word-boundary search over every tracked `.py`, `.md`, `.ipynb`, config
file, both notebook generators' string templates, `tests/`, `models/`; an AST
import graph over all 196 `.py` files (2,652 imports, lazy and `TYPE_CHECKING`
tagged); `git log -S`. `scratch/` is gitignored and has no hits. "Unused" means
unused in this repo.

### chemistry/

- **Recipe layer (W1).** `chem_recipe.py` (274 lines) is imported only by three
  `tests/legacy/` files via `from scripts.chem_recipe import ...`; no `scripts/`
  exists and `tests/legacy` is not in `testpaths`. `CHEM_DB` has 47 literal keys
  but 39 unique; eight are duplicated and three differ (`Na2MoO4·2H2O` is
  `CT_MoO4` at line 94, `CT_Mo7O24` at line 132, the later wins). Four
  compound tables (`COMPOUND_DB`, `CHEM_DB`, `Chemical` defaults,
  `SALT_DISSOCIATION_MAP`); 16 of 25 shared compounds disagree on molar mass by
  up to 0.006 g/mol. `SolutionRecipe`, `AqueousTotals*`, `SALT_DISSOCIATION_MAP`,
  `ION_TO_ENGINE_KEY`, `normalize_ion_label`, `ion_to_engine_key`,
  `map_user_ions_to_engine`, `validate_compound_ids`: no consumer.
- **`thermo_params.py` (W2), 765 lines.** Helpers (36-103) are used only by
  `EquilibriumSet`; `AcidDefinition`, `WaterDefinition`, `ThermoSnapshot`,
  `collect_thermo_params`, `validate_thermodynamics` have no consumer;
  `ThermodynamicConfig` is reached only from `models/vlmodels/adm1/bsm2.py:875`
  (`apply_to_cv`, which sets `cv._thermo_config`, which nothing reads) and one
  test fixture. `equilibria.py` imports its helpers at top level and
  `thermo_params.py` defers `EquilibriumSet` (lines 205, 290): the only real
  import cycle in these packages (hoisting it fails with `ImportError`).
- **`partition.py` (W3).** Deprecated `HenryPartition`/`RaoultPartition` aliases
  ("kept for one phase") still have about 50 call sites in 7 test files. `_M_WATER`
  is unused. Line 654 imports `IdealGasEOS` and never uses it; its "avoids
  circular" comment is false. `MultispeciesVLEPartition`'s docstring promises a
  non-ideal-EOS `NotImplementedError`, but there is no EOS parameter.
  `_kH_mol_L_atm_from_ref` duplicates `HenryEquilibrium._kH_mol_L_atm`.
- **`equilibria.py` (W4).** Presets `bsm2_diprotic_co2`, `bsm2_with_sulfide`,
  `adm1_full` have no caller. The engine passes only `pKas`, `dH_J_per_mol`,
  `T_ref_K`, `total_key`, `species_refs`. Line 510-511 says the `_HA`/`_A-`
  fallback "is removed in the PARTITION_MODEL phase"; it is not (already in
  OPEN_WORK).
- **`database.py` (W5).** `extend()` builds `ReactionSystem(existing + new)` and
  drops `solver`, `label` and the engine config. Reproduced: a
  `solver="newton_raphson"` system extended comes back `charge_balance`.
  Latent (no current caller), untested. The docstring names the deleted
  `ReactionSet`.
- **`registry.validate_compound_id` (W6).** `PHController` doses by writing
  `n_mol[chemical_id]`, so its ids are species ids, but the validator checks
  `COMPOUND_DB`. `NH3`, `NH4+`, `H+` warn; `HCl` (not a species) passes.
- **`__init__.py`.** 24 of 37 exported names are never imported via the root,
  including single-letter `g`, `L`, `kg`, `mg`, `mL`.

### kinetics/

- **W7.** Only `tests/standalone/test_kinetics.py` (11 tests) imports it.
  `YeastAcetateV1` is a scaffold; `YeastAcetateContoisO2V1` and
  `DEFAULT_CANONICAL_IDS` have no reference anywhere. `PyOMES/README.md` says it
  is "used by reaction models"; it is not.

### reactions/

- **W8.** `reaction_system.py` `_conservation_monitor` is write-only (the CV keeps
  its own); `_accuracy_monitor` reaches engines that never read it (see the
  OPEN_WORK entry on high-ionic-strength warnings). Docstring lines 42-47
  describe the deleted `SpeciationPropertySolver`. Not split; not changed here.
- **W9.** `plots.py`: no pytest coverage, `plot_vant_hoff` has no caller,
  matplotlib is not declared in any extra although three notebooks call
  `plot_speciation`.
- **W10.** `builder.py`: unused `_NH3`, `Dict` imports; hard-coded `_MW_O2 = 32.0`.
  `aerobic_growth` does not check that the yield is achievable: acetate at
  0.9 g/g returns without error, with CO2 consumed and O2 produced (elements
  still balance). Logged, not fixed here.
- **Redundant deferred imports** in `kinetic.py:102`, `equilibrium.py:254`
  (`_shared` is already imported at top level) and `:321` (`plots`).

### Cross-package

- **X1. Import structure.** Top-level module graph is acyclic; the package-level
  `chemistry <-> reactions` cycle exists only because `databases/` lives inside
  `chemistry/` (see the D7 note for what remains). Scratch experiment: hoisting
  every deferred import listed in Part A gave 2072 passed. Not tested:
  `chemical_equilibrium -> reactions` (3 lazy imports) and
  `reaction_system -> engines` (2); left alone.
- **X2. van 't Hoff and constants.** At least nine mass-action copies in
  `chemistry/`, `reactions/`, `thermo/` (OPEN_WORK lists three) plus four in
  `models/`. After this phase seven remain in the package. Extends the two
  OPEN_WORK entries; not fixed here (numerics-changing by 1 ulp).
- **X3. Monod.** Three implementations: `reactions/builder.py:298`,
  `templates/stirred_tank/factory.py:288`, `Monod` in
  `templates/stirred_tank/kinetics.py:167`. Measured on the real code paths
  (factory closure extracted by AST): identical bit patterns on 200,000
  log-uniform physical points and on a 144-point edge grid with non-negative
  parameters. They differ only for negative `mu_max` or negative `Ks`; the
  template version returns 0.0 there. The builder's O2 path
  (`Ko2_gL` set) against `DualSubstrateMonod(secondary_in_mol_L=False,
  secondary_MW=32.0)` is also bit-identical on 200,000 physical points; it
  differs only for `Ko2_gL = 0` with zero O2 (builder raises
  `ZeroDivisionError`, `DualSubstrateMonod` returns 0.0).
- **X4. `equilibria/`.** No runtime consumer, no test; `GasEOS` survives only as a
  `TYPE_CHECKING` annotation on the always-`None` `ThermoFramework.gas_eos`.
  The Peng-Robinson code gives plausible compressibility factors (CO2 at 20 atm
  0.897, N2 at 50 atm 0.987, CH4 at 50 atm 0.910; a plausibility check against
  approximate reference values, not a validation).
- **X5. Docs on deleted or historical code.** `ReactionSet` (4 places),
  `fermenter.*` (5), `scripts/chem_recipe.py` (2), `v8_split` (4),
  `SpeciationPropertySolver`, `fermenter_unit`, `chem_env` (4);
  `common_species.py:23` shows an import that raises `ImportError`.

## Verification

- Full suite at every checkpoint; baseline 2072 passed.
- **Fingerprints, frozen as committed tests** (per the OPEN_WORK note about engine
  fingerprints): rate-function values as hex bit patterns on a fixed grid for
  `Monod` and for the O2 path (`DualSubstrateMonod`), captured before the
  refactor; `Henry/Raoult/KspEquilibrium` outputs on a fixed grid, before and
  after.
- BSM2 sentinels (`test_bsm2_reference.py`) unchanged.
- Import guard: an AST test asserting the top-level module graph is acyclic.
- Before every delete or move, re-run the search over `.py`, `.ipynb`, `.md`,
  generators.
- `test_gas_constant_single_source.py` still passes after each deletion and move.

## Expected test count

2072 − 11 (kinetics) − 6 (alias tests) − 2 (validator, D2) + 1 (import guard)
+ 2 (rate-law fingerprints) + 8 (rate-law behaviour, D6) + 4 (gas EOS, D4)
+ 2 (`extend()`) = **about 2070**. The counts for new tests are estimates; the
subtractions are exact.

## Deferred (and where it is logged)

All logged in `OPEN_WORK.md` at checkpoint 17.

- van 't Hoff single implementation, constants sweep: existing entries, extended
  with this review's counts.
- `ReactionSystem` monitor plumbing: existing "no engine emits a high-ionic-strength
  warning" entry.
- A species-based recipe helper (weighed salt → species amounts, molar masses from
  `Species` / `ChemicalRegistry`), replacing the deleted recipe layer (D1).
- Gas EOS: separate partial pressure from fugacity; route gas pressure through
  `ThermoFramework.gas_eos` (D4 note).
- Mapping-based parameters for several limiting or inhibiting species, and
  composable inhibition (D5/D6 note).
- Package-level layering phase (D7 note).
- `ReactionBuilder.aerobic_growth`: check that the yield is achievable (W10).
- Molar-mass unification across `Chemical` and `Species`; relocating
  `EquilibriumSet` next to the Bisection engine; `plot_vant_hoff`.
