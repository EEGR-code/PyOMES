# chemistry / kinetics / reactions cleanup — design discussion

> **Status: design discussion, 2026-09-21.** No branch, no code yet. Written from
> a review of `PyOMES/chemistry/` (3,893 lines), `PyOMES/kinetics/` (308) and
> `PyOMES/reactions/` (2,806), plus the packages they touch (`thermo/`, the
> top-level `equilibria/`, `chemical_equilibrium/`, `core/`, `templates/`,
> `models/`). When picked up: follow `README.md`'s "How to start one", copy
> `PHASE_KICKOFF_TEMPLATE.md` to `CHEMISTRY_REACTIONS_KINETICS_CLEANUP_CHECKLIST.md`,
> branch, implement, ship.

## Commit discipline for this phase

**The repo owner runs every `git add`, `git commit` and `git push`.** At each
checkpoint: edit, run the checkpoint's sanity check, report, then stop for
review and commit.

## Goal

Apply the scrutiny of `chemical-equilibrium-engines-subfolder` to three more
packages: remove dead and duplicated code, remove import workarounds that no
longer work around anything, fix two defects the review found, and stop the
docs describing deleted code. Three parts, in this order:

- **Part A — hygiene.** Hoist deferred imports that are no longer load-bearing,
  remove unused imports/constants, fix false docstrings. Bit-identical.
- **Part B — deletions and consolidation.** Delete code with no consumer;
  give Monod one source. Bit-identical (proved by fingerprint, not assumed).
- **Part C — behaviour fixes.** `ChemistryDatabase.extend()` losing the solver;
  the `PHController` id validator. Behaviour changes, no numerics change to
  existing tests.

Parts A and B are pure refactor: anything that looks like a bug is noted, not
fixed. Part C is isolated in its own commits.

## Settled by the owner (2026-09-21)

- **S1.** There are no outside users of the package. No shims, aliases or
  deprecation periods; repo-internal callers are migrated.
- **S2.** `tests/legacy/` is deleted. GitHub history is the archive.
- **S3.** The three Monod implementations become one source. Bit-identity of the
  current three was measured (see X3 below).
- Git history before 2026-08-06 is not available and is not needed.

## Decisions (proposed defaults — confirm before starting)

| # | Decision | Proposed default | Why |
|---|---|---|---|
| D1 | The recipe layer (`recipe.py`, `chem_recipe.py`, `types.py`, salt/ion maps in `registry.py`) | Delete | No working consumer; parallel tables disagree (W1). Supersedes the open question in `STRONG_ION_INFERENCE_GENERALIZATION.md`. |
| D2 | `PHController` id validation | Drop constructor-time check; log a follow-up in `OPEN_WORK.md` to validate dosing ids at attach time | It checks the wrong table (W6). |
| D3 | `EquilibriumSet`'s unit-conversion kwargs (`Ka`, `lnKa`, `dH`+`dH_unit`, `T_ref`+`T_unit`, `Kw`) | Delete; keep `pKas`, `dH_J_per_mol`, `T_ref_K`, `total_key`, `species_refs` | No caller, no test; deleting them leaves nothing in `thermo_params.py` worth keeping. Alternative: fold the ~70 lines of helpers into `equilibria.py`. |
| D4 | The top-level `equilibria/` package (555 lines) and `ThermoFramework.gas_eos` | Delete both (removes 1 test) | No runtime consumer and no test. Restorable from git. |
| D5 | Where the single Monod lives | New `reactions/rate_laws.py` | `reactions/builder.py` must call it; importing from `templates/` would point the dependency upward. |
| D6 | The seven other rate laws in `templates/stirred_tank/kinetics.py` (Contois, Andrews, ContoisAndrews, Tessier, Moser, Blackman, DualSubstrateMonod) | Delete unless wanted as a library | No test, no notebook, no caller. |
| D7 | Move `chemistry/database.py` + `databases/` above `reactions` | Not in this phase | Module-level top imports are already acyclic; the package-level cycle is only geographic. |
| D8 | `MultispeciesPartitionModel` / `MultispeciesVLEPartition` | Keep; fix defects only | `PHENOMENA_PROTOCOL.md` uses them as its deliberate "excluded" example. Revisit when that resolves. |

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
- **Redundant deferred imports** in `kinetic.py:102`, `equilibrium.py:254`
  (`_shared` is already imported at top level) and `:321` (`plots`).

### Cross-package

- **X1. Import structure.** Top-level module graph is acyclic; the package-level
  `chemistry <-> reactions` cycle exists only because `databases/` lives inside
  `chemistry/`. Scratch experiment: hoisting every deferred import listed in
  Part A gave 2072 passed. Not tested: `chemical_equilibrium -> reactions` (3
  lazy imports) and `reaction_system -> engines` (2); left alone.
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
  template version returns 0.0 there. **Not measured:** the builder's
  optional O2 term (`Ko2_gL` set); the fingerprint must cover it before the
  refactor.
- **X4. `equilibria/`.** No runtime consumer, no test; `GasEOS` survives only as a
  `TYPE_CHECKING` annotation on the always-`None` `ThermoFramework.gas_eos`.
- **X5. Docs on deleted or historical code.** `ReactionSet` (4 places),
  `fermenter.*` (5), `scripts/chem_recipe.py` (2), `v8_split` (4),
  `SpeciationPropertySolver`, `fermenter_unit`, `chem_env` (4);
  `common_species.py:23` shows an import that raises `ImportError`.

## Verification

- Full suite at every checkpoint; baseline 2072 passed.
- **Fingerprints, frozen as committed tests** (per the OPEN_WORK note about engine
  fingerprints): Monod rate values as hex bit patterns on a fixed grid (with and
  without `Ko2_gL`), captured before the refactor; `Henry/Raoult/KspEquilibrium`
  outputs on a fixed grid, before and after.
- BSM2 sentinels (`test_bsm2_reference.py`) unchanged.
- Import guard: an AST test asserting the top-level module graph is acyclic.
- Before every delete, re-run the search over `.py`, `.ipynb`, `.md`, generators.
- `test_gas_constant_single_source.py` still passes after each deletion.

## Expected test count

2072 − 11 (kinetics) − 6 (alias tests) − 1 (`gas_eos` test) − 2 (validator, D2)
+ 1 (import guard) + 1 (Monod fingerprint) + 2 (`extend()`) = **2056**.
If D2 is retargeted instead of dropped: 2058. If D4 keeps `gas_eos`: +1.

## Deferred (and where it is logged)

- van 't Hoff single implementation, constants sweep: existing OPEN_WORK entries,
  extended with this review's counts at the end of the phase.
- `ReactionSystem` monitor plumbing: existing OPEN_WORK "no engine emits a
  high-ionic-strength warning" entry.
- Molar-mass unification across `Chemical` and `Species`; relocating
  `EquilibriumSet` next to the Bisection engine; `plot_vant_hoff`; moving
  `databases/`: new OPEN_WORK entries.
