# Phase Kickoff Checklist — chemistry-reactions-kinetics-cleanup

> Checklist for [`CHEMISTRY_REACTIONS_KINETICS_CLEANUP.md`](CHEMISTRY_REACTIONS_KINETICS_CLEANUP.md),
> the source of truth for goals, audit and decisions D1-D8 (all settled
> 2026-09-21); do not restate them here. See `README.md`'s "Branching and tagging
> convention". Modelled on
> [`../shipped/CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER_CHECKLIST.md`](../shipped/CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER_CHECKLIST.md).

**Working rules**

- Parts A and B are pure refactor: no behaviour or numerical change. Anything that
  looks like a bug is noted, not fixed. Part C changes behaviour deliberately.
- The repo owner runs every `git add`, `git commit`, `git push`. At each
  checkpoint: edit, run its sanity check, report, stop.
- No shims, aliases or deprecation periods (S1).
- Before deleting or moving anything, re-search `.py`, `.ipynb`, `.md`, and both
  notebook generators' string templates. If the plan is wrong, stop and report.
- Use `git mv` for moves. Leave `docs/dev/implementation/shipped/` untouched.
- Test-count arithmetic in the checkpoints is cumulative from the 2072 baseline;
  counts for new tests are estimates.

## Pre-flight

- [x] `git status -sb` clean (checked before branching, 2026-09-21)
- [x] `git log origin/main..main --oneline` empty
- [x] Decisions D1-D8 confirmed and recorded in the plan doc (2026-09-21)
- [x] Branch created off `main`: `chemistry-reactions-kinetics-cleanup`
- [x] This checklist and the plan doc committed as the first commit (`9fcf31a`;
      the updated versions are committed as the second)

## During

- [x] Plan doc exists in `docs/dev/implementation/upcoming/`
- [ ] Checkpoints tracked below, one commit each
- [ ] **If work stalls:** add a status banner to the top of the plan doc at once

### Checkpoint 1 inventory (record here)

- [ ] Suite baseline (expected 2072 passed); line counts (chemistry 3,893, kinetics
      308, reactions 2,806, equilibria 555); import-graph and consumer scripts rerun
      and results pasted; searches repeated at branch start.

### Checkpoints

**Part A — hygiene (bit-identical)**

- [ ] 1. Baseline and inventory (above).
- [ ] 2. Add the import guard test: top-level module graph of `PyOMES/` acyclic
      (function-level and `TYPE_CHECKING` imports excluded). It passes today.
      Sanity: 2073 tests.
- [ ] 3. Hoist the deferred imports that are not load-bearing: `partition.py`
      (six sites; `_resolve_species` reads `chemistry.common_species` directly),
      `database.py:85`, `equilibria.py:493,549`, `reaction_system.py:171,484`,
      `stoichiometry.py:175`, `kinetic.py:102`, `equilibrium.py:254,321`.
      Leave `thermo_params.py` (deleted in 7) and the `chemical_equilibrium`
      engine imports. Sanity: full suite; import guard green.
- [ ] 4. Remove unused imports and constants (`_M_WATER`, the `IdealGasEOS` import
      and its false comment in `partition.py`, `_NH3`, `Dict`, `warnings`,
      `field`/`Optional` in `equilibria.py`, `Species` in `databases/aqueous.py`,
      `Sequence` in `plots.py`). D8: make the `MultispeciesVLEPartition` docstring
      say it is ideal-only, and replace `_kH_mol_L_atm_from_ref` with one helper
      shared with `HenryEquilibrium` (leave the `101325` literal for the constants
      sweep). Fix false docstrings in files that survive (`ReactionSet`,
      `fermenter.*`, `common_species.py:23`, `equilibria.py:510`,
      `reaction_system.py:42-47`). Sanity: Henry/Raoult/Ksp fingerprint
      byte-identical; the 11 `Multispecies` tests pass.

**Part B — deletions, moves and consolidation (bit-identical)**

- [ ] 5. Delete `PyOMES/kinetics/` and `tests/standalone/test_kinetics.py`; remove
      the `kinetics/` row in `PyOMES/README.md`, the row in the root `README.md`
      test table, and the `kinetics/` line in `docs/architecture.md`. Sanity:
      2062 tests.
- [ ] 6. Delete `tests/legacy/` (15 files, not collected). Sanity: grep finds no
      reference in docs or CI; suite unchanged.
- [ ] 7. (D3) Delete `thermo_params.py`, the unit-conversion arguments of
      `EquilibriumSet.add`/`set_water`, and `_VALID_*`. Drop `thermo=` and the
      no-op `apply_to_cv` call from `build_bsm2_cv` (and its docstring); edit the
      fixture in `test_bsm2_reference.py`; prune `chemistry/__init__.py`. Sanity:
      BSM2 sentinels unchanged; suite unchanged.
- [ ] 8. Delete the three unused `EquilibriumSet` presets. Sanity: suite unchanged.
- [ ] 9. Migrate ~50 `HenryPartition`/`RaoultPartition` call sites (7 test files) to
      `HenryEquilibrium`/`RaoultEquilibrium`; delete the aliases and their 6
      tests; fix the `anaerobic_digestion.py` docstring that names the old type.
      Sanity: 2056 tests; Henry/Raoult/Ksp fingerprint byte-identical.
- [ ] 10. (D1) Delete `recipe.py`, `chem_recipe.py`, `types.py` and the ion/salt maps
       in `registry.py`; update `chemistry/__init__.py`; update the recipe-layer
       section of `STRONG_ION_INFERENCE_GENERALIZATION.md` to "resolved, see the
       cleanup plan". Leave `COMPOUND_DB`/`resolve_compound`/`validate_compound_id`
       for checkpoint 15.
- [ ] 11. (D4) Move `equilibria/` to `thermo/gas_eos.py` with `git mv`: `GasEOS`,
       `IdealGasEOS`, `PengRobinsonEOS`, `CriticalProperties`, `BIOGAS_SPECIES`,
       `BIOGAS_KIJ`. Drop `HenryIdealVLE` and the `equilibria/` package. Make
       `ThermoFramework.gas_eos` a typed `Optional[GasEOS]` (real import, no
       `TYPE_CHECKING`); export from `thermo/__init__.py`. Add tests pinning current
       behaviour: known Z values (CO2 20 atm ≈ 0.897, N2 50 atm ≈ 0.987, CH4 50 atm
       ≈ 0.910, 1 atm mixture ≈ ideal) and that `PengRobinsonEOS.partial_pressures_atm`
       returns fugacities while `IdealGasEOS`'s returns partial pressures. Fix the
       stale `fermenter.equilibria.vle` docstring, `thermo/liquid_phase_model.py`'s
       reference, the `PyOMES/README.md` row and the `docs/architecture.md` lines.
       Sanity: `test_gas_eos_none_by_default` still passes; about 2060 tests.
- [ ] 12. (D5, D6) **Before editing:** freeze the rate-function fingerprints as
       tests, for the `Monod` path and for the builder's `Ko2_gL` path (measured
       bit-identical to `DualSubstrateMonod(secondary_in_mol_L=False,
       secondary_MW=32.0)`), captured against the current code. Then create
       `reactions/rate_laws.py` with `GrowthKinetics`, the shared μ → mol/h
       conversion, `Monod`, `DualSubstrateMonod` and the six other laws (unchanged);
       make the factory fallback and `ReactionBuilder.monod_aerobic_growth`
       delegate to it; delete `templates/stirred_tank/kinetics.py`; repoint
       `templates/stirred_tank/__init__.py`, the builder docstring and any notebook
       imports. Add behaviour-level tests for all eight laws (values at fixed points,
       `make_rate_fn` outputs), pinning `Andrews` and `ContoisAndrews`. Sanity:
       fingerprints bit-identical; note the one edge-case change (zero `Ko2_gL` with
       zero O2 returns 0.0); about 2070 tests.
- [ ] 13. (D7) Move `chemistry/database.py` and `chemistry/databases/` to
       `PyOMES/databases/` with `git mv` (name provisional; confirm at start).
       Repoint imports in `templates/stirred_tank/factory.py`, `models/vlmodels/adm1`,
       four test files, 8 notebooks (cell sources), `raw_construction.py`, the
       tutorials' READMEs; drop the `ChemistryDatabase` export from
       `chemistry/__init__.py`; add the row to `PyOMES/README.md`. Sanity: suite
       green; import guard green.

**Part C — behaviour (isolated commits)**

- [ ] 14. Fix `ChemistryDatabase.extend()` to preserve solver, label and engine
       config; fix the `ReactionSet` docstring. Sanity: +2 tests; stock databases
       unaffected.
- [ ] 15. (D2) Drop the constructor-time validation in `PHController`; delete
       `registry.py`; remove the two tests in `TestCVPHControllerRegistryValidation`;
       add the documentation warning to `PHController`'s `chemical_id` /
       `base_chemical_id` parameters (ids are not validated; must be a
       strong-corrector alias or a species in a declared equilibrium, else the dose
       accumulates as inert); add one line to
       `PHCONTROLLER_CORRECTOR_VALIDATION.md` saying the old check was removed and
       that note is now the replacement. Sanity: −2 tests.
- [ ] 16. Shrink `chemistry.__all__` to names still used; add a `plots` extra
       (matplotlib) to `pyproject.toml` and to `all`.
- [ ] 17. Docs sweep and logging: extend the OPEN_WORK van 't Hoff and constants
       entries with this review's counts; add OPEN_WORK entries for every item in the
       plan doc's "Deferred" list. Full suite; record the final count (about 2070).

## Shipping

- [ ] Full test suite green on the branch (about 2070; see the plan doc)
- [ ] `git checkout main`
- [ ] `git merge --no-ff chemistry-reactions-kinetics-cleanup -m "Merge chemistry-reactions-kinetics-cleanup: <summary>"`
- [ ] `git tag chemistry-reactions-kinetics-cleanup-shipped <commit-hash>`
- [ ] `git push && git push --tags`
- [ ] `git branch -d chemistry-reactions-kinetics-cleanup` and
      `git push origin --delete chemistry-reactions-kinetics-cleanup`
- [ ] Move the plan doc and this checklist to `docs/dev/implementation/shipped/`; add
      "Shipped" banners
- [ ] Update `upcoming/README.md`'s "Recently shipped" list
