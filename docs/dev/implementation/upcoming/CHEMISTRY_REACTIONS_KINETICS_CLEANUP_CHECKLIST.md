# Phase Kickoff Checklist — chemistry-reactions-kinetics-cleanup

> Checklist for [`CHEMISTRY_REACTIONS_KINETICS_CLEANUP.md`](CHEMISTRY_REACTIONS_KINETICS_CLEANUP.md),
> the source of truth for goals, audit and decisions D1-D8; do not restate them
> here. See `README.md`'s "Branching and tagging convention". Modelled on
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

## Pre-flight

- [ ] `git status -sb` clean
- [ ] `git log origin/main..main --oneline` empty
- [ ] Decisions D1-D8 confirmed or amended in the plan doc
- [ ] Branch created off `main`: `git checkout -b chemistry-reactions-kinetics-cleanup`
- [ ] This checklist and the plan doc committed as the first commit

## During

- [ ] Plan doc exists in `docs/dev/implementation/upcoming/`
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
      Sanity: +1 test.
- [ ] 3. Hoist the deferred imports that are not load-bearing: `partition.py`
      (six sites; `_resolve_species` reads `chemistry.common_species` directly),
      `database.py:85`, `equilibria.py:493,549`, `reaction_system.py:171,484`,
      `stoichiometry.py:175`, `kinetic.py:102`, `equilibrium.py:254,321`.
      Leave `thermo_params.py` (deleted in 7) and the `chemical_equilibrium`
      engine imports. Sanity: full suite; import guard green.
- [ ] 4. Remove unused imports and constants (`_M_WATER`, the `IdealGasEOS` import,
      `_NH3`, `Dict`, `warnings`, `field`/`Optional` in `equilibria.py`, `Species` in
      `databases/aqueous.py`, `Sequence` in `plots.py`); fix false docstrings in
      files that survive (`ReactionSet`, `fermenter.*`, `common_species.py:23`,
      `equilibria.py:510`, `reaction_system.py:42-47`, the `Multispecies` docstring).

**Part B — deletions and consolidation (bit-identical)**

- [ ] 5. Delete `PyOMES/kinetics/` and `tests/standalone/test_kinetics.py`; remove
      the `kinetics/` row in `PyOMES/README.md`, the row in the root `README.md`
      test table, and the `kinetics/` line in `docs/architecture.md`. Sanity: −11 tests.
- [ ] 6. Delete `tests/legacy/` (15 files, not collected). Sanity: grep finds no
      reference in docs or CI; suite unchanged.
- [ ] 7. Delete `thermo_params.py` (D3: also remove the unit kwargs and helpers from
      `EquilibriumSet`; delete `_VALID_*`). Drop `thermo=` and the no-op
      `apply_to_cv` call from `build_bsm2_cv` (and its docstring); edit the fixture
      in `test_bsm2_reference.py`; prune `chemistry/__init__.py`. Sanity: BSM2
      sentinels unchanged; suite unchanged.
- [ ] 8. Delete the three unused `EquilibriumSet` presets. Sanity: suite unchanged.
- [ ] 9. Migrate ~50 `HenryPartition`/`RaoultPartition` call sites (7 test files) to
      `HenryEquilibrium`/`RaoultEquilibrium`; delete the aliases and their 6
      tests; fix the `anaerobic_digestion.py` docstring that names the old type.
      Sanity: −6 tests; Henry/Raoult/Ksp fingerprint byte-identical.
- [ ] 10. (D1) Delete `recipe.py`, `chem_recipe.py`, `types.py` and the ion/salt maps
       in `registry.py`; update `chemistry/__init__.py`; update the recipe-layer
       section of `STRONG_ION_INFERENCE_GENERALIZATION.md` to "resolved". Leave
       `COMPOUND_DB`/`resolve_compound`/`validate_compound_id` for 14.
- [ ] 11. (D4) Delete `PyOMES/equilibria/`, the `ThermoFramework.gas_eos` field, its
       `TYPE_CHECKING` import and its test; fix `thermo/liquid_phase_model.py`'s
       docstring reference, `PyOMES/README.md` and `docs/architecture.md` lines.
       Sanity: −1 test.
- [ ] 12. (S3, D5, D6) **Before editing:** freeze the Monod fingerprint as a test,
       covering `Ko2_gL=None` **and** `Ko2_gL` set (the O2 path is not yet
       measured). Then create `reactions/rate_laws.py` with `Monod` (and
       `GrowthKinetics`); make the factory fallback and
       `ReactionBuilder.monod_aerobic_growth` call it, keeping the O2 factor in the
       builder and the operation order unchanged; delete the seven unused laws
       (D6); repoint `templates/stirred_tank/__init__.py` and any notebook
       imports. Sanity: fingerprint bit-identical; +1 test.

**Part C — behaviour (isolated commits)**

- [ ] 13. Fix `ChemistryDatabase.extend()` to preserve solver, label and engine
       config; fix the `ReactionSet` docstring. Sanity: +2 tests; note that stock
       databases are unaffected.
- [ ] 14. (D2) Drop the constructor-time validation in `PHController`; delete
       `registry.py`; remove the two tests in `TestCVPHControllerRegistryValidation`;
       add the OPEN_WORK follow-up on validating dosing ids. Sanity: −2 tests.
- [ ] 15. Shrink `chemistry.__all__` to names still used; add a `plots` extra
       (matplotlib) to `pyproject.toml` and to `all`.
- [ ] 16. Docs sweep and logging: extend the OPEN_WORK van 't Hoff and constants
       entries with this review's counts; add entries for the deferred items in the
       plan doc. Full suite; record the final count.

## Shipping

- [ ] Full test suite green on the branch (expected 2056; see the plan doc)
- [ ] `git checkout main`
- [ ] `git merge --no-ff chemistry-reactions-kinetics-cleanup -m "Merge chemistry-reactions-kinetics-cleanup: <summary>"`
- [ ] `git tag chemistry-reactions-kinetics-cleanup-shipped <commit-hash>`
- [ ] `git push && git push --tags`
- [ ] `git branch -d chemistry-reactions-kinetics-cleanup` and
      `git push origin --delete chemistry-reactions-kinetics-cleanup`
- [ ] Move the plan doc and this checklist to `docs/dev/implementation/shipped/`; add
      "Shipped" banners
- [ ] Update `upcoming/README.md`'s "Recently shipped" list
