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
- [x] Checkpoints tracked below, one commit each
- [ ] **If work stalls:** add a status banner to the top of the plan doc at once

### Checkpoint 1 inventory

- [x] Suite baseline (expected 2072 passed); line counts (chemistry 3,893, kinetics
      308, reactions 2,806, equilibria 555); import-graph and consumer scripts rerun
      and results pasted; searches repeated at branch start.

Recorded 2026-09-21 at `9896752` (working tree clean; branch level with
`origin/chemistry-reactions-kinetics-cleanup`). Searches ran over every tracked
file: 196 `.py`, 32 `.ipynb`, 96 `.md`, 3 `.yml`, 2 `.toml`, 2 `.json`, `.cff`. That
covers `tests/`, `models/`, `.github/` and both notebook generators' string templates.
The gitignored `scratch/` and `notes/` are excluded, and so are the plan doc and this
checklist. The only dynamic import is `tests/run_tests.py:209`, over test module
names; it does not name any candidate. **No decision D1-D8 is affected.**

**Baseline:** `python -m pytest` (configured `testpaths`) = **2072 passed**, 0 failed,
0 skipped, 334 warnings, 2m08s. Line counts (tracked `.py`) all match the plan:
chemistry 3,893, kinetics 308, reactions 2,806, equilibria 555; also
`thermo_params.py` 765 and `chem_recipe.py` 274.

**Import graph** (AST, 196 files, 2,652 imports, tagged top / lazy / `TYPE_CHECKING`):

| Graph | Cycles |
|---|---|
| Module level, top-level imports only, all of `PyOMES/` | **none** |
| Module level, plus lazy imports | 3: `chemistry.equilibria` <-> `chemistry.thermo_params` (lazy at `thermo_params.py:205,290`; the one the plan expects); and two outside these packages, see below |
| Module level, plus `TYPE_CHECKING` | adds one 11-module group (the engines, `reaction_system`, `plots`, `partition`, ...) that closes only through `plots.py:17-18`; not a runtime cycle |
| Package level, top-level imports only | `chemistry` <-> `reactions` (the nine top-level `reactions` imports in the three `databases/` files), and `control` <-> `core` |

Not in the plan and outside these packages, left alone: lazy cycles
`core.control_volume` <-> `core.solvers` (`solvers.py:345`) and `control.param_path`
<-> `core` <-> `core.simulation` (`param_path.py:177`; three in `simulation.py`), and
the package-level `control` <-> `core`. None affects the import guard (checkpoint 2).

**Consumers** (word-boundary counts over the same corpus, 166 top-level names):

| Package | Names | Not used outside own file | Used only inside own package | Used only by tests | Used elsewhere |
|---|---|---|---|---|---|
| chemistry | 124 | 33 | 31 | 17 | 43 |
| kinetics | 8 | 1 | 1 | 4 | 2 (false matches) |
| reactions | 34 | 5 | 8 | 2 | 19 |

The two kinetics "used elsewhere" are false matches: `Environment` (a generic word in
`reactions/__init__.py` and one notebook; nothing imports it from `PyOMES.kinetics`)
and `KineticModel` (only the `PyOMES/README.md` sentence). Tests-only in chemistry
includes the two aliases, `CHEM_DB`, `recipe_to_totals`, `recipe_g_L_to_mol_L` and
the `Multispecies*` classes.

**Delete-candidate searches** (all confirm the plan):

- **`kinetics/`:** only `tests/standalone/test_kinetics.py` (3 import lines, 11 tests)
  imports `PyOMES.kinetics`. `templates/stirred_tank/__init__.py:30` `from .kinetics`
  is the template's own module. Non-code mentions: `PyOMES/README.md:19`,
  `docs/architecture.md:413`, `README.md:153`, `upcoming/README.md:213`.
- **`thermo_params`:** `ThermodynamicConfig` is used only at `bsm2.py:768,872,874-875`
  and `test_bsm2_reference.py:100,103`; nothing reads `_thermo_config`. No caller
  passes `Ka`, `lnKa`, `dH`, `dH_unit`, `T_ref`, `T_unit` or `Kw` to `EquilibriumSet`
  except `ThermodynamicConfig.add_acid`/`set_water` in `thermo_params.py` itself
  (`equilibria.py:291` names it in a docstring). The three presets have no caller;
  `bsm2_diprotic_co2` is referenced only by `adm1_full`.
- **Recipe layer (`recipe.py`, `chem_recipe.py`, `types.py`, `registry.py`):** every
  name appears only in `chemistry/`, `tests/legacy/` and docs, except
  `validate_compound_id` (`control/cv_loops.py:110-118` and one test,
  `test_simulation.py:3070`). `CHEM_DB` has 47 literal keys, 39 unique.
- **`HenryPartition`/`RaoultPartition`:** 51 call sites (37 + 14) in 7 test files; 6
  alias tests (3 + 3 in `test_equilibrium_constraint.py`).
- **`equilibria/`:** imported only by `partition.py:654` (lazy, unused) and
  `thermo/framework.py:23` (`TYPE_CHECKING`); no test imports it.
  `thermo/gas_eos.py` does not exist and the names do not clash.
- **`tests/legacy/`:** 15 files; `testpaths` excludes it; no CI reference.
- **Monod:** `reactions/builder.py:298`, `templates/stirred_tank/factory.py:288` and
  `kinetics.py:167` are where the plan says; `-k Multispecies` selects 11 tests.

**Additions for later checkpoints** (nothing blocks a decision):

- **Checkpoint 5:** `README.md:153` lists `test_kinetics.py` as one item in a table
  cell, not a row. `upcoming/README.md:213` is a dated prose line naming
  `PyOMES/kinetics/`.
- **Checkpoint 6:** root `README.md:90` and `:137` name `tests/legacy/`, so "no
  reference in docs" needs those two lines edited. No CI reference.
- **Checkpoint 9:** `test_equilibrium_constraint_dual_role.py` and
  `test_transfer_models.py` mention the aliases in prose only (no code use).
- **Checkpoint 10:** `PyOMES/README.md:13` advertises "solution recipe builders";
  `tests/data/gas_equilibrated_pH_standards.json:8` mentions "CHEM_DB keys" in a
  prose note.

### Checkpoints

**Part A — hygiene (bit-identical)**

- [x] 1. Baseline and inventory (above).
- [x] 2. Add the import guard test: top-level module graph of `PyOMES/` acyclic
      (function-level and `TYPE_CHECKING` imports excluded). It passes today.
      Sanity: 2073 tests.
      _Notes: done 2026-09-22. New `tests/standalone/test_import_graph_acyclic.py`,
      one test (`test_top_level_import_graph_is_acyclic`). Its own AST walker mirrors
      the checkpoint-1 scratch script: an `ast.NodeVisitor` skips `FunctionDef`/
      `AsyncFunctionDef` bodies (depth-tracked) and `if TYPE_CHECKING:` bodies (the
      `else` branch is kept, since it runs), so only edges that execute at import
      time are collected; relative imports are resolved against each file's own
      dotted module name, and `from pkg import name` resolves to the submodule
      `pkg.name` when one exists, else to `pkg`. A DFS cycle finder (white/gray/black
      colouring) walks the resulting graph. Before asserting on the real graph, the
      test checks the finder itself against two tiny synthetic graphs (one with a
      3-cycle, one without), so a clean result on `PyOMES/` isn't from an untested
      detector. Sanity: full suite **2073 passed**, 0 failed (2072 + this one test);
      import guard green._
- [x] 3. Hoist the deferred imports that are not load-bearing: `partition.py`
      (six sites; `_resolve_species` reads `chemistry.common_species` directly),
      `database.py:85`, `equilibria.py:493,549`, `reaction_system.py:171,484`,
      `stoichiometry.py:175`, `kinetic.py:102`, `equilibrium.py:254,321`.
      Leave `thermo_params.py` (deleted in 7) and the `chemical_equilibrium`
      engine imports. Sanity: full suite; import guard green.
      _Notes: done 2026-09-22, 7 files (23 insertions, 28 deletions). Before
      editing, confirmed by AST scan that `python -c "import PyOMES"` and the
      checkpoint-2 guard both stay meaningful checks: the only lazy imports left
      in these 7 files afterwards are `equilibria.py`'s two deferred `import
      warnings` (stdlib, not in this checkpoint's list) and
      `reaction_system.py`'s two `chemical_equilibrium.engines.*` imports (the
      ones explicitly left alone). `partition.py`'s six sites: `_resolve_species`
      no longer goes through `reactions.stoichiometry._get_common_species` at
      all — it now builds the same `{v.id: v for v in vars(...) if
      isinstance(v, Species)}` lookup directly against the hoisted
      `from . import common_species`, same dict-building logic, same error
      message, so still bit-identical; the other five (two duplicate
      `StoichiometryEntry`, `_parse_stoichiometry`, `vant_hoff_log_K`, the
      already-dead `IdealGasEOS` with its false "avoids circular" comment
      preserved verbatim, since checkpoint 4 removes it) moved to the top import
      block unchanged. `database.py:85`: hoisted; the `TYPE_CHECKING`-only copy
      of the same `ReactionSystem` import (line 29) is now redundant with the
      real one and was dropped so the name isn't bound twice. `equilibria.py`:
      the two `.common_species` imports (5 names, 6 names, 1 overlap-free union)
      merged into one top-level import of all 11 names. `reaction_system.py`:
      both hoisted (`chemistry.species_check`, `.plots`); confirmed `.plots`
      doesn't import matplotlib at its own module level (it's deferred inside
      each plotting function, per its own docstring), so this doesn't force a
      hard matplotlib dependency onto `reaction_system`. `stoichiometry.py:175`:
      hoisted (`_get_common_species`'s own lazy import of `.common_species`,
      separate from partition.py's now-removed use of that function).
      `kinetic.py:102` / `equilibrium.py:254`: the redundant `_infer_elements`
      import merged into the existing top-level `._shared` import already in
      each file, local line deleted. `equilibrium.py:321`: `.plots` hoisted too
      (same matplotlib-safety reasoning); the module-level `plot_vant_hoff`
      import doesn't collide with the class's own `plot_vant_hoff` method (a
      bare name inside the method resolves against the module scope, not the
      class body).
      Traced the resulting load order by hand before editing (`PyOMES/__init__.py`
      reaches `chemistry` before `reactions`, so `chemistry/__init__.py`'s own
      sequence — pulled through by the newly-top-level imports — is what actually
      runs first): every module newly reached this way (`reactions.stoichiometry`,
      `.equilibrium`, `.kinetic`, `.blackbox`, `.environment`, `._shared`,
      `.protocols`, `.builder`, `.reaction_system`, `.plots`,
      `chemistry.species_check`, `PyOMES.equilibria.vle`) bottoms out at leaf
      modules (`chemistry.species`, `chemistry.common_species`, `PyOMES.units`,
      stdlib) with no edge back into `chemistry.partition`/`.database`/`.equilibria`
      — matching the plan's own X1 scratch-experiment finding. Verified
      empirically too: `python -c "import PyOMES"` succeeds, the checkpoint-2
      guard test passes, and the full suite is unchanged at **2073 passed**,
      0 failed (pure code motion, no behaviour change)._
- [x] 4. Remove unused imports and constants (`_M_WATER`, the `IdealGasEOS` import
      and its false comment in `partition.py`, `_NH3`, `Dict`, `warnings`,
      `field`/`Optional` in `equilibria.py`, `Species` in `databases/aqueous.py`,
      `Sequence` in `plots.py`). D8: make the `MultispeciesVLEPartition` docstring
      say it is ideal-only, and replace `_kH_mol_L_atm_from_ref` with one helper
      shared with `HenryEquilibrium` (leave the `101325` literal for the constants
      sweep). Fix false docstrings in files that survive (`ReactionSet`,
      `fermenter.*`, `common_species.py:23`, `equilibria.py:510`,
      `reaction_system.py:42-47`). Sanity: Henry/Raoult/Ksp fingerprint
      byte-identical; the 11 `Multispecies` tests pass.
      _Notes: done 2026-09-22, 8 files. Each of the 7 unused imports/constants was
      confirmed genuinely unused first (an AST + word-boundary scan: only its own
      import/def line, no other reference in that file) before removing it;
      `_CO2`/`_H2O` in `builder.py` stay (both used), only `_NH3` goes. D8: the
      `_kH_mol_L_atm_from_ref` module-level helper (already parametrised, used by
      `MultispeciesVLEPartition`) is now the one implementation;
      `HenryEquilibrium._kH_mol_L_atm` delegates to it instead of repeating the
      same formula — same operations in the same order, so bit-identical by
      inspection, then confirmed by fingerprint (below). The class and method
      docstrings and the removed "avoids circular"/"§CP7" comments no longer
      imply an `IdealGasEOS` object is used (there is no EOS parameter on this
      class at all). **False docstrings:** `ReactionSet` — 4 places, matching the
      checkpoint-1 audit's count exactly: `database.py`'s two (real, factual
      errors — the field's type is `ReactionSystem`) fixed by renaming;
      `reaction_system.py`'s two (historical "replaces the deleted ReactionSet"
      framing) rewritten to describe current design in the present tense, keeping
      the reasoning and dropping the label, along with a `checkpoint 3 of the
      state-unification phase` / `C2` / `C3+` passage that also falsely described
      the deleted `SpeciationPropertySolver` as still owning the engine, and an
      `EQUILIBRIUM_CONSTRAINT_UNIFICATION CP3` pointer (that doc is shipped/
      historical, so dropped per the same test the previous phase used: keep a
      pointer only if canonical, current, too long to summarise, and at a stable
      path — this failed all four). `fermenter.*` — searched the surviving files
      (`chemistry/{species,common_species,species_check,compounds,partition,
      equilibria,database}.py`, `chemistry/databases/`, all of `reactions/`,
      `thermo/`, `templates/stirred_tank/`): none found; the 5 the checkpoint-1
      audit counted are all in `registry.py`/`thermo_params.py`, which don't
      survive (deleted at checkpoints 15/7), so nothing to do here — they leave
      with those files.
      `common_species.py:23`: its doctest imported `StoichiometryEntry` from
      `PyOMES.chemistry`, which doesn't re-export it (confirmed: absent from
      `chemistry/__init__.py`'s `__all__`) and would raise `ImportError`;
      repointed to `PyOMES.reactions.stoichiometry`, which does. `equilibria.py`
      (VFA entries): dropped the false "removed in the PARTITION_MODEL phase"
      claim — confirmed live via a repo search (`chemical_equilibrium/engines/
      bisection/{acid_base,engine}.py` still emit the `{name}_HA`/`{name}_A-`
      keys for entries without `species_refs`, and `S_ac`/`S_pro`/`S_bu`/`S_va`
      right below this comment are exactly such entries) — reworded to state
      what happens now, without inventing a deprecation status nothing tracks.
      **Verification:** a fingerprint of `HenryEquilibrium` (`_kH_mol_L_atm`,
      `log_K`, `partition_ratio`, `equilibrium_a_moles`), `RaoultEquilibrium`,
      `MultispeciesVLEPartition.equilibrium_all_a_moles`, and `KspEquilibrium`
      (via `vant_hoff_log_K`), 4,880 values over a grid of H_ref/dlnH/T_ref/T_K/
      capacities/alphas, run against the current code and against the
      pre-checkpoint-4 `partition.py` (loaded from `git show HEAD:...` as an
      isolated module, same SHA-256 both times): byte-identical. The 11
      `Multispecies` tests pass, and so do `test_partition_model.py`,
      `test_equilibrium_constraint.py`, `test_equilibrium_classification.py`,
      `test_chemistry_database.py` and `test_bsm2_reference.py` individually.
      Full suite unchanged at **2073 passed**, 0 failed._

**Part B — deletions, moves and consolidation (bit-identical)**

- [x] 5. Delete `PyOMES/kinetics/` and `tests/standalone/test_kinetics.py`; remove
      the `kinetics/` row in `PyOMES/README.md`, the row in the root `README.md`
      test table, and the `kinetics/` line in `docs/architecture.md`. Sanity:
      2062 tests.
      _Notes: done 2026-09-22. Fresh repo-wide search first (all tracked file
      types): the only live references were `PyOMES/README.md`'s row,
      `README.md`'s test-table cell, and `docs/architecture.md`'s two-line tree
      entry (row + its `core.py, mapping.py, ...` continuation, removed
      together). Deleted the 6 files under `PyOMES/kinetics/` and
      `tests/standalone/test_kinetics.py` (7 files, matching the checkpoint-1
      inventory exactly). **Plan correction to the checkpoint text:** the root
      `README.md` test table has no separate `kinetics` row — `test_kinetics.py`
      was one of four items in the "Reaction stoichiometry & kinetics" cell (the
      row's other three tests stay; the row's own name is a generic category
      label, not a reference to the deleted package, so it stays too) — flagged
      in the checkpoint-1 note, handled as a cell edit, not a row deletion.
      **Left alone (checked, not an oversight):** `docs/dev/implementation/
      shipped/STIRRED_TANK_TEMPLATE_MIGRATION.md` (shipped/ is off-limits per
      the working rules); this phase's own plan and checklist docs, which
      correctly describe the audit as it was before deletion; `tests/legacy/`'s
      unrelated `v8_split_11.kinetics.core` path (a different, fictional
      package, deleted itself at checkpoint 6); `upcoming/README.md:213`, a
      dated 2026-09-15 changelog entry explaining a past decision, not a
      current-behaviour claim. Verified: `python -c "import PyOMES"` succeeds;
      the import guard test passes; a follow-up search for
      `PyOMES.kinetics`/`PyOMES/kinetics` found only those left-alone,
      historical/plan-doc hits. Full suite **2062 passed**, 0 failed (2073 − 11,
      exactly as predicted)._
- [x] 6. Delete `tests/legacy/` (15 files, not collected). Sanity: grep finds no
      reference in docs or CI; suite unchanged.
      _Notes: done 2026-09-22. Fresh repo-wide search first: 15 files, confirmed
      not in `testpaths` (`pyproject.toml:18`), no `.github/workflows/tests.yml`
      reference, no `tests/run_tests.py` reference, no notebook-generator
      reference. Deleted all 15. Fixed the two root `README.md` mentions flagged
      at checkpoint 1 (`:90`'s repository-layout line dropped `legacy/` from the
      `tests/` list; `:137`'s sentence about non-default test dirs dropped the
      `tests/legacy/` clause, keeping the `tests/performance/` one). **Left
      alone:** `docs/dev/implementation/shipped/EQUILIBRIUM_RESULT.md` (shipped/
      is off-limits) and this phase's own plan/checklist docs, which correctly
      describe the audit as it was before deletion. Follow-up search: zero
      `tests/legacy` or `tests.legacy` hits outside those. Full suite **2062
      passed**, 0 failed — unchanged, since nothing in `tests/legacy/` was ever
      collected._
- [x] 7. (D3) Delete `thermo_params.py`, the unit-conversion arguments of
      `EquilibriumSet.add`/`set_water`, and `_VALID_*`. Drop `thermo=` and the
      no-op `apply_to_cv` call from `build_bsm2_cv` (and its docstring); edit the
      fixture in `test_bsm2_reference.py`; prune `chemistry/__init__.py`. Sanity:
      BSM2 sentinels unchanged; suite unchanged.
      _Notes: done 2026-09-22, 8 files. Deleted `chemistry/thermo_params.py`
      (765 lines) entirely. Before touching `EquilibriumSet`, a fresh repo-wide
      search (`.py`, `.ipynb`) for `Ka=`, `lnKa=`, `dH_unit=`, `T_unit=`, bare
      `dH=`/`T_ref=`/`Kw=` confirmed **zero real callers** pass any of the six
      deleted arguments (`Ka`, `lnKa`, `dH`+`dH_unit`, `T_ref`+`T_unit`, `Kw`) —
      every match was either inside `equilibria.py` itself, an unrelated local
      variable of the same name elsewhere (the bisection engine, an HPLC
      column model), or a different class's constructor
      (`HenryPartition(..., T_ref=...)`), so removing them is bit-identical for
      every existing call site, not just narrower. `EquilibriumSet.add` keeps
      `pKas`, `n_active`, `correction`, `dH_J_per_mol`, `T_ref_K`, `total_key`,
      `species_refs`; `set_water` keeps `pKw`, `correction`, `dH_J_per_mol`,
      `T_ref_K` (also dropped its unused `import warnings`, dead even before
      this edit — never called). `_VALID_K_FORMATS`, `_VALID_DH_UNITS`,
      `_VALID_T_UNITS` and the `_convert_*`/`_Ka_to_pKa`/`_lnKa_to_pKa` helpers
      went with the deleted arguments; `_VALID_CORRECTIONS` (unrelated to
      units — validates `correction`, which stays) and `_R_J` (the van 't Hoff
      gas constant, used by `EquilibriumDef.pKas_at_T`/`WaterDef.Kw_at_T`,
      which stay) moved into `equilibria.py` itself — `_R_J` now `from
      ..units import R_J_PER_MOL_K as _R_J`, the single source, same as
      `thermo_params.py` already did. Rewrote the module docstring's and
      `add`'s doctest examples to use `dH_J_per_mol` instead of `dH`+`dH_unit`
      (same J/mol values, e.g. 7.646/14.9 kJ/mol → 7646.0/14900.0, matching
      `bsm2_default`'s own CO2 entry). Confirmed `PyOMES.chemical_equilibrium.
      engines.bisection.engine._add_water`'s two `eq_set.set_water(...)` calls
      (the only external caller) already use only kept kwargs (`pKw`,
      `correction`, `dH_J_per_mol`, `T_ref_K`). `chemistry/__init__.py`: removed
      the `thermo_params` import block and its 4 `__all__` entries
      (`ThermodynamicConfig`, `validate_thermodynamics`,
      `collect_thermo_params`, `ThermoSnapshot`; `AcidDefinition`/
      `WaterDefinition` were imported but never exported, an existing
      inconsistency that disappears with the file). `models/vlmodels/adm1/
      bsm2.py`: removed `build_bsm2_cv`'s `thermo=` parameter, its docstring
      paragraph, and the 13-line no-op `apply_to_cv` block (both branches —
      confirmed by a fresh search that `cv._thermo_config` has no reader
      anywhere in the repo, only shipped/historical docs). `test_bsm2_reference.py`:
      fixture drops the `ThermodynamicConfig` import and `thermo=` kwarg,
      calling `build_bsm2_cv(rxn_set)` with defaults.
      **Extra fixes found by the fresh search, not explicitly in this
      checkpoint's text:** `docs/architecture.md`'s chemistry/ tree line still
      listed `thermo_params.py` as a live file — removed.
      `docs/dev/implementation/upcoming/STRONG_ION_INFERENCE_GENERALIZATION.md`
      (a live open-question doc, not this phase's own) cited
      `ThermodynamicConfig.compute_CT_cation_from_charge_balance` as
      "unaffected" by a proposed fix — confirmed zero other callers repo-wide,
      then reworded to record that it was deleted here, with no replacement
      yet. `tests/standalone/test_import_graph_acyclic.py`'s own docstring
      named `chemistry.equilibria <-> chemistry.thermo_params` as "the one
      real cycle" — a fresh module-graph rescan (top+lazy, all of `PyOMES/`)
      confirmed that cycle is gone (only the same two out-of-scope cycles from
      checkpoint 1 remain: `core.control_volume <-> core.solvers` and
      `control.param_path <-> core <-> core.simulation`), so the docstring now
      cites the still-current `core.control_volume` example instead.
      `OPEN_WORK.md:346`'s mention is a dated, past-tense record of the
      2026-07-02 gas-constant re-baseline and is left as historical.
      **Verification:** `python -c "import PyOMES"` succeeds; the import guard
      test passes; `test_bsm2_reference.py` (BSM2 sentinels),
      `test_equilibrium_constraint.py`, `test_equilibrium_classification.py`,
      `test_chemistry_database.py` — 103 passed together. Full suite unchanged
      at **2062 passed**, 0 failed._
- [x] 8. Delete the three unused `EquilibriumSet` presets. Sanity: suite unchanged.
      _Notes: done 2026-09-22. Fresh search confirmed zero external callers of
      `bsm2_diprotic_co2`, `bsm2_with_sulfide`, `adm1_full` (only `adm1_full`'s
      own now-deleted call to `bsm2_diprotic_co2`, and `shipped/`/plan-doc
      mentions). Deleted all three (43 lines). Their deletion left
      `H3PO4`, `H2PO4_minus`, `HPO4_2minus`, `PO4_3minus`, `HSO4_minus`,
      `SO4_2minus` unused in this file (only `adm1_full` referenced them), so
      trimmed the `common_species` import to the five names `bsm2_default`
      still uses. One file changed. Verified: `python -c "import PyOMES"`
      succeeds; import guard green; an AST unused-import check on the file
      found nothing left over. Full suite unchanged at **2062 passed**,
      0 failed._
- [x] 9. Migrate ~50 `HenryPartition`/`RaoultPartition` call sites (7 test files) to
      `HenryEquilibrium`/`RaoultEquilibrium`; delete the aliases and their 6
      tests; fix the `anaerobic_digestion.py` docstring that names the old type.
      Sanity: 2056 tests; Henry/Raoult/Ksp fingerprint byte-identical.
      _Notes: done 2026-09-22, 9 files. Fresh search found 9 test files
      mentioning the alias names, not 7 — but 2 of them
      (`test_equilibrium_constraint_dual_role.py`, `test_transfer_models.py`,
      already flagged at checkpoint 1) are prose-only, no code use, confirmed
      again here; the 7 with real constructor calls total exactly **51**
      (matching the checkpoint-1 count), not ~105 raw textual occurrences (most
      of the difference is this codebase's style of a fresh `from PyOMES.chemistry
      import HenryPartition` inside every test method, plus the alias tests'
      own dual imports). Deleted `TestHen
      ryPartitionAlias`/
      `TestRaoultPartitionAlias` (6 tests) from `test_equilibrium_constraint.py`
      first — including the module docstring line describing them and the
      `import warnings`, which nothing else in that file used — so the
      remaining rename couldn't collide with a test that intentionally imports
      both names. Then a word-boundary rename (`HenryPartition`→
      `HenryEquilibrium`, `RaoultPartition`→`RaoultEquilibrium`) across the
      other 6 files: every `from PyOMES.chemistry import HenryPartition` line,
      every constructor call, every `-> HenryPartition` return-type annotation
      on local `_hp()` fixtures, and every docstring/comment describing the
      test file's own current subject (e.g. "Tests for PartitionModel protocol
      and HenryPartition (C2)." → "...HenryEquilibrium (C2)."). Checked for
      duplicate-name imports and re-parsed all 6 files: clean. Deleted the two
      `HenryPartition`/`RaoultPartition` alias functions from `partition.py`
      (`Any` and `import warnings` then became unused there too — only used by
      the deleted functions — removed both); dropped both names from
      `chemistry/__init__.py`'s import and `__all__`. Fixed
      `anaerobic_digestion.py`'s docstring per the checkpoint text (the old
      type name inside an otherwise-historical "Previously these were..."
      sentence). **Left alone, deliberately:** the same old-name mention in
      `test_equilibrium_constraint_dual_role.py`/`test_transfer_models.py`
      (genuinely historical "before this phase"/"previously" framing, not
      named in this checkpoint's text, unlike `anaerobic_digestion.py`); every
      `docs/dev/ideas/*.md` and `shipped/*.md` mention (design proposals and
      shipped history, both off-limits); `upcoming/README.md`'s four dated
      changelog entries. **Verification:** `python -c "import PyOMES"`
      succeeds; import guard green; a 4,880-value fingerprint of
      `HenryEquilibrium`/`RaoultEquilibrium`/`KspEquilibrium`/
      `MultispeciesVLEPartition` against the pre-checkpoint `partition.py`
      (from `git show HEAD:...`) is byte-identical (same SHA-256 as
      checkpoint 4's) — expected, since the aliases were pure pass-through
      wrappers and nothing in the fingerprint ever called them. The 7 migrated
      files: 261 passed, no `DeprecationWarning`s left. Full suite **2056
      passed**, 0 failed, exactly as predicted (2062 − 6)._
- [x] 10. (D1) Delete `recipe.py`, `chem_recipe.py`, `types.py` and the ion/salt maps
       in `registry.py`; update `chemistry/__init__.py`; update the recipe-layer
       section of `STRONG_ION_INFERENCE_GENERALIZATION.md` to "resolved, see the
       cleanup plan". Leave `COMPOUND_DB`/`resolve_compound`/`validate_compound_id`
       for checkpoint 15.
       _Notes: done 2026-09-22, 8 files, 744 lines removed. Fresh repo-wide search
       (`.py`, `.ipynb`) before deleting: zero consumers anywhere of
       `SolutionRecipe`, `AqueousTotals`/`AqueousTotalsUser`/`AqueousEquilibrium`,
       `SALT_DISSOCIATION_MAP`, `ION_TO_ENGINE_KEY`, `normalize_ion_label`,
       `ion_to_engine_key`, `map_user_ions_to_engine`, `ChemSpec`, `CHEM_DB`,
       `recipe_to_totals`, `recipe_g_L_to_mol_L`, outside the files being deleted
       and `chemistry/__init__.py`'s own import/export lines — matching the audit.
       Also deleted `validate_compound_ids` (plural): the checkpoint text names
       only 3 survivors (`COMPOUND_DB`, `resolve_compound`,
       `validate_compound_id` singular), and the plural form had zero callers,
       matching its "no consumer" classification in the original W1 audit
       (distinct from the singular, which `control/cv_loops.py` and one test
       still use and which stays). `registry.py`'s module docstring described
       three purposes (strong-ion inference, ion-label normalization, mapping to
       engine keys) — all three now deleted — so rewrote it to describe what
       actually remains (the compound → molar-mass/composition lookup); also
       fixed two runtime message strings (`resolve_compound`'s `KeyError`,
       `validate_compound_id`'s `UserWarning`) that still said
       `fermenter.chemistry.registry.COMPOUND_DB` — the pre-rename project name,
       missed at checkpoint 4 because I'd read `registry.py` as not surviving at
       all, when only its ion/salt-map portion doesn't; no test asserts the old
       string. `chemistry/__init__.py`'s own docstring described the package as
       providing "*typed* inputs/outputs" — `types.py`'s own stated purpose,
       word for word — so rewrote it to describe the broader surviving scope
       (species, partition models, equilibrium sets, the compound registry).
       `STRONG_ION_INFERENCE_GENERALIZATION.md`: appended a dated resolution
       note to the "Adjacent, out of scope" section (matching that doc's own
       established pattern of dated update notes; the design-discussion prose
       above it is left as the historical record). Two further doc fixes flagged
       at checkpoint 1: `PyOMES/README.md`'s chemistry/ row dropped "solution
       recipe builders"; the orphaned `tests/data/gas_equilibrated_pH_standards.json`
       (zero Python readers, confirmed again here) had its "CHEM_DB keys" note
       reworded to not cite the deleted name. **Left alone:** the same file's
       mention of `run_test_NIST_buffer_standards_with_recipe.py`, a script that
       already doesn't exist in the repo — a separate, pre-existing dangling
       reference, not caused by D1. **Verification:** JSON re-parses;
       `python -c "import PyOMES"` succeeds; import guard green; an AST
       unused-import check on `registry.py`/`chemistry/__init__.py` found
       nothing real (the `__init__.py` "hits" are all names used only via
       `__all__`, the normal pattern for a package init); `validate_compound_id`/
       `PHController` tests (21) and the full suite (**2056 passed**, 0 failed)
       both unchanged — nothing tested the deleted code._
- [x] 11. (D4) Move `equilibria/` to `thermo/gas_eos.py` with `git mv`: `GasEOS`,
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
       _Notes: done 2026-09-22. `git mv PyOMES/equilibria/vle.py
       PyOMES/thermo/gas_eos.py` (clean rename, git tracked it); `GasEOS` and
       `IdealGasEOS` moved as-is, `HenryIdealVLE` dropped (fresh search first:
       zero consumers anywhere, confirming the audit). `peng_robinson.py`'s
       content (`CriticalProperties`, `BIOGAS_SPECIES`, `BIOGAS_KIJ`,
       `PengRobinsonEOS` and its helpers) merged into the same file by hand
       (two files into one has no single `git mv` target, so this half loses
       rename tracking; `vle.py`'s larger content keeps it) — code unchanged
       except one shared `_R_L_ATM_PER_MOL_K` import (`as R`) replacing the
       two files' separate aliases of the same constant, and a merged module
       docstring. Then `git rm` the two-file remainder (`peng_robinson.py`,
       `__init__.py`); the package directory carries no tracked files now.
       `ThermoFramework.gas_eos` is `Optional[GasEOS]` with a real top-level
       import (`gas_eos.py` only depends on `..units`, a leaf, so no cycle);
       the now-empty `if TYPE_CHECKING:` block and its now-unused
       `TYPE_CHECKING` import removed. `thermo/__init__.py` exports all six
       names. **New tests** (`tests/standalone/test_gas_eos.py`, 6 — this
       module had zero test coverage before, per the checkpoint-1 audit):
       4 Z-factor pins at fixed, reproducible `(n, V, T=308.15 K)` inputs
       (chosen by bisection to hit the stated pressures exactly), each
       asserting both the precise computed value (`rel=1e-9`, a regression
       pin) and the checkpoint's approximate figure (`abs=0.01`, confirming
       the pin is in the right ballpark) — CO2 20 atm: Z=0.8946 (target
       ≈0.897); N2 50 atm: Z=0.9898 (≈0.987); CH4 50 atm: Z=0.9111 (≈0.910);
       a 60/40 CH4/CO2 mix at 1 atm: Z=0.9971 (≈1, the ideal-gas limit) — all
       close to, not exactly, the audit's approximate figures (a different
       plausibility script, most likely a different exact molar quantity for
       the same target pressure); 2 tests contrasting
       `PengRobinsonEOS.partial_pressures_atm` (fugacities) against
       `IdealGasEOS.partial_pressures_atm` (partial pressures): at ~50 atm
       CH4 they differ by ~17% (fugacity 45.50 vs partial pressure 54.88);
       at ~1 atm (the same mixture) they agree to within 1%, the expected
       φ→1 low-pressure convergence. **Docs fixed:** `peng_robinson.py`'s two
       `fermenter.equilibria.vle.GasEOS`-style docstring cross-references
       (old project name and path) resolved by merging into one file (no
       cross-reference needed); `thermo/liquid_phase_model.py:4`'s path
       reference; `PyOMES/README.md` (dropped the `equilibria/` row, folded
       gas EOS into the `thermo/` row); `docs/architecture.md` (the
       "Equilibrium pathways" section's item 2, presenting
       `ProcessCoupledEquilibrator`/`HenryEquilibriumInterface` as a live
       second pathway, was already false before this checkpoint — that
       package was CUFermenter-era dead code, per `equilibria/__init__.py`'s
       own docstring — removed rather than re-pathed; the tree's `equilibria/`
       entry removed, `thermo/` entry gains a `gas_eos.py` line).
       **`OPEN_WORK.md` updated too** (not explicitly in this checkpoint's
       text, but its own tracked "`docs/architecture.md` still describes
       deleted CUFermenter-era code" entry cites exactly these two spots by
       line number): added a dated partial-fix note, and corrected the "No
       longer exist" bullet's now-stale carve-out ("except `vle.py` and
       `peng_robinson.py`" — both are gone too now) and the gas-constant
       entry's `equilibria/vle.py` path mention. **Verification:**
       `python -c "import PyOMES"` succeeds; import guard green; an AST
       unused-import check on all four touched `thermo/` files found nothing
       real (the `__init__.py` "hits" are names used only via `__all__`, as
       in earlier checkpoints); `test_thermo_framework.py` (16, including
       `test_gas_eos_none_by_default`) and the new `test_gas_eos.py` (6) all
       pass. Full suite **2062 passed**, 0 failed (2056 + 6 new tests, exactly
       as predicted)._
- [x] 12. (D5, D6) **Before editing:** freeze the rate-function fingerprints as
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
       _Notes: done 2026-09-22, 7 files (6 edited + 1 new test file). Before any
       edit, captured a 2000-point grid fingerprint (SHA-256 of packed doubles) of
       `ReactionBuilder.monod_aerobic_growth`'s rate_fn, both the plain-Monod
       path and the `Ko2_gL` path, against the pre-refactor code. `git mv
       templates/stirred_tank/kinetics.py reactions/rate_laws.py` (clean rename;
       this alone satisfies "delete templates/stirred_tank/kinetics.py" — move
       and delete are the same operation here). All 8 classes moved unchanged
       except the header docstring (reworded reactor-agnostic, per D5: "public,
       reactor-agnostic") and 5 imports (`field`, `asdict`, `Any`, `Dict`,
       `Optional`) that were already completely unused before this move — a
       pre-existing issue, not caused by D5/D6, but cheap to fix while already
       rewriting this exact file's header; verified pre-existing via
       `git show HEAD:...` before touching. `reactions/__init__.py` now exports
       all 9 names (the protocol + 8 laws) — a new, natural home alongside
       `ReactionBuilder`/`ReactionSystem`; `rate_laws.py` has zero internal
       PyOMES imports (pure stdlib), so no cycle risk either direction.
       `templates/stirred_tank/__init__.py` re-export kept (still part of that
       package's public fluent-builder API, not a deprecated shim), repointed to
       `PyOMES.reactions.rate_laws` (absolute, matching this exact file's own
       existing style for cross-package refs). `StirredTankBuilder.substrate`'s
       docstring ("the builder docstring") now shows the canonical import path,
       noting the templates re-export. Zero notebooks reference any of this
       (fresh search: only `rate_laws.py`'s own docstring matched).
       **Delegation:** `ReactionBuilder.monod_aerobic_growth`'s inline `_rate_fn`
       closure (formerly `builder.py:298`, the audit's third Monod
       implementation) replaced with `Monod(...)`/`DualSubstrateMonod(
       secondary_in_mol_L=False, secondary_MW=32.0, ...)` + `.make_rate_fn(...)`,
       chosen by whether `Ko2_gL is None`; verified against the captured
       fingerprint: byte-for-byte `==` match (not `approx`) on both paths, and
       the new edge case (`Ko2_gL=0.0`, O2=0) returns `0.0` where the old
       closure raised `ZeroDivisionError`, exactly as D5 specifies. The
       factory's default-Monod fallback closure (`factory.py:288`, the second
       implementation) replaced with `Monod(...).make_rate_fn(...)` the same
       way; verified separately (formula comparison plus an empirical 5,000-point
       check against the pre-edit closure, `==` on every point). **New tests**
       (`tests/standalone/test_rate_laws.py`, 36 — more than the plan's ~10
       estimate, since this code had zero coverage before this phase and
       D5/D6 both ask for thorough pinning): fixed-point `mu()` values for all
       8 laws (hand-verified by simple arithmetic, e.g. Monod at S=Ks gives
       μ_max/2); `Andrews`/`ContoisAndrews` checked against a frozen Haldane-form
       reference function (μ = μ_max·r/(Ks+r+r²/Ki)) across 3-2 parameter sets ×
       6-7 values each, plus a peak-below-μ_max check and a
       biomass-independence/scale-invariance check — the explicit "reference for
       a later composable-inhibition redesign" the D5/D6 note asks for;
       `make_rate_fn` extensive-rate (mol/h) pins for `Monod` and
       `DualSubstrateMonod`; the frozen-reference fingerprint (embedded copy of
       the pre-refactor closure, not just a one-time script) comparing the
       *current* `ReactionBuilder.monod_aerobic_growth` against it on a 500-point
       log-uniform grid for both paths, `==` exact, plus the edge-case test
       (current code returns 0.0; the frozen reference still raises
       `ZeroDivisionError`, proving the test would catch a regression either
       way). **Verification:** `python -c "import PyOMES"` succeeds; import
       guard green; an AST unused-import check across all 6 touched files found
       nothing left over (`factory.py`'s `SimulationConfig` and
       `templates/stirred_tank/builder.py`'s `Callable` were flagged too but
       confirmed pre-existing and unrelated via `git show HEAD:...`, left alone);
       `test_builder.py`/`test_configs.py` (86, StirredTank templates, exercises
       the factory's real `create_volume` path) and the new
       `test_rate_laws.py` (36) all pass. Full suite **2098 passed**, 0 failed
       (2062 + 36 new tests)._
- [x] 13. (D7) Move `chemistry/database.py` and `chemistry/databases/` to
       `PyOMES/databases/` with `git mv` (name provisional; confirm at start).
       Repoint imports in `templates/stirred_tank/factory.py`, `models/vlmodels/adm1`,
       four test files, 8 notebooks (cell sources), `raw_construction.py`, the
       tutorials' READMEs; drop the `ChemistryDatabase` export from
       `chemistry/__init__.py`; add the row to `PyOMES/README.md`. Sanity: suite
       green; import guard green.
       _Notes: done 2026-09-22, 20 files (5 renamed, 15 edited). Name confirmed:
       `PyOMES/databases/` — no conflict (fresh check: the path didn't exist),
       matches the plan's target layout, no reason to deviate. `git mv`
       `chemistry/database.py` and the 4 files under `chemistry/databases/`
       (clean renames); `chemistry/databases/`'s own directory is now empty of
       tracked files. Internal imports fixed for the new depth: `database.py`'s
       `..reactions.reaction_system` import is unchanged (same depth, `chemistry/`
       and `databases/` are both direct `PyOMES/` children), but its
       `TYPE_CHECKING` imports of `.species`/`.partition` become
       `..chemistry.species`/`..chemistry.partition` (now cross-package); the
       three stock-database files' `..common_species`/`..species`/`..partition`
       become `..chemistry.*` and their `...reactions.*`/`...thermo.*` (three
       dots, from the old depth-2 `chemistry.databases.X`) become `..reactions.*`/
       `..thermo.*` (two dots, from the new depth-1 `databases.X`); their
       sibling imports (`.database`, `.aqueous`, `.bioprocess_basic`) are
       unchanged. Docstring usage examples in all 5 files updated to the new
       import path. `chemistry/__init__.py`: dropped `from .database import
       ChemistryDatabase` and its `__all__` entry.
       **External consumers** (fresh search first, all file types):
       `templates/stirred_tank/factory.py`'s two import lines and one docstring
       cross-reference; the 4 test files (`test_chemistry_database.py` — by far
       the most, including 2 combined `from PyOMES.chemistry import
       ChemistryDatabase, Species`-style lines split into two clean imports each
       — `test_cv_advance.py`, `test_equilibrium_constraint_dual_role.py`,
       `test_nr_tableau_gas_liquid.py`); all 8 notebooks (2 each in
       `Example1_mtp_well.ipynb` and `chemistry_database.ipynb`, 1 each in the
       other 6, matching the plan's count exactly) — edited as raw bytes (not
       text-mode) after a first attempt corrupted line endings across three
       CRLF-normalized notebooks (a text-mode read translates CRLF→LF, and
       writing back non-translating produced 300+ line diffs instead of the
       expected 1-2; caught by reviewing `git diff --stat` before it went
       further, reverted with `git checkout --`, redone in binary mode: exactly
       9 changed lines across all 8 notebooks, verified JSON-valid and
       source-only, no saved output touched); `raw_construction.py`'s one
       docstring path mention (not an import — the file imports nothing from
       this package). **`docs/tutorials/reactions/README.md`** (the one
       tutorials' README that mentions the topic) only names the `ChemistryDatabase`
       class, no module path — needed no change. **`models/vlmodels/adm1`**:
       fresh search (whole directory, all files) found no real import at all —
       only a path-free comment in `base.py` mentioning "ChemistryDatabase
       rollout" in prose; the plan's expectation of an import to repoint didn't
       hold here, noted rather than treated as a blocker (everything else about
       D7 and the target layout is unaffected). Added the `databases/` row to
       `PyOMES/README.md`, and, since already editing the adjacent `reactions/`
       row, folded in its own missing mention of `rate_laws.py` (added in
       checkpoint 12, never documented there). `docs/architecture.md` has no
       `chemistry.database`/`ChemistryDatabase` mention at all — nothing to fix.
       **Verification:** `python -c "import PyOMES"` and direct imports of all
       4 `databases.*` modules succeed; import guard green; an AST unused-import
       check found nothing real (TYPE_CHECKING-only annotation imports in
       `database.py`, and `factory.py`'s pre-existing, already-confirmed-unrelated
       `SimulationConfig`, both false positives). The specifically affected
       files — `test_chemistry_database.py`, `test_cv_advance.py`,
       `test_equilibrium_constraint_dual_role.py`, `test_nr_tableau_gas_liquid.py`,
       `test_builder.py`, `test_configs.py` — 218 passed. Notebooks were not
       executed (not asked for by this checkpoint's sanity check), only verified
       structurally: JSON-valid, source-only diff, correct new paths. Full suite
       unchanged at **2098 passed**, 0 failed._

**Part C — behaviour (isolated commits)**

- [x] 14. Fix `ChemistryDatabase.extend()` to preserve solver, label and engine
       config; fix the `ReactionSet` docstring. Sanity: +2 tests; stock databases
       unaffected.
       _Notes: done 2026-09-22, 2 files. `extend()` now reads the existing
       system's `label`, `._solver` and `._engine_config` (there is no public
       `solver` accessor; `tests/standalone/test_nr_speciation_engine.py`
       already reads `._solver` directly, so this matches existing
       within-repo practice) and passes them to the merged
       `ReactionSystem` — `label=`/`solver=` at construction (the only
       place `solver` can be set), then `.configure_engine(**engine_config)`
       for `use_activity`/`activity_model`. When this database has no
       existing reactions (`self.reactions is None`), the merged system
       falls back to `ReactionSystem`'s own defaults, since there is nothing
       to inherit — documented in the docstring and covered by its own test.
       `ReactionSet` docstring: already fixed at checkpoint 4 (both mentions
       in this file were renamed to `ReactionSystem` then); a fresh check
       found none left, so nothing to do here — the checkpoint text
       restates a completed item. **New tests** (3, not the estimated 2):
       `test_extend_preserves_solver_and_label`,
       `test_extend_preserves_engine_config`, and
       `test_extend_with_no_existing_reactions_uses_defaults` (the fallback
       path, worth pinning since I added it to the docstring). **Stock
       databases unaffected, verified two ways:** all 36 tests in
       `test_chemistry_database.py` pass unchanged, including every
       `TestStockDatabases`/`TestStockDatabasePartitionModels` test; and
       directly inspecting `AQUEOUS_DEFAULT`/`BIOPROCESS_BASIC`/`AD_BASIC`
       (built via a real `.extend()` chain) after the fix shows all three
       still carry the plain defaults (`label=""`, `solver="charge_balance"`,
       `engine_config={"use_activity": False, "activity_model": "davies"}`)
       — expected, since nothing in the aqueous → bioprocess_basic → AD_BASIC
       chain ever sets `label=`/`solver=`/`configure_engine()` to anything
       else, so the fix has nothing non-default to propagate for them; the
       audit's reproduction case (`solver="newton_raphson"` dropped to
       `"charge_balance"`) had no current caller. **Verification:**
       `python -c "import PyOMES"` succeeds; import guard green; no real
       unused imports (same TYPE_CHECKING-only false positives as
       checkpoint 13). Full suite **2101 passed**, 0 failed (2098 + 3)._
       
- [x] 15. (D2) Drop the constructor-time validation in `PHController`; delete
       `registry.py`; remove the two tests in `TestCVPHControllerRegistryValidation`;
       add the documentation warning to `PHController`'s `chemical_id` /
       `base_chemical_id` parameters (ids are not validated; must be a
       strong-corrector alias or a species in a declared equilibrium, else the dose
       accumulates as inert); add one line to
       `PHCONTROLLER_CORRECTOR_VALIDATION.md` saying the old check was removed and
       that note is now the replacement. Sanity: −2 tests.
       _Notes: done 2026-09-22, 6 files. Deleted `PHController.__post_init__`
       entirely — its whole body was the two `validate_compound_id` calls, so
       nothing else was lost. Documentation warning added to both
       `chemical_id` and `base_chemical_id` in the class's Attributes
       docstring, naming the real failure mode (silent inert dose) and
       pointing at `PHCONTROLLER_CORRECTOR_VALIDATION.md` for the proper
       check. Fresh repo-wide search (all file types) before deleting
       `registry.py`: its only consumer was the just-removed
       `__post_init__` call and `chemistry/__init__.py`'s export — `COMPOUND_DB`,
       `resolve_compound`, `_COMPOUND_ALIASES` had zero other callers anywhere,
       confirming the whole file (not just `validate_compound_id`) was safe to
       delete, exactly as checkpoint 1's inventory anticipated ("Leave
       COMPOUND_DB/resolve_compound/validate_compound_id for checkpoint 15").
       `chemistry/__init__.py`: dropped the import/exports; its module
       docstring said "the compound registry" ambiguously (registry.py's
       deleted `COMPOUND_DB` vs. `compounds.py`'s still-live
       `ChemicalRegistry`) — reworded to name `ChemicalRegistry` explicitly
       so the still-true half of that sentence doesn't read as newly false.
       Removed `TestCVPHControllerRegistryValidation` (2 tests) from
       `test_simulation.py`; a repo-wide check found no other test asserting
       a construction-time warning from `PHController`. Added the one line
       to `PHCONTROLLER_CORRECTOR_VALIDATION.md` (as a dated update note
       under its existing Status callout, matching that doc's own style) —
       and, since `STRONG_ION_INFERENCE_GENERALIZATION.md`'s checkpoint-10
       note explicitly said `COMPOUND_DB`/`resolve_compound`/
       `validate_compound_id` were "left for checkpoint 15 to decide on,"
       added a matching follow-up note there too, closing that loop.
       **Verified:** `python -c "import PyOMES"` succeeds; import guard
       green; constructing `PHController(chemical_id="NotARealCompound_XYZ")`
       under `warnings.simplefilter("error")` no longer raises (previously
       would have, on the old code path); no real unused imports (`cv_loops.py`
       clean; `chemistry/__init__.py`'s hits are the same `__all__`-only
       false positives as before). `test_simulation.py` (309 tests) passes.
       Full suite **2099 passed**, 0 failed (2101 − 2, exactly as predicted)._
- [x] 16. Shrink `chemistry.__all__` to names still used; add a `plots` extra
       (matplotlib) to `pyproject.toml` and to `all`.
       _Notes: done 2026-09-22, 3 files. Of the 14 names `chemistry/__init__.py`
       exported (well down from the original 37, since D1/D2/D4/D7 already
       deleted most of the churn), a fresh check of every real
       `from PyOMES.chemistry import ...` line in the repo (all file types)
       found 4 never used that way: `Chemical`, `ChemicalRegistry` (always
       imported from `.compounds` directly — even `PyOMES/__init__.py` itself
       does this), `EquilibriumSet` (always from `.equilibria`), and
       `EquilibriumDef` (never imported at all outside its own definition
       file; every other hit was a docstring/comment naming the class, not an
       import). Removed all 4 from both the import statements and `__all__` —
       dropping only the `__all__` entry would have left a dead, undocumented
       side door (`from PyOMES.chemistry import Chemical` would still have
       worked, just not via `__all__`); no `from PyOMES.chemistry import *`
       exists anywhere to make that distinction matter either way. Left the
       module docstring's "compound database"/"equilibrium sets" wording
       alone — it describes what the subpackage's submodules provide, which
       is still accurate, not specifically the root `__all__`. **`plots`
       extra:** added `plots = ["matplotlib"]` to `pyproject.toml` and
       `"matplotlib"` to `all`; left `test` alone (`plots.py` has no pytest
       coverage per the checkpoint-1 audit, so CI doesn't need it). Also
       added a `plots` row to root `README.md`'s "Optional Features" table
       and its `pip install -e ".[plots]"` example line, alongside the
       `pyproject.toml` change the checkpoint named, for discoverability.
       **Verification:** `python -c "import PyOMES"` succeeds; import guard
       green; `from PyOMES.chemistry import Chemical` now raises
       `ImportError` (confirmed); `pyproject.toml` re-parses as valid TOML
       with the expected extras. Full suite unchanged at **2099 passed**,
       0 failed (no tests added, removed, or behaviour-changed)._
- [x] 17. Docs sweep and logging: extend the OPEN_WORK van 't Hoff and constants
       entries with this review's counts; add OPEN_WORK entries for every item in the
       plan doc's "Deferred" list. Full suite; record the final count (about 2070).
       _Notes: done 2026-09-22, `OPEN_WORK.md` only, no code changes. Extended
       3 existing entries: the van 't Hoff copy entry now reports the plan's
       "at least nine down to seven" estimate alongside what this review could
       concretely verify (thermo_params.py deletion removed 2 at checkpoint 7;
       partition.py's D8 delegation collapsed one pair into one implementation
       at checkpoint 4) and an independent recount of 8 distinct
       `math.exp`-shaped temperature-correction implementations across
       `chemistry/`, `reactions/`, `thermo/` — close to, not exactly matching,
       the plan's "seven", left unreconciled rather than forced to agree; the
       constants-sweep entry notes which files the 2026-09-20 survey counted
       literals in have since moved or been deleted (thermo_params.py,
       recipe.py/chem_recipe.py, kinetics.py→rate_laws.py,
       equilibria/→gas_eos.py, database.py/databases/→PyOMES/databases/), that
       this phase was pure-refactor for those literals so aggregate totals
       should be roughly unchanged, plus a spot-check of 18 occurrences of
       `101325`/`298.15` across the moved/new files; the monitor-plumbing
       entry notes `ReactionSystem._conservation_monitor` is also write-only
       (confirmed via `grep -n "_conservation_monitor\." reaction_system.py`
       returning empty), with `core/control_volume.py` both setting it via
       `attach_conservation_monitor()` and separately reading its own
       `cv._conservation_monitor` directly — the actual source of the
       `ConservationWarning`s seen throughout this phase's test runs. Added
       6 new entries, one per remaining "Deferred" bullet (the plan's own
       last bullet — molar-mass unification, `EquilibriumSet` location,
       `plot_vant_hoff` — kept bundled as one entry, matching how the plan
       groups them): a species-based replacement for the deleted recipe
       layer (D1), superseding the recipe-layer question in
       `STRONG_ION_INFERENCE_GENERALIZATION.md`; the `partial_pressures_atm`
       ideal-vs-fugacity split and `ThermoFramework.gas_eos` being read by
       nothing in production, with the 7 hard-coded ideal-gas-law call sites
       named (D4); mapping-based parameters and composable inhibition for
       the rate laws, including why the Haldane form (`Andrews`,
       `ContoisAndrews`) can't be a plain multiplicative wrapper (D5/D6);
       the remaining package-level `chemistry`<->`reactions` cycle via
       `partition.py`'s three `*Equilibrium` classes needing
       `StoichiometryEntry`/`vant_hoff_log_K` (D7, module-level graph still
       proven acyclic by `test_import_graph_acyclic.py`); the W10
       yield-achievability gap in `ReactionBuilder.aerobic_growth`; and the
       three bundled loose ends. **Verification:** full suite re-run in the
       background, **2099 passed**, 0 failed, 166 warnings, 116.64s —
       unchanged from checkpoint 16 as expected (docs-only checkpoint); the
       plan's "about 2070" was an estimate made before the exact
       checkpoint-16 count (2099) was known, so 2099 stands as the final
       count for this phase, not 2070._
- [x] 18. Move `_ATOMIC_WEIGHTS` from `chemistry/species.py` to `PyOMES/units.py`;
       add a design note for three ambient species-resolution fallbacks found
       along the way.
       _Notes: done 2026-09-22 at `21ff78d`. Surfaced during an exploratory,
       conversational review of whether `common_species.py`/`compounds.py`/
       `_ATOMIC_WEIGHTS` belong in `PyOMES/databases/` instead of `chemistry/`
       (not part of the original D1-D8 audit; additive scope). `_ATOMIC_WEIGHTS`
       is a fixed IUPAC 2021 physical-constants table (18 entries), not
       domain-specific reaction/species data — `units.py`'s own docstring
       already states its purpose as "the single, authoritative definition of
       shared constants," matching `R_J_PER_MOL_K` etc. Moved as
       `ATOMIC_WEIGHTS` (public, no leading underscore, matching `units.py`'s
       naming convention); `species.py`'s `__post_init__` and two docstring
       pointers updated to match; `common_species.py`'s docstring pointer fixed
       too (it named `species.py` as the table's home). Confirmed zero other
       consumers repo-wide before moving. Separately, the same review found
       `reactions/stoichiometry.py`, `chemistry/partition.py` and
       `core/control_volume.py` all resolve unrecognized species ids by
       scanning `common_species.py`'s entire module namespace via `vars()`,
       not from anything the model itself declared — logged as
       [`EXPLICIT_SPECIES_RESOLUTION.md`](EXPLICIT_SPECIES_RESOLUTION.md)
       (design note only, no code change; not part of this phase).
       **Verification:** `test_species.py` + `test_species_check.py`
       (25 passed); `test_import_graph_acyclic.py` green;
       `Species(id="CO2", atoms={"C":1,"O":2}).MW == 44.009` confirmed by hand.
       Full-suite confirmation folded into checkpoint 19 below (no tests
       added/removed/changed by this checkpoint on its own)._
- [x] 19. Move `chemistry/compounds.py` to `PyOMES/compounds.py`.
       _Notes: done 2026-09-22 at `444e651`. Same review: `ChemicalRegistry`/
       `Chemical` have zero imports of `Species` or anything else in
       `chemistry/` (confirmed by reading the file), aren't exported from
       `chemistry/__init__.py`, and their real consumers —
       `stream_adapter.py`'s `FeedState` and the stirred-tank template's
       default organism/substrate composition lookup
       (`OrganismConfig.resolve()`/`SubstrateConfig.resolve()`, confirmed
       load-bearing via `test_configs.py:258-259`) — already reached past
       `chemistry/` directly, the same way `PyOMES/__init__.py` does. Moved
       via `git mv` to sit alongside `units.py`/`config.py`/`stream_adapter.py`;
       updated all 9 import sites (`PyOMES/__init__.py`, `stream_adapter.py`,
       `templates/stirred_tank/{factory,configs}.py` ×4,
       `tests/run_tests.py`, `tests/standalone/{conftest,test_compounds,
       test_feed_state}.py`) plus doc references (`chemistry/__init__.py`
       docstring, `PyOMES/README.md` table, `OPEN_WORK.md`'s
       molar-mass-unification entry). Does not affect the D7
       `chemistry`<->`reactions` cycle discussion — `compounds.py` has no
       edge into either package. **Verification:** full standalone suite,
       **2064 passed**, 0 failed (covers checkpoint 18 above too; neither
       checkpoint added/removed a test)._
- [x] 20. Remove `ChemicalRegistry.IDs` and `FeedState`/`stream_adapter.py`
       (no consumers found anywhere in the repo).
       _Notes: done 2026-09-22/23 at `74dcba5` + `856e507` (landed as two
       commits — the first `git add` only staged the two file deletions and
       missed the other 8 edited files; caught by re-checking `git status`
       and `git show --stat` after push, finished in the follow-up commit).
       `ChemicalRegistry.IDs` existed only for bioSTEAM-shape compatibility
       (`self.chemicals.IDs`); a repo-wide search found exactly one caller,
       its own test. `FeedState` had no consumer anywhere outside its own
       test file and fixtures — not in `templates/`, `models/`, or any
       tutorial notebook; its one prior justification in
       `docs/dev/implementation/shipped/CUFERMENTER_SUNSET.md` (a
       `strong_ions.py` caller) no longer holds, since `strong_ions.py` was
       deleted in the `chemical-equilibrium-engines-subfolder` phase. Deleted
       `stream_adapter.py` and `test_feed_state.py` outright (`git rm`);
       removed the now-unused `simple_feed`/`rich_feed` conftest fixtures,
       the `.IDs` property and its test, the `FeedState` import/export in
       `PyOMES/__init__.py`, the dead `_simple_feed`/`_rich_feed` helpers in
       `tests/run_tests.py`, and stale doc references (`PyOMES/README.md`,
       root `README.md`, `docs/architecture.md`). Left
       `tests/run_tests.py`'s separate, pre-existing
       `create_standalone_fermenter` import failure alone — logged in
       `OPEN_WORK.md` since 2026-09-18, unrelated to this checkpoint, not
       part of this phase. **Verification:** full standalone suite,
       **2051 passed** (2064 − 13: the 12 `test_feed_state.py` tests plus
       `test_ids_property`), 0 failed. Full configured suite (`testpaths` =
       standalone + validation) re-confirmed 2026-09-23: **2086 passed**,
       0 failed, 166 warnings, 163.76s — this phase's final count._

## Shipping

- [x] Full test suite green on the branch — **2086 passed**, 0 failed,
      confirmed 2026-09-23 (checkpoint 20's number above; supersedes the
      "about 2070" estimate from checkpoint 17, same reconciliation as that
      checkpoint already did once for its own estimate)
- [ ] `git checkout main`
- [ ] `git merge --no-ff chemistry-reactions-kinetics-cleanup -m "Merge chemistry-reactions-kinetics-cleanup: <summary>"`
- [ ] `git tag chemistry-reactions-kinetics-cleanup-shipped <commit-hash>`
- [ ] `git push && git push --tags`
- [ ] `git branch -d chemistry-reactions-kinetics-cleanup` and
      `git push origin --delete chemistry-reactions-kinetics-cleanup`
- [ ] Move the plan doc and this checklist to `docs/dev/implementation/shipped/`; add
      "Shipped" banners
- [ ] Update `upcoming/README.md`'s "Recently shipped" list
