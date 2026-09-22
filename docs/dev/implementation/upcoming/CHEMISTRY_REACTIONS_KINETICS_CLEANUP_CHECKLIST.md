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
      own dual imports). Deleted `TestHenryPartitionAlias`/
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
