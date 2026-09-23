# Phase Kickoff Checklist — equilibrium-set-relocation

> Checklist for [`EQUILIBRIUM_SET_RELOCATION.md`](EQUILIBRIUM_SET_RELOCATION.md),
> the source of truth for motivation and design; do not restate it here. Where
> this checklist and the note disagree, this checklist wins: the note was written
> from a conversational audit and the "Re-verification" section below corrects
> it. See [`README.md`](README.md)'s "Branching and tagging convention". Modelled
> on
> [`../shipped/CHEMISTRY_REACTIONS_KINETICS_CLEANUP_CHECKLIST.md`](../shipped/CHEMISTRY_REACTIONS_KINETICS_CLEANUP_CHECKLIST.md).

**Working rules**

- Checkpoint 1 is a pure move: no behaviour change. Checkpoint 2 deletes one
  method and rewrites two tests; nothing else changes behaviour.
- The repo owner runs every `git` command (branch, add, commit, push). At each
  checkpoint: edit, run its sanity check, report, stop, and hand over the
  commands.
- No shims, aliases or re-exports (decision 1 below).
- Move the file with a plain filesystem move (no git commands from Claude). Leave `docs/dev/implementation/shipped/` untouched,
  and do not edit `docs/dev/ideas/` or the historical comment in
  `tests/standalone/test_bsm2_reference.py:190`; they describe past states.
- Docs and docstrings describe current behaviour only: no phase or checkpoint
  labels, no pointers to this checklist or the design note.
- Anything that looks like a bug or dead code beyond scope is noted, not fixed.

**Decisions** (settled 2026-09-23, before starting)

1. No re-export of `EquilibriumSet`/`EquilibriumDef` from
   `chemical_equilibrium/__init__.py`; the deep path stays the only path.
2. The `EquilibriumSet` import in `BisectionChemicalEquilibriumEngine.from_reactions()`
   stays lazy, beside the `reactions.equilibrium` import it already shares a
   block with.
3. The move and the `bsm2_default()` deletion ship together as two separate
   checkpoints in this one checklist.

## Pre-flight

- [x] `git status -sb` clean, `main` level with `origin/main` (checked
      2026-09-23)
- [x] Full-suite baseline on `main`: **2086 passed**, 0 failed, 166 warnings,
      4m19s (`python -m pytest -p no:cacheprovider`)
- [ ] Branch created off `main`: `equilibrium-set-relocation`
- [ ] This checklist committed on that branch as the first commit

## During

- [x] Plan doc exists in `docs/dev/implementation/upcoming/`
      ([`EQUILIBRIUM_SET_RELOCATION.md`](EQUILIBRIUM_SET_RELOCATION.md))
- [ ] Checkpoints tracked below, one commit each
- [ ] **If work stalls:** add a status banner to the top of the plan doc at once

### Re-verification of the design note (2026-09-23, on `main`)

Searched every tracked file type (`.py`, `.ipynb`, `.md`, `.toml`, `.yml`,
`.json`), plus the gitignored `notes/`, `scratch/` and `PyOMES.egg-info/` with
plain `grep`, for `EquilibriumSet`, `EquilibriumDef`, `WaterDef`,
`bsm2_default` and `chemistry.equilibria`/`chemistry/equilibria`. `.github/`
has no hit. `notes/` and `PyOMES.egg-info/` are gitignored and untracked, so
nothing to edit there. The install is editable, so the move needs no reinstall.

**Confirmed as the note says**

- `engine.py:202` is the lazy import; its docstring refs are at lines 18, 92,
  156. `acid_base.py:solve_from_equilibrium_set` takes the set duck-typed, with
  no import.
- The only import sites of `chemistry.equilibria` are `engine.py:202` and
  `tests/standalone/test_cv_compute_interface.py:278,294`.
- `bsm2_default()` is *invoked* only by its own module docstring
  (`equilibria.py:23`) and `test_cv_compute_interface.py:279,295`. `bsm2.py`
  never calls it: its `_BSM2_*` constants are independent literals and the
  "matches bit-for-bit" mentions (`bsm2.py:472,734`) are comments only.
- `chemistry/` has zero references to `chemical_equilibrium`. Every
  `chemical_equilibrium` reference elsewhere in `PyOMES/` is a docstring or an
  in-function import (`reaction_system.py:261,273`, `core/solvers.py:653`), so
  the new edge cannot close a cycle. `tests/standalone/test_import_graph_acyclic.py`
  confirms it after the move.
- `chemistry/__init__.py` does not export the names; its docstring clause is
  "and acid-base equilibrium sets".

**Discrepancies (all folded into the checkpoints below)**

- "Exactly three files" is wrong. Additional references to `bsm2_default()` or
  the class names that the deletion would leave false:
  `acid_base.py:1083,1105,1130,1145` (comments and a docstring),
  `test_speciation.py:491-494` (class docstring), `bsm2.py:472,734`
  (comments), `OPEN_WORK.md:273-277,374`. The
  `01_bisection_engine_basics.ipynb` notebook names the classes in prose only
  (via `engine._equilibrium_set`); no import path, no change needed.
- `OPEN_WORK.md:273-277` says the BSM2 model still uses the deprecated
  `{name}_HA`/`{name}_A-` fallback through its four VFA rows. Checked by
  building the engine from `_build_bsm2_equilibrium_reactions()`: every entry,
  VFAs included, carries `species_refs` (`S_ac` -> `['S_ac', 'S_ac-']`, and so
  on). The fallback in `_compute_species_eq` is reachable only through a
  hand-built `EquilibriumSet` with no `species_refs`. The code stays; the
  wording is corrected, and the dead path is logged in `OPEN_WORK.md`.
- After the deletion, the `common_species` import in `equilibria.py` (`CO2`,
  `HCO3_minus`, `CO3_2minus`, `NH4_plus`, `NH3`) is unused: `bsm2_default()` is
  its only user. Checkpoint 1 rewrites it as the note says; checkpoint 2
  removes it, so the moved file ends with no import from `chemistry/` and the
  `chemical_equilibrium` -> `chemistry` edge is transient.
- Stale line numbers: `OPEN_WORK.md:273` cites `equilibria.py:510-511` (the
  file has 499 lines); `architecture.md` lists `bisection/` at 387-388, not
  386. Its `chemistry/` listing also names files that no longer exist
  (`recipe.py`, `registry.py`, `types.py`, `chem_recipe.py`, `compounds.py`);
  only the `equilibria.py` entry is in scope here.
- Unrelated, left alone: `upcoming/README.md:9` has a stray `+-*` line.

### Checkpoints

- [x] 1. Move the module (a plain filesystem move, since the owner runs git;
      `git add` records it as a rename)
      `PyOMES/chemistry/equilibria.py` ->
      `PyOMES/chemical_equilibrium/engines/bisection/equilibria.py`; rewrite its
      two imports to `from ....units import R_J_PER_MOL_K as _R_J` and
      `from ....chemistry.common_species import CO2, HCO3_minus, CO3_2minus,
      NH4_plus, NH3`. `engine.py`: lazy import becomes
      `from .equilibria import EquilibriumSet`; repoint the three docstring refs
      (lines 18, 92, 156) to
      `PyOMES.chemical_equilibrium.engines.bisection.equilibria.EquilibriumSet`.
      `test_cv_compute_interface.py:278,294`: update the import path.
      `chemistry/__init__.py`: drop "and acid-base equilibrium sets" from the
      module docstring. `docs/architecture.md`: remove `equilibria.py` from the
      `chemistry/` listing and add it to the `bisection/` entry. Nothing
      else moves; `bsm2_default()` is still present.
      Sanity: `python -c "import PyOMES"`; import-graph guard green; both
      `tests/standalone/test_cv_compute_interface.py` and
      `test_bsm2_reference.py` pass; full suite **2086 passed**; repo-wide
      search finds no remaining `chemistry.equilibria` outside shipped docs,
      `docs/dev/ideas/`, the design note and `OPEN_WORK.md`.
      _Notes: done 2026-09-23. Sanity: `import PyOMES` and the new deep import
      work; the four files `test_import_graph_acyclic.py`,
      `test_cv_compute_interface.py`, `test_bsm2_reference.py` and
      `test_speciation.py` give 75 passed; full suite **2086 passed**, 0 failed,
      166 warnings, 3m04s (unchanged from the baseline). No `.py`, `.ipynb`,
      `.toml`, `.yml` or `.json` file outside `docs/dev/` still names
      `chemistry.equilibria`. The moved file keeps the `....` relative-import
      style its neighbours (`engine.py`, `acid_base.py`) already use, even
      though the package also uses absolute `from PyOMES...` imports elsewhere
      (`control/`, `templates/`, `chemistry/partition.py`); the mixed style is
      logged in checkpoint 2's `OPEN_WORK.md` edits, not changed here. A stale
      `chemistry/__pycache__/equilibria.cpython-312.pyc` is left behind
      (gitignored, and never imported without its source)._
- [ ] 2. Delete `EquilibriumSet.bsm2_default()`.
      `equilibria.py`: remove the method and its "Factory presets" header, the
      "Load a preset and modify (approach 2)" block from the module docstring
      (relabel the remaining example), and the now-unused `common_species`
      import. `test_cv_compute_interface.py`: replace the two
      `bsm2_default()` call sites with a small synthetic `EquilibriumSet`
      (one helper beside `_make_speciation_reactions()`; entries with
      `species_refs` for a carbonate and an ammonium ladder, entries without for
      two acids), drop BSM2 branding from the test names and docstrings, and
      make the first test also assert named species so it cannot pass
      vacuously. Reword the comments and docstrings that name `bsm2_default()`:
      `acid_base.py:1083,1105,1130,1145`, `test_speciation.py:491-494`,
      `bsm2.py:472,734` (keep the statement that the constants are Rosen &
      Jeppsson 2006 values; drop the claim of matching a preset that no longer
      exists). `OPEN_WORK.md`: correct the bullets at ~273-277 and ~374 (fallback
      is not reached by BSM2 and there is no preset), add a short entry for the
      now-unreachable-from-repo-code fallback, add one for the mixed
      relative/absolute import style (15 lines of 3+ dot imports, all in
      `engines/bisection/` and `engines/nr/`; absolute `from PyOMES...` used in
      `control/`, `templates/`, `chemistry/partition.py`), close out the `EquilibriumSet`
      location bullet (~681-684) and fix the "Three small loose ends" heading
      count.
      Sanity: repo-wide search for `bsm2_default` finds only historical docs
      (shipped, `docs/dev/ideas/`, `test_bsm2_reference.py:190`, this
      checklist and the design note); both named test files pass; full suite
      **2086 passed** (two tests replaced by two); `equilibria.py` ~465 lines
      with no import from `chemistry/`.

## Shipping

- [ ] Full test suite green on the branch (expected 2086 passed)
- [ ] `git checkout main`
- [ ] `git merge --no-ff equilibrium-set-relocation -m "Merge equilibrium-set-relocation: <summary>"`
- [ ] `git tag equilibrium-set-relocation-shipped <commit-hash>`
- [ ] `git push && git push --tags`
- [ ] `git branch -d equilibrium-set-relocation` and
      `git push origin --delete equilibrium-set-relocation`
- [ ] Move the design note and this checklist to `docs/dev/implementation/shipped/`;
      add "Shipped" banners
- [ ] Update `upcoming/README.md`: remove the "Design discussions" entry, add a
      "Recently shipped" entry
