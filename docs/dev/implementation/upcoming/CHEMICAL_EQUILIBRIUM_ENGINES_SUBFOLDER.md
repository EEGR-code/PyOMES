# chemical_equilibrium engines/ subfolder — design discussion

> Status: pre-phase design discussion, 2026-09-20. No branch, no checklist,
> no code yet. Written up from a planning conversation that audited
> `PyOMES/chemical_equilibrium/` (16 flat `.py` files, ~5,970 lines). When this
> is picked up: follow `README.md`'s "How to start one" — write a checklist
> file from `PHASE_KICKOFF_TEMPLATE.md`, branch, implement, ship.
>
> **Scope note (2026-09-20):** widened the same day to include retiring the
> `activity_models.py` / `sit.py` compatibility layer in favour of `PyOMES.thermo`
> (see "Part A" below), done *before* the engine move so those imports are
> not rewritten twice. Also widened the same day to remove the orphaned
> `strong_ions.py` (see "Part B" below); the `engines/` move is now "Part C".
>
> **Scope note (2026-09-20, during checkpoint 7):** widened once more with
> "Part D" — unify every copy of the gas constant onto the one definition in
> `PyOMES/units.py`. Unlike Parts A–C, **Part D deliberately changes numbers**
> (the copies are not all the same value), so it runs last and is isolated in
> its own commits.

## Commit discipline for this phase

**The repo owner runs every `git add`, `git commit`, and `git push` for this
phase personally — an assistant executing this plan should not run those
commands.** At each checkpoint: make the file edits, report exactly what
changed (and run the checkpoint's sanity check), then stop and wait for the
owner to review, stage, and commit before moving to the next checkpoint.

## Goal

Declutter `PyOMES/chemical_equilibrium/` and make the engine boundaries
visible in the directory tree. Three parts, in this order:

- **Part A — retire the `thermo` compatibility layer.** `sit.py` and
  `activity_models.py` largely re-export models that now live in
  `PyOMES/thermo/`. Move the small amount of real code they hold into
  `thermo`, repoint importers at `PyOMES.thermo`, and delete both files.
- **Part B — remove `strong_ions.py`.** It has no production callers, holds
  a small hand-maintained salt-dissociation lookup users must trust, and
  belongs neither in this package nor (after
  `STRONG_ION_INFERENCE_GENERALIZATION.md`) in a form the engines need.
  Delete it and its test.
- **Part C — `engines/` subfolder.** Today the shared infrastructure and the
  three engine implementations sit side by side, and the Bisection engine is
  called just `engine.py`, so nothing in a file listing says which files
  belong to which engine. Move the engine-specific files under `engines/`.

- **Part D — one source of truth for the gas constant.** `PyOMES/units.py`
  already defines it (and says it is the authoritative place), but about ten
  other definitions and hard-coded literals remain, at three different
  precisions. Repoint them all at `units.py`, derive the L·atm value from the
  J value so there is exactly one numeric literal, and add a guard test so
  new copies fail CI.

**Parts A–C are pure file moves, import rewrites and dead-code removal — no
behaviour change.** Part D is the one deliberate numerical change (see its
audit and Decisions 11–15 below).

## Audit: what is shared vs. engine-specific

| Group | Files |
|---|---|
| Shared (stay at top level) | `protocols.py`, `activity.py`, `activity_dispatch.py`, `numerical_gradient.py` |
| Thermo compatibility layer (Part A: dissolved into `thermo/`) | `activity_models.py`, `sit.py` |
| Orphaned helper (Part B: deleted) | `strong_ions.py` |
| Orphaned Bisection-only entry points (checkpoint 9b: deleted) | `api.py`, `factory.py` |
| Bisection engine | `engine.py`, `acid_base.py` (~2,040 lines) |
| NR engine | `nr_engine.py`, `nr_tableau.py`, `nr_solver.py` (~2,520 lines) |
| PHREEQC engine | `phreeqc_engine.py` (327 lines) |

The engines' internals do not import each other. The engine-specific files
depend only on `protocols` and `activity_models`, so the split follows an
existing seam.

Notes:

- `acid_base.py` calls itself the "Core acid-base speciation solver (Level 1
  foundation)", but its only consumers are `engine.py` and the package
  `__init__` re-export of `solve_acid_base`. It is Bisection's back end.
- `numerical_gradient.py` wraps *any* engine and imports only `protocols`, so
  it stays shared.
- `api.py` / `factory.py` were originally planned to stay at the top level as
  "the package's public entry points". A re-check during Part C (2026-09-20)
  found they are hard-wired to the Bisection engine (`SpeciationFactory`
  builds only that engine; `SpeciationEngineAdapter` is typed to it and calls
  its legacy keyword interface) and have **no callers anywhere in the repo**
  (no test, notebook, model or tutorial). Engine selection already lives in
  `ReactionSystem.engine`. They are deleted in checkpoint 9b (Decision 16).

### Thermo compatibility layer (Part A)

- **`sit.py` is a pure shim**: it re-exports `SITLiquidModel` and the legacy
  `SITActivityModel` alias from `thermo`. Nothing in the package imports it;
  the only importers are two tests in `test_liquid_phase_model.py` that test
  the alias itself.
- **`activity_models.py` is a mix:**

  | Content | Kind | Users |
  |---|---|---|
  | Water-property re-exports (`debye_huckel_A`, etc.) | shim | `activity.py` only |
  | `IdealActivityModel` / `DaviesActivityModel` aliases | shim | package `__init__.py`; `DaviesActivityModel` in 3 validation notebooks, `tests/validation/speciation/_generate_notebooks.py`, `test_saturation_index.py`, `test_liquid_phase_model.py` |
  | `ActivityModel` protocol (per-ion `gamma(z, I, T_K)`) | **real code** | `acid_base.py` type hint |
  | `make_activity_model(use_activity, name)` factory | **real code** | `engine.py`, `nr_engine.py`, `acid_base.py`, `test_nr_tableau_gas_liquid.py`, `OPEN_WORK.md` |

  The real code cannot simply be deleted: `thermo` has no string-to-model
  factory (`ThermoFramework` takes model instances), and the engines' `use_activity`
  / `activity_model="davies"` kwargs need one.
- **Dead code in `activity.py`:** `davies_log10_gamma` and `davies_gamma` have
  no callers anywhere in the repo and duplicate `DaviesLiquidModel`. Their
  module docstring ("scaffolded for future use") is stale.

### `strong_ions.py` (Part B)

`strong_ions_from_feed_molL(feed)` infers `CT_*` strong-ion totals from a
BioSTEAM stream, a `FeedState`, or a plain dict, expanding neutral salts via
`SALT_DISSOCIATION_MAP` in `chemistry/registry.py`. Findings (repo-wide search
of `.py` and `.md` files; outside users cannot be seen):

- **No production callers.** Only the package `__init__.py` re-export, its
  own test (`test_strong_ions.py`), and a comment in `registry.py` reference
  it. The engines take `strong_ions=` as a plain dict and never call it.
- **A one-consumer lookup.** `SALT_DISSOCIATION_MAP` (exported from
  `PyOMES.chemistry`) has this function as its only user. `chem_recipe.py`'s
  `ChemSpec` registry holds the same salt-to-ion information again for about
  30 compounds — two parallel tables for one job.
- **Silent-drop behaviour.** The output is a fixed set of 11 `CT_*` keys,
  narrower than the NR engine's own allowlist (no Fe, Cu, MoO4, `S_cat` /
  `S_an`). An unrecognised salt contributes nothing, without a warning.
  `_mol_L_from_object` returns `0.0` on any exception, so a bad stream or a
  missing `F_vol` silently drops ions. Bare ids such as `"K"` / `"Na"` are
  treated as the ion, which is ambiguous against BioSTEAM element-style ids.
  This is the same silent-drop pattern that `STRONG_ION_INFERENCE_GENERALIZATION.md`
  documents.
- **Stale docs.** The docstring cites `fermenter.stream_adapter.FeedState`;
  `FeedState` now lives in `PyOMES/stream_adapter.py`.

Rationale for deleting rather than relocating: a model knows which ions it
is charge-balancing and can state them (`strong_ions=` today, or species in
`n_mol` once Phase 1 of the strong-ion note lands). If BioSTEAM coupling later
needs neutral-salt expansion, it belongs in `PyOMES/stream_adapter.py`, emitting
species ids rather than `CT_*` keys and warning on unmapped species — not in
the equilibrium package.

### Gas-constant definitions (Part D)

Audited 2026-09-20 (`.py` and `.ipynb`; `scratch/` excluded). **The root
exists:** `PyOMES/units.py` defines `R_J_PER_MOL_K = 8.31446261815324` (CODATA
2018) and `R_L_ATM_PER_MOL_K = 0.08205736608096`, and its docstring says it
exists because "the codebase historically contained repeated definitions of
common constants (e.g., the ideal gas constant in L·atm/(mol·K))". A first
consolidation of the J value was done 2026-07-02
(`equilibrium-constraint-unification` CP1, BSM2 sentinels re-baselined; see the
comment in `test_bsm2_reference.py`); it missed the copies below and never
touched the L·atm value.

| Definition | Value | vs root | Notes |
|---|---|---|---|
| `core/phases.py:28` `R_L_ATM_MOL_K` | 0.0820574 | +4.1e-7 | Different *name* and value from the root. The most-imported copy: 25 files (7 notebooks, both generators, ~10 tests, `nr_solver`, `boundaries`, `solvers`, `gas_liquid_link`, `templates/stirred_tank/factory`) |
| `chemistry/partition.py:51` `_R_L_ATM_MOL_K` | 0.0820574 | +4.1e-7 | Comment says "must match core.phases" — the duplication is acknowledged in the code |
| `equilibria/peng_robinson.py:38` `R` | 0.0820574 | +4.1e-7 | |
| `control/cv_loops.py:912` `_R_L_ATM_PER_MOL_K` | 0.08205736608095958 | +5e-15 | A private copy of the precise value |
| `control/cv_loops.py:1129` `_R_UNIV` | 8.31446261815324 | identical | Bit-identical private copy |
| `models/vlmodels/adm1/{base:224, bsm2:328, bsm2_direct:272}` `_R_J` | 8.31446 | −3.2e-7 | Van 't Hoff `Ka(T)`; the same rounding the 2026-07-02 re-baseline removed from `thermo_params.py` / `framework.py` |
| `reactions/plots.py:18` `_R_GAS` | 8.314 | −5.6e-5 | Van 't Hoff plot only |
| Literals in tests | 0.0820574 (`test_gas_liquid_link:895`, `test_nr_gas_liquid_cp2` ×2, `test_partition_model` ×2, `test_precipitation_gas_liquid_cp5:29`); 8.31446261815324 (`test_equilibrium_constants:22`, added in checkpoint 7 as a frozen reference) | | |
| `docs/tutorials/D2C_workshop/Example1_mtp_well.ipynb:202` `R_LA` | 0.08205 | −9e-5 | |

Already correct (import from `units`): `thermo/framework.py`,
`thermo/equilibrium_constants.py`, `chemistry/thermo_params.py`,
`chemistry/partition.py` (J only), `reactions/equilibrium.py`,
`equilibria/vle.py`, `models/vlmodels/headspace.py`.

Two findings matter beyond tidiness:

- **Production currently uses two values of R for the same physics.**
  `equilibria/vle.py` and `models/vlmodels/headspace.py` use the precise
  0.08205736608096; `core/phases.py` (partial pressure, gas moles, boundaries,
  solvers) and `partition.py` use 0.0820574, 4.1e-7 higher. Every gas-liquid
  calculation therefore mixes them.
- **The root is not itself single-sourced.** `R_L_ATM_PER_MOL_K` is a second
  literal, not derived from `R_J_PER_MOL_K`. Three "precise" versions exist
  (root `0.08205736608096`, `R_J / 101.325` = `0.08205736608095968`, and the
  `cv_loops` copy `0.08205736608095958`); they differ by up to 4e-15.

Other fundamental-looking constants (`273.15`, `101325`, `9.81`, `96485`,
Boltzmann, Avogadro) were not audited in depth; a quick count found 28
literal occurrences of them across 16 package files. Out of scope for this phase (only
`R` is unified), but the Part D guard test is the natural place to extend later.

## Target layout

```
chemical_equilibrium/
    __init__.py            (public re-exports, unchanged)
    protocols.py
    activity.py  activity_dispatch.py
    numerical_gradient.py
    engines/
        __init__.py
        bisection/   __init__.py, engine.py, acid_base.py
        nr/          __init__.py, engine.py, tableau.py, solver.py
        phreeqc.py

thermo/
    liquid_phase_model.py  (gains: ActivityModel protocol, make_activity_model —
                            or a small thermo/factory.py; decide at kickoff)
```

The top level drops from 16 files to 5 (including `__init__.py`), all of them
engine-agnostic. `activity_models.py`, `sit.py`, `strong_ions.py`, `api.py` and
`factory.py` are gone. `nr_engine.py` becomes
`engines/nr/engine.py`, so the `nr_` prefixes become redundant. PHREEQC is a
flat file because it is a single 327-line module; it can become a package if it
grows. The optional `phreeqpython` dependency stays confined to one file.

## Decisions

1. **Layout:** the full split above (chosen over renaming in place, or an
   `engines/` folder holding only the three engine classes).
2. **No back-compat shims for the engine modules.** The old import paths
   (`chemical_equilibrium.nr_engine`, `.engine`, `.acid_base`, ...) stop
   working when the files move; every in-repo caller is updated in the phase.
   This is a deliberate breaking change for any outside code importing those
   paths (see Risks).
3. **`_vant_hoff` duplication:** resolve it *before* the move, so the shared
   helper does not end up copied into separate `engines/` subfolders.
   `acid_base._vant_hoff_K` and `nr_tableau._vant_hoff_log_K` overlap, and
   `notes/review-notes.md` already flags this. **The single copy lives in
   `PyOMES/thermo/`** (a small module such as `thermo/equilibrium_constants.py`;
   exact filename at kickoff): temperature-correcting an equilibrium constant is
   thermodynamics, and `chemical_equilibrium` already depends on `thermo`.
4. **Part A goes first.** Retiring `activity_models.py` / `sit.py` removes
   imports from `engine.py`, `nr_engine.py` and `acid_base.py` before those
   files move, so they are not rewritten twice.
5. **`ActivityModel` / `make_activity_model` move into `thermo/`.** The
   dependency direction stays correct (`chemical_equilibrium` already depends
   on `thermo`, not the reverse). `liquid_phase_model.py` already documents
   `ActivityModel` as its sibling protocol, with a stale `PyOMES/speciation/`
   path that this fixes.
6. **Dead Davies helpers are deleted** from `activity.py`, not migrated.
7. **`strong_ions.py` is deleted, not relocated** (Part B), and its test with
   it. It has no production callers; see the audit above.
8. **`SALT_DISSOCIATION_MAP` and `chem_recipe.py` are left in place.** Once
   `strong_ions.py` is gone the map has no consumer, but whether to keep,
   merge or remove the salt tables is the recipe-layer question that
   `STRONG_ION_INFERENCE_GENERALIZATION.md` already defers to a separate note.
9. **`IdealActivityModel` / `DaviesActivityModel` are dropped from the package
   `__init__.py`** and `__all__`, with no alias or deprecation period. Callers
   use `IdealLiquidModel` / `DaviesLiquidModel` from `PyOMES.thermo`.
10. **`strong_ions_from_feed_molL` is dropped from the package `__init__.py`**
    and `__all__` along with the module, with no stub or deprecation period.
    Subject to the `.ipynb` search in checkpoint 5.
11. **`PyOMES/units.py` stays the single source of truth for `R`** (Part D).
    It already documents itself as that, and 8 modules import from it, so no new
    `constants.py`. Fix its stale `fermenter.units` docstring name. *Proposed —
    confirm at kickoff of Part D.*
12. **Exactly one numeric literal for `R`.** `R_J_PER_MOL_K` stays the literal;
    `R_L_ATM_PER_MOL_K` becomes `R_J_PER_MOL_K / (PA_PER_ATM / L_PER_M3)`
    (`0.08205736608095968`, 3.9e-15 from today's root value), so the two can
    never disagree. A unit test pins the relation. *Proposed.*
13. **Rename, do not alias, `core.phases.R_L_ATM_MOL_K`.** Importers switch to
    `PyOMES.units.R_L_ATM_PER_MOL_K` (25 files, including 7 notebooks and both
    generators). This follows Decision 2 (no shims), and the two names
    (`..._MOL_K` vs `..._PER_MOL_K`) are themselves part of how the split
    happened. The alternative — keep `core.phases.R_L_ATM_MOL_K` as a
    re-export of the `units` value — is single-valued and far smaller, but
    keeps two names. *Proposed: rename; say so at kickoff if you prefer the
    smaller alias.*
14. **The ADM1 / BSM2 `_R_J = 8.31446` copies are unified too**, with the BSM2
    sentinels re-baselined and the before/after recorded, exactly as the
    2026-07-02 consolidation did for `thermo_params.py` / `framework.py`.
    **Needs your call:** these are benchmark-reference implementations, so if the
    rounded value is deliberate (to match a published ADM1/BSM2 specification),
    keep it and replace the bare literal with a commented, named constant instead.
15. **A guard test fails CI on any new `R` literal** outside `units.py`
    (a source scan of `PyOMES/` and `models/` for the gas-constant literal
    patterns, with an explicit allowlist). The allowlist is added in
    checkpoint 13 listing the known remaining copies, and emptied in
    checkpoint 14.
16. **`api.py` and `factory.py` are deleted, not moved or generalised**
    (checkpoint 9b). Neither has any caller in the repo, the package has no
    outside users yet (so the exported names `SpeciationFactory` and
    `SpeciationEngineAdapter` carry little weight), and both are
    Bisection-only: the adapter's typed API (`CT_TIC`/`CT_NH_T`/`CT_P` totals,
    a pH scan with `n_scan`/`pH_min`/`tol`) does not fit NR (totals by master
    species plus `strong_ions`, built from declared reactions) or PHREEQC (a
    `component_map`). A general factory would also duplicate
    `ReactionSystem.engine`. Restorable from commit `2e5554a`
    (`git show 2e5554a:PyOMES/chemical_equilibrium/api.py`, and likewise
    `factory.py`). Consequences outside this phase's scope are logged in the
    checklist and in `STRONG_ION_INFERENCE_GENERALIZATION.md`'s
    recipe-layer section.

## Proposed checkpoints

1. Inventory every import site (`.py` and `.ipynb`) of the modules being moved
   or deleted. Record counts in the checklist. Baseline: full test suite green.

**Part A — thermo consolidation**

2. Add `ActivityModel` and `make_activity_model` to `thermo/`; export from
   `PyOMES.thermo`. Update the stale `PyOMES/speciation/` docstring paths in
   `liquid_phase_model.py`. `activity_models.py` still exists at this point.
3. Repoint importers at `PyOMES.thermo`: `engine.py`, `nr_engine.py`,
   `acid_base.py`, `activity.py` (`thermo.water_properties`),
   `test_nr_tableau_gas_liquid.py`, and the `DaviesActivityModel` users in
   `tests/validation/speciation/` (tests, generator, 3 notebooks —
   `DaviesLiquidModel`). Remove `IdealActivityModel` / `DaviesActivityModel`
   from `__init__.py` and `__all__` (Decision 9).
4. Delete `sit.py` and `activity_models.py`; delete the two alias tests in
   `test_liquid_phase_model.py`; repoint the factory test there. Delete
   `davies_log10_gamma` / `davies_gamma` and fix `activity.py`'s docstring.
   Run `test_liquid_phase_model.py`, `test_thermo_framework.py`,
   `test_speciation.py`, `tests/validation/speciation/`.

**Part B — remove `strong_ions.py`**

5. Search `.ipynb` files too (the earlier audit covered `.py` and `.md`) for
   `strong_ions_from_feed_molL` / `PyOMES.chemical_equilibrium.strong_ions`.
   If clear, delete `strong_ions.py` and `tests/standalone/test_strong_ions.py`;
   remove the `__init__.py` import and `__all__` entry (Decision 10);
   drop the stale `strong_ions_from_feed_molL` comment in
   `chemistry/registry.py`. `FeedState` and its fixtures in `conftest.py` /
   `test_feed_state.py` stay.
6. Update the "Adjacent, out of scope" section of
   `STRONG_ION_INFERENCE_GENERALIZATION.md`: `strong_ions.py` was removed and
   `SALT_DISSOCIATION_MAP` now has no consumer. Run `test_feed_state.py` and the
   speciation tests.

**Part C — engines/ subfolder**

7. De-duplicate `_vant_hoff_K` / `_vant_hoff_log_K` into a single helper in
   `PyOMES/thermo/` (Decision 3); own commit, own tests. Repoint `acid_base.py`
   and `nr_tableau.py` at it.
8. Move the NR files (`engines/nr/`); rewrite relative imports. A file in
   `engines/nr/` (or `engines/bisection/`) is two packages deeper than
   before, so `..units` becomes `....units` and `..core.phases` becomes
   `....core.phases` (four dots); a sibling-of-engines import such as
   `.protocols` becomes `...protocols`. (The flat `engines/phreeqc.py` is
   only one deeper: `..units` becomes `...units`.) Run `test_nr_*` and
   `test_equilibrium_classification.py`.
9. Move the Bisection files (`engines/bisection/`). Run `test_speciation*.py`
   and `test_bisection_chemical_equilibrium_engine_alias.py`.
9b. *Added during Part C.* Delete `api.py` and `factory.py` and their two
    package-level exports (Decision 16); reword the two `protocols.py`
    docstrings that cite them. Log the knock-on effects outside this phase's
    scope (checklist and `STRONG_ION_INFERENCE_GENERALIZATION.md`). Bit-identical
    apart from the removed names.
10. Move `phreeqc_engine.py` to `engines/phreeqc.py`. Run
   `tests/validation/speciation/`.
10b. *Added during Part C; unrelated to the layout.* Surfaced while
    verifying checkpoint 10: the two notebooks that guard the PHREEQC import
    used different flag names (`_HAVE_PHREEQC` in the ArXiv tutorial,
    `HAS_PHREEQC` in validation notebook 10), and the ArXiv guard did not
    actually detect a missing `phreeqpython` (the engine imports it lazily, so
    importing the module never fails). Standardise on `HAS_PHREEQC` and add the
    presence check to the ArXiv guard, in the generator template and the
    notebook together.
11. Update `__init__.py` re-exports, then external callers: `PyOMES/reactions/
   reaction_system.py`, `models/vlmodels/adm1/{base,bsm2}.py`, the notebook
   generators, and tests. No shims are left at the old paths (Decision 2), so
   finish with a repo-wide search (`.py`, `.ipynb`, `.md`) for each old module
   path to confirm nothing still points at it.
12. Update docstring cross-references (about 15 files use
    `PyOMES.chemical_equilibrium.engine...`-style paths) and any current docs,
    including `OPEN_WORK.md` (cites `activity_models.make_activity_model`) and
    `docs/architecture.md` (lists `activity_models.py`, `sit.py`). Leave
    historical docs under `docs/dev/implementation/shipped/` and
    `docs/dev/ideas/` alone — they record what was true when written.

**Part D — gas-constant unification** (runs after Part C so Parts A–C stay
bit-identical and verifiable; Part D changes numbers)

13. *Value-preserving.* Re-verify the audit table above with a fresh search.
    Replace the one bit-identical private copy (`cv_loops._R_UNIV`) with an
    import from `units`, and make `test_equilibrium_constants.py` import `R`
    instead of holding a literal. Add the guard test (Decision 15) with an
    allowlist of every copy that remains. Results must be bit-identical
    (fingerprint, as in checkpoint 7).
14. *Numerics-changing, own commit.* Derive `R_L_ATM_PER_MOL_K` from
    `R_J_PER_MOL_K` (Decision 12); repoint `core/phases.py` (rename per
    Decision 13, all 25 files), `partition.py`, `peng_robinson.py`,
    `cv_loops.py`, `plots.py`, and — per Decision 14 — the three ADM1/BSM2
    files; replace the test literals and the `Example1_mtp_well.ipynb`
    `R_LA` with imports; empty the guard allowlist. Expect: every
    `p = nRT/V` shifts by about 4e-7 relative, and the BSM2 sentinels
    (`RTOL_SENTINEL = 1e-9`) and any test with a tight tolerance on a gas
    quantity will move. Before changing anything, record a fingerprint of
    gas-liquid results (partial pressures, gas moles, BSM2 final state); after,
    record the shift and re-baseline with a dated before/after comment in
    `test_bsm2_reference.py`, as the 2026-07-02 note does. The shift should
    match the 4.1e-7 (or 3.2e-7 for ADM1 `Ka(T)`) prediction and nothing else.
15. Full suite green, then ship per the kickoff template.

Each move checkpoint should be a `git mv` plus import rewrites only, so history
follows the files.

## Risks

- **Wide import churn.** Tests, notebook generators, `models/vlmodels/adm1/*`
  and `reactions/reaction_system.py` all import the old module paths directly.
  Missed sites fail at import time, so the test suite catches them, but
  notebooks (`.ipynb`) and lazy in-function imports (e.g.
  `reaction_system.py`, `nr_tableau.py:435`) only fail when executed.
- **`ArXiv_preprint` and `tests/validation/speciation` generators** embed import
  paths inside string templates, so a plain search-and-replace of `import`
  lines can miss them. If `NOTEBOOK_GENERATOR_REMOVAL.md` ships first, this
  shrinks; sequence the two phases accordingly.
- **Deliberate breaking changes, no deprecation period.** Decisions 2, 9 and
  10 remove import paths and package exports outright: the old engine module
  paths, `IdealActivityModel` / `DaviesActivityModel`, and
  `strong_ions_from_feed_molL`. Any outside code using them breaks at import
  time. Acceptable for this fork, but note it in release notes or the
  changelog if one is kept.
- **Saved notebooks.** The three `tests/validation/speciation/` notebooks have
  outputs embedded; only the import line changes, but they are still edited
  files and need re-running to confirm.
- **Phase-doc paths.** Existing docs already cite the stale `src/speciation/`
  layout; do not "fix" those as part of this phase.

## Open questions

All resolved 2026-09-20; the answers are recorded in Decisions above. Kept
here as a record.

1. **Shims at the old paths?** *Resolved: no shims* (Decision 2).
2. **`IdealActivityModel` / `DaviesActivityModel` in the package `__init__.py`?**
   *Resolved: drop outright* (Decision 9).
3. **Where does the shared `_vant_hoff` helper live?** *Resolved: `PyOMES/thermo/`*
   (Decision 3).
4. **Should `api.py` / `factory.py` move under `engines/bisection/`?**
   *Originally resolved (default, not put to a vote): no*, on the grounds that
   they were public entry points and the natural place for engine selection.
   *Reopened and re-resolved 2026-09-20 (Decision 16): delete them.* They have
   no callers, the package has no outside users yet, and a general factory
   would duplicate `ReactionSystem.engine`.
5. **Interaction with `NOTEBOOK_GENERATOR_REMOVAL.md`.** *Resolved (default):*
   Part A edits `tests/validation/speciation/_generate_notebooks.py` and its
   notebooks. If the generator has been retired first, only the notebooks need
   editing; otherwise edit both in the same checkpoint so they do not drift.
6. **Interaction with `STRONG_ION_INFERENCE_GENERALIZATION.md`.** *Resolved
   (default):* that note edits `nr_engine.py`, `nr_solver.py`, `engine.py` and
   `acid_base.py`. Land this phase first if both are imminent, since it is the
   smaller logical change; whichever lands second rebases onto the moved paths.
7. **`strong_ions_from_feed_molL` in the package `__init__.py`?** *Resolved:
   drop outright* (Decision 10), subject to the `.ipynb` search in checkpoint 5.
