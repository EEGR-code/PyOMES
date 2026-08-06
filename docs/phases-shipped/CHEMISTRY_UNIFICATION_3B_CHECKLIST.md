# Chemistry Unification 3b — Checklist

> **Status: Shipped 2026-06-03** — all 12 active checkpoints landed
> across 9 commits on the `chemistry-unification-3b` branch. 945
> standalone tests passing (913 baseline + 32 net new). BSM2 reference
> sentinels re-baselined twice: once for CO₂ key rename (C5, fixing the
> ghost-key CO₂ transfer bug) and once for `_refresh_derived` writeback
> list update (C7). Tag: `chemistry-unification-3b-shipped`.

Working checklist for `chemistry-unification-3b`. Conceptual framing
is in
[CHEMISTRY_UNIFICATION.md](CHEMISTRY_UNIFICATION.md) under "Phase E"
and "Leakage That Won't Go Away Automatically → Leak 2". The branch
scope and resolved design knobs are in
[CHEMISTRY_UNIFICATION_PLAN.md](CHEMISTRY_UNIFICATION_PLAN.md) under
"Phase 3b mandate". This file is the implementation plan: file-level
edits, ordered checkpoints, and sanity checks.

---

## Goal

Three tightly coupled pieces:

1. **Unified Species-ID emission.** The speciation engine emits
   species using their declared `Species.id` directly — no
   `_CANONICAL_NAMES` table, no synthesised `CO2aq`/`_HA`/`_A-`
   decorators for properly-declared species. CO₂ dissolution adopts
   phase-agnostic IDs: `n_mol["CO2"]` in the liquid phase (replacing
   `"CO2aq"`), matching PHREEQC convention and eliminating the naming
   asymmetry that blocked `derive_speciation_keys` from being a
   one-liner.

2. **BSM2/ADM1 builder migration.** Builders declare cross-phase
   partition reactions (gas ⇌ liquid) and drop explicit
   `speciation_keys=` kwargs and the `speciation_correction()` builder
   method. `derive_speciation_keys` auto-populates the mapping from
   the declared reactions; `speciation_ladders` auto-populates the
   alpha denominator from the connected single-phase equilibria.

3. **ChemistryDatabase packaging.** `ChemistryDatabase` frozen
   dataclass with `.extend()` composition. `ThermoFramework` absorbs
   engine activity configuration and `ThermodynamicConfig`'s
   temperature-correction helpers. Three stock database modules.
   `ControlVolume.chemistry_db` optional field. `ThermoFramework`
   mismatch warning on link construction. Q7 demo.

---

## Scope: Narrow option

This phase ships the **cross-phase partition wiring** and the
**CO₂/NH₄⁺ Species-ID emission fix**. VFA acid-base equilibrium
reactions (pKa split of `S_ac`, `S_pro`, `S_bu`, `S_va`, `H2S`) are
**not** added to the stock databases in this phase. Their alpha
remains 1.0 (no change from today). The generic `_HA`/`_A-` emission
path in `_compute_species_eq` is retained as a **deprecated legacy
fallback** for `EquilibriumDef` entries without `species_refs` (the
VFA rows in `bsm2_default()`). That fallback is removed in the
VFA Full-Closure phase (see § VFA Full-Closure Plan below).

---

## Resolved design decisions

- **Phase-agnostic Species IDs.** `n_mol["CO2"]` replaces
  `n_mol["CO2aq"]` everywhere. Gas and liquid dissolved CO₂ share the
  same `Species(id="CO2", ...)` object. `derive_speciation_keys`
  trivially yields `speciation_keys["CO2"] = "CO2"`.

- **`EquilibriumDef.species_refs` carries the ladder.** New optional
  field `species_refs: Tuple[Species, ...]` on `EquilibriumDef`,
  ordered most-protonated to least. When set, `_compute_species_eq`
  emits by `sp.id`. When absent (legacy VFA rows in `bsm2_default()`),
  the existing `_HA`/`_A-` path fires as before.

- **`_CANONICAL_NAMES` deleted.** CO₂, NH₄⁺, phosphate, bisulfate
  entries in `bsm2_default()` all gain `species_refs`, so the table
  has no remaining consumers. Deleted alongside the canonical-emission
  branch in `_compute_species_eq`.

- **`bsm2_default()` partially migrated.** CO₂ and NH₄⁺ entries gain
  Species references; VFA entries (`S_ac`, `S_pro`, `S_bu`, `S_va`)
  remain string-based and tagged with a deprecation note pointing to
  the VFA Full-Closure plan.

- **Ladder derivation at CV construction.** `derive_speciation_keys`
  extended to build `self.speciation_ladders: Dict[str, List[str]]`
  from the `ReactionSystem`'s single-phase equilibria. `_alpha_for`
  reads from `self.speciation_ladders`; `_CANONICAL_NAMES` import
  deleted from `gas_liquid_link.py`.

- **Cross-phase partition reactions only.** Builders declare gas ⇌
  liquid partition reactions (no `log_K`, no acid-base pKa). VFA
  cross-phase reactions are declared so `speciation_keys` is
  auto-wired; their alpha remains 1.0 because no acid-base equilibrium
  is declared for them.

- **ThermoFramework scope.** Absorbs `use_activity` / `activity_model`
  from the engine constructor and `pKa_at_T` / `Kw_at_T` /
  `context_kwargs_at_current_T` from `ThermodynamicConfig`. No Peng–
  Robinson fugacity surface (trigger-gated,
  [PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md](PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md)).

---

## Checkpoints

Ordered so each leaves the suite runnable. C1–C4 form the emission
core; C5 is the sentinel-bearing CO₂ rename; C6 closes the ladder gap;
C7 completes the builder migration; C8–C11 are the packaging layer;
C12 is the demo and final green.

### C1 — `EquilibriumDef` gains `species_refs`

- [x] [`src/chemistry/equilibria.py`](../../src/chemistry/equilibria.py)
  — add `species_refs: Tuple[Any, ...] = ()` to `EquilibriumDef`
  (typed `Tuple[Species, ...]` once import is clean). No behaviour
  change — field is optional and ignored at this checkpoint.
- [x] Verify: existing tests unaffected (`python -m pytest tests/standalone -q`).

### C2 — Engine helpers carry `species_refs`

- [x] [`src/speciation/engine.py`](../../src/speciation/engine.py)
  `_add_acid` — extract product `StoichiometryEntry` items (exclude
  `H+` and solvents). Order by charge descending (most protonated
  first): insert the acid form (reactant) at position 0. Build
  `species_refs` tuple and pass to `eq_set.add(...)` via a new
  keyword. Update `EquilibriumSet.add()` to accept and store
  `species_refs`.
- [x] Verify: `from_reactions`-built BSM2 engine emits same keys as
  before (no change yet — `_compute_species_eq` still uses
  `_CANONICAL_NAMES`).

### C3 — `_compute_species_eq` emits by Species ID; `_CANONICAL_NAMES` deleted

- [x] [`src/speciation/acid_base.py`](../../src/speciation/acid_base.py)
  — in `_compute_species_eq`: add branch: if `eq_def.species_refs`
  is non-empty, emit `out[sp.id]` for each species in `species_refs`
  (same positional logic as current canonical branch). Keep the
  existing generic `_HA`/`_A-` fallback for entries without
  `species_refs`.
- [x] Delete `_CANONICAL_NAMES` constant and the branch that consults
  it. Update the function's module docstring.
- [x] Verify: tests that assert on canonical names (`"CO2aq"`, `"NH4+"`,
  `"NH3"`) are not yet green — they will pass after C4 migrates
  `bsm2_default()`. Tests that use `from_reactions` path may already
  emit different keys if the bsm2 reactions carry species_refs.
  Run `test_speciation.py` to audit which assertions flip.

### C4 — `bsm2_default()` partially migrated; CO₂/NH₄⁺ entries gain `species_refs`

- [x] [`src/chemistry/equilibria.py`](../../src/chemistry/equilibria.py)
  `bsm2_default()` and related presets (`bsm2_diprotic_co2`,
  `bsm2_with_sulfide`, `adm1_full`) — pass `species_refs` for
  CO₂ (using `common_species.CO2, HCO3_minus, CO3_2minus`),
  NH₄⁺ (using `NH4_plus, NH3`), phosphate, and bisulfate. VFA
  entries (`S_ac`, `S_pro`, `S_bu`, `S_va`, `H2S`) remain
  string-based; add a comment:
  ```python
  # DEPRECATED: string-based entry; will be replaced by
  # Species-declared equilibrium reactions in the VFA
  # Full-Closure phase (see CHEMISTRY_UNIFICATION_3B_CHECKLIST.md).
  ```
- [x] After this checkpoint, `_CANONICAL_NAMES` has zero consumers in
  the `from_reactions` path. The `_HA`/`_A-` generic path fires
  only for the string-based VFA rows.
- [x] `test_speciation.py` assertions on `"CO2aq"` now fail (emission
  is `"CO2"` via `CO2.id`). Assertions on `"NH4+"`, `"NH3"` pass
  (Species ids match). Assertions on `"AceticAcid_HA"`,
  `"AceticAcid_A-"` pass (generic fallback still active).

### C5 — CO₂ rename: `n_mol["CO2aq"]` → `n_mol["CO2"]`

Sentinel-bearing step. Run `test_bsm2_reference.py` sentinels before
and after; expect zero numerical drift (structural rename only).

- [x] [`models/vlmodels/adm1/bsm2.py`](../../models/vlmodels/adm1/bsm2.py)
  — BSM2 liquid initial conditions: `"CO2aq"` → `"CO2"`.
- [x] [`models/vlmodels/adm1/base.py`](../../models/vlmodels/adm1/base.py)
  — any `"CO2aq"` references in ADM1 builder.
- [x] [`tests/standalone/test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py)
  — initial condition dict and sentinel dict: `"CO2aq"` → `"CO2"`.
  Update comment block to record the rename rationale.
- [x] [`tests/standalone/test_cv_advance.py`](../../tests/standalone/test_cv_advance.py)
  — `n_mol` setup sites: `"CO2aq"` → `"CO2"`.
- [x] [`tests/standalone/test_gas_liquid_link.py`](../../tests/standalone/test_gas_liquid_link.py)
  — `n_liq` setup sites: `"CO2aq"` → `"CO2"`.
- [x] [`tests/standalone/test_speciation.py`](../../tests/standalone/test_speciation.py)
  — assertions on `"CO2aq"` → `"CO2"`.
- [x] Grep for any remaining `"CO2aq"` in `src/` and `tests/standalone/`;
  update doc strings in `gas_liquid_link.py` and `common_species.py`
  that mention `"CO2aq"`.
- [x] Run BSM2 sentinels; confirm zero drift; update comment.

### C6 — Ladder derivation; `_CANONICAL_NAMES` import removed from link

- [x] [`src/core/gas_liquid_link.py`](../../src/core/gas_liquid_link.py)
  — add `speciation_ladders: Dict[str, List[str]]` field (default
  empty dict).
- [x] `derive_speciation_keys` — after populating `speciation_keys`
  from cross-phase reactions, also populate `speciation_ladders`:
  for each `(gas_id, liq_id)` pair in the derived mapping, walk
  `reaction_system.single_phase_equilibria` to find all Species
  in the same acid-base chain as `liq_id`. Store as
  `self.speciation_ladders[gas_id] = [sp.id for sp in ladder]`.
  If no single-phase equilibria are available for a species
  (e.g. VFAs under Narrow scope), store a one-element ladder
  `[liq_id]` so alpha naturally evaluates to 1.0 (graceful
  degradation).
- [x] `_alpha_for` — replace `_CANONICAL_NAMES` lookup with
  `self.speciation_ladders.get(species)`. Delete `from ..speciation
  .acid_base import _CANONICAL_NAMES` import.
- [x] [`tests/standalone/test_gas_liquid_link.py`](../../tests/standalone/test_gas_liquid_link.py)
  — update `derive_speciation_keys` tests to also assert on
  `speciation_ladders`; update `test_cv_init_hook_runs_derive`.
- [x] Verify: CO₂ alpha in BSM2 reference test is numerically
  unchanged (ladder `["CO2", "HCO3-", "CO3--"]` gives identical
  result to former `_CANONICAL_NAMES`-based sum).

### C7 — Builder migration: cross-phase partition reactions declared

- [x] [`models/vlmodels/adm1/bsm2.py`](../../models/vlmodels/adm1/bsm2.py)
  — replace `TransferConfig(speciation_keys={"CO2": "CO2aq"})` with
  a declared cross-phase equilibrium reaction in the reaction builder
  for CO₂. `TransferConfig(species={})` — no `speciation_keys` kwarg.
- [x] [`models/vlmodels/adm1/base.py`](../../models/vlmodels/adm1/base.py)
  — replace `.speciation_correction("CO2", "CO2aq")` and the
  VFA/H2S `.speciation_correction(sp, f"{sp}_HA")` calls with
  declared cross-phase partition reactions in `build_adm1_reactions`.
- [x] [`models/vlmodels/fermenter/config/builder.py`](../../models/vlmodels/fermenter/config/builder.py)
  — delete `speciation_correction()` method body (becomes a stub
  that raises `NotImplementedError` with a migration message, then
  deleted after confirming no callers remain).
- [x] [`models/vlmodels/fermenter/config/configs.py`](../../models/vlmodels/fermenter/config/configs.py)
  — delete `speciation_keys` field from `TransferConfig`.
- [x] [`models/vlmodels/fermenter/config/factory.py`](../../models/vlmodels/fermenter/config/factory.py)
  — remove `speciation_keys=...` kwarg from `KineticGasLiquidLink`
  constructor call.
- [x] `KineticGasLiquidLink.__init__` — remove `speciation_keys`
  constructor parameter (field still exists, populated only via
  `derive_speciation_keys`).
- [x] Update `test_gas_liquid_link.py` tests that pass explicit
  `speciation_keys=` to use reaction declarations instead.
- [x] Run BSM2 reference test: zero numerical drift expected (the
  cross-phase reaction declaration is structural wiring only).

### C8 — `ThermoFramework` + `ThermodynamicConfig` absorption

- [x] New module
  [`src/thermo/framework.py`](../../src/thermo/framework.py) —
  `ThermoFramework` frozen dataclass:
  ```python
  @dataclass(frozen=True)
  class ThermoFramework:
      activity_model: str = "davies"
      use_activity: bool = False
      standard_T_K: float = 298.15
      standard_P_atm: float = 1.0
  ```
  Add `pKa_at_T(pKa_ref, dH, T_K)` and `Kw_at_T(T_K)` as static/
  class methods (migrated from `ThermodynamicConfig`).
- [x] New `src/thermo/__init__.py` exposing `ThermoFramework`.
- [x] [`src/chemistry/thermo_params.py`](../../src/chemistry/thermo_params.py)
  — migrate `ThermodynamicConfig.pKa_at_T`, `Kw_at_T`,
  `context_kwargs_at_current_T` to delegate to `ThermoFramework`.
  Mark `ThermodynamicConfig` as deprecated with a docstring note.
  Do **not** delete yet — callers in `make_bsm2_chem_env_fn` still
  reference it; clean up when those callers are removed.
- [x] [`src/speciation/engine.py`](../../src/speciation/engine.py)
  — `SpeciationEngine.__init__` accepts an optional
  `thermo: Optional[ThermoFramework] = None`; when set, reads
  `use_activity` and `activity_model` from it (overriding explicit
  kwargs). Backwards-compatible — explicit kwargs still work.
- [x] Add `src/thermo/` to package `__init__` exports.
- [x] Tests: `test_thermo_framework.py` — construction, `pKa_at_T`,
  `Kw_at_T`, `dataclasses.replace()` override pattern.

### C9 — `ChemistryDatabase` + stock databases + `ControlVolume.chemistry_db`

- [x] New module
  [`src/chemistry/database.py`](../../src/chemistry/database.py) —
  `ChemistryDatabase` frozen dataclass:
  ```python
  @dataclass(frozen=True)
  class ChemistryDatabase:
      thermo: ThermoFramework
      species: Dict[str, Species]
      reactions: ReactionSet

      def extend(self, species=None, reactions=None, thermo=None) -> "ChemistryDatabase":
          ...
  ```
- [x] New package
  [`src/chemistry/databases/`](../../src/chemistry/databases/) —
  `__init__.py`, then three modules:
  - `aqueous.py` — water dissociation + carbonate system + NH₄⁺/NH₃.
    Exports `AQUEOUS_DEFAULT: ChemistryDatabase`.
  - `bioprocess_basic.py` — extends `AQUEOUS_DEFAULT` with O₂
    dissolution placeholder, phosphate, bisulfate. Exports
    `BIOPROCESS_BASIC: ChemistryDatabase`.
  - `anaerobic_digestion.py` — extends `BIOPROCESS_BASIC` with CO₂,
    CH₄, H₂ cross-phase partition declarations (no VFA acid-base
    equilibria — Narrow scope). Exports `AD_BASIC: ChemistryDatabase`.
    Include a `# TODO VFA Full-Closure` comment at the VFA placeholder.
- [x] [`src/core/control_volume.py`](../../src/core/control_volume.py)
  — add optional `chemistry_db: Optional[ChemistryDatabase] = None`
  field. When set, extract `thermo` for passing to the speciation
  engine at construction.
- [x] Tests: `test_chemistry_database.py` — construction, `.extend()`,
  `species` merge, `thermo` override, import from stock modules.

### C10 — `ThermoFramework` mismatch warning

- [x] [`src/core/links.py`](../../src/core/links.py) (or wherever
  `AdvectiveLink`/`DiffusiveLink` are defined) — at construction,
  if both CVs have `chemistry_db` set and their `thermo` instances
  differ, emit:
  ```
  UserWarning: AdvectiveLink connects CVs with different ThermoFrameworks.
  Expected if CVs represent physically distinct phases; potentially a
  bug if they represent zones of the same phase.
  ```
- [x] Tests: one test each for same-thermo (no warning) and
  different-thermo (warning).

### C11 — Demo

- [x] [`demos/model_api/chemistry/chemistry_database.py`](../../demos/model_api/chemistry/chemistry_database.py)
  — demonstrates the create → extend → import lifecycle using
  `AD_BASIC` from `anaerobic_digestion.py`. Shows `.extend()` adding
  a custom species. Verified runnable.

### C12 — Documentation, test cleanup, final green suite

- [x] [`docs/phases-upcoming/CHEMISTRY_UNIFICATION_PLAN.md`](CHEMISTRY_UNIFICATION_PLAN.md)
  — update Phase 3b section status to "shipped"; move to
  `phases-shipped/` alongside the design doc.
- [x] [`docs/architecture.md`](../architecture.md) — describe
  `ChemistryDatabase` and `ThermoFramework` as the user-facing
  chemistry config surface.
- [x] Remaining `"CO2aq"` references in docs/comments: update.
- [x] `test_speciation.py` `test_vfa_emits_generic` — retarget to
  document the **legacy fallback** behaviour; rename to
  `test_vfa_generic_emission_legacy_fallback`. Add new test
  `test_species_based_emission` that uses a `from_reactions`-built
  VFA equilibrium reaction and asserts Species IDs are emitted.
- [x] Verify 913 + N tests green; N ≥ new tests added in C8/C9/C10/C11.

---

## VFA Full-Closure Plan (deferred, trigger-gated)

> **Not implemented in this phase.** Documented here as a
> self-contained checkpoint to be pulled when a model-accuracy
> requirement justifies the BSM2/ADM1 kinetics update.
>
> **Trigger condition:** a concrete use case where VFA gas-transfer
> with α < 1 materially affects simulation output (e.g. high-pH
> headspace stripping of acetic acid at pH 7.5, or H₂S/HS⁻
> equilibrium in a sulfate-reduction scenario).

### What ships

- **VFA acid-base equilibrium reactions** in
  `src/chemistry/databases/anaerobic_digestion.py`. For each VFA
  (acetate, propionate, butyrate, valerate) and H₂S:
  ```python
  # acetate example — add to AD_BASIC stock database
  AceticAcid = Species(id="AceticAcid", ...)  # protonated form
  Acetate    = Species(id="Acetate-",   ...)  # deprotonated
  rxn_acetate = EquilibriumReaction(
      stoichiometry=[
          StoichiometryEntry(AceticAcid, "liquid", -1),
          StoichiometryEntry(Acetate,    "liquid", +1),
          StoichiometryEntry(H_plus,     "liquid", +1),
      ],
      log_K=-4.76, T_ref_K=298.15,
      dH_J_per_mol=None,
  )
  ```
- **Remove deprecated string-based VFA rows** from
  `bsm2_default()`. The `EquilibriumSet.bsm2_default()` preset
  either adopts Species-based declarations or is superseded by
  the stock database path.
- **Delete the generic `_HA`/`_A-` fallback** in
  `_compute_species_eq`. Every `EquilibriumDef` in production paths
  has `species_refs` by this point.
- **n_mol split for VFAs.** After the speciation engine writes the
  protonated fraction to `n_mol["AceticAcid"]` and the deprotonated
  fraction to `n_mol["Acetate-"]`, `n_mol["AceticAcid"]` changes
  meaning from *total acetate* to *protonated fraction only*.
- **BSM2/ADM1 kinetics update.** Every kinetic reaction that currently
  consumes or produces `"AceticAcid"` as a total (e.g.
  aceticlastic methanogenesis) must be rewritten to use the total:
  `env.S("AceticAcid") + env.S("Acetate-")` (or a helper
  `env.total("AceticAcid")`). Equivalent updates for propionate,
  butyrate, valerate, H₂S/HS⁻.
- **Speciation-ladder derivation** for VFAs becomes non-trivial:
  `speciation_ladders["AceticAcid"] = ["AceticAcid", "Acetate-"]`,
  enabling real alpha values at neutral pH.
- **BSM2 reference sentinels re-baseline.** Expect non-trivial
  drift at neutral pH (α_acetate ≈ 0.001 at pH 7.3, pKa 4.76) —
  this is the intended physical correction.

### Pre-work needed before triggering

1. Decide on Species IDs for all VFA deprotonated forms and ensure
   consistency with any downstream consumers (e.g. `S_ac` ADM1
   notation vs `Acetate-` IUPAC-style).
2. Audit all kinetic reactions in `build_bsm2_reactions()` and
   `build_adm1_reactions()` that reference VFA species as totals.
3. Decide whether `env.total(species_id)` is the right helper or
   whether kinetics declare stoichiometry against the total species
   directly.
4. Plan sentinel re-baseline scope and RTOL adjustment.
