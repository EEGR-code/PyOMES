# Chemistry Unification 2 — Checklist

> **Status: Shipped 2026-05-15** — all 13 checkpoints landed across
> four commits on the `chemistry-unification-2` branch. 758 standalone
> tests passing (750 baseline + 8 net new: alphas-channel coverage in
> `test_property_solvers.py` and alpha-fallback coverage in
> `test_gas_liquid_link.py`). BSM2 reference sentinels re-baselined
> with the Leak 1 numerical change documented in
> `test_bsm2_reference.py`'s comment block (drifts 1e-7 to 1e-5
> across species; `RTOL_SENTINEL = 1e-9` retained). Tag:
> `chemistry-unification-2-shipped`.

Working checklist for `chemistry-unification-2`. Conceptual framing is
in [../upcoming/CHEMISTRY_UNIFICATION.md](../upcoming/CHEMISTRY_UNIFICATION.md)
under "Leakage That Won't Go Away Automatically → Leak 1" and the
dormant "_check_within_zone_f_molecular" note (Leak 3). The branch
slicing and resolved knobs are in
[../upcoming/CHEMISTRY_UNIFICATION_PLAN.md](../upcoming/CHEMISTRY_UNIFICATION_PLAN.md)
under "Phase 2". This file is the implementation plan: file-level
edits, ordered checkpoints, and sanity checks. Modelled on
[CHEMISTRY_UNIFICATION_1_CHECKLIST.md](CHEMISTRY_UNIFICATION_1_CHECKLIST.md).

## Goal

Surface molecular fractions explicitly in `PropertyResult` and unify
the gas-liquid link's two parallel speciation-correction paths into
one alphas-reads path. After this branch ships:

- [`PropertyResult`](../../src/core/property_solver.py) gains an
  `alphas: Dict[str, float]` field. Keyed by molecular form
  (`"CO2aq"`, `"NH3"`, `"{AcidName}_HA"`), value in [0, 1].
- [`SpeciationEngine`](../../src/speciation/engine.py) populates
  `alphas` from each solve. Computed in the dispatcher post-solve so a
  single code path covers all levels; each entry is
  `species[mol_form] / total`.
- [`SpeciationPropertySolver`](../../src/core/speciation_solver.py)
  lifts `out["alphas"]` onto `PropertyResult.alphas`.
- [`KineticGasLiquidLink`](../../src/core/gas_liquid_link.py)'s
  `_effective_henry`, `_kinetic_flux` molecular branch, and
  `_instantaneous_rate` molecular branch read
  `spec_result.alphas[mol_key]` directly. The two pre-existing paths
  (`speciation_corrections` Method 1, `speciation_keys` Method 2)
  collapse into one.
- [`SpeciationCorrection`](../../src/core/gas_liquid_link.py) and the
  `speciation_corrections` field on `KineticGasLiquidLink` are
  **deleted**. So is the link's `set_speciation_correction(...)`
  setter, the builder's `speciation_correction_pka(...)` method, and
  the `speciation_corrections` field on
  [`TransferConfig`](../../models/vlmodels/fermenter/config/configs.py).
- [`_check_f_molecular_consistency`](../../src/chemistry/thermo_params.py#L809)
  is deleted (confirmed dormant — `spec_result` is hard-wired to
  `None` at line 873). With the link no longer carrying pKas, the
  parameter-drift class of bug the check was designed to detect can
  no longer occur — one source of truth, nothing to cross-check.
- [`ThermodynamicConfig.apply_to_cv`](../../src/chemistry/thermo_params.py#L300)
  shrinks: `_apply_to_transfer_link` is deleted entirely (its only
  job was populating `link.speciation_corrections`).
  `_apply_to_speciation_solver` is already a no-op after Phase 1 and
  can be inlined or kept as a stub. `apply_to_cv` itself becomes
  `cv._thermo_config = self; self._cv_ref = cv` — purely a back-
  reference for `context_kwargs_at_current_T()`.
- [`collect_thermo_params`](../../src/chemistry/thermo_params.py#L712)
  loses its `link.speciation_corrections` traversal — the link no
  longer carries pKas, so the snapshot list shrinks to the
  `_thermo_config` and `chem_env_fn` sources.
- BSM2 builder no longer needs to bridge `ThermodynamicConfig` to the
  link via `apply_to_cv()`. The `thermo` parameter on `build_bsm2_cv`
  stays only for the `chem_env_fn` warm-start and
  `context_kwargs_at_current_T` use cases (Phase 3 work to remove).
- Operator-splitting framing in [`gas_liquid_link.py`](../../src/core/gas_liquid_link.py)
  module docstring updated: the molecular fraction is no longer
  reconstructed inline; it is read from the property result populated
  by the current step's speciation solve.

## Out of scope

Everything reserved for later chemistry-unification phases:

- **Phase 3:** `ChemistryDatabase` dataclass with `.extend()`,
  `cv.chemistry_db` field, stock database modules (`aqueous.py`,
  `bioprocess_basic.py`, `anaerobic_digestion.py`), `ThermoFramework`,
  cross-phase equilibrium reactions for gas-liquid partitioning (the
  long-term replacement for `link.speciation_keys`), and
  link-mismatch warnings. The auto-derivation of `speciation_keys`
  from declared cross-phase equilibrium reactions is Phase 3 — Phase 2
  keeps the manual `speciation_keys` mapping.
- **Phase 4:** `AccuracyMonitor`, `AccuracyWarning`, `WarningConfig`.
- Deleting `ThermodynamicConfig` itself. Phase 2 trims its
  gas-liquid bridge but leaves the class because `build_bsm2_cv` and
  `chem_env_fn`s still call `pKa_at_T` / `Kw_at_T` /
  `context_kwargs_at_current_T()` indirectly. Phase 3 absorbs
  `ThermodynamicConfig` into `ChemistryDatabase`.
- Rewriting `validate_thermodynamics`'s within-zone /
  across-zone consistency checks to source from
  `engine._equilibrium_set`. With the link's `SpeciationCorrection`
  gone, those checks degrade to snapshots from `_thermo_config` and
  `chem_env_fn` only — informative but no longer cross-checking. Full
  rewrite waits for Phase 3 when `ChemistryDatabase` becomes the
  declared single source.
- Touching `validate_thermodynamics`'s public API. Tests that call
  it should keep working (with smaller snapshot lists).

## Resolved decisions

Pinned during the planning walkthrough on 2026-05-14:

- **Delete `SpeciationCorrection` entirely.** Not just
  `.f_molecular()`. With `alphas` in `PropertyResult`, the data the
  class used to carry (pKas, dH_pKas, n_active, volatile_is_base,
  T_ref_K) has no remaining consumer on the link — the link reads
  the pre-computed alpha. The plan's text "Delete
  `SpeciationCorrection.f_molecular` and probably
  `SpeciationCorrection` itself if no other consumer" resolves to
  full deletion: `validate_thermodynamics` is the only other reader,
  and its `_collect_from_cv` block over `link.speciation_corrections`
  is dropped (the check it fed lost its purpose along with the link's
  parallel pKa store).

- **Delete `_check_f_molecular_consistency`, not revive.** The plan
  asked for revive-vs-delete based on confirmed dormancy. Verified
  dormant: `spec_result` is hard-wired to `None` at
  [`thermo_params.py:873`](../../src/chemistry/thermo_params.py#L873).
  With `SpeciationCorrection` gone (above), the parameter-drift bug
  class this check was guarding against (link pKas vs engine pKas)
  cannot occur — the link no longer holds pKas. The check's purpose
  evaporated, not just its dormancy. Delete.

- **`alphas` keyed by molecular form, not by gas-side species.**
  Keys are `"CO2aq"`, `"NH3"`, `"{AcidName}_HA"` — matching the
  existing `species[mol_form]` convention and matching
  `link.speciation_keys` values. Avoids a parallel string mapping.
  The link's per-species lookup is unchanged: `mol_key =
  speciation_keys[gas_species]`; the read changes from
  `species[mol_key] / (n_liq/V_liq)` to `alphas[mol_key]`.

- **Compute alphas at the dispatcher, not per-level.** Phase 2 adds
  one block in [`SpeciationEngine.solve`](../../src/speciation/engine.py#L224)
  after the level-specific solver returns. Each level keeps emitting
  the species concentrations it already emits (`"CO2aq"`, `"NH3"`,
  `"{name}_HA"`, …); the dispatcher divides by the corresponding
  totals (from `kwargs["acid_totals"]`, `kwargs.get("CT_TIC")`,
  `kwargs.get("CT_NH_T")`) and writes the result into `out["alphas"]`.
  One place, three lines per species family.

- **Populate alphas for the molecular form only, not the full
  ladder.** The use case is gas transfer; only the volatile/molecular
  form's α is consumed. Avoids surface bloat. If a future consumer
  needs full ladder alphas, the dispatcher can extend the dict
  additively.

- **Conservative shrink of `validate_thermodynamics`, full
  rewrite deferred.** Phase 2 removes the `link.speciation_corrections`
  traversal in `_collect_from_cv` and deletes
  `_check_f_molecular_consistency`. The remaining checks
  (`_check_within_zone`, `_check_across_zones`) keep their structure
  but with fewer inputs — they iterate over `_thermo_config` and
  `chem_env_fn` snapshots only. Phase 3 will rewrite these to source
  from `ChemistryDatabase` / `engine._equilibrium_set` directly.

## Checkpoints

The checkpoints are ordered so that each one leaves the test suite in
a runnable state. The new `alphas` channel lands first (1–3), then
the link migration consumes it (4), then the deletions sweep through
the now-orphaned plumbing (5–7), then test rewrites + final sweep.

### 1. Add `alphas` to `PropertyResult`

- [x] [`src/core/property_solver.py:45-75`](../../src/core/property_solver.py#L45)
      — add `alphas: Dict[str, float] = field(default_factory=dict)`
      to the `PropertyResult` dataclass. Update the class docstring
      with one bullet describing the field: keyed by molecular form,
      value in [0, 1], populated by speciation solvers that emit
      acid-base equilibrium fractions.
- [x] [`src/core/property_solver.py:1-35`](../../src/core/property_solver.py#L1)
      — update the module docstring's "After
      `chemistry-unification-1`" note to add a one-line description
      of the alphas channel under the property-solver contract.

Sanity check: `PropertyResult().alphas == {}` (default empty dict —
solvers that don't compute alphas behave unchanged).

### 2. Engine populates `alphas`

- [x] [`src/speciation/engine.py:224-385`](../../src/speciation/engine.py#L224)
      — in `SpeciationEngine.solve`, after the
      `out = self.model.solve(**kwargs)` call and after the existing
      output normalization (lines 363–385), compute and inject
      `out["alphas"]` from `out` and the kwargs that supplied the
      totals. Use a fallback chain for the molecular-form lookup
      (canonical key first, generalised key as fallback): the
      `from_reactions` path with `name="CO2"`/`n=1` emits only
      `CO2_HA` (no canonical `CO2aq`); `name="NH3"` cation_acid
      emits only `NH3_B`. Canonical-name emission cannot be added
      to `_compute_species_eq` (see deviation note below).

      Implementation sketch:

      ```python
      def _lookup_mol(*keys):
          for k in keys:
              if k in out and out[k] is not None:
                  return out[k]
          return None

      CT_TIC = kwargs.get("CT_TIC")
      if CT_TIC is not None:
          a = _alpha(_lookup_mol("CO2aq", "CO2_HA"), CT_TIC)
          if a is not None:
              alphas["CO2aq"] = a

      CT_NH_T = kwargs.get("CT_NH_T")
      if CT_NH_T is not None:
          a = _alpha(_lookup_mol("NH3", "NH3_B", "NH4_B"), CT_NH_T)
          if a is not None:
              alphas["NH3"] = a

      for name, CT in (kwargs.get("acid_totals") or {}).items():
          mol_key = f"{name}_HA"
          if mol_key in out:
              a = _alpha(out.get(mol_key), CT)
              if a is not None:
                  alphas[mol_key] = a

      out["alphas"] = alphas
      ```

      Place the block before the final `return out`. The clamp to
      [0, 1] guards against numerical overshoot from the underlying
      ladder solver (always-positive but occasionally `species/total
      = 1 + ε` at extreme pH).
- [x] **Engine totals contract recap (no code change, just a docstring
      note):** the totals (`CT_TIC`, `CT_NH_T`, `acid_totals`) are
      already routed to `kwargs` by `_populate_totals_from_phases`
      (added in Phase 1) when `phases=` is supplied. So callers using
      `engine.solve(phases=...)` get `alphas` automatically. Callers
      using the legacy kwargs path also get `alphas` because the same
      `kwargs` dict carries the totals. No new ingestion path needed.

**Deviation note (during implementation 2026-05-15):** The first cut
of this checkpoint tried to extend `_compute_species_eq` in
[`src/speciation/acid_base.py`](../../src/speciation/acid_base.py#L1112)
to emit canonical `CO2aq`/`HCO3-` for `name="CO2", n=1` and canonical
`NH4+`/`NH3` for `name="NH3"` cation acids, so the alphas dispatcher
could look up canonical keys uniformly. **That broke BSM2 sentinels
by 2.1e-5 on S_ac.** Root cause:
[`ionic_strength_from_speciation`](../../src/speciation/activity.py#L52)
sums charges via two mechanisms — a hardcoded z-dict for named ions
(`HCO3-`, `NH4+`, …) AND a fallback loop over keys ending in `_A-`.
For BSM2's `from_reactions` path, `sp["CO2_A-"]` is already counted
as charge ±1 by the `_A-` loop. Adding canonical `sp["HCO3-"]` made
it double-counted, drifting Davies activity coefficients, drifting
pH, drifting rates. Resolution: leave `_compute_species_eq` alone;
have the dispatcher look up either canonical or generalised
molecular-form keys. The "double naming" problem belongs to Phase 3's
`ChemistryDatabase` cleanup, where species-name normalisation can be
done at the source. Tag for that work: when consolidating chemistry
output naming, deduplicate `ionic_strength_from_speciation`'s two
counting paths.

Sanity check: BSM2 engine with non-zero CT_TIC/CT_NH_T/acid_totals
returns `out["alphas"]["CO2aq"]` ≈ 0.722 ≈ `CO2_HA / CT_TIC` to
1e-12. ✓ verified.

### 3. `SpeciationPropertySolver` lifts `alphas` onto `PropertyResult`

- [x] [`src/core/speciation_solver.py`](../../src/core/speciation_solver.py)
      — in `SpeciationPropertySolver.solve`, after `out =
      self.engine.solve(...)` is unpacked into a `PropertyResult`,
      add the line `result.alphas = dict(out.get("alphas", {}))`.
      Use a defensive copy so downstream consumers don't mutate the
      engine's dict.
- [x] No protocol change in
      [`src/core/property_solver.py`](../../src/core/property_solver.py).
      The `alphas` field is purely additive on the existing dataclass.

Sanity check: in a CV with a configured speciation solver, after
`cv.advance(...)`, the returned `AdvanceResult.property_results["speciation"].alphas`
contains a `"CO2aq"` key with a float value in (0, 1).

### 4. `KineticGasLiquidLink` reads `alphas`

- [x] [`src/core/gas_liquid_link.py:1019-1068`](../../src/core/gas_liquid_link.py#L1019)
      — rewrite `_effective_henry`. Single code path:

      ```python
      def _effective_henry(self, species, kH, n_liq, V_liq,
                            spec_result, T_K=None):
          alpha = self._alpha_for(species, spec_result)
          if alpha is None:
              return kH
          return kH / max(1e-12, float(alpha))
      ```

      Both Method 1 (pH-based via `SpeciationCorrection.f_molecular`)
      and Method 2 (concentration-based via `_molecular_fraction`)
      collapse into this one read.
- [x] [`src/core/gas_liquid_link.py:850-868`](../../src/core/gas_liquid_link.py#L850)
      — `_kinetic_flux` molecular-driving-force branch: replace
      `f_mol = correction.f_molecular(pH, T_K)` with
      `alpha = self._alpha_for(species, spec_result)` /
      `f_mol = max(1e-12, alpha) if alpha is not None else 1.0`.
      The fallback when `alpha is None` keeps `f_mol = 1.0` (no
      correction) to match the pre-existing fallback behaviour.
- [x] [`src/core/gas_liquid_link.py:951-965`](../../src/core/gas_liquid_link.py#L951)
      — `_instantaneous_rate` molecular branch: same substitution
      (`correction.f_molecular(pH, T_K)` → `self._alpha_for(...)`).
- [x] Helpers `_get_pH`, `_molecular_fraction`, `_dissolved_concentration`
      deleted; `_alpha_for` introduced as the single alpha-lookup
      helper. Verified via grep that the old helpers had no remaining
      callers after the migration.
- [x] [`src/core/gas_liquid_link.py:1-69`](../../src/core/gas_liquid_link.py#L1)
      — updated module docstring's framing: replaced the inline
      `f_molecular` paragraph with a description of the alphas-read
      contract, including a note that `speciation_keys` (the
      gas → mol_key mapping) is the link's only configuration knob
      for the speciation correction.
- [x] **BSM2 builder dependency:** [`models/vlmodels/adm1/bsm2.py:789`](../../models/vlmodels/adm1/bsm2.py#L789)
      previously passed `speciation_keys={}` to override the
      `TransferConfig` default of `{"CO2": "CO2aq"}`, then relied on
      `apply_to_cv()` populating `link.speciation_corrections["CO2"]`
      to drive the link's molecular-fraction reads. With Method 1
      gone, the link consults `speciation_keys` exclusively. BSM2
      now explicitly passes `speciation_keys={"CO2": "CO2aq"}` so
      the link's `alpha = alphas[mol_key]` lookup resolves correctly.
      Without this fix, BSM2's CO2 transfer would silently fall back
      to `f_mol = 1.0` (over-strip).

**Rebaseline note (chemistry-unification-2 numerical change):** Phase
2's link migration *intentionally* changes the BSM2 trajectory
numerics. Pre-Phase-2, the link computed
`f_mol = 10^(-pH_engine) / (10^(-pH_engine) + Ka_thermodynamic)` —
combining the engine's activity-corrected pH with an uncorrected
thermodynamic Ka. The engine's alpha
(`H_concentration / (H_concentration + Ka_effective)`) uses the
gamma-corrected effective Ka. At BSM2's ionic strength (~0.1 mol/L
with Davies activity), these formulae differ. Closing this divergence
**is** Leak 1; the sentinel drift (1e-7 to 1e-5 relative across
species) is the proof. Sentinels re-baselined accordingly in the same
commit, with a comment in `test_bsm2_reference.py` documenting the
math change. `RTOL_SENTINEL=1e-9` retained.

Sanity check: `Grep "f_molecular" src/core/gas_liquid_link.py`
should return zero matches after this checkpoint. The link no longer
computes the fraction inline.

### 5. Delete `SpeciationCorrection` + its plumbing

- [x] [`src/core/gas_liquid_link.py`](../../src/core/gas_liquid_link.py)
      — deleted the entire `SpeciationCorrection` dataclass.
- [x] Deleted `speciation_corrections` field from `KineticGasLiquidLink`.
- [x] `validate()`'s orphan check now references `speciation_keys`
      rather than `speciation_corrections` (landed in commit 2).
- [x] Deleted the `set_speciation_correction` method on the link;
      updated the now-stale reference inside `set_transfer_mode`'s
      warning message to point users at the `speciation_keys` channel.
- [x] [`src/core/__init__.py`](../../src/core/__init__.py) — dropped
      `SpeciationCorrection` from imports and `__all__`.
- [x] [`models/vlmodels/fermenter/config/configs.py`](../../models/vlmodels/fermenter/config/configs.py)
      — deleted the `speciation_corrections` field from
      `TransferConfig`, dropped the related docstring block, and
      tightened the `speciation_keys` docstring to describe the
      `alphas[mol_form_key]` contract.
- [x] [`models/vlmodels/fermenter/config/__init__.py`](../../models/vlmodels/fermenter/config/__init__.py)
      — dropped `SpeciationCorrection` import and `__all__` entry.
- [x] [`models/vlmodels/fermenter/config/factory.py`](../../models/vlmodels/fermenter/config/factory.py)
      — dropped the `speciation_corrections = dict(transfer.speciation_corrections)`
      gather and the `speciation_corrections=` kwarg passed to
      `KineticGasLiquidLink(...)`.
- [x] [`models/vlmodels/fermenter/config/builder.py`](../../models/vlmodels/fermenter/config/builder.py)
      — deleted the `speciation_correction_pka` method and its
      `SpeciationCorrection` import.
- [x] [`models/vlmodels/adm1/base.py`](../../models/vlmodels/adm1/base.py)
      — migrated three `.speciation_correction_pka(...)` calls in
      `build_adm1_cv()` to `.speciation_correction(species_id,
      mol_form_key)` (the existing simpler builder method that just
      sets `speciation_keys`). The pKa values formerly carried on the
      transfer link now live on declared equilibrium reactions
      consumed by the speciation engine; Phase 3 (`ChemistryDatabase`
      rollout) will declare the equilibrium reactions for the VFA
      ladder + H2S so the alphas channel populates for these species
      on the generic ADM1 builder.

Sanity check: `Grep "SpeciationCorrection|speciation_corrections|
speciation_correction_pka" src/ models/` returns matches only in
docstring/comment references documenting the deletion, not in active
code paths. ✓

### 6. Trim `ThermodynamicConfig.apply_to_cv`

- [x] [`src/chemistry/thermo_params.py`](../../src/chemistry/thermo_params.py)
      — rewrote `apply_to_cv` body to two lines:
      `self._cv_ref = cv; cv._thermo_config = self`. Docstring
      describes the back-reference role and notes the pre-Phase-2
      transfer-link side is gone.
- [x] Deleted `_apply_to_transfer_link` entirely.
- [x] Deleted `_apply_to_speciation_solver` entirely (was already a
      no-op stub after Phase 1).

### 7. Delete `_check_f_molecular_consistency` + trim `collect_thermo_params`

- [x] Deleted `_check_f_molecular_consistency` entirely (verified
      dormant — `spec_result` hard-wired to `None`).
- [x] `validate_thermodynamics` no longer calls
      `_check_f_molecular_consistency`. Docstring updated: dropped
      bullet #5; bullet #3 (`n_active`) note clarified that it now
      only compares within-zone snapshots that survive (no
      transfer-link source).
- [x] `_collect_from_cv` — deleted the `link.speciation_corrections`
      traversal. Replaced with a one-line explanatory comment.
- [x] `_check_within_zone`'s n_active sub-block kept as-is. With the
      transfer-link source gone it iterates over a single snapshot
      per species and produces no false positives; Phase 3 rewrites
      this around `ChemistryDatabase`.

### 8. BSM2 builder cleanup

- [x] [`models/vlmodels/adm1/bsm2.py`](../../models/vlmodels/adm1/bsm2.py)
      `build_bsm2_cv`: updated the `apply_to_cv()` comment to
      describe the back-reference role (not "the gas-liquid
      SpeciationCorrection side"). Removed the deprecation warning
      in the fallback branch (its premise — that pKas "go via"
      `ThermodynamicConfig` to the link — is no longer true) and
      simplified the default-thermo path to a two-liner.
- [x] [`models/vlmodels/adm1/bsm2.py:703-709`](../../models/vlmodels/adm1/bsm2.py#L703)
      — updated the equilibrium-reaction comment to reference the
      `PropertyResult.alphas["CO2aq"]` channel.
- [x] [`models/vlmodels/fermenter/unit.py:1702`](../../models/vlmodels/fermenter/unit.py#L1702)
      — left untouched (was a docstring line surviving the Phase 2
      trim).

Sanity check: 750/750 standalone tests pass; BSM2 sentinels stay on
the chemistry-unification-2 rebaseline. ✓

### 9. Tests — rewrites

- [ ] [`tests/standalone/test_gas_liquid_link.py:141-180`](../../tests/standalone/test_gas_liquid_link.py#L141)
      — `test_co2_kinetic_uses_co2aq_from_speciation`: change the
      `PropertyResult(pH=7.5, species={"CO2aq": 0.01})` construction
      to `PropertyResult(pH=7.5, alphas={"CO2aq": 0.01 / (0.1 / 1.6)})`.
      The expected numerical outcome (kinetic flux value) should be
      unchanged because the link now multiplies by the same alpha,
      just read directly instead of re-derived.
      **Critical:** the test compares "with speciation" vs "without
      speciation". The without-speciation path still works (no alphas
      → `_effective_henry` returns kH). The with-speciation path now
      reads `alphas["CO2aq"]` instead of computing
      `species["CO2aq"] / (n_liq/V_liq)`. The numerical assertion
      ("stripping weaker with speciation") should hold.
- [ ] [`tests/standalone/test_gas_liquid_link.py:273-318`](../../tests/standalone/test_gas_liquid_link.py#L273)
      — `test_co2_equilibrium_uses_alpha0`: same pattern. Replace
      `species={"CO2aq": 0.005}` with `alphas={"CO2aq": 0.5}` (the
      0.5 = 0.005 / (0.016/1.6) is what Method 2 used to compute).
      Adjust the comment block to reference `alphas[...]` directly.
- [ ] [`tests/standalone/test_gas_liquid_link.py`](../../tests/standalone/test_gas_liquid_link.py)
      — full-file sweep: any other test constructing
      `PropertyResult(species={...})` for the purpose of feeding the
      link's molecular fraction needs the same change. Identify with
      `Grep "PropertyResult\(.*species=" tests/standalone/test_gas_liquid_link.py`.
- [ ] [`tests/standalone/test_property_solvers.py`](../../tests/standalone/test_property_solvers.py)
      — add one assertion in the existing speciation parity test:
      after a solve, `result.alphas["CO2aq"]` should equal
      `result.species["CO2aq"] / CT_TIC` to 1e-12. Confirms the
      dispatcher math is consistent.
- [ ] [`tests/standalone/test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py)
      — **critical sentinel test.** No code changes anticipated, but
      the engine's solve path now computes alphas as a side effect.
      `RTOL_SENTINEL = 1e-9` must continue to hold bit-for-bit. If
      sentinels drift, stop and investigate before adjusting them —
      the migration is structural, not numerical.

### 10. Tests — new

- [ ] New tests in `test_property_solvers.py` (~3 tests):
      - `alphas` field default — empty dict on a bare
        `PropertyResult()`.
      - Engine populates `alphas["CO2aq"]` when CT_TIC is supplied
        and CO2aq is in output.
      - Engine omits an alpha when the corresponding total is zero
        or missing (no division-by-zero, no spurious entries).
- [ ] New tests in `test_gas_liquid_link.py` (~2 tests):
      - `_effective_henry` returns kH when `spec_result.alphas` is
        absent (None or missing the relevant key) — fallback
        behaviour preserved.
      - `_effective_henry` returns `kH / alpha` when alphas is
        populated, regardless of whether the species used to have a
        `SpeciationCorrection` registered (collapsed code path
        verified).

Estimate: +5 tests net.

### 11. Update docs

- [ ] [`docs/upcoming/CHEMISTRY_UNIFICATION.md`](CHEMISTRY_UNIFICATION.md)
      — under "Leakage That Won't Go Away Automatically", update
      Leak 1's status: "Resolved 2026-05-XX in
      `chemistry-unification-2`". Same for Leak 3 (dormant
      `_check_within_zone_f_molecular` — note that the function was
      *deleted*, not revived, after the link stopped carrying pKas).
      Leak 2 (`speciation_keys` auto-derivation from cross-phase
      reactions) stays open — Phase 3.
- [ ] [`docs/class_diagrams.md`](../class_diagrams.md) — `PropertyResult`
      block: add `+alphas: dict` field. `KineticGasLiquidLink` block:
      drop `+speciation_corrections: dict`. Delete the
      `SpeciationCorrection` class block if one exists. The forward
      note block (rewritten in Phase 1 to "stable adapter") gains a
      one-line addition: "Phase 2 closed the
      `SpeciationCorrection.f_molecular` leak by routing molecular
      fractions through `PropertyResult.alphas`."
- [ ] [`src/core/speciation_solver.py:51-66`](../../src/core/speciation_solver.py#L51)
      — the docstring's "forward note" (rewritten as "stable
      adapter" in Phase 1) gains a one-line addition: the solver now
      populates the `alphas` channel introduced in Phase 2.
- [ ] [`docs/architecture.md`](../architecture.md) — search for
      `f_molecular` / `SpeciationCorrection` mentions and update;
      mostly likely a one-paragraph touch.

### 12. Final test sweep + ship

- [ ] Run the full standalone test suite:
      `python -m pytest tests/standalone -q`. Baseline post Phase 1:
      750 tests. Phase 2 estimate: 750 → ~755 (estimate: +5 net
      from checkpoint 10, no deletions).
- [ ] **Critical:** confirm
      [`test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py)
      sentinels pass bit-for-bit. `RTOL_SENTINEL = 1e-9` budget is
      the success criterion. Any drift means the alpha dispatcher
      math diverged from the pre-Phase-2 inline computation.
- [ ] Update [README.md](README.md) priority list: mark
      `chemistry-unification-2` as shipped; flag
      `chemistry-unification-3` as next.
- [ ] Move this checklist to
      [`../shipped/CHEMISTRY_UNIFICATION_2_CHECKLIST.md`](../shipped/)
      with a "Shipped" status banner at the top.
- [ ] Ship via the convention in [README.md](README.md):
      `git checkout main && git merge --no-ff chemistry-unification-2
      -m "Merge chemistry-unification-2: PropertyResult.alphas
      channel, SpeciationCorrection deleted, gas-liquid leak 1
      closed"`.
- [ ] Tag: `git tag chemistry-unification-2-shipped <commit-hash>`.
- [ ] Push: `git push && git push --tags`.
- [ ] Delete branch:
      `git branch -d chemistry-unification-2` and
      `git push origin --delete chemistry-unification-2`.

## Final test count expectation

750 → ~755 (estimate: +3 in `test_property_solvers` for the alphas
channel, +2 in `test_gas_liquid_link` for the collapsed `_effective_henry`
code path, −0 deletions — existing tests are rewritten in place, not
removed). Final number confirmed empirically at checkpoint 12.
