# LAYER1_GAP_CLOSURE — Phase 2 of 2

> **Status:** Shipped 2026-07-03. Tag `layer1-gap-closure-shipped` (after
>   merge). CP1–CP7 all landed, plus one addition beyond the original scope:
>   the `src/speciation/` → `src/chemical_equilibrium/` module rename (§18's
>   own deferred decision), agreed during CP6 to execute at CP7 rather than
>   defer again. Real findings surfaced and fixed along the way, not just
>   mechanical execution: a species-id collision bug in `nr_solver.py`'s
>   internal dicts (gas and liquid secondaries sharing a bare id, e.g. both
>   `"CO2"`, silently overwrote each other — fixed via `SecondaryEntry.c_key`,
>   CP2); `step_internal_transfer()` had *no* scope filter at all before this
>   phase, not a partial one as the plan assumed (CP3); `WaterVapourBoundary`'s
>   only real consumer never tracked liquid-phase H2O, so retiring it required
>   also seeding a real liquid water pool, not just swapping the mechanism
>   (CP4); a folded-vs-unfolded precipitation comparison surfaced a genuine
>   physical effect (CO₂ stripping promotes CaCO₃ scaling) worth locking in
>   as a regression rather than a bug to paper over (CP5); a naive rename
>   script corrupted several "old name → new name" mapping statements into
>   nonsensical self-references, caught before committing (CP6).
> **Depends on:** [EQUILIBRIUM_CONSTRAINT_UNIFICATION.md](../phases-shipped/EQUILIBRIUM_CONSTRAINT_UNIFICATION.md)
>   (Phase 1) — shipped 2026-07-02, tag `equilibrium-constraint-unification-shipped`.
>   This phase consumes its auto-classification path, `*Equilibrium` sibling
>   types, and `activity_for_entry()` helper.
> **Design authority:** [MASS_EXCHANGE_ARCHITECTURE.md](../design/MASS_EXCHANGE_ARCHITECTURE.md)
>   §10.2 (Layer 1 algebraic solve), §10.4 (status table — updated, "⚠ gap"
>   rows now closed), §14 (constraint families, `NRTableau` folding limits);
>   [THERMODYNAMIC_MODEL_ARCHITECTURE.md](../design/THERMODYNAMIC_MODEL_ARCHITECTURE.md)
>   §6 (gas-liquid-solid coupling), §8.4/§8.5 (deferred Jacobian primitives —
>   §8.4's `DifferentiableLiquidModel` shipped here as standalone groundwork,
>   not yet wired into the inner NR loop), §8.6 (water activity / osmotic
>   coefficient — new tracked follow-up surfaced during CP4, not fixed here);
>   [CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md](../design/CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md)
>   §18 (rename — shipped here, CP6/CP7).
> **Tag:** `layer1-gap-closure-shipped`.

---

## Background and motivation

`MASS_EXCHANGE_ARCHITECTURE.md` §10.4 has stood as "⚠ gap" since the doc's
first draft: `EquilibriumTransferModel` runs `cv.advance()`'s
`step_internal_transfer()` *after* the speciation engine's `solve()`, as a
sequential correction, rather than as rows in the same `g(y,z) = 0` system.
For CO₂/HCO₃⁻ where gas-liquid and acid-base equilibria are tightly coupled
(§10.5: high kLa, > ~500 h⁻¹), this sequential splitting introduces O(dt)
error per step. This phase closes that gap for gas-liquid (Henry, and Raoult
as its own checkpoint), building on Phase 1's declaration-side unification.

Precipitation (Ksp) is not a separate gap to close here — it's already folded
into `NRSpeciationEngine.solve()`'s single call, via a nested active-set loop.
This phase's job re: precipitation is narrower: verify it continues working
unchanged alongside newly-folded gas-liquid rows (CP5), not to change its
numerics.

---

## Key design decisions

**Mass-action rows draw non-ideality corrections directly from `ThermoFramework`
primitives, not by wrapping `PartitionModel`.** Per the wrapping litmus test
(`MASS_EXCHANGE_ARCHITECTURE.md` §14.2): `MultispeciesVLEPartition` internally
*solves* the coupled VLE, so wrapping it inside an outer Newton row would nest
a solver inside a solver. Instead, folded rows call
`GasEOS.partial_pressures_atm()` and `LiquidPhaseModel.gamma_all()` directly at
the current Newton trial point (both pure evaluate functions), using Phase 1's
CP4 dispatch helper, and let the tableau's own iteration do the fixed-point
search.

**`NRTableau`'s one-master-per-connected-component limitation is confirmed
safe for realistic gas-liquid folding, not safe for general multi-ion Ksp.**
Verified by reading `build_tableau()` directly (2026-07-01 investigation):
Henry-type entries either attach onto an existing single-component ladder
(CO₂, NH₃, H₂S) or form a trivial standalone singleton (inert gases: O₂, CH₄,
N₂, H₂) — neither merges two previously independent multi-species components.
This phase folds gas-liquid generally; it does **not** attempt to generalize
`NRTableau`'s master selection for multi-ion minerals (that's
`MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md`'s territory) — single-
driving-ion Ksp keeps using the existing nested active-set loop.

**Ideal-gas case is the primary target; non-ideal EOS coupling is explicitly
deferred within this phase.** `∂φ_i/∂y_j` (the `GasEOS` Jacobian, §8.5 of
`THERMODYNAMIC_MODEL_ARCHITECTURE.md`) doesn't exist yet. For `IdealGasEOS`
(`φ_i = 1` identically) this term is trivially zero, which covers VLsim's
near-term bioprocess use case (≤ ~1 atm). Folding with `PengRobinsonEOS`
requires that Jacobian primitive first — tracked as a follow-on, not blocking
this phase.

**Raoult/evaporation folds in this phase, as its own distinct checkpoint** (per
explicit scoping decision) — not merged into the Henry checkpoint, since it
also retires `WaterVapourBoundary`/`VentWaterLoss` as a side effect, a
separate migration from the core Henry-folding work.

**The `SpeciationEngine` → `ChemicalEquilibriumEngine` rename (§18) closes this
phase**, once CP1-CP5 demonstrate the engine actually resolves gas-liquid
(and continues to resolve solid-liquid) simultaneously with acid-base — the
point at which the new name stops being aspirational.

---

## Scope

**In scope:**

- `build_tableau()` changes: stop skipping gas-liquid-classified reactions
  (Phase 1's classifier); add gas-phase species as graph nodes using the
  activity-dispatch helper; attach onto existing components or form
  singletons per the confirmed-safe cases above. Raise a clear
  `ConfigurationError` (not silent misbehavior) if a declared gas-liquid
  reaction would require merging two independent multi-species components —
  explicitly out of this phase's capability, point the error at
  `MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md`.
- Residual/Jacobian assembly for folded mass-action rows: `φ_i`/`γ_i` read
  from `GasEOS`/`LiquidPhaseModel` at the current NR trial point (per the
  wrapping litmus test), not via `MultispeciesVLEPartition` or any other
  solve-oriented `PartitionModel`. `IdealGasEOS` path only; `PengRobinsonEOS`
  path explicitly out of scope pending the `GasEOS` Jacobian primitive.
- `algebraic_species()` expansion: engine reports folded gas-liquid species as
  algebraic, alongside already-included solid-liquid species.
- `step_internal_transfer()` scope-filter fix (§13.5 of
  `MASS_EXCHANGE_ARCHITECTURE.md`, currently unimplemented): skip **both**
  `KineticTransferModel` and `EquilibriumTransferModel` for engine-owned
  species, not just the equilibrium one.
- Raoult/evaporation folding as its own checkpoint: fold `H2O` gas-liquid
  equilibrium the same way as Henry; retire `WaterVapourBoundary` and
  `VentWaterLoss`.
- Precipitation regression verification: confirm the existing nested
  active-set loop is unaffected by newly-folded gas rows in the same
  `solve()` call; explicit regression tests covering a CV with both
  precipitation and folded gas-liquid active simultaneously.
- The rename: `SpeciationEngine` → `ChemicalEquilibriumEngine`,
  `NRSpeciationEngine` → `NRChemicalEquilibriumEngine`,
  `PHREEQCEngine` → `PHREEQCChemicalEquilibriumEngine`,
  `NumericalGradientEngine` → `NumericalGradientEquilibriumEngine`,
  `SpeciationEngineProtocol` → `ChemicalEquilibriumEngineProtocol` — full
  mapping per §18 of `CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md`. Module
  rename (`src/speciation/` → `src/chemical_equilibrium/`) decided at this
  checkpoint, not before (§18 explicitly deferred this).
- Doc closeout: `MASS_EXCHANGE_ARCHITECTURE.md` §10.4 status table updated
  (Henry/Raoult row: gap closed; Ksp row: unchanged, nested-loop noted as
  accepted); `CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md` status header and
  §18 marked shipped; both phase docs moved to `phases-shipped/`.

**Out of scope (tracked elsewhere):**

- `PengRobinsonEOS` wiring / non-ideal gas Jacobian — §8.5 of
  `THERMODYNAMIC_MODEL_ARCHITECTURE.md`, separate future phase.
- Multi-ion Ksp as one flat Jacobian — `MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md`,
  separate standalone pathway; this phase only guarantees single-ion Ksp
  keeps working, not that it becomes a flat fold.
- Active-set Jacobian discontinuity at precipitation on/off boundaries (§7 of
  `THERMODYNAMIC_MODEL_ARCHITECTURE.md`) — remains an open design gate for
  the white-box BDF path; this phase doesn't need precipitation to be
  differentiable across the switch, only correct at converged points.
- Isotherm family (Langmuir/Freundlich) — still deferred, no concrete use
  case.
- `PropertyCalculator`'s `phase.speciation` → `engine.current_result`
  migration gap (noted in `CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md` §8) —
  harmless under the current explicit-Euler path; opportunistic fix here if
  convenient, not required for this phase's correctness.
- Surface complexation / site conservation (§12 Q9) — deferred, unrelated to
  this phase.
- Water activity / osmotic coefficient (`a_H2O = 1` assumed everywhere,
  model-invariant — Davies/SIT have no term for a neutral solvent) — found
  during CP4, now tracked at §8.6 of `THERMODYNAMIC_MODEL_ARCHITECTURE.md`.
  Relevant to high-solids/dry AD and concentrated fermentation broths;
  requires water becoming a real tracked NR master (deferred separately,
  see CP4's design discussion) plus a Pitzer/SIT osmotic-coefficient
  method, neither attempted in this phase.

---

## Checkpoints

### CP1 — `build_tableau()` gas-liquid folding

- [x] Stop filtering gas-liquid-classified reactions out of `build_tableau()`'s
      graph (`src/speciation/nr_tableau.py`); add gas-phase species as nodes
      using Phase 1's classifier output.
- [x] Confirm attach-onto-existing-component / form-singleton logic handles
      both cases correctly (CO₂/NH₃/H₂S onto existing ladders; O₂/CH₄/N₂/H₂ as
      singletons) — regression tests for each.
- [x] Add explicit detection + `ConfigurationError` for the case where a
      gas-liquid reaction would bridge two independent multi-species
      components (out of this phase's scope) — test with a constructed
      pathological case to confirm the error fires rather than silently
      producing wrong results.

### CP2 — Mass-action row residual/Jacobian assembly (ideal-gas case)

- [x] Wire Phase 1's CP4 activity-dispatch helper into the NR residual
      evaluation for gas-liquid rows: `g_i = log(φ_i·y_i·P) − log(γ_i·kH_i·x_i)`
      (or equivalent log-linear form matching the existing tableau
      convention), evaluated at the current Newton trial point.
- [x] `IdealGasEOS` path only (`φ_i = 1`, `∂φ_i/∂y_j = 0` — no Jacobian
      contribution needed from the gas side for this checkpoint).
- [x] Liquid-side Jacobian contribution (`∂γ_i/∂C_j`) — implement the
      `DifferentiableLiquidModel` extension from §8.4 of
      `THERMODYNAMIC_MODEL_ARCHITECTURE.md` for Davies/SIT
      (`∂γ_i/∂C_j = (∂γ_i/∂I)(z_j²/2)`, "cheap analytically" per that section).
- [x] Tests: a CV with CO₂ gas-liquid + carbonate ladder converges to the same
      equilibrium as today's SNIA sequential path at low kLa (consistency
      check), and diverges from it in the expected direction at high kLa
      (demonstrating the gap this phase closes — compare against a
      known-accurate reference, e.g. a very fine SNIA sub-stepping as ground
      truth).

### CP3 — `algebraic_species()` + `step_internal_transfer()` scope-filter fix

- [x] Expand `algebraic_species()` on the engine to include newly-folded
      gas-liquid species alongside already-included solid-liquid species.
      **Verified it already did — it unions bare `species_id` across all
      secondaries regardless of `phase`, so no code change was needed.**
- [x] Fix `step_internal_transfer()` (§13.5 gap) to skip **both**
      `KineticTransferModel` and `EquilibriumTransferModel` for engine-owned
      species. **Found no filter existed at all before this phase (not a
      partial one, as this bullet's premise assumed) — added one, scoped to
      `NRChemicalEquilibriumEngine.gas_liquid_species()` (species with an
      actual folded gas row), deliberately narrower than the full
      `algebraic_species()` set to avoid disabling legitimate
      `transfer_models` usage for species with ordinary acid-base chemistry
      but no folded Henry/Raoult row.**
- [x] Tests: a CV with a folded species also carrying a (now-redundant)
      registered `KineticTransferModel` for the same species confirms the
      kinetic model is skipped, not double-applied.

### CP4 — Raoult / evaporation folding

- [x] Apply the same folding pattern (CP1/CP2) to `H2O` via a
      `RaoultEquilibrium` instance (Phase 1's type, satisfying
      `EquilibriumConstraint`).
- [x] Retire `WaterVapourBoundary` and `VentWaterLoss` per §6.2 of
      `MASS_EXCHANGE_ARCHITECTURE.md` — confirm no other call sites depend on
      them before deletion (repeat the CP0-style dead-code audit pattern from
      `EQUILIBRIUM_RESULT.md` if any ambiguity).
- [x] Tests: liquid ↔ gas H2O mass balance closes automatically without the
      retired boundary objects; existing water-vapor-saturation test coverage
      migrated to the folded path.

### CP5 — Precipitation regression verification

- [x] Regression tests: a CV with both precipitation (existing nested
      active-set) and newly-folded gas-liquid rows active simultaneously in
      the same `solve()` call — confirm both converge correctly and
      independently (no cross-contamination between the outer active-set loop
      and the newly-folded gas rows' Newton unknowns).
- [x] Confirm `NR_PRECIPITATION_CV_INTEGRATION.md`'s `SolidPhase` writeback (if
      shipped by this point) is unaffected by the gas-liquid folding changes
      in CP1/CP2. **It hadn't shipped (status header still read "Upcoming —
      not yet started") — documented as a no-op finding, not silently
      skipped.**
- [x] No changes to precipitation's numerics themselves — this checkpoint is
      verification only.

### CP6 — Rename (`SpeciationEngine` → `ChemicalEquilibriumEngine`)

- [x] Apply the full identifier mapping from §18 of
      `CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md`: `SpeciationEngine` →
      `ChemicalEquilibriumEngine`, `NRSpeciationEngine` →
      `NRChemicalEquilibriumEngine`, `PHREEQCEngine` →
      `PHREEQCChemicalEquilibriumEngine`, `NumericalGradientEngine` →
      `NumericalGradientEquilibriumEngine`, `SpeciationEngineProtocol` →
      `ChemicalEquilibriumEngineProtocol`. `algebraic_species()`,
      `GrayBoxEngineProtocol`, `WhiteBoxEngineProtocol` unchanged (already
      generic names).
- [x] Decide module rename (`src/speciation/` → `src/chemical_equilibrium/`)
      now, per §18's deferred decision — either execute it or explicitly
      record why it's kept as-is with an updated docstring. **Decided here:
      keep the module path as `src/speciation/` for this checkpoint (lower
      risk, narrower diff — the class rename alone already touched 80
      files), with the module rename itself moved explicitly into CP7's
      scope rather than deferred indefinitely. Executed at CP7.**
- [x] Backward-compat aliases for one phase (old names as deprecated
      re-exports) given the size of the rename's blast radius, matching the
      pattern used for `DaviesActivityModel = DaviesLiquidModel` etc. in
      `THERMODYNAMIC_MODEL_ARCHITECTURE.md`. **Removed 2026-07-03**, same day:
      a grep across `src/`, `tests/`, `demos/`, `docs/` for the old names
      turned up nothing outside the five alias-definition lines themselves, so
      the "one phase" window closed immediately rather than carrying dead
      aliases forward. See §18 of `CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md`.
- [x] Mechanical sweep across `src/`, `tests/`, `demos/`, `docs/` for the old
      names; full suite green.

### CP7 — Doc closeout + phase move

- [x] Update `MASS_EXCHANGE_ARCHITECTURE.md` §10.4: Henry/Raoult row —
      "⚠ gap" → "✅ closed"; Ksp row — unchanged, with a note that it remains
      a nested active-set solve, not a flat Jacobian, by accepted design
      (§14.3).
- [x] Update `CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md` status header and
      §18 ("planned rename") to "shipped."
- [x] Move this file to `docs/phases-shipped/`.
      (`EQUILIBRIUM_CONSTRAINT_UNIFICATION.md` already moved there at its own
      CP5 close, 2026-07-02 — it shipped and merged independently rather than
      waiting for this phase, per its own closing checkpoint.)
- [x] Tag `layer1-gap-closure-shipped`, merge `--no-ff` into `main`, delete
      the feature branch.

---

## Files changed

Paths below are as of when each checkpoint's changes actually landed
(`src/speciation/...`); CP7 subsequently renamed the whole module to
`src/chemical_equilibrium/...` — see that row.

| File | Nature of change |
|---|---|
| `src/speciation/nr_tableau.py` | CP1: gas-liquid graph folding (`ConfigurationError`, `SecondaryEntry.phase`, `_derive_gas_secondary`/`_derive_solvent_gas_secondary`); CP6: rename |
| `src/speciation/nr_solver.py` | CP2: volume-aware residual/Jacobian assembly, `SecondaryEntry.c_key` collision fix; CP6: rename |
| `src/speciation/engine.py`, `nr_engine.py` | CP2: `EquilibriumResult.partial_pressures_atm` population, `retain_jacobian` guard; CP3: `gas_liquid_species()`; CP6: rename |
| `src/thermo/liquid_phase_model.py`, `sit_liquid_model.py` | CP2: `DifferentiableLiquidModel` extension (Davies/SIT `∂γ_i/∂C_j`) |
| `src/core/control_volume.py` | CP3: `step_internal_transfer()` scope-filter fix (`_gas_liquid_engine_owned_species()`) |
| `models/vlmodels/adm1/base.py` | CP4: `WaterVapourBoundary` → `transfer_models=EquilibriumTransferModel(RaoultEquilibrium())` + liquid H2O seeding |
| `src/core/boundaries.py`, `src/core/__init__.py` | CP4: `WaterVapourBoundary`/`VentWaterLoss` deletion |
| `src/speciation/` → `src/chemical_equilibrium/` (whole module, all files) | CP7: module rename, all internal/external import paths updated |
| `tests/` | New tests per checkpoint (`test_nr_tableau_gas_liquid.py`, `test_nr_gas_liquid_cp2.py`, `test_step_internal_transfer_scope_filter.py`, `test_raoult_h2o_fold_cp4.py`, `test_precipitation_gas_liquid_cp5.py`); mechanical rename sweep CP6/CP7 |
| `docs/design/MASS_EXCHANGE_ARCHITECTURE.md`, `THERMODYNAMIC_MODEL_ARCHITECTURE.md`, `CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md` | CP7: status closeout; §8.6 water-activity follow-up (CP4 finding) |

---

## Trigger conditions

[EQUILIBRIUM_CONSTRAINT_UNIFICATION.md](../phases-shipped/EQUILIBRIUM_CONSTRAINT_UNIFICATION.md)
has shipped and tagged (2026-07-02) — start whenever you're ready. No other
external trigger; this is the agreed second half of Layer 1 gap closure.

## How to start

Follow the standard convention in [README.md](README.md): branch
`layer1-gap-closure` off `main`, work through CP1→CP7 in order, one commit
per checkpoint, full suite run at each. Tag and merge on completion.
