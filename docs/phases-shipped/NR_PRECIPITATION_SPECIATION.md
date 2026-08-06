# NR_PRECIPITATION_SPECIATION — Phase 1 of 2

> **Status:** Shipped 2026-06-23 — commit `16db5c7`, tag `nr-precipitation-speciation-shipped`.
> **Follows:** `nr-speciation-engine` (shipped 2026-06-23).
> **Precedes:** [NR_PRECIPITATION_CV_INTEGRATION.md](../phases-upcoming/NR_PRECIPITATION_CV_INTEGRATION.md) (Phase 2).

---

## Background and motivation

The NR speciation engine shipped in `nr-speciation-engine` solves dissolved
aqueous equilibria for arbitrary reaction networks.  The saturation-index
notebook (`05_saturation_index.ipynb`) demonstrates that the engine can flag
supersaturation (SI > 0) but cannot enforce it: the engine runs in
*dissolved-only* mode and returns concentrations as if precipitation were
suppressed.

The next natural step is **precipitation equilibrium**: when SI > 0 the solver
should find the split between dissolved and solid mineral such that IAP = Ksp.
This is required for:

- Accurate dissolved-Ca²⁺ and dissolved-carbonate predictions in hard water,
  lime-softening, and struvite systems.
- Starting simulations from a known solid feed (e.g. add CaCO₃ solid to a
  tank; determine how much dissolves).
- Providing a thermodynamically correct foundation for future kinetic
  crystallisation models.

This document covers **Phase 1: the speciation layer only**.  The engine's
`solve()` method gains an outer active-set precipitation loop; results are
returned in a `"minerals"` key on the output dict.  No changes to
`ControlVolume`, `SolidPhase`, or the ODE state occur here — those are Phase 2.

---

## Key design decisions (resolved in design discussion 2026-06-23)

### 1. No new class for mineral declaration

Minerals are declared as `EquilibriumReaction` objects with one
`StoichiometryEntry(phase="solid")` entry (the solid species) and one or more
`StoichiometryEntry(phase="liquid")` entries (the dissolved ionic products).
`ReactionSystem` routes any reaction containing a solid-phase entry into a new
`precipitation_equilibria` bucket.  No `MineralPhase` type is introduced.

Example — calcite dissolution:

```python
CaCO3 = Species(id="CaCO3", atoms={"Ca":1,"C":1,"O":3}, charge=0, MW=100.086)

calcite = EquilibriumReaction(
    stoichiometry=[
        StoichiometryEntry(species=CaCO3,         phase="solid",  coefficient=-1),
        StoichiometryEntry(species=Ca_plus_plus,  phase="liquid", coefficient=+1),
        StoichiometryEntry(species=CO3_2minus,    phase="liquid", coefficient=+1),
    ],
    log_K=-8.48,   # Ksp at 25 °C; a(CaCO3 solid) = 1 by convention
    label="calcite_dissolution",
)
```

The `log_K` here is the Ksp: because the solid is a pure phase its activity is
1, so the equilibrium expression contains only dissolved species.
`dH_J_per_mol` is supported for Van't Hoff correction (same as for dissolved
reactions).

### 2. No InertComponent — Ca²⁺ stays outside the inner NR system

Ca²⁺ (and other precipitating ions with no dissolved complexation) remains a
strong ion.  It does not become a master unknown in the NR tableau.  The outer
active-set loop adjusts its effective contribution to the charge balance:

```
CT_Ca_eff  = strong_ions["CT_Ca"]   −  ξ_calcite
CT_CO2_eff = totals["CO2"]          −  ξ_calcite
```

and passes these adjusted values to the unchanged inner `solve_nr()`.

If Ca²⁺ complexation reactions are declared in future (e.g. CaHCO₃⁺), Ca²⁺
would naturally enter the NR tableau as a master species via the BFS step
already in `build_tableau()` — the inner solver handles it without any new
apparatus, because the log-linear expansion of CaHCO₃⁺ in terms of
x_Ca + x_CO2 + x_H is structurally identical to how any multi-hop secondary
is handled today.  The only change required at that point is the
`element_stoichiometry` generalisation described in §4 below.

### 3. Outer active-set algorithm (MINTEQ/PHREEQC-style)

```
ξ = {} (empty — no active minerals)
active_set = {}

loop (outer active-set changes):
    compute effective totals / strong ions from ξ
    
    inner ξ-convergence loop (fixed active set):
        run solve_nr() with effective totals / strong ions
        compute f_j(ξ) = log(IAP_j) - log_Ksp_j  for each active mineral j
        if ||f||∞ < tol_prec:  break  (ξ converged)
        Newton step on ξ (finite-difference Jacobian)
        clamp each ξ_j ≥ 0 (physical constraint)

    active-set update:
        remove minerals with ξ_j < 0 → set ξ_j = 0
        add most-supersaturated inactive mineral if SI > 0

    if active set unchanged:  break (outer converged)
```

The finite-difference Jacobian for the outer Newton step costs 1 additional
`solve_nr()` call per active mineral per outer iteration.  For 1–3 minerals
this is negligible.

IAP is computed from the `solve_nr()` output concentrations with Davies
activity corrections applied:

- For dissolved species in the NR output (e.g. CO₃²⁻):
  `a_j = γ_j × out[species_id]`  where `γ_j = davies.gamma(z_j, I)`.
- For strong-ion precipitating species (e.g. Ca²⁺):
  `a_Ca = γ_Ca × CT_Ca_eff`  where `CT_Ca_eff = CT_Ca − ξ_calcite`.

Species charges are available from the `StoichiometryEntry.species.charge`
fields of the dissolution reaction.

Van't Hoff correction to log_Ksp is applied when `rxn.dH_J_per_mol` is set,
using the same `_vant_hoff_log_K` helper already in `nr_tableau.py`.

### 4. element_stoichiometry fix (cross-component mass balances)

The current `_residual_and_jacobian` assigns each secondary species to exactly
one component and only counts it toward that component's mass balance.  This is
correct for single-component systems but wrong for cross-component species like
CaHCO₃⁺, which contains both Ca and C and should contribute to both the Ca
and the TIC mass balances.

The fix: add `element_stoichiometry: Dict[str, float]` to `SecondaryEntry`.
This maps each **non-H⁺ master ID** to the count of that master's "formula
unit" in the secondary species.  For simple secondaries in existing systems,
`element_stoichiometry` equals `{master_id: nu[master_id]}` — so the fix is
backward-compatible.

The accumulation loop in `_residual_and_jacobian` is changed from:

```python
# OLD: secondary counted only for the component it was assigned to
for row_idx, comp in enumerate(components):
    for sp_id in comp.species_ids:
        R[row_idx] += c[sp_id]
```

to:

```python
# NEW: secondary counted for every component according to element stoichiometry
for sec in secondaries:
    c_j = concentrations[sec.species_id]
    for row_idx, comp in enumerate(components):
        coeff = sec.element_stoichiometry.get(comp.master_id, 0.0)
        if coeff != 0.0:
            R[row_idx] += coeff * c_j
```

`element_stoichiometry` is computed in `build_tableau()` during the BFS
derivation step: for a secondary species with `nu = {m1: ν1, m2: ν2, "H+": νH}`,
the element stoichiometry is `{m1: ν1, m2: ν2}` (H⁺ is excluded because H⁺
closes the charge balance, not a mass balance).

This is the standard representation used in PHREEQC and all major geochemical
codes: each species carries a "formula vector" decomposing it in terms of
master species.  The VLsim component-grouping shortcut is a special case of
this general structure.

---

## Scope

**In scope:**

- `EquilibriumReaction` with `phase="solid"` entries recognised and routed to
  `ReactionSystem.precipitation_equilibria`.
- `element_stoichiometry` field on `SecondaryEntry`; mass balance accumulation
  generalised.
- Outer active-set loop and finite-difference outer-Newton in
  `NRSpeciationEngine.solve()`.
- `solve()` output dict gains `"minerals"` key:
  `{"calcite": {"xi_mol_L": float, "SI": float}}`.
- `Ca_plus_plus` and `Mg_plus_plus` added to `common_species.py`.
- Tests: `TestPrecipitationEquilibria` class in
  `tests/standalone/test_nr_speciation_engine.py`.
- Demo notebook `05_precipitation_equilibrium.ipynb`.
- README table updated.
- Docstring on `NRSpeciationEngine` and `NRSpeciationEngine.solve()` reference
  Phase 2 for CV/SolidPhase integration.

**Out of scope (deferred to Phase 2):**

- `SolidPhase` writeback from ξ.
- `ControlVolume._read_from_phases` summing liquid + solid contributions.
- Lazy-creation / consistency of phase types (GasPhase, LiquidPhase,
  SolidPhase) on ControlVolume.
- Multiple co-precipitating minerals sharing a common component (e.g. calcite +
  aragonite both consuming CT_CO2) — tested only with independent minerals in
  Phase 1.
- Kinetic crystallisation rates (separate future phase, if needed).
- Redox equilibria (pe as a master) — deferred per NR design doc.

---

## Files changed

| File | Nature of change |
|---|---|
| `src/chemistry/common_species.py` | Add `Ca_plus_plus`, `Mg_plus_plus` |
| `src/reactions/reaction_system.py` | Add `precipitation_equilibria` bucket; populate when any stoichiometry entry has `phase="solid"` |
| `src/speciation/nr_tableau.py` | Add `element_stoichiometry` to `SecondaryEntry`; compute it in `build_tableau()` |
| `src/speciation/nr_solver.py` | Generalise mass balance accumulation in `_residual_and_jacobian` to use `element_stoichiometry` |
| `src/speciation/nr_engine.py` | Outer active-set loop; `from_reactions()` accepts `precipitation_reactions`; `solve()` returns `"minerals"` key; docstring Phase 2 reference |
| `tests/standalone/test_nr_speciation_engine.py` | New `TestPrecipitationEquilibria` class; new test for `element_stoichiometry` correctness |
| `demos/model_api/chemistry/speciation/_generate_notebooks.py` | Add `05_precipitation_equilibrium.ipynb`; update README table |

---

## Test cases (TestPrecipitationEquilibria)

1. **`test_calcite_precipitation_from_supersaturated`** — CT_Ca=0.002,
   CT_CO2=0.010, CT_Na=0.005 (known SI > 0 from notebook 03); after solve
   `xi > 0` and `log(a_Ca × a_CO₃) ≈ −8.48 ± 0.01`.

2. **`test_calcite_mass_balance`** — total Ca conserved:
   `out["minerals"]["calcite"]["xi_mol_L"] + CT_Ca_dissolved ≈ CT_Ca_input`.

3. **`test_undersaturated_no_precipitation`** — low CT_Ca + CT_CO2 at low pH;
   SI < 0; `xi == 0.0`.

4. **`test_dissolution_from_solid`** — declare n_mol_CaCO3_solid > 0, no Ca in
   liquid; after solve some dissolves; IAP ≈ Ksp; dissolved_Ca + xi ≈ initial
   solid.

5. **`test_complete_dissolution_high_ksp`** — NH₄Cl-like dissolution reaction
   with log_Ksp = +5; starting solid amount n_mol_solid = 0.01; after solve
   xi ≈ n_mol_solid (all dissolves).

6. **`test_element_stoichiometry_computed`** — build tableau for a
   cross-component reaction (Ca²⁺ + HCO₃⁻ → CaHCO₃⁺); verify
   `sec.element_stoichiometry == {"CO2": 1.0, "Ca++": 1.0}`.

7. **`test_element_stoichiometry_single_component_unchanged`** — existing
   carbonate + ammonia system; verify backward-compatible element_stoichiometry
   (HCO₃⁻: `{"CO2": 1.0}`, NH₄⁺: `{"NH3": 1.0}`).

---

## Demo notebook: 05_precipitation_equilibrium.ipynb

Sections:

1. **Background** — equilibrium precipitation vs SI assessment; Ksp as the
   governing constraint; comparison with the dissolved-only (SI > 0) result
   from notebook 03.
2. **Declaration** — `EquilibriumReaction` with `phase="solid"` entry; Ca²⁺
   as a strong ion with `strong_ions={"CT_Ca": ...}`.
3. **Single-point solve** — CT_Ca=0.002, CT_CO2=0.010, CT_Na=0.005; compare
   dissolved [Ca²⁺] and [CO₃²⁻] with and without precipitation; verify
   IAP = Ksp after precipitation solve.
4. **Dissolved Ca²⁺ depletion curve** — sweep CT_Ca at fixed CT_CO2 and pH;
   show xi as a function of CT_Ca; show [Ca²⁺]_dissolved saturating at Ksp-set
   ceiling.
5. **Dissolution from solid** — start with CaCO₃ solid = 0.001 mol/L, no Ca
   in liquid; sweep CT_CO2 from 0 to 0.020 mol/L; show fraction dissolved as a
   function of CT_CO2.
6. **Scope note** — dissolved-only assumption for Ca²⁺ (no complexation);
   pointer to Phase 2 for CV/SolidPhase integration.

---

## How to start

1. Create branch `nr-precipitation-speciation` off `main`.
2. Work through the files in the order listed in *Files changed*.
3. Run `pytest tests/standalone/test_nr_speciation_engine.py -v` after each
   file to catch regressions early.
4. Generate and execute notebooks with
   `python demos/model_api/chemistry/speciation/_generate_notebooks.py`
   followed by `jupyter nbconvert --execute`.
5. When green: merge `--no-ff`, tag `nr-precipitation-speciation-shipped`,
   delete branch.
