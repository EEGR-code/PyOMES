# Polyprotic Acid Chains in SpeciationEngine

**Date:** 2026-06-21  
**Status:** Shipped (1502/0 tests)

---

## Background

`EquilibriumReaction` models a single deprotonation step: one reactant acid,
one conjugate base, one H⁺.  A three-pKa system like phosphate therefore
requires three separate objects:

```python
# H3PO4 ⇌ H2PO4- + H+  (pKa 2.15)
# H2PO4- ⇌ HPO4-- + H+  (pKa 7.20)
# HPO4-- ⇌ PO4--- + H+  (pKa 12.35)
```

Before this fix, `SpeciationEngine.from_reactions` registered each step as an
**independent monoprotic acid** in the `EquilibriumSet`, giving three separate
`EquilibriumDef` entries.  This caused two bugs:

1. **Species overwriting** — `_compute_species_eq` iterates the entries in
   order and assigns `out[species_id] = CT * alpha`.  The HPO4-- step had
   `CT = 0` (because `n_mol["HPO4--"]` starts at zero), so it overwrote the
   correct HPO4-- value that had just been written by the H2PO4- step with
   zero.

2. **Progressive phosphate drain** — because `n_mol["HPO4--"]` was always
   forced to zero, the H2PO4- step used only `n_mol["H2PO4-"]` as its total,
   which was already reduced by the previous write-back, causing H2PO4- to
   decrease each step and total phosphate to drift away from its conserved
   value.

Together these produced spurious `ConservationWarning` on both the
`charge_cum` and `element_step_H` checks after `cv.equilibrate_to_pH()`.

---

## Fix

### Convention: `total_id` marks polyprotic steps

Set `total_id="<ladder_name>"` identically on every step that belongs to the
same acid:

```python
# bioprocess_basic.py
EquilibriumReaction(..., log_K=-2.15, total_id="phosphate", label="eq_phosphate_1")
EquilibriumReaction(..., log_K=-7.20, total_id="phosphate", label="eq_phosphate_2")
EquilibriumReaction(..., log_K=-12.35, total_id="phosphate", label="eq_phosphate_3")
```

### Engine: merge same-`total_id` groups

`from_reactions` now does a two-pass approach:

1. Classify each reaction and compute its resolved name
   (`total_id` if set, else the acid species id).
2. Group by resolved name.  Groups with a single reaction are handled by the
   existing `_add_acid` path.  Groups with multiple reactions are handed to
   the new `_add_polyprotic_acid` helper.

`_add_polyprotic_acid` sorts the group by ascending pKa, assembles the full
species chain from most-protonated to least-protonated (H3PO4 → H2PO4- →
HPO4-- → PO4---), and calls `eq_set.add(name, pKas=(2.15, 7.20, 12.35),
species_refs=(...))` **once**.  The resulting single `EquilibriumDef` is
identical to what `EquilibriumSet.bsm2_default()` uses for phosphate, so the
solver's polyprotic alpha computation conserves total phosphate by
construction.

### Result

After the fix, `_compute_species_eq` processes phosphate as one entry and
emits all four species (H3PO4, H2PO4-, HPO4--, PO4---) from a single
CT-conserving alpha expansion.  No overwriting, no drain.  Between
consecutive `advance()` steps the total phosphate in `n_mol` is stable to
machine precision and the charge residual is ≈ 1 × 10⁻⁸ mol (threshold
1 × 10⁻⁶ mol).

---

## Files changed

| File | Change |
|---|---|
| `src/chemistry/databases/bioprocess_basic.py` | `total_id="phosphate"` on all three phosphate `EquilibriumReaction` objects |
| `src/speciation/engine.py` | `from_reactions` groups by resolved name; new `_add_polyprotic_acid()` helper |
| `src/speciation/engine.py` | `_refresh_derived` writeback list extended to include H3PO4, H2PO4-, HPO4--, PO4--- |
| `src/core/control_volume.py` | `_collect_species_registry()` supplemental scan of `common_species` catalog (catches Cl-, Na+, K+, H2PO4-, HPO4--, PO4---) |
| `src/reactions/builder.py` | `monod_aerobic_growth()` — optional `Ko2_gL` O₂ Monod term |
| `tests/standalone/test_cv_advance.py` | `test_advance_after_equilibrate_no_conservation_warning` regression guard |

---

## Rule for new acid systems

Any acid that requires more than one deprotonation step and is defined using
individual `EquilibriumReaction` objects **must** carry a shared `total_id`
on every step.  Omitting `total_id` is correct only for genuinely independent
monoprotic buffers.

Single-step acids (H2S/HS-, HSO4-/SO4--, organic VFAs) do not need
`total_id`; they naturally produce a single `EquilibriumDef`.
