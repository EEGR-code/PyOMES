# Implementation Plan: v0.12.6 (Pre-UAT) and v0.13.0 (Post-UAT)

**Codebase:** fermenter v0.12.5  
**Working directory:** `fermenter_codebase_v12_5/`  
**Current state:** 690 tests, 7 warnings  
**Date:** 2026-03-24  

This document captures architectural issues found during deep code review and
the implementation plan for addressing them. Items G.1–G.6 are pre-UAT
protective fixes (warnings and freezing). Items H.1–H.5 are post-UAT
structural fixes that require decomposing `fermenter_unit.py`.

---

## Context: Issues Found

### Configuration vs Reality Mismatches

| # | Issue | Severity | Root Cause |
|---|-------|----------|------------|
| 1 | Speciation engine `T_C` frozen at construction — pKa temperature corrections use stale temperature if `self.T` is changed | HIGH | `SpeciationEngine` stores `T_C` once; no setter or per-call override in the fermenter's direct speciation path |
| 2 | `V_liq_L` frozen at construction — never updated during simulation despite pH dosing adding liquid volume | HIGH | `self.V_liq_L` set once at line 709; no dosing volume tracking |
| 3 | HPLC derived quantities (`_u_cm_min`, `_dz`, K_D_eff) frozen at construction — mutating `flow_rate_mL_min` etc. after construction gives stale simulation | MEDIUM | Plain attributes with no setter protection; derived quantities computed once |
| 4 | `acid_pKas` is a mutable dict shared by reference between fermenter and all speciation calls — silent mutation possible | MEDIUM | `self.acid_pKas = dict(acid_pKas)` shared; `SpeciationPropertySolver` captures a separate copy |
| 5 | `strong_kwargs` accumulated during a run but not persisted between runs — user cannot inspect post-run strong ion state | MEDIUM | Local dict in `_calc_ODE_with_headspace`, not stored on `self` |

### Silent Failure Points

| # | Issue | Severity | Root Cause |
|---|-------|----------|------------|
| 6 | Module-level import of `speciate_multi_acids_TIC` in `fermenter_unit.py` — unused but creates hard dependency on legacy adapter | MEDIUM | Line 28 compatibility import |
| 7 | HPLC species dataclasses not frozen — `sp.K_D = 0.5` after column construction silently corrupts state | LOW | `@dataclass` without `frozen=True` |
| 8 | `freeze_speciation=True` in `ScipyODESolver` has no warning — user may not realise pH is frozen within ODE steps | LOW | No diagnostic at construction |

### Architectural Inconsistencies

| # | Issue | Severity | Root Cause |
|---|-------|----------|------------|
| 9 | Two parallel speciation paths (direct at t=0, CV-based for t>0) with different capabilities — initial pH may be inconsistent with ODE loop pH | HIGH | Path A (line 1686) lacks `dH_` context that Path B (CV advance) receives |
| 10 | Three separate sources of truth for Henry constants (`self.kH_CO2`, `KineticGasLiquidLink.henry`, `HenryEquilibriumInterface`) | MEDIUM | Historical layering of pre-CV and CV-based code paths |
| 11 | `self.T` is a bare attribute with no change propagation — changing it updates nothing except Henry constants (on next `_update_henry_constants()` call) | MEDIUM | No property/setter; no observer pattern |

---

## Pre-UAT Items (v0.12.6)

These items prevent silent incorrect results, are self-contained with no
regression risk, and add warnings that make existing limitations visible.
None require restructuring `fermenter_unit.py`.

### G.1: Freeze HPLC column against post-construction mutation (45 min)
*Addresses issues #3, #7*

**Files:** `src/fermenter/models/hplc_column.py`, `tests/standalone/test_hplc_column.py`

**G.1a: Freeze species dataclasses.**
Change `@dataclass` to `@dataclass(frozen=True)` on `LangmuirSpecies` and
`PartitionSpecies`. The column constructor already copies the species list,
so no downstream code should be mutating species objects. Any code that does
will get a clear `FrozenInstanceError`.

**G.1b: Make species list a tuple.**
Change `self.species = list(species or [])` to `self.species = tuple(species)`
(after the empty-list validation). All iteration and indexing works identically
on tuples.

**G.1c: Freeze physical parameters via `__setattr__`.**
Add a `_FROZEN_ATTRS` frozenset and override `__setattr__` to reject changes
to physical parameters after construction:

```python
_FROZEN_ATTRS = frozenset({
    "length_cm", "diameter_cm", "void_fraction", "flow_rate_mL_min",
    "bulk_density_g_mL", "D_ax_default", "T_K", "n_cells",
    "mobile_phase_pH", "K_D_reference_pH", "advection_scheme",
    "dispersion_scheme", "intraparticle_porosity",
})

def __setattr__(self, name, value):
    if getattr(self, "_frozen", False) and name in self._FROZEN_ATTRS:
        raise AttributeError(
            f"Cannot modify '{name}' after construction. "
            f"Create a new HPLCColumn with the desired parameters.")
    super().__setattr__(name, value)
```

Set `self._frozen = True` as the last line of `__init__`.

**Tests to add:**
- `test_species_frozen` — `sp.K_D = 0.5` raises `FrozenInstanceError`
- `test_column_attribute_frozen` — `col.flow_rate_mL_min = 1.0` raises `AttributeError`
- `test_non_frozen_attribute_allowed` — `col.label = "new"` succeeds
- `test_species_list_is_tuple` — `isinstance(col.species, tuple)`

---

### G.2: Warn on speciation engine temperature drift (30 min)
*Addresses issue #1 — mitigation only*

**Files:** `src/fermenter/fermenter_unit.py`

**G.2a:** Record construction temperature:
```python
self._construction_T_K = float(T)
```
(in `__init__`, after `self.T = float(T)`)

**G.2b:** At start of `_calc_ODE_with_headspace`, check for drift:
```python
if abs(float(self.T) - self._construction_T_K) > 0.1:
    warnings.warn(
        f"Reactor temperature ({self.T:.1f} K) differs from construction "
        f"({self._construction_T_K:.1f} K). Henry constants will be "
        f"recalculated, but speciation pKa temperature corrections still "
        f"use the original temperature. For accurate temperature-dependent "
        f"speciation, reconstruct the fermenter or use the ADM1/BSM2 path "
        f"with ThermodynamicConfig.",
        RuntimeWarning, stacklevel=2,
    )
```

**Tests to add:**
- `test_temperature_drift_warns` — change `self.T`, call simulate, assert `RuntimeWarning`
- `test_no_warn_at_construction_temperature` — no warning when T unchanged

---

### G.3: Document V_liq_L limitation and add dosing volume diagnostic (30 min)
*Addresses issue #2 — mitigation only*

**Files:** `src/fermenter/fermenter_unit.py`

**G.3a:** Add docstring note to `_calc_ODE_with_headspace`:
```
Note: Liquid volume (V_liq_L) is assumed constant during the simulation.
This is physically correct for batch fermentation and scenarios where
dosing volumes are negligible relative to total liquid volume.
```

**G.3b:** At end of `_calc_ODE_with_headspace`, after the ODE loop, estimate
total dosed volume and warn if >1% of liquid volume:

```python
# Estimate dosed liquid volume from ledger
dosed_mol = self._ledger.dosed_total_mol if hasattr(self._ledger, 'dosed_total_mol') else {}
total_dosed_mol = sum(abs(v) for v in dosed_mol.values())
# Rough estimate: 1 mol dosing solution ≈ 0.1 L (concentrated acid/base)
est_dosed_L = total_dosed_mol * 0.1
if est_dosed_L > 0.01 * self.V_liq_L:
    warnings.warn(
        f"Estimated dosing volume (~{est_dosed_L:.3f} L) exceeds 1% of "
        f"liquid volume ({self.V_liq_L:.1f} L). Concentrations may be "
        f"slightly inaccurate because V_liq_L is held constant during "
        f"simulation.",
        RuntimeWarning,
    )
```

**Tests to add:**
- Verify warning fires with large dosing (may need a targeted integration test)

---

### G.4: Remove unused module-level import (5 min)
*Addresses issue #6*

**File:** `src/fermenter/fermenter_unit.py`

Check if `speciate_multi_acids_TIC` is referenced anywhere in the file beyond
the import line. If not, remove line 28:
```python
from .speciation.legacy_adapter import speciate_multi_acids_TIC  # legacy solver exposed for compatibility
```

The function is already re-exported from `fermenter.speciation.__init__` —
that is the correct import path. Verify no tests import it from
`fermenter.fermenter_unit`.

---

### G.5: Add `freeze_speciation` warning (10 min)
*Addresses issue #8*

**File:** `src/fermenter/core/solvers.py`

In `ScipyODESolver.__init__`, when `freeze_speciation=True`:
```python
if freeze_speciation:
    import warnings
    warnings.warn(
        "freeze_speciation=True: speciation will be computed once per "
        "ODE step and reused for all internal substeps. This improves "
        "performance but may reduce accuracy for pH-sensitive systems.",
        UserWarning, stacklevel=2,
    )
```

**Tests to add:**
- `test_freeze_speciation_warns` — assert `UserWarning` emitted at construction

---

### G.6: Update CHANGELOG and run full test suite (30 min)

Add v0.12.6 entry covering G.1–G.5. Verify 690+ tests pass.

---

## Post-UAT Items (v0.13.0 / P3.3)

These items are interconnected and stem from the same root cause:
`fermenter_unit.py` maintains parallel state (plain attributes, a
SpeciationEngine, a ControlVolume, and local dicts) rather than deriving
everything from a single authoritative ControlVolume. They should be planned
as a coordinated refactoring sprint after UAT feedback confirms which
simulation paths users exercise most.

### Dependency Graph

```
H.3 (unified speciation) ──┐
                            ├── H.4 (Henry consolidation)
H.2 (dynamic V_liq_L) ─────┤
                            ├── H.5 (pKa/strong ion ownership)
H.1 (speciation T) ────────┘
        │
        └── All depend on P3.3 (CV as persistent state owner)
```

### H.1: Unify speciation temperature with phase temperature (2–3 hours)
*Addresses issues #1, #9, #11*

**Approach:** Add optional `T_K` kwarg to `SpeciationEngine.solve()`. If
provided, override the model's `T_C` for that call. The
`SpeciationPropertySolver` already reads T_K from the liquid phase — extend
it to pass `T_K` through to the engine. In the fermenter's direct speciation
path (line 1686), also pass `T_K=float(self.T)`.

**Key insight:** The Level 1 and Level 2 models already compute
`T_K = 273.15 + self.T_C` inside `solve()`. The change is making this
overridable per call rather than frozen at construction.

**Dependency:** Should be done alongside H.3 to avoid fixing the direct
path only to have it removed.

---

### H.2: Make `V_liq_L` dynamic (3–4 hours)
*Addresses issue #2*

**Approach:** Replace `self.V_liq_L` with a property that reads from the
CV's liquid phase. Dosing logic updates the CV's liquid volume when liquid
is added. This requires the CV to exist for the entire simulation lifetime,
not just the ODE loop — which is the core of P3.3.

```python
@property
def V_liq_L(self):
    if hasattr(self, '_cv') and self._cv is not None:
        liq = self._cv.phases.get("liquid")
        if liq is not None:
            return float(liq.V_L)
    return self._V_liq_L_init
```

**Dependency:** Requires P3.3 (CV as persistent state owner).

---

### H.3: Unify the two speciation paths (2–3 hours)
*Addresses issues #9, #10*

**Approach:** Eliminate Path A (direct `self.speciation.solve()` at t=0).
Construct the CV before the initial speciation. Initial speciation becomes:

```python
_cv = ControlVolume(...)
_advance_result = _cv.advance(dt_h=0, context=initial_context)
```

This ensures temperature corrections, pKa sources, and context kwargs are
identical for all timesteps including t=0.

**Dependency:** Requires P3.3 (CV lifetime management).

---

### H.4: Consolidate Henry constant sources (1–2 hours)
*Addresses issue #10*

**Approach:** Remove `self.kH_CO2`, `self.henry_CO2_mol_per_L_atm` etc.
from the fermenter object. All Henry constant lookups go through the
`KineticGasLiquidLink.henry` dict (already temperature-corrected per
timestep). The fermenter's `_update_henry_constants()` method becomes
unnecessary once the CV path is the only path.

**Dependency:** Requires H.3.

---

### H.5: `acid_pKas` ownership and `strong_kwargs` persistence (1 hour)
*Addresses issues #4, #5*

**Approach:**
- `acid_pKas` should be owned by `ThermodynamicConfig` (ADM1) or
  `SpeciationPropertySolver` (fermenter), not as a mutable dict on the
  fermenter object.
- `strong_kwargs` persistence: if inter-run accumulation is desired, store
  strong ion state on the CV's liquid phase. If non-persistence is correct,
  document explicitly.

**Dependency:** Benefits from P3.3 but can be done independently.

---

## Summary Table

| Item | Effort | Timing | Issues | Risk |
|------|--------|--------|--------|------|
| G.1 | 45 min | Pre-UAT | #3, #7 | None — HPLC is self-contained |
| G.2 | 30 min | Pre-UAT | #1 | None — warning only |
| G.3 | 30 min | Pre-UAT | #2 | None — warning only |
| G.4 | 5 min | Pre-UAT | #6 | Minimal — verify no external users |
| G.5 | 10 min | Pre-UAT | #8 | None — warning only |
| G.6 | 30 min | Pre-UAT | — | Test + changelog |
| H.1 | 2–3 hrs | Post-UAT | #1, #9, #11 | Medium — touches engine API |
| H.2 | 3–4 hrs | Post-UAT | #2 | High — requires P3.3 |
| H.3 | 2–3 hrs | Post-UAT | #9, #10 | High — requires P3.3 |
| H.4 | 1–2 hrs | Post-UAT | #10 | Medium — requires H.3 |
| H.5 | 1 hr | Post-UAT | #4, #5 | Low — mostly documentation |

Pre-UAT total: ~2.5 hours  
Post-UAT total: ~10–13 hours (overlaps with P3.3)

---

## Completed Work (for reference)

| Version | What was done | Tests |
|---------|--------------|-------|
| v0.12.3 | HPLC model, numerics, silent failure fixes, dH corrections, print→logging, pyproject.toml, README | 584 |
| v0.12.4 | legacy_adapter/biosteam_wrapper clarification, solvers/ cleanup, `__version__`, CHANGELOG | 653 |
| v0.12.5 | Phase invariants on core objects, HPLC validation, live V_liq reads in ADM1, config warnings, equilibria exception narrowing, class docstrings, redundant fallback removal | 690 |
