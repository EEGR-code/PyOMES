# EQUILIBRIUM_RESULT

> **Status:** Shipped 2026-07-01. Tag `equilibrium-result-shipped`, merged
>   to `main`. All CP0-CP7 checkpoints complete. 1793 → 1800 tests.
> **Depends on:** [THERMODYNAMIC_MODEL_ARCHITECTURE.md](THERMODYNAMIC_MODEL_ARCHITECTURE.md)
>   shipped 2026-07-01. Design authority:
>   [CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md](../design/CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md)
>   §6 (`EquilibriumResult`), §7 (z-update strategies), §8 (design decisions),
>   §18 (planned rename — out of scope here, tracked separately).
> **Companion branch:** `legacy-cleanup`, tag `legacy-cleanup-shipped`
>   (CP0's dead-code deletion, merged to `main` separately before this
>   branch's own work began).

---

## Background and motivation

`SpeciationEngine.solve()` and `NRSpeciationEngine.solve()` currently return a
mutable `dict` and, when called with `phases=`, write derived species back to
`phase.n_mol` as an inseparable side effect of the same call (via
`_refresh_derived`). Every production call site that passes `phases=` discards
the return value entirely and relies solely on that side effect.

This is incorrect inside adaptive/implicit ODE solvers: the solver evaluates
`compute_rhs()` / the RHS closure multiple times per accepted step, including
on rejected trial steps, and each call re-triggers the write. Separating
*solve* (pure, returns an immutable result) from *commit* (`apply_to_phases()`,
called explicitly only when a step is accepted) removes this inconsistency and
unblocks the PHREEQC solver-mode phase (§5.2, §17 Phase 6 in the architecture
doc), which depends on `solve() → EquilibriumResult` being in place first.

Investigating the current call graph for this phase also surfaced dead code
left over from `CUFERMENTER_SUNSET` (`fermenter_unit.py` deleted, but
`src/equilibria/engine.py`, `src/equilibria/coupled.py`, and
`src/core/gl_equilibrium.py`'s `HenryEquilibriumInterface` were not). CP0
audits this before the main migration touches any of it.

---

## Key design decisions

Decisions recorded this session (2026-07-01), superseding/resolving the "not
yet scoped" note in the architecture doc's header for the `solve()` return
type:

**No `Mapping` shim — clean break to attribute-only access.** Confirmed
despite the migration touching more call sites than §8 originally scoped
(`numerical_gradient.py`, `speciation/api.py`, ~20 test files, in addition to
`cv.advance()`/`ReactionSystem`). One deliberate, narrow concession:
`EquilibriumResult.to_dict()` is added as an explicit opt-in convenience for
call sites that need a dict snapshot (e.g. `speciation/api.py`'s
`raw=` field) — this is not a `Mapping` protocol implementation and does not
support `result["pH"]` / `result.get(...)`.

**Atomic cutover — `solve()` stops auto-writing entirely, in the same phase.**
All 7 production call sites that currently rely on `solve(phases=...)`
auto-writing (`control_volume.py` ×2, `solvers.py` ×3, `simulation.py` ×1) are
migrated to the explicit `result = engine.solve(...); result.apply_to_phases(...)`
pattern within this phase, not deferred. A partial migration would silently
break pH/derived-species updates on whichever call sites were missed.

**Field inventory (from tracing `acid_base.py` and `nr_solver.py` output
keys) — named fields, not a generic `extra` dump, for the common core:**
`pH`, `pH_conc`, `logH`, `aH`, `gamma_H`, `gamma_OH`, `ionic_strength`,
`charge_residual`, `n_iter` (verify whether ever populated — treat as
always-`None` if not), `species_mol_L`, `partial_pressures_atm`,
`saturation_indices`, `alphas`, `extra` (catch-all).

**`species_mol_L` split rule differs by engine — must preserve today's exact
writeback behavior:**
- `SpeciationEngine`: only the fixed canonical species tuple already used by
  the current writeback block (engine.py:420-430) goes into `species_mol_L`.
  Generic polyprotic-ladder keys (`{name}_HA`, `{name}_A-`, `{name}_BH+`,
  `{name}_B`) are real dict keys today but were never written back — they go
  into `extra`, not `species_mol_L`, so `apply_to_phases()` doesn't change
  writeback behavior.
- `NRSpeciationEngine`: no split needed — `species_mol_L` is the full
  masters + secondaries + `H2O` set, matching today's `_writeback` exactly.
- `PHREEQCEngine`: `species_mol_L` is every `out` key except the `{pH, logH,
  IonicStrength}` meta keys.

**Precipitation output:** mineral `SI` → `saturation_indices`. Mineral `ξ`
(mol/L precipitated) has no phase-writeback target yet (deferred to
[NR_PRECIPITATION_CV_INTEGRATION.md](NR_PRECIPITATION_CV_INTEGRATION.md)) —
goes into `extra["minerals_xi_mol_L"]`, diagnostic only.

---

## Scope

**In scope:**

- `EquilibriumResult` frozen dataclass + `apply_to_phases()` in
  `src/speciation/protocols.py`.
- `SpeciationEngine`, `NRSpeciationEngine`, `PHREEQCEngine` all construct and
  return `EquilibriumResult` from `solve()`.
- `numerical_gradient.py`, `speciation/api.py` updated to consume
  `EquilibriumResult` (conditional on CP0 not deleting `speciation/api.py`
  outright).
- All 7 production call sites cut over to explicit
  `solve()` + `apply_to_phases()`.
- Mechanical test-suite migration off dict-style access.
- CP0 audit (and, if approved, deletion) of dead code found during this
  phase's investigation.

**Out of scope (tracked elsewhere):**

- §18 rename (`SpeciationEngine` → `ChemicalEquilibriumEngine` etc.) —
  accompanies Layer 1 gap closure, not this phase.
- PHREEQC solver mode (§5.2, architecture doc §17 Phase 6) — depends on this
  phase shipping first, not part of it.
- `NRChemicalEquilibriumEngine` white-box `SplitJacobianCapable`/
  `ResidualCapable` completeness — already shipped (Phase 4, prior session);
  this phase only needs to confirm those methods don't read the now-changed
  `solve()` return value.
- Mineral `ξ` writeback to `SolidPhase.n_mol` — `NR_PRECIPITATION_CV_INTEGRATION.md`,
  separate phase, not yet started.

---

## Checkpoints

### CP0 — Legacy dead-code audit (CUFERMENTER_SUNSET follow-up) — SHIPPED 2026-07-01

- [x] Enumerate the full `src/equilibria/` package (`engine.py`, `coupled.py`,
      `vle.py`, `interfaces.py`, `__init__.py` if present); classify each
      file/class live/dead by tracing the full import closure from
      `src/core/__init__.py` and every test/demo entry point.
      `vle.py`/`peng_robinson.py` are live (used by `chemistry/partition.py`,
      `thermo/framework.py`); `engine.py`/`coupled.py`/`factory.py`/
      `interfaces.py` were dead.
- [x] Trace `src/core/gl_equilibrium.py` (`HenryEquilibriumInterface`):
      confirmed zero live constructors/callers across `src/`, `tests/`,
      `demos/` beyond the `core/__init__.py` re-export; it already raised a
      `DeprecationWarning` at `__init__` pointing to `EquilibriumTransferModel`
      (added when TRANSFER_MODEL shipped 2026-06-19).
- [x] Trace `src/speciation/api.py` (`SpeciationEngineAdapter`): confirmed
      zero live callers, but self-contained and not self-deprecated —
      classified "keep."
- [x] Extended the search to `demos/` — no hits for any candidate.
- [x] Checked `docs/` — none of the candidates are documented as public API
      in TRANSFER_MODEL.md or CUFERMENTER_SUNSET.md.
- [x] Cross-checked against the SIMULATION_CLASS-phase "Legacy
      CUFermentationSpeciation island" memory (`run_batch`, `loops.py`,
      `system.py`, `types.py`) — that island's files (`loops.py`, `system.py`,
      `types.py`) no longer exist in the repo at all (removed in a later,
      undocumented cleanup); `equilibria/engine.py`/`coupled.py` are a
      **separate, unrelated leftover**, not part of that island.
      `factory.py::CoupledEquilibriumFactory.create_process_for_fermenter`
      imports `vlmodels.fermenter.unit` (confirmed via repo-wide search: does
      not exist anywhere) and expects a `Fermenter`-shaped object deleted in
      CUFERMENTER_SUNSET — direct evidence of unreachability, not just absence
      of callers.
- [x] Classification table produced and presented (see chat transcript
      2026-07-01): `engine.py`/`coupled.py`/`factory.py`/`interfaces.py`/
      `gl_equilibrium.py` → dead, delete now. `speciation/api.py` → unused but
      not broken, keep (feeds CP5). `vle.py`/`peng_robinson.py` → live, not a
      candidate.
- [x] User approved "Delete now" via AskUserQuestion.
- [x] Deleted on branch `legacy-cleanup` (commit `ee26f40`); trimmed
      `src/core/__init__.py` (import + `__all__` entry) and
      `src/equilibria/__init__.py` (re-exports → docstring-only). Merged
      `--no-ff` into `main` (`f644c80`), tagged `legacy-cleanup-shipped`,
      branch deleted.
- [x] Full suite run after deletion: 1793 passed, 27 skipped, 0 failures —
      unchanged from baseline, confirming zero test coverage depended on the
      deleted files.
- [x] `speciation/api.py` verdict ("keep") fed forward into CP5.

### CP1 — `EquilibriumResult` + `apply_to_phases()` — SHIPPED 2026-07-01 (commit `2045503`)

- [x] Add frozen dataclass `EquilibriumResult` to `src/speciation/protocols.py`
      with fields: `pH`, `pH_conc`, `logH`, `aH`, `gamma_H`, `gamma_OH`,
      `ionic_strength`, `charge_residual`, `n_iter: Optional[int] = None`,
      `species_mol_L: Dict[str, float]`,
      `partial_pressures_atm: Dict[str, float] = field(default_factory=dict)`,
      `saturation_indices: Dict[str, float] = field(default_factory=dict)`,
      `alphas: Dict[str, float] = field(default_factory=dict)`,
      `extra: Dict[str, Any] = field(default_factory=dict)`.
- [x] Add `apply_to_phases(self, phases, liquid_key="liquid") -> None`: guard
      `liq = phases.get(liquid_key)`, `V_L = getattr(liq, "V_L", 0.0)`, skip
      if `V_L <= 0` or no `_refresh_derived`; build
      `{sp: conc * V_L for sp, conc in species_mol_L.items()}` and call
      `liq._refresh_derived(...)`.
- [x] Add `to_dict()` convenience method (`dataclasses.asdict(self)`) — the
      one deliberate concession to dict-shaped consumers (`speciation/api.py`).
- [x] Update `SpeciationEngineProtocol.solve()` return annotation from
      `Dict[str, Any]` to `"EquilibriumResult"` (forward ref).
- [x] Delete the unused `SolveResult` TypedDict (protocols.py:65) — already
      dead scaffolding, superseded by `EquilibriumResult`.
- [x] Test coverage added to `tests/standalone/test_speciation_protocols.py`
      (existing "protocols.py tests" file — a dedicated
      `tests/speciation/test_equilibrium_result.py` would have duplicated
      that file's role; `TestSolveResult` was replaced with
      `TestEquilibriumResult` in place, since deleting `SolveResult` broke
      that section's import): construction (minimal + full), frozen-immutability,
      shared-mutable-default-dict guard, `to_dict()` round-trip,
      `apply_to_phases()` against a real `LiquidPhase` (V_L scaling,
      untouched-species preservation, no-op guards for missing liquid /
      no `_refresh_derived`).
- [x] Ran `tests/standalone/test_speciation_protocols.py`: 27 passed.

### CP2 — `SpeciationEngine.solve()` → `EquilibriumResult` — SHIPPED 2026-07-01 (commit `3864a2e`)

- [x] In `src/speciation/engine.py::solve()`, after the existing `out` dict
      is fully assembled (through the `alphas` block), construct
      `EquilibriumResult`: `species_mol_L` = the existing hardcoded canonical
      tuple (now `_CANONICAL_WRITEBACK_SPECIES` module constant) filtered to
      keys present in `out`; everything else not matching a named field goes
      into `extra`.
- [x] Delete the writeback block — superseded by `apply_to_phases()`.
- [x] Update `get_CO2aq_from_totals()`: checks both `species_mol_L` and
      `extra` for `"CO2aq"`/`"CO2"` — "CO2aq" isn't in the canonical tuple
      (only "CO2" is), so it can land in `extra` depending on solve path;
      checking only `species_mol_L` would have silently broken this method.
- [x] Delete the redundant `SpeciationResult` dataclass — was already dead
      (`as_dict` branch never triggered). Also removed its
      `speciation/__init__.py` re-export (replaced with `EquilibriumResult`)
      — deleting the class alone broke import collection repo-wide via
      `from .engine import SpeciationEngine, SpeciationResult`.
- [x] Update module docstring to describe `EquilibriumResult`.
- [x] `tests/standalone/test_speciation.py`: converted 9× `out["pH"]` /
      `out.get("pH", np.nan)` patterns and 6 named-variable variants
      (`out_ideal`, `out_davies`, `out1/out2`, `out_ws/out_no`,
      `out_cold/out_hot`, `out_a/out_b`) to `.pH` attribute access; added a
      `_species(result)` test helper (`{**species_mol_L, **extra}`) for the
      3 `TestCanonicalEmission` membership-check tests, since custom acids
      like `"AceticAcid"` aren't in the canonical writeback tuple and land
      in `extra`, not `species_mol_L`.
- [x] Ran `tests/standalone/test_speciation.py` + `test_speciation_protocols.py`:
      64 passed.
- [x] **Known transient state, by design**: production call sites
      (`control_volume.py`, `solvers.py`, `simulation.py`) still discard
      `solve()`'s return value and don't call `apply_to_phases()` yet (CP6).
      Until CP6 lands, `SpeciationEngine`-backed CV/simulation paths (the
      default engine — `ReactionSystem.engine` picks `SpeciationEngine`
      unless `solver="newton_raphson"`) do not receive derived-species
      writeback. Full-suite run intentionally deferred to CP6/CP7 — do not
      treat a full-suite run between CP2 and CP6 as a regression signal.

### CP3 — `NRSpeciationEngine.solve()` → `EquilibriumResult` — SHIPPED 2026-07-01 (commit `1a40d48`)

- [x] In `src/speciation/nr_engine.py::solve()`, keep popping the white-box
      internal keys (`_jacobian_matrix`, `_log_activities`,
      `_ionic_strength_final`) into `self._cached_*` **before** constructing
      `EquilibriumResult` (order matters — don't let them leak into
      `species_mol_L`/`extra`).
- [x] `species_mol_L` = `self._tableau.masters` + `self._tableau.secondaries`
      + `"H2O": _C_WATER_MOL_L` (mirrors `_writeback` exactly — no canonical
      restriction, unlike `SpeciationEngine`).
- [x] `saturation_indices` = `{label: info["SI"] for label, info in
      out.get("minerals", {}).items()}`; `extra["minerals_xi_mol_L"]` =
      `{label: info["xi_mol_L"] ...}`.
- [x] Delete the `_writeback` method and its call site — superseded by
      `apply_to_phases()`.
- [x] Verify `jacobian_dg_dz()`, `jacobian_dg_dy()`, `residual()`,
      `jacobian_dz_dy()` read only `self._cached_*`, never the `solve()`
      return value — confirm no hidden dependency before/after this change.
- [x] Enumerate `tests/` files using `NRSpeciationEngine`/`from_reactions(` +
      `.solve(`; convert dict access (including `minerals` assertions →
      `saturation_indices`/`extra`) to attribute access.
- [x] Run affected NR-engine test files.

### CP4 — `PHREEQCEngine.solve()` → `EquilibriumResult` — SHIPPED 2026-07-01 (commit `bcebf82`)

- [x] In `src/speciation/phreeqc_engine.py::solve()`, after
      `out = self._build_output(sol)`, construct `EquilibriumResult`:
      `species_mol_L` = all `out` keys except `{"pH", "logH",
      "IonicStrength"}`; `ionic_strength` = `out["IonicStrength"]`.
- [x] Delete the `_writeback` method and its call site — superseded by
      `apply_to_phases()`.
- [x] Update module docstring example and `solve()` docstring's call-pattern
      examples to attribute form.
- [x] Locate and update PHREEQC test file(s) (likely skip-marked when
      `phreeqpython` isn't installed — update source regardless).

### CP5 — Wrapper updates — SHIPPED 2026-07-01 (commit `cd75dec`)

- [x] `numerical_gradient.py::jacobian_dz_dy()`: `out_plus.get(sp_id)` /
      `out_minus.get(sp_id)` → `out_plus.species_mol_L.get(sp_id)` /
      `out_minus.species_mol_L.get(sp_id)`.
- [x] `numerical_gradient.py::solve()` passthrough type hint → `EquilibriumResult`.
- [x] `speciation/api.py::SpeciationEngineAdapter.equilibrate()` — **only if
      CP0 recommends "keep as public API"**, otherwise this item is replaced
      by CP0's deletion: `out.get("pH")` → `out.pH`,
      `out.get("CO2")` → `out.species_mol_L.get("CO2")`,
      `out.get("IonicStrength")` → `out.ionic_strength`,
      `raw=dict(out)` → `raw=out.to_dict()`, `out.get("n_iter")` → `out.n_iter`.
- [x] Confirm `equilibria/engine.py` / `equilibria/coupled.py` are left
      untouched here regardless of CP0's outcome (CP0 owns their fate on its
      own branch).
- [x] Run tests covering `numerical_gradient.py` and (if kept) `speciation/api.py`.

### CP6 — Production call-site cutover (atomic — lands as one unit) — SHIPPED 2026-07-01 (commit `99cd071`)

- [x] `control_volume.py:576` (`advance()`, step 1): capture `result`, call
      `result.apply_to_phases(self.phases)` before step 1b (property
      calculators) and step 3 (`compute_reaction_rates`).
- [x] `control_volume.py:651` (`compute_rhs()`): commit before the trailing
      `return self.compute_reaction_rates(t_h)`.
- [x] `solvers.py:240` (snapshot step solver, "2b. Speciation"): commit to
      `{"liquid": snap_liq, "gas": snap_gas}` before step 2c (transfer) and
      step 2d (`compute_reaction_rates` on `snap_liq_cv`).
- [x] `solvers.py:315` ("6. FINAL SPECIATION on live state"): commit before
      `_monitor_pH_post_step(cv)`.
- [x] `solvers.py:613` (ScipyODESolver RHS closure): commit to
      `{"liquid": _temp_liq, "gas": _temp_gas}` — confirmed unconditional
      commit on every RHS eval (including rejected trial steps) is safe:
      these are scratch snapshots rebuilt each call and discarded except
      when `sv.write_to_cv(...)` persists the accepted final state.
- [x] `solvers.py:791` ("Final speciation on live state" post-integration):
      commit before `_monitor_pH_post_step(cv)`. Used `eq_result` as the
      variable name (not `result`) to avoid shadowing the pre-existing
      scipy `solve_ivp` result object in scope at that point.
- [x] `simulation.py:1236` (`_initial_solve`): commit, keeping the existing
      `try/except (TypeError, ValueError, AttributeError, KeyError)`
      wrapping around both the `solve()` and `apply_to_phases()` calls.
- [x] Re-ran the `engine.solve(` inventory grep against the final diff —
      confirmed exactly 7 production call sites, no 8th missed.
- [x] Traced each site forward to its next `phase.n_mol`/`phase.pH` read —
      all correctly ordered.
- [x] **Found and fixed one casualty**: `test_simulation.py`'s
      `TestC9EndToEndPHController.FakeSpeciation` test double predated the
      new contract — it wrote directly to `phase.n_mol` inside `solve()`
      and returned `None` (implicit), which broke under
      `result.apply_to_phases()` in `_initial_solve` (`AttributeError`,
      caught by the existing broad except there, surfacing as a NaN-pH
      test failure rather than a clear error). Fixed to return a real
      `EquilibriumResult`. Full suite: 1798 passed, 27 skipped, 0
      failures — first clean full run since CP2.

### CP7 — Test suite migration + closeout — SHIPPED 2026-07-01 (commit `337e23e`)

- [x] Re-ran repo-wide grep for dict-style access; the only hits outside
      already-migrated files were `tests/legacy/*` (excluded from
      `testpaths`, out of scope) and two pre-existing `_HighIEngine` test
      doubles in `test_accuracy_monitor.py` that return plain dicts but are
      never routed through `apply_to_phases()` in their own tests (testing
      `AccuracyMonitor`'s warning mechanism, not the solve() contract) —
      left as-is, not a live bug, noted here for visibility rather than
      speculatively "fixed".
- [x] Converted all real dict-access sites (done incrementally in
      CP2-CP6, verified complete here).
- [x] Added `TestApplyToPhasesWriteback` in `test_speciation.py` (NR and
      PHREEQC already had equivalent coverage via `TestPhaseWriteback` /
      `TestPhasesWriteback`, updated in CP3/CP4) — asserts
      `apply_to_phases()` writeback exactly matches `species_mol_L * V_L`,
      and that `solve()` itself does not auto-write (sentinel check).
- [x] Added `test_repeated_advance_updates_pH_each_step` in
      `test_cv_advance.py` — three `cv.advance()` calls with a feed
      injected between steps 1-2; asserts pH shifts by step 3 (once
      speciation re-solves on the fed totals per the pre-step-state
      operator-splitting invariant). Calibrated threshold to the actual
      observed shift (~0.02 pH units) after an initial too-strict 0.05
      threshold false-failed on first run.
- [x] Full suite run: 1793 → 1800 passing, 27 skipped, 0 failures.
- [x] Updated `docs/design/CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md`:
      header banner, §8 resolution note, §17 phase table (split Phase 1
      into 1 + 1b to reflect that `EquilibriumResult` shipped separately
      from the rest of Phase 1's original scope).
- [x] Moved this file to `docs/phases-shipped/EQUILIBRIUM_RESULT.md`.
- [x] Tag `equilibrium-result-shipped`, merge `--no-ff` into `main`,
      delete the `equilibrium-result` branch.

---

## Files changed

| File | Nature of change | Outcome |
|---|---|---|
| `src/equilibria/engine.py`, `coupled.py`, `factory.py`, `interfaces.py`, `src/core/gl_equilibrium.py` | CP0: audit | Deleted (`legacy-cleanup-shipped`) — confirmed dead, unreachable since `CUFERMENTER_SUNSET` |
| `src/speciation/protocols.py` | CP1: add `EquilibriumResult`, `apply_to_phases()`, `to_dict()`; delete `SolveResult` | Shipped |
| `src/speciation/engine.py`, `__init__.py` | CP2: `solve()` returns `EquilibriumResult`; delete writeback block + `SpeciationResult` (+ its re-export) | Shipped |
| `src/speciation/nr_engine.py` | CP3: `solve()` returns `EquilibriumResult`; delete `_writeback` | Shipped |
| `src/speciation/phreeqc_engine.py` | CP4: `solve()` returns `EquilibriumResult`; delete `_writeback` | Shipped |
| `src/speciation/numerical_gradient.py` | CP5: attribute access in `jacobian_dz_dy()` | Shipped |
| `src/speciation/api.py` | CP5: attribute access, `to_dict()` | Shipped — CP0 kept it (unused but not broken) |
| `src/core/control_volume.py` | CP6: 2 call sites, explicit `apply_to_phases()` | Shipped |
| `src/core/solvers.py` | CP6: 4 call sites, explicit `apply_to_phases()` | Shipped |
| `src/core/simulation.py` | CP6: 1 call site, explicit `apply_to_phases()` | Shipped |
| `src/core/property_calculator.py` | CP6: docstring wording update | Shipped |
| `tests/` | CP1-CP7: new tests + dict-access migration across `test_speciation.py`, `test_speciation_protocols.py`, `test_nr_speciation_engine.py`, `test_nr_whitebox.py`, `test_phreeqc_engine.py`, `test_numerical_gradient_engine.py`, `test_simulation.py`, `test_cv_advance.py` | Shipped |
| `docs/design/CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md` | CP7: mark shipped | Shipped |

---

## How it went

Ran CP0 → CP7 in order as planned, one commit per checkpoint (plus a couple
of doc-only commits marking progress). No checkpoint required re-scoping.
Two things surfaced during execution that weren't fully anticipated in the
plan:

1. A test double (`test_simulation.py`'s `FakeSpeciation`) predated the new
   `solve()` contract and broke under CP6's cutover — fixed as part of that
   checkpoint's commit, not spun out separately, since it was a direct
   consequence of the cutover itself.
2. The `species_mol_L` vs `extra` split for `SpeciationEngine` (canonical
   tuple only) turned out to matter for more call sites than expected —
   `get_CO2aq_from_totals()`, `speciation/api.py`'s CO2 lookup, and a
   cross-engine test comparison all needed to check both buckets, since
   `"CO2aq"` isn't in the canonical tuple but `"CO2"` is.

Final: 1793 → 1800 tests, 27 skipped (unchanged — `phreeqpython` not
installed in this environment), 0 failures.
