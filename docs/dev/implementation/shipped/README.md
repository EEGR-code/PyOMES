# Phases — Shipped

Historical record of completed refactor phases. Each entry's planning
docs and implementation logs are kept here indefinitely as the
authoritative record of what was decided and what was done. Status
banners at the top of each doc indicate when the work landed.

These docs are **frozen**: once a phase ships, its plan and checklist
do not get updated except for status banners or factual corrections.
Forward-looking ideas that came out of the work go to
[../upcoming/](../upcoming/) as separate planning notes.

## Bug fixes

### ScipyODESolver — two latent bugs fixed (2026-06-18)

Two bugs in [src/core/solvers.py](../../src/core/solvers.py) were preventing
`ScipyODESolver` from running at all when used with a CV that has gas-liquid
transfer (e.g. a fermenter with `KineticGasLiquidLink`):

- **`spec_result` NameError** — the ODE RHS closure `f()` referenced
  `spec_result` before it was assigned. The speciation engine writes its
  results back into `_temp_liq.n_mol` directly via `_refresh_derived`; no
  result object is returned. Fix: `spec_result = None` initialised alongside
  `pH = None` at the top of the closure.

- **Stale `property_results` kwarg** — `link.compute_flux()` was called with
  `property_results={...}` which was removed from
  `KineticGasLiquidLink.compute_flux()` during an earlier phase. The kwarg was
  removed; the call is now
  `link.compute_flux(gas, liq, dt_h, instantaneous=True)`.

Both were introduced when `ScipyODESolver` was wired up but the full
chemistry-enabled path was never exercised. 1458/1458 tests pass after the fix.

## Demo and documentation improvements

### Aerobic fermentation dynamic simulation (2026-06-18)

Added Section 6 to
[demos/aerobic_fermentation_stoichiometry.ipynb](../../demos/aerobic_fermentation_stoichiometry.ipynb):
a fully dynamic batch simulation complementing the existing CHNOSP
stoichiometry analysis.

**Setup**
- Initial medium: 30 g/L glucose, 30 g/L KH₂PO₄, 5 g/L NH₄Cl,
  0.5 g/L K₂SO₄, 0.4 g/L MgSO₄·7H₂O
- Standard `FermenterBuilder` + `Simulation` architecture; 5 L vessel,
  20 % headspace, T = 310.15 K, air sparge at 1 vvm, kLa(O₂) = 200 h⁻¹
- Full acid-base equilibria via `BIOPROCESS_BASIC`: phosphate ladder
  (pKa 2.15/7.20/12.35), bisulfate (pKa 1.99), ammonium (pKa 9.25),
  carbonate (pKa 6.35)
- Simple Monod kinetics (μ_max = 0.5 h⁻¹, K_S = 0.05 g/L), Y_XS = 0.4 g/g;
  two biomass compositions — Roels extended and Upcraft

**Key implementation notes**
- Kinetic reaction built via `KineticReaction` + `StoichiometryEntry` using
  shared objects from `common_species` (`NH4_plus`, `CO2`, `H2O`) — required
  to avoid `SpeciesConflictError` when combining user-defined kinetics with
  `BIOPROCESS_BASIC` equilibrium reactions that already reference the same
  species with canonical charge/MW values.
- NH₄⁺ used as the nitrogen source species (not NH₃): at fermentation
  pH 4.5–7, NH₃ is < 0.1 % of total ammoniacal nitrogen, so consuming from
  the NH₃ pool would clamp to zero each step and break the N mass balance.
- Adaptive DOP853 via `ScipyODESolver(rtol=1e-4, atol=1e-9, max_step=0.1)`
  instead of fixed-step Euler, to handle rapid O₂ drawdown at peak growth
  without needing an O₂ co-limitation Monod term.

**Outputs**
- Per-composition concentration report (initial / final; all species in
  mmol/L; total phosphate CT_P; pH at start and end)
- 6-panel time-series plot saved to
  `demos/aerobic_fermentation/simulation_timeseries.png`

### FBA demo improvements (shipped 2026-06-01)

Demo-quality work; no `src/` or `models/` changes. Closed the FBA
stoichiometry verification item carried forward from the 2026-05-28
exploration session.

Key outcomes:
- **`fba_toy.py`**: attribution corrected — Orth, Thiele & Palsson
  (2010) Box 1 is the S-matrix formalism, not a toy reaction list;
  "Future work" disclaimer removed; misleading print statement fixed.
- **`ecoli_core.json`**: network extended from 14 to 18 reactions.
  PFL alone does not fix anaerobic growth in the lumped model (NADH
  imbalance remains); the minimal fix is PFL + ADH together. Added
  PFL, ADH, EX_for_e, EX_etoh_e reactions and for_e / etoh_e
  metabolites. "UNVERIFIED / FUTURE WORK" framing replaced with
  clear scope documentation.
- **`fba_ecoli_core.py`**: attribution updated; Formate + Ethanol
  added to trajectory output; anaerobic-behaviour note corrected.
- **`demos/model_api/README.md`** and **`demos/README.md`**: stale
  "~14-reaction / future replacement" and "attribution warnings"
  descriptions updated.

- **[FBA_DEMOS.md](FBA_DEMOS.md)** — full record of changes and
  the stoichiometry findings that drove them.

## Phases shipped

### Simulation class — `simulation-class` (shipped 2026-05-27)

Third and final of the three sequenced phases. New `Simulation`
class is the single orchestration pathway for all CV-based models:
subsumes `MultiCVSystem` (deleted), absorbs `run_batch`'s time loop,
and promotes controllers and profiles to first-class members with
structured per-step records in the new `BatchResult`. Pattern 1
controller protocol (single `compute(state, dt_h) -> ControlAction`;
legacy `Commands` / `Actuator` deleted); RunContext-based lifecycle
gating across CVs / Phases / links / `_LockableList` wrappers;
Pattern B unchecked-setter dispatch for orchestrator-mediated
mutation. Final standalone suite 977/0. Legacy
`CUFermentationSpeciation` island (run_batch, legacy `loops.py` /
`system.py`, `vlmodels/fermenter/types.py`) retained for its
out-of-scope-but-still-alive consumers.

- **[SIMULATION_CLASS.md](SIMULATION_CLASS.md)** — design note
  with 12 resolved decisions (incl. RunContext mechanism, Pattern 1
  collapse, Pattern B unchecked-setter dispatch). Status banner
  at top.
- **[SIMULATION_CLASS_CHECKLIST.md](SIMULATION_CLASS_CHECKLIST.md)**
  — implementation log with per-checkpoint commit / test deltas
  (C0 design lock-in through C14 deletion sweep + C15 doc-shipping).
  Tag: `simulation-class-shipped`.

### CV refactor — Phases 1 through 5 (shipped 2026-04-29)

Made `ControlVolume.advance()` run property solvers first, threaded
`property_results` through reactions and phase interfaces, eliminated
the `_last_properties` cache, made `KineticGasLiquidLink` a dual
`CVLink` / `PhaseInterface`, and reduced `GasLiquidVolume` to a thin
wrapper around a single `ControlVolume` with both phases inside. 740
standalone tests passing throughout.

- **[CV_UPDATE.md](CV_UPDATE.md)** — original plan covering all five
  phases. Status banner at the top confirms the shipping date and
  test count. Body preserved as the historical plan.
- **[PHASE5_CHECKLIST.md](PHASE5_CHECKLIST.md)** — implementation log
  for Phase 5 (the structural collapse to single-CV). All 52
  checkpoints ticked.
- **[ORDERING.md](ORDERING.md)** — design rationale: feed-before-
  reactions argument and the correct advance sequence.
- **[ORDERING_CRITIQUE.md](ORDERING_CRITIQUE.md)** — three independent
  expert reviews of ORDERING.md, surfacing nuances around
  operator-splitting and EulerSnapshotSolver alignment.

### Solver promotion — Phase 6 (shipped 2026-05-01)

Promoted `EulerSnapshotSolver` and `ScipyODESolver` from being
`GasLiquidVolume`-attached strategies to being usable directly via
`cv.advance(solver=...)` on any `ControlVolume`.  Boundaries became a
CV-level concept; `GasLiquidAdvanceResult` was dropped (the unified
`AdvanceResult` gained `transfer_record` and `boundary_records`);
`MultiCVSystem.advance_all` exposes per-CV solver choice through a
`solvers={cv_key: solver}` parameter; the remaining `context` →
`chem_env` naming holdovers (`context_fn`, `add_context`,
`apply_context`) were renamed.  747 standalone tests passing
(740 baseline + 7 new tests).

- **[SOLVER_PROMOTION.md](SOLVER_PROMOTION.md)** — original design
  note, kept as a historical record.
- **[PHASE6_CHECKLIST.md](PHASE6_CHECKLIST.md)** — implementation log.
  All 70 checkpoints ticked.

### BSM2 reference test (shipped 2026-05-12)

Pre-chemistry-unification regression safety net for the BSM2
reference model ([`models/vlmodels/adm1/bsm2.py`](../../models/vlmodels/adm1/bsm2.py)).
Added [`tests/standalone/test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py)
with three hardcoded end-of-trajectory sentinel assertions (liquid
concentrations, gas mole counts, pH at `rel_tol=1e-9`) plus three
invariants (pH physically valid, concentrations non-negative,
speciation populated). Also added a "Forward note" docstring on
`SpeciationPropertySolver` cross-referencing
[CHEMISTRY_UNIFICATION.md](../upcoming/CHEMISTRY_UNIFICATION.md)
and the class-diagrams deferred-refactor block, so the in-code class
itself signals its scheduled absorption into the unified `Reaction`
framework. 714 standalone tests passing (708 baseline + 6 new).

- **[BSM2_REFERENCE_TEST.md](BSM2_REFERENCE_TEST.md)** — checklist
  with goal, scope, resolved decisions, and final test count.

### GLV removal — Phase 7 (shipped 2026-05-12)

Deleted `GasLiquidVolume` entirely after Phase 6 reduced it to a
near-transparent shim.  Fermenters are now plain `ControlVolume`s with
`{"gas", "liquid"}` phases plus a `KineticGasLiquidLink` as an
internal interface; callers construct one via
`FermenterFactory.create_volume(...)` (which inlines the link wiring)
or `FermenterBuilder().build()`.  Variables named `glv` were swept to
`cv` across `models/`, `systems/`, and `tests/`; the `glv.transfer_link`
accessor was replaced by an `internal_interfaces` lookup against
`KineticGasLiquidLink` (option 6b); the `pH` / `ionic_strength` /
`_last_result` accessors were dropped (callers read
`result.properties["speciation"].pH` off the `AdvanceResult` instead).
The `ADM1` and `BSM2` builders were renamed
(`build_adm1_glv` → `build_adm1_cv`, `build_bsm2_glv` → `build_bsm2_cv`),
and `ThermodynamicConfig.apply_to_glv` → `apply_to_cv`.  712 standalone
tests passing (744 pre-deletion, minus the deleted
`test_gas_liquid_volume.py` whose coverage was already duplicated in
`test_cv_advance.py`'s `TestCVWithKineticGasLiquidLink`).

- **[GLV_REMOVAL.md](GLV_REMOVAL.md)** — design note, kept as a
  historical record.
- **[PHASE7_CHECKLIST.md](PHASE7_CHECKLIST.md)** — implementation log.
  All 5 checkpoints ticked.

## Cross-references

- Solver characterisation post-refactor:
  [../solvers.md](../solvers.md).
- Architecture overview reflecting the post-Phase-7 shape:
  [../architecture.md](../architecture.md).
- Class diagrams reflecting the post-Phase-7 shape:
  [../class_diagrams.md](../class_diagrams.md).
- Open follow-on work after the CV refactor:
  [../upcoming/CHEMISTRY_UNIFICATION.md](../upcoming/CHEMISTRY_UNIFICATION.md),
  [../upcoming/RUN_HISTORY.md](../upcoming/RUN_HISTORY.md),
  [../upcoming/CONTAINER_LAYERING.md](../upcoming/CONTAINER_LAYERING.md).
