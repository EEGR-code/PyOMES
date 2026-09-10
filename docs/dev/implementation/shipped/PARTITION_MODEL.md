# PARTITION_MODEL — Phase Checklist

> **Status:** Shipped 2026-06-16. Tag: `partition-model-shipped`. Suite: 1428/0.

---

## Background and motivation

The original trigger for this phase was labelled "VFA Full-Closure": correct
the gas-liquid alpha factor for VFA species at non-trivial pH.  Investigation
revealed two things that reframe the problem:

1. **VFA stripping is negligible regardless of alpha.** VFA Henry constants
   (~5 000 mol/L/atm for acetic acid) are ~160 000× larger than CO₂.  Even
   with alpha = 1.0, essentially zero VFA transfers to the gas phase at
   operating pH.  The alpha error is real but numerically inert for VFAs.

2. **H₂S is the genuine accuracy gap.** H₂S has kH ≈ 0.10 mol/L/atm and
   pKa = 7.0, which lands squarely in the AD operating pH range.  The alpha
   correction changes the predicted equilibrium gas-phase fraction from ~9%
   (alpha = 1.0) to ~5% (alpha = 0.5 at pH 7) — a 2× error that propagates
   directly into dissolved sulfide inhibition kinetics.

3. **Henry constants live in the wrong layer.** The `_HENRY_PARAMS` table in
   `models/vlmodels/fermenter/config/factory.py` is a private fermenter-builder
   artefact.  Users cannot inspect, extend, or override it through the same
   declaration surface as the rest of the chemistry.  The `vfa_volatility` flag
   in `build_adm1_reactions()` / `build_adm1_cv()` compounds this: it couples
   species identity ("these four VFAs") to thermodynamic behaviour ("are
   volatile") inside a model-specific builder, rather than letting the user
   declare volatility as part of the chemistry definition.

This phase addresses both the accuracy gap (H₂S) and the architectural gap
(Henry constants as user-declared chemistry).  VFA gas transfer is out of
scope: the standard BSM2 and ADM1 models do not include it, and the numerical
impact of the alpha error is negligible in any case.

---

## Goal

1. Introduce a `PartitionModel` protocol in the chemistry layer so that
   phase-partition relationships (Henry's law, Langmuir, Raoult, etc.) are
   declared explicitly by the user alongside reactions and pKas, rather than
   looked up from a private builder table.

2. Add `partition_models: Dict[str, PartitionModel]` to `ChemistryDatabase`
   and populate the stock databases with `HenryPartition` objects for the
   species they already cover.

3. Fix the H₂S alpha correction by declaring the H₂S ⇌ HS⁻ equilibrium with
   `species_refs` in `AD_BASIC`, giving the speciation engine a real two-element
   ladder and enabling correct alpha computation in `_alpha_for`.

4. Update `KineticGasLiquidLink` to read `PartitionModel` objects from the
   chemistry database rather than carrying its own Henry constant dicts.

5. Remove the `vfa_volatility` flag from the ADM1 builders, making those
   builders faithful to the published BSM2/ADM1 models.

6. Delete `_HENRY_PARAMS` and `henry_mol_L_atm()` from `factory.py`.

---

## Scope

**In scope:**
- `PartitionModel` protocol + `HenryPartition` concrete implementation
- `ChemistryDatabase.partition_models` field and `.extend()` support
- H₂S / HS⁻ species declaration and equilibrium in `AD_BASIC`
- Henry constants for CO₂, CH₄, H₂, O₂, N₂, NH₃, H₂S moved to stock databases
- `KineticGasLiquidLink` reads partition models from chemistry database
- `vfa_volatility` flag removed from `build_adm1_reactions()` / `build_adm1_cv()`
- BSM2 reference sentinel re-baseline if H₂S alpha correction affects any
  existing sentinel (only if BSM2 runs with `sulfate_reduction=True`)

**Out of scope:**
- VFA acid-base equilibria (VFA stripping is negligible; no kinetics update needed)
- Non-linear partition models (`LangmuirPartition`, `FreundlichPartition`) — the
  `PartitionModel` protocol is designed to accommodate them, but only
  `HenryPartition` ships in this phase
- `KineticSolidLiquidLink` or any new link type
- `_HA` / `_A-` legacy fallback removal (separate clean-up)
- HPC_CHECKPOINTING (separate trigger-gated phase)

---

## Resolved design decisions

**`PartitionModel` is a protocol, not a base class.**
Consistent with the rest of the codebase (SystemSolver, ControllerBase, etc.).
`HenryPartition` is a frozen dataclass that satisfies the protocol.

**`beta()` returns `Optional[float]`.**
`HenryPartition.beta()` returns a finite float; the exponential analytical step
solution in `_kinetic_flux` is exact for linear (Henry) partitions.  Future
non-linear models (Langmuir) return `None` from `beta()` to signal that the ODE
instantaneous-rate path must be used.  The backend already has both paths;
`None` from `beta()` routes to the existing `instantaneous=True` branch.

**Capacity parameters are phase-agnostic.**
`PartitionModel` methods take `capacity_a` and `capacity_b` (generic floats)
rather than `V_liq` / `V_gas`.  `HenryPartition` interprets them as volumes
(L); a future `LangmuirPartition` would interpret `capacity_b` as solid mass
(kg).  This keeps the protocol usable for gas-liquid, solid-liquid, and
gas-solid partition without schema changes.

**Volatility is implicit in the database.**
Any species present in `chemistry_db.partition_models` is automatically
gas-transferable.  The `FermenterBuilder` reads `partition_models` at build
time rather than requiring explicit `transfer_species()` calls for species
whose Henry constants are already in the database.  kLa and transfer mode
(kinetic vs. equilibrium) remain builder-level operational parameters — they
are not thermodynamic properties and do not belong in the database.

**`_HENRY_PARAMS` and `henry_mol_L_atm()` are deleted, not deprecated.**
The package is pre-release; no backward-compatibility shim is needed.
All callers are internal and are updated in the same phase.

**H₂S equilibrium goes in `AD_BASIC`, not `BIOPROCESS_BASIC`.**
H₂S / HS⁻ chemistry is specific to anaerobic digestion contexts.
`BIOPROCESS_BASIC` is used by aerobic models where sulfide speciation is
irrelevant.

**`vfa_volatility` is removed, not deprecated.**
The flag is a non-standard extension that produces silently wrong results.
Removing it makes `build_adm1_reactions()` and `build_adm1_cv()` accurate
to the published ADM1 specification.  Users who genuinely want VFA gas
transfer declare VFA `HenryPartition` objects and equilibria explicitly in
their own database extension.

---

## Checkpoints

Ordered so each checkpoint leaves the test suite runnable.

### C1 — `HS_minus` Species + H₂S equilibrium in `AD_BASIC`

- [ ] `src/chemistry/common_species.py` — add `HS_minus = Species(id="HS-", ...)`
      (charge=-1, atoms={"H":1,"S":1}).
- [ ] `src/chemistry/databases/anaerobic_digestion.py` — add monoprotic
      `EquilibriumReaction` for H₂S ⇌ HS⁻ + H⁺ with `log_K = -7.0`,
      `species_refs = (H2S_sp, HS_minus)`, `correction="van_t_hoff"`,
      `dH_J_per_mol = 20000.0`.  Wire using `Species` objects imported from
      `common_species` (or declared locally) so the ladder is two-element.
- [ ] `tests/standalone/test_chemistry_database.py` — assert `AD_BASIC`
      single-phase equilibria include an H₂S entry; assert `species_refs` is
      `(H2S_sp, HS_minus)`; assert `HS-` is present in `AD_BASIC.species`.
- [ ] Verify suite green (additive change only).

### C2 — `PartitionModel` protocol + `HenryPartition`

- [ ] New module `src/chemistry/partition.py`:
  ```python
  class PartitionModel(Protocol):
      def equilibrium_a_moles(
          self, n_total: float, capacity_a: float, capacity_b: float,
          T_K: float, alpha: float = 1.0,
      ) -> float: ...

      def beta(
          self, capacity_a: float, capacity_b: float,
          T_K: float, alpha: float = 1.0,
      ) -> Optional[float]:
          """Return dimensionless partition coefficient for the analytical
          exponential step solution.  Return None for non-linear models that
          require the ODE instantaneous-rate path instead."""
          ...

  @dataclass(frozen=True)
  class HenryPartition:
      H_ref: float        # mol/(m³·Pa), Sander convention
      dlnH: float         # K, van't Hoff d(ln kH)/d(1/T)
      T_ref: float = 298.15

      def beta(self, capacity_a, capacity_b, T_K, alpha=1.0) -> float:
          kH = self._kH_mol_L_atm(T_K)
          return (kH / max(1e-12, alpha)) * R_L_ATM_MOL_K * T_K * capacity_a / capacity_b

      def equilibrium_a_moles(self, n_total, capacity_a, capacity_b, T_K, alpha=1.0):
          b = self.beta(capacity_a, capacity_b, T_K, alpha)
          return b * n_total / (1.0 + b)

      def _kH_mol_L_atm(self, T_K):
          H_ref_mol_L_atm = (self.H_ref / 1000.0) * 101325.0
          return H_ref_mol_L_atm * math.exp(self.dlnH * (1.0 / T_K - 1.0 / self.T_ref))
  ```
- [ ] Export `PartitionModel` and `HenryPartition` from `src/chemistry/__init__.py`.
- [ ] `tests/standalone/test_partition_model.py` — construction, `beta()`,
      `equilibrium_a_moles()`, temperature dependence, alpha correction,
      `dataclasses.replace()` override pattern.
- [ ] Verify suite green (new module, no existing code touched).

### C3 — `ChemistryDatabase` gains `partition_models`

- [ ] `src/chemistry/database.py` — add `partition_models: Dict[str, PartitionModel]`
      field (default empty dict, frozen).
- [ ] `ChemistryDatabase.extend()` — merge `partition_models` dicts (child
      overrides parent for the same species key).
- [ ] `tests/standalone/test_chemistry_database.py` — extend existing tests:
      assert `partition_models` merges correctly; assert child entry overrides
      parent entry for same species key.
- [ ] Verify suite green (additive field, no existing behaviour changed).

### C4 — Stock databases populated with `HenryPartition` objects

Moves values from `_HENRY_PARAMS` in `factory.py` into the databases where
they belong.  `_HENRY_PARAMS` is not deleted yet (C6 does that).

- [ ] `src/chemistry/databases/bioprocess_basic.py` — add `partition_models`:
      `O₂`, `N₂` `HenryPartition` objects (values from `_HENRY_PARAMS`).
- [ ] `src/chemistry/databases/anaerobic_digestion.py` — add `partition_models`:
      `CO₂`, `CH₄`, `H₂`, `NH₃`, `H₂S` `HenryPartition` objects (values from
      `_HENRY_PARAMS`).  `H₂S` entry here matches the equilibrium added in C1.
- [ ] `tests/standalone/test_chemistry_database.py` — assert each stock database
      exposes the expected `partition_models` keys; spot-check one kH value at
      298.15 K against the Sander compilation reference.
- [ ] Verify suite green (database fields populated; link still reads from its
      own `henry` dict at this checkpoint).

### C5 — `KineticGasLiquidLink` reads `PartitionModel`

Breaking change: replaces the link's `henry` / `henry_params` dicts with
`partition_models`.  Tests that construct links with explicit `henry=` dicts
must be updated.

- [ ] `src/core/gas_liquid_link.py`:
  - Replace `henry: Dict[str, float]` and `henry_params` fields with
    `partition_models: Dict[str, PartitionModel]` (default empty dict).
  - `compute_flow()` iterates over `partition_models` keys (was `henry` keys).
  - `_kinetic_flux()`: call `model.beta(V_liq, V_gas, T_K, alpha)`.  If
    `beta()` returns `None`, fall back to `_instantaneous_rate()` path
    (set `instantaneous=True` for that species).
  - `_equilibrium_flux()`: call `model.equilibrium_a_moles(n_total, V_liq,
    V_gas, T_K, alpha)`.
  - `_instantaneous_rate()`: derive `C_star` from the model via
    `model.equilibrium_a_moles(n_total, V_liq, V_gas, T_K, alpha) / V_liq`.
  - `_effective_henry()`: remove (logic absorbed into `PartitionModel`).
  - Temperature recompute loop in `compute_flow()`: remove (temperature
    dependence moves into `HenryPartition._kH_mol_L_atm()`).
  - `set_henry()`, `set_kLa_with_co2_ratio()` mutators: update signatures;
    `set_henry()` now takes a `PartitionModel` not a float.
  - `validate()`: update to check `partition_models` keys instead of `henry`.
- [ ] `tests/standalone/test_gas_liquid_link.py` — update all fixtures to
      pass `partition_models={"CO2": HenryPartition(...), ...}` instead of
      `henry={...}`.  Assert H₂S alpha correction produces correct equilibrium
      fraction (target: ~5% in gas at pH 7, vs ~9% with alpha=1.0).
- [ ] Run BSM2 reference sentinels; expect zero drift (H₂S alpha correction
      only affects models using `sulfate_reduction=True`, not standard BSM2).

### C6 — Builder updated; `_HENRY_PARAMS` deleted

- [ ] `models/vlmodels/fermenter/config/factory.py`:
  - Delete `_HENRY_PARAMS` dict and `henry_mol_L_atm()` function.
  - Henry constant construction loop (lines 173–184) reads from
    `chemistry_db.partition_models` instead.  For each species in
    `chemistry_db.partition_models`, add its `HenryPartition` to the link's
    `partition_models`.  Explicit `henry_mol_L_atm` kwarg on
    `transfer_species()` is replaced by an explicit `PartitionModel` override.
- [ ] `models/vlmodels/fermenter/config/builder.py`:
  - `transfer_species()` `henry_mol_L_atm` parameter replaced by optional
    `partition_model: Optional[PartitionModel]`.  When `None`, the species
    is expected to be present in `chemistry_db.partition_models`.
  - Species present in `chemistry_db.partition_models` are automatically
    registered for transfer at build time without requiring an explicit
    `transfer_species()` call.
- [ ] Verify suite green.

### C7 — `vfa_volatility` flag removed; ADM1 builders cleaned

- [ ] `models/vlmodels/adm1/base.py`:
  - Remove `vfa_volatility` parameter from `build_adm1_reactions()` and
    `build_adm1_cv()`.
  - Remove the `if vfa_volatility:` block that appended cross-phase
    `EquilibriumReaction` entries for `AceticAcid`, `Propionate`,
    `Butyrate`, `Valerate`.
  - `sulfate_reduction=True` no longer calls `b.transfer_species("H2S")` —
    H₂S is now in `AD_BASIC.partition_models` and auto-registered.
  - Update any callers of `build_adm1_cv()` / `build_adm1_reactions()` that
    pass `vfa_volatility=True` (grep across `models/` and `tests/`).
- [ ] `models/vlmodels/adm1/bsm2_direct.py` — check for `vfa_volatility`;
      remove if present.
- [ ] `tests/standalone/` — grep for `vfa_volatility`; remove or update.
- [ ] Run full suite; re-baseline ADM1+sulfate_reduction sentinels if the
      H₂S alpha correction changes their steady-state values.

### C8 — Demo + documentation

- [ ] `demos/model_api/chemistry/partition_model.py` — demonstrates the
      explicit declaration pattern:
      1. Build a simple batch vessel (2 L liquid, 0.5 L headspace, buffered
         to pH 7) containing dissolved H₂S.
      2. Show predicted equilibrium gas-phase fraction with alpha = 1.0
         (wrong: ~9%) vs. the correct alpha-corrected result (~5%) by
         comparing a database that includes the H₂S equilibrium against one
         that does not.
      3. Show a user extending `AD_BASIC` with a custom species
         (`HenryPartition`) via `database.extend(partition_models={...})`.
- [x] `docs/OPEN_WORK.md` — VFA Full-Closure entry replaced with
      PARTITION_MODEL; marked as ready-to-implement (done 2026-06-12).
- [ ] `src/chemistry/__init__.py` — ensure `PartitionModel`, `HenryPartition`
      are in the public API.
- [ ] Verify 1390 + N tests green (N ≥ new tests in C1/C2/C3/C4/C5).

---

## What does NOT change

- `KineticGasLiquidLink` control flow (exponential approach, clamping,
  operator-splitting, `_alpha_for`, `_liquid_total`, `speciation_ladders`).
- BSM2 model (`bsm2.py`) — no VFA gas transfer in the standard model; H₂S
  alpha correction does not apply unless `sulfate_reduction=True` is used.
- The `_HA` / `_A-` legacy fallback in `_compute_species_eq` — left for a
  later clean-up phase.
- `ControlVolume`, `Simulation`, `SystemSolver` — no changes needed.
- VFA acid-base equilibria — not declared; VFA stripping remains at effective
  zero (correct for the standard models).
