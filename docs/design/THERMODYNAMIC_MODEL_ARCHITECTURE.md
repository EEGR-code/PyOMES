# Thermodynamic Model Architecture

> **Status:** Shipped 2026-07-01 — all 7 checkpoints implemented; 1793 tests
> passing. See commit `7a842d3`.
>
> This document captures the design vision that emerged from the mass-exchange
> architecture review. It covers the thermodynamic layer that sits beneath the
> `PartitionModel` / `PhaseInterface` machinery: activity models, equations of
> state, and the `ThermoFramework` that coordinates them.

---

## 1. Motivation

Three gaps motivate this design:

**1. `ActivityModel` is speciation-internal.**
`src/chemical_equilibrium/activity_models.py` defines an `ActivityModel` protocol and
implements Davies and SIT models. These correct equilibrium constants for
ion activity within the speciation engine. But partition models
(`HenryPartition`, the forthcoming `RaoultPartition`, `KspPartition`) currently
apply no activity correction to the dissolved-phase driving force. A CO₂ Henry
constant is not corrected for ionic strength; a Ksp driving force is not
corrected for ion activity. Thermodynamic consistency requires the same activity
model to be used throughout a CV — by both the speciation engine and the
partition models.

**2. `GasEOS` is not connected to the partition machinery.**
`src/equilibria/vle.py` defines a `GasEOS` protocol and `IdealGasEOS`.
`src/equilibria/peng_robinson.py` implements `PengRobinsonEOS`. These exist and
are correctly designed — `partial_pressures_atm` takes the full gas composition
and returns fugacities. But nothing in `src/chemistry/partition.py` or
`src/core/` uses them. `HenryPartition` assumes ideal gas. The PR implementation
is unreachable from the partition machinery.

**3. `ThermoFramework` holds strings, not objects.**
`src/thermo/framework.py` is correctly positioned as a CV-level thermodynamic
coordinator, but it stores `activity_model: str` and `use_activity: bool` rather
than an `ActivityModel` (or `LiquidPhaseModel`) instance. The actual model object
is constructed inside the speciation engine. A partition model that wants the
same activity model has no way to obtain it from `ThermoFramework`.

**The goal** is a coherent thermodynamic layer in which:
- Both phases (gas, liquid) have a first-class EOS protocol
- `ThermoFramework` holds instances of both
- The speciation engine and all partition models are constructed from the same
  `ThermoFramework`, ensuring thermodynamic consistency structurally

---

## 2. LiquidPhaseModel — liquid-side EOS

### 2.1 Symmetry with GasEOS

Gas-liquid equilibrium is governed by equality of fugacities:

```
f_i^gas  =  f_i^liq

φ_i(T, P, y) × y_i × P  =  γ_i(T, x) × x_i × f_i^ref
     ↑ GasEOS                    ↑ LiquidPhaseModel
```

`GasEOS` computes fugacity coefficients φ_i from the gas-phase composition.
`LiquidPhaseModel` computes activity coefficients γ_i from the liquid-phase
composition. Both are non-ideality corrections; both take the full phase
composition as input; both return per-species values. The protocols are
structurally symmetric.

### 2.2 Protocol definition

```python
class LiquidPhaseModel(Protocol):
    """Liquid-phase non-ideality model (liquid-side EOS).

    Symmetric with GasEOS: takes full phase composition, returns per-species
    activity coefficients.
    """
    name: str

    def gamma_all(
        self,
        x_mol: Dict[str, float],   # species concentrations (mol/L) or mole fractions
        T_K: float,
        *,
        charge: Dict[str, int],    # z_i per species — needed by ionic models
    ) -> Dict[str, float]:
        """Activity coefficients γ_i for all species at the given composition."""
        ...
```

The `charge` dict is passed as a keyword argument so non-ionic models (NRTL,
UNIQUAC) can ignore it without a distinct protocol. Species absent from the
return dict are treated as having γ_i = 1.0.

### 2.3 Implementations

**`IdealLiquidModel`** — γ_i = 1 for all species. The default; correct for
dilute non-ionic systems (CH₄, H₂, O₂ at bioprocess concentrations).

**`DaviesLiquidModel`** — Davies equation with temperature-dependent A(T).
Compresses `x_mol` to ionic strength `I = 0.5 Σ z_i² C_i` then applies:

```
log₁₀(γ_i) = −A(T) z_i² ( √I/(1+√I) − 0.3 I )
```

This is the existing `DaviesActivityModel` re-expressed under the new protocol.
The ionic-strength proxy is the Debye-Hückel insight: for dilute electrolytes, I
captures all the relevant composition information. Non-ionic species (z=0) get
γ_i = 1.0 automatically.

**`SITLiquidModel`** — Specific Ion Interaction theory. Higher accuracy at
elevated ionic strength (> ~0.5 M). Also ionic-strength-based but with
species-pair interaction coefficients ε(i,j). The existing `SITActivityModel`
re-expressed under the new protocol.

**`NRTLLiquidModel`** (future) — Non-Random Two Liquid. Uses the full `x_mol`
vector for all species. Required for organic liquid mixtures or concentrated
non-electrolyte systems where ionic strength is not the relevant variable.

### 2.4 Relationship to existing `ActivityModel`

The existing `ActivityModel` protocol (`gamma(z, I_molL, T_K) → float`) is
single-species and ionic-strength-based. It is retained internally in the
speciation engine where per-ion calls at known ionic strength remain the
natural interface. `LiquidPhaseModel` is the outward-facing protocol used by
`ThermoFramework` and partition models. The Davies and SIT implementations
satisfy both protocols.

### 2.5 Module home

`LiquidPhaseModel` and its implementations move to `src/thermo/`. The existing
`src/chemical_equilibrium/activity_models.py` retains `ActivityModel` for backward
compatibility but imports from `thermo/`. `GasEOS` stays in `src/equilibria/`.

---

## 3. ThermoFramework update

### 3.1 Current state

```python
@dataclass(frozen=True)
class ThermoFramework:
    activity_model: str = "davies"    # string name
    use_activity: bool = False
    standard_T_K: float = 298.15
    standard_P_atm: float = 1.0
```

The string is resolved to an `ActivityModel` object inside the speciation engine
via `make_activity_model()`. Nothing else can obtain the live object.

### 3.2 Proposed state

```python
@dataclass(frozen=True)
class ThermoFramework:
    liquid_activity: LiquidPhaseModel = field(default_factory=IdealLiquidModel)
    gas_eos: GasEOS                   = field(default_factory=IdealGasEOS)
    standard_T_K: float = 298.15
    standard_P_atm: float = 1.0
    # T-correction helpers (pKa_at_T, Kw_at_T) unchanged
```

**Common pre-built instances** (replaces the existing `THERMO_DAVIES` / `THERMO_IDEAL`):

```python
THERMO_IDEAL  = ThermoFramework(liquid_activity=IdealLiquidModel(),
                                gas_eos=IdealGasEOS())

THERMO_DAVIES = ThermoFramework(liquid_activity=DaviesLiquidModel(),
                                gas_eos=IdealGasEOS())

THERMO_DAVIES_PR = ThermoFramework(liquid_activity=DaviesLiquidModel(),
                                   gas_eos=PengRobinsonEOS(BIOGAS_SPECIES))
```

### 3.3 Role in the CV

One `ThermoFramework` instance is attached to each `ControlVolume`. Both the
speciation engine and every partition model registered on that CV are
constructed from it. Thermodynamic consistency is enforced structurally —
it is impossible to have the speciation engine using Davies while a partition
model uses ideal, because both draw from the same `ThermoFramework`.

```
ControlVolume(thermo=THERMO_DAVIES)
  ├── reaction_system.engine    ← uses thermo.liquid_activity (Davies)
  │                               produces: I, pH → written to phase.speciation
  └── internal_interfaces
        "CO2" → HenryPartition(..., thermo=cv.thermo)  ← uses thermo.liquid_activity
        "H2O" → RaoultPartition(thermo=cv.thermo)
        "CH4" → HenryPartition(..., thermo=THERMO_IDEAL)  ← override: CH4 is non-ionic
```

Per-interface overrides are permitted when a species is known to be non-ionic
(γ_i = 1) and the default `cv.thermo` would apply an unnecessary correction.
This is an efficiency choice, not a thermodynamic inconsistency, because ideal
(γ=1) is correct for non-ionic species regardless of ionic strength.

---

## 4. PartitionModel hierarchy

### 4.1 Existing single-species PartitionModel

The current `PartitionModel` protocol is species-local: it computes the
equilibrium distribution of **one species** given its total moles, the two
phase capacities, temperature, and an ionisation fraction `alpha`:

```python
class PartitionModel(Protocol):
    def equilibrium_a_moles(self, n_total, capacity_a, capacity_b,
                             T_K, alpha=1.0) -> float: ...
    def partition_ratio(self, capacity_a, capacity_b,
                        T_K, alpha=1.0) -> Optional[float]: ...
```

`partition_ratio() → float` enables the analytical exponential step (exact for
linear models). `partition_ratio() → None` signals non-linearity; the ODE
instantaneous-rate path is used instead. This capability-detection mechanism
already exists and is retained unchanged.

This protocol is correct and sufficient for Henry's law (ideal or
activity-corrected), `RaoultPartition`, and single-species Ksp. It covers the
majority of bioprocess partition phenomena.

### 4.2 MultispeciesPartitionModel

A parallel protocol for phenomena where the equilibrium of species i depends on
the concentrations of other species j in the same phase:

```python
class MultispeciesPartitionModel(Protocol):
    """Phase-global partition model for composition-dependent equilibria.

    Used when mixing rules, competitive binding, or non-ideal activity
    models require the full phase state to compute any single species'
    equilibrium distribution.
    """

    def equilibrium_all_a_moles(
        self,
        phase_a: Phase,
        phase_b: Phase,
        T_K: float,
    ) -> Dict[str, float]:
        """Equilibrium moles in phase 'a' for all coupled species.

        Receives the full state of both phases. Returns a dict mapping
        species_id → n_a_eq. Species not in the return dict are not
        handled by this model.
        """
        ...
```

There is no `partition_ratio` concept here: the composition coupling means no
single-species analytical step exists.

### 4.3 Trigger condition

| Model type | Use `PartitionModel` | Use `MultispeciesPartitionModel` |
|---|---|---|
| Henry's law (ideal gas) | ✓ | — |
| Henry's law + Davies activity | ✓ (uses I from phase.speciation) | — |
| Raoult's law (H₂O) | ✓ | — |
| Simple Ksp (single driving ion) | ✓ | — |
| Henry + Peng-Robinson gas | ✓ if species are non-interacting | Use multi if k_ij ≠ 0 |
| Competitive Langmuir adsorption | — | ✓ (shared sites couple all adsorbates) |
| NRTL liquid activity | — | ✓ (γ_i depends on all x_j) |
| Ksp with multi-ion stoichiometry | — | ✓ (IP couples multiple ions) |
| PR gas + NRTL liquid (full VLE) | — | ✓ |

**Practical rule:** use `PartitionModel` when species are thermodynamically
independent at operating conditions. Use `MultispeciesPartitionModel` when mixing
rules, competitive binding, or non-ideal composition-dependent activity links
them. For bioprocess at 1 atm with dilute solutions, `PartitionModel` covers
essentially all cases.

### 4.4 MultispeciesVLEPartition

The first `MultispeciesPartitionModel` implementation, combining a
`LiquidPhaseModel` and a `GasEOS`:

```python
@dataclass(frozen=True)
class MultispeciesVLEPartition:
    """VLE partition using composition-dependent EOS on both sides.

    Uses ThermoFramework.gas_eos for gas fugacities and
    ThermoFramework.liquid_activity for liquid activity corrections.
    Warranted when gas-phase mixing rules couple species (PR at high pressure)
    or when liquid non-ideality is composition-dependent (NRTL).
    """
    thermo: ThermoFramework
    kH_ref: Dict[str, float]    # Henry solubility per species
    dlnH:   Dict[str, float]    # d(ln kH)/d(1/T) per species

    def equilibrium_all_a_moles(self, phase_a, phase_b, T_K):
        # 1. Get full gas fugacities from gas_eos (accounts for mixing rules)
        # 2. Get activity coefficients from liquid_activity (accounts for composition)
        # 3. Solve coupled VLE: φ_i × y_i × P = γ_i × kH_i(T) × x_i for all i
        ...
```

### 4.5 RaoultPartition

Raoult's law applies to the solvent (water) in an ideal solution:

```
p_w = x_w × P_sat(T)
```

Henry's law applies to dilute solutes. The naming distinction matters
physically: Raoult and Henry are the two limits of the same underlying fugacity
equality, applicable in opposite concentration regimes.

```python
@dataclass(frozen=True)
class RaoultPartition:
    """Raoult's law gas-liquid partition for the solvent (H₂O).

    p_w = x_w × P_sat(T).  Correct for the solvent in an ideal or
    near-ideal solution; use HenryPartition for dilute solutes.
    """
    T_ref: float = 298.15
    P_sat_ref_atm: float = 0.03169  # P_sat(H₂O) at 25 °C
    dlnPsat: float = ...            # d(ln P_sat)/d(1/T) from Clausius-Clapeyron
```

`RaoultPartition` implements the single-species `PartitionModel` protocol —
Raoult's law is species-local (x_w appears alone on the liquid side) and linear,
so `partition_ratio() → float` is available and the analytical step applies.

---

## 5. Per-species EOS flexibility

### 5.1 What the architecture supports

The `PhaseInterface`-per-species registration pattern already supports distinct
partition models per species:

```python
ControlVolume(
    thermo=THERMO_DAVIES,
    transfer_models={
        "H2O":  EquilibriumTransferModel(RaoultPartition()),
        "CO2":  KineticTransferModel(HenryPartition(kH_ref=..., dlnH=...)),
        "NH3":  KineticTransferModel(HenryPartition(kH_ref=..., dlnH=...)),
        "CH4":  KineticTransferModel(HenryPartition(kH_ref=..., dlnH=...,
                                     thermo=THERMO_IDEAL)),  # non-ionic: no activity
    }
)
```

Per-species Henry constants (`kH_ref`, `dlnH`) and Setschenow salting-out
coefficients live in each `PartitionModel` instance. The EOS family
(`liquid_activity`, `gas_eos`) lives in `ThermoFramework` and is shared. The
two levels of variation are cleanly separated.

### 5.2 Quantities that couple "independent" species

Even with per-species `PartitionModel` instances, four quantities couple the
equilibrium of different species and prevent fully independent solution:

**Ionic strength I** — `I = 0.5 Σ z_i² C_i`. Every Davies-corrected Henry
constant depends on what every ion is doing. Shifting CO₂ from gas to liquid
changes carbonate speciation, changes I, changes γ(Ca²⁺), changes the Ksp
driving force.

**Total gas pressure P** — `P = Σ p_j`. The Raoult condition for H₂O
(`p_w = x_w × P_sat`) and every Henry condition (`p_i = y_i × P`) share the
same P. Each species' equilibrium condition references the same total.

**Gas mole fractions y_j** — PR-EOS mixing rules: `φ_i^PR(T, P, y)` depends
on all `y_j` through the van der Waals combining rule. The fugacity of CO₂ is
affected by the CH₄ mole fraction even if CH₄ uses ideal-gas partitioning,
because CH₄ appears in the mole fraction denominator.

**Mass balance** — total moles of each element are conserved across all phases.
Moving NH₃ to the gas phase changes the liquid nitrogen balance, which changes
the NH₄⁺/NH₃ ratio, which changes pH, which changes the CO₂/HCO₃⁻ split.

### 5.3 SNIA: when per-species independence holds

Under SNIA (the current default), coupling is broken by sequencing:

```
Step 1 — Speciation engine: computes I, pH, alpha(NH₃), alpha(CO₂)
          using current gas composition as fixed (pre-step state)

Step 4 — Partition models run sequentially:
          "CO2": uses I from step 1, current y_CO₂, P from current gas
          "NH3": uses alpha from step 1, updated P (after CO₂ ran)
          "H2O": uses updated P (after CO₂ and NH₃ ran)
```

The ordering within step 4 introduces a within-step O(dt) splitting error
between species — on top of the SNIA O(dt) error between steps 1 and 4.

For bioprocess timescales (dilute solutions, kLa ≪ 500 h⁻¹, 1 atm), both
errors are negligible and per-species `PartitionModel` instances with SNIA
ordering are both correct and efficient.

### 5.4 Simultaneous solve: one EOS family per phase

When all instantaneous equilibria are gathered into a single NR system (the
Layer 1 coupled solve, see §6), the Jacobian for mixed EOS choices becomes:

```
∂r_CO₂/∂y_CH₄ :  through φ_CO₂^PR(y) mixing rules — non-zero if PR
∂r_CO₂/∂C_Ca  :  through γ_CO₂(I) ionic strength — non-zero if Davies
∂r_Ksp/∂C_CO₂ :  through γ(Ca²⁺)γ(CO₃²⁻) — non-zero if Davies
```

Mixed model families (some species PR, some ideal) are not thermodynamically
incorrect, but they complicate Jacobian assembly: the PR block spans all gas
species regardless of which are labelled "PR" because mole fractions always sum
to 1.

The practical recommendation (used in Aspen, PHREEQC, and HYSYS) is:
**one EOS per phase; per-species variation through parameters, not model family**.
Species-specific constants — `kH_ref`, `dlnH`, Setschenow coefficient `k_s`,
binary interaction parameter `k_ij` — capture per-species thermodynamic
variation within a common EOS framework. Setting `k_ij = 0` for non-interacting
PR pairs recovers ideal-gas behaviour for those species without changing the
model family.

`ThermoFramework` holding one `LiquidPhaseModel` and one `GasEOS` per CV is
therefore not just a cleanliness preference — it is the natural consistency unit
for a simultaneous solve.

---

## 6. Gas-liquid-solid coupling

### 6.1 The three constraint sets

A CV with gas, liquid, and solid phases has three simultaneous equilibrium
conditions, all coupled through the liquid composition:

```
Acid-base (Layer 1, NR)         already solved simultaneously in NRChemicalEquilibriumEngine
  + ionic activity corrections via LiquidPhaseModel

Gas-liquid (Layer 1 gap)        currently solved after speciation (SNIA)
  ← depends on I from acid-base solve
  → changes liquid composition → changes I → couples back to acid-base

Liquid-solid / Ksp (wrong layer) currently inside NRChemicalEquilibriumEngine, not PhaseInterface
  ← depends on I from acid-base solve
  → changes ionic species → changes I → couples to everything
```

The MASS_EXCHANGE Layer 1 gap (`EquilibriumTransferModel` runs after speciation
rather than inside the algebraic solve) becomes more significant when activity
corrections are large: the Ksp condition for CaCO₃ depends on
`γ(Ca²⁺) × γ(CO₃²⁻)` computed from `I`, which depends on the gas-liquid CO₂
partition, which was computed using `I` from before the partition ran.

### 6.2 Severity by system type

| System | Coupling severity | SNIA adequate? |
|---|---|---|
| Dilute fermentation broth, 1 atm | I < 50 mM, kLa < 200 h⁻¹ | Yes |
| Saline media (~500 mM) | Setschenow correction ~10–20% | Usually yes |
| Seawater-scale, CO₂ system | Activity correction > 30% | Borderline |
| Concentrated slurry + Ksp + high kLa | Strong three-way coupling | No — needs Layer 1 |
| High-pressure biogas (> 5 atm) | PR φ_CO₂ departs from ideal | Depends on accuracy target |

### 6.3 Layer 1 with activity-corrected EOS

When the Layer 1 gap is closed (gas-liquid and Ksp constraints folded into the
NR tableau alongside acid-base), the activity-corrected driving forces are
evaluated inside the NR iteration. This is the correct treatment: I is updated
at each Newton step as species redistribute across phases, and convergence
guarantees global thermodynamic consistency.

`NRChemicalEquilibriumEngine` already does exactly this for the acid-base + Davies
coupling — activity corrections are inside the NR loop. Extending the tableau to
include gas-liquid and Ksp rows inherits this property without additional
structural work, provided the `LiquidPhaseModel` exposes its Jacobian
contribution `∂γ_i/∂C_j` (or equivalently `∂γ_i/∂I × ∂I/∂C_j`). For the
ionic-strength-based models (Davies, SIT) this derivative is straightforward;
for NRTL it requires the full `∂γ_i/∂x_j` matrix.

This connection is the primary reason to clean up the thermodynamic layer before
closing the Layer 1 gap — the gap closure work is easier and more correct when
`LiquidPhaseModel` is a first-class protocol rather than a string embedded in the
speciation engine.

---

## 7. Module locations

| Component | Location | Status | Notes |
|---|---|---|---|
| `LiquidPhaseModel` protocol | `src/thermo/liquid_phase_model.py` | ✓ shipped | Symmetric with `GasEOS`; `@runtime_checkable` |
| `IdealLiquidModel` | `src/thermo/liquid_phase_model.py` | ✓ shipped | Replaces implicit ideal default; absent → γ=1.0 |
| `DaviesLiquidModel` | `src/thermo/liquid_phase_model.py` | ✓ shipped | Satisfies `LiquidPhaseModel` and `ActivityModel` |
| `SITLiquidModel` | `src/thermo/sit_liquid_model.py` | ✓ shipped | Satisfies both protocols; `compute_gammas()` for BC |
| Water property helpers | `src/thermo/water_properties.py` | ✓ shipped | Moved from `src/chemical_equilibrium/`; re-exported for BC |
| `DaviesActivityModel` / `SITActivityModel` | `src/chemical_equilibrium/activity_models.py` / `sit.py` | ✓ aliases | `= DaviesLiquidModel` / `= SITLiquidModel`; BC |
| `ThermoFramework` | `src/thermo/framework.py` | ✓ updated | `liquid_activity: LiquidPhaseModel`; `@property` BC shims |
| `THERMO_IDEAL` / `THERMO_DAVIES` | `src/thermo/framework.py` | ✓ shipped | Pre-built instances |
| `GasEOS` / `IdealGasEOS` | `src/equilibria/vle.py` | unchanged | Gas-side counterpart |
| `PengRobinsonEOS` | `src/equilibria/peng_robinson.py` | unchanged | Non-ideal gas; not yet wired to `ThermoFramework` |
| `PartitionModel` protocol | `src/chemistry/partition.py` | unchanged | — |
| `HenryPartition` | `src/chemistry/partition.py` | ✓ updated | Optional `thermo` field; `kH_eff = kH / γ_i` |
| `MultispeciesPartitionModel` protocol | `src/chemistry/partition.py` | ✓ shipped | `equilibrium_all_a_moles()` for coupled VLE |
| `MultispeciesVLEPartition` | `src/chemistry/partition.py` | ✓ shipped | IdealGasEOS path; PR path deferred |
| `RaoultPartition` | `src/chemistry/partition.py` | ✓ shipped | Clausius-Clapeyron P_sat for H₂O |

---

## 8. Design decisions — resolved 2026-06-30

### 8.1 LiquidPhaseModel base signature — Option A

`gamma_all(x_mol, T_K, *, charge) → Dict[str, float]` is the sole base method.
Davies/SIT compress `x_mol` to ionic strength internally; NRTL uses the full
composition. One protocol, one call pattern for all callers — symmetric with
`GasEOS.partial_pressures_atm(n_gas_mol: Dict, T_K, V_L) → Dict`.

**`x_mol` unit convention:** concentrations in **mol/L** (VLsim's native unit).
Implementations are responsible for internal conversion (e.g., Davies converts
mol/L to molality before applying the Debye-Hückel equation). Callers always
pass mol/L and do not convert.

**Backward compatibility:** `ActivityModel.gamma(z, I_molL, T_K) → float` is
retained as an internal-only protocol within `src/chemical_equilibrium/`. `DaviesLiquidModel`
and `SITLiquidModel` satisfy both protocols: `LiquidPhaseModel` for external
callers (ThermoFramework, partition models) and `ActivityModel` for the NR
engine's inner loop, which does not need to change.

**Future models:** NRTL, eNRTL, UNIQUAC, Pitzer, COSMO-RS, and ML-based models
all work correctly with `gamma_all`. Non-electrolyte models ignore the `charge`
kwarg. A pressure parameter is absent; if high-pressure liquid-phase corrections
become in-scope, add `*, P_atm: float = 1.0` as a defaulted kwarg without
breaking existing implementations.

### 8.2 GasEOS — one per ThermoFramework (per CV)

Confirmed correct. PR with `k_ij = 0` for non-interacting pairs gives
φ_i ≈ 0.98–1.0 at 1 atm — negligible departure from ideal. At pressures where
PR matters (> 3–5 atm biogas), all gas species should use PR uniformly, which
is exactly what one EOS per CV provides. Species-specific variation comes through
parameters (`kH_ref`, `dlnH`, `k_ij`), not separate EOS instances.

### 8.3 EquilibriumScope vs algebraic_species() — algebraic_species() only

`EquilibriumScope` is dropped. `algebraic_species() → FrozenSet[str]` on
`ChemicalEquilibriumEngineProtocol` is the sole ownership mechanism.
`cv.advance()` uses it to skip partition models for species the engine already
owns:

```python
engine_owned = engine.algebraic_species()
for species_id, model in self._partition_models.items():
    if species_id in engine_owned:
        continue
    model.run(...)
```

`MASS_EXCHANGE_ARCHITECTURE.md` §12 Q8 and the `EquilibriumScope` proposal in
§13.3 are superseded by this decision. If phenomenon-category introspection is
needed in future, a separate optional `handled_phenomena() → FrozenSet[str]`
mixin can be added without replacing `algebraic_species()`.

### 8.4 Jacobian contribution from LiquidPhaseModel — deferred

`LiquidPhaseModel` does not expose a Jacobian method in this phase. The
`∂γ_i/∂C_j` contribution is only needed when gas-liquid and Ksp constraints are
folded into the NR tableau (Layer 1 gap closure). At that point, add a
`DifferentiableLiquidModel` extension protocol — parallel to `SplitJacobianCapable`
on the engine — exposing `jacobian_dgamma_dx(x_mol, T_K, *, charge) → ndarray`.
For Davies/SIT: `∂γ_i/∂C_j = (∂γ_i/∂I)(z_j²/2)` (cheap analytically).
For NRTL: full `∂γ_i/∂x_j` matrix required.

### 8.5 Jacobian contribution from GasEOS — deferred, not yet recorded elsewhere

The gas-side counterpart to §8.4, noted in the 2026-07-01 constraint-family
discussion (`MASS_EXCHANGE_ARCHITECTURE.md` §14) but not previously written
down. `GasEOS` does not expose a Jacobian method. `∂φ_i/∂y_j` is needed for a
true white-box fold of non-ideal gas-liquid rows — trivially zero for
`IdealGasEOS` (`φ_i = 1` identically), real but standard for `PengRobinsonEOS`
(closed-form cubic-EOS fugacity-coefficient derivatives from the mixing
rules). Not urgent: VLsim's near-term bioprocess use case stays close to
ideal-gas at ~1 atm, where this term vanishes. Add a `DifferentiableGasEOS`
extension protocol, parallel to `DifferentiableLiquidModel` in §8.4, when
`PengRobinsonEOS` is actually wired into folded gas-liquid rows.

### 8.6 Water activity (osmotic coefficient) — deferred, follow-up identified 2026-07-02

Found during `LAYER1_GAP_CLOSURE` CP4 (the Raoult/H2O fold): `a_H2O = 1` is
assumed everywhere in this codebase, for every reaction that references
water as a participant — not a new assumption CP4 introduced, but a
pre-existing one (`_SOLVENT_IDS` excludes H2O from `NRTableau`'s graph
entirely, structurally, before any activity model is consulted) that CP4's
Raoult fold made newly visible rather than newly true.

**This assumption is model-invariant, not just unimplemented.** Verified
directly: with `activity_model="davies"` active, H2O still never appears in
`tableau.masters`/`secondaries`/`master_charges` — the exclusion happens at
graph-construction time, upstream of where Davies/SIT could act. And even if
H2O *were* pushed through `DaviesLiquidModel.gamma(z=0, I)`, it returns
exactly `1.0` at every `I` tested (0.01–6.0 mol/L) — not because Davies
correctly evaluates water's activity as negligibly different from 1 at high
`I` (it isn't, physically), but because Davies/SIT are `z²`-driven ionic-strength theories
with no term for a neutral solute or solvent at all. Turning on activity
corrections for solutes does not partially correct water's own activity;
it's an entirely separate, unimplemented quantity.

**Where the approximation breaks down:** roughly the same regime where
Davies itself stops being valid (`I ≳ 0.5 mol/kg`) — concentrated
electrolytes (seawater `a_H2O ≈ 0.98`, brines lower), and, directly relevant
to this codebase's AD use case, **high-solids/dry anaerobic digestion**
(> ~15–20% TS), a recognized methanogen-inhibition factor in the AD
literature that this codebase cannot currently represent even in principle.
Concentrated non-ionic broths (high-gravity fermentation) show the same
effect via ordinary Raoult's-law crowding, independent of ionic strength.

**What a correct treatment needs:** not a `γ_H2O` correction (wrong
quantity) but an **osmotic coefficient** `φ`, the standard Pitzer/SIT
quantity for solvent activity: `ln(a_w) = −φ·M_w·Σm_i` (approximately). In a
self-consistent Pitzer/SIT formalism, `φ` is derivable from the *same*
ion-interaction parameters already used for solute `γ_i`, via the
Gibbs-Duhem relation — so this is not an unrelated model, just the half of
SIT/Pitzer (`SITLiquidModel` in particular) that was never implemented,
only the solute-`γ` half. Fixing it requires two coupled pieces, neither
attempted here: (1) water becoming a real tracked NR master with a finite
total (the "Option A" design explicitly deferred in CP4 — see
`LAYER1_GAP_CLOSURE.md`'s CP4 discussion — since without a tracked total
there's no `x_H2O`/`I` to feed an osmotic correction into anyway), and (2)
an osmotic-coefficient method on `LiquidPhaseModel` (or a new
`DifferentiableLiquidModel`-style extension) implementing the Gibbs-Duhem
consistency with whatever `γ_i` parameters are already in use.

---

## 9. Cross-references

- [`MASS_EXCHANGE_ARCHITECTURE.md`](MASS_EXCHANGE_ARCHITECTURE.md) —
  transport architecture; `PartitionModel` as the bridge; Layer 1 gap (§10.2);
  `RaoultPartition` (§6.2); open questions §12 Q3 and Q8 addressed here
- [`SOLVER_ARCHITECTURE.md`](SOLVER_ARCHITECTURE.md) —
  Layer 1/2/3 solver coupling; the simultaneous solve that §6.3 requires
- [`CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md`](CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md) —
  engine/solver split; black/gray/white-box engine protocols; `EquilibriumResult`;
  z-update strategies; `GrayBoxEngineProtocol.jacobian_dz_dy()` relates to §8.4;
  `SpeciationEngine` → `ChemicalEquilibriumEngine` rename (§18)
- [`src/thermo/framework.py`](../../src/thermo/framework.py) —
  `ThermoFramework` (shipped); `liquid_activity: LiquidPhaseModel` + `gas_eos`
- [`src/thermo/liquid_phase_model.py`](../../src/thermo/liquid_phase_model.py) —
  `LiquidPhaseModel`, `IdealLiquidModel`, `DaviesLiquidModel` (all shipped)
- [`src/thermo/sit_liquid_model.py`](../../src/thermo/sit_liquid_model.py) —
  `SITLiquidModel`, `SIT_EPSILON`, `ION_CHARGES` (shipped)
- [`src/chemical_equilibrium/activity_models.py`](../../src/chemical_equilibrium/activity_models.py) —
  `ActivityModel` protocol; backward-compat aliases `DaviesActivityModel = DaviesLiquidModel`
- [`src/equilibria/vle.py`](../../src/equilibria/vle.py) —
  `GasEOS` protocol and `IdealGasEOS`; gas-side counterpart to `LiquidPhaseModel`
- [`src/equilibria/peng_robinson.py`](../../src/equilibria/peng_robinson.py) —
  `PengRobinsonEOS`; not yet wired to `ThermoFramework` (deferred)
- [`src/chemistry/partition.py`](../../src/chemistry/partition.py) —
  `HenryPartition` (activity-corrected), `RaoultPartition`, `MultispeciesVLEPartition`
  (all shipped); `MultispeciesPartitionModel` protocol
