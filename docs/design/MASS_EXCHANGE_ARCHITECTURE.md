# Mass Exchange Architecture — Design Notes

*Status: exploration / pre-implementation. This document records a design discussion
and is intended to be refined iteratively before any code changes are made.*

---

## 1. Motivation

This document captures a design exploration of how mass exchange and transport
phenomena are (and should be) conceptualised in VLsim, prompted by the question:

> Is there a clean separation between physical transport phenomena (advection,
> diffusion, gas absorption, solid precipitation, etc.) and chemical phenomena
> (equilibrium reactions, kinetically driven reactions)?

The central finding is that the current architecture is largely well-structured
but has several areas where the conceptual map and the type system diverge — most
notably around precipitation, evaporation, and the rigid topology assumption
baked into the `PhaseInterface` / `CVLink` split.

---

## 2. Conceptual Taxonomy: Two Axes

All mass exchange phenomena can be located on two orthogonal axes:

### Axis 1 — Domain

| Domain | Meaning | Changes chemical identity? |
|---|---|---|
| **Transport** | Moves species between locations or phases | No — species identity preserved |
| **Chemistry** | Transforms species via bond formation/breaking | Yes |

### Axis 2 — Topology (where are the two endpoints?)

| Topology | Current protocol |
|---|---|
| Same CV, different phases | `PhaseInterface` |
| Different CVs | `CVLink` / `InterzonalFlow` (see §4) |
| External system boundary | `ExternalBoundary` |

These two axes are **independent**. The driving force (thermodynamic, spatial,
electrical) is a property of the phenomenon, not of the topology. This matters
— see §7.

---

## 3. The `PartitionModel` as the Bridge

`PartitionModel` is the only object that straddles both domains. It encodes a
chemical equilibrium criterion (what the thermodynamics requires) but is consumed
by transport code (how and how fast the system moves toward it):

```
Chemistry domain:
  EquilibriumReaction  ──►  HenryPartition / KspPartition / LangmuirPartition
                             (encodes WHAT the equilibrium is)

Transport domain:
  TransferModel (Kinetic or Equilibrium)  ──►  PhaseInterface
                             (encodes HOW FAST equilibrium is approached)
```

The invariant:
- A `PartitionModel` **never moves moles**.
- A `TransferModel` **never knows chemistry**.
- A `PhaseInterface` **wires them together**.

The thermodynamic layer beneath `PartitionModel` — activity models, equations of
state, and the `ThermoFramework` that coordinates them — is detailed in
[`THERMODYNAMIC_MODEL_ARCHITECTURE.md`](THERMODYNAMIC_MODEL_ARCHITECTURE.md).

This separation is already correctly implemented for gas-liquid transfer
(`KineticTransferModel + HenryPartition`) and is the design target for all other
phase-transfer phenomena.

---

## 4. Terminology Recommendations

The current name `CVLink` describes the *graph structure* of the code, not the
physical phenomenon. Recommended rename:

| Current | Recommended | Rationale |
|---|---|---|
| `CVLink` (protocol) | `InterzonalFlow` | "interzonal" is standard in compartmental modelling; "flow" encodes what is computed |
| `AdvectiveLink` | `AdvectiveFlow` | consistent naming |
| `DiffusiveLink` | `DispersiveFlow` | more accurate: the coefficient is an effective dispersion coefficient (turbulent eddies + molecular diffusion), not purely Fickian |

The `kLa` parameter on `DiffusiveLink` / `DispersiveFlow` is also misleading
(borrowed from gas-liquid notation); `k_disp` or `D_eff` would be clearer.

---

## 5. Transport Phenomena Map: Existing vs Preferred

### 5.1 Inter-zone transport (InterzonalFlow)

All produce `{species: flux_mol_per_h}` via `compute_flow(cvs, dt_h)`.

| Phenomenon | Driving force | Formula | Existing | Status |
|---|---|---|---|---|
| **Advection** | Bulk flow rate Q | `Q × C_i(src)` | `AdvectiveLink` | ✅ |
| **Dispersion** | Concentration gradient | `k_i × V_eff × ΔC_i` | `DiffusiveLink` | ✅ (rename) |
| **Gravitational settling** | Particle buoyancy | `v_s_i × A × C_i` | — | ❌ not yet |
| **Level-triggered overflow** | Liquid level > weir | `Q(h) × C_i` | — | ❌ not yet |
| **Electromigration** | Electric potential gradient | `z_i × u_i × C_i × ΔV/L` | — | ❌ not yet |
| **Nernst-Planck (full)** | ΔC + ΔV coupled | `-D_i × ΔC_i/L − u_i × C_i × ΔV/L` | — | ❌ not yet |
| **Thermodynamic (Henry's law)** | Chemical potential | `kLa × V × (C* − C)` via `PartitionModel` | — | ❌ **gap** (see §7) |

### 5.2 External boundary transport (ExternalBoundary)

| Phenomenon | Existing class | Status |
|---|---|---|
| Liquid feed (advective in) | `LiquidFeed` | ✅ |
| Gas feed (advective in) | `GasFeed` | ✅ |
| Liquid drain (advective out) | `LiquidDrain` | ✅ |
| Membrane permeation | `MembraneGasBoundary` | ✅ |
| Pressure relief (threshold) | `PressureReliefVent` | ✅ |
| Proportional gas outlet | `ProportionalGasOutlet` | ✅ |
| Water vapour saturation | `WaterVapourBoundary` | ⚠ misplaced (see §5.3) |

### 5.3 Cross-phase intra-CV transport (PhaseInterface)

All share the pattern: `TransferModel(PartitionModel)`.

| Phenomenon | Partition model | Transfer model | Existing | Status |
|---|---|---|---|---|
| Gas ↔ liquid (Henry, kinetic) | `HenryPartition` | `KineticTransferModel` | `transfer_models=` kwarg | ✅ preferred path |
| Gas ↔ liquid (Henry, legacy) | hardcoded | — | `KineticGasLiquidLink` | ⚠ retire |
| Gas ↔ liquid (equilibrium) | `HenryPartition` | `EquilibriumTransferModel` | `transfer_models=` kwarg | ✅ |
| Evaporation / condensation | `RaoultPartition` | `EquilibriumTransferModel` | `WaterVapourBoundary` | ⚠ wrong layer |
| Liquid ↔ solid (Ksp, equil.) | `KspPartition` | `EquilibriumTransferModel` | inside `NRChemicalEquilibriumEngine` | ⚠ wrong layer |
| Liquid ↔ solid (Ksp, kinetic) | `KspPartition` | `KineticTransferModel` | — | ❌ not possible yet |
| Liquid ↔ solid (adsorption) | `LangmuirPartition` / `FreundlichPartition` | `KineticTransferModel` | — | ❌ not yet |

### 5.4 Chemistry (changes species identity)

| Phenomenon | Existing | Status |
|---|---|---|
| Kinetic reactions | `KineticReaction` → `ReactionSystem` | ✅ |
| Acid-base / complexation equilibria | `EquilibriumReaction` → `BisectionChemicalEquilibriumEngine` | ✅ |
| Multi-component speciation | `EquilibriumReaction` → `NRChemicalEquilibriumEngine` | ✅ |
| Phase-equilibrium criterion | `PartitionModel` (Henry, future: Ksp, Langmuir) | ✅ design; partial impl |

---

## 6. The Two Misplaced Phenomena

### 6.1 Precipitation / dissolution (Ksp) — inside the speciation engine

The `NRChemicalEquilibriumEngine` currently acts as both:
1. Solver of the aqueous equilibrium constraint — **correct**
2. Arbiter of how much solid forms / dissolves — **correct for the criterion**
3. Holder of the solid-phase mole amounts — **wrong layer**

The speciation engine should report the saturation index (SI = log(IP/Ksp)) but
not move moles to `SolidPhase`. The mass transfer belongs in a `PhaseInterface`.

**Preferred design (as originally proposed — superseded, see §14):**
```python
# KspPartition implements PartitionModel
# Signals: IP > Ksp → precipitate; IP < Ksp + solid present → dissolve
cv = ControlVolume(
    phases={"liquid": liq, "solid": SolidPhase(...)},
    transfer_models={
        "CaCO3": EquilibriumTransferModel(
            KspPartition(Ksp=3.36e-9, stoichiometry={"Ca++": 1, "CO3--": 1})
        ),
        # Kinetic dissolution is then a one-word change:
        "CaCO3": KineticTransferModel(
            KspPartition(Ksp=3.36e-9, stoichiometry={"Ca++": 1, "CO3--": 1}),
            k_transfer=k_diss,
        ),
    }
)
```

> **Superseded 2026-07-01.** Moving Ksp out to `PhaseInterface` is no longer
> recommended. `EquilibriumReaction` (log-K) and Ksp are the same mass-action
> representation — the shipped `NR_PRECIPITATION_SPECIATION` already encodes
> Ksp as `log_K`. Extracting it while folding the mathematically identical
> Henry's-law case *in* would be inconsistent. See §14 for the resolution and
> the real (different) limitation this uncovered.

### 6.2 Evaporation / condensation — as ExternalBoundary

`WaterVapourBoundary` is an `ExternalBoundary` that adjusts gas-phase H₂O moles
to maintain saturation. Conceptually it is an **equilibrium gas-liquid
PhaseInterface** where the driving force is Raoult's law. Being an
`ExternalBoundary`:
- Does not decrement liquid-phase H₂O (water "appears from nowhere")
- Required `VentWaterLoss` as a patch, introducing an ordering dependency
- Treats the liquid as an infinite reservoir rather than a tracked phase

**Preferred design:**
```python
cv = ControlVolume(
    transfer_models={
        "H2O": EquilibriumTransferModel(RaoultPartition())
        # Raoult's law: P_H2O = x_H2O × P_sat(T)
        # Automatically closes the liquid ↔ gas H2O mass balance
    }
)
```

This would retire both `WaterVapourBoundary` and `VentWaterLoss`.

---

## 7. The Topology Rigidity Problem

### The assumption

The current architecture assumes:
- Thermodynamic driving forces → **intra-CV** (`PhaseInterface`)
- Spatial driving forces → **inter-CV** (`InterzonalFlow`)

This is a *convention*, not a physical law. The same Henry's law physics applies
whether gas and liquid are two phases in one CV or two separate CVs.

### Reasonable use cases that break this assumption

**Case A: Gas headspace CV + bulk liquid CV**
```
cv_headspace: phases = {"gas": GasPhase}
cv_bulk:      phases = {"liquid": LiquidPhase}
connection:   Henry's law transfer  ← currently no InterzonalFlow implementation for this
```
Motivation: model the headspace as a spatially distinct compartment; connect the
same headspace to multiple liquid CVs with different kLa values; model
stratification in the gas phase.

**Case B: Bulk liquid CV + biofilm CV**
```
cv_bulk:    phases = {"liquid": LiquidPhase}
cv_biofilm: phases = {"liquid": LiquidPhase, "solid": SolidPhase}
```
The outer surface connection (bulk → biofilm pore water) is `DispersiveFlow`
— already supported. Internal biofilm chemistry is handled by `cv_biofilm`'s
own reaction system. This use case is **already handled**.

**Case C: Liquid height tracking in fed-batch**

In a single-CV design, `V_gas = V_reactor − V_liquid` is automatic. In a
multi-CV design (separate gas and liquid CVs), a **geometric constraint** would
be needed to tie the two volumes together — this concept does not currently
exist. The single-CV design handles fed-batch geometry more naturally; separate
CVs for this purpose add complexity without benefit.

### The deeper insight

A transport link connects **two (CV, phase) pairs**. Whether those pairs are in
the same CV or different CVs should be secondary to the driving force model:

```
TransportLink
  source: (cv_key, phase_key)
  sink:   (cv_key, phase_key)
  model:  TransferModel(PartitionModel)   ← same regardless of topology
```

Under this framing:
- `PhaseInterface` = special case where both CV keys are identical
- `InterzonalFlow` = special case where CV keys differ

The driving force model is **orthogonal to topology**.

---

## 8. Electrokinetic Transport

Electrokinetics introduces a third class of driving force for `InterzonalFlow`:
the electric potential gradient. The Nernst-Planck equation shows that
electrokinetic transport is the superposition of existing and new terms:

```
J_i  =  −D_i × ΔC_i / L          (Fick diffusion → DispersiveFlow)
      −  u_i × C_i  × ΔV / L     (electromigration → new)
      +  C_i × v                  (advection, incl. electro-osmosis → AdvectiveFlow)
```

where `u_i = z_i × D_i × F / (RT)` is the signed ionic mobility.

### What is already present

| Requirement | Status | Location |
|---|---|---|
| `Species.charge` (`z_i`) | ✅ first-class field | `Species.charge` on every species |
| Charge balance tracking | ✅ per-step per-CV | `ConservationMonitor._charge_residual()` |
| Charge residual flagging | ✅ | `AccuracyMonitor.check_charge_residual()` |
| Ionic strength | ✅ computed | `liq.speciation["IonicStrength"]` |

### What is new

**Ionic mobility** (`u_i`, m²·mol/(J·s)) — the transport parameter specific to
electromigration. Currently no concept of per-species mobility exists. It is
linked to diffusivity via the Einstein relation but is most simply stored
explicitly on the link.

**Electric potential** (`φ`, V) as a field. Two implementation depths:
- *Imposed field* (simple): `delta_V_V: float` parameter on the link — sufficient
  for electrodialysis, electrokinetic remediation with applied voltage
- *Computed field* (advanced): `phi_V` as a zone state variable solved via the
  Poisson equation — requires a new coupled field solver

**New `InterzonalFlow` implementations:**
```python
ElectromigrationFlow:
    ionic_mobility: Dict[str, float]   # u_i per species (explicit)
    charge: Dict[str, int]             # z_i (from Species registry)
    delta_V_V: float                   # φ_source − φ_sink (V)
    gap_m: float                       # zone separation (m)
    # flux_i = z_i × u_i × C_i × delta_V_V / gap_m

NernstPlanckFlow:   # combined; one D_i serves both terms via Einstein relation
    diffusivity: Dict[str, float]      # D_i (m²/h)
    charge: Dict[str, int]
    delta_V_V: float
    gap_m: float
    T_K: float = 298.15
```

### Charge conservation coupling

When electromigration moves Na⁺ from zone A to zone B, zone A loses positive
charge and zone B gains it — deliberately and correctly. The current
`ConservationMonitor` treats per-zone charge drift as an error. A small protocol
extension is needed:

```python
# Orchestrator passes acknowledged inter-zone charge transfer
monitor.check_step(phases, acknowledged_charge_transfer=charge_transferred_mol)
```

At system boundaries, electrode (Faradaic) reactions close the electrical circuit
and are modelled as `KineticReaction` objects at the boundary CVs.

---

## 9. Full Updated Transport Map

```
MASS EXCHANGE PHENOMENA
│
├── TRANSPORT (species identity preserved)
│   │
│   ├── [A] InterzonalFlow — between (CV, phase) pairs in DIFFERENT CVs
│   │     │
│   │     ├── Spatial driving force
│   │     │     AdvectiveFlow        Q × C_i               [✅ exists]
│   │     │     DispersiveFlow       k_i × V_eff × ΔC_i    [✅ exists, rename]
│   │     │     SettlingFlow         v_s_i × A × C_i       [❌ not yet]
│   │     │     OverflowFlow         Q(level) × C_i        [❌ not yet]
│   │     │
│   │     ├── Thermodynamic driving force (GAP — see §7)
│   │     │     KineticInterzonalTransfer(HenryPartition)   gas CV ↔ liquid CV
│   │     │     KineticInterzonalTransfer(KspPartition)     liquid CV ↔ solid CV
│   │     │     EquilibriumInterzonalTransfer(...)          any phase pair
│   │     │
│   │     └── Electric potential driving force
│   │           ElectromigrationFlow  z × u × C × ΔV/L    [❌ not yet]
│   │           NernstPlanckFlow      D + ΔV coupled       [❌ not yet]
│   │
│   ├── [B] ExternalBoundary — one endpoint is an external reservoir
│   │     LiquidFeed, GasFeed, LiquidDrain                 [✅]
│   │     MembraneGasBoundary                               [✅]
│   │     PressureReliefVent, ProportionalGasOutlet         [✅]
│   │     WaterVapourBoundary                               [⚠ retire → PhaseInterface]
│   │
│   └── [C] PhaseInterface — between (CV, phase) pairs in the SAME CV
│         All use TransferModel(PartitionModel)
│         Gas ↔ Liquid   KineticTransferModel(HenryPartition)        [✅]
│         Gas ↔ Liquid   EquilibriumTransferModel(HenryPartition)     [✅]
│         Evaporation     EquilibriumTransferModel(RaoultPartition) [❌ not yet]
│         Liquid ↔ Solid  EquilibriumTransferModel(KspPartition)      [❌ not yet]
│         Liquid ↔ Solid  KineticTransferModel(KspPartition)          [❌ not yet]
│         Liquid ↔ Solid  KineticTransferModel(LangmuirPartition)     [❌ not yet]
│
└── CHEMISTRY (species identity changes)
      KineticReaction    → ReactionSystem                   [✅]
      EquilibriumReaction → BisectionChemicalEquilibriumEngine                [✅]
      Phase criterion    → PartitionModel                   [✅ design; partial impl]
```

---

## 10. Solver Coupling: When Separation Works and When It Doesn't

### 10.1 The timescale argument

Operator splitting (treating each phenomenon sequentially within a timestep) is
valid when phenomena have **well-separated timescales**. The DAE formulation
documented in `SOLVER_ARCHITECTURE.md` codifies this:

```
dy/dt = f(y, z)     ← differential variables: conserved totals, changed by kinetics
0     = g(y, z)     ← algebraic variables: equilibrium distributions, instantaneous
```

This separation is valid when equilibration timescales (microseconds for
acid-base, sub-seconds for Henry gas-liquid) are fast relative to kinetic
timescales (minutes to hours for bioprocess reactions). When this holds, the
stiff fast-equilibrium dynamics disappear from the ODE and only the slow kinetics
remain. When it does *not* hold — or when multiple "instantaneous" phenomena are
coupled to each other rather than just to slow kinetics — sequential splitting
introduces errors of order O(dt) per step.

### 10.2 Three nested layers of simultaneous solve

The architecture supports (or is designed to support) three distinct levels of
simultaneous solve, each handling a different coupling scope:

```
Layer 3 — Multi-CV ODE (MonolithicODE)
  ┌─────────────────────────────────────────────────┐
  │  dy_A/dt = f_A(y_A, z_A) + J_{A←B}(y_A, y_B)  │  ← inter-zone flux in combined RHS
  │  dy_B/dt = f_B(y_B, z_B) + J_{B←A}(y_A, y_B)  │
  │                                                 │
  │  Layer 2 — Per-CV DAE (compute_rhs)             │
  │  ┌──────────────────────────────────────────┐   │
  │  │  dy/dt  = f(y, z)    ← kinetics         │   │
  │  │  0      = g(y, z)    ← equilibria       │   │
  │  │  evaluated simultaneously at (y, z)     │   │
  │  │                                         │   │
  │  │  Layer 1 — Algebraic (speciation)       │   │
  │  │  ┌───────────────────────────────────┐  │   │
  │  │  │  0 = g(y, z)                     │  │   │
  │  │  │  all equilibria solved together  │  │   │
  │  │  │  via Newton-Raphson given y      │  │   │
  │  │  └───────────────────────────────────┘  │   │
  │  └──────────────────────────────────────────┘   │
  └─────────────────────────────────────────────────┘
```

**Layer 1 — Algebraic solve (speciation engine)**

All instantaneous equilibrium constraints within a CV are solved simultaneously
in a single Newton-Raphson system. This is not sequential — all acid-base
equilibria, complexation, and (eventually) precipitation constraints are rows of
the same g(y, z) = 0 system.

What *must* also live here: `EquilibriumTransferModel`. If gas-liquid
equilibrium is truly instantaneous, the Henry's law partition constraint is
another row of g(y, z) = 0 — coupled to the acid-base constraints because
`z_liq,CO2` appears in both. The current implementation runs `EquilibriumTransferModel`
*after* the speciation engine as a sequential correction. For CO₂/HCO₃⁻ where
gas-liquid and acid-base equilibria are coupled, they must be solved together for
consistency.

**Layer 2 — Per-CV DAE (`compute_rhs`)**

`compute_rhs()` evaluates kinetic reactions and kinetic `PhaseInterface` transfer
at the **same (y, z) snapshot** — eliminating the previous sequential splitting
between reactions and transfer. The DAE extension (white-box engine exposing
`∂g/∂z`, `∂g/∂y`) allows an implicit stiff integrator (BDF, Radau) to treat z
as algebraic variables in the combined system, removing step-size restriction
from the equilibrium timescale.

**Layer 3 — Multi-CV ODE (MonolithicODE)**

The `MonolithicODE` solver assembles a combined RHS across all CVs, making
inter-zone transport terms part of the simultaneous integration rather than
sequential pre-corrections. This resolves temporal lag errors when zone coupling
is tight (transfer rate comparable to reaction rate).

### 10.3 The reactive transport dilemma

The hardest coupling case is `InterzonalFlow` of species in **local chemical
equilibrium** within each zone. Consider CT (total dissolved inorganic carbon =
CO₂ + HCO₃⁻ + CO₃²⁻) diffusing between two liquid zones:

- CO₂(aq), HCO₃⁻, CO₃²⁻ have different diffusivities
- But they are in instantaneous equilibrium within each zone
- They cannot diffuse independently while maintaining local equilibrium

If `DispersiveFlow` carries CO₂(aq) and HCO₃⁻ separately with their own
coefficients, the zones leave equilibrium after each transport step — the
speciation re-solve then re-distributes, but the re-distribution does not
conserve what the transport step intended.

**The correct treatment** is to transport conserved totals (CT, the sum) with a
single effective diffusivity, then let each zone's speciation re-establish the
equilibrium split. The effective diffusivity of CT is:

```
D_eff,CT = Σ_i  (∂z_i/∂CT) × D_i
```

where `∂z_i/∂CT` is the fraction of a marginal unit of CT that appears as
species i at local equilibrium — exactly the column for CT in the gray-box
Jacobian `∂z/∂y`. This is the **direct substitution approach** to reactive
transport; it maintains thermodynamic consistency at the cost of computing the
Jacobian at each evaluation.

**Current architecture position:** `DispersiveFlow` operates on individual
species without awareness of which are in local equilibrium. For non-equilibrium
species this is exact. For equilibrium species it is a first-order-accurate
approximation. The gray-box engine protocol exposes what is needed for the D_eff
correction, but `DispersiveFlow` does not yet consume it.

### 10.4 What each combination requires

| Combination | Required solver path | Current status |
|---|---|---|
| Instantaneous equilibria (acid-base, complexation) | Layer 1 — simultaneous NR | ✅ |
| Instantaneous equilibria + gas-liquid transfer (Henry, Raoult) | Layer 1 — both in algebraic solve | ✅ closed by `LAYER1_GAP_CLOSURE` (tag `layer1-gap-closure-shipped`) — folded into the NR tableau alongside acid-base, ideal-gas case; `PengRobinsonEOS` non-ideal case remains §8.5-deferred |
| Instantaneous equilibria + solid-liquid transfer (Ksp) | Layer 1 — both in algebraic solve | ✅ closed — but via a nested active-set loop around the inner NR solve (`NRSpeciationEngine._solve_with_precipitation`), not a single flat Jacobian; accepted by design for single-driving-ion Ksp (§14.3). Multi-ion Ksp as one flat Jacobian is `MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md`'s territory, not this gap |
| Kinetic reactions + kinetic phase transfer | Layer 2 — `compute_rhs` same snapshot | ✅ |
| Stiff kinetics + fast equilibria (DAE) | Layer 2 — white-box engine + BDF/Radau | ⚠ designed; not yet implemented |
| Tight inter-zone kinetic coupling | Layer 3 — MonolithicODE | ✅ closed — `MonolithicODESolver` shipped 2026-06-12 (tag `monolithic-ode-shipped`); its `_build_rhs()` folds inter-CV transport rates directly into the same `dydt` as per-CV chemistry (`system_solver.py`), exactly the "combined RHS across all CVs" this row describes. This status was stale — corrected 2026-07-06, verified directly against the shipped code |
| Inter-zone transport of equilibrium species | Layer 3 + gray-box D_eff correction | ❌ approximate only — confirmed still open: `_build_rhs()`'s transport term moves each species independently at its own `kLa`/`Q`, no gray-box Jacobian/D_eff correction present. Distinct from the row above — MonolithicODESolver closes general tight coupling but not this species-in-local-equilibrium case |
| Instantaneous equilibrium inter-zone transfer | Cross-CV algebraic solve | ❌ no mechanism; rarely needed in practice |

### 10.5 Practical significance

For bioprocess timescales the timescale separation is generally sufficient and
operator splitting is accurate. The cases where it materially breaks down are:

1. **Very high kLa (> ~500 h⁻¹) with CO₂ / NH₃**: gas-liquid transfer timescale
   approaches the pH response timescale. Folding gas-liquid rows into the NR
   tableau's algebraic solve resolves this — closed by `LAYER1_GAP_CLOSURE`
   CP2, demonstrated concretely by a realistic AD biogas-sparging scenario
   where the old sequential (SNIA) path overshoots the simultaneous-fold
   answer by >20% in one macro-timestep.

2. **Tight inter-zone coupling at coarse timesteps**: one zone feeding another at
   a rate comparable to the reaction rate. MonolithicODE (Layer 3) resolves this.

3. **Reactive diffusion of equilibrium species at large pH gradients**: the D_eff
   correction matters when equilibrium species diffusivities differ significantly
   (CO₂ vs HCO₃⁻ differs ~1.5×) and zones are at substantially different pH. At
   similar pH the effect is small and the current approximation is adequate.

---

## 11. Priority Gaps

Ordered by design impact:

| Priority | Gap | What's needed |
|---|---|---|
| 1 | ~~Precipitation out of speciation engine~~ — **superseded 2026-07-01, see §14** | Ksp stays in-engine as a mass-action row alongside acid-base and (once folded) Henry; `KspPartition`/`PhaseInterface` remains available only for genuinely kinetic dissolution |
| 2 | Thermodynamic inter-zone transfer | `KineticInterzonalTransfer(PartitionModel)` — reuses existing models, new topology |
| 3 | ~~`EquilibriumTransferModel` folded into algebraic solve~~ — **closed 2026-07-03, `LAYER1_GAP_CLOSURE` CP1/CP2** | Gas-liquid (Henry) rows now fold into `g(y, z) = 0` alongside acid-base (§10.2); ideal-gas case only, `PengRobinsonEOS` remains §8.5-deferred |
| 4 | ~~Evaporation as PhaseInterface~~ — **closed 2026-07-03, `LAYER1_GAP_CLOSURE` CP4** | `RaoultEquilibrium` + `transfer_models=EquilibriumTransferModel(RaoultEquilibrium())`; `WaterVapourBoundary`/`VentWaterLoss` retired |
| 5 | Adsorption | `LangmuirPartition` + `FreundlichPartition` + wire to `SolidPhase` |
| 6 | Gravitational settling | `SettlingFlow` — new `InterzonalFlow` implementation |
| 7 | Electrokinetics (imposed field) | `ElectromigrationFlow` + `ConservationMonitor` extension |
| 8 | Reactive inter-zone transport (D_eff) | `DispersiveFlow` consumes gray-box `∂z/∂y` for effective diffusivity of equilibrium species |
| 9 | Geometric volume constraints | `VolumeConstraint` for multi-CV designs where volumes are coupled |
| 10 | `cv.advance()` engine-scope contract — **partially closed 2026-07-03, `LAYER1_GAP_CLOSURE` CP3** | `step_internal_transfer()` scope-filter fix shipped (`NRChemicalEquilibriumEngine.gas_liquid_species()`, a narrower pragmatic mechanism than the full `EquilibriumScope` declaration sketched in §13); construction-time overlap validation remains open |

---

## 12. Open Questions

1. **Unified transport link**: Should `PhaseInterface` and `InterzonalFlow` be
   unified under a single protocol that takes `(source_cv, source_phase,
   sink_cv, sink_phase)` with `same CV` as a degenerate case? Or is the
   execution-context difference (intra-step vs inter-step) a sufficient reason to
   keep them separate?

2. **Execution order for thermodynamic inter-zone transfer**: When a `KineticInterzonalTransfer(HenryPartition)` connects two CVs, the `alpha` correction (ionisation fraction) needed for the Henry's law driving force comes from the speciation engine. The orchestrator runs speciation *within* each `cv.advance()` call — so the alpha used for inter-zone flux would be one step stale. Is this acceptable, or does the orchestrator need to trigger a pre-flux speciation solve?

3. **KspPartition signature**: Precipitation involves multiple ions with a
   stoichiometric product (`[Ca²⁺]¹[CO₃²⁻]¹ = Ksp`). The `PartitionModel`
   protocol currently takes a single `n_total` and capacities. Does `KspPartition`
   need a different protocol, or can it read individual ion concentrations from
   the phase directly?
   > **Addressed in [`THERMODYNAMIC_MODEL_ARCHITECTURE.md`](THERMODYNAMIC_MODEL_ARCHITECTURE.md) §4.**
   > A `MultispeciesPartitionModel` protocol that receives the full phase state is
   > the recommended parallel track. `KspPartition` for multi-ion minerals is one
   > of its primary use cases. Resolve §8.1 of that document before implementation.
   > **Further addressed in §14** — multi-ion Ksp folding is also gated on
   > `NRTableau`'s one-master-per-component limitation, not just the
   > `PartitionModel` signature.

4. **Solid phase volume**: `SolidPhase` currently has `n_mol` but its volume is
   not dynamic. When solids precipitate, does `V_solid` grow and displace liquid?
   For dilute systems this is negligible; for dense slurries it matters.

5. **Electrokinetics — imposed vs computed potential**: Is the imposed-field
   (`delta_V_V` as a link parameter) sufficient for all near-term use cases, or
   is there a bioprocess application that requires a spatially resolved potential
   field (Poisson equation coupled to ion distribution)?

6. **Engine scope validation**: Should scope-overlap errors be raised at CV
   construction time (strict; catches misconfiguration early) or at solve time
   (flexible; allows runtime engine switching)? Construction-time validation is
   safer for the common case but conflicts with the ability to hot-swap engines
   for benchmarking.

7. **z-staleness across SIA iterations**: In SIA, `cv.advance()` is called
   multiple times per timestep. On the second iteration the phase state includes
   source-term corrections from the first iteration. Should the engine re-solve
   at the start of each iteration (correct z, more expensive) or reuse the z
   from the first pass (approximate, consistent with current SNIA ordering)?

8. **`EquilibriumScope` vs `algebraic_species()` — resolved 2026-06-30.**
   `algebraic_species() → FrozenSet[str]` on `ChemicalEquilibriumEngineProtocol`
   is the sole ownership mechanism. `EquilibriumScope` is dropped. `cv.advance()`
   uses `engine.algebraic_species()` to skip partition models for species the
   engine already owns. The `EquilibriumScope` proposal in §13.3 is superseded.
   See `THERMODYNAMIC_MODEL_ARCHITECTURE.md` §8.3 for full rationale.

9. **Surface complexation / reactive adsorption — deferred, not in scope for
   Layer 1 gap closure.** (2026-07-01 design discussion, not yet implemented.)

   When an adsorbed species undergoes further reaction on the surface itself
   (e.g. surface protonation `≡SOH ⇌ ≡SO⁻ + H⁺`, or two adsorbed species
   forming a surface complex), a plain `PartitionModel`/`MultispeciesPartitionModel`
   isotherm is insufficient — the surface needs to be treated as a phase with
   its own conserved site-total and its own mass-action reaction ladder,
   structurally identical to the CO₂ gas-liquid + carbonate-ladder coupling
   that motivates Layer 1 folding in the first place (§10.2): a species crosses
   a phase boundary via a partition/isotherm step, then reacts further within
   the destination phase, and the two must be solved together.

   `SolidPhase` can plausibly serve as the container for adsorbed surface
   species without needing a fourth `Phase` subclass — in many real systems
   (e.g. ferrihydrite surface complexation) site density scales with the
   amount of precipitated mineral present, so co-locating them is arguably
   correct physics, not just container reuse. Two things are still needed
   before this can be folded into a tableau, neither of which is in scope now:

   - **Species-role-based activity dispatch**, not just phase-membership-based.
     A bulk precipitate species has activity fixed at 1 by convention; an
     adsorbed surface species' activity analog is fractional site occupancy
     `θ_i` (possibly with an electrostatic/Boltzmann correction for charged
     sites). Both could live in the same `SolidPhase.n_mol` dict but need
     different activity treatment — extending the per-species activity-model
     override already used for gas/liquid species (`THERMODYNAMIC_MODEL_ARCHITECTURE.md`
     §3.3) to solid-phase species as well.
   - **Site-conservation as a component-total row**, potentially referencing
     another species' current amount (`total_sites = site_density × n_mol[mineral]`)
     rather than a fixed constant — the same "total assembled from another
     species' contribution" pattern already used by
     `NR_PRECIPITATION_CV_INTEGRATION.md`'s liquid+solid component totals,
     just one more instance of it.

   Explicitly out of scope for the current Layer 1 gap closure / declaration-API
   generalization work. Recorded here so it isn't rediscovered from scratch when
   `BisectionChemicalEquilibriumEngine`'s backend framework (tableau compilation, phase
   representation) is next revisited.

---

## 13. Variable-Scope Engines and the cv.advance() Contract

### 13.1 The engine scope problem

A key challenge arises with engines whose scope of what they resolve is
configurable — PHREEQC being the primary example. PHREEQC can operate in three
distinct modes:

| Mode | Engine resolves | Left for PhaseInterface |
|---|---|---|
| **Aqueous only** | Acid-base, complexation, redox | Henry's law, Ksp, adsorption |
| **Aqueous + gas phase** | All of above + Henry's law for all gas species | Kinetic kLa only |
| **Aqueous + gas + minerals** | All of above + Ksp equilibria for all minerals | Nothing (engine owns all equilibria) |

The `cv.advance()` method currently calls `engine.solve(phases=self.phases, T_K=...)`
then, later, `step_internal_transfer(dt_h)` which runs **all registered
`PhaseInterface` objects unconditionally**. This is safe only when the engine is in
aqueous-only mode. When PHREEQC is in Mode 2 or 3, the engine writes equilibrated
values to both `gas.n_mol` and `liquid.n_mol` in step 1. If `EquilibriumTransferModel`
for the same species is then run in step 4, two competing equilibrium models are
applied sequentially — one overwriting the other's result.

If the two models used identical thermodynamic parameters the re-partition would
compute zero correction. In practice they use different databases, different
activity corrections, and different temperature models — so the second pass
introduces a residual. This is not idempotent; it is **thermodynamically
inconsistent**.

For `KineticTransferModel` in Mode 2: the engine drives the species to
equilibrium; the subsequent kinetic transfer sees (C* − C) ≈ 0 and computes near-
zero flux. This is approximately idempotent only if both the engine and the
`PartitionModel` use the same Henry constant. Any discrepancy produces a small
spurious kinetic flux at each timestep, systematically biasing the result.

### 13.2 The cv.advance() implicit contract

Reading the sequential body of `cv.advance()` (lines ~540–629 of
`src/core/control_volume.py`), five assumptions are currently implicit rather than
documented:

```
Step 1  engine.solve(phases, T_K)

        IMPLICIT: engine scope = AQUEOUS ONLY
        IMPLICIT: post-condition is z = h(y) for single-phase chemistry only
        No assertion that engine did not modify gas.n_mol or solid.n_mol

Step 2  apply external_source_terms
Step 2b apply boundary fluxes

        IMPLICIT: z is now stale (engine was not re-called after source terms)
        The code comment reads: "operator-splitting contract: feeds in step 2
        do not shift the same-step pH" — this is the SNIA assumption,
        documented as a known approximation but not flagged as strategy-
        dependent behaviour

Step 3  compute_reaction_rates(t_h)

        IMPLICIT: uses z from step 1 (pre-feed state), not re-solved
        Code comment: "compute_rhs() would re-run speciation on the post-feed
        totals and break that invariant, so advance() uses compute_reaction_rates
        directly" — this is a deliberate inconsistency traded for performance

Step 4  step_internal_transfer(dt_h)

        IMPLICIT: engine did NOT handle any of these transfer phenomena
        IMPLICIT: all registered PhaseInterface objects are still needed
        No filtering by engine scope; no overlap check
```

Three things need to be explicit:

1. **Engine equilibrium scope** — what the engine resolved and what it did not
2. **Operator-splitting strategy** — the "pre-feed speciation for reactions" choice
   is the SNIA assumption embedded in the method body; SIA and Monolithic paths
   have different z-staleness semantics
3. **Non-overlap invariant** — nothing validates that PhaseInterface and engine
   scope don't cover the same phenomena

### 13.3 Resolution A: engine scope declaration

The engine declares at construction time what it resolved. The CV uses this to
filter `step_internal_transfer`:

```python
@dataclass(frozen=True)
class EquilibriumScope:
    """Phenomena an engine resolves internally beyond aqueous equilibria."""
    gas_liquid: frozenset = frozenset()   # species whose Henry's law is engine-owned
    mineral:    frozenset = frozenset()   # minerals whose Ksp is engine-owned

# Engine declares at construction
engine = PHREEQCChemicalEquilibriumEngine(
    component_map={"CO2": "C", "H2S": "S(-2)"},
    scope=EquilibriumScope(
        gas_liquid=frozenset({"CO2", "CH4", "H2", "H2S", "NH3"}),
        mineral=frozenset({"Calcite", "FeS"}),
    )
)

# step_internal_transfer becomes scope-aware
for interface in self.internal_interfaces:
    if engine_scope.owns(interface):
        continue     # engine already resolved this; skip to avoid double-handling
    interface.compute_flux(...)
```

**Advantage:** flexible — engine scope varies per instance without changing
registered `PhaseInterface` objects.
**Risk:** scope declaration and actual engine behaviour can diverge silently. A
PHREEQC instance configured without a gas phase but declared with
`gas_liquid={...}` would silently skip PhaseInterface objects that should run.

### 13.4 Resolution B: strict separation with construction-time validation

The engine scope is always `AQUEOUS_ONLY`. Users configuring PHREEQC to handle
gas-liquid simply do not register `EquilibriumTransferModel` for those species.
The CV validates at construction that no interface claims phenomena that would
overlap with the engine's scope:

```python
def _validate_no_scope_overlap(self):
    engine_scope = getattr(self.reaction_system.engine, "scope", EquilibriumScope())
    for interface in self.internal_interfaces:
        species_claimed = _species_claimed_by_interface(interface)
        if overlap := species_claimed & engine_scope.gas_liquid:
            raise ConfigurationError(
                f"PhaseInterface claims {overlap} already resolved by engine scope."
                f" Remove the interface or set engine scope to AQUEOUS_ONLY."
            )
```

**Advantage:** errors at construction, not at runtime; impossible to have silent
double-handling.
**Trade-off:** stricter configuration; harder to benchmark engine vs PhaseInterface
side-by-side.

A practical middle ground: implement Resolution A (scope declaration + filtering)
but add a **construction-time warning** whenever a registered `PhaseInterface`
would be skipped due to scope overlap, prompting the user to either remove the
interface or acknowledge the skip explicitly.

### 13.5 cv.advance() z-staleness and solver strategy

The "pre-feed speciation used for reactions" ordering is the SNIA assumption
embedded in `cv.advance()`. Its behaviour changes across the three solver
strategies:

| Strategy | z used in step 3 reactions | z staleness | Error order |
|---|---|---|---|
| **SNIA** | step-1 speciation (pre-feed) | always one substep stale | O(dt) |
| **SIA** | step-1 of current iteration | converges across iterations | O(dt²) at convergence |
| **Monolithic / `compute_rhs()`** | engine re-called at each RHS evaluation | current state | correct |

The current `cv.advance()` hardcodes the SNIA path. The Monolithic path goes
through `cv.compute_rhs()` instead, which re-runs the engine at each call (engine
writeback is part of `compute_rhs`). SIA — iterating `(transport → cv.advance())`
until convergence — is not yet implemented as a `SystemSolver` variant.

The key implication for engine scope: in the Monolithic path, PHREEQC in Mode 2
handles gas-liquid equilibrium **inside each `compute_rhs()` call**. The
gas-liquid partition is then an algebraic constraint evaluated by PHREEQC at each
RHS evaluation, not a rate term in the ODE. This means:

- A `KineticTransferModel` for species in `engine.scope.gas_liquid` is
  **also wrong in Monolithic mode** (not just `EquilibriumTransferModel`), because
  the engine enforces the equilibrium constraint instantaneously — there is no
  kinetic approach to equilibrium for those species.
- The scope-filter in `step_internal_transfer` is insufficient for Monolithic: the
  filter must also block `KineticTransferModel` when the engine scope covers those
  species.
- Construction-time validation (Resolution B) catches this earlier and more
  reliably than runtime filtering.

### 13.6 Proposed refined cv.advance() contract

```
cv.advance(dt_h, t_h, external_source_terms, solver, *, splitting_note="snia")

  PRECONDITIONS:
    - phase.n_mol contains conserved totals (y variables)
    - engine.scope (or default AQUEOUS_ONLY) declares what engine resolves
    - no PhaseInterface is registered for phenomena in engine.scope [validated at construction]

  STEP 1 — Algebraic solve
    engine.solve(phases, T_K)
    POST: z variables populated; phenomena in engine.scope are at equilibrium

  STEP 2 — External mass input
    apply external_source_terms
    apply boundary fluxes
    NOTE: z is now stale (y changed; engine not re-called)
          "SNIA contract": reactions in step 3 use z from step 1

  STEP 3 — Kinetic chemistry
    compute_reaction_rates(t_h)      [reads z from step 1, not re-run]
    NOTE: large feeds (step 2) shift y significantly → z from step 1 may be
          inaccurate for that step; corrected by SIA iteration or Monolithic path

  STEP 4 — Kinetic inter-phase transfer (scope-filtered)
    for each PhaseInterface NOT in engine.scope:
        interface.compute_flux(state_a, state_b, dt_h)
    [interfaces IN engine.scope are skipped — engine already resolved them]

  STEP 5 — Conservation audit
    conservation_monitor.check_step(phases)

  POSTCONDITIONS:
    - phase.n_mol updated with net changes from all active sub-steps
    - phenomena in engine.scope are at equilibrium (from step 1, now stale)
    - kinetic phenomena advanced by dt_h
    - z variables reflect step-1 speciation (SNIA); stale relative to step-2 y
```

---

## 14. Constraint-family taxonomy and folding limits (2026-07-01 design discussion)

> **§14.1 shipped 2026-07-02** as
> [`EQUILIBRIUM_CONSTRAINT_UNIFICATION.md`](../phases-shipped/EQUILIBRIUM_CONSTRAINT_UNIFICATION.md)
> (tag `equilibrium-constraint-unification-shipped`). Two corrections against
> what's written below, settled during implementation: (1) `EquilibriumConstraint`
> is attribute-based (`stoichiometry`, `log_K`, `dH_J_per_mol`, `T_ref_K` as
> plain attributes, not `log_K(T_K)`/`dH_J_per_mol()` methods) — a separate
> `vant_hoff_log_K()` free function does the temperature correction, since
> `EquilibriumReaction` already stores these as attributes and a method-based
> protocol would have forced it to change for no benefit; (2) `RaoultEquilibrium`
> keeps `RaoultPartition`'s actual shipped field names (`P_sat_ref`, `dH_vap`,
> `T_ref`, `C_water_mol_L`), not the placeholder `P_sat_ref_atm`/`dlnPsat` used
> below — `RaoultPartition` had already shipped (via
> `THERMODYNAMIC_MODEL_ARCHITECTURE.md`) by the time this section was written.
> §14.2/§14.3's folding-limit analysis is unaffected by either correction and
> remains accurate; the actual folding described there (non-ideality
> corrections wired into the Newton residual, Ksp/Henry rows merged into one
> Jacobian) is still unimplemented — this phase only unified the *declaration*
> surface (`classify_equilibrium_constraint()`, one object satisfying both
> `PartitionModel` and `EquilibriumConstraint`) and built
> `activity_for_entry()` as standalone groundwork; it is Phase 2's job
> (`LAYER1_GAP_CLOSURE.md`).

*Status: design discussion, not yet a checkpoint plan.* Captures the resolution
of the §6.1/§11 Priority 1 contradiction and a concrete limit on how far
folding can go without further work.

### 14.1 Two constraint families, not one `PartitionModel` wrapping decision

`EquilibriumReaction` (log-K, mass-action) and the `PartitionModel` family
(Henry, Ksp, Raoult) are not two competing representations of chemistry —
Henry's law and Ksp are *already* the same mass-action structure as an
acid-base reaction, just with domain-conventional parameters (`H_ref`/`dlnH`
map directly onto `log_K`/van't Hoff slope; the shipped
`NR_PRECIPITATION_SPECIATION` already sets `log_K = Ksp`). This is also the
standard architecture in geochemical equilibrium software (PHREEQC, MINTEQ,
EQ3/6): aqueous complexation, gas exchange, and mineral precipitation all
share one log-K/master-species formalism, differentiated only by which
non-ideality correction applies to which species.

Two families result, not one:

1. **Mass-action rows** — acid-base, Henry, Ksp, Raoult. Not one type with
   alternate constructors, but sibling types sharing a structural
   `EquilibriumConstraint` Protocol (`stoichiometry`, `log_K(T_K)`,
   `dH_J_per_mol()`), each with `__init__` keywords matching its own standard
   domain terminology rather than collapsing into a generic `log_K` at
   construction time: `EquilibriumReaction` (acid-base, needs no change — its
   native parameterization already *is* `log_K`), `HenryEquilibrium` (`H_ref`,
   `dlnH` — renamed from `HenryPartition`), `KspEquilibrium` (`Ksp`), and
   `RaoultEquilibrium` (`P_sat_ref_atm`, `dlnPsat`). Naming and the full
   constructor shape are decided in
   [`EQUILIBRIUM_CONSTRAINT_UNIFICATION.md`](../phases-upcoming/EQUILIBRIUM_CONSTRAINT_UNIFICATION.md)
   (Phase 1). The `*Equilibrium` types additionally satisfy `PartitionModel`
   unchanged (their existing role feeding `TransferModel`/kinetic use) — one
   constructed instance serves both the reaction-list/tableau role and the
   kinetic-transfer role, fixing the `HenryPartition`/reaction
   double-declaration bug that motivated this taxonomy in the first place.
   Non-ideality corrections (`γ_i` from `LiquidPhaseModel`, `φ_i` from
   `GasEOS`) are still read directly from the shared `ThermoFramework` when
   assembling the Newton residual, not via these types' own methods — see
   §14.2 for why.
2. **Isotherms** — Langmuir, Freundlich. Not mass-action-reformattable (no
   log-K exists for a saturable-capacity curve). Remain a distinct constraint
   type, safe to build on `PartitionModel`/`MultispeciesPartitionModel`.

### 14.2 Wrapping litmus test

A type may satisfy `EquilibriumConstraint` alongside its existing
`PartitionModel` role — or a constraint type may otherwise delegate to a
`PartitionModel`-family object — only if the relevant method is a pure
evaluate-at-a-point function, never one with an internal iteration or
convergence loop. `HenryEquilibrium`/`RaoultEquilibrium`/single-ion
`KspEquilibrium` all satisfy this (each has `partition_ratio() → float`
today — the codebase's own signal of linearity), which is why they can safely
dual-satisfy both protocols on one object per §14.1.
`MultispeciesVLEPartition.equilibrium_all_a_moles`
(§4.4 of `THERMODYNAMIC_MODEL_ARCHITECTURE.md`) explicitly "solves" the
coupled VLE internally; wrapping it inside an outer Newton row would nest a
solver inside a solver — the same anti-pattern as §13.1's PHREEQC dual-write
and the precipitation active-set loop below. The fix: mass-action rows for
non-ideal gas-liquid equilibria call `GasEOS.partial_pressures_atm()` and
`LiquidPhaseModel.gamma_all()` directly (both already pure evaluate
functions), letting the *outer* tableau's own Newton iteration do the fixed-point
search instead of nesting `MultispeciesVLEPartition`'s internal one.

Competitive Langmuir (`θ_i = K_i·C_i / (1 + Σ_j K_j·C_j)`) is directly
evaluable at any trial concentration vector — no internal solve — so wrapping
is safe there, provided future isotherm implementations keep this property.

### 14.3 Precipitation resolution and its real limit

Ksp stays folded into the engine, evaluated within the same `solve()` call and
timestep as acid-base and (once implemented) Henry rows — no sequential SNIA
correction. But "folded into the same tableau" is not, today, literally "one
flat Jacobian" for every mineral:

`NRTableau.build_tableau()` selects **one master species per connected
reaction-graph component** ([`src/speciation/nr_tableau.py`](../../src/speciation/nr_tableau.py))
and explicitly skips cross-phase reactions (`rxn.is_cross_phase`) — which
includes precipitation reactions. The shipped `NR_PRECIPITATION_SPECIATION`
therefore never feeds precipitation into this graph at all; it runs a
separate outer active-set loop around the inner acid-base NR instead. That
outer loop is why `THERMODYNAMIC_MODEL_ARCHITECTURE.md` §7 flags "active-set
discontinuities" in `∂z/∂y` at precipitation on/off boundaries — the
discontinuity is the Jacobian-level symptom of the nested-loop structure, not
a separate issue.

**Practical consequence:**
- **Gas-liquid (Henry) folding is unaffected by this limit.** A gas species
  either attaches onto an existing single-component ladder (CO₂, NH₃, H₂S —
  already have their own acid-base component) or forms a trivial standalone
  singleton (inert gases: O₂, CH₄, N₂, H₂). Neither merges two previously
  independent multi-species components, so realistic bioprocess/AD gas-liquid
  folding does not require any change to master-species selection.
- **Multi-ion Ksp folding as one flat Jacobian is gated on generalizing
  `NRTableau`'s master selection** — exactly the capability
  `MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md` Phase 1 is building,
  for unrelated reasons (Fe/Ca/phosphate/citrate networks). Single-driving-ion
  Ksp (the common bioprocess case — CaCO₃, struvite with one limiting ion) can
  continue using the existing nested active-set approach; it works and is
  tested, it's just not a single Jacobian, which matters only once the
  white-box BDF path (§10.4, deferred) is reached.
- **Convergence target, not immediate merge**: per 2026-07-01 discussion, the
  standalone `MultiComponentEquilibriumEngine`
  (`MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md`) keeps its existing
  build-separately-first rule unchanged. Its eventual integration decision
  (already gated on its own acceptance criteria) should target implementing
  the same constraint-row protocol described here, so it becomes a pluggable
  engine backend rather than a fourth independent precipitation pathway
  alongside the two already reconciled above.

## Cross-references

- [`THERMODYNAMIC_MODEL_ARCHITECTURE.md`](THERMODYNAMIC_MODEL_ARCHITECTURE.md) —
  thermodynamic layer beneath `PartitionModel`: `LiquidPhaseModel`, `GasEOS`,
  `ThermoFramework` update, `MultispeciesPartitionModel`, per-species EOS
  flexibility, and gas-liquid-solid coupling
- [`SOLVER_ARCHITECTURE.md`](SOLVER_ARCHITECTURE.md) —
  companion doc; Layer 1/2/3 solver coupling; SNIA/SIA/Monolithic strategies;
  identified extensions (SIA SystemSolver, reactive D_eff)
- [`CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md`](CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md) —
  engine/solver split; `ChemicalEquilibriumSystem`; black/gray/white-box
  protocols; `EquilibriumResult`; z-update strategies; `algebraic_species()`;
  `EquilibriumScope` overlap question (§12 Q8)
