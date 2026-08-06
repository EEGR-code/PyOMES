# TransferModel — Design Direction

> **Status:** Design decision record, 2026-06-18. Shipped 2026-06-19,
> tag `transfer-model-shipped`. `KineticTransferModel`, `EquilibriumTransferModel`,
> and `HenryPartition` are all live. Material now covered in
> `docs/design/MASS_EXCHANGE_ARCHITECTURE.md` §3 and §5.

---

## Motivation

The current gas-liquid transfer API has several entangled concerns that
make it harder to use and harder to generalise:

- **Verbose at the call site.** The common intra-CV case requires explicit
  construction of a `KineticGasLiquidLink` with four phase-routing keys,
  then passing it into `ControlVolume` as `internal_interfaces=[link]`.
  This exposes plumbing that the user shouldn't need to see.

- **Thermodynamics and transport are split across two dicts.**
  `partition_models` (chemistry) and `kLa` (equipment) are separate dicts
  on the link, keyed by the same species IDs. They describe two aspects of
  the same per-species phenomenon but are not composed together.

- **Transport mode is a separate bookkeeping set.**
  Whether a species is kinetically limited or at instantaneous equilibrium
  is recorded in `equilibrium_species` — a set maintained alongside the
  partition-model and kLa dicts. This is fragile and not self-documenting.

- **The analytical exponential solution is a hidden numerical fix.**
  `_kinetic_flux` uses an exact analytical step solution that is
  unconditionally stable regardless of `k_transfer × dt_h`. This exists
  solely to avoid Euler instability at high kLa. It couples the flux
  calculation to the step size, violating the principle that physics and
  numerics should be separate. The correct response to an oversized Euler
  step is a warning, not a silent correction.

- **`kLa` is gas-liquid specific by name and convention.**
  The underlying physics applies equally to solid-liquid adsorption,
  membrane transport, and any two-phase partition. The naming and
  architecture should be phase-agnostic.

---

## Design decisions

### 1. New `TransferModel` hierarchy

Introduce two concrete types that bundle the per-species physics for
inter-phase transfer:

```python
@dataclass
class KineticTransferModel:
    partition_model: PartitionModel   # thermodynamics — driving force
    k_transfer: float                 # transport rate coefficient (h⁻¹ or equivalent)

@dataclass
class EquilibriumTransferModel:
    partition_model: PartitionModel   # thermodynamics — equilibrium state
    # no k_transfer — transfer is instantaneous (algebraic constraint)
```

The transport mode is **implicit in the type**, not encoded as a flag or
inferred from the presence/absence of a key. The link/backend routes to the
correct calculation by `isinstance` check.

`k_transfer` is a mutable scalar on `KineticTransferModel`, directly
addressable by controllers (e.g. via ParamPath) without needing to reach
through the link. Changing agitation speed → update
`cv.transfer_models["CO2"].k_transfer`.

### 2. `PartitionModel` — single required method

The protocol is simplified to one required method:

```python
class PartitionModel(Protocol):
    def equilibrium_a_moles(
        self,
        n_total: float,
        capacity_a: float,
        capacity_b: float,
        T_K: float,
        alpha: float = 1.0,
    ) -> float:
        """Moles in phase 'a' at equilibrium given conserved total."""
        ...

    def partition_ratio(
        self,
        capacity_a: float,
        capacity_b: float,
        T_K: float,
        alpha: float = 1.0,
    ) -> Optional[float]:
        """Optional. Constant dimensionless ratio for linear models.

        When provided, enables performance-optimised calculation paths.
        Return None for concentration-dependent (non-linear) models.
        """
        ...
```

`equilibrium_a_moles` is the fundamental contract — it answers "where does
equilibrium lie?" for any model, linear or not. `partition_ratio` is an
optional performance hint; returning `None` causes a fallback to
`equilibrium_a_moles`. Neither method carries transport parameters or
knows anything about the solver.

`HenryPartition` is the canonical **reference implementation** for linear
gas-liquid partitioning. It is not privileged in the architecture; it
ships as a concrete starting point that users can copy and adapt.

### 3. `ControlVolume` gets a `transfer_models` kwarg

```python
cv = ControlVolume(
    phases={
        "gas":    gas_phase,    # GasPhase instance
        "liquid": liquid_phase, # LiquidPhase instance
    },
    transfer_models={
        "CO2": KineticTransferModel(HenryPartition(H_ref=3.4e-4, dlnH=2400.0), k_transfer=200.0),
        "CH4": KineticTransferModel(HenryPartition(H_ref=1.4e-5, dlnH=1600.0), k_transfer=150.0),
        "N2":  EquilibriumTransferModel(HenryPartition(H_ref=6.5e-6, dlnH=1300.0)),
    },
)
```

The CV inspects `phases` for a `GasPhase`/`LiquidPhase` pair (or
`LiquidPhase`/`SolidPhase` for adsorption, etc.) using the typed phase
classes that already exist, and creates the internal link automatically.
No import of `KineticGasLiquidLink` is required at the call site.

The same API generalises to solid-liquid without any structural change:

```python
cv = ControlVolume(
    phases={"liquid": liquid_phase, "solid": resin_phase},
    transfer_models={
        "Protein": KineticTransferModel(LangmuirPartition(q_max=0.12, K=8.5), k_transfer=0.3),
    },
)
```

**Edge case:** when a CV has more than two phases (e.g. gas + liquid +
solid), automatic phase-pair inference is ambiguous. In this case a
`phase_pair=("gas", "liquid")` kwarg specifies which pair
`transfer_models` applies to. This is a named exception, not the default.

### 4. How `TransferModel` types slot into the solver

See [DAE_ODE_STRUCTURE.md](DAE_ODE_STRUCTURE.md) for the full `compute_rhs()` contract and
the y/z variable split.

- **`KineticTransferModel`** contributes to `f_kinetic(y, z)`. It supplies
  the following term to the combined RHS:

  ```
  dn_liq/dt = k_transfer × (C_eq − C_liq)
  ```

  where `C_eq` is derived from `partition_model.equilibrium_a_moles()` at
  the current `(y, z)` state. This is a true instantaneous rate, returned
  inside `compute_rhs()` and applied once by the caller's forward step.

- **`EquilibriumTransferModel`** contributes to `solve_algebraic(y)`. It
  defines `z` — the per-phase distribution of a species — as an algebraic
  variable. It contributes **nothing to `dy/dt`**.

The two-type hierarchy makes this structural distinction explicit: the type
alone determines whether a species participates in the kinetic path or the
algebraic path, with no flag or separate bookkeeping set required.

### 5. `equilibrium_species` set is removed

The `equilibrium_species` set on `KineticGasLiquidLink` is superseded.
Transport mode is now encoded in the `TransferModel` type, and equilibrium
species contribute to `z` (the algebraic solve) rather than to `dy/dt`.
There is no separate set to maintain, and no risk of inconsistency between
a species appearing in `kLa` and also in `equilibrium_species`.

---

## Inter-CV transfer

`transfer_models` on `ControlVolume` is for **intra-CV** phase partitioning
only. Inter-CV transfer (linking two separate `ControlVolume` objects)
remains explicit at the `Simulation` level.

Inter-CV links reuse the same `TransferModel` types for per-species
physics but wrap them with explicit phase routing:

```python
simulation.add_link(
    CVLink(
        from_cv="reactor", from_phase="gas",
        to_cv="scrubber",  to_phase="liquid",
        transfer_models={
            "CO2": KineticTransferModel(HenryPartition(...), k_transfer=200.0),
        },
    )
)
```

`CVLink` is a new type. `KineticGasLiquidLink` may persist internally as
its implementation substrate, or be retired separately. The `TransferModel`
objects are identical whether used intra-CV or inter-CV.

---

## Boundaries

`Boundary` objects (feeds, vents, membranes) are **unchanged** and outside
the scope of this design. They represent material entering or leaving the
system from outside, driven by external conditions. They are not phase
partition models and should not be forced into the `TransferModel`
abstraction.

---

## What stays the same

| Component | Status |
|---|---|
| `PartitionModel` protocol | Simplified (single required method); otherwise unchanged |
| `HenryPartition` | Kept as reference implementation |
| `GasPhase`, `LiquidPhase`, `SolidPhase` | Unchanged — type inference relies on these |
| `Boundary` abstraction | Unchanged |
| `KineticGasLiquidLink` | Retained as internal implementation substrate; hidden from public API in the intra-CV case |
| `speciation_ladders` logic on the link | Unchanged internally |

---

## Deferred

- `CVLink` implementation details and how it integrates with `Simulation`
- Migration path and deprecation of the direct `KineticGasLiquidLink` public API
- Multiple phase pairs within one CV (gas + liquid + solid triple-phase case)
- Optional Jacobian methods on `PartitionModel` for non-linear implicit integration
- Stability warning thresholds and messaging for fixed-step Euler + high `k_transfer`
- Naming of the new types (`KineticTransferModel` / `EquilibriumTransferModel` vs shorter aliases for notebook use)
