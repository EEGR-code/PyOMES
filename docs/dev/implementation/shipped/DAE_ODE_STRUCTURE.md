# DAE/ODE Structure — Design Direction

> **Status:** Design decision record, 2026-06-18. Implemented in Phase A
> (CV_COMPUTE_INTERFACE, shipped 2026-06-11, tag `cv-compute-interface-shipped`).
> `compute_rhs()`, the y/z split, and the algebraic writeback are all live.
> Key concepts now documented in `docs/design/SOLVER_ARCHITECTURE.md`.

---

## Differential vs algebraic variables

The architecture enforces a clean DAE separation throughout `cv.advance()`:

```
dy/dt = f(y, z)     ← differential equations  (kinetic terms)
0     = g(y, z)     ← algebraic equations     (constraints)
```

- **`y` — differential variables**: conserved totals that kinetics change
  over time. Total moles of CO2 across both phases, total carbonate, total
  acetate, etc. Only the forward step (`y += dy/dt × dt`) modifies these.

- **`z` — algebraic variables**: the current equilibrium distribution of
  those totals, derived from `y` by solving the algebraic constraints.
  pH, speciation fractions, per-phase distributions of equilibrium species.
  Recomputed from `y` at each evaluation; never modified by the forward step.

`z` values appear freely in kinetic rate laws — that is their purpose. For
example, acetoclastic methanogenesis uses `[Ac⁻]` as its substrate
concentration. `[Ac⁻]` is an algebraic variable derived from total acetate
(`y`) and pH (`z`). The kinetic rate depends on `z`, but the returned
derivative is `d(y_total_acetate)/dt`. After the forward step reduces total
acetate, the next algebraic solve re-derives `[Ac⁻]`; the acid-base
equilibrium re-establishes implicitly.

### Timescale justification

This separation is valid when the equilibration timescale (microseconds for
acid-base, sub-second for Henry gas-liquid equilibria) is fast relative to
the kinetic timescale (minutes to hours for bioprocess reactions). When this
holds, the stiff fast-equilibrium dynamics disappear from the ODE and only
the slow kinetic dynamics remain.

---

## `compute_rhs()` contract

`compute_rhs()` encapsulates both the algebraic resolution and the kinetic
rate evaluation in a single call:

```
compute_rhs(t_h):
    z = solve_algebraic(y)      # speciation + equilibrium phase distributions
                                # z is derived from y; y is NOT mutated
    rates = f_kinetic(y, z)     # KineticReaction + kinetic transfer models
                                # all evaluated simultaneously at (y, z)
    return rates                # dy/dt only — no side-effects on y
```

Consequences:

- **No sequential splitting between reactions and kinetic transfer.** Both
  are evaluated at the same `(y, z)` snapshot and returned as a combined
  RHS. The caller applies a single forward step. The previous
  operator-splitting pattern — advancing reactions then advancing transfer
  separately — is eliminated.
- **ODE solvers call `compute_rhs()` at tentative states correctly.** Each
  call re-derives `z` from the current tentative `y`. This is not
  redundant — at a different `y`, `z` is genuinely different.

### `advance()` simplification

`advance()` reduces to two operations:

```python
rhs = cv.compute_rhs(t_h)
for sp, rate in rhs:
    phase.n_mol[sp] += rate * dt_h
```

### Writeback side-effect

The algebraic solve writes derived species (pH, HCO3⁻, per-phase
distributions of equilibrium species) back to `phase.n_mol` as a
side-effect so that rate laws can read them without extra plumbing. This
writeback is not a mutation of the differential state `y`. `compute_rhs()`
must not modify conserved totals — only the caller's forward step does that.

---

## Algebraic equilibrium path

Equilibrium models define `z` — the equilibrium distribution of a species
across phases — as an algebraic variable resolved during `solve_algebraic`
inside `compute_rhs()`. Given the conserved total `y_total` (moles summed
across both phases):

```
z_liq = partition_model.equilibrium_a_moles(y_total, V_liq, V_gas, T_K, alpha)
z_gas = y_total − z_liq
```

There is no `Δ` applied to the differential state. `y_total` is unchanged
by equilibrium transfer; it only changes when a kinetic reaction produces
or consumes the species. The per-phase values are written to `phase.n_mol`
during the algebraic solve so downstream rate laws can read them.

Equilibrium models contribute **nothing to `dy/dt`**. They participate only
in `solve_algebraic`, not in `f_kinetic`. This structural distinction is
what the two-type `TransferModel` hierarchy (see
[TRANSFER_MODEL.md](TRANSFER_MODEL.md)) makes explicit.

---

## Removal of the analytical exponential workaround

The previous `_kinetic_flux` used an exact analytical step solution that is
unconditionally stable regardless of `k_transfer × dt_h`. This existed
solely to avoid Euler instability at high kLa. It coupled the flux
calculation to the step size, violating the principle that physics and
numerics should be separate.

Under the new structure, kinetic transfer contributes a true instantaneous
rate returned inside `compute_rhs()` and applied once by the caller's
forward step. If a fixed-step Euler solver is used and
`k_transfer × dt_h` exceeds a stability threshold, the system warns rather
than silently compensates.
