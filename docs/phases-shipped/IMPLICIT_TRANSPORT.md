# IMPLICIT_TRANSPORT — Phase D

## Status

**Shipped 2026-06-12.** Tag: `implicit-transport-shipped`. 1339→1363 tests.
See [`SOLVER_ARCHITECTURE.md`](SOLVER_ARCHITECTURE.md) for context.

## Goal

Deliver `ImplicitTransportSystemSolver`: an IMEX (implicit-explicit) Axis 2
solver that treats inter-CV transport implicitly, making the system
unconditionally stable for transport stiffness. No new dependencies beyond
`scipy` (already required). Per-CV physics (Axis 1) stays explicit Euler
via `cv.advance()`.

## When this solver is warranted

**Trigger**: `MultirateSystemSolver` (Phase C) requires M > ~10 substeps
at the target `dt_h` to maintain CFL stability, and that M-fold overhead
dominates runtime.

The CFL number arises when inter-CV flow rates are fast relative to CV
volume: `τ_circ = V_cv / Q_circ`. Observed ranges across model families:

| Model | τ_circ | dt_h at CFL=1 |
|---|---|---|
| CSTR multi-tank (Q_recirc small) | 60–600 min | 60–600 min (no issue) |
| Compartmental fermenter, gentle 2-zone | 30–60 min | 30–60 min (manageable) |
| Compartmental fermenter, CFD-reduced | < 0.5 min | < 30 s (M ≈ 100+) |
| HPLC column, 20 cells | 0.05–0.5 min | 3–30 s (M ≈ 200+) |

For the bottom two families, explicit subcycling is impractical.
`ImplicitTransportSystemSolver` is appropriate from M ≈ 10 upward.

## IMEX scheme

Operator split: transport implicit, per-CV chemistry explicit.

At each macro step:

```
Step 1: compute per-CV explicit chemistry (cv.advance for each CV)
Step 2: assemble transport matrix A and source vector b from link fluxes
Step 3: solve (I + dt·A)·n* = n_after_chem + dt·b
Step 4: write n* back to cv phases
Step 5: fire controller events due in [t, t+dt]; apply ZOH-cached actions
```

The IMEX structure means the chemistry computation (Step 1) uses the
pre-step species concentrations; the transport system (Steps 2–4) uses
chemistry-updated concentrations as its initial state. This is the
Lie–Trotter IMEX split. Strang variant (transport half-step bracketing
chemistry) is possible as a `strang=True` option but is not the default.

## State vector and matrix assembly

The state vector follows the packing convention from Phase C:

```
n ∈ R^(N_cv × N_phase × N_species)
cv0_phaseA_sp0, cv0_phaseA_sp1, ..., cv0_phaseB_sp0, ..., cv1_phaseA_sp0, ...
```

Transport matrix **A** is assembled from all active links. Each link
contributes to **A** in the species subspace it acts on:

**AdvectiveLink** (Q L h⁻¹ from CV i to CV j, phase p, species s):

```
A[i,s, i,s] += Q / V_i
A[j,s, i,s] -= Q / V_i    (inflow to j from i)
```

(Rows: species entering/leaving each CV. Convention: positive diagonal = loss
from source.)

**DiffusiveLink** (exchange_rate E L h⁻¹ between CV i and CV j, phase p):

```
A[i,s, i,s] += E / V_i
A[i,s, j,s] -= E / V_j    (gain at i from j; coefficient of n_j is E/V_j)
A[j,s, j,s] += E / V_j
A[j,s, i,s] -= E / V_i    (gain at j from i; coefficient of n_i is E/V_i)
```

Note: an earlier draft of this document had the off-diagonal subscripts swapped
(A[i,j] -= E/V_i and A[j,i] -= E/V_j). The correct derivation from
dn_i/dt = −(E/V_i)·n_i + (E/V_j)·n_j gives the formulas above.

Both contribute symmetric off-diagonal entries, making **A** symmetric for
purely diffusive networks. Purely advective networks produce non-symmetric
**A** (upwind convention).

## Linear system solve

```
(I + dt·A)·n* = rhs
```

`scipy.sparse.linalg.spsolve` on a CSR-format sparse matrix. For networks
up to a few thousand state entries, this is one LU factorisation per solve.

**LU caching**: **A** changes only when:
- The link topology changes (links added/removed)
- Phase volumes change (CV volume tracking added in the future)
- `dt_h` changes between macro steps

For the common fixed-topology, fixed-volume, fixed-dt_h case, `(I + dt·A)`
is factored once and its LU decomposition cached for the entire run. Cache
invalidation checks at the start of each macro step compare current topology
hash and `dt_h` against the cached values. Recomputing an LU for a
50-node HPLC system takes ~1 ms; checking the hash takes ~1 µs.

```python
class ImplicitTransportSystemSolver:
    _cached_dt_h: Optional[float] = None
    _cached_lu: Optional[SuperLU] = None
    _cached_topo_hash: Optional[int] = None
```

## Controller events

Periodic controllers fire at T_c boundaries. Between T_c events, ZOH holds
the last `ControlAction` output. The IMEX solver fires controller events
only at the end of each macro step — it does not subcycle, so no sub-step
event handling is needed. If T_c < dt_h, the user should choose
`MonolithicODESolver` (Phase E) instead, which uses ODE-event stepping.

## Comparison to `MultirateSystemSolver`

| Property | MultirateSystemSolver | ImplicitTransportSystemSolver |
|---|---|---|
| Transport stability | CFL-bounded (needs M substeps) | Unconditionally stable |
| Per-step cost | O(M · link evaluations) | O(1 LU factorisation + backsolve) |
| Chemistry order | First-order Euler, no sub-steps | First-order Euler, no sub-steps |
| Controller events | Handled in sub-step loop | Macro-step boundaries only |
| T_c < dt_h support | Yes | No (use Phase E) |
| Sparse assembly | Not needed | Yes (scipy.sparse) |
| Extra dependencies | None | None (scipy already required) |

For M ≈ 10+, `ImplicitTransportSystemSolver` is typically faster.
For M < 5, `MultirateSystemSolver` is simpler and roughly equivalent in cost.

## Implementation checklist

- [x] `ImplicitTransportSystemSolver` in `src/core/system_solver.py`
- [x] `_assemble_transport_matrix(sim) → scipy.sparse.csr_matrix`
- [x] `_assemble_rhs` — RHS is `_pack_state(sim)` after chemistry; b=0 for current link types
- [x] `_pack_state(sim) → np.ndarray` (shared with Phase C utilities)
- [x] `_unpack_state(vec, sim) → None` (shared with Phase C utilities)
- [x] LU cache with topology-hash + dt_h invalidation
- [x] `DiffusiveLink` contributes symmetric blocks (Phase D corrects doc typo)
- [x] Controller events at macro-step boundaries via `sim._invoke_controllers`
      (EventScheduler not needed; Phase D does not subcycle)
- [x] ZOH semantics: chemistry step reads previous controller action; transport
      step reads post-chemistry state
- [x] Export from `vlsim.core`
- [x] Tests: 2-CV advective — single-step exact, mass conservation, non-negative
      at 100000× CFL
- [x] Tests: 2-CV diffusive — single-step exact, equilibrium, mass conservation
- [x] Tests: HPLC 10-cell chain — stable at dt_h=1h with Q=1000 (CFL=1000);
      convergence within 2% of ExplicitEuler at small dt_h
- [x] Tests: LU cache hit on second call; invalidated when `dt_h` changes
- [x] Tests: LU cache invalidated when link removed; invalidated when link added
- [x] Tests: controller fires at every macro-step boundary; same count as ExplicitEuler
- [x] Tests: single-CV no-link path runs without error (pure chemistry passthrough)

## Cross-references

- [`SOLVER_ARCHITECTURE.md`](SOLVER_ARCHITECTURE.md) — two-axis context;
  transport stiffness in the stiffness taxonomy
- [`SYSTEM_SOLVER_PROTOCOL.md`](SYSTEM_SOLVER_PROTOCOL.md) — Phase C prerequisite;
  `SystemSolver` protocol, `DiffusiveLink`, packing convention, `EventScheduler`
- [`MONOLITHIC_ODE.md`](MONOLITHIC_ODE.md) — Phase E; alternative for T_c < dt_h
- [`../../src/core/links.py`](../../src/core/links.py)
- [`../../src/core/simulation.py`](../../src/core/simulation.py)
