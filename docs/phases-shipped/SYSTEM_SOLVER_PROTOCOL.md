# SYSTEM_SOLVER_PROTOCOL — Phase C

## Status

**Shipped 2026-06-12.** Tag: `system-solver-protocol-shipped`.
See [`SOLVER_ARCHITECTURE.md`](SOLVER_ARCHITECTURE.md) for context.

## Goal

Define the `SystemSolver` protocol (the Axis 2 extension point), wire
`Simulation._step()` dispatch, implement the event scheduling infrastructure
shared by all system solvers, and deliver three concrete solvers:
`ExplicitEulerSystemSolver` (validates wiring), `StrangSplittingSystemSolver`,
and `MultirateSystemSolver`. Add `DiffusiveLink` as a prerequisite for HPLC
and compartmental-fermenter models.

## `SystemSolver` protocol

```python
class SystemSolver(Protocol):
    def advance_system(
        self,
        sim: "Simulation",
        dt_h: float,
        t_h: float,
    ) -> Dict[str, AdvanceResult]: ...
```

`advance_system` is fully responsible for advancing all CVs, applying all
links, firing controller events, and updating controller states for one macro
timestep. It receives the full `Simulation` object and reads CVs and links
through the standard `sim.cvs` and link-iteration interfaces. No geometry-
specific assumptions.

The per-CV physics strategy (Axis 1) is retrieved from `Simulation.cv_solver`
(default for all CVs) or `Simulation.cv_solvers` (per-CV overrides) and
forwarded to each `cv.advance(solver=per_cv_solver)` call.

## `Simulation` dispatch wiring

```python
# Simulation.__init__ gains:
system_solver: Optional[SystemSolver] = None

# Simulation._step() becomes:
def _step(self, dt_h, t_h):
    if self._system_solver is not None:
        return self._system_solver.advance_system(self, dt_h, t_h)
    # existing explicit Euler body — unchanged when no system_solver set
    ...
```

The existing sequential body remains the default when `system_solver=None`.
`ExplicitEulerSystemSolver` wraps that body exactly; passing it explicitly
should produce identical results to the default path.

## `SystemEnv` type

Defined here; shared with Phase B's `state_rates` signature:

```python
@dataclass(frozen=True)
class SystemEnv:
    cvs: Dict[str, "ControlVolume"]   # read-only view of current CV states
    t_h: float

    def liquid_pH(self, cv_key: str) -> Optional[float]:
        liq = self.cvs[cv_key].phases.get("liquid")
        return liq.pH if liq is not None and "H+" in liq.n_mol else None
```

Passed to `controller.state_rates(env, t_h)` at each RHS evaluation.

## Event scheduling infrastructure

`EventScheduler` manages the priority queue of controller update events:

```python
class EventScheduler:
    def __init__(self, controllers, t_start_h: float): ...

    def next_event_h(self) -> Optional[float]:
        """Smallest scheduled fire time, or None if no periodic controllers."""

    def controllers_due(self, t_h: float) -> List[Controller]:
        """Controllers whose next fire time <= t_h + tolerance."""

    def advance_clocks(self, t_h: float) -> None:
        """Reschedule fired controllers to t_h + T_c."""
```

System solvers that take ODE sub-steps (Phases D, E) pass event times to the
integrator as stopping points. Simple solvers (below) use `EventScheduler`
only at macro step boundaries for correctness.

**ZOH semantics**: between T_c events, each controller's `ControlAction`
output is cached and re-applied unchanged. `compute()` is called only at
event times, not at every RHS evaluation (unless `update_period_h` is `None`).

## State vector packing convention

Shared across Phases C, D, and E. Defined once in `_pack_state` /
`_unpack_state` utilities:

```
for each cv (sim.cvs, insertion order):
    for each phase_key (alphabetical):
        for each species_key (alphabetical within phase.n_mol):
            append n_mol[species_key]
```

Phase-aware by construction: each link specifies `source_phase_key` and
`target_phase_key`. Any model expressed as a CV graph with consistent string
keys for CVs, phases, and species packs and unpacks automatically.

## Concrete solver implementations

### `ExplicitEulerSystemSolver`

Wraps the current `Simulation._step()` sequential body exactly. Ships first
to validate the `SystemSolver` dispatch wiring. Must produce identical results
to the non-`SystemSolver` path on all existing multi-CV tests. Zero behaviour
change; zero risk.

### `StrangSplittingSystemSolver`

Second-order Strang splitting of the link-flux operator from per-CV physics:

```
links(dt/2) → [cv.advance(dt) for each CV] → links(dt/2)
```

Cost: one extra link-flux evaluation per macro step — essentially free.
Same CFL stability bound as explicit Euler; halves the splitting error.
Periodic controllers fire at the half-step midpoint and at the endpoint.

### `MultirateSystemSolver`

Adaptive subcycling: M link evaluations at `dt_h / M` between one full
per-CV solve per macro step. M is chosen each macro step from the CFL number:

```python
M = ceil(dt_h / min_over_links(V_source_phase / Q_link))
```

CFL-free for transport. Cost: O(M · K · S) link evaluations vs O(N · S³)
for per-CV chemistry. Most valuable when M is large (HPLC-like models,
fast-circulation fermenters). For chemistry-dominated CVs, M-fold link
evaluations dominate; verify that the cost is justified per-run.

Controller events: if T_c < `dt_h / M`, periodic controllers still fire at
their correct T_c boundaries within the sub-stepping window.

## `DiffusiveLink`

```python
class DiffusiveLink:
    """Symmetric concentration-gradient exchange between two CVs.

    flux = exchange_rate × (C_j − C_i)   [symmetric; same magnitude both ways]
    exchange_rate : float   [L h⁻¹; equivalent to E_ij or D_ax · A / Δz]
    """
    def __init__(self, cv_i_key, cv_j_key, phase_key, exchange_rate_L_per_h): ...
```

Both HPLC axial dispersion and compartmental-fermenter turbulent backmixing
require symmetric reversible exchange that `AdvectiveLink` cannot represent.
`DiffusiveLink` contributes symmetric off-diagonal entries to the transport
matrix **A** (Phase D) and is handled correctly by all Phase C solvers via
its `compute_flux` method.

`source_phase_key` / `target_phase_key` attributes must be confirmed or added
to the `CVLink` base class as a prerequisite for Phase D matrix assembly.

## Implementation checklist

- [ ] `SystemSolver` protocol in `src/core/system_solver.py`
- [ ] `SystemEnv` dataclass
- [ ] `EventScheduler` class with priority queue
- [ ] `_pack_state(sim) → np.ndarray` utility
- [ ] `_unpack_state(vec, sim) → None` utility
- [ ] `Simulation.__init__` gains `system_solver` parameter
- [ ] `Simulation._step()` dispatch wired
- [ ] `ExplicitEulerSystemSolver` — identical behaviour to current code verified
- [ ] `StrangSplittingSystemSolver`
- [ ] `MultirateSystemSolver` with adaptive M and CFL estimate
- [ ] `DiffusiveLink` in `src/core/links.py`
- [ ] `source_phase_key` / `target_phase_key` confirmed or added to `CVLink` base
- [ ] Export `SystemSolver`, all three solvers, `DiffusiveLink`, `EventScheduler`
      from `vlsim.core`
- [ ] Tests: `ExplicitEulerSystemSolver` identical to current code on all
      existing multi-CV tests
- [ ] Tests: `StrangSplittingSystemSolver` halves splitting error vs explicit
      Euler on a 2-CV case with known analytical solution
- [ ] Tests: `MultirateSystemSolver` stable at `dt_h` that causes explicit
      Euler CFL instability on a deliberately fast-link case
- [ ] Tests: `DiffusiveLink` produces symmetric exchange; coexists with
      `AdvectiveLink` on the same CV boundary
- [ ] Tests: periodic controller events fire at correct T_c boundaries under
      both `StrangSplittingSystemSolver` and `MultirateSystemSolver`
- [ ] Tests: ZOH — control action frozen between events; integral still evolves

## Cross-references

- [`SOLVER_ARCHITECTURE.md`](SOLVER_ARCHITECTURE.md) — two-axis context
- [`CV_COMPUTE_INTERFACE.md`](CV_COMPUTE_INTERFACE.md) — Phase A prerequisite
- [`CONTROLLER_STATE_PROTOCOL.md`](CONTROLLER_STATE_PROTOCOL.md) — Phase B prerequisite; `SystemEnv`
- [`IMPLICIT_TRANSPORT.md`](IMPLICIT_TRANSPORT.md) — Phase D; uses `EventScheduler` + packing convention
- [`MONOLITHIC_ODE.md`](MONOLITHIC_ODE.md) — Phase E; uses `EventScheduler` + packing convention
- [`../../src/core/simulation.py`](../../src/core/simulation.py)
- [`../../src/core/links.py`](../../src/core/links.py)
