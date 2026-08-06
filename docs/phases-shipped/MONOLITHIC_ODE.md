# MONOLITHIC_ODE — Phase E

## Status

**Shipped 2026-06-12.** Tag: `monolithic-ode-shipped`. 1363→1390 tests.
See [`../phases-shipped/SOLVER_ARCHITECTURE.md`](../phases-shipped/SOLVER_ARCHITECTURE.md) for context.

## Goal

Deliver `MonolithicODESolver`: an Axis 2 solver that co-integrates all CV
species concentrations, controller integral states, and inter-CV transport
in a single `scipy.solve_ivp` call. This correctly handles tight feedback
loops (T_c < dt_h), multi-rate digital controllers, and kinetically stiff
per-CV systems without the splitting error that accumulates when controllers
update only once per macro step.

## When this solver is warranted

**Primary trigger**: a controller with `update_period_h` < `dt_h` — the
common macro-step ZOH firing is too coarse, and integral state accuracy
matters (tight pH/DO control, observer-based state estimation).

**Secondary trigger**: per-CV kinetic stiffness that benefits from adaptive
sub-stepping at the system level rather than inside each `cv.advance()`.

For slow regulatory control (T_c >> dt_h), the existing default behaviour
is correct and requires no change. The `MonolithicODESolver` adds overhead
without benefit for those cases.

## System state vector

Extends the Phase C packing convention with controller differential state:

```
y ∈ R^(N_cv × N_phase × N_species + Σ_c |differential_state(c)|)

[   CV species block (Phase C packing order)   | controller 0 | controller 1 | ... ]
```

Controller indices are appended in the order they appear in
`sim.controllers`, with keys ordered alphabetically within each controller's
`differential_state()` dict. The same ordering is used for write-back via
`set_state`.

Only controllers with `len(differential_state()) > 0` contribute entries.
Controllers without differential state are not part of the state vector but
may still have their `compute()` called at event times.

## RHS assembler

```python
def _rhs(t_h: float, y: np.ndarray) -> np.ndarray:
    # 1. Unpack y → cv phases + controller states
    _unpack_extended_state(y, sim, controllers_with_state)

    # 2. Per-CV chemistry rates (Axis 1)
    dydt_chem = np.zeros_like(y)
    for cv in sim.cvs.values():
        rates = cv.compute_rhs(t_h)
        _write_rates(dydt_chem, rates, cv_offset_map)

    # 3. Inter-CV transport rates (link fluxes as explicit derivatives)
    dydt_transport = np.zeros_like(y)
    for link in sim.links:
        flux = link.compute_flux()
        _add_flux(dydt_transport, flux, cv_offset_map)

    # 4. Controller state rates
    env = SystemEnv(cvs=dict(sim.cvs), t_h=t_h)
    dydt_ctrl = np.zeros(n_ctrl_states)
    for c, offset in ctrl_offset_map.items():
        rates = c.state_rates(env, t_h)
        for key_i, key in enumerate(sorted(rates)):
            dydt_ctrl[offset + key_i] = rates[key]

    return np.concatenate([dydt_chem + dydt_transport, dydt_ctrl])
```

Transport appears as explicit derivatives here (not as a separate implicit
step). This is appropriate when the solver (`solve_ivp` with `RK45` or
`LSODA`) can take sub-steps shorter than the CFL timescale. For extremely
transport-stiff systems, `ImplicitTransportSystemSolver` (Phase D) is the
better choice.

## Controller event handling

`scipy.solve_ivp(events=...)` is used to locate T_c fire times exactly:

```python
event_fns = []
for c in periodic_controllers:
    def make_event(T_c, t_start):
        def event(t, y):
            return ((t - t_start) % T_c)   # zero at every T_c multiple
        event.terminal = True
        event.direction = -1
        return event
    event_fns.append(make_event(c.update_period_h, t_h))
```

The integrator is called in a loop:

```python
t_now = t_h
while t_now < t_h + dt_h:
    t_end = min(t_h + dt_h, scheduler.next_event_h())
    sol = solve_ivp(_rhs, [t_now, t_end], y_current,
                    method=method, events=event_fns, dense_output=False)
    y_current = sol.y[:, -1]
    t_now = sol.t[-1]

    # fire controllers that are due
    for c in scheduler.controllers_due(t_now):
        _unpack_extended_state(y_current, sim, controllers_with_state)
        env = SystemEnv(cvs=dict(sim.cvs), t_h=t_now)
        action = c.compute(...)
        _apply_action(sim, action)
        # refresh ZOH cache for this controller
        _zoh_cache[c] = action
    scheduler.advance_clocks(t_now)
```

**ZOH cache**: between T_c events, each controller's `ControlAction` is
held in `_zoh_cache`. The action is applied at the start of each sub-step.
This implements zero-order hold without multiple `compute()` calls.

**Multi-rate hierarchies**: the `EventScheduler` priority queue holds all
controllers' next fire times. Events interleave naturally. At a T_c2 (slow
outer loop) event, the outer controller fires and may update a setpoint that
the inner controller reads via `SystemEnv`. No special nesting logic; the
composition falls out of independent `update_period_h` declarations.

## ODE solver selection

```python
method: str = "RK45"          # default; accurate for smooth problems
                               # swap to "LSODA" for stiff per-CV kinetics
jac_sparsity: Optional[...] = None  # populated if CV Jacobians available
```

`LSODA` auto-selects stiff/non-stiff internally — useful when some CVs are
kinetically stiff and others are not. Provide `jac_sparsity` from the block-
diagonal structure of the system Jacobian to guide finite-difference
approximation.

If all CVs in the system implement `compute_jacobian(t_h)` (Phase A) and
return non-None, the full analytical block-diagonal Jacobian is assembled:

```
J_full = block_diag(*[cv.compute_jacobian(t_h) for cv in sim.cvs.values()])
# link flux Jacobians: analytical d(flux)/d(n_mol) if links expose Jacobian,
# else omit (FD fallback for off-diagonal transport blocks)
```

## Comparison to other Axis 2 solvers

| Property | ExplicitEuler | MultirateSystemSolver | ImplicitTransport | MonolithicODE |
|---|---|---|---|---|
| T_c < dt_h accuracy | No (fires once) | Partial | No (fires at end) | Yes |
| Transport stability | CFL-bounded | CFL-free (subcycles) | Unconditionally stable | Solver-adaptive |
| Co-integrated ctrl state | No | No | No | Yes |
| Multi-rate event handling | No | Partial | No | Yes |
| LSODA stiff per-CV | No | No | No | Yes (select `method`) |
| Extra dependencies | None | None | scipy.sparse | None |
| Per-step cost | Lowest | Medium | Low (cached LU) | Highest |

For tight control loops: `MonolithicODESolver`.
For transport-stiff HPLC/compartmental models with no tight controllers:
`ImplicitTransportSystemSolver` is cheaper.

## Implementation checklist

- [ ] `MonolithicODESolver` in `src/core/system_solver.py`
- [ ] `_pack_extended_state(sim, controllers) → np.ndarray`
- [ ] `_unpack_extended_state(y, sim, controllers) → None`
- [ ] `_rhs(t_h, y)` assembler (chem + transport + controller rates)
- [ ] `_make_event_fn(controller, t_start)` returning terminal event fn
- [ ] ZOH cache management (initialise, update at T_c, apply at each sub-step)
- [ ] Outer loop: run `solve_ivp` in segments between T_c events
- [ ] `EventScheduler` integration (imported from Phase C)
- [ ] `method` parameter (default `"RK45"`; user-selectable)
- [ ] `jac_sparsity` wiring when all `cv.compute_jacobian()` return non-None
- [ ] Export from `vlsim.core`
- [ ] Tests: PI controller — integral at t=T_c matches `solve_ivp` of
      `di/dt = e(t)` alone (isolated ODE reference)
- [ ] Tests: T_c=0.01h, dt_h=0.1h — 10 controller updates per macro step,
      integral matches continuous reference to within RK45 tolerance
- [ ] Tests: two-rate hierarchy (T_c1=0.02h, T_c2=0.1h) — both controllers
      fire at correct independent boundaries; no missed events
- [ ] Tests: outer controller setpoint update propagates to inner controller
      state_rates at correct event time
- [ ] Tests: MonolithicODESolver with no periodic controllers — identical
      results to ExplicitEulerSystemSolver (regression)
- [ ] Tests: LSODA method stable on a kinetically stiff single-CV case where
      RK45 requires many sub-steps
- [ ] Tests: ZOH — control action frozen between T_c events; integral
      continues to evolve

## Cross-references

- [`SOLVER_ARCHITECTURE.md`](SOLVER_ARCHITECTURE.md) — two-axis context;
  controller hybrid dynamics; multi-rate digital controllers
- [`SYSTEM_SOLVER_PROTOCOL.md`](SYSTEM_SOLVER_PROTOCOL.md) — Phase C prerequisite;
  `SystemSolver` protocol, `EventScheduler`, packing convention, `SystemEnv`
- [`CONTROLLER_STATE_PROTOCOL.md`](CONTROLLER_STATE_PROTOCOL.md) — Phase B;
  `differential_state`, `state_rates`, `set_state`, `update_period_h`
- [`CV_COMPUTE_INTERFACE.md`](CV_COMPUTE_INTERFACE.md) — Phase A;
  `compute_rhs`, `compute_jacobian`
- [`IMPLICIT_TRANSPORT.md`](IMPLICIT_TRANSPORT.md) — Phase D; alternative for
  transport-stiff systems without tight control loops
- [`../../src/core/simulation.py`](../../src/core/simulation.py)
- [`../../src/core/links.py`](../../src/core/links.py)
