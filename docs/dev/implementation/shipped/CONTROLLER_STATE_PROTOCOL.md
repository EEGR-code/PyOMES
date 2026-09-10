# CONTROLLER_STATE_PROTOCOL — Phase B

## Status

**Shipped 2026-06-12.** Tag `controller-state-protocol-shipped` (after merge).
Independent of Phase A; required by Phase C/E. See
[`SOLVER_ARCHITECTURE.md`](SOLVER_ARCHITECTURE.md) for context.

## Goal

Extend the Pattern 1 controller protocol with optional differential-state and
update-period interfaces, enabling the Phase E system integrator to
co-integrate controller integral states alongside CV species concentrations
and to schedule discrete controller update events correctly. All additions are
optional; existing controllers require no changes.

## Protocol additions

The four new members default to no-ops / empty returns. A controller that
implements none of them behaves exactly as today.

### `update_period_h: Optional[float]`

```python
@property
def update_period_h(self) -> Optional[float]:
    return None   # default: fires at every system integrator RHS evaluation
```

`None` — continuous mode. `compute()` is called at every ODE sub-step RHS
evaluation. Appropriate for purely proportional controllers or smooth
lookup-based outputs.

`float` — discrete mode. `compute()` is called only at T_c boundaries
(multiples of `update_period_h`). Control output is held constant (zero-order
hold, ZOH) between events. The system integrator treats each T_c boundary as
an event: ODE integration halts, controller fires, output is refreshed,
integration resumes.

### `differential_state() → Dict[str, float]`

```python
def differential_state(self) -> Dict[str, float]:
    return {}   # default: no differential state
```

Returns current values of integral accumulators, observer states, or any other
controller state governed by a differential equation. Keys are arbitrary string
labels consistent with `state_rates`. The system integrator packs these values
into the system state vector alongside CV species concentrations when
`len(differential_state()) > 0`.

### `state_rates(env, t_h) → Dict[str, float]`

```python
def state_rates(self, env: "SystemEnv", t_h: float) -> Dict[str, float]:
    return {}   # default: no differential state
```

Returns `d/dt` of each entry in `differential_state()` evaluated at the current
system environment `env` and time `t_h`. This is the ODE right-hand side for
the controller's state. The system integrator calls this at every RHS
evaluation and packs the result alongside CV rate vectors.

### `set_state(state: Dict[str, float]) → None`

```python
def set_state(self, state: Dict[str, float]) -> None:
    pass   # default: no-op
```

Write-back from the solver after each accepted ODE step. The system integrator
calls this to update the controller's state to the solver's accepted values
before calling `compute()` to produce the next control action.

## Backward compatibility

Controllers that do not implement `differential_state`, `state_rates`, or
`set_state` return empty dicts / pass from the defaults. The system integrator
checks `len(c.differential_state()) > 0` before including a controller in the
extended state vector. The existing `compute(state) -> ControlAction` signature
is unchanged.

## Example: PI pH controller

```python
class PIpHController:
    def __init__(self, setpoint: float, Kp: float, Ki: float, T_c_h: float):
        self._setpoint = setpoint
        self._Kp, self._Ki = Kp, Ki
        self._integral = 0.0
        self._T_c_h = T_c_h

    @property
    def update_period_h(self) -> float:
        return self._T_c_h           # discrete: fire at T_c boundaries

    def differential_state(self) -> Dict[str, float]:
        return {"integral": self._integral}

    def state_rates(self, env, t_h: float) -> Dict[str, float]:
        pH = env.liquid_pH("fermenter")
        return {"integral": self._setpoint - pH}   # di/dt = e(t)

    def set_state(self, state: Dict[str, float]) -> None:
        self._integral = state["integral"]          # solver write-back

    def compute(self, state) -> ControlAction:
        pH = state.cv("fermenter")["liquid"].pH
        error = self._setpoint - pH
        u = self._Kp * error + self._Ki * self._integral
        return ControlAction(param_path="fermenter.base_feed_rate_mol_h",
                             value=max(0.0, u))
```

`_integral` is managed by the system integrator between T_c events via
`state_rates` + `set_state`, not accumulated manually. `compute()` reads
the current solver-managed value at the moment the controller fires.

## `SystemEnv` type

`state_rates` receives a `SystemEnv` object providing read-only access to
current CV state (all phases, current speciation outputs, property calculator
results) and current time `t_h`. `SystemEnv` is defined in Phase C alongside
the `SystemSolver` protocol.

## Hierarchical controllers

For a two-tier hierarchy (fast inner loop at T_c1, slow outer loop at T_c2),
both controllers declare their `update_period_h`. The Phase C `EventScheduler`
maintains a joint priority queue. At T_c2 events the outer controller fires
and potentially updates a setpoint read by the inner controller's `state_rates`.
No special-casing in the protocol — the composition emerges from independent
`update_period_h` declarations.

## Implementation checklist

- [x] Add `update_period_h` property to `ControllerBase` mixin (default `None`)
- [x] Add `differential_state()` to `ControllerBase` (default `{}`)
- [x] Add `state_rates(env, t_h)` to `ControllerBase` (default `{}`)
- [x] Add `set_state(state)` to `ControllerBase` (default no-op)
- [x] Define `SystemEnv` read-only type in `src/core/system_env.py`; includes `liquid_pH()` stub
- [x] Concrete built-in controllers require no changes — all duck-typed via `getattr`
- [x] Tests: PI controller — integral tracked by Euler loop matches manual at n=10/100/1000 steps
- [x] Tests: `update_period_h = None` (default) confirmed via `ControllerBase` and `getattr`
- [x] Tests: `update_period_h = T_c` readable from discrete controllers
- [x] Tests: backward compatibility — all six built-in controllers unchanged; `getattr` pattern verified

**Implementation note:** `Controller(Protocol)` retains only `compute()` as its required member
(preserving all existing `isinstance` checks). The four Phase B defaults live on
`ControllerBase`, an opt-in mixin. Phase C uses `getattr(ctrl, 'differential_state', None)`
etc. for duck-typed backward compat with structural controllers.

## Cross-references

- [`SOLVER_ARCHITECTURE.md`](SOLVER_ARCHITECTURE.md) — two-axis context; controller hybrid dynamics
- [`SYSTEM_SOLVER_PROTOCOL.md`](SYSTEM_SOLVER_PROTOCOL.md) — Phase C; `EventScheduler` + `SystemEnv` definition
- [`MONOLITHIC_ODE.md`](MONOLITHIC_ODE.md) — Phase E; primary consumer of this protocol
- [`../../src/core/simulation.py`](../../src/core/simulation.py) — current controller dispatch
