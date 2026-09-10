"""Controller protocol — the CV-native control loop interface (C7).

A :class:`Controller` reads a :class:`~PyOMES.core.snapshot.CVSnapshot`
(single-CV controllers) or
:class:`~PyOMES.core.snapshot.SimulationSnapshot` (multi-CV
controllers) and returns a structured
:class:`~PyOMES.control.actions.ControlAction` describing the
mutation it wants the orchestrator to apply.

Pattern 1 collapse (SIMULATION_CLASS.md decision 2, resolved
2026-05-26): the legacy two-step ``Controller.compute(state, cmd) ->
cmd`` plus ``Actuator.rates(state, cmd) -> dict`` is replaced by a
single ``compute(state, dt_h) -> ControlAction`` method. The
intermediate ``Commands`` dataclass and the separate ``Actuator``
protocol are deleted; controller-internal state (setpoint,
integral, derivative, sampling-hold buffer) lives on ``self``.

Concrete controllers in :mod:`PyOMES.control.loops` migrate to this
shape at C8. The orchestrator (:class:`~PyOMES.core.Simulation`)
invokes ``compute`` per step at C9. Until then the protocol is
declared but not consumed.

Phase B (CONTROLLER_STATE_PROTOCOL, 2026-06-11) adds :class:`ControllerBase`,
an optional mixin that provides default implementations of the four new
optional members required by the Phase C/E system integrator
(``update_period_h``, ``differential_state``, ``state_rates``,
``set_state``). Controllers that do not implement or inherit these members
are handled by the system integrator via ``getattr`` with safe defaults.
"""
from __future__ import annotations

from typing import Dict, Optional, Protocol, Union, runtime_checkable

from PyOMES.core.snapshot import CVSnapshot, SimulationSnapshot
from PyOMES.core.system_env import SystemEnv
from PyOMES.control.actions import ControlAction


@runtime_checkable
class Controller(Protocol):
    """CV-native controller protocol (single required method).

    Any object with a ``compute(state, dt_h) -> ControlAction`` method
    satisfies this protocol.

    **Optional duck-typed members** consumed by :class:`~PyOMES.core.Simulation`
    via ``getattr``:

    * ``sample_period_h: float`` or ``sample_period_s: float`` — macro-step
      gating. When set, ``compute()`` fires only at sampling instants; the
      last :class:`ControlAction` is reused between samples (zero-order hold).
    * ``target_cv_key: Optional[str]`` — disambiguates which CV snapshot is
      passed to ``compute()`` in multi-CV simulations.
    * ``reset() -> None`` — called by :meth:`Simulation.run` at run entry.
      Controllers without ``reset`` are tolerated.

    **Phase B optional members** consumed by the Phase C/E system integrator
    via ``getattr`` (default values shown). Inherit from :class:`ControllerBase`
    to receive these defaults without re-implementing them:

    * ``update_period_h: Optional[float]`` — ``None`` for continuous mode;
      positive float for discrete ZOH mode (fire at T_c boundaries only).
    * ``differential_state() -> Dict[str, float]`` — current values of
      ODE-governed controller states. Empty dict → excluded from the system
      state vector.
    * ``state_rates(env, t_h) -> Dict[str, float]`` — ODE RHS for the
      controller's differential state at (:class:`SystemEnv` ``env``, ``t_h``).
    * ``set_state(state: Dict[str, float]) -> None`` — solver write-back after
      each accepted ODE step.
    """

    def compute(
        self,
        state: Union[CVSnapshot, SimulationSnapshot],
        dt_h: float,
    ) -> ControlAction: ...


class ControllerBase:
    """Optional mixin providing Phase B default implementations.

    Inherit from this class to receive no-op / empty-return defaults for all
    four Phase B protocol members. Controllers that only implement
    ``compute()`` and do not need co-integration by the Phase C/E system
    integrator do **not** need to inherit from this class — the system
    integrator retrieves these members via ``getattr`` with safe defaults.

    Example::

        class PIpHController(ControllerBase):
            def __init__(self, setpoint, Kp, Ki, T_c_h):
                self._setpoint = setpoint
                self._Kp, self._Ki = Kp, Ki
                self._integral = 0.0
                self._T_c_h = T_c_h

            @property
            def update_period_h(self) -> float:
                return self._T_c_h

            def differential_state(self):
                return {"integral": self._integral}

            def state_rates(self, env, t_h):
                pH = env.liquid_pH("fermenter")
                return {"integral": self._setpoint - pH}

            def set_state(self, state):
                self._integral = state["integral"]

            def compute(self, state, dt_h):
                ...
    """

    @property
    def update_period_h(self) -> Optional[float]:
        """Discrete update period (hours), or ``None`` for continuous mode."""
        return None

    def differential_state(self) -> Dict[str, float]:
        """Current values of ODE-governed controller states. Empty → no state."""
        return {}

    def state_rates(
        self,
        env: SystemEnv,
        t_h: float,
    ) -> Dict[str, float]:
        """ODE RHS for controller differential state at (*env*, *t_h*)."""
        return {}

    def set_state(self, state: Dict[str, float]) -> None:
        """Write solver-accepted state values back onto the controller."""
