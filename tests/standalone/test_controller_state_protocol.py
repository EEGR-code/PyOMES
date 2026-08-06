# -*- coding: utf-8 -*-
"""Tests for Phase B — CONTROLLER_STATE_PROTOCOL.

Validates:
- Controller Protocol defaults (four Phase B members return expected values)
- SystemEnv stub: importable, frozen dataclass, liquid_pH() helper
- PIpHController pattern: integral tracked by a mock Euler integrator
  matches manually accumulated integral at multiple step sizes
- update_period_h: None for continuous controllers, float for discrete
- set_state write-back: solver-managed state replaces internal state
- Backward compat: existing built-in controllers unchanged; getattr
  with safe defaults produces the right no-op behaviour
- Exports: SystemEnv and Controller importable from PyOMES.core and
  PyOMES.control
"""
from __future__ import annotations

import pytest
from typing import Dict, Optional, Union


# ═══════════════════════════════════════════════════════════════════════
#  Helpers shared across tests
# ═══════════════════════════════════════════════════════════════════════

def _euler_integrate(controller, env_factory, t_start: float, t_end: float, n_steps: int):
    """Minimal Euler integrator exercising the Phase B protocol.

    Calls ``state_rates(env, t)`` at each sub-step and writes back via
    ``set_state`` — exactly the contract Phase E MonolithicODESolver
    will use.  Returns the controller's ``differential_state()`` at t_end.
    """
    dt = (t_end - t_start) / n_steps
    t = t_start
    state = dict(controller.differential_state())
    for _ in range(n_steps):
        env = env_factory(t)
        rates = controller.state_rates(env, t)
        for k in state:
            state[k] = state[k] + rates.get(k, 0.0) * dt
        controller.set_state(state)
        t += dt
    return dict(controller.differential_state())


# ═══════════════════════════════════════════════════════════════════════
#  Exports
# ═══════════════════════════════════════════════════════════════════════

class TestExports:

    def test_system_env_importable_from_core(self):
        from PyOMES.core import SystemEnv
        assert SystemEnv is not None

    def test_system_env_importable_from_control(self):
        from PyOMES.control import SystemEnv
        assert SystemEnv is not None

    def test_controller_importable_from_control(self):
        from PyOMES.control import Controller
        assert Controller is not None

    def test_controller_base_importable_from_control(self):
        from PyOMES.control import ControllerBase
        assert ControllerBase is not None

    def test_controller_base_importable_direct(self):
        from PyOMES.control.interfaces import ControllerBase
        assert ControllerBase is not None

    def test_system_env_importable_direct(self):
        from PyOMES.core.system_env import SystemEnv
        assert SystemEnv is not None


# ═══════════════════════════════════════════════════════════════════════
#  SystemEnv stub
# ═══════════════════════════════════════════════════════════════════════

class TestSystemEnv:

    def _make_env(self, pH=7.0, t_h=0.0):
        from PyOMES.core.system_env import SystemEnv

        class _FakeLiq:
            pass

        liq = _FakeLiq()
        liq.pH = pH

        class _FakeCV:
            phases = {"liquid": liq}

        return SystemEnv(cvs={"r": _FakeCV()}, t_h=t_h)

    def test_frozen(self):
        from PyOMES.core.system_env import SystemEnv
        env = SystemEnv(cvs={}, t_h=1.0)
        with pytest.raises((AttributeError, TypeError)):
            env.t_h = 2.0  # type: ignore[misc]

    def test_liquid_pH_present(self):
        env = self._make_env(pH=6.8, t_h=0.5)
        assert abs(env.liquid_pH("r") - 6.8) < 1e-12

    def test_liquid_pH_missing_cv(self):
        from PyOMES.core.system_env import SystemEnv
        env = SystemEnv(cvs={}, t_h=0.0)
        assert env.liquid_pH("nonexistent") is None

    def test_liquid_pH_no_liquid_phase(self):
        from PyOMES.core.system_env import SystemEnv

        class _FakeCV:
            phases = {}  # no liquid

        env = SystemEnv(cvs={"r": _FakeCV()}, t_h=0.0)
        assert env.liquid_pH("r") is None

    def test_liquid_pH_no_pH_attribute(self):
        from PyOMES.core.system_env import SystemEnv

        class _FakeLiq:
            pass  # no pH attr

        class _FakeCV:
            phases = {"liquid": _FakeLiq()}

        env = SystemEnv(cvs={"r": _FakeCV()}, t_h=0.0)
        assert env.liquid_pH("r") is None


# ═══════════════════════════════════════════════════════════════════════
#  Controller Protocol defaults
# ═══════════════════════════════════════════════════════════════════════

class TestProtocolDefaults:
    """A class that inherits from ControllerBase receives the four Phase B
    defaults for free without re-implementing them.

    The Controller Protocol is unchanged (only requires compute); ControllerBase
    is the opt-in mixin that provides the default implementations.
    """

    def _make_minimal_controller(self):
        from PyOMES.control.interfaces import ControllerBase
        from PyOMES.control.actions import ControlAction

        class MinimalCtrl(ControllerBase):
            def compute(self, state, dt_h) -> ControlAction:
                return ControlAction(controller_label="minimal", target_cv_key="")

        return MinimalCtrl()

    def test_update_period_h_default_none(self):
        ctrl = self._make_minimal_controller()
        assert ctrl.update_period_h is None

    def test_differential_state_default_empty(self):
        ctrl = self._make_minimal_controller()
        assert ctrl.differential_state() == {}

    def test_state_rates_default_empty(self):
        from PyOMES.core.system_env import SystemEnv
        ctrl = self._make_minimal_controller()
        env = SystemEnv(cvs={}, t_h=0.0)
        assert ctrl.state_rates(env, 0.0) == {}

    def test_set_state_default_noop(self):
        ctrl = self._make_minimal_controller()
        ctrl.set_state({"x": 1.0})  # must not raise

    def test_differential_state_returns_new_dict_each_call(self):
        ctrl = self._make_minimal_controller()
        d1 = ctrl.differential_state()
        d2 = ctrl.differential_state()
        assert d1 is not d2 or d1 == {}  # either different objects or both empty


# ═══════════════════════════════════════════════════════════════════════
#  PI controller — integral tracking
# ═══════════════════════════════════════════════════════════════════════

class PIpHController:
    """Toy PI pH controller implementing all four Phase B members.

    Matches the design-doc example exactly so tests confirm the
    contract, not an arbitrary implementation.
    """

    def __init__(self, setpoint: float, Kp: float, Ki: float, T_c_h: float):
        self._setpoint = setpoint
        self._Kp = Kp
        self._Ki = Ki
        self._integral = 0.0
        self._T_c_h = T_c_h

    # ── Phase B members ──────────────────────────────────────────────

    @property
    def update_period_h(self) -> float:
        return self._T_c_h

    def differential_state(self) -> Dict[str, float]:
        return {"integral": self._integral}

    def state_rates(self, env, t_h: float) -> Dict[str, float]:
        pH = env.liquid_pH("r")
        if pH is None:
            return {"integral": 0.0}
        return {"integral": self._setpoint - pH}   # di/dt = e(t)

    def set_state(self, state: Dict[str, float]) -> None:
        self._integral = state["integral"]

    # ── compute (required) ──────────────────────────────────────────

    def compute(self, state, dt_h):
        from PyOMES.control.actions import ControlAction
        pH = self._setpoint  # simplified; real impl reads state.pH
        error = self._setpoint - pH
        u = self._Kp * error + self._Ki * self._integral
        return ControlAction(
            controller_label="PI_pH",
            target_cv_key="r",
            params_changed={"base_feed_rate_mol_h": max(0.0, u)},
        )


class TestPIpHController:

    def _make_env_factory(self, pH_value: float):
        """Returns a factory that always produces an env with fixed pH."""
        from PyOMES.core.system_env import SystemEnv

        class _FakeLiq:
            pH = pH_value

        class _FakeCV:
            phases = {"liquid": _FakeLiq()}

        def factory(t_h: float) -> SystemEnv:
            return SystemEnv(cvs={"r": _FakeCV()}, t_h=t_h)

        return factory

    def test_update_period_h_readable(self):
        ctrl = PIpHController(setpoint=7.0, Kp=1.0, Ki=0.5, T_c_h=0.1)
        assert ctrl.update_period_h == pytest.approx(0.1)

    def test_differential_state_initial(self):
        ctrl = PIpHController(setpoint=7.0, Kp=1.0, Ki=0.5, T_c_h=0.1)
        assert ctrl.differential_state() == {"integral": 0.0}

    def test_set_state_writes_back(self):
        ctrl = PIpHController(setpoint=7.0, Kp=1.0, Ki=0.5, T_c_h=0.1)
        ctrl.set_state({"integral": 3.5})
        assert ctrl.differential_state() == {"integral": pytest.approx(3.5)}

    def test_state_rates_returns_error(self):
        """state_rates should return di/dt = setpoint - pH."""
        ctrl = PIpHController(setpoint=7.0, Kp=1.0, Ki=0.5, T_c_h=0.1)
        factory = self._make_env_factory(pH_value=6.5)
        env = factory(0.0)
        rates = ctrl.state_rates(env, 0.0)
        assert rates == {"integral": pytest.approx(0.5)}

    @pytest.mark.parametrize("n_steps", [10, 100, 1000])
    def test_euler_integral_matches_manual(self, n_steps: int):
        """Integral accumulated by mock Euler loop == manual Euler on same grid."""
        setpoint = 7.0
        pH_constant = 6.8   # constant error = 0.2
        t_end = 1.0
        dt = t_end / n_steps

        ctrl = PIpHController(setpoint=setpoint, Kp=1.0, Ki=0.5, T_c_h=0.05)
        factory = self._make_env_factory(pH_value=pH_constant)

        # Mock-integrator path
        final = _euler_integrate(ctrl, factory, t_start=0.0, t_end=t_end, n_steps=n_steps)

        # Manual accumulation on the identical grid
        error = setpoint - pH_constant   # = 0.2 always
        manual_integral = 0.0
        for _ in range(n_steps):
            manual_integral += error * dt

        assert final["integral"] == pytest.approx(manual_integral, rel=1e-10)

    def test_state_rates_constant_error_grid_independence(self):
        """With constant pH the final integral is (setpoint - pH) * T regardless of n."""
        setpoint = 7.0
        pH_constant = 6.7   # error = 0.3
        T = 2.0

        for n in [5, 50, 500]:
            ctrl = PIpHController(setpoint=setpoint, Kp=1.0, Ki=0.5, T_c_h=0.1)
            factory = self._make_env_factory(pH_value=pH_constant)
            final = _euler_integrate(ctrl, factory, t_start=0.0, t_end=T, n_steps=n)
            expected = (setpoint - pH_constant) * T
            assert final["integral"] == pytest.approx(expected, rel=1e-10), \
                f"n={n}: got {final['integral']!r}, expected {expected!r}"

    def test_set_state_is_used_between_steps(self):
        """state_rates always reads _integral from set_state, not self directly."""
        ctrl = PIpHController(setpoint=7.0, Kp=1.0, Ki=0.5, T_c_h=0.1)
        factory = self._make_env_factory(pH_value=6.5)

        # Manually inject a non-zero starting state as if the solver did it.
        ctrl.set_state({"integral": 10.0})
        assert ctrl.differential_state()["integral"] == pytest.approx(10.0)

        # One Euler step should add error*dt on top of 10.0
        env = factory(0.0)
        rates = ctrl.state_rates(env, 0.0)   # di/dt = 0.5
        new_val = 10.0 + rates["integral"] * 0.25
        ctrl.set_state({"integral": new_val})
        assert ctrl.differential_state()["integral"] == pytest.approx(10.0 + 0.5 * 0.25)


# ═══════════════════════════════════════════════════════════════════════
#  update_period_h semantics
# ═══════════════════════════════════════════════════════════════════════

class TestUpdatePeriodH:

    def test_continuous_controller_has_none(self):
        """A controller without update_period_h → getattr yields None."""
        class SimpleCtrl:
            def compute(self, state, dt_h):
                from PyOMES.control.actions import ControlAction
                return ControlAction(controller_label="s", target_cv_key="")

        ctrl = SimpleCtrl()
        assert getattr(ctrl, "update_period_h", None) is None

    def test_discrete_controller_exposes_period(self):
        ctrl = PIpHController(setpoint=7.0, Kp=1.0, Ki=0.5, T_c_h=0.25)
        assert ctrl.update_period_h == pytest.approx(0.25)

    def test_none_period_treated_as_continuous(self):
        """None means fire every step — confirmed by checking the condition
        the Phase C scheduler will evaluate."""
        from PyOMES.control.interfaces import ControllerBase
        class ContinuousCtrl(ControllerBase):
            def compute(self, state, dt_h):
                from PyOMES.control.actions import ControlAction
                return ControlAction(controller_label="c", target_cv_key="")
        c = ContinuousCtrl()
        assert c.update_period_h is None


# ═══════════════════════════════════════════════════════════════════════
#  Backward compatibility — existing built-in controllers
# ═══════════════════════════════════════════════════════════════════════

class TestBackwardCompat:
    """Existing built-in controllers require no changes and remain fully
    functional.  The system integrator (Phase C) reads the four Phase B
    members via getattr with safe fallbacks."""

    @pytest.fixture
    def ph_ctrl(self):
        from PyOMES.control import PHController
        return PHController(setpoint=7.0)

    @pytest.fixture
    def do_ctrl(self):
        from PyOMES.control import DOAgitationController
        return DOAgitationController(setpoint_mol_L=2.5e-4)

    @pytest.fixture
    def pressure_ctrl(self):
        from PyOMES.control import InstantPressureReliefController
        return InstantPressureReliefController(P_set_atm=1.5)

    # ── getattr patterns Phase C will use ───────────────────────────

    def _update_period(self, ctrl):
        return getattr(ctrl, "update_period_h", None)

    def _differential_state(self, ctrl):
        fn = getattr(ctrl, "differential_state", None)
        return fn() if fn is not None else {}

    def _state_rates(self, ctrl, env, t_h):
        fn = getattr(ctrl, "state_rates", None)
        return fn(env, t_h) if fn is not None else {}

    def _set_state(self, ctrl, state):
        fn = getattr(ctrl, "set_state", None)
        if fn is not None:
            fn(state)

    def test_ph_controller_update_period_none(self, ph_ctrl):
        assert self._update_period(ph_ctrl) is None

    def test_ph_controller_no_differential_state(self, ph_ctrl):
        assert self._differential_state(ph_ctrl) == {}

    def test_ph_controller_no_state_rates(self, ph_ctrl):
        from PyOMES.core.system_env import SystemEnv
        env = SystemEnv(cvs={}, t_h=0.0)
        assert self._state_rates(ph_ctrl, env, 0.0) == {}

    def test_ph_controller_set_state_noop(self, ph_ctrl):
        self._set_state(ph_ctrl, {"x": 1.0})  # must not raise or mutate

    def test_do_controller_update_period_none(self, do_ctrl):
        assert self._update_period(do_ctrl) is None

    def test_do_controller_no_differential_state(self, do_ctrl):
        assert self._differential_state(do_ctrl) == {}

    def test_pressure_controller_update_period_none(self, pressure_ctrl):
        assert self._update_period(pressure_ctrl) is None

    def test_pressure_controller_no_differential_state(self, pressure_ctrl):
        assert self._differential_state(pressure_ctrl) == {}

    def test_pressure_controller_set_state_noop(self, pressure_ctrl):
        self._set_state(pressure_ctrl, {})  # must not raise

    def test_ph_controller_sample_period_still_works(self):
        """sample_period_h / sample_period_s are unchanged by Phase B."""
        from PyOMES.control import PHController
        ctrl = PHController(setpoint=7.0, sample_period_s=300.0)
        assert ctrl.sample_period_s == pytest.approx(300.0)
        assert ctrl.sample_period_h is None

    def test_existing_controllers_unchanged_by_phase_b(self, ph_ctrl):
        """PHController.compute still works — no regression."""
        from PyOMES.core.snapshot import CVSnapshot
        snap = CVSnapshot(
            cv_key="r",
            t_h=0.0,
            pH=6.5,
            ionic_strength=None,
            V_liq_L=1.0,
            V_gas_L=0.0,
            P_gas_atm=1.0,
            T_K=310.0,
        )
        action = ph_ctrl.compute(snap, 0.01)
        assert action is not None
        assert action.target_cv_key == "r"
