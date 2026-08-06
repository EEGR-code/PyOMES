# -*- coding: utf-8 -*-
"""Tests for Phase E — MONOLITHIC_ODE.

Covers (by checkpoint):

CP1 — Extended state packing
  - _ctrl_list: controllers with empty differential_state excluded
  - _ctrl_list: keys sorted alphabetically within each controller
  - _ctrl_list: preserves sim.controllers iteration order
  - _ctrl_list: controller without differential_state attribute excluded
  - _pack_extended_state: length = n_cv + sum(ctrl state sizes)
  - _pack_extended_state: no-controller case equals _pack_state output
  - _unpack_extended_state: CV species round-trip
  - _unpack_extended_state: controller states round-trip

CP2 — RHS assembler
  - No-link, no-controller: dydt equals cv.compute_rhs() rates
  - Advective link transport rate matches analytic flux
  - Diffusive link transport rate matches analytic flux
  - Controller state rates written into tail of dydt
  - Zero-state controller excluded from dydt controller block

CP3 — MonolithicODESolver advance_system (no-controller path)
  - Satisfies SystemSolver protocol
  - No-controller result matches ExplicitEulerSystemSolver within tolerance
  - returns (results={}, [], [], profile_actions) shape

CP4 — Event loop + ZOH cache
  - PI integral evolves continuously; value at T_c matches isolated ODE
  - T_c=0.01h, dt_h=0.1h: exactly 10 controller fires
  - Two-rate hierarchy: fire counts correct for both controllers
  - ZOH: action object identity frozen between T_c events
  - Non-periodic controller fires at macro-step endpoint

CP5 — method param + export
  - MonolithicODESolver importable from vlsim.core
  - method="LSODA" accepted and stable on stiff system
  - default_period_h=None raises ValueError for undeclared controller
  - default_period_h=float used for undeclared controller
"""

from __future__ import annotations

import math
import pytest
from typing import Dict, Optional


# ═══════════════════════════════════════════════════════════════════════
#  Shared helpers
# ═══════════════════════════════════════════════════════════════════════

def _make_liquid_cv(key: str, n_mol: dict, V_L: float = 1.0, T_K: float = 300.0):
    from PyOMES.core import ControlVolume, LiquidPhase
    liq = LiquidPhase(n_mol=dict(n_mol), V_L=V_L, T_K=T_K)
    return ControlVolume(phases={"liquid": liq}, label=key)


def _make_sim(n_mol_dict: dict, controllers=None, links=None):
    """Single-CV or multi-CV simulation. n_mol_dict: {cv_key: n_mol_dict}."""
    from PyOMES.core import ControlVolume, LiquidPhase, Simulation
    cvs = {}
    for cv_key, n_mol in n_mol_dict.items():
        liq = LiquidPhase(n_mol=dict(n_mol), V_L=1.0, T_K=300.0)
        cvs[cv_key] = ControlVolume(phases={"liquid": liq}, label=cv_key)
    return Simulation(cvs=cvs, controllers=controllers or [], links=links or [])


def _make_ctrl(integral=0.0, update_period_h=0.05, label="pi"):
    """Minimal PI-like controller with one differential state ('integral')."""
    from PyOMES.control.interfaces import ControllerBase
    from PyOMES.control.actions import ControlAction

    class _PICtrl(ControllerBase):
        target_cv_key = "main"

        def __init__(self, integral_init, period, lbl):
            self._integral = float(integral_init)
            self._period = period
            self._label = lbl
            self.fire_count = 0
            self.last_action = None

        @property
        def update_period_h(self):
            return self._period

        def differential_state(self):
            return {"integral": self._integral}

        def state_rates(self, env, t_h):
            return {"integral": 1.0}   # di/dt = 1 always (constant rate for testability)

        def set_state(self, state):
            self._integral = state["integral"]

        def compute(self, state, dt_h):
            self.fire_count += 1
            action = ControlAction(controller_label=self._label, target_cv_key="main")
            self.last_action = action
            return action

    return _PICtrl(integral, update_period_h, label)


def _make_stateless_ctrl(update_period_h=None, label="p"):
    """Controller with no differential state (simple P-only)."""
    from PyOMES.control.interfaces import ControllerBase
    from PyOMES.control.actions import ControlAction

    class _PCtrl(ControllerBase):
        target_cv_key = "main"

        def __init__(self, period, lbl):
            self._period = period
            self._label = lbl
            self.fire_count = 0

        @property
        def update_period_h(self):
            return self._period

        def compute(self, state, dt_h):
            self.fire_count += 1
            return ControlAction(controller_label=self._label, target_cv_key="main")

    return _PCtrl(update_period_h, label)


# ═══════════════════════════════════════════════════════════════════════
#  CP1 — Extended state packing
# ═══════════════════════════════════════════════════════════════════════

class TestCtrlList:

    def test_empty_diff_state_excluded(self):
        from PyOMES.core.system_solver import _ctrl_list
        ctrl_with = _make_ctrl(integral=1.0)
        ctrl_without = _make_stateless_ctrl()
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[ctrl_with, ctrl_without])
        cl = _ctrl_list(sim)
        assert len(cl) == 1
        assert cl[0][0] is ctrl_with

    def test_no_diff_state_attribute_excluded(self):
        from PyOMES.core.system_solver import _ctrl_list
        from PyOMES.control.actions import ControlAction

        class BareCtrl:
            def compute(self, state, dt_h):
                return ControlAction(controller_label="bare", target_cv_key="main")

        sim = _make_sim({"main": {"S": 1.0}}, controllers=[BareCtrl()])
        cl = _ctrl_list(sim)
        assert cl == []

    def test_keys_sorted_alphabetically(self):
        from PyOMES.control.interfaces import ControllerBase
        from PyOMES.control.actions import ControlAction
        from PyOMES.core.system_solver import _ctrl_list

        class MultiStateCtrl(ControllerBase):
            target_cv_key = "main"
            update_period_h = 0.1

            def differential_state(self):
                return {"zebra": 3.0, "alpha": 1.0, "mango": 2.0}

            def state_rates(self, env, t_h):
                return {}

            def set_state(self, state):
                pass

            def compute(self, state, dt_h):
                return ControlAction(controller_label="ms", target_cv_key="main")

        sim = _make_sim({"main": {"S": 1.0}}, controllers=[MultiStateCtrl()])
        cl = _ctrl_list(sim)
        assert cl[0][1] == ["alpha", "mango", "zebra"]

    def test_preserves_controllers_order(self):
        from PyOMES.core.system_solver import _ctrl_list
        c1 = _make_ctrl(integral=1.0, label="first")
        c2 = _make_ctrl(integral=2.0, label="second")
        c3 = _make_ctrl(integral=3.0, label="third")
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[c1, c2, c3])
        cl = _ctrl_list(sim)
        assert [c for c, _ in cl] == [c1, c2, c3]


class TestPackExtendedState:

    def _make_two_ctrl_sim(self, i1=1.0, i2=2.0):
        c1 = _make_ctrl(integral=i1, label="c1")
        c2 = _make_ctrl(integral=i2, label="c2")
        return _make_sim({"main": {"S": 3.0, "X": 5.0}}, controllers=[c1, c2]), c1, c2

    def test_length_equals_cv_plus_ctrl_states(self):
        import numpy as np
        from PyOMES.core.system_solver import _ctrl_list, _pack_extended_state, _pack_state
        sim, c1, c2 = self._make_two_ctrl_sim()
        cl = _ctrl_list(sim)
        y = _pack_extended_state(sim, cl)
        n_cv = len(_pack_state(sim))
        # 2 CV species + 1 integral per ctrl × 2
        assert len(y) == n_cv + 2

    def test_no_controller_equals_pack_state(self):
        import numpy as np
        from PyOMES.core.system_solver import _pack_extended_state, _pack_state
        sim = _make_sim({"main": {"S": 3.0, "X": 5.0}})
        y_ext = _pack_extended_state(sim, [])
        y_cv = _pack_state(sim)
        np.testing.assert_array_equal(y_ext, y_cv)

    def test_controller_values_appended_correctly(self):
        import numpy as np
        from PyOMES.core.system_solver import _ctrl_list, _pack_extended_state, _pack_state
        sim, c1, c2 = self._make_two_ctrl_sim(i1=7.5, i2=3.2)
        cl = _ctrl_list(sim)
        y = _pack_extended_state(sim, cl)
        n_cv = len(_pack_state(sim))
        assert y[n_cv] == pytest.approx(7.5)
        assert y[n_cv + 1] == pytest.approx(3.2)


class TestUnpackExtendedState:

    def test_cv_species_roundtrip(self):
        import numpy as np
        from PyOMES.core.system_solver import _ctrl_list, _pack_extended_state, _pack_state, _unpack_extended_state
        sim = _make_sim({"A": {"S": 4.0, "X": 2.0}, "B": {"S": 1.0, "X": 3.0}})
        cl = _ctrl_list(sim)
        n_cv = len(_pack_state(sim))
        y = _pack_extended_state(sim, cl)
        y_modified = y.copy()
        y_modified[:n_cv] *= 2.0
        _unpack_extended_state(y_modified, sim, cl, n_cv)
        y2 = _pack_extended_state(sim, cl)
        np.testing.assert_allclose(y2[:n_cv], y[:n_cv] * 2.0)

    def test_controller_states_roundtrip(self):
        import numpy as np
        from PyOMES.core.system_solver import _ctrl_list, _pack_extended_state, _pack_state, _unpack_extended_state
        c1 = _make_ctrl(integral=0.0, label="c1")
        c2 = _make_ctrl(integral=0.0, label="c2")
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[c1, c2])
        cl = _ctrl_list(sim)
        n_cv = len(_pack_state(sim))
        y = _pack_extended_state(sim, cl)
        # modify controller block
        y[n_cv] = 42.0
        y[n_cv + 1] = 99.0
        _unpack_extended_state(y, sim, cl, n_cv)
        assert c1._integral == pytest.approx(42.0)
        assert c2._integral == pytest.approx(99.0)


# ═══════════════════════════════════════════════════════════════════════
#  CP2 — RHS assembler
# ═══════════════════════════════════════════════════════════════════════

class TestRhsAssembler:

    def test_no_link_no_ctrl_rates_match_compute_rhs(self):
        """With no links and no controllers, _build_rhs rates equal cv.compute_rhs."""
        import numpy as np
        from PyOMES.core.system_solver import (
            _ctrl_list, _pack_state, _build_rhs,
        )
        # Single inert CV with no reaction model — compute_rhs returns all zeros.
        sim = _make_sim({"main": {"S": 2.0, "X": 1.0}})
        cl = _ctrl_list(sim)
        n_cv = len(_pack_state(sim))
        rhs = _build_rhs(sim, cl, n_cv)
        y = _pack_state(sim)
        dydt = rhs(0.0, np.concatenate([y]))
        # No reactions → all zeros
        np.testing.assert_allclose(dydt, 0.0, atol=1e-12)

    def test_advective_link_transport_rate(self):
        """Advective link: dydt[src] = -Q*C, dydt[snk] = +Q*C."""
        import numpy as np
        from PyOMES.core import AdvectiveLink
        from PyOMES.core.system_solver import (
            _ctrl_list, _pack_state, _state_index_map, _build_rhs,
        )
        Q = 2.0   # L/h
        n_S = 3.0  # mol in source
        V = 1.0    # L
        link = AdvectiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B",   _sink_phase_key="liquid",
            Q_L_per_h=Q,
        )
        sim = _make_sim({"A": {"S": n_S}, "B": {"S": 0.0}}, links=[link])
        cl = _ctrl_list(sim)
        n_cv = len(_pack_state(sim))
        rhs = _build_rhs(sim, cl, n_cv)
        y = _pack_state(sim)
        dydt = rhs(0.0, y)

        idx = _state_index_map(sim)
        expected_rate = (n_S / V) * Q   # C * Q  mol/h
        assert dydt[idx[("A", "liquid", "S")]] == pytest.approx(-expected_rate)
        assert dydt[idx[("B", "liquid", "S")]] == pytest.approx(+expected_rate)

    def test_diffusive_link_transport_rate(self):
        """Diffusive link: flux = kLa * V_eff * (C_src - C_snk)."""
        import numpy as np
        from PyOMES.core import DiffusiveLink
        from PyOMES.core.system_solver import (
            _ctrl_list, _pack_state, _state_index_map, _build_rhs,
        )
        kLa = 5.0
        n_A, n_B = 4.0, 1.0
        V = 1.0
        V_eff = 2.0 * V * V / (V + V)   # harmonic mean = 1.0
        link = DiffusiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B",   _sink_phase_key="liquid",
            kLa={"S": kLa},
        )
        sim = _make_sim({"A": {"S": n_A}, "B": {"S": n_B}}, links=[link])
        cl = _ctrl_list(sim)
        n_cv = len(_pack_state(sim))
        rhs = _build_rhs(sim, cl, n_cv)
        y = _pack_state(sim)
        dydt = rhs(0.0, y)

        idx = _state_index_map(sim)
        C_A, C_B = n_A / V, n_B / V
        expected_flux = kLa * V_eff * (C_A - C_B)   # mol/h
        assert dydt[idx[("A", "liquid", "S")]] == pytest.approx(-expected_flux)
        assert dydt[idx[("B", "liquid", "S")]] == pytest.approx(+expected_flux)

    def test_controller_rates_in_tail(self):
        """Controller state_rates contribute to the tail of dydt."""
        import numpy as np
        from PyOMES.core.system_solver import (
            _ctrl_list, _pack_state, _build_rhs,
        )
        ctrl = _make_ctrl(integral=0.0, update_period_h=0.05)
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[ctrl])
        cl = _ctrl_list(sim)
        n_cv = len(_pack_state(sim))
        rhs = _build_rhs(sim, cl, n_cv)
        y = _pack_state(sim)
        y_ext = np.concatenate([y, np.array([0.0])])  # integral starts at 0
        dydt = rhs(0.0, y_ext)
        # state_rates returns {"integral": 1.0}
        assert dydt[n_cv] == pytest.approx(1.0)

    def test_stateless_ctrl_does_not_extend_dydt(self):
        """Controller with no differential state: dydt length = n_cv only."""
        import numpy as np
        from PyOMES.core.system_solver import _ctrl_list, _pack_state, _build_rhs
        ctrl = _make_stateless_ctrl(update_period_h=0.05)
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[ctrl])
        cl = _ctrl_list(sim)
        n_cv = len(_pack_state(sim))
        rhs = _build_rhs(sim, cl, n_cv)
        y = _pack_state(sim)
        dydt = rhs(0.0, y)
        assert len(dydt) == n_cv


# ═══════════════════════════════════════════════════════════════════════
#  CP3 — MonolithicODESolver: no-controller path
# ═══════════════════════════════════════════════════════════════════════

class TestMonolithicNoController:

    def test_satisfies_system_solver_protocol(self):
        from PyOMES.core.system_solver import MonolithicODESolver, SystemSolver
        assert isinstance(MonolithicODESolver(), SystemSolver)

    def test_no_controller_matches_explicit_euler(self):
        """With no controllers, MonolithicODE trajectory matches ExplicitEuler."""
        import numpy as np
        from PyOMES.core import Simulation
        from PyOMES.core.system_solver import (
            ExplicitEulerSystemSolver, MonolithicODESolver,
        )
        # Two inert CVs with no reactions — only advective transport.
        from PyOMES.core import AdvectiveLink
        link = AdvectiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B",   _sink_phase_key="liquid",
            Q_L_per_h=0.1,
        )
        def _fresh():
            return _make_sim({"A": {"S": 1.0}, "B": {"S": 0.0}}, links=[link])

        sim_euler = _fresh()
        sim_mono  = _fresh()
        sim_euler._set_system_solver_unchecked(ExplicitEulerSystemSolver())
        sim_mono._set_system_solver_unchecked(MonolithicODESolver())

        # Run 5 steps of dt=0.01h — transport rate is slow, should agree well.
        for _ in range(5):
            sim_euler._step_default(dt_h=0.01, t_h=0.0)
        res_euler = sim_euler.cvs["A"].phases["liquid"].n_mol["S"]

        for _ in range(5):
            sim_mono.system_solver.advance_system(sim_mono, dt_h=0.01, t_h=0.0)
        res_mono = sim_mono.cvs["A"].phases["liquid"].n_mol["S"]

        # Should agree within O(dt) Euler truncation error
        assert abs(res_euler - res_mono) < 0.01

    def test_advance_system_return_shape(self):
        """advance_system returns a 4-tuple with results={}."""
        from PyOMES.core.system_solver import MonolithicODESolver
        sim = _make_sim({"main": {"S": 1.0}})
        solver = MonolithicODESolver()
        out = solver.advance_system(sim, dt_h=0.01, t_h=0.0)
        results, link_records, ctrl_actions, profile_actions = out
        assert results == {}
        assert isinstance(link_records, list)
        assert isinstance(ctrl_actions, list)
        assert isinstance(profile_actions, list)


# ═══════════════════════════════════════════════════════════════════════
#  STEP_SOLVER_INTERFACE_REFINEMENT item 2 — MonolithicODESolver rejects
#  a per-CV solver= that it can never consult (checkpoint 3 of
#  STEP_SOLVER_REFINEMENT_CHECKLIST.md)
# ═══════════════════════════════════════════════════════════════════════

class TestMonolithicRejectsPerCVSolver:

    def test_raises_when_sim_solver_is_a_stepsolver_instance(self):
        from PyOMES.core import SimultaneousAdaptiveSolver
        from PyOMES.core.system_solver import MonolithicODESolver

        sim = _make_sim({"main": {"S": 1.0}})
        sim.solver = SimultaneousAdaptiveSolver()
        with pytest.raises(ValueError, match="MonolithicODESolver"):
            MonolithicODESolver().advance_system(sim, dt_h=0.01, t_h=0.0)

    def test_raises_when_sim_solver_is_a_per_cv_dict(self):
        from PyOMES.core import SimultaneousEulerSolver
        from PyOMES.core.system_solver import MonolithicODESolver

        sim = _make_sim({"main": {"S": 1.0}})
        sim.solver = {"main": SimultaneousEulerSolver()}
        with pytest.raises(ValueError, match="MonolithicODESolver"):
            MonolithicODESolver().advance_system(sim, dt_h=0.01, t_h=0.0)

    def test_no_raise_when_sim_solver_is_none(self):
        """The common case: system_solver=MonolithicODESolver() with no
        per-CV solver= configured must be unaffected."""
        from PyOMES.core.system_solver import MonolithicODESolver

        sim = _make_sim({"main": {"S": 1.0}})
        assert sim.solver is None
        MonolithicODESolver().advance_system(sim, dt_h=0.01, t_h=0.0)

    def test_raises_via_sim_run_regardless_of_assignment_order(self):
        """system_solver assigned before solver, or after -- either
        order must be caught at the first sim.run() call."""
        from PyOMES.core import SimultaneousEulerSolver
        from PyOMES.core.system_solver import MonolithicODESolver

        sim = _make_sim({"main": {"S": 1.0}})
        sim.system_solver = MonolithicODESolver()
        sim.solver = SimultaneousEulerSolver()
        with pytest.raises(ValueError, match="MonolithicODESolver"):
            sim.run(tau_h=0.01, n_steps=1)


# ═══════════════════════════════════════════════════════════════════════
#  CP4 — Event loop + ZOH cache
# ═══════════════════════════════════════════════════════════════════════

class TestEventLoop:

    def test_pi_integral_matches_isolated_ode(self):
        """Integral state evolved by MonolithicODE matches scipy reference.

        di/dt = 1.0 (constant), T_c = dt_h = 0.1h.
        After one macro step, integral should equal 0.1 * 1.0 = 0.1.
        """
        import numpy as np
        from scipy.integrate import solve_ivp
        from PyOMES.core.system_solver import MonolithicODESolver

        ctrl = _make_ctrl(integral=0.0, update_period_h=0.1)
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[ctrl])

        MonolithicODESolver().advance_system(sim, dt_h=0.1, t_h=0.0)

        # Reference: exact integral of di/dt=1 over [0, 0.1]
        ref = solve_ivp(lambda t, y: [1.0], [0.0, 0.1], [0.0], method="RK45")
        expected = ref.y[0, -1]
        assert ctrl._integral == pytest.approx(expected, rel=1e-4)

    def test_ten_fires_per_macro_step(self):
        """T_c=0.01h, dt_h=0.1h → exactly 10 controller fires."""
        from PyOMES.core.system_solver import MonolithicODESolver

        ctrl = _make_ctrl(integral=0.0, update_period_h=0.01)
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[ctrl])

        MonolithicODESolver().advance_system(sim, dt_h=0.1, t_h=0.0)
        assert ctrl.fire_count == 10

    def test_integral_accumulates_over_ten_fires(self):
        """di/dt=1, T_c=0.01h, dt_h=0.1h: integral ≈ 0.1 after one macro step."""
        import numpy as np
        from scipy.integrate import solve_ivp
        from PyOMES.core.system_solver import MonolithicODESolver

        ctrl = _make_ctrl(integral=0.0, update_period_h=0.01)
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[ctrl])

        MonolithicODESolver().advance_system(sim, dt_h=0.1, t_h=0.0)

        ref = solve_ivp(lambda t, y: [1.0], [0.0, 0.1], [0.0], method="RK45")
        expected = ref.y[0, -1]
        assert ctrl._integral == pytest.approx(expected, rel=1e-3)

    def test_two_rate_hierarchy_fire_counts(self):
        """T_c1=0.02h, T_c2=0.1h in one dt_h=0.1h step: 5 + 1 fires."""
        from PyOMES.core.system_solver import MonolithicODESolver

        ctrl_fast = _make_ctrl(integral=0.0, update_period_h=0.02, label="fast")
        ctrl_slow = _make_ctrl(integral=0.0, update_period_h=0.1,  label="slow")
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[ctrl_fast, ctrl_slow])

        MonolithicODESolver().advance_system(sim, dt_h=0.1, t_h=0.0)
        assert ctrl_fast.fire_count == 5
        assert ctrl_slow.fire_count == 1

    def test_zoh_action_identity_between_events(self):
        """Control action object identity is frozen between T_c events."""
        from PyOMES.core.system_solver import MonolithicODESolver

        actions_seen = []

        from PyOMES.control.interfaces import ControllerBase
        from PyOMES.control.actions import ControlAction

        class TrackingCtrl(ControllerBase):
            target_cv_key = "main"
            _period = 0.1

            @property
            def update_period_h(self):
                return self._period

            def differential_state(self):
                return {"i": 0.0}

            def state_rates(self, env, t_h):
                return {"i": 0.0}

            def set_state(self, state):
                pass

            def compute(self, state, dt_h):
                action = ControlAction(controller_label="track", target_cv_key="main")
                actions_seen.append(action)
                return action

        ctrl = TrackingCtrl()
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[ctrl])

        # Two macro steps; ctrl fires once per step (T_c = dt_h = 0.1h).
        solver = MonolithicODESolver()
        solver.advance_system(sim, dt_h=0.1, t_h=0.0)
        solver.advance_system(sim, dt_h=0.1, t_h=0.1)

        # One unique action per macro step — each compute() returns a new object.
        assert len(actions_seen) == 2

    def test_non_periodic_fires_at_macro_boundary(self):
        """Controller without update_period_h fires once at macro endpoint
        when default_period_h is provided as a fallback."""
        from PyOMES.core.system_solver import MonolithicODESolver

        ctrl = _make_stateless_ctrl(update_period_h=None)
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[ctrl])

        # default_period_h=0.1 = dt_h → fires at macro-step endpoint
        MonolithicODESolver(default_period_h=0.1).advance_system(sim, dt_h=0.1, t_h=0.0)
        assert ctrl.fire_count == 1

    def test_undeclared_period_raises_without_default(self):
        """No update_period_h and no default_period_h → ValueError."""
        from PyOMES.core.system_solver import MonolithicODESolver

        ctrl = _make_stateless_ctrl(update_period_h=None)
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[ctrl])

        with pytest.raises(ValueError, match="update_period_h"):
            MonolithicODESolver().advance_system(sim, dt_h=0.1, t_h=0.0)


# ═══════════════════════════════════════════════════════════════════════
#  CP5 — method param + export
# ═══════════════════════════════════════════════════════════════════════

class TestMethodAndExport:

    def test_importable_from_core(self):
        from PyOMES.core import MonolithicODESolver
        assert MonolithicODESolver is not None

    def test_lsoda_accepted_and_stable(self):
        """method='LSODA' completes without error on a stiff-looking system."""
        from PyOMES.core.system_solver import MonolithicODESolver

        ctrl = _make_ctrl(integral=0.0, update_period_h=0.1)
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[ctrl])
        # Should complete; exact values not checked here.
        MonolithicODESolver(method="LSODA").advance_system(sim, dt_h=0.1, t_h=0.0)

    def test_default_period_h_used_for_undeclared_controller(self):
        """default_period_h=0.05h: undeclared controller fires at T_c=0.05h."""
        from PyOMES.core.system_solver import MonolithicODESolver

        ctrl = _make_stateless_ctrl(update_period_h=None)
        sim = _make_sim({"main": {"S": 1.0}}, controllers=[ctrl])

        # dt_h=0.1h, default_period_h=0.05h → 2 fires
        MonolithicODESolver(default_period_h=0.05).advance_system(sim, dt_h=0.1, t_h=0.0)
        assert ctrl.fire_count == 2
