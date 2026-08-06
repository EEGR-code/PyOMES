# -*- coding: utf-8 -*-
"""Tests for Phase D — IMPLICIT_TRANSPORT.

Covers:
- ImplicitTransportSystemSolver importable from vlsim.core
- SystemSolver protocol satisfied
- _assemble_transport_matrix: shape, AdvectiveLink entries, DiffusiveLink entries
- 2-CV advective: mass conserved; non-negative at CFL-violating dt_h
- 2-CV diffusive: mass conserved; equilibrium reached at equal concentrations
- Analytical implicit step for 1-way advection (single step exact check)
- LU cache hit on second call (same object identity)
- LU cache invalidated when dt_h changes
- LU cache invalidated when a link is removed
- HPLC 10-cell chain: stable at dt_h >> CFL dt (Q=1000, V=1, dt=1h)
- HPLC 10-cell chain: convergence at small dt_h — matches ExplicitEuler
- Controller fires at every macro-step boundary under ImplicitTransport
- ZOH: same controller count per step as ExplicitEuler
"""

from __future__ import annotations

import math
import pytest
from typing import Dict


# ═══════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════

def _make_liquid_cv(key: str, n_mol: dict, V_L: float = 1.0, T_K: float = 300.0):
    from PyOMES.core import ControlVolume, LiquidPhase
    liq = LiquidPhase(n_mol=dict(n_mol), V_L=V_L, T_K=T_K)
    return ControlVolume(phases={"liquid": liq}, label=key)


def _make_two_cv_sim(n_A: dict, n_B: dict, V_L: float = 1.0, links=None):
    from PyOMES.core import Simulation
    cv_A = _make_liquid_cv("A", n_A, V_L=V_L)
    cv_B = _make_liquid_cv("B", n_B, V_L=V_L)
    return Simulation(cvs={"A": cv_A, "B": cv_B}, links=links or [])


def _make_chain_sim(n_cells: int, Q: float, V_L: float = 1.0, init_S: float = 1.0):
    """n_cells linked in series by advective links at Q L/h; bolus in cell 0."""
    from PyOMES.core import ControlVolume, LiquidPhase, Simulation, AdvectiveLink
    cvs = {}
    links = []
    for k in range(n_cells):
        n_mol = {"S": init_S if k == 0 else 0.0}
        liq = LiquidPhase(n_mol=n_mol, V_L=V_L, T_K=300.0)
        cvs[f"cv{k}"] = ControlVolume(phases={"liquid": liq}, label=f"cv{k}")
    for k in range(n_cells - 1):
        links.append(AdvectiveLink(
            _source_cv_key=f"cv{k}",
            _source_phase_key="liquid",
            _sink_cv_key=f"cv{k+1}",
            _sink_phase_key="liquid",
            Q_L_per_h=Q,
        ))
    return Simulation(cvs=cvs, links=links)


# ═══════════════════════════════════════════════════════════════════════
#  Exports
# ═══════════════════════════════════════════════════════════════════════

class TestExports:

    def test_importable_from_core(self):
        from PyOMES.core import ImplicitTransportSystemSolver
        assert ImplicitTransportSystemSolver is not None

    def test_protocol_satisfied(self):
        from PyOMES.core import ImplicitTransportSystemSolver
        from PyOMES.core.system_solver import SystemSolver
        assert isinstance(ImplicitTransportSystemSolver(), SystemSolver)


# ═══════════════════════════════════════════════════════════════════════
#  Matrix assembly unit tests
# ═══════════════════════════════════════════════════════════════════════

class TestAssembleTransportMatrix:

    def _solver(self):
        from PyOMES.core import ImplicitTransportSystemSolver
        return ImplicitTransportSystemSolver()

    def test_empty_links_returns_zero_matrix(self):
        import numpy as np
        sim = _make_two_cv_sim({"S": 1.0}, {"S": 0.0})
        A = self._solver()._assemble_transport_matrix(sim)
        assert A.shape == (2, 2)
        assert A.nnz == 0

    def test_advective_link_diagonal_and_off_diagonal(self):
        import numpy as np
        from PyOMES.core import AdvectiveLink
        from PyOMES.core.system_solver import _state_index_map
        link = AdvectiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B", _sink_phase_key="liquid",
            Q_L_per_h=2.0,
        )
        sim = _make_two_cv_sim({"S": 1.0}, {"S": 0.0}, V_L=4.0, links=[link])
        A = self._solver()._assemble_transport_matrix(sim)
        idx = _state_index_map(sim)
        i = idx[("A", "liquid", "S")]
        j = idx[("B", "liquid", "S")]
        A_dense = A.toarray()
        # A[i, i] = Q/V = 2/4 = 0.5
        assert A_dense[i, i] == pytest.approx(0.5)
        # A[j, i] = -Q/V = -0.5  (sink gains)
        assert A_dense[j, i] == pytest.approx(-0.5)
        # Off-diagonals that should be zero
        assert A_dense[i, j] == pytest.approx(0.0)
        assert A_dense[j, j] == pytest.approx(0.0)

    def test_diffusive_link_symmetric_blocks(self):
        import numpy as np
        from PyOMES.core import DiffusiveLink
        from PyOMES.core.system_solver import _state_index_map
        link = DiffusiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B", _sink_phase_key="liquid",
            kLa={"S": 1.0}, V_eff_L=2.0,
        )
        sim = _make_two_cv_sim({"S": 1.0}, {"S": 0.0}, V_L=2.0, links=[link])
        A = self._solver()._assemble_transport_matrix(sim)
        idx = _state_index_map(sim)
        i = idx[("A", "liquid", "S")]
        j = idx[("B", "liquid", "S")]
        A_dense = A.toarray()
        # E = kLa * V_eff = 1.0 * 2.0 = 2.0, V_i = V_j = 2.0
        # A[i,i] = E/V_i = 2/2 = 1.0; A[j,j] = E/V_j = 1.0
        # A[i,j] = -E/V_j = -1.0; A[j,i] = -E/V_i = -1.0
        assert A_dense[i, i] == pytest.approx(1.0)
        assert A_dense[j, j] == pytest.approx(1.0)
        assert A_dense[i, j] == pytest.approx(-1.0)
        assert A_dense[j, i] == pytest.approx(-1.0)

    def test_diffusive_link_asymmetric_volumes(self):
        """Verify correct subscript assignment when V_i ≠ V_j."""
        import numpy as np
        from PyOMES.core import ControlVolume, LiquidPhase, Simulation, DiffusiveLink
        from PyOMES.core.system_solver import _state_index_map
        # V_A = 1, V_B = 3, kLa=1, V_eff_L=1 → E=1
        liq_A = LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=300.0)
        liq_B = LiquidPhase(n_mol={"S": 0.0}, V_L=3.0, T_K=300.0)
        cv_A = ControlVolume(phases={"liquid": liq_A}, label="A")
        cv_B = ControlVolume(phases={"liquid": liq_B}, label="B")
        link = DiffusiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B", _sink_phase_key="liquid",
            kLa={"S": 1.0}, V_eff_L=1.0,
        )
        sim = Simulation(cvs={"A": cv_A, "B": cv_B}, links=[link])
        A = self._solver()._assemble_transport_matrix(sim)
        idx = _state_index_map(sim)
        i = idx[("A", "liquid", "S")]
        j = idx[("B", "liquid", "S")]
        A_dense = A.toarray()
        # E=1, V_i=1, V_j=3
        # A[i,i]=1/1=1; A[j,j]=1/3; A[i,j]=-1/3; A[j,i]=-1/1=-1
        assert A_dense[i, i] == pytest.approx(1.0)
        assert A_dense[j, j] == pytest.approx(1.0 / 3.0)
        assert A_dense[i, j] == pytest.approx(-1.0 / 3.0)
        assert A_dense[j, i] == pytest.approx(-1.0)

    def test_matrix_shape_ten_cells(self):
        """10-cell chain: 10 state entries, (10×10) matrix."""
        sim = _make_chain_sim(10, Q=1.0)
        A = self._solver()._assemble_transport_matrix(sim)
        assert A.shape == (10, 10)


# ═══════════════════════════════════════════════════════════════════════
#  Analytical checks
# ═══════════════════════════════════════════════════════════════════════

class TestAnalytical:
    """Two-CV systems with exact or checkable solutions."""

    def _run(self, sim, dt_h, n_steps):
        from PyOMES.core import ImplicitTransportSystemSolver
        sim._system_solver = ImplicitTransportSystemSolver()
        result = sim.run(tau_h=dt_h * n_steps, n_steps=n_steps)
        return result

    # ── Advective ──────────────────────────────────────────────────────

    def test_advective_single_step_exact(self):
        """One implicit step of 1-way advective link matches closed-form."""
        from PyOMES.core import AdvectiveLink
        Q, V, dt_h = 2.0, 1.0, 0.5
        link = AdvectiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B", _sink_phase_key="liquid",
            Q_L_per_h=Q,
        )
        sim = _make_two_cv_sim({"S": 1.0}, {"S": 0.0}, V_L=V, links=[link])
        self._run(sim, dt_h=dt_h, n_steps=1)
        n_A = sim.cvs["A"].phases["liquid"].n_mol["S"]
        n_B = sim.cvs["B"].phases["liquid"].n_mol["S"]
        rate = Q / V
        # n*_A = 1 / (1 + dt*rate)
        n_A_expected = 1.0 / (1.0 + dt_h * rate)
        # n*_B = dt*rate * n*_A
        n_B_expected = dt_h * rate * n_A_expected
        assert n_A == pytest.approx(n_A_expected, rel=1e-10)
        assert n_B == pytest.approx(n_B_expected, rel=1e-10)

    def test_advective_mass_conserved(self):
        """Total moles conserved across any number of advective steps."""
        from PyOMES.core import AdvectiveLink
        link = AdvectiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B", _sink_phase_key="liquid",
            Q_L_per_h=1.0,
        )
        sim = _make_two_cv_sim({"S": 3.0}, {"S": 1.0}, V_L=1.0, links=[link])
        result = self._run(sim, dt_h=0.2, n_steps=10)
        A_final = float(result.phase_mol["A"]["liquid"]["S"][-1])
        B_final = float(result.phase_mol["B"]["liquid"]["S"][-1])
        assert A_final + B_final == pytest.approx(4.0, abs=1e-10)

    def test_advective_non_negative_at_huge_dt(self):
        """Implicit solve stays non-negative even at dt_h >> τ_CFL."""
        from PyOMES.core import AdvectiveLink
        # Q=1000, V=1 → τ_CFL = 0.001 h; dt_h = 100 h → 100000× CFL
        link = AdvectiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B", _sink_phase_key="liquid",
            Q_L_per_h=1000.0,
        )
        sim = _make_two_cv_sim({"S": 5.0}, {"S": 0.0}, V_L=1.0, links=[link])
        result = self._run(sim, dt_h=100.0, n_steps=3)
        for cv_key in ("A", "B"):
            vals = result.phase_mol[cv_key]["liquid"]["S"]
            assert all(v >= -1e-12 for v in vals), f"{cv_key} went negative"

    # ── Diffusive ──────────────────────────────────────────────────────

    def test_diffusive_equilibrium_equal_volumes(self):
        """Symmetric diffusion converges to equal concentrations."""
        from PyOMES.core import DiffusiveLink
        link = DiffusiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B", _sink_phase_key="liquid",
            kLa={"S": 5.0}, V_eff_L=1.0,
        )
        sim = _make_two_cv_sim({"S": 4.0}, {"S": 0.0}, V_L=1.0, links=[link])
        result = self._run(sim, dt_h=0.1, n_steps=200)
        A_final = float(result.phase_mol["A"]["liquid"]["S"][-1])
        B_final = float(result.phase_mol["B"]["liquid"]["S"][-1])
        assert A_final + B_final == pytest.approx(4.0, abs=1e-8)
        assert A_final == pytest.approx(2.0, abs=0.01)
        assert B_final == pytest.approx(2.0, abs=0.01)

    def test_diffusive_mass_conserved(self):
        """Mass is conserved under diffusive link for arbitrary dt_h."""
        from PyOMES.core import DiffusiveLink
        link = DiffusiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B", _sink_phase_key="liquid",
            kLa={"S": 2.0}, V_eff_L=1.0,
        )
        sim = _make_two_cv_sim({"S": 6.0}, {"S": 2.0}, V_L=1.0, links=[link])
        result = self._run(sim, dt_h=0.5, n_steps=20)
        A_final = float(result.phase_mol["A"]["liquid"]["S"][-1])
        B_final = float(result.phase_mol["B"]["liquid"]["S"][-1])
        assert A_final + B_final == pytest.approx(8.0, abs=1e-8)

    def test_diffusive_single_step_exact(self):
        """Exact check for one implicit diffusive step with equal volumes."""
        from PyOMES.core import DiffusiveLink
        kla, V_eff, V, dt_h = 1.0, 1.0, 1.0, 0.5
        E = kla * V_eff  # = 1.0
        link = DiffusiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B", _sink_phase_key="liquid",
            kLa={"S": kla}, V_eff_L=V_eff,
        )
        sim = _make_two_cv_sim({"S": 2.0}, {"S": 0.0}, V_L=V, links=[link])
        self._run(sim, dt_h=dt_h, n_steps=1)
        n_A = sim.cvs["A"].phases["liquid"].n_mol["S"]
        n_B = sim.cvs["B"].phases["liquid"].n_mol["S"]
        # E=1, V_i=V_j=1 → A = [[1,-1],[-1,1]]
        # M = I + dt*A = [[1+dt,-dt],[-dt,1+dt]] with dt=0.5
        # M = [[1.5, -0.5], [-0.5, 1.5]]
        # rhs = [2, 0]; solve Mx=rhs:
        # det = 1.5^2 - 0.5^2 = 2.25 - 0.25 = 2.0
        # x_A = (1.5*2 - (-0.5)*0) / 2.0 = 3.0/2.0 = 1.5
        # x_B = (1.5*0 - (-0.5)*2) / 2.0 = 1.0/2.0 = 0.5
        assert n_A == pytest.approx(1.5, rel=1e-10)
        assert n_B == pytest.approx(0.5, rel=1e-10)


# ═══════════════════════════════════════════════════════════════════════
#  LU cache
# ═══════════════════════════════════════════════════════════════════════

class TestLUCache:

    def _make_solver_and_sim(self, Q=1.0, dt_h=0.1):
        from PyOMES.core import ImplicitTransportSystemSolver, AdvectiveLink
        link = AdvectiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B", _sink_phase_key="liquid",
            Q_L_per_h=Q,
        )
        sim = _make_two_cv_sim({"S": 1.0}, {"S": 0.0}, links=[link])
        solver = ImplicitTransportSystemSolver()
        return solver, sim

    def test_cache_hit_on_second_call(self):
        """Second _get_lu call with same dt_h returns the same SuperLU object."""
        solver, sim = self._make_solver_and_sim()
        lu1 = solver._get_lu(sim, dt_h=0.1)
        lu2 = solver._get_lu(sim, dt_h=0.1)
        assert lu1 is lu2

    def test_cache_miss_on_dt_change(self):
        """Changing dt_h invalidates the cache."""
        solver, sim = self._make_solver_and_sim()
        lu1 = solver._get_lu(sim, dt_h=0.1)
        lu2 = solver._get_lu(sim, dt_h=0.2)
        assert lu1 is not lu2

    def test_cache_miss_on_link_removal(self):
        """Removing a link changes the topology hash → cache invalidated."""
        from PyOMES.core import ImplicitTransportSystemSolver, AdvectiveLink
        link = AdvectiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B", _sink_phase_key="liquid",
            Q_L_per_h=1.0,
        )
        sim = _make_two_cv_sim({"S": 1.0}, {"S": 0.0}, links=[link])
        solver = ImplicitTransportSystemSolver()
        lu1 = solver._get_lu(sim, dt_h=0.1)
        sim.links.clear()
        lu2 = solver._get_lu(sim, dt_h=0.1)
        assert lu1 is not lu2

    def test_cache_miss_on_link_added(self):
        """Adding a new link also invalidates the cache."""
        from PyOMES.core import ImplicitTransportSystemSolver, AdvectiveLink
        sim = _make_two_cv_sim({"S": 1.0}, {"S": 0.0})
        solver = ImplicitTransportSystemSolver()
        lu1 = solver._get_lu(sim, dt_h=0.1)
        link = AdvectiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B", _sink_phase_key="liquid",
            Q_L_per_h=1.0,
        )
        sim.links.append(link)
        lu2 = solver._get_lu(sim, dt_h=0.1)
        assert lu1 is not lu2

    def test_cache_miss_on_rate_change(self):
        """Changing Q on a link changes the hash."""
        from PyOMES.core import ImplicitTransportSystemSolver, AdvectiveLink
        link = AdvectiveLink(
            _source_cv_key="A", _source_phase_key="liquid",
            _sink_cv_key="B", _sink_phase_key="liquid",
            Q_L_per_h=1.0,
        )
        sim = _make_two_cv_sim({"S": 1.0}, {"S": 0.0}, links=[link])
        solver = ImplicitTransportSystemSolver()
        lu1 = solver._get_lu(sim, dt_h=0.1)
        link.Q_L_per_h = 5.0
        lu2 = solver._get_lu(sim, dt_h=0.1)
        assert lu1 is not lu2


# ═══════════════════════════════════════════════════════════════════════
#  HPLC 10-cell chain
# ═══════════════════════════════════════════════════════════════════════

class TestHPLCChain:
    """10 cells in series with fast advective links.

    τ_CFL = V/Q = 1/1000 = 0.001 h.  Explicit Euler at dt_h=1.0 would be
    wildly unstable (CFL number = 1000).  ImplicitTransport must be stable.
    """

    def _run_chain(self, system_solver, Q, dt_h, n_steps, n_cells=10, V_L=1.0):
        sim = _make_chain_sim(n_cells, Q=Q, V_L=V_L)
        sim._system_solver = system_solver
        result = sim.run(tau_h=dt_h * n_steps, n_steps=n_steps)
        return result

    def test_stable_at_cfl_violating_dt(self):
        """At Q=1000, V=1, dt=1h (CFL=1000), implicit stays non-negative."""
        from PyOMES.core import ImplicitTransportSystemSolver
        result = self._run_chain(
            ImplicitTransportSystemSolver(), Q=1000.0, dt_h=1.0, n_steps=5
        )
        total_mass = sum(
            float(result.phase_mol[f"cv{k}"]["liquid"]["S"][-1])
            for k in range(10)
        )
        for k in range(10):
            vals = result.phase_mol[f"cv{k}"]["liquid"]["S"]
            assert all(v >= -1e-10 for v in vals), f"cv{k} went negative"
        assert total_mass == pytest.approx(1.0, abs=1e-6)

    def test_mass_conserved_all_steps(self):
        """At every recorded step, total chain mass equals initial mass."""
        from PyOMES.core import ImplicitTransportSystemSolver
        result = self._run_chain(
            ImplicitTransportSystemSolver(), Q=1000.0, dt_h=1.0, n_steps=10
        )
        import numpy as np
        for step_idx in range(11):  # n_steps+1 recorded points
            total = sum(
                float(result.phase_mol[f"cv{k}"]["liquid"]["S"][step_idx])
                for k in range(10)
            )
            assert total == pytest.approx(1.0, abs=1e-6), f"step {step_idx}: mass={total}"

    def test_convergence_at_small_dt(self):
        """At small dt_h where explicit is stable, both solvers agree closely."""
        from PyOMES.core import ImplicitTransportSystemSolver, ExplicitEulerSystemSolver
        # Q=0.1, V=1 → τ_CFL = 10h; dt=0.01h is well within CFL
        Q, dt_h, n_steps = 0.1, 0.01, 50
        result_implicit = self._run_chain(
            ImplicitTransportSystemSolver(), Q=Q, dt_h=dt_h, n_steps=n_steps
        )
        result_explicit = self._run_chain(
            ExplicitEulerSystemSolver(), Q=Q, dt_h=dt_h, n_steps=n_steps
        )
        for k in range(10):
            implicit_val = float(result_implicit.phase_mol[f"cv{k}"]["liquid"]["S"][-1])
            explicit_val = float(result_explicit.phase_mol[f"cv{k}"]["liquid"]["S"][-1])
            # First-order agreement within 2% for small dt
            assert abs(implicit_val - explicit_val) < 0.02, (
                f"cv{k}: implicit={implicit_val:.6f}, explicit={explicit_val:.6f}"
            )


# ═══════════════════════════════════════════════════════════════════════
#  Controller events
# ═══════════════════════════════════════════════════════════════════════

class TestControllerEvents:

    def _make_counting_ctrl(self):
        from PyOMES.control.interfaces import ControllerBase
        from PyOMES.control.actions import ControlAction

        class Counter(ControllerBase):
            def __init__(self):
                self.count = 0

            def compute(self, state, dt_h):
                self.count += 1
                return ControlAction(controller_label="counter", target_cv_key="main")

        return Counter()

    def _run_with_ctrl(self, system_solver, ctrl, n_steps=5):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv}, controllers=[ctrl], system_solver=system_solver)
        sim.run(tau_h=0.1 * n_steps, n_steps=n_steps)
        return ctrl.count

    def test_controller_fires_every_step(self):
        """Without update_period_h, controller fires once per macro step."""
        from PyOMES.core import ImplicitTransportSystemSolver
        ctrl = self._make_counting_ctrl()
        n = self._run_with_ctrl(ImplicitTransportSystemSolver(), ctrl, n_steps=7)
        assert n == 7

    def test_same_count_as_explicit_euler(self):
        """ImplicitTransport fires controllers the same number of times as ExplicitEuler."""
        from PyOMES.core import ImplicitTransportSystemSolver, ExplicitEulerSystemSolver
        ctrl_i = self._make_counting_ctrl()
        ctrl_e = self._make_counting_ctrl()
        n_implicit = self._run_with_ctrl(ImplicitTransportSystemSolver(), ctrl_i, n_steps=6)
        n_explicit = self._run_with_ctrl(ExplicitEulerSystemSolver(), ctrl_e, n_steps=6)
        assert n_implicit == n_explicit == 6

    def test_no_chemistry_single_cv_runs_without_error(self):
        """Single CV with no links: solver degrades to pure-chemistry path."""
        from PyOMES.core import ImplicitTransportSystemSolver, Simulation
        cv = _make_liquid_cv("main", {"S": 2.0})
        sim = Simulation(cvs={"main": cv}, system_solver=ImplicitTransportSystemSolver())
        result = sim.run(tau_h=0.5, n_steps=5)
        # No transport, no chemistry → S should be unchanged
        final = float(result.phase_mol["main"]["liquid"]["S"][-1])
        assert final == pytest.approx(2.0, rel=1e-10)
