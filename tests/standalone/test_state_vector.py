# -*- coding: utf-8 -*-
"""Tests for src/core/state_vector.py (STEP_SOLVER_INTERFACE_REFINEMENT.md
item 4 -- unified state-vector packing). Checkpoint 7 of
STEP_SOLVER_REFINEMENT_CHECKLIST.md.

Covers the StateVector class directly (pack/unpack roundtrip,
exclude_species, ctrl_list extension, index_map) plus a regression
check that the five thin-wrapper functions in system_solver.py
(_pack_state/_unpack_state/_state_index_map/_pack_extended_state/
_unpack_extended_state) still produce identical output to what they
did before delegating to StateVector.
"""

from __future__ import annotations

import numpy as np
import pytest


def _make_sim(n_mol_by_cv):
    """{cv_key: {phase_key: {species: mol}}} -> Simulation."""
    from PyOMES.core import ControlVolume, LiquidPhase, GasPhase, Simulation
    cvs = {}
    for cv_key, phases in n_mol_by_cv.items():
        phase_objs = {}
        for pk, n_mol in phases.items():
            cls = GasPhase if pk == "gas" else LiquidPhase
            phase_objs[pk] = cls(n_mol=dict(n_mol), V_L=1.0, T_K=300.0)
        cvs[cv_key] = ControlVolume(phases=phase_objs, label=cv_key)
    return Simulation(cvs=cvs)


class TestStateVectorPackUnpack:

    def test_pack_unpack_roundtrip(self):
        from PyOMES.core.state_vector import StateVector
        sim = _make_sim({"A": {"liquid": {"S": 1.0, "X": 2.0}}})
        sv = StateVector(sim.cvs)
        y = sv.pack()
        assert y.shape == (2,)
        y2 = y * 2.0
        sv.unpack(y2)
        assert sim.cvs["A"].phases["liquid"].n_mol["S"] == pytest.approx(2.0)
        assert sim.cvs["A"].phases["liquid"].n_mol["X"] == pytest.approx(4.0)

    def test_ordering_is_cv_then_phase_then_species_alphabetical(self):
        from PyOMES.core.state_vector import StateVector
        sim = _make_sim({
            "B": {"liquid": {"Z": 1.0, "A": 2.0}},
            "A": {"liquid": {"Y": 3.0}},
        })
        sv = StateVector(sim.cvs)
        # dict insertion order: B, then A (as passed to _make_sim)
        assert sv.index_map == {
            ("B", "liquid", "A"): 0,
            ("B", "liquid", "Z"): 1,
            ("A", "liquid", "Y"): 2,
        }

    def test_multi_phase_alphabetical_ordering(self):
        from PyOMES.core.state_vector import StateVector
        sim = _make_sim({"A": {"gas": {"O2": 1.0}, "liquid": {"S": 2.0}}})
        sv = StateVector(sim.cvs)
        assert list(sv.index_map) == [("A", "gas", "O2"), ("A", "liquid", "S")]

    def test_exclude_species(self):
        from PyOMES.core.state_vector import StateVector
        sim = _make_sim({"A": {"liquid": {"S": 1.0, "H+": 1e-7}}})
        sv = StateVector(sim.cvs, exclude_species=frozenset({"H+"}))
        assert sv.n_cv == 1
        assert ("A", "liquid", "H+") not in sv.index_map

    def test_unpack_without_floor_allows_negative(self):
        from PyOMES.core.state_vector import StateVector
        sim = _make_sim({"A": {"liquid": {"S": 1.0}}})
        sv = StateVector(sim.cvs)
        sv.unpack(np.array([-5.0]), floor=False)
        assert sim.cvs["A"].phases["liquid"].n_mol["S"] == pytest.approx(-5.0)

    def test_unpack_with_floor_clamps_at_zero(self):
        from PyOMES.core.state_vector import StateVector
        sim = _make_sim({"A": {"liquid": {"S": 1.0}}})
        sv = StateVector(sim.cvs)
        sv.unpack(np.array([-5.0]), floor=True)
        assert sim.cvs["A"].phases["liquid"].n_mol["S"] == 0.0


class TestStateVectorControllerExtension:

    def _make_ctrl(self, value=0.0):
        from PyOMES.control.interfaces import ControllerBase

        class _Ctrl(ControllerBase):
            def __init__(self):
                self.integral = value

            def differential_state(self):
                return {"integral": self.integral}

            def set_state(self, state):
                self.integral = state["integral"]

        return _Ctrl()

    def test_pack_appends_controller_state_after_cv_block(self):
        from PyOMES.core.state_vector import StateVector
        sim = _make_sim({"A": {"liquid": {"S": 1.0}}})
        ctrl = self._make_ctrl(value=42.0)
        sv = StateVector(sim.cvs, ctrl_list=[(ctrl, ["integral"])])
        y = sv.pack()
        assert sv.n_cv == 1
        assert sv.n_total == 2
        assert y[0] == pytest.approx(1.0)
        assert y[1] == pytest.approx(42.0)

    def test_unpack_writes_back_to_controller(self):
        from PyOMES.core.state_vector import StateVector
        sim = _make_sim({"A": {"liquid": {"S": 1.0}}})
        ctrl = self._make_ctrl(value=0.0)
        sv = StateVector(sim.cvs, ctrl_list=[(ctrl, ["integral"])])
        sv.unpack(np.array([9.0, 99.0]))
        assert sim.cvs["A"].phases["liquid"].n_mol["S"] == pytest.approx(9.0)
        assert ctrl.integral == pytest.approx(99.0)


class TestThinWrapperRegression:
    """The five system_solver.py functions now delegate to StateVector
    internally; confirm their public signatures/output are unaffected."""

    def test_pack_state_unpack_state_roundtrip(self):
        from PyOMES.core.system_solver import _pack_state, _unpack_state
        sim = _make_sim({"A": {"liquid": {"S": 1.0, "X": 2.0}}})
        y = _pack_state(sim)
        _unpack_state(y * 3.0, sim)
        assert sim.cvs["A"].phases["liquid"].n_mol["S"] == pytest.approx(3.0)
        assert sim.cvs["A"].phases["liquid"].n_mol["X"] == pytest.approx(6.0)

    def test_state_index_map(self):
        from PyOMES.core.system_solver import _state_index_map
        sim = _make_sim({"A": {"liquid": {"S": 1.0}}})
        idx = _state_index_map(sim)
        assert idx == {("A", "liquid", "S"): 0}

    def test_pack_extended_state_with_empty_ctrl_list(self):
        from PyOMES.core.system_solver import _pack_extended_state, _pack_state
        sim = _make_sim({"A": {"liquid": {"S": 1.0}}})
        assert _pack_extended_state(sim, []).tolist() == _pack_state(sim).tolist()

    def test_unpack_extended_state_ignores_stale_n_cv_arg(self):
        """n_cv is accepted for signature compatibility but StateVector
        recomputes its own offset -- a wrong/stale n_cv must not corrupt
        the result."""
        from PyOMES.core.system_solver import _unpack_extended_state
        sim = _make_sim({"A": {"liquid": {"S": 1.0}}})
        _unpack_extended_state(np.array([7.0]), sim, [], n_cv=999)
        assert sim.cvs["A"].phases["liquid"].n_mol["S"] == pytest.approx(7.0)
