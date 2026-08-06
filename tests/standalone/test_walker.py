# -*- coding: utf-8 -*-
"""Integration tests for the reflective walker (C3).

Validates:
- String path dispatch: AttrSegment traversal, ListSelector (type/index),
  dict-leaf write, descriptor-leaf write, _set_unchecked fallback
- Typed-form dispatch: _DictItemRef, MutableScalar class-level, MutableDict class-level
- ParamPathError raised for unknown attribute, unknown class, empty collection
- ParamPathError raised for malformed path strings (propagated from parser)
- End-to-end: Simulation.run() with new-format params_changed strings
"""

import pytest

from PyOMES.chemistry import HenryPartition


def _hp(kH: float, dlnH: float = 0.0) -> HenryPartition:
    return HenryPartition(H_ref=kH * 1000.0 / 101325.0, dlnH=dlnH)


# ══════════════════════════════════════════════════════════════════════
#  Fixtures
# ══════════════════════════════════════════════════════════════════════

def _make_gas_liquid_cv():
    """CV with GasPhase + LiquidPhase + KineticGasLiquidLink + GasFeed."""
    from PyOMES.core import (
        ControlVolume, GasPhase, LiquidPhase, KineticGasLiquidLink,
    )
    from PyOMES.core.boundaries import GasFeed
    gas = GasPhase(n_mol={"O2": 0.5, "CO2": 0.05, "N2": 1.8}, V_L=0.3, T_K=308.15)
    liq = LiquidPhase(n_mol={"O2": 1e-4, "CO2": 1e-3}, V_L=1.0, T_K=308.15)
    link = KineticGasLiquidLink(
        gas_cv_key="main", gas_phase_key="gas",
        liquid_cv_key="main", liquid_phase_key="liquid",
        partition_models={"O2": _hp(1.3e-3), "CO2": _hp(3.4e-2)},
        kLa={"O2": 100.0, "CO2": 90.0},
    )
    feed = GasFeed(vvm_min=1.0, y={"O2": 0.21, "N2": 0.79})
    return ControlVolume(
        phases={"gas": gas, "liquid": liq},
        internal_interfaces=[link],
        boundaries=[feed],
        label="main",
    ), link, feed


def _make_running_sim(cv):
    from PyOMES.core import Simulation
    sim = Simulation(cvs={"main": cv})
    sim._context.is_running = True
    return sim


# ══════════════════════════════════════════════════════════════════════
#  String path — AttrSegment traversal
# ══════════════════════════════════════════════════════════════════════

class TestAttrTraversal:

    def test_phases_liquid_T_K(self):
        cv, link, feed = _make_gas_liquid_cv()
        sim = _make_running_sim(cv)
        try:
            sim._apply_param_change(cv, "phases.liquid.T_K", 310.0)
            assert cv.phases["liquid"].T_K == pytest.approx(310.0)
        finally:
            sim._context.is_running = False

    def test_phases_gas_T_K(self):
        cv, link, feed = _make_gas_liquid_cv()
        sim = _make_running_sim(cv)
        try:
            sim._apply_param_change(cv, "phases.gas.T_K", 312.0)
            assert cv.phases["gas"].T_K == pytest.approx(312.0)
        finally:
            sim._context.is_running = False


# ══════════════════════════════════════════════════════════════════════
#  String path — ListSelector (type-based)
# ══════════════════════════════════════════════════════════════════════

class TestListSelector:

    def test_kLa_O2_via_full_path(self):
        cv, link, feed = _make_gas_liquid_cv()
        sim = _make_running_sim(cv)
        try:
            sim._apply_param_change(
                cv,
                "internal_interfaces[KineticGasLiquidLink].kLa.O2",
                250.0,
            )
            assert link.kLa["O2"] == pytest.approx(250.0)
        finally:
            sim._context.is_running = False

    def test_kLa_CO2_via_full_path(self):
        cv, link, feed = _make_gas_liquid_cv()
        sim = _make_running_sim(cv)
        try:
            sim._apply_param_change(
                cv,
                "internal_interfaces[KineticGasLiquidLink].kLa.CO2",
                200.0,
            )
            assert link.kLa["CO2"] == pytest.approx(200.0)
        finally:
            sim._context.is_running = False

    def test_boundaries_gas_feed_vvm_min(self):
        cv, link, feed = _make_gas_liquid_cv()
        sim = _make_running_sim(cv)
        try:
            sim._apply_param_change(cv, "boundaries[GasFeed].vvm_min", 3.0)
            assert feed.vvm_min == pytest.approx(3.0)
        finally:
            sim._context.is_running = False

    def test_boundaries_gas_feed_y_O2(self):
        cv, link, feed = _make_gas_liquid_cv()
        sim = _make_running_sim(cv)
        try:
            sim._apply_param_change(cv, "boundaries[GasFeed].y.O2", 0.40)
            assert feed.y["O2"] == pytest.approx(0.40)
        finally:
            sim._context.is_running = False

    def test_index_selector_zero(self):
        """boundaries[0] selects the first boundary by numeric index."""
        cv, link, feed = _make_gas_liquid_cv()
        sim = _make_running_sim(cv)
        try:
            sim._apply_param_change(cv, "boundaries[0].vvm_min", 5.0)
            assert feed.vvm_min == pytest.approx(5.0)
        finally:
            sim._context.is_running = False


# ══════════════════════════════════════════════════════════════════════
#  Typed form dispatch
# ══════════════════════════════════════════════════════════════════════

class TestTypedFormDispatch:

    def test_dict_item_ref_kLa_O2(self):
        """Typed _DictItemRef path using KineticGasLiquidLink.kLa (MutableDict in C4)."""
        from PyOMES.core import KineticGasLiquidLink
        from PyOMES.control.descriptors import _DictItemRef
        cv, link, feed = _make_gas_liquid_cv()
        sim = _make_running_sim(cv)
        try:
            path = KineticGasLiquidLink.kLa["O2"]
            assert isinstance(path, _DictItemRef)
            sim._apply_param_change(cv, path, 300.0)
            assert link.kLa["O2"] == pytest.approx(300.0)
        finally:
            sim._context.is_running = False

    def test_dict_item_ref_gas_feed_y(self):
        """Typed _DictItemRef path using GasFeed.y."""
        from PyOMES.core.boundaries import GasFeed
        from PyOMES.control.descriptors import _DictItemRef
        cv, link, feed = _make_gas_liquid_cv()
        sim = _make_running_sim(cv)
        try:
            path = GasFeed.y["O2"]
            assert isinstance(path, _DictItemRef)
            sim._apply_param_change(cv, path, 0.40)
            assert feed.y["O2"] == pytest.approx(0.40)
        finally:
            sim._context.is_running = False

    def test_mutable_scalar_T_K(self):
        from PyOMES.core.phases import GasPhase
        from PyOMES.control.descriptors import MutableScalar
        # GasPhase.T_K is NOT yet a MutableScalar in C3 (migration is C4).
        # Test with a synthetic host instead to verify the dispatch code path.
        from PyOMES.control.descriptors import MutableScalar

        class _Host:
            T_K = MutableScalar(float, positive=True)
            _context = None
            def __init__(self, T_K):
                self.T_K = T_K

        host = _Host(T_K=300.0)
        # Use the walker directly by patching a CV search.
        from PyOMES.core import Simulation, ControlVolume, LiquidPhase
        from PyOMES.core.lifecycle import RunContext
        liq = LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=300.0)
        cv = ControlVolume(phases={"liquid": liq})
        cv.boundaries = [host]  # put _Host in boundaries

        sim = Simulation(cvs={"main": cv})
        sim._context.is_running = True
        try:
            path = _Host.T_K  # class-level → MutableScalar descriptor
            assert isinstance(path, MutableScalar)
            sim._apply_param_change(cv, path, 320.0)
            assert host.T_K == pytest.approx(320.0)
        finally:
            sim._context.is_running = False


# ══════════════════════════════════════════════════════════════════════
#  ParamPathError cases
# ══════════════════════════════════════════════════════════════════════

class TestParamPathErrors:

    def test_malformed_path_raises(self):
        from PyOMES.control.param_path import ParamPathError
        cv, link, feed = _make_gas_liquid_cv()
        sim = _make_running_sim(cv)
        try:
            with pytest.raises(ParamPathError):
                sim._apply_param_change(cv, "frobnicate.the.widget", 42.0)
        finally:
            sim._context.is_running = False

    def test_unknown_class_in_selector_raises(self):
        from PyOMES.control.param_path import ParamPathError
        cv, link, feed = _make_gas_liquid_cv()
        sim = _make_running_sim(cv)
        try:
            with pytest.raises(ParamPathError):
                sim._apply_param_change(
                    cv,
                    "internal_interfaces[NoSuchClass].kLa.O2",
                    100.0,
                )
        finally:
            sim._context.is_running = False

    def test_class_not_in_cv_raises(self):
        """ListSelector matching succeeds but no instance of that class exists."""
        from PyOMES.control.param_path import ParamPathError
        # CV with no KineticGasLiquidLink
        from PyOMES.core import ControlVolume, LiquidPhase, Simulation
        liq = LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=300.0)
        cv = ControlVolume(phases={"liquid": liq})
        sim = Simulation(cvs={"main": cv})
        sim._context.is_running = True
        try:
            with pytest.raises(ParamPathError):
                sim._apply_param_change(
                    cv,
                    "internal_interfaces[KineticGasLiquidLink].kLa.O2",
                    100.0,
                )
        finally:
            sim._context.is_running = False

    def test_unknown_leaf_attribute_raises(self):
        from PyOMES.control.param_path import ParamPathError
        cv, link, feed = _make_gas_liquid_cv()
        sim = _make_running_sim(cv)
        try:
            with pytest.raises(ParamPathError):
                # vvm_max doesn't exist and has no _set_vvm_max_unchecked
                sim._apply_param_change(
                    cv, "boundaries[GasFeed].vvm_max", 99.0,
                )
        finally:
            sim._context.is_running = False

    def test_invalid_path_type_raises(self):
        from PyOMES.control.param_path import ParamPathError
        cv, link, feed = _make_gas_liquid_cv()
        sim = _make_running_sim(cv)
        try:
            with pytest.raises(ParamPathError):
                sim._apply_param_change(cv, 12345, 1.0)
        finally:
            sim._context.is_running = False


# ══════════════════════════════════════════════════════════════════════
#  End-to-end: Simulation.run() with new-format params_changed
# ══════════════════════════════════════════════════════════════════════

class TestEndToEnd:

    def test_kLa_updated_end_to_end(self):
        """Simulation.run() dispatches kLa via the new path and link updates."""
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction

        cv, link, feed = _make_gas_liquid_cv()

        class SetkLa:
            def compute(self, state, dt_h):
                return ControlAction(
                    target_cv_key=state.cv_key,
                    params_changed={
                        "internal_interfaces[KineticGasLiquidLink].kLa.O2": 500.0,
                    },
                )

        sim = Simulation(cvs={"main": cv}, controllers=[SetkLa()])
        sim.run(tau_h=0.01, n_steps=1)
        assert link.kLa["O2"] == pytest.approx(500.0)

    def test_vvm_and_y_updated_end_to_end(self):
        """Simulation.run() dispatches boundaries[GasFeed] paths correctly."""
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction

        cv, link, feed = _make_gas_liquid_cv()

        class SetFeed:
            def compute(self, state, dt_h):
                return ControlAction(
                    target_cv_key=state.cv_key,
                    params_changed={
                        "boundaries[GasFeed].vvm_min": 3.5,
                        "boundaries[GasFeed].y.O2": 0.30,
                    },
                )

        sim = Simulation(cvs={"main": cv}, controllers=[SetFeed()])
        sim.run(tau_h=0.01, n_steps=1)
        assert feed.vvm_min == pytest.approx(3.5)
        assert feed.y["O2"] == pytest.approx(0.30)

    def test_unknown_path_raises_during_run(self):
        """ParamPathError propagates out of Simulation.run() (Q3)."""
        from PyOMES.core import Simulation, ControlVolume, LiquidPhase
        from PyOMES.control.actions import ControlAction
        from PyOMES.control.param_path import ParamPathError

        class BadPath:
            def compute(self, state, dt_h):
                return ControlAction(
                    target_cv_key=state.cv_key,
                    params_changed={"kLa.O2": 100.0},  # old format → error
                )

        liq = LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=300.0)
        cv = ControlVolume(phases={"liquid": liq})
        sim = Simulation(cvs={"main": cv}, controllers=[BadPath()])
        with pytest.raises(ParamPathError):
            sim.run(tau_h=0.01, n_steps=1)
