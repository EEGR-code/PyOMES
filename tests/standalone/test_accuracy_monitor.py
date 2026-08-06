# -*- coding: utf-8 -*-
"""Tests for ``PyOMES.monitoring.accuracy`` and ``PyOMES.config``.

Covers:
- ``WarningConfig`` defaults + classmethod presets.
- ``VLSIM_WARNINGS`` env-var loading (silent / verbose /
  production / empty / unknown).
- ``AccuracyMonitor`` per-check methods (positive / negative /
  no-op input handling for each of six checks).
- Throttle modes (once / first_N / always / silent).
- Integration: ``ControlVolume.__init__`` attaches a monitor;
  the speciation property solver and the snapshot solver call
  into it.
"""

from __future__ import annotations

import os
import warnings

import pytest

import PyOMES
from PyOMES import AccuracyWarning, WarningConfig
from PyOMES.config import _warning_config_from_env
from PyOMES.monitoring.accuracy import (
    AccuracyMonitor,
    _summary_counter,
    print_accuracy_summary,
    reset_accuracy_summary,
)


# ── Fixtures ──────────────────────────────────────────────────────


@pytest.fixture(autouse=True)
def _restore_warning_config():
    """Reset thresholds, throttle, and the module-level summary
    counter to their factory defaults between tests so cross-test
    state doesn't leak."""
    snapshot = PyOMES.config.warnings
    PyOMES.config.warnings = WarningConfig()
    reset_accuracy_summary()
    try:
        yield
    finally:
        PyOMES.config.warnings = snapshot
        reset_accuracy_summary()


def _aw(buf):
    """Count AccuracyWarning instances in a captured warning buffer."""
    return sum(1 for w in buf if issubclass(w.category, AccuracyWarning))


# ════════════════════════════════════════════════════════════════════
#  Config + env-var
# ════════════════════════════════════════════════════════════════════


class TestWarningConfig:

    def test_default_thresholds(self):
        cfg = WarningConfig()
        assert cfg.pH_change_threshold == pytest.approx(0.3)
        assert cfg.newton_iters_threshold == 15
        assert cfg.charge_residual_threshold == pytest.approx(1e-8)
        assert cfg.dt_over_tau_min_threshold == pytest.approx(0.1)
        assert cfg.scipy_rejection_threshold == pytest.approx(0.3)
        assert cfg.ionic_strength_ideal_threshold == pytest.approx(0.10)
        assert cfg.ionic_strength_davies_threshold == pytest.approx(0.50)
        assert cfg.throttle == "once"
        assert cfg.first_N == 10

    def test_silent_preset(self):
        assert WarningConfig.silent().throttle == "silent"

    def test_verbose_preset(self):
        assert WarningConfig.verbose().throttle == "always"

    def test_production_preset(self):
        prod = WarningConfig.production()
        assert prod.throttle == "first_N"
        assert prod.first_N == 3


class TestEnvVar:

    def test_silent(self, monkeypatch):
        monkeypatch.setenv("VLSIM_WARNINGS", "silent")
        assert _warning_config_from_env().throttle == "silent"

    def test_verbose_case_insensitive(self, monkeypatch):
        monkeypatch.setenv("VLSIM_WARNINGS", "VERBOSE")
        assert _warning_config_from_env().throttle == "always"

    def test_production_whitespace_trimmed(self, monkeypatch):
        monkeypatch.setenv("VLSIM_WARNINGS", "  production  ")
        cfg = _warning_config_from_env()
        assert cfg.throttle == "first_N"
        assert cfg.first_N == 3

    def test_empty_string_falls_through(self, monkeypatch):
        monkeypatch.setenv("VLSIM_WARNINGS", "")
        assert _warning_config_from_env().throttle == "once"  # default

    def test_unset_uses_default(self, monkeypatch):
        monkeypatch.delenv("VLSIM_WARNINGS", raising=False)
        assert _warning_config_from_env().throttle == "once"

    def test_unknown_value_warns_and_defaults(self, monkeypatch):
        monkeypatch.setenv("VLSIM_WARNINGS", "bogus")
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            cfg = _warning_config_from_env()
        assert cfg.throttle == "once"  # fell back to default
        msgs = [str(w.message) for w in buf if issubclass(w.category, UserWarning)]
        assert any("VLSIM_WARNINGS" in m and "bogus" in m for m in msgs)


# ════════════════════════════════════════════════════════════════════
#  Per-check methods
# ════════════════════════════════════════════════════════════════════


class TestCheckPHJump:

    def test_first_call_no_warn(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.pH_change_threshold = 0.01
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_pH_jump(7.0)
        assert _aw(buf) == 0
        assert m.last_pH == pytest.approx(7.0)
        assert m.step_count == 1

    def test_jump_above_threshold_warns(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.pH_change_threshold = 0.1
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_pH_jump(7.0)
            m.check_pH_jump(7.5)
        assert _aw(buf) == 1

    def test_none_and_nan_no_op(self):
        PyOMES.config.warnings.throttle = "always"
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_pH_jump(None)
            m.check_pH_jump(float("nan"))
        assert _aw(buf) == 0


class TestCheckNewtonIters:

    def test_above_threshold_warns(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.newton_iters_threshold = 15
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_newton_iters(20)
        assert _aw(buf) == 1

    def test_at_threshold_no_warn(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.newton_iters_threshold = 15
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_newton_iters(15)
        assert _aw(buf) == 0

    def test_none_negative_nonint_no_op(self):
        PyOMES.config.warnings.throttle = "always"
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_newton_iters(None)
            m.check_newton_iters(-5)
            m.check_newton_iters("not-a-number")
        assert _aw(buf) == 0


class TestCheckChargeResidual:

    def test_above_threshold_warns(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.charge_residual_threshold = 1e-8
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_charge_residual(1e-6)
            m.check_charge_residual(-1e-6)  # absolute value
        assert _aw(buf) == 2

    def test_below_threshold_no_warn(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.charge_residual_threshold = 1e-6
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_charge_residual(1e-9)
        assert _aw(buf) == 0

    def test_none_and_nan_no_op(self):
        PyOMES.config.warnings.throttle = "always"
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_charge_residual(None)
            m.check_charge_residual(float("nan"))
        assert _aw(buf) == 0


class TestCheckNegativeMole:

    def test_negative_species_warns(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.negative_mole_tolerance = 1e-9
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_negative_mole({"S": -1e-3, "X": 5.0})
        assert _aw(buf) == 1

    def test_within_tolerance_no_warn(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.negative_mole_tolerance = 1e-6
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_negative_mole({"S": -1e-9, "X": 5.0})
        assert _aw(buf) == 0

    def test_all_nonnegative_no_warn(self):
        PyOMES.config.warnings.throttle = "always"
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_negative_mole({"S": 1.0, "X": 0.0})
        assert _aw(buf) == 0

    def test_message_names_worst_offender(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.negative_mole_tolerance = 1e-9
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_negative_mole({"S": -1e-3, "T": -5.0, "X": 5.0})
        msgs = [str(w.message) for w in buf if issubclass(w.category, AccuracyWarning)]
        assert len(msgs) == 1
        assert "'T'" in msgs[0]  # T is the most negative


class TestCheckClampInvoked:

    def test_changed_rate_warns(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.clamp_invoked_tolerance = 1e-9
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_clamp_invoked({"S": -12.3}, {"S": -8.1}, dt_h=0.5)
        assert _aw(buf) == 1

    def test_unchanged_rate_no_warn(self):
        PyOMES.config.warnings.throttle = "always"
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_clamp_invoked({"S": -12.3}, {"S": -12.3}, dt_h=0.5)
        assert _aw(buf) == 0

    def test_within_tolerance_no_warn(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.clamp_invoked_tolerance = 1e-3
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_clamp_invoked({"S": -12.3}, {"S": -12.3000001}, dt_h=0.5)
        assert _aw(buf) == 0

    def test_message_names_species_and_rates(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.clamp_invoked_tolerance = 1e-9
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_clamp_invoked({"S_lcfa": -12.3}, {"S_lcfa": -8.1}, dt_h=0.5)
        msgs = [str(w.message) for w in buf if issubclass(w.category, AccuracyWarning)]
        assert len(msgs) == 1
        assert "S_lcfa" in msgs[0]


class TestCheckDtVsTauMin:

    def test_ratio_above_threshold_warns(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.dt_over_tau_min_threshold = 0.1
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_dt_vs_tau_min(dt_h=1.0, tau_min_h=5.0)  # ratio 0.2
        assert _aw(buf) == 1

    def test_tau_none_no_op(self):
        PyOMES.config.warnings.throttle = "always"
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_dt_vs_tau_min(dt_h=1.0, tau_min_h=None)
        assert _aw(buf) == 0

    def test_negative_or_zero_no_op(self):
        PyOMES.config.warnings.throttle = "always"
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_dt_vs_tau_min(dt_h=-1.0, tau_min_h=5.0)
            m.check_dt_vs_tau_min(dt_h=1.0, tau_min_h=0.0)
        assert _aw(buf) == 0


class TestCheckScipyRejections:

    def test_high_rejection_ratio_warns(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.scipy_rejection_threshold = 0.3
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_scipy_rejections(n_accepted=10, n_attempted=20)
        assert _aw(buf) == 1

    def test_low_rejection_no_warn(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.scipy_rejection_threshold = 0.3
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_scipy_rejections(n_accepted=10, n_attempted=11)
        assert _aw(buf) == 0

    def test_trivially_short_run_no_op(self):
        PyOMES.config.warnings.throttle = "always"
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_scipy_rejections(n_accepted=1, n_attempted=1)
            m.check_scipy_rejections(n_accepted=0, n_attempted=5)
        assert _aw(buf) == 0


class TestCheckIonicStrength:

    def test_davies_above_threshold(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.ionic_strength_davies_threshold = 0.5
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_ionic_strength(0.6, activity_model="davies", use_activity=True)
        assert _aw(buf) == 1
        assert _summary_counter["ionic_strength_davies"] == 1

    def test_ideal_assumption_uses_ideal_threshold(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.ionic_strength_ideal_threshold = 0.10
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            # use_activity=False forces ideal regime regardless of model name
            m.check_ionic_strength(0.15, activity_model="davies", use_activity=False)
        assert _aw(buf) == 1
        assert _summary_counter["ionic_strength_ideal"] == 1

    def test_unknown_activity_model_no_op(self):
        PyOMES.config.warnings.throttle = "always"
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_ionic_strength(10.0, activity_model="some_future_model", use_activity=True)
        assert _aw(buf) == 0

    def test_none_no_op(self):
        PyOMES.config.warnings.throttle = "always"
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m.check_ionic_strength(None, activity_model="davies", use_activity=True)
        assert _aw(buf) == 0


# ════════════════════════════════════════════════════════════════════
#  Throttle modes
# ════════════════════════════════════════════════════════════════════


class TestThrottle:

    def test_once_emits_one_per_monitor_per_category(self):
        PyOMES.config.warnings.throttle = "once"
        PyOMES.config.warnings.pH_change_threshold = 0.0
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            for v in (7.0, 7.1, 7.2, 7.3, 7.4):
                m.check_pH_jump(v)
        assert _aw(buf) == 1

    def test_once_is_per_monitor_not_global(self):
        PyOMES.config.warnings.throttle = "once"
        PyOMES.config.warnings.pH_change_threshold = 0.0
        m1, m2 = AccuracyMonitor(), AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m1.check_pH_jump(7.0); m1.check_pH_jump(7.5)  # 1 warn
            m1.check_pH_jump(8.0)                          # suppressed by once
            m2.check_pH_jump(7.0); m2.check_pH_jump(7.5)  # 1 warn
            m2.check_pH_jump(8.0)                          # suppressed
        assert _aw(buf) == 2  # one per monitor

    def test_first_N_caps_globally(self):
        # first_N is a global cap via the module-level summary counter.
        PyOMES.config.warnings.throttle = "first_N"
        PyOMES.config.warnings.first_N = 2
        PyOMES.config.warnings.pH_change_threshold = 0.0
        m1, m2 = AccuracyMonitor(), AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            m1.check_pH_jump(7.0); m1.check_pH_jump(7.5)  # 1 warn
            m2.check_pH_jump(7.0); m2.check_pH_jump(7.5)  # 1 warn — total 2, cap hit
            m1.check_pH_jump(8.0)                          # suppressed by cap
            m2.check_pH_jump(8.0)                          # suppressed by cap
        assert _aw(buf) == 2

    def test_always_emits_unconditionally(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.pH_change_threshold = 0.0
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            for v in (7.0, 7.1, 7.2, 7.3):
                m.check_pH_jump(v)
        # First call no-ops (no prior); next 3 each cross the 0.0 threshold.
        assert _aw(buf) == 3

    def test_silent_suppresses_all(self):
        PyOMES.config.warnings.throttle = "silent"
        PyOMES.config.warnings.pH_change_threshold = 0.0
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            for v in (7.0, 7.5, 8.0):
                m.check_pH_jump(v)
        assert _aw(buf) == 0


class TestSummaryHelpers:

    def test_summary_counter_increments_on_emission(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.pH_change_threshold = 0.0
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            m.check_pH_jump(7.0); m.check_pH_jump(7.1); m.check_pH_jump(7.2)
        assert _summary_counter["pH_change"] == 2

    def test_reset_clears_monitor_state_only(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.pH_change_threshold = 0.0
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            m.check_pH_jump(7.0); m.check_pH_jump(7.1)
        m.reset()
        assert m.last_pH is None
        assert m.step_count == 0
        assert m._warned_categories == set()
        # Counter persists — separate concept.
        assert _summary_counter["pH_change"] == 1

    def test_reset_accuracy_summary_clears_counter(self):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.pH_change_threshold = 0.0
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            m.check_pH_jump(7.0); m.check_pH_jump(7.1)
        assert _summary_counter["pH_change"] == 1
        reset_accuracy_summary()
        assert _summary_counter["pH_change"] == 0

    def test_print_summary_no_op_when_empty(self, capsys):
        print_accuracy_summary()
        captured = capsys.readouterr()
        assert captured.out == ""

    def test_print_summary_includes_categories(self, capsys):
        PyOMES.config.warnings.throttle = "always"
        PyOMES.config.warnings.pH_change_threshold = 0.0
        m = AccuracyMonitor()
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            m.check_pH_jump(7.0); m.check_pH_jump(7.5)
        print_accuracy_summary()
        captured = capsys.readouterr()
        assert "AccuracyWarning summary" in captured.out
        assert "pH_change=1" in captured.out


# ════════════════════════════════════════════════════════════════════
#  Integration with ControlVolume + solvers
# ════════════════════════════════════════════════════════════════════


class TestCVIntegration:

    def test_control_volume_attaches_monitor(self):
        from PyOMES.core import ControlVolume, GasPhase, LiquidPhase
        gas = GasPhase({"O2": 1.0, "N2": 3.0}, V_L=200.0, T_K=305.15)
        liq = LiquidPhase({"CO2": 0.01}, V_L=800.0, T_K=305.15)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})
        assert isinstance(cv._accuracy_monitor, AccuracyMonitor)
        assert cv._tau_min_h_estimate is None  # deferred — see checklist

    def test_monitor_attached_to_reaction_system_engine(self):
        """state-unification C4d: CV.__init__ calls
        cv.reaction_system.attach_monitor(self._accuracy_monitor),
        which propagates to the engine. The engine-side accuracy
        checks then route to the CV's monitor without any
        property_solvers list."""
        from PyOMES.core import ControlVolume, LiquidPhase
        from PyOMES.reactions import ReactionSystem

        class _StubEngine:
            activity_model = "davies"
            use_activity = True
            def solve(self, **kw):
                return {"IonicStrength": 0.1, "pH": 7.0, "logH": -7.0, "alphas": {}}

        liq = LiquidPhase({"H+": 1e-7}, V_L=1.0, T_K=298.15)
        system = ReactionSystem([], label="stub_system")
        system.attach_engine(_StubEngine())
        cv = ControlVolume(phases={"liquid": liq}, reaction_system=system)
        assert system._engine._accuracy_monitor is cv._accuracy_monitor


class TestEngineIonicStrength:
    """state-unification C4d: the ionic-strength regime check fires
    when the engine is attached to a ReactionSystem with a monitor.
    The legacy engine.py RuntimeWarning was replaced by AccuracyWarning
    in chemistry-unification-4; this surface survives the C4d
    property-solver collapse since it lives on the engine directly.
    """

    def test_high_ionic_strength_emits_accuracy_warning(self):
        from PyOMES.core import LiquidPhase
        from PyOMES.reactions import ReactionSystem

        class _HighIEngine:
            activity_model = "davies"
            use_activity = True
            _accuracy_monitor = None  # set via attach_monitor below
            def solve(self, **kw):
                # Mirror SpeciationPropertySolver's monitor hook
                out = {"IonicStrength": 0.8, "pH": 7.0, "logH": -7.0, "alphas": {}}
                if self._accuracy_monitor is not None:
                    self._accuracy_monitor.check_ionic_strength(
                        out["IonicStrength"],
                        activity_model=self.activity_model,
                        use_activity=self.use_activity,
                    )
                return out

        PyOMES.config.warnings.throttle = "always"
        system = ReactionSystem([], label="stub")
        system.attach_engine(_HighIEngine())
        system.attach_monitor(AccuracyMonitor())
        liq = LiquidPhase({"H+": 1e-7}, V_L=1.0, T_K=298.15)
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            system._engine.solve(phases={"liquid": liq}, T_K=298.15)
        assert _aw(buf) == 1
        assert _summary_counter["ionic_strength_davies"] == 1

    def test_standalone_engine_no_monitor_silent(self):
        """An engine called directly (no monitor) is silent — the
        regime warning only fires when a monitor is attached."""
        from PyOMES.core import LiquidPhase

        class _HighIEngine:
            activity_model = "davies"
            use_activity = True
            _accuracy_monitor = None
            def solve(self, **kw):
                return {"IonicStrength": 0.8, "pH": 7.0, "logH": -7.0, "alphas": {}}

        engine = _HighIEngine()
        liq = LiquidPhase({"H+": 1e-7}, V_L=1.0, T_K=298.15)
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            engine.solve(phases={"liquid": liq}, T_K=298.15)
        assert _aw(buf) == 0
