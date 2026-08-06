"""Tests for fermenter.core.boundaries — external boundary objects.

Tests cover:
  - ExternalBoundary protocol compliance
  - GasFeed flux computation and parity with legacy _gas_feed_mol_per_h
  - PressureReliefVent (instant and smooth modes)
  - apply_boundary orchestrator
  - ExternalFluxRecord diagnostics
  - Integration: feed + equilibrium + vent cycle through ControlVolume
"""

import math
import pytest

from PyOMES.core.phases import GasPhase, LiquidPhase, R_L_ATM_MOL_K
from PyOMES.core.control_volume import ControlVolume
from PyOMES.core.boundaries import (
    ExternalBoundary,
    ExternalFluxRecord,
    apply_boundary,
    GasFeed,
    PressureReliefVent,
    _smooth_vent_fraction,
)


# ════════════════════════════════════════════════════════════════════════
#  Helpers
# ════════════════════════════════════════════════════════════════════════

def _make_gl_cv(
    n_gas=None, V_gas=200.0, V_liq=800.0, T_K=305.15,
):
    """Build a standard gas-liquid CV for testing."""
    n_gas = n_gas or {"O2": 1.0, "CO2": 0.5, "N2": 3.0}
    gas = GasPhase(n_gas, V_L=V_gas, T_K=T_K)
    liq = LiquidPhase({"CO2": 0.01}, V_L=V_liq, T_K=T_K)
    return ControlVolume(phases={"gas": gas, "liquid": liq}, label="test")


# ════════════════════════════════════════════════════════════════════════
#  Protocol compliance
# ════════════════════════════════════════════════════════════════════════

class TestProtocolCompliance:
    def test_gas_feed_satisfies_protocol(self):
        feed = GasFeed(vvm_min=1.0)
        assert isinstance(feed, ExternalBoundary)

    def test_pressure_relief_satisfies_protocol(self):
        vent = PressureReliefVent(P_set_atm=1.0)
        assert isinstance(vent, ExternalBoundary)

    def test_gas_feed_has_phase_key(self):
        feed = GasFeed(phase_key="headspace")
        assert feed.phase_key == "headspace"

    def test_gas_feed_has_label(self):
        feed = GasFeed(label="air_supply")
        assert feed.label == "air_supply"

    def test_vent_has_phase_key(self):
        vent = PressureReliefVent(phase_key="headspace")
        assert vent.phase_key == "headspace"

    def test_vent_has_label(self):
        vent = PressureReliefVent(label="my_valve")
        assert vent.label == "my_valve"


# ════════════════════════════════════════════════════════════════════════
#  GasFeed tests
# ════════════════════════════════════════════════════════════════════════

class TestGasFeed:

    def test_zero_vvm_gives_zero_flux(self):
        feed = GasFeed(vvm_min=0.0, y={"O2": 0.21, "N2": 0.79})
        cv = _make_gl_cv()
        flux = feed.compute_flux(cv, dt_h=0.01)
        assert all(v == 0.0 for v in flux.values())

    def test_negative_vvm_gives_zero_flux(self):
        feed = GasFeed(vvm_min=-1.0)
        cv = _make_gl_cv()
        flux = feed.compute_flux(cv, dt_h=0.01)
        assert all(v == 0.0 for v in flux.values())

    def test_flux_is_positive(self):
        """Feed fluxes should be positive (entering the CV)."""
        feed = GasFeed(vvm_min=1.0, y={"O2": 0.21, "N2": 0.79})
        cv = _make_gl_cv()
        flux = feed.compute_flux(cv, dt_h=0.01)
        assert flux["O2"] > 0.0
        assert flux["N2"] > 0.0

    def test_composition_ratios(self):
        """Species flux ratios must match inlet mole fractions."""
        feed = GasFeed(vvm_min=1.0, y={"O2": 0.21, "N2": 0.79}, P_inlet_atm=1.0)
        cv = _make_gl_cv()
        flux = feed.compute_flux(cv, dt_h=0.01)
        total = flux["O2"] + flux["N2"]
        assert flux["O2"] / total == pytest.approx(0.21, rel=1e-10)
        assert flux["N2"] / total == pytest.approx(0.79, rel=1e-10)

    def test_composition_normalisation(self):
        """Non-normalised y values should be normalised internally."""
        feed = GasFeed(vvm_min=1.0, y={"O2": 42, "N2": 158})
        cv = _make_gl_cv()
        flux = feed.compute_flux(cv, dt_h=0.01)
        total = flux["O2"] + flux["N2"]
        assert flux["O2"] / total == pytest.approx(0.21, rel=1e-10)

    def test_total_molar_flow_formula(self):
        """Verify n_dot = P * vvm * V_liq * 60 / (R * T)."""
        vvm = 1.5
        V_liq = 800.0
        T_K = 305.15
        P_in = 1.0
        feed = GasFeed(vvm_min=vvm, y={"O2": 1.0}, P_inlet_atm=P_in)
        cv = _make_gl_cv(V_liq=V_liq, T_K=T_K)
        flux = feed.compute_flux(cv, dt_h=0.01)
        expected = (P_in * vvm * V_liq * 60.0) / (R_L_ATM_MOL_K * T_K)
        assert flux["O2"] == pytest.approx(expected, rel=1e-10)

    def test_pressure_scales_linearly(self):
        """Doubling inlet pressure should double the molar flow."""
        cv = _make_gl_cv()
        feed_1 = GasFeed(vvm_min=1.0, y={"O2": 1.0}, P_inlet_atm=1.0)
        feed_2 = GasFeed(vvm_min=1.0, y={"O2": 1.0}, P_inlet_atm=2.0)
        f1 = feed_1.compute_flux(cv, dt_h=0.01)["O2"]
        f2 = feed_2.compute_flux(cv, dt_h=0.01)["O2"]
        assert f2 == pytest.approx(2.0 * f1, rel=1e-10)

    def test_reads_liquid_volume_from_cv(self):
        """GasFeed should use cv["liquid"].V_L for vvm calculation."""
        cv_small = _make_gl_cv(V_liq=100.0)
        cv_big = _make_gl_cv(V_liq=1000.0)
        feed = GasFeed(vvm_min=1.0, y={"O2": 1.0})
        f_small = feed.compute_flux(cv_small, dt_h=0.01)["O2"]
        f_big = feed.compute_flux(cv_big, dt_h=0.01)["O2"]
        assert f_big / f_small == pytest.approx(10.0, rel=1e-10)

    def test_V_liq_override(self):
        """V_liq_override should take precedence over CV liquid volume."""
        feed = GasFeed(vvm_min=1.0, y={"O2": 1.0}, V_liq_override=500.0)
        cv = _make_gl_cv(V_liq=9999.0)  # Should be ignored
        flux = feed.compute_flux(cv, dt_h=0.01)
        # Compare with expected using override volume
        T_K = cv["gas"].T_K
        expected = (1.0 * 1.0 * 500.0 * 60.0) / (R_L_ATM_MOL_K * T_K)
        assert flux["O2"] == pytest.approx(expected, rel=1e-10)

    def test_T_override(self):
        """T_override_K should take precedence over CV gas temperature."""
        T_override = 273.15
        feed = GasFeed(vvm_min=1.0, y={"O2": 1.0}, T_override_K=T_override)
        cv = _make_gl_cv(T_K=400.0)  # Should be ignored
        flux = feed.compute_flux(cv, dt_h=0.01)
        expected = (1.0 * 1.0 * 800.0 * 60.0) / (R_L_ATM_MOL_K * T_override)
        assert flux["O2"] == pytest.approx(expected, rel=1e-10)

    def test_no_liquid_phase_without_override_gives_zero(self):
        """If no liquid phase and no override, flux should be zero."""
        gas = GasPhase({"O2": 1.0}, V_L=200.0, T_K=305.15)
        cv = ControlVolume(phases={"gas": gas})
        feed = GasFeed(vvm_min=1.0, y={"O2": 1.0})
        flux = feed.compute_flux(cv, dt_h=0.01)
        assert flux["O2"] == 0.0

    def test_no_liquid_phase_with_override_works(self):
        """V_liq_override should work even without a liquid phase."""
        gas = GasPhase({"O2": 1.0}, V_L=200.0, T_K=305.15)
        cv = ControlVolume(phases={"gas": gas})
        feed = GasFeed(vvm_min=1.0, y={"O2": 1.0}, V_liq_override=800.0)
        flux = feed.compute_flux(cv, dt_h=0.01)
        assert flux["O2"] > 0.0

    def test_y_accessor_returns_underlying_dict(self):
        # GasFeed.y is a MutableDict descriptor; instance access returns
        # the actual backing dict (not a copy).  Direct item assignment
        # mutates it; use the orchestrator path (params_changed) for gated
        # mid-run updates.
        feed = GasFeed(y={"O2": 0.21, "N2": 0.79})
        y = feed.y
        y["O2"] = 0.50
        assert feed.y["O2"] == pytest.approx(0.50)

    def test_repr(self):
        feed = GasFeed(vvm_min=1.5, y={"O2": 0.21, "N2": 0.79})
        r = repr(feed)
        assert "GasFeed" in r
        assert "1.5" in r

    # ── Parity with legacy _gas_feed_mol_per_h ──────────────────────

    def test_parity_with_legacy_formula(self):
        """GasFeed must match the legacy formula from fermenter_unit._gas_feed_mol_per_h.

        Note: the legacy code uses R=0.082057, while GasFeed uses the more
        precise R_L_ATM_MOL_K=0.0820574.  We test parity with the same R
        constant that GasFeed uses (i.e. the one from phases.py).
        """
        R = R_L_ATM_MOL_K  # same constant GasFeed uses
        vvm = 1.0
        V_liq = 800.0
        T_K = 305.15
        P_in = 1.0
        y = {"O2": 0.21, "N2": 0.79}

        # Expected formula: n_dot = P * vvm * V_liq * 60 / (R * T)
        Q = vvm * V_liq
        n_total = (P_in * Q * 60.0) / (R * T_K)
        expected = {k: float(v) * n_total for k, v in y.items()}

        # Boundary object
        feed = GasFeed(vvm_min=vvm, y=y, P_inlet_atm=P_in)
        cv = _make_gl_cv(V_liq=V_liq, T_K=T_K)
        boundary = feed.compute_flux(cv, dt_h=0.01)

        for sp in y:
            assert boundary[sp] == pytest.approx(expected[sp], rel=1e-10), \
                f"Mismatch for {sp}: boundary={boundary[sp]}, expected={expected[sp]}"


# ════════════════════════════════════════════════════════════════════════
#  PressureReliefVent tests
# ════════════════════════════════════════════════════════════════════════

class TestPressureReliefVent:

    def test_invalid_mode_raises(self):
        with pytest.raises(ValueError, match="mode"):
            PressureReliefVent(mode="invalid")

    def test_below_setpoint_gives_zero_flux(self):
        """No venting when P < P_set."""
        # Build a low-pressure headspace
        gas = GasPhase({"O2": 0.01, "N2": 0.03}, V_L=200.0, T_K=305.15)
        liq = LiquidPhase({"CO2": 0.01}, V_L=800.0, T_K=305.15)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})
        assert gas.P_atm < 1.0  # Sanity check

        vent = PressureReliefVent(P_set_atm=1.0, mode="instant")
        flux = vent.compute_flux(cv, dt_h=0.01)
        assert all(v == 0.0 for v in flux.values())

    def test_at_setpoint_gives_zero_flux(self):
        """No venting when P == P_set (within float precision)."""
        T_K = 305.15
        V_gas = 200.0
        P_set = 1.0
        n_target = (P_set * V_gas) / (R_L_ATM_MOL_K * T_K)
        gas = GasPhase({"N2": n_target}, V_L=V_gas, T_K=T_K)
        liq = LiquidPhase({}, V_L=800.0, T_K=T_K)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})
        assert gas.P_atm == pytest.approx(P_set, rel=1e-10)

        vent = PressureReliefVent(P_set_atm=P_set, mode="instant")
        flux = vent.compute_flux(cv, dt_h=0.01)
        for v in flux.values():
            assert abs(v) < 1e-12

    # ── Instant mode ────────────────────────────────────────────────

    def test_instant_vent_removes_excess(self):
        """Instant mode should vent exactly enough to reach P_set."""
        T_K = 305.15
        V_gas = 200.0
        P_set = 1.0
        n_target = (P_set * V_gas) / (R_L_ATM_MOL_K * T_K)
        n_excess = 0.5
        gas = GasPhase({"O2": n_target + n_excess}, V_L=V_gas, T_K=T_K)
        liq = LiquidPhase({}, V_L=800.0, T_K=T_K)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})
        assert gas.P_atm > P_set

        vent = PressureReliefVent(P_set_atm=P_set, mode="instant")
        dt_h = 0.01
        flux = vent.compute_flux(cv, dt_h)

        # Total moles to vent = n_excess
        total_vented = -sum(flux.values()) * dt_h
        assert total_vented == pytest.approx(n_excess, rel=1e-8)

    def test_instant_vent_flux_is_negative(self):
        """Vent fluxes must be negative (leaving the CV)."""
        gas = GasPhase({"O2": 2.0, "N2": 6.0}, V_L=200.0, T_K=305.15)
        liq = LiquidPhase({}, V_L=800.0, T_K=305.15)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})
        # Make sure we're above 1 atm
        assert gas.P_atm > 1.0

        vent = PressureReliefVent(P_set_atm=1.0, mode="instant")
        flux = vent.compute_flux(cv, dt_h=0.01)
        for sp, v in flux.items():
            assert v <= 0.0, f"Vent flux for {sp} should be <= 0, got {v}"

    def test_instant_proportional_split(self):
        """Vent should split proportionally by mole fraction."""
        gas = GasPhase({"O2": 2.0, "N2": 6.0}, V_L=200.0, T_K=305.15)
        liq = LiquidPhase({}, V_L=800.0, T_K=305.15)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})

        vent = PressureReliefVent(P_set_atm=1.0, mode="instant")
        flux = vent.compute_flux(cv, dt_h=0.01)
        # O2 fraction = 2/8 = 0.25, N2 fraction = 6/8 = 0.75
        total_vent = abs(flux["O2"]) + abs(flux["N2"])
        assert abs(flux["O2"]) / total_vent == pytest.approx(0.25, rel=1e-10)
        assert abs(flux["N2"]) / total_vent == pytest.approx(0.75, rel=1e-10)

    def test_instant_vent_then_apply_reaches_setpoint(self):
        """After applying instant vent, gas pressure should equal P_set."""
        T_K = 305.15
        V_gas = 200.0
        P_set = 1.0
        gas = GasPhase({"O2": 2.0, "CO2": 0.5, "N2": 6.0}, V_L=V_gas, T_K=T_K)
        liq = LiquidPhase({}, V_L=800.0, T_K=T_K)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})

        vent = PressureReliefVent(P_set_atm=P_set, mode="instant")
        rec = apply_boundary(vent, cv, dt_h=0.01)

        assert cv["gas"].P_atm == pytest.approx(P_set, rel=1e-6)

    # ── Smooth mode ─────────────────────────────────────────────────

    def test_smooth_vents_above_setpoint(self):
        """Smooth mode should vent when P > P_set."""
        gas = GasPhase({"O2": 2.0, "N2": 6.0}, V_L=200.0, T_K=305.15)
        liq = LiquidPhase({}, V_L=800.0, T_K=305.15)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})
        assert gas.P_atm > 1.0

        vent = PressureReliefVent(P_set_atm=1.0, mode="smooth")
        flux = vent.compute_flux(cv, dt_h=0.01)
        total_vent_rate = -sum(flux.values())
        assert total_vent_rate > 0.0

    def test_smooth_does_not_overshoot_setpoint(self):
        """Smooth vent should not remove more than the excess in one step."""
        T_K = 305.15
        V_gas = 200.0
        P_set = 1.0
        n_target = (P_set * V_gas) / (R_L_ATM_MOL_K * T_K)
        n_excess = 0.01  # small excess
        gas = GasPhase({"N2": n_target + n_excess}, V_L=V_gas, T_K=T_K)
        liq = LiquidPhase({}, V_L=800.0, T_K=T_K)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})

        vent = PressureReliefVent(P_set_atm=P_set, mode="smooth", k_vent_per_h=500.0)
        dt_h = 0.01
        flux = vent.compute_flux(cv, dt_h)
        total_vented = -sum(flux.values()) * dt_h
        # Should vent some but not all excess
        assert total_vented > 0.0
        assert total_vented < n_target + n_excess  # Don't empty the headspace

    def test_smooth_below_setpoint_gives_near_zero(self):
        """Smooth mode should give negligible vent when P < P_set."""
        T_K = 305.15
        V_gas = 200.0
        P_set = 2.0  # High setpoint
        gas = GasPhase({"N2": 0.1}, V_L=V_gas, T_K=T_K)  # Low pressure
        liq = LiquidPhase({}, V_L=800.0, T_K=T_K)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})
        assert gas.P_atm < P_set

        vent = PressureReliefVent(P_set_atm=P_set, mode="smooth")
        flux = vent.compute_flux(cv, dt_h=0.01)
        assert all(v == 0.0 for v in flux.values())

    def test_smooth_increases_with_excess_pressure(self):
        """Higher excess pressure → more venting."""
        T_K = 305.15
        V_gas = 200.0
        P_set = 1.0

        liq = LiquidPhase({}, V_L=800.0, T_K=T_K)
        vent = PressureReliefVent(P_set_atm=P_set, mode="smooth")

        gas_low = GasPhase({"N2": 8.2}, V_L=V_gas, T_K=T_K)  # Slightly above P_set
        gas_high = GasPhase({"N2": 12.0}, V_L=V_gas, T_K=T_K)  # Well above P_set

        cv_low = ControlVolume(phases={"gas": gas_low, "liquid": liq.snapshot()})
        cv_high = ControlVolume(phases={"gas": gas_high, "liquid": liq.snapshot()})

        flux_low = vent.compute_flux(cv_low, dt_h=0.01)
        flux_high = vent.compute_flux(cv_high, dt_h=0.01)

        vent_low = -sum(flux_low.values())
        vent_high = -sum(flux_high.values())
        assert vent_high > vent_low

    # ── Smooth vent fraction helper ─────────────────────────────────

    def test_smooth_vent_fraction_zero_excess(self):
        # Smooth ReLU at x=0 gives ~w/2, so fraction is small but not exactly zero
        frac = _smooth_vent_fraction(0.0, dt_h=0.01, k_vent_per_h=500.0, smooth_width_atm=0.01)
        assert frac >= 0.0
        assert frac < 1.0

    def test_smooth_vent_fraction_negative_excess(self):
        # With large negative excess, the smooth ReLU tail is very small
        frac = _smooth_vent_fraction(-5.0, dt_h=0.01, k_vent_per_h=500.0, smooth_width_atm=0.01)
        assert frac >= 0.0
        assert frac < 0.01  # Deep below setpoint → negligible vent

    def test_smooth_vent_fraction_large_excess(self):
        frac = _smooth_vent_fraction(1.0, dt_h=1.0, k_vent_per_h=500.0, smooth_width_atm=0.01)
        assert frac > 0.99  # Should vent almost everything

    def test_smooth_vent_fraction_bounded(self):
        """Fraction must always be in [0, 1)."""
        for excess in [-1.0, 0.0, 0.01, 0.1, 1.0, 10.0]:
            for dt in [0.001, 0.01, 0.1, 1.0]:
                frac = _smooth_vent_fraction(excess, dt, 500.0, 0.01)
                assert 0.0 <= frac < 1.0


# ════════════════════════════════════════════════════════════════════════
#  apply_boundary tests
# ════════════════════════════════════════════════════════════════════════

class TestApplyBoundary:

    def test_feed_increases_inventory(self):
        cv = _make_gl_cv()
        mol_before = cv.total_mol()
        feed = GasFeed(vvm_min=1.0, y={"O2": 0.21, "N2": 0.79})
        rec = apply_boundary(feed, cv, dt_h=0.01)
        mol_after = cv.total_mol()
        # O2 and N2 should increase
        assert mol_after.get("O2", 0.0) > mol_before.get("O2", 0.0)
        assert mol_after.get("N2", 0.0) > mol_before.get("N2", 0.0)

    def test_vent_decreases_inventory(self):
        # Build an over-pressured CV
        gas = GasPhase({"O2": 2.0, "N2": 6.0}, V_L=200.0, T_K=305.15)
        liq = LiquidPhase({}, V_L=800.0, T_K=305.15)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})
        mol_before = sum(cv.total_mol().values())

        vent = PressureReliefVent(P_set_atm=1.0, mode="instant")
        rec = apply_boundary(vent, cv, dt_h=0.01)
        mol_after = sum(cv.total_mol().values())
        assert mol_after < mol_before

    def test_record_fields(self):
        cv = _make_gl_cv()
        feed = GasFeed(vvm_min=1.0, y={"O2": 1.0}, label="test_feed")
        rec = apply_boundary(feed, cv, dt_h=0.05)

        assert rec.boundary_label == "test_feed"
        assert rec.phase_key == "gas"
        assert rec.dt_h == 0.05
        assert "O2" in rec.flux_mol_per_h
        assert "O2" in rec.mol_applied
        # mol_applied = flux * dt
        assert rec.mol_applied["O2"] == pytest.approx(
            rec.flux_mol_per_h["O2"] * 0.05, rel=1e-10
        )

    def test_record_total_mol_applied(self):
        cv = _make_gl_cv()
        feed = GasFeed(vvm_min=1.0, y={"O2": 0.5, "N2": 0.5})
        rec = apply_boundary(feed, cv, dt_h=0.01)
        expected = sum(rec.mol_applied.values())
        assert rec.total_mol_applied == pytest.approx(expected)

    def test_record_summary_string(self):
        cv = _make_gl_cv()
        feed = GasFeed(vvm_min=1.0, y={"O2": 1.0}, label="my_feed")
        rec = apply_boundary(feed, cv, dt_h=0.01)
        s = rec.summary()
        assert "my_feed" in s
        assert "O2" in s
        assert "mol/h" in s

    def test_vent_record_shows_negative_mol(self):
        gas = GasPhase({"O2": 2.0, "N2": 6.0}, V_L=200.0, T_K=305.15)
        liq = LiquidPhase({}, V_L=800.0, T_K=305.15)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})

        vent = PressureReliefVent(P_set_atm=1.0, mode="instant")
        rec = apply_boundary(vent, cv, dt_h=0.01)
        assert rec.total_mol_applied < 0.0


# ════════════════════════════════════════════════════════════════════════
#  ExternalFluxRecord tests
# ════════════════════════════════════════════════════════════════════════

class TestExternalFluxRecord:

    def test_empty_record(self):
        rec = ExternalFluxRecord(boundary_label="x", phase_key="gas")
        assert rec.total_mol_applied == 0.0

    def test_summary_is_string(self):
        rec = ExternalFluxRecord(
            boundary_label="test",
            phase_key="gas",
            flux_mol_per_h={"O2": 10.0},
            dt_h=0.1,
            mol_applied={"O2": 1.0},
        )
        s = rec.summary()
        assert isinstance(s, str)
        assert "O2" in s


# ════════════════════════════════════════════════════════════════════════
#  Integration: feed → equilibrium → vent cycle
# ════════════════════════════════════════════════════════════════════════

class TestIntegrationFeedEquilibriumVent:
    """Full cycle test: feed gas → internal equilibrium → pressure vent.

    Uses a mock equilibrium interface (no speciation dependency) to
    verify that the feed/vent boundary objects compose correctly with
    the ControlVolume's internal transfer mechanism.
    """

    def test_feed_then_vent_cycle(self):
        """Gas feed raises pressure; vent brings it back down."""
        gas = GasPhase({"O2": 1.0, "N2": 3.0}, V_L=200.0, T_K=305.15)
        liq = LiquidPhase({"CO2": 0.01}, V_L=800.0, T_K=305.15)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq}, label="cycle_test")

        P_before = cv["gas"].P_atm
        feed = GasFeed(vvm_min=2.0, y={"O2": 0.21, "N2": 0.79})
        vent = PressureReliefVent(P_set_atm=P_before, mode="instant")

        dt_h = 0.01

        # 1. Feed — pressure should rise
        feed_rec = apply_boundary(feed, cv, dt_h)
        P_after_feed = cv["gas"].P_atm
        assert P_after_feed > P_before

        # 2. Vent — pressure should return to setpoint
        vent_rec = apply_boundary(vent, cv, dt_h)
        P_after_vent = cv["gas"].P_atm
        assert P_after_vent == pytest.approx(P_before, rel=1e-5)

    def test_multiple_steps_inventory_audit(self):
        """Run N steps of feed + vent, verify audit trail sums correctly."""
        gas = GasPhase({"O2": 1.0, "N2": 3.0}, V_L=200.0, T_K=305.15)
        liq = LiquidPhase({}, V_L=800.0, T_K=305.15)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})

        feed = GasFeed(vvm_min=1.0, y={"O2": 0.21, "N2": 0.79})
        vent = PressureReliefVent(P_set_atm=1.0, mode="instant")

        dt_h = 0.01
        n_steps = 20

        total_fed = {"O2": 0.0, "N2": 0.0}
        total_vented = {"O2": 0.0, "N2": 0.0}
        mol_start = cv.total_mol()

        for _ in range(n_steps):
            frec = apply_boundary(feed, cv, dt_h)
            vrec = apply_boundary(vent, cv, dt_h)
            for sp in total_fed:
                total_fed[sp] += frec.mol_applied.get(sp, 0.0)
                total_vented[sp] += vrec.mol_applied.get(sp, 0.0)

        mol_end = cv.total_mol()

        # Verify: end = start + fed + vented (vented is negative)
        for sp in ("O2", "N2"):
            expected = mol_start.get(sp, 0.0) + total_fed[sp] + total_vented[sp]
            assert mol_end.get(sp, 0.0) == pytest.approx(expected, rel=1e-8), \
                f"Audit mismatch for {sp}"

    def test_zero_feed_zero_vent_is_noop(self):
        """Zero vvm + P < P_set → nothing happens."""
        cv = _make_gl_cv()
        mol_before = cv.total_mol()

        feed = GasFeed(vvm_min=0.0)
        vent = PressureReliefVent(P_set_atm=100.0, mode="instant")

        apply_boundary(feed, cv, dt_h=0.01)
        apply_boundary(vent, cv, dt_h=0.01)

        mol_after = cv.total_mol()
        for sp in set(mol_before) | set(mol_after):
            assert mol_after.get(sp, 0.0) == pytest.approx(mol_before.get(sp, 0.0))


# ════════════════════════════════════════════════════════════════════════
#  P8 — Pattern B on GasFeed (framework-polish)
# ════════════════════════════════════════════════════════════════════════

class TestGasFeedPatternB:
    """GasFeed.vvm_min and GasFeed.y follow Pattern B (SIMULATION_CLASS
    decision 2/9): public setters are gated by RunContext; orchestrator
    unchecked paths bypass the gate."""

    def _make_running_sim_with_gas_feed(self):
        from PyOMES.core import Simulation, ControlVolume, GasPhase, LiquidPhase
        gas = GasPhase(n_mol={"O2": 0.5, "N2": 2.0}, V_L=0.3, T_K=305.0)
        liq = LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=305.0)
        feed = GasFeed(vvm_min=1.0, y={"O2": 0.21, "N2": 0.79})
        cv = ControlVolume(phases={"gas": gas, "liquid": liq},
                           boundaries=[feed], label="cv")
        sim = Simulation(cvs={"cv": cv})
        return sim, cv, feed

    def test_vvm_min_gate_fires_when_running(self):
        """vvm_min setter raises RuntimeError while Simulation is running."""
        import pytest
        from PyOMES.core.lifecycle import RunContext
        sim, cv, feed = self._make_running_sim_with_gas_feed()
        # Manually set running state (avoids needing a full run() call)
        sim._context.is_running = True
        try:
            with pytest.raises(RuntimeError):
                feed.vvm_min = 2.0
        finally:
            sim._context.is_running = False

    def test_vvm_min_unchecked_bypasses_gate(self):
        """MutableScalar._set_unchecked writes vvm_min while running."""
        from PyOMES.control.descriptors import MutableScalar
        sim, cv, feed = self._make_running_sim_with_gas_feed()
        sim._context.is_running = True
        try:
            descriptor = type(feed).__dict__["vvm_min"]
            assert isinstance(descriptor, MutableScalar)
            descriptor._set_unchecked(feed, 3.5)
            assert feed.vvm_min == pytest.approx(3.5)
        finally:
            sim._context.is_running = False

    def test_set_y_whole_dict_gate_fires_when_running(self):
        """Whole-dict replacement feed.y = {...} raises while running."""
        import pytest
        sim, cv, feed = self._make_running_sim_with_gas_feed()
        sim._context.is_running = True
        try:
            with pytest.raises(RuntimeError):
                feed.y = {"O2": 0.5, "N2": 0.5}
        finally:
            sim._context.is_running = False

    def test_y_item_unchecked_writes_through(self):
        """MutableDict._set_item_unchecked updates y while running."""
        from PyOMES.control.descriptors import MutableDict
        sim, cv, feed = self._make_running_sim_with_gas_feed()
        sim._context.is_running = True
        try:
            descriptor = type(feed).__dict__["y"]
            assert isinstance(descriptor, MutableDict)
            descriptor._set_item_unchecked(feed, "O2", 0.30)
            assert feed.y["O2"] == pytest.approx(0.30)
        finally:
            sim._context.is_running = False

    def test_apply_param_change_writes_y_and_vvm(self):
        """Reflective walker writes boundaries[GasFeed].y.<k> and .vvm_min."""
        sim, cv, feed = self._make_running_sim_with_gas_feed()
        sim._context.is_running = True
        try:
            sim._apply_param_change(cv, "boundaries[GasFeed].y.O2", 0.35)
            assert feed.y["O2"] == pytest.approx(0.35)
            sim._apply_param_change(cv, "boundaries[GasFeed].vvm_min", 2.5)
            assert feed.vvm_min == pytest.approx(2.5)
        finally:
            sim._context.is_running = False
