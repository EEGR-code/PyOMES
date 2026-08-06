"""Tests for fermenter.properties.viscosity — viscosity model framework."""

import pytest
import math
import numpy as np
from PyOMES.properties.viscosity import (
    BrothState,
    ViscosityResult,
    ViscosityModel,
    wrap_viscosity_model,
    WaterViscosity,
    ArrheniusBiomassViscosity,
    PowerLawBiomassViscosity,
    mu_water_Pa_s,
    _CallableViscosityAdapter,
)


class TestBrothState:
    def test_defaults(self):
        bs = BrothState()
        assert bs.T_K == 305.15
        assert bs.X_g_L == 0.0
        assert bs.pH is None

    def test_custom(self):
        bs = BrothState(T_K=310.0, X_g_L=50.0, pH=6.5)
        assert bs.T_K == 310.0
        assert bs.X_g_L == 50.0
        assert bs.pH == 6.5

    def test_frozen(self):
        bs = BrothState()
        with pytest.raises(AttributeError):
            bs.T_K = 999.0


class TestViscosityResult:
    def test_basic(self):
        vr = ViscosityResult(mu_Pa_s=1e-3)
        assert vr.mu_Pa_s == 1e-3
        assert vr.mu_rel is None
        assert vr.extras == {}


class TestMuWater:
    def test_25C(self):
        mu = mu_water_Pa_s(298.15)
        assert mu == pytest.approx(8.9e-4, rel=0.05)

    def test_decreases_with_temperature(self):
        mu_20 = mu_water_Pa_s(293.15)
        mu_50 = mu_water_Pa_s(323.15)
        assert mu_20 > mu_50

    def test_positive(self):
        for T in [273.15, 298.15, 310.15, 373.15]:
            assert mu_water_Pa_s(T) > 0


class TestWaterViscosity:
    def test_returns_result(self):
        model = WaterViscosity()
        r = model.compute(BrothState(T_K=298.15))
        assert isinstance(r, ViscosityResult)
        assert r.mu_Pa_s > 0
        assert r.mu_rel == 1.0

    def test_satisfies_protocol(self):
        assert isinstance(WaterViscosity(), ViscosityModel)


class TestArrheniusBiomassViscosity:
    def test_zero_biomass_equals_water(self):
        model = ArrheniusBiomassViscosity(k_X=0.005)
        r = model.compute(BrothState(T_K=298.15, X_g_L=0.0))
        mu_w = mu_water_Pa_s(298.15)
        assert r.mu_Pa_s == pytest.approx(mu_w, rel=1e-6)
        assert r.mu_rel == pytest.approx(1.0, rel=1e-6)

    def test_increases_with_biomass(self):
        model = ArrheniusBiomassViscosity(k_X=0.005)
        r_low = model.compute(BrothState(T_K=305.15, X_g_L=10.0))
        r_high = model.compute(BrothState(T_K=305.15, X_g_L=100.0))
        assert r_high.mu_Pa_s > r_low.mu_Pa_s

    def test_known_correction(self):
        model = ArrheniusBiomassViscosity(k_X=0.01)
        r = model.compute(BrothState(T_K=305.15, X_g_L=50.0))
        expected_rel = math.exp(0.01 * 50.0)
        assert r.mu_rel == pytest.approx(expected_rel, rel=1e-6)

    def test_extras_populated(self):
        model = ArrheniusBiomassViscosity(k_X=0.005)
        r = model.compute(BrothState(T_K=305.15, X_g_L=20.0))
        assert "mu_water_Pa_s" in r.extras
        assert "k_X" in r.extras
        assert r.extras["k_X"] == 0.005

    def test_satisfies_protocol(self):
        assert isinstance(ArrheniusBiomassViscosity(), ViscosityModel)


class TestPowerLawBiomassViscosity:
    def test_zero_biomass_equals_water(self):
        model = PowerLawBiomassViscosity(k=0.001, n=1.5)
        r = model.compute(BrothState(T_K=298.15, X_g_L=0.0))
        mu_w = mu_water_Pa_s(298.15)
        assert r.mu_Pa_s == pytest.approx(mu_w, rel=1e-6)

    def test_increases_with_biomass(self):
        model = PowerLawBiomassViscosity(k=0.001, n=1.5)
        r_low = model.compute(BrothState(T_K=305.15, X_g_L=10.0))
        r_high = model.compute(BrothState(T_K=305.15, X_g_L=100.0))
        assert r_high.mu_Pa_s > r_low.mu_Pa_s

    def test_satisfies_protocol(self):
        assert isinstance(PowerLawBiomassViscosity(), ViscosityModel)


class TestWrapViscosityModel:
    def test_none_returns_none(self):
        assert wrap_viscosity_model(None) is None

    def test_protocol_object_returned_as_is(self):
        model = WaterViscosity()
        assert wrap_viscosity_model(model) is model

    def test_callable_wrapped(self):
        fn = lambda state: 1.5e-3
        wrapped = wrap_viscosity_model(fn)
        assert isinstance(wrapped, _CallableViscosityAdapter)
        r = wrapped.compute(BrothState())
        assert r.mu_Pa_s == pytest.approx(1.5e-3)

    def test_invalid_type_raises(self):
        with pytest.raises(TypeError, match="viscosity_model"):
            wrap_viscosity_model(42)


class TestCallableAdapter:
    def test_receives_broth_state(self):
        """The callable should receive a BrothState and can read its fields."""
        captured = {}
        def fn(state):
            captured["T_K"] = state.T_K
            captured["X"] = state.X_g_L
            return 1e-3
        adapter = _CallableViscosityAdapter(fn)
        adapter.compute(BrothState(T_K=310.0, X_g_L=25.0))
        assert captured["T_K"] == 310.0
        assert captured["X"] == 25.0

    def test_ml_model_pattern(self):
        """Simulate wrapping a scikit-learn-style predict method."""
        class FakeMLModel:
            def predict(self, features):
                return [features[0] * 1e-3]  # trivial prediction

        ml = FakeMLModel()
        fn = lambda state: ml.predict([state.X_g_L])[0]
        adapter = wrap_viscosity_model(fn)
        r = adapter.compute(BrothState(X_g_L=1.0))
        assert r.mu_Pa_s == pytest.approx(1e-3)
