"""Tests for fermenter.kinetics — model interface and parameter handling."""

import pytest
import numpy as np
from PyOMES.kinetics.core import ParameterSpec, ParameterSet, Environment
from PyOMES.kinetics.models import YeastAcetateV1
from PyOMES.kinetics.mapping import StateMapping


@pytest.mark.unit
class TestParameterSet:
    """Test ParameterSet construction and validation."""

    def test_from_schema_defaults(self):
        schema = (
            ParameterSpec("mu_max", default=0.5),
            ParameterSpec("Ks", default=0.01),
        )
        ps = ParameterSet.from_schema(schema)
        assert ps["mu_max"] == 0.5
        assert ps["Ks"] == 0.01

    def test_overrides(self):
        schema = (ParameterSpec("mu_max", default=0.5),)
        ps = ParameterSet.from_schema(schema, overrides={"mu_max": 0.8})
        assert ps["mu_max"] == 0.8

    def test_unknown_param_raises(self):
        schema = (ParameterSpec("mu_max", default=0.5),)
        with pytest.raises(KeyError, match="Unknown parameter"):
            ParameterSet.from_schema(schema, overrides={"bad_param": 1.0})

    def test_bounds_violation_raises(self):
        schema = (ParameterSpec("mu_max", default=0.5, lower=0.0),)
        with pytest.raises(ValueError, match="lower bound"):
            ParameterSet.from_schema(schema, overrides={"mu_max": -1.0})


@pytest.mark.unit
class TestEnvironment:
    """Test Environment dataclass."""

    def test_construction(self):
        env = Environment(T_K=305.15, V_L=1.0, pH=7.0, I_molal=0.0,
                         DO_mol_L=1e-4, CO2aq_mol_L=1e-3, P_atm=1.0)
        assert env.T_K == 305.15
        assert env.DO_mol_L == 1e-4


@pytest.mark.unit
class TestYeastAcetateV1:
    """Test the plug-in kinetic model."""

    def test_state_ids(self):
        model = YeastAcetateV1()
        assert len(model.state_ids) > 0
        assert "X_gDW" in model.state_ids

    def test_rhs_shape(self):
        model = YeastAcetateV1()
        params = ParameterSet.from_schema(model.param_schema)
        env = Environment(T_K=305.15, V_L=1.0, pH=6.0, I_molal=0.01,
                         DO_mol_L=1e-4, CO2aq_mol_L=1e-3, P_atm=1.0)
        y = np.array([0.1, 0.01, 1e-4, 1e-3])
        dydt, outputs = model.rhs(t_h=0.0, y=y, env=env, params=params)
        assert dydt.shape == y.shape
        assert isinstance(outputs, dict)

    def test_growth_positive(self):
        """With substrate and O2 available, biomass growth should be positive."""
        model = YeastAcetateV1()
        params = ParameterSet.from_schema(model.param_schema)
        env = Environment(T_K=305.15, V_L=1.0, pH=6.0, I_molal=0.01,
                         DO_mol_L=1e-4, CO2aq_mol_L=1e-3, P_atm=1.0)
        y = np.array([0.1, 0.01, 1e-4, 1e-3])
        dydt, _ = model.rhs(t_h=0.0, y=y, env=env, params=params)
        assert dydt[0] > 0, "dX/dt should be positive with available substrate"

    def test_no_growth_without_substrate(self):
        """With no substrate, growth should be negligible."""
        model = YeastAcetateV1()
        params = ParameterSet.from_schema(model.param_schema)
        env = Environment(T_K=305.15, V_L=1.0, pH=6.0, I_molal=0.01,
                         DO_mol_L=1e-4, CO2aq_mol_L=1e-3, P_atm=1.0)
        y = np.array([0.1, 0.0, 1e-4, 1e-3])  # zero acetate
        dydt, _ = model.rhs(t_h=0.0, y=y, env=env, params=params)
        assert dydt[0] == pytest.approx(0.0, abs=1e-15)


@pytest.mark.unit
class TestStateMapping:
    """Test StateMapping between canonical and model state IDs."""

    def test_auto_mapping(self):
        mapping = StateMapping.auto(("X_g_L", "SAc_g_L"))
        assert mapping.state_ids == ("X_g_L", "SAc_g_L")

    def test_explicit_mapping(self):
        mapping = StateMapping(
            state_ids=("X_gDW", "Acetate_mol"),
            canonical_to_state_id={"X_g_L": "X_gDW", "SAc_g_L": "Acetate_mol"},
        )
        assert mapping.canonical_to_state_id["X_g_L"] == "X_gDW"
