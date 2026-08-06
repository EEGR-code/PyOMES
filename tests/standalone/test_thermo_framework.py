# -*- coding: utf-8 -*-
"""Tests for ThermoFramework."""
from __future__ import annotations

import dataclasses
import math
import pytest


class TestThermoFrameworkConstruction:
    def test_default_liquid_activity_is_ideal(self):
        from PyOMES.thermo import ThermoFramework, IdealLiquidModel
        tf = ThermoFramework()
        assert isinstance(tf.liquid_activity, IdealLiquidModel)

    def test_default_use_activity_false(self):
        from PyOMES.thermo import ThermoFramework
        tf = ThermoFramework()
        assert tf.use_activity is False

    def test_default_activity_model_name(self):
        from PyOMES.thermo import ThermoFramework
        tf = ThermoFramework()
        assert tf.activity_model == "ideal"

    def test_default_standard_conditions(self):
        from PyOMES.thermo import ThermoFramework
        tf = ThermoFramework()
        assert tf.standard_T_K == pytest.approx(298.15)
        assert tf.standard_P_atm == pytest.approx(1.0)

    def test_davies_liquid_activity(self):
        from PyOMES.thermo import ThermoFramework, DaviesLiquidModel
        tf = ThermoFramework(liquid_activity=DaviesLiquidModel())
        assert isinstance(tf.liquid_activity, DaviesLiquidModel)
        assert tf.use_activity is True
        assert tf.activity_model == "davies"

    def test_sit_liquid_activity(self):
        from PyOMES.thermo import ThermoFramework
        from PyOMES.thermo import SITLiquidModel
        tf = ThermoFramework(liquid_activity=SITLiquidModel())
        assert tf.use_activity is True
        assert tf.activity_model == "sit"

    def test_frozen(self):
        from PyOMES.thermo import ThermoFramework
        tf = ThermoFramework()
        with pytest.raises((dataclasses.FrozenInstanceError, AttributeError)):
            tf.standard_T_K = 310.0  # type: ignore

    def test_replace(self):
        from PyOMES.thermo import ThermoFramework, DaviesLiquidModel
        tf = ThermoFramework()
        tf2 = dataclasses.replace(tf, liquid_activity=DaviesLiquidModel())
        assert tf2.use_activity is True
        assert tf.use_activity is False  # original unchanged

    def test_gas_eos_none_by_default(self):
        from PyOMES.thermo import ThermoFramework
        assert ThermoFramework().gas_eos is None


class TestThermoFrameworkPKaAtT:
    def test_no_correction_returns_ref(self):
        from PyOMES.thermo import ThermoFramework
        tf = ThermoFramework()
        pKa = tf.pKa_at_T(pKa_ref=6.35, dH_J_per_mol=0.0, T_K=308.15)
        assert pKa == pytest.approx(6.35, abs=1e-12)

    def test_same_T_returns_ref(self):
        from PyOMES.thermo import ThermoFramework
        tf = ThermoFramework()
        pKa = tf.pKa_at_T(pKa_ref=6.35, dH_J_per_mol=7646.0,
                           T_K=298.15, T_ref_K=298.15)
        assert pKa == pytest.approx(6.35, abs=1e-10)

    def test_van_t_hoff_correction(self):
        """CO₂ pKa₁ shifts from 6.35 at 25°C to ~6.24 at 35°C (BSM2 values)."""
        from PyOMES.thermo import ThermoFramework
        tf = ThermoFramework()
        pKa_35 = tf.pKa_at_T(pKa_ref=6.35, dH_J_per_mol=7646.0,
                              T_K=308.15, T_ref_K=298.15)
        assert 6.20 < pKa_35 < 6.35

    def test_negative_dH_shifts_up(self):
        from PyOMES.thermo import ThermoFramework
        tf = ThermoFramework()
        pKa_hot = tf.pKa_at_T(pKa_ref=7.0, dH_J_per_mol=-5000.0,
                               T_K=318.15, T_ref_K=298.15)
        assert pKa_hot > 7.0


class TestThermoFrameworkKwAtT:
    def test_at_ref_T(self):
        from PyOMES.thermo import ThermoFramework
        tf = ThermoFramework()
        Kw = tf.Kw_at_T(T_K=298.15)
        assert Kw == pytest.approx(1e-14, rel=1e-6)

    def test_higher_T_increases_Kw(self):
        from PyOMES.thermo import ThermoFramework
        tf = ThermoFramework()
        Kw_35 = tf.Kw_at_T(T_K=308.15)
        assert Kw_35 > 1e-14


class TestChemicalEquilibriumEngineThermo:
    def test_thermo_overrides_explicit_kwargs(self):
        from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
        from PyOMES.thermo import ThermoFramework, DaviesLiquidModel
        tf = ThermoFramework(liquid_activity=DaviesLiquidModel())
        eng = BisectionChemicalEquilibriumEngine(
            use_activity=False,       # overridden by thermo
            activity_model="ideal",   # overridden by thermo
            thermo=tf,
        )
        assert eng.use_activity is True
        assert eng.activity_model == "davies"

    def test_no_thermo_uses_explicit(self):
        from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
        eng = BisectionChemicalEquilibriumEngine(use_activity=True, activity_model="ideal")
        assert eng.use_activity is True
        assert eng.activity_model == "ideal"

    def test_thermo_liquid_activity_stored(self):
        from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
        from PyOMES.thermo import ThermoFramework, DaviesLiquidModel
        model = DaviesLiquidModel()
        tf = ThermoFramework(liquid_activity=model)
        eng = BisectionChemicalEquilibriumEngine(thermo=tf)
        assert eng._liquid_activity is model


class TestThermoPresets:
    def test_thermo_ideal(self):
        from PyOMES.thermo.framework import THERMO_IDEAL
        from PyOMES.thermo import IdealLiquidModel
        assert isinstance(THERMO_IDEAL.liquid_activity, IdealLiquidModel)
        assert THERMO_IDEAL.use_activity is False
        assert THERMO_IDEAL.activity_model == "ideal"

    def test_thermo_davies(self):
        from PyOMES.thermo.framework import THERMO_DAVIES
        from PyOMES.thermo import DaviesLiquidModel
        assert isinstance(THERMO_DAVIES.liquid_activity, DaviesLiquidModel)
        assert THERMO_DAVIES.use_activity is True
        assert THERMO_DAVIES.activity_model == "davies"
