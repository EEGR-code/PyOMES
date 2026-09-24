# -*- coding: utf-8 -*-
"""Tests for PyOMES.thermo.gas_eos: IdealGasEOS and PengRobinsonEOS.

Pins current behaviour after moving ``PyOMES/equilibria/`` (``vle.py`` +
``peng_robinson.py``) into this one file (checkpoint 11, decision D4),
which previously had no test coverage at all:

- Compressibility factors (Z) at biogas-relevant pressures, cross-checked
  against the plausibility figures from the checkpoint-1 audit (a
  plausibility check against approximate reference values, not a
  validation).
- ``PengRobinsonEOS.partial_pressures_atm`` returns fugacities
  (``f_i = y_i * phi_i * P``); ``IdealGasEOS``'s returns plain partial
  pressures (``y_i * P``). The two differ at non-ideal (high) pressure
  and converge at low pressure.
"""
from __future__ import annotations

import pytest


class TestPengRobinsonCompressibilityFactors:
    """(n, V, T) chosen so the mixture sits at the stated total pressure."""

    T_K = 308.15

    def _eos(self):
        from PyOMES.thermo import PengRobinsonEOS, BIOGAS_SPECIES
        return PengRobinsonEOS(BIOGAS_SPECIES)

    def test_co2_20atm(self):
        eos = self._eos()
        n_gas = {"CO2": 0.884102110278784}
        V_L = 1.0
        P = eos.pressure_mixture_atm(n_gas, T_K=self.T_K, V_L=V_L)
        assert P == pytest.approx(20.0, rel=1e-9)
        Z, _ = eos.compressibility_factors(n_gas, T_K=self.T_K, V_L=V_L)
        assert Z == pytest.approx(0.894638994442354, rel=1e-9)
        assert Z == pytest.approx(0.897, abs=0.01)

    def test_n2_50atm(self):
        eos = self._eos()
        n_gas = {"N2": 1.997687212732119}
        V_L = 1.0
        P = eos.pressure_mixture_atm(n_gas, T_K=self.T_K, V_L=V_L)
        assert P == pytest.approx(50.0, rel=1e-9)
        Z, _ = eos.compressibility_factors(n_gas, T_K=self.T_K, V_L=V_L)
        assert Z == pytest.approx(0.9898349174524119, rel=1e-9)
        assert Z == pytest.approx(0.987, abs=0.01)

    def test_ch4_50atm(self):
        eos = self._eos()
        n_gas = {"CH4": 2.170219515757866}
        V_L = 1.0
        P = eos.pressure_mixture_atm(n_gas, T_K=self.T_K, V_L=V_L)
        assert P == pytest.approx(50.0, rel=1e-9)
        Z, _ = eos.compressibility_factors(n_gas, T_K=self.T_K, V_L=V_L)
        assert Z == pytest.approx(0.9111431092351556, rel=1e-9)
        assert Z == pytest.approx(0.910, abs=0.01)

    def test_1atm_biogas_mixture_is_near_ideal(self):
        eos = self._eos()
        n_gas = {"CH4": 0.6, "CO2": 0.4}
        V_L = 25.213137354330655
        P = eos.pressure_mixture_atm(n_gas, T_K=self.T_K, V_L=V_L)
        assert P == pytest.approx(1.0, rel=1e-6)
        Z, _ = eos.compressibility_factors(n_gas, T_K=self.T_K, V_L=V_L)
        assert Z == pytest.approx(0.9971193518650184, rel=1e-9)
        assert Z == pytest.approx(1.0, abs=0.01)


class TestPartialPressuresVsFugacities:
    """IdealGasEOS returns partial pressures; PengRobinsonEOS returns
    fugacities. They differ at high pressure and converge at low pressure."""

    T_K = 308.15

    def test_differ_at_high_pressure(self):
        from PyOMES.thermo import PengRobinsonEOS, IdealGasEOS, BIOGAS_SPECIES
        pr = PengRobinsonEOS(BIOGAS_SPECIES)
        ideal = IdealGasEOS()
        n_gas = {"CH4": 2.170219515757866}  # ~50 atm CH4 (see above)
        V_L = 1.0
        fugacity = pr.partial_pressures_atm(n_gas, T_K=self.T_K, V_L=V_L)["CH4"]
        partial_p = ideal.partial_pressures_atm(n_gas, T_K=self.T_K, V_L=V_L)["CH4"]
        assert fugacity == pytest.approx(45.495234588722575, rel=1e-9)
        assert partial_p == pytest.approx(54.87612153701265, rel=1e-9)
        # phi < 1 here (attractive-dominated non-ideality): fugacity < partial pressure.
        assert fugacity < partial_p
        assert fugacity / partial_p == pytest.approx(0.8291, abs=1e-3)

    def test_converge_at_low_pressure(self):
        from PyOMES.thermo import PengRobinsonEOS, IdealGasEOS, BIOGAS_SPECIES
        pr = PengRobinsonEOS(BIOGAS_SPECIES)
        ideal = IdealGasEOS()
        n_gas = {"CH4": 0.6, "CO2": 0.4}
        V_L = 25.213137354330655  # ~1 atm total (see above)
        fug = pr.partial_pressures_atm(n_gas, T_K=self.T_K, V_L=V_L)
        pp = ideal.partial_pressures_atm(n_gas, T_K=self.T_K, V_L=V_L)
        for sp in ("CH4", "CO2"):
            assert fug[sp] == pytest.approx(pp[sp], rel=1e-2)
