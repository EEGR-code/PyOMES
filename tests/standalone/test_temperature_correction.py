"""Tests for PyOMES.thermo.temperature_correction and arrhenius_factor.

The ``_legacy_*`` functions below write the van 't Hoff formulas out in full,
as independent references, so ``vant_hoff_K`` and ``vant_hoff_log_K`` are
pinned to *exact* (bit-for-bit) equality, not just approximate agreement.
"""

import itertools
import math
from types import SimpleNamespace

import numpy as np
import pytest

from PyOMES.reactions import arrhenius_factor
from PyOMES.reactions.equilibrium.constraint import constraint_log_K_at
from PyOMES.thermo.temperature_correction import (
    clausius_clapeyron,
    henry_constant,
    ln_correction,
    vant_hoff_K,
    vant_hoff_log_K,
)
from PyOMES.units import R_J_PER_MOL_K

_LOG10_E = np.log10(np.e)


def _legacy_vant_hoff_K(K_ref, dH_J_per_mol, T_K, T_ref_K=298.15):
    K_ref = float(K_ref)
    dH_J_per_mol = float(dH_J_per_mol)
    T_K = float(T_K)
    T_ref_K = float(T_ref_K)
    if not np.isfinite(K_ref) or K_ref <= 0.0:
        return float(K_ref)
    if not np.isfinite(dH_J_per_mol) or abs(dH_J_per_mol) < 1e-30:
        return float(K_ref)
    if not np.isfinite(T_K) or T_K <= 0.0:
        return float(K_ref)
    return float(K_ref * np.exp(-(dH_J_per_mol / R_J_PER_MOL_K) * (1.0 / T_K - 1.0 / T_ref_K)))


def _legacy_vant_hoff_log_K(log_K_ref, dH_J_per_mol, T_K, T_ref_K):
    if dH_J_per_mol is None or abs(dH_J_per_mol) < 1e-30:
        return float(log_K_ref)
    if abs(T_K - T_ref_K) < 1e-10:
        return float(log_K_ref)
    delta_ln_K = -(float(dH_J_per_mol) / R_J_PER_MOL_K) * (1.0 / float(T_K) - 1.0 / float(T_ref_K))
    return float(log_K_ref) + delta_ln_K * _LOG10_E


_K_REFS = [1e-14, 4.47e-7, 5.6e-11, 1.0, 3.3e-9]
_LOG_KS = [-14.0, -6.35, -10.33, -8.48, 0.0, 2.5]
_DHS = [55830.0, 9160.0, -12100.0, -22000.0, 1.0e5, 0.0]
_TEMPS = [273.15, 278.15, 288.15, 298.15, 310.15, 323.15, 373.15]
_T_REFS = [298.15, 293.15]


class TestBitIdenticalToLegacy:
    def test_vant_hoff_K_matches_legacy_exactly(self):
        for K, dH, T, Tr in itertools.product(_K_REFS, _DHS, _TEMPS, _T_REFS):
            assert vant_hoff_K(K, dH, T, Tr) == _legacy_vant_hoff_K(K, dH, T, Tr)

    def test_vant_hoff_K_default_reference_temperature(self):
        assert vant_hoff_K(1e-14, 55830.0, 310.15) == _legacy_vant_hoff_K(1e-14, 55830.0, 310.15)

    def test_vant_hoff_log_K_matches_legacy_exactly(self):
        for lk, dH, T, Tr in itertools.product(_LOG_KS, _DHS + [None], _TEMPS, _T_REFS):
            assert vant_hoff_log_K(lk, dH, T, Tr) == _legacy_vant_hoff_log_K(lk, dH, T, Tr)


class TestVantHoffK:
    def test_zero_enthalpy_is_identity(self):
        assert vant_hoff_K(1e-14, 0.0, 330.0) == 1e-14

    def test_reference_temperature_is_identity(self):
        assert vant_hoff_K(1e-14, 55830.0, 298.15) == 1e-14

    def test_endothermic_increases_K_with_temperature(self):
        assert vant_hoff_K(1e-14, 55830.0, 323.15) > 1e-14
        assert vant_hoff_K(1e-14, 55830.0, 278.15) < 1e-14

    def test_exothermic_decreases_K_with_temperature(self):
        assert vant_hoff_K(1e-3, -20000.0, 323.15) < 1e-3

    @pytest.mark.parametrize("K_ref", [0.0, -1.0, float("nan"), float("inf")])
    def test_bad_K_ref_returned_unchanged(self, K_ref):
        out = vant_hoff_K(K_ref, 55830.0, 310.15)
        assert out == K_ref or (math.isnan(out) and math.isnan(K_ref))

    @pytest.mark.parametrize("dH", [float("nan"), float("inf")])
    def test_non_finite_enthalpy_returns_K_ref(self, dH):
        assert vant_hoff_K(1e-14, dH, 310.15) == 1e-14

    @pytest.mark.parametrize("T_K", [0.0, -5.0, float("nan")])
    def test_bad_temperature_returns_K_ref(self, T_K):
        assert vant_hoff_K(1e-14, 55830.0, T_K) == 1e-14


class TestVantHoffLogK:
    def test_none_enthalpy_is_identity(self):
        assert vant_hoff_log_K(-8.48, None, 330.0, 298.15) == -8.48

    def test_zero_enthalpy_is_identity(self):
        assert vant_hoff_log_K(-8.48, 0.0, 330.0, 298.15) == -8.48

    def test_reference_temperature_is_identity(self):
        assert vant_hoff_log_K(-8.48, -12100.0, 298.15, 298.15) == -8.48

    def test_matches_manual_calculation(self):
        dH, T = 55830.0, 323.15
        expected = -14.0 + (-(dH / R_J_PER_MOL_K) * (1.0 / T - 1.0 / 298.15)) / math.log(10.0)
        assert vant_hoff_log_K(-14.0, dH, T, 298.15) == pytest.approx(expected, rel=1e-12)

    def test_endothermic_increases_log_K_with_temperature(self):
        assert vant_hoff_log_K(-14.0, 55830.0, 323.15, 298.15) > -14.0


class TestConsistencyBetweenForms:
    def test_K_and_log_K_forms_agree(self):
        for K, dH, T in itertools.product(_K_REFS, _DHS, _TEMPS):
            via_K = math.log10(vant_hoff_K(K, dH, T, 298.15))
            via_log = vant_hoff_log_K(math.log10(K), dH, T, 298.15)
            assert via_K == pytest.approx(via_log, abs=1e-12)

class TestLnCorrection:
    def test_sign_and_zero(self):
        E = 55830.0 / R_J_PER_MOL_K
        assert ln_correction(E, 298.15, 298.15) == 0.0
        assert ln_correction(E, 323.15, 298.15) > 0.0
        assert ln_correction(-E, 323.15, 298.15) < 0.0

    def test_accepts_arrays(self):
        T = np.array([280.0, 298.15, 330.0])
        out = ln_correction(2000.0, T, 298.15)
        assert out.shape == (3,)
        assert out[1] == 0.0
        assert [ln_correction(2000.0, float(t), 298.15) for t in T] == list(out)


# Each of these is the expression a copy used before it called the kernel,
# written out in full, so the wrappers are pinned to exact equality with it.
_TEMPS_FINE = [273.15, 290.0, 298.15, 298.16, 308.15, 330.0, 373.15]


class TestHenryConstant:
    def test_matches_sander_form_exactly(self):
        for kH, dlnH, T in itertools.product((3.3e-4, 1.3e-5, 1.0), (-2400.0, 0.0, 1700.0, 2400.0), _TEMPS_FINE):
            assert henry_constant(kH, dlnH, T, 298.15) == kH * math.exp(dlnH * (1.0 / T - 1.0 / 298.15))

    def test_more_soluble_when_cold_for_positive_dlnH(self):
        assert henry_constant(3.3e-4, 2400.0, 283.15, 298.15) > 3.3e-4


class TestClausiusClapeyron:
    def test_matches_both_written_forms_exactly(self):
        for T in _TEMPS_FINE:
            # (ΔH_vap/R)·(1/T_ref − 1/T), with ΔH_vap/R given in kelvin ...
            assert clausius_clapeyron(0.0313, 5290.0, T, 298.15) == 0.0313 * math.exp(5290.0 * (1.0 / 298.15 - 1.0 / T))
            # ... and −ΔH_vap/R·(1/T − 1/T_ref), with ΔH_vap in J/mol.
            assert clausius_clapeyron(0.03169, 44011.0 / R_J_PER_MOL_K, T, 298.15) == (
                0.03169 * math.exp(-44011.0 / R_J_PER_MOL_K * (1.0 / T - 1.0 / 298.15)))

    def test_vapour_pressure_rises_with_temperature(self):
        assert clausius_clapeyron(0.0313, 5290.0, 308.15, 298.15) > 0.0313


class TestArrheniusFactor:
    def test_matches_written_form_exactly(self):
        for Ea_R, T in itertools.product((0.0, 2000.0, 7000.0), _TEMPS_FINE):
            assert arrhenius_factor(Ea_R, T, 308.15) == math.exp(Ea_R * (1.0 / 308.15 - 1.0 / T))

    def test_one_at_reference_and_faster_when_warmer(self):
        assert arrhenius_factor(7000.0, 308.15, 308.15) == 1.0
        assert arrhenius_factor(7000.0, 318.15, 308.15) > 1.0


class TestConstraintLogKAt:
    def _c(self, log_K, dH, T_ref=298.15):
        return SimpleNamespace(log_K=log_K, dH_J_per_mol=dH, T_ref_K=T_ref, stoichiometry=())

    def test_equals_vant_hoff_log_K_exactly(self):
        for lk, dH, T, Tr in itertools.product(_LOG_KS, _DHS + [None], _TEMPS, _T_REFS):
            assert constraint_log_K_at(self._c(lk, dH, Tr), T) == vant_hoff_log_K(lk, dH, T, Tr)

    def test_no_enthalpy_is_identity(self):
        assert constraint_log_K_at(self._c(-8.48, None), 330.0) == -8.48
