"""Tests for PyOMES.thermo.equilibrium_constants — shared van 't Hoff helpers.

The two wrappers replace ``acid_base._vant_hoff_K`` and
``nr_tableau._vant_hoff_log_K``. The ``_legacy_*`` functions below are those
originals, kept verbatim as references so the refactor is pinned to *exact*
(bit-for-bit) equality, not just approximate agreement.
"""

import itertools
import math

import numpy as np
import pytest

from PyOMES.thermo.equilibrium_constants import (
    vant_hoff_K,
    vant_hoff_delta_ln_K,
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

    def test_delta_ln_K_sign_and_zero(self):
        assert vant_hoff_delta_ln_K(55830.0, 298.15, 298.15) == 0.0
        assert vant_hoff_delta_ln_K(55830.0, 323.15, 298.15) > 0.0
        assert vant_hoff_delta_ln_K(-55830.0, 323.15, 298.15) < 0.0
