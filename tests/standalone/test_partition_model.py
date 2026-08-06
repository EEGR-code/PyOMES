# -*- coding: utf-8 -*-
"""Tests for PartitionModel protocol and HenryPartition (C2)."""
from __future__ import annotations

import dataclasses
import math
import pytest

# R in L·atm/(mol·K) — must match partition.py and core.phases.R_L_ATM_MOL_K
_R = 0.0820574
_T_REF = 298.15

# A representative H2S HenryPartition: kH ≈ 0.10 mol/(L·atm) at 298.15 K.
# H_ref (mol/m³/Pa) back-calculated: 0.10 * 1000 / 101325 ≈ 9.869e-4
_H2S_H_REF = 0.10 * 1000.0 / 101325.0
_H2S_DLN_H = 2100.0  # K (illustrative van't Hoff)


class TestHenryPartitionConstruction:
    def test_basic(self):
        from PyOMES.chemistry import HenryPartition
        hp = HenryPartition(H_ref=_H2S_H_REF, dlnH=_H2S_DLN_H)
        assert hp.H_ref == pytest.approx(_H2S_H_REF)
        assert hp.dlnH == _H2S_DLN_H
        assert hp.T_ref == _T_REF

    def test_custom_t_ref(self):
        from PyOMES.chemistry import HenryPartition
        hp = HenryPartition(H_ref=1e-4, dlnH=0.0, T_ref=310.0)
        assert hp.T_ref == 310.0

    def test_frozen(self):
        from PyOMES.chemistry import HenryPartition
        hp = HenryPartition(H_ref=1e-4, dlnH=0.0)
        with pytest.raises((dataclasses.FrozenInstanceError, AttributeError)):
            hp.H_ref = 99.0  # type: ignore

    def test_dataclasses_replace(self):
        from PyOMES.chemistry import HenryPartition
        hp = HenryPartition(H_ref=_H2S_H_REF, dlnH=_H2S_DLN_H)
        hp2 = dataclasses.replace(hp, dlnH=0.0)
        assert hp2.dlnH == 0.0
        assert hp2.H_ref == hp.H_ref


class TestHenryPartitionKH:
    def test_kH_at_ref_temperature(self):
        from PyOMES.chemistry import HenryPartition
        hp = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0)
        # With dlnH=0 the exp factor is 1; kH = H_ref/1000 * 101325
        expected = (_H2S_H_REF / 1000.0) * 101325.0
        assert hp._kH_mol_L_atm(_T_REF) == pytest.approx(expected, rel=1e-9)

    def test_kH_matches_construction_value(self):
        from PyOMES.chemistry import HenryPartition
        # Build so that kH at T_ref exactly equals 0.10 mol/L/atm
        hp = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0)
        assert hp._kH_mol_L_atm(_T_REF) == pytest.approx(0.10, rel=1e-6)

    def test_kH_temperature_dependence(self):
        from PyOMES.chemistry import HenryPartition
        hp = HenryPartition(H_ref=_H2S_H_REF, dlnH=_H2S_DLN_H)
        kH_ref = hp._kH_mol_L_atm(_T_REF)
        kH_hot = hp._kH_mol_L_atm(308.15)
        # Positive dlnH → more soluble at lower T → kH decreases as T rises
        assert kH_hot < kH_ref

    def test_kH_zero_dlnh_is_temperature_invariant(self):
        from PyOMES.chemistry import HenryPartition
        hp = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0)
        assert hp._kH_mol_L_atm(280.0) == pytest.approx(hp._kH_mol_L_atm(320.0))


class TestHenryPartitionPartitionRatio:
    def _hp(self):
        from PyOMES.chemistry import HenryPartition
        return HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0)

    def test_partition_ratio_formula(self):
        hp = self._hp()
        V_liq, V_gas = 2.0, 0.5
        T = _T_REF
        kH = hp._kH_mol_L_atm(T)
        expected = kH * _R * T * V_liq / V_gas
        assert hp.partition_ratio(V_liq, V_gas, T) == pytest.approx(expected, rel=1e-9)

    def test_partition_ratio_alpha_halved_doubles_ratio(self):
        hp = self._hp()
        b1 = hp.partition_ratio(1.0, 1.0, _T_REF, alpha=1.0)
        b2 = hp.partition_ratio(1.0, 1.0, _T_REF, alpha=0.5)
        assert b2 == pytest.approx(2.0 * b1, rel=1e-9)

    def test_partition_ratio_returns_float_not_none(self):
        # HenryPartition is a linear model — partition_ratio() must return a float
        hp = self._hp()
        result = hp.partition_ratio(1.0, 1.0, _T_REF)
        assert result is not None
        assert isinstance(result, float)

    def test_partition_ratio_larger_liquid_volume(self):
        hp = self._hp()
        b_small = hp.partition_ratio(1.0, 1.0, _T_REF)
        b_large = hp.partition_ratio(2.0, 1.0, _T_REF)
        assert b_large == pytest.approx(2.0 * b_small, rel=1e-9)

    def test_partition_ratio_alpha_near_zero_clamped(self):
        hp = self._hp()
        # alpha=0 would divide by zero; implementation clamps to 1e-12
        b = hp.partition_ratio(1.0, 1.0, _T_REF, alpha=0.0)
        assert math.isfinite(b)
        assert b > 0.0


class TestHenryPartitionEquilibriumMoles:
    def _hp(self):
        from PyOMES.chemistry import HenryPartition
        return HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0)

    def test_equilibrium_fraction(self):
        hp = self._hp()
        n_total = 1.0
        V_liq, V_gas = 1.0, 1.0
        b = hp.partition_ratio(V_liq, V_gas, _T_REF)
        expected = b * n_total / (1.0 + b)
        assert hp.equilibrium_a_moles(n_total, V_liq, V_gas, _T_REF) == pytest.approx(expected)

    def test_total_conservation(self):
        hp = self._hp()
        V_liq, V_gas, n_total = 2.0, 0.5, 3.0
        n_liq = hp.equilibrium_a_moles(n_total, V_liq, V_gas, _T_REF)
        assert 0.0 < n_liq < n_total

    def test_more_soluble_at_lower_temperature(self):
        from PyOMES.chemistry import HenryPartition
        hp = HenryPartition(H_ref=_H2S_H_REF, dlnH=_H2S_DLN_H)
        n_total, V_liq, V_gas = 1.0, 1.0, 1.0
        n_liq_cold = hp.equilibrium_a_moles(n_total, V_liq, V_gas, 280.0)
        n_liq_hot  = hp.equilibrium_a_moles(n_total, V_liq, V_gas, 320.0)
        assert n_liq_cold > n_liq_hot

    def test_alpha_correction_keeps_more_in_liquid(self):
        from PyOMES.chemistry import HenryPartition
        hp = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0)
        n_total, V_liq, V_gas = 1.0, 1.0, 1.0
        n_liq_no_alpha  = hp.equilibrium_a_moles(n_total, V_liq, V_gas, _T_REF, alpha=1.0)
        n_liq_with_alpha = hp.equilibrium_a_moles(n_total, V_liq, V_gas, _T_REF, alpha=0.5)
        # At alpha=0.5 (half ionised), effective Henry is higher → more stays in liquid
        assert n_liq_with_alpha > n_liq_no_alpha

    def test_high_solubility_mostly_liquid(self):
        from PyOMES.chemistry import HenryPartition
        # kH = 100 mol/L/atm → extremely soluble
        H_ref_high = 100.0 * 1000.0 / 101325.0
        hp = HenryPartition(H_ref=H_ref_high, dlnH=0.0)
        n_liq = hp.equilibrium_a_moles(1.0, 1.0, 1.0, _T_REF)
        assert n_liq > 0.99

    def test_low_solubility_mostly_gas(self):
        from PyOMES.chemistry import HenryPartition
        # kH = 1e-6 mol/L/atm → almost insoluble
        H_ref_low = 1e-6 * 1000.0 / 101325.0
        hp = HenryPartition(H_ref=H_ref_low, dlnH=0.0)
        n_liq = hp.equilibrium_a_moles(1.0, 1.0, 1.0, _T_REF)
        assert n_liq < 0.01


class TestPartitionModelProtocol:
    def test_henry_partition_satisfies_protocol(self):
        from PyOMES.chemistry import HenryPartition, PartitionModel
        hp = HenryPartition(H_ref=1e-4, dlnH=0.0)
        assert callable(hp.partition_ratio)
        assert callable(hp.equilibrium_a_moles)


class TestHenryPartitionActivityCorrection:
    """CP4: activity correction via thermo.liquid_activity."""

    _H2S_H_REF = _H2S_H_REF  # from module-level constant

    def _davies_thermo(self):
        from PyOMES.thermo import ThermoFramework, DaviesLiquidModel
        return ThermoFramework(liquid_activity=DaviesLiquidModel())

    def test_no_thermo_gamma_one(self):
        from PyOMES.chemistry import HenryPartition
        hp = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0)
        r_no_thermo = hp.partition_ratio(1.0, 1.0, _T_REF)
        r_with_empty = hp.partition_ratio(1.0, 1.0, _T_REF,
                                          species_id="H2S", x_mol={}, charge={})
        assert r_no_thermo == pytest.approx(r_with_empty)

    def test_ideal_thermo_no_change(self):
        from PyOMES.chemistry import HenryPartition
        from PyOMES.thermo import ThermoFramework
        thermo = ThermoFramework()  # IdealLiquidModel
        hp = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0, thermo=thermo)
        x = {"H2S": 0.01, "H+": 1e-7, "HS-": 1e-7}
        ch = {"H2S": 0, "H+": 1, "HS-": -1}
        r_ideal = hp.partition_ratio(1.0, 1.0, _T_REF,
                                     species_id="H2S", x_mol=x, charge=ch)
        r_base  = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0).partition_ratio(1.0, 1.0, _T_REF)
        assert r_ideal == pytest.approx(r_base)

    def test_davies_thermo_neutral_species_no_change(self):
        """Neutral species (z=0) get γ=1.0 from Davies → no correction."""
        from PyOMES.chemistry import HenryPartition
        thermo = self._davies_thermo()
        hp_with = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0, thermo=thermo)
        hp_bare = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0)
        x = {"H2S": 0.01, "H+": 1e-7, "HS-": 1e-7}
        ch = {"H2S": 0, "H+": 1, "HS-": -1}
        r_with = hp_with.partition_ratio(1.0, 1.0, _T_REF,
                                         species_id="H2S", x_mol=x, charge=ch)
        r_bare = hp_bare.partition_ratio(1.0, 1.0, _T_REF)
        assert r_with == pytest.approx(r_bare)

    def test_davies_thermo_ionic_reduces_partition(self):
        """Ionic species with γ<1 → effective kH / γ > kH → more stays liquid."""
        from PyOMES.chemistry import HenryPartition
        thermo = self._davies_thermo()
        # HS- is ionic (z=-1); at I>0, Davies gives γ<1
        hp_with = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0, thermo=thermo)
        hp_bare = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0)
        x = {"HS-": 0.05, "Na+": 0.05}
        ch = {"HS-": -1, "Na+": 1}
        r_with = hp_with.partition_ratio(1.0, 1.0, _T_REF,
                                         species_id="HS-", x_mol=x, charge=ch)
        r_bare = hp_bare.partition_ratio(1.0, 1.0, _T_REF)
        assert r_with > r_bare  # γ<1 → divide by smaller → larger ratio

    def test_activity_correction_equilibrium_moles(self):
        """equilibrium_a_moles passes activity kwargs through."""
        from PyOMES.chemistry import HenryPartition
        thermo = self._davies_thermo()
        hp = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0, thermo=thermo)
        x = {"HS-": 0.05, "Na+": 0.05}
        ch = {"HS-": -1, "Na+": 1}
        n_liq = hp.equilibrium_a_moles(1.0, 1.0, 1.0, _T_REF,
                                        species_id="HS-", x_mol=x, charge=ch)
        n_liq_bare = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0).equilibrium_a_moles(
            1.0, 1.0, 1.0, _T_REF)
        assert n_liq > n_liq_bare

    def test_thermo_field_excluded_from_equality(self):
        """thermo field has compare=False so equal H_ref/dlnH/T_ref instances compare equal."""
        from PyOMES.chemistry import HenryPartition
        from PyOMES.thermo import ThermoFramework, DaviesLiquidModel
        hp1 = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0,
                              thermo=ThermoFramework(liquid_activity=DaviesLiquidModel()))
        hp2 = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0)
        assert hp1 == hp2


class TestRaoultPartition:
    _R = 0.0820574  # L·atm/(mol·K)
    _C_W = 55.51     # mol/L pure water

    def test_construction_defaults(self):
        from PyOMES.chemistry import RaoultPartition
        rp = RaoultPartition()
        assert rp.P_sat_ref == pytest.approx(0.03169, rel=1e-3)
        assert rp.T_ref == pytest.approx(298.15)

    def test_psat_at_ref_returns_ref(self):
        from PyOMES.chemistry import RaoultPartition
        rp = RaoultPartition()
        assert rp.P_sat(298.15) == pytest.approx(rp.P_sat_ref, rel=1e-9)

    def test_psat_increases_with_temperature(self):
        from PyOMES.chemistry import RaoultPartition
        rp = RaoultPartition()
        assert rp.P_sat(373.15) > rp.P_sat(298.15)

    def test_psat_at_100c_near_one_atm(self):
        from PyOMES.chemistry import RaoultPartition
        rp = RaoultPartition()
        # Clausius-Clapeyron gives ~0.94 atm at 100°C (empirically ~1 atm)
        assert 0.8 < rp.P_sat(373.15) < 1.5

    def test_partition_ratio_formula(self):
        from PyOMES.chemistry import RaoultPartition
        rp = RaoultPartition()
        T = 298.15
        V_liq, V_gas = 2.0, 0.1
        Ps = rp.P_sat(T)
        expected = (rp.C_water_mol_L * V_liq * self._R * T) / (Ps * V_gas)
        assert rp.partition_ratio(V_liq, V_gas, T) == pytest.approx(expected, rel=1e-9)

    def test_partition_ratio_large_at_low_T(self):
        """Water is mostly liquid at low T (low P_sat → large ratio)."""
        from PyOMES.chemistry import RaoultPartition
        rp = RaoultPartition()
        r_cold = rp.partition_ratio(1.0, 1.0, 280.0)
        r_hot  = rp.partition_ratio(1.0, 1.0, 370.0)
        assert r_cold > r_hot

    def test_equilibrium_moles_mostly_liquid_at_ambient(self):
        from PyOMES.chemistry import RaoultPartition
        rp = RaoultPartition()
        n_liq = rp.equilibrium_a_moles(1.0, 1.0, 0.001, 298.15)
        assert n_liq > 0.99

    def test_alpha_ignored(self):
        """alpha has no physical meaning for water — partition_ratio is invariant."""
        from PyOMES.chemistry import RaoultPartition
        rp = RaoultPartition()
        r1 = rp.partition_ratio(1.0, 1.0, 298.15, alpha=1.0)
        r2 = rp.partition_ratio(1.0, 1.0, 298.15, alpha=0.5)
        assert r1 == pytest.approx(r2)

    def test_satisfies_partition_model_protocol(self):
        from PyOMES.chemistry import RaoultPartition, PartitionModel
        rp = RaoultPartition()
        assert callable(rp.partition_ratio)
        assert callable(rp.equilibrium_a_moles)

    def test_partition_ratio_returns_float(self):
        from PyOMES.chemistry import RaoultPartition
        r = RaoultPartition().partition_ratio(1.0, 1.0, 298.15)
        assert isinstance(r, float)
        assert r > 0.0

    def test_frozen(self):
        from PyOMES.chemistry import RaoultPartition
        import dataclasses
        rp = RaoultPartition()
        with pytest.raises((dataclasses.FrozenInstanceError, AttributeError)):
            rp.P_sat_ref = 0.1  # type: ignore


# ── Module-level constants for MultispeciesVLE tests ─────────────────────────
_H_H2S = 0.10 * 1000.0 / 101325.0   # 0.10 mol/(L·atm) → Sander mol/(m³·Pa)
_H_CO2 = 3.4e-4 * 1000.0 / 101325.0  # 3.4e-4 mol/(L·atm) → Sander units


class TestMultispeciesPartitionModelProtocol:
    def test_multispecies_vle_satisfies_protocol(self):
        from PyOMES.chemistry import MultispeciesVLEPartition, MultispeciesPartitionModel
        vle = MultispeciesVLEPartition(kH_ref={"H2S": _H_H2S})
        assert isinstance(vle, MultispeciesPartitionModel)

    def test_protocol_has_equilibrium_all_a_moles(self):
        from PyOMES.chemistry import MultispeciesVLEPartition
        vle = MultispeciesVLEPartition(kH_ref={"H2S": _H_H2S})
        assert callable(vle.equilibrium_all_a_moles)


class TestMultispeciesVLEPartition:
    def _vle(self):
        from PyOMES.chemistry import MultispeciesVLEPartition
        return MultispeciesVLEPartition(
            kH_ref={"H2S": _H_H2S, "CO2": _H_CO2},
            dlnH={"H2S": 2100.0, "CO2": 2400.0},
        )

    def test_construction(self):
        vle = self._vle()
        assert "H2S" in vle.kH_ref
        assert "CO2" in vle.kH_ref
        assert vle.T_ref == pytest.approx(298.15)

    def test_returns_dict_for_known_species(self):
        vle = self._vle()
        result = vle.equilibrium_all_a_moles(
            {"H2S": 0.01, "CO2": 0.02}, capacity_a=1.0, capacity_b=0.1, T_K=308.15
        )
        assert "H2S" in result
        assert "CO2" in result

    def test_unknown_species_absent(self):
        vle = self._vle()
        result = vle.equilibrium_all_a_moles(
            {"H2S": 0.01, "N2": 0.1}, capacity_a=1.0, capacity_b=0.1, T_K=308.15
        )
        assert "H2S" in result
        assert "N2" not in result

    def test_moles_conserved_per_species(self):
        vle = self._vle()
        n_in = {"H2S": 0.01, "CO2": 0.02}
        result = vle.equilibrium_all_a_moles(n_in, capacity_a=1.0, capacity_b=0.5, T_K=308.15)
        for sp, n_tot in n_in.items():
            assert 0.0 < result[sp] <= n_tot

    def test_more_liquid_at_larger_liquid_volume(self):
        vle = self._vle()
        r_small = vle.equilibrium_all_a_moles({"H2S": 0.01}, 1.0, 1.0, 308.15)
        r_large = vle.equilibrium_all_a_moles({"H2S": 0.01}, 2.0, 1.0, 308.15)
        assert r_large["H2S"] > r_small["H2S"]

    def test_agrees_with_henry_partition(self):
        """IdealGasEOS path must match per-species HenryPartition exactly."""
        from PyOMES.chemistry import HenryPartition, MultispeciesVLEPartition
        vle = MultispeciesVLEPartition(kH_ref={"H2S": _H_H2S}, dlnH={"H2S": 2100.0})
        hp  = HenryPartition(H_ref=_H_H2S, dlnH=2100.0)

        n_tot, V_liq, V_gas, T = 0.05, 2.0, 0.5, 308.15
        n_liq_multi = vle.equilibrium_all_a_moles({"H2S": n_tot}, V_liq, V_gas, T)["H2S"]
        n_liq_henry = hp.equilibrium_a_moles(n_tot, V_liq, V_gas, T)
        assert n_liq_multi == pytest.approx(n_liq_henry, rel=1e-9)

    def test_more_soluble_species_has_more_in_liquid(self):
        """CO2 (lower kH) stays more in gas than H2S (higher kH) at same conditions."""
        vle = self._vle()
        result = vle.equilibrium_all_a_moles(
            {"H2S": 0.01, "CO2": 0.01}, capacity_a=1.0, capacity_b=0.5, T_K=298.15
        )
        assert result["H2S"] > result["CO2"]

    def test_zero_total_gives_zero(self):
        vle = self._vle()
        result = vle.equilibrium_all_a_moles({"H2S": 0.0}, 1.0, 1.0, 298.15)
        assert result["H2S"] == pytest.approx(0.0)

    def test_frozen(self):
        import dataclasses
        from PyOMES.chemistry import MultispeciesVLEPartition
        vle = MultispeciesVLEPartition(kH_ref={"H2S": _H_H2S})
        with pytest.raises((dataclasses.FrozenInstanceError, AttributeError)):
            vle.T_ref = 310.0  # type: ignore
