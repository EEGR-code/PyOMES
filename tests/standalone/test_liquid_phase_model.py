# -*- coding: utf-8 -*-
"""Tests for LiquidPhaseModel protocol and implementations."""
from __future__ import annotations

import pytest


class TestLiquidPhaseModelProtocol:
    def test_ideal_satisfies_protocol(self):
        from PyOMES.thermo import IdealLiquidModel, LiquidPhaseModel
        model = IdealLiquidModel()
        assert isinstance(model, LiquidPhaseModel)

    def test_ideal_has_name(self):
        from PyOMES.thermo import IdealLiquidModel
        assert IdealLiquidModel().name == "ideal"

    def test_ideal_gamma_all_returns_dict(self):
        from PyOMES.thermo import IdealLiquidModel
        result = IdealLiquidModel().gamma_all(
            {"CO2": 1.0, "HCO3-": 0.5},
            T_K=298.15,
            charge={"CO2": 0, "HCO3-": -1},
        )
        assert isinstance(result, dict)

    def test_absent_means_gamma_one(self):
        from PyOMES.thermo import IdealLiquidModel
        result = IdealLiquidModel().gamma_all(
            {"H+": 1e-7, "OH-": 1e-7},
            T_K=298.15,
            charge={"H+": 1, "OH-": -1},
        )
        assert result.get("H+", 1.0) == pytest.approx(1.0)
        assert result.get("OH-", 1.0) == pytest.approx(1.0)
        assert result.get("HCO3-", 1.0) == pytest.approx(1.0)

    def test_ideal_temperature_invariant(self):
        from PyOMES.thermo import IdealLiquidModel
        model = IdealLiquidModel()
        r_cold = model.gamma_all({"Ca++": 1e-3}, 280.0, charge={"Ca++": 2})
        r_hot  = model.gamma_all({"Ca++": 1e-3}, 370.0, charge={"Ca++": 2})
        assert r_cold == r_hot

    def test_ideal_composition_invariant(self):
        from PyOMES.thermo import IdealLiquidModel
        model = IdealLiquidModel()
        r_dilute = model.gamma_all({"NaCl": 0.001}, 298.15, charge={"NaCl": 0})
        r_conc   = model.gamma_all({"NaCl": 5.0},   298.15, charge={"NaCl": 0})
        assert r_dilute == r_conc

    def test_ideal_ignores_charge(self):
        from PyOMES.thermo import IdealLiquidModel
        model = IdealLiquidModel()
        r_with_charge    = model.gamma_all({"Ca++": 0.01}, 298.15, charge={"Ca++": 2})
        r_without_charge = model.gamma_all({"Ca++": 0.01}, 298.15, charge={})
        assert r_with_charge == r_without_charge

    def test_ideal_empty_composition(self):
        from PyOMES.thermo import IdealLiquidModel
        result = IdealLiquidModel().gamma_all({}, 298.15, charge={})
        assert result == {}

    def test_protocol_export_from_thermo(self):
        from PyOMES.thermo import LiquidPhaseModel, IdealLiquidModel
        assert LiquidPhaseModel is not None
        assert IdealLiquidModel is not None


class TestDaviesLiquidModel:
    def test_satisfies_liquid_phase_model_protocol(self):
        from PyOMES.thermo import DaviesLiquidModel, LiquidPhaseModel
        model = DaviesLiquidModel()
        assert isinstance(model, LiquidPhaseModel)

    def test_name(self):
        from PyOMES.thermo import DaviesLiquidModel
        assert DaviesLiquidModel().name == "davies"

    def test_frozen(self):
        import dataclasses
        from PyOMES.thermo import DaviesLiquidModel
        with pytest.raises((dataclasses.FrozenInstanceError, AttributeError)):
            DaviesLiquidModel().name = "other"  # type: ignore

    def test_neutral_species_absent_from_result(self):
        from PyOMES.thermo import DaviesLiquidModel
        model = DaviesLiquidModel()
        result = model.gamma_all(
            {"CO2": 1.0, "H+": 1e-3, "HCO3-": 1e-3},
            T_K=298.15,
            charge={"CO2": 0, "H+": 1, "HCO3-": -1},
        )
        assert "CO2" not in result
        assert result.get("CO2", 1.0) == pytest.approx(1.0)

    def test_ionic_species_present_and_below_one(self):
        from PyOMES.thermo import DaviesLiquidModel
        model = DaviesLiquidModel()
        result = model.gamma_all(
            {"H+": 0.01, "Cl-": 0.01},
            T_K=298.15,
            charge={"H+": 1, "Cl-": -1},
        )
        assert "H+" in result
        assert "Cl-" in result
        assert result["H+"] < 1.0   # Davies suppresses activity
        assert result["Cl-"] < 1.0

    def test_zero_concentration_ionic_excluded_at_zero_I(self):
        from PyOMES.thermo import DaviesLiquidModel
        model = DaviesLiquidModel()
        result = model.gamma_all(
            {"H+": 0.0, "Cl-": 0.0},
            T_K=298.15,
            charge={"H+": 1, "Cl-": -1},
        )
        # I=0 → gamma=1.0 (default) for all species; ionic species may be
        # included but must satisfy gamma = 1.0
        for v in result.values():
            assert v == pytest.approx(1.0)

    def test_higher_charge_larger_correction(self):
        from PyOMES.thermo import DaviesLiquidModel
        model = DaviesLiquidModel()
        # Same I, compare z=1 vs z=2
        result1 = model.gamma_all(
            {"Na+": 0.05, "Cl-": 0.05},
            T_K=298.15,
            charge={"Na+": 1, "Cl-": -1},
        )
        result2 = model.gamma_all(
            {"Ca++": 0.05, "SO4--": 0.05},
            T_K=298.15,
            charge={"Ca++": 2, "SO4--": -2},
        )
        # z=2 gets stronger suppression
        assert result2["Ca++"] < result1["Na+"]

    def test_temperature_affects_gamma(self):
        from PyOMES.thermo import DaviesLiquidModel
        model = DaviesLiquidModel()
        x = {"H+": 0.05, "Cl-": 0.05}
        ch = {"H+": 1, "Cl-": -1}
        g_25 = model.gamma_all(x, T_K=298.15, charge=ch)["H+"]
        g_50 = model.gamma_all(x, T_K=323.15, charge=ch)["H+"]
        assert g_25 != pytest.approx(g_50, rel=1e-3)

    # ── ActivityModel (per-ion) compatibility ──────────────────────────

    def test_gamma_at_zero_ionic_strength_is_one(self):
        from PyOMES.thermo import DaviesLiquidModel
        g = DaviesLiquidModel().gamma(1, 0.0, T_K=298.15)
        assert g == pytest.approx(1.0)

    def test_gamma_neutral_ion_is_one(self):
        from PyOMES.thermo import DaviesLiquidModel
        g = DaviesLiquidModel().gamma(0, 0.5, T_K=298.15)
        assert g == pytest.approx(1.0)

    def test_gamma_monovalent_decreases_with_I(self):
        from PyOMES.thermo import DaviesLiquidModel
        model = DaviesLiquidModel()
        g_low  = model.gamma(1, 0.01, T_K=298.15)
        g_high = model.gamma(1, 0.10, T_K=298.15)
        assert g_low > g_high  # more suppression at higher I (Davies regime)

    def test_gamma_all_consistency_with_gamma(self):
        """gamma_all at known composition must match gamma(z, I)."""
        from PyOMES.thermo import DaviesLiquidModel
        model = DaviesLiquidModel()
        x = {"H+": 0.05, "Cl-": 0.05}
        ch = {"H+": 1, "Cl-": -1}
        gammas = model.gamma_all(x, T_K=298.15, charge=ch)
        # I = 0.5*(1²*0.05 + 1²*0.05) = 0.05 mol/L
        expected = model.gamma(1, 0.05, T_K=298.15)
        assert gammas["H+"] == pytest.approx(expected, rel=1e-9)

    # ── Backward compat alias ──────────────────────────────────────────

    def test_davies_activity_model_alias(self):
        from PyOMES.chemical_equilibrium.activity_models import DaviesActivityModel, DaviesLiquidModel
        assert DaviesActivityModel is DaviesLiquidModel

    def test_make_activity_model_returns_davies_liquid_model(self):
        from PyOMES.chemical_equilibrium.activity_models import make_activity_model, DaviesLiquidModel
        m = make_activity_model(True, "davies")
        assert isinstance(m, DaviesLiquidModel)

    def test_water_helpers_still_importable_from_speciation(self):
        from PyOMES.chemical_equilibrium.activity_models import (
            debye_huckel_A,
            ionic_strength_molal_from_molar,
            water_density_kg_per_m3,
        )
        assert debye_huckel_A(298.15) == pytest.approx(0.509, rel=0.02)


class TestSITLiquidModel:
    def test_satisfies_liquid_phase_model_protocol(self):
        from PyOMES.thermo import SITLiquidModel, LiquidPhaseModel
        assert isinstance(SITLiquidModel(), LiquidPhaseModel)

    def test_name(self):
        from PyOMES.thermo import SITLiquidModel
        assert SITLiquidModel().name == "sit"

    def test_gamma_all_returns_ionic_species_only(self):
        from PyOMES.thermo import SITLiquidModel
        model = SITLiquidModel()
        result = model.gamma_all(
            {"Na+": 0.1, "Cl-": 0.1, "CO2": 0.01},
            T_K=298.15,
            charge={"Na+": 1, "Cl-": -1, "CO2": 0},
        )
        assert "Na+" in result
        assert "Cl-" in result
        assert "CO2" not in result

    def test_gamma_all_values_finite(self):
        from PyOMES.thermo import SITLiquidModel
        model = SITLiquidModel()
        result = model.gamma_all(
            {"Na+": 0.5, "Cl-": 0.5},
            T_K=298.15,
            charge={"Na+": 1, "Cl-": -1},
        )
        for v in result.values():
            assert v > 0.0
            import math
            assert math.isfinite(v)

    def test_gamma_all_empty_at_zero_composition(self):
        from PyOMES.thermo import SITLiquidModel
        result = SITLiquidModel().gamma_all({}, 298.15, charge={})
        assert result == {}

    def test_gamma_at_zero_I_is_one(self):
        from PyOMES.thermo import SITLiquidModel
        assert SITLiquidModel().gamma(1, 0.0, T_K=298.15) == pytest.approx(1.0)

    def test_gamma_neutral_is_one(self):
        from PyOMES.thermo import SITLiquidModel
        assert SITLiquidModel().gamma(0, 0.5, T_K=298.15) == pytest.approx(1.0)

    def test_compute_gammas_returns_ion_charges_keys(self):
        from PyOMES.thermo import SITLiquidModel
        from PyOMES.thermo.sit_liquid_model import ION_CHARGES
        model = SITLiquidModel()
        result = model.compute_gammas({"Na+": 0.1, "Cl-": 0.1}, 0.1, 298.15)
        for ion in ION_CHARGES:
            assert ion in result

    def test_compute_gammas_values_positive(self):
        from PyOMES.thermo import SITLiquidModel
        model = SITLiquidModel()
        result = model.compute_gammas({"Na+": 0.5, "Cl-": 0.5}, 0.5, 298.15)
        assert all(v > 0.0 for v in result.values())

    def test_sit_activity_model_alias(self):
        from PyOMES.chemical_equilibrium.sit import SITActivityModel, SITLiquidModel
        assert SITActivityModel is SITLiquidModel

    def test_make_activity_model_returns_sit(self):
        from PyOMES.chemical_equilibrium.activity_models import make_activity_model
        from PyOMES.thermo import SITLiquidModel
        m = make_activity_model(True, "sit")
        assert isinstance(m, SITLiquidModel)

    def test_sit_epsilon_importable_from_speciation(self):
        from PyOMES.chemical_equilibrium.sit import SIT_EPSILON, ION_CHARGES
        assert isinstance(SIT_EPSILON, dict)
        assert ("Na+", "Cl-") in SIT_EPSILON
        assert "H+" in ION_CHARGES
