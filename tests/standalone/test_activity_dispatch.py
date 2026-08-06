# -*- coding: utf-8 -*-
"""Tests for CP4: activity_for_entry() per-entry activity dispatch helper.

Standalone verification — activity_for_entry is not wired into
build_tableau()'s Newton residual yet (Phase 2's job); these tests only
prove the phase-based dispatch matches calling gamma_all()/
partial_pressures_atm() directly.
"""
from __future__ import annotations

import pytest


def _entry(species_id, phase, charge=0, atoms=None, coefficient=-1.0):
    from PyOMES.chemistry.species import Species
    from PyOMES.reactions.stoichiometry import StoichiometryEntry
    sp = Species(id=species_id, atoms=atoms or {"H": 1}, charge=charge)
    return StoichiometryEntry(species=sp, phase=phase, coefficient=coefficient)


class TestSolidEntry:
    def test_returns_one_unconditionally(self):
        from PyOMES.chemical_equilibrium.activity_dispatch import activity_for_entry
        entry = _entry("CaCO3", "solid")
        # Garbage phases/thermo — solid dispatch must not touch either.
        result = activity_for_entry(entry, phases={}, thermo=None, T_K=298.15)
        assert result == 1.0

    def test_ignores_charge_kwarg(self):
        from PyOMES.chemical_equilibrium.activity_dispatch import activity_for_entry
        entry = _entry("CaCO3", "solid", charge=2)
        result = activity_for_entry(
            entry, phases={}, thermo=None, T_K=298.15, charge={"CaCO3": 2},
        )
        assert result == 1.0


class TestLiquidEntry:
    def _liquid_phase(self, n_mol):
        from PyOMES.core.phases import LiquidPhase
        return LiquidPhase(n_mol=n_mol, V_L=1.0, T_K=298.15)

    def test_ideal_liquid_model_gives_gamma_one(self):
        from PyOMES.thermo import ThermoFramework
        from PyOMES.chemical_equilibrium.activity_dispatch import activity_for_entry
        entry = _entry("HS-", "liquid", charge=-1)
        phases = {"liquid": self._liquid_phase({"HS-": 0.05, "Na+": 0.05})}
        thermo = ThermoFramework()  # IdealLiquidModel
        result = activity_for_entry(entry, phases, thermo, T_K=298.15)
        assert result == pytest.approx(1.0)

    def test_matches_gamma_all_directly_with_explicit_charge_map(self):
        """Davies with a full charge map — must exactly match calling
        gamma_all() with the same x_mol/charge/T_K."""
        from PyOMES.thermo import ThermoFramework, DaviesLiquidModel
        from PyOMES.chemical_equilibrium.activity_dispatch import activity_for_entry

        x_mol = {"HS-": 0.05, "Na+": 0.05}
        charge_map = {"HS-": -1, "Na+": +1}
        thermo = ThermoFramework(liquid_activity=DaviesLiquidModel())
        entry = _entry("HS-", "liquid", charge=-1)
        phases = {"liquid": self._liquid_phase(x_mol)}

        result = activity_for_entry(
            entry, phases, thermo, T_K=298.15, charge=charge_map,
        )
        expected = thermo.liquid_activity.gamma_all(
            x_mol, 298.15, charge=charge_map,
        )["HS-"]
        assert result == pytest.approx(expected, rel=1e-12)
        assert result != pytest.approx(1.0)  # sanity: correction is non-trivial

    def test_fallback_charge_uses_only_this_species(self):
        """Without an explicit charge= map, the fallback treats every
        other species in x_mol as neutral (charge 0) — documented
        under-estimate of ionic strength when multiple ions are present."""
        from PyOMES.thermo import ThermoFramework, DaviesLiquidModel
        from PyOMES.chemical_equilibrium.activity_dispatch import activity_for_entry

        x_mol = {"HS-": 0.05, "Na+": 0.05}
        thermo = ThermoFramework(liquid_activity=DaviesLiquidModel())
        entry = _entry("HS-", "liquid", charge=-1)
        phases = {"liquid": self._liquid_phase(x_mol)}

        result_fallback = activity_for_entry(entry, phases, thermo, T_K=298.15)
        result_full_charge = activity_for_entry(
            entry, phases, thermo, T_K=298.15, charge={"HS-": -1, "Na+": +1},
        )
        expected_fallback = thermo.liquid_activity.gamma_all(
            x_mol, 298.15, charge={"HS-": -1},
        )["HS-"]
        assert result_fallback == pytest.approx(expected_fallback, rel=1e-12)
        # Fallback under-counts ionic strength (missing Na+'s contribution)
        # relative to the full-charge-map case — a real, distinct path.
        assert result_fallback != pytest.approx(result_full_charge)

    def test_neutral_species_gamma_one_under_davies(self):
        from PyOMES.thermo import ThermoFramework, DaviesLiquidModel
        from PyOMES.chemical_equilibrium.activity_dispatch import activity_for_entry
        thermo = ThermoFramework(liquid_activity=DaviesLiquidModel())
        entry = _entry("H2S", "liquid", charge=0)
        phases = {"liquid": self._liquid_phase({"H2S": 0.01, "H+": 1e-7, "HS-": 1e-7})}
        result = activity_for_entry(
            entry, phases, thermo, T_K=298.15,
            charge={"H2S": 0, "H+": 1, "HS-": -1},
        )
        assert result == pytest.approx(1.0)

    def test_missing_liquid_phase_key_raises_keyerror(self):
        from PyOMES.thermo import ThermoFramework
        from PyOMES.chemical_equilibrium.activity_dispatch import activity_for_entry
        entry = _entry("HS-", "liquid", charge=-1)
        with pytest.raises(KeyError):
            activity_for_entry(entry, phases={}, thermo=ThermoFramework(), T_K=298.15)


class TestGasEntry:
    def _gas_phase(self, n_mol, V_L=0.4):
        from PyOMES.core.phases import GasPhase
        return GasPhase(n_mol=n_mol, V_L=V_L, T_K=305.15)

    def test_matches_partial_pressures_atm_directly_with_explicit_eos(self):
        from PyOMES.equilibria.vle import IdealGasEOS
        from PyOMES.thermo import ThermoFramework
        from PyOMES.chemical_equilibrium.activity_dispatch import activity_for_entry

        n_mol = {"CO2": 0.01, "O2": 0.05}
        gas = self._gas_phase(n_mol)
        thermo = ThermoFramework(gas_eos=IdealGasEOS())
        entry = _entry("CO2", "gas", charge=0)
        phases = {"gas": gas}

        result = activity_for_entry(entry, phases, thermo, T_K=305.15)
        expected = IdealGasEOS().partial_pressures_atm(
            n_mol, T_K=305.15, V_L=gas.V_L,
        )["CO2"]
        assert result == pytest.approx(expected, rel=1e-12)

    def test_falls_back_to_ideal_gas_eos_when_thermo_gas_eos_is_none(self):
        from PyOMES.equilibria.vle import IdealGasEOS
        from PyOMES.thermo import ThermoFramework
        from PyOMES.chemical_equilibrium.activity_dispatch import activity_for_entry

        n_mol = {"CO2": 0.01}
        gas = self._gas_phase(n_mol)
        thermo = ThermoFramework()  # gas_eos=None by default
        entry = _entry("CO2", "gas", charge=0)
        phases = {"gas": gas}

        result = activity_for_entry(entry, phases, thermo, T_K=305.15)
        expected = IdealGasEOS().partial_pressures_atm(
            n_mol, T_K=305.15, V_L=gas.V_L,
        )["CO2"]
        assert result == pytest.approx(expected, rel=1e-12)

    def test_species_absent_from_gas_phase_returns_zero(self):
        from PyOMES.thermo import ThermoFramework
        from PyOMES.chemical_equilibrium.activity_dispatch import activity_for_entry
        gas = self._gas_phase({"O2": 0.05})
        thermo = ThermoFramework()
        entry = _entry("CO2", "gas", charge=0)
        result = activity_for_entry(entry, {"gas": gas}, thermo, T_K=305.15)
        assert result == 0.0

    def test_missing_gas_phase_key_raises_keyerror(self):
        from PyOMES.thermo import ThermoFramework
        from PyOMES.chemical_equilibrium.activity_dispatch import activity_for_entry
        entry = _entry("CO2", "gas", charge=0)
        with pytest.raises(KeyError):
            activity_for_entry(entry, phases={}, thermo=ThermoFramework(), T_K=305.15)


class TestInvalidPhase:
    def test_unrecognised_phase_raises_value_error(self):
        from PyOMES.thermo import ThermoFramework
        from PyOMES.chemical_equilibrium.activity_dispatch import activity_for_entry
        entry = _entry("X", "plasma", charge=0)
        with pytest.raises(ValueError):
            activity_for_entry(entry, phases={}, thermo=ThermoFramework(), T_K=298.15)
