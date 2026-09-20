# -*- coding: utf-8 -*-
"""Tighter coverage for 03_saturation_index.ipynb.

`tests/standalone/test_nr_speciation_engine.py::TestPrecipitationEquilibria`
covers the *precipitating* active-set engine (SI driven to ~0 once calcite
forms). This file covers the notebook's own scenario: a *non-precipitating*
(dissolved-only) engine, where SI is reported as a diagnostic that can sit
persistently above zero -- a different code path (no `precipitation_reactions`
on the engine at all) that the existing precipitation tests don't exercise.

Reference point is the notebook's own worked example: CT_Ca=0.002 mol/L,
CT_CO2=0.010 mol/L, CT_Na=0.005 mol/L, T=298.15 K, calcite Ksp=10**-8.48.
"""
from __future__ import annotations

import numpy as np
import pytest

from PyOMES.chemistry.common_species import (
    H2O, H_plus, OH_minus,
    CO2, HCO3_minus, CO3_2minus,
)
from PyOMES.reactions.equilibrium import EquilibriumReaction
from PyOMES.reactions.stoichiometry import StoichiometryEntry
from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
from PyOMES.thermo import DaviesLiquidModel

T_K = 298.15
KSP_CALCITE = 10 ** (-8.48)
CT_CA = 0.002
CT_CO2 = 0.010
CT_NA = 0.005


def _e(species, phase, coeff):
    return StoichiometryEntry(species=species, phase=phase, coefficient=coeff)


def _make_dissolved_only_engine():
    """No precipitation_reactions -- Ca2+ is a pure spectator strong ion,
    so SI can be computed and reported without being forced to zero."""
    water = EquilibriumReaction(
        stoichiometry=[_e(H2O, "liquid", -1), _e(H_plus, "liquid", +1), _e(OH_minus, "liquid", +1)],
        log_K=-14.0, label="water",
    )
    co2_first = EquilibriumReaction(
        stoichiometry=[_e(CO2, "liquid", -1), _e(H2O, "liquid", -1), _e(HCO3_minus, "liquid", +1), _e(H_plus, "liquid", +1)],
        log_K=-6.35, total_id="CO2", label="co2_first",
    )
    co2_second = EquilibriumReaction(
        stoichiometry=[_e(HCO3_minus, "liquid", -1), _e(CO3_2minus, "liquid", +1), _e(H_plus, "liquid", +1)],
        log_K=-10.33, total_id="CO2", label="co2_second",
    )
    return NRChemicalEquilibriumEngine.from_reactions(
        [water, co2_first, co2_second], use_activity=True, activity_model="davies",
    )


class TestDaviesGamma:
    def test_gamma_at_notebook_ionic_strength(self):
        davies = DaviesLiquidModel()
        # I = 0.01101 mol/L is the notebook's own solved ionic strength at
        # this recipe -- verified independently below in TestSaturationIndex.
        gamma = davies.gamma(2, 0.01101, T_K=T_K)
        assert gamma == pytest.approx(0.6493, abs=0.001)


class TestSaturationIndex:
    def _solve(self):
        engine = _make_dissolved_only_engine()
        return engine.solve(totals={"CO2": CT_CO2}, strong_ions={"CT_Ca": CT_CA, "CT_Na": CT_NA}, T_K=T_K)

    def test_ph_and_ionic_strength(self):
        out = self._solve()
        assert float(out.pH) == pytest.approx(7.2518, abs=0.001)
        assert float(out.ionic_strength) == pytest.approx(0.01101, abs=1e-4)

    def test_si_ideal(self):
        out = self._solve()
        c_co3 = float(out.species_mol_L["CO3--"])
        iap_ideal = CT_CA * c_co3
        si_ideal = np.log10(iap_ideal / KSP_CALCITE)
        assert si_ideal == pytest.approx(0.7967, abs=0.01)
        assert si_ideal > 0  # supersaturated

    def test_si_davies(self):
        out = self._solve()
        davies = DaviesLiquidModel()
        c_co3 = float(out.species_mol_L["CO3--"])
        gamma = davies.gamma(2, float(out.ionic_strength), T_K=T_K)
        iap_davies = (gamma * CT_CA) * (gamma * c_co3)
        si_davies = np.log10(iap_davies / KSP_CALCITE)
        assert si_davies == pytest.approx(0.4217, abs=0.01)
        assert si_davies > 0  # still supersaturated, but less so than ideal
