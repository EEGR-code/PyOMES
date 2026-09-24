# -*- coding: utf-8 -*-
"""NR-vs-PHREEQC agreement, closing a real gap in 06_phreeqc_benchmark.ipynb's coverage.

`tests/standalone/test_phreeqc_engine.py` only checks PHREEQCChemicalEquilibriumEngine's
internal self-consistency (pH range, logH == -pH, warmstart-vs-fresh agreement) --
nothing in the pytest suite compares NR's and PHREEQC's *solved values* against
each other for the same chemistry, which is exactly what the notebook does.

Uses PHREEQC's own log K values (not PyOMES's textbook defaults) to isolate
solver-numerical agreement from thermodynamic-data differences, per the
notebook's own stated approach. Tolerances are the notebook's own printed
per-section maxima (pure carbonate ~1.1e-3, carbonate+NH3 ~1.5e-3), not the
notebook's Section 2/6 numbers -- those are dominated by a documented,
model-difference (missing Na+/Ca2+ ion pairs in the NR engine), not solver
disagreement, and are explicitly excluded here as out of scope.
"""
from __future__ import annotations

import pytest

from PyOMES.chemistry.common_species import (
    H2O, H_plus, OH_minus,
    CO2, HCO3_minus, CO3_2minus,
    NH3, NH4_plus,
)
from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
from PyOMES.chemical_equilibrium.engines.nr.engine import NRChemicalEquilibriumEngine
from PyOMES.chemical_equilibrium.engines.phreeqc import PHREEQCChemicalEquilibriumEngine

T_K = 298.15
PH_TOL_CARBONATE = 2e-3
PH_TOL_CARBONATE_NH3 = 2e-3


def _e(species, phase, coeff):
    return StoichiometryEntry(species=species, phase=phase, coefficient=coeff)


def _make_nr_carbonate_engine():
    # PHREEQC's own phreeqc.dat log K values, per notebook 06.
    water = EquilibriumReaction(
        stoichiometry=[_e(H2O, "liquid", -1), _e(H_plus, "liquid", +1), _e(OH_minus, "liquid", +1)],
        log_K=-13.998, label="water",
    )
    co2_first = EquilibriumReaction(
        stoichiometry=[_e(CO2, "liquid", -1), _e(H2O, "liquid", -1), _e(HCO3_minus, "liquid", +1), _e(H_plus, "liquid", +1)],
        log_K=-6.352, total_id="CO2", label="co2_first",
    )
    co2_second = EquilibriumReaction(
        stoichiometry=[_e(HCO3_minus, "liquid", -1), _e(CO3_2minus, "liquid", +1), _e(H_plus, "liquid", +1)],
        log_K=-10.329, total_id="CO2", label="co2_second",
    )
    return NRChemicalEquilibriumEngine.from_reactions(
        [water, co2_first, co2_second], use_activity=True, activity_model="davies",
    )


def _make_nr_carbonate_nh3_engine():
    water = EquilibriumReaction(
        stoichiometry=[_e(H2O, "liquid", -1), _e(H_plus, "liquid", +1), _e(OH_minus, "liquid", +1)],
        log_K=-13.998, label="water",
    )
    co2_first = EquilibriumReaction(
        stoichiometry=[_e(CO2, "liquid", -1), _e(H2O, "liquid", -1), _e(HCO3_minus, "liquid", +1), _e(H_plus, "liquid", +1)],
        log_K=-6.352, total_id="CO2", label="co2_first",
    )
    co2_second = EquilibriumReaction(
        stoichiometry=[_e(HCO3_minus, "liquid", -1), _e(CO3_2minus, "liquid", +1), _e(H_plus, "liquid", +1)],
        log_K=-10.329, total_id="CO2", label="co2_second",
    )
    nh3_rxn = EquilibriumReaction(
        stoichiometry=[_e(NH4_plus, "liquid", -1), _e(NH3, "liquid", +1), _e(H_plus, "liquid", +1)],
        log_K=-9.252, total_id="NH3", label="nh3",
    )
    return NRChemicalEquilibriumEngine.from_reactions(
        [water, co2_first, co2_second, nh3_rxn], use_activity=True, activity_model="davies",
    )


@pytest.fixture(scope="module")
def phreeqpython_mod():
    return pytest.importorskip("phreeqpython")


@pytest.fixture(scope="module")
def phreeqc_carbonate_engine(phreeqpython_mod):
    return PHREEQCChemicalEquilibriumEngine(
        {"CO2": 1.0}, component_map={"CO2": "C(4)"}, T_C=25.0, use_warmstart=False,
    )


@pytest.fixture(scope="module")
def phreeqc_carbonate_nh3_engine(phreeqpython_mod):
    return PHREEQCChemicalEquilibriumEngine(
        {"CO2": 1.0, "NH3": 0.5}, component_map={"CO2": "C(4)", "NH3": "N(-3)"}, T_C=25.0, use_warmstart=False,
    )


class TestCarbonateAgreement:
    @pytest.mark.parametrize("CT_CO2", [0.001, 0.010, 0.050, 0.100])
    def test_ph_agrees_with_phreeqc(self, CT_CO2, phreeqc_carbonate_engine):
        nr = _make_nr_carbonate_engine()
        nr_out = nr.solve(totals={"CO2": CT_CO2}, strong_ions={}, T_K=T_K)
        pq_out = phreeqc_carbonate_engine.solve(totals={"CO2": CT_CO2})
        assert abs(float(nr_out.pH) - float(pq_out.pH)) < PH_TOL_CARBONATE


class TestCarbonateNH3Agreement:
    @pytest.mark.parametrize("CT_NH3", [0.0, 0.003, 0.006, 0.010])
    def test_ph_agrees_with_phreeqc(self, CT_NH3, phreeqc_carbonate_nh3_engine):
        nr = _make_nr_carbonate_nh3_engine()
        nr_out = nr.solve(totals={"CO2": 0.010, "NH3": CT_NH3}, strong_ions={}, T_K=T_K)
        pq_out = phreeqc_carbonate_nh3_engine.solve(totals={"CO2": 0.010, "NH3": CT_NH3})
        assert abs(float(nr_out.pH) - float(pq_out.pH)) < PH_TOL_CARBONATE_NH3
