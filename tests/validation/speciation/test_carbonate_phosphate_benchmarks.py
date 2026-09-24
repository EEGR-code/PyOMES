# -*- coding: utf-8 -*-
"""Tighter NR-engine benchmarks for 01_single_component_benchmarks.ipynb.

The existing `tests/standalone/test_speciation.py::TestCarbonateSystem` and
`TestPhosphateSystem` only exercise `BisectionChemicalEquilibriumEngine` at
loose tolerances (+/-0.3 pH for carbonate, a bare range check for phosphate).
This file closes the gap the notebook actually demonstrates: NRChemicalEquilibriumEngine
against closed-form analytical formulas, at the notebook's own tight
tolerance (<0.01 pH for carbonate and for H3PO4 itself).

For the amphoteric phosphate salts (NaH2PO4, Na2HPO4, Na3PO4), the notebook
is explicit that NR is the *correct* reference and the textbook
half-sum-of-pKa formula is only a rough approximation (off by up to 1.3 pH
units for Na3PO4) -- so those cases are locked in as regression checks
against the notebook's own printed NR values, not against the formula.
"""
from __future__ import annotations

import numpy as np
import pytest

from PyOMES.chemistry.common_species import (
    H2O, H_plus, OH_minus,
    CO2, HCO3_minus, CO3_2minus,
    H3PO4, H2PO4_minus, HPO4_2minus, PO4_3minus,
)
from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
from PyOMES.chemical_equilibrium.engines.nr.engine import NRChemicalEquilibriumEngine


def _e(species, phase, coeff):
    return StoichiometryEntry(species=species, phase=phase, coefficient=coeff)


def _make_carbonate_engine():
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
    return NRChemicalEquilibriumEngine.from_reactions([water, co2_first, co2_second])


def _make_phosphate_engine():
    water = EquilibriumReaction(
        stoichiometry=[_e(H2O, "liquid", -1), _e(H_plus, "liquid", +1), _e(OH_minus, "liquid", +1)],
        log_K=-14.0, label="water",
    )
    p1 = EquilibriumReaction(
        stoichiometry=[_e(H3PO4, "liquid", -1), _e(H2PO4_minus, "liquid", +1), _e(H_plus, "liquid", +1)],
        log_K=-2.15, total_id="H3PO4", label="p1",
    )
    p2 = EquilibriumReaction(
        stoichiometry=[_e(H2PO4_minus, "liquid", -1), _e(HPO4_2minus, "liquid", +1), _e(H_plus, "liquid", +1)],
        log_K=-7.20, total_id="H3PO4", label="p2",
    )
    p3 = EquilibriumReaction(
        stoichiometry=[_e(HPO4_2minus, "liquid", -1), _e(PO4_3minus, "liquid", +1), _e(H_plus, "liquid", +1)],
        log_K=-12.35, total_id="H3PO4", label="p3",
    )
    return NRChemicalEquilibriumEngine.from_reactions([water, p1, p2, p3])


# =============================================================================
#  Carbonate: NR vs. diprotic-quadratic analytical formula
# =============================================================================

class TestCarbonateAnalyticalAgreement:
    """NR must match the closed-form diprotic pH formula to < 0.01 pH units.

    Formula neglects CO3-- and OH- (valid across the CT range tested, per
    the notebook): CT*Ka1 = H^2 + Ka1*H  =>  H = (-Ka1 + sqrt(Ka1^2+4*Ka1*CT))/2.
    """

    PH_TOL = 0.01
    Ka1 = 10 ** (-6.35)

    @staticmethod
    def _analytical_ph(CT, Ka1):
        disc = Ka1 ** 2 + 4 * Ka1 * CT
        H = (-Ka1 + disc ** 0.5) / 2
        return -np.log10(H)

    @pytest.mark.parametrize("CT", [1e-4, 1e-3, 1e-2, 1e-1, 5e-1])
    def test_ph_matches_analytical(self, CT):
        engine = _make_carbonate_engine()
        out = engine.solve(totals={"CO2": CT}, strong_ions={})
        ph_analytical = self._analytical_ph(CT, self.Ka1)
        assert abs(float(out.pH) - ph_analytical) < self.PH_TOL

    def test_charge_residual_near_zero(self):
        engine = _make_carbonate_engine()
        out = engine.solve(totals={"CO2": 1e-2}, strong_ions={})
        assert abs(float(out.charge_residual)) < 1e-8


# =============================================================================
#  Phosphate: H3PO4 vs. analytical quadratic (tight), amphoteric salts (regression)
# =============================================================================

class TestPhosphateAnalyticalAgreement:
    PH_TOL = 0.01
    Ka1 = 10 ** (-2.15)

    def test_h3po4_matches_analytical(self):
        CT_p = 0.01
        disc = self.Ka1 ** 2 + 4 * self.Ka1 * CT_p
        H = (-self.Ka1 + disc ** 0.5) / 2
        ph_analytical = -np.log10(H)

        engine = _make_phosphate_engine()
        out = engine.solve(totals={"H3PO4": CT_p}, strong_ions={})
        assert abs(float(out.pH) - ph_analytical) < self.PH_TOL


class TestPhosphateAmphotericRegression:
    """NR is the correct reference for amphoteric salts (notebook 01, Section 1.2) --
    the half-sum-of-pKa approximation diverges from it by up to 1.3 pH units.
    These lock in the notebook's own printed NR values rather than the formula.
    """

    CT_p = 0.01
    PH_TOL = 0.01

    @pytest.mark.parametrize(
        "CT_Na, expected_ph",
        [
            (0.01, 4.7918),   # NaH2PO4
            (0.02, 9.5189),   # Na2HPO4
            (0.03, 11.8746),  # Na3PO4
        ],
    )
    def test_amphoteric_salt_ph_regression(self, CT_Na, expected_ph):
        engine = _make_phosphate_engine()
        out = engine.solve(totals={"H3PO4": self.CT_p}, strong_ions={"CT_Na": CT_Na})
        assert abs(float(out.pH) - expected_ph) < self.PH_TOL

    @pytest.mark.parametrize("CT_Na", [0.0, 0.01, 0.02, 0.03])
    def test_total_phosphorus_conserved(self, CT_Na):
        engine = _make_phosphate_engine()
        out = engine.solve(totals={"H3PO4": self.CT_p}, strong_ions={"CT_Na": CT_Na})
        total_P = sum(out.species_mol_L[sp] for sp in ("H3PO4", "H2PO4-", "HPO4--", "PO4---"))
        assert total_P == pytest.approx(self.CT_p, rel=1e-6)
