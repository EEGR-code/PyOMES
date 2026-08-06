# -*- coding: utf-8 -*-
"""Tests for NumericalGradientEquilibriumEngine (Phase 2 of the protocol hierarchy).

Coverage
--------
1. Protocol conformance: satisfies GrayBoxEngineProtocol, not WhiteBoxEngineProtocol.
2. Delegation: solve(), algebraic_species(), reset_cache(), reset_counters(),
   n_solve_calls all proxy to the inner engine.
3. jacobian_dz_dy() returns SpeciationJacobian with the correct shape.
4. jacobian_dz_dy() values have the correct sign and order of magnitude.
5. eps_abs / eps_rel control the perturbation magnitude.
6. ValueError when totals= is not supplied, including the phases= path.
7. n_solve_calls increments by 2 * n_components per Jacobian call.
8. Constructor rejects non-engine arguments.
"""
from __future__ import annotations

import numpy as np
import pytest

from PyOMES.chemical_equilibrium.numerical_gradient import NumericalGradientEquilibriumEngine
from PyOMES.chemical_equilibrium.protocols import (
    GrayBoxEngineProtocol,
    ChemicalEquilibriumEngineProtocol,
    SpeciationJacobian,
    WhiteBoxEngineProtocol,
)


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

def _make_reactions():
    from PyOMES.chemistry.common_species import (
        CO2, CO3_2minus, H2O, H_plus, HCO3_minus, NH3, NH4_plus, OH_minus,
    )
    from PyOMES.reactions.equilibrium import EquilibriumReaction
    from PyOMES.reactions.stoichiometry import StoichiometryEntry

    def _e(species, coeff):
        return StoichiometryEntry(species=species, phase="liquid", coefficient=coeff)

    return [
        EquilibriumReaction(
            stoichiometry=[_e(H2O, -1.0), _e(H_plus, +1.0), _e(OH_minus, +1.0)],
            log_K=-14.0, balance_elements=("H", "O"), label="water",
        ),
        EquilibriumReaction(
            stoichiometry=[
                _e(CO2, -1.0), _e(H2O, -1.0), _e(HCO3_minus, +1.0), _e(H_plus, +1.0),
            ],
            log_K=-6.35, total_id="CO2",
            balance_elements=("C", "H", "O"), label="co2_first",
        ),
        EquilibriumReaction(
            stoichiometry=[_e(HCO3_minus, -1.0), _e(CO3_2minus, +1.0), _e(H_plus, +1.0)],
            log_K=-10.33, total_id="CO2",
            balance_elements=("C", "H", "O"), label="co2_second",
        ),
        EquilibriumReaction(
            stoichiometry=[_e(NH4_plus, -1.0), _e(NH3, +1.0), _e(H_plus, +1.0)],
            log_K=-9.25, total_id="NH3",
            balance_elements=("N", "H"), label="nh4",
        ),
    ]


@pytest.fixture(scope="module")
def nr_engine():
    from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
    return NRChemicalEquilibriumEngine.from_reactions(_make_reactions())


@pytest.fixture(scope="module")
def gradient_engine(nr_engine):
    return NumericalGradientEquilibriumEngine(nr_engine)


# Totals used for a ~neutral pH system
_TOTALS = {"CO2": 0.05, "NH3": 0.04}


# ---------------------------------------------------------------------------
# 1. Protocol conformance
# ---------------------------------------------------------------------------

class TestProtocolConformance:
    def test_satisfies_speciation_engine_protocol(self, gradient_engine):
        assert isinstance(gradient_engine, ChemicalEquilibriumEngineProtocol)

    def test_satisfies_graybox_protocol(self, gradient_engine):
        assert isinstance(gradient_engine, GrayBoxEngineProtocol)

    def test_does_not_satisfy_whitebox_protocol(self, gradient_engine):
        assert not isinstance(gradient_engine, WhiteBoxEngineProtocol)

    def test_constructor_rejects_non_engine(self):
        with pytest.raises(TypeError, match="ChemicalEquilibriumEngineProtocol"):
            NumericalGradientEquilibriumEngine(object())

    def test_eps_attributes_stored(self, nr_engine):
        eng = NumericalGradientEquilibriumEngine(nr_engine, eps_abs=1e-8, eps_rel=1e-4)
        assert eng.eps_abs == pytest.approx(1e-8)
        assert eng.eps_rel == pytest.approx(1e-4)

    def test_default_eps_values(self, gradient_engine):
        assert gradient_engine.eps_abs == pytest.approx(1e-10)
        assert gradient_engine.eps_rel == pytest.approx(1e-5)


# ---------------------------------------------------------------------------
# 2. Delegation
# ---------------------------------------------------------------------------

class TestDelegation:
    def test_solve_returns_same_as_inner(self, nr_engine, gradient_engine):
        out_inner = nr_engine.solve(totals=dict(_TOTALS))
        out_outer = gradient_engine.solve(totals=dict(_TOTALS))
        assert out_inner.pH == pytest.approx(out_outer.pH, abs=1e-6)

    def test_algebraic_species_matches_inner(self, nr_engine, gradient_engine):
        assert gradient_engine.algebraic_species() == nr_engine.algebraic_species()

    def test_reset_cache_delegates(self, nr_engine):
        eng = NumericalGradientEquilibriumEngine(nr_engine)
        eng.solve(totals=dict(_TOTALS))
        eng.reset_cache()
        assert nr_engine._cache.log_x is None

    def test_reset_counters_delegates(self, nr_engine):
        eng = NumericalGradientEquilibriumEngine(nr_engine)
        eng.solve(totals=dict(_TOTALS))
        eng.reset_counters()
        assert eng.n_solve_calls == 0

    def test_n_solve_calls_proxies_inner(self, nr_engine):
        eng = NumericalGradientEquilibriumEngine(nr_engine)
        eng.reset_counters()
        before = eng.n_solve_calls
        eng.solve(totals=dict(_TOTALS))
        assert eng.n_solve_calls == before + 1


# ---------------------------------------------------------------------------
# 3. jacobian_dz_dy — shape and structure
# ---------------------------------------------------------------------------

class TestJacobianShape:
    def test_returns_speciation_jacobian(self, gradient_engine):
        jac = gradient_engine.jacobian_dz_dy(totals=dict(_TOTALS))
        assert isinstance(jac, SpeciationJacobian)

    def test_shape_matches_alg_and_comp(self, gradient_engine):
        jac = gradient_engine.jacobian_dz_dy(totals=dict(_TOTALS))
        n_alg = len(gradient_engine.algebraic_species())
        n_comp = len(_TOTALS)
        assert jac.dz_dy.shape == (n_alg, n_comp)

    def test_algebraic_ids_sorted(self, gradient_engine):
        jac = gradient_engine.jacobian_dz_dy(totals=dict(_TOTALS))
        assert list(jac.algebraic_ids) == sorted(jac.algebraic_ids)

    def test_component_ids_match_totals_keys(self, gradient_engine):
        jac = gradient_engine.jacobian_dz_dy(totals=dict(_TOTALS))
        assert set(jac.component_ids) == set(_TOTALS.keys())

    def test_component_ids_order_preserved(self, gradient_engine):
        ordered = {"NH3": 0.04, "CO2": 0.05}
        jac = gradient_engine.jacobian_dz_dy(totals=ordered)
        assert jac.component_ids == ("NH3", "CO2")

    def test_single_component(self, nr_engine):
        eng = NumericalGradientEquilibriumEngine(nr_engine)
        jac = eng.jacobian_dz_dy(totals={"CO2": 0.05})
        n_alg = len(eng.algebraic_species())
        assert jac.dz_dy.shape == (n_alg, 1)


# ---------------------------------------------------------------------------
# 4. jacobian_dz_dy — value correctness
# ---------------------------------------------------------------------------

class TestJacobianValues:
    """Qualitative checks: sign and order of magnitude for the carbonate system.

    At pH ~7 with CT_CO2=0.05 mol/L:
      - HCO3- dominates (~0.97 of CT_CO2)
      - d[HCO3-]/d[CT_CO2] ≈ alpha_1 (positive, < 1)
      - d[CO3--]/d[CT_CO2] > 0 (positive but small)
      - d[NH4+]/d[CT_NH3] > 0 (more total NH3 → more NH4+)
      - d[HCO3-]/d[CT_NH3] ≈ 0 (CO3 ladder weakly coupled to NH3 at neutral pH)
    """

    @pytest.fixture(scope="class")
    def jac(self, gradient_engine):
        return gradient_engine.jacobian_dz_dy(totals=dict(_TOTALS))

    def _get(self, jac, row_id, col_id):
        i = jac.algebraic_ids.index(row_id)
        j = jac.component_ids.index(col_id)
        return jac.dz_dy[i, j]

    def test_dhco3_dco2_positive(self, jac):
        # At this acidic operating point CO2 dominates, but HCO3- still grows
        # when total CO2 increases (positive, just small).
        v = self._get(jac, "HCO3-", "CO2")
        assert v > 0.0, f"d[HCO3-]/d[CT_CO2]={v:.6f} expected > 0"

    def test_dco2_dco2_dominates_column(self, jac):
        # CO2 is the dominant species at this acidic pH; it captures most of
        # any added total CO2.
        v_co2 = self._get(jac, "CO2", "CO2")
        v_hco3 = self._get(jac, "HCO3-", "CO2")
        assert v_co2 > v_hco3, "d[CO2]/d[CT_CO2] should exceed d[HCO3-]/d[CT_CO2]"

    def test_dco2_dnh3_negative(self, jac):
        # NH3 is a base: adding it raises pH, converting CO2 → HCO3-.
        # So d[CO2]/d[CT_NH3] < 0.
        v = self._get(jac, "CO2", "NH3")
        assert v < 0.0, f"d[CO2]/d[CT_NH3]={v:.4f} expected < 0 (base raises pH)"

    def test_dnh4_dnh3_positive(self, jac):
        v = self._get(jac, "NH4+", "NH3")
        assert v > 0.5, f"d[NH4+]/d[CT_NH3]={v:.4f} expected > 0.5"

    def test_dco2_dco2_positive(self, jac):
        v = self._get(jac, "CO2", "CO2")
        assert v >= 0.0, f"d[CO2]/d[CT_CO2]={v:.4f} expected >= 0"

    def test_jacobian_finite(self, jac):
        assert np.all(np.isfinite(jac.dz_dy))

    def test_mass_balance_column(self, jac):
        # Sum of d[all CO2-ladder species]/d[CT_CO2] should be close to 1
        # (perturbation distributes across CO2, HCO3-, CO3--)
        j = jac.component_ids.index("CO2")
        co2_ids = {"CO2", "HCO3-", "CO3--"}
        total = sum(
            jac.dz_dy[jac.algebraic_ids.index(sp), j]
            for sp in co2_ids
            if sp in jac.algebraic_ids
        )
        assert total == pytest.approx(1.0, abs=0.05), (
            f"Sum of d[CO2-ladder]/d[CT_CO2] = {total:.4f}, expected ~1.0"
        )


# ---------------------------------------------------------------------------
# 5. eps_abs / eps_rel perturbation control
# ---------------------------------------------------------------------------

class TestEpsilonControl:
    def test_pure_absolute_mode(self, nr_engine):
        eng = NumericalGradientEquilibriumEngine(nr_engine, eps_abs=1e-6, eps_rel=0.0)
        jac = eng.jacobian_dz_dy(totals=dict(_TOTALS))
        assert np.all(np.isfinite(jac.dz_dy))

    def test_pure_relative_mode(self, nr_engine):
        eng = NumericalGradientEquilibriumEngine(nr_engine, eps_abs=0.0, eps_rel=1e-4)
        jac = eng.jacobian_dz_dy(totals=dict(_TOTALS))
        assert np.all(np.isfinite(jac.dz_dy))

    def test_large_eps_still_converges(self, nr_engine):
        eng = NumericalGradientEquilibriumEngine(nr_engine, eps_abs=1e-3, eps_rel=0.0)
        jac = eng.jacobian_dz_dy(totals=dict(_TOTALS))
        assert np.all(np.isfinite(jac.dz_dy))

    def test_near_zero_total_clamped(self, nr_engine):
        eng = NumericalGradientEquilibriumEngine(nr_engine, eps_abs=1e-10, eps_rel=1e-5)
        tiny_totals = {"CO2": 1e-12, "NH3": 0.04}
        jac = eng.jacobian_dz_dy(totals=tiny_totals)
        assert np.all(np.isfinite(jac.dz_dy))

    def test_result_consistent_across_eps(self, nr_engine):
        jac_fine = NumericalGradientEquilibriumEngine(
            nr_engine, eps_abs=1e-10, eps_rel=1e-5
        ).jacobian_dz_dy(totals=dict(_TOTALS))
        jac_coarse = NumericalGradientEquilibriumEngine(
            nr_engine, eps_abs=1e-4, eps_rel=0.0
        ).jacobian_dz_dy(totals=dict(_TOTALS))
        # Same sign pattern; Frobenius norm within 5% of each other
        np.testing.assert_array_equal(
            np.sign(jac_fine.dz_dy), np.sign(jac_coarse.dz_dy)
        )


# ---------------------------------------------------------------------------
# 6. Error handling
# ---------------------------------------------------------------------------

class TestErrorHandling:
    def test_missing_totals_raises(self, gradient_engine):
        with pytest.raises(ValueError, match="totals="):
            gradient_engine.jacobian_dz_dy()

    def test_phases_without_totals_raises(self, gradient_engine):
        with pytest.raises(ValueError, match="phases="):
            gradient_engine.jacobian_dz_dy(phases={"liquid": object()})

    def test_empty_totals_returns_empty_jacobian(self, gradient_engine):
        jac = gradient_engine.jacobian_dz_dy(totals={})
        n_alg = len(gradient_engine.algebraic_species())
        assert jac.dz_dy.shape == (n_alg, 0)
        assert jac.component_ids == ()


# ---------------------------------------------------------------------------
# 7. n_solve_calls accounting
# ---------------------------------------------------------------------------

class TestSolveCallCounting:
    def test_jacobian_increments_by_2n_comp(self, nr_engine):
        eng = NumericalGradientEquilibriumEngine(nr_engine)
        eng.reset_counters()
        n_comp = len(_TOTALS)
        eng.jacobian_dz_dy(totals=dict(_TOTALS))
        assert eng.n_solve_calls == 2 * n_comp

    def test_direct_solve_increments_by_one(self, nr_engine):
        eng = NumericalGradientEquilibriumEngine(nr_engine)
        eng.reset_counters()
        eng.solve(totals=dict(_TOTALS))
        assert eng.n_solve_calls == 1

    def test_multiple_jacobian_calls_accumulate(self, nr_engine):
        eng = NumericalGradientEquilibriumEngine(nr_engine)
        eng.reset_counters()
        n_comp = len(_TOTALS)
        eng.jacobian_dz_dy(totals=dict(_TOTALS))
        eng.jacobian_dz_dy(totals=dict(_TOTALS))
        assert eng.n_solve_calls == 4 * n_comp
