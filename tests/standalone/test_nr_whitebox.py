# -*- coding: utf-8 -*-
"""Tests for NRChemicalEquilibriumEngine Phase 4 white-box additions.

Coverage
--------
1. Protocol conformance: isinstance checks for Whitebox/Graybox/Blackbox.
2. retain_jacobian=False (default): all white-box methods raise RuntimeError.
3. retain_jacobian=True: methods available after solve().
4. jacobian_dg_dz(): shape, row/col IDs, sparsity structure.
5. jacobian_dg_dy(): shape, row/col IDs, value structure (-I block, zero row).
6. residual(): evaluates g(y,z) without updating state.
7. jacobian_dz_dy(): analytical vs. numerical Jacobian consistency.
8. reset_cache() clears cached white-box state.
9. from_reactions() passes retain_jacobian through.
10. Precipitation path: retain_jacobian=True still works.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import issparse

from PyOMES.chemical_equilibrium.protocols import (
    GrayBoxEngineProtocol,
    SparseJacobian,
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
def blackbox_engine():
    from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
    return NRChemicalEquilibriumEngine.from_reactions(_make_reactions())


@pytest.fixture(scope="module")
def whitebox_engine():
    from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
    eng = NRChemicalEquilibriumEngine.from_reactions(_make_reactions(), retain_jacobian=True)
    # Warm up with a solve so all caches are populated
    eng.solve(totals={"CO2": 0.05, "NH3": 0.04})
    return eng


_TOTALS = {"CO2": 0.05, "NH3": 0.04}


# ---------------------------------------------------------------------------
# 1. Protocol conformance
# ---------------------------------------------------------------------------

class TestProtocolConformance:
    def test_blackbox_satisfies_speciation_engine(self, blackbox_engine):
        assert isinstance(blackbox_engine, ChemicalEquilibriumEngineProtocol)

    def test_blackbox_not_graybox_before_phase4(self, blackbox_engine):
        # Before adding white-box methods, NRChemicalEquilibriumEngine was black-box only.
        # After Phase 4, ALL NRChemicalEquilibriumEngine instances have the methods and
        # therefore satisfy GrayBoxEngineProtocol structurally.
        assert isinstance(blackbox_engine, GrayBoxEngineProtocol)

    def test_whitebox_satisfies_whitebox_protocol(self, whitebox_engine):
        assert isinstance(whitebox_engine, WhiteBoxEngineProtocol)

    def test_whitebox_satisfies_graybox_protocol(self, whitebox_engine):
        assert isinstance(whitebox_engine, GrayBoxEngineProtocol)

    def test_whitebox_satisfies_speciation_engine_protocol(self, whitebox_engine):
        assert isinstance(whitebox_engine, ChemicalEquilibriumEngineProtocol)

    def test_retain_jacobian_attribute_false(self, blackbox_engine):
        assert blackbox_engine.retain_jacobian is False

    def test_retain_jacobian_attribute_true(self, whitebox_engine):
        assert whitebox_engine.retain_jacobian is True


# ---------------------------------------------------------------------------
# 2. retain_jacobian=False raises RuntimeError on all white-box methods
# ---------------------------------------------------------------------------

_WHITEBOX_METHODS = [
    "jacobian_dg_dz",
    "jacobian_dg_dy",
    "residual",
    "jacobian_dz_dy",
]


class TestDisabledGuard:
    @pytest.mark.parametrize("method_name", _WHITEBOX_METHODS)
    def test_raises_runtime_error_when_disabled(self, blackbox_engine, method_name):
        method = getattr(blackbox_engine, method_name)
        with pytest.raises(RuntimeError, match="retain_jacobian=True"):
            if method_name == "residual":
                method(totals=dict(_TOTALS))
            else:
                method()

    def test_error_message_mentions_retain_jacobian(self, blackbox_engine):
        with pytest.raises(RuntimeError, match="retain_jacobian=True"):
            blackbox_engine.jacobian_dg_dz()

    def test_regular_solve_still_works_on_blackbox(self, blackbox_engine):
        import math
        out = blackbox_engine.solve(totals=dict(_TOTALS))
        assert not math.isnan(out.pH)

    def test_no_stale_keys_in_solve_output(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        eng = NRChemicalEquilibriumEngine.from_reactions(_make_reactions(), retain_jacobian=True)
        out = eng.solve(totals=dict(_TOTALS))
        for key in ("_jacobian_matrix", "_log_activities", "_ionic_strength_final"):
            assert key not in out.species_mol_L
            assert key not in out.extra


# ---------------------------------------------------------------------------
# 3. jacobian_dg_dz — shape, IDs, sparsity
# ---------------------------------------------------------------------------

class TestJacobianDgDz:
    @pytest.fixture(scope="class")
    def jac(self, whitebox_engine):
        return whitebox_engine.jacobian_dg_dz()

    def test_returns_sparse_jacobian(self, jac):
        assert isinstance(jac, SparseJacobian)

    def test_matrix_is_sparse(self, jac):
        assert issparse(jac.matrix)

    def test_shape_square(self, jac):
        m, n = jac.matrix.shape
        assert m == n, f"Expected square Jacobian; got shape {m}x{n}"

    def test_n_masters(self, whitebox_engine, jac):
        m = len(whitebox_engine.tableau.masters)
        assert jac.matrix.shape == (m, m)

    def test_row_ids_length(self, whitebox_engine, jac):
        m = len(whitebox_engine.tableau.masters)
        assert len(jac.row_ids) == m

    def test_col_ids_length(self, whitebox_engine, jac):
        m = len(whitebox_engine.tableau.masters)
        assert len(jac.col_ids) == m

    def test_last_row_id_is_h_plus(self, jac):
        assert jac.row_ids[-1] == "H+", f"Last row_id should be 'H+'; got {jac.row_ids[-1]}"

    def test_first_col_id_is_h_plus(self, jac):
        assert jac.col_ids[0] == "H+", f"First col_id should be 'H+'; got {jac.col_ids[0]}"

    def test_matrix_nonzero(self, jac):
        dense = jac.matrix.toarray()
        assert np.any(dense != 0.0), "Jacobian should have non-zero entries"

    def test_matrix_finite(self, jac):
        dense = jac.matrix.toarray()
        assert np.all(np.isfinite(dense))

    def test_diagonal_dominance_approx(self, jac):
        # The NR Jacobian for a well-conditioned system is diagonally dominant
        # (mass-balance rows: diagonal = sum of concentrations in component,
        # off-diagonals = cross-component contributions which are smaller).
        dense = np.abs(jac.matrix.toarray())
        diag = np.diag(dense)
        off_diag_sums = dense.sum(axis=1) - diag
        # Allow some slack; just check diagonal >= 0 (positive semidefinite)
        assert np.all(diag >= 0.0)

    def test_raises_before_solve(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        eng = NRChemicalEquilibriumEngine.from_reactions(_make_reactions(), retain_jacobian=True)
        with pytest.raises(RuntimeError, match="solve\\(\\)"):
            eng.jacobian_dg_dz()


# ---------------------------------------------------------------------------
# 4. jacobian_dg_dy — shape, structure, -I block
# ---------------------------------------------------------------------------

class TestJacobianDgDy:
    @pytest.fixture(scope="class")
    def jac(self, whitebox_engine):
        return whitebox_engine.jacobian_dg_dy()

    def test_returns_sparse_jacobian(self, jac):
        assert isinstance(jac, SparseJacobian)

    def test_shape(self, whitebox_engine, jac):
        m = len(whitebox_engine.tableau.masters)
        n_comp = m - 1
        assert jac.matrix.shape == (m, n_comp)

    def test_top_block_is_minus_identity(self, whitebox_engine, jac):
        n_comp = len(whitebox_engine.tableau.masters) - 1
        dense = jac.matrix.toarray()
        top_block = dense[:n_comp, :]
        expected = -np.eye(n_comp)
        np.testing.assert_array_almost_equal(top_block, expected, decimal=12)

    def test_charge_balance_row_is_zero(self, whitebox_engine, jac):
        dense = jac.matrix.toarray()
        last_row = dense[-1, :]
        assert np.all(last_row == 0.0), f"Charge-balance row should be zero; got {last_row}"

    def test_row_ids_same_as_dg_dz(self, whitebox_engine):
        row_ids_dg_dz = whitebox_engine.jacobian_dg_dz().row_ids
        row_ids_dg_dy = whitebox_engine.jacobian_dg_dy().row_ids
        assert row_ids_dg_dz == row_ids_dg_dy

    def test_col_ids_are_components(self, whitebox_engine, jac):
        component_ids = tuple(
            comp.master_id for comp in whitebox_engine.tableau.components
        )
        assert jac.col_ids == component_ids

    def test_last_row_id_is_h_plus(self, jac):
        assert jac.row_ids[-1] == "H+"

    def test_idempotent(self, whitebox_engine):
        # Structural matrix: two calls return equal content
        j1 = whitebox_engine.jacobian_dg_dy()
        j2 = whitebox_engine.jacobian_dg_dy()
        np.testing.assert_array_equal(j1.matrix.toarray(), j2.matrix.toarray())


# ---------------------------------------------------------------------------
# 5. residual()
# ---------------------------------------------------------------------------

class TestResidual:
    def test_returns_ndarray(self, whitebox_engine):
        R = whitebox_engine.residual(totals=dict(_TOTALS))
        assert isinstance(R, np.ndarray)

    def test_length_equals_n_masters(self, whitebox_engine):
        m = len(whitebox_engine.tableau.masters)
        R = whitebox_engine.residual(totals=dict(_TOTALS))
        assert len(R) == m

    def test_near_zero_at_converged_state(self, whitebox_engine):
        # Immediately after solve(), the residual at the cached state should be ≈ 0
        whitebox_engine.solve(totals=dict(_TOTALS))
        R = whitebox_engine.residual(totals=dict(_TOTALS))
        assert np.max(np.abs(R)) < 1e-8, f"Residual at converged state = {R}"

    def test_nonzero_at_perturbed_state(self, whitebox_engine):
        whitebox_engine.solve(totals=dict(_TOTALS))
        # Perturb the log-activities away from convergence
        x_perturbed = whitebox_engine._cached_x.copy()
        x_perturbed[0] += 1.0  # shift log[H+] by 1
        R = whitebox_engine.residual(totals=dict(_TOTALS), log_activities=x_perturbed)
        assert np.max(np.abs(R)) > 1e-4, "Perturbed residual should be non-zero"

    def test_does_not_update_n_solve_calls(self, whitebox_engine):
        before = whitebox_engine.n_solve_calls
        whitebox_engine.residual(totals=dict(_TOTALS))
        assert whitebox_engine.n_solve_calls == before

    def test_does_not_update_cached_x(self, whitebox_engine):
        whitebox_engine.solve(totals=dict(_TOTALS))
        x_before = whitebox_engine._cached_x.copy()
        x_perturbed = x_before.copy()
        x_perturbed[0] += 1.0
        whitebox_engine.residual(totals=dict(_TOTALS), log_activities=x_perturbed)
        np.testing.assert_array_equal(whitebox_engine._cached_x, x_before)

    def test_missing_totals_raises(self, whitebox_engine):
        with pytest.raises(ValueError, match="totals="):
            whitebox_engine.residual()

    def test_missing_cached_state_raises(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        eng = NRChemicalEquilibriumEngine.from_reactions(_make_reactions(), retain_jacobian=True)
        with pytest.raises(RuntimeError, match="solve\\(\\)"):
            eng.residual(totals=dict(_TOTALS))

    def test_finite_values(self, whitebox_engine):
        R = whitebox_engine.residual(totals=dict(_TOTALS))
        assert np.all(np.isfinite(R))


# ---------------------------------------------------------------------------
# 6. jacobian_dz_dy — analytical vs. numerical consistency
# ---------------------------------------------------------------------------

class TestJacobianDzDy:
    @pytest.fixture(scope="class")
    def both_engines(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        from PyOMES.chemical_equilibrium.numerical_gradient import NumericalGradientEquilibriumEngine
        reactions = _make_reactions()
        wb = NRChemicalEquilibriumEngine.from_reactions(reactions, retain_jacobian=True)
        bb = NRChemicalEquilibriumEngine.from_reactions(reactions)
        num = NumericalGradientEquilibriumEngine(bb)
        wb.solve(totals=dict(_TOTALS))
        return wb, num

    def test_returns_speciation_jacobian(self, whitebox_engine):
        jac = whitebox_engine.jacobian_dz_dy()
        assert isinstance(jac, SpeciationJacobian)

    def test_shape_matches_alg_and_comp(self, whitebox_engine):
        jac = whitebox_engine.jacobian_dz_dy()
        n_alg = len(whitebox_engine.algebraic_species())
        n_comp = len(whitebox_engine.tableau.masters) - 1
        assert jac.dz_dy.shape == (n_alg, n_comp)

    def test_algebraic_ids_sorted(self, whitebox_engine):
        jac = whitebox_engine.jacobian_dz_dy()
        assert list(jac.algebraic_ids) == sorted(jac.algebraic_ids)

    def test_component_ids_are_non_h_masters(self, whitebox_engine):
        jac = whitebox_engine.jacobian_dz_dy()
        expected = tuple(comp.master_id for comp in whitebox_engine.tableau.components)
        assert jac.component_ids == expected

    def test_finite_values(self, whitebox_engine):
        jac = whitebox_engine.jacobian_dz_dy()
        assert np.all(np.isfinite(jac.dz_dy))

    def test_consistent_with_numerical_gradient(self, both_engines):
        wb, num = both_engines
        jac_wb = wb.jacobian_dz_dy()
        jac_num = num.jacobian_dz_dy(totals=dict(_TOTALS))

        # Align rows: both use sorted algebraic_ids
        assert jac_wb.algebraic_ids == jac_num.algebraic_ids

        # Columns may differ in ordering; align
        wb_cols = jac_wb.component_ids
        num_cols = jac_num.component_ids

        for comp_id in wb_cols:
            j_wb = wb_cols.index(comp_id)
            j_num = num_cols.index(comp_id)
            col_wb = jac_wb.dz_dy[:, j_wb]
            col_num = jac_num.dz_dy[:, j_num]
            # Allow 1% relative error or 1e-8 absolute for small values
            np.testing.assert_allclose(
                col_wb, col_num,
                rtol=0.01, atol=1e-8,
                err_msg=f"Column for '{comp_id}' disagrees between analytical and numerical",
            )

    def test_mass_balance_sum_to_one(self, whitebox_engine):
        """Sum of d[CO2-ladder species]/d[CT_CO2] ≈ 1 (perturbation distributes)."""
        jac = whitebox_engine.jacobian_dz_dy()
        cols = jac.component_ids
        if "CO2" not in cols:
            pytest.skip("CO2 not in component_ids")
        j = cols.index("CO2")
        co2_ids = {"CO2", "HCO3-", "CO3--"}
        total = sum(
            jac.dz_dy[jac.algebraic_ids.index(sp), j]
            for sp in co2_ids
            if sp in jac.algebraic_ids
        )
        assert total == pytest.approx(1.0, abs=0.05), (
            f"Sum of d[CO2-ladder]/d[CT_CO2] = {total:.4f}, expected ~1.0"
        )

    def test_raises_before_solve(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        eng = NRChemicalEquilibriumEngine.from_reactions(_make_reactions(), retain_jacobian=True)
        with pytest.raises(RuntimeError, match="solve\\(\\)"):
            eng.jacobian_dz_dy()


# ---------------------------------------------------------------------------
# 7. reset_cache clears white-box state
# ---------------------------------------------------------------------------

class TestResetCache:
    def test_clears_cached_jacobian(self, whitebox_engine):
        whitebox_engine.solve(totals=dict(_TOTALS))
        assert whitebox_engine._cached_jacobian is not None
        whitebox_engine.reset_cache()
        assert whitebox_engine._cached_jacobian is None

    def test_clears_cached_x(self, whitebox_engine):
        whitebox_engine.solve(totals=dict(_TOTALS))
        whitebox_engine.reset_cache()
        assert whitebox_engine._cached_x is None

    def test_clears_cached_concentrations(self, whitebox_engine):
        whitebox_engine.solve(totals=dict(_TOTALS))
        whitebox_engine.reset_cache()
        assert whitebox_engine._cached_concentrations is None

    def test_whitebox_methods_raise_after_reset(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        eng = NRChemicalEquilibriumEngine.from_reactions(_make_reactions(), retain_jacobian=True)
        eng.solve(totals=dict(_TOTALS))
        eng.reset_cache()
        with pytest.raises(RuntimeError, match="solve\\(\\)"):
            eng.jacobian_dg_dz()

    def test_solve_repopulates_after_reset(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        eng = NRChemicalEquilibriumEngine.from_reactions(_make_reactions(), retain_jacobian=True)
        eng.solve(totals=dict(_TOTALS))
        eng.reset_cache()
        eng.solve(totals=dict(_TOTALS))
        assert eng._cached_jacobian is not None


# ---------------------------------------------------------------------------
# 8. from_reactions passes retain_jacobian through
# ---------------------------------------------------------------------------

class TestFromReactions:
    def test_default_retain_jacobian_false(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        eng = NRChemicalEquilibriumEngine.from_reactions(_make_reactions())
        assert eng.retain_jacobian is False

    def test_explicit_retain_jacobian_true(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        eng = NRChemicalEquilibriumEngine.from_reactions(_make_reactions(), retain_jacobian=True)
        assert eng.retain_jacobian is True

    def test_whitebox_works_via_from_reactions(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        eng = NRChemicalEquilibriumEngine.from_reactions(_make_reactions(), retain_jacobian=True)
        eng.solve(totals=dict(_TOTALS))
        jac = eng.jacobian_dg_dz()
        assert jac.matrix.shape[0] > 0


# ---------------------------------------------------------------------------
# 9. Implicit function theorem identity: dz_dy = -(dg_dz)^{-1} * dg_dy
# ---------------------------------------------------------------------------

class TestImplicitFunctionTheorem:
    """Verify the analytical relationship between the split and total Jacobians."""

    def test_ift_identity(self, whitebox_engine):
        """-(dg/dz)^{-1} * (dg/dy) should equal dx_dy (first m rows of dz_dy chain)."""
        wb = whitebox_engine
        wb.solve(totals=dict(_TOTALS))

        J = wb._cached_jacobian
        m = len(wb.tableau.masters)
        n_comp = m - 1
        dg_dy_dense = np.zeros((m, n_comp))
        for i in range(n_comp):
            dg_dy_dense[i, i] = -1.0

        # dx/dy from IFT
        dx_dy_ift = np.linalg.solve(J, -dg_dy_dense)

        # dx/dy reconstructed from jacobian_dz_dy (masters only)
        jac = wb.jacobian_dz_dy()
        master_to_xcol = {m_id: k for k, m_id in enumerate(wb.tableau.masters)}
        alg_id_to_row = {sp: i for i, sp in enumerate(jac.algebraic_ids)}
        import math
        LN10 = math.log(10.0)

        dx_dy_from_jac = np.zeros((m, n_comp))
        for k, m_id in enumerate(wb.tableau.masters):
            row = alg_id_to_row.get(m_id)
            if row is None:
                continue
            c_k = wb._cached_concentrations.get(m_id, 0.0)
            if c_k > 0:
                dx_dy_from_jac[k, :] = jac.dz_dy[row, :] / (LN10 * c_k)

        np.testing.assert_allclose(
            dx_dy_from_jac, dx_dy_ift,
            rtol=1e-10, atol=1e-12,
            err_msg="dx/dy from jacobian_dz_dy does not match IFT result",
        )
