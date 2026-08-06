# -*- coding: utf-8 -*-
"""Tests for src/speciation/protocols.py — Phase 1 of the protocol hierarchy.

Coverage
--------
1. Both existing engines pass isinstance(engine, ChemicalEquilibriumEngineProtocol).
2. Neither existing engine passes isinstance(engine, GrayBoxEngineProtocol).
3. SpeciationJacobian and SparseJacobian are immutable (frozen dataclasses).
4. EquilibriumResult: construction, immutability, apply_to_phases(), to_dict().
5. WhiteBoxEngineProtocol composes all capability protocols correctly.
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pytest
import scipy.sparse

from PyOMES.chemical_equilibrium.protocols import (
    EquilibriumResult,
    GrayBoxEngineProtocol,
    ResidualCapable,
    SolutionJacobianCapable,
    SparseJacobian,
    ChemicalEquilibriumEngineProtocol,
    SpeciationJacobian,
    SplitJacobianCapable,
    WhiteBoxEngineProtocol,
)


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

def _make_reactions():
    """Minimal carbonate + ammonia reaction set (reused from NR engine tests)."""
    from PyOMES.chemistry.common_species import (
        CO2, CO3_2minus, H2O, H_plus, HCO3_minus, NH3, NH4_plus, OH_minus,
    )
    from PyOMES.reactions.equilibrium import EquilibriumReaction
    from PyOMES.reactions.stoichiometry import StoichiometryEntry

    def _e(species, coeff):
        return StoichiometryEntry(species=species, phase="liquid", coefficient=coeff)

    water = EquilibriumReaction(
        stoichiometry=[_e(H2O, -1.0), _e(H_plus, +1.0), _e(OH_minus, +1.0)],
        log_K=-14.0,
        balance_elements=("H", "O"),
        label="water",
    )
    co2_first = EquilibriumReaction(
        stoichiometry=[
            _e(CO2, -1.0), _e(H2O, -1.0), _e(HCO3_minus, +1.0), _e(H_plus, +1.0),
        ],
        log_K=-6.35,
        total_id="CO2",
        balance_elements=("C", "H", "O"),
        label="co2_first",
    )
    co2_second = EquilibriumReaction(
        stoichiometry=[_e(HCO3_minus, -1.0), _e(CO3_2minus, +1.0), _e(H_plus, +1.0)],
        log_K=-10.33,
        total_id="CO2",
        balance_elements=("C", "H", "O"),
        label="co2_second",
    )
    nh4 = EquilibriumReaction(
        stoichiometry=[_e(NH4_plus, -1.0), _e(NH3, +1.0), _e(H_plus, +1.0)],
        log_K=-9.25,
        total_id="NH3",
        balance_elements=("N", "H"),
        label="nh4",
    )
    return [water, co2_first, co2_second, nh4]


@pytest.fixture(scope="module")
def speciation_engine():
    from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
    return BisectionChemicalEquilibriumEngine.from_reactions(_make_reactions())


@pytest.fixture(scope="module")
def nr_engine():
    from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
    return NRChemicalEquilibriumEngine.from_reactions(_make_reactions())


# ---------------------------------------------------------------------------
# 1. Both existing engines satisfy ChemicalEquilibriumEngineProtocol
# ---------------------------------------------------------------------------

class TestChemicalEquilibriumEngineProtocol:
    def test_speciation_engine_isinstance(self, speciation_engine):
        assert isinstance(speciation_engine, ChemicalEquilibriumEngineProtocol)

    def test_nr_engine_isinstance(self, nr_engine):
        assert isinstance(nr_engine, ChemicalEquilibriumEngineProtocol)

    def test_arbitrary_object_is_not(self):
        assert not isinstance(object(), ChemicalEquilibriumEngineProtocol)

    def test_missing_one_method_is_not(self):
        class AlmostEngine:
            n_solve_calls = 0
            def solve(self, **kwargs): ...
            def algebraic_species(self): ...
            def reset_cache(self): ...
            # reset_counters intentionally omitted

        assert not isinstance(AlmostEngine(), ChemicalEquilibriumEngineProtocol)


# ---------------------------------------------------------------------------
# 2. Neither existing engine satisfies GrayBoxEngineProtocol
# ---------------------------------------------------------------------------

class TestGrayBoxProtocol:
    def test_speciation_engine_not_graybox(self, speciation_engine):
        assert not isinstance(speciation_engine, GrayBoxEngineProtocol)

    def test_nr_engine_is_graybox_after_phase4(self, nr_engine):
        # Phase 4 added jacobian_dz_dy() (and the split-Jacobian methods) to
        # NRChemicalEquilibriumEngine so it now structurally satisfies GrayBoxEngineProtocol.
        assert isinstance(nr_engine, GrayBoxEngineProtocol)

    def test_graybox_requires_jacobian_dz_dy(self):
        class FakeGrayBox:
            n_solve_calls = 0
            def solve(self, **kwargs): ...
            def algebraic_species(self): ...
            def reset_cache(self): ...
            def reset_counters(self): ...
            def jacobian_dz_dy(self, **kwargs): ...

        assert isinstance(FakeGrayBox(), GrayBoxEngineProtocol)

    def test_whitebox_implies_graybox(self):
        class FakeWhiteBox:
            n_solve_calls = 0
            def solve(self, **kwargs): ...
            def algebraic_species(self): ...
            def reset_cache(self): ...
            def reset_counters(self): ...
            def jacobian_dz_dy(self, **kwargs): ...
            def residual(self, **kwargs): ...
            def jacobian_dg_dz(self, **kwargs): ...
            def jacobian_dg_dy(self, **kwargs): ...

        engine = FakeWhiteBox()
        assert isinstance(engine, GrayBoxEngineProtocol)
        assert isinstance(engine, WhiteBoxEngineProtocol)


# ---------------------------------------------------------------------------
# 3. SpeciationJacobian and SparseJacobian are frozen (immutable)
# ---------------------------------------------------------------------------

class TestDataContainersImmutability:
    def test_speciation_jacobian_is_frozen(self):
        j = SpeciationJacobian(
            dz_dy=np.eye(2),
            algebraic_ids=("OH-", "HCO3-"),
            component_ids=("CO2", "NH3"),
        )
        with pytest.raises((dataclasses.FrozenInstanceError, AttributeError)):
            j.dz_dy = np.zeros((2, 2))

    def test_speciation_jacobian_field_assignment_frozen(self):
        j = SpeciationJacobian(
            dz_dy=np.zeros((3, 2)),
            algebraic_ids=("H+", "OH-", "HCO3-"),
            component_ids=("CO2", "NH3"),
        )
        with pytest.raises((dataclasses.FrozenInstanceError, AttributeError)):
            j.algebraic_ids = ("x",)

    def test_sparse_jacobian_is_frozen(self):
        m = scipy.sparse.eye(3, format="csr")
        sj = SparseJacobian(
            matrix=m,
            row_ids=("r0", "r1", "r2"),
            col_ids=("c0", "c1", "c2"),
        )
        with pytest.raises((dataclasses.FrozenInstanceError, AttributeError)):
            sj.matrix = scipy.sparse.eye(3, format="csr")

    def test_sparse_jacobian_col_ids_frozen(self):
        m = scipy.sparse.csr_matrix((2, 2))
        sj = SparseJacobian(matrix=m, row_ids=("a", "b"), col_ids=("x", "y"))
        with pytest.raises((dataclasses.FrozenInstanceError, AttributeError)):
            sj.col_ids = ("z", "w")

    def test_speciation_jacobian_stores_values(self):
        arr = np.array([[1.0, 2.0], [3.0, 4.0]])
        j = SpeciationJacobian(
            dz_dy=arr,
            algebraic_ids=("OH-", "HCO3-"),
            component_ids=("CO2", "NH3"),
        )
        np.testing.assert_array_equal(j.dz_dy, arr)
        assert j.algebraic_ids == ("OH-", "HCO3-")
        assert j.component_ids == ("CO2", "NH3")

    def test_sparse_jacobian_stores_values(self):
        m = scipy.sparse.eye(4, format="csr")
        ids = ("a", "b", "c", "d")
        sj = SparseJacobian(matrix=m, row_ids=ids, col_ids=ids)
        assert sj.matrix.shape == (4, 4)
        assert sj.row_ids == ids
        assert sj.col_ids == ids


# ---------------------------------------------------------------------------
# 4. EquilibriumResult: construction, immutability, apply_to_phases(), to_dict()
# ---------------------------------------------------------------------------

class TestEquilibriumResult:
    def test_construction_minimal(self):
        result = EquilibriumResult(pH=7.0)
        assert result.pH == 7.0
        assert result.species_mol_L == {}
        assert result.extra == {}

    def test_construction_full(self):
        result = EquilibriumResult(
            pH=6.8,
            pH_conc=6.8,
            logH=-6.8,
            aH=1.58e-7,
            gamma_H=1.0,
            gamma_OH=1.0,
            ionic_strength=0.01,
            charge_residual=1e-13,
            n_iter=None,
            species_mol_L={"CO2": 0.001, "HCO3-": 0.0005, "H+": 1.58e-7},
            partial_pressures_atm={"CO2": 0.05},
            saturation_indices={"calcite": -0.3},
            alphas={"CO2aq": 0.6},
            extra={"minerals_xi_mol_L": {"calcite": 0.0}},
        )
        assert result.ionic_strength == pytest.approx(0.01)
        assert result.species_mol_L["CO2"] == pytest.approx(0.001)
        assert result.saturation_indices["calcite"] == pytest.approx(-0.3)

    def test_is_frozen(self):
        result = EquilibriumResult(pH=7.0)
        with pytest.raises((dataclasses.FrozenInstanceError, AttributeError)):
            result.pH = 8.0

    def test_species_mol_L_default_is_independent_dict(self):
        # Guards against a shared-mutable-default bug: two instances must
        # not alias the same dict.
        a = EquilibriumResult(pH=7.0)
        b = EquilibriumResult(pH=7.0)
        a.species_mol_L["CO2"] = 1.0
        assert "CO2" not in b.species_mol_L

    def test_to_dict_round_trips_all_fields(self):
        result = EquilibriumResult(
            pH=7.0,
            ionic_strength=0.02,
            species_mol_L={"CO2": 0.001},
        )
        d = result.to_dict()
        assert d["pH"] == 7.0
        assert d["ionic_strength"] == pytest.approx(0.02)
        assert d["species_mol_L"] == {"CO2": 0.001}
        assert set(d.keys()) == {f.name for f in dataclasses.fields(EquilibriumResult)}

    def test_apply_to_phases_writes_n_mol_scaled_by_V_L(self):
        from PyOMES.core.phases import LiquidPhase

        liq = LiquidPhase(n_mol={"H+": 1e-7}, V_L=2.0, T_K=298.15)
        result = EquilibriumResult(
            pH=7.0,
            species_mol_L={"H+": 5e-8, "OH-": 2e-7, "CO2": 0.001},
        )
        result.apply_to_phases({"liquid": liq})

        assert liq.n_mol["H+"] == pytest.approx(1e-7)   # 5e-8 mol/L * 2.0 L
        assert liq.n_mol["OH-"] == pytest.approx(4e-7)  # 2e-7 mol/L * 2.0 L
        assert liq.n_mol["CO2"] == pytest.approx(0.002)  # 0.001 mol/L * 2.0 L

    def test_apply_to_phases_leaves_untouched_species_alone(self):
        from PyOMES.core.phases import LiquidPhase

        liq = LiquidPhase(n_mol={"H+": 1e-7, "Na+": 0.05}, V_L=1.0, T_K=298.15)
        result = EquilibriumResult(pH=7.0, species_mol_L={"H+": 1e-7})
        result.apply_to_phases({"liquid": liq})

        assert liq.n_mol["Na+"] == pytest.approx(0.05)

    def test_apply_to_phases_noop_when_liquid_missing(self):
        result = EquilibriumResult(pH=7.0, species_mol_L={"H+": 1e-7})
        # Must not raise — liquid_key absent from phases dict.
        result.apply_to_phases({})

    def test_apply_to_phases_noop_when_V_L_zero(self):
        from PyOMES.core.phases import GasPhase

        # GasPhase has no _refresh_derived / V_L in the LiquidPhase sense —
        # exercises the hasattr guard rather than the V_L guard.
        gas = GasPhase(n_mol={"O2": 1.0}, V_L=1.0, T_K=298.15)
        result = EquilibriumResult(pH=7.0, species_mol_L={"O2": 5.0})
        # Must not raise even though "liquid" key holds a phase without
        # _refresh_derived.
        result.apply_to_phases({"liquid": gas})
        assert gas.n_mol["O2"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# 5. Capability mixin protocols work independently
# ---------------------------------------------------------------------------

class TestCapabilityMixins:
    def test_residual_capable_protocol(self):
        class HasResidual:
            def residual(self, **kwargs): ...

        assert isinstance(HasResidual(), ResidualCapable)
        assert not isinstance(object(), ResidualCapable)

    def test_solution_jacobian_capable_protocol(self):
        class HasJacobian:
            def jacobian_dz_dy(self, **kwargs): ...

        assert isinstance(HasJacobian(), SolutionJacobianCapable)
        assert not isinstance(object(), SolutionJacobianCapable)

    def test_split_jacobian_capable_protocol(self):
        class HasSplitJacobians:
            def jacobian_dg_dz(self, **kwargs): ...
            def jacobian_dg_dy(self, **kwargs): ...

        assert isinstance(HasSplitJacobians(), SplitJacobianCapable)

    def test_split_jacobian_requires_both_methods(self):
        class OnlyDgDz:
            def jacobian_dg_dz(self, **kwargs): ...

        assert not isinstance(OnlyDgDz(), SplitJacobianCapable)
