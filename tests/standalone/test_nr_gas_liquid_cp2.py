# -*- coding: utf-8 -*-
"""Tests for CP2 of LAYER1_GAP_CLOSURE: mass-action row residual/Jacobian
assembly (ideal-gas case).

Where CP1 (``test_nr_tableau_gas_liquid.py``) only proved tableau
*topology*, this file proves the volume-aware *numerics* CP2 wires in:

* ``TestInertGasMatchesPartitionModel`` — a gas-liquid species with no
  acid-base coupling (O2) must solve to exactly the same split as
  ``HenryEquilibrium.equilibrium_a_moles()``'s own closed-form answer,
  since with no pH coupling the simultaneous fold has nothing to
  simultaneously solve for.
* ``TestCoupledMassConservation`` — CO2 + carbonate ladder + Henry: total
  carbon (gas CO2 + liquid CO2 + HCO3- + CO3--) must be exactly conserved
  by the simultaneous solve.
* ``TestSNIAConsistencyAndDivergence`` — the plan's own test requirement:
  fine-substepped SNIA (many small sequential speciation+transfer steps,
  continuously updating the pH-driven alpha correction) converges to the
  same answer as the simultaneous fold; a single coarse macro-step (what
  ``cv.advance()`` does today) diverges from it, because it uses the pH
  from *before* the transfer to decide how much should transfer, never
  re-checking that the transfer itself would have shifted that pH. The
  SNIA reference reproduces ``KineticGasLiquidLink._kinetic_flux``'s own
  exact exponential-step formula (verified against that source), not an
  approximation of it.
* ``TestDifferentiableLiquidModel`` — the §8.4 Davies/SIT Jacobian
  extension, checked against finite differences.

Design note on the divergence test (see PR discussion): the coarse/fine
gap is not a clean function of kLa alone — it also depends on β (Henry
favorability), which is itself pH-dependent, so β can dominate over kLa.
This file ships one concrete, realistic scenario (AD biogas sparging into
lightly-buffered digestate) rather than a kLa sweep implying strict
monotonicity that the underlying physics doesn't actually have.
"""
from __future__ import annotations

import math

import pytest


# ═══════════════════════════════════════════════════════════════════════════
#  Shared chemistry builders
# ═══════════════════════════════════════════════════════════════════════════

def _water_rxn():
    from PyOMES.chemistry.common_species import H2O, H_plus, OH_minus
    from PyOMES.reactions.equilibrium import EquilibriumReaction
    from PyOMES.reactions.stoichiometry import StoichiometryEntry as E
    return EquilibriumReaction(
        stoichiometry=[
            E(species=H2O, phase="liquid", coefficient=-1.0),
            E(species=H_plus, phase="liquid", coefficient=+1.0),
            E(species=OH_minus, phase="liquid", coefficient=+1.0),
        ],
        log_K=-14.0, balance_elements=("H", "O"), label="water",
    )


def _carbonate_ladder():
    from PyOMES.chemistry.common_species import H2O, H_plus, CO2, HCO3_minus, CO3_2minus
    from PyOMES.reactions.equilibrium import EquilibriumReaction
    from PyOMES.reactions.stoichiometry import StoichiometryEntry as E
    co2_first = EquilibriumReaction(
        stoichiometry=[
            E(species=CO2, phase="liquid", coefficient=-1.0),
            E(species=H2O, phase="liquid", coefficient=-1.0),
            E(species=HCO3_minus, phase="liquid", coefficient=+1.0),
            E(species=H_plus, phase="liquid", coefficient=+1.0),
        ],
        log_K=-6.35, total_id="CO2", balance_elements=("C", "H", "O"), label="co2_first",
    )
    co2_second = EquilibriumReaction(
        stoichiometry=[
            E(species=HCO3_minus, phase="liquid", coefficient=-1.0),
            E(species=CO3_2minus, phase="liquid", coefficient=+1.0),
            E(species=H_plus, phase="liquid", coefficient=+1.0),
        ],
        log_K=-10.33, total_id="CO2", balance_elements=("C", "H", "O"), label="co2_second",
    )
    return [co2_first, co2_second]


def _co2_henry(H_ref=3.4e-4, dlnH=2400.0):
    from PyOMES.chemistry import HenryEquilibrium
    return HenryEquilibrium(
        H_ref=H_ref, dlnH=dlnH, gas_species="CO2", liquid_species="CO2",
        label="henry_CO2",
    )


# ═══════════════════════════════════════════════════════════════════════════
#  Inert gas (no acid-base coupling) matches PartitionModel exactly
# ═══════════════════════════════════════════════════════════════════════════

class TestInertGasMatchesPartitionModel:
    """O2: no acid-base ladder, so the folded solve has nothing to
    simultaneously couple with — it must reduce to exactly the same
    answer as HenryEquilibrium's own closed-form PartitionModel role."""

    def test_o2_matches_equilibrium_a_moles(self):
        from PyOMES.chemistry.species import Species
        from PyOMES.chemistry import HenryEquilibrium
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

        O2 = Species(id="O2", atoms={"O": 2}, charge=0)
        henry = HenryEquilibrium(H_ref=1.3e-5, dlnH=1500.0,
                                  gas_species=O2, liquid_species=O2, label="henry_O2")
        engine = NRChemicalEquilibriumEngine.from_reactions([_water_rxn(), henry], T_K=298.15)

        V_liq, V_gas, T_K = 1.0, 0.5, 298.15
        n_total_O2 = 0.01
        out = engine.solve(
            totals={"O2": n_total_O2 / V_liq}, strong_ions={},
            V_liq_L=V_liq, V_gas_L=V_gas,
        )

        n_liq_expected = henry.equilibrium_a_moles(n_total_O2, V_liq, V_gas, T_K)
        C_liq_expected = n_liq_expected / V_liq
        n_gas_expected = n_total_O2 - n_liq_expected
        p_gas_expected = n_gas_expected * 0.0820574 * T_K / V_gas

        assert out.species_mol_L["O2"] == pytest.approx(C_liq_expected, rel=1e-8)
        assert out.partial_pressures_atm["O2"] == pytest.approx(p_gas_expected, rel=1e-8)


# ═══════════════════════════════════════════════════════════════════════════
#  Coupled CO2 + carbonate: exact mass conservation
# ═══════════════════════════════════════════════════════════════════════════

class TestCoupledMassConservation:

    def test_total_carbon_conserved_across_phases(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

        engine = NRChemicalEquilibriumEngine.from_reactions(
            [_water_rxn()] + _carbonate_ladder() + [_co2_henry()], T_K=298.15,
        )
        V_liq, V_gas, T_K = 1.0, 0.2, 298.15
        n_total_C = 0.05
        out = engine.solve(
            totals={"CO2": n_total_C / V_liq}, strong_ions={},
            V_liq_L=V_liq, V_gas_L=V_gas,
        )

        C_liq_total = (
            out.species_mol_L["CO2"] + out.species_mol_L["HCO3-"]
            + out.species_mol_L["CO3--"]
        )
        n_liq = C_liq_total * V_liq
        n_gas = out.partial_pressures_atm["CO2"] * V_gas / (0.0820574 * T_K)
        assert (n_liq + n_gas) == pytest.approx(n_total_C, abs=1e-9)

    def test_pH_in_plausible_range(self):
        """Sanity check: 0.05 mol/L total carbon with no added base should
        give a mildly acidic pH (CO2 is a weak acid), not something
        pathological."""
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

        engine = NRChemicalEquilibriumEngine.from_reactions(
            [_water_rxn()] + _carbonate_ladder() + [_co2_henry()], T_K=298.15,
        )
        out = engine.solve(
            totals={"CO2": 0.05}, strong_ions={}, V_liq_L=1.0, V_gas_L=0.2,
        )
        assert 3.0 < out.pH < 7.0


# ═══════════════════════════════════════════════════════════════════════════
#  SNIA reference — reproduces KineticGasLiquidLink._kinetic_flux exactly
# ═══════════════════════════════════════════════════════════════════════════

def _ab_only_engine(T_K=308.15):
    from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
    return NRChemicalEquilibriumEngine.from_reactions(
        [_water_rxn()] + _carbonate_ladder(), T_K=T_K,
    )


def _alpha_from_speciation(ab_engine, n_liq_CT_mol, V_liq, strong_ions):
    """Molecular CO2 fraction of total dissolved carbon — matches
    KineticGasLiquidLink._alpha_for's own definition exactly."""
    out = ab_engine.solve(
        totals={"CO2": n_liq_CT_mol / V_liq}, strong_ions=strong_ions,
    )
    CT = out.species_mol_L["CO2"] + out.species_mol_L["HCO3-"] + out.species_mol_L["CO3--"]
    alpha = out.species_mol_L["CO2"] / CT if CT > 0 else 1.0
    return alpha, out.pH


def _snia_run(
    ab_engine, henry, n_liq0, n_gas0, V_liq, V_gas, T_K, strong_ions,
    kLa, dt_total_h, n_substeps,
):
    """Sequential (SNIA) integration over dt_total_h, split into
    n_substeps equal sub-steps. Each sub-step: (1) re-solve acid-base
    speciation on the *current* liquid total to get a fresh alpha, (2)
    apply KineticGasLiquidLink._kinetic_flux's exact analytical
    exponential-step formula (total-IC basis) for that sub-step's dt.

    n_substeps=1 reproduces exactly what a single cv.advance() call does
    today (one alpha evaluation per macro step); large n_substeps is the
    "very fine SNIA sub-stepping" ground truth the plan's CP2 test bullet
    calls for.
    """
    dt = dt_total_h / n_substeps
    n_liq, n_gas = n_liq0, n_gas0
    for _ in range(n_substeps):
        alpha, _ = _alpha_from_speciation(ab_engine, n_liq, V_liq, strong_ions)
        n_total = n_liq + n_gas
        beta = henry.partition_ratio(V_liq, V_gas, T_K, alpha=alpha)
        n_liq_eq = beta * n_total / (1.0 + beta)
        B = kLa * (1.0 + beta)
        delta = n_liq_eq - n_liq
        exponent = B * dt
        frac = 1.0 if exponent > 50.0 else 1.0 - math.exp(-exponent)
        n_liq = n_liq + delta * frac
        n_gas = n_total - n_liq
    return n_liq, n_gas


def _simultaneous_n_liq(sim_engine, n_total_C, V_liq, V_gas, T_K, strong_ions):
    out = sim_engine.solve(
        totals={"CO2": n_total_C / V_liq}, strong_ions=strong_ions,
        V_liq_L=V_liq, V_gas_L=V_gas,
    )
    CT = out.species_mol_L["CO2"] + out.species_mol_L["HCO3-"] + out.species_mol_L["CO3--"]
    return CT * V_liq


# AD biogas-sparging scenario: a pulse of CO2-rich biogas injected into
# lightly-buffered digestate, 90% of total carbon starting in the gas
# phase. Realistic AD-scale volumes/temperature/timestep/kLa.
_V_LIQ = 2.0
_V_GAS = 0.5
_T_K = 308.15
_N_TOTAL_C = 0.15
_STRONG_IONS = {"CT_Na": 0.01}
_N_LIQ0 = _N_TOTAL_C * 0.10
_N_GAS0 = _N_TOTAL_C * 0.90
_DT_H = 0.02       # one realistic cv.advance() macro step
_KLA = 100.0        # realistic sparged-AD kLa (h^-1)


class TestSNIAConsistencyAndDivergence:

    @pytest.fixture(scope="class")
    def engines(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        ab_engine = _ab_only_engine(T_K=_T_K)
        sim_engine = NRChemicalEquilibriumEngine.from_reactions(
            [_water_rxn()] + _carbonate_ladder() + [_co2_henry()], T_K=_T_K,
        )
        return ab_engine, sim_engine

    @pytest.fixture(scope="class")
    def henry(self):
        return _co2_henry()

    def test_fine_substepped_snia_matches_simultaneous_fold(self, engines, henry):
        """Consistency check (plan's own CP2 test bullet): many small
        sequential steps, continuously refreshing the pH/alpha
        correction, converge to the same equilibrium as the one-shot
        simultaneous fold."""
        ab_engine, sim_engine = engines
        n_liq_sim = _simultaneous_n_liq(
            sim_engine, _N_TOTAL_C, _V_LIQ, _V_GAS, _T_K, _STRONG_IONS,
        )
        n_liq_fine, _ = _snia_run(
            ab_engine, henry, _N_LIQ0, _N_GAS0, _V_LIQ, _V_GAS, _T_K,
            _STRONG_IONS, kLa=_KLA, dt_total_h=1.0, n_substeps=2000,
        )
        # 1.0 h total elapsed (many multiples of the equilibration
        # timescale at kLa=100) lets fine sub-stepping actually reach
        # equilibrium, unlike the single 0.02 h macro-step below.
        assert n_liq_fine == pytest.approx(n_liq_sim, rel=1e-3)

    def test_coarse_single_step_diverges_from_simultaneous_fold(self, engines, henry):
        """Divergence demonstration: a single coarse macro-step (what
        cv.advance() does today, one alpha evaluation before the
        transfer) overshoots the true simultaneous-fold answer by >20%
        for this realistic AD-sparging scenario, because it never
        re-checks that dissolving that much CO2 would itself drop the pH
        (and hence the Henry solubility) partway through the step."""
        ab_engine, sim_engine = engines
        n_liq_sim = _simultaneous_n_liq(
            sim_engine, _N_TOTAL_C, _V_LIQ, _V_GAS, _T_K, _STRONG_IONS,
        )
        n_liq_coarse, _ = _snia_run(
            ab_engine, henry, _N_LIQ0, _N_GAS0, _V_LIQ, _V_GAS, _T_K,
            _STRONG_IONS, kLa=_KLA, dt_total_h=_DT_H, n_substeps=1,
        )
        rel_err = (n_liq_coarse - n_liq_sim) / n_liq_sim
        assert rel_err > 0.20, (
            f"expected coarse SNIA to overshoot the simultaneous fold by "
            f">20% in this scenario; got {rel_err:.1%}"
        )

    def test_fine_substepping_within_same_macro_step_tracks_simultaneous_direction(
        self, engines, henry,
    ):
        """Within the *same* single 0.02 h macro-step budget (not run to
        full 1.0 h convergence), refining resolution moves the answer
        toward the simultaneous fold's direction rather than away from
        it — confirming the coarse step's error is a discretization
        artifact of the stale one-shot alpha, not a difference in the
        underlying physical model."""
        ab_engine, sim_engine = engines
        n_liq_sim = _simultaneous_n_liq(
            sim_engine, _N_TOTAL_C, _V_LIQ, _V_GAS, _T_K, _STRONG_IONS,
        )
        n_liq_coarse, _ = _snia_run(
            ab_engine, henry, _N_LIQ0, _N_GAS0, _V_LIQ, _V_GAS, _T_K,
            _STRONG_IONS, kLa=_KLA, dt_total_h=_DT_H, n_substeps=1,
        )
        n_liq_finer, _ = _snia_run(
            ab_engine, henry, _N_LIQ0, _N_GAS0, _V_LIQ, _V_GAS, _T_K,
            _STRONG_IONS, kLa=_KLA, dt_total_h=_DT_H, n_substeps=200,
        )
        assert abs(n_liq_finer - n_liq_sim) < abs(n_liq_coarse - n_liq_sim)


# ═══════════════════════════════════════════════════════════════════════════
#  DifferentiableLiquidModel — §8.4 Jacobian extension
# ═══════════════════════════════════════════════════════════════════════════

class TestDifferentiableLiquidModel:

    def _finite_difference(self, model, x_mol, charge, T_K=298.15, eps=1e-7):
        import numpy as np
        species_ids = sorted(x_mol)
        n = len(species_ids)
        g0 = model.gamma_all(x_mol, T_K, charge=charge)
        Jfd = np.zeros((n, n))
        for j, sp_j in enumerate(species_ids):
            x2 = dict(x_mol)
            x2[sp_j] = x2[sp_j] + eps
            g1 = model.gamma_all(x2, T_K, charge=charge)
            for i, sp_i in enumerate(species_ids):
                Jfd[i, j] = (g1.get(sp_i, 1.0) - g0.get(sp_i, 1.0)) / eps
        return Jfd

    def test_davies_satisfies_differentiable_liquid_model(self):
        from PyOMES.thermo.liquid_phase_model import (
            DaviesLiquidModel, DifferentiableLiquidModel,
        )
        assert isinstance(DaviesLiquidModel(), DifferentiableLiquidModel)

    def test_sit_satisfies_differentiable_liquid_model(self):
        from PyOMES.thermo.sit_liquid_model import SITLiquidModel
        from PyOMES.thermo.liquid_phase_model import DifferentiableLiquidModel
        assert isinstance(SITLiquidModel(), DifferentiableLiquidModel)

    def test_davies_jacobian_matches_finite_difference(self):
        from PyOMES.thermo.liquid_phase_model import DaviesLiquidModel
        model = DaviesLiquidModel()
        x_mol = {"Na+": 0.05, "Cl-": 0.05, "Ca++": 0.01, "CO2": 0.02}
        charge = {"Na+": 1, "Cl-": -1, "Ca++": 2, "CO2": 0}
        J = model.jacobian_dgamma_dx(x_mol, 298.15, charge=charge)
        Jfd = self._finite_difference(model, x_mol, charge)
        assert J == pytest.approx(Jfd, abs=1e-4)

    def test_sit_jacobian_dh_term_only_disagrees_with_full_fd(self):
        """SIT's jacobian intentionally differentiates only the DH term
        (not the ion-pair epsilon cross-terms — see the method's own
        docstring) — this test documents that limitation rather than
        asserting a false equivalence with the full finite difference."""
        from PyOMES.thermo.sit_liquid_model import SITLiquidModel
        model = SITLiquidModel()
        x_mol = {"Na+": 0.05, "Cl-": 0.05, "Ca++": 0.01, "CO2": 0.02}
        charge = {"Na+": 1, "Cl-": -1, "Ca++": 2, "CO2": 0}
        J = model.jacobian_dgamma_dx(x_mol, 298.15, charge=charge)
        Jfd = self._finite_difference(model, x_mol, charge)
        assert not (J == pytest.approx(Jfd, abs=1e-2))

    def test_neutral_species_have_zero_jacobian_row_and_column(self):
        from PyOMES.thermo.liquid_phase_model import DaviesLiquidModel
        model = DaviesLiquidModel()
        x_mol = {"Na+": 0.05, "Cl-": 0.05, "CO2": 0.02}
        charge = {"Na+": 1, "Cl-": -1, "CO2": 0}
        J = model.jacobian_dgamma_dx(x_mol, 298.15, charge=charge)
        species_ids = sorted(x_mol)
        co2_idx = species_ids.index("CO2")
        assert (J[co2_idx, :] == 0).all()
        assert (J[:, co2_idx] == 0).all()

    def test_dilute_limit_returns_zero_matrix(self):
        from PyOMES.thermo.liquid_phase_model import DaviesLiquidModel
        model = DaviesLiquidModel()
        x_mol = {"Na+": 0.0, "Cl-": 0.0}
        charge = {"Na+": 1, "Cl-": -1}
        J = model.jacobian_dgamma_dx(x_mol, 298.15, charge=charge)
        assert (J == 0).all()
