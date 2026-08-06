# -*- coding: utf-8 -*-
"""Tests for CP4 of LAYER1_GAP_CLOSURE: Raoult / evaporation folding.

Two independent pieces, per the design discussion that shaped this
checkpoint:

* ``TestRaoultTableauFold`` — the standalone NR-tableau capability
  (matching CP1/CP2's pattern, applied to H2O via RaoultEquilibrium).
  Unlike CO2/NH3/H2S, water has no tracked "total" component — its
  liquid-side activity is fixed at 1 by the universal pure-solvent
  convention (the same one that already excludes it from the acid-base
  graph), so the fold produces a T-dependent *constant* relation
  (``nu={}``) with no coupling to any master and no mass-balance
  contribution. Water mass conservation is *not* handled by this fold —
  see the next class.
* ``TestWaterVapourBoundaryRetirement`` — the actual mass-conservation
  fix, and the reason WaterVapourBoundary/VentWaterLoss are retired:
  ``transfer_models={"H2O": EquilibriumTransferModel(RaoultEquilibrium())}``
  properly moves mass between the *tracked* liquid and gas H2O pools
  (§6.2 of MASS_EXCHANGE_ARCHITECTURE.md), unlike WaterVapourBoundary
  (which "conjured" gas-phase water from nowhere). This is engine-
  agnostic — it works whether the CV's reaction_system uses
  NRChemicalEquilibriumEngine or the old BisectionChemicalEquilibriumEngine — which is why it, not
  the NR-tableau fold, is what closes the mass balance for
  models/vlmodels/adm1/base.py's build_adm1_cv() (which still uses the
  old engine).
"""
from __future__ import annotations

import math

import pytest


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


# ═══════════════════════════════════════════════════════════════════════════
#  NR-tableau Raoult fold — standalone capability
# ═══════════════════════════════════════════════════════════════════════════

class TestRaoultTableauFold:

    @pytest.fixture(scope="class")
    def engine(self):
        from PyOMES.chemistry import RaoultEquilibrium
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        return NRChemicalEquilibriumEngine.from_reactions(
            [_water_rxn(), RaoultEquilibrium()], T_K=298.15,
        )

    def test_h2o_does_not_become_a_master(self, engine):
        """Water has no tracked total — folding it must not create a new
        component/master, unlike CO2/NH3/H2S in CP1."""
        assert engine.tableau.masters == ["H+"]
        assert engine.tableau.components == []

    def test_h2o_gas_secondary_present_with_empty_nu(self, engine):
        gas_secs = [s for s in engine.tableau.secondaries if s.phase == "gas"]
        assert len(gas_secs) == 1
        sec = gas_secs[0]
        assert sec.species_id == "H2O"
        assert sec.nu == {}
        assert sec.element_stoichiometry == {}

    def test_partial_pressure_matches_raoult_P_sat(self, engine):
        """The fold's value must exactly match RaoultEquilibrium's own
        P_sat(T) — the whole point of the constant-relation design."""
        from PyOMES.chemistry import RaoultEquilibrium
        raoult = RaoultEquilibrium()
        out = engine.solve(totals={}, strong_ions={}, V_liq_L=1.0, V_gas_L=0.2)
        assert out.partial_pressures_atm["H2O"] == pytest.approx(
            raoult.P_sat(298.15), rel=1e-9,
        )

    def test_partial_pressure_temperature_dependence_matches_raoult(self):
        """Van't Hoff correction via the tableau must agree with
        RaoultEquilibrium.P_sat() at a non-reference temperature too."""
        from PyOMES.chemistry import RaoultEquilibrium
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

        raoult = RaoultEquilibrium()
        engine = NRChemicalEquilibriumEngine.from_reactions(
            [_water_rxn(), raoult], T_K=308.15,
        )
        out = engine.solve(totals={}, strong_ions={}, V_liq_L=1.0, V_gas_L=0.2)
        assert out.partial_pressures_atm["H2O"] == pytest.approx(
            raoult.P_sat(308.15), rel=1e-9,
        )

    def test_h2o_fold_does_not_perturb_ph(self, engine):
        """A constant relation with nu={} must not affect the H+ solve at
        all — pH should be neutral (no other totals declared)."""
        out = engine.solve(totals={}, strong_ions={}, V_liq_L=1.0, V_gas_L=0.2)
        assert out.pH == pytest.approx(7.0, abs=1e-6)

    def test_h2o_fold_coexists_with_ordinary_acid_base_chemistry(self):
        """A CO2 ladder alongside the Raoult H2O fold: CO2 gets a real
        component/master (finite total), H2O does not — both coexist in
        the same tableau without interfering."""
        from PyOMES.chemistry.common_species import (
            H2O, H_plus, CO2, HCO3_minus, CO3_2minus,
        )
        from PyOMES.chemistry import RaoultEquilibrium
        from PyOMES.reactions.equilibrium import EquilibriumReaction
        from PyOMES.reactions.stoichiometry import StoichiometryEntry as E
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

        co2_first = EquilibriumReaction(
            stoichiometry=[
                E(species=CO2, phase="liquid", coefficient=-1.0),
                E(species=H2O, phase="liquid", coefficient=-1.0),
                E(species=HCO3_minus, phase="liquid", coefficient=+1.0),
                E(species=H_plus, phase="liquid", coefficient=+1.0),
            ],
            log_K=-6.35, total_id="CO2", balance_elements=("C", "H", "O"),
            label="co2_first",
        )
        co2_second = EquilibriumReaction(
            stoichiometry=[
                E(species=HCO3_minus, phase="liquid", coefficient=-1.0),
                E(species=CO3_2minus, phase="liquid", coefficient=+1.0),
                E(species=H_plus, phase="liquid", coefficient=+1.0),
            ],
            log_K=-10.33, total_id="CO2", balance_elements=("C", "H", "O"),
            label="co2_second",
        )
        engine = NRChemicalEquilibriumEngine.from_reactions(
            [_water_rxn(), co2_first, co2_second, RaoultEquilibrium()],
            T_K=298.15,
        )
        assert set(engine.tableau.masters) == {"H+", "CO2"}
        out = engine.solve(
            totals={"CO2": 0.05}, strong_ions={}, V_liq_L=1.0, V_gas_L=0.2,
        )
        assert out.species_mol_L["CO2"] > 0.0
        assert out.partial_pressures_atm["H2O"] == pytest.approx(0.03169, rel=1e-6)


# ═══════════════════════════════════════════════════════════════════════════
#  WaterVapourBoundary/VentWaterLoss retirement
# ═══════════════════════════════════════════════════════════════════════════

class TestRetiredClassesGone:

    def test_water_vapour_boundary_not_importable(self):
        import PyOMES.core as core
        assert not hasattr(core, "WaterVapourBoundary")

    def test_vent_water_loss_not_importable(self):
        import PyOMES.core as core
        assert not hasattr(core, "VentWaterLoss")

    def test_water_vapour_boundary_not_in_boundaries_module(self):
        import PyOMES.core.boundaries as boundaries
        assert not hasattr(boundaries, "WaterVapourBoundary")
        assert not hasattr(boundaries, "VentWaterLoss")


class TestWaterVapourBoundaryRetirement:
    """The replacement mechanism: transfer_models=EquilibriumTransferModel(
    RaoultEquilibrium()) closes the liquid <-> gas H2O mass balance
    automatically, without either retired boundary object. Engine-
    agnostic — no NRChemicalEquilibriumEngine/reaction_system involved at all,
    matching how build_adm1_cv() uses it (that model still uses the old
    BisectionChemicalEquilibriumEngine)."""

    def _make_cv(self, *, liq_h2o_mol, gas_h2o_mol=0.0, V_liq=1.6, V_gas=0.4, T_K=308.15):
        from PyOMES.chemistry import RaoultEquilibrium
        from PyOMES.core.phases import GasPhase, LiquidPhase
        from PyOMES.core.control_volume import ControlVolume
        from PyOMES.core.transfer_models import EquilibriumTransferModel

        gas = GasPhase(n_mol={"H2O": gas_h2o_mol}, V_L=V_gas, T_K=T_K)
        liq = LiquidPhase(n_mol={"H2O": liq_h2o_mol}, V_L=V_liq, T_K=T_K)
        return ControlVolume(
            phases={"gas": gas, "liquid": liq},
            transfer_models={"H2O": EquilibriumTransferModel(RaoultEquilibrium())},
            label="water_only",
        )

    def test_liquid_decrements_as_gas_gains_water(self):
        """Unlike WaterVapourBoundary (which never touched the liquid
        phase), the replacement must draw the gas-phase gain from the
        liquid phase's own tracked pool."""
        cv = self._make_cv(liq_h2o_mol=88.8, gas_h2o_mol=0.0)
        n_liq_before = cv.phases["liquid"].n_mol["H2O"]
        n_gas_before = cv.phases["gas"].n_mol["H2O"]

        cv.step_internal_transfer(dt_h=0.1)

        n_liq_after = cv.phases["liquid"].n_mol["H2O"]
        n_gas_after = cv.phases["gas"].n_mol["H2O"]
        assert n_gas_after > n_gas_before
        assert n_liq_after < n_liq_before

    def test_total_water_conserved_across_transfer(self):
        cv = self._make_cv(liq_h2o_mol=88.8, gas_h2o_mol=0.0)
        total_before = cv.phases["liquid"].n_mol["H2O"] + cv.phases["gas"].n_mol["H2O"]

        cv.step_internal_transfer(dt_h=0.1)

        total_after = cv.phases["liquid"].n_mol["H2O"] + cv.phases["gas"].n_mol["H2O"]
        assert total_after == pytest.approx(total_before, rel=1e-9)

    def test_converges_toward_saturation_partial_pressure(self):
        """Repeated steps should drive gas-phase H2O partial pressure
        toward RaoultEquilibrium's P_sat(T), matching
        WaterVapourBoundary's own stated behaviour (it targeted the same
        physical saturation state, just via a different, non-conserving
        mechanism)."""
        from PyOMES.chemistry import RaoultEquilibrium
        cv = self._make_cv(liq_h2o_mol=88.8, gas_h2o_mol=0.0, T_K=308.15)
        for _ in range(50):
            cv.step_internal_transfer(dt_h=0.05)
        p_final = cv.phases["gas"].p_atm.get("H2O", 0.0)
        expected = RaoultEquilibrium().P_sat(308.15)
        assert p_final == pytest.approx(expected, rel=1e-3)

    def test_empty_liquid_pool_transfers_nothing(self):
        """Sanity check on the mass-conservation fix: with no liquid
        water tracked at all (the old, unseeded default), the transfer
        correctly moves ~nothing — this is exactly the silent-regression
        risk seeding guards against, made explicit as a regression test."""
        cv = self._make_cv(liq_h2o_mol=0.0, gas_h2o_mol=0.0)
        cv.step_internal_transfer(dt_h=0.1)
        assert cv.phases["gas"].n_mol["H2O"] == pytest.approx(0.0, abs=1e-9)
