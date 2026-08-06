# -*- coding: utf-8 -*-
"""Tests for MembraneGasBoundary (Stage C)."""

import pytest
from PyOMES.core.phases import GasPhase, LiquidPhase, R_L_ATM_MOL_K
from PyOMES.core.control_volume import ControlVolume
from PyOMES.core.boundaries import MembraneGasBoundary, ExternalBoundary, apply_boundary
from PyOMES.chemistry import HenryPartition


def _hp(kH: float, dlnH: float = 0.0) -> HenryPartition:
    return HenryPartition(H_ref=kH * 1000.0 / 101325.0, dlnH=dlnH)


def _make_gas_cv(n_mol, V_L=1e-4, T_K=310.15):
    return ControlVolume(phases={"gas": GasPhase(n_mol=dict(n_mol), V_L=V_L, T_K=T_K)})

def _make_membrane(**kw):
    d = dict(permeability={"O2": 0.03, "CO2": 0.15, "N2": 0.015}, area_m2=3.14e-5)
    d.update(kw)
    return MembraneGasBoundary(**d)

def _make_well_plate_cv(boundaries=None, reaction_system=None, kLa=None, depleted_o2=False):
    """Well plate fermenter CV. If depleted_o2=True, start with no O2 in headspace."""
    from PyOMES.core.control_volume import ControlVolume
    from PyOMES.core.gas_liquid_link import KineticGasLiquidLink
    V_gas, V_liq, T_K = 1e-4, 2e-4, 310.15
    n_N2 = 0.7808 * V_gas / (R_L_ATM_MOL_K * T_K)
    if depleted_o2:
        n_O2 = 0.0  # depleted headspace
    else:
        n_O2 = 0.2095 * V_gas / (R_L_ATM_MOL_K * T_K)
    lkw = dict(gas_cv_key="gas", gas_phase_key="gas", liquid_cv_key="liquid", liquid_phase_key="liquid",
               partition_models={"O2": _hp(1.3e-3), "CO2": _hp(3.4e-2), "N2": _hp(6.5e-4)})
    if kLa: lkw["kLa"] = kLa; lkw["equilibrium_species"] = {"N2"}
    else: lkw["equilibrium_species"] = {"O2", "CO2", "N2"}
    return ControlVolume(
        phases={
            "gas": GasPhase(n_mol={"O2": n_O2, "N2": n_N2}, V_L=V_gas, T_K=T_K),
            "liquid": LiquidPhase(n_mol={"O2": 0.0, "CO2": 0.0}, V_L=V_liq, T_K=T_K),
        },
        internal_interfaces=[KineticGasLiquidLink(**lkw)],
        boundaries=boundaries, reaction_system=reaction_system)


class TestO2Inward:
    def test_o2_permeates_inward(self):
        flux = _make_membrane().compute_flux(_make_gas_cv({"N2": 1e-5}), 0.01)
        assert flux.get("O2", 0.0) > 0.0
    def test_o2_flux_proportional_to_dp(self):
        cv1 = _make_gas_cv({"O2": 0.0, "N2": 1e-5})
        cv2 = _make_gas_cv({"O2": 0.1*1e-4/(R_L_ATM_MOL_K*310.15), "N2": 1e-5})
        m = _make_membrane()
        assert m.compute_flux(cv1, 0.01)["O2"] > m.compute_flux(cv2, 0.01)["O2"] > 0

class TestCO2Outward:
    def test_co2_permeates_outward(self):
        cv = _make_gas_cv({"CO2": 0.05*1e-4/(R_L_ATM_MOL_K*310.15), "N2": 1e-5})
        assert _make_membrane().compute_flux(cv, 0.01).get("CO2", 0.0) < 0.0

class TestZeroFlux:
    def test_no_flux_at_equilibrium(self):
        V, T = 1e-4, 310.15
        atm = {"O2": 0.2095, "CO2": 0.0004, "N2": 0.7808}
        n = {s: p*V/(R_L_ATM_MOL_K*T) for s, p in atm.items()}
        flux = _make_membrane().compute_flux(_make_gas_cv(n, V, T), 0.01)
        for s in atm: assert abs(flux.get(s, 0.0)) < 1e-15

class TestScaling:
    def test_scales_with_area(self):
        cv = _make_gas_cv({"N2": 1e-5})
        f1 = MembraneGasBoundary(permeability={"O2": 0.03}, area_m2=1e-5).compute_flux(cv, 0.01)["O2"]
        f2 = MembraneGasBoundary(permeability={"O2": 0.03}, area_m2=2e-5).compute_flux(cv, 0.01)["O2"]
        assert f2 == pytest.approx(2.0 * f1, rel=1e-12)
    def test_scales_with_permeability(self):
        cv = _make_gas_cv({"N2": 1e-5})
        f1 = MembraneGasBoundary(permeability={"O2": 0.01}, area_m2=1e-5).compute_flux(cv, 0.01)["O2"]
        f2 = MembraneGasBoundary(permeability={"O2": 0.02}, area_m2=1e-5).compute_flux(cv, 0.01)["O2"]
        assert f2 == pytest.approx(2.0 * f1, rel=1e-12)

class TestEdgeCases:
    def test_unlisted_species(self):
        cv = _make_gas_cv({"O2": 1e-6, "Ar": 1e-6})
        assert "Ar" not in MembraneGasBoundary(permeability={"O2": 0.03}, area_m2=1e-5).compute_flux(cv, 0.01)
    def test_missing_gas_phase(self):
        cv = ControlVolume(phases={"liquid": LiquidPhase(n_mol={}, V_L=1.0, T_K=300.0)})
        assert _make_membrane().compute_flux(cv, 0.01) == {}
    def test_zero_area(self):
        assert MembraneGasBoundary(permeability={"O2": 0.03}, area_m2=0.0).compute_flux(_make_gas_cv({"N2": 1e-5}), 0.01) == {}

class TestAtmosphere:
    def test_default_is_standard_air(self):
        assert _make_membrane().external_atmosphere["O2"] == pytest.approx(0.2095)
    def test_custom_atmosphere(self):
        m = MembraneGasBoundary(permeability={"CO2": 0.15}, area_m2=1e-5, external_atmosphere={"CO2": 0.05})
        assert m.external_atmosphere["CO2"] == pytest.approx(0.05)
    def test_co2_incubator_reduces_outflow(self):
        cv = _make_gas_cv({"CO2": 0.05*1e-4/(R_L_ATM_MOL_K*310.15), "N2": 1e-6})
        f_std = MembraneGasBoundary(permeability={"CO2": 0.15}, area_m2=1e-5).compute_flux(cv, 0.01).get("CO2", 0)
        f_inc = MembraneGasBoundary(permeability={"CO2": 0.15}, area_m2=1e-5, external_atmosphere={"CO2": 0.05}).compute_flux(cv, 0.01).get("CO2", 0)
        assert f_std < 0.0
        assert abs(f_inc) < abs(f_std)

class TestMutators:
    def test_set_external_atmosphere(self):
        m = _make_membrane(); m.set_external_atmosphere({"O2": 0.95})
        assert m.external_atmosphere["O2"] == pytest.approx(0.95)
    def test_set_permeability(self):
        m = _make_membrane(); m.set_permeability("O2", 0.1)
        assert m.permeability["O2"] == pytest.approx(0.1)

class TestProtocol:
    def test_satisfies_external_boundary(self):
        assert isinstance(_make_membrane(), ExternalBoundary)
    def test_phase_key(self):
        assert _make_membrane().phase_key == "gas"
    def test_label(self):
        assert _make_membrane(label="x").label == "x"
    def test_repr(self):
        assert "MembraneGasBoundary" in repr(_make_membrane())

class TestApplyBoundary:
    def test_adds_o2(self):
        cv = _make_gas_cv({"N2": 1e-5})
        rec = apply_boundary(_make_membrane(), cv, 0.01)
        assert rec.mol_applied.get("O2", 0.0) > 0.0

class TestWithFermenterCV:
    def test_membrane_on_depleted_cv_adds_o2(self):
        """Fermenter CV with depleted O₂ headspace: membrane should add O₂."""
        cv = _make_well_plate_cv(boundaries=[_make_membrane()], depleted_o2=True)
        tb = cv.total_mol()
        cv.advance(dt_h=0.01)
        ta = cv.total_mol()
        assert ta.get("O2", 0.0) > tb.get("O2", 0.0)

class TestWellPlateIntegration:
    def test_o2_steady_state(self):
        """Well plate with O₂ consumption reaches steady state via membrane supply."""
        from PyOMES.reactions import KineticReaction, StoichiometryEntry
        from PyOMES.chemistry import Species
        O2 = Species(id="O2", atoms={"O": 2})
        Waste = Species(id="Waste", atoms={"O": 2})
        rxn = KineticReaction([StoichiometryEntry(species=O2, phase="liquid", coefficient=-1.0),
                               StoichiometryEntry(species=Waste, phase="liquid", coefficient=+1.0)],
                              rate_fn=lambda env: 0.5*env.S("O2")*env.V_L, balance_elements=("O",))
        cv = _make_well_plate_cv(boundaries=[_make_membrane()], reaction_system=rxn)
        hist = []
        for _ in range(2000):
            cv.advance(dt_h=0.001)
            hist.append(cv.phases["gas"].n_mol.get("O2", 0.0))
        assert hist[-1] > 0.0, "O2 should not be fully depleted"
        tail = hist[-100:]
        assert (max(tail)-min(tail))/max(sum(tail)/len(tail), 1e-30) < 0.05, "Should reach steady state"

    def test_co2_accumulates(self):
        """CO₂ produced by metabolism accumulates in headspace and vents through membrane."""
        from PyOMES.reactions import KineticReaction, StoichiometryEntry
        from PyOMES.core.control_volume import ControlVolume
        from PyOMES.core.gas_liquid_link import KineticGasLiquidLink
        from PyOMES.chemistry import Species
        A = Species(id="A", atoms={"C": 1})
        CO2 = Species(id="CO2", atoms={"C": 1, "O": 2})
        rxn = KineticReaction([StoichiometryEntry(species=A, phase="liquid", coefficient=-1.0),
                               StoichiometryEntry(species=CO2, phase="liquid", coefficient=+1.0)],
                              rate_fn=lambda env: 1e-6, balance_elements=("C",))
        V_g, V_l, T = 1e-4, 2e-4, 310.15
        n_N2 = 0.78*V_g/(R_L_ATM_MOL_K*T)
        link = KineticGasLiquidLink(gas_cv_key="gas", gas_phase_key="gas",
                                     liquid_cv_key="liquid", liquid_phase_key="liquid",
                                     partition_models={"CO2": _hp(3.4e-2), "N2": _hp(6.5e-4)}, kLa={"CO2": 500.0},
                                     equilibrium_species={"N2"})
        membrane = MembraneGasBoundary(permeability={"CO2": 0.15}, area_m2=3.14e-5)
        cv = ControlVolume(
            phases={
                "gas": GasPhase(n_mol={"N2": n_N2, "CO2": 0.0}, V_L=V_g, T_K=T),
                "liquid": LiquidPhase(n_mol={"A": 1.0, "CO2": 0.0}, V_L=V_l, T_K=T),
            },
            internal_interfaces=[link],
            reaction_system=rxn,
            boundaries=[membrane])
        for _ in range(500): cv.advance(dt_h=0.001)
        assert cv.phases["gas"].n_mol.get("CO2", 0.0) > 0.0
