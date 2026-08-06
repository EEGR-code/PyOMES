"""Tests for fermenter.core — phases, interfaces, and ControlVolume."""

import pytest
import math
from PyOMES.core.phases import GasPhase, LiquidPhase, SolidPhase, R_L_ATM_MOL_K, Phase
from PyOMES.core.interfaces import PhaseInterface, TransferDiagnostics
from PyOMES.core.control_volume import ControlVolume


# ════════════════════════════════════════════════════════════════════════
#  GasPhase tests
# ════════════════════════════════════════════════════════════════════════

class TestGasPhase:
    def test_construction(self):
        g = GasPhase({"O2": 1.0, "N2": 3.0}, V_L=100.0, T_K=300.0)
        assert g.n_mol["O2"] == 1.0
        assert g.n_mol["N2"] == 3.0
        assert g.V_L == 100.0
        assert g.T_K == 300.0

    def test_total_mol(self):
        g = GasPhase({"O2": 1.0, "CO2": 0.5, "N2": 3.0}, V_L=100.0, T_K=300.0)
        assert g.n_total == pytest.approx(4.5)

    def test_pressure_ideal_gas(self):
        g = GasPhase({"O2": 1.0}, V_L=100.0, T_K=300.0)
        expected = 1.0 * R_L_ATM_MOL_K * 300.0 / 100.0
        assert g.P_atm == pytest.approx(expected, rel=1e-6)

    def test_mole_fractions(self):
        g = GasPhase({"O2": 1.0, "N2": 3.0}, V_L=100.0, T_K=300.0)
        y = g.y
        assert y["O2"] == pytest.approx(0.25)
        assert y["N2"] == pytest.approx(0.75)

    def test_partial_pressures_sum_to_total(self):
        g = GasPhase({"O2": 1.0, "CO2": 0.5, "N2": 3.0}, V_L=100.0, T_K=300.0)
        p = g.p_atm
        assert sum(p.values()) == pytest.approx(g.P_atm, rel=1e-10)

    def test_apply_flux_adds_moles(self):
        g = GasPhase({"O2": 1.0}, V_L=100.0, T_K=300.0)
        g.apply_flux({"O2": 10.0}, dt_h=0.1)  # +1 mol
        assert g.n_mol["O2"] == pytest.approx(2.0)

    def test_apply_flux_removes_moles(self):
        g = GasPhase({"O2": 2.0}, V_L=100.0, T_K=300.0)
        g.apply_flux({"O2": -10.0}, dt_h=0.1)  # -1 mol
        assert g.n_mol["O2"] == pytest.approx(1.0)

    def test_apply_flux_floors_at_zero(self):
        g = GasPhase({"O2": 0.5}, V_L=100.0, T_K=300.0)
        g.apply_flux({"O2": -100.0}, dt_h=1.0)  # would go to -99.5
        assert g.n_mol["O2"] == 0.0

    def test_apply_flux_new_species(self):
        g = GasPhase({"O2": 1.0}, V_L=100.0, T_K=300.0)
        g.apply_flux({"CO2": 5.0}, dt_h=0.2)  # +1 mol CO2
        assert g.n_mol["CO2"] == pytest.approx(1.0)

    def test_snapshot_independence(self):
        g = GasPhase({"O2": 1.0}, V_L=100.0, T_K=300.0)
        snap = g.snapshot()
        g.apply_flux({"O2": 100.0}, dt_h=1.0)
        assert snap.n_mol["O2"] == 1.0  # unchanged
        assert g.n_mol["O2"] == 101.0

    def test_zero_volume_raises(self):
        with pytest.raises(ValueError, match="must be positive"):
            GasPhase({"O2": 1.0}, V_L=0.0, T_K=300.0)

    def test_negative_volume_raises(self):
        with pytest.raises(ValueError, match="must be positive"):
            GasPhase({"O2": 1.0}, V_L=-1.0, T_K=300.0)

    def test_zero_temperature_raises(self):
        with pytest.raises(ValueError, match="must be positive"):
            GasPhase({"O2": 1.0}, V_L=1.0, T_K=0.0)

    def test_negative_temperature_raises(self):
        with pytest.raises(ValueError, match="must be positive"):
            GasPhase({"O2": 1.0}, V_L=1.0, T_K=-10.0)

    def test_volume_setter_rejects_zero(self):
        g = GasPhase({"O2": 1.0}, V_L=1.0, T_K=300.0)
        with pytest.raises(ValueError, match="must be positive"):
            g.V_L = 0.0

    def test_temperature_setter_rejects_negative(self):
        g = GasPhase({"O2": 1.0}, V_L=1.0, T_K=300.0)
        with pytest.raises(ValueError, match="must be positive"):
            g.T_K = -5.0

    def test_satisfies_protocol(self):
        g = GasPhase({"O2": 1.0}, V_L=100.0, T_K=300.0)
        assert isinstance(g, Phase)

    def test_properties_empty_by_default(self):
        g = GasPhase({"O2": 1.0}, V_L=1.0, T_K=300.0)
        assert g.properties == {}

    def test_properties_accepts_initial_values(self):
        g = GasPhase({"O2": 1.0}, V_L=1.0, T_K=300.0, properties={"density_kg_L": 1.2})
        assert g.properties["density_kg_L"] == pytest.approx(1.2)

    def test_snapshot_copies_properties(self):
        g = GasPhase({"O2": 1.0}, V_L=1.0, T_K=300.0, properties={"density_kg_L": 1.2})
        snap = g.snapshot()
        assert snap.properties["density_kg_L"] == pytest.approx(1.2)
        snap.properties["density_kg_L"] = 99.0
        assert g.properties["density_kg_L"] == pytest.approx(1.2)


# ════════════════════════════════════════════════════════════════════════
#  LiquidPhase tests
# ════════════════════════════════════════════════════════════════════════

class TestLiquidPhase:
    def test_construction(self):
        liq = LiquidPhase({"AceticAcid": 0.01}, V_L=1.0, T_K=305.15)
        assert liq.n_mol["AceticAcid"] == 0.01

    def test_zero_volume_raises(self):
        with pytest.raises(ValueError, match="must be positive"):
            LiquidPhase({"A": 1.0}, V_L=0.0, T_K=300.0)

    def test_negative_volume_raises(self):
        with pytest.raises(ValueError, match="must be positive"):
            LiquidPhase({"A": 1.0}, V_L=-1.0, T_K=300.0)

    def test_zero_temperature_raises(self):
        with pytest.raises(ValueError, match="must be positive"):
            LiquidPhase({"A": 1.0}, V_L=1.0, T_K=0.0)

    def test_volume_setter_rejects_zero(self):
        liq = LiquidPhase({"A": 1.0}, V_L=1.0, T_K=300.0)
        with pytest.raises(ValueError, match="must be positive"):
            liq.V_L = 0.0

    def test_temperature_setter_rejects_negative(self):
        liq = LiquidPhase({"A": 1.0}, V_L=1.0, T_K=300.0)
        with pytest.raises(ValueError, match="must be positive"):
            liq.T_K = -5.0

    def test_concentrations_mol_L(self):
        liq = LiquidPhase({"AceticAcid": 0.5}, V_L=2.0, T_K=305.15)
        assert liq.concentrations_mol_L["AceticAcid"] == pytest.approx(0.25)

    def test_concentrations_g_L(self):
        liq = LiquidPhase(
            {"AceticAcid": 0.5}, V_L=2.0, T_K=305.15,
            mw={"AceticAcid": 60.052},
        )
        assert liq.concentrations_g_L["AceticAcid"] == pytest.approx(0.5 * 60.052 / 2.0)

    def test_concentrations_g_L_missing_mw(self):
        liq = LiquidPhase({"X": 0.1}, V_L=1.0, T_K=305.15, mw={})
        assert "X" not in liq.concentrations_g_L

    def test_pH_from_n_mol_H_plus(self):
        """pH is derived from n_mol["H+"] / V_L (state-unification C3).
        For pH 6.5, [H+] = 10^-6.5 mol/L. Storage is in moles, so
        n_mol["H+"] = [H+] * V_L."""
        H_mol_L = 10 ** -6.5
        liq = LiquidPhase(
            {"CO2": 0.01, "H+": H_mol_L * 1.0},
            V_L=1.0, T_K=305.15,
        )
        assert liq.pH == pytest.approx(6.5)

    def test_pH_raises_when_no_H_plus_in_n_mol(self):
        """Phase without speciation has no n_mol['H+']; pH access raises
        (state-unification C3)."""
        liq = LiquidPhase({"CO2": 0.01}, V_L=1.0, T_K=305.15)
        with pytest.raises(ValueError, match=r"'H\+'"):
            _ = liq.pH

    def test_apply_flux(self):
        liq = LiquidPhase({"CO2": 1.0}, V_L=1.0, T_K=305.15)
        liq.apply_flux({"CO2": -5.0}, dt_h=0.1)  # -0.5 mol
        assert liq.n_mol["CO2"] == pytest.approx(0.5)

    def test_snapshot_independence(self):
        liq = LiquidPhase({"CO2": 1.0}, V_L=1.0, T_K=305.15,
                           speciation={"pH": 7.0})
        snap = liq.snapshot()
        liq.speciation["pH"] = 3.0
        assert snap.speciation["pH"] == 7.0

    def test_satisfies_protocol(self):
        liq = LiquidPhase({"CO2": 0.1}, V_L=1.0, T_K=305.15)
        assert isinstance(liq, Phase)


# ════════════════════════════════════════════════════════════════════════
#  SolidPhase tests
# ════════════════════════════════════════════════════════════════════════

class TestSolidPhase:
    def test_construction(self):
        s = SolidPhase({"CaCO3": 0.5}, V_L=0.01, T_K=305.15,
                        mw={"CaCO3": 100.09}, surface_area_m2=2.0)
        assert s.n_mol["CaCO3"] == 0.5
        assert s.surface_area_m2 == 2.0

    def test_mass_g(self):
        s = SolidPhase({"CaCO3": 0.5}, V_L=0.01, T_K=305.15,
                        mw={"CaCO3": 100.09})
        assert s.mass_g["CaCO3"] == pytest.approx(50.045)

    def test_apply_flux(self):
        s = SolidPhase({"CaCO3": 1.0}, V_L=0.01, T_K=305.15)
        s.apply_flux({"CaCO3": -2.0}, dt_h=0.25)  # -0.5 mol
        assert s.n_mol["CaCO3"] == pytest.approx(0.5)

    def test_satisfies_protocol(self):
        s = SolidPhase({"CaCO3": 0.1}, V_L=0.01, T_K=305.15)
        assert isinstance(s, Phase)

    def test_properties_empty_by_default(self):
        s = SolidPhase({"CaCO3": 0.1}, V_L=0.01, T_K=305.15)
        assert s.properties == {}

    def test_snapshot_copies_properties(self):
        s = SolidPhase({"CaCO3": 0.1}, V_L=0.01, T_K=305.15, properties={"porosity_eff": 0.4})
        snap = s.snapshot()
        assert snap.properties["porosity_eff"] == pytest.approx(0.4)
        snap.properties["porosity_eff"] = 99.0
        assert s.properties["porosity_eff"] == pytest.approx(0.4)


# ════════════════════════════════════════════════════════════════════════
#  Mock interface for testing ControlVolume
# ════════════════════════════════════════════════════════════════════════

class MockGLInterface:
    """Simple mock: transfers a fixed rate of CO2 from gas to liquid."""

    def __init__(self, rate: float = 1.0, rate_mol_per_h: float = None):
        self.rate = float(rate_mol_per_h if rate_mol_per_h is not None else rate)

    @property
    def phase_a_key(self) -> str:
        return "gas"

    @property
    def phase_b_key(self) -> str:
        return "liquid"

    def compute_flux(self, state_a, state_b, dt_h, property_results=None):
        # Fixed-rate CO2 transfer from gas to liquid
        return {"CO2": self.rate}


class MockSLInterface:
    """Simple mock: dissolves solid into liquid at a fixed rate."""

    def __init__(self, rate: float = 0.1, rate_mol_per_h: float = None):
        self.rate = float(rate_mol_per_h if rate_mol_per_h is not None else rate)

    @property
    def phase_a_key(self) -> str:
        return "solid"

    @property
    def phase_b_key(self) -> str:
        return "liquid"

    def compute_flux(self, state_a, state_b, dt_h, property_results=None):
        return {"CaCO3": self.rate}


# ════════════════════════════════════════════════════════════════════════
#  ControlVolume tests
# ════════════════════════════════════════════════════════════════════════

class TestControlVolume:
    def test_construction_gas_liquid(self):
        gas = GasPhase({"O2": 1.0, "CO2": 0.5, "N2": 3.0}, V_L=200.0, T_K=305.15)
        liq = LiquidPhase({"CO2": 0.01}, V_L=800.0, T_K=305.15)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq}, label="fermenter")
        assert "gas" in cv
        assert "liquid" in cv
        assert cv.label == "fermenter"

    def test_getitem(self):
        gas = GasPhase({"O2": 1.0}, V_L=100.0, T_K=300.0)
        cv = ControlVolume(phases={"gas": gas})
        assert cv["gas"] is gas

    def test_total_mol_aggregates(self):
        gas = GasPhase({"CO2": 0.5}, V_L=100.0, T_K=300.0)
        liq = LiquidPhase({"CO2": 0.3}, V_L=800.0, T_K=300.0)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})
        totals = cv.total_mol()
        assert totals["CO2"] == pytest.approx(0.8)

    def test_invalid_interface_key_raises(self):
        gas = GasPhase({"O2": 1.0}, V_L=100.0, T_K=300.0)
        with pytest.raises(KeyError, match="phase_b_key"):
            ControlVolume(
                phases={"gas": gas},
                internal_interfaces=[MockGLInterface()],
            )

    def test_phase_keys(self):
        gas = GasPhase({"O2": 1.0}, V_L=100.0, T_K=300.0)
        liq = LiquidPhase({"CO2": 0.01}, V_L=800.0, T_K=300.0)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})
        assert set(cv.phase_keys) == {"gas", "liquid"}


def _make_gl_cv(rate=1.0):
    """Helper: build a gas-liquid CV with a mock CO2 transfer interface."""
    gas = GasPhase({"O2": 1.0, "CO2": 2.0, "N2": 3.0}, V_L=200.0, T_K=305.15)
    liq = LiquidPhase({"CO2": 0.1}, V_L=800.0, T_K=305.15)
    iface = MockGLInterface(rate_mol_per_h=rate)
    return ControlVolume(
        phases={"gas": gas, "liquid": liq},
        internal_interfaces=[iface],
        label="test_gl",
    )


class TestControlVolumeInternalTransfer:
    """Tests for step_internal_transfer with mass balance diagnostics."""

    def test_transfer_conserves_mass(self):
        cv = _make_gl_cv(rate=10.0)
        diag = cv.step_internal_transfer(dt_h=0.01)
        assert diag.is_conserved()

    def test_transfer_moves_material(self):
        cv = _make_gl_cv(rate=10.0)
        co2_gas_before = cv["gas"].n_mol["CO2"]
        co2_liq_before = cv["liquid"].n_mol["CO2"]
        diag = cv.step_internal_transfer(dt_h=0.1)
        # 10 mol/h * 0.1 h = 1 mol transferred
        assert cv["gas"].n_mol["CO2"] == pytest.approx(co2_gas_before - 1.0)
        assert cv["liquid"].n_mol["CO2"] == pytest.approx(co2_liq_before + 1.0)

    def test_residual_is_zero(self):
        cv = _make_gl_cv(rate=5.0)
        diag = cv.step_internal_transfer(dt_h=0.02)
        assert diag.max_abs_residual == pytest.approx(0.0, abs=1e-15)

    def test_untransferred_species_unchanged(self):
        cv = _make_gl_cv(rate=5.0)
        o2_before = cv["gas"].n_mol["O2"]
        n2_before = cv["gas"].n_mol["N2"]
        cv.step_internal_transfer(dt_h=0.01)
        assert cv["gas"].n_mol["O2"] == o2_before
        assert cv["gas"].n_mol["N2"] == n2_before

    def test_diagnostics_fluxes_recorded(self):
        cv = _make_gl_cv(rate=3.0)
        diag = cv.step_internal_transfer(dt_h=0.01)
        assert len(diag.fluxes) == 1
        assert diag.fluxes[0]["CO2"] == pytest.approx(3.0)

    def test_no_interfaces_is_noop(self):
        gas = GasPhase({"O2": 1.0}, V_L=100.0, T_K=300.0)
        cv = ControlVolume(phases={"gas": gas})
        diag = cv.step_internal_transfer(dt_h=0.1)
        assert diag.is_conserved()
        assert diag.total_residual == 0.0
        assert len(diag.fluxes) == 0

    def test_floor_at_zero_breaks_conservation(self):
        """When a flux would drive a phase negative, the floor at 0
        means total moles decrease.  This is expected and detectable."""
        gas = GasPhase({"CO2": 0.01}, V_L=100.0, T_K=300.0)
        liq = LiquidPhase({"CO2": 0.0}, V_L=800.0, T_K=300.0)
        iface = MockGLInterface(rate_mol_per_h=100.0)  # huge rate
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            internal_interfaces=[iface],
        )
        diag = cv.step_internal_transfer(dt_h=1.0)
        # Gas wanted to lose 100 mol but only had 0.01 → floored to 0
        # Liquid gained 100 mol
        # Total CO2 before: 0.01, after: 0 + 100 = 100 → NOT conserved
        assert not diag.is_conserved()

    def test_summary_string(self):
        cv = _make_gl_cv(rate=1.0)
        diag = cv.step_internal_transfer(dt_h=0.01)
        s = diag.summary()
        assert "Conserved" in s
        assert "Total mol" in s


class TestControlVolumeThreePhase:
    """Tests for a three-phase (gas-liquid-solid) CV."""

    def test_three_phase_construction(self):
        gas = GasPhase({"O2": 1.0, "CO2": 0.5}, V_L=200.0, T_K=305.15)
        liq = LiquidPhase({"CO2": 0.01, "CaCO3": 0.0}, V_L=800.0, T_K=305.15)
        sol = SolidPhase({"CaCO3": 1.0}, V_L=0.5, T_K=305.15)
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq, "solid": sol},
            internal_interfaces=[MockGLInterface(rate=0.5), MockSLInterface(rate=0.1)],
            label="three_phase",
        )
        assert len(cv.phase_keys) == 3
        assert len(cv.internal_interfaces) == 2

    def test_three_phase_transfer_conserves(self):
        gas = GasPhase({"O2": 1.0, "CO2": 2.0}, V_L=200.0, T_K=305.15)
        liq = LiquidPhase({"CO2": 0.1, "CaCO3": 0.0}, V_L=800.0, T_K=305.15)
        sol = SolidPhase({"CaCO3": 5.0}, V_L=0.5, T_K=305.15)
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq, "solid": sol},
            internal_interfaces=[MockGLInterface(rate=1.0), MockSLInterface(rate=0.5)],
        )
        diag = cv.step_internal_transfer(dt_h=0.01)
        assert diag.is_conserved()

    def test_total_mol_aggregates_three_phases(self):
        gas = GasPhase({"CO2": 1.0}, V_L=100.0, T_K=300.0)
        liq = LiquidPhase({"CO2": 0.5, "CaCO3": 0.1}, V_L=800.0, T_K=300.0)
        sol = SolidPhase({"CaCO3": 2.0}, V_L=0.5, T_K=300.0)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq, "solid": sol})
        totals = cv.total_mol()
        assert totals["CO2"] == pytest.approx(1.5)
        assert totals["CaCO3"] == pytest.approx(2.1)


class TestControlVolumeExternalFlux:
    """Tests for apply_external_flux."""

    def test_external_flux_changes_inventory(self):
        gas = GasPhase({"O2": 1.0}, V_L=100.0, T_K=300.0)
        cv = ControlVolume(phases={"gas": gas})
        mol_before = cv.total_mol()
        cv.apply_external_flux("gas", {"O2": 10.0}, dt_h=0.5)  # +5 mol
        assert cv["gas"].n_mol["O2"] == pytest.approx(6.0)

    def test_external_flux_not_tracked_by_internal_diag(self):
        """External fluxes are intentional inventory changes, not
        conservation violations."""
        gas = GasPhase({"O2": 1.0}, V_L=100.0, T_K=300.0)
        cv = ControlVolume(phases={"gas": gas})
        cv.apply_external_flux("gas", {"O2": 100.0}, dt_h=1.0)
        # Internal transfer with no interfaces should still report conserved
        diag = cv.step_internal_transfer(dt_h=0.01)
        assert diag.is_conserved()


class TestControlVolumeSnapshot:
    def test_snapshot_is_independent(self):
        gas = GasPhase({"O2": 1.0, "CO2": 0.5}, V_L=200.0, T_K=305.15)
        liq = LiquidPhase({"CO2": 0.1}, V_L=800.0, T_K=305.15)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq}, label="orig")
        snap = cv.snapshot()

        cv.apply_external_flux("gas", {"O2": 100.0}, dt_h=1.0)
        assert snap["gas"].n_mol["O2"] == 1.0  # unchanged


class TestTransferDiagnostics:
    def test_empty_diagnostics(self):
        d = TransferDiagnostics()
        assert d.is_conserved()
        assert d.total_residual == 0.0
        assert d.max_abs_residual == 0.0

    def test_nonzero_residual(self):
        d = TransferDiagnostics(
            species_before={"CO2": 1.0},
            species_after={"CO2": 1.001},
            residual={"CO2": 0.001},
            total_mol_before=1.0,
            total_mol_after=1.001,
        )
        assert not d.is_conserved(tol_mol=1e-6)
        assert d.is_conserved(tol_mol=0.01)
        assert d.max_abs_residual == pytest.approx(0.001)
