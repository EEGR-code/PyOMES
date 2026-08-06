# -*- coding: utf-8 -*-
"""Tests for KineticGasLiquidLink (Stage A).

Validates:
- O₂ absorption (kinetic): positive flux, scales with kLa
- CO₂ stripping (kinetic): negative flux when liquid supersaturated
- Mixed modes (kinetic O₂, equilibrium N₂)
- Equilibrium mode reaches exact Henry partition
- CO₂ kinetic mode uses molecular CO₂(aq) from speciation
- CO₂ kinetic mode falls back to total dissolved without speciation
- CO₂ equilibrium mode uses α₀-corrected effective Henry constant
- Conservation: gas + liquid total unchanged
- Safety clamping: doesn't exceed available inventory
- Zero kLa gives zero flux
- Missing CVs returns empty
- set_kLa and set_kLa_with_co2_ratio mutators
- Protocol compliance (CVLink)
- Integration: kinetic absorption approaches equilibrium over many steps
"""

import numpy as np
import pytest

from PyOMES.core.phases import GasPhase, LiquidPhase, R_L_ATM_MOL_K
from PyOMES.core.control_volume import ControlVolume
from PyOMES.core.gas_liquid_link import KineticGasLiquidLink
from PyOMES.chemistry import HenryPartition


def _hp(kH_mol_L_atm: float, dlnH: float = 0.0) -> HenryPartition:
    """Build a HenryPartition from a mol/L/atm kH for test fixtures."""
    return HenryPartition(H_ref=kH_mol_L_atm * 1000.0 / 101325.0, dlnH=dlnH)


# ═══════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════

def _make_two_cv_system(n_gas, n_liq, V_gas=0.4, V_liq=1.6, T_K=305.15):
    """Create a gas CV + liquid CV pair as a dict."""
    gas = GasPhase(n_mol=dict(n_gas), V_L=V_gas, T_K=T_K)
    liq = LiquidPhase(n_mol=dict(n_liq), V_L=V_liq, T_K=T_K)
    cv_gas = ControlVolume(phases={"gas": gas}, label="gas")
    cv_liq = ControlVolume(phases={"liquid": liq}, label="liquid")
    return {"gas": cv_gas, "liquid": cv_liq}


def _henry_equilibrium_conc(kH, n_gas_species, V_gas, T_K):
    """Compute Henry equilibrium dissolved concentration from gas moles."""
    p = n_gas_species * R_L_ATM_MOL_K * T_K / V_gas
    return kH * p


# ═══════════════════════════════════════════════════════════════════════
#  Kinetic mode: O₂ absorption
# ═══════════════════════════════════════════════════════════════════════

class TestKineticO2Absorption:

    def test_o2_absorbs_from_gas_to_liquid(self):
        """O₂ should transfer from gas to liquid (positive flux)."""
        cvs = _make_two_cv_system(
            n_gas={"O2": 0.5}, n_liq={"O2": 0.0},
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3)}, kLa={"O2": 150.0},
        )
        flow = link.compute_flow(cvs, dt_h=0.01)
        assert "O2" in flow
        assert flow["O2"] > 0.0  # gas → liquid

    def test_flux_scales_with_kla(self):
        """Doubling kLa should roughly double the flux (in the linear regime)."""
        # Use small kLa × dt so the exponential is in its linear regime
        cvs1 = _make_two_cv_system(n_gas={"O2": 0.5}, n_liq={"O2": 0.0})
        cvs2 = _make_two_cv_system(n_gas={"O2": 0.5}, n_liq={"O2": 0.0})

        link1 = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3)}, kLa={"O2": 1.0},  # kLa*dt=0.001 (deeply linear)
        )
        link2 = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3)}, kLa={"O2": 2.0},  # kLa*dt=0.002 (deeply linear)
        )
        f1 = link1.compute_flow(cvs1, dt_h=0.001)["O2"]
        f2 = link2.compute_flow(cvs2, dt_h=0.001)["O2"]
        assert f2 == pytest.approx(2.0 * f1, rel=0.01)

    def test_no_flux_at_equilibrium(self):
        """At Henry equilibrium, kinetic flux should be zero."""
        V_gas, V_liq, T_K = 0.4, 1.6, 305.15
        kH = 1.3e-3
        n_O2_gas = 0.5
        C_star = _henry_equilibrium_conc(kH, n_O2_gas, V_gas, T_K)
        n_O2_liq = C_star * V_liq

        cvs = _make_two_cv_system(
            n_gas={"O2": n_O2_gas}, n_liq={"O2": n_O2_liq},
            V_gas=V_gas, V_liq=V_liq, T_K=T_K,
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(kH)}, kLa={"O2": 150.0},
        )
        flow = link.compute_flow(cvs, dt_h=0.01)
        assert flow.get("O2", 0.0) == pytest.approx(0.0, abs=1e-12)

    def test_zero_kla_gives_zero_flux(self):
        cvs = _make_two_cv_system(n_gas={"O2": 0.5}, n_liq={"O2": 0.0})
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3)}, kLa={"O2": 0.0},
        )
        flow = link.compute_flow(cvs, dt_h=0.01)
        assert flow == {}


# ═══════════════════════════════════════════════════════════════════════
#  Kinetic mode: CO₂ stripping
# ═══════════════════════════════════════════════════════════════════════

class TestKineticCO2Stripping:

    def test_co2_strips_from_liquid_to_gas(self):
        """CO₂ supersaturated in liquid should strip (negative flux)."""
        cvs = _make_two_cv_system(
            n_gas={"CO2": 0.001},  # very little gas CO₂
            n_liq={"CO2": 0.1},   # lots of dissolved CO₂
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"CO2": _hp(3.4e-2)}, kLa={"CO2": 135.0},
        )
        flow = link.compute_flow(cvs, dt_h=0.01)
        assert "CO2" in flow
        assert flow["CO2"] < 0.0  # liquid → gas

    def test_co2_kinetic_uses_alpha_from_speciation(self):
        """When speciation is available, driving force should use the
        molecular fraction read from ``PropertyResult.alphas``.

        Exercises the alphas channel introduced in
        ``chemistry-unification-2`` — the link no longer reconstructs
        the molecular fraction inline from pH + pKas (Method 1) or by
        dividing ``species[mol_key] / (n_liq / V_liq)`` (Method 2);
        both paths collapsed to a single ``alphas[mol_key]`` read.
        """
        # state-unification C4: alpha is computed inline from
        # liq_phase.n_mol. To get alpha=0.16, populate n_mol with
        # CO2aq=0.016 mol, HCO3-=0.084 mol (sum=0.1 = total TIC).
        cvs = _make_two_cv_system(
            n_gas={"CO2": 0.001},
            n_liq={"CO2": 0.016, "HCO3-": 0.084},
            V_liq=1.6,
        )

        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"CO2": _hp(3.4e-2)}, kLa={"CO2": 135.0},
        )
        gas_phase = cvs["gas"].phases["gas"]
        liq_phase = cvs["liquid"].phases["liquid"]
        flow_with_spec = link.compute_flux(
            gas_phase, liq_phase, dt_h=0.01,
        )

        # Without speciation: α defaults to 1.0, full TIC drives stripping.
        # n_liq stored under the legacy "CO2" key so the link's
        # pre-speciation-state fallback finds it.
        cvs2 = _make_two_cv_system(
            n_gas={"CO2": 0.001}, n_liq={"CO2": 0.1}, V_liq=1.6,
        )
        flow_without_spec = link.compute_flux(
            cvs2["gas"].phases["gas"],
            cvs2["liquid"].phases["liquid"],
            dt_h=0.01,
        )

        # Stripping should be weaker with speciation (only molecular
        # CO2aq participates in transfer, not total TIC).
        assert abs(flow_with_spec["CO2"]) < abs(flow_without_spec["CO2"])


# ═══════════════════════════════════════════════════════════════════════
#  Equilibrium mode
# ═══════════════════════════════════════════════════════════════════════

class TestEquilibriumMode:

    def test_equilibrium_reaches_henry_partition(self):
        """Equilibrium mode should partition to exactly C = kH × p."""
        V_gas, V_liq, T_K = 0.4, 1.6, 305.15
        kH_N2 = 6.5e-4
        n_total = 0.5  # all in gas initially

        cvs = _make_two_cv_system(
            n_gas={"N2": n_total}, n_liq={"N2": 0.0},
            V_gas=V_gas, V_liq=V_liq, T_K=T_K,
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"N2": _hp(kH_N2)},
            equilibrium_species={"N2"},
        )
        flow = link.compute_flow(cvs, dt_h=1.0)

        # Apply the flux manually and check Henry condition
        n_transferred = flow["N2"] * 1.0
        n_gas_eq = n_total - n_transferred
        n_liq_eq = n_transferred

        p_eq = n_gas_eq * R_L_ATM_MOL_K * T_K / V_gas
        C_eq = n_liq_eq / V_liq

        assert C_eq == pytest.approx(kH_N2 * p_eq, rel=1e-10)

    def test_equilibrium_already_at_equilibrium_gives_zero(self):
        """If already at Henry equilibrium, flux should be zero."""
        V_gas, V_liq, T_K = 0.4, 1.6, 305.15
        kH = 6.5e-4
        # Set up at equilibrium
        n_total = 0.5
        beta = kH * R_L_ATM_MOL_K * T_K * V_liq / V_gas
        n_liq = beta * n_total / (1.0 + beta)
        n_gas = n_total - n_liq

        cvs = _make_two_cv_system(
            n_gas={"N2": n_gas}, n_liq={"N2": n_liq},
            V_gas=V_gas, V_liq=V_liq, T_K=T_K,
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"N2": _hp(kH)}, equilibrium_species={"N2"},
        )
        flow = link.compute_flow(cvs, dt_h=1.0)
        assert flow.get("N2", 0.0) == pytest.approx(0.0, abs=1e-12)


# ═══════════════════════════════════════════════════════════════════════
#  Mixed modes
# ═══════════════════════════════════════════════════════════════════════

class TestMixedModes:

    def test_kinetic_o2_equilibrium_n2(self):
        """O₂ kinetic and N₂ equilibrium in the same link."""
        cvs = _make_two_cv_system(
            n_gas={"O2": 0.3, "N2": 2.0},
            n_liq={"O2": 0.0, "N2": 0.0},
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3), "N2": _hp(6.5e-4)},
            kLa={"O2": 150.0},
            equilibrium_species={"N2"},
        )
        flow = link.compute_flow(cvs, dt_h=0.01)

        # Both should transfer (gas→liquid)
        assert flow.get("O2", 0.0) > 0.0
        assert flow.get("N2", 0.0) > 0.0

        # N₂ equilibrium flux should be much larger than O₂ kinetic flux
        # (equilibrium reaches partition in one step, kinetic is gradual)
        assert abs(flow["N2"]) > abs(flow["O2"])


# ═══════════════════════════════════════════════════════════════════════
#  CO₂ equilibrium with speciation correction
# ═══════════════════════════════════════════════════════════════════════

class TestCO2EquilibriumWithSpeciation:

    def test_co2_equilibrium_uses_alpha0(self):
        """CO₂ equilibrium should use α₀-corrected Henry constant.

        Reads α₀ from ``PropertyResult.alphas["CO2aq"]`` populated by
        the speciation engine (chemistry-unification-2).
        """
        V_gas, V_liq, T_K = 0.4, 1.6, 305.15
        kH_CO2 = 3.4e-2

        # state-unification C4: alpha computed inline from n_mol.
        # For α₀=0.5, populate CO2aq=0.008 mol, HCO3-=0.008 mol
        # (sum=0.016, ratio=0.5).
        cvs = _make_two_cv_system(
            n_gas={"CO2": 0.05},
            n_liq={"CO2": 0.008, "HCO3-": 0.008},
            V_gas=V_gas, V_liq=V_liq, T_K=T_K,
        )

        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"CO2": _hp(kH_CO2)},
            equilibrium_species={"CO2"},
        )
        gas_phase = cvs["gas"].phases["gas"]
        liq_phase = cvs["liquid"].phases["liquid"]
        flow_with_spec = link.compute_flux(
            gas_phase, liq_phase, dt_h=1.0,
        )

        # Without speciation (α₀=1): kH_eff = kH. n_liq stored under
        # the legacy "CO2" key so the link's pre-speciation-state
        # fallback finds it.
        cvs2 = _make_two_cv_system(
            n_gas={"CO2": 0.05}, n_liq={"CO2": 0.016},
            V_gas=V_gas, V_liq=V_liq, T_K=T_K,
        )
        flow_without = link.compute_flux(
            cvs2["gas"].phases["gas"],
            cvs2["liquid"].phases["liquid"],
            dt_h=1.0,
        )

        # With α₀ < 1, effective kH is larger → more dissolves → larger flux
        assert abs(flow_with_spec.get("CO2", 0.0)) > abs(flow_without.get("CO2", 0.0))


# ═══════════════════════════════════════════════════════════════════════
#  Alpha correction via _alpha_for
# ═══════════════════════════════════════════════════════════════════════

class TestAlphaCorrection:
    """Verify alpha lookup via speciation_ladders and _alpha_for."""

    def test_alpha_for_returns_none_without_ladder(self):
        """With no speciation_ladders entry, _alpha_for returns None
        (link treats alpha=1 — no correction).
        """
        from PyOMES.core import LiquidPhase
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"CO2": _hp(3.4e-2)},
            equilibrium_species={"CO2"},
        )
        liq = LiquidPhase(n_mol={"O2": 0.001}, V_L=1.6, T_K=305.15)
        assert link._alpha_for("CO2", liq) is None

    def test_alpha_for_reads_ladder_ratio(self):
        """With ladder set, _alpha_for returns n_mol[mol_key] / ladder_sum."""
        from PyOMES.core import LiquidPhase
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"CO2": _hp(3.4e-2)},
            equilibrium_species={"CO2"},
        )
        link.speciation_ladders["CO2"] = ["CO2", "HCO3-", "CO3--"]

        for alpha in (0.5, 0.1, 0.01, 0.99):
            liq = LiquidPhase(
                n_mol={"CO2": alpha, "HCO3-": 1.0 - alpha, "CO3--": 0.0},
                V_L=1.6, T_K=305.15,
            )
            assert link._alpha_for("CO2", liq) == pytest.approx(alpha, abs=1e-15)


# ═══════════════════════════════════════════════════════════════════════
#  Conservation
# ═══════════════════════════════════════════════════════════════════════

class TestConservation:

    def test_kinetic_conserves_total_moles(self):
        """Gas + liquid total must be unchanged by a kinetic flux application."""
        from PyOMES.core import Simulation

        cvs_dict = _make_two_cv_system(
            n_gas={"O2": 0.5, "CO2": 0.1, "N2": 2.0},
            n_liq={"O2": 0.001, "CO2": 0.01, "N2": 0.001},
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3), "CO2": _hp(3.4e-2), "N2": _hp(6.5e-4)},
            kLa={"O2": 150.0, "CO2": 135.0, "N2": 150.0},
        )
        sys = Simulation(cvs=cvs_dict, links=[link])
        mol_before = sys.total_mol()
        sys.run(tau_h=0.01, n_steps=1)
        mol_after = sys.total_mol()

        for sp in ("O2", "CO2", "N2"):
            assert mol_after[sp] == pytest.approx(mol_before[sp], abs=1e-12), \
                f"{sp} not conserved"

    def test_equilibrium_conserves_total_moles(self):
        from PyOMES.core import Simulation

        cvs_dict = _make_two_cv_system(
            n_gas={"N2": 2.0}, n_liq={"N2": 0.0},
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"N2": _hp(6.5e-4)}, equilibrium_species={"N2"},
        )
        sys = Simulation(cvs=cvs_dict, links=[link])
        mol_before = sys.total_mol()
        sys.run(tau_h=1.0, n_steps=1)
        mol_after = sys.total_mol()
        assert mol_after["N2"] == pytest.approx(mol_before["N2"], abs=1e-12)


# ═══════════════════════════════════════════════════════════════════════
#  Safety clamping
# ═══════════════════════════════════════════════════════════════════════

class TestSafetyClamping:

    def test_absorption_clamped_to_gas_inventory(self):
        """When equilibrium amount exceeds gas inventory, flux is clamped."""
        # Very high kH: equilibrium wants more dissolved than total gas available
        cvs = _make_two_cv_system(
            n_gas={"O2": 0.001}, n_liq={"O2": 0.0},
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.0)},   # very high solubility
            kLa={"O2": 1e6},     # very fast transfer
        )
        flow = link.compute_flow(cvs, dt_h=0.1)
        # Max flux = 0.001 / 0.1 = 0.01 mol/h (gas inventory limit)
        assert flow["O2"] <= 0.01 + 1e-12

    def test_stripping_clamped_to_liquid_inventory(self):
        """Stripping flux should not exceed available liquid moles / dt."""
        cvs = _make_two_cv_system(
            n_gas={"CO2": 1e-6},   # near-zero gas
            n_liq={"CO2": 0.001},  # some dissolved
            V_liq=1.0,
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"CO2": _hp(3.4e-2)}, kLa={"CO2": 1e6},
        )
        flow = link.compute_flow(cvs, dt_h=0.1)
        # CO₂ should strip (negative). Max = -0.001 / 0.1 = -0.01 mol/h
        if "CO2" in flow:
            assert flow["CO2"] >= -0.01 - 1e-12

    def test_exponential_form_prevents_overshoot(self):
        """High kLa should not cause oscillation or overshoot of equilibrium."""
        from PyOMES.core import Simulation

        V_gas, V_liq, T_K = 0.4, 1.6, 305.15
        kH = 1.3e-3
        n_total = 0.5

        cvs_dict = _make_two_cv_system(
            n_gas={"O2": n_total}, n_liq={"O2": 0.0},
            V_gas=V_gas, V_liq=V_liq, T_K=T_K,
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(kH)}, kLa={"O2": 500.0},  # kLa*dt=5.0
        )
        sys = Simulation(cvs=cvs_dict, links=[link])

        # Run 20 steps; the per-step liquid n_mol time series comes
        # straight off the BatchResult — convergence should be
        # monotonic, not oscillate.
        result = sys.run(tau_h=0.2, n_steps=20)
        n_liq_series = result.liquid_mol["liquid"]["O2"]
        prev_n_liq = 0.0
        for step, n_liq in enumerate(n_liq_series):
            assert float(n_liq) >= prev_n_liq - 1e-15, (
                f"Step {step}: n_liq decreased"
            )
            prev_n_liq = float(n_liq)


# ═══════════════════════════════════════════════════════════════════════
#  Edge cases
# ═══════════════════════════════════════════════════════════════════════

class TestEdgeCases:

    def test_missing_gas_cv_returns_empty(self):
        link = KineticGasLiquidLink(
            gas_cv_key="missing", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3)}, kLa={"O2": 150.0},
        )
        assert link.compute_flow({}, dt_h=0.01) == {}

    def test_missing_liquid_cv_returns_empty(self):
        cvs = _make_two_cv_system(n_gas={"O2": 0.5}, n_liq={})
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="missing", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3)}, kLa={"O2": 150.0},
        )
        assert link.compute_flow(cvs, dt_h=0.01) == {}

    def test_species_not_in_henry_not_transferred(self):
        """Species without a Henry constant should not be transferred."""
        cvs = _make_two_cv_system(
            n_gas={"O2": 0.5, "Ar": 0.1}, n_liq={},
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3)}, kLa={"O2": 150.0},  # no Ar
        )
        flow = link.compute_flow(cvs, dt_h=0.01)
        assert "Ar" not in flow


# ═══════════════════════════════════════════════════════════════════════
#  Mutators
# ═══════════════════════════════════════════════════════════════════════

class TestMutators:

    def test_set_kla(self):
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3)}, kLa={"O2": 100.0},
        )
        link.set_kLa("O2", 300.0)
        assert link.kLa["O2"] == 300.0

    def test_set_kla_with_co2_ratio(self):
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3), "CO2": _hp(3.4e-2)},
            kLa={"O2": 100.0, "CO2": 90.0},
        )
        link.set_kLa_with_co2_ratio(200.0, co2_ratio=0.9)
        assert link.kLa["O2"] == 200.0
        assert link.kLa["CO2"] == pytest.approx(180.0, abs=1e-10)


# ═══════════════════════════════════════════════════════════════════════
#  Protocol compliance
# ═══════════════════════════════════════════════════════════════════════

class TestProtocol:

    def test_satisfies_cvlink(self):
        from PyOMES.core.links import CVLink
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3)},
        )
        assert isinstance(link, CVLink)

    def test_source_sink_properties(self):
        link = KineticGasLiquidLink(
            gas_cv_key="hs", gas_phase_key="gas",
            liquid_cv_key="bulk", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3)},
        )
        assert link.source_cv_key == "hs"
        assert link.source_phase_key == "gas"
        assert link.sink_cv_key == "bulk"
        assert link.sink_phase_key == "liquid"

    def test_label_auto(self):
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3)},
        )
        assert "gas" in link.label and "liquid" in link.label

    def test_label_custom(self):
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3)}, _label="my_transfer",
        )
        assert link.label == "my_transfer"

    def test_repr(self):
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3), "N2": _hp(6.5e-4)},
            kLa={"O2": 150.0},
            equilibrium_species={"N2"},
        )
        r = repr(link)
        assert "O2" in r and "N2" in r and "eq" in r


# ═══════════════════════════════════════════════════════════════════════
#  Integration: kinetic absorption approaches equilibrium
# ═══════════════════════════════════════════════════════════════════════

class TestIntegrationKineticApproach:

    def test_kinetic_o2_approaches_equilibrium(self):
        """Over many steps, kinetic O₂ absorption should approach Henry equilibrium."""
        from PyOMES.core import Simulation

        V_gas, V_liq, T_K = 0.4, 1.6, 305.15
        kH = 1.3e-3
        n_O2_total = 0.5

        cvs_dict = _make_two_cv_system(
            n_gas={"O2": n_O2_total}, n_liq={"O2": 0.0},
            V_gas=V_gas, V_liq=V_liq, T_K=T_K,
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"O2": _hp(kH)}, kLa={"O2": 500.0},  # high kLa → fast approach
        )
        sys = Simulation(cvs=cvs_dict, links=[link])

        sys.run(tau_h=5.0, n_steps=500)

        # Should be near equilibrium
        n_gas_final = sys["gas"]["gas"].n_mol["O2"]
        n_liq_final = sys["liquid"]["liquid"].n_mol["O2"]
        p_final = n_gas_final * R_L_ATM_MOL_K * T_K / V_gas
        C_final = n_liq_final / V_liq

        assert C_final == pytest.approx(kH * p_final, rel=0.01)

        # Total must be conserved
        assert (n_gas_final + n_liq_final) == pytest.approx(n_O2_total, abs=1e-10)

    def test_kinetic_co2_stripping_approaches_equilibrium(self):
        """CO₂ supersaturated in liquid should strip toward equilibrium."""
        from PyOMES.core import Simulation

        V_gas, V_liq, T_K = 0.4, 1.6, 305.15
        kH = 3.4e-2
        # Start with all CO₂ in liquid (supersaturated)
        n_CO2_total = 0.05

        cvs_dict = _make_two_cv_system(
            n_gas={"CO2": 0.0}, n_liq={"CO2": n_CO2_total},
            V_gas=V_gas, V_liq=V_liq, T_K=T_K,
        )
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"CO2": _hp(kH)}, kLa={"CO2": 500.0},
        )
        sys = Simulation(cvs=cvs_dict, links=[link])

        sys.run(tau_h=5.0, n_steps=500)

        n_gas_final = sys["gas"]["gas"].n_mol.get("CO2", 0.0)
        n_liq_final = sys["liquid"]["liquid"].n_mol.get("CO2", 0.0)

        # Some CO₂ should have moved to gas phase
        assert n_gas_final > 0.001

        # Total conserved
        assert (n_gas_final + n_liq_final) == pytest.approx(n_CO2_total, abs=1e-10)


# ═══════════════════════════════════════════════════════════════════════
#  derive_speciation_keys (chemistry-unification-3)
# ═══════════════════════════════════════════════════════════════════════

class TestDeriveSpeciationKeys:
    """Validate that cross-phase equilibrium reactions populate
    :attr:`KineticGasLiquidLink.speciation_keys` correctly.
    """

    def _make_link(self):
        return KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"CO2": _hp(0.034)},
            speciation_keys={},  # start empty for clarity
        )

    def test_populates_from_cross_phase_reaction(self):
        """A cross-phase equilibrium reaction declaring
        ``CO2(gas) ⇌ CO2(liquid)`` records the gas → liquid mapping
        in ``speciation_keys``.
        """
        from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry, ReactionSystem
        from PyOMES.chemistry.common_species import CO2

        rxn = EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=CO2, phase="gas", coefficient=-1.0),
                StoichiometryEntry(species=CO2, phase="liquid", coefficient=+1.0),
            ],
            balance_elements=("C", "O"),
            label="CO2 partition",
        )
        rxn_set = ReactionSystem([rxn])
        link = self._make_link()
        link.derive_speciation_keys(rxn_set)
        assert link.speciation_keys == {"CO2": "CO2"}

    def test_skips_single_phase_equilibrium(self):
        """Single-phase equilibrium reactions must not populate
        ``speciation_keys`` via the cross-phase path.  The identity
        entry for CO2 still appears because Option A seeds it from
        ``partition_models``.
        """
        from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry, ReactionSystem
        from PyOMES.chemistry.common_species import H_plus, OH_minus, H2O

        rxn = EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=H2O, phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=H_plus, phase="liquid", coefficient=+1.0),
                StoichiometryEntry(species=OH_minus, phase="liquid", coefficient=+1.0),
            ],
            log_K=-14.0,
            balance_elements=("H", "O"),
        )
        link = self._make_link()
        link.derive_speciation_keys(ReactionSystem([rxn]))
        # CO2 seeded from partition_models (Option A identity seed); the water
        # reaction contributed nothing via the cross-phase path.
        assert link.speciation_keys == {"CO2": "CO2"}

    def test_skips_reactions_with_mismatched_phase_keys(self):
        """A cross-phase reaction whose phase keys don't match the
        link's ``gas_phase_key`` / ``liquid_phase_key`` is ignored —
        defensive for multi-CV systems with multiple links.
        """
        from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry, ReactionSystem
        from PyOMES.chemistry.common_species import CO2

        rxn = EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=CO2, phase="headspace", coefficient=-1.0),
                StoichiometryEntry(species=CO2, phase="bulk", coefficient=+1.0),
            ],
            balance_elements=("C", "O"),
        )
        link = self._make_link()  # gas_phase_key="gas", liquid_phase_key="liquid"
        link.derive_speciation_keys(ReactionSystem([rxn]))
        # CO2 seeded from partition_models (Option A); mismatched reaction
        # contributed nothing via the cross-phase path.
        assert link.speciation_keys == {"CO2": "CO2"}

    def test_skips_kinetic_reactions(self):
        """Kinetic reactions must be ignored even if their
        stoichiometry happens to span phases (the engine's
        cross-phase semantics are scoped to equilibrium reactions).
        """
        from PyOMES.reactions import KineticReaction, StoichiometryEntry, ReactionSystem
        from PyOMES.chemistry.common_species import CO2

        rxn = KineticReaction(
            stoichiometry=[
                StoichiometryEntry(species=CO2, phase="gas", coefficient=-1.0),
                StoichiometryEntry(species=CO2, phase="liquid", coefficient=+1.0),
            ],
            rate_fn=lambda env: 0.0,
            balance_elements=("C", "O"),
        )
        link = self._make_link()
        link.derive_speciation_keys(ReactionSystem([rxn]))
        # CO2 seeded from partition_models (Option A); kinetic reaction
        # contributed nothing via the cross-phase path.
        assert link.speciation_keys == {"CO2": "CO2"}

    def test_additive_behavior_preserves_existing_entries(self):
        """``derive_speciation_keys`` adds to the existing dict; it
        does not overwrite pre-existing entries.  Option A also seeds
        any partition_models species not yet present.
        """
        from PyOMES.reactions import ReactionSystem
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"CO2": _hp(0.034), "NH3": _hp(58.0)},
            speciation_keys={"CO2": "CO2"},  # explicit pre-existing
        )
        # No cross-phase reactions declared; Option A seeds NH3 from partition_models.
        link.derive_speciation_keys(ReactionSystem([]))
        assert link.speciation_keys == {"CO2": "CO2", "NH3": "NH3"}

    def test_cv_init_hook_runs_derive(self):
        """``ControlVolume.__init__`` calls ``derive_speciation_keys``
        on every internal :class:`KineticGasLiquidLink` after the
        reaction model is attached.
        """
        from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry, ReactionSystem
        from PyOMES.chemistry.common_species import CO2

        rxn = EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=CO2, phase="gas", coefficient=-1.0),
                StoichiometryEntry(species=CO2, phase="liquid", coefficient=+1.0),
            ],
            balance_elements=("C", "O"),
        )
        link = KineticGasLiquidLink(
            gas_cv_key="cv", gas_phase_key="gas",
            liquid_cv_key="cv", liquid_phase_key="liquid",
            partition_models={"CO2": _hp(0.034)}, kLa={"CO2": 100.0},
            speciation_keys={},  # empty; the hook should populate
        )
        gas = GasPhase(n_mol={"CO2": 0.001}, V_L=0.4, T_K=305.15)
        liq = LiquidPhase(n_mol={"CO2": 0.01}, V_L=1.6, T_K=305.15)
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            internal_interfaces=[link],
            reaction_system=ReactionSystem([rxn]),
            label="test",
        )
        # Hook fires during __init__, populating keys and ladders.
        iface = cv.internal_interfaces[0]
        assert iface.speciation_keys == {"CO2": "CO2"}
        # Ladder for CO2 is a one-element list (no single-phase equilibria
        # declared for CO2 in this fixture — graceful degradation).
        assert iface.speciation_ladders.get("CO2") == ["CO2"]

    def test_partition_models_seed_speciation_keys_without_cross_phase_reaction(self):
        """Option A: a species in ``partition_models`` gets an identity
        mapping in ``speciation_keys`` even when no cross-phase
        EquilibriumReaction is declared.  This allows users to omit the
        routing-only boilerplate reaction for Henry-partitioned species.
        """
        from PyOMES.reactions import ReactionSystem
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={
                "CO2": _hp(0.034),
                "O2":  _hp(0.0013),
                "NH3": _hp(58.0),
            },
            speciation_keys={},  # start empty — no explicit declarations
        )
        link.derive_speciation_keys(ReactionSystem([]))  # no reactions at all
        # Every partition_models species gets an identity entry.
        assert link.speciation_keys == {"CO2": "CO2", "O2": "O2", "NH3": "NH3"}


# ═══════════════════════════════════════════════════════════════════════
#  H₂S alpha correction
# ═══════════════════════════════════════════════════════════════════════

class TestH2SAlphaCorrection:
    """Verify that the H₂S ionisation correction materially changes the
    gas-phase fraction.

    At pH ≈ pKa = 7.0: alpha₀ ≈ 0.5 (equal H₂S and HS⁻).
    System: V_liq=2 L, V_gas=0.5 L, T=298.15 K, 1 mol total H₂S.

    Without correction (alpha=1.0): ~9% in gas.
    With correction (alpha=0.5): ~5% in gas.
    """

    _kH_H2S = 0.10   # mol/L/atm at 298.15 K (Sander 2015)
    _V_liq = 2.0
    _V_gas = 0.5
    _T_K = 298.15

    def _link(self, alpha_override=None):
        """Build a link; set speciation_ladders if alpha_override given."""
        kH = self._kH_H2S
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"H2S": HenryPartition(
                H_ref=kH * 1000.0 / 101325.0, dlnH=0.0,
            )},
            equilibrium_species={"H2S"},
        )
        if alpha_override is not None:
            link.speciation_keys["H2S"] = "H2S"
            link.speciation_ladders["H2S"] = ["H2S", "HS-"]
        return link, alpha_override

    def _gas_fraction(self, alpha=1.0):
        R = 0.0820574
        T, V_liq, V_gas = self._T_K, self._V_liq, self._V_gas
        kH = self._kH_H2S
        beta = (kH / alpha) * R * T * V_liq / V_gas
        return 1.0 / (1.0 + beta)

    def test_alpha1_gives_approx_9pct_in_gas(self):
        """Without ionisation correction, ~9% of H₂S inventory is in gas."""
        frac = self._gas_fraction(alpha=1.0)
        assert frac == pytest.approx(0.09, abs=0.02)

    def test_alpha05_gives_approx_5pct_in_gas(self):
        """With alpha=0.5 (pH=pKa), ~5% of H₂S inventory is in gas."""
        frac = self._gas_fraction(alpha=0.5)
        assert frac == pytest.approx(0.05, abs=0.02)

    def test_alpha_correction_reduces_gas_fraction(self):
        """Model-driven: alpha=0.5 must give strictly less gas than alpha=1."""
        from PyOMES.core import LiquidPhase, GasPhase
        V_liq, V_gas, T_K = self._V_liq, self._V_gas, self._T_K
        n_total = 1.0

        # Without correction (no speciation ladder → alpha=1)
        hp = HenryPartition(H_ref=self._kH_H2S * 1000.0 / 101325.0, dlnH=0.0)
        n_liq_alpha1 = hp.equilibrium_a_moles(n_total, V_liq, V_gas, T_K, alpha=1.0)

        # With correction (pH=pKa → alpha=0.5)
        n_liq_alpha05 = hp.equilibrium_a_moles(n_total, V_liq, V_gas, T_K, alpha=0.5)

        gas_frac_alpha1 = (n_total - n_liq_alpha1) / n_total
        gas_frac_alpha05 = (n_total - n_liq_alpha05) / n_total

        assert gas_frac_alpha05 < gas_frac_alpha1
        assert gas_frac_alpha1 == pytest.approx(0.09, abs=0.02)
        assert gas_frac_alpha05 == pytest.approx(0.05, abs=0.02)
