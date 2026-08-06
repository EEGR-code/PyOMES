# -*- coding: utf-8 -*-
"""Tests for the TransferModel API.

Covers:
- KineticTransferModel and EquilibriumTransferModel construction and validation
- ControlVolume.transfer_models kwarg: auto phase-pair detection, link creation,
  kLa / equilibrium_species / molecular_driving_force routing
- transfer_basis="molecular" routing
- Error cases: bad transfer_basis, no gas/liquid pair
- Coexistence of transfer_models and internal_interfaces
- snapshot() round-trips correctly without double-adding the link
- partition_ratio() rename on HenryPartition and PartitionModel protocol
"""

import math
import warnings
import pytest

from PyOMES.core import (
    ControlVolume,
    EquilibriumTransferModel,
    GasPhase,
    KineticTransferModel,
    LiquidPhase,
    SolidPhase,
)
from PyOMES.core.gas_liquid_link import KineticGasLiquidLink
from PyOMES.chemistry import HenryEquilibrium


# ── Fixtures ──────────────────────────────────────────────────────────────────

def _henry(kH_mol_L_atm: float) -> HenryEquilibrium:
    return HenryEquilibrium(H_ref=kH_mol_L_atm * 1000.0 / 101325.0, dlnH=0.0)


def _simple_phases():
    gas = GasPhase({"O2": 0.01, "CO2": 0.005, "N2": 0.04}, V_L=0.4, T_K=305.15)
    liq = LiquidPhase({"O2": 1e-4, "CO2": 1e-4}, V_L=1.6, T_K=305.15)
    return gas, liq


# ── KineticTransferModel ──────────────────────────────────────────────────────

class TestKineticTransferModel:
    def test_construction_defaults(self):
        m = KineticTransferModel(_henry(1e-3), k_transfer=150.0)
        assert m.k_transfer == 150.0
        assert m.transfer_basis == "total"

    def test_construction_molecular_basis(self):
        m = KineticTransferModel(_henry(1e-3), k_transfer=200.0, transfer_basis="molecular")
        assert m.transfer_basis == "molecular"

    def test_k_transfer_coerced_to_float(self):
        m = KineticTransferModel(_henry(1e-3), k_transfer=100)
        assert isinstance(m.k_transfer, float)

    def test_invalid_transfer_basis_raises(self):
        with pytest.raises(ValueError, match="transfer_basis"):
            KineticTransferModel(_henry(1e-3), k_transfer=100.0, transfer_basis="volumetric")

    def test_partition_model_stored(self):
        hp = _henry(1.3e-3)
        m = KineticTransferModel(hp, k_transfer=50.0)
        assert m.partition_model is hp


# ── EquilibriumTransferModel ──────────────────────────────────────────────────

class TestEquilibriumTransferModel:
    def test_construction(self):
        hp = _henry(6.4e-4)
        m = EquilibriumTransferModel(hp)
        assert m.partition_model is hp

    def test_has_no_k_transfer(self):
        m = EquilibriumTransferModel(_henry(6.4e-4))
        assert not hasattr(m, "k_transfer")


# ── HenryPartition.partition_ratio rename ────────────────────────────────────

class TestPartitionRatioRename:
    def test_partition_ratio_exists(self):
        hp = _henry(1.3e-3)
        assert hasattr(hp, "partition_ratio")

    def test_beta_removed(self):
        hp = _henry(1.3e-3)
        assert not hasattr(hp, "beta")

    def test_partition_ratio_returns_float(self):
        hp = _henry(1.3e-3)
        r = hp.partition_ratio(1.6, 0.4, 305.15, alpha=1.0)
        assert isinstance(r, float)
        assert r > 0.0

    def test_equilibrium_a_moles_uses_partition_ratio(self):
        hp = _henry(1.3e-3)
        n_liq = hp.equilibrium_a_moles(1.0, 1.6, 0.4, 305.15)
        r = hp.partition_ratio(1.6, 0.4, 305.15)
        assert math.isclose(n_liq, r / (1.0 + r), rel_tol=1e-10)

    def test_alpha_correction(self):
        hp = _henry(3.4e-2)
        r_alpha = hp.partition_ratio(1.6, 0.4, 305.15, alpha=0.5)
        r_full = hp.partition_ratio(1.6, 0.4, 305.15, alpha=1.0)
        # lower alpha → higher ratio (more partitioned to liquid)
        assert r_alpha > r_full


# ── ControlVolume.transfer_models kwarg ──────────────────────────────────────

class TestControlVolumeTransferModels:
    def test_basic_construction_creates_link(self):
        gas, liq = _simple_phases()
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            transfer_models={
                "O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0),
            },
        )
        assert len(cv.internal_interfaces) == 1
        assert isinstance(cv.internal_interfaces[0], KineticGasLiquidLink)

    def test_transfer_models_dict_stored(self):
        gas, liq = _simple_phases()
        tm = {"O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0)}
        cv = ControlVolume(phases={"gas": gas, "liquid": liq}, transfer_models=tm)
        assert "O2" in cv.transfer_models
        assert cv.transfer_models["O2"] is tm["O2"]

    def test_kinetic_species_in_kLa(self):
        gas, liq = _simple_phases()
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            transfer_models={
                "O2":  KineticTransferModel(_henry(1.3e-3), k_transfer=150.0),
                "CO2": KineticTransferModel(_henry(2.94e-2), k_transfer=135.0),
            },
        )
        link = cv.internal_interfaces[0]
        assert link.kLa["O2"] == pytest.approx(150.0)
        assert link.kLa["CO2"] == pytest.approx(135.0)

    def test_equilibrium_species_routed_correctly(self):
        gas, liq = _simple_phases()
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            transfer_models={
                "O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0),
                "N2": EquilibriumTransferModel(_henry(5.3e-4)),
            },
        )
        link = cv.internal_interfaces[0]
        assert "N2" in link.equilibrium_species
        assert "O2" not in link.equilibrium_species

    def test_molecular_basis_routed_to_molecular_driving_force(self):
        gas, liq = _simple_phases()
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            transfer_models={
                "CO2": KineticTransferModel(
                    _henry(2.94e-2), k_transfer=135.0, transfer_basis="molecular"
                ),
            },
        )
        link = cv.internal_interfaces[0]
        assert "CO2" in link.molecular_driving_force

    def test_total_basis_not_in_molecular_driving_force(self):
        gas, liq = _simple_phases()
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            transfer_models={
                "O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0),
            },
        )
        link = cv.internal_interfaces[0]
        assert "O2" not in link.molecular_driving_force

    def test_no_transfer_models_no_link_appended(self):
        gas, liq = _simple_phases()
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})
        assert cv.internal_interfaces == []
        assert cv.transfer_models == {}

    def test_coexists_with_explicit_internal_interfaces(self):
        """transfer_models link appended after explicit internal_interfaces."""
        gas, liq = _simple_phases()
        from PyOMES.core.interfaces import PhaseInterface

        class DummyIface:
            phase_a_key = "gas"
            phase_b_key = "liquid"
            def compute_flux(self, a, b, dt_h):
                return {}

        dummy = DummyIface()
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            internal_interfaces=[dummy],
            transfer_models={"O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0)},
        )
        assert len(cv.internal_interfaces) == 2
        assert cv.internal_interfaces[0] is dummy
        assert isinstance(cv.internal_interfaces[1], KineticGasLiquidLink)

    def test_missing_gas_phase_raises(self):
        liq = LiquidPhase({"O2": 1e-4}, V_L=1.6, T_K=305.15)
        with pytest.raises(ValueError, match="GasPhase"):
            ControlVolume(
                phases={"liquid": liq},
                transfer_models={"O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0)},
            )

    def test_missing_liquid_phase_raises(self):
        gas = GasPhase({"O2": 0.01}, V_L=0.4, T_K=305.15)
        with pytest.raises(ValueError, match="LiquidPhase"):
            ControlVolume(
                phases={"gas": gas},
                transfer_models={"O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0)},
            )

    def test_phase_pair_override(self):
        gas = GasPhase({"O2": 0.01}, V_L=0.4, T_K=305.15)
        liq = LiquidPhase({"O2": 1e-4}, V_L=1.6, T_K=305.15)
        # Non-standard keys — auto-detection would find them, but explicit override works too
        cv = ControlVolume(
            phases={"headspace": gas, "broth": liq},
            transfer_models={"O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0)},
            phase_pair=("headspace", "broth"),
        )
        link = cv.internal_interfaces[0]
        assert link.gas_phase_key == "headspace"
        assert link.liquid_phase_key == "broth"

    def test_deprecation_warning_not_emitted_via_transfer_models(self):
        """Internal factory construction must not surface DeprecationWarning to users."""
        gas, liq = _simple_phases()
        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            # Should NOT raise — the factory suppresses the warning internally
            ControlVolume(
                phases={"gas": gas, "liquid": liq},
                transfer_models={"O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0)},
            )


# ── Snapshot round-trip ───────────────────────────────────────────────────────

class TestSnapshotWithTransferModels:
    def test_snapshot_has_one_link(self):
        gas, liq = _simple_phases()
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            transfer_models={"O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0)},
        )
        snap = cv.snapshot()
        assert len(snap.internal_interfaces) == 1
        assert isinstance(snap.internal_interfaces[0], KineticGasLiquidLink)

    def test_snapshot_transfer_models_dict_preserved(self):
        gas, liq = _simple_phases()
        tm = {"O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0)}
        cv = ControlVolume(phases={"gas": gas, "liquid": liq}, transfer_models=tm)
        snap = cv.snapshot()
        assert "O2" in snap.transfer_models

    def test_snapshot_explicit_interfaces_preserved(self):
        gas, liq = _simple_phases()

        class DummyIface:
            phase_a_key = "gas"
            phase_b_key = "liquid"
            def compute_flux(self, a, b, dt_h):
                return {}

        dummy = DummyIface()
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            internal_interfaces=[dummy],
            transfer_models={"O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0)},
        )
        snap = cv.snapshot()
        assert len(snap.internal_interfaces) == 2
        assert snap.internal_interfaces[0] is dummy


# ── Physics sanity: transfer actually moves moles ────────────────────────────

class TestTransferPhysics:
    def test_kinetic_transfer_moves_o2_gas_to_liquid(self):
        """Gas-to-liquid O2 absorption reduces gas moles and increases liquid moles."""
        gas = GasPhase({"O2": 0.02}, V_L=0.4, T_K=305.15)
        liq = LiquidPhase({"O2": 0.0}, V_L=1.6, T_K=305.15)
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            transfer_models={"O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0)},
        )
        n_gas_before = cv.phases["gas"].n_mol["O2"]
        cv.step_internal_transfer(dt_h=0.01)
        n_gas_after = cv.phases["gas"].n_mol["O2"]
        n_liq_after = cv.phases["liquid"].n_mol.get("O2", 0.0)
        assert n_gas_after < n_gas_before
        assert n_liq_after > 0.0

    def test_equilibrium_transfer_reaches_partition(self):
        """Equilibrium species reaches near-Henry partition in one step."""
        hp = _henry(5.3e-4)
        gas = GasPhase({"N2": 0.05}, V_L=0.4, T_K=305.15)
        liq = LiquidPhase({"N2": 0.0}, V_L=1.6, T_K=305.15)
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            transfer_models={"N2": EquilibriumTransferModel(hp)},
        )
        n_total = 0.05
        cv.step_internal_transfer(dt_h=0.01)
        n_liq = cv.phases["liquid"].n_mol.get("N2", 0.0)
        n_liq_eq = hp.equilibrium_a_moles(n_total, 1.6, 0.4, 305.15)
        assert abs(n_liq - n_liq_eq) < 1e-10

    def test_conservation_across_transfer(self):
        """Total moles conserved across a kinetic transfer step."""
        gas = GasPhase({"O2": 0.02}, V_L=0.4, T_K=305.15)
        liq = LiquidPhase({"O2": 0.001}, V_L=1.6, T_K=305.15)
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            transfer_models={"O2": KineticTransferModel(_henry(1.3e-3), k_transfer=150.0)},
        )
        n_before = cv.total_mol().get("O2", 0.0)
        diag = cv.step_internal_transfer(dt_h=0.01)
        n_after = cv.total_mol().get("O2", 0.0)
        assert abs(n_after - n_before) < 1e-12
        assert diag.is_conserved()
