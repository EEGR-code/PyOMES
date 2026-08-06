# -*- coding: utf-8 -*-
"""Tests for configuration dataclasses (Stage D).

Validates:
- Default construction and derived properties
- Validation rejects invalid inputs
- to_dict / from_dict round-trip serialisation
- TransferConfig convenience constructors (default_kinetic, default_equilibrium)
- OrganismConfig.resolve() and SubstrateConfig.resolve() from registry
- Explicit atoms/MW override registry lookup
- TransferMode enum string parsing
- SpeciesTransferConfig from raw dict in TransferConfig
"""

import pytest
import copy

from vlmodels.fermenter.config import (
    TransferMode,
    VesselConfig,
    GasFeedConfig,
    SpeciesTransferConfig,
    TransferConfig,
    ChemistryConfig,
    OrganismConfig,
    SubstrateConfig,
    SimulationConfig,
)


# ═══════════════════════════════════════════════════════════════════════
#  VesselConfig
# ═══════════════════════════════════════════════════════════════════════

class TestVesselConfig:

    def test_defaults(self):
        v = VesselConfig()
        assert v.V_total_L == 2.0
        assert v.headspace_frac == 0.20
        assert v.T_K == 305.15

    def test_derived_volumes(self):
        v = VesselConfig(V_total_L=100.0, headspace_frac=0.25)
        assert v.V_headspace_L == pytest.approx(25.0)
        assert v.V_liquid_L == pytest.approx(75.0)

    def test_yN2_computed_from_balance(self):
        v = VesselConfig(yO2_init=0.21, yCO2_init=0.04)
        assert v.yN2_init == pytest.approx(0.75, abs=1e-10)

    def test_yN2_explicit(self):
        v = VesselConfig(yO2_init=0.21, yCO2_init=0.04, yN2_init=0.5)
        assert v.yN2_init == 0.5

    def test_invalid_volume_raises(self):
        with pytest.raises(ValueError, match="V_total_L"):
            VesselConfig(V_total_L=-1.0)

    def test_invalid_headspace_frac_raises(self):
        with pytest.raises(ValueError, match="headspace_frac"):
            VesselConfig(headspace_frac=0.0)
        with pytest.raises(ValueError, match="headspace_frac"):
            VesselConfig(headspace_frac=1.0)

    def test_invalid_temperature_raises(self):
        with pytest.raises(ValueError, match="T_K"):
            VesselConfig(T_K=-10)

    def test_invalid_pressure_raises(self):
        with pytest.raises(ValueError, match="P_init_atm"):
            VesselConfig(P_init_atm=0)

    def test_round_trip(self):
        v1 = VesselConfig(V_total_L=500, headspace_frac=0.15, T_K=310.0)
        d = v1.to_dict()
        v2 = VesselConfig.from_dict(d)
        assert v2.V_total_L == v1.V_total_L
        assert v2.headspace_frac == v1.headspace_frac
        assert v2.T_K == v1.T_K

    def test_from_dict_ignores_extra_keys(self):
        d = {"V_total_L": 100, "headspace_frac": 0.3, "T_K": 300, "extra": "ignored"}
        v = VesselConfig.from_dict(d)
        assert v.V_total_L == 100


# ═══════════════════════════════════════════════════════════════════════
#  GasFeedConfig
# ═══════════════════════════════════════════════════════════════════════

class TestGasFeedConfig:

    def test_defaults(self):
        g = GasFeedConfig()
        assert g.vvm_min == 1.0
        assert g.P_inlet_atm == 1.0

    def test_composition_normalised(self):
        g = GasFeedConfig(composition={"O2": 1, "N2": 3})
        assert g.composition["O2"] == pytest.approx(0.25)
        assert g.composition["N2"] == pytest.approx(0.75)

    def test_zero_vvm_for_no_sparging(self):
        g = GasFeedConfig(vvm_min=0.0)
        assert g.vvm_min == 0.0

    def test_negative_vvm_raises(self):
        with pytest.raises(ValueError, match="vvm_min"):
            GasFeedConfig(vvm_min=-1.0)

    def test_invalid_pressure_raises(self):
        with pytest.raises(ValueError, match="P_inlet_atm"):
            GasFeedConfig(P_inlet_atm=0.0)

    def test_round_trip(self):
        g1 = GasFeedConfig(vvm_min=2.0, composition={"O2": 0.5, "N2": 0.5})
        d = g1.to_dict()
        g2 = GasFeedConfig.from_dict(d)
        assert g2.vvm_min == g1.vvm_min
        assert g2.composition["O2"] == pytest.approx(g1.composition["O2"])


# ═══════════════════════════════════════════════════════════════════════
#  SpeciesTransferConfig
# ═══════════════════════════════════════════════════════════════════════

class TestSpeciesTransferConfig:

    def test_default_equilibrium(self):
        s = SpeciesTransferConfig()
        assert s.mode == TransferMode.EQUILIBRIUM
        assert s.kLa_per_h == 0.0

    def test_string_mode_parsing(self):
        s = SpeciesTransferConfig(mode="kinetic")
        assert s.mode == TransferMode.KINETIC

    def test_case_insensitive_mode(self):
        s = SpeciesTransferConfig(mode="EQUILIBRIUM")
        assert s.mode == TransferMode.EQUILIBRIUM

    def test_invalid_mode_raises(self):
        with pytest.raises(ValueError):
            SpeciesTransferConfig(mode="invalid")

    def test_negative_kla_raises(self):
        with pytest.raises(ValueError, match="kLa_per_h"):
            SpeciesTransferConfig(kLa_per_h=-1.0)

    def test_to_dict_serialises_mode_as_string(self):
        s = SpeciesTransferConfig(mode=TransferMode.KINETIC, kLa_per_h=100.0)
        d = s.to_dict()
        assert d["mode"] == "kinetic"
        assert d["kLa_per_h"] == 100.0


# ═══════════════════════════════════════════════════════════════════════
#  TransferConfig
# ═══════════════════════════════════════════════════════════════════════

class TestTransferConfig:

    def test_default_kinetic_factory(self):
        t = TransferConfig.default_kinetic(kLa_O2=200.0)
        assert t.species["O2"].mode == TransferMode.KINETIC
        assert t.species["O2"].kLa_per_h == 200.0
        assert t.species["CO2"].mode == TransferMode.KINETIC
        assert t.species["CO2"].kLa_per_h == pytest.approx(180.0)
        assert t.species["N2"].mode == TransferMode.EQUILIBRIUM

    def test_default_equilibrium_factory(self):
        t = TransferConfig.default_equilibrium()
        for sp in ("O2", "CO2", "N2"):
            assert t.species[sp].mode == TransferMode.EQUILIBRIUM

    def test_custom_co2_ratio(self):
        t = TransferConfig.default_kinetic(kLa_O2=100.0, kLa_CO2_ratio=0.8)
        assert t.species["CO2"].kLa_per_h == pytest.approx(80.0)
        assert t.kLa_CO2_ratio == 0.8

    def test_accepts_raw_dicts(self):
        t = TransferConfig(
            species={
                "O2": {"mode": "kinetic", "kLa_per_h": 150.0},
                "N2": {"mode": "equilibrium"},
            }
        )
        assert isinstance(t.species["O2"], SpeciesTransferConfig)
        assert t.species["O2"].mode == TransferMode.KINETIC

    def test_invalid_species_entry_raises(self):
        with pytest.raises(TypeError, match="SpeciesTransferConfig"):
            TransferConfig(species={"O2": 42})

    def test_round_trip(self):
        t1 = TransferConfig.default_kinetic(kLa_O2=250.0)
        d = t1.to_dict()
        t2 = TransferConfig.from_dict(d)
        assert t2.species["O2"].kLa_per_h == 250.0
        assert t2.species["CO2"].mode == TransferMode.KINETIC
        assert t2.kLa_CO2_ratio == t1.kLa_CO2_ratio


# ═══════════════════════════════════════════════════════════════════════
#  ChemistryConfig
# ═══════════════════════════════════════════════════════════════════════

class TestChemistryConfig:

    def test_defaults(self):
        c = ChemistryConfig()
        assert c.use_activity is False
        assert "AceticAcid" in c.acid_pKas

    def test_round_trip(self):
        c1 = ChemistryConfig(use_activity=True)
        d = c1.to_dict()
        c2 = ChemistryConfig.from_dict(d)
        assert c2.use_activity is True


# ═══════════════════════════════════════════════════════════════════════
#  OrganismConfig
# ═══════════════════════════════════════════════════════════════════════

class TestOrganismConfig:

    def test_defaults(self):
        o = OrganismConfig()
        assert o.organism_id == "Yeast"
        assert o.balance_basis == "CHO"
        assert o.atoms is None  # resolved later

    def test_explicit_atoms_and_mw(self):
        o = OrganismConfig(
            organism_id="Custom",
            atoms={"C": 1, "H": 2, "O": 1},
            MW=30.0,
            balance_basis="CHO",
        )
        assert o.atoms["C"] == 1
        assert o.MW == 30.0

    def test_invalid_balance_raises(self):
        with pytest.raises(ValueError, match="balance_basis"):
            OrganismConfig(balance_basis="XYZ")

    def test_normalises_balance_case(self):
        o = OrganismConfig(balance_basis="chno")
        assert o.balance_basis == "CHNO"

    def test_invalid_mw_raises(self):
        with pytest.raises(ValueError, match="MW"):
            OrganismConfig(MW=-5.0)

    def test_resolve_from_default_registry(self):
        o = OrganismConfig(organism_id="Yeast")
        resolved = o.resolve()
        assert resolved.atoms is not None
        assert resolved.atoms["C"] == pytest.approx(1.0)
        assert resolved.MW is not None
        assert resolved.MW > 0

    def test_resolve_keeps_explicit_atoms(self):
        o = OrganismConfig(
            organism_id="Yeast",
            atoms={"C": 99, "H": 99},
            MW=999.0,
        )
        resolved = o.resolve()
        assert resolved.atoms["C"] == 99
        assert resolved.MW == 999.0

    def test_resolve_unknown_organism_without_atoms_raises(self):
        o = OrganismConfig(organism_id="UnknownBug")
        with pytest.raises(KeyError):
            o.resolve()

    def test_round_trip(self):
        o1 = OrganismConfig(organism_id="E_coli", balance_basis="CHNO",
                            atoms={"C": 1, "H": 1.77, "O": 0.49, "N": 0.24}, MW=23.7)
        d = o1.to_dict()
        o2 = OrganismConfig.from_dict(d)
        assert o2.organism_id == "E_coli"
        assert o2.balance_basis == "CHNO"
        assert o2.atoms["N"] == pytest.approx(0.24)


# ═══════════════════════════════════════════════════════════════════════
#  SubstrateConfig
# ═══════════════════════════════════════════════════════════════════════

class TestSubstrateConfig:

    def test_defaults(self):
        s = SubstrateConfig()
        assert s.substrate_id == "AceticAcid"
        assert s.mu_max == 0.5

    def test_explicit_atoms_and_mw(self):
        s = SubstrateConfig(
            substrate_id="Custom",
            atoms={"C": 6, "H": 12, "O": 6},
            MW=180.0,
            mu_max=0.8,
            Ks=0.02,
            yield_gX_gS=0.5,
        )
        assert s.MW == 180.0

    def test_invalid_mu_max_raises(self):
        with pytest.raises(ValueError, match="mu_max"):
            SubstrateConfig(mu_max=0.0)

    def test_invalid_yield_raises(self):
        with pytest.raises(ValueError, match="yield_gX_gS"):
            SubstrateConfig(yield_gX_gS=-0.1)

    def test_invalid_mw_raises(self):
        with pytest.raises(ValueError, match="MW"):
            SubstrateConfig(MW=-5.0)

    def test_negative_ks_raises(self):
        with pytest.raises(ValueError, match="Ks"):
            SubstrateConfig(Ks=-0.001)

    def test_resolve_from_registry(self):
        s = SubstrateConfig(substrate_id="AceticAcid")
        resolved = s.resolve()
        assert resolved.atoms is not None
        assert resolved.atoms["C"] == pytest.approx(2.0)
        assert resolved.MW == pytest.approx(60.052)
        # Kinetics should be preserved
        assert resolved.mu_max == s.mu_max
        assert resolved.yield_gX_gS == s.yield_gX_gS

    def test_resolve_keeps_explicit(self):
        s = SubstrateConfig(
            substrate_id="AceticAcid",
            atoms={"C": 99},
            MW=999.0,
        )
        resolved = s.resolve()
        assert resolved.atoms["C"] == 99
        assert resolved.MW == 999.0

    def test_resolve_unknown_without_atoms_raises(self):
        s = SubstrateConfig(substrate_id="UnknownSubstrate")
        with pytest.raises(KeyError):
            s.resolve()

    def test_round_trip(self):
        s1 = SubstrateConfig(substrate_id="Glucose", atoms={"C": 6, "H": 12, "O": 6},
                             MW=180.0, mu_max=0.8, Ks=0.02, yield_gX_gS=0.5)
        d = s1.to_dict()
        s2 = SubstrateConfig.from_dict(d)
        assert s2.substrate_id == "Glucose"
        assert s2.MW == 180.0
        assert s2.mu_max == 0.8


# ═══════════════════════════════════════════════════════════════════════
#  SimulationConfig
# ═══════════════════════════════════════════════════════════════════════

class TestSimulationConfig:

    def test_defaults(self):
        s = SimulationConfig()
        assert s.tau == 5.0
        assert s.n_steps == 1000
        assert s.t_lag == 0.0

    def test_dt_property(self):
        s = SimulationConfig(tau=10.0, n_steps=500)
        assert s.dt_h == pytest.approx(0.02)

    def test_invalid_tau_raises(self):
        with pytest.raises(ValueError, match="tau"):
            SimulationConfig(tau=0)

    def test_invalid_nsteps_raises(self):
        with pytest.raises(ValueError, match="n_steps"):
            SimulationConfig(n_steps=0)

    def test_negative_lag_raises(self):
        with pytest.raises(ValueError, match="t_lag"):
            SimulationConfig(t_lag=-1.0)

    def test_round_trip(self):
        s1 = SimulationConfig(tau=24, n_steps=2000, t_lag=0.5)
        d = s1.to_dict()
        s2 = SimulationConfig.from_dict(d)
        assert s2.tau == 24
        assert s2.n_steps == 2000
        assert s2.t_lag == 0.5


# ═══════════════════════════════════════════════════════════════════════
#  TransferMode enum
# ═══════════════════════════════════════════════════════════════════════

class TestTransferMode:

    def test_values(self):
        assert TransferMode.KINETIC.value == "kinetic"
        assert TransferMode.EQUILIBRIUM.value == "equilibrium"
        assert TransferMode.NONE.value == "none"

    def test_string_constructible(self):
        assert TransferMode("kinetic") == TransferMode.KINETIC
        assert TransferMode("equilibrium") == TransferMode.EQUILIBRIUM

    def test_is_string(self):
        # str(Enum) subclass → can be used as string
        assert TransferMode.KINETIC == "kinetic"
