# -*- coding: utf-8 -*-
"""Tests for configuration dataclasses (Stage D).

Validates:
- Default construction and derived properties
- Validation rejects invalid inputs
- to_dict / from_dict round-trip serialisation
- TransferConfig convenience constructors (kinetic, equilibrium)
- OrganismConfig / SubstrateConfig take an id, an id with atoms, or a Species
- TransferMode enum string parsing
- SpeciesTransferConfig from raw dict in TransferConfig
"""

import pytest
import copy

from PyOMES.chemistry import Species

from PyOMES.templates.stirred_tank import (
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

    def test_gas_composition_defaults_to_empty(self):
        assert VesselConfig().gas_composition == {}

    def test_gas_composition_kept_as_given(self):
        v = VesselConfig(gas_composition={"O2": 0.21, "Ar": 0.01})
        assert v.gas_composition == {"O2": 0.21, "Ar": 0.01}

    def test_negative_gas_fraction_raises(self):
        with pytest.raises(ValueError, match="gas_composition"):
            VesselConfig(gas_composition={"O2": -0.1})

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
        g = GasFeedConfig(composition={"N2": 1.0})
        assert g.vvm_min == 1.0
        assert g.P_inlet_atm == 1.0

    def test_feed_without_composition_raises(self):
        with pytest.raises(ValueError, match="composition"):
            GasFeedConfig()

    def test_zero_vvm_needs_no_composition(self):
        assert GasFeedConfig(vvm_min=0.0).composition == {}

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

    def test_kinetic_factory(self):
        t = TransferConfig.kinetic({"O2": 200.0, "CO2": 180.0}, equilibrium=["N2"])
        assert list(t.species) == ["O2", "CO2", "N2"]
        assert t.species["O2"].mode == TransferMode.KINETIC
        assert t.species["O2"].kLa_per_h == 200.0
        assert t.species["CO2"].kLa_per_h == 180.0
        assert t.species["N2"].mode == TransferMode.EQUILIBRIUM

    def test_equilibrium_factory(self):
        t = TransferConfig.equilibrium(["CH4", "H2"])
        assert list(t.species) == ["CH4", "H2"]
        for sp in ("CH4", "H2"):
            assert t.species[sp].mode == TransferMode.EQUILIBRIUM

    def test_default_is_no_transfer(self):
        assert TransferConfig().species == {}

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
        t1 = TransferConfig.kinetic({"O2": 250.0, "CO2": 225.0}, equilibrium=["N2"])
        d = t1.to_dict()
        t2 = TransferConfig.from_dict(d)
        assert t2.species["O2"].kLa_per_h == 250.0
        assert t2.species["CO2"].mode == TransferMode.KINETIC
        assert t2.species["N2"].mode == TransferMode.EQUILIBRIUM


# ═══════════════════════════════════════════════════════════════════════
#  ChemistryConfig
# ═══════════════════════════════════════════════════════════════════════

class TestChemistryConfig:

    def test_defaults(self):
        c = ChemistryConfig()
        assert c.activity_model == "ideal"

    def test_round_trip(self):
        c1 = ChemistryConfig(activity_model="davies")
        d = c1.to_dict()
        c2 = ChemistryConfig.from_dict(d)
        assert c2.activity_model == "davies"


# ═══════════════════════════════════════════════════════════════════════
#  OrganismConfig
# ═══════════════════════════════════════════════════════════════════════

class TestOrganismConfig:

    def test_organism_is_required(self):
        with pytest.raises(TypeError):
            OrganismConfig()

    def test_defaults(self):
        o = OrganismConfig("Yeast")
        assert o.organism == "Yeast"
        assert o.balance_basis == "CHO"
        assert o.atoms is None and o.MW is None
        assert o.n_source_id is None
        assert o.overwrite is False

    def test_explicit_atoms_and_mw(self):
        o = OrganismConfig(
            organism="Custom",
            atoms={"C": 1, "H": 2, "O": 1},
            MW=30.0,
            balance_basis="CHO",
        )
        assert o.atoms["C"] == 1
        assert o.MW == 30.0

    def test_species_accepted(self):
        sp = Species(id="E_coli", atoms={"C": 1, "H": 1.77, "O": 0.49, "N": 0.24})
        assert OrganismConfig(sp).organism is sp

    def test_species_with_atoms_raises(self):
        sp = Species(id="E_coli", atoms={"C": 1})
        with pytest.raises(ValueError, match="is a Species"):
            OrganismConfig(sp, atoms={"C": 1})

    def test_mw_without_atoms_raises(self):
        with pytest.raises(ValueError, match="needs atoms="):
            OrganismConfig("Yeast", MW=24.6)

    def test_non_id_raises(self):
        with pytest.raises(TypeError, match="species id"):
            OrganismConfig(42)

    def test_invalid_balance_raises(self):
        with pytest.raises(ValueError, match="balance_basis"):
            OrganismConfig("Yeast", balance_basis="XYZ")

    def test_normalises_balance_case(self):
        o = OrganismConfig("Yeast", balance_basis="chno", n_source_id="NH3")
        assert o.balance_basis == "CHNO"

    def test_chno_without_n_source_raises(self):
        with pytest.raises(ValueError, match="n_source_id"):
            OrganismConfig("Yeast", balance_basis="CHNO")

    def test_invalid_mw_raises(self):
        with pytest.raises(ValueError, match="MW"):
            OrganismConfig("Custom", atoms={"C": 1}, MW=-5.0)

    def test_round_trip(self):
        o1 = OrganismConfig(organism="E_coli", balance_basis="CHNO", n_source_id="NH3",
                            atoms={"C": 1, "H": 1.77, "O": 0.49, "N": 0.24}, MW=23.7,
                            o2_id="O2_aq", overwrite=True)
        d = o1.to_dict()
        o2 = OrganismConfig.from_dict(d)
        assert o2.organism == "E_coli"
        assert o2.balance_basis == "CHNO"
        assert o2.atoms["N"] == pytest.approx(0.24)
        assert (o2.n_source_id, o2.o2_id, o2.overwrite) == ("NH3", "O2_aq", True)

    def test_round_trip_keeps_species(self):
        sp = Species(id="E_coli", atoms={"C": 1, "H": 1.77, "O": 0.49})
        o2 = OrganismConfig.from_dict(OrganismConfig(sp).to_dict())
        assert o2.organism is sp


# ═══════════════════════════════════════════════════════════════════════
#  SubstrateConfig
# ═══════════════════════════════════════════════════════════════════════

class TestSubstrateConfig:

    def test_substrate_is_required(self):
        with pytest.raises(TypeError):
            SubstrateConfig()

    def test_defaults(self):
        s = SubstrateConfig("AceticAcid")
        assert s.substrate == "AceticAcid"
        assert s.mu_max == 0.5
        assert s.overwrite is False

    def test_explicit_atoms_and_mw(self):
        s = SubstrateConfig(
            substrate="Custom",
            atoms={"C": 6, "H": 12, "O": 6},
            MW=180.0,
            mu_max=0.8,
            Ks=0.02,
            yield_gX_gS=0.5,
        )
        assert s.MW == 180.0

    def test_species_accepted(self):
        sp = Species(id="Glucose", atoms={"C": 6, "H": 12, "O": 6})
        assert SubstrateConfig(sp).substrate is sp

    def test_mw_without_atoms_raises(self):
        with pytest.raises(ValueError, match="needs atoms="):
            SubstrateConfig("Glucose", MW=180.0)

    def test_invalid_mu_max_raises(self):
        with pytest.raises(ValueError, match="mu_max"):
            SubstrateConfig("AceticAcid", mu_max=0.0)

    def test_invalid_yield_raises(self):
        with pytest.raises(ValueError, match="yield_gX_gS"):
            SubstrateConfig("AceticAcid", yield_gX_gS=-0.1)

    def test_invalid_mw_raises(self):
        with pytest.raises(ValueError, match="MW"):
            SubstrateConfig("Custom", atoms={"C": 1}, MW=-5.0)

    def test_negative_ks_raises(self):
        with pytest.raises(ValueError, match="Ks"):
            SubstrateConfig("AceticAcid", Ks=-0.001)

    def test_round_trip(self):
        s1 = SubstrateConfig(substrate="Glucose", atoms={"C": 6, "H": 12, "O": 6},
                             MW=180.0, mu_max=0.8, Ks=0.02, yield_gX_gS=0.5)
        d = s1.to_dict()
        s2 = SubstrateConfig.from_dict(d)
        assert s2.substrate == "Glucose"
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
