# -*- coding: utf-8 -*-
"""Tests for FermenterBuilder fluent API (Stage F).

Validates:
- Minimal build (defaults only)
- Full chain: vessel → gas_feed → transfer → chemistry → organism → substrate → build
- build() returns ControlVolume
- build_simulation_and_run() returns BatchResult
- Multiple substrates via repeated .substrate() calls
- .controller() accumulates controllers
- .reaction_system() overrides organism/substrates
- .no_gas_feed() disables gas feed
- .transfer_kinetic() and .transfer_equilibrium() convenience methods
- .transfer() with custom TransferConfig
- .label() sets the label
- .get_configs() returns config dataclasses without building
- Warning when substrates set without organism
- Builder is chainable (every method returns self)
- repr shows current state
- Identical output to direct FermenterFactory.create_volume call
"""

import warnings
import numpy as np
import pytest

from vlmodels.fermenter.config import (
    FermenterBuilder,
    FermenterFactory,
    VesselConfig,
    GasFeedConfig,
    TransferConfig,
    ChemistryConfig,
    OrganismConfig,
    SubstrateConfig,
)
from PyOMES.core.control_volume import ControlVolume
from PyOMES.core.recorder import BatchResult
from PyOMES.core.gas_liquid_link import KineticGasLiquidLink


def _gl_link(cv):
    """Find the gas-liquid transfer link in a fermenter CV."""
    return next(
        i for i in cv.internal_interfaces
        if isinstance(i, KineticGasLiquidLink)
    )


# ═══════════════════════════════════════════════════════════════════════
#  Basic construction
# ═══════════════════════════════════════════════════════════════════════

class TestBuildBasic:

    def test_minimal_build(self):
        """Building with no calls should produce a ControlVolume with defaults."""
        cv = FermenterBuilder().build()
        assert isinstance(cv, ControlVolume)

    def test_full_chain(self):
        cv = (
            FermenterBuilder()
            .vessel(V_total_L=2000, T_K=305.15)
            .gas_feed(vvm_min=1.0, composition={"O2": 0.21, "N2": 0.79})
            .transfer_kinetic(kLa_O2=150.0)
            .chemistry()
            .organism("Yeast")
            .substrate("AceticAcid", mu_max=0.5, Ks=5e-3, yield_gX_gS=0.36)
            .label("my_fermenter")
            .build()
        )
        assert isinstance(cv, ControlVolume)
        assert cv.label == "my_fermenter"

    def test_build_returns_cv(self):
        cv = FermenterBuilder().vessel(V_total_L=100).build()
        assert isinstance(cv, ControlVolume)

    def test_build_simulation_and_run_returns_batch_result(self):
        result = (
            FermenterBuilder()
            .vessel(V_total_L=100, T_K=305.15)
            .transfer_equilibrium()
            .build_simulation_and_run(tau_h=0.1, n_steps=5)
        )
        assert isinstance(result, BatchResult)
        assert len(result.t_h) == 6


# ═══════════════════════════════════════════════════════════════════════
#  Vessel
# ═══════════════════════════════════════════════════════════════════════

class TestVessel:

    def test_vessel_params_forwarded(self):
        cv = (
            FermenterBuilder()
            .vessel(V_total_L=500, headspace_frac=0.30, T_K=310.0)
            .build()
        )
        assert cv.phases["gas"].V_L == pytest.approx(150.0)
        assert cv.phases["liquid"].V_L == pytest.approx(350.0)
        assert cv.phases["gas"].T_K == pytest.approx(310.0)

    def test_default_vessel(self):
        cv = FermenterBuilder().build()
        # Default VesselConfig: V_total=2.0, headspace=0.20
        assert cv.phases["gas"].V_L == pytest.approx(0.4)
        assert cv.phases["liquid"].V_L == pytest.approx(1.6)


# ═══════════════════════════════════════════════════════════════════════
#  Gas feed
# ═══════════════════════════════════════════════════════════════════════

class TestGasFeed:

    def test_gas_feed_creates_boundary(self):
        cv = (
            FermenterBuilder()
            .gas_feed(vvm_min=1.0, composition={"O2": 0.21, "N2": 0.79})
            .build()
        )
        assert len(cv.boundaries) == 1

    def test_no_gas_feed(self):
        cv = FermenterBuilder().no_gas_feed().build()
        assert len(cv.boundaries) == 0

    def test_no_gas_feed_by_default(self):
        cv = FermenterBuilder().build()
        assert len(cv.boundaries) == 0


# ═══════════════════════════════════════════════════════════════════════
#  Transfer
# ═══════════════════════════════════════════════════════════════════════

class TestTransfer:

    def test_kinetic(self):
        cv = FermenterBuilder().transfer_kinetic(kLa_O2=200.0).build()
        link = _gl_link(cv)
        assert link.kLa.get("O2", 0.0) == pytest.approx(200.0)
        assert "N2" in link.equilibrium_species

    def test_equilibrium(self):
        cv = FermenterBuilder().transfer_equilibrium().build()
        link = _gl_link(cv)
        assert "O2" in link.equilibrium_species
        assert "CO2" in link.equilibrium_species

    def test_custom_transfer(self):
        cfg = TransferConfig.default_kinetic(kLa_O2=999.0, kLa_CO2_ratio=0.5)
        cv = FermenterBuilder().transfer(cfg).build()
        assert _gl_link(cv).kLa["O2"] == pytest.approx(999.0)
        assert _gl_link(cv).kLa["CO2"] == pytest.approx(499.5)

    def test_default_transfer_is_equilibrium(self):
        cv = FermenterBuilder().build()
        link = _gl_link(cv)
        assert "O2" in link.equilibrium_species


# ═══════════════════════════════════════════════════════════════════════
#  Chemistry
# ═══════════════════════════════════════════════════════════════════════

class TestChemistry:

    def test_chemistry_forwarded(self):
        cv = (
            FermenterBuilder()
            .chemistry(use_activity=True, activity_model="davies")
            .build()
        )
        # state-unification C4d: chemistry config now configures the
        # lazy-build engine defaults on cv.reaction_system; the legacy
        # property_solvers list is gone. Without declared equilibria
        # the engine stays unbuilt (no chemistry to solve).
        assert cv.reaction_system is None or cv.reaction_system.engine is None

# test_custom_pkas was removed in chemistry-unification-1.
# `FermenterBuilder.chemistry()` no longer accepts `acid_pKas` — pKa
# values come from declared equilibrium reactions consumed by
# `BisectionChemicalEquilibriumEngine.from_reactions()`.


# ═══════════════════════════════════════════════════════════════════════
#  Organism and substrates
# ═══════════════════════════════════════════════════════════════════════

class TestOrganismSubstrates:

    def test_single_substrate(self):
        cv = (
            FermenterBuilder()
            .organism("Yeast")
            .substrate("AceticAcid", yield_gX_gS=0.36)
            .build()
        )
        assert cv.reaction_system is not None

    def test_multiple_substrates(self):
        from PyOMES.reactions import ReactionSystem
        cv = (
            FermenterBuilder()
            .organism("Yeast")
            .substrate("AceticAcid", yield_gX_gS=0.36)
            .substrate("PropionicAcid", yield_gX_gS=0.54)
            .build()
        )
        assert isinstance(cv.reaction_system, ReactionSystem)
        assert len(cv.reaction_system.reactions) == 2

    def test_explicit_organism_atoms(self):
        cv = (
            FermenterBuilder()
            .organism("Custom", atoms={"C": 1, "H": 2, "O": 1}, MW=30.0)
            .substrate("AceticAcid", yield_gX_gS=0.36)
            .build()
        )
        assert cv.reaction_system is not None

    def test_no_organism_no_reaction(self):
        cv = FermenterBuilder().build()
        assert cv.reaction_system is None

    def test_substrates_without_organism_warns(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cv = (
                FermenterBuilder()
                .substrate("AceticAcid")
                .build()
            )
            assert any("organism" in str(warning.message).lower() for warning in w)
        assert cv.reaction_system is None

    def test_chno_balance(self):
        cv = (
            FermenterBuilder()
            .organism("Yeast", balance_basis="CHNO")
            .substrate("AceticAcid", yield_gX_gS=0.36)
            .build()
        )
        assert cv.reaction_system is not None


# ═══════════════════════════════════════════════════════════════════════
#  Custom reaction model
# ═══════════════════════════════════════════════════════════════════════

class TestCustomReactionModel:

    def test_reaction_system_overrides(self):
        class MockRxn:
            def compute_rates(self, env):
                return {}

        mock = MockRxn()
        cv = (
            FermenterBuilder()
            .organism("Yeast")
            .substrate("AceticAcid")
            .reaction_system(mock)
            .build()
        )
        assert cv.reaction_system is mock


# ═══════════════════════════════════════════════════════════════════════
#  Controllers
# ═══════════════════════════════════════════════════════════════════════

class TestControllers:

    def test_controller_accumulates(self):
        b = (
            FermenterBuilder()
            .controller("ctrl_a")
            .controller("ctrl_b")
        )
        assert len(b._controllers) == 2


# ═══════════════════════════════════════════════════════════════════════
#  Label
# ═══════════════════════════════════════════════════════════════════════

class TestLabel:

    def test_label_set(self):
        cv = FermenterBuilder().label("my_well_plate").build()
        assert cv.label == "my_well_plate"


# ═══════════════════════════════════════════════════════════════════════
#  Context
# ═══════════════════════════════════════════════════════════════════════

# TestChemEnv deleted in state-unification C4e: the chem_env /
# chem_env_fn paths are removed entirely. Strong ions live in
# phase.n_mol (PHREEQC convention); pH derives from n_mol["H+"]
# populated by cv.reaction_system.engine.

# ═══════════════════════════════════════════════════════════════════════
#  get_configs
# ═══════════════════════════════════════════════════════════════════════

class TestGetConfigs:

    def test_returns_config_dict(self):
        cfgs = (
            FermenterBuilder()
            .vessel(V_total_L=500, T_K=310.0)
            .gas_feed(vvm_min=2.0)
            .transfer_kinetic(kLa_O2=200.0)
            .chemistry()
            .organism("Yeast")
            .substrate("AceticAcid", yield_gX_gS=0.36)
            .get_configs()
        )
        assert isinstance(cfgs["vessel"], VesselConfig)
        assert cfgs["vessel"].V_total_L == 500
        assert isinstance(cfgs["gas_feed"], GasFeedConfig)
        assert cfgs["gas_feed"].vvm_min == 2.0
        assert isinstance(cfgs["transfer"], TransferConfig)
        assert isinstance(cfgs["chemistry"], ChemistryConfig)
        assert isinstance(cfgs["organism"], OrganismConfig)
        assert cfgs["organism"].organism_id == "Yeast"
        assert len(cfgs["substrates"]) == 1
        assert isinstance(cfgs["substrates"][0], SubstrateConfig)


# ═══════════════════════════════════════════════════════════════════════
#  Chaining
# ═══════════════════════════════════════════════════════════════════════

class TestChaining:

    def test_every_method_returns_self(self):
        b = FermenterBuilder()
        assert b.vessel() is b
        assert b.gas_feed() is b
        assert b.no_gas_feed() is b
        assert b.transfer_kinetic() is b
        assert b.transfer_equilibrium() is b
        assert b.transfer(TransferConfig.default_equilibrium()) is b
        assert b.chemistry() is b
        assert b.organism() is b
        assert b.substrate() is b
        assert b.controller("x") is b
        assert b.reaction_system(None) is b
        assert b.label("x") is b
        # state-unification C4e: chem_env / chem_env_fn methods
        # removed from FermenterBuilder.


# ═══════════════════════════════════════════════════════════════════════
#  Repr
# ═══════════════════════════════════════════════════════════════════════

class TestRepr:

    def test_repr_empty(self):
        r = repr(FermenterBuilder())
        assert "empty" in r

    def test_repr_with_state(self):
        b = (
            FermenterBuilder()
            .vessel(V_total_L=2000)
            .organism("Yeast")
            .substrate("AceticAcid")
        )
        r = repr(b)
        assert "2000" in r
        assert "Yeast" in r
        assert "AceticAcid" in r


# ═══════════════════════════════════════════════════════════════════════
#  Parity with direct factory call
# ═══════════════════════════════════════════════════════════════════════

class TestParityWithFactory:

    def test_builder_matches_factory(self):
        """Builder and factory should produce equivalent ControlVolumes."""
        # Builder path
        cv_builder = (
            FermenterBuilder()
            .vessel(V_total_L=100, headspace_frac=0.25, T_K=305.15)
            .gas_feed(vvm_min=1.0, composition={"O2": 0.21, "N2": 0.79})
            .transfer_kinetic(kLa_O2=150.0, kLa_CO2_ratio=0.9)
            .chemistry()
            .build()
        )

        # Factory path
        cv_factory = FermenterFactory.create_volume(
            vessel=VesselConfig(V_total_L=100, headspace_frac=0.25, T_K=305.15),
            gas_feed=GasFeedConfig(vvm_min=1.0, composition={"O2": 0.21, "N2": 0.79}),
            transfer=TransferConfig.default_kinetic(kLa_O2=150.0, kLa_CO2_ratio=0.9),
            chemistry=ChemistryConfig(),
        )

        # Compare phase inventories
        for sp in cv_builder.phases["gas"].n_mol:
            assert cv_builder.phases["gas"].n_mol[sp] == pytest.approx(
                cv_factory.phases["gas"].n_mol[sp], abs=1e-12
            )
        for sp in cv_builder.phases["liquid"].n_mol:
            assert cv_builder.phases["liquid"].n_mol[sp] == pytest.approx(
                cv_factory.phases["liquid"].n_mol[sp], abs=1e-12
            )

        # Compare transfer link kLa
        assert _gl_link(cv_builder).kLa["O2"] == _gl_link(cv_factory).kLa["O2"]
        assert _gl_link(cv_builder).kLa["CO2"] == _gl_link(cv_factory).kLa["CO2"]
