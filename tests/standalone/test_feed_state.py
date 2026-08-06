"""Tests for fermenter.stream_adapter — FeedState."""

import pytest
from PyOMES.stream_adapter import FeedState
from PyOMES.chemistry.compounds import ChemicalRegistry


class TestFeedStateConstruction:
    """Test all FeedState construction paths."""

    def test_from_mass_concentrations(self, registry):
        feed = FeedState.from_mass_concentrations(
            {"AceticAcid": 60.052},  # exactly 1 mol/L
            registry=registry,
        )
        assert feed.mol_L("AceticAcid") == pytest.approx(1.0, rel=1e-6)
        assert feed.g_L("AceticAcid") == pytest.approx(60.052, rel=1e-6)

    def test_from_molar_direct(self, registry):
        feed = FeedState(
            concentrations_mol_L={"AceticAcid": 0.5},
            registry=registry,
        )
        assert feed.mol_L("AceticAcid") == pytest.approx(0.5)
        assert feed.g_L("AceticAcid") == pytest.approx(0.5 * 60.052, rel=1e-4)

    def test_from_mixed_concentrations(self, registry):
        feed = FeedState.from_mixed_concentrations(
            mass_g_L={"AceticAcid": 1.0},
            molar_mol_L={"NH3": 0.01},
            registry=registry,
        )
        assert feed.mol_L("NH3") == pytest.approx(0.01)
        assert feed.g_L("AceticAcid") == pytest.approx(1.0, rel=1e-3)

    def test_zero_concentration_excluded(self, registry):
        feed = FeedState.from_mass_concentrations(
            {"AceticAcid": 1.0, "Yeast": 0.0},
            registry=registry,
        )
        assert feed.has("AceticAcid")
        assert not feed.has("Yeast")

    def test_unknown_compound_raises(self, registry):
        with pytest.raises(KeyError, match="not found"):
            FeedState.from_mass_concentrations(
                {"FakeCompound": 1.0}, registry=registry,
            )

    def test_default_registry_created_automatically(self):
        feed = FeedState.from_mass_concentrations({"AceticAcid": 1.0})
        assert feed.registry is not None
        assert feed.g_L("AceticAcid") == pytest.approx(1.0, rel=1e-3)


class TestFeedStateAccessors:
    """Test mol_L, g_L, has, species_ids."""

    def test_mol_L_missing_species_returns_zero(self, simple_feed):
        assert simple_feed.mol_L("CO2") == 0.0

    def test_g_L_missing_species_returns_zero(self, simple_feed):
        assert simple_feed.g_L("CO2") == 0.0

    def test_has_returns_false_for_absent(self, simple_feed):
        assert not simple_feed.has("CO2")

    def test_species_ids(self, simple_feed):
        ids = simple_feed.species_ids
        assert "AceticAcid" in ids
        assert "Yeast" in ids

    def test_mol_to_g_roundtrip(self, registry):
        """g/L → mol/L → g/L should roundtrip exactly."""
        feed = FeedState.from_mass_concentrations(
            {"AceticAcid": 5.0, "NH3": 0.5},
            registry=registry,
        )
        assert feed.g_L("AceticAcid") == pytest.approx(5.0, rel=1e-6)
        assert feed.g_L("NH3") == pytest.approx(0.5, rel=1e-6)


class TestFeedStateMixedAdditive:
    """If a species appears in both g/L and mol/L inputs, they should add."""

    def test_mixed_adds_contributions(self, registry):
        mw = registry["NH3"].MW
        feed = FeedState.from_mixed_concentrations(
            mass_g_L={"NH3": mw},        # = 1.0 mol/L
            molar_mol_L={"NH3": 0.5},    # + 0.5 mol/L
            registry=registry,
        )
        assert feed.mol_L("NH3") == pytest.approx(1.5, rel=1e-6)


