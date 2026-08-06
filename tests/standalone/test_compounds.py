"""Tests for fermenter.chemistry.compounds — standalone chemical registry."""

import pytest
from PyOMES.chemistry.compounds import Chemical, ChemicalRegistry


class TestChemical:
    """Tests for the Chemical dataclass."""

    def test_basic_construction(self):
        c = Chemical("TestAcid", MW=60.0, atoms={"C": 2, "H": 4, "O": 2})
        assert c.ID == "TestAcid"
        assert c.MW == 60.0
        assert c.atoms == {"C": 2, "H": 4, "O": 2}

    def test_mw_must_be_positive(self):
        with pytest.raises(ValueError, match="MW must be > 0"):
            Chemical("Bad", MW=0.0)
        with pytest.raises(ValueError, match="MW must be > 0"):
            Chemical("Bad", MW=-1.0)

    def test_default_phase(self):
        c = Chemical("X", MW=1.0)
        assert c.phase == "l"

    def test_frozen(self):
        c = Chemical("X", MW=1.0)
        with pytest.raises(AttributeError):
            c.MW = 2.0


class TestChemicalRegistry:
    """Tests for the ChemicalRegistry."""

    def test_default_registry_not_empty(self, registry):
        assert len(registry) > 0

    def test_lookup_by_id(self, registry):
        ac = registry["AceticAcid"]
        assert ac.MW == pytest.approx(60.052)
        assert ac.atoms["C"] == 2

    def test_missing_key_raises(self, registry):
        with pytest.raises(KeyError, match="not found"):
            registry["NonExistentCompound"]

    def test_contains(self, registry):
        assert "AceticAcid" in registry
        assert "FakeChemical" not in registry

    def test_ids_property(self, registry):
        ids = registry.IDs
        assert isinstance(ids, tuple)
        assert "AceticAcid" in ids
        assert "O2" in ids
        assert "Yeast" in ids

    def test_register_custom_compound(self, registry):
        custom = Chemical("MyBug", MW=25.0, atoms={"C": 1, "H": 1.8, "O": 0.5})
        registry.register(custom)
        assert "MyBug" in registry
        assert registry["MyBug"].MW == 25.0

    def test_register_overwrites(self, registry):
        old_mw = registry["Yeast"].MW
        registry.register(Chemical("Yeast", MW=99.0, atoms={"C": 1}))
        assert registry["Yeast"].MW == 99.0
        # Restore for other tests
        registry.register(Chemical("Yeast", MW=old_mw,
                                   atoms={"C": 1, "H": 1.61, "O": 0.56, "N": 0.16}))

    def test_known_mw_water(self, registry):
        assert registry["Water"].MW == pytest.approx(18.015, rel=1e-3)

    def test_known_mw_co2(self, registry):
        assert registry["CO2"].MW == pytest.approx(44.009, rel=1e-3)

    def test_known_mw_o2(self, registry):
        assert registry["O2"].MW == pytest.approx(31.998, rel=1e-3)

    def test_known_mw_n2(self, registry):
        assert registry["N2"].MW == pytest.approx(28.014, rel=1e-3)

    def test_known_mw_nh3(self, registry):
        assert registry["NH3"].MW == pytest.approx(17.031, rel=1e-3)

    def test_known_mw_h3po4(self, registry):
        assert registry["H3PO4"].MW == pytest.approx(97.994, rel=1e-3)

    def test_known_mw_koh(self, registry):
        assert registry["KOH"].MW == pytest.approx(56.106, rel=1e-3)

    def test_yeast_has_nitrogen(self, registry):
        y = registry["Yeast"]
        assert y.atoms.get("N", 0) > 0, "Default Yeast should have N in formula"

    def test_yeast_cho_has_no_nitrogen(self, registry):
        y = registry["Yeast_CHO"]
        assert y.atoms.get("N", 0) == 0, "Yeast_CHO should have no N"
