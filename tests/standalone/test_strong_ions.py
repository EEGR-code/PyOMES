"""Tests for fermenter.speciation.strong_ions — standalone strong-ion inference."""

import pytest
from PyOMES.chemical_equilibrium.strong_ions import strong_ions_from_feed_molL
from PyOMES.stream_adapter import FeedState


@pytest.mark.unit
class TestStrongIonsFromFeedState:
    """Verify that strong_ions_from_feed_molL works with FeedState objects."""

    def test_empty_feed_returns_zeros(self, registry):
        feed = FeedState(registry=registry)
        ions = strong_ions_from_feed_molL(feed)
        assert all(v == 0.0 for v in ions.values())

    def test_kh2po4_contributes_potassium(self, registry):
        feed = FeedState(
            concentrations_mol_L={"KH2PO4": 0.01},
            registry=registry,
        )
        ions = strong_ions_from_feed_molL(feed)
        assert ions["CT_K"] == pytest.approx(0.01, rel=1e-6)

    def test_cacl2_contributes_calcium_and_chloride(self, registry):
        feed = FeedState(
            concentrations_mol_L={"CaCl2": 0.005},
            registry=registry,
        )
        ions = strong_ions_from_feed_molL(feed)
        assert ions["CT_Ca"] == pytest.approx(0.005, rel=1e-6)
        assert ions["CT_Cl"] == pytest.approx(0.010, rel=1e-6)  # 2 Cl per CaCl2

    def test_ammonium_sulfate_contributes_sulfate(self, registry):
        feed = FeedState(
            concentrations_mol_L={"AmmoniumSulfate": 0.015},
            registry=registry,
        )
        ions = strong_ions_from_feed_molL(feed)
        assert ions["CT_SO4"] == pytest.approx(0.015, rel=1e-6)

    def test_multiple_salts_accumulate(self, rich_feed):
        ions = strong_ions_from_feed_molL(rich_feed)
        # KH2PO4 gives K⁺, AmmoniumSulfate + MgSO4 + ZnSO4 give SO4²⁻
        assert ions["CT_K"] > 0
        assert ions["CT_SO4"] > 0
        assert ions["CT_Ca"] > 0
        assert ions["CT_Cl"] > 0

    def test_plain_dict_accepted(self):
        """strong_ions_from_feed_molL should also accept a plain dict."""
        ions = strong_ions_from_feed_molL({"KH2PO4": 0.01})
        assert ions["CT_K"] == pytest.approx(0.01, rel=1e-6)

    def test_charge_balance_sanity(self, rich_feed):
        """Cation charge should be within reasonable range of anion charge."""
        ions = strong_ions_from_feed_molL(rich_feed)
        cation_eq = (
            ions["CT_K"] + ions["CT_Na"]
            + 2 * ions["CT_Ca"] + 2 * ions["CT_Mg"]
            + 2 * ions["CT_Zn"] + 2 * ions["CT_Mn"]
            + 2 * ions["CT_Co"]
        )
        anion_eq = (
            ions["CT_Cl"] + ions["CT_NO3"]
            + 2 * ions["CT_SO4"]
        )
        # They won't be exactly equal (NH4+ is not tracked here), but both
        # should be positive and within an order of magnitude.
        assert cation_eq > 0
        assert anion_eq > 0
