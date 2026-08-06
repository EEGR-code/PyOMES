# -*- coding: utf-8 -*-
"""Tests for the string stoichiometry parser."""
import pytest

from PyOMES.chemistry import Species
from PyOMES.reactions.stoichiometry import StoichiometryEntry, _parse_stoichiometry
from PyOMES.reactions import EquilibriumReaction, KineticReaction


# ── Shared local-species fixtures ─────────────────────────────────────────────

ACETIC_ACID   = Species(id="AceticAcid", atoms={"C": 2, "H": 4, "O": 2}, charge=0,  MW=60.052)
ACETATE_MINUS = Species(id="Acetate-",   atoms={"C": 2, "H": 3, "O": 2}, charge=-1, MW=59.044)

_LOCAL = {"AceticAcid": ACETIC_ACID, "Acetate-": ACETATE_MINUS}


class TestStringStoichiometry:
    """Tests for _parse_stoichiometry and its integration into reaction classes."""

    # ── 1. Standard equilibrium string — only common_species IDs ─────────────

    def test_water_dissociation_no_species_kwarg(self):
        """H2O <-> H+ + OH- resolves from common_species with no species kwarg."""
        entries = _parse_stoichiometry("H2O,aq <-> H+,aq + OH-,aq", None)
        assert len(entries) == 3
        by_id = {e.species.id: e for e in entries}
        assert by_id["H2O"].coefficient == pytest.approx(-1.0)
        assert by_id["H+"].coefficient  == pytest.approx(+1.0)
        assert by_id["OH-"].coefficient == pytest.approx(+1.0)

    def test_common_species_through_equilibrium_reaction(self):
        """EquilibriumReaction accepts a string using only common_species IDs."""
        rxn = EquilibriumReaction(
            "H2O,aq <-> H+,aq + OH-,aq",
            log_K=-14.0,
            balance_elements=("H", "O"),
        )
        assert rxn.log_K == pytest.approx(-14.0)
        ids = {e.species.id for e in rxn.stoichiometry}
        assert ids == {"H2O", "H+", "OH-"}

    # ── 2. String with locally declared species via species kwarg ─────────────

    def test_local_species_via_kwarg(self):
        """AceticAcid dissociation resolves when caller supplies species dict."""
        entries = _parse_stoichiometry(
            "AceticAcid,aq <-> Acetate-,aq + H+,aq",
            _LOCAL,
        )
        by_id = {e.species.id: e for e in entries}
        assert by_id["AceticAcid"].coefficient == pytest.approx(-1.0)
        assert by_id["Acetate-"].coefficient   == pytest.approx(+1.0)
        assert by_id["H+"].coefficient         == pytest.approx(+1.0)

    def test_local_species_through_equilibrium_reaction(self):
        """EquilibriumReaction string with local species kwarg constructs correctly."""
        rxn = EquilibriumReaction(
            "AceticAcid,aq <-> Acetate-,aq + H+,aq",
            species=_LOCAL,
            log_K=-4.756,
            balance_elements=("C", "H", "O"),
        )
        assert rxn.log_K == pytest.approx(-4.756)
        assert set(rxn.species_ids) == {"AceticAcid", "Acetate-", "H+"}

    # ── 3. <-> signs: reactants negative, products positive ───────────────────

    def test_equilibrium_arrow_sign_convention(self):
        """Reactants get negative coefficients; products get positive."""
        entries = _parse_stoichiometry(
            "CO2,aq + H2O,aq <-> HCO3-,aq + H+,aq", None
        )
        by_id = {e.species.id: e for e in entries}
        assert by_id["CO2"].coefficient   == pytest.approx(-1.0)
        assert by_id["H2O"].coefficient   == pytest.approx(-1.0)
        assert by_id["HCO3-"].coefficient == pytest.approx(+1.0)
        assert by_id["H+"].coefficient    == pytest.approx(+1.0)

    # ── 4. -> parsed correctly for kinetic reaction ───────────────────────────

    def test_kinetic_arrow_parsed(self):
        """'->' arrow produces correct sign convention for KineticReaction."""
        rxn = KineticReaction(
            "CO2,g -> CO2,aq",
            rate_fn=lambda env: 0.0,
            balance_elements=("C", "O"),
        )
        by_phase = {e.phase: e for e in rxn.stoichiometry}
        assert by_phase["gas"].coefficient    == pytest.approx(-1.0)
        assert by_phase["liquid"].coefficient == pytest.approx(+1.0)

    def test_kinetic_arrow_sign_convention(self):
        """Direct parser: '->' gives reactant negative, product positive."""
        entries = _parse_stoichiometry("CO2,g -> CO2,aq", None, reaction_type="kinetic")
        by_phase = {e.phase: e for e in entries}
        assert by_phase["gas"].coefficient    == pytest.approx(-1.0)
        assert by_phase["liquid"].coefficient == pytest.approx(+1.0)

    # ── 5. All phase aliases resolve correctly ────────────────────────────────

    def test_phase_aq_resolves_to_liquid(self):
        entries = _parse_stoichiometry("H2O,aq <-> H+,aq + OH-,aq", None)
        for e in entries:
            assert e.phase == "liquid"

    def test_phase_l_resolves_to_liquid(self):
        entries = _parse_stoichiometry("H2O,l <-> H+,l + OH-,l", None)
        for e in entries:
            assert e.phase == "liquid"

    def test_phase_g_resolves_to_gas(self):
        entries = _parse_stoichiometry("CO2,g -> CO2,aq", None)
        gas_entries = [e for e in entries if e.phase == "gas"]
        assert len(gas_entries) == 1
        assert gas_entries[0].species.id == "CO2"

    def test_phase_s_resolves_to_solid(self):
        # Use NH3 which is in common_species; solid phase is unusual but valid
        entries = _parse_stoichiometry("NH3,s -> NH3,aq", None)
        phases = {e.phase for e in entries}
        assert "solid" in phases
        assert "liquid" in phases

    # ── 6. Implicit coefficient 1.0 ──────────────────────────────────────────

    def test_implicit_coefficient_is_one(self):
        entries = _parse_stoichiometry("H2O,aq <-> H+,aq + OH-,aq", None)
        for e in entries:
            assert abs(e.coefficient) == pytest.approx(1.0)

    # ── 7. Float coefficient ──────────────────────────────────────────────────

    def test_float_coefficient(self):
        entries = _parse_stoichiometry(
            "2.0 H2O,aq <-> 2.0 H+,aq + 2.0 OH-,aq", None
        )
        by_id = {e.species.id: e for e in entries}
        assert by_id["H2O"].coefficient == pytest.approx(-2.0)
        assert by_id["H+"].coefficient  == pytest.approx(+2.0)
        assert by_id["OH-"].coefficient == pytest.approx(+2.0)

    def test_fractional_float_coefficient(self):
        entries = _parse_stoichiometry(
            "1.5 H2O,aq <-> 1.5 H+,aq + 1.5 OH-,aq", None
        )
        by_id = {e.species.id: e for e in entries}
        assert by_id["H2O"].coefficient == pytest.approx(-1.5)
        assert by_id["H+"].coefficient  == pytest.approx(+1.5)

    # ── 8. Arrow direction mismatch raises ValueError ─────────────────────────

    def test_equilibrium_arrow_in_kinetic_raises(self):
        with pytest.raises(ValueError, match="<->"):
            _parse_stoichiometry(
                "H2O,aq <-> H+,aq + OH-,aq", None, reaction_type="kinetic"
            )

    def test_kinetic_arrow_in_equilibrium_raises(self):
        with pytest.raises(ValueError, match="->"):
            _parse_stoichiometry(
                "CO2,g -> CO2,aq", None, reaction_type="equilibrium"
            )

    def test_equilibrium_reaction_rejects_kinetic_arrow(self):
        with pytest.raises(ValueError, match="EquilibriumReaction"):
            EquilibriumReaction(
                "CO2,g -> CO2,aq",
                log_K=-1.5,
                balance_elements=("C", "O"),
            )

    def test_kinetic_reaction_rejects_equilibrium_arrow(self):
        with pytest.raises(ValueError, match="KineticReaction"):
            KineticReaction(
                "H2O,aq <-> H+,aq + OH-,aq",
                rate_fn=lambda env: 0.0,
                balance_elements=("H", "O"),
            )

    # ── 9. Unknown species ID raises ValueError listing available names ────────

    def test_unknown_species_raises(self):
        with pytest.raises(ValueError, match="Bogus"):
            _parse_stoichiometry("Bogus,aq <-> H+,aq + OH-,aq", None)

    def test_unknown_species_message_lists_common(self):
        with pytest.raises(ValueError, match="H2O"):
            # H2O is in common — message should list it
            _parse_stoichiometry("NotARealSpecies,aq <-> H+,aq", None)

    def test_unknown_species_suggests_species_kwarg(self):
        with pytest.raises(ValueError, match="species"):
            # no species kwarg provided → message should hint at it
            _parse_stoichiometry("AceticAcid,aq <-> Acetate-,aq + H+,aq", None)

    # ── 10. Missing phase suffix raises ValueError ────────────────────────────

    def test_missing_phase_suffix_raises(self):
        with pytest.raises(ValueError, match="phase suffix"):
            _parse_stoichiometry("H2O <-> H+ + OH-", None)

    def test_missing_phase_suffix_names_offending_term(self):
        with pytest.raises(ValueError, match="H2O"):
            _parse_stoichiometry("H2O <-> H+,aq + OH-,aq", None)

    # ── 11. Unknown phase suffix raises ValueError ────────────────────────────

    def test_unknown_phase_suffix_raises(self):
        with pytest.raises(ValueError, match="plasma"):
            _parse_stoichiometry("H2O,plasma <-> H+,aq + OH-,aq", None)

    def test_unknown_phase_lists_valid_options(self):
        with pytest.raises(ValueError, match="aq"):
            _parse_stoichiometry("H2O,xyz <-> H+,aq + OH-,aq", None)

    # ── 12. No-species-kwarg still works for all-common-species equation ──────

    def test_carbonate_equilibrium_no_species_kwarg(self):
        """CO2 + H2O <-> HCO3- + H+ needs no species kwarg."""
        rxn = EquilibriumReaction(
            "CO2,aq + H2O,aq <-> HCO3-,aq + H+,aq",
            log_K=-6.35,
            balance_elements=("C", "H", "O"),
        )
        assert set(rxn.species_ids) == {"CO2", "H2O", "HCO3-", "H+"}

    # ── 13. No-arrow string raises ────────────────────────────────────────────

    def test_no_arrow_raises(self):
        with pytest.raises(ValueError, match="No arrow"):
            _parse_stoichiometry("H2O,aq + H+,aq + OH-,aq", None)

    # ── 14. Both arrows raises ────────────────────────────────────────────────

    def test_both_arrows_raises(self):
        with pytest.raises(ValueError, match="both"):
            _parse_stoichiometry("A,aq <-> B,aq -> C,aq", None)

    # ── 15. Multiple terms each side ─────────────────────────────────────────

    def test_multiple_terms_each_side(self):
        """CO2 + H2O <-> HCO3- + H+ — two reactants, two products."""
        entries = _parse_stoichiometry(
            "CO2,aq + H2O,aq <-> HCO3-,aq + H+,aq", None
        )
        reactants = [e for e in entries if e.coefficient < 0]
        products  = [e for e in entries if e.coefficient > 0]
        assert len(reactants) == 2
        assert len(products)  == 2
