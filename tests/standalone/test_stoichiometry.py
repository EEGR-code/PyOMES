# -*- coding: utf-8 -*-
"""Tests for the string stoichiometry parser."""
import pytest

from PyOMES.chemistry import Species
from PyOMES.reactions.stoichiometry import StoichiometryEntry, _parse_stoichiometry
from PyOMES.reactions import EquilibriumReaction, KineticReaction


# ── Species defined by the test, as a user would ─────────────────────────────

H2O        = Species(id="H2O",   atoms={"H": 2, "O": 1})
H_PLUS     = Species(id="H+",    atoms={"H": 1},                 charge=+1)
OH_MINUS   = Species(id="OH-",   atoms={"O": 1, "H": 1},         charge=-1)
CO2        = Species(id="CO2",   atoms={"C": 1, "O": 2})
HCO3_MINUS = Species(id="HCO3-", atoms={"H": 1, "C": 1, "O": 3}, charge=-1)
NH3        = Species(id="NH3",   atoms={"N": 1, "H": 3})

ACETIC_ACID   = Species(id="AceticAcid", atoms={"C": 2, "H": 4, "O": 2}, charge=0,  MW=60.052)
ACETATE_MINUS = Species(id="Acetate-",   atoms={"C": 2, "H": 3, "O": 2}, charge=-1, MW=59.044)


def _sp(*species):
    """``{id: Species}`` for exactly the species given."""
    return {s.id: s for s in species}


_WATER     = _sp(H2O, H_PLUS, OH_MINUS)
_CARBONATE = _sp(CO2, H2O, HCO3_MINUS, H_PLUS)
_LOCAL     = _sp(ACETIC_ACID, ACETATE_MINUS, H_PLUS)


class TestStringStoichiometry:
    """Tests for _parse_stoichiometry and its integration into reaction classes."""

    # ── 1. Standard equilibrium string — user-defined species ───────────────

    def test_water_dissociation(self):
        """H2O <-> H+ + OH- resolves from the species passed."""
        entries = _parse_stoichiometry("H2O,aq <-> H+,aq + OH-,aq", _WATER)
        assert len(entries) == 3
        by_id = {e.species.id: e for e in entries}
        assert by_id["H2O"].coefficient == pytest.approx(-1.0)
        assert by_id["H+"].coefficient  == pytest.approx(+1.0)
        assert by_id["OH-"].coefficient == pytest.approx(+1.0)
        assert by_id["H2O"].species is H2O

    def test_string_through_equilibrium_reaction(self):
        """EquilibriumReaction accepts a string with the species passed."""
        rxn = EquilibriumReaction(
            "H2O,aq <-> H+,aq + OH-,aq",
            species=_WATER,
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
            "CO2,aq + H2O,aq <-> HCO3-,aq + H+,aq", _CARBONATE
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
            species=_sp(CO2),
            rate_fn=lambda env: 0.0,
            balance_elements=("C", "O"),
        )
        by_phase = {e.phase: e for e in rxn.stoichiometry}
        assert by_phase["gas"].coefficient    == pytest.approx(-1.0)
        assert by_phase["liquid"].coefficient == pytest.approx(+1.0)

    def test_kinetic_arrow_sign_convention(self):
        """Direct parser: '->' gives reactant negative, product positive."""
        entries = _parse_stoichiometry("CO2,g -> CO2,aq", _sp(CO2), reaction_type="kinetic")
        by_phase = {e.phase: e for e in entries}
        assert by_phase["gas"].coefficient    == pytest.approx(-1.0)
        assert by_phase["liquid"].coefficient == pytest.approx(+1.0)

    # ── 5. All phase aliases resolve correctly ────────────────────────────────

    def test_phase_aq_resolves_to_liquid(self):
        entries = _parse_stoichiometry("H2O,aq <-> H+,aq + OH-,aq", _WATER)
        for e in entries:
            assert e.phase == "liquid"

    def test_phase_l_resolves_to_liquid(self):
        entries = _parse_stoichiometry("H2O,l <-> H+,l + OH-,l", _WATER)
        for e in entries:
            assert e.phase == "liquid"

    def test_phase_g_resolves_to_gas(self):
        entries = _parse_stoichiometry("CO2,g -> CO2,aq", _sp(CO2))
        gas_entries = [e for e in entries if e.phase == "gas"]
        assert len(gas_entries) == 1
        assert gas_entries[0].species.id == "CO2"

    def test_phase_s_resolves_to_solid(self):
        # Solid NH3 is unusual but valid for the parser
        entries = _parse_stoichiometry("NH3,s -> NH3,aq", _sp(NH3))
        phases = {e.phase for e in entries}
        assert "solid" in phases
        assert "liquid" in phases

    # ── 6. Implicit coefficient 1.0 ──────────────────────────────────────────

    def test_implicit_coefficient_is_one(self):
        entries = _parse_stoichiometry("H2O,aq <-> H+,aq + OH-,aq", _WATER)
        for e in entries:
            assert abs(e.coefficient) == pytest.approx(1.0)

    # ── 7. Float coefficient ──────────────────────────────────────────────────

    def test_float_coefficient(self):
        entries = _parse_stoichiometry(
            "2.0 H2O,aq <-> 2.0 H+,aq + 2.0 OH-,aq", _WATER
        )
        by_id = {e.species.id: e for e in entries}
        assert by_id["H2O"].coefficient == pytest.approx(-2.0)
        assert by_id["H+"].coefficient  == pytest.approx(+2.0)
        assert by_id["OH-"].coefficient == pytest.approx(+2.0)

    def test_fractional_float_coefficient(self):
        entries = _parse_stoichiometry(
            "1.5 H2O,aq <-> 1.5 H+,aq + 1.5 OH-,aq", _WATER
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
            _parse_stoichiometry("Bogus,aq <-> H+,aq + OH-,aq", _WATER)

    def test_unknown_species_message_lists_species_passed(self):
        with pytest.raises(ValueError, match=r"available: .*H2O"):
            _parse_stoichiometry("NotARealSpecies,aq <-> H+,aq", _WATER)

    def test_unknown_species_suggests_species_kwarg(self):
        with pytest.raises(ValueError, match="no species were passed.*species="):
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

    # ── 12. Ids are looked up only in the species passed ──────────────────────

    def test_carbonate_equilibrium(self):
        """CO2 + H2O <-> HCO3- + H+ resolves from the species passed."""
        rxn = EquilibriumReaction(
            "CO2,aq + H2O,aq <-> HCO3-,aq + H+,aq",
            species=_CARBONATE,
            log_K=-6.35,
            balance_elements=("C", "H", "O"),
        )
        assert set(rxn.species_ids) == {"CO2", "H2O", "HCO3-", "H+"}

    def test_no_species_kwarg_raises(self):
        """Without species=, even water's ids are unknown."""
        with pytest.raises(ValueError, match="'H2O' is not among the species passed"):
            EquilibriumReaction(
                "H2O,aq <-> H+,aq + OH-,aq",
                log_K=-14.0,
                balance_elements=("H", "O"),
            )

    def test_id_missing_from_species_passed_raises(self):
        """An id the caller did not pass raises, whatever else defines it."""
        with pytest.raises(ValueError, match="'OH-' is not among the species passed"):
            _parse_stoichiometry("H2O,aq <-> H+,aq + OH-,aq", _LOCAL | _sp(H2O))

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
            "CO2,aq + H2O,aq <-> HCO3-,aq + H+,aq", _CARBONATE
        )
        reactants = [e for e in entries if e.coefficient < 0]
        products  = [e for e in entries if e.coefficient > 0]
        assert len(reactants) == 2
        assert len(products)  == 2
