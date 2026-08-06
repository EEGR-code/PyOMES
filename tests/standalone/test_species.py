# -*- coding: utf-8 -*-
"""Tests for the Species type and common_species declarations."""

import pytest


# ── Species data class ────────────────────────────────────────────────

class TestSpecies:

    def test_default_charge_and_MW(self):
        from PyOMES.chemistry import Species
        s = Species(id="CH4", atoms={"C": 1, "H": 4})
        assert s.charge == 0
        # MW auto-computed: 12.011 + 4*1.008 = 16.043
        assert s.MW == pytest.approx(16.043, rel=1e-4)

    def test_explicit_fields(self):
        from PyOMES.chemistry import Species
        s = Species(id="HCO3-", atoms={"H": 1, "C": 1, "O": 3},
                    charge=-1, MW=61.016)
        assert s.id == "HCO3-"
        assert s.charge == -1
        assert s.MW == pytest.approx(61.016)
        assert dict(s.atoms) == {"H": 1, "C": 1, "O": 3}

    def test_value_equality(self):
        from PyOMES.chemistry import Species
        a = Species(id="CO2", atoms={"C": 1, "O": 2}, charge=0, MW=44.01)
        b = Species(id="CO2", atoms={"C": 1, "O": 2}, charge=0, MW=44.01)
        assert a == b
        assert hash(a) == hash(b)

    def test_inequality_on_atoms(self):
        from PyOMES.chemistry import Species
        a = Species(id="X", atoms={"C": 1})
        b = Species(id="X", atoms={"C": 2})
        assert a != b

    def test_inequality_on_charge(self):
        from PyOMES.chemistry import Species
        a = Species(id="X", atoms={"C": 1}, charge=0)
        b = Species(id="X", atoms={"C": 1}, charge=-1)
        assert a != b

    def test_set_membership(self):
        from PyOMES.chemistry import Species
        a = Species(id="A", atoms={"C": 1})
        b = Species(id="B", atoms={"C": 1})
        s = {a, b}
        assert len(s) == 2

    def test_atoms_is_immutable(self):
        from PyOMES.chemistry import Species
        s = Species(id="X", atoms={"C": 1})
        with pytest.raises(TypeError):
            s.atoms["C"] = 2  # MappingProxyType blocks assignment

    def test_frozen_blocks_attribute_assignment(self):
        from PyOMES.chemistry import Species
        s = Species(id="X", atoms={"C": 1})
        with pytest.raises(Exception):
            s.id = "Y"

    def test_repr_includes_id(self):
        from PyOMES.chemistry import Species
        s = Species(id="CO2", atoms={"C": 1, "O": 2})
        assert "CO2" in repr(s)


# ── common_species ───────────────────────────────────────────────────

class TestCommonSpecies:

    def test_water_atoms(self):
        from PyOMES.chemistry.common_species import H2O
        assert dict(H2O.atoms) == {"H": 2, "O": 1}
        assert H2O.charge == 0

    def test_protons_and_hydroxide_charges(self):
        from PyOMES.chemistry.common_species import H_plus, OH_minus
        assert H_plus.charge == 1
        assert OH_minus.charge == -1

    def test_carbonate_double_negative(self):
        from PyOMES.chemistry.common_species import CO3_2minus
        assert CO3_2minus.charge == -2

    def test_ammonium_positive(self):
        from PyOMES.chemistry.common_species import NH4_plus
        assert NH4_plus.charge == 1
        assert dict(NH4_plus.atoms) == {"N": 1, "H": 4}

    def test_identity_preserved_across_imports(self):
        from PyOMES.chemistry.common_species import CO2 as A
        from PyOMES.chemistry.common_species import CO2 as B
        assert A is B


# ── validate_balance with charge ──────────────────────────────────────

class TestChargeBalance:

    def test_balanced_dissociation_passes(self):
        from PyOMES.reactions import StoichiometryEntry, validate_balance
        from PyOMES.chemistry.common_species import H_plus, OH_minus, H2O
        entries = [
            StoichiometryEntry(species=H2O,      phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=H_plus,   phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=OH_minus, phase="liquid", coefficient=+1.0),
        ]
        validate_balance(entries, elements=("H", "O"), check_charge=True)

    def test_imbalanced_charge_raises(self):
        from PyOMES.reactions import StoichiometryEntry, StoichiometryError, validate_balance
        from PyOMES.chemistry import Species
        from PyOMES.chemistry.common_species import H_plus, H2O
        # H2O -> H+ + neutral OH placeholder: elements balance, charge does not.
        entries = [
            StoichiometryEntry(species=H2O,    phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=H_plus, phase="liquid", coefficient=+1.0),
            StoichiometryEntry(
                species=Species(id="OH-neutral", atoms={"O": 1, "H": 1}, charge=0),
                phase="liquid", coefficient=+1.0),
        ]
        validate_balance(entries, elements=("H", "O"))  # element-only passes
        with pytest.raises(StoichiometryError) as exc:
            validate_balance(entries, elements=("H", "O"), check_charge=True)
        assert exc.value.element == "charge"
