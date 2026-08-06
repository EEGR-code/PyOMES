# -*- coding: utf-8 -*-
"""Tests for ``check_species_consistency``."""

import warnings

import pytest


def _rxn_with(*species):
    """Wrap Species objects in a stub reaction with a ``stoichiometry``
    attribute, suitable for ``check_species_consistency``."""
    class _Entry:
        def __init__(self, sp):
            self.species = sp

    class _Rxn:
        def __init__(self, ss):
            self.stoichiometry = [_Entry(s) for s in ss]

    return _Rxn(list(species))


class TestCheckSpeciesConsistency:

    def test_clean_imports_no_warning(self):
        from PyOMES.chemistry import check_species_consistency
        from PyOMES.chemistry.common_species import H2O, CO2, H_plus, OH_minus
        with warnings.catch_warnings():
            warnings.simplefilter("error")  # any warning would fail
            check_species_consistency([
                _rxn_with(H2O, H_plus, OH_minus),
                _rxn_with(CO2, H2O, H_plus),
            ])

    def test_hard_conflict_atoms_raises(self):
        from PyOMES.chemistry import Species, check_species_consistency, SpeciesConflictError
        a = Species(id="X", atoms={"C": 1})
        b = Species(id="X", atoms={"C": 2})
        with pytest.raises(SpeciesConflictError) as exc:
            check_species_consistency([_rxn_with(a), _rxn_with(b)])
        assert exc.value.species_id == "X"

    def test_hard_conflict_charge_raises(self):
        from PyOMES.chemistry import Species, check_species_consistency, SpeciesConflictError
        a = Species(id="X", atoms={"O": 1, "H": 1}, charge=0)
        b = Species(id="X", atoms={"O": 1, "H": 1}, charge=-1)
        with pytest.raises(SpeciesConflictError):
            check_species_consistency([_rxn_with(a), _rxn_with(b)])

    def test_hard_conflict_MW_raises(self):
        from PyOMES.chemistry import Species, check_species_consistency, SpeciesConflictError
        a = Species(id="X", atoms={"C": 1}, MW=12.0)
        b = Species(id="X", atoms={"C": 1}, MW=13.0)
        with pytest.raises(SpeciesConflictError):
            check_species_consistency([_rxn_with(a), _rxn_with(b)])

    def test_soft_conflict_warns(self):
        from PyOMES.chemistry import Species, check_species_consistency
        a = Species(id="X", atoms={"C": 1})
        b = Species(id="X", atoms={"C": 1})  # equal data, distinct object
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            check_species_consistency([_rxn_with(a), _rxn_with(b)])
            soft = [x for x in w if issubclass(x.category, UserWarning)]
            assert len(soft) == 1
            assert "X" in str(soft[0].message)

    def test_soft_conflict_raise_mode(self):
        from PyOMES.chemistry import Species, check_species_consistency, SpeciesConflictError
        a = Species(id="X", atoms={"C": 1})
        b = Species(id="X", atoms={"C": 1})
        with pytest.raises(SpeciesConflictError):
            check_species_consistency(
                [_rxn_with(a), _rxn_with(b)],
                soft_conflicts="raise",
            )

    def test_soft_conflict_ignore_mode(self):
        from PyOMES.chemistry import Species, check_species_consistency
        a = Species(id="X", atoms={"C": 1})
        b = Species(id="X", atoms={"C": 1})
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            check_species_consistency(
                [_rxn_with(a), _rxn_with(b)],
                soft_conflicts="ignore",
            )

    def test_invalid_soft_conflicts_value_raises(self):
        from PyOMES.chemistry import check_species_consistency
        with pytest.raises(ValueError):
            check_species_consistency([], soft_conflicts="bogus")

    def test_single_object_referenced_many_times_silent(self):
        """One Species shared by 50 entries is the happy path."""
        from PyOMES.chemistry import Species, check_species_consistency
        sp = Species(id="X", atoms={"C": 1})
        rxns = [_rxn_with(sp) for _ in range(50)]
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            check_species_consistency(rxns)
