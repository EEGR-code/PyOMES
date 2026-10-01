# -*- coding: utf-8 -*-
"""ADM1 and BSM2 name only species they declare.

ADM1's species lookup raises on an id that is neither a declared species nor
one of its organisms. BSM2 declares its strong-ion lumps ``S_cat`` / ``S_an``
and passes its species table to the CV, so the conservation monitor counts
their charge.
"""
import pytest


class TestADM1Species:

    def test_declared_species_resolve(self):
        from vlmodels.adm1.base import _get_species
        assert _get_species("Glucose").atoms["C"] == 6

    def test_organisms_resolve_to_biomass(self):
        from vlmodels.adm1.base import _get_species, ORG, CHON_ORGS, BIO_CHON, BIO_CHO
        for org in ORG.values():
            expected = BIO_CHON if org in CHON_ORGS else BIO_CHO
            assert dict(_get_species(org).atoms) == expected

    def test_unknown_id_raises(self):
        from vlmodels.adm1.base import _get_species
        with pytest.raises(KeyError, match="'Glucoes' is not an ADM1 species"):
            _get_species("Glucoes")


class TestBSM2Species:

    @pytest.fixture(scope="class")
    def cv(self):
        from vlmodels.adm1.bsm2 import (
            build_bsm2_reactions, build_bsm2_cv, seed_bsm2_strong_ions,
        )
        cv = build_bsm2_cv(build_bsm2_reactions())
        seed_bsm2_strong_ions(cv)
        return cv

    def test_strong_ion_lumps_are_declared(self):
        from vlmodels.adm1.bsm2 import SPECIES
        assert (SPECIES["S_cat"].charge, dict(SPECIES["S_cat"].atoms)) == (+1, {})
        assert (SPECIES["S_an"].charge, dict(SPECIES["S_an"].atoms)) == (-1, {})

    def test_cv_holds_the_model_species(self, cv):
        from vlmodels.adm1.bsm2 import SPECIES
        for sp_id, sp in SPECIES.items():
            assert cv.species[sp_id] is sp

    def test_monitor_counts_the_lumps(self, cv):
        registry = cv._conservation_monitor._species_registry
        assert registry["S_cat"].charge == +1 and registry["S_an"].charge == -1
