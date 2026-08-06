# -*- coding: utf-8 -*-
"""Tests for the ChemicalEquilibriumEngine -> BisectionChemicalEquilibriumEngine rename.

The unqualified name ``ChemicalEquilibriumEngine`` read as "the" default/
canonical engine, when it is actually one of three peer implementations of
``ChemicalEquilibriumEngineProtocol`` (the others being
``NRChemicalEquilibriumEngine`` and ``PHREEQCChemicalEquilibriumEngine``) —
specifically the original, simplest one (1-D bisection over the charge
balance). Renamed to ``BisectionChemicalEquilibriumEngine``, with
``ChemicalEquilibriumEngine`` kept as a deprecated alias for one phase.

These tests exercise the alias directly — unlike the rest of the suite,
which was migrated to the new name.
"""
from __future__ import annotations

import warnings

import pytest


def _water_reaction():
    from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
    from PyOMES.chemistry.common_species import H2O, H_plus, OH_minus
    return EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=H2O, phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=H_plus, phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=OH_minus, phase="liquid", coefficient=+1.0),
        ],
        log_K=-14.0, label="water",
    )


class TestDeprecatedAlias:
    def test_direct_construction_emits_deprecation_warning(self):
        from PyOMES.chemical_equilibrium import ChemicalEquilibriumEngine
        with pytest.warns(DeprecationWarning):
            ChemicalEquilibriumEngine()

    def test_direct_construction_is_a_bisection_engine(self):
        from PyOMES.chemical_equilibrium import (
            BisectionChemicalEquilibriumEngine, ChemicalEquilibriumEngine,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            eng = ChemicalEquilibriumEngine()
        assert isinstance(eng, BisectionChemicalEquilibriumEngine)

    def test_from_reactions_emits_deprecation_warning(self):
        from PyOMES.chemical_equilibrium import ChemicalEquilibriumEngine
        with pytest.warns(DeprecationWarning):
            ChemicalEquilibriumEngine.from_reactions([_water_reaction()], T_K=298.15)

    def test_from_reactions_returns_working_engine(self):
        """from_reactions is a cls-based classmethod — the alias subclass
        must not break its construction path."""
        from PyOMES.chemical_equilibrium import ChemicalEquilibriumEngine
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            eng = ChemicalEquilibriumEngine.from_reactions(
                [_water_reaction()], T_K=298.15,
            )
        assert eng._equilibrium_set is not None

    def test_new_name_construction_emits_no_warning(self):
        from PyOMES.chemical_equilibrium import BisectionChemicalEquilibriumEngine
        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            BisectionChemicalEquilibriumEngine()  # must not raise

    def test_new_name_from_reactions_emits_no_warning(self):
        from PyOMES.chemical_equilibrium import BisectionChemicalEquilibriumEngine
        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            BisectionChemicalEquilibriumEngine.from_reactions(
                [_water_reaction()], T_K=298.15,
            )  # must not raise

    def test_alias_satisfies_engine_protocol(self):
        from PyOMES.chemical_equilibrium.protocols import ChemicalEquilibriumEngineProtocol
        from PyOMES.chemical_equilibrium import ChemicalEquilibriumEngine
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            eng = ChemicalEquilibriumEngine()
        assert isinstance(eng, ChemicalEquilibriumEngineProtocol)

    def test_solve_behavior_identical_via_alias_and_direct(self):
        """Same reactions/inputs -> identical EquilibriumResult regardless
        of which name constructed the engine (pure rename, no behavior
        change)."""
        from PyOMES.chemical_equilibrium import (
            BisectionChemicalEquilibriumEngine, ChemicalEquilibriumEngine,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            eng_alias = ChemicalEquilibriumEngine.from_reactions(
                [_water_reaction()], T_K=298.15,
            )
        eng_direct = BisectionChemicalEquilibriumEngine.from_reactions(
            [_water_reaction()], T_K=298.15,
        )
        out_alias = eng_alias.solve(CT_TIC=0.0)
        out_direct = eng_direct.solve(CT_TIC=0.0)
        assert out_alias.pH == pytest.approx(out_direct.pH)
