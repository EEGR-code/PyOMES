# -*- coding: utf-8 -*-
"""Tests for the ``activity_model`` argument: names, model objects and bad input.

The same argument is taken by both chemical-equilibrium engines,
``ReactionSystem.configure_engine``, ``StirredTankBuilder.chemistry`` and
``ChemistryConfig``; it is a name (``"ideal"``, ``"davies"``, ``"sit"``) or an
activity model object.
"""
from __future__ import annotations

import json

import pytest

from PyOMES.chemical_equilibrium.engines.bisection.engine import BisectionChemicalEquilibriumEngine
from PyOMES.chemical_equilibrium.engines.nr.engine import NRChemicalEquilibriumEngine
from PyOMES.reactions import EquilibriumReaction, ReactionSystem, StoichiometryEntry
from PyOMES.templates.stirred_tank import ChemistryConfig
from PyOMES.thermo import DaviesLiquidModel, IdealLiquidModel, SITLiquidModel, make_activity_model


def _carbonate_reactions():
    """Water and the two carbonate dissociations, for building either engine."""
    from PyOMES.chemistry.common_species import (
        CO2, CO3_2minus, H2O, H_plus, HCO3_minus, OH_minus,
    )

    def e(species, coeff):
        return StoichiometryEntry(species=species, phase="liquid", coefficient=coeff)

    return [
        EquilibriumReaction([e(H2O, -1), e(H_plus, 1), e(OH_minus, 1)],
                            log_K=-14.0, balance_elements=("H", "O"), label="water"),
        EquilibriumReaction([e(CO2, -1), e(H2O, -1), e(HCO3_minus, 1), e(H_plus, 1)],
                            log_K=-6.35, total_id="CO2", balance_elements=("C", "H", "O"),
                            label="co2_first"),
        EquilibriumReaction([e(HCO3_minus, -1), e(CO3_2minus, 1), e(H_plus, 1)],
                            log_K=-10.33, total_id="CO2", balance_elements=("C", "H", "O"),
                            label="co2_second"),
    ]


def _bisection(activity_model):
    return BisectionChemicalEquilibriumEngine.from_reactions(
        _carbonate_reactions(), activity_model=activity_model)


def _nr(activity_model):
    return NRChemicalEquilibriumEngine.from_reactions(
        _carbonate_reactions(), activity_model=activity_model)


def _solve_bisection(engine):
    return engine.solve(CT_TIC=0.02, CT_Na=0.3, CT_Cl=0.25)


def _solve_nr(engine):
    return engine.solve(totals={"CO2": 0.02}, strong_ions={"CT_Na": 0.3, "CT_Cl": 0.25})


class _NamelessDavies:
    """A user-written model with ``gamma`` only and no ``name`` label.

    Forwards to the Davies equation, so an engine using it should reproduce
    Davies results exactly.
    """

    def __init__(self):
        self._davies = DaviesLiquidModel()

    def gamma(self, z, I_molL, *, T_K):
        return self._davies.gamma(z, I_molL, T_K=T_K)


class TestMakeActivityModel:
    @pytest.mark.parametrize("name, cls", [
        ("ideal", IdealLiquidModel), ("davies", DaviesLiquidModel), ("sit", SITLiquidModel),
        ("Davies", DaviesLiquidModel), (" SIT ", SITLiquidModel),
    ])
    def test_name_builds_model(self, name, cls):
        assert type(make_activity_model(name)) is cls

    def test_default_is_ideal(self):
        assert isinstance(make_activity_model(), IdealLiquidModel)

    def test_model_object_returned_unchanged(self):
        model = SITLiquidModel(epsilon={("Na+", "Cl-"): 0.05})
        assert make_activity_model(model) is model

    def test_unknown_name_raises(self):
        with pytest.raises(ValueError, match="pitzer"):
            make_activity_model("pitzer")

    def test_non_model_object_raises(self):
        with pytest.raises(TypeError, match="gamma"):
            make_activity_model(object())


@pytest.mark.parametrize("build", [_bisection, _nr], ids=["bisection", "nr"])
class TestEngineArgument:
    @pytest.mark.parametrize("name, cls", [
        ("ideal", IdealLiquidModel), ("davies", DaviesLiquidModel), ("sit", SITLiquidModel),
    ])
    def test_name_resolved_at_construction(self, build, name, cls):
        assert type(build(name).activity_model) is cls

    def test_default_is_ideal(self, build):
        engine = (BisectionChemicalEquilibriumEngine() if build is _bisection
                  else NRChemicalEquilibriumEngine.from_reactions(_carbonate_reactions()))
        assert isinstance(engine.activity_model, IdealLiquidModel)

    def test_model_object_stored_as_given(self, build):
        model = DaviesLiquidModel()
        assert build(model).activity_model is model

    def test_unknown_name_fails_at_construction(self, build):
        with pytest.raises(ValueError):
            build("pitzer")


class TestNamelessModelIsNotTreatedAsIdeal:
    """A model object with no ``name`` goes through the ionic-strength loop.

    The ideal fast path is chosen by type (``IdealLiquidModel``), so this model
    must give Davies results, not ideal ones.
    """

    def test_bisection(self):
        nameless = _solve_bisection(_bisection(_NamelessDavies()))
        davies = _solve_bisection(_bisection("davies"))
        ideal = _solve_bisection(_bisection("ideal"))
        assert nameless.pH == davies.pH
        assert nameless.pH != pytest.approx(ideal.pH, abs=1e-3)

    def test_nr(self):
        nameless = _solve_nr(_nr(_NamelessDavies()))
        davies = _solve_nr(_nr("davies"))
        ideal = _solve_nr(_nr("ideal"))
        assert nameless.pH == davies.pH
        assert nameless.pH != pytest.approx(ideal.pH, abs=1e-3)


class TestReactionSystemConfigureEngine:
    def test_name_reaches_engine(self):
        rs = ReactionSystem(_carbonate_reactions())
        rs.configure_engine(activity_model="sit")
        assert isinstance(rs.engine.activity_model, SITLiquidModel)

    def test_model_object_reaches_engine(self):
        model = SITLiquidModel(epsilon={("Na+", "Cl-"): 0.05})
        rs = ReactionSystem(_carbonate_reactions(), solver="newton_raphson")
        rs.configure_engine(activity_model=model)
        assert rs.engine.activity_model is model

    def test_unknown_name_fails_at_the_call(self):
        rs = ReactionSystem(_carbonate_reactions())
        with pytest.raises(ValueError):
            rs.configure_engine(activity_model="pitzer")


class TestChemistryConfig:
    def test_unknown_name_fails_at_construction(self):
        with pytest.raises(ValueError):
            ChemistryConfig(activity_model="pitzer")

    def test_named_model_round_trips_through_json(self):
        c = ChemistryConfig(activity_model="davies")
        c2 = ChemistryConfig.from_dict(json.loads(json.dumps(c.to_dict())))
        assert c2.activity_model == "davies"

    def test_model_object_round_trips_as_the_object(self):
        model = SITLiquidModel(epsilon={("Na+", "Cl-"): 0.05})
        c2 = ChemistryConfig.from_dict(ChemistryConfig(activity_model=model).to_dict())
        assert c2.activity_model is model
