# -*- coding: utf-8 -*-
"""What each engine does with a negative component total.

A total within roundoff of zero (at most ``ROUNDOFF_RTOL`` of the sum of the
magnitudes it was summed from) reads as zero. A larger negative total, or any
negative amount of a strong ion, raises ``NegativeTotalError`` naming the
component; nothing is changed to make it fit. Run on the user-defined model of
``test_user_defined_model.py`` under the Bisection and NR engines.
"""
import importlib.util
import math
from pathlib import Path

import pytest

from PyOMES.chemical_equilibrium import NegativeTotalError
from PyOMES.chemical_equilibrium.totals import ROUNDOFF_RTOL, checked_total
from PyOMES.chemistry import Species
from PyOMES.core import ControlVolume, Simulation
from PyOMES.reactions import ReactionSystem


def _load_user_defined_model():
    path = Path(__file__).with_name("test_user_defined_model.py")
    spec = importlib.util.spec_from_file_location("_negative_totals_user_defined_model", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_UDM = _load_user_defined_model()
NA_PLUS = Species(id="Na+", atoms={"Na": 1}, charge=+1)
ENGINES = ("bisection", "nr")


def _cv(engine):
    base = _UDM._build_cv()
    reaction_system = base.reaction_system
    if engine == "nr":
        reaction_system = ReactionSystem(
            list(reaction_system.reactions), label="user_defined", solver="newton_raphson",
        )
    return ControlVolume(phases=dict(base.phases), reaction_system=reaction_system,
                         label="negative_totals", species=[NA_PLUS])


def _solve(cv):
    liquid = cv.phases["liquid"]
    result = cv.reaction_system.engine.solve(phases=cv.phases, T_K=float(liquid.T_K))
    result.apply_to_phases(cv.phases)


class TestCheckedTotal:
    def test_a_non_negative_total_is_returned_as_summed(self):
        assert checked_total("X", [1.0, 2.0], engine="E") == 3.0
        assert checked_total("X", [], engine="E") == 0.0

    def test_a_negative_total_within_roundoff_reads_as_zero(self):
        assert checked_total("X", [1.0, -1.0 - 0.5 * ROUNDOFF_RTOL], engine="E") == 0.0

    def test_a_negative_total_beyond_roundoff_raises(self):
        with pytest.raises(NegativeTotalError, match="'X'"):
            checked_total("X", [1.0, -1.0 - 10.0 * ROUNDOFF_RTOL], engine="E")

    def test_a_single_negative_amount_raises_however_small(self):
        with pytest.raises(NegativeTotalError):
            checked_total("Na+", [-1e-30], engine="E")

    def test_it_is_not_a_value_error(self):
        """Code that recovers from a ValueError in a solve must not swallow it."""
        assert not issubclass(NegativeTotalError, ValueError)


@pytest.mark.parametrize("engine", ENGINES)
class TestEngines:
    def test_a_negative_acid_total_raises_naming_it(self, engine):
        cv = _cv(engine)
        cv.phases["liquid"].n_mol["AceticAcid"] = -1e-6
        with pytest.raises(NegativeTotalError, match="AceticAcid"):
            _solve(cv)

    def test_a_total_negative_within_roundoff_solves_as_zero(self, engine):
        cv = _cv(engine)
        liquid = cv.phases["liquid"].n_mol
        liquid["AceticAcid"] = 1e-3
        liquid["Acetate-"] = -1e-3 * (1.0 + 0.1 * ROUNDOFF_RTOL)
        _solve(cv)
        total = liquid.get("AceticAcid", 0.0) + liquid.get("Acetate-", 0.0)
        assert abs(total) <= 1e-40
        assert math.isfinite(float(cv.phases["liquid"].pH))

    def test_a_negative_strong_ion_raises(self, engine):
        cv = _cv(engine)
        cv.phases["liquid"].n_mol["Na+"] = -1e-9
        with pytest.raises(NegativeTotalError, match="Na\\+"):
            _solve(cv)

    def test_advance_names_the_step(self, engine):
        cv = _cv(engine)
        cv.phases["liquid"].n_mol["AceticAcid"] = -1e-6
        with pytest.raises(NegativeTotalError) as info:
            cv.advance(dt_h=0.05, t_h=1.5)
        notes = " ".join(getattr(info.value, "__notes__", []))
        assert "'negative_totals'" in notes and "0.05 h" in notes and "t_h = 1.5" in notes

    def test_the_first_solve_of_a_run_raises_rather_than_warning(self, engine):
        cv = _cv(engine)
        cv.phases["liquid"].n_mol["AceticAcid"] = -1e-6
        with pytest.raises(NegativeTotalError) as info:
            Simulation(cvs={"main": cv}, label="t").run(tau_h=0.1, n_steps=1)
        assert "first solve" in " ".join(getattr(info.value, "__notes__", []))
