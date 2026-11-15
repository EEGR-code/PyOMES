# -*- coding: utf-8 -*-
"""What each engine does with a negative component total.

A total negative within tolerance of zero (the larger of
``PyOMES.config.negative_total_tolerance_mol`` and ``ROUNDOFF_RTOL`` of the sum of
the magnitudes it was summed from) reads as zero. A larger negative total raises
``NegativeTotalError`` naming the component; nothing is changed to make it fit. Run on the user-defined model of
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

    def test_a_negative_total_within_roundoff_of_large_amounts_reads_as_zero(self):
        """Beyond the absolute tolerance but within ROUNDOFF_RTOL of the forms."""
        big = 1e6
        assert checked_total("X", [big, -big * (1.0 + 0.5 * ROUNDOFF_RTOL)], engine="E") == 0.0

    def test_a_negative_total_beyond_roundoff_of_large_amounts_raises(self):
        big = 1e6
        with pytest.raises(NegativeTotalError, match="'X'"):
            checked_total("X", [big, -big * (1.0 + 10.0 * ROUNDOFF_RTOL)], engine="E")

    def test_a_negative_total_within_the_absolute_tolerance_reads_as_zero(self):
        """An integrator taking a species to zero leaves trial totals within its
        tolerance of zero: -4e-18 mol of a total whose forms are 1.5e-17 mol."""
        assert checked_total("X", [5.7e-18, -9.7e-18], engine="E") == 0.0
        assert checked_total("Na+", [-0.5e-9], engine="E") == 0.0

    def test_a_negative_total_beyond_the_absolute_tolerance_raises(self):
        with pytest.raises(NegativeTotalError, match="negative_total_tolerance_mol"):
            checked_total("Na+", [-2e-9], engine="E")

    def test_the_absolute_tolerance_is_configurable(self):
        import PyOMES
        saved = PyOMES.config.negative_total_tolerance_mol
        try:
            PyOMES.config.negative_total_tolerance_mol = 1e-6
            assert checked_total("Na+", [-5e-7], engine="E") == 0.0
            PyOMES.config.negative_total_tolerance_mol = 0.0
            with pytest.raises(NegativeTotalError):
                checked_total("Na+", [-1e-30], engine="E")
        finally:
            PyOMES.config.negative_total_tolerance_mol = saved

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
        cv.phases["liquid"].n_mol["Na+"] = -1e-6
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


class TestAdaptiveSolverRetries:
    """An engine raising during an adaptive step makes the solver integrate the
    step again from its start with smaller maximum steps."""

    @staticmethod
    def _engine_raising(cv, calls_to_fail):
        """Make the engine raise NegativeTotalError on the given call numbers."""
        engine = cv.reaction_system.engine
        real = engine.solve
        count = {"n": 0}

        def solve(*args, **kwargs):
            count["n"] += 1
            if count["n"] in calls_to_fail:
                raise NegativeTotalError("Bisection engine: the total of 'X' is -1e-6 mol.")
            return real(*args, **kwargs)

        engine.solve = solve
        return count

    def _step(self, cv, solver):
        cv.advance(dt_h=0.0)
        cv.advance(dt_h=0.05, t_h=0.0, solver=solver)

    @pytest.mark.parametrize("method", ["BDF", "Radau", "LSODA", "RK45"])
    def test_a_failed_trial_is_retried_and_the_step_completes(self, method):
        import warnings
        from PyOMES import AccuracyWarning
        from PyOMES.core import SimultaneousAdaptiveSolver

        reference = _cv("bisection")
        self._step(reference, SimultaneousAdaptiveSolver(method=method))

        cv = _cv("bisection")
        cv.advance(dt_h=0.0)
        self._engine_raising(cv, calls_to_fail={3})
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cv.advance(dt_h=0.05, t_h=0.0, solver=SimultaneousAdaptiveSolver(method=method))
        retries = [str(w.message) for w in caught
                   if issubclass(w.category, AccuracyWarning) and "again from its start" in str(w.message)]
        assert len(retries) == 1
        assert "'X'" in retries[0] and "attempt 1 of" in retries[0]
        for sp in ("AceticAcid", "O2", "CO2"):
            assert cv.phases["liquid"].n_mol[sp] == pytest.approx(
                reference.phases["liquid"].n_mol[sp], rel=1e-4)

    def test_a_total_that_stays_negative_raises_after_the_retries(self):
        import warnings
        from PyOMES.core import SimultaneousAdaptiveSolver

        cv = _cv("bisection")
        cv.advance(dt_h=0.0)
        before = dict(cv.phases["liquid"].n_mol)
        self._engine_raising(cv, calls_to_fail=set(range(3, 10_000)))
        solver = SimultaneousAdaptiveSolver(method="BDF")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(NegativeTotalError) as info:
                cv.advance(dt_h=0.05, t_h=0.0, solver=solver)
        notes = " ".join(getattr(info.value, "__notes__", []))
        assert f"{solver.NEGATIVE_TOTAL_RETRIES} times" in notes
        assert "the model, not the step size" in notes
        for sp in ("AceticAcid", "O2", "CO2"):
            assert cv.phases["liquid"].n_mol[sp] == before[sp]
