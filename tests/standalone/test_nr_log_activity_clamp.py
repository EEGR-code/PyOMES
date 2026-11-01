# -*- coding: utf-8 -*-
"""Regression tests for the Newton solver's log-activity clamp.

``_compute_concentrations`` clamps each log10-activity to
``[-max_log_activity, +max_log_activity]`` before exponentiating, and
``solve_nr`` pins the master of an absent component at ``-max_log_activity``.
The clamp is written with Python's ``min``/``max`` on floats, so these tests
state the behaviour any replacement has to keep: an array version would give
``inf`` and a NumPy warning on overflow, and reversing the arguments to
``min``/``max`` would turn NaN into a limit.
"""
from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

from PyOMES.chemical_equilibrium.engines.nr import solver
from PyOMES.chemical_equilibrium.engines.nr.engine import NRChemicalEquilibriumEngine
from PyOMES.chemical_equilibrium.engines.nr.solver import (
    NRSolverCache,
    _compute_concentrations,
    solve_nr,
)
from PyOMES.chemical_equilibrium.engines.nr.tableau import NRTableau, build_tableau
from PyOMES.databases.aqueous import CO2, CO3_2minus, H2O, H_plus, HCO3_minus, NH3, NH4_plus, OH_minus
from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
from PyOMES.thermo import IdealLiquidModel


def _e(species, coeff):
    return StoichiometryEntry(species=species, phase="liquid", coefficient=coeff)


def _reactions():
    """Water, a two-step carbonate ladder and ammonia: masters H+, CO2, NH3."""
    return [
        EquilibriumReaction(
            stoichiometry=[_e(H2O, -1), _e(H_plus, +1), _e(OH_minus, +1)],
            log_K=-14.0, label="water",
        ),
        EquilibriumReaction(
            stoichiometry=[_e(CO2, -1), _e(H2O, -1), _e(HCO3_minus, +1), _e(H_plus, +1)],
            log_K=-6.35, total_id="CO2", label="co2_first",
        ),
        EquilibriumReaction(
            stoichiometry=[_e(HCO3_minus, -1), _e(CO3_2minus, +1), _e(H_plus, +1)],
            log_K=-10.33, total_id="CO2", label="co2_second",
        ),
        EquilibriumReaction(
            stoichiometry=[_e(NH4_plus, -1), _e(NH3, +1), _e(H_plus, +1)],
            log_K=-9.25, total_id="NH3", label="nh4",
        ),
    ]


@pytest.fixture(scope="module")
def tableau() -> NRTableau:
    return build_tableau(_reactions(), T_K=298.15)


def _x(tableau, **log_activities):
    """A log-activity vector with each named master set, the others at -3."""
    x = np.full(len(tableau.masters), -3.0)
    for name, value in log_activities.items():
        x[tableau.masters.index({"H": "H+"}.get(name, name))] = value
    return x


def _raw_log_activities(tableau, x):
    """Every species' log10-activity at *x*, before any clamp."""
    return [float(value) for value in x] + [
        sec.log_K_prime + sum(nu * x[tableau.masters.index(k)] for k, nu in sec.nu.items())
        for sec in tableau.secondaries
    ]


class TestExtremeValuesAreClamped:
    """Raw values as extreme as a cold-start Newton step produces."""

    @pytest.mark.parametrize("raw, expected", [(1545.0, 1e50), (-539.0, 1e-50)])
    def test_master_is_held_at_the_limit(self, tableau, raw, expected):
        c = _compute_concentrations(tableau, _x(tableau, CO2=raw), {})
        assert c["CO2"] == pytest.approx(expected, rel=1e-12, abs=0)

    def test_secondary_is_held_at_the_limit(self, tableau):
        # CO3-- = CO2 - 2 H+ + log K', so x[H+] = -800 sends it far above the limit.
        c = _compute_concentrations(tableau, _x(tableau, H=-800.0), {})
        assert c["CO3--"] == pytest.approx(1e50, rel=1e-12, abs=0)
        assert c["H+"] == pytest.approx(1e-50, rel=1e-12, abs=0)

    def test_extreme_values_stay_finite_and_raise_no_warning(self, tableau):
        x = _x(tableau, H=1545.0, CO2=-539.0, NH3=1545.0)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            c = _compute_concentrations(tableau, x, {})
        assert all(math.isfinite(value) for value in c.values())

    def test_activity_coefficient_divides_the_clamped_activity(self, tableau):
        c = _compute_concentrations(tableau, _x(tableau, CO2=1545.0), {"CO2": 0.5})
        assert c["CO2"] == pytest.approx(2e50, rel=1e-12, abs=0)

    def test_values_inside_the_limits_are_untouched(self, tableau):
        c = _compute_concentrations(tableau, _x(tableau, CO2=-49.0), {})
        assert c["CO2"] == 10.0 ** -49.0

    def test_limit_argument_sets_the_clamp(self, tableau):
        c = _compute_concentrations(
            tableau, _x(tableau, CO2=1545.0, NH3=-539.0), {}, max_log_activity=20.0,
        )
        assert c["CO2"] == pytest.approx(1e20, rel=1e-12, abs=0)
        assert c["NH3"] == pytest.approx(1e-20, rel=1e-12, abs=0)

    def test_limit_too_large_for_a_float_overflows(self, tableau):
        """Why the limit needs an upper bound: 10.0 ** 309 is not a float."""
        with pytest.raises(OverflowError):
            _compute_concentrations(
                tableau, _x(tableau, CO2=1545.0), {}, max_log_activity=1000.0,
            )

    def test_full_solve_survives_clamped_trial_steps(self, tableau, monkeypatch):
        """A cold start in strong acid clamps on trial steps, then converges cleanly."""
        raw_peaks = []

        def spy(tableau_, x, gammas, max_log_activity=50.0):
            raw_peaks.append(max(_raw_log_activities(tableau_, x)))
            return _compute_concentrations(tableau_, x, gammas, max_log_activity=max_log_activity)

        monkeypatch.setattr(solver, "_compute_concentrations", spy)
        totals = {"CO2": 1e-3, "NH3": 1e-3}
        out = solve_nr(tableau, totals, {"CT_Cl": 1.0}, activity_model=IdealLiquidModel())

        assert max(raw_peaks) > 50.0, "this case no longer clamps, so it tests nothing"
        assert raw_peaks[-1] <= 50.0
        assert out["pH"] == pytest.approx(0.0, abs=1e-2)
        assert out["CO2"] + out["HCO3-"] + out["CO3--"] == pytest.approx(1e-3, rel=1e-8)
        assert out["NH3"] + out["NH4+"] == pytest.approx(1e-3, rel=1e-8)


class TestNaNStaysNaN:
    """A failed solve has to stay visibly failed, not become a plausible limit."""

    def test_nan_master_stays_nan(self, tableau):
        c = _compute_concentrations(tableau, _x(tableau, CO2=float("nan")), {})
        assert math.isnan(c["CO2"])

    def test_nan_master_reaches_its_secondaries(self, tableau):
        c = _compute_concentrations(tableau, _x(tableau, CO2=float("nan")), {})
        assert math.isnan(c["HCO3-"])
        assert math.isnan(c["CO3--"])

    def test_nan_leaves_other_components_alone(self, tableau):
        c = _compute_concentrations(tableau, _x(tableau, CO2=float("nan")), {})
        assert c["NH3"] == pytest.approx(1e-3, rel=1e-12, abs=0)


class TestZeroTotalComponentIsPinned:
    """An absent component's master sits at -max_log_activity, not at -infinity."""

    @staticmethod
    def _solve(tableau, totals, max_log_activity=50.0):
        cache = NRSolverCache()
        out = solve_nr(
            tableau, totals, {}, activity_model=IdealLiquidModel(), cache=cache,
            max_log_activity=max_log_activity,
        )
        return out, cache

    def test_master_is_pinned_at_the_limit(self, tableau):
        out, cache = self._solve(tableau, {"CO2": 0.05, "NH3": 0.0})
        assert cache.log_x[tableau.masters.index("NH3")] == -50.0
        assert out["NH3"] == pytest.approx(1e-50, rel=1e-12, abs=0)

    def test_pin_follows_the_limit(self, tableau):
        out, cache = self._solve(tableau, {"CO2": 0.05, "NH3": 0.0}, max_log_activity=100.0)
        assert cache.log_x[tableau.masters.index("NH3")] == -100.0
        assert out["NH3"] == pytest.approx(1e-100, rel=1e-12, abs=0)

    def test_missing_total_is_pinned_like_a_zero_one(self, tableau):
        _, cache = self._solve(tableau, {"CO2": 0.05})
        assert cache.log_x[tableau.masters.index("NH3")] == -50.0

    def test_total_below_the_threshold_is_pinned(self, tableau):
        _, cache = self._solve(tableau, {"CO2": 0.05, "NH3": 1e-21})
        assert cache.log_x[tableau.masters.index("NH3")] == -50.0

    def test_secondary_of_a_pinned_master_is_negligible(self, tableau):
        out, _ = self._solve(tableau, {"CO2": 0.05, "NH3": 0.0})
        assert out["NH4+"] < 1e-30

    def test_pinned_component_does_not_move_the_answer(self, tableau):
        """Leaving NH3 out altogether would give the same carbonate result."""
        reactions = [r for r in _reactions() if r.label != "nh4"]
        without = solve_nr(
            build_tableau(reactions, T_K=298.15), {"CO2": 0.05}, {},
            activity_model=IdealLiquidModel(),
        )
        pinned, _ = self._solve(tableau, {"CO2": 0.05, "NH3": 0.0})
        assert pinned["pH"] == pytest.approx(without["pH"], abs=1e-12)

    def test_component_recovers_when_its_total_returns(self, tableau):
        """A warmstart holding the pin must not stall the next, non-zero solve."""
        cache = NRSolverCache()
        solve_nr(tableau, {"CO2": 0.05, "NH3": 0.0}, {}, activity_model=IdealLiquidModel(), cache=cache)
        out = solve_nr(tableau, {"CO2": 0.05, "NH3": 0.04}, {}, activity_model=IdealLiquidModel(), cache=cache)
        assert out["NH3"] + out["NH4+"] == pytest.approx(0.04, rel=1e-8)

    def test_engine_reports_the_pinned_value(self):
        engine = NRChemicalEquilibriumEngine.from_reactions(_reactions())
        out = engine.solve(totals={"CO2": 0.05, "NH3": 0.0})
        assert out.species_mol_L["NH3"] == pytest.approx(1e-50, rel=1e-12, abs=0)


class TestUnsafeLimitsAreRejected:
    """The engine refuses a limit that cannot work, before any solve runs."""

    @pytest.mark.parametrize(
        "limit", [0.0, -1.0, -50.0, 150.5, 300.0, 1000.0, float("inf"), float("-inf"), float("nan")],
    )
    def test_constructor_rejects(self, tableau, limit):
        with pytest.raises(ValueError, match="max_log_activity must be positive and at most 150"):
            NRChemicalEquilibriumEngine(tableau, max_log_activity=limit)

    def test_from_reactions_rejects(self):
        with pytest.raises(ValueError, match="max_log_activity"):
            NRChemicalEquilibriumEngine.from_reactions(_reactions(), max_log_activity=1000.0)

    def test_message_names_the_value(self, tableau):
        with pytest.raises(ValueError, match="got 1000.0"):
            NRChemicalEquilibriumEngine(tableau, max_log_activity=1000.0)

    @pytest.mark.parametrize("limit", [1e-3, 20.0, 50.0, 100.0, 150.0, 150])
    def test_constructor_accepts(self, tableau, limit):
        engine = NRChemicalEquilibriumEngine(tableau, max_log_activity=limit)
        assert engine.max_log_activity == float(limit)

    def test_default_limit_is_accepted(self, tableau):
        assert NRChemicalEquilibriumEngine(tableau).max_log_activity == 50.0

    def test_largest_accepted_limit_clamps_without_numeric_warnings(self, tableau, monkeypatch):
        """Trial steps pass the ceiling here; a higher limit overflows the line search."""
        raw_peaks = []

        def spy(tableau_, x, gammas, max_log_activity=50.0):
            raw_peaks.append(max(_raw_log_activities(tableau_, x)))
            return _compute_concentrations(tableau_, x, gammas, max_log_activity=max_log_activity)

        monkeypatch.setattr(solver, "_compute_concentrations", spy)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            out = solve_nr(
                tableau, {"CO2": 1e-6, "NH3": 1e-3}, {"CT_Na": 1.0},
                activity_model=IdealLiquidModel(), max_log_activity=150.0,
            )
        assert max(raw_peaks) > 150.0, "this case no longer clamps, so it tests nothing"
        assert out["pH"] == pytest.approx(14.0, abs=1e-2)
