# -*- coding: utf-8 -*-
"""Regression tests for the Newton solver's log-activity clamp.

``_compute_concentrations`` clamps each log10-activity to
``[-max_log_activity, +max_log_activity]`` before exponentiating, and
``solve_nr`` pins the master of an absent component at ``-max_log_activity``.
The clamp is written with Python's ``min``/``max`` on floats, so these tests
state the behaviour any replacement has to keep: an array version would give
``inf`` and a NumPy warning on overflow, and reversing the arguments to
``min``/``max`` would turn NaN into a limit. They also cover the engine's limit
on ``max_log_activity`` and its report of when the limit or the pin affects a
final answer.
"""
from __future__ import annotations

import math
import warnings
from types import SimpleNamespace

import numpy as np
import pytest

from PyOMES.chemical_equilibrium.engines.nr import engine as engine_module
from PyOMES.chemical_equilibrium.engines.nr import solver
from PyOMES.chemical_equilibrium.engines.nr.engine import NRChemicalEquilibriumEngine
from PyOMES.chemical_equilibrium.engines.nr.solver import (
    NRSolverCache,
    _compute_concentrations,
    _limit_effects,
    solve_nr,
)
from PyOMES.chemical_equilibrium.engines.nr.tableau import NRTableau, SecondaryEntry, build_tableau
from PyOMES.chemistry.species import Species
from PyOMES.databases.aqueous import CO2, CO3_2minus, H2O, H_plus, HCO3_minus, NH3, NH4_plus, OH_minus
from PyOMES.databases.bioprocess_basic import Ca_2plus
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


CLAMPED = "clamped_log_activities"
ABSENT = "absent_component_species"


def _forged_tableau(*secondaries, masters=("H+", "A")):
    """The parts of a tableau that _limit_effects reads."""
    return SimpleNamespace(masters=list(masters), secondaries=list(secondaries))


def _pin_on(tableau, master):
    return [(0, tableau.masters.index(master), -50.0)]


class TestLimitEffectsDetection:
    """_limit_effects on a forged final state, so each rule is tested on its own."""

    def test_high_master_is_reported_with_its_raw_value(self):
        clamped, absent = _limit_effects(
            _forged_tableau(), np.array([-7.0, 60.5]), {"H+": 1e-7, "A": 1e50}, [], 50.0, 1e-10,
        )
        assert clamped == {"A": 60.5}
        assert absent == {}

    def test_high_secondary_is_reported_with_its_raw_value(self):
        sec = SecondaryEntry("S", 0, {"A": 1.0}, 5.0)
        clamped, _ = _limit_effects(
            _forged_tableau(sec), np.array([-7.0, 50.0]),
            {"H+": 1e-7, "A": 1e50, "S": 1e50}, [], 50.0, 1e-10,
        )
        assert clamped == {"S": 55.0}

    def test_values_at_the_limit_are_not_reported(self):
        sec = SecondaryEntry("S", 0, {"A": 1.0}, 0.0)
        clamped, _ = _limit_effects(
            _forged_tableau(sec), np.array([-7.0, -50.0]),
            {"H+": 1e-7, "A": 1e-50, "S": 1e-50}, [], 50.0, 1e-10,
        )
        assert clamped == {}

    def test_low_values_of_an_active_component_are_reported(self):
        sec = SecondaryEntry("S", 0, {"A": 1.0}, -3.0)
        clamped, _ = _limit_effects(
            _forged_tableau(sec), np.array([-7.0, -60.0]),
            {"H+": 1e-7, "A": 1e-50, "S": 1e-50}, [], 50.0, 1e-10,
        )
        assert clamped == {"A": -60.0, "S": -63.0}

    def test_low_values_due_to_the_pin_are_not_reported(self):
        sec = SecondaryEntry("S", 0, {"A": 1.0}, -3.0)
        tab = _forged_tableau(sec)
        clamped, absent = _limit_effects(
            tab, np.array([-7.0, -50.0000001]),
            {"H+": 1e-7, "A": 1e-50, "S": 1e-50}, _pin_on(tab, "A"), 50.0, 1e-10,
        )
        assert clamped == {}
        assert absent == {}

    def test_low_value_is_reported_when_the_pin_is_not_what_lowers_it(self):
        """A pinned master with a negative coefficient pushes the species up."""
        sec = SecondaryEntry("S", 0, {"A": -1.0}, -200.0)
        tab = _forged_tableau(sec)
        clamped, _ = _limit_effects(
            tab, np.array([-7.0, -50.0]),
            {"H+": 1e-7, "A": 1e-50, "S": 1e-50}, _pin_on(tab, "A"), 50.0, 1e-10,
        )
        assert clamped == {"S": -150.0}

    def test_high_value_of_a_pinned_component_is_reported(self):
        sec = SecondaryEntry("S", 0, {"A": -1.0}, 5.0)
        tab = _forged_tableau(sec)
        clamped, _ = _limit_effects(
            tab, np.array([-7.0, -50.0]),
            {"H+": 1e-7, "A": 1e-50, "S": 1e50}, _pin_on(tab, "A"), 50.0, 1e-10,
        )
        assert clamped == {"S": 55.0}

    def test_nan_is_not_reported_as_a_limit_hit(self):
        clamped, absent = _limit_effects(
            _forged_tableau(), np.array([-7.0, float("nan")]),
            {"H+": 1e-7, "A": float("nan")}, [], 50.0, 1e-10,
        )
        assert clamped == {}
        assert absent == {}

    def test_absent_report_is_strictly_above_the_tolerance(self):
        sec = SecondaryEntry("S", 0, {"A": 1.0}, 0.0)
        tab = _forged_tableau(sec)
        x = np.array([-7.0, -5.0])
        at_tol = _limit_effects(tab, x, {"H+": 1e-7, "A": 1e-10, "S": 1e-10},
                                _pin_on(tab, "A"), 50.0, 1e-10)
        above = _limit_effects(tab, x, {"H+": 1e-7, "A": 2e-10, "S": 3e-10},
                               _pin_on(tab, "A"), 50.0, 1e-10)
        assert at_tol[1] == {}
        assert above[1] == {"A": 2e-10, "S": 3e-10}

    def test_species_of_an_active_component_are_never_absent(self):
        sec = SecondaryEntry("S", 0, {"A": 1.0}, 0.0)
        _, absent = _limit_effects(
            _forged_tableau(sec), np.array([-7.0, -1.0]),
            {"H+": 1e-7, "A": 0.1, "S": 0.1}, [], 50.0, 1e-10,
        )
        assert absent == {}

    def test_gas_secondary_is_keyed_apart_from_its_liquid_parent(self):
        sec = SecondaryEntry("A", 0, {"A": 1.0}, 0.0, phase="gas")
        clamped, _ = _limit_effects(
            _forged_tableau(sec), np.array([-7.0, -60.0]),
            {"H+": 1e-7, "A": 1e-50, "A:gas": 1e-50}, [], 50.0, 1e-10,
        )
        assert clamped == {"A": -60.0, "A:gas": -60.0}


class TestEngineReportsLimitEffects:
    """What a caller of NRChemicalEquilibriumEngine.solve() sees."""

    # 1 mol/L strong acid: CO3-- has a raw log10-activity of about -19.7.
    ACID = {"totals": {"CO2": 1e-3, "NH3": 1e-3}, "strong_ions": {"CT_Cl": 1.0}}
    # Strong base with no carbon: at a limit of 20, CO3-- rises above the
    # tolerance from a zero CO2 total.
    BASE = {"totals": {"CO2": 0.0, "NH3": 0.1}, "strong_ions": {"CT_Na": 1.0}}

    @staticmethod
    def _engine(**kwargs):
        return NRChemicalEquilibriumEngine.from_reactions(_reactions(), **kwargs)

    def test_ordinary_solve_reports_nothing(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            out = self._engine().solve(totals={"CO2": 0.05, "NH3": 0.04})
        assert out.extra == {}

    def test_clamping_only_on_trial_steps_reports_nothing(self):
        """The strong-acid cold start clamps on trial steps: see the test above."""
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            out = self._engine(use_warmstart=False).solve(**self.ACID)
        assert out.extra == {}

    def test_absent_component_below_the_limit_reports_nothing(self, tableau):
        """At high pH, NH4+ of an absent NH3 falls below the limit: no false alarm."""
        cache = NRSolverCache()
        solve_nr(
            tableau, {"CO2": 0.05, "NH3": 0.0}, {"CT_Na": 1.0},
            activity_model=IdealLiquidModel(), cache=cache,
        )
        raw = _raw_log_activities(tableau, cache.log_x)
        assert min(raw) < -50.0, "this case no longer reaches the limit, so it tests nothing"

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            out = self._engine().solve(
                totals={"CO2": 0.05, "NH3": 0.0}, strong_ions={"CT_Na": 1.0},
            )
        assert out.extra == {}

    def test_clamped_trace_species_is_reported(self):
        with pytest.warns(RuntimeWarning, match=r"\['CO3--'\] are held at the limit"):
            out = self._engine(max_log_activity=15.0).solve(**self.ACID)
        assert set(out.extra) == {CLAMPED}
        assert set(out.extra[CLAMPED]) == {"CO3--"}
        assert out.extra[CLAMPED]["CO3--"] == pytest.approx(-19.7, abs=0.1)
        assert out.species_mol_L["CO3--"] == pytest.approx(1e-15, rel=1e-12, abs=0)

    def test_absent_component_species_is_reported(self):
        with pytest.warns(RuntimeWarning, match=r"\['CO3--'\] belong to absent components"):
            out = self._engine(max_log_activity=20.0).solve(**self.BASE)
        assert set(out.extra) == {ABSENT}
        assert set(out.extra[ABSENT]) == {"CO3--"}
        assert out.extra[ABSENT]["CO3--"] == out.species_mol_L["CO3--"]
        assert out.extra[ABSENT]["CO3--"] > 1e-10

    def test_default_limit_keeps_the_same_species_below_the_tolerance(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            out = self._engine().solve(**self.BASE)
        assert out.extra == {}
        assert out.species_mol_L["CO3--"] < 1e-10

    def test_report_does_not_change_the_answer(self, tableau):
        engine = self._engine(max_log_activity=15.0, use_warmstart=False)
        with pytest.warns(RuntimeWarning):
            reported = engine.solve(**self.ACID)
        plain = solve_nr(
            tableau, self.ACID["totals"], self.ACID["strong_ions"],
            activity_model=IdealLiquidModel(), max_log_activity=15.0,
        )
        assert reported.pH == plain["pH"]
        assert reported.species_mol_L["CO3--"] == plain["CO3--"]

    def test_warns_once_per_species_but_reports_every_solve(self):
        engine = self._engine(max_log_activity=15.0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            results = [engine.solve(**self.ACID) for _ in range(3)]
            engine.reset_cache()
            results.append(engine.solve(**self.ACID))
        assert len(caught) == 1
        assert all("CO3--" in result.extra[CLAMPED] for result in results)

    def test_warns_again_for_a_new_species(self):
        """H+ near pH 14 is beyond a limit of 13, which also holds the reported pH at 13."""
        engine = self._engine(max_log_activity=13.0)
        with pytest.warns(RuntimeWarning, match="CO3--"):
            engine.solve(**self.ACID)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            out = engine.solve(totals={"CO2": 1e-3, "NH3": 1e-3}, strong_ions={"CT_Na": 1.0})
        assert len(caught) == 1
        assert "'H+'" in str(caught[0].message)
        assert "CO3--" not in str(caught[0].message)
        assert out.extra[CLAMPED]["H+"] == pytest.approx(-14.0, abs=0.1)
        assert out.pH == pytest.approx(13.0, abs=1e-9)

    def test_each_engine_warns_for_itself(self):
        for _ in range(2):
            with pytest.warns(RuntimeWarning, match="held at the limit"):
                self._engine(max_log_activity=15.0).solve(**self.ACID)

    def test_warning_names_the_limit(self):
        with pytest.warns(RuntimeWarning, match="max_log_activity=15"):
            self._engine(max_log_activity=15.0).solve(**self.ACID)

    def test_private_keys_do_not_reach_the_result(self):
        engine = self._engine(max_log_activity=15.0, retain_jacobian=True)
        with pytest.warns(RuntimeWarning):
            out = engine.solve(**self.ACID)
        assert not [key for key in {**out.species_mol_L, **out.extra} if key.startswith("_")]

    def test_solve_nr_reports_only_when_asked(self, tableau):
        args = (tableau, self.ACID["totals"], self.ACID["strong_ions"])
        kwargs = {"activity_model": IdealLiquidModel(), "max_log_activity": 15.0}
        assert "_clamped_log_activities" not in solve_nr(*args, **kwargs)
        asked = solve_nr(*args, report_limit_effects=True, **kwargs)
        assert "CO3--" in asked["_clamped_log_activities"]

    def test_final_evaluation_only(self, monkeypatch):
        """The check stays out of the Newton iterations and the line search."""
        checks = []
        evaluations = []
        original_check = solver._limit_effects

        def count_check(*args, **kwargs):
            checks.append(len(evaluations))
            return original_check(*args, **kwargs)

        def count_evaluation(*args, **kwargs):
            evaluations.append(1)
            return _compute_concentrations(*args, **kwargs)

        monkeypatch.setattr(solver, "_limit_effects", count_check)
        monkeypatch.setattr(solver, "_compute_concentrations", count_evaluation)
        self._engine(use_warmstart=False).solve(**self.ACID)
        assert len(checks) == 1
        assert checks[0] == len(evaluations) > 1


class TestPrecipitationReportsOnce:
    """The active-set loop runs solve_nr many times; only the final solve reports."""

    @staticmethod
    def _engine(**kwargs):
        calcite = EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(
                    species=Species(id="CaCO3", atoms={"Ca": 1, "C": 1, "O": 3}, charge=0),
                    phase="solid", coefficient=-1,
                ),
                _e(Ca_2plus, +1),
                _e(CO3_2minus, +1),
            ],
            log_K=-8.48, label="calcite",
        )
        return NRChemicalEquilibriumEngine.from_reactions(_reactions() + [calcite], **kwargs)

    def test_only_the_final_solve_is_checked(self, monkeypatch):
        calls = []
        original_solve = solver.solve_nr

        def count_solve(*args, **kwargs):
            calls.append(kwargs["report_limit_effects"])
            return original_solve(*args, **kwargs)

        monkeypatch.setattr(engine_module, "solve_nr", count_solve)
        out = self._engine().solve(
            totals={"CO2": 0.010, "NH3": 1e-3}, strong_ions={"CT_Ca": 0.002, "CT_Na": 0.005},
        )
        assert out.extra["minerals_xi_mol_L"]["calcite"] > 0.0
        assert len(calls) > 1, "this case no longer iterates, so it tests nothing"
        assert calls.count(True) == 1
        assert calls[-1] is True

    def test_report_sits_beside_the_mineral_amounts(self):
        with pytest.warns(RuntimeWarning, match="held at the limit"):
            out = self._engine(max_log_activity=15.0).solve(
                totals={"CO2": 1e-3, "NH3": 1e-3}, strong_ions={"CT_Cl": 1.0, "CT_Ca": 1e-3},
            )
        assert out.extra["minerals_xi_mol_L"]["calcite"] == 0.0
        assert "CO3--" in out.extra[CLAMPED]

    def test_ordinary_precipitation_reports_nothing_new(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            out = self._engine().solve(
                totals={"CO2": 0.010, "NH3": 1e-3}, strong_ions={"CT_Ca": 0.002, "CT_Na": 0.005},
            )
        assert set(out.extra) == {"minerals_xi_mol_L"}
