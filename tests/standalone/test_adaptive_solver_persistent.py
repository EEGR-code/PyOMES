# -*- coding: utf-8 -*-
"""SimultaneousAdaptiveSolver keeping one integrator across macro steps.

Each restart trigger has a test here that fails if the trigger is removed.
"""
from __future__ import annotations

import pickle
import warnings

import numpy as np
import pytest

from PyOMES.chemistry import Species
from PyOMES.control.actions import ControlAction, ProfileRecord
from PyOMES.control.descriptors import MutableDict, MutableScalar, parameter_generation
from PyOMES.core import ControlVolume, LiquidPhase, Simulation
from PyOMES.core.links import AdvectiveLink
from PyOMES.core import solvers as solvers_mod
from PyOMES.core.solvers import SimultaneousAdaptiveSolver
from PyOMES.databases.aqueous import H2O, H_plus, OH_minus
from PyOMES.reactions import (
    EquilibriumReaction, KineticReaction, ReactionSystem, StoichiometryEntry,
)

ACID = Species(id="AcOH", atoms={"C": 2, "H": 4, "O": 2}, charge=0)
BASE = Species(id="Ac-", atoms={"C": 2, "H": 3, "O": 2}, charge=-1)
PRODUCT = Species(id="P", atoms={"C": 2, "H": 4, "O": 2}, charge=0)
SODIUM = Species(id="Na+", atoms={"Na": 1}, charge=1)
DT_H = 0.05


def _entry(species, coefficient):
    return StoichiometryEntry(species=species, phase="liquid", coefficient=coefficient)


def _make_cv(extra_mol=None):
    system = ReactionSystem(
        [
            KineticReaction(
                stoichiometry=[_entry(ACID, -1), _entry(PRODUCT, +1)],
                rate_fn=lambda env: 2.0 * env.S("AcOH") * env.V_L,
                label="decay",
            ),
            EquilibriumReaction(
                stoichiometry=[_entry(H2O, -1), _entry(H_plus, +1), _entry(OH_minus, +1)],
                log_K=-14.0, label="water",
            ),
            EquilibriumReaction(
                stoichiometry=[_entry(ACID, -1), _entry(BASE, +1), _entry(H_plus, +1)],
                log_K=-4.76, label="acetic",
            ),
        ],
        solver="newton_raphson",
    )
    n_mol = {"AcOH": 0.05, "Ac-": 0.0, "P": 0.0, "Na+": 0.01, "H+": 1e-5,
             "OH-": 0.0, "H2O": 55.5}
    n_mol.update(extra_mol or {})
    return ControlVolume(
        phases={"liquid": LiquidPhase(n_mol=n_mol, V_L=1.0, T_K=298.15)},
        reaction_system=system, label="weak_acid", species=[SODIUM, PRODUCT],
    )


def _solver(**kwargs):
    kwargs.setdefault("method", "BDF")
    kwargs.setdefault("persistent", True)
    return SimultaneousAdaptiveSolver(**kwargs)


@pytest.fixture
def steppers(monkeypatch):
    """Records each integrator made."""
    made = []
    real = solvers_mod._make_stepper

    def counting(*args, **kwargs):
        stepper = real(*args, **kwargs)
        made.append(stepper)
        return stepper

    monkeypatch.setattr(solvers_mod, "_make_stepper", counting)
    return made


@pytest.fixture
def builds(monkeypatch):
    """Records each restart: a new right-hand side built from the current state."""
    made = []
    real = SimultaneousAdaptiveSolver._build_integrand

    def counting(self, *args, **kwargs):
        integrand = real(self, *args, **kwargs)
        made.append(integrand)
        return integrand

    monkeypatch.setattr(SimultaneousAdaptiveSolver, "_build_integrand", counting)
    return made


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def _advance(cv, solver, n=1, **kwargs):
    """Step *cv*, keeping time running on from its last step."""
    dt_h = kwargs.get("dt_h", DT_H)
    for _ in range(n):
        t_h = cv.__dict__.setdefault("_test_t_h", 0.0)
        solver.solve_step(cv, dt_h=dt_h, t_h=t_h,
                          external_source_terms=kwargs.get("sources"))
        cv._test_t_h = t_h + dt_h


class TestKeepsOneIntegrator:
    def test_nothing_changed_means_one_integrator(self, builds):
        cv, solver = _make_cv(), _solver()
        _advance(cv, solver, n=5)
        assert len(builds) == 1

    def test_persistent_is_the_default(self, steppers):
        cv, solver = _make_cv(), SimultaneousAdaptiveSolver(method="BDF")
        assert solver.persistent is True
        _advance(cv, solver, n=2)
        assert steppers

    def test_not_persistent_builds_no_stepper(self, steppers):
        cv, solver = _make_cv(), _solver(persistent=False)
        _advance(cv, solver, n=3)
        assert steppers == []

    def test_freeze_speciation_restarts_each_step(self, steppers):
        cv, solver = _make_cv(), _solver(freeze_speciation=True)
        _advance(cv, solver, n=3)
        assert steppers == []

    @pytest.mark.parametrize("method", ["BDF", "Radau", "LSODA", "RK45", "RK23", "DOP853"])
    def test_agrees_with_restarting_within_tolerance(self, method):
        a, b = _make_cv(), _make_cv()
        _advance(a, _solver(method=method), n=20)
        _advance(b, _solver(method=method, persistent=False), n=20)
        for sp in ("P", "AcOH"):
            assert a.phases["liquid"].n_mol[sp] == pytest.approx(
                b.phases["liquid"].n_mol[sp], rel=1e-4,
            )

    def test_unknown_method_is_refused(self):
        cv = _make_cv()
        with pytest.raises(ValueError, match="Unknown method"):
            _advance(cv, _solver(method="nope"))

    def test_no_step_in_a_macro_step_still_moves_the_output(self):
        """Macro steps much shorter than the integrator's own, so some take no step."""
        cv, solver = _make_cv(), _solver()
        values = []
        for _ in range(30):
            _advance(cv, solver, dt_h=1e-4)
            values.append(cv.phases["liquid"].n_mol["P"])
        assert all(later > earlier for earlier, later in zip(values, values[1:]))


class TestTotalsAreKept:
    """State entries of resolved species are only parts of a total, so a flux can
    take one negative. Flooring it would create material."""

    @pytest.mark.parametrize("method", ["BDF", "Radau", "RK45"])
    def test_rate_law_on_a_resolved_species_conserves_carbon(self, method):
        cv, solver = _make_cv(), _solver(method=method)
        liquid = cv.phases["liquid"].n_mol
        _advance(cv, solver, dt_h=DT_H)
        total = sum(liquid[sp] for sp in ("AcOH", "Ac-", "P"))
        _advance(cv, solver, n=40)
        assert sum(liquid[sp] for sp in ("AcOH", "Ac-", "P")) == pytest.approx(total, rel=1e-5)
        assert all(liquid[sp] >= 0.0 for sp in ("AcOH", "Ac-", "P"))

    def test_drain_matches_washout(self):
        """With no reaction, a drain must take every total down as exp(-D t)."""
        from PyOMES.core.boundaries import LiquidDrain

        cv = _make_cv()
        cv.reaction_system._kinetic_reactions[0].rate_fn = lambda env: 0.0
        cv.boundaries.append(LiquidDrain(Q_L_per_h=2.0))
        solver = _solver()
        liquid = cv.phases["liquid"].n_mol
        _advance(cv, solver, dt_h=DT_H)
        start = liquid["AcOH"] + liquid["Ac-"]
        _advance(cv, solver, n=20)
        assert liquid["AcOH"] + liquid["Ac-"] == pytest.approx(
            start * np.exp(-2.0 * 20 * DT_H), rel=1e-4,
        )

    def test_a_species_the_engine_does_not_resolve_is_still_floored(self):
        """Na+ is a strong ion: its state entry is the amount itself."""
        cv, solver = _make_cv(), _solver()
        integrand = solver._build_integrand(cv, DT_H, None, keep_totals=True)
        y = -np.ones_like(integrand.y0)
        floored = integrand.floor(y)
        by_name = dict(zip((key[2] for key in integrand.sv._layout), floored))
        assert by_name["Na+"] == 0.0
        assert by_name["P"] == 0.0
        assert by_name["AcOH"] == -1.0


class TestRestartTriggers:
    def _two_steps(self, change, builds, **kwargs):
        cv, solver = _make_cv(), _solver()
        _advance(cv, solver)
        assert len(builds) == 1
        change(cv, solver)
        _advance(cv, solver, **kwargs)
        return len(builds)

    def test_amount_changed(self, builds):
        def dose(cv, solver):
            cv.phases["liquid"].n_mol["Na+"] += 1e-3
        assert self._two_steps(dose, builds) == 2

    def test_a_very_small_amount_change_counts(self, builds):
        def nudge(cv, solver):
            n_mol = cv.phases["liquid"].n_mol
            n_mol["H2O"] = np.nextafter(n_mol["H2O"], 100.0)
        assert self._two_steps(nudge, builds) == 2

    def test_temperature_changed(self, builds):
        def warm(cv, solver):
            cv.phases["liquid"].T_K = 303.15
        assert self._two_steps(warm, builds) == 2

    def test_volume_changed(self, builds):
        def grow(cv, solver):
            cv.phases["liquid"].V_L = 1.1
        assert self._two_steps(grow, builds) == 2

    def test_species_added(self, builds):
        def add(cv, solver):
            cv.phases["liquid"].n_mol["Na+x"] = 0.0
        assert self._two_steps(add, builds) == 2

    def test_dt_changed(self, builds):
        assert self._two_steps(lambda cv, solver: None, builds, dt_h=DT_H / 2) == 2

    def test_sources_changed(self, builds):
        sources = {"liquid": {"P": 1e-3}}
        assert self._two_steps(lambda cv, solver: None, builds, sources=sources) == 2

    def test_same_sources_do_not_restart(self, builds):
        cv, solver = _make_cv(), _solver()
        sources = {"liquid": {"P": 1e-3}}
        _advance(cv, solver, n=3, sources=sources)
        assert len(builds) == 1

    def test_engine_replaced(self, builds):
        def swap(cv, solver):
            old = cv.reaction_system.engine
            cv.reaction_system.attach_engine(type(old)(old.tableau, T_C=old.T_C))
        assert self._two_steps(swap, builds) == 2

    def test_species_added_by_the_closing_solve(self, builds):
        """The snapshot is taken after the closing solve, so it cannot see this. The
        state vector was built without the species the solve then writes."""
        cv = _make_cv()
        liquid = cv.phases["liquid"].n_mol
        del liquid["Ac-"], liquid["OH-"]
        solver = _solver()
        _advance(cv, solver)
        assert "Ac-" in liquid
        _advance(cv, solver)
        assert len(builds) == 2

    def test_reaction_system_replaced_with_the_same_engine(self, builds):
        """The engine check cannot see this: only the rate laws differ."""
        def swap(cv, solver):
            system = _make_cv().reaction_system
            system.attach_engine(cv.reaction_system.engine)
            cv._set_reaction_system_unchecked(system)
        assert self._two_steps(swap, builds) == 2

    def test_time_jump_restarts(self, builds):
        cv, solver = _make_cv(), _solver()
        _advance(cv, solver, n=2)
        cv._test_t_h += 5.0
        _advance(cv, solver)
        assert len(builds) == 2

    def test_time_following_on_does_not_restart(self, builds):
        """Simulation's grid comes from linspace, so it differs by rounding."""
        cv, solver = _make_cv(), _solver()
        for t_h in np.linspace(0.0, 1.0, 21)[:-1]:
            solver.solve_step(cv, dt_h=0.05, t_h=float(t_h))
        assert len(builds) == 1

    def test_explicit_restart(self, builds):
        assert self._two_steps(lambda cv, solver: solver.restart(), builds) == 2

    def test_restart_starts_from_the_changed_state(self):
        """Dosing inert product mid-run must appear in the result, not be overwritten."""
        cv, solver = _make_cv(), _solver()
        _advance(cv, solver, n=2)
        before = cv.phases["liquid"].n_mol["P"]
        cv.phases["liquid"].n_mol["P"] = before + 0.5
        _advance(cv, solver)
        assert cv.phases["liquid"].n_mol["P"] > before + 0.5


class _Knob:
    """A parameter holder of the kind boundaries and links are."""
    k = MutableScalar(float)
    table = MutableDict(value_type=float)

    def __init__(self):
        self.k = 1.0
        self.table = {"O2": 1.0}


class TestRateLawsReadSimulationTime:
    @staticmethod
    def _cv_recording_time(seen):
        cv = _make_cv()
        reaction = cv.reaction_system._kinetic_reactions[0]
        rate = reaction.rate_fn

        def timed(env):
            seen.append(env.t_h)
            return rate(env)

        reaction.rate_fn = timed
        return cv

    def test_time_runs_across_macro_steps(self):
        seen = []
        cv = self._cv_recording_time(seen)
        Simulation(cvs={"main": cv}, solver=_solver(), label="t").run(tau_h=1.0, n_steps=10)
        # Simulation passes the end of each step as t_h, so time starts at 0.1.
        assert min(seen) == pytest.approx(0.1)
        assert max(seen) >= 1.0

    def test_a_start_time_is_honoured(self):
        seen = []
        cv = self._cv_recording_time(seen)
        Simulation(cvs={"main": cv}, solver=_solver(), label="t").run(
            tau_h=1.0, n_steps=10, start_t_h=100.0,
        )
        assert min(seen) == pytest.approx(100.1)
        assert max(seen) >= 101.0

    def test_not_persistent_keeps_its_time_from_zero(self):
        seen = []
        cv = self._cv_recording_time(seen)
        Simulation(cvs={"main": cv}, solver=_solver(persistent=False), label="t").run(
            tau_h=1.0, n_steps=10,
        )
        assert min(seen) == 0.0
        assert max(seen) <= 0.1 + 1e-12

    def test_jacobian_columns_read_the_same_time(self):
        seen = []
        cv = self._cv_recording_time(seen)
        Simulation(cvs={"main": cv}, solver=_solver(use_engine_jacobian=True), label="t").run(
            tau_h=1.0, n_steps=10,
        )
        assert min(seen) >= 0.1 - 1e-12


class TestStopsAtTheStepEndAfterARestart:
    """A stepper past the macro-step end has to interpolate back, and may have reached
    trial states a floor acts on. Worth it only once the state has stopped changing."""

    def test_a_new_run_lands_exactly_on_the_step_end(self, steppers):
        cv, solver = _make_cv(), _solver()
        _advance(cv, solver)
        assert steppers[-1].t == DT_H

    def test_it_runs_freely_after_a_step_with_no_change(self, steppers):
        cv, solver = _make_cv(), _solver()
        _advance(cv, solver)
        _advance(cv, solver)
        assert len(steppers) == 2
        assert steppers[-1].t_bound > 1e6
        _advance(cv, solver, n=3)
        assert len(steppers) == 2

    def test_a_change_goes_back_to_landing_on_the_step_end(self, steppers):
        cv, solver = _make_cv(), _solver()
        _advance(cv, solver, n=3)
        cv.phases["liquid"].n_mol["Na+"] += 1e-3
        _advance(cv, solver)
        assert steppers[-1].t == DT_H
        assert steppers[-1].t_bound == DT_H

    def test_every_step_changing_never_runs_freely(self, steppers):
        cv, solver = _make_cv(), _solver()
        for _ in range(4):
            _advance(cv, solver)
            cv.phases["liquid"].n_mol["Na+"] += 1e-3
        assert len(steppers) == 4
        assert all(stepper.t_bound == DT_H for stepper in steppers)


class TestFailure:
    @staticmethod
    def _failing(stepper):
        def fail():
            stepper.status = "failed"
            return "forced"
        stepper.step = fail

    def _carbon(self, cv):
        liquid = cv.phases["liquid"].n_mol
        return sum(liquid[sp] for sp in ("AcOH", "Ac-", "P"))

    def test_a_continuing_run_retries_from_a_new_integrator(self, builds, steppers):
        cv, solver = _make_cv(), _solver()
        _advance(cv, solver, n=3)
        assert len(builds) == 1
        total = self._carbon(cv)
        run = solver._runs[cv]
        self._failing(run.stepper)
        run.t_end = run.stepper.t + 1.0  # the next macro step has to step
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            _advance(cv, solver)
        assert len(builds) == 2
        assert solver._runs[cv].stepper.t == DT_H
        assert self._carbon(cv) == pytest.approx(total, rel=1e-6)

    def test_a_second_failure_restores_the_state_and_warns(self, monkeypatch, steppers):
        cv, solver = _make_cv(), _solver()
        _advance(cv, solver, n=3)
        total = self._carbon(cv)
        made = solvers_mod._make_stepper

        def failing(*args, **kwargs):
            stepper = made(*args, **kwargs)
            self._failing(stepper)
            return stepper

        monkeypatch.setattr(solvers_mod, "_make_stepper", failing)
        run = solver._runs[cv]
        self._failing(run.stepper)
        run.t_end = run.stepper.t + 1.0
        with pytest.warns(RuntimeWarning, match="forced"):
            _advance(cv, solver)
        assert self._carbon(cv) == pytest.approx(total, rel=1e-12)
        assert cv not in solver._runs

    def test_a_new_run_that_fails_does_not_retry(self, monkeypatch, builds):
        cv, solver = _make_cv(), _solver()
        made = solvers_mod._make_stepper

        def failing(*args, **kwargs):
            stepper = made(*args, **kwargs)
            self._failing(stepper)
            return stepper

        monkeypatch.setattr(solvers_mod, "_make_stepper", failing)
        with pytest.warns(RuntimeWarning, match="forced"):
            _advance(cv, solver)
        assert len(builds) == 1


class TestAccuracyMonitor:
    class _Spy:
        def __init__(self):
            self.calls = []

        def check_scipy_rejections(self, n_accepted, n_attempted):
            self.calls.append((n_accepted, n_attempted))

        def check_pH_jump(self, pH):
            pass

    def test_it_sees_each_macro_steps_own_counts(self):
        cv, solver = _make_cv(), _solver(method="RK45")
        cv._accuracy_monitor = spy = self._Spy()
        _advance(cv, solver, n=12)
        assert spy.calls
        assert all(n_accepted > 0 and n_attempted > 0 for n_accepted, n_attempted in spy.calls)
        total = sum(n_attempted for _, n_attempted in spy.calls)
        assert total <= solver._runs[cv].stepper.nfev + 100

    def test_a_macro_step_with_no_integrator_step_is_not_judged(self):
        cv, solver = _make_cv(), _solver(method="RK45")
        cv._accuracy_monitor = spy = self._Spy()
        for _ in range(30):
            _advance(cv, solver, dt_h=1e-4)
        assert len(spy.calls) < 30


class TestParameterGeneration:
    def _changes(self, write):
        before = parameter_generation()
        write()
        return parameter_generation() - before

    def test_new_objects_do_not_count(self):
        assert self._changes(_Knob) == 0
        assert self._changes(lambda: LiquidPhase(n_mol={"X": 1.0}, V_L=1.0, T_K=300.0)) == 0

    def test_scalar_change_counts_and_same_value_does_not(self):
        knob = _Knob()
        assert self._changes(lambda: setattr(knob, "k", 2.0)) == 1
        assert self._changes(lambda: setattr(knob, "k", 2.0)) == 0

    def test_dict_replacement_counts_and_same_value_does_not(self):
        knob = _Knob()
        assert self._changes(lambda: setattr(knob, "table", {"O2": 2.0})) == 1
        assert self._changes(lambda: setattr(knob, "table", {"O2": 2.0})) == 0

    def test_dict_item_counts_and_same_value_does_not(self):
        knob = _Knob()
        write = lambda value: type(knob).__dict__["table"]._set_item_unchecked(knob, "O2", value)
        assert self._changes(lambda: write(3.0)) == 1
        assert self._changes(lambda: write(3.0)) == 0
        assert self._changes(lambda: write(3.0) or
                             type(knob).__dict__["table"]._set_item_unchecked(knob, "CO2", 1.0)) == 1

    def test_plain_dict_leaf_counts_and_same_value_does_not(self):
        """`_write_leaf` writes `link.kLa["O2"]` into a plain dict."""
        sim = Simulation(cvs={"main": _make_cv()}, label="t")
        table = {"O2": 1.0}
        assert self._changes(lambda: sim._write_leaf(table, "O2", 2.0, None)) == 1
        assert self._changes(lambda: sim._write_leaf(table, "O2", 2.0, None)) == 0

    def test_unchecked_setter_fallback_counts(self):
        class Target:
            def _set_rate_unchecked(self, value):
                self.rate = value

        sim = Simulation(cvs={"main": _make_cv()}, label="t")
        assert self._changes(lambda: sim._write_leaf(Target(), "rate", 2.0, None)) == 1


class TestParameterTriggers:
    def test_a_parameter_change_restarts(self, builds):
        cv, solver, knob = _make_cv(), _solver(), _Knob()
        _advance(cv, solver)
        knob.k = 5.0
        _advance(cv, solver)
        assert len(builds) == 2

    def test_rewriting_a_parameter_unchanged_does_not_restart(self, builds):
        cv, solver, knob = _make_cv(), _solver(), _Knob()
        _advance(cv, solver)
        knob.k = 1.0
        _advance(cv, solver)
        assert len(builds) == 1

    def test_scratch_phases_of_another_cv_do_not_restart(self, builds):
        """A restart in one CV builds phases; that must not restart the other."""
        solver = _solver()
        a, b = _make_cv(), _make_cv()
        _advance(a, solver)
        _advance(b, solver)
        for _ in range(3):
            b.phases["liquid"].n_mol["Na+"] += 1e-4  # b restarts every step
            _advance(b, solver)
            _advance(a, solver)
        assert len([s for s in builds]) == 5  # a once, b four times

    def _run(self, **kwargs):
        cv = _make_cv()
        sim = Simulation(cvs={"main": cv}, solver=_solver(), label="t", **kwargs)
        sim.run(tau_h=0.1, n_steps=5)

    @staticmethod
    def _profile(targets_at):
        class Profile:
            def apply(self, t_h, sim):
                return ProfileRecord(profile_label="p", t_h=t_h, targets=targets_at(t_h))
        return Profile()

    def test_profile_with_a_changing_report_restarts(self, builds):
        self._run(profiles=[self._profile(lambda t_h: {"plain.attribute": t_h})])
        assert len(builds) == 5

    def test_profile_with_a_constant_report_does_not_restart(self, builds):
        self._run(profiles=[self._profile(lambda t_h: {"plain.attribute": 1.0})])
        assert len(builds) == 1

    def test_controller_doing_nothing_does_not_restart(self, builds):
        class Idle:
            label = "idle"

            def compute(self, state, dt_h):
                return ControlAction(controller_label="idle", target_cv_key="main",
                                     t_h=0.0, dt_h=dt_h)
        self._run(controllers=[Idle()])
        assert len(builds) == 1


class TestAgainstExactAnswers:
    @pytest.mark.parametrize("method", ["BDF", "Radau", "LSODA", "RK45", "DOP853"])
    def test_time_dependent_rate_law(self, method):
        """dA/dt = -k t A has a closed form; Simulation's t_h leads by one step."""
        k, tau_h, n_steps = 3.0, 1.0, 10
        dt_h = tau_h / n_steps
        cv = _make_cv({"Na+": 0.01})
        species_a = Species(id="A", atoms={"C": 1}, charge=0)
        species_b = Species(id="B", atoms={"C": 1}, charge=0)
        cv.phases["liquid"].n_mol.update({"A": 1.0, "B": 0.0})
        system = cv.reaction_system
        system._kinetic_reactions[0] = KineticReaction(
            stoichiometry=[_entry(species_a, -1), _entry(species_b, +1)],
            rate_fn=lambda env: k * env.t_h * env.S("A") * env.V_L,
            label="timed",
        )
        solver = _solver(method=method, rtol=1e-9, atol=1e-11)
        Simulation(cvs={"main": cv}, solver=solver, label="t").run(tau_h=tau_h, n_steps=n_steps)
        exact = np.exp(-k * ((tau_h + dt_h) ** 2 - dt_h ** 2) / 2.0)
        assert cv.phases["liquid"].n_mol["A"] == pytest.approx(exact, rel=1e-6)

    def test_linked_cvs_agree_with_restarting(self):
        """Links change the state every macro step, so every step restarts."""
        def run(**kwargs):
            cvs = {"a": _make_cv(), "b": _make_cv({"AcOH": 0.01})}
            links = [AdvectiveLink(
                _source_cv_key=x, _source_phase_key="liquid",
                _sink_cv_key=y, _sink_phase_key="liquid", Q_L_per_h=0.5,
            ) for x, y in (("a", "b"), ("b", "a"))]
            sim = Simulation(cvs=cvs, links=links, solver=_solver(**kwargs), label="t")
            sim.run(tau_h=0.5, n_steps=10)
            return cvs

        persistent, restarting = run(), run(persistent=False)
        for key in ("a", "b"):
            for sp in ("P", "AcOH", "Ac-"):
                assert persistent[key].phases["liquid"].n_mol[sp] == pytest.approx(
                    restarting[key].phases["liquid"].n_mol[sp], rel=1e-5,
                )

    def test_resuming_with_a_pickled_solver_agrees_with_an_unbroken_run(self):
        """A checkpoint pickles the solver, which restarts once when loaded."""
        unbroken = _make_cv()
        Simulation(cvs={"main": unbroken}, solver=_solver(), label="t").run(tau_h=1.0, n_steps=20)

        cv = _make_cv()
        sim = Simulation(cvs={"main": cv}, solver=_solver(), label="t")
        sim.run(tau_h=0.5, n_steps=10)
        sim.solver = pickle.loads(pickle.dumps(sim.solver))
        sim.run(tau_h=0.5, n_steps=10, start_t_h=0.5)
        assert cv.phases["liquid"].n_mol["P"] == pytest.approx(
            unbroken.phases["liquid"].n_mol["P"], rel=1e-5,
        )


class TestScope:
    def test_each_cv_has_its_own_integrator(self, builds):
        solver = _solver()
        a, b = _make_cv(), _make_cv({"AcOH": 0.02})
        for _ in range(4):
            _advance(a, solver)
            _advance(b, solver)
        assert len(builds) == 2

        alone = _make_cv({"AcOH": 0.02})
        _advance(alone, _solver(), n=4)
        assert b.phases["liquid"].n_mol["P"] == pytest.approx(
            alone.phases["liquid"].n_mol["P"], rel=1e-12,
        )

    def test_pickled_solver_restarts_and_runs(self, builds):
        cv, solver = _make_cv(), _solver()
        _advance(cv, solver, n=2)
        loaded = pickle.loads(pickle.dumps(solver))
        assert loaded.persistent is True
        _advance(cv, loaded)
        assert len(builds) == 2
