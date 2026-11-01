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
from PyOMES.core import ControlVolume, LiquidPhase
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
    kwargs.setdefault("restart_each_step", False)
    return SimultaneousAdaptiveSolver(**kwargs)


@pytest.fixture
def builds(monkeypatch):
    """Records each new integrator."""
    made = []
    real = solvers_mod._make_stepper

    def counting(*args, **kwargs):
        stepper = real(*args, **kwargs)
        made.append(stepper)
        return stepper

    monkeypatch.setattr(solvers_mod, "_make_stepper", counting)
    return made


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def _advance(cv, solver, n=1, **kwargs):
    for _ in range(n):
        solver.solve_step(cv, dt_h=kwargs.get("dt_h", DT_H), t_h=0.0,
                          external_source_terms=kwargs.get("sources"))


class TestKeepsOneIntegrator:
    def test_nothing_changed_means_one_integrator(self, builds):
        cv, solver = _make_cv(), _solver()
        _advance(cv, solver, n=5)
        assert len(builds) == 1

    def test_restart_each_step_builds_no_stepper(self, builds):
        cv, solver = _make_cv(), _solver(restart_each_step=True)
        _advance(cv, solver, n=3)
        assert builds == []

    def test_freeze_speciation_restarts_each_step(self, builds):
        cv, solver = _make_cv(), _solver(freeze_speciation=True)
        _advance(cv, solver, n=3)
        assert builds == []

    @pytest.mark.parametrize("method", ["BDF", "Radau", "LSODA", "RK45", "RK23", "DOP853"])
    def test_agrees_with_restarting_within_tolerance(self, method):
        a, b = _make_cv(), _make_cv()
        _advance(a, _solver(method=method), n=20)
        _advance(b, _solver(method=method, restart_each_step=True), n=20)
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
            solver.solve_step(cv, dt_h=1e-4, t_h=0.0)
            values.append(cv.phases["liquid"].n_mol["P"])
        assert all(later > earlier for earlier, later in zip(values, values[1:]))


class TestTotalsAreKept:
    """State entries of resolved species are only parts of a total, so a flux can
    take one negative. Flooring it would create material."""

    @pytest.mark.parametrize("method", ["BDF", "Radau", "RK45"])
    def test_rate_law_on_a_resolved_species_conserves_carbon(self, method):
        cv, solver = _make_cv(), _solver(method=method)
        liquid = cv.phases["liquid"].n_mol
        solver.solve_step(cv, dt_h=DT_H, t_h=0.0)
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
        solver.solve_step(cv, dt_h=DT_H, t_h=0.0)
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
        assert loaded.restart_each_step is False
        _advance(cv, loaded)
        assert len(builds) == 2

    def test_failure_restores_the_state_and_warns(self, builds, monkeypatch):
        cv, solver = _make_cv(), _solver()
        _advance(cv, solver)
        stepper = builds[0]

        def fail():
            stepper.status = "failed"
            return "forced"

        monkeypatch.setattr(stepper, "step", fail)
        solver_run = solver._runs[cv]
        solver_run.t_end = stepper.t + 1.0  # the next macro step has to step
        total = sum(cv.phases["liquid"].n_mol[sp] for sp in ("AcOH", "Ac-", "P"))
        with pytest.warns(RuntimeWarning, match="forced"):
            _advance(cv, solver)
        after = sum(cv.phases["liquid"].n_mol[sp] for sp in ("AcOH", "Ac-", "P"))
        assert after == pytest.approx(total, rel=1e-12)
        assert cv not in solver._runs
