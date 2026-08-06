# -*- coding: utf-8 -*-
"""Tests for Phase 5 of the speciation protocol hierarchy:
SimultaneousAdaptiveSolver use_engine_jacobian flag.

Coverage
--------
1. Constructor: flag stored, default False, ValueError for freeze+jac incompatibility.
2. Repr includes use_engine_jacobian.
3. UserWarning when GrayBox engine present and flag is False (BDF/Radau/LSODA).
4. No warning when flag is True, no speciation, or explicit method.
5. Step succeeds with use_engine_jacobian=True and produces same result as without.
6. jac_callable returns a finite (n, n) matrix.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Helpers shared by multiple test classes
# ---------------------------------------------------------------------------

def _make_reactions():
    """Minimal carbonate + ammonia equilibrium reactions."""
    from PyOMES.chemistry.common_species import (
        CO2, CO3_2minus, H2O, H_plus, HCO3_minus, NH3, NH4_plus, OH_minus,
    )
    from PyOMES.reactions.equilibrium import EquilibriumReaction
    from PyOMES.reactions.stoichiometry import StoichiometryEntry

    def _e(species, coeff):
        return StoichiometryEntry(species=species, phase="liquid", coefficient=coeff)

    water = EquilibriumReaction(
        stoichiometry=[_e(H2O, -1.0), _e(H_plus, +1.0), _e(OH_minus, +1.0)],
        log_K=-14.0,
        balance_elements=("H", "O"),
        label="water",
    )
    co2_first = EquilibriumReaction(
        stoichiometry=[
            _e(CO2, -1.0), _e(H2O, -1.0), _e(HCO3_minus, +1.0), _e(H_plus, +1.0),
        ],
        log_K=-6.35,
        total_id="CO2",
        balance_elements=("C", "H", "O"),
        label="co2_first",
    )
    co2_second = EquilibriumReaction(
        stoichiometry=[_e(HCO3_minus, -1.0), _e(CO3_2minus, +1.0), _e(H_plus, +1.0)],
        log_K=-10.33,
        total_id="CO2",
        balance_elements=("C", "H", "O"),
        label="co2_second",
    )
    nh4 = EquilibriumReaction(
        stoichiometry=[_e(NH4_plus, -1.0), _e(NH3, +1.0), _e(H_plus, +1.0)],
        log_K=-9.25,
        total_id="NH3",
        balance_elements=("N", "H"),
        label="nh4",
    )
    return [water, co2_first, co2_second, nh4]


def _make_cv_with_nr_engine(retain_jacobian: bool = False):
    """CV with NRChemicalEquilibriumEngine, carbonate + ammonia chemistry."""
    from PyOMES.core import ControlVolume, GasPhase, LiquidPhase
    from PyOMES.reactions import ReactionSystem
    from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

    reactions = _make_reactions()
    gas = GasPhase(n_mol={"CO2": 0.01, "N2": 0.5}, V_L=0.4, T_K=308.15)
    liquid = LiquidPhase(n_mol={"CO2": 0.005, "NH3": 0.04}, V_L=1.6, T_K=308.15)

    rxn_system = ReactionSystem(reactions, solver="newton_raphson")
    if retain_jacobian:
        nr_engine = NRChemicalEquilibriumEngine.from_reactions(reactions, retain_jacobian=True)
        rxn_system.attach_engine(nr_engine)

    return ControlVolume(
        phases={"gas": gas, "liquid": liquid},
        reaction_system=rxn_system,
    )


def _make_cv_with_blackbox_engine():
    """CV with the default BisectionChemicalEquilibriumEngine (charge-balance, not GrayBox)."""
    from PyOMES.core import ControlVolume, GasPhase, LiquidPhase
    from PyOMES.reactions import ReactionSystem

    reactions = _make_reactions()
    gas = GasPhase(n_mol={"CO2": 0.01, "N2": 0.5}, V_L=0.4, T_K=308.15)
    liquid = LiquidPhase(n_mol={"CO2": 0.005, "NH3": 0.04}, V_L=1.6, T_K=308.15)
    rxn_system = ReactionSystem(reactions)  # default solver="charge_balance"
    return ControlVolume(
        phases={"gas": gas, "liquid": liquid},
        reaction_system=rxn_system,
    )


def _make_cv_no_speciation():
    """CV without a reaction system."""
    from PyOMES.core import ControlVolume, GasPhase, LiquidPhase

    gas = GasPhase(n_mol={"CO2": 0.01, "N2": 0.5}, V_L=0.4, T_K=308.15)
    liquid = LiquidPhase(n_mol={"CO2": 0.005}, V_L=1.6, T_K=308.15)
    return ControlVolume(phases={"gas": gas, "liquid": liquid})


# ---------------------------------------------------------------------------
# 1. Constructor
# ---------------------------------------------------------------------------

class TestConstructor:
    def test_default_flag_is_false(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        solver = SimultaneousAdaptiveSolver()
        assert solver.use_engine_jacobian is False

    def test_flag_stored_when_true(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        solver = SimultaneousAdaptiveSolver(use_engine_jacobian=True)
        assert solver.use_engine_jacobian is True

    def test_freeze_and_jac_incompatible(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        with pytest.raises(ValueError, match="incompatible"):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                SimultaneousAdaptiveSolver(freeze_speciation=True, use_engine_jacobian=True)

    def test_freeze_alone_no_raise(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            solver = SimultaneousAdaptiveSolver(freeze_speciation=True)
        assert solver.freeze_speciation is True

    def test_use_engine_jacobian_alone_no_raise(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        solver = SimultaneousAdaptiveSolver(use_engine_jacobian=True)
        assert solver.use_engine_jacobian is True

    def test_other_params_still_work(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        solver = SimultaneousAdaptiveSolver(method="BDF", rtol=1e-4, atol=1e-7,
                                max_step=0.5, use_engine_jacobian=True)
        assert solver.method == "BDF"
        assert solver.rtol == pytest.approx(1e-4)
        assert solver.atol == pytest.approx(1e-7)
        assert solver.max_step == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# 2. Repr
# ---------------------------------------------------------------------------

class TestRepr:
    def test_repr_includes_flag_false(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        r = repr(SimultaneousAdaptiveSolver())
        assert "use_engine_jacobian=False" in r

    def test_repr_includes_flag_true(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        r = repr(SimultaneousAdaptiveSolver(use_engine_jacobian=True))
        assert "use_engine_jacobian=True" in r

    def test_repr_includes_method(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        r = repr(SimultaneousAdaptiveSolver(method="Radau"))
        assert "Radau" in r


# ---------------------------------------------------------------------------
# 3. UserWarning when GrayBox + flag=False + implicit method
# ---------------------------------------------------------------------------

class TestGrayBoxWarning:
    def _run_step(self, cv, solver, catch_exc=True):
        """Call solve_step, optionally suppressing all non-UserWarning exceptions."""
        try:
            solver.solve_step(cv, dt_h=0.01, t_h=0.0)
        except Exception:
            if not catch_exc:
                raise

    def test_warning_bdf_graybox_flag_false(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        cv = _make_cv_with_nr_engine()
        solver = SimultaneousAdaptiveSolver(method="BDF")
        with pytest.warns(UserWarning, match="use_engine_jacobian"):
            self._run_step(cv, solver)

    def test_warning_radau_graybox_flag_false(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        cv = _make_cv_with_nr_engine()
        solver = SimultaneousAdaptiveSolver(method="Radau")
        with pytest.warns(UserWarning, match="GrayBoxEngineProtocol"):
            self._run_step(cv, solver)

    def test_warning_lsoda_graybox_flag_false(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        cv = _make_cv_with_nr_engine()
        solver = SimultaneousAdaptiveSolver(method="LSODA")
        with pytest.warns(UserWarning, match="use_engine_jacobian"):
            self._run_step(cv, solver)

    def test_warning_message_exact_text(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        cv = _make_cv_with_nr_engine()
        solver = SimultaneousAdaptiveSolver(method="BDF")
        with pytest.warns(UserWarning) as record:
            self._run_step(cv, solver)
        messages = [str(w.message) for w in record]
        assert any("GrayBoxEngineProtocol" in m and "use_engine_jacobian=False" in m
                   for m in messages)

    def test_no_warning_when_flag_true(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        cv = _make_cv_with_nr_engine(retain_jacobian=True)
        solver = SimultaneousAdaptiveSolver(method="BDF", use_engine_jacobian=True)
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            # Should not raise UserWarning about use_engine_jacobian
            try:
                solver.solve_step(cv, dt_h=0.01, t_h=0.0)
            except UserWarning as e:
                if "use_engine_jacobian" in str(e):
                    pytest.fail(f"Unexpected UserWarning: {e}")

    def test_no_warning_for_blackbox_engine(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        cv = _make_cv_with_blackbox_engine()
        solver = SimultaneousAdaptiveSolver(method="BDF")
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            try:
                solver.solve_step(cv, dt_h=0.01, t_h=0.0)
            except UserWarning as e:
                if "use_engine_jacobian" in str(e):
                    pytest.fail(f"Unexpected jac UserWarning: {e}")

    def test_no_warning_no_speciation(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        cv = _make_cv_no_speciation()
        solver = SimultaneousAdaptiveSolver(method="BDF")
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            try:
                solver.solve_step(cv, dt_h=0.01, t_h=0.0)
            except UserWarning as e:
                if "use_engine_jacobian" in str(e):
                    pytest.fail(f"Unexpected jac UserWarning: {e}")

    def test_no_warning_for_explicit_method(self):
        """DOP853 doesn't support jac= so no advisory warning."""
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        cv = _make_cv_with_nr_engine()
        solver = SimultaneousAdaptiveSolver(method="DOP853")
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            try:
                solver.solve_step(cv, dt_h=0.01, t_h=0.0)
            except UserWarning as e:
                if "use_engine_jacobian" in str(e):
                    pytest.fail(f"Unexpected jac UserWarning for DOP853: {e}")

    def test_no_warning_for_rk45(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        cv = _make_cv_with_nr_engine()
        solver = SimultaneousAdaptiveSolver(method="RK45")
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            try:
                solver.solve_step(cv, dt_h=0.01, t_h=0.0)
            except UserWarning as e:
                if "use_engine_jacobian" in str(e):
                    pytest.fail(f"Unexpected jac UserWarning for RK45: {e}")


# ---------------------------------------------------------------------------
# 4. Analytical Jacobian functional tests
# ---------------------------------------------------------------------------

class TestAnalyticalJacobian:
    def test_step_succeeds_with_bdf_and_engine_jacobian(self):
        """BDF + use_engine_jacobian=True should complete without error."""
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        cv = _make_cv_with_nr_engine(retain_jacobian=True)
        solver = SimultaneousAdaptiveSolver(method="BDF", use_engine_jacobian=True)
        # No exception should be raised
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            solver.solve_step(cv, dt_h=0.01, t_h=0.0)

    def test_step_succeeds_with_radau_and_engine_jacobian(self):
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        cv = _make_cv_with_nr_engine(retain_jacobian=True)
        solver = SimultaneousAdaptiveSolver(method="Radau", use_engine_jacobian=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            solver.solve_step(cv, dt_h=0.01, t_h=0.0)

    def test_step_result_close_to_without_jacobian(self):
        """BDF results with and without engine Jacobian should be numerically close."""
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver

        cv1 = _make_cv_with_nr_engine(retain_jacobian=True)
        cv2 = _make_cv_with_nr_engine(retain_jacobian=True)

        solver_no_jac = SimultaneousAdaptiveSolver(method="BDF", rtol=1e-8, atol=1e-11)
        solver_jac = SimultaneousAdaptiveSolver(method="BDF", rtol=1e-8, atol=1e-11,
                                    use_engine_jacobian=True)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            solver_no_jac.solve_step(cv1, dt_h=0.05, t_h=0.0)
            solver_jac.solve_step(cv2, dt_h=0.05, t_h=0.0)

        for phase_key in ("gas", "liquid"):
            for sp in cv1.phases[phase_key].n_mol:
                if sp in cv2.phases[phase_key].n_mol:
                    v1 = cv1.phases[phase_key].n_mol[sp]
                    v2 = cv2.phases[phase_key].n_mol[sp]
                    assert abs(v1 - v2) <= 1e-6 * max(abs(v1), abs(v2), 1e-20), (
                        f"{phase_key}.{sp}: without_jac={v1}, with_jac={v2}"
                    )

    def test_jac_callable_returns_finite_matrix(self):
        """jac_callable returns a finite (n_ode, n_ode) numpy array."""
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver
        import numpy as np

        cv = _make_cv_with_nr_engine(retain_jacobian=True)
        solver = SimultaneousAdaptiveSolver(method="BDF", use_engine_jacobian=True)

        # Capture the jac callable by running one step and introspecting via monkeypatch
        captured = {}

        from scipy.integrate import solve_ivp as real_solve_ivp

        def patched_solve_ivp(f, t_span, y0, **kwargs):
            captured["jac"] = kwargs.get("jac")
            captured["y0"] = y0
            captured["f"] = f
            return real_solve_ivp(f, t_span, y0, **kwargs)

        import PyOMES.core.solvers as solvers_mod
        orig = getattr(solvers_mod, "_real_solve_ivp", None)

        # Patch scipy.integrate.solve_ivp inside the solver's import namespace
        import scipy.integrate
        orig_ivp = scipy.integrate.solve_ivp
        scipy.integrate.solve_ivp = patched_solve_ivp

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                solver.solve_step(cv, dt_h=0.01, t_h=0.0)
        finally:
            scipy.integrate.solve_ivp = orig_ivp

        assert "jac" in captured, "jac= kwarg not passed to solve_ivp"
        jac_fn = captured["jac"]
        assert jac_fn is not None, "jac= was None for use_engine_jacobian=True"

        y0 = captured["y0"]
        J = jac_fn(0.0, y0)

        assert isinstance(J, np.ndarray), f"jac_callable returned {type(J)}"
        assert J.shape == (len(y0), len(y0)), f"wrong shape: {J.shape}"
        assert np.all(np.isfinite(J)), "Jacobian contains non-finite values"

    def test_jac_not_passed_for_explicit_method(self):
        """DOP853 does not use jac= so it should not be passed."""
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver

        cv = _make_cv_with_nr_engine(retain_jacobian=True)
        solver = SimultaneousAdaptiveSolver(method="DOP853", use_engine_jacobian=True)

        captured = {}
        import scipy.integrate
        orig_ivp = scipy.integrate.solve_ivp

        def patched(f, t_span, y0, **kwargs):
            captured["jac"] = kwargs.get("jac")
            return orig_ivp(f, t_span, y0, **kwargs)

        scipy.integrate.solve_ivp = patched
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                solver.solve_step(cv, dt_h=0.01, t_h=0.0)
        finally:
            scipy.integrate.solve_ivp = orig_ivp

        # jac should be None because DOP853 is not in _JAC_METHODS
        assert captured.get("jac") is None

    def test_jac_fallback_when_retain_jacobian_false(self):
        """When engine Jacobian raises RuntimeError (retain_jacobian=False),
        jac_callable falls back to full FD and still returns a finite matrix."""
        from PyOMES.core.solvers import SimultaneousAdaptiveSolver

        cv = _make_cv_with_nr_engine(retain_jacobian=False)  # jac raises RuntimeError
        solver = SimultaneousAdaptiveSolver(method="BDF", use_engine_jacobian=True)

        captured = {}
        import scipy.integrate
        orig_ivp = scipy.integrate.solve_ivp

        def patched(f, t_span, y0, **kwargs):
            captured["jac"] = kwargs.get("jac")
            captured["y0"] = y0
            return orig_ivp(f, t_span, y0, **kwargs)

        scipy.integrate.solve_ivp = patched
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                solver.solve_step(cv, dt_h=0.01, t_h=0.0)
        finally:
            scipy.integrate.solve_ivp = orig_ivp

        assert captured.get("jac") is not None
        y0 = captured["y0"]
        J = captured["jac"](0.0, y0)
        import numpy as np
        assert np.all(np.isfinite(J))
