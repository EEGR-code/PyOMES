# -*- coding: utf-8 -*-
"""Tests for PyOMES/core/clamping.py (STEP_SOLVER_INTERFACE_REFINEMENT.md
item 7 -- shared, swappable clamp_fn module). Checkpoint 5 of
STEP_SOLVER_REFINEMENT_CHECKLIST.md.

Covers:
- proportional_clamp / floor_clamp as standalone pure functions,
  including "arbitrary subset of species" (not whole-dict ownership)
- clamp_fn wired into SequentialAdvanceSolver and SimultaneousEulerSolver
  as a swappable constructor parameter
- clamp_fn=None genuinely disables clamping (negative n_mol reachable)
- a bespoke per-species composite clamp_fn (the design note's own
  bsm2_style_clamp example) built out of proportional_clamp as an
  ingredient
"""

from __future__ import annotations

import pytest


# ═══════════════════════════════════════════════════════════════════════
#  proportional_clamp / floor_clamp as pure functions
# ═══════════════════════════════════════════════════════════════════════

class TestProportionalClamp:

    def test_no_op_when_nothing_would_go_negative(self):
        from PyOMES.core.clamping import proportional_clamp
        deltas = {"S": -1.0, "X": 2.0}
        current = {"S": 100.0, "X": 5.0}
        out = proportional_clamp(deltas, current, dt_h=0.1)
        assert out == deltas

    def test_scales_removal_to_exactly_zero(self):
        from PyOMES.core.clamping import proportional_clamp
        # current=1.0 mol, rate=-20 mol/h, dt=0.1h -> naive delta = -2.0 mol
        deltas = {"S": -20.0}
        current = {"S": 1.0}
        out = proportional_clamp(deltas, current, dt_h=0.1)
        n_after = current["S"] + out["S"] * 0.1
        assert n_after == pytest.approx(0.0, abs=1e-12)

    def test_operates_on_arbitrary_species_subset(self):
        """clamp_fn must not assume whole-dict ownership -- current_mol
        can carry species the deltas dict doesn't mention."""
        from PyOMES.core.clamping import proportional_clamp
        deltas = {"S": -20.0}
        current = {"S": 1.0, "X": 999.0, "Unrelated": 0.0}
        out = proportional_clamp(deltas, current, dt_h=0.1)
        assert set(out.keys()) == {"S"}

    def test_positive_rate_never_clamped(self):
        from PyOMES.core.clamping import proportional_clamp
        deltas = {"S": 50.0}
        current = {"S": 0.0}
        out = proportional_clamp(deltas, current, dt_h=1.0)
        assert out["S"] == 50.0


class TestFloorClamp:

    def test_floors_at_zero_by_default(self):
        from PyOMES.core.clamping import floor_clamp
        deltas = {"S": -20.0}
        current = {"S": 1.0}
        out = floor_clamp(deltas, current, dt_h=0.1)
        n_after = current["S"] + out["S"] * 0.1
        assert n_after == pytest.approx(0.0, abs=1e-12)

    def test_floors_at_eps_when_given(self):
        from PyOMES.core.clamping import floor_clamp
        deltas = {"S": -20.0}
        current = {"S": 1.0}
        out = floor_clamp(deltas, current, dt_h=0.1, eps=1e-9)
        n_after = current["S"] + out["S"] * 0.1
        assert n_after == pytest.approx(1e-9, abs=1e-12)


class TestFloorNonnegative:
    """floor_nonnegative(y) -- a distinct concern from clamp_fn: raw
    ODE state-array flooring for SimultaneousAdaptiveSolver (checkpoint
    7b), not a (deltas, current_mol, dt_h) triple."""

    def test_floors_negative_elements(self):
        import numpy as np
        from PyOMES.core.clamping import floor_nonnegative
        y = np.array([-5.0, 0.0, 3.0, -1e-12])
        out = floor_nonnegative(y)
        assert (out >= 0.0).all()
        assert out[2] == pytest.approx(3.0)

    def test_no_op_on_already_nonnegative(self):
        import numpy as np
        from PyOMES.core.clamping import floor_nonnegative
        y = np.array([0.0, 1.0, 2.5])
        out = floor_nonnegative(y)
        assert list(out) == list(y)


# ═══════════════════════════════════════════════════════════════════════
#  clamp_fn wired into the two discrete-step solvers
# ═══════════════════════════════════════════════════════════════════════

def _liquid_cv_with_fast_removal(k_removal=1000.0):
    """A liquid CV with one kinetic reaction that removes S at a rate
    aggressive enough to overdraw inventory in one dt_h -- the scenario
    clamp_fn is meant to guard against."""
    from PyOMES.core import ControlVolume, LiquidPhase
    from PyOMES.reactions import KineticReaction, ReactionSystem, StoichiometryEntry
    from PyOMES.chemistry import Species

    S = Species(id="S", atoms={"C": 1})
    Waste = Species(id="Waste", atoms={"C": 1})
    liq = LiquidPhase(n_mol={"S": 1.0, "Waste": 0.0}, V_L=1.0, T_K=300.0)
    rxn = KineticReaction(
        stoichiometry=[
            StoichiometryEntry(species=S, phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=Waste, phase="liquid", coefficient=+1.0),
        ],
        rate_fn=lambda env: k_removal,  # constant extensive rate, mol/h
        balance_elements=("C",),
        label="fast_removal",
    )
    rs = ReactionSystem(reactions=[rxn])
    return ControlVolume(phases={"liquid": liq}, reaction_system=rs, label="cv")


class TestClampFnOnSequentialAdvanceSolver:

    def test_default_clamp_prevents_negative_inventory(self):
        from PyOMES.core import SequentialAdvanceSolver
        cv = _liquid_cv_with_fast_removal()
        cv.advance(dt_h=0.1, solver=SequentialAdvanceSolver())
        assert cv.phases["liquid"].n_mol["S"] >= 0.0

    def test_clamp_fn_none_allows_negative_inventory(self):
        from PyOMES.core import SequentialAdvanceSolver
        cv = _liquid_cv_with_fast_removal()
        cv.advance(dt_h=0.1, solver=SequentialAdvanceSolver(clamp_fn=None))
        assert cv.phases["liquid"].n_mol["S"] < 0.0

    def test_default_solver_none_matches_explicit_default_clamp_fn(self):
        """solver=None (SequentialAdvanceSolver's own default clamp_fn)
        must be numerically identical to the pre-checkpoint-5 hardcoded
        floor -- proportional_clamp reduces to an exact floor for a
        single independently-clamped species."""
        cv_default = _liquid_cv_with_fast_removal()
        cv_explicit = _liquid_cv_with_fast_removal()
        from PyOMES.core import SequentialAdvanceSolver
        cv_default.advance(dt_h=0.1)
        cv_explicit.advance(dt_h=0.1, solver=SequentialAdvanceSolver())
        assert cv_default.phases["liquid"].n_mol["S"] == pytest.approx(
            cv_explicit.phases["liquid"].n_mol["S"]
        )


class TestClampFnOnSimultaneousEulerSolver:

    def test_default_clamp_prevents_negative_inventory(self):
        from PyOMES.core import SimultaneousEulerSolver
        cv = _liquid_cv_with_fast_removal()
        cv.advance(dt_h=0.1, solver=SimultaneousEulerSolver())
        assert cv.phases["liquid"].n_mol["S"] >= 0.0

    def test_clamp_fn_none_allows_negative_inventory(self):
        from PyOMES.core import SimultaneousEulerSolver
        cv = _liquid_cv_with_fast_removal()
        cv.advance(dt_h=0.1, solver=SimultaneousEulerSolver(clamp_fn=None))
        assert cv.phases["liquid"].n_mol["S"] < 0.0

    def test_floor_clamp_swappable_in(self):
        from PyOMES.core import SimultaneousEulerSolver
        from PyOMES.core.clamping import floor_clamp
        cv = _liquid_cv_with_fast_removal()
        cv.advance(dt_h=0.1, solver=SimultaneousEulerSolver(clamp_fn=floor_clamp))
        assert cv.phases["liquid"].n_mol["S"] >= 0.0


# ═══════════════════════════════════════════════════════════════════════
#  Bespoke composite clamp_fn (design note's bsm2_style_clamp example)
# ═══════════════════════════════════════════════════════════════════════

class TestCompositeClampFn:

    def test_bsm2_style_clamp_composite(self):
        """A clamp_fn built out of proportional_clamp as an ingredient,
        with per-species overrides -- the shared module's default must
        stay usable as a building block, not a monolithic
        take-it-or-leave-it function (per the source note's explicit
        requirement)."""
        from PyOMES.core.clamping import proportional_clamp

        NEVER_CLAMP = {"NeverClamped"}
        FLOOR_AT_EPSILON = {"FloorSpecies"}
        EPS = 1e-9

        def bsm2_style_clamp(deltas, current_mol, dt_h):
            result = dict(deltas)
            for sp in NEVER_CLAMP:
                if sp in deltas:
                    result[sp] = deltas[sp]
            for sp in FLOOR_AT_EPSILON:
                if sp in deltas:
                    n_current = current_mol.get(sp, 0.0)
                    delta_mol = deltas[sp] * dt_h
                    n_after = n_current + delta_mol
                    if n_after < EPS and delta_mol < 0.0:
                        max_removal = n_current - EPS
                        scale = max(0.0, min(1.0, max_removal / abs(delta_mol)))
                        result[sp] = deltas[sp] * scale
            remaining = {
                k: v for k, v in deltas.items()
                if k not in NEVER_CLAMP and k not in FLOOR_AT_EPSILON
            }
            result.update(proportional_clamp(remaining, current_mol, dt_h))
            return result

        deltas = {"NeverClamped": -100.0, "FloorSpecies": -100.0, "Normal": -100.0}
        current = {"NeverClamped": 1.0, "FloorSpecies": 1.0, "Normal": 1.0}
        out = bsm2_style_clamp(deltas, current, dt_h=1.0)

        # NeverClamped passes through unscaled, even though it would go negative.
        assert out["NeverClamped"] == -100.0
        # FloorSpecies lands at EPS, not exactly zero.
        n_after_floor = current["FloorSpecies"] + out["FloorSpecies"] * 1.0
        assert n_after_floor == pytest.approx(EPS, abs=1e-12)
        # Normal uses the shared proportional_clamp default -> lands at exactly zero.
        n_after_normal = current["Normal"] + out["Normal"] * 1.0
        assert n_after_normal == pytest.approx(0.0, abs=1e-12)


# ═══════════════════════════════════════════════════════════════════════
#  Diagnostics wired into the solvers (STEP_SOLVER_INTERFACE_REFINEMENT.md
#  item 8) -- checkpoint 6 of STEP_SOLVER_REFINEMENT_CHECKLIST.md
# ═══════════════════════════════════════════════════════════════════════

@pytest.fixture(autouse=True)
def _reset_warning_state():
    import PyOMES
    from PyOMES import WarningConfig
    from PyOMES.monitoring.accuracy import reset_accuracy_summary
    snapshot = PyOMES.config.warnings
    PyOMES.config.warnings = WarningConfig(throttle="always")
    reset_accuracy_summary()
    try:
        yield
    finally:
        PyOMES.config.warnings = snapshot
        reset_accuracy_summary()


class TestNegativeMoleDiagnostic:

    def test_clamp_fn_none_emits_negative_mole_warning(self):
        import warnings
        from PyOMES import AccuracyWarning
        from PyOMES.core import SequentialAdvanceSolver
        cv = _liquid_cv_with_fast_removal()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cv.advance(dt_h=0.1, solver=SequentialAdvanceSolver(clamp_fn=None))
        categories = [w.category for w in caught]
        assert any(issubclass(c, AccuracyWarning) for c in categories)

    def test_default_clamp_emits_no_negative_mole_warning(self):
        import warnings
        from PyOMES import AccuracyWarning
        from PyOMES.core import SequentialAdvanceSolver
        cv = _liquid_cv_with_fast_removal()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cv.advance(dt_h=0.1, solver=SequentialAdvanceSolver())
        neg_mole = [
            w for w in caught
            if issubclass(w.category, AccuracyWarning)
            and "went negative" in str(w.message)
        ]
        assert neg_mole == []


class TestClampInvokedDiagnostic:

    def test_default_clamp_emits_clamp_invoked_on_aggressive_kinetics(self):
        import warnings
        from PyOMES import AccuracyWarning
        from PyOMES.core import SimultaneousEulerSolver
        cv = _liquid_cv_with_fast_removal()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cv.advance(dt_h=0.1, solver=SimultaneousEulerSolver())
        clamp_invoked = [
            w for w in caught
            if issubclass(w.category, AccuracyWarning)
            and "clamp_fn scaled" in str(w.message)
        ]
        assert len(clamp_invoked) >= 1

    def test_mild_kinetics_emits_no_clamp_invoked(self):
        import warnings
        from PyOMES import AccuracyWarning
        from PyOMES.core import SimultaneousEulerSolver
        cv = _liquid_cv_with_fast_removal(k_removal=0.01)  # far from overdrawing
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cv.advance(dt_h=0.1, solver=SimultaneousEulerSolver())
        clamp_invoked = [
            w for w in caught
            if issubclass(w.category, AccuracyWarning)
            and "clamp_fn scaled" in str(w.message)
        ]
        assert clamp_invoked == []
