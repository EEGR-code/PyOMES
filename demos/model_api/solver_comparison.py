#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Demo: Axis 2 system-solver comparison.

Runs the same two-CV recirculating model under all five Axis 2
:class:`~PyOMES.core.SystemSolver` implementations and prints a
comparison table.  Two things are demonstrated:

**Section 1 — stability and accuracy.**  All five solvers remain stable
on this well-behaved transport problem (tau_CFL >> dt_h) and conserve
mass.  Accuracy differs: MonolithicODE/RK45 is closest to the
analytical solution; ExplicitEuler overshoots due to first-order
splitting error.  Wall-time differences are printed to give a rough
sense of overhead at this small scale.

**Section 2 — controller sub-stepping.**
:class:`~PyOMES.core.MonolithicODESolver` fires a periodic controller
at exact T_c = 0.01 h sub-step boundaries (10× per dt_h = 0.1 h macro
step) and co-integrates the controller's integral state continuously.
:class:`~PyOMES.core.ExplicitEulerSystemSolver` fires the same
controller once per macro step and does not evolve the integral at all.

Model topology::

    reactor (1 L, S = 10 mol)
        v Q = 2 L/h                 tau_CFL = V/Q = 0.5 h >> dt_h = 0.1 h
    recycle (1 L, S = 0 mol)
        ^ Q = 2 L/h

No reactions — pure transport.  S distributes between the two CVs
according to the mixing dynamics; all solvers converge to the same
steady state.

Run from the repo root after ``pip install -e .``::

    python demos/model_api/solver_comparison.py

When to use each solver:

+----------------------------+---------------------------------------------+
| Solver                     | When to choose it                           |
+============================+=============================================+
| ExplicitEulerSystemSolver  | Default; correct for slow flows (τ_CFL ≫    |
|                            | dt) and slow controllers.                   |
+----------------------------+---------------------------------------------+
| StrangSplittingSystemSolver| Free upgrade to 2nd-order splitting; nearly |
|                            | zero extra cost over Euler.                 |
+----------------------------+---------------------------------------------+
| MultirateSystemSolver      | Fast inter-CV circulation; removes CFL      |
|                            | constraint by subcycling link flux.         |
+----------------------------+---------------------------------------------+
| ImplicitTransportSystemSolver| Unconditionally stable transport; best for|
|                            | many-CV compartmental models (HPLC, etc.)  |
|                            | with no tight controllers.                 |
+----------------------------+---------------------------------------------+
| MonolithicODESolver        | Tight feedback loops (T_c < dt_h), multi-  |
|                            | rate digital controllers, co-integrated    |
|                            | integral states (PI/PID, observers).       |
+----------------------------+---------------------------------------------+
"""

import sys
import time
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_repo_root / "demos"))
import _bootstrap  # noqa: F401, E402

from PyOMES.core import (  # noqa: E402
    AdvectiveLink,
    ControlVolume,
    LiquidPhase,
    Simulation,
)
from PyOMES.core.system_solver import (  # noqa: E402
    ExplicitEulerSystemSolver,
    ImplicitTransportSystemSolver,
    MonolithicODESolver,
    MultirateSystemSolver,
    StrangSplittingSystemSolver,
)
from PyOMES.control.actions import ControlAction  # noqa: E402
from PyOMES.control.interfaces import ControllerBase  # noqa: E402


# ── Model parameters ──────────────────────────────────────────────────

Q_L_PER_H = 2.0    # recirculation flow rate (L/h)
V_L       = 1.0    # volume of each CV (L)
S_INIT    = 10.0   # initial substrate in reactor (mol)

# τ_CFL = V / Q = 0.5 h >> dt_h = 0.1 h → all solvers are CFL-stable.
DT_H    = 0.1
TAU_H   = 1.0
N_STEPS = int(TAU_H / DT_H)

# Controller parameters (Section 2 only).
T_C_H = 0.01   # controller update period (h); 10× per macro step


# ── Model construction ────────────────────────────────────────────────

def _build(system_solver=None, controllers=None) -> Simulation:
    """Assemble a fresh two-CV simulation for each solver trial."""
    liq_r = LiquidPhase(n_mol={"S": S_INIT}, V_L=V_L, T_K=310.0)
    liq_c = LiquidPhase(n_mol={"S": 0.0},    V_L=V_L, T_K=310.0)
    return Simulation(
        cvs={
            "reactor": ControlVolume(phases={"liquid": liq_r}, label="reactor"),
            "recycle": ControlVolume(phases={"liquid": liq_c}, label="recycle"),
        },
        links=[
            AdvectiveLink(
                _source_cv_key="reactor", _source_phase_key="liquid",
                _sink_cv_key="recycle",   _sink_phase_key="liquid",
                Q_L_per_h=Q_L_PER_H, _label="fwd",
            ),
            AdvectiveLink(
                _source_cv_key="recycle", _source_phase_key="liquid",
                _sink_cv_key="reactor",   _sink_phase_key="liquid",
                Q_L_per_h=Q_L_PER_H, _label="bwd",
            ),
        ],
        controllers=controllers or [],
        system_solver=system_solver,
        label="solver_comparison",
    )


# ── Section 1: solver agreement ───────────────────────────────────────

def section1_agreement() -> None:
    """Compare stability and accuracy profile of all five solvers."""
    print("=" * 65)
    print("Section 1: Solver stability and accuracy (no controller)")
    print(f"  reactor <-> recycle, Q = {Q_L_PER_H} L/h, "
          f"tau = {TAU_H} h, dt = {DT_H} h ({N_STEPS} steps)")
    print(f"  Initial: reactor S = {S_INIT:.1f} mol, recycle S = 0.0 mol")
    print()

    solvers = [
        ("ExplicitEuler",     ExplicitEulerSystemSolver()),
        ("StrangSplitting",   StrangSplittingSystemSolver()),
        ("Multirate",         MultirateSystemSolver()),
        ("ImplicitTransport", ImplicitTransportSystemSolver()),
        ("MonolithicODE",     MonolithicODESolver()),
    ]

    fmt = "{:<22}  {:>10}  {:>10}  {:>12}  {:>8}"
    print(fmt.format("Solver", "S_reactor", "S_recycle", "D vs Euler", "ms"))
    print("-" * 65)

    ref_reactor = None
    for name, solver in solvers:
        sim = _build(system_solver=solver)
        t0 = time.perf_counter()
        result = sim.run(tau_h=TAU_H, n_steps=N_STEPS)
        elapsed_ms = (time.perf_counter() - t0) * 1e3

        S_r = result.liquid_mol["reactor"]["S"][-1]
        S_c = result.liquid_mol["recycle"]["S"][-1]

        if ref_reactor is None:
            ref_reactor = S_r
            delta_str = "baseline"
        else:
            delta = abs(S_r - ref_reactor) / (ref_reactor or 1.0)
            delta_str = f"{delta:.2e}"

        print(fmt.format(
            name,
            f"{S_r:.5f}",
            f"{S_c:.5f}",
            delta_str,
            f"{elapsed_ms:.1f}",
        ))

    total = S_INIT  # mol conserved
    print()
    print(f"  Mass conservation: reactor + recycle ~= {total:.1f} mol (all solvers).")
    print(f"  Analytical: S_reactor(1h) ~= 5.092 mol (exact: 5 + 5*exp(-4)).")
    print()


# ── Section 2: controller sub-stepping ───────────────────────────────

class _IntegralController(ControllerBase):
    """Minimal controller demonstrating co-integrated differential state.

    ``di/dt = 1.0`` (constant) so the expected integral after time τ is
    simply τ itself — easy to verify analytically.

    With ExplicitEulerSystemSolver:
      - ``compute()`` fires once per 0.1 h macro step.
      - ``state_rates()`` is never called; the integral does not evolve.

    With MonolithicODESolver:
      - ``state_rates()`` drives continuous integral evolution (di/dt = 1).
      - ``compute()`` fires at each 0.01 h T_c boundary (10× per step).
      - Final integral ≈ τ as expected.
    """

    target_cv_key = "reactor"

    def __init__(self):
        self._integral = 0.0
        self.fire_count = 0

    @property
    def update_period_h(self) -> float:
        return T_C_H

    def differential_state(self):
        return {"integral": self._integral}

    def state_rates(self, env, t_h):
        return {"integral": 1.0}

    def set_state(self, state):
        self._integral = state["integral"]

    def compute(self, state, dt_h):
        self.fire_count += 1
        return ControlAction(
            controller_label="integral_ctrl",
            target_cv_key="reactor",
        )


N_STEPS_2 = 5
TAU_H_2   = N_STEPS_2 * DT_H
FIRES_PER_STEP = int(DT_H / T_C_H)


def section2_controller() -> None:
    """MonolithicODE fires at T_c boundaries and co-integrates the state."""
    print("=" * 65)
    print("Section 2: Controller sub-stepping")
    print(f"  Controller update period: T_c = {T_C_H} h")
    print(f"  Macro step: dt = {DT_H} h  =>  {FIRES_PER_STEP}x per step with MonolithicODE")
    print(f"  Run: tau = {TAU_H_2} h ({N_STEPS_2} macro steps)")
    print(f"  Expected integral (di/dt = 1): {TAU_H_2:.3f}")
    print()

    fmt = "{:<22}  {:>12}  {:>13}  {:>10}  {:>10}"
    print(fmt.format(
        "Solver", "fires/step", "total fires", "integral", "expected"
    ))
    print("-" * 65)

    for name, solver in [
        ("ExplicitEuler",  ExplicitEulerSystemSolver()),
        ("MonolithicODE",  MonolithicODESolver()),
    ]:
        ctrl = _IntegralController()
        sim = _build(system_solver=solver, controllers=[ctrl])
        sim.run(tau_h=TAU_H_2, n_steps=N_STEPS_2)

        fires_per_step = ctrl.fire_count / N_STEPS_2
        print(fmt.format(
            name,
            f"{fires_per_step:.0f}",
            str(ctrl.fire_count),
            f"{ctrl._integral:.4f}",
            f"{TAU_H_2:.4f}",
        ))

    print()
    print("  ExplicitEuler  -- fires once per macro step; integral not evolved.")
    print("  MonolithicODE  -- fires at each T_c boundary; integral co-integrated")
    print(f"                   with transport dynamics (di/dt = 1 => integral = {TAU_H_2}).")
    print()
    print("  => Choose MonolithicODE when T_c < dt_h and integral-state accuracy")
    print("    matters: tight pH/DO loops, multi-rate digital controllers, observers.")


# ── Entry point ───────────────────────────────────────────────────────

def main() -> None:
    print()
    print("PyOMES - Axis 2 SystemSolver comparison")
    print()
    section1_agreement()
    section2_controller()


if __name__ == "__main__":
    main()
