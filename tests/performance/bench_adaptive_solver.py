# -*- coding: utf-8 -*-
"""Speed and accuracy benchmark for ``SimultaneousAdaptiveSolver``.

Run by hand from the repo root; pytest does not collect it::

    python tests/performance/bench_adaptive_solver.py check
    python tests/performance/bench_adaptive_solver.py check --exact
    python tests/performance/bench_adaptive_solver.py check --cases iron_bdf bsm2_bdf
    python tests/performance/bench_adaptive_solver.py check --solver-kwargs "{\\"rtol\\": 1e-7}"
    python tests/performance/bench_adaptive_solver.py baseline
    python tests/performance/bench_adaptive_solver.py reference

``check`` compares a run against ``baselines/adaptive_solver.json``. It gates on
numbers that do not depend on the machine: engine solves and error against the stored
reference, or, with ``--exact``, a fingerprint of every recorded value. CPU time is
reported and never gated. ``baseline`` and ``reference`` rewrite their parts of that
file; use them only when a change to the stored numbers is intended.

The reference is Radau at ``REFERENCE_TOLERANCES``, on each case's own macro grid,
with one integrator across macro steps. A restarting reference is not used: restarts
do not conserve totals as well (total Fe drifted 1.6e-6 on iron and 9.8e-5 with a feed
and drain, against about 1e-7), and that bias would be charged to whatever is compared
with it. Persistent Radau, BDF and LSODA agreed to within ``ERROR_FLOOR`` on the
feed-and-drain and weak-acid cases, so smaller errors are treated as equal. Each timed
run builds its model first and times only ``Simulation.run``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import sys
import time
import warnings
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

REPO = Path(__file__).resolve().parents[2]
for _path in (REPO, REPO / "models", REPO / "tests" / "standalone",
              REPO / "tests" / "validation" / "speciation"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import numpy as np

from PyOMES.chemical_equilibrium.engines.nr.engine import NRChemicalEquilibriumEngine
from PyOMES.chemistry import Species
from PyOMES.control import PHController
from PyOMES.control.cv_profiles import TemperatureRamp
from PyOMES.core import ControlVolume, LiquidPhase, Simulation
from PyOMES.core.boundaries import LiquidDrain, LiquidFeed
from PyOMES.core.solvers import SimultaneousAdaptiveSolver
from PyOMES.databases.aqueous import H2O, H_plus, OH_minus
from PyOMES.reactions import (
    EquilibriumReaction, KineticReaction, ReactionSystem, StoichiometryEntry,
)

BASELINE_FILE = Path(__file__).with_name("baselines") / "adaptive_solver.json"
REFERENCE_TOLERANCES = {
    "method": "Radau", "rtol": 1e-9, "atol": 1e-11, "persistent": True,
}
REFERENCE_POINTS = 50  # stored per trajectory, at most
ERROR_FLOOR = 5e-7     # errors below this are not told apart from the reference's own


# ---------------------------------------------------------------- models

def _iron_cv():
    import test_iron_oxidation as iron
    return iron._make_cv(pH_target=6.5)


def _iron_davies_cv():
    cv = _iron_cv()
    system = cv.reaction_system
    old = system.engine
    system.attach_engine(NRChemicalEquilibriumEngine(
        old.tableau, activity_model="davies", T_C=old.T_C,
    ))
    cv.equilibrate_to_pH({"K+": 1, "OH-": 1}, 6.5)
    return cv


def _iron_feed_drain_cv():
    cv = _iron_cv()
    cv.boundaries.append(LiquidFeed(
        Q_L_per_h=0.5,
        feed_conc_mol_L={"Fe2+": 8.78e-4, "K+": 0.2204, "H2PO4-": 0.2204},
    ))
    cv.boundaries.append(LiquidDrain(Q_L_per_h=0.5))
    return cv


def _bsm2_cv():
    from vlmodels.adm1.bsm2 import build_bsm2_cv, build_bsm2_reactions
    import test_bsm2_reference as reference
    cv = build_bsm2_cv(build_bsm2_reactions(verbose=False))
    V_L = cv.phases["liquid"].V_L
    for species_id, conc in reference.INITIAL_LIQUID_CONC.items():
        cv.phases["liquid"].n_mol[species_id] = conc * V_L
    return cv


def _weak_acid_cv():
    """A non-stiff model for the explicit methods: no gas link."""
    acid = Species(id="AcOH", atoms={"C": 2, "H": 4, "O": 2}, charge=0)
    base = Species(id="Ac-", atoms={"C": 2, "H": 3, "O": 2}, charge=-1)
    product = Species(id="P", atoms={"C": 2, "H": 4, "O": 2}, charge=0)
    sodium = Species(id="Na+", atoms={"Na": 1}, charge=1)

    def entry(species, coefficient):
        return StoichiometryEntry(species=species, phase="liquid", coefficient=coefficient)

    system = ReactionSystem(
        [
            KineticReaction(
                stoichiometry=[entry(acid, -1), entry(product, +1)],
                rate_fn=lambda env: 2.0 * env.S("AcOH") * env.V_L,
                label="decay",
            ),
            EquilibriumReaction(
                stoichiometry=[entry(H2O, -1), entry(H_plus, +1), entry(OH_minus, +1)],
                log_K=-14.0, label="water",
            ),
            EquilibriumReaction(
                stoichiometry=[entry(acid, -1), entry(base, +1), entry(H_plus, +1)],
                log_K=-4.76, label="acetic",
            ),
        ],
        solver="newton_raphson",
    )
    liquid = LiquidPhase(
        n_mol={"AcOH": 0.05, "Ac-": 0.0, "P": 0.0, "Na+": 0.01, "H+": 1e-5,
               "OH-": 0.0, "H2O": 55.5},
        V_L=1.0, T_K=298.15,
    )
    return ControlVolume(
        phases={"liquid": liquid}, reaction_system=system, label="weak_acid",
        species=[sodium, product],
    )


def _ph_controller(**kwargs):
    return PHController(
        setpoint=6.5, acid_dose="H3PO4", base_dose={"K+": 1, "OH-": 1},
        Kp=1.0, Ki=1.0, max_add_molL_hr=10, **kwargs,
    )


IRON_SPECIES = ("Fe2+", "O2", "Fe(OH)2+")

CASES = {
    "iron_bdf": dict(build=_iron_cv, method="BDF", jac=True, tau_h=2.0, n_steps=200,
                     track=IRON_SPECIES),
    "iron_radau": dict(build=_iron_cv, method="Radau", tau_h=2.0, n_steps=200,
                       track=IRON_SPECIES),
    "iron_lsoda": dict(build=_iron_cv, method="LSODA", tau_h=2.0, n_steps=200,
                       track=IRON_SPECIES),
    "iron_davies_bdf": dict(build=_iron_davies_cv, method="BDF", jac=True, tau_h=2.0,
                            n_steps=200, track=IRON_SPECIES),
    "iron_feed_drain_bdf": dict(build=_iron_feed_drain_cv, method="BDF", jac=True,
                                tau_h=2.0, n_steps=200, track=IRON_SPECIES),
    # A dose fires every macro step: the case with the most restarts.
    "iron_ph_control_bdf": dict(build=_iron_cv, method="BDF", jac=True, tau_h=2.0,
                                n_steps=200, track=IRON_SPECIES,
                                controllers=lambda: [_ph_controller()]),
    # A dose decided every 20 macro steps and held between them.
    "iron_ph_sampled_bdf": dict(build=_iron_cv, method="BDF", jac=True, tau_h=2.0,
                                n_steps=200, track=IRON_SPECIES,
                                controllers=lambda: [_ph_controller(
                                    sample_period_h=0.2, deadband=0.005)]),
    # pH stays inside the band, so the controller runs and never doses.
    "iron_ph_idle_bdf": dict(build=_iron_cv, method="BDF", jac=True, tau_h=2.0,
                             n_steps=200, track=IRON_SPECIES,
                             controllers=lambda: [_ph_controller(
                                 sample_period_h=0.2, deadband=0.05)]),
    # A parameter changes and no amount does.
    "iron_temperature_ramp_bdf": dict(build=_iron_cv, method="BDF", jac=True, tau_h=2.0,
                                      n_steps=200, track=IRON_SPECIES,
                                      profiles=lambda: [TemperatureRamp(
                                          target_cv_key="main",
                                          waypoints=[(0.0, 298.15), (2.0, 303.15)])]),
    "bsm2_bdf": dict(build=_bsm2_cv, method="BDF", tau_h=1.0, n_steps=100,
                     track=("S_ac", "S_h2", "S_ch4", "CO2")),
    "weak_acid_rk45": dict(build=_weak_acid_cv, method="RK45", tau_h=2.0, n_steps=100,
                           track=("AcOH", "P")),
    "weak_acid_rk23": dict(build=_weak_acid_cv, method="RK23", tau_h=2.0, n_steps=100,
                           track=("AcOH", "P")),
    "weak_acid_dop853": dict(build=_weak_acid_cv, method="DOP853", tau_h=2.0, n_steps=100,
                             track=("AcOH", "P")),
}


# ---------------------------------------------------------------- running

def _prepare(case, solver_kwargs):
    kwargs = {"method": case["method"]}
    if case.get("jac"):
        kwargs["use_engine_jacobian"] = True
    kwargs.update(solver_kwargs)
    cv = case["build"]()
    sim = Simulation(
        cvs={"main": cv},
        solver=SimultaneousAdaptiveSolver(**kwargs),
        controllers=case["controllers"]() if "controllers" in case else None,
        profiles=case["profiles"]() if "profiles" in case else None,
        label="bench",
    )
    return cv, sim


def _run_once(case, solver_kwargs):
    cv, sim = _prepare(case, solver_kwargs)
    start = time.process_time()
    result = sim.run(tau_h=case["tau_h"], n_steps=case["n_steps"])
    return cv, result, time.process_time() - start


def _fingerprint(result):
    """SHA-256 over every recorded amount, pH and the time grid."""
    digest = hashlib.sha256(result.t_h.tobytes())
    phases = result.phase_mol["main"]
    for phase_key in sorted(phases):
        for species_id in sorted(phases[phase_key]):
            digest.update(f"{phase_key}.{species_id}".encode())
            digest.update(phases[phase_key][species_id].tobytes())
    digest.update(result.pH["main"].tobytes())
    return digest.hexdigest()


def _stride(case):
    return max(1, case["n_steps"] // REFERENCE_POINTS)


def _trajectories(case, result):
    stride = _stride(case)
    out = {sp: result.liquid_mol["main"][sp][::stride].tolist() for sp in case["track"]}
    out["pH"] = result.pH["main"][::stride].tolist()
    return out


def _errors(case, result, reference):
    """Largest error per tracked species, relative to that species' largest value;
    absolute for pH."""
    now = _trajectories(case, result)
    errors = {}
    for key, ref in reference.items():
        ref_arr, now_arr = np.array(ref, dtype=float), np.array(now[key], dtype=float)
        largest = float(np.nanmax(np.abs(now_arr - ref_arr)))
        scale = 1.0 if key == "pH" else float(np.nanmax(np.abs(ref_arr)))
        errors[key] = largest / scale if scale > 0.0 else largest
    return errors


def _measure(name, solver_kwargs, repeats, reference):
    case = CASES[name]
    # The first run in a process pays for imports and warm-up, so it is not timed.
    _run_once(case, solver_kwargs)
    cpu, fingerprints = [], set()
    for _ in range(repeats):
        cv, result, seconds = _run_once(case, solver_kwargs)
        cpu.append(seconds)
        fingerprints.add(_fingerprint(result))
    return {
        "engine_solves": cv.reaction_system.engine.n_solve_calls,
        "fingerprint": fingerprints.pop() if len(fingerprints) == 1 else "not repeatable",
        "cpu_s_median": statistics.median(cpu),
        "cpu_s_min": min(cpu),
        "cpu_s_max": max(cpu),
        "errors": _errors(case, result, reference) if reference else {},
    }


# ---------------------------------------------------------------- commands

def _load():
    if BASELINE_FILE.exists():
        return json.loads(BASELINE_FILE.read_text(encoding="utf-8"))
    return {"reference": {}, "baseline": {}}


def _save(data):
    BASELINE_FILE.parent.mkdir(exist_ok=True)
    BASELINE_FILE.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")


def _command_reference(names, data):
    for name in names:
        start = time.perf_counter()
        _, result, _ = _run_once(CASES[name], REFERENCE_TOLERANCES)
        data["reference"][name] = _trajectories(CASES[name], result)
        _save(data)
        print(f"{name:28s} reference done in {time.perf_counter() - start:.1f} s", flush=True)


def _command_baseline(names, data, repeats):
    for name in names:
        if name not in data["reference"]:
            sys.exit(f"{name}: no reference stored; run the reference command first.")
        row = _measure(name, {}, repeats, data["reference"][name])
        data["baseline"][name] = row
        _save(data)
        print(f"{name:28s} solves={row['engine_solves']:7d} "
              f"cpu={row['cpu_s_median']:7.3f} s "
              f"worst error={max(row['errors'].values()):.2e}", flush=True)


def _command_check(names, data, repeats, solver_kwargs, exact):
    failures = []
    for name in names:
        base = data["baseline"].get(name)
        if base is None:
            sys.exit(f"{name}: no baseline stored; run the baseline command first.")
        row = _measure(name, solver_kwargs, repeats, data["reference"][name])
        worst = max(row["errors"].values())
        worst_base = max(base["errors"].values())
        same = row["fingerprint"] == base["fingerprint"]
        print(f"{name:28s} solves {base['engine_solves']:7d} -> {row['engine_solves']:7d}  "
              f"cpu {base['cpu_s_median']:7.3f} -> {row['cpu_s_median']:7.3f} s "
              f"(x{base['cpu_s_median'] / row['cpu_s_median']:5.1f})  "
              f"error {worst_base:.2e} -> {worst:.2e}  "
              f"{'identical' if same else 'differs'}", flush=True)
        if exact:
            if not same:
                failures.append(f"{name}: results are not bit-identical")
        else:
            if row["engine_solves"] > base["engine_solves"]:
                failures.append(f"{name}: more engine solves")
            if worst > max(worst_base, ERROR_FLOOR):
                failures.append(f"{name}: larger error")
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("command", choices=("check", "baseline", "reference"))
    parser.add_argument("--cases", nargs="+", choices=sorted(CASES), default=sorted(CASES))
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--exact", action="store_true",
                        help="check: require bit-identical results.")
    parser.add_argument("--solver-kwargs", type=json.loads, default={},
                        help="check: JSON of extra SimultaneousAdaptiveSolver arguments.")
    args = parser.parse_args()

    warnings.simplefilter("ignore")
    data = _load()
    if args.command == "reference":
        _command_reference(args.cases, data)
    elif args.command == "baseline":
        _command_baseline(args.cases, data, args.repeats)
    else:
        failures = _command_check(
            args.cases, data, args.repeats, args.solver_kwargs, args.exact,
        )
        for failure in failures:
            print("FAILED", failure)
        sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
