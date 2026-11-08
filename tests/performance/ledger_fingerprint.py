# -*- coding: utf-8 -*-
"""Fingerprint of fixed runs, for the conservation-ledger phase.

Run by hand from the repo root; pytest does not collect it::

    py tests/performance/ledger_fingerprint.py capture before.json
    py tests/performance/ledger_fingerprint.py capture after.json --cases udm_bisection_plain d2c
    py tests/performance/ledger_fingerprint.py compare before.json after.json
    py tests/performance/ledger_fingerprint.py list

``capture`` builds each case, runs it through ``Simulation.run`` and writes, per
case, every phase's amount of every species at every recorded time, the pH series,
the warnings raised (by category, with the first message of each) and any error.
Floats are written as JSON numbers, which round-trip exactly, so two captures of
the same code compare bit for bit. ``compare`` reports, per case, whether the two
captures are identical and, if not, the largest absolute and relative difference
per phase and species, and the change in warnings.

The cases are those listed in
``docs/dev/implementation/upcoming/CONSERVATION_LEDGER_CHECKLIST.md`` ("Measurement").
They are built here, not imported from the tutorials, so a tutorial edit does not
move the fingerprint; the user-defined model and BSM2's initial state are taken from
their tests, which the checklist names as the cases. Output files are not committed.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import traceback
import warnings
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

REPO = Path(__file__).resolve().parents[2]
for _path in (REPO, REPO / "models", REPO / "tests" / "standalone"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import numpy as np  # noqa: E402

from PyOMES.chemistry import Species  # noqa: E402
from PyOMES.core import (  # noqa: E402
    AdvectiveLink, ControlVolume, GasPhase, Simulation,
)
from PyOMES.core.boundaries import LiquidFeed, PressureReliefVent  # noqa: E402
from PyOMES.core.solvers import (  # noqa: E402
    SimultaneousAdaptiveSolver, SimultaneousEulerSolver,
)
from PyOMES.core import system_solver as _ss  # noqa: E402
from PyOMES.control.cv_loops import PHController  # noqa: E402
from PyOMES.reactions import ReactionSystem  # noqa: E402

R_L_ATM = 0.082057366080960


# ---------------------------------------------------------------- cases

def _udm_cv(engine="bisection", label="user_defined"):
    """The user-defined model of ``test_user_defined_model.py``, with either engine."""
    import test_user_defined_model as udm
    cv = udm._build_cv()
    if engine == "nr":
        reactions = list(cv.reaction_system.reactions)
        cv = ControlVolume(
            phases=dict(cv.phases),
            reaction_system=ReactionSystem(reactions, label="user_defined",
                                           solver="newton_raphson"),
            label=label,
        )
    return cv


def _udm_feed(cv):
    cv.boundaries.append(LiquidFeed(
        Q_L_per_h=0.1, feed_conc_mol_L={"AceticAcid": 1e-2, "O2": 2e-3},
        label="feed",
    ))
    return cv


def _udm_vent(cv):
    """A headspace that starts above the vent's set point."""
    V_gas, T_K = 0.5, 298.15
    n_total = 1.5 * V_gas / (R_L_ATM * T_K)
    gas = GasPhase({"O2": 0.2 * n_total, "CO2": 0.8 * n_total}, V_L=V_gas, T_K=T_K)
    cv = ControlVolume(
        phases={"gas": gas, **cv.phases},
        reaction_system=cv.reaction_system,
        label=cv.label,
    )
    cv.boundaries.append(PressureReliefVent(P_set_atm=1.1, mode="instant"))
    return cv


NA_PLUS = Species(id="Na+", atoms={"Na": 1}, charge=+1)


def _udm_dose(cv):
    """Base dosed by a pH controller; returns (cv, controllers)."""
    cv = ControlVolume(
        phases=dict(cv.phases), reaction_system=cv.reaction_system,
        label=cv.label, species=[NA_PLUS],
    )
    ctrl = PHController(setpoint=4.5, Kp=0.5, Ki=0.0,
                        base_dose={"Na+": 1, "OH-": 1}, max_add_molL_hr=0.05)
    return cv, [ctrl]


def _sim(cvs, *, controllers=None, solver=None, links=None, system_solver=None):
    return Simulation(cvs=cvs, controllers=controllers or [], solver=solver,
                      links=links or [], system_solver=system_solver,
                      label="ledger_fingerprint")


def case_udm(engine, variant):
    def build():
        cv = _udm_cv(engine)
        controllers = None
        if variant == "feed":
            cv = _udm_feed(cv)
        elif variant == "vent":
            cv = _udm_vent(cv)
        elif variant == "dose":
            cv, controllers = _udm_dose(cv)
        return _sim({"main": cv}, controllers=controllers), 1.0, 20
    return build


def case_udm_solver(solver_factory, variant="plain"):
    def build():
        cv = _udm_cv("bisection")
        if variant == "feed":
            cv = _udm_feed(cv)
        return _sim({"main": cv}, solver=solver_factory()), 1.0, 20
    return build


def case_d2c():
    import importlib.util
    path = REPO / "docs" / "tutorials" / "D2C_workshop" / "raw_construction.py"
    spec = importlib.util.spec_from_file_location("_ledger_raw_construction", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.build(), 1.0, 200


#: ADM1 starting state (mol/L) for this fingerprint only. The repo has no canonical
#: one; these are round, moderate values so that every process runs.
ADM1_INITIAL_MOL_L = {
    "Glucose": 5e-3, "AminoAcid": 2e-3, "LCFA": 2e-4, "AceticAcid": 1e-2,
    "Propionate": 2e-3, "Butyrate": 2e-3, "Valerate": 1e-3, "H2": 1e-7,
    "CO2": 2e-2, "NH3": 2e-2, "SO4": 2e-3, "X_composite": 2e-2,
}
ADM1_ORGANISM_MOL_L = 5e-3


def case_adm1():
    from vlmodels.adm1.base import (
        build_adm1_cv, build_adm1_reactions, get_organism_list, seed_adm1_strong_ions,
    )
    cv = build_adm1_cv(build_adm1_reactions(), V_total_L=2.0, T_K=308.15)
    seed_adm1_strong_ions(cv, CT_Na=0.100, CT_Cl=0.030)
    liquid = cv.phases["liquid"]
    for sp, conc in ADM1_INITIAL_MOL_L.items():
        liquid.n_mol[sp] = conc * liquid.V_L
    for org in get_organism_list():
        liquid.n_mol[org] = ADM1_ORGANISM_MOL_L * liquid.V_L
    return _sim({"main": cv}), 0.2, 8


def case_bsm2():
    from vlmodels.adm1.bsm2 import build_bsm2_cv, build_bsm2_reactions
    from test_bsm2_reference import DT_H, INITIAL_LIQUID_CONC, N_STEPS
    cv = build_bsm2_cv(build_bsm2_reactions(verbose=False))
    V = cv.phases["liquid"].V_L
    for sp, conc in INITIAL_LIQUID_CONC.items():
        cv.phases["liquid"].n_mol[sp] = conc * V
    return _sim({"main": cv}), N_STEPS * DT_H, N_STEPS


def _tank(fed):
    """The stirred-tank template as ``fed_batch_fermenter.py`` builds it.

    The batch case is the same tank without the feed: ``batch_fermenter.py``
    seeds no substrate or organism, so nothing in it grows.
    """
    from PyOMES.control.cv_loops import DOAgitationController
    from PyOMES.databases.anaerobic_digestion import AD_BASIC
    from PyOMES.databases.bioprocess_basic import AIR
    from PyOMES.templates.stirred_tank import StirredTankBuilder

    V_total, headspace = 2000.0, 0.20
    builder = (
        StirredTankBuilder()
        .vessel(V_total_L=V_total, headspace_frac=headspace, T_K=305.15)
        .initial_gas(AIR)
        .gas_feed(vvm_min=1.0, composition={"O2": 0.21, "N2": 0.79})
        .transfer_kinetic({"O2": 150.0, "CO2": 135.0}, equilibrium=["N2"])
        .chemistry(chemistry_db=AD_BASIC)
        .organism("Yeast", balance_basis="CHO")
        .substrate("AceticAcid", mu_max=0.5, Ks=5e-3, yield_gX_gS=0.36)
        .label("fed_batch" if fed else "batch")
    )
    cv = builder.build()
    V_liq = cv.phases["liquid"].V_L
    cv.boundaries.append(PressureReliefVent(P_set_atm=1.10, mode="instant"))
    if fed:
        cv.boundaries.append(LiquidFeed(
            Q_L_per_h=5.0, feed_conc_mol_L={"AceticAcid": 50.0 / 60.052},
            label="substrate_feed",
        ))
    controllers = [
        PHController(setpoint=5.5, Kp=0.5, Ki=0.2, acid_dose="H3PO4",
                     base_dose={"Na+": 1, "OH-": 1}, max_add_molL_hr=0.05),
        DOAgitationController(setpoint_mol_L=1e-4, Kp=1e5, Ki=2e4),
    ]
    cv.phases["liquid"].n_mol["AceticAcid"] = (1.0 / 60.052) * V_liq
    cv.phases["liquid"].n_mol["Yeast"] = (0.1 / AD_BASIC.species["Yeast"].MW) * V_liq
    return _sim({"main": cv}, controllers=controllers), 2.0, 100


def case_link(system_solver_factory):
    def build():
        cvs = {"a": _udm_cv(label="a"), "b": _udm_cv(label="b")}
        cvs["b"].phases["liquid"].n_mol["AceticAcid"] *= 2.0
        links = [AdvectiveLink(
            _source_cv_key=x, _source_phase_key="liquid",
            _sink_cv_key=y, _sink_phase_key="liquid", Q_L_per_h=0.5,
        ) for x, y in (("a", "b"), ("b", "a"))]
        system_solver = system_solver_factory() if system_solver_factory else None
        return _sim(cvs, links=links, system_solver=system_solver), 1.0, 20
    return build


CASES = {}
for _engine in ("bisection", "nr"):
    for _variant in ("plain", "feed", "vent", "dose"):
        CASES[f"udm_{_engine}_{_variant}"] = case_udm(_engine, _variant)
CASES.update({
    "udm_euler": case_udm_solver(SimultaneousEulerSolver),
    "udm_adaptive": case_udm_solver(SimultaneousAdaptiveSolver),
    "udm_adaptive_feed": case_udm_solver(SimultaneousAdaptiveSolver, "feed"),
    "udm_adaptive_feed_restart": case_udm_solver(
        lambda: SimultaneousAdaptiveSolver(persistent=False), "feed"),
    "d2c": case_d2c,
    "adm1": case_adm1,
    "bsm2": case_bsm2,
    "tank_batch": lambda: _tank(fed=False),
    "tank_fed_batch": lambda: _tank(fed=True),
    "link_default": case_link(None),
    "link_explicit_euler": case_link(_ss.ExplicitEulerSystemSolver),
    "link_strang": case_link(_ss.StrangSplittingSystemSolver),
    "link_multirate": case_link(_ss.MultirateSystemSolver),
    "link_implicit": case_link(_ss.ImplicitTransportSystemSolver),
    "link_monolithic": case_link(_ss.MonolithicODESolver),
})


# ---------------------------------------------------------------- capture

def _floats(values):
    """A list of floats with NaN and infinities as strings (JSON has neither)."""
    out = []
    for v in np.asarray(values, dtype=float).tolist():
        out.append(v if math.isfinite(v) else repr(v))
    return out


def _run_case(name):
    record = {"warnings": {}, "error": None}
    started = time.process_time()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            sim, tau_h, n_steps = CASES[name]()
            result = sim.run(tau_h=tau_h, n_steps=n_steps)
            record["t_h"] = _floats(result.t_h)
            record["phase_mol"] = {
                cv_key: {
                    pk: {sp: _floats(series) for sp, series in sorted(phases.items())}
                    for pk, phases in sorted(cv_phases.items())
                }
                for cv_key, cv_phases in sorted(result.phase_mol.items())
            }
            record["pH"] = {k: _floats(v) for k, v in sorted(result.pH.items())}
        except Exception:  # recorded, so one failing case does not stop the rest
            record["error"] = traceback.format_exc(limit=6)
    for w in caught:
        key = w.category.__name__
        entry = record["warnings"].setdefault(key, {"count": 0, "first": str(w.message)[:300]})
        entry["count"] += 1
    record["cpu_s"] = round(time.process_time() - started, 2)
    return record


def command_capture(path, names):
    out = {}
    for name in names:
        print(f"{name} ...", end=" ", flush=True)
        out[name] = _run_case(name)
        status = "ERROR" if out[name]["error"] else "ok"
        print(f"{status} ({out[name]['cpu_s']} s CPU)", flush=True)
    Path(path).write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"wrote {path}")


# ---------------------------------------------------------------- compare

def _as_array(values):
    return np.array([float(v) for v in values])


def _compare_case(a, b):
    lines = []
    if a.get("error") or b.get("error"):
        if (a.get("error") or "") != (b.get("error") or ""):
            lines.append("  error differs:")
            for tag, rec in (("A", a), ("B", b)):
                last = (rec.get("error") or "none").strip().splitlines()[-1:]
                lines.append(f"    {tag}: {last[0] if last else 'none'}")
        return lines
    for field in ("t_h",):
        if a[field] != b[field]:
            lines.append(f"  {field} differs")
    for cv_key in sorted(set(a["phase_mol"]) | set(b["phase_mol"])):
        pa, pb = a["phase_mol"].get(cv_key, {}), b["phase_mol"].get(cv_key, {})
        for pk in sorted(set(pa) | set(pb)):
            sa, sb = pa.get(pk, {}), pb.get(pk, {})
            for sp in sorted(set(sa) | set(sb)):
                if sp not in sa or sp not in sb:
                    lines.append(f"  {cv_key}.{pk}.{sp}: only in {'B' if sp not in sa else 'A'}")
                    continue
                if sa[sp] == sb[sp]:
                    continue
                xa, xb = _as_array(sa[sp]), _as_array(sb[sp])
                if xa.shape != xb.shape:
                    lines.append(f"  {cv_key}.{pk}.{sp}: length {xa.shape} vs {xb.shape}")
                    continue
                diff = np.abs(xa - xb)
                scale = np.maximum(np.abs(xa), 1e-300)
                lines.append(
                    f"  {cv_key}.{pk}.{sp}: max abs {np.nanmax(diff):.3e}, "
                    f"max rel {np.nanmax(diff / scale):.3e}, final {xa[-1]:.6e} -> {xb[-1]:.6e}"
                )
    for cv_key in sorted(set(a["pH"]) | set(b["pH"])):
        if a["pH"].get(cv_key) != b["pH"].get(cv_key):
            xa, xb = _as_array(a["pH"].get(cv_key, [])), _as_array(b["pH"].get(cv_key, []))
            if xa.shape == xb.shape:
                lines.append(f"  {cv_key}.pH: max abs {np.nanmax(np.abs(xa - xb)):.3e}, "
                             f"final {xa[-1]:.6f} -> {xb[-1]:.6f}")
            else:
                lines.append(f"  {cv_key}.pH: length differs")
    wa, wb = a["warnings"], b["warnings"]
    for cat in sorted(set(wa) | set(wb)):
        ca = wa.get(cat, {}).get("count", 0)
        cb = wb.get(cat, {}).get("count", 0)
        if ca != cb:
            lines.append(f"  warnings {cat}: {ca} -> {cb}")
    return lines


def command_compare(path_a, path_b):
    a = json.loads(Path(path_a).read_text(encoding="utf-8"))
    b = json.loads(Path(path_b).read_text(encoding="utf-8"))
    identical = True
    for name in sorted(set(a) | set(b)):
        if name not in a or name not in b:
            print(f"{name}: only in {'B' if name not in a else 'A'}")
            identical = False
            continue
        lines = _compare_case(a[name], b[name])
        if lines:
            identical = False
            print(f"{name}: differs")
            print("\n".join(lines))
        else:
            print(f"{name}: identical")
    print("ALL IDENTICAL" if identical else "DIFFERENCES FOUND")
    return 0 if identical else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    cap = sub.add_parser("capture")
    cap.add_argument("out")
    cap.add_argument("--cases", nargs="+", choices=sorted(CASES), default=list(CASES))
    cmp_ = sub.add_parser("compare")
    cmp_.add_argument("a")
    cmp_.add_argument("b")
    sub.add_parser("list")
    args = parser.parse_args()
    if args.command == "capture":
        command_capture(args.out, args.cases)
    elif args.command == "compare":
        sys.exit(command_compare(args.a, args.b))
    else:
        print("\n".join(CASES))


if __name__ == "__main__":
    main()
