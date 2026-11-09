# -*- coding: utf-8 -*-
"""Record every clamp and floor that changes an amount, for the conservation-ledger phase.

Run by hand from the repo root; pytest does not collect it.

As a script runner, on the cases of ``ledger_fingerprint.py`` or on Python scripts::

    py tests/performance/clamp_audit.py cases
    py tests/performance/clamp_audit.py --out audit.json cases --only d2c bsm2
    py tests/performance/clamp_audit.py scripts docs/tutorials/templates/fed_batch_fermenter.py
    py tests/performance/clamp_audit.py --no-clamp --no-trial-floor cases
    py tests/performance/clamp_audit.py --no-clamp scripts docs/tutorials/templates/cstr_fermenter.py --log-dir out

(The options before the subcommand apply to the whole run; ``--log-dir`` writes each
script's printed output to a file, so a run with and one without ``--no-clamp`` can be
compared.)

On notebooks (each runs in its own kernel; a first cell installs the recorder and a
last cell writes its report)::

    py tests/performance/clamp_audit.py notebooks docs/tutorials/templates/batch_fermenter.ipynb

As a pytest plugin, per test::

    $env:PYTHONPATH = "tests/performance"
    py -m pytest -p clamp_audit -p no:cacheprovider -q --clamp-audit-out audit.json

What is recorded, per run (a case, a script, a notebook or a test), per guard:

``clamp_fn``
    A step solver's clamp (``SequentialAdvanceSolver``, ``SimultaneousEulerSolver``),
    whatever function was passed. Records the removal it cut, in mol.
``phase_floor``
    ``apply_flux(..., clamp=True)`` on a gas, liquid or solid phase setting a negative
    result to zero. Records the material it added.
``trial_floor``, ``final_floor``, ``jacobian_floor``, ``unpack_floor``
    ``floor_nonnegative`` in the adaptive solver: on a trial state in the right-hand
    side, on the state written back after a step, in the Jacobian, and in
    ``StateVector.unpack(floor=True)`` on the retry and restore paths. Records what a
    floored vector gained. The trial and Jacobian floors change only what rates are
    evaluated at; the final and unpack floors change the amounts.
``integrate_reactions``
    Any call of ``ControlVolume._integrate_reactions`` (it has no caller in the repo).
``nonfinite_state``
    Not a guard: a NaN or infinite value in a state the adaptive solver hands to a
    floor (a trial state, a written-back state). Records the number of such values.

Each guard's entry holds the number of calls that changed something, the number of
values changed, the total amount (mol, positive) and the species most affected.

``--no-clamp`` makes the step solvers' default ``clamp_fn`` ``None``; ``--no-trial-floor``
evaluates the adaptive solver's right-hand side and Jacobian at the unfloored trial
state, and its ``trial_floor`` and ``jacobian_floor`` entries then count the negative
values that were left as they are. Both are for measuring what depends on them
(checkpoint 2a); the other floors stay as they are and are still recorded. Under
``persistent=True`` the engine-resolved species are not floored and are not counted.
In the pytest plugin, work done by a module-scoped fixture is counted against the
first test that uses it. This does not see ``if`` guards, ``np.clip``,
or the guards inside rate laws, links and engines; ``bench_adaptive_solver.py clamps``
counts some of those.
"""
from __future__ import annotations

import argparse
import json
import os
import runpy
import sys
import traceback
import warnings
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _path in (REPO, REPO / "models", REPO / "tests" / "standalone", HERE):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import numpy as np  # noqa: E402

GUARDS = ("clamp_fn", "phase_floor", "trial_floor", "final_floor", "jacobian_floor",
          "unpack_floor", "integrate_reactions", "nonfinite_state")


class _Recorder:
    def __init__(self):
        self.reports = {}
        self.errors = {}
        self.run = "(none)"

    def start(self, run):
        self.run = run
        self.reports.setdefault(run, {})

    def add(self, guard, species_amounts, n_values=None):
        """species_amounts: {species or index: positive amount}."""
        if not species_amounts:
            return
        entry = self.reports.setdefault(self.run, {}).setdefault(
            guard, {"calls": 0, "values": 0, "amount_mol": 0.0, "by_species": defaultdict(float)})
        entry["calls"] += 1
        entry["values"] += n_values if n_values is not None else len(species_amounts)
        for sp, amount in species_amounts.items():
            entry["amount_mol"] += amount
            entry["by_species"][str(sp)] += amount

    def as_json(self):
        out = {}
        for run, guards in self.reports.items():
            out[run] = {}
            for guard, e in guards.items():
                top = sorted(e["by_species"].items(), key=lambda kv: -kv[1])[:5]
                out[run][guard] = {"calls": e["calls"], "values": e["values"],
                                   "amount_mol": e["amount_mol"], "top_species": top}
        return out


RECORDER = _Recorder()
_INSTALLED = False


def _wrap_clamp(fn):
    if fn is None or getattr(fn, "_clamp_audit", False):
        return fn

    def recorded(deltas, current_mol, dt_h, *args, **kwargs):
        out = fn(deltas, current_mol, dt_h, *args, **kwargs)
        cut = {sp: abs(float(r) - float(out.get(sp, 0.0))) * float(dt_h)
               for sp, r in deltas.items() if float(out.get(sp, 0.0)) != float(r)}
        RECORDER.add("clamp_fn", cut)
        return out

    recorded._clamp_audit = True
    return recorded


def install(*, no_clamp=False, no_trial_floor=False):
    """Patch the guards to record into ``RECORDER``. Idempotent."""
    global _INSTALLED
    if _INSTALLED:
        return
    _INSTALLED = True

    from PyOMES.core import clamping, control_volume, phases, solvers, state_vector  # noqa: F401

    # Step solvers' clamp_fn (a default argument, so wrapped after __init__).
    for cls in (solvers.SequentialAdvanceSolver, solvers.SimultaneousEulerSolver):
        original_init = cls.__init__

        def init(self, *args, _original=original_init, **kwargs):
            if no_clamp and not args and "clamp_fn" not in kwargs:
                kwargs["clamp_fn"] = None
            _original(self, *args, **kwargs)
            self.clamp_fn = _wrap_clamp(self.clamp_fn)

        cls.__init__ = init

    # Phase floors.
    for cls in (phases.GasPhase, phases.LiquidPhase, phases.SolidPhase):
        original_apply = cls.apply_flux

        def apply_flux(self, flux_mol_per_h, dt_h, *, clamp=True, _original=original_apply):
            if clamp:
                dt = float(dt_h)
                added = {}
                for sp, rate in flux_mol_per_h.items():
                    new_value = self._n_mol.get(sp, 0.0) + float(rate) * dt
                    if new_value < 0.0:
                        added[sp] = -new_value
                RECORDER.add("phase_floor", added)
            return _original(self, flux_mol_per_h, dt_h, clamp=clamp)

        cls.apply_flux = apply_flux

    # floor_nonnegative, labelled by where it is called from.
    real_floor = clamping.floor_nonnegative

    def _label():
        names = []
        frame = sys._getframe(2)
        for _ in range(4):
            if frame is None:
                break
            names.append(frame.f_code.co_name)
            frame = frame.f_back
        # The caller is floor_state (the persistent integrator's wrapper) or the
        # solver code itself: f (right-hand side), jac_callable, solve_step (the
        # write-back after a step) or StateVector.unpack.
        if names[:1] == ["unpack"]:
            return "unpack_floor"
        if names[:1] == ["jac_callable"]:
            return "jacobian_floor"
        if "solve_step" in names[:3] and "f" not in names[:2]:
            return "final_floor"
        return "trial_floor"

    def _context():
        """(indices the caller restores after flooring, index -> species name)."""
        kept, names = (), {}
        frame = sys._getframe(2)
        if frame is not None and frame.f_code.co_name == "floor_state":
            kept = frame.f_locals.get("kept", ()) or ()
        for _ in range(6):
            if frame is None:
                break
            sv = frame.f_locals.get("sv")
            index_map = getattr(sv, "index_map", None)
            if callable(index_map):
                index_map = index_map()
            if isinstance(index_map, dict):
                names = {i: f"{pk}.{sp}" for (_cv, pk, sp), i in index_map.items()}
                break
            layout = frame.f_locals.get("phase_layout")
            offsets = frame.f_locals.get("phase_offset")
            if isinstance(layout, dict) and isinstance(offsets, dict):
                names = {offsets[pk] + j: f"{pk}.{sp}"
                         for pk, species in layout.items() if pk in offsets
                         for j, sp in enumerate(species)}
                break
            frame = frame.f_back
        return set(kept), names

    def floor_nonnegative(y):
        label = _label()
        y_arr = np.asarray(y, dtype=float)
        nonfinite = ~np.isfinite(y_arr)
        if nonfinite.any():
            RECORDER.add("nonfinite_state", {label: float(nonfinite.sum())},
                         n_values=int(nonfinite.sum()))
        negative = y_arr < 0.0
        if negative.any():
            kept, names = _context()
            gained = {names.get(int(i), f"index {int(i)}"): float(-y_arr[i])
                      for i in np.flatnonzero(negative) if int(i) not in kept}
            RECORDER.add(label, gained)
        if no_trial_floor and label in ("trial_floor", "jacobian_floor"):
            return np.array(y_arr, copy=True)
        return real_floor(y)

    clamping.floor_nonnegative = floor_nonnegative
    solvers.floor_nonnegative = floor_nonnegative

    # ControlVolume._integrate_reactions (no caller expected).
    original_integrate = control_volume.ControlVolume._integrate_reactions

    def integrate_reactions(self, *args, **kwargs):
        RECORDER.add("integrate_reactions", {"(call)": 1.0})
        return original_integrate(self, *args, **kwargs)

    control_volume.ControlVolume._integrate_reactions = integrate_reactions


# ---------------------------------------------------------------- runners

def _run(name, body):
    RECORDER.start(name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            body()
            return None
        except BaseException as exc:  # a script's sys.exit or error is recorded, not raised
            if isinstance(exc, KeyboardInterrupt):
                raise
            RECORDER.errors[name] = traceback.format_exc(limit=4)
            return exc


def run_cases(names):
    import ledger_fingerprint as lf
    for name in names:
        def body(name=name):
            sim, tau_h, n_steps = lf.CASES[name]()
            sim.run(tau_h=tau_h, n_steps=n_steps)
        print(f"{name} ...", flush=True)
        _run(name, body)


def run_scripts(paths, log_dir=None, suffix=""):
    import contextlib
    for name in paths:
        path = Path(name).resolve()
        print(f"{name} ...", flush=True)
        cwd = os.getcwd()

        def body(path=path):
            log = None
            if log_dir:
                Path(log_dir).mkdir(parents=True, exist_ok=True)
                log = open(Path(log_dir) / f"{path.stem}{suffix}.txt", "w", encoding="utf-8")
            os.chdir(path.parent)
            try:
                with contextlib.ExitStack() as stack:
                    if log is not None:
                        stack.enter_context(contextlib.redirect_stdout(log))
                        stack.enter_context(contextlib.redirect_stderr(log))
                    runpy.run_path(str(path), run_name="__main__")
            finally:
                os.chdir(cwd)
                if log is not None:
                    log.close()

        _run(str(name), body)


def run_notebooks(paths, out_dir, flags):
    import nbclient
    import nbformat
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for path in paths:
        path = Path(path)
        report = (out_dir / (path.stem + ".clamp_audit.json")).resolve()
        nb = nbformat.read(str(path), as_version=4)
        setup = (
            "import sys\n"
            f"sys.path.insert(0, {str(HERE)!r})\n"
            "import clamp_audit\n"
            f"clamp_audit.install(no_clamp={flags['no_clamp']}, "
            f"no_trial_floor={flags['no_trial_floor']})\n"
            f"clamp_audit.RECORDER.start({str(path)!r})\n"
        )
        teardown = (
            "import json, clamp_audit\n"
            f"open({str(report)!r}, 'w').write(json.dumps(clamp_audit.RECORDER.as_json()))\n"
        )
        nb.cells.insert(0, nbformat.v4.new_code_cell(setup))
        nb.cells.append(nbformat.v4.new_code_cell(teardown))
        print(f"{path} ...", flush=True)
        try:
            nbclient.NotebookClient(nb, timeout=3600,
                                    resources={"metadata": {"path": str(path.parent)}}).execute()
            RECORDER.reports.update(json.loads(report.read_text()))
        except Exception:
            RECORDER.errors[str(path)] = traceback.format_exc(limit=4)


def _print(reports):
    for run, guards in reports.items():
        fired = {g: e for g, e in guards.items() if g != "error"}
        error = guards.get("error")
        if not fired and not error:
            print(f"{run}: nothing fired")
            continue
        print(f"{run}:")
        for guard in GUARDS:
            if guard in fired:
                e = fired[guard]
                top = ", ".join(f"{sp} {amt:.3g}" for sp, amt in e["top_species"])
                print(f"  {guard:<20} calls {e['calls']:>6}  values {e['values']:>7}  "
                      f"total {e['amount_mol']:.4g} mol  [{top}]")
        if error:
            print("  ERROR: " + error.strip().splitlines()[-1])


def main():
    # Library messages print arrows; a redirected Windows console is cp1252.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", help="write the report as JSON")
    parser.add_argument("--no-clamp", action="store_true")
    parser.add_argument("--no-trial-floor", action="store_true")
    sub = parser.add_subparsers(dest="command", required=True)
    c = sub.add_parser("cases")
    c.add_argument("--only", nargs="+")
    s = sub.add_parser("scripts")
    s.add_argument("paths", nargs="+")
    s.add_argument("--log-dir", help="write each script's printed output to "
                   "<log-dir>/<name>.txt (with .no-clamp / .no-trial-floor added)")
    n = sub.add_parser("notebooks")
    n.add_argument("paths", nargs="+")
    n.add_argument("--report-dir", default=os.environ.get("TEMP", "."))
    args = parser.parse_args()

    flags = {"no_clamp": args.no_clamp, "no_trial_floor": args.no_trial_floor}
    if args.command == "notebooks":
        run_notebooks(args.paths, args.report_dir, flags)
        reports = dict(RECORDER.reports)
    else:
        install(**flags)
        if args.command == "cases":
            import ledger_fingerprint as lf
            run_cases(args.only or list(lf.CASES))
        else:
            suffix = (".no-clamp" if args.no_clamp else "") + (
                ".no-trial-floor" if args.no_trial_floor else "")
            run_scripts(args.paths, args.log_dir, suffix)
        reports = RECORDER.as_json()
    for run, error in RECORDER.errors.items():
        reports.setdefault(run, {})["error"] = error
    _print(reports)
    if args.out:
        Path(args.out).write_text(json.dumps(reports, indent=1), encoding="utf-8")
        print(f"wrote {args.out}")


# ---------------------------------------------------------------- pytest plugin

def pytest_addoption(parser):
    parser.addoption("--clamp-audit-out", default=None,
                     help="write the clamp audit's per-test report here (JSON)")
    parser.addoption("--clamp-audit-no-clamp", action="store_true")
    parser.addoption("--clamp-audit-no-trial-floor", action="store_true")


def pytest_configure(config):
    install(no_clamp=config.getoption("--clamp-audit-no-clamp"),
            no_trial_floor=config.getoption("--clamp-audit-no-trial-floor"))


def pytest_runtest_setup(item):
    RECORDER.start(item.nodeid)


def pytest_unconfigure(config):
    reports = {run: g for run, g in RECORDER.as_json().items() if g}
    out = config.getoption("--clamp-audit-out")
    if out:
        Path(out).write_text(json.dumps(reports, indent=1), encoding="utf-8")
    print(f"\nclamp audit: {len(reports)} tests where a clamp or floor fired"
          + (f"; report in {out}" if out else ""))


if __name__ == "__main__":
    main()
