# -*- coding: utf-8 -*-
"""Mass balance closure on the user-defined model, under each engine.

The standard (``docs/dev/implementation/upcoming/MASS_BALANCE_CLOSURE.md``): every
element and the charge are conserved to roundoff, with what crosses the control
volume's boundary accounted for. Each test runs the user-defined model of
``test_user_defined_model.py`` for 1 h in 20 steps through ``Simulation`` and checks,
for one element or for charge, that the change in the total equals what crossed the
boundary, to 1e-12 relative. The baseline is taken before the first solve, so the
initial speciation is checked too.

What crossed is recorded here, not taken from the package: a boundary is wrapped so
that what its ``compute_flux`` returns is kept, and a controller's dose is kept as it
reaches ``ControlVolume.apply_external_flux``. The tests therefore do not depend on
the package's own record of what crossed, and can later check it.

Four cases: plain; a ``LiquidFeed``; a headspace that starts above a
``PressureReliefVent``'s set point; a base dose (Na+ with OH-) by a ``PHController``.
Tests that fail today are strict expected failures whose reason gives today's number
and its cause.
"""
import importlib.util
from pathlib import Path

import pytest

from PyOMES.chemistry import Species
from PyOMES.control.cv_loops import PHController
from PyOMES.core import ControlVolume, GasPhase, Simulation
from PyOMES.core.boundaries import LiquidFeed, PressureReliefVent
from PyOMES.reactions import ReactionSystem
from PyOMES.units import R_L_ATM_PER_MOL_K


def _load_user_defined_model():
    """The model of ``test_user_defined_model.py``, loaded from its file (test
    modules are not importable by name here)."""
    path = Path(__file__).with_name("test_user_defined_model.py")
    spec = importlib.util.spec_from_file_location("_closure_user_defined_model", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_UDM = _load_user_defined_model()
USER_SPECIES, _build_cv = _UDM.USER_SPECIES, _UDM._build_cv

NA_PLUS = Species(id="Na+", atoms={"Na": 1}, charge=+1)
SPECIES = dict(USER_SPECIES, **{NA_PLUS.id: NA_PLUS})

RTOL = 1e-12
TAU_H, N_STEPS = 1.0, 20


# ── What crossed the boundary ────────────────────────────────────────────────

def _add(acc, amounts_mol):
    """Add the element and charge content of {species: mol} to ``acc``."""
    for sp_id, n in amounts_mol.items():
        species = SPECIES[sp_id]
        for element, count in species.atoms.items():
            acc[element] = acc.get(element, 0.0) + n * count
        acc["charge"] = acc.get("charge", 0.0) + n * species.charge


class _Recording:
    """A boundary that keeps what the boundary it wraps returned."""

    def __init__(self, inner, crossed):
        self._inner = inner
        self._crossed = crossed

    def compute_flux(self, cv, dt_h, instantaneous=False):
        flux = self._inner.compute_flux(cv, dt_h, instantaneous=instantaneous)
        _add(self._crossed, {sp: rate * dt_h for sp, rate in flux.items()})
        return flux

    def __getattr__(self, name):
        return getattr(self._inner, name)


def _record_doses(cv, crossed):
    apply = cv.apply_external_flux

    def recorded(phase_key, flux_mol_per_h, dt_h):
        _add(crossed, {sp: rate * dt_h for sp, rate in flux_mol_per_h.items()})
        return apply(phase_key, flux_mol_per_h, dt_h)

    cv.apply_external_flux = recorded


# ── The cases ────────────────────────────────────────────────────────────────

def _build(engine, case, crossed):
    base = _build_cv()
    reaction_system = base.reaction_system
    if engine == "nr":
        reaction_system = ReactionSystem(
            list(reaction_system.reactions), label="user_defined", solver="newton_raphson",
        )
    phases = dict(base.phases)
    if case == "vent":
        V_gas, T_K = 0.5, 298.15
        n_gas = 1.5 * V_gas / (R_L_ATM_PER_MOL_K * T_K)  # starts at 1.5 atm
        phases["gas"] = GasPhase({"O2": 0.2 * n_gas, "CO2": 0.8 * n_gas}, V_L=V_gas, T_K=T_K)
    cv = ControlVolume(
        phases=phases, reaction_system=reaction_system, label="closure",
        species=[NA_PLUS] if case == "dose" else None,
    )
    controllers = []
    if case == "feed":
        cv.boundaries.append(_Recording(LiquidFeed(
            Q_L_per_h=0.1, feed_conc_mol_L={"AceticAcid": 1e-2, "O2": 2e-3},
        ), crossed))
    elif case == "vent":
        cv.boundaries.append(_Recording(PressureReliefVent(P_set_atm=1.1, mode="instant"), crossed))
    elif case == "dose":
        controllers.append(PHController(
            setpoint=4.5, Kp=0.5, Ki=0.0, base_dose={"Na+": 1, "OH-": 1}, max_add_molL_hr=0.05,
        ))
        _record_doses(cv, crossed)
    return cv, controllers


def _inventory(cv):
    """Element totals, charge, and the ionic inventory sum(|z| n) over every phase."""
    totals = {"charge": 0.0, "ions": 0.0}
    for phase in cv.phases.values():
        for sp_id, n in phase.n_mol.items():
            species = SPECIES[sp_id]
            for element, count in species.atoms.items():
                totals[element] = totals.get(element, 0.0) + n * count
            totals["charge"] += n * species.charge
            totals["ions"] += abs(n * species.charge)
    return totals


_RUNS = {}


def _run(engine, case):
    """(before, after, crossed), run once per engine and case."""
    key = (engine, case)
    if key not in _RUNS:
        crossed = {}
        cv, controllers = _build(engine, case, crossed)
        before = _inventory(cv)
        Simulation(cvs={"main": cv}, controllers=controllers, label="closure").run(
            tau_h=TAU_H, n_steps=N_STEPS,
        )
        _RUNS[key] = (before, _inventory(cv), crossed)
    return _RUNS[key]


# ── Expected failures, with today's numbers ──────────────────────────────────

def _xfail(reason):
    return pytest.mark.xfail(strict=True, reason=reason)


_BIS_WATER = _xfail(
    "The Bisection engine does not book water: the H2O that CO2 + H2O <-> HCO3- + H+ "
    "consumes is not taken from n_mol, so each HCO3- formed adds 2 H and 1 O "
    "(+4.7e-6 mol H, +4.3e-8 relative, over the first solve and 1 h)."
)
_BIS_DOSE = _xfail(
    "The Bisection engine forces the liquid neutral and does not book water: the OH- "
    "dosed with Na+ is replaced by the pH the charge balance gives, and the water it "
    "forms is not added, so close to 2 H and 1 O are lost per OH- (5.0e-3 mol dosed; "
    "-9.9e-3 mol H, -4.95e-3 mol O, -8.9e-5 relative, over 1 h)."
)
_NR_ELEMENTS = _xfail(
    "The NR engine sets liquid water to 55.509 mol/L x V_L on every solve (55.5 mol "
    "is given, so the first solve adds 9 mmol) and conserves component totals only to "
    "its residual tolerance (H and O +1.5e-4 relative, +6.4e-5 with the dose; C +1.5e-8 "
    "to +2.4e-8; over the first solve and 1 h)."
)
_NR_CHARGE = _xfail(
    "The NR engine conserves component totals, and so charge, only to its residual "
    "tolerance (1e-10 mol/L): about -4e-14 mol of charge over 1 h, -4e-13 with "
    "the dose, 3e-11 to 4e-11 of the ionic inventory."
)

CASES = ("plain", "feed", "vent", "dose")


def _params(engine):
    out = []
    for case in CASES:
        for quantity in ("C", "H", "O", "Na", "charge"):
            if quantity == "Na" and case != "dose":
                continue
            marks = []
            if engine == "bisection" and quantity in ("H", "O"):
                marks = [_BIS_DOSE if case == "dose" else _BIS_WATER]
            elif engine == "nr" and quantity in ("C", "H", "O"):
                marks = [_NR_ELEMENTS]
            elif engine == "nr" and quantity == "charge":
                marks = [_NR_CHARGE]
            out.append(pytest.param(case, quantity, marks=marks, id=f"{case}-{quantity}"))
    return out


def _check(engine, case, quantity):
    before, after, crossed = _run(engine, case)
    change = after.get(quantity, 0.0) - before.get(quantity, 0.0)
    imbalance = change - crossed.get(quantity, 0.0)
    if quantity == "charge":
        scale = max(before["ions"], after["ions"])
    else:
        scale = max(abs(before.get(quantity, 0.0)), abs(after.get(quantity, 0.0)))
    assert abs(imbalance) <= RTOL * scale, (
        f"{engine}, {case}: {quantity} changed by {change:+.6e} mol and "
        f"{crossed.get(quantity, 0.0):+.6e} mol crossed the boundary, an imbalance of "
        f"{imbalance:+.3e} mol ({imbalance / scale:+.2e} of {scale:.4g})"
    )


@pytest.mark.parametrize("case,quantity", _params("bisection"))
def test_bisection_closes(case, quantity):
    _check("bisection", case, quantity)


@pytest.mark.parametrize("case,quantity", _params("nr"))
def test_nr_closes(case, quantity):
    _check("nr", case, quantity)


def test_the_boundaries_move_material():
    """The feed, vent and dose cases are not trivially closed: something crosses."""
    for case, quantity in (("feed", "C"), ("vent", "C"), ("dose", "Na")):
        _, _, crossed = _run("bisection", case)
        assert abs(crossed.get(quantity, 0.0)) > 1e-6, (case, quantity, crossed)
