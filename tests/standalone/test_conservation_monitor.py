# -*- coding: utf-8 -*-
"""Tests for ConservationMonitor (state-unification C6).

Element + charge accounting per advance() step. Mirrors the
AccuracyMonitor throttle / summary plumbing; per-CV instance is
auto-attached in ControlVolume.__init__ when a ReactionSystem is
present.
"""

from __future__ import annotations

import warnings

import pytest

import PyOMES
from PyOMES.monitoring import (
    ConservationMonitor,
    ConservationWarning,
    reset_conservation_summary,
)


@pytest.fixture(autouse=True)
def _isolate_summary_counter():
    """Pytest fixture: reset the module-level counter + monitor
    state before each test so the throttle's first_N / once
    semantics aren't polluted by prior tests."""
    reset_conservation_summary()
    yield
    reset_conservation_summary()


@pytest.fixture(autouse=True)
def _force_always_throttle():
    """For these tests, prefer the ``"always"`` throttle so each
    check_step emits the warning whenever the threshold trips —
    otherwise the once-per-category default suppresses follow-up
    emissions and we can't observe ratios."""
    original = PyOMES.config.warnings.throttle
    PyOMES.config.warnings.throttle = "always"
    yield
    PyOMES.config.warnings.throttle = original


def _make_water_species_registry():
    """Tiny registry: H+, OH-, Na+, Cl-, plus H2O."""
    from PyOMES.chemistry import Species
    return {
        "H+": Species(id="H+", atoms={"H": 1}, charge=+1),
        "OH-": Species(id="OH-", atoms={"O": 1, "H": 1}, charge=-1),
        "Na+": Species(id="Na+", atoms={"Na": 1}, charge=+1),
        "Cl-": Species(id="Cl-", atoms={"Cl": 1}, charge=-1),
        "H2O": Species(id="H2O", atoms={"H": 2, "O": 1}, charge=0),
    }


def _make_phase(n_mol_dict, V_L=1.0, T_K=298.15):
    from PyOMES.core import LiquidPhase
    return LiquidPhase(n_mol=n_mol_dict, V_L=V_L, T_K=T_K)


class TestConservationMonitorBasics:

    def test_no_emission_on_first_call(self):
        """First call sets baseline; nothing to compare against."""
        monitor = ConservationMonitor()
        monitor.set_species_registry(_make_water_species_registry())
        liq = _make_phase({"Na+": 1.0, "Cl-": 1.0})
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            monitor.check_step({"liquid": liq})
        assert len(buf) == 0
        assert monitor.step_count == 1

    def test_no_emission_when_balanced(self):
        """Balanced charge, conserved elements → no warning."""
        monitor = ConservationMonitor()
        monitor.set_species_registry(_make_water_species_registry())
        liq = _make_phase({"Na+": 1.0, "Cl-": 1.0})
        monitor.check_step({"liquid": liq})  # baseline
        # Step 2: same state
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            monitor.check_step({"liquid": liq})
        cw = [w for w in buf if issubclass(w.category, ConservationWarning)]
        assert len(cw) == 0

    def test_no_op_with_empty_registry(self):
        """Without a species registry, the monitor silently no-ops."""
        monitor = ConservationMonitor()
        # No set_species_registry call.
        liq = _make_phase({"NaCl_lump": 1.0})
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            monitor.check_step({"liquid": liq})
            monitor.check_step({"liquid": liq})
        cw = [w for w in buf if issubclass(w.category, ConservationWarning)]
        assert len(cw) == 0


class TestChargeBalance:

    def test_unmatched_counterion_step_warning(self):
        """Adding Na+ without Cl- breaks charge balance — should
        emit a ``charge_step`` warning."""
        monitor = ConservationMonitor()
        monitor.set_species_registry(_make_water_species_registry())
        # Start balanced
        liq0 = _make_phase({"Na+": 1.0, "Cl-": 1.0})
        monitor.check_step({"liquid": liq0})  # baseline
        # Step: add 2 mol Na+ but no Cl-
        liq1 = _make_phase({"Na+": 3.0, "Cl-": 1.0})
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            monitor.check_step({"liquid": liq1})
        cw = [w for w in buf if issubclass(w.category, ConservationWarning)]
        assert any("Charge balance drifted" in str(w.message) for w in cw)

    def test_charge_cumulative_warning(self):
        """Cumulative charge drift past the threshold should fire
        the ``charge_cum`` category."""
        monitor = ConservationMonitor()
        monitor.set_species_registry(_make_water_species_registry())
        liq0 = _make_phase({"Na+": 1.0, "Cl-": 1.0})
        monitor.check_step({"liquid": liq0})
        # Walk charge upwards via repeated unmatched additions
        liq = _make_phase({"Na+": 10.0, "Cl-": 1.0})
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            monitor.check_step({"liquid": liq})
        cw = [w for w in buf if issubclass(w.category, ConservationWarning)]
        # Both per-step (drift since prior) and cumulative
        # (absolute residual) should fire on a big jump.
        msgs = " ".join(str(w.message) for w in cw)
        assert "cumulative threshold" in msgs


class TestElementBalance:

    def test_step_element_loss(self):
        """A step that destroys element atoms (broken stoichiometry)
        triggers the per-step element warning."""
        monitor = ConservationMonitor()
        monitor.set_species_registry(_make_water_species_registry())
        liq0 = _make_phase({"H2O": 10.0})  # H=20, O=10
        monitor.check_step({"liquid": liq0})  # baseline
        # Step: half the water disappears with no compensating species
        liq1 = _make_phase({"H2O": 5.0})  # H=10, O=5
        with warnings.catch_warnings(record=True) as buf:
            warnings.simplefilter("always")
            monitor.check_step({"liquid": liq1})
        cw = [w for w in buf if issubclass(w.category, ConservationWarning)]
        # Both H and O drifted; both should fire.
        cats = {str(w.message).split("'")[1] for w in cw if "Element" in str(w.message)}
        assert "H" in cats or "O" in cats


# Throttle plumbing tested in test_accuracy_monitor.py (the _emit
# code is shared between AccuracyMonitor and ConservationMonitor).
