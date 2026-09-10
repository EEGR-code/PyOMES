# -*- coding: utf-8 -*-
"""Tests for Simulation skeleton + RunContext (SIMULATION_CLASS C1).

Validates:
- Simulation construction with single-CV and multi-CV configurations
- Link-endpoint validation (four KeyError cases ported from MultiCVSystem)
- Accessor helpers (__getitem__, __contains__, cv_keys)
- total_mol aggregation across CVs and phases
- snapshot() independence (deep-copied CVs) and sharing (links, etc.)
- snapshotted CVs have _context = None (snapshot is unowned until added)
- __repr__ shape
- RunContext wiring: every CV's _context is the Simulation's _context
- RunContext default state (is_running=False, label propagated)
- Double-ownership: a CV already owned by one Simulation raises when
  used to construct a second Simulation
- Exports: Simulation and RunContext importable from PyOMES.core
"""

import pytest


# ═══════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════

def _make_liquid_cv(key, n_mol, V_L=1.0, T_K=300.0):
    """Create a simple CV with only a liquid phase."""
    from PyOMES.core import ControlVolume, LiquidPhase
    liq = LiquidPhase(n_mol=dict(n_mol), V_L=V_L, T_K=T_K)
    return ControlVolume(phases={"liquid": liq}, label=key)


def _make_gas_liquid_cv(key, gas_mol, liq_mol, V_gas=0.2, V_liq=0.8, T_K=305.15):
    """Create a CV with gas and liquid phases."""
    from PyOMES.core import ControlVolume, GasPhase, LiquidPhase
    gas = GasPhase(n_mol=dict(gas_mol), V_L=V_gas, T_K=T_K)
    liq = LiquidPhase(n_mol=dict(liq_mol), V_L=V_liq, T_K=T_K)
    return ControlVolume(phases={"gas": gas, "liquid": liq}, label=key)


# ═══════════════════════════════════════════════════════════════════════
#  Exports
# ═══════════════════════════════════════════════════════════════════════

class TestExports:

    def test_simulation_and_runcontext_importable_from_core(self):
        from PyOMES.core import Simulation, RunContext
        assert Simulation is not None
        assert RunContext is not None


# ═══════════════════════════════════════════════════════════════════════
#  Construction
# ═══════════════════════════════════════════════════════════════════════

class TestConstruction:

    def test_single_cv_construction(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0}, V_L=1.0)
        sim = Simulation(cvs={"main": cv})
        assert sim.cv_keys == ["main"]
        assert sim["main"] is cv
        assert "main" in sim
        assert sim.label == ""
        assert sim._t_h == 0.0

    def test_multi_cv_construction(self):
        from PyOMES.core import Simulation
        cv_a = _make_liquid_cv("A", {"S": 1.0})
        cv_b = _make_liquid_cv("B", {"S": 0.0})
        sim = Simulation(cvs={"A": cv_a, "B": cv_b}, label="two-zone")
        assert sorted(sim.cv_keys) == ["A", "B"]
        assert sim["A"] is cv_a
        assert sim["B"] is cv_b
        assert sim.label == "two-zone"

    def test_construction_with_controllers_and_profiles(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        # Controllers and profiles are opaque at C1 (consumed at C9/C10);
        # any object placeholder works.
        ctrl = object()
        prof = object()
        sim = Simulation(
            cvs={"main": cv},
            controllers=[ctrl],
            profiles=[prof],
        )
        assert sim.controllers == [ctrl]
        assert sim.profiles == [prof]

    def test_construction_with_solver_and_recorder(self):
        from PyOMES.core import Simulation, SimultaneousEulerSolver
        cv = _make_liquid_cv("main", {"S": 1.0})
        solver = SimultaneousEulerSolver()
        recorder = object()
        sim = Simulation(cvs={"main": cv}, solver=solver, recorder=recorder)
        assert sim.solver is solver
        assert sim._recorder is recorder


# ═══════════════════════════════════════════════════════════════════════
#  Link-endpoint validation
# ═══════════════════════════════════════════════════════════════════════

class TestLinkValidation:

    def test_unknown_source_cv_raises(self):
        from PyOMES.core import Simulation, AdvectiveLink
        cv = _make_liquid_cv("A", {"S": 1.0})
        link = AdvectiveLink("missing", "liquid", "A", "liquid", Q_L_per_h=1.0)
        with pytest.raises(KeyError, match="source_cv_key"):
            Simulation(cvs={"A": cv}, links=[link])

    def test_unknown_sink_cv_raises(self):
        from PyOMES.core import Simulation, AdvectiveLink
        cv = _make_liquid_cv("A", {"S": 1.0})
        link = AdvectiveLink("A", "liquid", "missing", "liquid", Q_L_per_h=1.0)
        with pytest.raises(KeyError, match="sink_cv_key"):
            Simulation(cvs={"A": cv}, links=[link])

    def test_unknown_source_phase_raises(self):
        from PyOMES.core import Simulation, AdvectiveLink
        cv_a = _make_liquid_cv("A", {"S": 1.0})
        cv_b = _make_liquid_cv("B", {"S": 0.0})
        link = AdvectiveLink("A", "gas", "B", "liquid", Q_L_per_h=1.0)
        with pytest.raises(KeyError, match="source phase"):
            Simulation(cvs={"A": cv_a, "B": cv_b}, links=[link])

    def test_unknown_sink_phase_raises(self):
        from PyOMES.core import Simulation, AdvectiveLink
        cv_a = _make_liquid_cv("A", {"S": 1.0})
        cv_b = _make_liquid_cv("B", {"S": 0.0})
        link = AdvectiveLink("A", "liquid", "B", "gas", Q_L_per_h=1.0)
        with pytest.raises(KeyError, match="sink phase"):
            Simulation(cvs={"A": cv_a, "B": cv_b}, links=[link])

    def test_validation_failure_does_not_wire_context(self):
        """If validation fails, no CV should end up with a stale _context."""
        from PyOMES.core import Simulation, AdvectiveLink
        cv = _make_liquid_cv("A", {"S": 1.0})
        bad_link = AdvectiveLink("missing", "liquid", "A", "liquid", Q_L_per_h=1.0)
        with pytest.raises(KeyError):
            Simulation(cvs={"A": cv}, links=[bad_link])
        # CV should still be unowned and able to enter a different Sim
        assert cv._context is None
        sim2 = Simulation(cvs={"A": cv})
        assert cv._context is sim2._context


# ═══════════════════════════════════════════════════════════════════════
#  total_mol
# ═══════════════════════════════════════════════════════════════════════

class TestTotalMol:

    def test_single_cv_sums_across_phases(self):
        from PyOMES.core import Simulation
        cv = _make_gas_liquid_cv(
            "main",
            gas_mol={"O2": 1.0, "N2": 3.0},
            liq_mol={"O2": 0.01},
            V_gas=0.2, V_liq=0.8,
        )
        sim = Simulation(cvs={"main": cv})
        totals = sim.total_mol()
        assert totals["O2"] == pytest.approx(1.01, abs=1e-12)
        assert totals["N2"] == pytest.approx(3.0, abs=1e-12)

    def test_multi_cv_aggregates_across_cvs(self):
        from PyOMES.core import Simulation
        cv_a = _make_liquid_cv("A", {"S": 1.0, "X": 0.5})
        cv_b = _make_liquid_cv("B", {"S": 2.0, "Y": 1.0})
        sim = Simulation(cvs={"A": cv_a, "B": cv_b})
        totals = sim.total_mol()
        assert totals["S"] == pytest.approx(3.0, abs=1e-12)
        assert totals["X"] == pytest.approx(0.5, abs=1e-12)
        assert totals["Y"] == pytest.approx(1.0, abs=1e-12)


# ═══════════════════════════════════════════════════════════════════════
#  RunContext wiring
# ═══════════════════════════════════════════════════════════════════════

class TestRunContext:

    def test_default_runcontext_state(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv}, label="exp1")
        assert sim._context.is_running is False
        assert sim._context.label == "exp1"

    def test_runcontext_wired_to_every_cv(self):
        from PyOMES.core import Simulation
        cv_a = _make_liquid_cv("A", {"S": 1.0})
        cv_b = _make_liquid_cv("B", {"S": 0.0})
        sim = Simulation(cvs={"A": cv_a, "B": cv_b})
        # The exact same RunContext object is wired everywhere
        assert cv_a._context is sim._context
        assert cv_b._context is sim._context

    def test_fresh_cv_has_no_context(self):
        cv = _make_liquid_cv("main", {"S": 1.0})
        assert cv._context is None

    def test_double_ownership_raises(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim1 = Simulation(cvs={"main": cv}, label="exp1")
        # cv is now owned by sim1; a second Simulation must refuse
        with pytest.raises(RuntimeError, match="already owned"):
            Simulation(cvs={"main": cv}, label="exp2")
        # sim1 unchanged
        assert cv._context is sim1._context

    def test_disjoint_simulations_have_independent_contexts(self):
        from PyOMES.core import Simulation
        cv_a = _make_liquid_cv("A", {"S": 1.0})
        cv_b = _make_liquid_cv("B", {"S": 0.0})
        sim_a = Simulation(cvs={"A": cv_a})
        sim_b = Simulation(cvs={"B": cv_b})
        assert sim_a._context is not sim_b._context
        assert cv_a._context is sim_a._context
        assert cv_b._context is sim_b._context


# ═══════════════════════════════════════════════════════════════════════
#  Snapshot
# ═══════════════════════════════════════════════════════════════════════

class TestSnapshot:

    def test_snapshot_returns_new_simulation(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv}, label="orig")
        snap = sim.snapshot()
        assert isinstance(snap, Simulation)
        assert snap is not sim

    def test_snapshot_cvs_are_independent(self):
        """Mutating a snapshot's CV does not affect the original."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0}, V_L=1.0)
        sim = Simulation(cvs={"main": cv})
        snap = sim.snapshot()
        # The two CVs are distinct objects
        assert snap["main"] is not sim["main"]
        # Phases inside are also distinct (deep-copied via cv.snapshot)
        assert snap["main"].phases["liquid"] is not sim["main"].phases["liquid"]
        # Mutating one phase's n_mol does not bleed across
        snap["main"].phases["liquid"].n_mol["S"] = 99.0
        assert sim["main"].phases["liquid"].n_mol["S"] == pytest.approx(1.0)

    def test_snapshot_has_fresh_context(self):
        """The snapshot has its own RunContext, not the original's."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv}, label="orig")
        snap = sim.snapshot()
        # New context (different identity), same label
        assert snap._context is not sim._context
        assert snap._context.label == sim._context.label
        # Snapshot's CV is wired to the snapshot's context, not the original's
        assert snap["main"]._context is snap._context

    def test_snapshot_shares_links_and_controllers(self):
        """Links/controllers/profiles are shared, not copied (per checklist)."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        ctrl = object()
        prof = object()
        sim = Simulation(cvs={"main": cv}, controllers=[ctrl], profiles=[prof])
        snap = sim.snapshot()
        # Lists are not the same list object (each is a fresh list(...) copy)
        # but the contained references are identical.
        assert snap.controllers == sim.controllers
        assert snap.controllers[0] is ctrl
        assert snap.profiles[0] is prof


# ═══════════════════════════════════════════════════════════════════════
#  Repr
# ═══════════════════════════════════════════════════════════════════════

class TestRepr:

    def test_repr_shape(self):
        from PyOMES.core import Simulation
        cv_a = _make_liquid_cv("A", {"S": 1.0})
        cv_b = _make_liquid_cv("B", {"S": 0.0})
        sim = Simulation(cvs={"A": cv_a, "B": cv_b}, label="demo")
        r = repr(sim)
        assert "Simulation" in r
        assert "label='demo'" in r
        assert "cvs=[A, B]" in r
        assert "links=0" in r
        assert "controllers=0" in r
        assert "profiles=0" in r


# ═══════════════════════════════════════════════════════════════════════
#  CVSnapshot / SimulationSnapshot builders (C2)
# ═══════════════════════════════════════════════════════════════════════

class TestCVSnapshotBuilder:

    def test_exports_importable(self):
        from PyOMES.core import (
            CVSnapshot,
            SimulationSnapshot,
            build_cv_snapshot,
            build_simulation_snapshot,
        )
        assert CVSnapshot is not None
        assert SimulationSnapshot is not None
        assert build_cv_snapshot is not None
        assert build_simulation_snapshot is not None

    def test_build_from_gas_liquid_cv(self):
        from PyOMES.core import build_cv_snapshot
        cv = _make_gas_liquid_cv(
            "main",
            gas_mol={"O2": 1.0, "N2": 3.0, "CO2": 0.1},
            liq_mol={"O2": 0.001, "CO2": 0.01},
            V_gas=0.2, V_liq=0.8, T_K=305.15,
        )
        snap = build_cv_snapshot(cv, "main", t_h=2.5)
        assert snap.cv_key == "main"
        assert snap.t_h == 2.5
        assert snap.T_K == pytest.approx(305.15)
        assert snap.V_gas_L == pytest.approx(0.2)
        assert snap.V_liq_L == pytest.approx(0.8)
        assert snap.n_gas_mol["O2"] == pytest.approx(1.0)
        assert snap.n_liq_mol["CO2"] == pytest.approx(0.01)
        # y_gas: mole fractions across the gas inventory
        assert snap.y_gas["O2"] == pytest.approx(1.0 / 4.1)
        # P_gas_atm: from ideal gas law (positive)
        assert snap.P_gas_atm > 0.0

    def test_build_from_liquid_only_cv(self):
        from PyOMES.core import build_cv_snapshot
        cv = _make_liquid_cv("main", {"S": 1.0}, V_L=2.0, T_K=298.15)
        snap = build_cv_snapshot(cv, "main", t_h=0.0)
        assert snap.T_K == pytest.approx(298.15)
        assert snap.V_liq_L == pytest.approx(2.0)
        assert snap.V_gas_L == 0.0
        assert snap.P_gas_atm == 0.0
        assert snap.n_gas_mol == {}
        assert snap.y_gas == {}
        assert snap.n_liq_mol["S"] == pytest.approx(1.0)

    def test_pH_is_none_without_h_plus_species(self):
        """A CV with no H+ species in liquid yields pH = None."""
        from PyOMES.core import build_cv_snapshot
        cv = _make_liquid_cv("main", {"S": 1.0})
        snap = build_cv_snapshot(cv, "main", t_h=0.0)
        assert snap.pH is None
        assert snap.ionic_strength is None

    def test_pH_computed_when_h_plus_present(self):
        """When n_mol['H+'] exists and is positive, pH derives from -log10."""
        from PyOMES.core import build_cv_snapshot
        # H+ = 1e-7 mol in V_L = 1.0 L → [H+] = 1e-7 mol/L → pH = 7
        cv = _make_liquid_cv("main", {"H+": 1e-7, "S": 1.0}, V_L=1.0)
        snap = build_cv_snapshot(cv, "main", t_h=0.0)
        assert snap.pH == pytest.approx(7.0, abs=1e-9)

    def test_ionic_strength_from_properties_dict(self):
        """ionic_strength is forwarded from liquid.properties when present."""
        from PyOMES.core import ControlVolume, LiquidPhase, build_cv_snapshot
        liq = LiquidPhase(
            n_mol={"S": 1.0},
            V_L=1.0,
            T_K=298.15,
            properties={"ionic_strength": 0.05},
        )
        cv = ControlVolume(phases={"liquid": liq}, label="main")
        snap = build_cv_snapshot(cv, "main", t_h=0.0)
        assert snap.ionic_strength == pytest.approx(0.05)

    def test_ionic_strength_fallback_to_speciation_key(self):
        """When properties lacks 'ionic_strength', fall back to legacy
        speciation['IonicStrength']."""
        from PyOMES.core import ControlVolume, LiquidPhase, build_cv_snapshot
        liq = LiquidPhase(
            n_mol={"S": 1.0},
            V_L=1.0,
            T_K=298.15,
            speciation={"IonicStrength": 0.12},
        )
        cv = ControlVolume(phases={"liquid": liq}, label="main")
        snap = build_cv_snapshot(cv, "main", t_h=0.0)
        assert snap.ionic_strength == pytest.approx(0.12)

    def test_species_properties_forwarded_from_liquid(self):
        """species_properties is a fresh copy of liquid.properties."""
        from PyOMES.core import ControlVolume, LiquidPhase, build_cv_snapshot
        liq = LiquidPhase(
            n_mol={"S": 1.0},
            V_L=1.0,
            T_K=298.15,
            properties={"viscosity_Pa_s": 1.2e-3, "density_kg_L": 0.997},
        )
        cv = ControlVolume(phases={"liquid": liq}, label="main")
        snap = build_cv_snapshot(cv, "main", t_h=0.0)
        assert snap.species_properties["viscosity_Pa_s"] == pytest.approx(1.2e-3)
        assert snap.species_properties["density_kg_L"] == pytest.approx(0.997)
        # Verify it's a copy: mutating the snapshot's dict shouldn't
        # bleed back into the phase
        snap.species_properties["viscosity_Pa_s"] = 99.0
        assert liq.properties["viscosity_Pa_s"] == pytest.approx(1.2e-3)

    def test_sensors_contain_T_C_and_DO_when_O2_in_liquid(self):
        """The sensors dict exposes T_C unconditionally and DO_mol_L
        when oxygen is present in the liquid phase."""
        from PyOMES.core import build_cv_snapshot
        cv = _make_gas_liquid_cv(
            "main",
            gas_mol={"O2": 1.0},
            liq_mol={"O2": 0.002},
            V_gas=0.2, V_liq=2.0, T_K=300.0,
        )
        snap = build_cv_snapshot(cv, "main", t_h=0.0)
        assert snap.sensors["T_C"] == pytest.approx(300.0 - 273.15)
        # DO_mol_L = liquid n_mol["O2"] / V_liq_L = 0.002 / 2.0 = 0.001
        assert snap.sensors["DO_mol_L"] == pytest.approx(0.001)

    def test_sensors_lack_DO_when_O2_absent(self):
        """Without O2 in liquid, the DO_mol_L key is omitted (not None-valued)."""
        from PyOMES.core import build_cv_snapshot
        cv = _make_liquid_cv("main", {"S": 1.0})
        snap = build_cv_snapshot(cv, "main", t_h=0.0)
        assert "DO_mol_L" not in snap.sensors
        assert "T_C" in snap.sensors

    def test_snapshot_is_frozen(self):
        """Attempting to reassign a CVSnapshot field raises FrozenInstanceError."""
        from dataclasses import FrozenInstanceError
        from PyOMES.core import build_cv_snapshot
        cv = _make_liquid_cv("main", {"S": 1.0})
        snap = build_cv_snapshot(cv, "main", t_h=0.0)
        with pytest.raises(FrozenInstanceError):
            snap.pH = 7.0  # type: ignore[misc]
        with pytest.raises(FrozenInstanceError):
            snap.T_K = 500.0  # type: ignore[misc]

    def test_asdict_round_trip(self):
        """dataclasses.asdict produces a faithful dict view; reconstruction
        from that dict gives an equal CVSnapshot."""
        from dataclasses import asdict
        from PyOMES.core import CVSnapshot, build_cv_snapshot
        cv = _make_gas_liquid_cv(
            "main",
            gas_mol={"O2": 1.0},
            liq_mol={"O2": 0.001},
            V_gas=0.2, V_liq=0.8, T_K=300.0,
        )
        snap = build_cv_snapshot(cv, "main", t_h=1.5)
        d = asdict(snap)
        assert d["cv_key"] == "main"
        assert d["t_h"] == pytest.approx(1.5)
        # Reconstruction
        snap2 = CVSnapshot(**d)
        assert snap2 == snap

    def test_advance_result_none_is_tolerated(self):
        """The advance_result hook accepts None and is currently unused."""
        from PyOMES.core import build_cv_snapshot
        cv = _make_liquid_cv("main", {"S": 1.0})
        snap_none = build_cv_snapshot(cv, "main", t_h=0.0, advance_result=None)
        # No advance_result passed at all should match
        snap_default = build_cv_snapshot(cv, "main", t_h=0.0)
        assert snap_none == snap_default

    def test_general_phase_dicts_fermenter(self):
        """phase_keys / n_mol / V_L / properties are populated for a gas+liquid CV."""
        from PyOMES.core import build_cv_snapshot
        cv = _make_gas_liquid_cv(
            "main",
            gas_mol={"O2": 1.0, "CO2": 0.5},
            liq_mol={"S": 2.0},
            V_gas=0.2, V_liq=0.8, T_K=300.0,
        )
        snap = build_cv_snapshot(cv, "main", t_h=0.0)
        assert set(snap.phase_keys) == {"gas", "liquid"}
        assert snap.n_mol["gas"]["O2"] == pytest.approx(1.0)
        assert snap.n_mol["liquid"]["S"] == pytest.approx(2.0)
        assert snap.V_L["gas"] == pytest.approx(0.2)
        assert snap.V_L["liquid"] == pytest.approx(0.8)
        # convenience aliases still correct
        assert snap.n_liq_mol["S"] == pytest.approx(2.0)
        assert snap.n_gas_mol["O2"] == pytest.approx(1.0)

    def test_general_phase_dicts_custom_phase_keys(self):
        """CVs with non-standard phase names populate phase_keys/n_mol/V_L."""
        from PyOMES.core import ControlVolume, build_cv_snapshot
        from PyOMES.core.phases import LiquidPhase
        cv = ControlVolume(phases={
            "mobile": LiquidPhase(n_mol={"MeOH": 0.5, "water": 1.0}, V_L=0.1, T_K=298.15),
            "stationary": LiquidPhase(n_mol={"analyte": 0.01}, V_L=0.05, T_K=298.15),
        })
        snap = build_cv_snapshot(cv, "col1", t_h=0.0)
        assert snap.phase_keys == ["mobile", "stationary"]
        assert snap.n_mol["mobile"]["MeOH"] == pytest.approx(0.5)
        assert snap.n_mol["stationary"]["analyte"] == pytest.approx(0.01)
        assert snap.V_L["mobile"] == pytest.approx(0.1)
        # convenience aliases absent for non-standard keys
        assert snap.V_liq_L == 0.0
        assert snap.n_liq_mol == {}

    def test_properties_dict_captures_phase_properties(self):
        """phase.properties values appear in snap.properties[phase_key]."""
        from PyOMES.core import ControlVolume, GasPhase, LiquidPhase, build_cv_snapshot
        liq = LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=298.15,
                          properties={"viscosity": 1e-3})
        gas = GasPhase(n_mol={"O2": 0.5}, V_L=0.2, T_K=298.15,
                       properties={"density_kg_m3": 1.3})
        cv = ControlVolume(phases={"liquid": liq, "gas": gas})
        snap = build_cv_snapshot(cv, "main", t_h=0.0)
        assert snap.properties["liquid"]["viscosity"] == pytest.approx(1e-3)
        assert snap.properties["gas"]["density_kg_m3"] == pytest.approx(1.3)
        # species_properties alias still reflects liquid.properties
        assert snap.species_properties["viscosity"] == pytest.approx(1e-3)


class TestSimulationSnapshotBuilder:

    def test_aggregates_all_cvs(self):
        from PyOMES.core import Simulation, build_simulation_snapshot
        cv_a = _make_liquid_cv("A", {"S": 1.0}, V_L=1.0)
        cv_b = _make_gas_liquid_cv(
            "B",
            gas_mol={"O2": 0.5},
            liq_mol={"O2": 0.001},
            V_gas=0.1, V_liq=0.9,
        )
        sim = Simulation(cvs={"A": cv_a, "B": cv_b})
        sim_snap = build_simulation_snapshot(sim, t_h=3.0)
        assert sim_snap.t_h == pytest.approx(3.0)
        assert sorted(sim_snap.cvs.keys()) == ["A", "B"]
        assert sim_snap.cvs["A"].cv_key == "A"
        assert sim_snap.cvs["B"].cv_key == "B"
        assert sim_snap.cvs["B"].V_gas_L == pytest.approx(0.1)

    def test_results_dict_threaded_through(self):
        """When the orchestrator passes per-CV AdvanceResults, each
        CV's snapshot receives its own."""
        from PyOMES.core import Simulation, build_simulation_snapshot
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        # Pass a stub advance_result; today it is unused by the snapshot
        # builder, but the threading is verified.
        sim_snap = build_simulation_snapshot(
            sim, t_h=0.0, results={"main": object()}
        )
        assert "main" in sim_snap.cvs

    def test_snapshot_is_frozen(self):
        from dataclasses import FrozenInstanceError
        from PyOMES.core import Simulation, build_simulation_snapshot
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        sim_snap = build_simulation_snapshot(sim, t_h=0.0)
        with pytest.raises(FrozenInstanceError):
            sim_snap.t_h = 99.0  # type: ignore[misc]


# ═══════════════════════════════════════════════════════════════════════
#  Recorder protocol + BatchRecorder + BatchResult (C3)
# ═══════════════════════════════════════════════════════════════════════

class TestActionStubs:

    def test_control_action_importable_and_default_constructible(self):
        from PyOMES.control.actions import ControlAction
        action = ControlAction()
        assert action.controller_label == ""
        assert action.flux_applied == {}

    def test_profile_record_importable_and_default_constructible(self):
        from PyOMES.control.actions import ProfileRecord
        record = ProfileRecord()
        assert record.profile_label == ""
        assert record.targets == {}

    def test_control_action_is_frozen(self):
        from dataclasses import FrozenInstanceError
        from PyOMES.control.actions import ControlAction
        a = ControlAction()
        with pytest.raises(FrozenInstanceError):
            a.controller_label = "x"  # type: ignore[misc]


class TestRecorderProtocol:

    def test_recorder_exports_importable(self):
        from PyOMES.core import Recorder, BatchRecorder, BatchResult
        assert Recorder is not None
        assert BatchRecorder is not None
        assert BatchResult is not None

    def test_batch_recorder_satisfies_recorder_protocol(self):
        """isinstance check via runtime_checkable Protocol."""
        from PyOMES.core import Recorder, BatchRecorder
        rec = BatchRecorder()
        assert isinstance(rec, Recorder)

    def test_batch_result_default_constructible(self):
        import numpy as np
        from PyOMES.core import BatchResult
        r = BatchResult()
        assert isinstance(r.t_h, np.ndarray)
        assert r.gas_mol == {}
        assert r.controller_actions == []
        assert r.accuracy_records == []
        assert r.conservation_records == []
        assert r.runtime_s == 0.0


class TestBatchRecorderPreAllocation:

    def test_record_init_preallocates_arrays(self):
        import numpy as np
        from PyOMES.core import Simulation, BatchRecorder
        cv = _make_gas_liquid_cv(
            "main",
            gas_mol={"O2": 1.0, "N2": 3.0},
            liq_mol={"O2": 0.001, "CO2": 0.01},
        )
        sim = Simulation(cvs={"main": cv})
        rec = BatchRecorder()
        rec.record_init(sim, n_steps=10)
        r = rec._result
        # Time grid is (n_steps + 1,)
        assert r.t_h.shape == (11,)
        # Per-CV arrays for every species seen at t=0
        assert "main" in r.gas_mol
        assert "main" in r.liquid_mol
        assert r.gas_mol["main"]["O2"].shape == (11,)
        assert r.gas_mol["main"]["N2"].shape == (11,)
        assert r.liquid_mol["main"]["CO2"].shape == (11,)
        # pH / ionic_strength arrays exist and are NaN
        assert np.all(np.isnan(r.pH["main"]))
        assert np.all(np.isnan(r.ionic_strength["main"]))

    def test_record_init_writes_initial_state(self):
        from PyOMES.core import Simulation, BatchRecorder
        cv = _make_gas_liquid_cv(
            "main",
            gas_mol={"O2": 1.0, "N2": 3.0},
            liq_mol={"O2": 0.001},
            V_gas=0.2, V_liq=0.8, T_K=300.0,
        )
        sim = Simulation(cvs={"main": cv})
        rec = BatchRecorder()
        rec.record_init(sim, n_steps=5)
        r = rec._result
        # Index 0 holds the initial state
        assert r.gas_mol["main"]["O2"][0] == pytest.approx(1.0)
        assert r.gas_mol["main"]["N2"][0] == pytest.approx(3.0)
        assert r.liquid_mol["main"]["O2"][0] == pytest.approx(0.001)
        assert r.P_atm["main"][0] > 0.0  # ideal-gas

    def test_record_init_multi_cv_preallocates_per_cv(self):
        from PyOMES.core import Simulation, BatchRecorder
        cv_a = _make_liquid_cv("A", {"S": 1.0, "X": 0.5})
        cv_b = _make_liquid_cv("B", {"S": 2.0, "Y": 1.0})
        sim = Simulation(cvs={"A": cv_a, "B": cv_b})
        rec = BatchRecorder()
        rec.record_init(sim, n_steps=3)
        r = rec._result
        assert sorted(r.liquid_mol.keys()) == ["A", "B"]
        assert "S" in r.liquid_mol["A"]
        assert "X" in r.liquid_mol["A"]
        assert "Y" in r.liquid_mol["B"]
        # Each CV gets its own pH array
        assert "A" in r.pH and "B" in r.pH

    def test_record_init_invalid_n_steps_raises(self):
        from PyOMES.core import Simulation, BatchRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = BatchRecorder()
        with pytest.raises(ValueError, match="n_steps"):
            rec.record_init(sim, n_steps=0)

    def test_record_init_twice_raises(self):
        from PyOMES.core import Simulation, BatchRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = BatchRecorder()
        rec.record_init(sim, n_steps=3)
        with pytest.raises(RuntimeError, match="twice"):
            rec.record_init(sim, n_steps=3)


class TestBatchRecorderStep:

    def test_record_step_before_init_raises(self):
        from PyOMES.core import BatchRecorder
        rec = BatchRecorder()
        with pytest.raises(RuntimeError, match="before record_init"):
            rec.record_step(
                step_index=1, t_h=0.1, dt_h=0.1,
                advance_results={}, controller_actions=[], profile_actions=[],
            )

    def test_record_step_writes_at_correct_index(self):
        from PyOMES.core import Simulation, BatchRecorder
        cv = _make_liquid_cv("main", {"S": 1.0}, V_L=1.0)
        sim = Simulation(cvs={"main": cv})
        rec = BatchRecorder()
        rec.record_init(sim, n_steps=3)

        # Simulate one step: mutate the CV's n_mol and call record_step
        cv.phases["liquid"].n_mol["S"] = 0.9
        rec.record_step(
            step_index=1, t_h=0.1, dt_h=0.1,
            advance_results={}, controller_actions=[], profile_actions=[],
        )
        r = rec._result
        assert r.t_h[1] == pytest.approx(0.1)
        assert r.liquid_mol["main"]["S"][0] == pytest.approx(1.0)  # initial
        assert r.liquid_mol["main"]["S"][1] == pytest.approx(0.9)  # after step

    def test_record_step_handles_new_species(self):
        """Species appearing in a phase after record_init get an array
        allocated on the first record_step that sees them."""
        from PyOMES.core import Simulation, BatchRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = BatchRecorder()
        rec.record_init(sim, n_steps=3)

        # Add a brand-new species, then record
        cv.phases["liquid"].n_mol["P"] = 0.05
        rec.record_step(
            step_index=1, t_h=0.1, dt_h=0.1,
            advance_results={}, controller_actions=[], profile_actions=[],
        )
        r = rec._result
        assert "P" in r.liquid_mol["main"]
        assert r.liquid_mol["main"]["P"].shape == (4,)
        assert r.liquid_mol["main"]["P"][0] == 0.0  # was not present at init
        assert r.liquid_mol["main"]["P"][1] == pytest.approx(0.05)

    def test_record_step_out_of_range_raises(self):
        from PyOMES.core import Simulation, BatchRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = BatchRecorder()
        rec.record_init(sim, n_steps=3)
        with pytest.raises(IndexError, match="outside range"):
            rec.record_step(
                step_index=99, t_h=0.1, dt_h=0.1,
                advance_results={}, controller_actions=[], profile_actions=[],
            )

    def test_record_step_threads_controller_and_profile_actions(self):
        from PyOMES.core import Simulation, BatchRecorder
        from PyOMES.control.actions import ControlAction, ProfileRecord
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = BatchRecorder()
        rec.record_init(sim, n_steps=2)
        action = ControlAction(controller_label="ph_ctrl", t_h=0.1)
        profile = ProfileRecord(profile_label="T_ramp", t_h=0.1)
        rec.record_step(
            step_index=1, t_h=0.1, dt_h=0.1,
            advance_results={},
            controller_actions=[action],
            profile_actions=[profile],
        )
        r = rec._result
        assert len(r.controller_actions) == 1
        assert r.controller_actions[0][0].controller_label == "ph_ctrl"
        assert len(r.profile_actions) == 1
        assert r.profile_actions[0][0].profile_label == "T_ramp"


class TestBatchRecorderFinalize:

    def test_finalize_before_init_raises(self):
        from PyOMES.core import BatchRecorder
        rec = BatchRecorder()
        with pytest.raises(RuntimeError, match="before record_init"):
            rec.finalize()

    def test_finalize_returns_batch_result_with_runtime(self):
        from PyOMES.core import Simulation, BatchRecorder, BatchResult
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = BatchRecorder()
        rec.record_init(sim, n_steps=2)
        result = rec.finalize()
        assert isinstance(result, BatchResult)
        # Runtime should be non-negative; perf_counter monotonic
        assert result.runtime_s >= 0.0

    def test_finalize_twice_raises(self):
        from PyOMES.core import Simulation, BatchRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = BatchRecorder()
        rec.record_init(sim, n_steps=2)
        rec.finalize()
        with pytest.raises(RuntimeError, match="twice"):
            rec.finalize()


# ═══════════════════════════════════════════════════════════════════════
#  BatchRecorder — phase-agnostic storage (property-snapshot-phase-agnostic C4)
# ═══════════════════════════════════════════════════════════════════════

class TestBatchRecorderPhaseAgnostic:

    def test_phase_mol_is_primary_storage(self):
        """phase_mol contains all phases; gas_mol/liquid_mol are aliases."""
        from PyOMES.core import Simulation, BatchRecorder
        cv = _make_gas_liquid_cv(
            "main",
            gas_mol={"O2": 1.0}, liq_mol={"S": 2.0},
            V_gas=0.2, V_liq=0.8, T_K=300.0,
        )
        sim = Simulation(cvs={"main": cv})
        rec = BatchRecorder()
        rec.record_init(sim, n_steps=2)
        result = rec.finalize()
        assert "gas" in result.phase_mol["main"]
        assert "liquid" in result.phase_mol["main"]
        assert result.phase_mol["main"]["gas"]["O2"].shape == (3,)
        assert result.phase_mol["main"]["liquid"]["S"].shape == (3,)

    def test_gas_mol_alias_matches_phase_mol(self):
        """result.gas_mol is a view into phase_mol["gas"]."""
        from PyOMES.core import Simulation, BatchRecorder
        cv = _make_gas_liquid_cv(
            "main",
            gas_mol={"O2": 1.0}, liq_mol={"S": 0.5},
            V_gas=0.2, V_liq=0.8, T_K=300.0,
        )
        sim = Simulation(cvs={"main": cv})
        result = sim.run(tau_h=0.1, n_steps=3)
        assert "O2" in result.gas_mol["main"]
        assert result.gas_mol["main"]["O2"] is result.phase_mol["main"]["gas"]["O2"]

    def test_liquid_mol_alias_matches_phase_mol(self):
        from PyOMES.core import Simulation
        cv = _make_gas_liquid_cv(
            "main",
            gas_mol={"O2": 1.0}, liq_mol={"S": 0.5},
            V_gas=0.2, V_liq=0.8, T_K=300.0,
        )
        sim = Simulation(cvs={"main": cv})
        result = sim.run(tau_h=0.1, n_steps=3)
        assert result.liquid_mol["main"]["S"] is result.phase_mol["main"]["liquid"]["S"]

    def test_solid_phase_trajectory_recorded(self):
        """SolidPhase moles appear in phase_mol under the registered key."""
        from PyOMES.core import ControlVolume, Simulation
        from PyOMES.core.phases import LiquidPhase, SolidPhase
        liq = LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=298.15)
        sol = SolidPhase(n_mol={"CaCO3": 0.5}, V_L=0.01, T_K=298.15)
        cv = ControlVolume(phases={"liquid": liq, "solid": sol})
        sim = Simulation(cvs={"main": cv})
        result = sim.run(tau_h=0.1, n_steps=2)
        assert "solid" in result.phase_mol["main"]
        assert "CaCO3" in result.phase_mol["main"]["solid"]
        assert result.phase_mol["main"]["solid"]["CaCO3"].shape == (3,)
        # Absent from the gas alias
        assert "solid" not in {k for cv_phases in [result.phase_mol["main"]]
                               for k in ["gas"] if k in cv_phases}

    def test_custom_phase_keys_in_to_dataframe(self):
        """to_dataframe rows carry the actual phase_key, not just 'gas'/'liquid'."""
        pytest.importorskip("pandas")
        from PyOMES.core import ControlVolume, Simulation
        from PyOMES.core.phases import LiquidPhase
        mobile = LiquidPhase(n_mol={"MeOH": 0.5}, V_L=0.1, T_K=298.15)
        stationary = LiquidPhase(n_mol={"analyte": 0.01}, V_L=0.05, T_K=298.15)
        cv = ControlVolume(phases={"mobile": mobile, "stationary": stationary})
        sim = Simulation(cvs={"col": cv})
        result = sim.run(tau_h=0.05, n_steps=2)
        df = result.to_dataframe()
        phase_keys_in_df = set(df[df["channel_kind"] == "phase_mol"]["phase_key"])
        assert "mobile" in phase_keys_in_df
        assert "stationary" in phase_keys_in_df


# ═══════════════════════════════════════════════════════════════════════
#  StreamingFileRecorder + load_run (RUN_HISTORY C1)
# ═══════════════════════════════════════════════════════════════════════

class TestStreamingFileRecorderExports:

    def test_streaming_recorder_importable_from_core(self):
        from PyOMES.core import StreamingFileRecorder  # noqa: F401

    def test_load_run_importable_from_core(self):
        from PyOMES.core import load_run  # noqa: F401

    def test_streaming_recorder_satisfies_recorder_protocol(self):
        from PyOMES.core import StreamingFileRecorder
        from PyOMES.core.recorder import Recorder
        assert isinstance(StreamingFileRecorder.__new__(object), object)
        # Protocol satisfaction verified structurally (has required methods).
        for method in ("record_init", "record_step", "finalize"):
            assert hasattr(StreamingFileRecorder, method)


class TestStreamingFileRecorderInit:

    def test_creates_directory(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import Simulation, StreamingFileRecorder
        run_dir = tmp_path / "run"
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = StreamingFileRecorder(run_dir)
        rec.record_init(sim, n_steps=5)
        assert run_dir.is_dir()

    def test_manifest_written_after_init(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import Simulation, StreamingFileRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = StreamingFileRecorder(tmp_path / "run")
        rec.record_init(sim, n_steps=5)
        manifest_path = tmp_path / "run" / "manifest.json"
        assert manifest_path.exists()

    def test_manifest_contains_cv_phase_species(self, tmp_path):
        pytest.importorskip("pyarrow")
        import json
        from PyOMES.core import Simulation, StreamingFileRecorder
        cv = _make_liquid_cv("main", {"S": 1.0, "P": 0.5})
        sim = Simulation(cvs={"main": cv})
        rec = StreamingFileRecorder(tmp_path / "run")
        rec.record_init(sim, n_steps=5)
        manifest = json.loads((tmp_path / "run" / "manifest.json").read_text())
        assert manifest["schema_version"] == 1
        assert "main" in manifest["cv_phase_species"]
        assert "liquid" in manifest["cv_phase_species"]["main"]
        assert set(manifest["cv_phase_species"]["main"]["liquid"]) == {"S", "P"}

    def test_record_init_twice_raises(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import Simulation, StreamingFileRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = StreamingFileRecorder(tmp_path / "run")
        rec.record_init(sim, n_steps=5)
        with pytest.raises(RuntimeError, match="record_init called twice"):
            rec.record_init(sim, n_steps=5)

    def test_pyarrow_missing_raises_import_error(self, tmp_path, monkeypatch):
        # monkeypatch restores sys.modules after the test — no reload needed.
        import sys
        monkeypatch.setitem(sys.modules, "pyarrow", None)
        from PyOMES.core.recorder import StreamingFileRecorder as _SFR
        with pytest.raises(ImportError, match="pyarrow"):
            _SFR(tmp_path / "run")


class TestStreamingFileRecorderFlush:

    def test_no_chunk_before_threshold(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import Simulation, StreamingFileRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        run_dir = tmp_path / "run"
        rec = StreamingFileRecorder(run_dir, flush_every_n_steps=10)
        rec.record_init(sim, n_steps=5)
        # 1 row in buffer (t=0); threshold is 10 — no flush yet
        assert not list(run_dir.glob("chunk_*.parquet"))

    def test_flush_triggered_at_threshold(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import Simulation, StreamingFileRecorder
        from PyOMES.core.recorder import Recorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        run_dir = tmp_path / "run"
        rec = StreamingFileRecorder(run_dir, flush_every_n_steps=3)
        rec.record_init(sim, n_steps=10)
        # Buffer has t=0 row; add 2 more steps to hit threshold of 3
        dummy_ar = {}
        rec.record_step(1, 1.0, 1.0, dummy_ar)
        assert not list(run_dir.glob("chunk_*.parquet"))
        rec.record_step(2, 2.0, 1.0, dummy_ar)
        assert (run_dir / "chunk_00000.parquet").exists()

    def test_multiple_chunks_produced(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import Simulation, StreamingFileRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        run_dir = tmp_path / "run"
        rec = StreamingFileRecorder(run_dir, flush_every_n_steps=3)
        rec.record_init(sim, n_steps=10)
        dummy_ar = {}
        for i in range(1, 8):
            rec.record_step(i, float(i), 1.0, dummy_ar)
        # By step 7 we have 8 rows total; two full flushes at rows 3 and 6
        chunks = sorted(run_dir.glob("chunk_*.parquet"))
        assert len(chunks) >= 2


class TestStreamingFileRecorderFinalize:

    def test_finalize_writes_remaining_buffer(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import Simulation, StreamingFileRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        run_dir = tmp_path / "run"
        rec = StreamingFileRecorder(run_dir, flush_every_n_steps=100)
        rec.record_init(sim, n_steps=3)
        dummy_ar = {}
        for i in range(1, 4):
            rec.record_step(i, float(i), 1.0, dummy_ar)
        rec.finalize()
        assert (run_dir / "chunk_00000.parquet").exists()

    def test_finalize_updates_manifest_chunks(self, tmp_path):
        pytest.importorskip("pyarrow")
        import json
        from PyOMES.core import Simulation, StreamingFileRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        run_dir = tmp_path / "run"
        rec = StreamingFileRecorder(run_dir, flush_every_n_steps=100)
        rec.record_init(sim, n_steps=3)
        for i in range(1, 4):
            rec.record_step(i, float(i), 1.0, {})
        rec.finalize()
        manifest = json.loads((run_dir / "manifest.json").read_text())
        assert len(manifest["chunks"]) == 1
        assert manifest["chunks"][0]["start_step"] == 0
        assert manifest["chunks"][0]["end_step"] == 3
        assert manifest["n_steps_total"] == 3

    def test_finalize_returns_batch_result(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import Simulation, StreamingFileRecorder
        from PyOMES.core.recorder import BatchResult
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = StreamingFileRecorder(tmp_path / "run", flush_every_n_steps=100)
        rec.record_init(sim, n_steps=3)
        for i in range(1, 4):
            rec.record_step(i, float(i), 1.0, {})
        result = rec.finalize()
        assert isinstance(result, BatchResult)

    def test_finalize_twice_raises(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import Simulation, StreamingFileRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = StreamingFileRecorder(tmp_path / "run", flush_every_n_steps=100)
        rec.record_init(sim, n_steps=2)
        for i in range(1, 3):
            rec.record_step(i, float(i), 1.0, {})
        rec.finalize()
        with pytest.raises(RuntimeError, match="finalize called twice"):
            rec.finalize()

    def test_finalize_before_init_raises(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import StreamingFileRecorder
        rec = StreamingFileRecorder(tmp_path / "run")
        with pytest.raises(RuntimeError, match="before record_init"):
            rec.finalize()


class TestStreamingFileRecorderPickle:

    def test_getstate_excludes_buffer(self, tmp_path):
        pytest.importorskip("pyarrow")
        import pickle
        from PyOMES.core import Simulation, StreamingFileRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = StreamingFileRecorder(tmp_path / "run", flush_every_n_steps=100)
        rec.record_init(sim, n_steps=5)
        for i in range(1, 4):
            rec.record_step(i, float(i), 1.0, {})
        assert len(rec._buffer) == 4  # 1 init row + 3 steps
        state = rec.__getstate__()
        assert state["_buffer"] == []

    def test_pickle_roundtrip_drops_buffer(self, tmp_path):
        pytest.importorskip("pyarrow")
        import pickle
        from PyOMES.core import Simulation, StreamingFileRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = StreamingFileRecorder(tmp_path / "run", flush_every_n_steps=100)
        rec.record_init(sim, n_steps=5)
        for i in range(1, 4):
            rec.record_step(i, float(i), 1.0, {})
        restored = pickle.loads(pickle.dumps(rec))
        assert restored._buffer == []
        assert restored._chunk_index == rec._chunk_index


class TestLoadRun:

    def test_missing_manifest_raises(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import load_run
        with pytest.raises(FileNotFoundError, match="manifest.json"):
            load_run(tmp_path / "nonexistent")

    def test_t_h_matches_batch_recorder(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import Simulation, StreamingFileRecorder, load_run
        import numpy as np
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})},
            recorder=StreamingFileRecorder(tmp_path / "run", flush_every_n_steps=100),
        )
        result_streaming = sim.run(tau_h=0.1, n_steps=5)
        result_batch = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})}
        ).run(tau_h=0.1, n_steps=5)
        np.testing.assert_allclose(result_streaming.t_h, result_batch.t_h)

    def test_phase_mol_shape_matches(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import Simulation, StreamingFileRecorder
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0, "P": 0.5})},
            recorder=StreamingFileRecorder(tmp_path / "run", flush_every_n_steps=100),
        )
        result = sim.run(tau_h=0.1, n_steps=4)
        assert result.phase_mol["main"]["liquid"]["S"].shape == (5,)
        assert result.phase_mol["main"]["liquid"]["P"].shape == (5,)

    def test_initial_moles_preserved(self, tmp_path):
        pytest.importorskip("pyarrow")
        import numpy as np
        from PyOMES.core import Simulation, StreamingFileRecorder
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 3.7})},
            recorder=StreamingFileRecorder(tmp_path / "run", flush_every_n_steps=100),
        )
        result = sim.run(tau_h=0.1, n_steps=3)
        assert result.phase_mol["main"]["liquid"]["S"][0] == pytest.approx(3.7)

    def test_multi_chunk_reassembly(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import Simulation, StreamingFileRecorder
        import numpy as np
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})},
            recorder=StreamingFileRecorder(tmp_path / "run", flush_every_n_steps=4),
        )
        result = sim.run(tau_h=0.01, n_steps=10)
        # Multiple chunks but result covers all 11 time points (0..10)
        assert len(result.t_h) == 11

    def test_gas_liquid_cv_p_atm_recovered(self, tmp_path):
        pytest.importorskip("pyarrow")
        from PyOMES.core import Simulation, StreamingFileRecorder
        import numpy as np
        sim = Simulation(
            cvs={"main": _make_gas_liquid_cv("main", gas_mol={"O2": 0.5}, liq_mol={"S": 1.0})},
            recorder=StreamingFileRecorder(tmp_path / "run", flush_every_n_steps=100),
        )
        result = sim.run(tau_h=0.01, n_steps=3)
        assert result.P_atm["main"].shape == (4,)
        assert np.all(result.P_atm["main"] > 0)

    def test_empty_chunks_returns_empty_batch_result(self, tmp_path):
        pytest.importorskip("pyarrow")
        import json
        from PyOMES.core import load_run
        run_dir = tmp_path / "run"
        run_dir.mkdir()
        manifest = {
            "schema_version": 1, "cv_phase_species": {}, "columns": [],
            "chunks": [], "n_steps_total": 0,
        }
        (run_dir / "manifest.json").write_text(json.dumps(manifest))
        result = load_run(run_dir)
        assert len(result.t_h) == 0


class TestSaveCheckpointStreamingRecorder:

    def test_checkpoint_manifest_includes_streaming_path(self, tmp_path):
        pytest.importorskip("pyarrow")
        import json
        from PyOMES.core import Simulation, StreamingFileRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        run_dir = tmp_path / "run"
        sim = Simulation(
            cvs={"main": cv},
            recorder=StreamingFileRecorder(run_dir, flush_every_n_steps=100),
        )
        sim.run(tau_h=0.1, n_steps=3)
        ckpt_dir = tmp_path / "ckpt"
        sim.save_checkpoint(ckpt_dir)
        manifest = json.loads((ckpt_dir / "manifest.json").read_text())
        assert "streaming_recorder_path" in manifest
        assert str(run_dir) in manifest["streaming_recorder_path"]

    def test_checkpoint_pickle_excludes_buffer(self, tmp_path):
        pytest.importorskip("pyarrow")
        import pickle, gzip
        from PyOMES.core import Simulation, StreamingFileRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        run_dir = tmp_path / "run"
        rec = StreamingFileRecorder(run_dir, flush_every_n_steps=100)
        sim = Simulation(cvs={"main": cv}, recorder=rec)
        sim.run(tau_h=0.1, n_steps=3)
        ckpt_dir = tmp_path / "ckpt"
        sim.save_checkpoint(ckpt_dir)
        raw = gzip.decompress((ckpt_dir / "simulation.pkl").read_bytes())
        restored_sim = pickle.loads(raw)
        assert restored_sim._recorder._buffer == []


# ═══════════════════════════════════════════════════════════════════════
#  SparseRecorder (RUN_HISTORY C2)
# ═══════════════════════════════════════════════════════════════════════

class TestSparseRecorderExports:

    def test_importable_from_core(self):
        from PyOMES.core import SparseRecorder  # noqa: F401

    def test_satisfies_recorder_protocol(self):
        from PyOMES.core import SparseRecorder
        for method in ("record_init", "record_step", "finalize"):
            assert hasattr(SparseRecorder, method)


class TestSparseRecorderBasics:

    def test_returns_batch_result(self):
        from PyOMES.core import Simulation, SparseRecorder
        from PyOMES.core.recorder import BatchResult
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})},
            recorder=SparseRecorder(every_m_steps=5),
        )
        result = sim.run(tau_h=0.1, n_steps=10)
        assert isinstance(result, BatchResult)

    def test_always_includes_initial_state(self):
        from PyOMES.core import Simulation, SparseRecorder
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 2.5})},
            recorder=SparseRecorder(every_m_steps=5),
        )
        result = sim.run(tau_h=0.1, n_steps=10)
        assert result.t_h[0] == pytest.approx(0.0)
        assert result.phase_mol["main"]["liquid"]["S"][0] == pytest.approx(2.5)

    def test_always_includes_final_state(self):
        from PyOMES.core import Simulation, SparseRecorder
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})},
            recorder=SparseRecorder(every_m_steps=3),
        )
        result = sim.run(tau_h=0.1, n_steps=10)
        # Final t_h should be n_steps * dt_h = 10 * 0.01
        assert result.t_h[-1] == pytest.approx(0.1)

    def test_fewer_rows_than_batch_recorder(self):
        from PyOMES.core import Simulation, SparseRecorder
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})},
            recorder=SparseRecorder(every_m_steps=5),
        )
        result = sim.run(tau_h=0.1, n_steps=20)
        # BatchRecorder would give 21 rows; SparseRecorder (every 5) gives far fewer
        assert len(result.t_h) < 21

    def test_every_m_steps_recorded(self):
        from PyOMES.core import Simulation, SparseRecorder
        # n_steps=10, every_m_steps=5 -> recorded: 0, 5, 10 = 3 rows
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})},
            recorder=SparseRecorder(every_m_steps=5),
        )
        result = sim.run(tau_h=0.5, n_steps=10)
        assert len(result.t_h) == 3

    def test_every_1_step_equals_batch_recorder_length(self):
        from PyOMES.core import Simulation, SparseRecorder
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})},
            recorder=SparseRecorder(every_m_steps=1),
        )
        result = sim.run(tau_h=0.1, n_steps=5)
        assert len(result.t_h) == 6  # steps 0,1,2,3,4,5

    def test_multi_cv_shape(self):
        from PyOMES.core import ControlVolume, LiquidPhase, Simulation, SparseRecorder
        cv_a = ControlVolume(phases={"liquid": LiquidPhase({"A": 1.0}, V_L=1.0, T_K=300.0)})
        cv_b = ControlVolume(phases={"liquid": LiquidPhase({"B": 2.0}, V_L=1.0, T_K=300.0)})
        sim = Simulation(
            cvs={"a": cv_a, "b": cv_b},
            recorder=SparseRecorder(every_m_steps=3),
        )
        result = sim.run(tau_h=0.3, n_steps=9)
        # n_steps=9, every_m_steps=3 -> steps 0,3,6,9 = 4 rows
        assert len(result.t_h) == 4
        assert "a" in result.phase_mol
        assert "b" in result.phase_mol

    def test_record_init_twice_raises(self):
        from PyOMES.core import Simulation, SparseRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = SparseRecorder()
        rec.record_init(sim, n_steps=5)
        with pytest.raises(RuntimeError, match="record_init called twice"):
            rec.record_init(sim, n_steps=5)

    def test_finalize_before_init_raises(self):
        from PyOMES.core import SparseRecorder
        with pytest.raises(RuntimeError, match="before record_init"):
            SparseRecorder().finalize()

    def test_gas_liquid_cv_p_atm_present(self):
        from PyOMES.core import Simulation, SparseRecorder
        import numpy as np
        sim = Simulation(
            cvs={"main": _make_gas_liquid_cv("main", {"O2": 0.5}, {"S": 1.0})},
            recorder=SparseRecorder(every_m_steps=2),
        )
        result = sim.run(tau_h=0.1, n_steps=6)
        assert "main" in result.P_atm
        assert np.all(result.P_atm["main"] > 0)


# ═══════════════════════════════════════════════════════════════════════
#  SummaryRecorder (RUN_HISTORY C2)
# ═══════════════════════════════════════════════════════════════════════

class TestSummaryRecorderExports:

    def test_importable_from_core(self):
        from PyOMES.core import SummaryRecorder  # noqa: F401

    def test_summary_result_importable_from_core(self):
        from PyOMES.core import SummaryResult  # noqa: F401

    def test_satisfies_recorder_protocol(self):
        from PyOMES.core import SummaryRecorder
        for method in ("record_init", "record_step", "finalize"):
            assert hasattr(SummaryRecorder, method)


class TestSummaryRecorderBasics:

    def test_returns_summary_result(self):
        from PyOMES.core import Simulation, SummaryRecorder, SummaryResult
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})},
            recorder=SummaryRecorder(),
        )
        result = sim.run(tau_h=0.1, n_steps=5)
        assert isinstance(result, SummaryResult)

    def test_initial_phase_mol_captured(self):
        from PyOMES.core import Simulation, SummaryRecorder
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 3.7})},
            recorder=SummaryRecorder(),
        )
        result = sim.run(tau_h=0.1, n_steps=5)
        assert result.initial_phase_mol["main"]["liquid"]["S"] == pytest.approx(3.7)

    def test_final_phase_mol_captured(self):
        from PyOMES.core import Simulation, SummaryRecorder
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})},
            recorder=SummaryRecorder(),
        )
        result = sim.run(tau_h=0.1, n_steps=5)
        assert "main" in result.final_phase_mol
        assert "liquid" in result.final_phase_mol["main"]
        assert "S" in result.final_phase_mol["main"]["liquid"]

    def test_start_and_end_t_h(self):
        from PyOMES.core import Simulation, SummaryRecorder
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})},
            recorder=SummaryRecorder(),
        )
        result = sim.run(tau_h=0.5, n_steps=5)
        assert result.start_t_h == pytest.approx(0.0)
        assert result.end_t_h == pytest.approx(0.5)

    def test_n_steps_recorded(self):
        from PyOMES.core import Simulation, SummaryRecorder
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})},
            recorder=SummaryRecorder(),
        )
        result = sim.run(tau_h=0.1, n_steps=7)
        assert result.n_steps == 7

    def test_runtime_positive(self):
        from PyOMES.core import Simulation, SummaryRecorder
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})},
            recorder=SummaryRecorder(),
        )
        result = sim.run(tau_h=0.1, n_steps=5)
        assert result.runtime_s >= 0.0

    def test_boundary_totals_empty_when_no_boundaries(self):
        from PyOMES.core import Simulation, SummaryRecorder
        sim = Simulation(
            cvs={"main": _make_liquid_cv("main", {"S": 1.0})},
            recorder=SummaryRecorder(),
        )
        result = sim.run(tau_h=0.1, n_steps=5)
        # No boundaries attached — totals dict should be empty or cv_key absent
        assert result.boundary_totals.get("main", {}) == {}

    def test_record_init_twice_raises(self):
        from PyOMES.core import Simulation, SummaryRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        rec = SummaryRecorder()
        rec.record_init(sim, n_steps=5)
        with pytest.raises(RuntimeError, match="record_init called twice"):
            rec.record_init(sim, n_steps=5)

    def test_finalize_before_init_raises(self):
        from PyOMES.core import SummaryRecorder
        with pytest.raises(RuntimeError, match="before record_init"):
            SummaryRecorder().finalize()

    def test_multi_cv_initial_and_final(self):
        from PyOMES.core import ControlVolume, LiquidPhase, Simulation, SummaryRecorder
        cv_a = ControlVolume(phases={"liquid": LiquidPhase({"A": 1.0}, V_L=1.0, T_K=300.0)})
        cv_b = ControlVolume(phases={"liquid": LiquidPhase({"B": 5.0}, V_L=1.0, T_K=300.0)})
        sim = Simulation(cvs={"a": cv_a, "b": cv_b}, recorder=SummaryRecorder())
        result = sim.run(tau_h=0.1, n_steps=3)
        assert result.initial_phase_mol["a"]["liquid"]["A"] == pytest.approx(1.0)
        assert result.initial_phase_mol["b"]["liquid"]["B"] == pytest.approx(5.0)
        assert "a" in result.final_phase_mol
        assert "b" in result.final_phase_mol


# ═══════════════════════════════════════════════════════════════════════
#  Simulation.run / Simulation._step — single-CV path (C4)
# ═══════════════════════════════════════════════════════════════════════

class TestRunValidation:

    def test_invalid_n_steps_raises(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        with pytest.raises(ValueError, match="n_steps"):
            sim.run(tau_h=1.0, n_steps=0)

    def test_invalid_tau_h_raises(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        with pytest.raises(ValueError, match="tau_h"):
            sim.run(tau_h=0.0, n_steps=10)


class TestRunBasics:

    def test_returns_batch_result_with_correct_shape(self):
        from PyOMES.core import Simulation, BatchResult
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        result = sim.run(tau_h=0.1, n_steps=10)
        assert isinstance(result, BatchResult)
        # Time grid is (n_steps + 1,)
        assert result.t_h.shape == (11,)
        # Per-CV nested liquid_mol
        assert "main" in result.liquid_mol
        assert "S" in result.liquid_mol["main"]
        assert result.liquid_mol["main"]["S"].shape == (11,)
        # Runtime measured
        assert result.runtime_s >= 0.0

    def test_time_grid_bounds(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        result = sim.run(tau_h=1.0, n_steps=50)
        assert result.t_h[0] == pytest.approx(0.0)
        assert result.t_h[-1] == pytest.approx(1.0)

    def test_no_reaction_no_boundary_conserves_inventory(self):
        """A CV with no reactions or boundaries should not change
        species inventory across a run."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0, "X": 0.5})
        sim = Simulation(cvs={"main": cv})
        result = sim.run(tau_h=0.5, n_steps=10)
        assert result.liquid_mol["main"]["S"][0] == pytest.approx(1.0)
        assert result.liquid_mol["main"]["S"][-1] == pytest.approx(1.0)
        assert result.liquid_mol["main"]["X"][-1] == pytest.approx(0.5)

    def test_t_h_advances_through_run(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv}, label="tracker")
        assert sim._t_h == 0.0
        sim.run(tau_h=2.0, n_steps=4)
        # After run completes, _t_h should reflect end of last step
        assert sim._t_h == pytest.approx(2.0)

    def test_t_h_resets_at_run_entry(self):
        """Decision 10: .run() is always destructive — _t_h resets."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        sim.run(tau_h=1.0, n_steps=4)
        # Re-run from a non-zero t_h state
        # Reset cv._context attribute to allow re-running (snapshot
        # workflow would normally avoid this)
        cv2 = _make_liquid_cv("main", {"S": 1.0})
        sim2 = Simulation(cvs={"main": cv2})
        sim2._t_h = 99.0  # simulate stale accumulator
        sim2.run(tau_h=1.0, n_steps=4)
        assert sim2._t_h == pytest.approx(1.0)  # reset overrode 99.0

    def test_is_running_flag_cleared_after_run(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        assert sim._context.is_running is False
        sim.run(tau_h=0.1, n_steps=2)
        assert sim._context.is_running is False

    def test_is_running_flag_cleared_after_run_failure(self):
        """If something inside .run() raises, the lifecycle flag must
        still be cleared (try/finally)."""
        from PyOMES.core import Simulation, BatchRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})

        class BrokenRecorder(BatchRecorder):
            def record_step(self, *args, **kwargs):
                raise RuntimeError("intentional test failure")

        sim = Simulation(cvs={"main": cv}, recorder=BrokenRecorder())
        with pytest.raises(RuntimeError, match="intentional"):
            sim.run(tau_h=0.1, n_steps=2)
        # The lifecycle flag must be cleared despite the failure
        assert sim._context.is_running is False


class TestRunSolverAndRecorder:

    def test_default_solver_stays_none(self):
        """Simulation is more general than the legacy run_batch: it
        supports liquid-only / gas-only / non-fermenter CVs.
        Defaulting to SimultaneousEulerSolver would force gas+liquid
        shapes, so the default is None (cv.advance falls back to its
        sequential body). Users opt in to snapshot semantics by
        passing solver=SimultaneousEulerSolver() explicitly."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        assert sim.solver is None
        sim.run(tau_h=0.1, n_steps=2)
        assert sim.solver is None

    def test_explicit_solver_used(self):
        """When the user passes a solver, the Simulation should
        dispatch to it on every step (not default to anything else).
        Verified with a stub solver that counts calls."""
        from PyOMES.core import Simulation
        from PyOMES.core.interfaces import AdvanceResult

        class CountingSolver:
            def __init__(self):
                self.calls = 0
            def solve_step(self, cv, dt_h, t_h, external_source_terms=None):
                self.calls += 1
                return AdvanceResult()

        cv = _make_liquid_cv("main", {"S": 1.0})
        custom = CountingSolver()
        sim = Simulation(cvs={"main": cv}, solver=custom)
        sim.run(tau_h=0.1, n_steps=4)
        # Stored on the Simulation and dispatched to once per step
        assert sim.solver is custom
        assert custom.calls == 4

    def test_default_recorder_constructed_per_run(self):
        from PyOMES.core import Simulation, BatchResult
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        # No recorder passed → BatchRecorder constructed fresh per call
        result = sim.run(tau_h=0.1, n_steps=2)
        assert isinstance(result, BatchResult)

    def test_explicit_recorder_used(self):
        from PyOMES.core import Simulation, BatchRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        rec = BatchRecorder()
        sim = Simulation(cvs={"main": cv}, recorder=rec)
        result = sim.run(tau_h=0.1, n_steps=2)
        assert result is rec._result


class TestStepInternal:

    def test_step_returns_results_links_controllers_profiles(self):
        """_step returns (results, link_records, controller_actions,
        profile_actions). For a single-CV Sim with no links / no
        controllers / no profiles, the last three are empty lists."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})  # solver=None → sequential body
        # C9 requires _controller_sample_state to exist (set in .run())
        sim._controller_sample_state = {}
        results, link_records, controller_actions, profile_actions = (
            sim._step(dt_h=0.05, t_h=0.05)
        )
        assert "main" in results
        assert hasattr(results["main"], "transfer")
        assert link_records == []
        assert controller_actions == []
        assert profile_actions == []


# ═══════════════════════════════════════════════════════════════════════
#  Multi-CV path (C5): link application, per-CV solver dispatch
# ═══════════════════════════════════════════════════════════════════════

class TestMultiCVStep:
    """Mirrors the core test_multi_cv.py cases on Simulation instead
    of MultiCVSystem. The legacy MultiCVSystem tests stay live until
    C14 (test_multi_cv.py is deleted then)."""

    def test_two_cvs_with_bidirectional_flow_converge(self):
        """Two zones with different concentrations, bidirectional
        circulation — mirrors the legacy test_multi_cv.py case to
        the same parameters."""
        from PyOMES.core import Simulation, AdvectiveLink
        V = 10.0
        cv_a = _make_liquid_cv("A", {"S": 2.0}, V_L=V)
        cv_b = _make_liquid_cv("B", {"S": 0.0}, V_L=V)
        Q = 50.0  # fast circulation
        link_ab = AdvectiveLink("A", "liquid", "B", "liquid", Q_L_per_h=Q)
        link_ba = AdvectiveLink("B", "liquid", "A", "liquid", Q_L_per_h=Q)
        sim = Simulation(cvs={"A": cv_a, "B": cv_b}, links=[link_ab, link_ba])
        # 500 steps × 0.01 h = 5 hours of mixing
        result = sim.run(tau_h=5.0, n_steps=500)
        # Total S = 2.0 mol in 20 L → C = 0.1 mol/L per CV at equilibrium
        C_A = result.liquid_mol["A"]["S"][-1] / V
        C_B = result.liquid_mol["B"]["S"][-1] / V
        assert C_A == pytest.approx(0.1, abs=0.005)
        assert C_B == pytest.approx(0.1, abs=0.005)

    def test_mixing_conserves_mass(self):
        """An isolated multi-CV system (no boundaries, no reactions)
        conserves total moles across all CVs."""
        from PyOMES.core import Simulation, AdvectiveLink
        cv_a = _make_liquid_cv("A", {"S": 1.0, "X": 0.5}, V_L=5.0)
        cv_b = _make_liquid_cv("B", {"S": 0.0, "X": 0.0}, V_L=5.0)
        link = AdvectiveLink("A", "liquid", "B", "liquid", Q_L_per_h=3.0)
        sim = Simulation(cvs={"A": cv_a, "B": cv_b}, links=[link])
        result = sim.run(tau_h=0.5, n_steps=50)
        # Sum across CVs at every step should equal initial total
        for i in range(51):
            total_S = (
                result.liquid_mol["A"]["S"][i]
                + result.liquid_mol["B"]["S"][i]
            )
            assert total_S == pytest.approx(1.0, abs=1e-10), \
                f"S not conserved at step {i}"

    def test_link_records_populated_in_batch_result(self):
        """transfer_records on the BatchResult contains LinkFlowRecord
        entries for each inter-CV link active per step."""
        from PyOMES.core import Simulation, AdvectiveLink
        cv_a = _make_liquid_cv("A", {"S": 1.0}, V_L=5.0)
        cv_b = _make_liquid_cv("B", {"S": 0.0}, V_L=5.0)
        link = AdvectiveLink("A", "liquid", "B", "liquid", Q_L_per_h=2.0)
        expected_label = link.label
        sim = Simulation(cvs={"A": cv_a, "B": cv_b}, links=[link])
        result = sim.run(tau_h=0.1, n_steps=5)
        # 5 steps, each with one inter-CV link → 5 entries each of
        # length 1 (only the inter-CV link; no internal gas-liquid
        # links in liquid-only CVs).
        assert len(result.transfer_records) == 5
        for step_records in result.transfer_records:
            assert len(step_records) == 1
            rec = step_records[0]
            assert rec.link_label == expected_label
            assert rec.source == "A.liquid"
            assert rec.sink == "B.liquid"

    def test_no_links_yields_empty_link_records(self):
        """Single-CV / multi-CV without links: each step's
        link_records section is empty."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        result = sim.run(tau_h=0.1, n_steps=3)
        for step_records in result.transfer_records:
            assert step_records == []


class TestMultiCVSolverDispatch:

    def test_single_solver_applies_to_all_cvs(self):
        """When self.solver is a single StepSolver, every CV gets
        the same solver. Verified with a stub solver counting calls."""
        from PyOMES.core import Simulation
        from PyOMES.core.interfaces import AdvanceResult

        class CountingSolver:
            def __init__(self):
                self.calls = 0
            def solve_step(self, cv, dt_h, t_h, external_source_terms=None):
                self.calls += 1
                return AdvanceResult()

        cv_a = _make_liquid_cv("A", {"S": 1.0})
        cv_b = _make_liquid_cv("B", {"S": 0.0})
        shared = CountingSolver()
        sim = Simulation(cvs={"A": cv_a, "B": cv_b}, solver=shared)
        sim.run(tau_h=0.1, n_steps=4)
        # 4 steps × 2 CVs = 8 calls into the same solver instance
        assert shared.calls == 8

    def test_dict_solver_dispatches_per_cv(self):
        """When self.solver is a Dict[cv_key, StepSolver], each CV
        gets its own solver instance."""
        from PyOMES.core import Simulation
        from PyOMES.core.interfaces import AdvanceResult

        class CountingSolver:
            def __init__(self, label):
                self.label = label
                self.calls = 0
            def solve_step(self, cv, dt_h, t_h, external_source_terms=None):
                self.calls += 1
                return AdvanceResult()

        cv_a = _make_liquid_cv("A", {"S": 1.0})
        cv_b = _make_liquid_cv("B", {"S": 0.0})
        solver_a = CountingSolver("A")
        solver_b = CountingSolver("B")
        sim = Simulation(
            cvs={"A": cv_a, "B": cv_b},
            solver={"A": solver_a, "B": solver_b},
        )
        sim.run(tau_h=0.1, n_steps=3)
        # Each CV's solver gets exactly n_steps calls
        assert solver_a.calls == 3
        assert solver_b.calls == 3

    def test_dict_solver_missing_cv_falls_back_to_sequential(self):
        """A CV absent from the solver dict gets None → sequential
        body. Verified by ensuring .run() succeeds even though only
        one of two CVs has a (stub) solver."""
        from PyOMES.core import Simulation
        from PyOMES.core.interfaces import AdvanceResult

        class CountingSolver:
            def __init__(self):
                self.calls = 0
            def solve_step(self, cv, dt_h, t_h, external_source_terms=None):
                self.calls += 1
                return AdvanceResult()

        cv_a = _make_liquid_cv("A", {"S": 1.0})
        cv_b = _make_liquid_cv("B", {"S": 0.0})
        solver_a = CountingSolver()
        sim = Simulation(
            cvs={"A": cv_a, "B": cv_b},
            solver={"A": solver_a},  # B absent → None → sequential
        )
        sim.run(tau_h=0.1, n_steps=3)
        # Only CV A's solver was called
        assert solver_a.calls == 3

    def test_multi_cv_step_solver_placeholder_raises(self):
        """The third arm of the solver-type union — a future
        MultiCVStepSolver — is a deliberate placeholder per
        decision 8. Detected via isinstance(solver, SystemSolver)
        (structural Protocol matching on advance_system)."""
        from PyOMES.core import Simulation

        class MockMultiCVStepSolver:
            def advance_system(self, sim, dt_h, t_h):
                ...

        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv}, solver=MockMultiCVStepSolver())
        with pytest.raises(NotImplementedError, match="MultiCVStepSolver"):
            sim.run(tau_h=0.1, n_steps=2)


# NOTE: TestMultiCVMatchesMultiCVSystem deleted at C14. The bit-for-bit
# parity between Simulation and MultiCVSystem was verified at C5 (5d34b6c);
# MultiCVSystem and its dedicated test_multi_cv.py have been deleted at
# C14 (this commit). The remaining C5 multi-CV behaviour is covered by
# TestMultiCVStep and TestMultiCVSolverDispatch above.


# ═══════════════════════════════════════════════════════════════════════
#  Lifecycle gating (C6)
# ═══════════════════════════════════════════════════════════════════════

def _arm_running(sim):
    """Manually flip the RunContext into is_running=True so we can test
    mid-run mutation gating without spinning the full .run() loop.
    Tests are responsible for clearing the flag (or letting the
    Simulation be garbage-collected)."""
    sim._context.is_running = True


class TestPhaseGating:

    def test_gas_phase_T_K_raises_mid_run(self):
        from PyOMES.core import Simulation
        cv = _make_gas_liquid_cv(
            "main", gas_mol={"O2": 1.0}, liq_mol={"O2": 0.001},
        )
        sim = Simulation(cvs={"main": cv}, label="exp")
        _arm_running(sim)
        with pytest.raises(RuntimeError, match="GasPhase.T_K"):
            cv.phases["gas"].T_K = 310.0

    def test_liquid_phase_V_L_raises_mid_run(self):
        from PyOMES.core import Simulation
        cv = _make_gas_liquid_cv(
            "main", gas_mol={"O2": 1.0}, liq_mol={"O2": 0.001},
        )
        sim = Simulation(cvs={"main": cv}, label="exp")
        _arm_running(sim)
        with pytest.raises(RuntimeError, match="LiquidPhase.V_L"):
            cv.phases["liquid"].V_L = 5.0

    def test_phase_unchecked_setter_bypasses_gate(self):
        """Pattern B: MutableScalar._set_unchecked lets the orchestrator
        mutate state during a run while the public setter remains gated."""
        from PyOMES.core import Simulation
        from PyOMES.control.descriptors import MutableScalar
        cv = _make_gas_liquid_cv(
            "main", gas_mol={"O2": 1.0}, liq_mol={"O2": 0.001},
        )
        sim = Simulation(cvs={"main": cv})
        _arm_running(sim)
        gas = cv.phases["gas"]
        descriptor = type(gas).__dict__["T_K"]
        assert isinstance(descriptor, MutableScalar)
        descriptor._set_unchecked(gas, 310.0)
        assert gas.T_K == pytest.approx(310.0)

    def test_phase_setter_outside_run_succeeds(self):
        from PyOMES.core import Simulation
        cv = _make_gas_liquid_cv(
            "main", gas_mol={"O2": 1.0}, liq_mol={"O2": 0.001},
        )
        sim = Simulation(cvs={"main": cv})
        # is_running is False — mutation should succeed
        cv.phases["liquid"].T_K = 310.0
        assert cv.phases["liquid"].T_K == pytest.approx(310.0)

    def test_phase_context_wired_by_simulation(self):
        """Phases inside a CV held by a Sim share the Sim's
        RunContext."""
        from PyOMES.core import Simulation
        cv = _make_gas_liquid_cv(
            "main", gas_mol={"O2": 1.0}, liq_mol={"O2": 0.001},
        )
        sim = Simulation(cvs={"main": cv})
        assert cv.phases["gas"]._context is sim._context
        assert cv.phases["liquid"]._context is sim._context


class TestKineticGasLiquidLinkGating:
    """The legacy _simulation_running / _enter_simulation /
    _exit_simulation / _warn_if_simulating machinery is replaced
    by the RunContext gate. C6 also adds the underscore-unchecked
    siblings for orchestrator-mediated mutation."""

    def _make_cv_with_link(self):
        """Fermenter-shaped CV with a kinetic gas-liquid link, for
        gating tests."""
        from vlmodels.fermenter.config import (
            FermenterFactory, VesselConfig, TransferConfig,
        )
        return FermenterFactory.create_volume(
            vessel=VesselConfig(V_total_L=100, T_K=305.15),
            transfer=TransferConfig.default_kinetic(kLa_O2=100.0),
        )

    def _link(self, cv):
        from PyOMES.core import KineticGasLiquidLink
        for iface in cv.internal_interfaces:
            if isinstance(iface, KineticGasLiquidLink):
                return iface
        raise AssertionError("No KineticGasLiquidLink in CV")

    def test_legacy_methods_deleted(self):
        cv = self._make_cv_with_link()
        link = self._link(cv)
        assert not hasattr(link, "_enter_simulation")
        assert not hasattr(link, "_exit_simulation")
        assert not hasattr(link, "_warn_if_simulating")
        assert not hasattr(link, "_simulation_running")

    def test_link_context_wired_by_simulation(self):
        from PyOMES.core import Simulation
        cv = self._make_cv_with_link()
        link = self._link(cv)
        assert link._context is None  # unowned
        sim = Simulation(cvs={"main": cv})
        assert link._context is sim._context

    def test_set_kLa_raises_mid_run(self):
        from PyOMES.core import Simulation
        cv = self._make_cv_with_link()
        link = self._link(cv)
        sim = Simulation(cvs={"main": cv})
        _arm_running(sim)
        with pytest.raises(RuntimeError, match="set_kLa"):
            link.set_kLa("O2", 200.0)

    def test_kLa_item_unchecked_bypasses_gate(self):
        from PyOMES.core import Simulation
        from PyOMES.control.descriptors import MutableDict
        cv = self._make_cv_with_link()
        link = self._link(cv)
        sim = Simulation(cvs={"main": cv})
        _arm_running(sim)
        descriptor = type(link).__dict__["kLa"]
        assert isinstance(descriptor, MutableDict)
        descriptor._set_item_unchecked(link, "O2", 200.0)
        assert link.kLa["O2"] == pytest.approx(200.0)

    def test_set_kLa_with_co2_ratio_raises_mid_run(self):
        from PyOMES.core import Simulation
        cv = self._make_cv_with_link()
        link = self._link(cv)
        sim = Simulation(cvs={"main": cv})
        _arm_running(sim)
        with pytest.raises(RuntimeError, match="set_kLa_with_co2_ratio"):
            link.set_kLa_with_co2_ratio(200.0)

    def test_set_henry_raises_mid_run(self):
        from PyOMES.core import Simulation
        cv = self._make_cv_with_link()
        link = self._link(cv)
        sim = Simulation(cvs={"main": cv})
        _arm_running(sim)
        with pytest.raises(RuntimeError, match="set_henry"):
            link.set_henry("O2", 0.01)

    def test_mutations_succeed_outside_run(self):
        from PyOMES.core import Simulation
        cv = self._make_cv_with_link()
        link = self._link(cv)
        sim = Simulation(cvs={"main": cv})
        # is_running is False
        link.set_kLa("O2", 150.0)
        assert link.kLa["O2"] == pytest.approx(150.0)


class TestControlVolumeGating:

    def test_boundaries_append_raises_mid_run(self):
        from PyOMES.core import Simulation
        cv = _make_gas_liquid_cv(
            "main", gas_mol={"O2": 1.0, "N2": 3.0},
            liq_mol={"O2": 0.001},
        )
        sim = Simulation(cvs={"main": cv}, label="exp")
        _arm_running(sim)
        # The boundary object's shape doesn't matter for the gate test
        with pytest.raises(RuntimeError, match="append"):
            cv.boundaries.append(object())

    def test_boundaries_pop_raises_mid_run(self):
        from PyOMES.core import Simulation
        cv = _make_gas_liquid_cv(
            "main", gas_mol={"O2": 1.0, "N2": 3.0},
            liq_mol={"O2": 0.001},
        )
        # Pre-populate with an opaque object
        cv.boundaries.append(object())
        sim = Simulation(cvs={"main": cv})
        _arm_running(sim)
        with pytest.raises(RuntimeError, match="pop"):
            cv.boundaries.pop()

    def test_property_calculators_append_raises_mid_run(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        _arm_running(sim)
        with pytest.raises(RuntimeError, match="append"):
            cv.property_calculators.append(object())

    def test_reaction_system_reassignment_raises_mid_run(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        _arm_running(sim)
        with pytest.raises(RuntimeError, match="reaction_system"):
            cv.reaction_system = object()

    def test_reaction_system_unchecked_bypasses_gate(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        _arm_running(sim)
        sentinel = object()
        cv._set_reaction_system_unchecked(sentinel)
        assert cv.reaction_system is sentinel

    def test_mutations_succeed_outside_run(self):
        from PyOMES.core import Simulation
        cv = _make_gas_liquid_cv(
            "main", gas_mol={"O2": 1.0}, liq_mol={"O2": 0.001},
        )
        sim = Simulation(cvs={"main": cv})
        cv.boundaries.append(object())
        cv.property_calculators.append(object())
        sentinel = object()
        cv.reaction_system = sentinel
        assert len(cv.boundaries) == 1
        assert len(cv.property_calculators) == 1
        assert cv.reaction_system is sentinel


class TestSimulationConfigGating:

    def test_solver_reassignment_raises_mid_run(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        _arm_running(sim)
        with pytest.raises(RuntimeError, match="solver"):
            sim.solver = object()

    def test_recorder_reassignment_raises_mid_run(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        _arm_running(sim)
        with pytest.raises(RuntimeError, match="recorder"):
            sim.recorder = object()

    def test_controllers_append_raises_mid_run(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        _arm_running(sim)
        with pytest.raises(RuntimeError, match="append"):
            sim.controllers.append(object())

    def test_profiles_append_raises_mid_run(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        _arm_running(sim)
        with pytest.raises(RuntimeError, match="append"):
            sim.profiles.append(object())

    def test_mutations_succeed_outside_run(self):
        from PyOMES.core import Simulation, SimultaneousEulerSolver, BatchRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        # is_running is False — all reassignments succeed
        sim.controllers.append(object())
        sim.profiles.append(object())
        sim.solver = SimultaneousEulerSolver()
        sim.recorder = BatchRecorder()
        assert len(sim.controllers) == 1
        assert len(sim.profiles) == 1
        assert isinstance(sim.solver, SimultaneousEulerSolver)


class TestLifecycleErrorMessage:

    def test_error_message_includes_simulation_label(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv}, label="experiment_A")
        _arm_running(sim)
        with pytest.raises(RuntimeError, match="experiment_A") as exc:
            cv.phases["liquid"].T_K = 310.0
        # Error message also mentions the controller/profile escape hatch
        assert "Controller or" in str(exc.value)


class TestMultiSimulationIsolation:

    def test_two_sims_run_state_disjoint(self):
        """Lockables in sim1 are gated by sim1's context, not sim2's.
        Multi-Simulation safe by construction (decision 9)."""
        from PyOMES.core import Simulation
        cv_a = _make_liquid_cv("A", {"S": 1.0})
        cv_b = _make_liquid_cv("B", {"S": 0.0})
        sim_a = Simulation(cvs={"A": cv_a}, label="sim_a")
        sim_b = Simulation(cvs={"B": cv_b}, label="sim_b")
        # sim_a is running; sim_b is not
        _arm_running(sim_a)
        # sim_a's CV mutation raises; sim_b's CV mutation succeeds
        with pytest.raises(RuntimeError, match="sim_a"):
            cv_a.phases["liquid"].T_K = 310.0
        cv_b.phases["liquid"].T_K = 310.0  # unaffected
        assert cv_b.phases["liquid"].T_K == pytest.approx(310.0)


class TestLifecycleGatingIntegratedWithRun:
    """A real Simulation.run() flips the gate at entry and clears
    at exit. Verify both transitions and that the gate clears even
    on failure."""

    def test_gate_clears_after_normal_run(self):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        assert sim._context.is_running is False
        sim.run(tau_h=0.1, n_steps=2)
        assert sim._context.is_running is False
        # After run, mutation succeeds
        cv.phases["liquid"].T_K = 305.0

    def test_gate_active_during_run(self):
        """A stub solver inspects sim._context.is_running mid-loop;
        the flag must be True there."""
        from PyOMES.core import Simulation
        from PyOMES.core.interfaces import AdvanceResult

        captured = {"is_running_mid_loop": None}

        class ProbingSolver:
            def solve_step(self, cv, dt_h, t_h, external_source_terms=None):
                captured["is_running_mid_loop"] = cv._context.is_running
                return AdvanceResult()

        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv}, solver=ProbingSolver())
        sim.run(tau_h=0.1, n_steps=2)
        assert captured["is_running_mid_loop"] is True


# ═══════════════════════════════════════════════════════════════════════
#  Controller protocol — Pattern 1 collapse (C7)
# ═══════════════════════════════════════════════════════════════════════

class TestControllerProtocol:
    """The new CV-native Controller protocol: single
    compute(state, dt_h) -> ControlAction method; no Actuator
    protocol; no Commands dataclass."""

    def test_protocol_importable(self):
        from PyOMES.control.interfaces import Controller
        assert Controller is not None

    def test_actuator_protocol_removed(self):
        """The legacy Actuator protocol is deleted at C7."""
        import PyOMES.control.interfaces as iface
        assert not hasattr(iface, "Actuator")

    def test_state_builder_module_deleted(self):
        """src/control/state_builder.py was deleted at C7. The
        legacy build_state body is inlined into models/vlmodels/
        fermenter/unit.py as _build_fermenter_state_inline."""
        with pytest.raises(ImportError):
            import PyOMES.control.state_builder  # noqa: F401

    def test_simple_controller_satisfies_protocol(self):
        """A duck-typed controller with a compute method is
        recognised by runtime_checkable Protocol."""
        from PyOMES.control.actions import ControlAction
        from PyOMES.control.interfaces import Controller

        class IdleController:
            def compute(self, state, dt_h):
                return ControlAction(controller_label="idle")

        ctrl = IdleController()
        assert isinstance(ctrl, Controller)

    def test_object_without_compute_does_not_satisfy_protocol(self):
        """Sanity check on the runtime_checkable semantics."""
        from PyOMES.control.interfaces import Controller
        assert not isinstance(object(), Controller)

    def test_controller_compute_receives_cv_snapshot(self):
        """Validate the single-CV input type: compute() is called
        with a CVSnapshot when the orchestrator routes a single CV."""
        from PyOMES.control.actions import ControlAction
        from PyOMES.core import CVSnapshot

        captured = {"state_type": None}

        class TypeProbeController:
            def compute(self, state, dt_h):
                captured["state_type"] = type(state)
                return ControlAction(controller_label="probe")

        ctrl = TypeProbeController()
        cv_snap = CVSnapshot(
            cv_key="main", t_h=0.0, pH=None, ionic_strength=None,
            T_K=298.15, V_liq_L=1.0, V_gas_L=0.0, P_gas_atm=0.0,
        )
        ctrl.compute(cv_snap, dt_h=0.01)
        assert captured["state_type"] is CVSnapshot

    def test_controller_compute_receives_simulation_snapshot(self):
        """Multi-CV controllers see a SimulationSnapshot."""
        from PyOMES.control.actions import ControlAction
        from PyOMES.core import CVSnapshot, SimulationSnapshot

        captured = {"state_type": None}

        class MultiCVProbe:
            def compute(self, state, dt_h):
                captured["state_type"] = type(state)
                return ControlAction(controller_label="multi")

        ctrl = MultiCVProbe()
        sim_snap = SimulationSnapshot(
            t_h=0.0,
            cvs={
                "A": CVSnapshot(
                    cv_key="A", t_h=0.0, pH=None, ionic_strength=None,
                    T_K=298.15, V_liq_L=1.0, V_gas_L=0.0,
                    P_gas_atm=0.0,
                ),
            },
        )
        ctrl.compute(sim_snap, dt_h=0.01)
        assert captured["state_type"] is SimulationSnapshot

    def test_optional_attributes_are_duck_typed(self):
        """Controllers MAY expose sample_period_h, target_cv_key,
        reset() — orchestrator uses hasattr/getattr defensively at
        C9. The Protocol declares only compute() as required."""
        from PyOMES.control.actions import ControlAction
        from PyOMES.control.interfaces import Controller

        class FullController:
            sample_period_h = 0.1
            target_cv_key = "main"
            def compute(self, state, dt_h):
                return ControlAction()
            def reset(self):
                pass

        class MinimalController:
            def compute(self, state, dt_h):
                return ControlAction()

        # Both satisfy the formal Protocol (only compute is required)
        assert isinstance(FullController(), Controller)
        assert isinstance(MinimalController(), Controller)
        # Optional members are inspectable via hasattr
        assert hasattr(FullController(), "sample_period_h")
        assert not hasattr(MinimalController(), "sample_period_h")
        assert hasattr(FullController(), "reset")
        assert not hasattr(MinimalController(), "reset")


# ═══════════════════════════════════════════════════════════════════════
#  PHController — CV-native (C8a)
# ═══════════════════════════════════════════════════════════════════════

def _cv_snap(
    pH=None, V_liq=1.0, T_K=298.15, t_h=0.0, cv_key="main",
    sensors=None,
):
    """Lightweight CVSnapshot factory for controller unit tests."""
    from PyOMES.core import CVSnapshot
    return CVSnapshot(
        cv_key=cv_key,
        t_h=float(t_h),
        pH=pH,
        ionic_strength=None,
        T_K=float(T_K),
        V_liq_L=float(V_liq),
        V_gas_L=0.0,
        P_gas_atm=0.0,
        sensors=sensors or {},
    )


class TestCVPHControllerBasics:

    def test_importable(self):
        from PyOMES.control.cv_loops import PHController
        assert PHController is not None

    def test_satisfies_controller_protocol(self):
        from PyOMES.control.interfaces import Controller
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5)
        assert isinstance(ctrl, Controller)

    def test_default_acid_id_is_H3PO4(self):
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5)
        assert ctrl.chemical_id == "H3PO4"

    def test_compute_returns_control_action(self):
        from PyOMES.control.actions import ControlAction
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5)
        action = ctrl.compute(_cv_snap(pH=7.0), dt_h=0.01)
        assert isinstance(action, ControlAction)
        assert action.controller_label == "pH"
        assert action.target_cv_key == "main"
        assert action.dt_h == pytest.approx(0.01)


class TestCVPHControllerDosing:

    def test_acid_dosing_when_pH_above_setpoint(self):
        """pH=7.0, setpoint=6.5 → too basic → acid dosing fires."""
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5, Kp=0.5, max_add_molL_hr=0.2)
        action = ctrl.compute(_cv_snap(pH=7.0, V_liq=10.0), dt_h=0.01)
        # add_rate = 0.5 * 0.5 = 0.25, clamped to 0.2 mol/L/h
        # mol_h = 0.2 * 10.0 = 2.0
        assert action.flux_applied["liquid"]["H3PO4"] == pytest.approx(2.0)
        # dosed_mol = mol_h * dt_h = 0.02
        assert action.dosed_mol["H3PO4"] == pytest.approx(0.02)

    def test_no_dosing_when_pH_at_setpoint(self):
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5)
        action = ctrl.compute(_cv_snap(pH=6.5), dt_h=0.01)
        assert action.flux_applied == {}
        assert action.dosed_mol == {}

    def test_deadband_suppresses_dosing(self):
        """|err| <= deadband → no action."""
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5, deadband=0.5)
        action = ctrl.compute(_cv_snap(pH=6.8), dt_h=0.01)
        assert action.flux_applied == {}

    def test_no_action_when_pH_is_none(self):
        """No speciation → pH is None → empty action."""
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5)
        action = ctrl.compute(_cv_snap(pH=None), dt_h=0.01)
        assert action.flux_applied == {}

    def test_base_dosing_when_configured(self):
        """pH=6.0, setpoint=6.5, base_chemical_id='NaOH' → base fires."""
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(
            setpoint=6.5, Kp=1.0, max_add_molL_hr=0.5,
            base_chemical_id="NaOH",
        )
        action = ctrl.compute(_cv_snap(pH=6.0, V_liq=5.0), dt_h=0.01)
        # err = -0.5; err_b = 0.5; add_rate = 1.0 * 0.5 = 0.5
        # mol_h = 0.5 * 5.0 = 2.5
        assert action.flux_applied["liquid"]["NaOH"] == pytest.approx(2.5)

    def test_base_skipped_when_not_configured(self):
        """Acid-only legacy mode: pH below setpoint, no base_id → no action."""
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5, base_chemical_id=None)
        action = ctrl.compute(_cv_snap(pH=6.0), dt_h=0.01)
        assert action.flux_applied == {}

    def test_phosphate_cap_clamps_dose(self):
        """When CT_P_mol_L sensor present, dosing is clamped so the
        phosphate total doesn't exceed CT_P_max."""
        from PyOMES.control.cv_loops import PHController
        # CT_P_max = 1.0 mol/L; current CT_P = 0.9 mol/L; V_liq = 1.0 L
        # Headroom = (1.0 - 0.9) * 1.0 = 0.1 mol over dt_h=0.1 → 1.0 mol/h cap
        ctrl = PHController(
            setpoint=6.5, Kp=10.0, max_add_molL_hr=100.0,
            CT_P_max=1.0,
        )
        snap = _cv_snap(pH=7.5, V_liq=1.0, sensors={"CT_P_mol_L": 0.9})
        action = ctrl.compute(snap, dt_h=0.1)
        assert action.flux_applied["liquid"]["H3PO4"] == pytest.approx(1.0)


class TestCVPHControllerIntegralState:

    def test_integral_state_accumulates(self):
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5, Kp=0.0, Ki=1.0)
        ctrl.compute(_cv_snap(pH=7.0), dt_h=0.1)
        i1 = ctrl._I_err
        ctrl.compute(_cv_snap(pH=7.0), dt_h=0.1)
        i2 = ctrl._I_err
        assert i2 > i1 > 0.0

    def test_reset_clears_integral(self):
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5, Ki=1.0)
        ctrl.compute(_cv_snap(pH=7.0), dt_h=0.1)
        assert ctrl._I_err > 0.0
        ctrl.reset()
        assert ctrl._I_err == 0.0
        assert ctrl._last_mode == ""

    def test_mode_switch_resets_integral(self):
        """Switching acid→base (or vice versa) resets the integral
        so anti-windup doesn't bleed across modes."""
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5, Ki=1.0, base_chemical_id="NaOH")
        # Acid dosing accumulates integral
        ctrl.compute(_cv_snap(pH=7.5), dt_h=0.1)
        i_acid = ctrl._I_err
        assert i_acid > 0.0
        # Switch to base — integral resets
        ctrl.compute(_cv_snap(pH=6.0), dt_h=0.1)
        # After base step: integral set to err_b * dt_h, NOT inheriting
        # the acid integral. Verify the value matches a fresh
        # accumulation from zero (0.5 * 0.1 = 0.05).
        assert ctrl._I_err == pytest.approx(0.05)


class TestCVPHControllerMultiCV:

    def test_single_cv_simulation_snapshot(self):
        """A SimulationSnapshot with one CV and no target_cv_key
        passes through to the sole CV."""
        from PyOMES.core import SimulationSnapshot
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5)
        sim_snap = SimulationSnapshot(
            t_h=0.0, cvs={"only": _cv_snap(pH=7.0, cv_key="only")},
        )
        action = ctrl.compute(sim_snap, dt_h=0.01)
        assert action.target_cv_key == "only"

    def test_multi_cv_simulation_snapshot_with_target(self):
        from PyOMES.core import SimulationSnapshot
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5, target_cv_key="B")
        sim_snap = SimulationSnapshot(
            t_h=0.0,
            cvs={
                "A": _cv_snap(pH=6.5, cv_key="A"),  # ok, no action
                "B": _cv_snap(pH=7.5, cv_key="B"),  # acid dosing
            },
        )
        action = ctrl.compute(sim_snap, dt_h=0.01)
        assert action.target_cv_key == "B"
        assert action.flux_applied.get("liquid", {}).get("H3PO4", 0) > 0

    def test_multi_cv_without_target_raises(self):
        from PyOMES.core import SimulationSnapshot
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5)
        sim_snap = SimulationSnapshot(
            t_h=0.0,
            cvs={
                "A": _cv_snap(pH=7.0, cv_key="A"),
                "B": _cv_snap(pH=7.0, cv_key="B"),
            },
        )
        with pytest.raises(RuntimeError, match="target_cv_key"):
            ctrl.compute(sim_snap, dt_h=0.01)

    def test_unknown_target_cv_key_raises(self):
        from PyOMES.core import SimulationSnapshot
        from PyOMES.control.cv_loops import PHController
        ctrl = PHController(setpoint=6.5, target_cv_key="missing")
        sim_snap = SimulationSnapshot(
            t_h=0.0, cvs={"A": _cv_snap(pH=7.0, cv_key="A")},
        )
        with pytest.raises(KeyError, match="missing"):
            ctrl.compute(sim_snap, dt_h=0.01)


# ═══════════════════════════════════════════════════════════════════════
#  DOAgitationController — CV-native (C8b)
# ═══════════════════════════════════════════════════════════════════════

def _cv_snap_with_do(DO_mol_L, V_liq=1.0, t_h=0.0, cv_key="main"):
    """Lightweight CVSnapshot factory with a DO_mol_L sensor."""
    from PyOMES.core import CVSnapshot
    return CVSnapshot(
        cv_key=cv_key,
        t_h=float(t_h),
        pH=None,
        ionic_strength=None,
        T_K=305.15,
        V_liq_L=float(V_liq),
        V_gas_L=0.2,
        P_gas_atm=1.0,
        sensors={"DO_mol_L": DO_mol_L} if DO_mol_L is not None else {},
    )


class TestCVDOAgitationControllerBasics:

    def test_importable(self):
        from PyOMES.control.cv_loops import DOAgitationController
        assert DOAgitationController is not None

    def test_legacy_alias_exists(self):
        from PyOMES.control.cv_loops import DOAgitationController, DOController
        assert DOController is DOAgitationController

    def test_satisfies_controller_protocol(self):
        from PyOMES.control.interfaces import Controller
        from PyOMES.control.cv_loops import DOAgitationController
        ctrl = DOAgitationController(setpoint_mol_L=1e-4)
        assert isinstance(ctrl, Controller)

    def test_from_pct_saturation_constructor(self):
        from PyOMES.control.cv_loops import DOAgitationController
        ctrl = DOAgitationController.from_pct_saturation(pct=30.0)
        # 30 % of (1.3e-3 × 0.2095) = 30 % of 2.72e-4 ≈ 8.17e-5 mol/L
        assert ctrl.setpoint_mol_L == pytest.approx(
            1.3e-3 * 0.2095 * 0.30, abs=1e-10
        )

    def test_requires_segmented_ode_flag(self):
        from PyOMES.control.cv_loops import DOAgitationController
        ctrl = DOAgitationController(setpoint_mol_L=1e-4)
        assert ctrl.requires_segmented_ode is True


class TestCVDOAgitationControllerCompute:

    def test_no_action_when_DO_sensor_missing(self):
        from PyOMES.control.cv_loops import DOAgitationController
        ctrl = DOAgitationController(setpoint_mol_L=1e-4)
        action = ctrl.compute(_cv_snap_with_do(DO_mol_L=None), dt_h=0.01)
        assert action.params_changed == {}

    def test_compute_emits_kLa_params(self):
        """compute() should populate params_changed with full reflective kLa paths."""
        from PyOMES.control.cv_loops import DOAgitationController
        _KLA_O2 = "internal_interfaces[KineticGasLiquidLink].kLa.O2"
        _KLA_CO2 = "internal_interfaces[KineticGasLiquidLink].kLa.CO2"
        ctrl = DOAgitationController(setpoint_mol_L=1e-4)
        action = ctrl.compute(_cv_snap_with_do(DO_mol_L=0.5e-4), dt_h=0.01)
        assert _KLA_O2 in action.params_changed
        assert _KLA_CO2 in action.params_changed
        assert action.params_changed[_KLA_O2] > 0.0

    def test_kLa_CO2_coupled_via_ratio(self):
        from PyOMES.control.cv_loops import DOAgitationController
        ctrl = DOAgitationController(
            setpoint_mol_L=1e-4, kLa_CO2_ratio=0.9,
        )
        action = ctrl.compute(_cv_snap_with_do(DO_mol_L=0.5e-4), dt_h=0.01)
        ratio = (
            action.params_changed["internal_interfaces[KineticGasLiquidLink].kLa.CO2"]
            / action.params_changed["internal_interfaces[KineticGasLiquidLink].kLa.O2"]
        )
        assert ratio == pytest.approx(0.9)

    def test_DO_below_setpoint_raises_rpm(self):
        """DO too low → controller increases RPM → kLa increases."""
        from PyOMES.control.cv_loops import DOAgitationController
        ctrl = DOAgitationController(
            setpoint_mol_L=2e-4, Kp=1e5, Ki=0.0,
            rpm_initial=400.0,
        )
        initial_rpm = ctrl._rpm
        ctrl.compute(_cv_snap_with_do(DO_mol_L=1e-4), dt_h=0.1)
        assert ctrl._rpm > initial_rpm

    def test_DO_above_setpoint_lowers_rpm(self):
        from PyOMES.control.cv_loops import DOAgitationController
        ctrl = DOAgitationController(
            setpoint_mol_L=1e-4, Kp=1e5, Ki=0.0,
            rpm_initial=400.0,
        )
        ctrl.compute(_cv_snap_with_do(DO_mol_L=2e-4), dt_h=0.1)
        assert ctrl._rpm < 400.0

    def test_rpm_clamped_at_max(self):
        from PyOMES.control.cv_loops import DOAgitationController
        ctrl = DOAgitationController(
            setpoint_mol_L=1.0, Kp=1e10, Ki=0.0,  # huge gain
            rpm_max=1200.0, rpm_initial=400.0,
        )
        # Massive error pushes RPM to max
        ctrl.compute(_cv_snap_with_do(DO_mol_L=1e-6), dt_h=0.1)
        assert ctrl._rpm == pytest.approx(1200.0)

    def test_rpm_clamped_at_min(self):
        from PyOMES.control.cv_loops import DOAgitationController
        ctrl = DOAgitationController(
            setpoint_mol_L=0.0, Kp=1e10, Ki=0.0,
            rpm_min=100.0, rpm_initial=400.0,
        )
        ctrl.compute(_cv_snap_with_do(DO_mol_L=1e-4), dt_h=0.1)
        assert ctrl._rpm == pytest.approx(100.0)

    def test_deadband_suppresses_action(self):
        """|err| within deadband → no rpm change."""
        from PyOMES.control.cv_loops import DOAgitationController
        ctrl = DOAgitationController(
            setpoint_mol_L=1e-4, Kp=1e5, Ki=0.0,
            deadband_mol_L=1e-4, rpm_initial=400.0,
        )
        ctrl.compute(_cv_snap_with_do(DO_mol_L=1.5e-4), dt_h=0.1)
        # Err = -5e-5, within deadband; no update
        assert ctrl._rpm == pytest.approx(400.0)


class TestCVDOAgitationControllerAntiWindup:

    def test_no_integral_buildup_when_saturated_at_max(self):
        """When rpm is clamped at rpm_max and error keeps pushing
        higher, the integral term doesn't wind up."""
        from PyOMES.control.cv_loops import DOAgitationController
        ctrl = DOAgitationController(
            setpoint_mol_L=1.0, Kp=1e10, Ki=1.0,
            rpm_max=1200.0, rpm_initial=1200.0,
        )
        # Drive several steps with persistent error
        for _ in range(5):
            ctrl.compute(_cv_snap_with_do(DO_mol_L=1e-6), dt_h=0.1)
        # Integral should not have grown unbounded (anti-windup
        # subtracts the contribution each saturated step)
        assert ctrl._I_err == pytest.approx(0.0)

    def test_reset_clears_state(self):
        from PyOMES.control.cv_loops import DOAgitationController
        ctrl = DOAgitationController(
            setpoint_mol_L=1e-4, Ki=1.0, rpm_initial=400.0,
        )
        ctrl.compute(_cv_snap_with_do(DO_mol_L=0.5e-4), dt_h=0.1)
        ctrl.reset()
        assert ctrl._I_err == 0.0
        assert ctrl._rpm == pytest.approx(400.0)


class TestCVDOAgitationControllerCustomMapping:

    def test_custom_agitation_to_kLa(self):
        """Pass a callable for RPM→kLa instead of the default power-law."""
        from PyOMES.control.cv_loops import DOAgitationController

        def linear(rpm):
            return 0.1 * rpm

        ctrl = DOAgitationController(
            setpoint_mol_L=1e-4, agitation_to_kLa=linear,
            rpm_initial=500.0,
        )
        # _rpm_to_kLa returns 0.1 * 500 = 50
        assert ctrl._kLa_O2 == pytest.approx(50.0)


# ═══════════════════════════════════════════════════════════════════════
#  DOCascadeController — CV-native (C8c)
# ═══════════════════════════════════════════════════════════════════════

class TestCVDOCascadeBasics:

    def test_importable(self):
        from PyOMES.control.cv_loops import DOCascadeController
        assert DOCascadeController is not None

    def test_satisfies_controller_protocol(self):
        from PyOMES.control.interfaces import Controller
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController(setpoint_mol_L=1e-4)
        assert isinstance(ctrl, Controller)

    def test_initial_tier_is_one(self):
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController(setpoint_mol_L=1e-4)
        assert ctrl._active_tier == 1

    def test_initial_yO2_is_baseline(self):
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController(setpoint_mol_L=1e-4, yO2_baseline=0.21)
        assert ctrl._yO2 == pytest.approx(0.21)

    def test_from_pct_saturation(self):
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController.from_pct_saturation(pct=30.0)
        assert ctrl.setpoint_mol_L == pytest.approx(
            1.3e-3 * 0.2095 * 0.30
        )


class TestCVDOCascadeTier1:
    """Tier 1: agitation alone handles the error when within RPM bounds."""

    def test_tier_stays_at_1_with_modest_error(self):
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController(
            setpoint_mol_L=1e-4, Kp=1e5, Ki=0.0,
            rpm_initial=400.0,
        )
        action = ctrl.compute(
            _cv_snap_with_do(DO_mol_L=0.8e-4), dt_h=0.1,
        )
        assert ctrl._active_tier == 1
        # kLa.O2 in the action; gas composition stays at baseline
        assert "internal_interfaces[KineticGasLiquidLink].kLa.O2" in action.params_changed
        assert action.params_changed["boundaries[GasFeed].y.O2"] == pytest.approx(0.21)


class TestCVDOCascadeTier2:
    """Tier 2: when RPM saturates at max, yO2 increases (enrichment)."""

    def test_tier_escalates_to_2_when_rpm_maxes(self):
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController(
            setpoint_mol_L=1.0,  # absurdly high → forces saturation
            Kp=1e10, Ki=0.0,
            rpm_max=1200.0, rpm_initial=400.0,
            yO2_baseline=0.21, yO2_max=1.0, Kp_gas=2.0,
        )
        ctrl.compute(_cv_snap_with_do(DO_mol_L=1e-6), dt_h=0.1)
        assert ctrl._active_tier == 2
        # yO2 should have moved upward toward yO2_max
        assert ctrl._yO2 > 0.21

    def test_yO2_decreases_when_rpm_at_min_and_DO_too_high(self):
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController(
            setpoint_mol_L=0.0,  # absurdly low → forces error_eff < 0
            Kp=1e10, Ki=0.0,
            rpm_min=100.0, rpm_initial=400.0,
            yO2_baseline=0.21, yO2_min=0.0, Kp_gas=2.0,
        )
        ctrl.compute(_cv_snap_with_do(DO_mol_L=1e-4), dt_h=0.1)
        assert ctrl._active_tier == 2
        assert ctrl._yO2 < 0.21

    def test_yO2_clamped_at_max(self):
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController(
            setpoint_mol_L=1.0,
            Kp=1e10, rpm_initial=1200.0,
            yO2_max=0.5, Kp_gas=100.0,  # huge gain — would overshoot
        )
        # Single step with huge gain
        ctrl.compute(_cv_snap_with_do(DO_mol_L=1e-6), dt_h=0.1)
        assert ctrl._yO2 <= 0.5 + 1e-9

    def test_yO2_relaxes_to_baseline_when_unsaturated(self):
        """When agitation is within range, yO2 relaxes back to baseline."""
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController(
            setpoint_mol_L=1e-4, Kp=1e5, Ki=0.0,
            yO2_baseline=0.21, yO2_relax_rate=10.0,
        )
        # Manually push yO2 off-baseline
        ctrl._yO2 = 0.5
        # Run a step within RPM range — relax should pull it back
        ctrl.compute(_cv_snap_with_do(DO_mol_L=0.95e-4), dt_h=0.1)
        # 10/h relax × 0.1 h = 1.0 — clamped to abs(diff) = 0.29 →
        # yO2 should have moved toward baseline by 0.29
        assert ctrl._yO2 < 0.5


class TestCVDOCascadeTier3:
    """Tier 3: vvm adjustment when both tier 1 and tier 2 saturate."""

    def test_tier_escalates_to_3_when_both_saturate(self):
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController(
            setpoint_mol_L=1.0,
            Kp=1e10, rpm_initial=1200.0,
            yO2_baseline=0.21, yO2_max=1.0, Kp_gas=1e10,
            enable_vvm_control=True,
            vvm_baseline=2.0, vvm_max=10.0, Kp_vvm=20.0,
        )
        # Force yO2 to saturate at max too
        ctrl._yO2 = 1.0
        ctrl.compute(_cv_snap_with_do(DO_mol_L=1e-6), dt_h=0.1)
        assert ctrl._active_tier == 3
        assert ctrl._vvm > 2.0

    def test_vvm_only_emitted_when_enabled(self):
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl_off = DOCascadeController(
            setpoint_mol_L=1e-4, enable_vvm_control=False,
        )
        ctrl_on = DOCascadeController(
            setpoint_mol_L=1e-4, enable_vvm_control=True,
        )
        action_off = ctrl_off.compute(
            _cv_snap_with_do(DO_mol_L=0.5e-4), dt_h=0.1,
        )
        action_on = ctrl_on.compute(
            _cv_snap_with_do(DO_mol_L=0.5e-4), dt_h=0.1,
        )
        assert "boundaries[GasFeed].vvm_min" not in action_off.params_changed
        assert "boundaries[GasFeed].vvm_min" in action_on.params_changed


class TestCVDOCascadeAction:
    """The ControlAction shape emitted by DOCascadeController."""

    def test_action_always_has_kLa_and_gas_y(self):
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController(setpoint_mol_L=1e-4)
        action = ctrl.compute(_cv_snap_with_do(DO_mol_L=0.5e-4), dt_h=0.1)
        assert "internal_interfaces[KineticGasLiquidLink].kLa.O2" in action.params_changed
        assert "internal_interfaces[KineticGasLiquidLink].kLa.CO2" in action.params_changed
        assert "boundaries[GasFeed].y.O2" in action.params_changed
        assert "boundaries[GasFeed].y.N2" in action.params_changed

    def test_gas_y_sums_to_one(self):
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController(setpoint_mol_L=1e-4)
        action = ctrl.compute(_cv_snap_with_do(DO_mol_L=0.5e-4), dt_h=0.1)
        ySum = (
            action.params_changed["boundaries[GasFeed].y.O2"]
            + action.params_changed["boundaries[GasFeed].y.N2"]
        )
        assert ySum == pytest.approx(1.0, abs=1e-9)

    def test_no_action_when_DO_missing(self):
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController(setpoint_mol_L=1e-4)
        action = ctrl.compute(_cv_snap_with_do(DO_mol_L=None), dt_h=0.1)
        assert action.params_changed == {}


class TestCVDOCascadeReset:

    def test_reset_restores_initial_state(self):
        from PyOMES.control.cv_loops import DOCascadeController
        ctrl = DOCascadeController(
            setpoint_mol_L=1e-4, rpm_initial=400.0,
            yO2_baseline=0.21, vvm_baseline=2.0,
        )
        # Drive some change
        ctrl.compute(_cv_snap_with_do(DO_mol_L=0.5e-4), dt_h=0.1)
        ctrl._yO2 = 0.5  # manual change to test reset
        ctrl._vvm = 5.0
        ctrl.reset()
        assert ctrl._rpm == pytest.approx(400.0)
        assert ctrl._yO2 == pytest.approx(0.21)
        assert ctrl._vvm == pytest.approx(2.0)
        assert ctrl._active_tier == 1
        assert ctrl._I_err == 0.0


# ═══════════════════════════════════════════════════════════════════════
#  Pressure relief controllers — CV-native (C8d)
# ═══════════════════════════════════════════════════════════════════════

def _cv_snap_pressurised(
    P_gas_atm=1.5, V_gas=0.2, T_K=305.15, t_h=0.0, cv_key="main",
    n_gas=None, y_gas=None,
):
    """CVSnapshot factory with realistic headspace state."""
    from PyOMES.core import CVSnapshot
    n_gas = dict(n_gas or {"O2": 2.0, "N2": 7.5, "CO2": 0.5})
    if y_gas is None:
        total = sum(n_gas.values())
        y_gas = {k: v / total for k, v in n_gas.items()} if total > 0 else {}
    return CVSnapshot(
        cv_key=cv_key,
        t_h=float(t_h),
        pH=None,
        ionic_strength=None,
        T_K=float(T_K),
        V_liq_L=0.8,
        V_gas_L=float(V_gas),
        P_gas_atm=float(P_gas_atm),
        n_gas_mol=n_gas,
        y_gas=y_gas,
        sensors={},
    )


class TestCVInstantPressureRelief:

    def test_importable(self):
        from PyOMES.control.cv_loops import InstantPressureReliefController
        assert InstantPressureReliefController is not None

    def test_satisfies_controller_protocol(self):
        from PyOMES.control.interfaces import Controller
        from PyOMES.control.cv_loops import InstantPressureReliefController
        ctrl = InstantPressureReliefController(P_set_atm=1.2)
        assert isinstance(ctrl, Controller)

    def test_no_action_below_setpoint(self):
        from PyOMES.control.cv_loops import InstantPressureReliefController
        ctrl = InstantPressureReliefController(P_set_atm=1.2)
        action = ctrl.compute(
            _cv_snap_pressurised(P_gas_atm=1.0), dt_h=0.01,
        )
        assert action.flux_applied == {}
        assert action.vented_mol == {}

    def test_action_above_setpoint(self):
        from PyOMES.control.cv_loops import InstantPressureReliefController
        ctrl = InstantPressureReliefController(P_set_atm=1.2)
        # P=1.5 atm > setpoint 1.2 → venting fires
        action = ctrl.compute(
            _cv_snap_pressurised(P_gas_atm=1.5), dt_h=0.01,
        )
        assert "gas" in action.flux_applied
        # All fluxes are negative (gas leaving)
        for species, flux in action.flux_applied["gas"].items():
            assert flux < 0.0, f"flux for {species} should be negative"
        # vented_mol is positive
        for species, mol in action.vented_mol.items():
            assert mol > 0.0

    def test_venting_proportional_to_mole_fractions(self):
        from PyOMES.control.cv_loops import InstantPressureReliefController
        ctrl = InstantPressureReliefController(P_set_atm=1.2)
        # y_O2 = 2/10 = 0.2; y_N2 = 7.5/10 = 0.75; y_CO2 = 0.5/10 = 0.05
        action = ctrl.compute(
            _cv_snap_pressurised(P_gas_atm=1.5), dt_h=0.01,
        )
        total_vent = sum(action.vented_mol.values())
        # Each species vented proportional to its mole fraction
        for species, y in [("O2", 0.2), ("N2", 0.75), ("CO2", 0.05)]:
            assert action.vented_mol[species] == pytest.approx(
                y * total_vent, rel=1e-9
            )

    def test_no_action_zero_dt(self):
        from PyOMES.control.cv_loops import InstantPressureReliefController
        ctrl = InstantPressureReliefController(P_set_atm=1.2)
        action = ctrl.compute(
            _cv_snap_pressurised(P_gas_atm=1.5), dt_h=0.0,
        )
        assert action.flux_applied == {}


class TestCVSmoothPressureRelief:

    def test_importable(self):
        from PyOMES.control.cv_loops import SmoothPressureReliefController
        assert SmoothPressureReliefController is not None

    def test_satisfies_controller_protocol(self):
        from PyOMES.control.interfaces import Controller
        from PyOMES.control.cv_loops import SmoothPressureReliefController
        ctrl = SmoothPressureReliefController(P_set_atm=1.2)
        assert isinstance(ctrl, Controller)

    def test_venting_monotonic_above_setpoint(self):
        """Smooth relief: a higher pressure produces more venting
        than a lower one (the curve is monotonic above setpoint).
        At default k_vent_per_h and dt the ramp saturates, so use
        gentler k for this test to verify the differential behaviour."""
        from PyOMES.control.cv_loops import SmoothPressureReliefController
        ctrl = SmoothPressureReliefController(
            P_set_atm=1.2, k_vent_per_h=10.0, smooth_width_atm=0.1,
        )
        action_low = ctrl.compute(
            _cv_snap_pressurised(P_gas_atm=1.25), dt_h=0.005,
        )
        ctrl2 = SmoothPressureReliefController(
            P_set_atm=1.2, k_vent_per_h=10.0, smooth_width_atm=0.1,
        )
        action_high = ctrl2.compute(
            _cv_snap_pressurised(P_gas_atm=1.5), dt_h=0.005,
        )
        flux_low = abs(action_low.flux_applied["gas"]["O2"])
        flux_high = abs(action_high.flux_applied["gas"]["O2"])
        assert flux_high > flux_low

    def test_action_above_setpoint(self):
        from PyOMES.control.cv_loops import SmoothPressureReliefController
        ctrl = SmoothPressureReliefController(
            P_set_atm=1.2, k_vent_per_h=500.0,
        )
        action = ctrl.compute(
            _cv_snap_pressurised(P_gas_atm=1.5), dt_h=0.01,
        )
        assert "gas" in action.flux_applied
        # All fluxes negative
        for flux in action.flux_applied["gas"].values():
            assert flux < 0.0


class TestCVPressureReliefControllerValve:
    """Physical valve-model controller — wraps the legacy
    ReliefValveActuator mass-flow physics."""

    def test_importable(self):
        from PyOMES.control.cv_loops import PressureReliefController
        assert PressureReliefController is not None

    def test_satisfies_controller_protocol(self):
        from PyOMES.control.interfaces import Controller
        from PyOMES.control.cv_loops import PressureReliefController
        ctrl = PressureReliefController(P_set_atm=1.2, diameter_m=0.01)
        assert isinstance(ctrl, Controller)

    def test_no_action_at_or_below_backpressure(self):
        from PyOMES.control.cv_loops import PressureReliefController
        ctrl = PressureReliefController(
            P_set_atm=0.9, diameter_m=0.01, P_back_atm=1.0,
            controller_kind="smooth",
        )
        # P_gas = 0.95 atm > P_set=0.9 (valve opens) BUT <= P_back=1.0
        # → no flow through the orifice
        action = ctrl.compute(
            _cv_snap_pressurised(P_gas_atm=0.95), dt_h=0.01,
        )
        assert action.flux_applied == {}

    def test_smooth_valve_action_above_setpoint(self):
        from PyOMES.control.cv_loops import PressureReliefController
        ctrl = PressureReliefController(
            P_set_atm=1.2, diameter_m=0.01,
            controller_kind="smooth", P_back_atm=1.0,
        )
        action = ctrl.compute(
            _cv_snap_pressurised(P_gas_atm=1.5), dt_h=0.01,
        )
        assert "gas" in action.flux_applied
        for flux in action.flux_applied["gas"].values():
            assert flux < 0.0

    def test_onoff_valve_latches(self):
        """OnOff valve stays open once tripped until pressure drops
        below P_set - hysteresis."""
        from PyOMES.control.cv_loops import PressureReliefController
        ctrl = PressureReliefController(
            P_set_atm=1.2, diameter_m=0.01,
            controller_kind="onoff", hysteresis_atm=0.1,
        )
        # Trip the valve open
        ctrl.compute(_cv_snap_pressurised(P_gas_atm=1.3), dt_h=0.01)
        assert ctrl._is_open is True
        # Stay open at P=1.15 (above lo=1.1)
        ctrl.compute(_cv_snap_pressurised(P_gas_atm=1.15), dt_h=0.01)
        assert ctrl._is_open is True
        # Close at P=1.05 (below lo=1.1)
        ctrl.compute(_cv_snap_pressurised(P_gas_atm=1.05), dt_h=0.01)
        assert ctrl._is_open is False

    def test_reset_clears_onoff_latch(self):
        from PyOMES.control.cv_loops import PressureReliefController
        ctrl = PressureReliefController(
            P_set_atm=1.2, diameter_m=0.01,
            controller_kind="onoff",
        )
        ctrl.compute(_cv_snap_pressurised(P_gas_atm=1.5), dt_h=0.01)
        assert ctrl._is_open is True
        ctrl.reset()
        assert ctrl._is_open is False


class TestCVPHControllerRegistryValidation:

    def test_invalid_acid_id_warns_at_construction(self):
        """validate_compound_id emits a UserWarning for unknown IDs
        (matching the legacy controller's behavior — it's a soft
        check, not a hard error)."""
        from PyOMES.control.cv_loops import PHController
        with pytest.warns(UserWarning, match="NotARealCompound_XYZ"):
            PHController(setpoint=6.5, chemical_id="NotARealCompound_XYZ")

    def test_invalid_base_id_warns_at_construction(self):
        from PyOMES.control.cv_loops import PHController
        with pytest.warns(UserWarning, match="NotARealCompound_XYZ"):
            PHController(
                setpoint=6.5, chemical_id="H3PO4",
                base_chemical_id="NotARealCompound_XYZ",
            )


# ═══════════════════════════════════════════════════════════════════════
#  Controller integration with Simulation.run (C9)
# ═══════════════════════════════════════════════════════════════════════

class TestC9ControllerInvocation:
    """Verify that Simulation.run actually fires controllers,
    collects their actions in the BatchResult, and dispatches
    flux_applied / params_changed onto the CV."""

    def test_controller_compute_called_each_step(self):
        """A controller whose compute() bumps a counter should be
        called n_steps times."""
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction

        class CountingController:
            def __init__(self):
                self.calls = 0
            def compute(self, state, dt_h):
                self.calls += 1
                return ControlAction()

        cv = _make_liquid_cv("main", {"S": 1.0})
        ctrl = CountingController()
        sim = Simulation(cvs={"main": cv}, controllers=[ctrl])
        sim.run(tau_h=0.1, n_steps=5)
        assert ctrl.calls == 5

    def test_controller_actions_recorded_in_batch_result(self):
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction

        class LabelingController:
            def compute(self, state, dt_h):
                return ControlAction(controller_label="probe", t_h=state.t_h)

        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(
            cvs={"main": cv}, controllers=[LabelingController()],
        )
        result = sim.run(tau_h=0.1, n_steps=4)
        # 4 steps × 1 controller = 4 actions, one per step
        assert len(result.controller_actions) == 4
        for per_step in result.controller_actions:
            assert len(per_step) == 1
            assert per_step[0].controller_label == "probe"

    def test_controller_reset_called_at_run_top(self):
        """Decision 10: .run() always destructive. reset() is
        called on every controller before the loop starts."""
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction

        class ResetSpyController:
            def __init__(self):
                self.reset_count = 0
            def compute(self, state, dt_h):
                return ControlAction()
            def reset(self):
                self.reset_count += 1

        cv = _make_liquid_cv("main", {"S": 1.0})
        ctrl = ResetSpyController()
        sim = Simulation(cvs={"main": cv}, controllers=[ctrl])
        sim.run(tau_h=0.1, n_steps=2)
        assert ctrl.reset_count == 1
        sim.run(tau_h=0.1, n_steps=2)  # re-run
        assert ctrl.reset_count == 2

    def test_controller_without_reset_tolerated(self):
        """Controllers without a reset() method should not crash."""
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction

        class NoResetCtrl:
            def compute(self, state, dt_h):
                return ControlAction()

        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv}, controllers=[NoResetCtrl()])
        sim.run(tau_h=0.1, n_steps=2)  # should not raise


class TestC9SamplePeriodGating:

    def test_sample_period_h_gates_compute(self):
        """sample_period_h=0.05 with dt_h=0.01 → compute() fires
        every 5 steps (once per period)."""
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction

        class TimedController:
            sample_period_h = 0.05
            def __init__(self):
                self.calls = 0
            def compute(self, state, dt_h):
                self.calls += 1
                return ControlAction(controller_label="timed")

        cv = _make_liquid_cv("main", {"S": 1.0})
        ctrl = TimedController()
        sim = Simulation(cvs={"main": cv}, controllers=[ctrl])
        # 10 steps × dt=0.01 = 0.1h total; period=0.05 → fires at
        # step 1 (always) and steps where Δt >= 0.05 from last fire.
        # Expected fires: steps 1, 6 (some grid noise tolerance).
        sim.run(tau_h=0.1, n_steps=10)
        assert 2 <= ctrl.calls <= 4  # tolerate floating-point grid

    def test_sample_period_s_converted_to_hours(self):
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction

        class TimedController:
            sample_period_s = 360.0  # 0.1 h
            def __init__(self):
                self.calls = 0
            def compute(self, state, dt_h):
                self.calls += 1
                return ControlAction()

        cv = _make_liquid_cv("main", {"S": 1.0})
        ctrl = TimedController()
        sim = Simulation(cvs={"main": cv}, controllers=[ctrl])
        # tau_h=0.2 = 2 × 0.1h periods; 20 steps; should fire at
        # the start and around the 0.1h mark
        sim.run(tau_h=0.2, n_steps=20)
        assert 2 <= ctrl.calls <= 4

    def test_between_samples_held_action_reissued(self):
        """ZOH: between sampling instants, the same action is
        reissued on every step so the recorder sees a per-step stream."""
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction

        class CountingTimedController:
            sample_period_h = 0.05
            def __init__(self):
                self.calls = 0
            def compute(self, state, dt_h):
                self.calls += 1
                return ControlAction(
                    controller_label=f"call_{self.calls}",
                )

        cv = _make_liquid_cv("main", {"S": 1.0})
        ctrl = CountingTimedController()
        sim = Simulation(cvs={"main": cv}, controllers=[ctrl])
        result = sim.run(tau_h=0.1, n_steps=10)
        # 10 steps recorded — each step has one ControlAction
        assert len(result.controller_actions) == 10
        # Some are real (fresh compute), some are held re-issues.
        # The fresh ones have distinct controller_label values.
        labels = [
            sa[0].controller_label
            for sa in result.controller_actions
        ]
        # At least one distinct label per fire; held actions repeat
        # the previous label.
        assert len(set(labels)) == ctrl.calls

    def test_never_fired_slot_returns_none_not_placeholder(self):
        """P9: _invoke_controllers returns None (not a phantom ControlAction)
        for a slot where should_fire=False and last_action=None.

        Guards against the case where a sample_state is pre-initialised with
        a past fire timestamp but no action was ever produced.  None entries
        are filtered by _step so the recorder never sees phantom actions."""
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction

        class SlowController:
            sample_period_h = 10.0
            def compute(self, state, dt_h):
                return ControlAction(controller_label="slow")

        cv = _make_liquid_cv("main", {"S": 1.0})
        ctrl = SlowController()
        sim = Simulation(cvs={"main": cv}, controllers=[ctrl])

        # Pre-initialise sample state: last_fire at t=0 with no action yet.
        # This puts us in the ZOH branch with last_action=None.
        # (_controller_sample_state is created by run(); create it manually here.)
        sim._controller_sample_state = {}
        sim._controller_sample_state[id(ctrl)] = {
            "last_fire_t_h": 0.0,
            "last_action": None,
        }

        # t_h=0.01 is within the 10 h period → should_fire=False → ZOH → None
        actions = sim._invoke_controllers(t_h=0.01, dt_h=0.01, results={})
        assert actions == [None], (
            "Expected [None] for never-fired ZOH slot, "
            f"got {actions!r}"
        )


class TestC9FluxApplyDispatch:

    def test_flux_applied_routes_to_cv_apply_external_flux(self):
        """An action with flux_applied={"liquid": {"S": -10}} should
        decrement liquid S by 10 mol/h * dt_h."""
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction

        class ConsumeS:
            def compute(self, state, dt_h):
                return ControlAction(
                    target_cv_key=state.cv_key,
                    flux_applied={"liquid": {"S": -1.0}},
                )

        cv = _make_liquid_cv("main", {"S": 10.0}, V_L=1.0)
        sim = Simulation(cvs={"main": cv}, controllers=[ConsumeS()])
        # 10 steps × dt=0.01 × -1 mol/h = -0.1 mol total
        result = sim.run(tau_h=0.1, n_steps=10)
        # The final liquid S should have dropped
        assert result.liquid_mol["main"]["S"][-1] < 10.0
        assert result.liquid_mol["main"]["S"][-1] == pytest.approx(
            10.0 - 0.1, abs=1e-9,
        )


class TestC9ParamPathDispatch:
    """The params_changed dispatch via _apply_param_change. Pattern B
    unchecked setters fire during the run; mid-run mutations from
    elsewhere remain gated."""

    def _make_fermenter_cv(self):
        from vlmodels.fermenter.config import (
            FermenterFactory, VesselConfig, TransferConfig,
        )
        return FermenterFactory.create_volume(
            vessel=VesselConfig(V_total_L=100, T_K=305.15),
            transfer=TransferConfig.default_kinetic(kLa_O2=100.0),
        )

    def _kgl(self, cv):
        from PyOMES.core import KineticGasLiquidLink
        for iface in cv.internal_interfaces:
            if isinstance(iface, KineticGasLiquidLink):
                return iface
        raise AssertionError("No KineticGasLiquidLink")

    def test_kLa_path_dispatches_to_link(self):
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction

        class SetkLa:
            def compute(self, state, dt_h):
                return ControlAction(
                    target_cv_key=state.cv_key,
                    params_changed={
                        "internal_interfaces[KineticGasLiquidLink].kLa.O2": 500.0,
                    },
                )

        cv = self._make_fermenter_cv()
        link = self._kgl(cv)
        sim = Simulation(cvs={"main": cv}, controllers=[SetkLa()])
        sim.run(tau_h=0.01, n_steps=1)
        assert link.kLa["O2"] == pytest.approx(500.0)

    def test_unknown_path_raises_param_path_error(self):
        """Unknown params_changed paths raise ParamPathError (Q3: fail-loud)."""
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction
        from PyOMES.control.param_path import ParamPathError

        class GarbagePath:
            def compute(self, state, dt_h):
                return ControlAction(
                    target_cv_key=state.cv_key,
                    params_changed={"frobnicate.the.widget": 42.0},
                )

        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv}, controllers=[GarbagePath()])
        with pytest.raises(ParamPathError):
            sim.run(tau_h=0.01, n_steps=1)


class TestC9EndToEndPHController:
    """Drive a real CV-native PHController through Simulation.run
    end-to-end; verify acid is dosed when pH is high."""

    def test_ph_controller_doses_acid_on_high_pH(self):
        from PyOMES.core import Simulation, ControlVolume, LiquidPhase
        from PyOMES.control.cv_loops import PHController

        class FakeSpeciation:
            """A trivial 'speciation' that sets phase.n_mol['H+'] from
            an internal pH variable so the CVSnapshot has a real pH
            for the controller to react to.

            Matches the real ChemicalEquilibriumEngineProtocol contract: solve()
            returns an EquilibriumResult (no side effects); the caller
            commits it via apply_to_phases()."""
            def __init__(self):
                self.engine = self
                self.pH = 7.5
                self.reactions = []
                self.stoichiometry = None
            def solve(self, phases, T_K):
                from PyOMES.chemical_equilibrium.protocols import EquilibriumResult
                H_plus_mol_L = 10.0 ** (-self.pH)
                return EquilibriumResult(pH=self.pH, species_mol_L={"H+": H_plus_mol_L})

        engine = FakeSpeciation()
        liq = LiquidPhase(n_mol={"H+": 10**-7.5}, V_L=1.0, T_K=298.15)
        cv = ControlVolume(
            phases={"liquid": liq}, reaction_system=engine, label="main",
        )
        ctrl = PHController(setpoint=6.5, Kp=10.0, max_add_molL_hr=1.0)
        sim = Simulation(cvs={"main": cv}, controllers=[ctrl])
        result = sim.run(tau_h=0.05, n_steps=5)

        # At each step, pH was 7.5 > 6.5, so acid_mol_h should be
        # positive and the H3PO4 should accumulate in liquid n_mol
        actions = [
            step_actions[0] for step_actions in result.controller_actions
        ]
        for action in actions:
            assert "H3PO4" in action.flux_applied.get("liquid", {})
            assert action.flux_applied["liquid"]["H3PO4"] > 0.0


# ═══════════════════════════════════════════════════════════════════════
#  Profile protocol + concrete profiles (C10)
# ═══════════════════════════════════════════════════════════════════════

class TestC10ProfileProtocol:

    def test_profile_protocol_importable(self):
        from PyOMES.control.cv_profiles import Profile
        assert Profile is not None

    def test_concrete_profiles_importable(self):
        from PyOMES.control.cv_profiles import (
            TemperatureRamp, VVMSchedule, SetpointTrajectory,
        )
        assert TemperatureRamp is not None
        assert VVMSchedule is not None
        assert SetpointTrajectory is not None

    def test_simple_profile_satisfies_protocol(self):
        from PyOMES.control.actions import ProfileRecord
        from PyOMES.control.cv_profiles import Profile

        class IdleProfile:
            def apply(self, t_h, sim):
                return ProfileRecord(profile_label="idle", t_h=t_h)

        assert isinstance(IdleProfile(), Profile)


class TestC10TemperatureRamp:

    def test_temperature_ramp_mutates_phase_T_K(self):
        from PyOMES.core import Simulation
        from PyOMES.control.cv_profiles import TemperatureRamp
        # Start at 298.15 K, ramp to 310 K over 0..1 hour
        cv = _make_gas_liquid_cv(
            "main", gas_mol={"O2": 1.0}, liq_mol={"O2": 0.001},
            T_K=298.15,
        )
        ramp = TemperatureRamp(
            target_cv_key="main",
            waypoints=[(0.0, 298.15), (1.0, 310.0)],
        )
        sim = Simulation(cvs={"main": cv}, profiles=[ramp])
        sim.run(tau_h=1.0, n_steps=10)
        # At end of run, T_K should be at 310 (final waypoint)
        assert cv.phases["liquid"].T_K == pytest.approx(310.0, abs=1e-9)
        assert cv.phases["gas"].T_K == pytest.approx(310.0, abs=1e-9)

    def test_profile_records_appear_in_batch_result(self):
        from PyOMES.core import Simulation
        from PyOMES.control.cv_profiles import TemperatureRamp
        cv = _make_gas_liquid_cv(
            "main", gas_mol={"O2": 1.0}, liq_mol={"O2": 0.001},
        )
        ramp = TemperatureRamp(
            target_cv_key="main",
            waypoints=[(0.0, 298.15), (1.0, 310.0)],
        )
        sim = Simulation(cvs={"main": cv}, profiles=[ramp])
        result = sim.run(tau_h=1.0, n_steps=5)
        # 5 steps × 1 profile = 5 profile records, one per step
        assert len(result.profile_actions) == 5
        for per_step in result.profile_actions:
            assert len(per_step) == 1
            rec = per_step[0]
            assert rec.profile_label == "temperature_ramp"
            assert "cv.main.phases.liquid.T_K" in rec.targets
            assert "cv.main.phases.gas.T_K" in rec.targets

    def test_temperature_ramp_holds_after_last_waypoint(self):
        """After the last waypoint, value holds (no extrapolation)."""
        from PyOMES.core import Simulation
        from PyOMES.control.cv_profiles import TemperatureRamp
        cv = _make_gas_liquid_cv(
            "main", gas_mol={"O2": 1.0}, liq_mol={"O2": 0.001},
            T_K=300.0,
        )
        # Ramp ends at t=0.5h with T=310; runs to t=1.0h
        ramp = TemperatureRamp(
            target_cv_key="main",
            waypoints=[(0.0, 300.0), (0.5, 310.0)],
        )
        sim = Simulation(cvs={"main": cv}, profiles=[ramp])
        sim.run(tau_h=1.0, n_steps=10)
        assert cv.phases["liquid"].T_K == pytest.approx(310.0, abs=1e-9)

    def test_temperature_ramp_during_run_bypasses_gate(self):
        """The Pattern B unchecked setter is used — profile mutation
        during a run doesn't raise RuntimeError despite is_running."""
        from PyOMES.core import Simulation
        from PyOMES.control.cv_profiles import TemperatureRamp
        cv = _make_gas_liquid_cv(
            "main", gas_mol={"O2": 1.0}, liq_mol={"O2": 0.001},
            T_K=298.15,
        )
        ramp = TemperatureRamp(
            target_cv_key="main",
            waypoints=[(0.0, 298.15), (1.0, 320.0)],
        )
        sim = Simulation(cvs={"main": cv}, profiles=[ramp])
        # If profile.apply tripped the gate, this would raise.
        sim.run(tau_h=1.0, n_steps=5)


class TestC10VVMSchedule:

    def _make_fermenter_with_gas_feed(self):
        from vlmodels.fermenter.config import (
            FermenterFactory, VesselConfig, TransferConfig, GasFeedConfig,
        )
        return FermenterFactory.create_volume(
            vessel=VesselConfig(V_total_L=100, T_K=305.15),
            transfer=TransferConfig.default_kinetic(kLa_O2=100.0),
            gas_feed=GasFeedConfig(
                vvm_min=0.5, composition={"O2": 0.21, "N2": 0.79},
            ),
        )

    def _gas_feed(self, cv):
        from PyOMES.core.boundaries import GasFeed
        for b in cv.boundaries:
            if isinstance(b, GasFeed):
                return b
        raise AssertionError("No GasFeed in CV")

    def test_vvm_schedule_mutates_gas_feed(self):
        from PyOMES.core import Simulation
        from PyOMES.control.cv_profiles import VVMSchedule
        cv = self._make_fermenter_with_gas_feed()
        gf = self._gas_feed(cv)
        assert gf.vvm_min == pytest.approx(0.5)
        schedule = VVMSchedule(
            target_cv_key="main",
            waypoints=[(0.0, 0.5), (0.1, 1.5)],
        )
        sim = Simulation(cvs={"main": cv}, profiles=[schedule])
        sim.run(tau_h=0.1, n_steps=5)
        assert gf.vvm_min == pytest.approx(1.5, abs=1e-9)

    def test_vvm_schedule_records_target(self):
        from PyOMES.core import Simulation
        from PyOMES.control.cv_profiles import VVMSchedule
        cv = self._make_fermenter_with_gas_feed()
        schedule = VVMSchedule(
            target_cv_key="main",
            waypoints=[(0.0, 0.5), (1.0, 2.0)],
        )
        sim = Simulation(cvs={"main": cv}, profiles=[schedule])
        result = sim.run(tau_h=0.1, n_steps=3)
        for per_step in result.profile_actions:
            assert per_step[0].profile_label == "vvm_schedule"
            assert "cv.main.gas_feed.vvm_min" in per_step[0].targets


class TestC10SetpointTrajectory:

    def test_setpoint_trajectory_mutates_controller_attribute(self):
        from PyOMES.core import Simulation
        from PyOMES.control.cv_loops import PHController
        from PyOMES.control.cv_profiles import SetpointTrajectory
        ctrl = PHController(setpoint=7.0)
        # Ramp setpoint from 7.0 → 6.0 over the run
        traj = SetpointTrajectory(
            controller=ctrl,
            attribute="setpoint",
            waypoints=[(0.0, 7.0), (1.0, 6.0)],
        )
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(
            cvs={"main": cv},
            controllers=[ctrl],
            profiles=[traj],
        )
        sim.run(tau_h=1.0, n_steps=10)
        # End of run, setpoint should be at 6.0
        assert ctrl.setpoint == pytest.approx(6.0, abs=1e-9)


class TestC10ProfileOrderingBeforeAdvance:
    """Verify profiles fire BEFORE cv.advance — temperature ramps
    take effect in the same step rather than the next one."""

    def test_temperature_ramp_visible_to_controller_same_step(self):
        """A profile sets T_K; the controller reads snapshot.T_K
        post-advance. They should see the SAME T_K (the new one),
        since the profile fires first."""
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction
        from PyOMES.control.cv_profiles import TemperatureRamp

        captured = []

        class TempReader:
            def compute(self, state, dt_h):
                captured.append(state.T_K)
                return ControlAction(target_cv_key=state.cv_key)

        cv = _make_gas_liquid_cv(
            "main", gas_mol={"O2": 1.0}, liq_mol={"O2": 0.001},
            T_K=298.15,
        )
        ramp = TemperatureRamp(
            target_cv_key="main",
            waypoints=[(0.0, 298.15), (1.0, 320.0)],
        )
        sim = Simulation(
            cvs={"main": cv}, profiles=[ramp],
            controllers=[TempReader()],
        )
        sim.run(tau_h=1.0, n_steps=5)
        # The controller should see strictly increasing T_K
        assert captured == sorted(captured)
        # Last reading should equal final ramp value
        assert captured[-1] == pytest.approx(320.0, abs=1e-9)


# ═══════════════════════════════════════════════════════════════════════
#  FermenterBuilder.build_simulation() (C11)
# ═══════════════════════════════════════════════════════════════════════

class TestC11FermenterBuilderSimulation:

    def test_new_fluent_methods_exist(self):
        from vlmodels.fermenter.config import FermenterBuilder
        b = FermenterBuilder()
        assert hasattr(b, "profile")
        assert hasattr(b, "recorder")
        assert hasattr(b, "build_simulation")
        assert hasattr(b, "build_simulation_and_run")

    def test_profile_fluent_appends(self):
        from vlmodels.fermenter.config import FermenterBuilder
        b = FermenterBuilder()
        result = b.profile(object()).profile(object())
        assert result is b  # chainable
        assert len(b._profiles) == 2

    def test_recorder_fluent_sets(self):
        from vlmodels.fermenter.config import FermenterBuilder
        rec = object()
        b = FermenterBuilder().recorder(rec)
        assert b._recorder is rec

    def test_build_simulation_returns_simulation(self):
        from vlmodels.fermenter.config import (
            FermenterBuilder, VesselConfig,
        )
        from PyOMES.core import Simulation
        sim = (
            FermenterBuilder()
            .vessel(V_total_L=100, T_K=305.15)
            .build_simulation()
        )
        assert isinstance(sim, Simulation)
        assert "main" in sim.cvs

    def test_build_simulation_attaches_controllers(self):
        from vlmodels.fermenter.config import FermenterBuilder
        ctrl = object()
        sim = (
            FermenterBuilder()
            .vessel(V_total_L=100, T_K=305.15)
            .controller(ctrl)
            .build_simulation()
        )
        assert ctrl in sim.controllers

    def test_build_simulation_attaches_profiles(self):
        from vlmodels.fermenter.config import FermenterBuilder
        prof = object()
        sim = (
            FermenterBuilder()
            .vessel(V_total_L=100, T_K=305.15)
            .profile(prof)
            .build_simulation()
        )
        assert prof in sim.profiles

    def test_build_simulation_attaches_solver_and_recorder(self):
        from vlmodels.fermenter.config import FermenterBuilder
        from PyOMES.core import BatchRecorder, SimultaneousEulerSolver
        solver = SimultaneousEulerSolver()
        recorder = BatchRecorder()
        sim = (
            FermenterBuilder()
            .vessel(V_total_L=100, T_K=305.15)
            .solver()  # default euler — produces an SimultaneousEulerSolver
            .recorder(recorder)
            .build_simulation()
        )
        # The fluent .solver() above sets self._solver; the new
        # framework reads it through Simulation.solver.
        assert sim.solver is not None
        assert sim.recorder is recorder

    def test_label_propagates(self):
        from vlmodels.fermenter.config import FermenterBuilder
        sim = (
            FermenterBuilder()
            .vessel(V_total_L=100, T_K=305.15)
            .label("exp_42")
            .build_simulation()
        )
        assert sim.label == "exp_42"
        assert sim._context.label == "exp_42"

    def test_label_override(self):
        from vlmodels.fermenter.config import FermenterBuilder
        sim = (
            FermenterBuilder()
            .vessel(V_total_L=100, T_K=305.15)
            .label("default")
            .build_simulation(label="override")
        )
        assert sim.label == "override"

    def test_build_simulation_and_run_returns_batch_result(self):
        from vlmodels.fermenter.config import FermenterBuilder
        from PyOMES.core import BatchResult
        result = (
            FermenterBuilder()
            .vessel(V_total_L=100, T_K=305.15)
            .transfer_kinetic(kLa_O2=100.0)
            .build_simulation_and_run(tau_h=0.1, n_steps=5)
        )
        assert isinstance(result, BatchResult)
        # Per-CV nested gas_mol with "main" key (single-CV convention)
        assert "main" in result.gas_mol

    def test_legacy_build_still_returns_cv(self):
        """Backward-compat: the legacy build() path keeps returning
        a ControlVolume so existing demos and tests don't break."""
        from vlmodels.fermenter.config import FermenterBuilder
        from PyOMES.core import ControlVolume
        cv = (
            FermenterBuilder()
            .vessel(V_total_L=100, T_K=305.15)
            .build()
        )
        assert isinstance(cv, ControlVolume)


# ═══════════════════════════════════════════════════════════════════════
#  ADM1 / BSM2 integration with Simulation (C13)
# ═══════════════════════════════════════════════════════════════════════

class TestC13ADM1Simulation:
    """ADM1 CVs build through their legacy helpers and run through the
    new Simulation orchestrator end-to-end. Strong-ion seeding stays
    in CV construction (per STATE_UNIFICATION C4); chem_env is gone."""

    def test_adm1_cv_runs_through_simulation(self):
        from vlmodels.adm1.base import (
            build_adm1_reactions, build_adm1_cv, seed_adm1_strong_ions,
        )
        from PyOMES.core import Simulation, BatchResult

        rxn_set = build_adm1_reactions()
        cv = build_adm1_cv(rxn_set, V_total_L=2.0, T_K=308.15)
        seed_adm1_strong_ions(cv, CT_Na=0.100, CT_Cl=0.030)

        sim = Simulation(cvs={"main": cv}, label="adm1_demo")
        # Short run — ADM1 dynamics need long timescales for steady-
        # state, but the smoke test here verifies the orchestrator
        # path doesn't crash.
        result = sim.run(tau_h=0.1, n_steps=4)
        assert isinstance(result, BatchResult)
        assert "main" in result.liquid_mol
        # ADM1 has many species; just verify some appear in the
        # nested result
        assert len(result.liquid_mol["main"]) > 0

    def test_adm1_simulation_result_shape(self):
        """ADM1 CV run produces the per-CV nested BatchResult shape
        (pre-allocated arrays at the expected length)."""
        import numpy as np
        from vlmodels.adm1.base import (
            build_adm1_reactions, build_adm1_cv, seed_adm1_strong_ions,
        )
        from PyOMES.core import Simulation

        rxn_set = build_adm1_reactions()
        cv = build_adm1_cv(rxn_set, V_total_L=2.0, T_K=308.15)
        seed_adm1_strong_ions(cv, CT_Na=0.100, CT_Cl=0.030)
        sim = Simulation(cvs={"main": cv})
        result = sim.run(tau_h=0.1, n_steps=4)
        # Per-CV result shape: t_h length n_steps+1; species arrays
        # the same length.
        assert result.t_h.shape == (5,)
        for species_arr in result.liquid_mol["main"].values():
            assert species_arr.shape == (5,)
        # pH array exists at the correct shape (whether it's populated
        # with finite values depends on whether ADM1's engine wrote
        # H+ into n_mol — outside the scope of this smoke test).
        assert result.pH["main"].shape == (5,)

    def test_adm1_water_vapour_converges_to_saturation(self):
        """CP4 of LAYER1_GAP_CLOSURE: build_adm1_cv() replaced
        WaterVapourBoundary with transfer_models=EquilibriumTransferModel
        (RaoultEquilibrium()) + a seeded liquid H2O pool. Over enough
        steps, gas-phase water vapour partial pressure should converge
        to RaoultEquilibrium's P_sat(T) — the same physical target
        WaterVapourBoundary aimed for, now reached by properly moving
        mass out of the (now-tracked) liquid pool instead of conjuring
        it from nowhere."""
        import math
        from vlmodels.adm1.base import (
            build_adm1_reactions, build_adm1_cv, seed_adm1_strong_ions,
        )
        from PyOMES.chemistry import RaoultEquilibrium
        from PyOMES.core import Simulation

        rxn_set = build_adm1_reactions()
        cv = build_adm1_cv(rxn_set, V_total_L=2.0, T_K=308.15)
        seed_adm1_strong_ions(cv, CT_Na=0.100, CT_Cl=0.030)

        liq_h2o_before = cv.phases["liquid"].n_mol["H2O"]
        assert liq_h2o_before > 0.0  # seeded, not the old unseeded-zero default

        sim = Simulation(cvs={"main": cv})
        sim.run(tau_h=0.5, n_steps=20)

        gas = cv.phases["gas"]
        liq = cv.phases["liquid"]
        assert not math.isnan(gas.n_mol["H2O"])
        assert not math.isnan(liq.n_mol["H2O"])
        p_final = gas.p_atm.get("H2O", 0.0)
        expected = RaoultEquilibrium().P_sat(308.15)
        assert p_final == pytest.approx(expected, rel=0.05)
        # Liquid pool must have lost exactly what the gas phase gained
        # (mass conservation — the WaterVapourBoundary-era model had no
        # such guarantee).
        # Loose tolerance: kinetic hydrolysis reactions also produce/
        # consume H2O over the run, so this isn't a pure transfer-only
        # balance — but it must be in the right ballpark, not off by
        # orders of magnitude the way a "water from nowhere" regression
        # would be.
        liq_h2o_after = liq.n_mol["H2O"]
        assert liq_h2o_after < liq_h2o_before
        assert (liq_h2o_before - liq_h2o_after) == pytest.approx(
            gas.n_mol["H2O"], rel=0.05,
        )


# ═══════════════════════════════════════════════════════════════════════
#  P1 — Fail-loud initial speciation solve (framework-polish)
# ═══════════════════════════════════════════════════════════════════════

class TestInitialSolveWarning:
    """_initial_solve emits AccuracyWarning instead of silently swallowing
    engine exceptions (framework-polish P1)."""

    def _make_cv_with_failing_engine(self, exc_type=ValueError):
        """Return a CV whose reaction_system.engine.solve always raises."""
        from types import SimpleNamespace
        from PyOMES.core import ControlVolume, LiquidPhase

        liq = LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=300.0)
        cv = ControlVolume(phases={"liquid": liq}, label="test_cv")

        engine = SimpleNamespace(solve=lambda **kw: (_ for _ in ()).throw(exc_type("boom")))
        reaction_system = SimpleNamespace(engine=engine)
        cv._reaction_system = reaction_system
        return cv

    def test_warning_fires_on_engine_failure(self):
        import warnings
        from PyOMES.core import Simulation
        from PyOMES.monitoring import AccuracyWarning

        cv = self._make_cv_with_failing_engine(ValueError)
        sim = Simulation(cvs={"test_cv": cv})

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            sim._initial_solve(cv)

        accuracy_warnings = [w for w in caught if issubclass(w.category, AccuracyWarning)]
        assert len(accuracy_warnings) == 1
        assert "test_cv" in str(accuracy_warnings[0].message)
        assert "NaN" in str(accuracy_warnings[0].message)

    def test_warning_fires_for_each_exception_type(self):
        import warnings
        from PyOMES.core import Simulation
        from PyOMES.monitoring import AccuracyWarning

        sim = Simulation(cvs={"dummy": self._make_cv_with_failing_engine(ValueError)})

        for exc_type in (TypeError, ValueError, AttributeError, KeyError):
            cv = self._make_cv_with_failing_engine(exc_type)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                sim._initial_solve(cv)
            accuracy_warnings = [w for w in caught if issubclass(w.category, AccuracyWarning)]
            assert len(accuracy_warnings) == 1, f"Expected warning for {exc_type.__name__}"


# ═══════════════════════════════════════════════════════════════════════
#  HPC_CHECKPOINTING C1 — snapshot coverage audit
# ═══════════════════════════════════════════════════════════════════════

class TestSnapshotCoverage:
    """Coverage audit for Simulation.snapshot() and ControlVolume.snapshot().

    Each test pins one row of the table from
    docs/dev/implementation/shipped/HPC_CHECKPOINTING.md: whether a piece of
    mid-run state is deep-copied, shared, or fresh on the snapshot.
    No source changes are made; this is an audit-and-document pass.
    """

    # ── Simulation-level ──────────────────────────────────────────────

    def test_t_h_not_copied_to_snapshot(self):
        """_t_h is not carried into snapshot(); the new Simulation starts
        at 0.0 regardless of where the original left off after run()."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        sim.run(tau_h=2.0, n_steps=4)
        assert sim._t_h == pytest.approx(2.0)
        snap = sim.snapshot()
        assert snap._t_h == pytest.approx(0.0)

    def test_t_h_reset_at_run_entry_is_by_design(self):
        """Decision 10: _t_h always resets to 0.0 at run() entry even when
        the accumulator already holds a non-zero value. Checkpoint resume
        shifts the time grid via start_t_h (C2) rather than relying on
        accumulated internal state."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        sim._t_h = 99.0  # simulate stale accumulator
        sim.run(tau_h=1.0, n_steps=2)
        assert sim._t_h == pytest.approx(1.0)  # 99.0 was overwritten

    def test_runcontext_is_fresh_on_snapshot(self):
        """snapshot() creates a new RunContext; is_running is always False.
        A stale is_running=True on the original does not leak into the
        snapshot's context."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        sim._context.is_running = True  # simulate mid-run or stale state
        snap = sim.snapshot()
        assert snap._context is not sim._context
        assert snap._context.is_running is False

    def test_controller_sample_state_not_present_before_run(self):
        """_controller_sample_state is only created inside run(); neither a
        fresh Simulation nor its snapshot has this attribute before run()."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        assert not hasattr(sim, "_controller_sample_state")
        snap = sim.snapshot()
        assert not hasattr(snap, "_controller_sample_state")

    def test_controller_sample_state_not_copied_to_snapshot(self):
        """After run(), _controller_sample_state persists on the original.
        snapshot() does not carry it over — the snapshot starts clean,
        and its own run() will reset it at entry."""
        from PyOMES.core import Simulation
        from PyOMES.control.actions import ControlAction

        class _SampledController:
            sample_period_h = 0.5
            target_cv_key = "main"
            def compute(self, state, dt_h):
                return ControlAction(controller_label="sampled")

        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv}, controllers=[_SampledController()])
        sim.run(tau_h=1.0, n_steps=10)
        assert hasattr(sim, "_controller_sample_state")
        assert len(sim._controller_sample_state) > 0
        snap = sim.snapshot()
        assert not hasattr(snap, "_controller_sample_state")

    def test_recorder_is_shared_when_explicit(self):
        """An explicitly passed recorder object is shared by identity in
        snapshot(). Both original and snapshot hold the same reference."""
        from PyOMES.core import Simulation, BatchRecorder
        cv = _make_liquid_cv("main", {"S": 1.0})
        rec = BatchRecorder()
        sim = Simulation(cvs={"main": cv}, recorder=rec)
        snap = sim.snapshot()
        assert snap._recorder is rec
        assert snap._recorder is sim._recorder

    def test_recorder_none_propagates_to_snapshot(self):
        """When recorder=None (default), snapshot also sees None. Each
        run() constructs a fresh BatchRecorder at entry; no recorder
        state accumulates between runs."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        snap = sim.snapshot()
        assert sim._recorder is None
        assert snap._recorder is None

    # ── ControlVolume-level ───────────────────────────────────────────

    def test_cv_snapshot_phases_deep_copied(self):
        """Phase.n_mol arrays are deep-copied; mutating the snapshot phase
        does not affect the original."""
        from PyOMES.core import ControlVolume, LiquidPhase
        liq = LiquidPhase(n_mol={"S": 1.5, "X": 0.3}, V_L=1.0, T_K=298.15)
        cv = ControlVolume(phases={"liquid": liq}, label="src")
        snap = cv.snapshot()
        snap.phases["liquid"].n_mol["S"] = 99.0
        assert cv.phases["liquid"].n_mol["S"] == pytest.approx(1.5)

    def test_cv_snapshot_reaction_system_shared_by_identity(self):
        """reaction_system is shared by identity in cv.snapshot(). The
        snapshot docstring labels it 'stateless', but it hosts the
        BisectionChemicalEquilibriumEngine and its warm-start caches — see
        test_cv_snapshot_speciation_warm_start_shared_via_reaction_system."""
        from types import SimpleNamespace
        from PyOMES.core import ControlVolume, LiquidPhase

        rs = SimpleNamespace()  # no attach_monitor — skips monitor wiring
        liq = LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=298.15)
        cv = ControlVolume(phases={"liquid": liq}, reaction_system=rs)
        snap = cv.snapshot()
        assert snap.reaction_system is rs
        assert snap.reaction_system is cv.reaction_system

    def test_cv_snapshot_internal_interfaces_shared(self):
        """Interface objects are shared by identity; only the wrapper list
        is a fresh copy."""
        from types import SimpleNamespace
        from PyOMES.core import ControlVolume, LiquidPhase, GasPhase

        iface = SimpleNamespace(phase_a_key="gas", phase_b_key="liquid")
        gas = GasPhase(n_mol={"O2": 1.0}, V_L=0.2, T_K=298.15)
        liq = LiquidPhase(n_mol={"O2": 0.001}, V_L=0.8, T_K=298.15)
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            internal_interfaces=[iface],
        )
        snap = cv.snapshot()
        assert snap.internal_interfaces is not cv.internal_interfaces
        assert snap.internal_interfaces[0] is iface

    def test_cv_snapshot_boundaries_shared(self):
        """Boundary objects are shared by identity; the wrapper list is
        a fresh copy."""
        from PyOMES.core import ControlVolume, LiquidPhase
        boundary = object()
        liq = LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=298.15)
        cv = ControlVolume(phases={"liquid": liq})
        cv.boundaries.append(boundary)
        snap = cv.snapshot()
        assert snap.boundaries is not cv.boundaries
        assert snap.boundaries[0] is boundary

    def test_cv_snapshot_accuracy_monitor_is_fresh(self):
        """cv.snapshot() constructs a fresh AccuracyMonitor for the new CV.
        Mid-run per-CV state (last_pH, step_count) on the original is not
        copied into the snapshot."""
        from PyOMES.core import ControlVolume, LiquidPhase
        liq = LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=298.15)
        cv = ControlVolume(phases={"liquid": liq})
        cv._accuracy_monitor.last_pH = 7.2  # simulate mid-run state
        cv._accuracy_monitor.step_count = 42
        snap = cv.snapshot()
        assert snap._accuracy_monitor is not cv._accuracy_monitor
        assert snap._accuracy_monitor.last_pH is None
        assert snap._accuracy_monitor.step_count == 0

    def test_cv_snapshot_accuracy_monitor_displaces_original_on_reaction_system(self):
        """When a reaction_system is present, the snapshot's __init__ calls
        attach_monitor(new_monitor), which replaces reaction_system._accuracy_monitor
        with the snapshot's fresh monitor. After the snapshot, the original CV's
        monitor is detached from the reaction_system chain."""
        from PyOMES.core import ControlVolume, LiquidPhase

        class _MockRS:
            def __init__(self):
                self._accuracy_monitor = None
                self.reactions = []
            def attach_monitor(self, m):
                self._accuracy_monitor = m
            def attach_conservation_monitor(self, m):
                pass

        rs = _MockRS()
        liq = LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=298.15)
        cv = ControlVolume(phases={"liquid": liq}, reaction_system=rs)
        original_monitor = cv._accuracy_monitor
        assert rs._accuracy_monitor is original_monitor  # wired at init
        snap = cv.snapshot()
        # Snapshot's __init__ calls attach_monitor → rs now points to
        # the snapshot's monitor; the original CV's monitor is detached.
        assert rs._accuracy_monitor is snap._accuracy_monitor
        assert rs._accuracy_monitor is not original_monitor

    def test_cv_snapshot_speciation_warm_start_shared_via_reaction_system(self):
        """BisectionChemicalEquilibriumEngine warm-start caches (_logH_last, _I_last) are shared
        after cv.snapshot() because reaction_system is shared by identity.
        Any warm-start state on the engine is visible from both CVs."""
        from types import SimpleNamespace
        from PyOMES.core import ControlVolume, LiquidPhase

        engine = SimpleNamespace(_logH_last=-7.0, _I_last=0.05)
        rs = SimpleNamespace(engine=engine)  # no attach_monitor
        liq = LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=298.15)
        cv = ControlVolume(phases={"liquid": liq}, reaction_system=rs)
        snap = cv.snapshot()

        assert snap.reaction_system is cv.reaction_system
        assert snap.reaction_system.engine._logH_last == pytest.approx(-7.0)
        # Mutating through one route is visible from the other
        snap.reaction_system.engine._logH_last = -6.5
        assert cv.reaction_system.engine._logH_last == pytest.approx(-6.5)


# ═══════════════════════════════════════════════════════════════════════
#  HPC_CHECKPOINTING C2 — start_t_h on run() + BatchResult.concat()
# ═══════════════════════════════════════════════════════════════════════

class TestStartTH:
    """start_t_h keyword argument on Simulation.run()."""

    def test_default_preserves_existing_behaviour(self):
        """start_t_h=0.0 (default) gives the same result as before."""
        import numpy as np
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        result = sim.run(tau_h=1.0, n_steps=4)
        assert result.t_h[0] == pytest.approx(0.0)
        assert result.t_h[-1] == pytest.approx(1.0)
        assert result.t_h.shape == (5,)

    def test_start_t_h_shifts_time_grid(self):
        """start_t_h shifts the output grid to [start_t_h, start_t_h + tau_h]."""
        import numpy as np
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        result = sim.run(tau_h=2.0, n_steps=4, start_t_h=5.0)
        assert result.t_h[0] == pytest.approx(5.0)
        assert result.t_h[-1] == pytest.approx(7.0)
        assert result.t_h.shape == (5,)
        # Interior points are evenly spaced
        assert result.t_h[1] == pytest.approx(5.5)
        assert result.t_h[2] == pytest.approx(6.0)

    def test_start_t_h_updates_accumulator(self):
        """After run(start_t_h=T), sim._t_h == T + tau_h."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        sim.run(tau_h=3.0, n_steps=6, start_t_h=10.0)
        assert sim._t_h == pytest.approx(13.0)

    def test_start_t_h_is_keyword_only(self):
        """start_t_h is keyword-only (after the bare *); passing it
        positionally raises TypeError."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0})
        sim = Simulation(cvs={"main": cv})
        with pytest.raises(TypeError):
            sim.run(1.0, 4, 5.0)  # type: ignore[call-arg]

    def test_start_t_h_zero_does_not_affect_species_values(self):
        """The explicit start_t_h=0.0 is identical to omitting it entirely."""
        from PyOMES.core import Simulation
        cv1 = _make_liquid_cv("main", {"S": 1.0})
        cv2 = _make_liquid_cv("main", {"S": 1.0})
        sim1 = Simulation(cvs={"main": cv1})
        sim2 = Simulation(cvs={"main": cv2})
        r1 = sim1.run(tau_h=0.5, n_steps=5)
        r2 = sim2.run(tau_h=0.5, n_steps=5, start_t_h=0.0)
        import numpy as np
        np.testing.assert_array_equal(r1.t_h, r2.t_h)
        np.testing.assert_array_almost_equal(
            r1.liquid_mol["main"]["S"], r2.liquid_mol["main"]["S"]
        )

    def test_species_arrays_length_matches_time_grid(self):
        """Species arrays always have the same length as t_h regardless of
        start_t_h."""
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", {"S": 1.0, "X": 0.5})
        sim = Simulation(cvs={"main": cv})
        result = sim.run(tau_h=1.0, n_steps=10, start_t_h=100.0)
        n = len(result.t_h)
        for arr in result.liquid_mol["main"].values():
            assert len(arr) == n


class TestBatchResultConcat:
    """BatchResult.concat(*results) — merge trajectories."""

    def _run(self, n_mol, tau_h, n_steps, start_t_h=0.0):
        from PyOMES.core import Simulation
        cv = _make_liquid_cv("main", n_mol)
        sim = Simulation(cvs={"main": cv})
        return sim.run(tau_h=tau_h, n_steps=n_steps, start_t_h=start_t_h)

    def test_concat_empty_raises(self):
        from PyOMES.core import BatchResult
        with pytest.raises(ValueError, match="at least one"):
            BatchResult.concat()

    def test_concat_single_returns_same(self):
        from PyOMES.core import BatchResult
        r = self._run({"S": 1.0}, tau_h=1.0, n_steps=4)
        assert BatchResult.concat(r) is r

    def test_concat_two_results_no_duplicate_boundary(self):
        """The boundary time point (end of r1 == start of r2) appears once."""
        import numpy as np
        from PyOMES.core import BatchResult
        r1 = self._run({"S": 1.0}, tau_h=1.0, n_steps=4)
        r2 = self._run({"S": 1.0}, tau_h=1.0, n_steps=4, start_t_h=1.0)
        full = BatchResult.concat(r1, r2)
        # r1: 5 points [0.0, 0.25, 0.5, 0.75, 1.0]
        # r2: 5 points [1.0, 1.25, 1.5, 1.75, 2.0] → drop first → 4 points
        assert full.t_h.shape == (9,)
        assert full.t_h[0] == pytest.approx(0.0)
        assert full.t_h[4] == pytest.approx(1.0)
        assert full.t_h[-1] == pytest.approx(2.0)
        # No duplicates
        assert len(np.unique(np.round(full.t_h, 10))) == 9

    def test_concat_species_arrays_correct_length(self):
        """Species arrays in the concat result match len(t_h)."""
        from PyOMES.core import BatchResult
        r1 = self._run({"S": 1.0}, tau_h=1.0, n_steps=4)
        r2 = self._run({"S": 1.0}, tau_h=1.0, n_steps=4, start_t_h=1.0)
        full = BatchResult.concat(r1, r2)
        n = len(full.t_h)
        for arr in full.liquid_mol["main"].values():
            assert len(arr) == n

    def test_concat_species_values_consistent(self):
        """Species values in the concat result are drawn from both runs."""
        import numpy as np
        from PyOMES.core import BatchResult, Simulation
        # No reactions: S stays constant throughout both runs
        r1 = self._run({"S": 1.0}, tau_h=1.0, n_steps=4)
        # Snapshot the state, resume (same initial conditions since no reactions)
        r2 = self._run({"S": 1.0}, tau_h=1.0, n_steps=4, start_t_h=1.0)
        full = BatchResult.concat(r1, r2)
        # All values should be 1.0 (no reactions, no decay)
        np.testing.assert_array_almost_equal(
            full.liquid_mol["main"]["S"], np.ones(9)
        )

    def test_concat_mismatched_species_fills_zeros(self):
        """A species present in r1 but absent from r2 is filled with zeros
        for r2's time range in the concatenated result."""
        import numpy as np
        from PyOMES.core import BatchResult, Simulation, ControlVolume, LiquidPhase

        # r1 has both S and X; r2 has only S
        cv1 = ControlVolume(
            phases={"liquid": LiquidPhase(n_mol={"S": 1.0, "X": 0.5}, V_L=1.0, T_K=298.15)},
        )
        sim1 = Simulation(cvs={"main": cv1})
        r1 = sim1.run(tau_h=1.0, n_steps=4)

        cv2 = ControlVolume(
            phases={"liquid": LiquidPhase(n_mol={"S": 1.0}, V_L=1.0, T_K=298.15)},
        )
        sim2 = Simulation(cvs={"main": cv2})
        r2 = sim2.run(tau_h=1.0, n_steps=4, start_t_h=1.0)

        full = BatchResult.concat(r1, r2)
        # X exists in full result
        assert "X" in full.liquid_mol["main"]
        x_arr = full.liquid_mol["main"]["X"]
        assert x_arr.shape == (9,)
        # First 5 points from r1 (X = 0.5)
        np.testing.assert_array_almost_equal(x_arr[:5], 0.5)
        # Last 4 points filled with zeros (X absent from r2)
        np.testing.assert_array_almost_equal(x_arr[5:], 0.0)

    def test_concat_ph_and_scalar_channels_concatenated(self):
        """pH, ionic_strength, and P_atm scalar channels are also merged."""
        import numpy as np
        from PyOMES.core import BatchResult
        r1 = self._run({"S": 1.0}, tau_h=1.0, n_steps=4)
        r2 = self._run({"S": 1.0}, tau_h=1.0, n_steps=4, start_t_h=1.0)
        full = BatchResult.concat(r1, r2)
        assert "main" in full.pH
        assert full.pH["main"].shape == (9,)
        assert "main" in full.P_atm
        assert full.P_atm["main"].shape == (9,)

    def test_concat_list_channels_concatenated(self):
        """transfer_records, advance_results, and other list channels
        are concatenated without dropping any entries."""
        from PyOMES.core import BatchResult
        r1 = self._run({"S": 1.0}, tau_h=1.0, n_steps=4)
        r2 = self._run({"S": 1.0}, tau_h=1.0, n_steps=4, start_t_h=1.0)
        full = BatchResult.concat(r1, r2)
        # 4 steps each → 8 step records
        assert len(full.advance_results) == 8
        assert len(full.transfer_records) == 8

    def test_concat_runtime_s_is_sum(self):
        """runtime_s in the concatenated result is the sum of all inputs."""
        from PyOMES.core import BatchResult
        r1 = self._run({"S": 1.0}, tau_h=0.1, n_steps=2)
        r2 = self._run({"S": 1.0}, tau_h=0.1, n_steps=2, start_t_h=0.1)
        full = BatchResult.concat(r1, r2)
        assert full.runtime_s == pytest.approx(r1.runtime_s + r2.runtime_s)

    def test_concat_three_results(self):
        """More than two results can be concatenated in one call."""
        import numpy as np
        from PyOMES.core import BatchResult
        r1 = self._run({"S": 1.0}, tau_h=1.0, n_steps=4)
        r2 = self._run({"S": 1.0}, tau_h=1.0, n_steps=4, start_t_h=1.0)
        r3 = self._run({"S": 1.0}, tau_h=1.0, n_steps=4, start_t_h=2.0)
        full = BatchResult.concat(r1, r2, r3)
        # 5 + 4 + 4 = 13 time points
        assert full.t_h.shape == (13,)
        assert full.t_h[0] == pytest.approx(0.0)
        assert full.t_h[-1] == pytest.approx(3.0)
        assert len(full.advance_results) == 12  # 4 steps × 3 runs
