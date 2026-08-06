# -*- coding: utf-8 -*-
"""Tests for ControlVolume.advance() and related helpers (Stage 3).

Validates:
- advance() with no reaction model (equilibrium + properties only)
- advance() with a simple KineticReaction attached
- advance() with external_source_terms
- advance() with both reaction model and external source terms
- Mass balance conservation through advance()
- Non-volatile species untouched by equilibrium
- AdvanceResult fields populated correctly
- Property solver results accessible via AdvanceResult.properties
- species_concentration convenience accessor
- _build_reaction_environment correctness
- Reaction integration produces correct mole changes
"""

import numpy as np
import pytest

from PyOMES.chemistry import HenryPartition


def _hp(kH: float, dlnH: float = 0.0) -> HenryPartition:
    return HenryPartition(H_ref=kH * 1000.0 / 101325.0, dlnH=dlnH)


# ═══════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════

def _make_simple_cv():
    """CV with gas + liquid, no interfaces, no reactions."""
    from PyOMES.core import ControlVolume, GasPhase, LiquidPhase
    gas = GasPhase(n_mol={"O2": 0.5, "CO2": 0.1, "N2": 2.0}, V_L=0.4, T_K=305.15)
    liquid = LiquidPhase(n_mol={"CO2": 0.01, "O2": 0.001}, V_L=1.6, T_K=305.15)
    return ControlVolume(phases={"gas": gas, "liquid": liquid}, label="test")


def _make_cv_with_speciation():
    """CV with speciation. state-unification C4: equilibria attach
    via ReactionSystem (chem_env / SpeciationPropertySolver gone).
    """
    from PyOMES.core import ControlVolume, GasPhase, LiquidPhase
    from PyOMES.reactions import ReactionSystem

    gas = GasPhase(n_mol={"O2": 0.5, "CO2": 0.1, "N2": 2.0}, V_L=0.4, T_K=305.15)
    liquid = LiquidPhase(
        n_mol={"CO2": 0.005, "O2": 0.001, "AceticAcid": 0.005 * 1.6},
        V_L=1.6, T_K=305.15,
    )
    equilibria = _make_speciation_reactions()
    system = ReactionSystem(equilibria, label="test_spec_chem")

    return ControlVolume(
        phases={"gas": gas, "liquid": liquid},
        reaction_system=system,
        label="test_spec",
    )


def _make_speciation_reactions():
    """Build a small acetate + carbonate equilibrium reaction list
    for tests. The engine is built lazily by ReactionSystem.engine
    on first access (state-unification C4)."""
    from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
    from PyOMES.chemistry import Species
    from PyOMES.chemistry.common_species import (
        H_plus, OH_minus, H2O, CO2, HCO3_minus,
    )

    HAc = Species(id="AceticAcid", atoms={"C": 2, "H": 4, "O": 2}, charge=0, MW=60.052)
    Ac_minus = Species(id="AceticAcid-", atoms={"C": 2, "H": 3, "O": 2}, charge=-1)
    rxns = [
        EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=H2O,      phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=H_plus,   phase="liquid", coefficient=+1.0),
                StoichiometryEntry(species=OH_minus, phase="liquid", coefficient=+1.0),
            ],
            log_K=-14.0,
            balance_elements=("H", "O"), label="water"),
        EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=CO2, phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=H2O, phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=HCO3_minus, phase="liquid", coefficient=+1.0),
                StoichiometryEntry(species=H_plus, phase="liquid", coefficient=+1.0),
            ],
            log_K=-6.35,
            balance_elements=("C", "H", "O"), label="CO2"),
        EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=HAc, phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=Ac_minus, phase="liquid", coefficient=+1.0),
                StoichiometryEntry(species=H_plus, phase="liquid", coefficient=+1.0),
            ],
            log_K=-4.76,
            balance_elements=("C", "H", "O"), label="acetate"),
    ]
    return rxns


def _make_simple_reaction():
    """A → B in liquid, rate = 0.5 * [A] * V_L mol/h."""
    from PyOMES.reactions import KineticReaction, StoichiometryEntry
    from PyOMES.chemistry import Species
    A = Species(id="A", atoms={"C": 1})
    B = Species(id="B", atoms={"C": 1})
    return KineticReaction(
        stoichiometry=[
            StoichiometryEntry(species=A, phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=B, phase="liquid", coefficient=+1.0),
        ],
        rate_fn=lambda env: 0.5 * env.S("A") * env.V_L,
        balance_elements=("C",),
        label="A_to_B",
    )


def _make_cv_with_reaction():
    """CV with a simple A→B reaction in liquid, no equilibrium interfaces."""
    from PyOMES.core import ControlVolume, GasPhase, LiquidPhase

    gas = GasPhase(n_mol={"O2": 0.5}, V_L=0.4, T_K=305.15)
    liquid = LiquidPhase(n_mol={"A": 0.1, "B": 0.0}, V_L=1.0, T_K=305.15)
    rxn = _make_simple_reaction()

    return ControlVolume(
        phases={"gas": gas, "liquid": liquid},
        reaction_system=rxn,
        label="test_rxn",
    )


# ═══════════════════════════════════════════════════════════════════════
#  PropertyCalculator protocol (state-unification C5)
# ═══════════════════════════════════════════════════════════════════════

class TestPropertyCalculatorSeam:
    """The PropertyCalculator protocol provides a seam for scalar
    derived properties (viscosity, density, …). state-unification
    C5 introduces the protocol + cv.property_calculators list +
    invocation point in advance(). No concrete implementers yet —
    these tests validate the plumbing only."""

    def test_default_no_calculators(self):
        cv = _make_simple_cv()
        assert cv.property_calculators == []

    def test_protocol_structural_typing(self):
        from PyOMES.core import PropertyCalculator

        class FakeViscosity:
            key = "viscosity"
            phase_key = "liquid"
            def compute(self, phase, T_K, P_atm):
                return 0.001

        assert isinstance(FakeViscosity(), PropertyCalculator)

    def test_calculator_writes_phase_properties(self):
        """A calculator's compute() return lands on
        phase.properties[calc.key] after advance()."""
        class FakeViscosity:
            key = "viscosity"
            phase_key = "liquid"
            def compute(self, phase, T_K, P_atm):
                return 42.0

        cv = _make_cv_with_reaction()
        cv.property_calculators.append(FakeViscosity())
        cv.advance(0.01, 0.0)
        assert cv.phases["liquid"].properties["viscosity"] == pytest.approx(42.0)

    def test_calculator_visible_to_rate_law(self):
        """A rate law can read phase.properties via env.prop()."""
        captured = []

        class FakeViscosity:
            key = "viscosity"
            phase_key = "liquid"
            def compute(self, phase, T_K, P_atm):
                return 1.234

        from PyOMES.reactions import KineticReaction, StoichiometryEntry
        from PyOMES.chemistry import Species

        A = Species(id="A", atoms={"C": 1})
        B = Species(id="B", atoms={"C": 1})

        def rate_fn(env):
            captured.append(env.prop("viscosity"))
            return 0.0

        rxn = KineticReaction(
            stoichiometry=[
                StoichiometryEntry(species=A, phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=B, phase="liquid", coefficient=+1.0),
            ],
            rate_fn=rate_fn,
            balance_elements=("C",),
        )

        from PyOMES.core import ControlVolume, GasPhase, LiquidPhase
        cv = ControlVolume(
            phases={
                "gas": GasPhase(n_mol={"O2": 0.5}, V_L=0.4, T_K=305.15),
                "liquid": LiquidPhase(n_mol={"A": 0.1, "B": 0.0}, V_L=1.0, T_K=305.15),
            },
            reaction_system=rxn,
            property_calculators=[FakeViscosity()],
        )
        cv.advance(0.01, 0.0)
        # rate_fn captured the viscosity from env.prop()
        assert any(v == pytest.approx(1.234) for v in captured)

    def test_gas_phase_calculator_writes_gas_properties(self):
        """A calculator with phase_key='gas' writes to gas.properties."""
        from PyOMES.core import ControlVolume, GasPhase, LiquidPhase

        class GasDensityCalc:
            key = "density_kg_m3"
            phase_key = "gas"
            def compute(self, phase, T_K, P_atm):
                return 1.25

        cv = ControlVolume(
            phases={
                "gas": GasPhase(n_mol={"O2": 1.0}, V_L=0.5, T_K=298.15),
                "liquid": LiquidPhase(n_mol={"S": 0.1}, V_L=1.0, T_K=298.15),
            },
            property_calculators=[GasDensityCalc()],
        )
        cv.advance(0.01, 0.0)
        assert cv.phases["gas"].properties["density_kg_m3"] == pytest.approx(1.25)
        assert "density_kg_m3" not in cv.phases["liquid"].properties

    def test_absent_phase_calculator_is_skipped(self):
        """Calculator whose phase_key names an absent phase is silently skipped."""
        from PyOMES.core import ControlVolume, LiquidPhase

        class SolidCalc:
            key = "porosity"
            phase_key = "solid"
            def compute(self, phase, T_K, P_atm):
                raise AssertionError("should not be called")

        cv = ControlVolume(
            phases={"liquid": LiquidPhase(n_mol={"S": 0.1}, V_L=1.0, T_K=298.15)},
            property_calculators=[SolidCalc()],
        )
        cv.advance(0.01, 0.0)  # must not raise


# ═══════════════════════════════════════════════════════════════════════
#  AdvanceResult
# ═══════════════════════════════════════════════════════════════════════

class TestAdvanceResult:
    """state-unification C4: AdvanceResult.properties field removed.
    pH and derived chemistry state live on cv.phases[...] directly."""

    def test_default_construction(self):
        from PyOMES.core.interfaces import AdvanceResult
        r = AdvanceResult()
        assert r.transfer is None
        assert r.reaction_sources is None
        assert r.transfer_record is None
        assert r.boundary_records == []

    def test_with_values(self):
        from PyOMES.core.interfaces import AdvanceResult, TransferDiagnostics
        td = TransferDiagnostics()
        r = AdvanceResult(
            transfer=td,
            reaction_sources={"liquid": {"A": -1.0}},
        )
        assert r.transfer is td
        assert r.reaction_sources == {"liquid": {"A": -1.0}}
        assert r.reaction_sources["liquid"]["A"] == -1.0


# ═══════════════════════════════════════════════════════════════════════
#  advance() — no reaction model
# ═══════════════════════════════════════════════════════════════════════

class TestAdvanceNoReaction:

    def test_advance_returns_advance_result(self):
        from PyOMES.core.interfaces import AdvanceResult
        cv = _make_simple_cv()
        result = cv.advance(dt_h=0.01)
        assert isinstance(result, AdvanceResult)

    def test_advance_transfer_diagnostics_present(self):
        cv = _make_simple_cv()
        result = cv.advance(dt_h=0.01)
        assert result.transfer is not None
        # No interfaces → no-op transfer → conserved
        assert result.transfer.is_conserved()

    def test_advance_no_reaction_sources(self):
        cv = _make_simple_cv()
        result = cv.advance(dt_h=0.01)
        assert result.reaction_sources is None

    def test_advance_with_speciation(self):
        """state-unification C4: pH lives on phase.pH after advance,
        derived from n_mol["H+"] populated by the engine writeback."""
        cv = _make_cv_with_speciation()
        cv.advance(0.01, 0.0)
        liq = cv.phases["liquid"]
        assert "H+" in liq.n_mol
        pH = liq.pH
        assert 2.0 < pH < 10.0

    def test_advance_populates_H_plus_in_n_mol(self):
        cv = _make_cv_with_speciation()
        cv.advance(0.01, 0.0)
        assert "H+" in cv.phases["liquid"].n_mol

    def test_advance_speciation_reflects_pre_step_state(self):
        """Speciation runs on the *pre-step* state per the C4 sequential
        body: speciation → external sources → reactions → transfer.
        External sources change moles AFTER speciation runs, so
        injecting a large CO2 feed during the step should not move pH
        within the same advance() call."""
        cv = _make_cv_with_speciation()
        cv.advance(0.01, 0.0, external_source_terms={"liquid": {"CO2": 100.0}})
        pH_with_feed = cv.phases["liquid"].pH

        cv_clean = _make_cv_with_speciation()
        cv_clean.advance(0.01, 0.0)
        pH_without_feed = cv_clean.phases["liquid"].pH
        assert pH_with_feed == pytest.approx(pH_without_feed, abs=1e-9)

    def test_repeated_advance_updates_pH_each_step(self):
        """EQUILIBRIUM_RESULT CP7 regression: pH must update on every
        advance() call, not just the first.

        Guards against a missed ``result.apply_to_phases(...)`` call at
        any of the CP6 production call sites — if the commit were
        silently skipped, n_mol["H+"] would freeze at its first-step
        value and pH would stay constant across steps even as a large
        acid feed is injected between them.
        """
        cv = _make_cv_with_speciation()
        cv.advance(0.01, 0.0)
        pH_step1 = cv.phases["liquid"].pH

        cv.advance(0.01, 0.01, external_source_terms={"liquid": {"AceticAcid": 0.05 * 1.6}})
        pH_step2 = cv.phases["liquid"].pH

        cv.advance(0.01, 0.02)
        pH_step3 = cv.phases["liquid"].pH

        # A large acid feed between steps 1 and 2 must shift pH downward
        # by step 3 (once speciation re-solves on the fed totals) — this
        # can only happen if apply_to_phases() actually committed the
        # step-3 solve.
        assert pH_step3 < pH_step2 - 0.005, (
            f"pH did not update across steps (frozen writeback?): "
            f"step1={pH_step1}, step2={pH_step2}, step3={pH_step3}"
        )


# ═══════════════════════════════════════════════════════════════════════
#  advance() — with external source terms
# ═══════════════════════════════════════════════════════════════════════

class TestAdvanceExternalSources:

    def test_external_source_adds_moles(self):
        cv = _make_simple_cv()
        mol_before = cv.total_mol()
        # Add 0.05 mol/h of CO2 to liquid over dt=1.0 h → +0.05 mol
        cv.advance(dt_h=1.0, external_source_terms={"liquid": {"CO2": 0.05}})
        mol_after = cv.total_mol()
        delta_CO2 = mol_after.get("CO2", 0.0) - mol_before.get("CO2", 0.0)
        assert delta_CO2 == pytest.approx(0.05, abs=1e-12)

    def test_external_source_unknown_phase_ignored(self):
        """Source terms for a non-existent phase should be silently skipped."""
        cv = _make_simple_cv()
        mol_before = cv.total_mol()
        cv.advance(dt_h=1.0, external_source_terms={"solid": {"CaCO3": 1.0}})
        mol_after = cv.total_mol()
        # Nothing should have changed
        for sp in mol_before:
            assert mol_after.get(sp, 0.0) == pytest.approx(mol_before[sp], abs=1e-15)


# ═══════════════════════════════════════════════════════════════════════
#  advance() — with reaction model
# ═══════════════════════════════════════════════════════════════════════

class TestAdvanceWithReaction:

    def test_reaction_consumes_substrate(self):
        cv = _make_cv_with_reaction()
        A_before = cv["liquid"].n_mol.get("A", 0.0)
        assert A_before > 0.0
        cv.advance(dt_h=0.1)
        A_after = cv["liquid"].n_mol.get("A", 0.0)
        assert A_after < A_before

    def test_reaction_produces_product(self):
        cv = _make_cv_with_reaction()
        B_before = cv["liquid"].n_mol.get("B", 0.0)
        cv.advance(dt_h=0.1)
        B_after = cv["liquid"].n_mol.get("B", 0.0)
        assert B_after > B_before

    def test_reaction_conserves_carbon(self):
        """A→B with equal C atoms: total C (A+B) must be conserved."""
        cv = _make_cv_with_reaction()
        total_before = cv["liquid"].n_mol.get("A", 0.0) + cv["liquid"].n_mol.get("B", 0.0)
        cv.advance(dt_h=0.1)
        total_after = cv["liquid"].n_mol.get("A", 0.0) + cv["liquid"].n_mol.get("B", 0.0)
        assert total_after == pytest.approx(total_before, abs=1e-12)

    def test_reaction_sources_in_result(self):
        cv = _make_cv_with_reaction()
        result = cv.advance(dt_h=0.1)
        assert result.reaction_sources is not None
        assert "liquid" in result.reaction_sources
        # A should be consumed (negative rate)
        assert result.reaction_sources["liquid"]["A"] < 0.0
        # B should be produced (positive rate)
        assert result.reaction_sources["liquid"]["B"] > 0.0

    def test_reaction_plus_external_both_applied(self):
        """Both reaction and external source terms should be applied."""
        cv = _make_cv_with_reaction()
        # Add some extra B from external source
        result = cv.advance(dt_h=0.1, external_source_terms={"liquid": {"B": 1.0}})
        # B should have gained from both reaction and external
        B_final = cv["liquid"].n_mol.get("B", 0.0)
        # External adds 1.0 mol/h * 0.1 h = 0.1 mol, plus some from reaction
        assert B_final > 0.1

    def test_reaction_non_reactive_species_untouched(self):
        """Species not in the reaction (O2 in gas) must be unchanged."""
        cv = _make_cv_with_reaction()
        O2_before = cv["gas"].n_mol.get("O2", 0.0)
        cv.advance(dt_h=0.1)
        O2_after = cv["gas"].n_mol.get("O2", 0.0)
        assert O2_after == pytest.approx(O2_before, abs=1e-15)


# ═══════════════════════════════════════════════════════════════════════
#  _build_reaction_environment
# ═══════════════════════════════════════════════════════════════════════

class TestBuildReactionEnvironment:
    """state-unification C4: _build_reaction_environment(t_h).
    Reads pH from phase.pH (raises if no H+ in n_mol; the
    builder catches that and returns env.pH=None). Concentrations
    come from phase.n_mol; t_h is passed directly.
    """

    def test_concentrations_from_liquid(self):
        cv = _make_cv_with_reaction()
        env = cv._build_reaction_environment()
        # A = 0.1 mol in 1.0 L → 0.1 mol/L
        assert env.S("A") == pytest.approx(0.1, abs=1e-12)

    def test_temperature_from_phase(self):
        cv = _make_cv_with_reaction()
        env = cv._build_reaction_environment()
        assert env.T_K == pytest.approx(305.15, abs=0.01)

    def test_volume_from_liquid(self):
        cv = _make_cv_with_reaction()
        env = cv._build_reaction_environment()
        assert env.V_L == pytest.approx(1.0, abs=1e-12)

    def test_pH_from_n_mol_H_plus(self):
        """pH derives from phase.n_mol["H+"] / V_L."""
        cv = _make_cv_with_reaction()
        # Inject H+ matching pH 6.5 directly into n_mol
        H_mol_L = 10 ** -6.5
        cv.phases["liquid"]._n_mol["H+"] = H_mol_L * cv.phases["liquid"].V_L
        env = cv._build_reaction_environment()
        assert env.pH == pytest.approx(6.5)

    def test_pH_none_when_no_H_plus(self):
        """No H+ in n_mol → env.pH is None (the builder catches the
        ValueError that phase.pH raises)."""
        cv = _make_cv_with_reaction()
        env = cv._build_reaction_environment()
        assert env.pH is None

    def test_t_h_passed_through(self):
        cv = _make_cv_with_reaction()
        env = cv._build_reaction_environment(2.5)
        assert env.t_h == pytest.approx(2.5)

    def test_gas_only_cv_uses_gas_volume(self):
        """P2: gas-only CV reads V_L from gas phase, not the 1.0 placeholder."""
        from PyOMES.core import ControlVolume, GasPhase
        gas = GasPhase(n_mol={"O2": 0.5, "N2": 2.0}, V_L=3.0, T_K=310.0)
        cv = ControlVolume(phases={"gas": gas}, label="gas_only")
        env = cv._build_reaction_environment()
        assert env.V_L == pytest.approx(3.0)
        assert env.T_K == pytest.approx(310.0)
        assert env.S("O2") == pytest.approx(0.0)  # no liquid — concentrations empty

    def test_no_phases_raises(self):
        """P2: CV with no phases raises ValueError — no volume to build env from."""
        from PyOMES.core import ControlVolume
        cv = ControlVolume(phases={}, label="empty_cv")
        with pytest.raises(ValueError, match="empty_cv"):
            cv._build_reaction_environment()


# ═══════════════════════════════════════════════════════════════════════
#  species_concentration accessor
# ═══════════════════════════════════════════════════════════════════════

class TestSpeciesConcentration:

    def test_returns_correct_concentration(self):
        cv = _make_cv_with_reaction()
        # A = 0.1 mol in 1.0 L → 0.1 mol/L
        assert cv.species_concentration("liquid", "A") == pytest.approx(0.1, abs=1e-12)

    def test_missing_species_returns_default(self):
        cv = _make_cv_with_reaction()
        assert cv.species_concentration("liquid", "Z") == 0.0

    def test_missing_phase_returns_default(self):
        cv = _make_cv_with_reaction()
        assert cv.species_concentration("solid", "A") == 0.0

    def test_custom_default(self):
        cv = _make_cv_with_reaction()
        assert cv.species_concentration("solid", "A", default=-1.0) == -1.0


# ═══════════════════════════════════════════════════════════════════════
#  P3 — species-vector-map cache (framework-polish)
# ═══════════════════════════════════════════════════════════════════════

class TestSpeciesVectorMapCache:
    """_build_species_vector_map is memoised per reaction_system identity;
    compute_rates should be called only once across multiple advance() steps
    (framework-polish P3)."""

    def test_single_compute_rates_call_across_steps(self):
        """compute_rates is called only once for map discovery, regardless
        of how many advance() steps are taken."""
        from unittest.mock import patch
        cv = _make_cv_with_reaction()

        call_count = {"n": 0}
        original_compute_rates = cv.reaction_system.compute_rates

        def counting_compute_rates(env):
            call_count["n"] += 1
            return original_compute_rates(env)

        with patch.object(cv.reaction_system, "compute_rates", side_effect=counting_compute_rates):
            # First call builds and caches the map
            cv._build_species_vector_map()
            # Subsequent calls should hit the cache, not call compute_rates again
            cv._build_species_vector_map()
            cv._build_species_vector_map()

        assert call_count["n"] == 1, (
            f"Expected 1 compute_rates call for map discovery, got {call_count['n']}"
        )

    def test_cache_cleared_on_reaction_system_reassignment(self):
        """Assigning a new reaction_system invalidates the cache."""
        cv = _make_cv_with_reaction()

        # Warm the cache
        cv._build_species_vector_map()
        assert cv._species_vector_map_cache is not None

        # Reassign via the unchecked path (simulates orchestrator behaviour)
        new_rs = cv.reaction_system  # same object — just testing invalidation
        cv._set_reaction_system_unchecked(new_rs)
        assert cv._species_vector_map_cache is None


# ═══════════════════════════════════════════════════════════════════════
#  Mass balance through advance()
# ═══════════════════════════════════════════════════════════════════════

class TestAdvanceMassBalance:

    def test_no_reaction_no_external_conserves(self):
        """Without reactions or external terms, advance must conserve total moles."""
        cv = _make_simple_cv()
        mol_before = cv.total_mol()
        cv.advance(dt_h=0.01)
        mol_after = cv.total_mol()
        for sp in set(mol_before) | set(mol_after):
            assert mol_after.get(sp, 0.0) == pytest.approx(
                mol_before.get(sp, 0.0), abs=1e-12
            ), f"Species {sp} not conserved"

    def test_reaction_conserves_within_reaction(self):
        """A→B reaction: A+B total moles conserved (C balanced)."""
        cv = _make_cv_with_reaction()
        A0 = cv["liquid"].n_mol["A"]
        B0 = cv["liquid"].n_mol.get("B", 0.0)
        cv.advance(dt_h=0.5)
        A1 = cv["liquid"].n_mol["A"]
        B1 = cv["liquid"].n_mol["B"]
        assert (A1 + B1) == pytest.approx(A0 + B0, abs=1e-12)

    def test_external_source_accounts_correctly(self):
        """External terms should change total by exactly the applied amount."""
        cv = _make_simple_cv()
        mol_before = cv.total_mol()
        # Add 0.1 mol/h * 0.5 h = 0.05 mol of new species "X" to liquid
        cv.advance(dt_h=0.5, external_source_terms={"liquid": {"X": 0.1}})
        mol_after = cv.total_mol()
        assert mol_after.get("X", 0.0) == pytest.approx(0.05, abs=1e-12)
        # Existing species should be unchanged (no interfaces)
        for sp in mol_before:
            assert mol_after.get(sp, 0.0) == pytest.approx(
                mol_before.get(sp, 0.0), abs=1e-12
            )


# ═══════════════════════════════════════════════════════════════════════
#  snapshot with advance-related attributes
# ═══════════════════════════════════════════════════════════════════════

class TestSnapshotWithAdvance:

    def test_snapshot_preserves_reaction_system(self):
        cv = _make_cv_with_reaction()
        cv2 = cv.snapshot()
        assert cv2.reaction_system is cv.reaction_system

    def test_snapshot_phases_independent(self):
        cv = _make_cv_with_reaction()
        cv2 = cv.snapshot()
        cv2["liquid"].n_mol["A"] = 999.0
        assert cv["liquid"].n_mol["A"] != 999.0


# ═══════════════════════════════════════════════════════════════════════
#  repr includes reaction_system
# ═══════════════════════════════════════════════════════════════════════

class TestRepr:

    def test_repr_shows_reaction_system_true(self):
        cv = _make_cv_with_reaction()
        assert "reaction_system=True" in repr(cv)

    def test_repr_shows_reaction_system_false(self):
        cv = _make_simple_cv()
        assert "reaction_system=False" in repr(cv)


# ═══════════════════════════════════════════════════════════════════════
#  Single CV with both gas + liquid phases + KineticGasLiquidLink
#  (Phase 3 enabling — see CV_UPDATE.md)
# ═══════════════════════════════════════════════════════════════════════

class TestCVWithKineticGasLiquidLink:
    """End-to-end: a single ControlVolume holding both phases with
    KineticGasLiquidLink as an internal PhaseInterface."""

    def _make_cv(self):
        from PyOMES.core import ControlVolume, GasPhase, LiquidPhase
        from PyOMES.core.gas_liquid_link import KineticGasLiquidLink
        from PyOMES.reactions import ReactionSystem

        gas = GasPhase(
            n_mol={"O2": 0.5, "CO2": 0.001, "N2": 2.0},
            V_L=0.4, T_K=305.15,
        )
        liquid = LiquidPhase(
            n_mol={"O2": 0.0, "CO2": 0.0, "N2": 0.0},
            V_L=1.6, T_K=305.15,
        )
        link = KineticGasLiquidLink(
            gas_cv_key="cv", gas_phase_key="gas",
            liquid_cv_key="cv", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3), "CO2": _hp(3.4e-2), "N2": _hp(6.5e-4)},
            kLa={"O2": 150.0, "CO2": 135.0},
            equilibrium_species={"N2"},
        )
        # Declared equilibria for the speciation engine. state-unification
        # C4: chemistry is defined by EquilibriumReactions; the engine
        # builds lazily off cv.reaction_system.
        system = ReactionSystem(_make_speciation_reactions(),
                                label="single_cv_gl_chem")
        return ControlVolume(
            phases={"gas": gas, "liquid": liquid},
            internal_interfaces=[link],
            reaction_system=system,
            label="single_cv_gl_transfer",
        )

    def test_construction_is_valid(self):
        cv = self._make_cv()
        assert "gas" in cv.phases and "liquid" in cv.phases
        assert len(cv.internal_interfaces) == 1
        assert cv.reaction_system is not None

    def test_advance_runs_to_completion(self):
        from PyOMES.core.interfaces import AdvanceResult
        cv = self._make_cv()
        result = cv.advance(0.01, 0.0)
        assert isinstance(result, AdvanceResult)

    def test_o2_transfers_gas_to_liquid(self):
        cv = self._make_cv()
        n_O2_gas_before = cv["gas"].n_mol["O2"]
        n_O2_liq_before = cv["liquid"].n_mol["O2"]

        cv.advance(0.01, 0.0)

        n_O2_gas_after = cv["gas"].n_mol["O2"]
        n_O2_liq_after = cv["liquid"].n_mol["O2"]
        # Gas should lose O2, liquid should gain O2
        assert n_O2_gas_after < n_O2_gas_before
        assert n_O2_liq_after > n_O2_liq_before

    def test_total_O2_conserved(self):
        """Internal transfer must conserve total moles of each species."""
        cv = self._make_cv()
        n_O2_total_before = cv["gas"].n_mol["O2"] + cv["liquid"].n_mol["O2"]
        cv.advance(0.01, 0.0)
        n_O2_total_after = cv["gas"].n_mol["O2"] + cv["liquid"].n_mol["O2"]
        assert n_O2_total_after == pytest.approx(n_O2_total_before, abs=1e-12)

    def test_transfer_diagnostics_conserved(self):
        cv = self._make_cv()
        result = cv.advance(0.01, 0.0)
        assert result.transfer is not None
        assert result.transfer.is_conserved()

    def test_speciation_runs_with_co2(self):
        """state-unification C4: pH after advance lives on
        cv.phases["liquid"].pH, derived from n_mol["H+"] populated
        by the engine writeback."""
        cv = self._make_cv()
        # Seed some dissolved CO2(aq) so speciation has something to work with
        cv["liquid"].n_mol["CO2"] = 0.005
        cv.advance(0.01, 0.0)
        liq = cv.phases["liquid"]
        assert "H+" in liq.n_mol
        assert 2.0 < liq.pH < 12.0

    def test_kinetic_approaches_equilibrium(self):
        """Over many steps, kinetic O₂ should approach Henry equilibrium.

        Ported from the deleted ``test_gas_liquid_volume.py`` (Phase 7 C5)."""
        from PyOMES.core import ControlVolume, GasPhase, LiquidPhase
        from PyOMES.core.gas_liquid_link import KineticGasLiquidLink
        from PyOMES.core.phases import R_L_ATM_MOL_K

        V_gas, V_liq, T_K = 0.4, 1.6, 305.15
        gas = GasPhase(
            n_mol={"O2": 0.5, "CO2": 0.01, "N2": 2.0},
            V_L=V_gas, T_K=T_K,
        )
        liquid = LiquidPhase(
            n_mol={"O2": 0.0, "CO2": 0.005, "N2": 0.0},
            V_L=V_liq, T_K=T_K,
        )
        link = KineticGasLiquidLink(
            gas_cv_key="cv", gas_phase_key="gas",
            liquid_cv_key="cv", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3), "CO2": _hp(3.4e-2), "N2": _hp(6.5e-4)},
            kLa={"O2": 500.0, "CO2": 450.0},
            equilibrium_species={"N2"},
        )
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liquid},
            internal_interfaces=[link],
        )

        for _ in range(500):
            cv.advance(dt_h=0.01)

        kH_O2 = 1.3e-3
        n_gas = cv.phases["gas"].n_mol["O2"]
        n_liq = cv.phases["liquid"].n_mol.get("O2", 0.0)
        p_eq = n_gas * R_L_ATM_MOL_K * T_K / V_gas
        C_eq = n_liq / V_liq
        assert C_eq == pytest.approx(kH_O2 * p_eq, rel=0.02)


# ═══════════════════════════════════════════════════════════════════════
#  CV-level boundaries (Phase 6, Checkpoint 1)
# ═══════════════════════════════════════════════════════════════════════

class TestCVBoundaries:
    """``boundaries`` is a ControlVolume-level concept (Phase 6)."""

    def test_liquid_feed_boundary_applies_without_gl_link(self):
        """A bare CV with a LiquidFeed boundary applies the feed flux
        through ``cv.advance(...)`` — no gas-liquid transfer link needed."""
        from PyOMES.core import ControlVolume, LiquidPhase, LiquidFeed

        liquid = LiquidPhase(n_mol={"AceticAcid": 0.0}, V_L=1.0, T_K=305.15)
        feed = LiquidFeed(
            Q_L_per_h=2.0,
            feed_conc_mol_L={"AceticAcid": 0.5},
        )
        cv = ControlVolume(
            phases={"liquid": liquid},
            boundaries=[feed],
            label="cv_with_liquid_feed",
        )

        cv.advance(dt_h=0.5)
        # 2 L/h × 0.5 mol/L × 0.5 h = 0.5 mol added
        assert cv["liquid"].n_mol["AceticAcid"] == pytest.approx(0.5, rel=1e-9)

    def test_no_boundaries_no_change(self):
        """Regression guard: a CV with no boundaries advances exactly as
        before — no extra fluxes from the new pass."""
        cv = _make_simple_cv()
        before = dict(cv.total_mol())
        cv.advance(dt_h=0.01)
        after = dict(cv.total_mol())
        for sp, n in before.items():
            assert after.get(sp, 0.0) == pytest.approx(n, abs=1e-12)

    def test_boundaries_list_is_mutable(self):
        """The boundaries list can be appended to after construction."""
        from PyOMES.core import LiquidFeed

        cv = _make_simple_cv()
        assert cv.boundaries == []
        feed = LiquidFeed(Q_L_per_h=0.0, feed_conc_mol_L={})
        cv.boundaries.append(feed)
        assert len(cv.boundaries) == 1


# ═══════════════════════════════════════════════════════════════════════
#  cv.advance(solver=...) dispatch (Phase 6, Checkpoint 3)
# ═══════════════════════════════════════════════════════════════════════

class TestCVAdvanceSolverDispatch:
    """`cv.advance(solver=...)` short-circuits the sequential body."""

    def _make_gl_cv(self, label="cv_with_link"):
        from PyOMES.core import (
            ControlVolume, GasPhase, LiquidPhase, KineticGasLiquidLink,
        )

        gas = GasPhase(
            n_mol={"O2": 0.5, "CO2": 0.001, "N2": 2.0},
            V_L=0.4, T_K=305.15,
        )
        liquid = LiquidPhase(
            n_mol={"O2": 0.0, "CO2": 0.0, "N2": 0.0},
            V_L=1.6, T_K=305.15,
        )
        link = KineticGasLiquidLink(
            gas_cv_key="cv", gas_phase_key="gas",
            liquid_cv_key="cv", liquid_phase_key="liquid",
            partition_models={"O2": _hp(1.3e-3), "CO2": _hp(3.4e-2), "N2": _hp(6.5e-4)},
            kLa={"O2": 150.0, "CO2": 135.0},
            equilibrium_species={"N2"},
        )
        # state-unification C4d: speciation via cv.reaction_system; no
        # SpeciationPropertySolver wrapper. CV has no equilibria
        # declared, so the engine lazy-builds and remains a no-op.
        return ControlVolume(
            phases={"gas": gas, "liquid": liquid},
            internal_interfaces=[link],
            label=label,
        )

    def test_bare_cv_with_euler_snapshot_solver_advances(self):
        """A bare CV under the snapshot solver advances without GLV."""
        from PyOMES.core import SimultaneousEulerSolver
        cv = self._make_gl_cv()
        n_O2_gas_before = cv["gas"].n_mol["O2"]
        n_O2_liq_before = cv["liquid"].n_mol["O2"]

        result = cv.advance(0.01, 0.0, solver=SimultaneousEulerSolver())

        # O2 should transfer from gas to liquid, just as with GLV
        assert cv["gas"].n_mol["O2"] < n_O2_gas_before
        assert cv["liquid"].n_mol["O2"] > n_O2_liq_before
        # Result populated
        assert result is not None

    def test_snapshot_solver_fails_fast_on_missing_liquid_phase(self):
        """A CV without a 'liquid' phase under SimultaneousEulerSolver
        must raise ValueError with a useful message (the C2 fail-fast) --
        speciation/reactions are liquid-scoped, so 'liquid' is the one
        phase this solver still requires."""
        from PyOMES.core import (
            ControlVolume, GasPhase, SimultaneousEulerSolver,
        )
        gas = GasPhase(n_mol={"X": 1.0}, V_L=1.0, T_K=300.0)
        cv = ControlVolume(phases={"gas": gas}, label="gas_only")
        with pytest.raises(ValueError, match="liquid"):
            cv.advance(dt_h=0.01, solver=SimultaneousEulerSolver())

    def test_snapshot_solver_generalizes_to_liquid_only_cv(self):
        """STEP_SOLVER_INTERFACE_REFINEMENT.md item 6: a liquid-only CV
        (feeds + kinetic reactions, no gas phase) now runs under
        SimultaneousEulerSolver -- previously this raised ValueError
        because the solver required exactly {'gas', 'liquid'}."""
        from PyOMES.core import (
            ControlVolume, LiquidPhase, SimultaneousEulerSolver,
        )
        liq = LiquidPhase(n_mol={"X": 1.0}, V_L=1.0, T_K=300.0)
        cv = ControlVolume(phases={"liquid": liq}, label="liquid_only")
        result = cv.advance(dt_h=0.01, solver=SimultaneousEulerSolver())
        assert result is not None
        assert "gas" not in cv.phases

    def test_adaptive_solver_fails_fast_on_missing_liquid_phase(self):
        """Checkpoint 7b: SimultaneousAdaptiveSolver's gas/liquid
        generalization -- a CV without a 'liquid' phase still fails
        fast with a useful message."""
        from PyOMES.core import (
            ControlVolume, GasPhase, SimultaneousAdaptiveSolver,
        )
        gas = GasPhase(n_mol={"X": 1.0}, V_L=1.0, T_K=300.0)
        cv = ControlVolume(phases={"gas": gas}, label="gas_only")
        with pytest.raises(ValueError, match="liquid"):
            cv.advance(dt_h=0.01, solver=SimultaneousAdaptiveSolver())

    @pytest.mark.parametrize("method", ["DOP853", "BDF", "Radau"])
    def test_adaptive_solver_generalizes_to_liquid_only_cv(self, method):
        """Checkpoint 7b: a liquid-only CV (no gas phase) now runs under
        SimultaneousAdaptiveSolver for both explicit (DOP853) and
        implicit (BDF/Radau) methods -- previously this raised
        ValueError for every method."""
        from PyOMES.core import (
            ControlVolume, LiquidPhase, SimultaneousAdaptiveSolver,
        )
        liq = LiquidPhase(n_mol={"X": 1.0}, V_L=1.0, T_K=300.0)
        cv = ControlVolume(phases={"liquid": liq}, label="liquid_only")
        result = cv.advance(dt_h=0.01, solver=SimultaneousAdaptiveSolver(method=method))
        assert result is not None
        assert "gas" not in cv.phases
        assert cv.phases["liquid"].n_mol["X"] == pytest.approx(1.0)

    def test_default_is_sugar_for_sequential_advance_solver(self):
        """cv.advance(dt_h, t_h) and cv.advance(dt_h, t_h,
        solver=SequentialAdvanceSolver()) must be byte-identical calls
        (STEP_SOLVER_INTERFACE_REFINEMENT.md item 3) -- solver=None
        dispatches to a real SequentialAdvanceSolver() instance rather
        than running a separate inline code path."""
        from PyOMES.core import SequentialAdvanceSolver

        cv_default = self._make_gl_cv(label="cv_default")
        cv_explicit = self._make_gl_cv(label="cv_explicit")

        result_default = cv_default.advance(0.01, 0.0)
        result_explicit = cv_explicit.advance(
            0.01, 0.0, solver=SequentialAdvanceSolver()
        )

        for phase_key in ("gas", "liquid"):
            assert cv_default[phase_key].n_mol == pytest.approx(
                cv_explicit[phase_key].n_mol
            )
        assert result_default.reaction_sources == result_explicit.reaction_sources


class TestEquilibrateToPH:
    """Regression tests for ControlVolume.equilibrate_to_pH.

    Uses an ionic medium (NH4+/Cl-/K+/H2PO4-) so no dissolution reactions
    are needed in the database. Verifies that "NaOH" as a strong-corrector
    alias raises pH correctly without requiring an eq_NaOH reaction.
    """

    def _make_cv(self):
        from PyOMES.chemistry.common_species import (
            NH4_plus, Cl_minus, K_plus, H2PO4_minus, HPO4_2minus, PO4_3minus,
            H3PO4, Na_plus,
        )
        from PyOMES.chemistry.databases.bioprocess_basic import (
            BIOPROCESS_BASIC, NH4Cl, KH2PO4,
        )
        from PyOMES.reactions import ReactionSystem
        from PyOMES.core import ControlVolume, LiquidPhase

        V_L = 250e-6  # 250 µL
        # Typical MTP-well loadings (5 g/L NH4Cl, 30 g/L KH2PO4)
        n_NH4 = (5.0 / float(NH4Cl.MW)) * V_L
        n_KH2 = (30.0 / float(KH2PO4.MW)) * V_L

        liquid = LiquidPhase(
            n_mol={
                NH4_plus.id:    n_NH4,
                Cl_minus.id:    n_NH4,
                K_plus.id:      n_KH2,
                H2PO4_minus.id: n_KH2,
                HPO4_2minus.id: 0.0,
                H3PO4.id:       0.0,
                PO4_3minus.id:  0.0,
                Na_plus.id:     0.0,
                "NH3": 0.0, "HCO3-": 0.0, "CO3--": 0.0,
                "OH-": 0.0, "H+": 1e-7 * V_L,
            },
            V_L=V_L, T_K=305.15,
        )

        relevant = {"eq_water", "eq_CO2", "eq_NH4",
                    "eq_phosphate_1", "eq_phosphate_2", "eq_phosphate_3"}
        rxns = ReactionSystem(
            [r for r in BIOPROCESS_BASIC.reactions if r.label in relevant]
        )

        return ControlVolume(
            phases={"liquid": liquid},
            reaction_system=rxns,
            label="test_eq_pH",
        )

    def test_ionic_medium_gives_acidic_baseline(self):
        """H2PO4- loaded directly → pH around 3.5–5.5 (no dissolution artefacts)."""
        cv = self._make_cv()
        cv.advance(dt_h=0.0, t_h=0.0)
        pH = float(cv["liquid"].pH)
        assert 3.0 < pH < 6.0, f"Expected baseline pH 3–6, got {pH:.3f}"

    def test_naoh_raises_ph_to_setpoint(self):
        """equilibrate_to_pH('NaOH', 6.0) must reach pH 6 ± 0.01 without error."""
        cv = self._make_cv()
        n_added = cv.equilibrate_to_pH("NaOH", 6.0)
        pH_final = float(cv["liquid"].pH)
        assert n_added > 0, "Expected positive moles of NaOH added"
        assert abs(pH_final - 6.0) < 0.01, f"Expected pH ≈ 6.0, got {pH_final:.4f}"

    def test_naoh_adds_na_plus_to_n_mol(self):
        """Strong-corrector path must increase Na+ in n_mol, not NaOH."""
        from PyOMES.chemistry.common_species import Na_plus
        cv = self._make_cv()
        na_before = cv["liquid"].n_mol.get(Na_plus.id, 0.0)
        cv.equilibrate_to_pH("NaOH", 6.0)
        na_after = cv["liquid"].n_mol.get(Na_plus.id, 0.0)
        assert na_after > na_before, "Na+ should increase after NaOH correction"
        assert cv["liquid"].n_mol.get("NaOH", 0.0) == 0.0, "NaOH should not appear in n_mol"

    def test_unknown_corrector_raises(self):
        """A corrector not in equilibrium stoichiometry must raise ValueError."""
        cv = self._make_cv()
        with pytest.raises(ValueError, match="does not appear in any EquilibriumReaction"):
            cv.equilibrate_to_pH("NaOH_unknown_salt", 6.0)

    def test_advance_after_equilibrate_no_conservation_warning(self):
        """10 advance steps after equilibrate_to_pH must fire no ConservationWarning.

        Regression guard for the charge-residual bug where strong ions
        (Cl-, Na+, K+) were absent from the species registry and caused
        spurious charge-imbalance warnings.
        """
        import warnings
        from PyOMES.monitoring.conservation import ConservationWarning

        cv = self._make_cv()
        cv.equilibrate_to_pH("NaOH", 6.0)

        with warnings.catch_warnings():
            warnings.simplefilter("error", ConservationWarning)
            for i in range(10):
                cv.advance(dt_h=0.01, t_h=i * 0.01)
