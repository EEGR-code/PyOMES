# -*- coding: utf-8 -*-
"""Tests for Phase A — CV_COMPUTE_INTERFACE.

Covers:
- compute_rhs: rates match advance() deltas on the same state
- snapshot_state / restore_state: round-trip fidelity
- algebraic_species: correct frozenset for declared and legacy engines
- compute_differential_rhs: same kinetic rates as compute_rhs at consistent state
- compute_algebraic_residual: NotImplementedError stub
- compute_jacobian: returns None
"""

import pytest


# ═══════════════════════════════════════════════════════════════════════
#  Helpers shared across test groups
# ═══════════════════════════════════════════════════════════════════════

def _make_simple_reaction():
    """A → B in liquid at rate 0.5 * [A] * V_L mol/h."""
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
    """CV: liquid only, A→B reaction, no speciation, no interfaces."""
    from PyOMES.core import ControlVolume, LiquidPhase
    liquid = LiquidPhase(n_mol={"A": 0.1, "B": 0.0}, V_L=1.0, T_K=298.15)
    return ControlVolume(
        phases={"liquid": liquid},
        reaction_system=_make_simple_reaction(),
        label="rxn_only",
    )


def _make_cv_no_reaction():
    """CV with no reaction_system."""
    from PyOMES.core import ControlVolume, LiquidPhase
    liquid = LiquidPhase(n_mol={"X": 0.5}, V_L=1.0, T_K=298.15)
    return ControlVolume(phases={"liquid": liquid}, label="no_rxn")


def _make_speciation_reactions():
    """Minimal equilibrium reaction set: water + CO2 + acetate."""
    from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
    from PyOMES.chemistry import Species
    from PyOMES.chemistry.common_species import (
        H_plus, OH_minus, H2O, CO2, HCO3_minus,
    )
    HAc = Species(id="AceticAcid", atoms={"C": 2, "H": 4, "O": 2}, charge=0)
    Ac = Species(id="AceticAcid-", atoms={"C": 2, "H": 3, "O": 2}, charge=-1)
    return [
        EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=H2O,      phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=H_plus,   phase="liquid", coefficient=+1.0),
                StoichiometryEntry(species=OH_minus, phase="liquid", coefficient=+1.0),
            ],
            log_K=-14.0, label="water",
        ),
        EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=CO2,        phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=H2O,        phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=HCO3_minus, phase="liquid", coefficient=+1.0),
                StoichiometryEntry(species=H_plus,     phase="liquid", coefficient=+1.0),
            ],
            log_K=-6.35, label="CO2",
        ),
        EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=HAc,  phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=Ac,   phase="liquid", coefficient=+1.0),
                StoichiometryEntry(species=H_plus, phase="liquid", coefficient=+1.0),
            ],
            log_K=-4.76, label="acetate",
        ),
    ]


def _make_cv_with_speciation_and_reaction():
    """CV: liquid, acetate equilibria + A→B kinetic reaction."""
    from PyOMES.core import ControlVolume, LiquidPhase
    from PyOMES.reactions import ReactionSystem
    liquid = LiquidPhase(
        n_mol={"A": 0.1, "B": 0.0, "AceticAcid": 0.01, "H+": 1e-7},
        V_L=1.0, T_K=298.15,
    )
    rxns = _make_speciation_reactions() + [_make_simple_reaction()]
    system = ReactionSystem(rxns, label="spec_and_rxn")
    return ControlVolume(
        phases={"liquid": liquid},
        reaction_system=system,
        label="spec_rxn",
    )


# ═══════════════════════════════════════════════════════════════════════
#  compute_rhs
# ═══════════════════════════════════════════════════════════════════════

class TestComputeRhs:

    def test_no_reaction_system_returns_empty(self):
        cv = _make_cv_no_reaction()
        assert cv.compute_rhs(t_h=0.0) == {}

    def test_rates_match_advance_deltas(self):
        """compute_rhs rates × dt_h should equal the n_mol change from advance()."""
        dt_h = 0.01
        cv = _make_cv_with_reaction()
        # Record pre-step state
        n_before = dict(cv.phases["liquid"].n_mol)

        # Get rates from compute_rhs (does NOT mutate state)
        rhs = cv.compute_rhs(t_h=0.0)

        # State must be unchanged after compute_rhs
        assert cv.phases["liquid"].n_mol == n_before

        # Now advance the CV
        cv.advance(dt_h=dt_h, t_h=0.0)
        n_after = dict(cv.phases["liquid"].n_mol)

        # For each species with a rate, delta should match rate * dt_h
        for sp, rate in rhs.get("liquid", {}).items():
            expected = max(0.0, n_before.get(sp, 0.0) + rate * dt_h)
            assert n_after[sp] == pytest.approx(expected, rel=1e-9, abs=1e-12), \
                f"species {sp!r}: expected {expected}, got {n_after[sp]}"

    def test_returns_dict_of_dicts(self):
        cv = _make_cv_with_reaction()
        rhs = cv.compute_rhs(t_h=0.0)
        assert isinstance(rhs, dict)
        for pk, sp_rates in rhs.items():
            assert isinstance(pk, str)
            assert isinstance(sp_rates, dict)

    def test_rates_are_mol_per_h(self):
        """A→B at rate 0.5 * [A] * V_L: with [A]=0.1 mol/L, V_L=1, rate=0.05 mol/h."""
        cv = _make_cv_with_reaction()
        rhs = cv.compute_rhs(t_h=0.0)
        liq_rates = rhs.get("liquid", {})
        # rate_A = -0.5 * (0.1 / 1.0) * 1.0 = -0.05 mol/h
        assert liq_rates["A"] == pytest.approx(-0.05, rel=1e-9)
        assert liq_rates["B"] == pytest.approx(+0.05, rel=1e-9)

    def test_does_not_mutate_state(self):
        """compute_rhs must not mutate differential state."""
        cv = _make_cv_with_reaction()
        snap = {pk: dict(p.n_mol) for pk, p in cv.phases.items()}
        cv.compute_rhs(t_h=0.0)
        for pk, n_mol in snap.items():
            assert cv.phases[pk].n_mol == n_mol

    def test_with_speciation_and_kinetics(self):
        """compute_rhs on a CV with both speciation and kinetics returns kinetic rates."""
        cv = _make_cv_with_speciation_and_reaction()
        rhs = cv.compute_rhs(t_h=0.0)
        liq = rhs.get("liquid", {})
        # A→B kinetics must still be in the rates
        assert "A" in liq
        assert "B" in liq
        assert liq["A"] == pytest.approx(-liq["B"], rel=1e-9)


# ═══════════════════════════════════════════════════════════════════════
#  snapshot_state / restore_state
# ═══════════════════════════════════════════════════════════════════════

class TestSnapshotRestore:

    def test_snapshot_returns_deep_copy(self):
        cv = _make_cv_with_reaction()
        snap = cv.snapshot_state()
        # Mutate original; snapshot must be unaffected
        cv.phases["liquid"].n_mol["A"] = 999.0
        assert snap["liquid"]["A"] == pytest.approx(0.1)

    def test_restore_overwrites_state(self):
        cv = _make_cv_with_reaction()
        snap = cv.snapshot_state()
        cv.phases["liquid"].n_mol["A"] = 999.0
        cv.restore_state(snap)
        assert cv.phases["liquid"].n_mol["A"] == pytest.approx(0.1)

    def test_round_trip_fidelity(self):
        """snapshot → mutate → restore → state is bit-identical to original."""
        cv = _make_cv_with_reaction()
        snap = cv.snapshot_state()
        original = {pk: dict(p.n_mol) for pk, p in cv.phases.items()}
        # Mutate
        cv.advance(dt_h=0.1, t_h=0.0)
        # Restore
        cv.restore_state(snap)
        for pk, n_mol in original.items():
            for sp, val in n_mol.items():
                assert cv.phases[pk].n_mol[sp] == val

    def test_snapshot_covers_all_phases(self):
        from PyOMES.core import ControlVolume, GasPhase, LiquidPhase
        gas = GasPhase(n_mol={"O2": 1.0, "N2": 2.0}, V_L=0.5, T_K=298.15)
        liq = LiquidPhase(n_mol={"X": 0.3}, V_L=1.0, T_K=298.15)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq})
        snap = cv.snapshot_state()
        assert set(snap.keys()) == {"gas", "liquid"}
        assert snap["gas"]["O2"] == pytest.approx(1.0)
        assert snap["liquid"]["X"] == pytest.approx(0.3)

    def test_restore_clears_added_species(self):
        """Species added after snapshot must be cleared on restore."""
        cv = _make_cv_with_reaction()
        snap = cv.snapshot_state()
        cv.phases["liquid"].n_mol["NEW_SP"] = 5.0
        cv.restore_state(snap)
        assert "NEW_SP" not in cv.phases["liquid"].n_mol

    def test_snapshot_then_compute_rhs_then_restore(self):
        """Standard trial-state pattern: compute without committing."""
        cv = _make_cv_with_reaction()
        snap = cv.snapshot_state()
        rhs = cv.compute_rhs(t_h=0.0)
        cv.restore_state(snap)
        # State must be back to original (compute_rhs should not have mutated it)
        assert cv.phases["liquid"].n_mol["A"] == pytest.approx(0.1)


# ═══════════════════════════════════════════════════════════════════════
#  algebraic_species
# ═══════════════════════════════════════════════════════════════════════

class TestAlgebraicSpecies:

    def test_legacy_engine_returns_empty_frozenset(self):
        """Engine constructed without from_reactions has no _equilibrium_set."""
        from PyOMES.chemical_equilibrium import BisectionChemicalEquilibriumEngine
        engine = BisectionChemicalEquilibriumEngine()
        result = engine.algebraic_species()
        assert result == frozenset()
        assert isinstance(result, frozenset)

    def test_from_reactions_includes_water_ions(self):
        """Any engine built from reactions must include H+ and OH-."""
        from PyOMES.chemical_equilibrium import BisectionChemicalEquilibriumEngine
        engine = BisectionChemicalEquilibriumEngine.from_reactions(_make_speciation_reactions())
        alg = engine.algebraic_species()
        assert "H+" in alg
        assert "OH-" in alg

    def test_from_reactions_includes_ladder_species(self):
        """species_refs from each EquilibriumReaction appear in algebraic_species."""
        from PyOMES.chemical_equilibrium import BisectionChemicalEquilibriumEngine
        from PyOMES.chemistry.common_species import CO2, HCO3_minus
        engine = BisectionChemicalEquilibriumEngine.from_reactions(_make_speciation_reactions())
        alg = engine.algebraic_species()
        # CO2 ladder: CO2, HCO3-
        assert CO2.id in alg
        assert HCO3_minus.id in alg
        # Acetate ladder: AceticAcid, AceticAcid-
        assert "AceticAcid" in alg
        assert "AceticAcid-" in alg

    def test_bsm2_equilibrium_set_includes_carbonate_and_nh(self):
        """EquilibriumSet.bsm2_default() species_refs are all in algebraic_species."""
        from PyOMES.chemical_equilibrium import BisectionChemicalEquilibriumEngine
        from PyOMES.chemistry.equilibria import EquilibriumSet
        eq_set = EquilibriumSet.bsm2_default()
        engine = BisectionChemicalEquilibriumEngine()
        engine._equilibrium_set = eq_set
        alg = engine.algebraic_species()

        assert "H+" in alg
        assert "OH-" in alg

        for eq_def in eq_set:
            for sp in eq_def.species_refs:
                assert sp.id in alg, f"expected {sp.id!r} in algebraic_species"

    def test_bsm2_vfa_entries_not_included(self):
        """BSM2 VFA entries have no species_refs and must NOT be in the frozenset."""
        from PyOMES.chemical_equilibrium import BisectionChemicalEquilibriumEngine
        from PyOMES.chemistry.equilibria import EquilibriumSet
        eq_set = EquilibriumSet.bsm2_default()
        engine = BisectionChemicalEquilibriumEngine()
        engine._equilibrium_set = eq_set
        alg = engine.algebraic_species()
        for name in ("S_ac", "S_pro", "S_bu", "S_va"):
            assert name not in alg

    def test_returns_frozenset_type(self):
        from PyOMES.chemical_equilibrium import BisectionChemicalEquilibriumEngine
        engine = BisectionChemicalEquilibriumEngine.from_reactions(_make_speciation_reactions())
        assert isinstance(engine.algebraic_species(), frozenset)


# ═══════════════════════════════════════════════════════════════════════
#  compute_differential_rhs
# ═══════════════════════════════════════════════════════════════════════

class TestComputeDifferentialRhs:

    def test_no_reaction_system_returns_empty(self):
        cv = _make_cv_no_reaction()
        assert cv.compute_differential_rhs(0.0, {}, {}) == {}

    def test_consistent_state_matches_compute_rhs(self):
        """compute_differential_rhs at (y=current, z={}) gives same rates as compute_rhs."""
        cv = _make_cv_with_reaction()
        # y_state = current differential state (no algebraic species)
        y_state = {pk: dict(p.n_mol) for pk, p in cv.phases.items()}
        z_state = {}

        diff_rates = cv.compute_differential_rhs(0.0, y_state, z_state)
        rhs_rates = cv.compute_rhs(0.0)

        for pk in rhs_rates:
            for sp, rate in rhs_rates[pk].items():
                assert diff_rates.get(pk, {}).get(sp, 0.0) == pytest.approx(
                    rate, rel=1e-9, abs=1e-14
                ), f"mismatch for {pk}/{sp}"

    def test_does_not_touch_phase_n_mol(self):
        """compute_differential_rhs must never write to phase.n_mol."""
        cv = _make_cv_with_reaction()
        n_before = dict(cv.phases["liquid"].n_mol)
        y_state = {"liquid": dict(cv.phases["liquid"].n_mol)}
        cv.compute_differential_rhs(0.0, y_state, {})
        assert cv.phases["liquid"].n_mol == n_before

    def test_uses_y_state_concentrations_not_phase(self):
        """Rates are evaluated at y_state, not at the current phase.n_mol."""
        cv = _make_cv_with_reaction()
        # Set y_state to [A]=0.2 (double), phase.n_mol stays at 0.1
        y_state = {"liquid": {"A": 0.2, "B": 0.0}}
        rates = cv.compute_differential_rhs(0.0, y_state, {})
        # rate = 0.5 * (0.2 / 1.0) * 1.0 = 0.1 mol/h
        assert rates.get("liquid", {}).get("A", None) == pytest.approx(-0.1, rel=1e-9)

    def test_z_state_overrides_y_state_for_same_species(self):
        """z_state values override y_state for species present in both."""
        cv = _make_cv_with_reaction()
        # y_state has A=0.1; z_state overrides A=0.4
        y_state = {"liquid": {"A": 0.1, "B": 0.0}}
        z_state = {"liquid": {"A": 0.4}}
        rates = cv.compute_differential_rhs(0.0, y_state, z_state)
        # rate = 0.5 * (0.4 / 1.0) * 1.0 = 0.2 mol/h
        assert rates.get("liquid", {}).get("A", None) == pytest.approx(-0.2, rel=1e-9)

    def test_returns_dict_of_dicts(self):
        cv = _make_cv_with_reaction()
        y_state = {"liquid": dict(cv.phases["liquid"].n_mol)}
        result = cv.compute_differential_rhs(0.0, y_state, {})
        assert isinstance(result, dict)
        for pk, sp_rates in result.items():
            assert isinstance(sp_rates, dict)


# ═══════════════════════════════════════════════════════════════════════
#  compute_algebraic_residual (stub) and compute_jacobian (stub)
# ═══════════════════════════════════════════════════════════════════════

class TestStubs:

    def test_compute_algebraic_residual_raises(self):
        cv = _make_cv_with_reaction()
        with pytest.raises(NotImplementedError):
            cv.compute_algebraic_residual(0.0, {}, {})

    def test_compute_jacobian_returns_none(self):
        cv = _make_cv_with_reaction()
        assert cv.compute_jacobian(0.0) is None

    def test_compute_jacobian_no_reaction(self):
        cv = _make_cv_no_reaction()
        assert cv.compute_jacobian(0.0) is None


# ═══════════════════════════════════════════════════════════════════════
#  advance() regression — behaviour unchanged after refactor
# ═══════════════════════════════════════════════════════════════════════

class TestAdvanceRegression:

    def test_simple_reaction_euler_step(self):
        """A→B with dt_h=0.1: A decreases, B increases, sum conserved."""
        cv = _make_cv_with_reaction()
        n_A_before = cv.phases["liquid"].n_mol["A"]
        cv.advance(dt_h=0.1, t_h=0.0)
        n_A_after = cv.phases["liquid"].n_mol["A"]
        n_B_after = cv.phases["liquid"].n_mol["B"]
        assert n_A_after < n_A_before
        assert n_B_after > 0.0
        # Carbon conserved
        assert n_A_after + n_B_after == pytest.approx(n_A_before, rel=1e-9)

    def test_advance_result_has_reaction_sources(self):
        cv = _make_cv_with_reaction()
        result = cv.advance(dt_h=0.01, t_h=0.0)
        assert result.reaction_sources is not None
        assert "liquid" in result.reaction_sources
        liq = result.reaction_sources["liquid"]
        assert "A" in liq
        assert "B" in liq

    def test_no_reaction_system_advance_is_noop(self):
        cv = _make_cv_no_reaction()
        snap = cv.snapshot_state()
        cv.advance(dt_h=0.1, t_h=0.0)
        assert cv.phases["liquid"].n_mol == snap["liquid"]
