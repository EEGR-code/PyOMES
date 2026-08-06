# -*- coding: utf-8 -*-
"""Tests for the reactions framework (Stage 0).

NOTE: Tests are written to be compatible with both pytest and the
unittest-based runner in run_tests.py (which does not support
pytest.warns or helper methods on test classes).
"""
import math
import warnings
import numpy as np
import pytest

# ── Helper factories (module-level) ──────────────────────────────────

def _sp(sp_id, atoms=None, charge=0, MW=0.0):
    """Build a Species — convenience for tests that don't care about
    MW or persisted identity."""
    from PyOMES.chemistry import Species
    return Species(id=sp_id, atoms=dict(atoms or {}), charge=charge, MW=float(MW))


def _entry(sp_id, coeff, atoms=None, *, phase="liquid", charge=0):
    """Build a StoichiometryEntry with a fresh Species."""
    from PyOMES.reactions import StoichiometryEntry
    return StoichiometryEntry(species=_sp(sp_id, atoms, charge=charge),
                              phase=phase, coefficient=coeff)


def _make_simple_reaction():
    from PyOMES.reactions import KineticReaction
    return KineticReaction(
        stoichiometry=[
            _entry("A", -1.0, {"C": 1}),
            _entry("B", +1.0, {"C": 1}),
        ],
        rate_fn=lambda env: 0.1 * env.S("A") * env.V_L,
        balance_elements=("C",), label="A_to_B",
    )

def _make_reaction_system():
    from PyOMES.reactions import KineticReaction, ReactionSystem
    rxn1 = KineticReaction(
        [_entry("A", -1.0, {"C": 1}),
         _entry("C", +1.0, {"C": 1})],
        rate_fn=lambda env: 2.0 * env.S("A") * env.V_L, balance_elements=("C",))
    rxn2 = KineticReaction(
        [_entry("B", -1.0, {"C": 2}),
         _entry("D", +1.0, {"C": 2})],
        rate_fn=lambda env: 0.5 * env.S("B") * env.V_L, balance_elements=("C",))
    return ReactionSystem([rxn1, rxn2])

def _make_balanced_blackbox():
    from PyOMES.reactions import BlackBoxReactionModel, FluxEntry
    class MockModel:
        def solve(self, inputs):
            rate = inputs.get("A", 0.0) * 0.5
            return {"consume_A": -rate, "produce_B": rate}
    return BlackBoxReactionModel(
        external_model=MockModel(),
        flux_mapping={
            "consume_A": FluxEntry("A", "liquid", {"C": 1}, sign=+1),
            "produce_B": FluxEntry("B", "liquid", {"C": 1}, sign=+1),
        }, balance_elements=("C",))

# ── Stoichiometry validation ─────────────────────────────────────────

class TestStoichiometryValidation:
    def test_balanced_reaction_passes(self):
        from PyOMES.reactions import validate_balance
        entries = [_entry("A",-1.0,{"C":1,"H":4}),
                   _entry("B",+1.0,{"C":1,"H":4})]
        validate_balance(entries, elements=("C","H"))

    def test_imbalanced_carbon_raises(self):
        from PyOMES.reactions import StoichiometryError, validate_balance
        entries = [_entry("A",-1.0,{"C":2}),
                   _entry("B",+1.0,{"C":1})]
        with pytest.raises(StoichiometryError):
            validate_balance(entries, elements=("C",))

    def test_error_carries_element_and_residual(self):
        from PyOMES.reactions import StoichiometryError, validate_balance
        entries = [_entry("A",-1.0,{"C":3}),
                   _entry("B",+1.0,{"C":1})]
        try:
            validate_balance(entries, elements=("C",))
            assert False, "Should have raised"
        except StoichiometryError as e:
            assert e.element == "C"
            assert abs(e.residual - (-2.0)) < 1e-10

    def test_tolerance_allows_small_residuals(self):
        from PyOMES.reactions import validate_balance
        entries = [_entry("A",-1.0,{"C":1.0}),
                   _entry("B",+1.0,{"C":1.0+1e-15})]
        validate_balance(entries, elements=("C",), atol=1e-10)

    def test_multi_element_checks_all(self):
        from PyOMES.reactions import StoichiometryError, validate_balance
        entries = [_entry("A",-1.0,{"C":1,"H":6}),
                   _entry("B",+1.0,{"C":1,"H":4})]
        validate_balance(entries, elements=("C",))
        with pytest.raises(StoichiometryError):
            validate_balance(entries, elements=("C","H"))

    def test_missing_element_treated_as_zero(self):
        from PyOMES.reactions import StoichiometryError, validate_balance
        entries = [_entry("O2",-1.0,{"O":2}),
                   _entry("H2O",+1.0,{"H":2,"O":1})]
        with pytest.raises(StoichiometryError):
            validate_balance(entries, elements=("O",))

# ── KineticReaction class ────────────────────────────────────────────

class TestReaction:
    def test_construction_validates(self):
        rxn = _make_simple_reaction()
        assert rxn.label == "A_to_B"

    def test_imbalanced_construction_raises(self):
        from PyOMES.reactions import KineticReaction, StoichiometryError
        with pytest.raises(StoichiometryError):
            KineticReaction([_entry("A",-1.0,{"C":2}),
                             _entry("B",+1.0,{"C":1})],
                            rate_fn=lambda env: 0.0, balance_elements=("C",))

    def test_compute_rates_correct_signs(self):
        from PyOMES.reactions import ReactionEnvironment
        rxn = _make_simple_reaction()
        env = ReactionEnvironment(concentrations={"A":1.0,"B":0.0}, V_L=10.0)
        rates = rxn.compute_rates(env)
        assert rates["liquid"]["A"] == pytest.approx(-1.0, abs=1e-12)
        assert rates["liquid"]["B"] == pytest.approx(+1.0, abs=1e-12)

    def test_zero_concentration_gives_zero_rate(self):
        from PyOMES.reactions import ReactionEnvironment
        rxn = _make_simple_reaction()
        env = ReactionEnvironment(concentrations={"A":0.0}, V_L=10.0)
        rates = rxn.compute_rates(env)
        assert abs(rates["liquid"]["A"]) < 1e-15

    def test_species_ids_property(self):
        assert _make_simple_reaction().species_ids == ["A","B"]

    def test_phases_property(self):
        assert _make_simple_reaction().phases == ["liquid"]

# ── ReactionBuilder ──────────────────────────────────────────────────

class TestReactionBuilder:
    def test_acetic_acid_cho_constructs(self):
        from PyOMES.reactions import ReactionBuilder
        rxn = ReactionBuilder.aerobic_growth(
            substrate_id="AceticAcid", substrate_atoms={"C":2,"H":4,"O":2}, MW_substrate=60.052,
            biomass_id="Yeast_CHO", biomass_atoms={"C":1,"H":1.61,"O":0.56}, MW_biomass=24.626,
            yield_gX_gS=0.36, rate_fn=lambda env: 1.0, balance="CHO")
        assert "AceticAcid" in rxn.species_ids

    def test_propionic_acid_cho_constructs(self):
        from PyOMES.reactions import ReactionBuilder
        ReactionBuilder.aerobic_growth(
            substrate_id="PropionicAcid", substrate_atoms={"C":3,"H":6,"O":2}, MW_substrate=74.079,
            biomass_id="Yeast_CHO", biomass_atoms={"C":1,"H":1.61,"O":0.56}, MW_biomass=24.626,
            yield_gX_gS=0.486*1.11, rate_fn=lambda env: 1.0, balance="CHO")

    def test_butyric_acid_cho_constructs(self):
        from PyOMES.reactions import ReactionBuilder
        ReactionBuilder.aerobic_growth(
            substrate_id="ButyricAcid", substrate_atoms={"C":4,"H":8,"O":2}, MW_substrate=88.106,
            biomass_id="Yeast_CHO", biomass_atoms={"C":1,"H":1.61,"O":0.56}, MW_biomass=24.626,
            yield_gX_gS=0.545*1.05, rate_fn=lambda env: 1.0, balance="CHO")

    def test_chno_mode_includes_nitrogen(self):
        from PyOMES.reactions import ReactionBuilder
        rxn = ReactionBuilder.aerobic_growth(
            substrate_id="AceticAcid", substrate_atoms={"C":2,"H":4,"O":2}, MW_substrate=60.052,
            biomass_id="Yeast", biomass_atoms={"C":1,"H":1.61,"O":0.56,"N":0.16}, MW_biomass=26.868,
            yield_gX_gS=0.36, rate_fn=lambda env: 1.0, balance="CHNO")
        assert "NH3" in rxn.species_ids

    def test_cho_parity_with_existing_stoichiometry(self):
        """Builder must match _stoich_aerobic_CHO from fermenter_unit.py."""
        from PyOMES.reactions import ReactionBuilder
        Cs,Hs,Os = 2.0,4.0,2.0; Cx,Hx,Ox = 1.0,1.61,0.56
        MW_S,MW_X = 60.052,24.626; Y_gXgS = 0.9*0.400
        Y_mol = Y_gXgS*MW_S/MW_X
        nCO2_exp = Cs - Y_mol*Cx
        nH2O_exp = (Hs - Y_mol*Hx)/2.0
        nO2_exp  = (Y_mol*Ox + 2.0*nCO2_exp + nH2O_exp - Os)/2.0
        rxn = ReactionBuilder.aerobic_growth(
            substrate_id="AceticAcid", substrate_atoms={"C":Cs,"H":Hs,"O":Os}, MW_substrate=MW_S,
            biomass_id="Yeast_CHO", biomass_atoms={"C":Cx,"H":Hx,"O":Ox}, MW_biomass=MW_X,
            yield_gX_gS=Y_gXgS, rate_fn=lambda env: 1.0, balance="CHO")
        coeffs = {e.species.id: e.coefficient for e in rxn.stoichiometry}
        assert coeffs["AceticAcid"] == pytest.approx(-1.0, abs=1e-12)
        assert coeffs["CO2"] == pytest.approx(nCO2_exp, abs=1e-10)
        assert coeffs["H2O"] == pytest.approx(nH2O_exp, abs=1e-10)
        assert coeffs["O2"]  == pytest.approx(-nO2_exp, abs=1e-10)
        assert coeffs["Yeast_CHO"] == pytest.approx(Y_mol, abs=1e-10)

    def test_chno_parity_nitrogen_demand(self):
        from PyOMES.reactions import ReactionBuilder
        Cs,Hs,Os,Ns = 2.0,4.0,2.0,0.0; Cx,Hx,Ox,Nx = 1.0,1.61,0.56,0.16
        MW_S,MW_X = 60.052,26.868; Y_gXgS = 0.36
        Y_mol = Y_gXgS*MW_S/MW_X
        nNH3_exp = Y_mol*Nx - Ns
        rxn = ReactionBuilder.aerobic_growth(
            substrate_id="AceticAcid", substrate_atoms={"C":Cs,"H":Hs,"O":Os,"N":Ns}, MW_substrate=MW_S,
            biomass_id="Yeast", biomass_atoms={"C":Cx,"H":Hx,"O":Ox,"N":Nx}, MW_biomass=MW_X,
            yield_gX_gS=Y_gXgS, rate_fn=lambda env: 1.0, balance="CHNO")
        coeffs = {e.species.id: e.coefficient for e in rxn.stoichiometry}
        assert coeffs["NH3"] == pytest.approx(-nNH3_exp, abs=1e-10)

    def test_zero_yield_all_substrate_to_co2(self):
        from PyOMES.reactions import ReactionBuilder
        rxn = ReactionBuilder.aerobic_growth(
            substrate_id="AceticAcid", substrate_atoms={"C":2,"H":4,"O":2}, MW_substrate=60.052,
            biomass_id="X", biomass_atoms={"C":1,"H":1.61,"O":0.56}, MW_biomass=24.626,
            yield_gX_gS=0.0, rate_fn=lambda env: 1.0, balance="CHO")
        coeffs = {e.species.id: e.coefficient for e in rxn.stoichiometry}
        assert coeffs["CO2"] == pytest.approx(2.0, abs=1e-10)
        assert abs(coeffs["X"]) < 1e-12

    def test_invalid_mw_raises(self):
        from PyOMES.reactions import ReactionBuilder
        with pytest.raises(ValueError):
            ReactionBuilder.aerobic_growth(
                substrate_id="A", substrate_atoms={"C":1}, MW_substrate=0.0,
                biomass_id="X", biomass_atoms={"C":1}, MW_biomass=1.0,
                yield_gX_gS=0.5, rate_fn=lambda env: 0.0)

    def test_from_coefficients(self):
        from PyOMES.reactions import ReactionBuilder
        rxn = ReactionBuilder.from_coefficients(
            entries=[_entry("A",-1.0,{"C":2}),
                     _entry("B",+0.5,{"C":2}),
                     _entry("CO2",+1.0,{"C":1})],
            rate_fn=lambda env: 1.0, balance_elements=("C",))
        assert len(rxn.stoichiometry) == 3

# ── ReactionBuilder.monod_aerobic_growth ─────────────────────────────

class TestMonodAerobicGrowth:
    """Tests for the Monod convenience wrapper."""

    @pytest.fixture
    def substrate(self):
        from PyOMES.chemistry import Species
        return Species(id="AceticAcid", atoms={"C": 2, "H": 4, "O": 2}, charge=0)

    @pytest.fixture
    def biomass(self):
        from PyOMES.chemistry import Species
        return Species(id="Yeast", atoms={"C": 1, "H": 1.61, "O": 0.56}, charge=0, MW=24.626)

    def test_returns_kinetic_reaction(self, substrate, biomass):
        from PyOMES.reactions import ReactionBuilder, KineticReaction
        rxn = ReactionBuilder.monod_aerobic_growth(
            substrate=substrate, biomass=biomass,
            mu_max_per_h=0.5, Ks_gL=5e-3, yield_gX_gS=0.36,
        )
        assert isinstance(rxn, KineticReaction)

    def test_stoichiometry_matches_aerobic_growth(self, substrate, biomass):
        """monod_aerobic_growth must produce identical stoichiometry to aerobic_growth."""
        from PyOMES.reactions import ReactionBuilder
        rxn_monod = ReactionBuilder.monod_aerobic_growth(
            substrate=substrate, biomass=biomass,
            mu_max_per_h=0.5, Ks_gL=5e-3, yield_gX_gS=0.36,
        )
        rxn_direct = ReactionBuilder.aerobic_growth(
            substrate_id=substrate.id, substrate_atoms=substrate.atoms,
            MW_substrate=float(substrate.MW),
            biomass_id=biomass.id, biomass_atoms=biomass.atoms,
            MW_biomass=float(biomass.MW),
            yield_gX_gS=0.36, rate_fn=lambda env: 0.0,
        )
        coeffs_m = {e.species.id: e.coefficient for e in rxn_monod.stoichiometry}
        coeffs_d = {e.species.id: e.coefficient for e in rxn_direct.stoichiometry}
        for sp in coeffs_d:
            assert coeffs_m[sp] == pytest.approx(coeffs_d[sp], abs=1e-10)

    def test_default_label(self, substrate, biomass):
        from PyOMES.reactions import ReactionBuilder
        rxn = ReactionBuilder.monod_aerobic_growth(
            substrate=substrate, biomass=biomass,
            mu_max_per_h=0.5, Ks_gL=5e-3, yield_gX_gS=0.36,
        )
        assert rxn.label == "monod_growth_AceticAcid"

    def test_custom_label(self, substrate, biomass):
        from PyOMES.reactions import ReactionBuilder
        rxn = ReactionBuilder.monod_aerobic_growth(
            substrate=substrate, biomass=biomass,
            mu_max_per_h=0.5, Ks_gL=5e-3, yield_gX_gS=0.36,
            label="my_growth",
        )
        assert rxn.label == "my_growth"

    def test_rate_fn_zero_at_zero_substrate(self, substrate, biomass):
        from PyOMES.reactions import ReactionBuilder, ReactionEnvironment
        rxn = ReactionBuilder.monod_aerobic_growth(
            substrate=substrate, biomass=biomass,
            mu_max_per_h=0.5, Ks_gL=5e-3, yield_gX_gS=0.36,
        )
        env = ReactionEnvironment(concentrations={"AceticAcid": 0.0, "Yeast": 0.1}, V_L=1.0, T_K=305.0)
        assert rxn.rate_fn(env) == 0.0

    def test_rate_fn_zero_at_zero_biomass(self, substrate, biomass):
        from PyOMES.reactions import ReactionBuilder, ReactionEnvironment
        rxn = ReactionBuilder.monod_aerobic_growth(
            substrate=substrate, biomass=biomass,
            mu_max_per_h=0.5, Ks_gL=5e-3, yield_gX_gS=0.36,
        )
        env = ReactionEnvironment(concentrations={"AceticAcid": 0.1, "Yeast": 0.0}, V_L=1.0, T_K=305.0)
        assert rxn.rate_fn(env) == 0.0

    def test_rate_fn_saturates_at_high_substrate(self, substrate, biomass):
        """At S >> Ks the rate approaches mu_max/Y * X/MW_S * V."""
        from PyOMES.reactions import ReactionBuilder, ReactionEnvironment
        MW_S = float(substrate.MW)
        MW_X = float(biomass.MW)
        mu_max, Ks, Y = 0.5, 5e-3, 0.36
        S_sat_gL = 10.0        # >> Ks=0.005 g/L
        X_gL     = 1.0
        rxn = ReactionBuilder.monod_aerobic_growth(
            substrate=substrate, biomass=biomass,
            mu_max_per_h=mu_max, Ks_gL=Ks, yield_gX_gS=Y,
        )
        env = ReactionEnvironment(
            concentrations={substrate.id: S_sat_gL / MW_S, biomass.id: X_gL / MW_X},
            V_L=1.0, T_K=305.0,
        )
        r = rxn.rate_fn(env)
        r_expected = mu_max / Y * X_gL / MW_S   # mol/h at V=1 L
        assert r == pytest.approx(r_expected, rel=1e-2)

    def test_rate_fn_half_saturation(self, substrate, biomass):
        """At S == Ks the rate is exactly half the saturation rate."""
        from PyOMES.reactions import ReactionBuilder, ReactionEnvironment
        MW_S = float(substrate.MW)
        MW_X = float(biomass.MW)
        mu_max, Ks, Y = 0.5, 5e-3, 0.36
        X_gL  = 1.0
        S_gL  = Ks        # S == Ks → μ = μ_max / 2
        rxn = ReactionBuilder.monod_aerobic_growth(
            substrate=substrate, biomass=biomass,
            mu_max_per_h=mu_max, Ks_gL=Ks, yield_gX_gS=Y,
        )
        env = ReactionEnvironment(
            concentrations={substrate.id: S_gL / MW_S, biomass.id: X_gL / MW_X},
            V_L=2.0, T_K=305.0,
        )
        r = rxn.rate_fn(env)
        r_half_sat = (mu_max / 2) / Y * X_gL / MW_S * 2.0
        assert r == pytest.approx(r_half_sat, rel=1e-9)

    def test_chno_balance(self, substrate):
        """CHNO balance mode is forwarded correctly."""
        from PyOMES.chemistry import Species
        from PyOMES.reactions import ReactionBuilder
        biomass_n = Species(id="Yeast_N", atoms={"C":1,"H":1.61,"O":0.56,"N":0.16}, charge=0, MW=26.868)
        rxn = ReactionBuilder.monod_aerobic_growth(
            substrate=substrate, biomass=biomass_n,
            mu_max_per_h=0.5, Ks_gL=5e-3, yield_gX_gS=0.36,
            balance="CHNO",
        )
        assert "NH3" in rxn.species_ids


# ── ReactionSystem ───────────────────────────────────────────────────

class TestReactionSystem:
    def test_parallel_reactions_sum(self):
        from PyOMES.reactions import ReactionEnvironment
        rset = _make_reaction_system()
        env = ReactionEnvironment(concentrations={"A":1.0,"B":2.0}, V_L=5.0)
        rates = rset.compute_rates(env)
        assert rates["liquid"]["A"] == pytest.approx(-10.0, abs=1e-10)
        assert rates["liquid"]["C"] == pytest.approx(+10.0, abs=1e-10)
        assert rates["liquid"]["B"] == pytest.approx(-5.0, abs=1e-10)
        assert rates["liquid"]["D"] == pytest.approx(+5.0, abs=1e-10)

    def test_series_reactions_net_rate(self):
        from PyOMES.reactions import KineticReaction, ReactionSystem, ReactionEnvironment
        rxn1 = KineticReaction([_entry("A",-1.0,{"C":1}),
                                _entry("B",+1.0,{"C":1})],
                               rate_fn=lambda env: 3.0, balance_elements=("C",))
        rxn2 = KineticReaction([_entry("B",-1.0,{"C":1}),
                                _entry("C",+1.0,{"C":1})],
                               rate_fn=lambda env: 1.0, balance_elements=("C",))
        rates = ReactionSystem([rxn1,rxn2]).compute_rates(ReactionEnvironment())
        assert rates["liquid"]["B"] == pytest.approx(+2.0, abs=1e-12)

    def test_multiple_organisms_shared_substrate(self):
        from PyOMES.reactions import KineticReaction, ReactionSystem, ReactionEnvironment
        rxn1 = KineticReaction([_entry("S",-1.0,{"C":2}),
                                _entry("X1",+0.5,{"C":1}),
                                _entry("CO2",+1.5,{"C":1,"O":2})],
                               rate_fn=lambda env: 0.3*env.S("S")*env.X("X1"), balance_elements=("C",))
        rxn2 = KineticReaction([_entry("S",-1.0,{"C":2}),
                                _entry("X2",+0.8,{"C":1}),
                                _entry("CO2",+1.2,{"C":1,"O":2})],
                               rate_fn=lambda env: 0.2*env.S("S")*env.X("X2"), balance_elements=("C",))
        env = ReactionEnvironment(concentrations={"S":5.0,"X1":1.0,"X2":2.0})
        rates = ReactionSystem([rxn1,rxn2]).compute_rates(env)
        r1,r2 = 0.3*5.0*1.0, 0.2*5.0*2.0
        assert rates["liquid"]["S"] == pytest.approx(-r1-r2, abs=1e-10)
        assert rates["liquid"]["X1"] == pytest.approx(0.5*r1, abs=1e-10)
        assert rates["liquid"]["X2"] == pytest.approx(0.8*r2, abs=1e-10)
        assert rates["liquid"]["CO2"] == pytest.approx(1.5*r1+1.2*r2, abs=1e-10)

    def test_len_and_iter(self):
        rset = _make_reaction_system()
        assert len(rset) == 2
        assert len(list(rset)) == 2


# ── EquilibriumReaction declaration / partition ──────────────────────

def _entry_sp(species, coeff, *, phase="liquid"):
    """Wrap an existing Species in a StoichiometryEntry."""
    from PyOMES.reactions import StoichiometryEntry
    return StoichiometryEntry(species=species, phase=phase, coefficient=coeff)


class TestEquilibriumDeclaration:

    def test_kinetic_reaction_carries_no_equilibrium_attrs(self):
        from PyOMES.reactions import KineticReaction
        rxn = _make_simple_reaction()
        assert isinstance(rxn, KineticReaction)
        # KineticReaction has no log_K / dH_J_per_mol / total_id attrs —
        # the kind discriminator is the type itself.
        assert not hasattr(rxn, "log_K")
        assert not hasattr(rxn, "dH_J_per_mol")
        assert not hasattr(rxn, "total_id")

    def test_equilibrium_construction(self):
        from PyOMES.reactions import EquilibriumReaction
        from PyOMES.chemistry.common_species import H_plus, OH_minus, H2O
        rxn = EquilibriumReaction(
            stoichiometry=[
                _entry_sp(H2O, -1.0), _entry_sp(H_plus, +1.0), _entry_sp(OH_minus, +1.0),
            ],
            log_K=-14.0,
            balance_elements=("H", "O"), label="water",
        )
        assert isinstance(rxn, EquilibriumReaction)
        assert rxn.log_K == pytest.approx(-14.0)
        # EquilibriumReaction has no rate_fn — equilibrium is algebraic.
        assert not hasattr(rxn, "rate_fn")

    def test_equilibrium_has_no_compute_rates(self):
        """EquilibriumReaction does not implement compute_rates —
        equilibrium reactions are routed to BisectionChemicalEquilibriumEngine.
        """
        from PyOMES.reactions import EquilibriumReaction
        from PyOMES.chemistry.common_species import H_plus, OH_minus, H2O
        rxn = EquilibriumReaction(
            stoichiometry=[
                _entry_sp(H2O, -1.0), _entry_sp(H_plus, +1.0), _entry_sp(OH_minus, +1.0),
            ],
            log_K=-14.0,
            balance_elements=("H", "O"),
        )
        assert not hasattr(rxn, "compute_rates")

    def test_equilibrium_without_log_K_rejected(self):
        from PyOMES.reactions import EquilibriumReaction
        from PyOMES.chemistry.common_species import H_plus, OH_minus, H2O
        with pytest.raises(ValueError):
            EquilibriumReaction(
                stoichiometry=[
                    _entry_sp(H2O, -1.0), _entry_sp(H_plus, +1.0), _entry_sp(OH_minus, +1.0),
                ],
                balance_elements=("H", "O"),
            )

    def test_is_cross_phase_property(self):
        """``EquilibriumReaction.is_cross_phase`` reflects whether the
        stoichiometry spans more than one phase.
        """
        from PyOMES.reactions import EquilibriumReaction
        from PyOMES.chemistry.common_species import CO2

        single_phase = EquilibriumReaction(
            stoichiometry=[
                _entry_sp(CO2, -1.0),
                _entry_sp(CO2, +1.0),
            ],
            log_K=0.0,
            balance_elements=("C", "O"),
        )
        assert single_phase.is_cross_phase is False

        cross_phase = EquilibriumReaction(
            stoichiometry=[
                _entry_sp(CO2, -1.0, phase="gas"),
                _entry_sp(CO2, +1.0, phase="liquid"),
            ],
            # log_K omitted for cross-phase
            balance_elements=("C", "O"),
        )
        assert cross_phase.is_cross_phase is True

    def test_cross_phase_equilibrium_log_K_optional(self):
        """Cross-phase equilibrium reactions may omit ``log_K``
        because the partition constant lives on the consuming
        gas-liquid link (Henry's law), not on the reaction.
        """
        from PyOMES.reactions import EquilibriumReaction
        from PyOMES.chemistry.common_species import CO2
        rxn = EquilibriumReaction(
            stoichiometry=[
                _entry_sp(CO2, -1.0, phase="gas"),
                _entry_sp(CO2, +1.0, phase="liquid"),
            ],
            balance_elements=("C", "O"),
            label="CO2 partition",
        )
        assert rxn.log_K is None
        assert rxn.is_cross_phase is True


class TestReactionSystemBucketing:
    """ReactionSystem pre-buckets its reactions by type at construction.
    The public ``kinetic_reactions`` / ``single_phase_equilibria`` /
    ``cross_phase_equilibria`` / ``blackbox_models`` projections
    replace the deleted ``partition()`` API and the
    ``has_equilibrium`` warning.
    """

    def _kinetic_rxn(self, label="k"):
        from PyOMES.reactions import KineticReaction
        return KineticReaction(
            stoichiometry=[_entry("A", -1.0, {"C": 1}), _entry("B", +1.0, {"C": 1})],
            rate_fn=lambda env: 0.0, balance_elements=("C",), label=label,
        )

    def _equilibrium_rxn(self, label="eq"):
        from PyOMES.reactions import EquilibriumReaction
        from PyOMES.chemistry.common_species import H_plus, OH_minus, H2O
        return EquilibriumReaction(
            stoichiometry=[
                _entry_sp(H2O, -1.0), _entry_sp(H_plus, +1.0), _entry_sp(OH_minus, +1.0),
            ],
            log_K=-14.0,
            balance_elements=("H", "O"), label=label,
        )

    def test_mixed_set_pre_buckets(self):
        from PyOMES.reactions import ReactionSystem
        system = ReactionSystem([
            self._kinetic_rxn("k1"),
            self._equilibrium_rxn("e1"),
            self._kinetic_rxn("k2"),
        ], label="mix")
        assert [r.label for r in system.kinetic_reactions] == ["k1", "k2"]
        assert [r.label for r in system.single_phase_equilibria] == ["e1"]
        assert system.cross_phase_equilibria == []
        assert system.blackbox_models == []
        assert len(system) == 3

    def test_kinetic_only_system_empty_equilibrium_bucket(self):
        from PyOMES.reactions import ReactionSystem
        kinetic_only = ReactionSystem(
            [self._kinetic_rxn(), self._kinetic_rxn()]
        )
        assert len(kinetic_only.kinetic_reactions) == 2
        assert kinetic_only.single_phase_equilibria == []
        assert kinetic_only.cross_phase_equilibria == []

    def test_equilibria_silently_skipped_in_compute_rates(self):
        """ReactionSystem.compute_rates iterates only kinetic +
        black-box buckets. Equilibria are silently skipped (they
        have no rate) — no RuntimeError, no UserWarning."""
        from PyOMES.reactions import ReactionSystem, ReactionEnvironment
        system = ReactionSystem([
            self._kinetic_rxn("k1"),
            self._equilibrium_rxn("e1"),
        ])
        # No raise — equilibria are pre-bucketed out of the rate path.
        env = ReactionEnvironment(concentrations={"A": 1.0}, V_L=1.0)
        rates = system.compute_rates(env)
        # Only kinetic species in the result.
        assert set(rates.get("liquid", {}).keys()) <= {"A", "B"}


class TestReactionEnvironmentHasPH:

    def test_has_pH_true_when_set(self):
        from PyOMES.reactions import ReactionEnvironment
        assert ReactionEnvironment(pH=7.0).has_pH is True

    def test_has_pH_false_when_none(self):
        from PyOMES.reactions import ReactionEnvironment
        assert ReactionEnvironment().has_pH is False

# ── BlackBoxReactionModel ────────────────────────────────────────────

class TestBlackBoxReactionModel:
    def test_balanced_model_no_warning(self):
        from PyOMES.reactions import ReactionEnvironment, MassBalanceWarning
        model = _make_balanced_blackbox()
        env = ReactionEnvironment(concentrations={"A":2.0})
        with warnings.catch_warnings():
            warnings.simplefilter("error", MassBalanceWarning)
            rates = model.compute_rates(env)
        assert rates["liquid"]["A"] == pytest.approx(-1.0, abs=1e-10)
        assert rates["liquid"]["B"] == pytest.approx(+1.0, abs=1e-10)

    def test_imbalanced_model_warns(self):
        from PyOMES.reactions import BlackBoxReactionModel, FluxEntry, ReactionEnvironment, MassBalanceWarning
        class BadModel:
            def solve(self, inputs): return {"produce_B": 1.0}
        model = BlackBoxReactionModel(
            external_model=BadModel(),
            flux_mapping={"produce_B": FluxEntry("B","liquid",{"C":1},sign=+1)},
            balance_elements=("C",), on_imbalance="warn")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            model.compute_rates(ReactionEnvironment())
            mass_warns = [x for x in w if issubclass(x.category, MassBalanceWarning)]
            assert len(mass_warns) >= 1

    def test_imbalanced_model_raises(self):
        from PyOMES.reactions import BlackBoxReactionModel, FluxEntry, ReactionEnvironment, MassBalanceError
        class BadModel:
            def solve(self, inputs): return {"produce_B": 1.0}
        model = BlackBoxReactionModel(
            external_model=BadModel(),
            flux_mapping={"produce_B": FluxEntry("B","liquid",{"C":1},sign=+1)},
            balance_elements=("C",), on_imbalance="raise")
        with pytest.raises(MassBalanceError):
            model.compute_rates(ReactionEnvironment())

    def test_imbalanced_model_ignore(self):
        from PyOMES.reactions import BlackBoxReactionModel, FluxEntry, ReactionEnvironment
        class BadModel:
            def solve(self, inputs): return {"produce_B": 1.0}
        model = BlackBoxReactionModel(
            external_model=BadModel(),
            flux_mapping={"produce_B": FluxEntry("B","liquid",{"C":1},sign=+1)},
            balance_elements=("C",), on_imbalance="ignore")
        model.compute_rates(ReactionEnvironment())
        assert len(model.imbalance_history) == 1

    def test_sign_convention(self):
        from PyOMES.reactions import BlackBoxReactionModel, FluxEntry, ReactionEnvironment
        class SignModel:
            def solve(self, inputs): return {"flux_a": 5.0}
        model = BlackBoxReactionModel(
            external_model=SignModel(),
            flux_mapping={"flux_a": FluxEntry("A","liquid",{"C":1},sign=-1.0)},
            balance_elements=())
        rates = model.compute_rates(ReactionEnvironment())
        assert rates["liquid"]["A"] == pytest.approx(-5.0, abs=1e-12)

# ── ReactionEnvironment ──────────────────────────────────────────────

class TestReactionEnvironment:
    def test_S_accessor(self):
        from PyOMES.reactions import ReactionEnvironment
        env = ReactionEnvironment(concentrations={"Glc":5.0,"Ac":1.0})
        assert env.S("Glc") == 5.0
        assert env.S("Missing") == 0.0

    def test_X_accessor(self):
        from PyOMES.reactions import ReactionEnvironment
        env = ReactionEnvironment(concentrations={"Yeast":3.0})
        assert env.X("Yeast") == 3.0

    def test_prop_accessor(self):
        from PyOMES.reactions import ReactionEnvironment
        env = ReactionEnvironment(properties={"DO_mol_L":1.3e-4})
        assert abs(env.prop("DO_mol_L") - 1.3e-4) < 1e-10

    def test_frozen(self):
        from PyOMES.reactions import ReactionEnvironment
        env = ReactionEnvironment(T_K=300.0)
        try:
            env.T_K = 310.0
            assert False, "Should have raised"
        except AttributeError:
            pass

# ── Protocol compliance ──────────────────────────────────────────────

class TestProtocolCompliance:
    def test_kinetic_reaction_is_reaction_model(self):
        from PyOMES.reactions import ReactionModel
        assert isinstance(_make_simple_reaction(), ReactionModel)

    def test_reaction_system_is_reaction_model(self):
        from PyOMES.reactions import ReactionModel
        assert isinstance(_make_reaction_system(), ReactionModel)

    def test_blackbox_is_reaction_model(self):
        from PyOMES.reactions import BlackBoxReactionModel, ReactionModel
        model = BlackBoxReactionModel(
            external_model=type("M",(),{"solve":lambda self,i:{}})(), flux_mapping={})
        assert isinstance(model, ReactionModel)
