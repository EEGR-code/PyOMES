# -*- coding: utf-8 -*-
"""A model whose species and reactions are all defined by the user.

Every ``Species`` here is built in this file and every reaction is given
the species it names. Nothing is imported from a database module. This is
the baseline way to define a model; if it fails, the routes that take
species from a database do not matter.

The checks: the model uses only the ``Species`` objects built here, it
advances, and it conserves each element. The stirred-tank template is
checked the same way: built from species and Henry constants given here,
with no database.
"""
import math

import pytest

from PyOMES.chemistry import Species
from PyOMES.core import ControlVolume, LiquidPhase
from PyOMES.reactions import EquilibriumReaction, KineticReaction, ReactionSystem


# ── Species, built by the user ───────────────────────────────────────────────

H2O           = Species(id="H2O",        atoms={"H": 2, "O": 1})
H_PLUS        = Species(id="H+",         atoms={"H": 1},                 charge=+1)
OH_MINUS      = Species(id="OH-",        atoms={"O": 1, "H": 1},         charge=-1)
CO2           = Species(id="CO2",        atoms={"C": 1, "O": 2})
HCO3_MINUS    = Species(id="HCO3-",      atoms={"H": 1, "C": 1, "O": 3}, charge=-1)
O2            = Species(id="O2",         atoms={"O": 2})
ACETIC_ACID   = Species(id="AceticAcid", atoms={"C": 2, "H": 4, "O": 2})
ACETATE_MINUS = Species(id="Acetate-",   atoms={"C": 2, "H": 3, "O": 2}, charge=-1)

USER_SPECIES = {
    s.id: s for s in (
        H2O, H_PLUS, OH_MINUS, CO2, HCO3_MINUS, O2, ACETIC_ACID, ACETATE_MINUS,
    )
}

V_L = 1.0


def _build_reactions():
    return ReactionSystem([
        EquilibriumReaction(
            "H2O,aq <-> H+,aq + OH-,aq",
            species=USER_SPECIES, log_K=-14.0,
            balance_elements=("H", "O"), label="water",
        ),
        EquilibriumReaction(
            "CO2,aq + H2O,aq <-> HCO3-,aq + H+,aq",
            species=USER_SPECIES, log_K=-6.35,
            balance_elements=("C", "H", "O"), label="carbonate",
        ),
        EquilibriumReaction(
            "AceticAcid,aq <-> Acetate-,aq + H+,aq",
            species=USER_SPECIES, log_K=-4.76,
            balance_elements=("C", "H", "O"), label="acetate",
        ),
        KineticReaction(
            "AceticAcid,aq + 2 O2,aq -> 2 CO2,aq + 2 H2O,aq",
            species=USER_SPECIES,
            rate_fn=lambda env: 50.0 * env.S("AceticAcid") * env.S("O2") * env.V_L,
            balance_elements=("C", "H", "O"), label="acetate_oxidation",
        ),
    ], label="user_defined")


def _build_cv():
    liquid = LiquidPhase(
        n_mol={
            "H2O": 55.5 * V_L,
            "H+": 1e-7 * V_L,
            "OH-": 1e-7 * V_L,
            "AceticAcid": 1e-2 * V_L,
            "O2": 2e-3 * V_L,
            "CO2": 1e-3 * V_L,
        },
        V_L=V_L, T_K=298.15,
    )
    return ControlVolume(
        phases={"liquid": liquid},
        reaction_system=_build_reactions(),
        label="user_defined",
    )


def _element_totals(cv):
    """Element totals over the species built here (other ids are checked
    separately, in ``test_every_n_mol_id_is_a_species_built_here``)."""
    totals = {}
    for sp_id, n in cv.total_mol().items():
        if sp_id not in USER_SPECIES:
            continue
        for element, count in USER_SPECIES[sp_id].atoms.items():
            totals[element] = totals.get(element, 0.0) + n * count
    return totals


class TestUserDefinedModel:

    def test_reactions_use_only_the_species_built_here(self):
        rs = _build_reactions()
        for rxn in rs.reactions:
            for entry in rxn.stoichiometry:
                assert entry.species is USER_SPECIES[entry.species.id], (
                    f"{rxn.label}: {entry.species.id!r} is not the object built here"
                )

    def test_monitor_knows_only_the_species_built_here(self):
        cv = _build_cv()
        registry = cv._conservation_monitor._species_registry
        assert set(registry) <= set(USER_SPECIES)
        for sp_id, sp in registry.items():
            assert sp is USER_SPECIES[sp_id]

    def test_every_n_mol_id_is_a_species_built_here(self):
        cv = _build_cv()
        for _ in range(20):
            cv.advance(dt_h=0.05)
        assert set(cv.total_mol()) <= set(USER_SPECIES)

    def test_declared_conjugate_base_is_written_back(self):
        cv = _build_cv()
        cv.advance(dt_h=0.0)
        assert cv["liquid"].n_mol.get("Acetate-", 0.0) > 0.0

    def test_advances_to_a_finite_pH(self):
        cv = _build_cv()
        for _ in range(20):
            cv.advance(dt_h=0.05)
        pH = float(cv["liquid"].pH)
        assert math.isfinite(pH)
        assert 2.0 < pH < 7.0

    def test_acetate_is_oxidised(self):
        cv = _build_cv()
        acetate_0 = cv["liquid"].n_mol["AceticAcid"]
        for _ in range(20):
            cv.advance(dt_h=0.05)
        total = cv["liquid"].n_mol.get("AceticAcid", 0.0) + cv["liquid"].n_mol.get("Acetate-", 0.0)
        assert total < acetate_0

    _DRIFT = pytest.mark.xfail(strict=True, reason=(
        "The Bisection engine treats H2O as a fixed solvent: the H2O the "
        "declared CO2 + H2O <-> HCO3- + H+ consumes is not taken from n_mol, "
        "so each HCO3- formed adds 2 H and 1 O (+2.5e-8 relative over 1 h)."
    ))

    @pytest.mark.parametrize("element", [
        "C",
        pytest.param("H", marks=_DRIFT),
        pytest.param("O", marks=_DRIFT),
    ])
    def test_conserves_each_element(self, element):
        cv = _build_cv()
        cv.advance(dt_h=0.0)  # speciate first, so the baseline is at equilibrium
        before = _element_totals(cv)[element]
        for _ in range(20):
            cv.advance(dt_h=0.05)
        after = _element_totals(cv)[element]
        assert after == pytest.approx(before, rel=1e-9)


# ── The stirred-tank template, with no database ──────────────────────────────

N2 = Species(id="N2", atoms={"N": 2})
YEAST_ATOMS, YEAST_MW = {"C": 1, "H": 1.61, "O": 0.56}, 24.626


def _build_tank():
    from PyOMES.templates.stirred_tank import StirredTankBuilder, TransferConfig
    return (
        StirredTankBuilder()
        .vessel(V_total_L=10.0, T_K=305.15)
        .initial_gas({"O2": 0.21, "N2": 0.79})
        .no_gas_feed()
        .transfer(TransferConfig(species={}))
        .transfer_species("O2", henry_mol_L_atm=1.2e-3)
        .transfer_species("CO2", henry_mol_L_atm=3.3e-2)
        .transfer_species("N2", henry_mol_L_atm=6.1e-4)
        .chemistry(species=[O2, CO2, H2O, N2])
        .organism("Yeast", atoms=YEAST_ATOMS, MW=YEAST_MW)
        .substrate("AceticAcid", atoms=dict(ACETIC_ACID.atoms), MW=ACETIC_ACID.MW)
        .build()
    )


class TestUserDefinedStirredTank:
    """The template builds and runs from the user's species and Henry
    constants alone."""

    def test_growth_uses_the_gases_built_here(self):
        rxn = _build_tank().reaction_system
        by_id = {e.species.id: e.species for e in rxn.stoichiometry}
        assert by_id["O2"] is O2 and by_id["CO2"] is CO2 and by_id["H2O"] is H2O

    def test_runs_and_grows(self):
        cv = _build_tank()
        V = cv["liquid"].V_L
        cv["liquid"].n_mol["AceticAcid"] = 1e-2 * V
        cv["liquid"].n_mol["Yeast"] = 1e-3 * V
        for _ in range(10):
            cv.advance(dt_h=0.1)
        assert cv["liquid"].n_mol["Yeast"] > 1e-3 * V


class TestStirredTankWithOtherGasIds:
    """The template has no fixed gas ids: gases named differently by the
    model are used throughout (headspace, transfer, growth)."""

    def _build(self):
        from PyOMES.templates.stirred_tank import StirredTankBuilder, TransferConfig
        o2 = Species(id="O2_aq", atoms={"O": 2})
        co2 = Species(id="CO2_aq", atoms={"C": 1, "O": 2})
        n2 = Species(id="N2_g", atoms={"N": 2})
        h2o = Species(id="H2O_l", atoms={"H": 2, "O": 1})
        cv = (
            StirredTankBuilder()
            .vessel(V_total_L=10.0, T_K=305.15)
            .initial_gas({"O2_aq": 0.21, "N2_g": 0.79})
            .no_gas_feed()
            .transfer(TransferConfig(species={}))
            .transfer_species("O2_aq", henry_mol_L_atm=1.2e-3)
            .transfer_species("CO2_aq", henry_mol_L_atm=3.3e-2)
            .transfer_species("N2_g", henry_mol_L_atm=6.1e-4)
            .chemistry(species=[o2, co2, n2, h2o])
            .organism("Yeast", atoms=YEAST_ATOMS, MW=YEAST_MW,
                      o2_id="O2_aq", co2_id="CO2_aq", h2o_id="H2O_l")
            .substrate("AceticAcid", atoms=dict(ACETIC_ACID.atoms), MW=ACETIC_ACID.MW)
            .build()
        )
        return cv, (o2, co2, n2, h2o)

    def test_gas_phase_holds_only_the_declared_gases(self):
        cv, _ = self._build()
        assert set(cv["gas"].n_mol) == {"O2_aq", "CO2_aq", "N2_g"}

    def test_growth_uses_the_renamed_gases(self):
        cv, (o2, co2, _, h2o) = self._build()
        by_id = {e.species.id: e.species for e in cv.reaction_system.stoichiometry}
        assert by_id["O2_aq"] is o2 and by_id["CO2_aq"] is co2 and by_id["H2O_l"] is h2o
        assert not {"O2", "CO2", "H2O", "N2"} & set(by_id)
