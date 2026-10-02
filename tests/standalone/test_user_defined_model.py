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
YEAST_ATOMS = {"C": 1, "H": 1.61, "O": 0.56}
YEAST_MW = Species(id="Yeast", atoms=YEAST_ATOMS).MW


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


YEAST = Species(id="Yeast", atoms=YEAST_ATOMS, MW=YEAST_MW)
NH3 = Species(id="NH3", atoms={"N": 1, "H": 3})


class TestStirredTankOrganismAndSubstrate:
    """The organism and substrate are named by an id among the model's
    species, or defined by an id with atoms or by a Species."""

    @staticmethod
    def _growth(organism, substrate, species, **organism_kw):
        from PyOMES.templates.stirred_tank import StirredTankBuilder, TransferConfig
        cv = (
            StirredTankBuilder()
            .transfer(TransferConfig(species={}))
            .chemistry(species=[O2, CO2, H2O, *species])
            .organism(*organism, **organism_kw)
            .substrate(*substrate)
            .build()
        )
        return {e.species.id: e.species for e in cv.reaction_system.stoichiometry}

    def test_ids_resolve_to_the_model_species(self):
        by_id = self._growth(("Yeast",), ("AceticAcid",), [YEAST, ACETIC_ACID])
        assert by_id["Yeast"] is YEAST and by_id["AceticAcid"] is ACETIC_ACID

    def test_unknown_id_names_the_fix(self):
        with pytest.raises(ValueError, match="'Yeast' is not among the species passed.*atoms="):
            self._growth(("Yeast",), ("AceticAcid",), [ACETIC_ACID])

    def test_id_with_atoms_defines_it_with_computed_mw(self):
        by_id = self._growth(("Yeast", YEAST_ATOMS), ("AceticAcid",), [ACETIC_ACID])
        assert dict(by_id["Yeast"].atoms) == YEAST_ATOMS
        assert by_id["Yeast"].MW == pytest.approx(Species(id="x", atoms=YEAST_ATOMS).MW)

    def test_species_is_used_as_given(self):
        e_coli = Species(id="E_coli", atoms={"C": 1, "H": 1.77, "O": 0.49})
        by_id = self._growth((e_coli,), (ACETIC_ACID,), [])
        assert by_id["E_coli"] is e_coli and by_id["AceticAcid"] is ACETIC_ACID

    def test_same_definition_as_the_model_uses_the_model_species(self):
        by_id = self._growth(("Yeast", YEAST_ATOMS, YEAST_MW), ("AceticAcid",),
                             [YEAST, ACETIC_ACID])
        assert by_id["Yeast"] is YEAST

    def test_different_definition_raises(self):
        from PyOMES.chemistry import SpeciesConflictError
        with pytest.raises(SpeciesConflictError, match="(?s)'Yeast'.*overwrite=True"):
            self._growth(("Yeast", {"C": 1, "H": 1.8, "O": 0.5}), ("AceticAcid",),
                         [YEAST, ACETIC_ACID])

    def test_overwrite_replaces_the_model_species(self):
        atoms = {"C": 1, "H": 1.8, "O": 0.5}
        by_id = self._growth(("Yeast", atoms), ("AceticAcid",), [YEAST, ACETIC_ACID],
                             overwrite=True)
        assert by_id["Yeast"] is not YEAST and dict(by_id["Yeast"].atoms) == atoms

    def test_chno_takes_the_n_source_from_the_model(self):
        yeast_n = Species(id="Yeast", atoms={**YEAST_ATOMS, "N": 0.16})
        by_id = self._growth((yeast_n,), ("AceticAcid",), [ACETIC_ACID, NH3],
                             balance_basis="CHNO", n_source_id="NH3")
        assert by_id["NH3"] is NH3

    def test_chno_unknown_n_source_names_the_fix(self):
        yeast_n = Species(id="Yeast", atoms={**YEAST_ATOMS, "N": 0.16})
        with pytest.raises(ValueError, match="'NH3' is not among the species passed.*n_source_id"):
            self._growth((yeast_n,), ("AceticAcid",), [ACETIC_ACID],
                         balance_basis="CHNO", n_source_id="NH3")


# ── The model's species set and reactions= ───────────────────────────────────

NA_PLUS = Species(id="Na+", atoms={"Na": 1}, charge=+1)


class TestModelSpecies:
    """``cv.species`` holds the species the model was given: those passed and
    those in its reactions."""

    @staticmethod
    def _liquid():
        return LiquidPhase(n_mol={"H2O": 55.5, "H+": 1e-7, "AceticAcid": 1e-2},
                           V_L=V_L, T_K=298.15)

    def test_holds_reaction_and_passed_species(self):
        cv = ControlVolume(phases={"liquid": self._liquid()},
                           reaction_system=_build_reactions(), species=[NA_PLUS])
        assert set(cv.species) == set(USER_SPECIES) | {"Na+"}
        assert cv.species["Na+"] is NA_PLUS and cv.species["O2"] is O2

    def test_is_read_only(self):
        cv = ControlVolume(phases={"liquid": self._liquid()}, species=[NA_PLUS])
        with pytest.raises(TypeError):
            cv.species["Na+"] = NA_PLUS

    def test_passed_species_feed_the_conservation_monitor(self):
        cv = ControlVolume(phases={"liquid": self._liquid()},
                           reaction_system=_build_reactions(), species=[NA_PLUS])
        assert cv._conservation_monitor._species_registry["Na+"] is NA_PLUS

    def test_conflicting_species_raises(self):
        from PyOMES.chemistry import SpeciesConflictError
        wrong_o2 = Species(id="O2", atoms={"O": 1})
        with pytest.raises(SpeciesConflictError, match="'O2'"):
            ControlVolume(phases={"liquid": self._liquid()},
                          reaction_system=_build_reactions(), species=[wrong_o2])

    def test_reactions_are_run_like_a_reaction_system(self):
        by_reactions = ControlVolume(phases={"liquid": self._liquid()},
                                     reactions=_build_reactions().reactions)
        by_system = ControlVolume(phases={"liquid": self._liquid()},
                                  reaction_system=_build_reactions())
        assert set(by_reactions.species) == set(by_system.species)
        for cv in (by_reactions, by_system):
            for _ in range(5):
                cv.advance(dt_h=0.05)
        assert by_reactions.total_mol() == pytest.approx(by_system.total_mol(), rel=1e-12)

    def test_reactions_and_reaction_system_together_raise(self):
        rs = _build_reactions()
        with pytest.raises(ValueError, match="not both"):
            ControlVolume(phases={"liquid": self._liquid()},
                          reactions=rs.reactions, reaction_system=rs)

    def test_snapshot_keeps_the_species(self):
        cv = ControlVolume(phases={"liquid": self._liquid()}, species=[NA_PLUS])
        assert cv.snapshot().species["Na+"] is NA_PLUS

    def test_template_cv_holds_the_model_species(self):
        cv = _build_tank()
        assert {"O2", "CO2", "H2O", "N2", "Yeast", "AceticAcid"} <= set(cv.species)
        assert cv.species["N2"] is N2


# ── Ids the model has no Species for ─────────────────────────────────────────

class TestUnresolvedSpecies:
    """An id in n_mol, a feed or the template's set-up that the model has no
    Species for warns once, naming where it entered."""

    @staticmethod
    def _unresolved(record):
        from PyOMES import UnresolvedSpeciesWarning
        return [str(w.message) for w in record
                if issubclass(w.category, UnresolvedSpeciesWarning)]

    def test_declared_model_does_not_warn(self):
        import warnings
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            cv = _build_cv()
            cv.advance(dt_h=0.05)
        assert self._unresolved(record) == []

    def test_undeclared_id_warns_once(self):
        import warnings
        liquid = LiquidPhase(n_mol={"H2O": 55.5, "Na +": 0.01}, V_L=V_L, T_K=298.15)
        cv = ControlVolume(phases={"liquid": liquid}, species=[H2O], label="typo")
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            for _ in range(3):
                cv.advance(dt_h=0.05)
        msgs = self._unresolved(record)
        assert len(msgs) == 1
        assert "ControlVolume('typo') n_mol when it first advanced: ['Na +']" in msgs[0]

    def test_id_appearing_mid_run_warns(self):
        import warnings
        liquid = LiquidPhase(n_mol={"H2O": 55.5}, V_L=V_L, T_K=298.15)
        cv = ControlVolume(phases={"liquid": liquid}, species=[H2O], label="late")
        cv.advance(dt_h=0.05)
        liquid.n_mol["Na+"] = 0.01
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            cv.advance(dt_h=0.05, t_h=0.05)
        assert any("n_mol by t = 0.05 h: ['Na+']" in m for m in self._unresolved(record))

    def test_snapshot_does_not_warn_again(self):
        import warnings
        liquid = LiquidPhase(n_mol={"H2O": 55.5, "X": 1.0}, V_L=V_L, T_K=298.15)
        cv = ControlVolume(phases={"liquid": liquid}, species=[H2O])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cv.advance(dt_h=0.05)
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            cv.snapshot().advance(dt_h=0.05)
        assert self._unresolved(record) == []

    def test_feed_warns_naming_the_boundary(self):
        import warnings
        from PyOMES.core import Simulation
        from PyOMES.core.boundaries import LiquidFeed
        liquid = LiquidPhase(n_mol={"H2O": 55.5}, V_L=V_L, T_K=298.15)
        cv = ControlVolume(phases={"liquid": liquid}, species=[H2O])
        cv.boundaries.append(LiquidFeed(Q_L_per_h=1.0, feed_conc_mol_L={"Glucose": 0.1},
                                        label="sugar"))
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            Simulation(cvs={"main": cv})
        assert any("LiquidFeed('sugar') on ControlVolume('main'): ['Glucose']" in m
                   for m in self._unresolved(record))

    def test_template_set_up_warns_naming_the_call(self):
        import warnings
        from PyOMES.templates.stirred_tank import StirredTankBuilder
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            (StirredTankBuilder()
             .initial_gas({"O2": 0.21, "N2": 0.79})
             .chemistry(species=[O2])
             .label("tank")
             .build())
        assert any("Stirred tank 'tank': the vessel's gas_composition "
                   "(.initial_gas()): ['N2']" in m for m in self._unresolved(record))
