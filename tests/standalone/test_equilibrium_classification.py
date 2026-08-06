# -*- coding: utf-8 -*-
"""Tests for CP2: classify_equilibrium_constraint() and the single flat-list
ingestion path across BisectionChemicalEquilibriumEngine.from_reactions,
NRChemicalEquilibriumEngine.from_reactions, and build_tableau().

Covers:
- classify_equilibrium_constraint() on EquilibriumReaction (single-phase and
  cross-phase), HenryEquilibrium, RaoultEquilibrium, KspEquilibrium, and the
  empty-stoichiometry error case.
- BisectionChemicalEquilibriumEngine.from_reactions(): cross-phase items are exposed via
  cross_phase_constraints instead of silently dropped.
- NRChemicalEquilibriumEngine.from_reactions(): solid-liquid items are auto-classified
  from one flat list; the deprecated precipitation_reactions= kwarg still
  works and warns.
- build_tableau(): still filters gas-liquid/solid-liquid items out of the
  acid-base graph, now via the shared classifier.
"""
from __future__ import annotations

import warnings

import pytest


# ── Shared fixtures ───────────────────────────────────────────────────────────

def _water_reaction():
    from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
    from PyOMES.chemistry.common_species import H2O, H_plus, OH_minus
    return EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=H2O, phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=H_plus, phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=OH_minus, phase="liquid", coefficient=+1.0),
        ],
        log_K=-14.0, label="water",
    )


def _co2_acid_reaction():
    from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
    from PyOMES.chemistry.common_species import CO2, H2O, HCO3_minus, H_plus
    return EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=CO2, phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=H2O, phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=HCO3_minus, phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=H_plus, phase="liquid", coefficient=+1.0),
        ],
        log_K=-6.35, total_id="CO2", label="co2_first",
        balance_elements=("C", "H", "O"),
    )


def _co2_gas_liquid_declaration():
    """Cross-phase EquilibriumReaction: CO2(gas) <-> CO2(liquid), no log_K
    (a partition declaration, consumed by KineticGasLiquidLink)."""
    from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
    from PyOMES.chemistry.common_species import CO2
    return EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=CO2, phase="gas", coefficient=-1.0),
            StoichiometryEntry(species=CO2, phase="liquid", coefficient=+1.0),
        ],
        label="co2_partition",
    )


def _henry_co2():
    from PyOMES.chemistry import HenryEquilibrium
    return HenryEquilibrium(
        H_ref=3.4e-4, dlnH=2400.0, gas_species="CO2", liquid_species="CO2",
    )


def _calcite_reaction():
    """Calcite dissolution: CaCO3(s) <-> Ca++ + CO3--  log_K = -8.48 (Ksp at 25 C)."""
    from PyOMES.chemistry.common_species import Ca_plus_plus, CO3_2minus
    from PyOMES.chemistry.species import Species
    from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
    CaCO3 = Species(id="CaCO3", atoms={"Ca": 1, "C": 1, "O": 3}, charge=0, MW=100.086)
    return EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=CaCO3, phase="solid", coefficient=-1.0),
            StoichiometryEntry(species=Ca_plus_plus, phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=CO3_2minus, phase="liquid", coefficient=+1.0),
        ],
        log_K=-8.48, label="calcite",
    )


def _species(result):
    return {**result.species_mol_L, **result.extra}


# ── classify_equilibrium_constraint() ────────────────────────────────────────

class TestClassifyEquilibriumConstraint:
    def test_water_reaction_is_acid_base(self):
        from PyOMES.reactions.equilibrium import classify_equilibrium_constraint
        assert classify_equilibrium_constraint(_water_reaction()) == "acid_base"

    def test_single_phase_acid_reaction_is_acid_base(self):
        from PyOMES.reactions.equilibrium import classify_equilibrium_constraint
        assert classify_equilibrium_constraint(_co2_acid_reaction()) == "acid_base"

    def test_cross_phase_equilibrium_reaction_is_gas_liquid(self):
        from PyOMES.reactions.equilibrium import classify_equilibrium_constraint
        assert classify_equilibrium_constraint(_co2_gas_liquid_declaration()) == "gas_liquid"

    def test_henry_equilibrium_is_gas_liquid(self):
        from PyOMES.reactions.equilibrium import classify_equilibrium_constraint
        assert classify_equilibrium_constraint(_henry_co2()) == "gas_liquid"

    def test_raoult_equilibrium_is_gas_liquid(self):
        from PyOMES.chemistry import RaoultEquilibrium
        from PyOMES.reactions.equilibrium import classify_equilibrium_constraint
        assert classify_equilibrium_constraint(RaoultEquilibrium()) == "gas_liquid"

    def test_equilibrium_reaction_calcite_is_solid_liquid(self):
        from PyOMES.reactions.equilibrium import classify_equilibrium_constraint
        assert classify_equilibrium_constraint(_calcite_reaction()) == "solid_liquid"

    def test_ksp_equilibrium_is_solid_liquid(self):
        from PyOMES.chemistry import KspEquilibrium
        from PyOMES.chemistry.species import Species
        from PyOMES.reactions.stoichiometry import StoichiometryEntry
        from PyOMES.reactions.equilibrium import classify_equilibrium_constraint
        MineralX_solid = Species(id="MineralX(s)", atoms={"Mn": 1, "O": 1}, charge=0)
        MineralX_aq = Species(id="MineralX", atoms={"Mn": 1, "O": 1}, charge=0)
        ksp = KspEquilibrium(
            stoichiometry=[
                StoichiometryEntry(species=MineralX_solid, phase="solid", coefficient=-1.0),
                StoichiometryEntry(species=MineralX_aq, phase="liquid", coefficient=+1.0),
            ],
            Ksp=1.0e-4,
        )
        assert classify_equilibrium_constraint(ksp) == "solid_liquid"

    def test_solid_priority_over_gas_liquid(self):
        """A solid+gas+liquid mix still classifies as solid_liquid (solid wins)."""
        from PyOMES.chemistry.species import Species
        from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
        from PyOMES.reactions.equilibrium import classify_equilibrium_constraint
        MineralX_solid = Species(id="MineralX(s)", atoms={"Mn": 1}, charge=0)
        MineralX_gas = Species(id="MineralX(g)", atoms={"Mn": 1}, charge=0)
        MineralX_aq = Species(id="MineralX", atoms={"Mn": 1}, charge=0)
        rxn = EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=MineralX_solid, phase="solid", coefficient=-1.0),
                StoichiometryEntry(species=MineralX_gas, phase="gas", coefficient=+1.0),
                StoichiometryEntry(species=MineralX_aq, phase="liquid", coefficient=+1.0),
            ],
            log_K=-4.0,
            balance_elements=(),
        )
        assert classify_equilibrium_constraint(rxn) == "solid_liquid"

    def test_empty_stoichiometry_raises_value_error(self):
        from PyOMES.chemistry import HenryEquilibrium
        from PyOMES.reactions.equilibrium import classify_equilibrium_constraint
        hp = HenryEquilibrium(H_ref=3.4e-4, dlnH=2400.0)  # gas/liquid species unset
        assert hp.stoichiometry == ()
        with pytest.raises(ValueError):
            classify_equilibrium_constraint(hp)


# ── BisectionChemicalEquilibriumEngine.from_reactions(): cross-phase exposure ─────────────────

class TestChemicalEquilibriumEngineCrossPhaseExposure:
    def test_default_cross_phase_constraints_empty(self):
        from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
        eng = BisectionChemicalEquilibriumEngine()
        assert eng.cross_phase_constraints == ()

    def test_cross_phase_reaction_not_dropped_silently(self):
        from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
        partition_decl = _co2_gas_liquid_declaration()
        eng = BisectionChemicalEquilibriumEngine.from_reactions(
            [_water_reaction(), _co2_acid_reaction(), partition_decl], T_K=298.15,
        )
        assert len(eng.cross_phase_constraints) == 1
        assert eng.cross_phase_constraints[0] is partition_decl

    def test_henry_equilibrium_exposed_as_cross_phase(self):
        from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
        henry = _henry_co2()
        eng = BisectionChemicalEquilibriumEngine.from_reactions(
            [_water_reaction(), _co2_acid_reaction(), henry], T_K=298.15,
        )
        assert henry in eng.cross_phase_constraints

    def test_no_cross_phase_items_gives_empty_tuple(self):
        from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
        eng = BisectionChemicalEquilibriumEngine.from_reactions(
            [_water_reaction(), _co2_acid_reaction()], T_K=298.15,
        )
        assert eng.cross_phase_constraints == ()

    def test_acid_base_solve_unaffected_by_cross_phase_item_presence(self):
        from PyOMES.core.phases import LiquidPhase
        from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
        eng_plain = BisectionChemicalEquilibriumEngine.from_reactions(
            [_water_reaction(), _co2_acid_reaction()], T_K=298.15,
        )
        eng_mixed = BisectionChemicalEquilibriumEngine.from_reactions(
            [_water_reaction(), _co2_acid_reaction(), _co2_gas_liquid_declaration()],
            T_K=298.15,
        )
        liq = LiquidPhase(n_mol={"CO2": 0.01}, V_L=1.0, T_K=298.15)
        out_plain = _species(eng_plain.solve(phases={"liquid": liq}, strong_kwargs={}, T_K=298.15))
        out_mixed = _species(eng_mixed.solve(phases={"liquid": liq}, strong_kwargs={}, T_K=298.15))
        assert out_plain["CO2"] == pytest.approx(out_mixed["CO2"])
        assert out_plain["HCO3-"] == pytest.approx(out_mixed["HCO3-"])

    def test_non_conforming_item_raises_value_error(self):
        from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
        with pytest.raises(ValueError):
            BisectionChemicalEquilibriumEngine.from_reactions([_water_reaction(), object()], T_K=298.15)


# ── NRChemicalEquilibriumEngine.from_reactions(): auto-classified precipitation ──────

class TestNRChemicalEquilibriumEngineAutoPrecipitation:
    def _carbonate_reactions(self):
        from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
        from PyOMES.chemistry.common_species import (
            CO2, H2O, HCO3_minus, CO3_2minus, H_plus,
        )
        co2_second = EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=HCO3_minus, phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=CO3_2minus, phase="liquid", coefficient=+1.0),
                StoichiometryEntry(species=H_plus, phase="liquid", coefficient=+1.0),
            ],
            log_K=-10.33, total_id="CO2", label="co2_second",
        )
        return [_water_reaction(), _co2_acid_reaction(), co2_second]

    def test_auto_classifies_solid_liquid_from_flat_list(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        calcite = _calcite_reaction()
        engine = NRChemicalEquilibriumEngine.from_reactions(
            self._carbonate_reactions() + [calcite],
            use_activity=True, activity_model="davies",
        )
        assert calcite in engine._precipitation_reactions

    def test_solve_matches_deprecated_kwarg_path(self):
        """Same physics whether the mineral arrives via the flat list or the
        deprecated precipitation_reactions= kwarg."""
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        engine_flat = NRChemicalEquilibriumEngine.from_reactions(
            self._carbonate_reactions() + [_calcite_reaction()],
            use_activity=True, activity_model="davies",
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            engine_kwarg = NRChemicalEquilibriumEngine.from_reactions(
                self._carbonate_reactions(),
                precipitation_reactions=[_calcite_reaction()],
                use_activity=True, activity_model="davies",
            )
        kwargs = dict(totals={"CO2": 0.010}, strong_ions={"CT_Ca": 0.002, "CT_Na": 0.005})
        out_flat = engine_flat.solve(**kwargs)
        out_kwarg = engine_kwarg.solve(**kwargs)
        assert out_flat.extra["minerals_xi_mol_L"]["calcite"] == pytest.approx(
            out_kwarg.extra["minerals_xi_mol_L"]["calcite"]
        )

    def test_deprecated_kwarg_emits_warning(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        with pytest.warns(DeprecationWarning):
            NRChemicalEquilibriumEngine.from_reactions(
                self._carbonate_reactions(),
                precipitation_reactions=[_calcite_reaction()],
                use_activity=True, activity_model="davies",
            )

    def test_no_solid_liquid_items_gives_no_precipitation_reactions(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        engine = NRChemicalEquilibriumEngine.from_reactions(self._carbonate_reactions())
        assert engine._precipitation_reactions == []

    def test_deprecated_kwarg_merges_with_auto_detected(self):
        """A solid-liquid item in the flat list AND the deprecated kwarg both
        end up in _precipitation_reactions."""
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        calcite = _calcite_reaction()
        other_mineral = _calcite_reaction()
        other_mineral.label = "other_mineral"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            engine = NRChemicalEquilibriumEngine.from_reactions(
                self._carbonate_reactions() + [calcite],
                precipitation_reactions=[other_mineral],
                use_activity=True, activity_model="davies",
            )
        assert calcite in engine._precipitation_reactions
        assert other_mineral in engine._precipitation_reactions
        assert len(engine._precipitation_reactions) == 2


# ── build_tableau(): still excludes gas-liquid/solid-liquid from the graph ──

class TestBuildTableauSharedClassifierFilter:
    def test_henry_equilibrium_excluded_from_graph(self):
        from PyOMES.chemical_equilibrium.nr_tableau import build_tableau
        tableau = build_tableau([_water_reaction(), _co2_acid_reaction(), _henry_co2()])
        assert "CO2" in tableau.masters

    def test_ksp_equilibrium_excluded_from_graph(self):
        from PyOMES.chemistry import KspEquilibrium
        from PyOMES.chemistry.species import Species
        from PyOMES.reactions.stoichiometry import StoichiometryEntry
        from PyOMES.chemical_equilibrium.nr_tableau import build_tableau
        MineralX_solid = Species(id="MineralX(s)", atoms={"Mn": 1}, charge=0)
        MineralX_aq = Species(id="MineralX", atoms={"Mn": 1}, charge=0)
        ksp = KspEquilibrium(
            stoichiometry=[
                StoichiometryEntry(species=MineralX_solid, phase="solid", coefficient=-1.0),
                StoichiometryEntry(species=MineralX_aq, phase="liquid", coefficient=+1.0),
            ],
            Ksp=1.0e-4,
        )
        tableau = build_tableau([_water_reaction(), _co2_acid_reaction(), ksp])
        assert "MineralX" not in tableau.masters
        assert "MineralX" not in [sec.species_id for sec in tableau.secondaries]

    def test_cross_phase_equilibrium_reaction_excluded_from_graph(self):
        from PyOMES.chemical_equilibrium.nr_tableau import build_tableau
        tableau = build_tableau(
            [_water_reaction(), _co2_acid_reaction(), _co2_gas_liquid_declaration()]
        )
        assert "CO2" in tableau.masters

    def test_non_equilibrium_constraint_item_silently_skipped(self):
        from PyOMES.chemical_equilibrium.nr_tableau import build_tableau
        tableau = build_tableau([_water_reaction(), _co2_acid_reaction(), object()])
        assert "CO2" in tableau.masters
