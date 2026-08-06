# -*- coding: utf-8 -*-
"""Tests for CP3 of LAYER1_GAP_CLOSURE: step_internal_transfer() scope-filter fix.

Before this phase, no scope filter existed at all in
``step_internal_transfer()`` — a species folded into a speciation
engine's own simultaneous gas-liquid solve (CP1/CP2) would *also* have
any separately-declared ``transfer_models`` entry applied, double-moving
it between phases. This file exercises the fix:

* ``ControlVolume.gas_liquid_species()`` /
  ``_gas_liquid_engine_owned_species()`` correctly identify only species
  with a *folded* gas-liquid row — not every species the engine resolves
  (ordinary acid-base-only species must keep using ``transfer_models``
  normally; over-filtering against the full ``algebraic_species()`` set
  would silently break that).
* A folded species carrying a redundant ``KineticTransferModel`` is
  skipped (not double-applied) by ``step_internal_transfer()``.
* A folded species carrying a redundant ``EquilibriumTransferModel`` is
  likewise skipped.
* A non-folded species (ordinary acid-base only, or no acid-base
  chemistry at all) keeps transferring normally — the filter is scoped,
  not a blanket disable.
"""
from __future__ import annotations

import pytest


def _o2_species():
    from PyOMES.chemistry.species import Species
    return Species(id="O2", atoms={"O": 2}, charge=0)


def _n2_species():
    from PyOMES.chemistry.species import Species
    return Species(id="N2", atoms={"N": 2}, charge=0)


def _water_rxn():
    from PyOMES.chemistry.common_species import H2O, H_plus, OH_minus
    from PyOMES.reactions.equilibrium import EquilibriumReaction
    from PyOMES.reactions.stoichiometry import StoichiometryEntry as E
    return EquilibriumReaction(
        stoichiometry=[
            E(species=H2O, phase="liquid", coefficient=-1.0),
            E(species=H_plus, phase="liquid", coefficient=+1.0),
            E(species=OH_minus, phase="liquid", coefficient=+1.0),
        ],
        log_K=-14.0, balance_elements=("H", "O"), label="water",
    )


def _o2_henry():
    from PyOMES.chemistry import HenryEquilibrium
    O2 = _o2_species()
    return HenryEquilibrium(
        H_ref=1.3e-5, dlnH=1500.0, gas_species=O2, liquid_species=O2,
        label="henry_O2",
    )


def _n2_henry():
    from PyOMES.chemistry import HenryEquilibrium
    N2 = _n2_species()
    return HenryEquilibrium(H_ref=6.4e-6, dlnH=1300.0, label="henry_N2_partition_only")


def _make_cv(reactions, transfer_models, *, gas_n_mol, liq_n_mol, solver="newton_raphson"):
    from PyOMES.core.phases import GasPhase, LiquidPhase
    from PyOMES.core.control_volume import ControlVolume
    from PyOMES.reactions.reaction_system import ReactionSystem

    gas = GasPhase(n_mol=dict(gas_n_mol), V_L=0.5, T_K=298.15)
    liq = LiquidPhase(n_mol=dict(liq_n_mol), V_L=1.0, T_K=298.15)
    rs = ReactionSystem(reactions, solver=solver)
    return ControlVolume(
        phases={"gas": gas, "liquid": liq},
        transfer_models=transfer_models,
        reaction_system=rs,
        label="test",
    )


# ═══════════════════════════════════════════════════════════════════════════
#  gas_liquid_species() / _gas_liquid_engine_owned_species() identification
# ═══════════════════════════════════════════════════════════════════════════

class TestOwnedSpeciesIdentification:

    def test_folded_species_reported_owned(self):
        from PyOMES.core.transfer_models import KineticTransferModel
        cv = _make_cv(
            [_water_rxn(), _o2_henry()],
            {"O2": KineticTransferModel(_o2_henry(), k_transfer=100.0)},
            gas_n_mol={"O2": 0.02}, liq_n_mol={"O2": 0.001},
        )
        assert cv._gas_liquid_engine_owned_species() == frozenset({"O2"})

    def test_no_folded_species_when_no_henry_declared(self):
        """Ordinary acid-base-only chemistry (no gas-liquid Henry row)
        must report no owned species — over-filtering here would break
        legitimate transfer_models usage for CO2-style species."""
        from PyOMES.core.transfer_models import KineticTransferModel
        from PyOMES.chemistry.common_species import (
            H2O, H_plus, CO2, HCO3_minus,
        )
        from PyOMES.reactions.equilibrium import EquilibriumReaction
        from PyOMES.reactions.stoichiometry import StoichiometryEntry as E

        co2_first = EquilibriumReaction(
            stoichiometry=[
                E(species=CO2, phase="liquid", coefficient=-1.0),
                E(species=H2O, phase="liquid", coefficient=-1.0),
                E(species=HCO3_minus, phase="liquid", coefficient=+1.0),
                E(species=H_plus, phase="liquid", coefficient=+1.0),
            ],
            log_K=-6.35, total_id="CO2", balance_elements=("C", "H", "O"),
            label="co2_first",
        )
        cv = _make_cv(
            [_water_rxn(), co2_first],
            {"CO2": KineticTransferModel(_o2_henry(), k_transfer=150.0)},
            gas_n_mol={"CO2": 0.01}, liq_n_mol={"CO2": 0.01},
        )
        assert cv._gas_liquid_engine_owned_species() == frozenset()

    def test_no_reaction_system_reports_empty(self):
        from PyOMES.core.phases import GasPhase, LiquidPhase
        from PyOMES.core.control_volume import ControlVolume

        gas = GasPhase(n_mol={"O2": 0.02}, V_L=0.5, T_K=298.15)
        liq = LiquidPhase(n_mol={"O2": 0.001}, V_L=1.0, T_K=298.15)
        cv = ControlVolume(phases={"gas": gas, "liquid": liq}, label="bare")
        assert cv._gas_liquid_engine_owned_species() == frozenset()


# ═══════════════════════════════════════════════════════════════════════════
#  step_internal_transfer() skips folded species, keeps non-folded ones
# ═══════════════════════════════════════════════════════════════════════════

class TestScopeFilterAppliedInTransferStep:

    def test_kinetic_transfer_model_skipped_for_folded_species(self):
        """The plan's own CP3 test requirement: a folded species with a
        redundant KineticTransferModel is not double-applied."""
        from PyOMES.core.transfer_models import KineticTransferModel

        cv = _make_cv(
            [_water_rxn(), _o2_henry()],
            {"O2": KineticTransferModel(_o2_henry(), k_transfer=500.0)},
            gas_n_mol={"O2": 0.02}, liq_n_mol={"O2": 0.001},
        )
        n_gas_before = cv.phases["gas"].n_mol["O2"]
        n_liq_before = cv.phases["liquid"].n_mol["O2"]

        diag = cv.step_internal_transfer(dt_h=0.1)

        assert cv.phases["gas"].n_mol["O2"] == pytest.approx(n_gas_before)
        assert cv.phases["liquid"].n_mol["O2"] == pytest.approx(n_liq_before)
        assert diag.residual.get("O2", 0.0) == pytest.approx(0.0, abs=1e-15)

    def test_equilibrium_transfer_model_skipped_for_folded_species(self):
        from PyOMES.core.transfer_models import EquilibriumTransferModel

        cv = _make_cv(
            [_water_rxn(), _o2_henry()],
            {"O2": EquilibriumTransferModel(_o2_henry())},
            gas_n_mol={"O2": 0.02}, liq_n_mol={"O2": 0.001},
        )
        n_gas_before = cv.phases["gas"].n_mol["O2"]
        n_liq_before = cv.phases["liquid"].n_mol["O2"]

        cv.step_internal_transfer(dt_h=0.1)

        assert cv.phases["gas"].n_mol["O2"] == pytest.approx(n_gas_before)
        assert cv.phases["liquid"].n_mol["O2"] == pytest.approx(n_liq_before)

    def test_without_filter_transfer_would_have_moved_mass(self):
        """Sanity check that the scenario above is a real test — the
        SAME transfer model, run standalone (no engine ownership), does
        move O2 between phases."""
        from PyOMES.core.transfer_models import KineticTransferModel
        from PyOMES.core.phases import GasPhase, LiquidPhase
        from PyOMES.core.control_volume import ControlVolume

        gas = GasPhase(n_mol={"O2": 0.02}, V_L=0.5, T_K=298.15)
        liq = LiquidPhase(n_mol={"O2": 0.001}, V_L=1.0, T_K=298.15)
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            transfer_models={"O2": KineticTransferModel(_o2_henry(), k_transfer=500.0)},
            label="no_engine",
        )
        n_gas_before = cv.phases["gas"].n_mol["O2"]
        cv.step_internal_transfer(dt_h=0.1)
        assert cv.phases["gas"].n_mol["O2"] != pytest.approx(n_gas_before)

    def test_mixed_cv_skips_only_folded_species(self):
        """A CV with one folded species (O2) and one non-folded species
        with its own transfer model (N2, no acid-base coupling declared
        for it at all) skips only O2's transfer, not N2's."""
        from PyOMES.core.transfer_models import KineticTransferModel

        cv = _make_cv(
            [_water_rxn(), _o2_henry()],
            {
                "O2": KineticTransferModel(_o2_henry(), k_transfer=500.0),
                "N2": KineticTransferModel(_n2_henry(), k_transfer=500.0),
            },
            gas_n_mol={"O2": 0.02, "N2": 0.03},
            liq_n_mol={"O2": 0.001, "N2": 0.001},
        )
        n_gas_O2_before = cv.phases["gas"].n_mol["O2"]
        n_gas_N2_before = cv.phases["gas"].n_mol["N2"]

        cv.step_internal_transfer(dt_h=0.1)

        assert cv.phases["gas"].n_mol["O2"] == pytest.approx(n_gas_O2_before)
        assert cv.phases["gas"].n_mol["N2"] != pytest.approx(n_gas_N2_before)
