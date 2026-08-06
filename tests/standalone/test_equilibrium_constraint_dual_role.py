# -*- coding: utf-8 -*-
"""Tests for CP3: fixing the Henry double-declaration bug (one object, two roles).

Before this phase, a gas-liquid species needed two independently-
parameterized declarations: a ``HenryPartition`` in
``KineticGasLiquidLink``'s ``partition_models=`` dict (feeding the
kinetic/equilibrium transfer math), and a separate, hand-built
cross-phase ``EquilibriumReaction`` in the CV's ``ReactionSystem``
(feeding ``derive_speciation_keys``/the speciation engine). Nothing
validated the two agreed.

``HenryEquilibrium`` satisfies both ``PartitionModel`` and
``EquilibriumConstraint`` from the same ``H_ref``/``dlnH`` fields — so
constructing *one* instance and using it in both places makes
disagreement structurally impossible, rather than merely catching it
at validation time. These tests exercise that pattern end-to-end
through ``ReactionSystem`` + ``KineticGasLiquidLink``, and demonstrate
the class of bug the old two-object pattern permitted.
"""
from __future__ import annotations

import math

import pytest


def _co2_henry(H_ref=3.4e-4, dlnH=2400.0):
    from PyOMES.chemistry import HenryEquilibrium
    return HenryEquilibrium(
        H_ref=H_ref, dlnH=dlnH, gas_species="CO2", liquid_species="CO2",
        label="partition_CO2",
    )


class TestOneObjectTwoRoles:
    def test_same_instance_satisfies_both_protocols(self):
        from PyOMES.chemistry import PartitionModel
        from PyOMES.reactions.equilibrium import EquilibriumConstraint
        co2 = _co2_henry()
        assert isinstance(co2, PartitionModel)
        assert isinstance(co2, EquilibriumConstraint)

    def test_reaction_system_accepts_it_directly(self):
        """ReactionSystem no longer hard-requires EquilibriumReaction —
        a HenryEquilibrium is classified and bucketed like any other
        EquilibriumConstraint (CP3's fix to the CP1/CP2 blocker)."""
        from PyOMES.reactions import ReactionSystem
        co2 = _co2_henry()
        system = ReactionSystem([co2])
        assert co2 in system.cross_phase_equilibria
        assert co2 not in system.single_phase_equilibria
        assert co2 not in system.precipitation_equilibria

    def test_partition_models_and_reaction_list_share_one_instance(self):
        """The exact object passed to KineticGasLiquidLink's
        partition_models= is the same object in the CV's reaction list —
        no separate declaration."""
        from PyOMES.reactions import ReactionSystem
        from PyOMES.core.gas_liquid_link import KineticGasLiquidLink

        co2 = _co2_henry()
        system = ReactionSystem([co2])
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"CO2": co2},
            speciation_keys={},
        )
        link.derive_speciation_keys(system)

        assert link.partition_models["CO2"] is co2
        assert system.cross_phase_equilibria[0] is co2
        # Both roles trace back to the identical object — by construction,
        # not by cross-validation.
        assert link.partition_models["CO2"] is system.cross_phase_equilibria[0]

    def test_derive_speciation_keys_populates_from_henry_equilibrium(self):
        """derive_speciation_keys picks up the gas<->liquid routing from
        the HenryEquilibrium's own stoichiometry, exactly as it would from
        a hand-built cross-phase EquilibriumReaction."""
        from PyOMES.reactions import ReactionSystem
        from PyOMES.core.gas_liquid_link import KineticGasLiquidLink

        co2 = _co2_henry()
        link = KineticGasLiquidLink(
            gas_cv_key="gas", gas_phase_key="gas",
            liquid_cv_key="liquid", liquid_phase_key="liquid",
            partition_models={"CO2": co2},
            speciation_keys={},
        )
        link.derive_speciation_keys(ReactionSystem([co2]))
        assert link.speciation_keys == {"CO2": "CO2"}

    def test_cv_construction_end_to_end(self):
        """Full ControlVolume construction with the shared instance: the
        internal KineticGasLiquidLink's speciation_keys populate correctly
        via the CV's __init__ hook."""
        from PyOMES.core.phases import GasPhase, LiquidPhase
        from PyOMES.core.control_volume import ControlVolume
        from PyOMES.core.gas_liquid_link import KineticGasLiquidLink
        from PyOMES.reactions import ReactionSystem

        co2 = _co2_henry()
        link = KineticGasLiquidLink(
            gas_cv_key="cv", gas_phase_key="gas",
            liquid_cv_key="cv", liquid_phase_key="liquid",
            partition_models={"CO2": co2}, kLa={"CO2": 100.0},
            speciation_keys={},
        )
        gas = GasPhase(n_mol={"CO2": 0.001}, V_L=0.4, T_K=305.15)
        liq = LiquidPhase(n_mol={"CO2": 0.01}, V_L=1.6, T_K=305.15)
        cv = ControlVolume(
            phases={"gas": gas, "liquid": liq},
            internal_interfaces=[link],
            reaction_system=ReactionSystem([co2]),
            label="test",
        )
        iface = cv.internal_interfaces[0]
        assert iface.speciation_keys == {"CO2": "CO2"}
        assert iface.partition_models["CO2"] is co2


class TestOldBugStructurallyImpossible:
    """Regression coverage for the double-declaration bug the old pattern
    permitted: two independently-parameterized Henry declarations for the
    same species that nothing checked agreed with each other."""

    def test_old_pattern_could_silently_disagree(self):
        """Demonstrates the bug this phase fixes: a partition_models dict
        entry and a separately-built cross-phase EquilibriumReaction for
        the same species, with nothing to relate their numbers. This is
        what a caller could do before HenryEquilibrium existed (and can
        still do today with two independent HenryEquilibrium instances,
        or a HenryEquilibrium alongside a hand-built EquilibriumReaction) —
        it is not prevented, only made unnecessary."""
        from PyOMES.chemistry import HenryEquilibrium
        from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
        from PyOMES.chemistry.common_species import CO2

        kinetic_side = HenryEquilibrium(H_ref=3.4e-4, dlnH=2400.0)
        # A hand-built "routing" reaction with no relation whatsoever to
        # kinetic_side's H_ref/dlnH — could be anything, including
        # (as here) a nonsensical disagreement, and nothing raises.
        declared_side = EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=CO2, phase="gas", coefficient=-1.0),
                StoichiometryEntry(species=CO2, phase="liquid", coefficient=+1.0),
            ],
            log_K=999.0,  # wildly inconsistent with kinetic_side.log_K; unchecked
            balance_elements=("C", "O"),
        )
        assert declared_side.log_K != pytest.approx(kinetic_side.log_K)
        # Nothing in either object's construction, or in ReactionSystem,
        # cross-validates these two — that's the bug. The fix isn't a
        # check; it's making the two-object pattern unnecessary (below).

    def test_new_pattern_cannot_disagree_with_itself(self):
        """With one HenryEquilibrium instance used in both roles, the
        PartitionModel-role kH(T) and the EquilibriumConstraint-role
        log_K/dH_J_per_mol are both always read from the same H_ref/dlnH
        fields — there is no second set of numbers to fall out of sync."""
        from PyOMES.reactions.equilibrium import vant_hoff_log_K

        co2 = _co2_henry(H_ref=3.4e-4, dlnH=2400.0)
        for T_K in (280.0, 298.15, 320.0, 350.0):
            kH_direct = co2._kH_mol_L_atm(T_K)  # PartitionModel-role computation
            log_K_via_constraint = vant_hoff_log_K(co2, T_K)  # EquilibriumConstraint role
            assert log_K_via_constraint == pytest.approx(
                math.log10(kH_direct), rel=1e-9
            )

    def test_changing_the_one_instance_changes_both_roles_together(self):
        """There is exactly one place to edit — H_ref/dlnH on the shared
        instance — and both roles see the update immediately, since both
        read from the same object rather than independent copies."""
        from PyOMES.chemistry import HenryEquilibrium
        import dataclasses

        original = _co2_henry(H_ref=3.4e-4, dlnH=2400.0)
        updated = dataclasses.replace(original, H_ref=5.0e-4)
        # partition_ratio-driving field and log_K both shift together
        assert updated.H_ref != original.H_ref
        assert updated.log_K != pytest.approx(original.log_K)
        assert updated.log_K == pytest.approx(
            math.log10(updated._kH_mol_L_atm(updated.T_ref)), rel=1e-9
        )


class TestAnaerobicDigestionDatabaseFixed:
    """AD_BASIC previously declared CO2's partition twice (a HenryPartition
    in partition_models, a separate log_K-less placeholder EquilibriumReaction
    in the reaction list). CP3 collapses these into one shared instance."""

    def test_co2_partition_model_and_reaction_are_the_same_object(self):
        from PyOMES.chemistry.databases.anaerobic_digestion import (
            AD_BASIC, _CO2_HENRY,
        )
        assert AD_BASIC.partition_models["CO2"] is _CO2_HENRY
        assert _CO2_HENRY in AD_BASIC.reactions.reactions

    def test_co2_classified_as_gas_liquid(self):
        from PyOMES.chemistry.databases.anaerobic_digestion import _CO2_HENRY
        from PyOMES.reactions.equilibrium import classify_equilibrium_constraint
        assert classify_equilibrium_constraint(_CO2_HENRY) == "gas_liquid"
