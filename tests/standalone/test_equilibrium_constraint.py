# -*- coding: utf-8 -*-
"""Tests for EquilibriumConstraint protocol and the *Equilibrium sibling types (CP1).

Covers:
- EquilibriumReaction conformance to EquilibriumConstraint + vant_hoff_log_K
  correctness against manual van't Hoff math.
- HenryEquilibrium / RaoultEquilibrium / KspEquilibrium: isinstance of both
  PartitionModel and EquilibriumConstraint; log_K / dH_J_per_mol manual
  conversion checks.
- HenryPartition / RaoultPartition deprecated-alias behavior: DeprecationWarning
  + construct-identical-object.
- KspEquilibrium single-ion solubility-cap behavior and multi-ion
  NotImplementedError behavior.
"""
from __future__ import annotations

import math
import warnings

import pytest

from PyOMES.units import R_J_PER_MOL_K as _R_J_MOL_K
_T_REF = 298.15


# ── EquilibriumReaction conformance ──────────────────────────────────────────

class TestEquilibriumReactionConformance:
    def _rxn(self, **kwargs):
        from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
        from PyOMES.chemistry.common_species import H2O, H_plus, OH_minus
        return EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=H2O,      phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=H_plus,    phase="liquid", coefficient=+1.0),
                StoichiometryEntry(species=OH_minus,  phase="liquid", coefficient=+1.0),
            ],
            log_K=-14.0,
            label="water",
            **kwargs,
        )

    def test_isinstance_equilibrium_constraint(self):
        from PyOMES.reactions.equilibrium import EquilibriumConstraint
        rxn = self._rxn()
        assert isinstance(rxn, EquilibriumConstraint)

    def test_attributes_are_plain_values_not_methods(self):
        rxn = self._rxn()
        assert isinstance(rxn.log_K, float)
        assert rxn.log_K == pytest.approx(-14.0)
        assert rxn.T_ref_K == pytest.approx(298.15)
        assert rxn.dH_J_per_mol is None

    def test_vant_hoff_no_dH_returns_log_K_unchanged(self):
        from PyOMES.reactions.equilibrium import vant_hoff_log_K
        rxn = self._rxn()
        assert vant_hoff_log_K(rxn, 350.0) == pytest.approx(rxn.log_K)

    def test_vant_hoff_at_reference_temperature_returns_log_K(self):
        from PyOMES.reactions.equilibrium import vant_hoff_log_K
        rxn = self._rxn(dH_J_per_mol=55800.0)
        assert vant_hoff_log_K(rxn, rxn.T_ref_K) == pytest.approx(rxn.log_K)

    def test_vant_hoff_matches_manual_calculation(self):
        from PyOMES.reactions.equilibrium import vant_hoff_log_K
        dH = 55800.0  # J/mol, water autoionization (endothermic)
        rxn = self._rxn(dH_J_per_mol=dH)
        T_K = 323.15
        delta_ln_K = -(dH / _R_J_MOL_K) * (1.0 / T_K - 1.0 / _T_REF)
        expected = rxn.log_K + delta_ln_K / math.log(10.0)
        assert vant_hoff_log_K(rxn, T_K) == pytest.approx(expected, rel=1e-12)

    def test_vant_hoff_endothermic_increases_log_K_with_temperature(self):
        from PyOMES.reactions.equilibrium import vant_hoff_log_K
        rxn = self._rxn(dH_J_per_mol=55800.0)
        assert vant_hoff_log_K(rxn, 323.15) > vant_hoff_log_K(rxn, _T_REF)


# ── HenryEquilibrium ──────────────────────────────────────────────────────────

_H2S_H_REF = 0.10 * 1000.0 / 101325.0  # kH = 0.10 mol/(L·atm) at T_ref
_H2S_DLN_H = 2100.0


class TestHenryEquilibriumConformance:
    def _henry(self, **kwargs):
        from PyOMES.chemistry import HenryEquilibrium
        return HenryEquilibrium(
            H_ref=_H2S_H_REF, dlnH=_H2S_DLN_H,
            gas_species="H2S", liquid_species="H2S",
            **kwargs,
        )

    def test_isinstance_partition_model(self):
        from PyOMES.chemistry import PartitionModel
        assert isinstance(self._henry(), PartitionModel)

    def test_isinstance_equilibrium_constraint(self):
        from PyOMES.reactions.equilibrium import EquilibriumConstraint
        assert isinstance(self._henry(), EquilibriumConstraint)

    def test_T_ref_K_matches_T_ref(self):
        hp = self._henry()
        assert hp.T_ref_K == pytest.approx(hp.T_ref)

    def test_log_K_matches_manual_conversion(self):
        hp = self._henry()
        kH_mol_L_atm = (hp.H_ref / 1000.0) * 101325.0
        assert hp.log_K == pytest.approx(math.log10(kH_mol_L_atm), rel=1e-12)

    def test_log_K_matches_kH_at_T_ref(self):
        hp = self._henry()
        assert hp.log_K == pytest.approx(math.log10(hp._kH_mol_L_atm(hp.T_ref)), rel=1e-12)

    def test_dH_J_per_mol_matches_manual_conversion(self):
        hp = self._henry()
        assert hp.dH_J_per_mol == pytest.approx(-_H2S_DLN_H * _R_J_MOL_K, rel=1e-9)

    def test_stoichiometry_when_species_set(self):
        hp = self._henry()
        entries = hp.stoichiometry
        assert len(entries) == 2
        phases = {e.phase for e in entries}
        assert phases == {"gas", "liquid"}
        gas_entry = next(e for e in entries if e.phase == "gas")
        liq_entry = next(e for e in entries if e.phase == "liquid")
        assert gas_entry.coefficient == pytest.approx(-1.0)
        assert liq_entry.coefficient == pytest.approx(+1.0)
        assert gas_entry.species.id == "H2S"
        assert liq_entry.species.id == "H2S"

    def test_stoichiometry_empty_when_species_unset(self):
        from PyOMES.chemistry import HenryEquilibrium
        hp = HenryEquilibrium(H_ref=_H2S_H_REF, dlnH=0.0)
        assert hp.stoichiometry == ()

    def test_vant_hoff_matches_kH_temperature_dependence(self):
        """vant_hoff_log_K(hp, T) must agree with the native _kH_mol_L_atm(T) path."""
        from PyOMES.reactions.equilibrium import vant_hoff_log_K
        hp = self._henry()
        T_K = 315.0
        expected = math.log10(hp._kH_mol_L_atm(T_K))
        assert vant_hoff_log_K(hp, T_K) == pytest.approx(expected, rel=1e-9)


# ── RaoultEquilibrium ─────────────────────────────────────────────────────────

class TestRaoultEquilibriumConformance:
    def test_isinstance_partition_model(self):
        from PyOMES.chemistry import RaoultEquilibrium, PartitionModel
        assert isinstance(RaoultEquilibrium(), PartitionModel)

    def test_isinstance_equilibrium_constraint(self):
        from PyOMES.chemistry import RaoultEquilibrium
        from PyOMES.reactions.equilibrium import EquilibriumConstraint
        assert isinstance(RaoultEquilibrium(), EquilibriumConstraint)

    def test_default_species_are_water(self):
        from PyOMES.chemistry import RaoultEquilibrium
        rp = RaoultEquilibrium()
        assert rp.gas_species == "H2O"
        assert rp.liquid_species == "H2O"

    def test_T_ref_K_matches_T_ref(self):
        from PyOMES.chemistry import RaoultEquilibrium
        rp = RaoultEquilibrium()
        assert rp.T_ref_K == pytest.approx(rp.T_ref)

    def test_log_K_matches_manual_conversion(self):
        from PyOMES.chemistry import RaoultEquilibrium
        rp = RaoultEquilibrium()
        assert rp.log_K == pytest.approx(-math.log10(rp.P_sat_ref), rel=1e-12)

    def test_dH_J_per_mol_matches_manual_conversion(self):
        from PyOMES.chemistry import RaoultEquilibrium
        rp = RaoultEquilibrium()
        assert rp.dH_J_per_mol == pytest.approx(-rp.dH_vap, rel=1e-12)

    def test_stoichiometry_gas_liquid_water(self):
        from PyOMES.chemistry import RaoultEquilibrium
        rp = RaoultEquilibrium()
        entries = rp.stoichiometry
        assert len(entries) == 2
        gas_entry = next(e for e in entries if e.phase == "gas")
        liq_entry = next(e for e in entries if e.phase == "liquid")
        assert gas_entry.species.id == "H2O"
        assert liq_entry.species.id == "H2O"
        assert gas_entry.coefficient == pytest.approx(-1.0)
        assert liq_entry.coefficient == pytest.approx(+1.0)

    def test_stoichiometry_empty_when_species_unset(self):
        from PyOMES.chemistry import RaoultEquilibrium
        rp = RaoultEquilibrium(gas_species=None, liquid_species=None)
        assert rp.stoichiometry == ()


# ── KspEquilibrium ────────────────────────────────────────────────────────────

class TestKspEquilibriumSingleIon:
    def _ksp(self, **kwargs):
        from PyOMES.chemistry import KspEquilibrium
        from PyOMES.chemistry.species import Species
        MineralX_solid = Species(id="MineralX(s)", atoms={"Mn": 1, "O": 1}, charge=0)
        MineralX_aq = Species(id="MineralX", atoms={"Mn": 1, "O": 1}, charge=0)
        from PyOMES.reactions.stoichiometry import StoichiometryEntry
        return KspEquilibrium(
            stoichiometry=[
                StoichiometryEntry(species=MineralX_solid, phase="solid", coefficient=-1.0),
                StoichiometryEntry(species=MineralX_aq, phase="liquid", coefficient=+1.0),
            ],
            Ksp=1.0e-4,
            label="MineralX solubility",
            **kwargs,
        )

    def test_isinstance_partition_model(self):
        from PyOMES.chemistry import PartitionModel
        assert isinstance(self._ksp(), PartitionModel)

    def test_isinstance_equilibrium_constraint(self):
        from PyOMES.reactions.equilibrium import EquilibriumConstraint
        assert isinstance(self._ksp(), EquilibriumConstraint)

    def test_log_K_equals_log10_Ksp(self):
        ksp = self._ksp()
        assert ksp.log_K == pytest.approx(math.log10(1.0e-4), rel=1e-12)

    def test_partition_ratio_returns_none(self):
        ksp = self._ksp()
        assert ksp.partition_ratio(1.0, 1.0, _T_REF) is None

    def test_equilibrium_a_moles_below_cap_conserves_total(self):
        ksp = self._ksp()
        # capacity_a (liquid volume) large enough that the solubility cap
        # exceeds n_total — all of it stays dissolved.
        n_liq = ksp.equilibrium_a_moles(n_total=1e-6, capacity_a=10.0, capacity_b=1.0, T_K=_T_REF)
        assert n_liq == pytest.approx(1e-6)

    def test_equilibrium_a_moles_above_cap_is_capped(self):
        ksp = self._ksp()
        # Ksp(T_ref) * capacity_a = 1e-4 * 1.0 = 1e-4 mol cap; n_total exceeds it.
        n_liq = ksp.equilibrium_a_moles(n_total=1.0, capacity_a=1.0, capacity_b=1.0, T_K=_T_REF)
        assert n_liq == pytest.approx(1.0e-4, rel=1e-9)

    def test_equilibrium_a_moles_matches_formula(self):
        ksp = self._ksp(dH_J_per_mol=12000.0)
        from PyOMES.reactions.equilibrium import vant_hoff_log_K
        T_K = 310.0
        Ksp_T = 10.0 ** vant_hoff_log_K(ksp, T_K)
        n_liq = ksp.equilibrium_a_moles(n_total=1.0, capacity_a=2.0, capacity_b=1.0, T_K=T_K)
        assert n_liq == pytest.approx(min(1.0, Ksp_T * 2.0), rel=1e-9)


class TestKspEquilibriumMultiIon:
    def _ksp(self):
        from PyOMES.chemistry import KspEquilibrium
        from PyOMES.chemistry.common_species import Ca_plus_plus, CO3_2minus
        from PyOMES.chemistry.species import Species
        from PyOMES.reactions.stoichiometry import StoichiometryEntry
        CaCO3_solid = Species(id="CaCO3(s)", atoms={"Ca": 1, "C": 1, "O": 3}, charge=0)
        return KspEquilibrium(
            stoichiometry=[
                StoichiometryEntry(species=CaCO3_solid, phase="solid", coefficient=-1.0),
                StoichiometryEntry(species=Ca_plus_plus, phase="liquid", coefficient=+1.0),
                StoichiometryEntry(species=CO3_2minus, phase="liquid", coefficient=+1.0),
            ],
            Ksp=3.3e-9,
            label="calcite solubility",
        )

    def test_isinstance_equilibrium_constraint(self):
        """Multi-ion Ksp still structurally satisfies EquilibriumConstraint."""
        from PyOMES.reactions.equilibrium import EquilibriumConstraint
        assert isinstance(self._ksp(), EquilibriumConstraint)

    def test_isinstance_partition_model(self):
        """isinstance only checks method presence, not single-ion validity."""
        from PyOMES.chemistry import PartitionModel
        assert isinstance(self._ksp(), PartitionModel)

    def test_partition_ratio_raises_not_implemented(self):
        ksp = self._ksp()
        with pytest.raises(NotImplementedError):
            ksp.partition_ratio(1.0, 1.0, _T_REF)

    def test_equilibrium_a_moles_raises_not_implemented(self):
        ksp = self._ksp()
        with pytest.raises(NotImplementedError):
            ksp.equilibrium_a_moles(1.0, 1.0, 1.0, _T_REF)


# ── Deprecated aliases ────────────────────────────────────────────────────────

class TestHenryPartitionAlias:
    def test_emits_deprecation_warning(self):
        from PyOMES.chemistry import HenryPartition
        with pytest.warns(DeprecationWarning):
            HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0)

    def test_constructs_henry_equilibrium(self):
        from PyOMES.chemistry import HenryPartition, HenryEquilibrium
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            hp = HenryPartition(H_ref=_H2S_H_REF, dlnH=0.0)
        assert type(hp) is HenryEquilibrium

    def test_identical_to_direct_construction(self):
        from PyOMES.chemistry import HenryPartition, HenryEquilibrium
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            via_alias = HenryPartition(H_ref=_H2S_H_REF, dlnH=_H2S_DLN_H)
        direct = HenryEquilibrium(H_ref=_H2S_H_REF, dlnH=_H2S_DLN_H)
        assert via_alias == direct


class TestRaoultPartitionAlias:
    def test_emits_deprecation_warning(self):
        from PyOMES.chemistry import RaoultPartition
        with pytest.warns(DeprecationWarning):
            RaoultPartition()

    def test_constructs_raoult_equilibrium(self):
        from PyOMES.chemistry import RaoultPartition, RaoultEquilibrium
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            rp = RaoultPartition()
        assert type(rp) is RaoultEquilibrium

    def test_identical_to_direct_construction(self):
        from PyOMES.chemistry import RaoultPartition, RaoultEquilibrium
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            via_alias = RaoultPartition(dH_vap=45000.0)
        direct = RaoultEquilibrium(dH_vap=45000.0)
        assert via_alias == direct
