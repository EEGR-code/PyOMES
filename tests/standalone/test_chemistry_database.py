# -*- coding: utf-8 -*-
"""Tests for ChemistryDatabase and stock database modules (C9)."""
from __future__ import annotations

import pytest


class TestChemistryDatabaseConstruction:
    def test_minimal_construction(self):
        from PyOMES.chemistry import ChemistryDatabase
        from PyOMES.thermo import ThermoFramework
        db = ChemistryDatabase(thermo=ThermoFramework())
        assert db.thermo is not None
        assert db.species == {}
        assert db.reactions is None

    def test_with_species(self):
        from PyOMES.chemistry import ChemistryDatabase, Species
        from PyOMES.thermo import ThermoFramework
        sp = Species(id="A", atoms={"C": 1}, charge=0, MW=12.0)
        db = ChemistryDatabase(
            thermo=ThermoFramework(),
            species={"A": sp},
        )
        assert "A" in db.species
        assert db.species["A"] is sp

    def test_frozen(self):
        import dataclasses
        from PyOMES.chemistry import ChemistryDatabase
        from PyOMES.thermo import ThermoFramework
        db = ChemistryDatabase(thermo=ThermoFramework())
        with pytest.raises((dataclasses.FrozenInstanceError, AttributeError)):
            db.thermo = ThermoFramework()  # type: ignore


class TestChemistryDatabaseExtend:
    def _base(self):
        from PyOMES.chemistry import ChemistryDatabase, Species
        from PyOMES.thermo import ThermoFramework
        sp = Species(id="A", atoms={"C": 1}, charge=0, MW=12.0)
        return ChemistryDatabase(
            thermo=ThermoFramework(),
            species={"A": sp},
        )

    def test_extend_adds_species(self):
        from PyOMES.chemistry import Species
        from PyOMES.thermo import ThermoFramework
        base = self._base()
        sp_b = Species(id="B", atoms={"N": 1}, charge=0, MW=14.0)
        extended = base.extend(species={"B": sp_b})
        assert "A" in extended.species
        assert "B" in extended.species
        assert "B" not in base.species  # original unchanged

    def test_extend_overrides_thermo(self):
        from PyOMES.thermo import ThermoFramework
        base = self._base()
        from PyOMES.thermo import DaviesLiquidModel
        new_thermo = ThermoFramework(liquid_activity=DaviesLiquidModel())
        extended = base.extend(thermo=new_thermo)
        assert extended.thermo.use_activity is True
        assert base.thermo.use_activity is False

    def test_extend_without_args_copies(self):
        base = self._base()
        extended = base.extend()
        assert extended.species == base.species
        assert extended.thermo == base.thermo

    def test_extend_adds_reactions(self):
        from PyOMES.chemistry.common_species import CO2, H_plus, H2O, HCO3_minus
        from PyOMES.reactions.equilibrium import EquilibriumReaction
        from PyOMES.reactions.stoichiometry import StoichiometryEntry
        base = self._base()
        rxn = EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=CO2,        phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=H2O,        phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=HCO3_minus, phase="liquid", coefficient=+1.0),
                StoichiometryEntry(species=H_plus,     phase="liquid", coefficient=+1.0),
            ],
            log_K=-6.35, total_id="CO2",
            balance_elements=("C", "H", "O"),
        )
        extended = base.extend(reactions=[rxn])
        assert extended.reactions is not None
        rxn_ids = [r.label for r in extended.reactions]
        rxns = list(extended.reactions)
        assert len(rxns) == 1


class TestStockDatabases:
    def test_aqueous_default_imports(self):
        from PyOMES.chemistry.databases.aqueous import AQUEOUS_DEFAULT
        assert AQUEOUS_DEFAULT is not None
        assert "CO2" in AQUEOUS_DEFAULT.species
        assert "HCO3-" in AQUEOUS_DEFAULT.species
        assert "NH4+" in AQUEOUS_DEFAULT.species

    def test_aqueous_has_reactions(self):
        from PyOMES.chemistry.databases.aqueous import AQUEOUS_DEFAULT
        assert AQUEOUS_DEFAULT.reactions is not None
        labels = [r.label for r in AQUEOUS_DEFAULT.reactions]
        assert "eq_water" in labels
        assert "eq_CO2" in labels
        assert "eq_NH4" in labels

    def test_bioprocess_basic_imports(self):
        from PyOMES.chemistry.databases.bioprocess_basic import BIOPROCESS_BASIC
        assert "H3PO4" in BIOPROCESS_BASIC.species
        assert "SO4--" in BIOPROCESS_BASIC.species
        # still has aqueous species
        assert "CO2" in BIOPROCESS_BASIC.species

    def test_bioprocess_basic_has_phosphate_reactions(self):
        from PyOMES.chemistry.databases.bioprocess_basic import BIOPROCESS_BASIC
        labels = {r.label for r in BIOPROCESS_BASIC.reactions}
        assert "eq_phosphate_1" in labels
        assert "eq_bisulfate" in labels

    def test_ad_basic_imports(self):
        from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC
        assert AD_BASIC is not None
        labels = {r.label for r in AD_BASIC.reactions}
        assert "partition_CO2" in labels

    def test_composition_chain(self):
        """AQUEOUS ⊂ BIOPROCESS_BASIC ⊂ AD_BASIC."""
        from PyOMES.chemistry.databases.aqueous import AQUEOUS_DEFAULT
        from PyOMES.chemistry.databases.bioprocess_basic import BIOPROCESS_BASIC
        from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC
        aqueous_rxn_labels = {r.label for r in AQUEOUS_DEFAULT.reactions}
        bp_rxn_labels = {r.label for r in BIOPROCESS_BASIC.reactions}
        ad_rxn_labels = {r.label for r in AD_BASIC.reactions}
        assert aqueous_rxn_labels.issubset(bp_rxn_labels)
        assert bp_rxn_labels.issubset(ad_rxn_labels)

    def test_thermo_inherited(self):
        from PyOMES.chemistry.databases.aqueous import AQUEOUS_DEFAULT
        from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC
        assert AD_BASIC.thermo == AQUEOUS_DEFAULT.thermo

    def test_ad_basic_h2s_equilibrium(self):
        from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC
        from PyOMES.chemistry.common_species import H2S, HS_minus
        from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
        from PyOMES.reactions.equilibrium import classify_equilibrium_constraint

        single_phase = [
            r for r in AD_BASIC.reactions
            if classify_equilibrium_constraint(r) == "acid_base"
        ]
        h2s_rxns = [
            r for r in single_phase
            if any(e.species.id == "H2S" for e in r.stoichiometry)
        ]
        assert len(h2s_rxns) == 1, "expected exactly one H2S equilibrium in AD_BASIC"

        engine = BisectionChemicalEquilibriumEngine.from_reactions(AD_BASIC.reactions)
        eq_def = engine._equilibrium_set.get("H2S")
        assert eq_def.species_refs == (H2S, HS_minus)

    def test_ad_basic_hs_minus_in_species(self):
        from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC
        assert "HS-" in AD_BASIC.species


class TestChemistryDatabasePartitionModels:
    def _base(self):
        from PyOMES.chemistry import ChemistryDatabase
        from PyOMES.thermo import ThermoFramework
        return ChemistryDatabase(thermo=ThermoFramework())

    def _hp(self, h_ref=1e-4, dln_h=0.0):
        from PyOMES.chemistry import HenryPartition
        return HenryPartition(H_ref=h_ref, dlnH=dln_h)

    def test_default_partition_models_empty(self):
        db = self._base()
        assert db.partition_models == {}

    def test_extend_adds_partition_model(self):
        db = self._base()
        hp = self._hp()
        extended = db.extend(partition_models={"CO2": hp})
        assert "CO2" in extended.partition_models
        assert extended.partition_models["CO2"] is hp

    def test_extend_leaves_parent_unchanged(self):
        db = self._base()
        extended = db.extend(partition_models={"CO2": self._hp()})
        assert db.partition_models == {}

    def test_extend_merges_parent_and_child(self):
        from PyOMES.chemistry import ChemistryDatabase
        from PyOMES.thermo import ThermoFramework
        hp_co2 = self._hp(h_ref=3.3e-4)
        hp_h2s = self._hp(h_ref=9.9e-4)
        parent = ChemistryDatabase(
            thermo=ThermoFramework(),
            partition_models={"CO2": hp_co2},
        )
        child = parent.extend(partition_models={"H2S": hp_h2s})
        assert "CO2" in child.partition_models
        assert "H2S" in child.partition_models

    def test_child_overrides_parent_on_collision(self):
        from PyOMES.chemistry import ChemistryDatabase
        from PyOMES.thermo import ThermoFramework
        hp_parent = self._hp(h_ref=1e-4)
        hp_child  = self._hp(h_ref=9e-4)
        parent = ChemistryDatabase(
            thermo=ThermoFramework(),
            partition_models={"CO2": hp_parent},
        )
        child = parent.extend(partition_models={"CO2": hp_child})
        assert child.partition_models["CO2"] is hp_child
        assert parent.partition_models["CO2"] is hp_parent

    def test_extend_without_partition_models_copies(self):
        from PyOMES.chemistry import ChemistryDatabase
        from PyOMES.thermo import ThermoFramework
        hp = self._hp()
        db = ChemistryDatabase(thermo=ThermoFramework(), partition_models={"CO2": hp})
        extended = db.extend()
        assert extended.partition_models == db.partition_models


class TestStockDatabasePartitionModels:
    def test_bioprocess_basic_has_o2_n2(self):
        from PyOMES.chemistry.databases.bioprocess_basic import BIOPROCESS_BASIC
        assert "O2" in BIOPROCESS_BASIC.partition_models
        assert "N2" in BIOPROCESS_BASIC.partition_models

    def test_ad_basic_has_co2_ch4_h2_nh3_h2s(self):
        from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC
        for key in ("CO2", "CH4", "H2", "NH3", "H2S"):
            assert key in AD_BASIC.partition_models, f"missing {key!r} in AD_BASIC.partition_models"

    def test_ad_basic_inherits_o2_n2_from_bioprocess(self):
        from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC
        assert "O2" in AD_BASIC.partition_models
        assert "N2" in AD_BASIC.partition_models

    def test_bioprocess_basic_does_not_have_ad_species(self):
        from PyOMES.chemistry.databases.bioprocess_basic import BIOPROCESS_BASIC
        for key in ("CO2", "CH4", "H2", "NH3", "H2S"):
            assert key not in BIOPROCESS_BASIC.partition_models

    def test_co2_kH_at_298_matches_sander(self):
        from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC
        hp = AD_BASIC.partition_models["CO2"]
        kH = hp._kH_mol_L_atm(298.15)
        # Sander 2015: kH(CO2) ≈ 3.4e-2 mol/(L·atm) at 25°C
        assert kH == pytest.approx(0.034, rel=0.05)

    def test_h2s_kH_at_298_matches_sander(self):
        from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC
        hp = AD_BASIC.partition_models["H2S"]
        kH = hp._kH_mol_L_atm(298.15)
        # Sander 2015: kH(H2S) ≈ 0.10 mol/(L·atm) at 25°C
        assert kH == pytest.approx(0.10, rel=0.05)

    def test_o2_kH_at_298_matches_sander(self):
        from PyOMES.chemistry.databases.bioprocess_basic import BIOPROCESS_BASIC
        hp = BIOPROCESS_BASIC.partition_models["O2"]
        kH = hp._kH_mol_L_atm(298.15)
        # Sander 2015: kH(O2) ≈ 1.3e-3 mol/(L·atm) at 25°C
        assert kH == pytest.approx(1.3e-3, rel=0.05)


class TestThermoMismatchWarning:
    def _make_cv(self, thermo=None):
        from PyOMES.core.control_volume import ControlVolume
        from PyOMES.core.phases import LiquidPhase
        from PyOMES.chemistry import ChemistryDatabase
        from PyOMES.thermo import ThermoFramework
        liq = LiquidPhase(n_mol={}, V_L=1.0, T_K=298.15)
        db = ChemistryDatabase(thermo=thermo or ThermoFramework()) if thermo is not None \
            else ChemistryDatabase(thermo=ThermoFramework())
        return ControlVolume(phases={"liquid": liq}, chemistry_db=db)

    def test_no_warning_same_thermo(self):
        import warnings
        from PyOMES.core import Simulation
        from PyOMES.core.links import AdvectiveLink
        from PyOMES.thermo import ThermoFramework
        thermo = ThermoFramework()
        cv_a = self._make_cv(thermo)
        cv_b = self._make_cv(thermo)
        link = AdvectiveLink(
            _source_cv_key="a", _source_phase_key="liquid",
            _sink_cv_key="b",   _sink_phase_key="liquid",
            Q_L_per_h=10.0,
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            Simulation(cvs={"a": cv_a, "b": cv_b}, links=[link])
        mismatch_warns = [x for x in w if "ThermoFramework" in str(x.message)]
        assert len(mismatch_warns) == 0

    def test_warning_different_thermo(self):
        import warnings
        from PyOMES.core import Simulation
        from PyOMES.core.links import AdvectiveLink
        from PyOMES.thermo import ThermoFramework
        import dataclasses
        from PyOMES.thermo import DaviesLiquidModel, IdealLiquidModel
        t1 = ThermoFramework(liquid_activity=DaviesLiquidModel())
        t2 = dataclasses.replace(t1, liquid_activity=IdealLiquidModel())
        cv_a = self._make_cv(t1)
        cv_b = self._make_cv(t2)
        link = AdvectiveLink(
            _source_cv_key="a", _source_phase_key="liquid",
            _sink_cv_key="b",   _sink_phase_key="liquid",
            Q_L_per_h=10.0,
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            Simulation(cvs={"a": cv_a, "b": cv_b}, links=[link])
        mismatch_warns = [x for x in w if "ThermoFramework" in str(x.message)]
        assert len(mismatch_warns) == 1


class TestControlVolumeChemistryDb:
    def test_cv_accepts_chemistry_db(self):
        from PyOMES.core.control_volume import ControlVolume
        from PyOMES.core.phases import LiquidPhase
        from PyOMES.chemistry.databases.aqueous import AQUEOUS_DEFAULT
        liq = LiquidPhase(n_mol={"CO2": 0.01}, V_L=1.0, T_K=298.15)
        cv = ControlVolume(
            phases={"liquid": liq},
            chemistry_db=AQUEOUS_DEFAULT,
        )
        assert cv.chemistry_db is AQUEOUS_DEFAULT
        assert cv.chemistry_db.thermo is not None

    def test_cv_without_chemistry_db_defaults_none(self):
        from PyOMES.core.control_volume import ControlVolume
        from PyOMES.core.phases import LiquidPhase
        liq = LiquidPhase(n_mol={}, V_L=1.0, T_K=298.15)
        cv = ControlVolume(phases={"liquid": liq})
        assert cv.chemistry_db is None
