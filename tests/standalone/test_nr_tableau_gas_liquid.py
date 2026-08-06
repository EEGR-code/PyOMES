# -*- coding: utf-8 -*-
"""Tests for CP1 of LAYER1_GAP_CLOSURE: build_tableau() gas-liquid folding.

Topology (per the CP1/CP2 split agreed for this phase) — these tests
assert on tableau *structure* (masters, components, secondaries, their
``phase``/``gas_species_ids`` tags) using the bespoke water + carbonate +
ammonia + sulfide chemistry already established in
``test_nr_speciation_engine.py``, extended with Henry gas-liquid
declarations. CP2's volume-aware residual/Jacobian numerics (and the
consistency-vs-SNIA tests that exercise them) live in
``test_nr_gas_liquid_cp2.py`` — this file also keeps the *permanent*
solve-time contract (folded tableaus require explicit volumes) close to
the topology tests that motivate it.

Three cases per the plan's confirmed-safe folding limits
(``MASS_EXCHANGE_ARCHITECTURE.md`` §14.3):

* **Attach**: CO2, NH3, H2S — gas forms of species that already have their
  own single-component acid-base ladder attach onto that component.
* **Singleton**: O2, CH4, N2, H2 — inert gases with no acid-base ladder
  form a trivial standalone component.
* **Bridging error**: a constructed pathological reaction whose gas-phase
  entry is tied to liquid-phase entries from two independent existing
  components must raise ``ConfigurationError``, not silently misattach.
"""
from __future__ import annotations

import pytest


def _base_reactions():
    """Water + two-step carbonate + ammonia — same chemistry as
    test_nr_speciation_engine.py's _make_reactions(), reused here so the
    carbonate/ammonia components already exist before gas-liquid folding
    is attempted onto them."""
    from PyOMES.chemistry.common_species import (
        H2O, H_plus, OH_minus,
        CO2, HCO3_minus, CO3_2minus,
        NH3, NH4_plus,
    )
    from PyOMES.reactions.equilibrium import EquilibriumReaction
    from PyOMES.reactions.stoichiometry import StoichiometryEntry

    def _e(species, phase, coeff):
        return StoichiometryEntry(species=species, phase=phase, coefficient=coeff)

    water = EquilibriumReaction(
        stoichiometry=[
            _e(H2O, "liquid", -1.0),
            _e(H_plus, "liquid", +1.0),
            _e(OH_minus, "liquid", +1.0),
        ],
        log_K=-14.0,
        balance_elements=("H", "O"),
        label="water",
    )
    co2_first = EquilibriumReaction(
        stoichiometry=[
            _e(CO2, "liquid", -1.0),
            _e(H2O, "liquid", -1.0),
            _e(HCO3_minus, "liquid", +1.0),
            _e(H_plus, "liquid", +1.0),
        ],
        log_K=-6.35,
        total_id="CO2",
        balance_elements=("C", "H", "O"),
        label="co2_first",
    )
    co2_second = EquilibriumReaction(
        stoichiometry=[
            _e(HCO3_minus, "liquid", -1.0),
            _e(CO3_2minus, "liquid", +1.0),
            _e(H_plus, "liquid", +1.0),
        ],
        log_K=-10.33,
        total_id="CO2",
        balance_elements=("C", "H", "O"),
        label="co2_second",
    )
    nh4 = EquilibriumReaction(
        stoichiometry=[
            _e(NH4_plus, "liquid", -1.0),
            _e(NH3, "liquid", +1.0),
            _e(H_plus, "liquid", +1.0),
        ],
        log_K=-9.25,
        total_id="NH3",
        balance_elements=("N", "H"),
        label="nh4",
    )
    return [water, co2_first, co2_second, nh4]


def _h2s_ladder():
    """H2S <-> HS- + H+, pKa=7.0 — a third single-component acid-base
    ladder to exercise the CO2/NH3/H2S 'attach' case named in the plan."""
    from PyOMES.chemistry.common_species import H2S, HS_minus, H_plus
    from PyOMES.reactions.equilibrium import EquilibriumReaction
    from PyOMES.reactions.stoichiometry import StoichiometryEntry

    return EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=H2S, phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=HS_minus, phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=H_plus, phase="liquid", coefficient=+1.0),
        ],
        log_K=-7.0,
        balance_elements=("H", "S"),
        label="eq_H2S",
    )


def _co2_henry():
    from PyOMES.chemistry import HenryEquilibrium
    return HenryEquilibrium(
        H_ref=3.4e-4, dlnH=2400.0, gas_species="CO2", liquid_species="CO2",
        label="henry_CO2",
    )


def _nh3_henry():
    from PyOMES.chemistry import HenryEquilibrium
    return HenryEquilibrium(
        H_ref=5.9e-1, dlnH=4200.0, gas_species="NH3", liquid_species="NH3",
        label="henry_NH3",
    )


def _h2s_henry():
    from PyOMES.chemistry import HenryEquilibrium
    return HenryEquilibrium(
        H_ref=1.0e-3, dlnH=2100.0, gas_species="H2S", liquid_species="H2S",
        label="henry_H2S",
    )


def _inert_gas_species():
    """O2/CH4/N2/H2 — not in common_species; local Species objects with no
    acid-base chemistry, matching the plan's inert-gas singleton case."""
    from PyOMES.chemistry.species import Species
    return {
        "O2": Species(id="O2", atoms={"O": 2}, charge=0),
        "CH4": Species(id="CH4", atoms={"C": 1, "H": 4}, charge=0),
        "N2": Species(id="N2", atoms={"N": 2}, charge=0),
        "H2": Species(id="H2", atoms={"H": 2}, charge=0),
    }


def _inert_henry(species_id, species_obj, H_ref, dlnH):
    from PyOMES.chemistry import HenryEquilibrium
    return HenryEquilibrium(
        H_ref=H_ref, dlnH=dlnH,
        gas_species=species_obj, liquid_species=species_obj,
        label=f"henry_{species_id}",
    )


# ═══════════════════════════════════════════════════════════════════════════
#  Attach case: CO2 / NH3 / H2S onto their existing acid-base ladders
# ═══════════════════════════════════════════════════════════════════════════

class TestAttachOntoExistingComponent:

    @pytest.fixture(scope="class")
    def tableau(self):
        from PyOMES.chemical_equilibrium.nr_tableau import build_tableau
        reactions = _base_reactions() + [_h2s_ladder(), _co2_henry(), _nh3_henry(), _h2s_henry()]
        return build_tableau(reactions, T_K=298.15)

    def test_masters_unchanged_by_gas_folding(self, tableau):
        """Folding gas-liquid rows must not introduce new masters — CO2/
        NH3/H2S gas forms attach onto their existing liquid masters."""
        assert set(tableau.masters) == {"H+", "CO2", "NH3", "H2S"}

    def test_gas_secondaries_present_with_gas_phase_tag(self, tableau):
        gas_secs = {s.species_id: s for s in tableau.secondaries if s.phase == "gas"}
        assert set(gas_secs) == {"CO2", "NH3", "H2S"}
        for sec in gas_secs.values():
            assert sec.phase == "gas"

    def test_liquid_secondaries_still_tagged_liquid(self, tableau):
        liquid_secs = [s for s in tableau.secondaries if s.phase == "liquid"]
        liquid_ids = {s.species_id for s in liquid_secs}
        assert {"OH-", "HCO3-", "CO3--", "NH4+", "HS-"} <= liquid_ids
        assert all(s.phase == "liquid" for s in liquid_secs)

    def test_component_gas_species_ids(self, tableau):
        comp_by_master = {c.master_id: c for c in tableau.components}
        assert comp_by_master["CO2"].gas_species_ids == ("CO2",)
        assert comp_by_master["NH3"].gas_species_ids == ("NH3",)
        assert comp_by_master["H2S"].gas_species_ids == ("H2S",)

    def test_component_liquid_species_ids_unaffected(self, tableau):
        """species_ids (liquid membership) is unchanged by gas folding —
        gas ids live only in the new gas_species_ids field."""
        comp_by_master = {c.master_id: c for c in tableau.components}
        assert set(comp_by_master["CO2"].species_ids) == {"CO2", "HCO3-", "CO3--"}
        assert set(comp_by_master["NH3"].species_ids) == {"NH3", "NH4+"}
        assert set(comp_by_master["H2S"].species_ids) == {"H2S", "HS-"}

    def test_no_gas_species_ids_on_unrelated_component(self, tableau):
        """Components with no folded gas-liquid equilibrium keep an empty
        gas_species_ids — there are none here since all three components
        have a Henry declaration, so this is exercised in the singleton
        test class instead (mixed system)."""
        for comp in tableau.components:
            assert isinstance(comp.gas_species_ids, tuple)

    def test_gas_co2_log_linear_relation(self, tableau):
        """log(p_CO2,gas) = x_CO2,liquid - log10(kH) — derived Henry
        relation should have nu={CO2: +1.0} and log_K_prime = -log10(kH)."""
        import math
        co2_henry = _co2_henry()
        gas_co2 = next(s for s in tableau.secondaries
                        if s.species_id == "CO2" and s.phase == "gas")
        assert gas_co2.nu == pytest.approx({"CO2": 1.0})
        expected_log_K_prime = -math.log10(co2_henry._kH_mol_L_atm(298.15))
        assert gas_co2.log_K_prime == pytest.approx(expected_log_K_prime, rel=1e-9)


# ═══════════════════════════════════════════════════════════════════════════
#  Singleton case: O2 / CH4 / N2 / H2 — no acid-base ladder
# ═══════════════════════════════════════════════════════════════════════════

class TestSingletonComponent:

    @pytest.fixture(scope="class")
    def tableau(self):
        from PyOMES.chemical_equilibrium.nr_tableau import build_tableau
        inert = _inert_gas_species()
        henries = [
            _inert_henry("O2", inert["O2"], 1.3e-5, 1500.0),
            _inert_henry("CH4", inert["CH4"], 1.4e-5, 1600.0),
            _inert_henry("N2", inert["N2"], 6.4e-6, 1300.0),
            _inert_henry("H2", inert["H2"], 7.7e-6, 500.0),
        ]
        reactions = _base_reactions() + henries
        return build_tableau(reactions, T_K=298.15)

    def test_inert_gases_become_trivial_masters(self, tableau):
        """No acid-base ladder for O2/CH4/N2/H2 means each forms its own
        singleton component with itself as master."""
        assert {"O2", "CH4", "N2", "H2"} <= set(tableau.masters)

    def test_singleton_components_have_single_liquid_member(self, tableau):
        comp_by_master = {c.master_id: c for c in tableau.components}
        for sp_id in ("O2", "CH4", "N2", "H2"):
            comp = comp_by_master[sp_id]
            assert comp.species_ids == (sp_id,)
            assert comp.gas_species_ids == (sp_id,)

    def test_singleton_gas_secondaries_tagged_gas(self, tableau):
        gas_secs = {s.species_id for s in tableau.secondaries if s.phase == "gas"}
        assert {"O2", "CH4", "N2", "H2"} <= gas_secs

    def test_singleton_master_charge_zero(self, tableau):
        for sp_id in ("O2", "CH4", "N2", "H2"):
            assert tableau.master_charges[sp_id] == 0

    def test_carbonate_ammonia_components_unaffected(self, tableau):
        """Mixed system: CO2/NH3 components exist alongside the new
        inert-gas singletons, with no gas_species_ids of their own (no
        Henry declared for them in this fixture)."""
        comp_by_master = {c.master_id: c for c in tableau.components}
        assert comp_by_master["CO2"].gas_species_ids == ()
        assert comp_by_master["NH3"].gas_species_ids == ()
        assert set(comp_by_master["CO2"].species_ids) == {"CO2", "HCO3-", "CO3--"}


# ═══════════════════════════════════════════════════════════════════════════
#  Bridging case: ConfigurationError, not silent misattachment
# ═══════════════════════════════════════════════════════════════════════════

class TestBridgingRaisesConfigurationError:

    def _pathological_bridge_reaction(self):
        """A hand-built gas-liquid reaction whose gas-phase entry is tied
        to liquid entries from two *independent* existing components
        (CO2's ladder and NH3's ladder) — structurally analogous to
        carbamate formation (CO2 + NH3 -> gas product), which would
        require merging two multi-species components. balance_elements=()
        skips elemental validation since the placeholder gas species is
        deliberately not a real formula."""
        from PyOMES.chemistry.common_species import CO2, NH3
        from PyOMES.chemistry.species import Species
        from PyOMES.reactions.equilibrium import EquilibriumReaction
        from PyOMES.reactions.stoichiometry import StoichiometryEntry

        placeholder_gas = Species(id="BridgeGas", atoms={}, charge=0)
        return EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=placeholder_gas, phase="gas", coefficient=-1.0),
                StoichiometryEntry(species=CO2, phase="liquid", coefficient=+1.0),
                StoichiometryEntry(species=NH3, phase="liquid", coefficient=+1.0),
            ],
            log_K=1.0,
            balance_elements=(),
            label="pathological_bridge",
        )

    def test_bridging_reaction_raises_configuration_error(self):
        from PyOMES.chemical_equilibrium.nr_tableau import build_tableau, ConfigurationError

        reactions = _base_reactions() + [self._pathological_bridge_reaction()]
        with pytest.raises(ConfigurationError, match="bridge"):
            build_tableau(reactions, T_K=298.15)

    def test_configuration_error_is_a_value_error(self):
        """ConfigurationError subclasses ValueError (matching the existing
        StoichiometryError/SpeciesConflictError convention) so callers
        that broadly catch ValueError still see it, while callers that
        want to distinguish 'out of capability' from 'malformed input'
        can catch ConfigurationError specifically."""
        from PyOMES.chemical_equilibrium.nr_tableau import ConfigurationError
        assert issubclass(ConfigurationError, ValueError)

    def test_non_bridging_system_still_builds(self):
        """Sanity check: the same base chemistry without the pathological
        reaction builds fine (isolates the failure to the bridge itself)."""
        from PyOMES.chemical_equilibrium.nr_tableau import build_tableau
        tableau = build_tableau(_base_reactions() + [_co2_henry()], T_K=298.15)
        assert "CO2" in tableau.masters


# ═══════════════════════════════════════════════════════════════════════════
#  Routing-only declarations remain silently skipped (unchanged behaviour)
# ═══════════════════════════════════════════════════════════════════════════

class TestUnparameterizedGasLiquidStillSkipped:

    def test_log_k_none_cross_phase_reaction_skipped(self):
        """A hand-built cross-phase EquilibriumReaction with log_K=None
        (pure KineticGasLiquidLink routing declaration) must still be
        silently skipped, exactly as before CP1 — it carries no mass-
        action constant to fold in."""
        from PyOMES.chemistry.common_species import CO2
        from PyOMES.reactions.equilibrium import EquilibriumReaction
        from PyOMES.reactions.stoichiometry import StoichiometryEntry
        from PyOMES.chemical_equilibrium.nr_tableau import build_tableau

        routing_only = EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=CO2, phase="gas", coefficient=-1.0),
                StoichiometryEntry(species=CO2, phase="liquid", coefficient=+1.0),
            ],
            label="routing_only",
        )
        tableau = build_tableau(_base_reactions() + [routing_only], T_K=298.15)
        gas_secs = {s.species_id for s in tableau.secondaries if s.phase == "gas"}
        assert "CO2" not in gas_secs
        comp_by_master = {c.master_id: c for c in tableau.components}
        assert comp_by_master["CO2"].gas_species_ids == ()

    def test_partition_model_only_henry_skipped(self):
        """A HenryEquilibrium with no gas_species/liquid_species set has
        an empty stoichiometry (classify_equilibrium_constraint raises
        ValueError internally) and is silently skipped, as before."""
        from PyOMES.chemistry import HenryEquilibrium
        from PyOMES.chemical_equilibrium.nr_tableau import build_tableau

        partition_only = HenryEquilibrium(H_ref=3.4e-4, dlnH=2400.0)
        tableau = build_tableau(_base_reactions() + [partition_only], T_K=298.15)
        assert set(tableau.masters) == {"H+", "CO2", "NH3"}


# ═══════════════════════════════════════════════════════════════════════════
#  Solve-time contract: folded tableaus require explicit volumes (CP2)
# ═══════════════════════════════════════════════════════════════════════════
# CP1 only built tableau topology; CP2 (see test_nr_gas_liquid_cp2.py) wires
# the actual volume-aware residual/Jacobian math. What CP1's own test file
# still needs to cover: the *permanent* contract that solving a folded
# tableau without volumes is refused, not silently wrong.

class TestSolveRequiresVolumes:

    def test_solve_nr_raises_without_volumes(self):
        from PyOMES.chemical_equilibrium.nr_tableau import build_tableau
        from PyOMES.chemical_equilibrium.nr_solver import solve_nr
        from PyOMES.chemical_equilibrium.activity_models import make_activity_model

        tableau = build_tableau(_base_reactions() + [_co2_henry()], T_K=298.15)
        am = make_activity_model(False, "ideal")
        with pytest.raises(ValueError, match="V_liq_L"):
            solve_nr(tableau, {"CO2": 0.05, "NH3": 0.04}, {}, activity_model=am)

    def test_engine_solve_raises_without_volumes(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

        engine = NRChemicalEquilibriumEngine.from_reactions(
            _base_reactions() + [_co2_henry()], T_K=298.15,
        )
        assert "CO2" in engine.tableau.masters
        with pytest.raises(ValueError, match="V_liq_L"):
            engine.solve(totals={"CO2": 0.05, "NH3": 0.04}, strong_ions={})

    def test_engine_solve_succeeds_with_volumes(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

        engine = NRChemicalEquilibriumEngine.from_reactions(
            _base_reactions() + [_co2_henry()], T_K=298.15,
        )
        out = engine.solve(
            totals={"CO2": 0.05, "NH3": 0.04}, strong_ions={},
            V_liq_L=1.0, V_gas_L=0.1,
        )
        assert out.species_mol_L["CO2"] > 0.0
        assert "CO2" in out.partial_pressures_atm

    def test_unfolded_tableau_still_solves_without_volumes(self):
        """Sanity/regression check: a tableau with no gas-liquid folding
        (the pre-CP1 case) is completely unaffected by the volume
        requirement."""
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

        engine = NRChemicalEquilibriumEngine.from_reactions(_base_reactions(), T_K=298.15)
        out = engine.solve(totals={"CO2": 0.05, "NH3": 0.04}, strong_ions={})
        assert out.species_mol_L["CO2"] > 0.0
        assert out.partial_pressures_atm == {}

    def test_ad_basic_database_reaction_list_requires_volumes(self):
        """The real-world case that motivated this guard: AD_BASIC's
        _CO2_HENRY is a fully-parameterized HenryEquilibrium (log_K set)
        already shipped in a reaction list consumed elsewhere. If it were
        ever fed into NRChemicalEquilibriumEngine.from_reactions() directly, CP1's
        topology change folds it — this confirms solving it without
        volumes is refused rather than silently producing a contaminated
        result."""
        from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

        # AD_BASIC.reactions already includes a water-dissociation reaction
        # (inherited from BIOPROCESS_BASIC) alongside _CO2_HENRY.
        engine = NRChemicalEquilibriumEngine.from_reactions(
            list(AD_BASIC.reactions.reactions), T_K=298.15,
        )
        comp_by_master = {c.master_id: c for c in engine.tableau.components}
        assert comp_by_master["CO2"].gas_species_ids == ("CO2",)
        with pytest.raises(ValueError, match="V_liq_L"):
            engine.solve(totals={"CO2": 0.05}, strong_ions={})
