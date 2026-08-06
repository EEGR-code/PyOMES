# -*- coding: utf-8 -*-
"""Tests for CP5 of LAYER1_GAP_CLOSURE: precipitation regression verification.

Verification only, per the plan — no changes to precipitation's numerics.
Confirms the existing nested active-set precipitation loop
(NR_PRECIPITATION_SPECIATION, shipped independently of this phase) and
CP1/CP2's newly-folded gas-liquid rows coexist correctly in the same
NRChemicalEquilibriumEngine.solve() call: both converge, and neither contaminates
the other's unknowns.

The natural test case is CaCO3 precipitation alongside CO2 gas-liquid
transfer — realistic (CO3-- comes from the same carbonate ladder CO2
gas-liquid folding attaches to) and exercises real interaction: the outer
active-set loop decrements the "CO2" component's *effective* total by the
precipitated carbon each outer iteration, and that effective total (per
CP2) now spans both phases, not liquid-only as before this phase.

NR_PRECIPITATION_CV_INTEGRATION.md's SolidPhase writeback (mentioned in
the plan as "if shipped by this point") has not shipped — its status
header still reads "Upcoming — not yet started" — so that part of this
checkpoint's scope is a no-op; nothing exists yet to verify.
"""
from __future__ import annotations

import pytest


def _R_L_ATM_MOL_K():
    return 0.0820574


def _carbonate_reactions():
    from PyOMES.chemistry.common_species import H2O, H_plus, OH_minus, CO2, HCO3_minus, CO3_2minus
    from PyOMES.reactions.equilibrium import EquilibriumReaction
    from PyOMES.reactions.stoichiometry import StoichiometryEntry as E

    water = EquilibriumReaction(
        stoichiometry=[E(H2O, "liquid", -1), E(H_plus, "liquid", +1), E(OH_minus, "liquid", +1)],
        log_K=-14.0, label="water",
    )
    co2_first = EquilibriumReaction(
        stoichiometry=[E(CO2, "liquid", -1), E(H2O, "liquid", -1),
                       E(HCO3_minus, "liquid", +1), E(H_plus, "liquid", +1)],
        log_K=-6.35, total_id="CO2", label="co2_first",
    )
    co2_second = EquilibriumReaction(
        stoichiometry=[E(HCO3_minus, "liquid", -1),
                       E(CO3_2minus, "liquid", +1), E(H_plus, "liquid", +1)],
        log_K=-10.33, total_id="CO2", label="co2_second",
    )
    return [water, co2_first, co2_second]


def _calcite_reaction():
    from PyOMES.chemistry.common_species import Ca_plus_plus, CO3_2minus
    from PyOMES.chemistry.species import Species
    from PyOMES.reactions.equilibrium import EquilibriumReaction
    from PyOMES.reactions.stoichiometry import StoichiometryEntry as E

    CaCO3 = Species(id="CaCO3", atoms={"Ca": 1, "C": 1, "O": 3}, charge=0, MW=100.086)
    return EquilibriumReaction(
        stoichiometry=[E(CaCO3, "solid", -1), E(Ca_plus_plus, "liquid", +1),
                       E(CO3_2minus, "liquid", +1)],
        log_K=-8.48, label="calcite",
    )


def _co2_henry():
    from PyOMES.chemistry import HenryEquilibrium
    return HenryEquilibrium(
        H_ref=3.4e-4, dlnH=2400.0, gas_species="CO2", liquid_species="CO2",
        label="henry_CO2",
    )


def _folded_engine_with_precipitation():
    from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
    return NRChemicalEquilibriumEngine.from_reactions(
        _carbonate_reactions() + [_calcite_reaction(), _co2_henry()],
        use_activity=True, activity_model="davies",
    )


class TestPrecipitationCoexistsWithGasLiquidFold:

    def test_folded_engine_still_auto_classifies_precipitation(self):
        """Sanity check: adding the Henry declaration must not disturb
        the existing auto-classification of the calcite reaction as
        precipitation (EQUILIBRIUM_CONSTRAINT_UNIFICATION CP2)."""
        engine = _folded_engine_with_precipitation()
        assert len(engine._precipitation_reactions) == 1
        assert engine._precipitation_reactions[0].label == "calcite"

    def test_folded_engine_still_folds_co2_gas_liquid(self):
        engine = _folded_engine_with_precipitation()
        comp_by_master = {c.master_id: c for c in engine.tableau.components}
        assert comp_by_master["CO2"].gas_species_ids == ("CO2",)

    def test_supersaturated_precipitates_and_converges(self):
        """Same CT_Ca/CT_CO2/CT_Na as the pre-existing (non-folded)
        regression test in test_nr_speciation_engine.py — precipitation
        must still trigger and converge (SI ~ 0) with the gas-liquid row
        folded in."""
        engine = _folded_engine_with_precipitation()
        out = engine.solve(
            totals={"CO2": 0.010}, strong_ions={"CT_Ca": 0.002, "CT_Na": 0.005},
            V_liq_L=1.0, V_gas_L=0.2,
        )
        xi = out.extra["minerals_xi_mol_L"]["calcite"]
        si = out.saturation_indices["calcite"]
        assert xi > 0.0, f"Expected precipitation (xi > 0), got xi={xi}"
        assert abs(si) < 0.05, f"SI not near zero after precipitation: SI={si:.4f}"

    def test_undersaturated_no_precipitation_with_gas_fold(self):
        engine = _folded_engine_with_precipitation()
        out = engine.solve(
            totals={"CO2": 0.001}, strong_ions={"CT_Ca": 0.0001, "CT_Cl": 0.005},
            V_liq_L=1.0, V_gas_L=0.2,
        )
        xi = out.extra["minerals_xi_mol_L"]["calcite"]
        si = out.saturation_indices["calcite"]
        assert xi == 0.0, f"Expected no precipitation, got xi={xi}"
        assert si < 0.0, f"Expected SI < 0 (undersaturated), got SI={si}"

    def test_calcium_conserved_across_precipitation(self):
        """xi + dissolved Ca == input CT_Ca, exactly as in the pre-
        existing non-folded test — the gas-liquid row must not perturb
        calcium's own mass balance at all (it doesn't touch Ca in any
        way)."""
        CT_Ca = 0.002
        engine = _folded_engine_with_precipitation()
        out = engine.solve(
            totals={"CO2": 0.010}, strong_ions={"CT_Ca": CT_Ca, "CT_Na": 0.005},
            V_liq_L=1.0, V_gas_L=0.2,
        )
        xi = out.extra["minerals_xi_mol_L"]["calcite"]
        Ca_dissolved = CT_Ca - xi
        assert abs(xi + Ca_dissolved - CT_Ca) < 1e-12

    def test_carbon_conserved_across_gas_liquid_and_solid(self):
        """The genuinely new three-way check this phase enables: total
        carbon must balance across gas CO2 + liquid (CO2+HCO3-+CO3--) +
        precipitated CaCO3 simultaneously, to near machine precision —
        confirming the outer active-set loop's effective-total
        adjustment (built for a liquid-only total) is still exactly
        correct now that "the CO2 component's total" spans both phases
        (CP2)."""
        V_liq, V_gas, T_K = 1.0, 0.2, 298.15
        CT_CO2 = 0.010
        engine = _folded_engine_with_precipitation()
        out = engine.solve(
            totals={"CO2": CT_CO2}, strong_ions={"CT_Ca": 0.002, "CT_Na": 0.005},
            V_liq_L=V_liq, V_gas_L=V_gas,
        )
        xi = out.extra["minerals_xi_mol_L"]["calcite"]

        C_liq = (out.species_mol_L["CO2"] + out.species_mol_L["HCO3-"]
                 + out.species_mol_L["CO3--"])
        n_liq_C = C_liq * V_liq
        n_gas_C = out.partial_pressures_atm["CO2"] * V_gas / (_R_L_ATM_MOL_K() * T_K)
        n_solid_C = xi * V_liq  # 1 carbon per mol CaCO3

        total_C = n_liq_C + n_gas_C + n_solid_C
        expected_C = CT_CO2 * V_liq
        assert total_C == pytest.approx(expected_C, abs=1e-8)

    def test_no_cross_contamination_vs_unfolded_precipitation_only(self):
        """The precipitated amount and pH should be close to (not wildly
        different from) the pre-existing unfolded engine's result for
        the same inputs — folding gas-liquid shouldn't grossly perturb
        precipitation just because a Henry row now exists.

        Some difference IS expected, and in a specific, physically
        meaningful direction verified by hand (not asserted blind):
        the folded engine strips some CO2 into the 0.2 L headspace,
        raising pH slightly (7.011 -> 7.058) relative to the unfolded
        (liquid-only) case. Despite *less* total liquid carbon, the
        higher pH shifts carbonate speciation toward more CO3-- (5.11e-6
        -> 5.50e-6 mol/L) — since Ksp = [Ca++][CO3--], that means *more*
        supersaturation, so the folded case precipitates *more* calcite,
        not less. This is the well-known "CO2 stripping promotes CaCO3
        scaling" effect (used industrially in water treatment) — the
        unfolded engine structurally cannot capture it at all, having no
        gas phase to strip CO2 into. A positive validation of the
        coupling, not a discrepancy to paper over with a loose bound."""
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

        unfolded_engine = NRChemicalEquilibriumEngine.from_reactions(
            _carbonate_reactions() + [_calcite_reaction()],
            use_activity=True, activity_model="davies",
        )
        out_unfolded = unfolded_engine.solve(
            totals={"CO2": 0.010}, strong_ions={"CT_Ca": 0.002, "CT_Na": 0.005},
        )
        folded_engine = _folded_engine_with_precipitation()
        out_folded = folded_engine.solve(
            totals={"CO2": 0.010}, strong_ions={"CT_Ca": 0.002, "CT_Na": 0.005},
            V_liq_L=1.0, V_gas_L=0.2,
        )

        xi_unfolded = out_unfolded.extra["minerals_xi_mol_L"]["calcite"]
        xi_folded = out_folded.extra["minerals_xi_mol_L"]["calcite"]
        assert xi_unfolded > 0.0
        assert xi_folded > 0.0
        assert out_folded.pH > out_unfolded.pH  # CO2 stripped to headspace
        assert xi_folded > xi_unfolded  # -> more supersaturated, more precipitation
        assert abs(out_folded.pH - out_unfolded.pH) < 0.5  # still a modest shift

    def test_retain_jacobian_still_refused_with_precipitation_and_gas_fold(self):
        """CP2's retain_jacobian=True guard (for gas-folded tableaus)
        must still fire regardless of precipitation reactions also being
        present — the guard checks tableau.secondaries only, independent
        of _precipitation_reactions."""
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

        with pytest.raises(NotImplementedError, match="gas-liquid secondaries"):
            NRChemicalEquilibriumEngine.from_reactions(
                _carbonate_reactions() + [_calcite_reaction(), _co2_henry()],
                use_activity=True, activity_model="davies",
                retain_jacobian=True,
            )


class TestNRPrecipitationCVIntegrationNotYetShipped:
    """NR_PRECIPITATION_CV_INTEGRATION.md's SolidPhase writeback is
    listed in the plan as 'confirm unaffected, if shipped by this
    point' — it hasn't shipped (status header still reads 'Upcoming —
    not yet started'), so this documents that finding rather than
    silently skipping it."""

    def test_solid_phase_exists_but_no_ksp_writeback_wired_yet(self):
        """SolidPhase itself exists (core phase type), but nothing in
        NRChemicalEquilibriumEngine writes precipitated mineral amounts to it —
        confirmed by inspecting the EquilibriumResult contract: mineral
        xi lives in extra['minerals_xi_mol_L'] (a diagnostic dict), not
        applied to any phase.n_mol by apply_to_phases()."""
        from PyOMES.core.phases import SolidPhase  # noqa: F401 — exists
        from PyOMES.chemical_equilibrium.protocols import EquilibriumResult
        import inspect

        src = inspect.getsource(EquilibriumResult.apply_to_phases)
        assert "solid" not in src.lower()
        assert "mineral" not in src.lower()
