# -*- coding: utf-8 -*-
"""Tests for PHREEQCChemicalEquilibriumEngine (Phase 3 of the protocol hierarchy).

Structure
---------
TestNameTranslation     — phreeqc_to_vlsim and _strip_oxidation_state;
                          no phreeqpython dependency, always runs.
TestImportError         — PHREEQCChemicalEquilibriumEngine raises ImportError when phreeqpython
                          is absent; uses monkeypatch, always runs.
TestPHREEQCChemicalEquilibriumEngine       — Protocol conformance, solve(), algebraic_species(),
                          reset_cache(), reset_counters(), phases writeback.
                          Skipped automatically when phreeqpython is absent.
"""
from __future__ import annotations

import sys

import pytest

from PyOMES.chemical_equilibrium.phreeqc_engine import (
    PHREEQCChemicalEquilibriumEngine,
    _strip_oxidation_state,
    phreeqc_to_vlsim,
)


# ---------------------------------------------------------------------------
# 1. Name-translation utilities (no phreeqpython required)
# ---------------------------------------------------------------------------

class TestNameTranslation:
    """phreeqc_to_vlsim: exceptions, structural regex, and pass-through."""

    # Exceptions dict
    def test_fe2_exception(self):
        assert phreeqc_to_vlsim("Fe+2") == "Fe2+"

    def test_fe3_exception(self):
        assert phreeqc_to_vlsim("Fe+3") == "Fe3+"

    # Structural regex: X+N → X + "+" * N
    def test_ca_plus2(self):
        assert phreeqc_to_vlsim("Ca+2") == "Ca++"

    def test_mg_plus2(self):
        assert phreeqc_to_vlsim("Mg+2") == "Mg++"

    def test_ba_plus2(self):
        assert phreeqc_to_vlsim("Ba+2") == "Ba++"

    # Structural regex: X-N → X + "-" * N
    def test_so4_minus2(self):
        assert phreeqc_to_vlsim("SO4-2") == "SO4--"

    def test_co3_minus2(self):
        assert phreeqc_to_vlsim("CO3-2") == "CO3--"

    def test_hpo4_minus2(self):
        assert phreeqc_to_vlsim("HPO4-2") == "HPO4--"

    def test_po4_minus3(self):
        assert phreeqc_to_vlsim("PO4-3") == "PO4---"

    # Pass-through: singly charged (no digit suffix)
    def test_h_plus_passthrough(self):
        assert phreeqc_to_vlsim("H+") == "H+"

    def test_oh_minus_passthrough(self):
        assert phreeqc_to_vlsim("OH-") == "OH-"

    def test_na_plus_passthrough(self):
        assert phreeqc_to_vlsim("Na+") == "Na+"

    def test_cl_minus_passthrough(self):
        assert phreeqc_to_vlsim("Cl-") == "Cl-"

    def test_hco3_minus_passthrough(self):
        assert phreeqc_to_vlsim("HCO3-") == "HCO3-"

    def test_nh4_plus_passthrough(self):
        assert phreeqc_to_vlsim("NH4+") == "NH4+"

    # Neutral species (no charge suffix at all)
    def test_neutral_species_passthrough(self):
        assert phreeqc_to_vlsim("CO2") == "CO2"
        assert phreeqc_to_vlsim("NH3") == "NH3"
        assert phreeqc_to_vlsim("H2O") == "H2O"


class TestStripOxidationState:
    def test_strip_negative_state(self):
        assert _strip_oxidation_state("S(-2)") == "S"

    def test_strip_positive_state(self):
        assert _strip_oxidation_state("C(4)") == "C"

    def test_strip_negative_3(self):
        assert _strip_oxidation_state("N(-3)") == "N"

    def test_strip_positive_6(self):
        assert _strip_oxidation_state("S(6)") == "S"

    def test_bare_name_unchanged(self):
        assert _strip_oxidation_state("Ca") == "Ca"
        assert _strip_oxidation_state("Na") == "Na"
        assert _strip_oxidation_state("Fe") == "Fe"

    def test_species_with_charge_unchanged(self):
        # Species like "Ca+2" don't have () notation — should pass through
        assert _strip_oxidation_state("Ca+2") == "Ca+2"


# ---------------------------------------------------------------------------
# 2. ImportError when phreeqpython is absent (monkeypatch, always runs)
# ---------------------------------------------------------------------------

class TestImportError:
    def test_raises_import_error_without_phreeqpython(self, monkeypatch):
        """PHREEQCChemicalEquilibriumEngine.__init__ must raise ImportError with install hint."""
        # Simulate phreeqpython being absent regardless of actual installation
        monkeypatch.setitem(sys.modules, "phreeqpython", None)
        with pytest.raises(ImportError, match="phreeqpython"):
            PHREEQCChemicalEquilibriumEngine({"CO2": 1.0}, component_map={"CO2": "C"})

    def test_import_error_message_contains_pip(self, monkeypatch):
        monkeypatch.setitem(sys.modules, "phreeqpython", None)
        with pytest.raises(ImportError, match="pip install"):
            PHREEQCChemicalEquilibriumEngine({"CO2": 1.0}, component_map={"CO2": "C"})


# ---------------------------------------------------------------------------
# 3. PHREEQCChemicalEquilibriumEngine integration tests (require phreeqpython)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def phreeqpython_mod():
    return pytest.importorskip("phreeqpython")


@pytest.fixture(scope="module")
def carbonate_engine(phreeqpython_mod):
    """PHREEQCChemicalEquilibriumEngine for a pure carbonate system."""
    return PHREEQCChemicalEquilibriumEngine(
        {"CO2": 1.0},
        component_map={"CO2": "C(4)"},
        T_C=25.0,
        use_warmstart=True,
    )


@pytest.fixture(scope="module")
def carbonate_nh3_engine(phreeqpython_mod):
    """PHREEQCChemicalEquilibriumEngine for carbonate + ammonia."""
    return PHREEQCChemicalEquilibriumEngine(
        {"CO2": 1.0, "NH3": 0.5},
        component_map={"CO2": "C(4)", "NH3": "N(-3)"},
        T_C=25.0,
        use_warmstart=True,
    )


class TestProtocolConformance:
    def test_satisfies_speciation_engine_protocol(self, carbonate_engine):
        from PyOMES.chemical_equilibrium.protocols import ChemicalEquilibriumEngineProtocol
        assert isinstance(carbonate_engine, ChemicalEquilibriumEngineProtocol)

    def test_not_graybox(self, carbonate_engine):
        from PyOMES.chemical_equilibrium.protocols import GrayBoxEngineProtocol
        assert not isinstance(carbonate_engine, GrayBoxEngineProtocol)

    def test_not_whitebox(self, carbonate_engine):
        from PyOMES.chemical_equilibrium.protocols import WhiteBoxEngineProtocol
        assert not isinstance(carbonate_engine, WhiteBoxEngineProtocol)


class TestSolveOutput:
    def test_returns_ph_key(self, carbonate_engine):
        out = carbonate_engine.solve(totals={"CO2": 0.010})
        assert out.pH is not None

    def test_returns_logh_key(self, carbonate_engine):
        out = carbonate_engine.solve(totals={"CO2": 0.010})
        assert out.logH is not None

    def test_returns_ionic_strength_key(self, carbonate_engine):
        out = carbonate_engine.solve(totals={"CO2": 0.010})
        assert out.ionic_strength is not None

    def test_logh_equals_minus_ph(self, carbonate_engine):
        out = carbonate_engine.solve(totals={"CO2": 0.010})
        assert out.logH == pytest.approx(-out.pH, abs=1e-9)

    def test_ph_in_reasonable_range(self, carbonate_engine):
        out = carbonate_engine.solve(totals={"CO2": 0.010})
        assert 3.0 < out.pH < 9.0, f"pH={out.pH} out of expected range"

    def test_ionic_strength_positive(self, carbonate_engine):
        out = carbonate_engine.solve(totals={"CO2": 0.010})
        assert out.ionic_strength >= 0.0

    def test_species_concentrations_present(self, carbonate_engine):
        out = carbonate_engine.solve(totals={"CO2": 0.010})
        sp = out.species_mol_L
        assert "HCO3-" in sp or "CO3--" in sp or "CO2" in sp

    def test_vlsim_naming_applied(self, carbonate_engine):
        out = carbonate_engine.solve(totals={"CO2": 0.010})
        # PHREEQC uses "CO3-2"; species_map should convert to "CO3--"
        assert "CO3-2" not in out.species_mol_L, "PHREEQC raw name should be translated"
        assert "CO3--" in out.species_mol_L, "PyOMES name CO3-- should be present"

    def test_carbonate_nh3_ph_above_pure_co2(
        self, carbonate_engine, carbonate_nh3_engine
    ):
        out_co2 = carbonate_engine.solve(totals={"CO2": 0.010})
        out_both = carbonate_nh3_engine.solve(totals={"CO2": 0.010, "NH3": 0.010})
        # Adding NH3 (base) should raise pH
        assert out_both.pH > out_co2.pH, (
            f"NH3 should raise pH: {out_both.pH:.3f} vs {out_co2.pH:.3f}"
        )


class TestAlgebraicSpecies:
    def test_returns_frozenset(self, carbonate_engine):
        alg = carbonate_engine.algebraic_species()
        assert isinstance(alg, frozenset)

    def test_non_empty(self, carbonate_engine):
        assert len(carbonate_engine.algebraic_species()) > 0

    def test_contains_h_plus(self, carbonate_engine):
        alg = carbonate_engine.algebraic_species()
        assert "H+" in alg

    def test_contains_oh_minus(self, carbonate_engine):
        alg = carbonate_engine.algebraic_species()
        assert "OH-" in alg

    def test_contains_hco3_minus(self, carbonate_engine):
        alg = carbonate_engine.algebraic_species()
        assert "HCO3-" in alg

    def test_contains_co3(self, carbonate_engine):
        alg = carbonate_engine.algebraic_species()
        assert "CO3--" in alg, f"algebraic_species={alg}"

    def test_vlsim_names_only(self, carbonate_engine):
        alg = carbonate_engine.algebraic_species()
        for sp in alg:
            assert "CO3-2" != sp, "PHREEQC raw name CO3-2 should be translated"


class TestWarmstart:
    def test_same_ph_warmstart_vs_fresh(self, phreeqpython_mod):
        warm = PHREEQCChemicalEquilibriumEngine(
            {"CO2": 1.0}, component_map={"CO2": "C(4)"}, use_warmstart=True
        )
        cold = PHREEQCChemicalEquilibriumEngine(
            {"CO2": 1.0}, component_map={"CO2": "C(4)"}, use_warmstart=False
        )
        for ct in (0.001, 0.005, 0.010, 0.050):
            out_w = warm.solve(totals={"CO2": ct})
            out_c = cold.solve(totals={"CO2": ct})
            assert out_w.pH == pytest.approx(out_c.pH, abs=0.02), (
                f"warmstart and fresh disagree at CT_CO2={ct}: "
                f"{out_w.pH:.4f} vs {out_c.pH:.4f}"
            )

    def test_warmstart_flag_does_not_retain_solution(self, phreeqpython_mod):
        eng = PHREEQCChemicalEquilibriumEngine(
            {"CO2": 1.0}, component_map={"CO2": "C(4)"}, use_warmstart=True
        )
        assert eng.use_warmstart is True
        assert eng._sol is None
        eng.solve(totals={"CO2": 0.010})
        assert eng._sol is None


class TestCounters:
    def test_n_solve_calls_increments(self, phreeqpython_mod):
        eng = PHREEQCChemicalEquilibriumEngine(
            {"CO2": 1.0}, component_map={"CO2": "C(4)"}
        )
        assert eng.n_solve_calls == 0
        eng.solve(totals={"CO2": 0.010})
        assert eng.n_solve_calls == 1
        eng.solve(totals={"CO2": 0.020})
        assert eng.n_solve_calls == 2

    def test_reset_counters(self, phreeqpython_mod):
        eng = PHREEQCChemicalEquilibriumEngine(
            {"CO2": 1.0}, component_map={"CO2": "C(4)"}
        )
        eng.solve(totals={"CO2": 0.010})
        eng.reset_counters()
        assert eng.n_solve_calls == 0


class TestResetCache:
    def test_reset_cache_allows_continued_solving(self, phreeqpython_mod):
        eng = PHREEQCChemicalEquilibriumEngine(
            {"CO2": 1.0}, component_map={"CO2": "C(4)"}, use_warmstart=True
        )
        out_before = eng.solve(totals={"CO2": 0.010})
        eng.reset_cache()
        assert eng._sol is None
        out_after = eng.solve(totals={"CO2": 0.010})
        assert out_after.pH == pytest.approx(out_before.pH, abs=0.02)

    def test_reset_cache_on_fresh_engine_safe(self, phreeqpython_mod):
        eng = PHREEQCChemicalEquilibriumEngine(
            {"CO2": 1.0}, component_map={"CO2": "C(4)"}, use_warmstart=False
        )
        eng.reset_cache()  # should not raise even if _sol is None


class TestNullSpeciesMap:
    def test_no_translation_when_species_map_none(self, phreeqpython_mod):
        eng = PHREEQCChemicalEquilibriumEngine(
            {"CO2": 1.0},
            component_map={"CO2": "C(4)"},
            species_map=None,
        )
        out = eng.solve(totals={"CO2": 0.010})
        # With species_map=None, PHREEQC raw names should appear
        assert "CO3-2" in out.species_mol_L, "raw PHREEQC name should appear when species_map=None"
        assert "CO3--" not in out.species_mol_L


class TestPhasesWriteback:
    """Verify writeback to phase.n_mol via EquilibriumResult.apply_to_phases()."""

    def test_writeback_called(self, phreeqpython_mod):
        from unittest.mock import MagicMock

        eng = PHREEQCChemicalEquilibriumEngine(
            {"CO2": 1.0}, component_map={"CO2": "C(4)"}
        )
        liq = MagicMock()
        liq.V_L = 1.0
        liq.n_mol = {"CO2": 0.010}
        phases = {"liquid": liq}

        result = eng.solve(phases=phases)
        result.apply_to_phases(phases)  # solve() no longer auto-writes
        liq._refresh_derived.assert_called_once()

    def test_writeback_contains_hco3(self, phreeqpython_mod):
        from unittest.mock import MagicMock

        eng = PHREEQCChemicalEquilibriumEngine(
            {"CO2": 1.0}, component_map={"CO2": "C(4)"}
        )
        written = {}
        liq = MagicMock()
        liq.V_L = 1.0
        liq.n_mol = {"CO2": 0.010}
        liq._refresh_derived.side_effect = lambda d: written.update(d)
        phases = {"liquid": liq}

        result = eng.solve(phases=phases)
        result.apply_to_phases(phases)  # solve() no longer auto-writes
        assert "HCO3-" in written, f"HCO3- not in writeback dict: {list(written.keys())}"
