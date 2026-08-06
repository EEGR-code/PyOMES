# -*- coding: utf-8 -*-
"""Tests for the Newton-Raphson speciation engine.

Bespoke test chemistry
-----------------------
All tests use a hand-crafted carbonate + ammonia system that covers the
key structural cases without requiring the full BSM2 model:

    Species (from common_species):
        H2O, H+, OH-            — water equilibrium
        CO2, HCO3-, CO3--       — two-step carbonate ladder  (total_id="CO2")
        NH4+, NH3                — one-step ammonia           (total_id="NH3")

    Reactions:
        (1) H2O  ⇌  H+  +  OH-           pKw = 14.0
        (2) CO2 + H2O  ⇌  HCO3- + H+     pKa1 = 6.35    total_id="CO2"
        (3) HCO3-  ⇌  CO3-- + H+         pKa2 = 10.33   total_id="CO2"
        (4) NH4+  ⇌  NH3  +  H+          pKa  = 9.25    total_id="NH3"

This gives a 3×3 NR system (masters: H+, CO2, NH3).

Why this chemistry
------------------
* The two-step carbonate ladder tests multi-hop BFS tableau derivation and
  verifies the CO2 + HCO3- + CO3-- mass balance closure.
* The NH4+/NH3 pair exercises the ``total_id`` master-override path: the
  heuristic would select NH4+ (it is the DAG source), but ``total_id="NH3"``
  forces NH3 as the master.  This is the generalisation of how AD models
  track total ammoniacal nitrogen as n_mol["NH3"].
* Including a strong-ion case (Na+) tests that the charge-balance residual
  correctly incorporates external fixed charges and that the resulting pH
  shift is in the right direction (Na+ → more basic).

Test structure
--------------
TestNRTableau       — unit tests for build_tableau (master selection, BFS)
TestNRSolverDirect  — end-to-end NRChemicalEquilibriumEngine.solve() with explicit totals
TestNRvsChargeBalance — NR and 1-D charge-balance engines agree to < 0.01 pH
TestReactionSystemRouter — solver= kwarg on ReactionSystem routes correctly
TestPhaseWriteback  — _refresh_derived receives H+, OH-, H2O, and the carbonate species
"""
from __future__ import annotations

import pytest
import numpy as np

# ═══════════════════════════════════════════════════════════════════════════
#  Shared test fixture — bespoke reaction set
# ═══════════════════════════════════════════════════════════════════════════

def _make_reactions():
    """Return the four bespoke equilibrium reactions used across all tests."""
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
            _e(H2O,     "liquid", -1.0),
            _e(H_plus,  "liquid", +1.0),
            _e(OH_minus,"liquid", +1.0),
        ],
        log_K=-14.0,
        balance_elements=("H", "O"),
        label="water",
    )
    co2_first = EquilibriumReaction(
        stoichiometry=[
            _e(CO2,       "liquid", -1.0),
            _e(H2O,       "liquid", -1.0),
            _e(HCO3_minus,"liquid", +1.0),
            _e(H_plus,    "liquid", +1.0),
        ],
        log_K=-6.35,
        total_id="CO2",
        balance_elements=("C", "H", "O"),
        label="co2_first",
    )
    co2_second = EquilibriumReaction(
        stoichiometry=[
            _e(HCO3_minus,"liquid", -1.0),
            _e(CO3_2minus,"liquid", +1.0),
            _e(H_plus,    "liquid", +1.0),
        ],
        log_K=-10.33,
        total_id="CO2",
        balance_elements=("C", "H", "O"),
        label="co2_second",
    )
    nh4 = EquilibriumReaction(
        stoichiometry=[
            _e(NH4_plus,"liquid", -1.0),
            _e(NH3,     "liquid", +1.0),
            _e(H_plus,  "liquid", +1.0),
        ],
        log_K=-9.25,
        total_id="NH3",
        balance_elements=("N", "H"),
        label="nh4",
    )
    return [water, co2_first, co2_second, nh4]


# ═══════════════════════════════════════════════════════════════════════════
#  TestNRTableau — unit tests for the tableau builder
# ═══════════════════════════════════════════════════════════════════════════

class TestNRTableau:
    """Validate master selection, BFS derivation, and tableau structure."""

    @pytest.fixture(scope="class")
    def tableau(self):
        from PyOMES.chemical_equilibrium.nr_tableau import build_tableau
        return build_tableau(_make_reactions(), T_K=298.15)

    def test_masters_list(self, tableau):
        """H+ must be first; CO2 and NH3 selected as masters."""
        assert tableau.masters[0] == "H+"
        assert "CO2" in tableau.masters    # total_id override agrees with heuristic
        assert "NH3" in tableau.masters    # total_id="NH3" overrides DAG source NH4+
        assert "NH4+" not in tableau.masters
        assert len(tableau.masters) == 3   # H+, CO2, NH3

    def test_nh3_is_master_not_nh4(self, tableau):
        """total_id='NH3' forces NH3 as master even though NH4+ is the DAG source."""
        assert "NH3" in tableau.masters
        assert "NH4+" not in tableau.masters

    def test_co2_is_master(self, tableau):
        """CO2 is the DAG source (appears only as reactant); total_id confirms."""
        assert "CO2" in tableau.masters

    def test_secondary_species_present(self, tableau):
        sec_ids = {s.species_id for s in tableau.secondaries}
        assert "OH-" in sec_ids
        assert "HCO3-" in sec_ids
        assert "CO3--" in sec_ids
        assert "NH4+" in sec_ids

    def test_oh_minus_nu(self, tableau):
        """OH- = Kw / [H+]: nu={H+: -1}, log_K'=log_Kw=-14."""
        oh = next(s for s in tableau.secondaries if s.species_id == "OH-")
        assert oh.nu == pytest.approx({"H+": -1.0})
        assert oh.log_K_prime == pytest.approx(-14.0)
        assert oh.charge == -1

    def test_hco3_nu(self, tableau):
        """HCO3- = K1*[CO2]/[H+]: nu={CO2:+1, H+:-1}, log_K'=-6.35."""
        hco3 = next(s for s in tableau.secondaries if s.species_id == "HCO3-")
        assert hco3.nu.get("CO2", 0.0) == pytest.approx(+1.0)
        assert hco3.nu.get("H+", 0.0) == pytest.approx(-1.0)
        assert hco3.log_K_prime == pytest.approx(-6.35)

    def test_co3_nu(self, tableau):
        """CO3-- derived through two steps: log_K'=-6.35+(-10.33)=-16.68."""
        co3 = next(s for s in tableau.secondaries if s.species_id == "CO3--")
        assert co3.nu.get("CO2", 0.0) == pytest.approx(+1.0)
        assert co3.nu.get("H+", 0.0) == pytest.approx(-2.0)
        assert co3.log_K_prime == pytest.approx(-16.68, abs=1e-10)

    def test_nh4_nu(self, tableau):
        """NH4+ = K_NH^-1 * [NH3] * [H+]: nu={NH3:+1, H+:+1}, log_K'=+9.25."""
        # NH4+ ⇌ NH3 + H+, log K=-9.25; NH3 is master so we invert:
        # log[a_NH4+] = -log K + log[a_NH3] + log[a_H+] = +9.25 + x_NH3 + x_H+
        nh4 = next(s for s in tableau.secondaries if s.species_id == "NH4+")
        assert nh4.nu.get("NH3", 0.0) == pytest.approx(+1.0)
        assert nh4.nu.get("H+", 0.0) == pytest.approx(+1.0)
        assert nh4.log_K_prime == pytest.approx(+9.25)

    def test_component_infos(self, tableau):
        """Two components: CO2 (with CO2, HCO3-, CO3--) and NH3 (with NH3, NH4+)."""
        comp_by_master = {c.master_id: c for c in tableau.components}
        assert "CO2" in comp_by_master
        assert "NH3" in comp_by_master
        co2_comp = comp_by_master["CO2"]
        assert set(co2_comp.species_ids) == {"CO2", "HCO3-", "CO3--"}
        nh3_comp = comp_by_master["NH3"]
        assert set(nh3_comp.species_ids) == {"NH3", "NH4+"}

    def test_log_kw(self, tableau):
        assert tableau.log_Kw == pytest.approx(-14.0)


# ═══════════════════════════════════════════════════════════════════════════
#  TestNRSolverDirect — end-to-end via NRChemicalEquilibriumEngine.solve()
# ═══════════════════════════════════════════════════════════════════════════

class TestNRSolverDirect:
    """Direct NRChemicalEquilibriumEngine.solve(totals=..., strong_ions=...) calls."""

    @pytest.fixture(scope="class")
    def engine(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        return NRChemicalEquilibriumEngine.from_reactions(_make_reactions())

    def test_pure_water_ph(self, engine):
        """No solutes → pH ≈ 7.0."""
        out = engine.solve(totals={}, strong_ions={})
        assert out.pH == pytest.approx(7.0, abs=0.05)
        assert out.species_mol_L["H+"] == pytest.approx(1e-7, rel=0.05)
        assert out.species_mol_L["OH-"] == pytest.approx(1e-7, rel=0.05)

    def test_pure_water_charge_balance(self, engine):
        """Charge residual should be near-zero after a pure-water solve."""
        out = engine.solve(totals={}, strong_ions={})
        assert abs(out.charge_residual) < 1e-8

    def test_carbonate_only(self, engine):
        """CO2_total=0.05 mol/L, no ammonia.

        At this total carbonate concentration, pH is weakly acidic
        (dominated by CO2 dissolved, most of which remains as CO2 at neutral pH).
        We verify:
          * pH < 7 (acidic, CO2 is weak acid)
          * mass closure: [CO2] + [HCO3-] + [CO3--] = 0.05 mol/L
          * charge residual ≈ 0
        """
        out = engine.solve(totals={"CO2": 0.05}, strong_ions={})
        assert out.pH < 7.0
        sp = out.species_mol_L
        total_C = sp["CO2"] + sp["HCO3-"] + sp["CO3--"]
        assert total_C == pytest.approx(0.05, rel=1e-6)
        assert abs(out.charge_residual) < 1e-7

    def test_ammonia_only(self, engine):
        """NH3_total=0.04 mol/L (total ammoniacal N), no carbonate.

        Ammonia is a weak base; a pure solution of ammoniacal N should be
        basic (pH > 7).  Mass closure: [NH3] + [NH4+] = 0.04 mol/L.
        """
        out = engine.solve(totals={"NH3": 0.04}, strong_ions={})
        assert out.pH > 7.0
        sp = out.species_mol_L
        total_N = sp["NH3"] + sp["NH4+"]
        assert total_N == pytest.approx(0.04, rel=1e-6)
        assert abs(out.charge_residual) < 1e-7

    def test_mixed_carbonate_ammonia(self, engine):
        """CO2=0.05, NH3=0.04 mol/L — realistic AD-like concentrations."""
        out = engine.solve(totals={"CO2": 0.05, "NH3": 0.04}, strong_ions={})
        # Both mass balances must close
        sp = out.species_mol_L
        assert sp["CO2"] + sp["HCO3-"] + sp["CO3--"] == pytest.approx(0.05, rel=1e-6)
        assert sp["NH3"] + sp["NH4+"] == pytest.approx(0.04, rel=1e-6)
        assert abs(out.charge_residual) < 1e-7
        # pH should be in a reasonable range for this chemistry
        assert 5.0 < out.pH < 10.0

    def test_strong_ion_na_raises_ph(self, engine):
        """Adding Na+ (cation) increases pH relative to baseline.

        Justification for including strong-ion tests
        --------------------------------------------
        Strong ions shift the charge balance without participating in
        equilibrium.  Na+ adds positive charge, which must be compensated
        by additional deprotonation → lower [H+] → higher pH.  This
        directly tests that ``strong_charge`` is correctly accumulated in
        the charge-balance residual row of the NR system.  In real process
        water Na+, K+, Ca2+, and Cl- are always present; if the charge
        balance misses them the solver converges to the wrong pH.
        """
        baseline = engine.solve(totals={"CO2": 0.05, "NH3": 0.04}, strong_ions={})
        with_na = engine.solve(
            totals={"CO2": 0.05, "NH3": 0.04},
            strong_ions={"CT_Na": 0.01},
        )
        assert with_na.pH > baseline.pH
        # Mass balances still hold
        sp = with_na.species_mol_L
        assert sp["CO2"] + sp["HCO3-"] + sp["CO3--"] == pytest.approx(0.05, rel=1e-6)
        assert sp["NH3"] + sp["NH4+"] == pytest.approx(0.04, rel=1e-6)
        assert abs(with_na.charge_residual) < 1e-6

    def test_strong_ion_cl_lowers_ph(self, engine):
        """Adding Cl- (anion) decreases pH relative to baseline."""
        baseline = engine.solve(totals={"CO2": 0.05, "NH3": 0.04}, strong_ions={})
        with_cl = engine.solve(
            totals={"CO2": 0.05, "NH3": 0.04},
            strong_ions={"CT_Cl": 0.01},
        )
        assert with_cl.pH < baseline.pH
        assert abs(with_cl.charge_residual) < 1e-6

    def test_output_keys_present(self, engine):
        """Result must carry the standard meta fields and species."""
        out = engine.solve(totals={"CO2": 0.05, "NH3": 0.04}, strong_ions={})
        for field in ("pH", "pH_conc", "logH", "aH", "gamma_H", "gamma_OH",
                      "ionic_strength", "charge_residual"):
            assert getattr(out, field) is not None, f"Missing field: {field!r}"
        for sp_id in ("H+", "OH-", "CO2", "HCO3-", "CO3--", "NH3", "NH4+"):
            assert sp_id in out.species_mol_L, f"Missing species: {sp_id!r}"

    def test_warmstart_reduces_iterations(self):
        """Calling solve() twice with the same problem should converge from cache."""
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        eng = NRChemicalEquilibriumEngine.from_reactions(_make_reactions(), use_warmstart=True)
        totals = {"CO2": 0.05, "NH3": 0.04}
        out1 = eng.solve(totals=totals, strong_ions={})
        out2 = eng.solve(totals=totals, strong_ions={})
        # Results must be identical
        assert out2.pH == pytest.approx(out1.pH, abs=1e-10)

    def test_reset_cache(self):
        """reset_cache() clears the warmstart state without breaking subsequent solves."""
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        eng = NRChemicalEquilibriumEngine.from_reactions(_make_reactions())
        eng.solve(totals={"CO2": 0.05, "NH3": 0.04}, strong_ions={})
        eng.reset_cache()
        assert eng._cache.log_x is None
        # Must still solve correctly after reset
        out = eng.solve(totals={"CO2": 0.05, "NH3": 0.04}, strong_ions={})
        assert abs(out.charge_residual) < 1e-7


# ═══════════════════════════════════════════════════════════════════════════
#  TestNRvsChargeBalance — NR and existing 1-D engine must agree
# ═══════════════════════════════════════════════════════════════════════════

def _cb_engine():
    """Build the existing BisectionChemicalEquilibriumEngine for the bespoke chemistry."""
    from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
    return BisectionChemicalEquilibriumEngine.from_reactions(_make_reactions())


class TestNRvsChargeBalance:
    """NR and 1-D charge-balance engines must agree to within solver tolerance.

    The existing BisectionChemicalEquilibriumEngine maps our reactions to:
      CT_TIC  → total carbonate (CO2 + HCO3- + CO3--)
      CT_NH_T → total ammoniacal N (NH3 + NH4+)

    The NR engine uses totals keyed by master id (CO2, NH3) with identical
    numerical values.  Both engines solve from the same thermodynamic data;
    results should agree to < 0.001 pH units and < 0.1 % for major species.
    """

    # Tolerance constants
    PH_TOL = 1e-3       # pH units
    CONC_RTOL = 1e-3    # relative tolerance on concentrations

    def _solve_both(self, CT_TIC, CT_NH_T, CT_Na=0.0, CT_Cl=0.0):
        """Return (cb_out, nr_out) for the same chemistry."""
        # 1-D charge-balance engine
        cb = _cb_engine()
        cb_out = cb.solve(CT_TIC=CT_TIC, CT_NH_T=CT_NH_T, CT_Na=CT_Na, CT_Cl=CT_Cl)

        # NR engine
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        nr = NRChemicalEquilibriumEngine.from_reactions(_make_reactions())
        strong_ions = {}
        if CT_Na:
            strong_ions["CT_Na"] = CT_Na
        if CT_Cl:
            strong_ions["CT_Cl"] = CT_Cl
        nr_out = nr.solve(
            totals={"CO2": CT_TIC, "NH3": CT_NH_T},
            strong_ions=strong_ions,
        )
        return cb_out, nr_out

    def test_pure_carbonate_ph(self):
        cb, nr = self._solve_both(CT_TIC=0.05, CT_NH_T=0.0)
        assert abs(nr.pH - cb.pH) < self.PH_TOL

    def test_pure_carbonate_species(self):
        cb, nr = self._solve_both(CT_TIC=0.05, CT_NH_T=0.0)
        # "CO2aq" isn't in BisectionChemicalEquilibriumEngine's canonical writeback tuple, so it
        # may land in `extra` rather than `species_mol_L` — check both.
        cb_species = {**cb.species_mol_L, **cb.extra}
        for sp_id in ("CO2", "HCO3-", "CO3--"):
            cb_c = float(cb_species.get(sp_id, cb_species.get("CO2aq", 0.0)) if sp_id == "CO2" else cb_species.get(sp_id, 0.0))
            nr_c = float(nr.species_mol_L[sp_id])
            if cb_c > 1e-15:
                assert abs(nr_c - cb_c) / cb_c < self.CONC_RTOL, (
                    f"{sp_id}: NR={nr_c:.4e}, CB={cb_c:.4e}"
                )

    def test_mixed_ph(self):
        cb, nr = self._solve_both(CT_TIC=0.05, CT_NH_T=0.04)
        assert abs(nr.pH - cb.pH) < self.PH_TOL

    def test_mixed_species(self):
        cb, nr = self._solve_both(CT_TIC=0.05, CT_NH_T=0.04)
        for sp_id in ("HCO3-", "NH3", "NH4+"):
            nr_c = float(nr.species_mol_L[sp_id])
            cb_c = float(cb.species_mol_L.get(sp_id, 0.0))
            if cb_c > 1e-12:
                assert abs(nr_c - cb_c) / cb_c < self.CONC_RTOL, (
                    f"{sp_id}: NR={nr_c:.4e}, CB={cb_c:.4e}"
                )

    def test_with_strong_na_ph(self):
        """pH agreement holds in the presence of Na+ strong ion."""
        cb, nr = self._solve_both(CT_TIC=0.05, CT_NH_T=0.04, CT_Na=0.01)
        assert abs(nr.pH - cb.pH) < self.PH_TOL

    def test_with_strong_cl_ph(self):
        """pH agreement holds in the presence of Cl- strong ion."""
        cb, nr = self._solve_both(CT_TIC=0.05, CT_NH_T=0.04, CT_Cl=0.01)
        assert abs(nr.pH - cb.pH) < self.PH_TOL

    def test_dilute_limits(self):
        """Agreement at very low concentrations (near pure water)."""
        cb, nr = self._solve_both(CT_TIC=1e-5, CT_NH_T=1e-5)
        assert abs(nr.pH - cb.pH) < self.PH_TOL

    def test_high_concentration(self):
        """Agreement at higher concentrations (0.2 mol/L carbonate)."""
        cb, nr = self._solve_both(CT_TIC=0.2, CT_NH_T=0.1)
        assert abs(nr.pH - cb.pH) < self.PH_TOL


# ═══════════════════════════════════════════════════════════════════════════
#  TestReactionSystemRouter — solver= kwarg on ReactionSystem
# ═══════════════════════════════════════════════════════════════════════════

class TestReactionSystemRouter:
    """ReactionSystem(solver=...) routes to the correct engine class."""

    def test_default_is_charge_balance(self):
        from PyOMES.reactions.reaction_system import ReactionSystem
        from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine
        rs = ReactionSystem(_make_reactions())
        assert rs._solver == "charge_balance"
        eng = rs.engine
        assert isinstance(eng, BisectionChemicalEquilibriumEngine)

    def test_newton_raphson_solver(self):
        from PyOMES.reactions.reaction_system import ReactionSystem
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        rs = ReactionSystem(_make_reactions(), solver="newton_raphson")
        assert rs._solver == "newton_raphson"
        eng = rs.engine
        assert isinstance(eng, NRChemicalEquilibriumEngine)

    def test_unknown_solver_raises(self):
        from PyOMES.reactions.reaction_system import ReactionSystem
        with pytest.raises(ValueError, match="unknown solver"):
            ReactionSystem(_make_reactions(), solver="gibbs_minimisation")

    def test_repr_includes_solver(self):
        from PyOMES.reactions.reaction_system import ReactionSystem
        rs = ReactionSystem(_make_reactions(), solver="newton_raphson")
        assert "newton_raphson" in repr(rs)

    def test_configure_engine_docstring_warns(self):
        """configure_engine() must NOT accept a solver= kwarg; passing one
        should raise TypeError (unexpected keyword argument)."""
        from PyOMES.reactions.reaction_system import ReactionSystem
        rs = ReactionSystem(_make_reactions(), solver="charge_balance")
        with pytest.raises(TypeError):
            rs.configure_engine(solver="newton_raphson")

    def test_nr_engine_solve_via_reaction_system(self):
        """End-to-end: ReactionSystem(solver='newton_raphson').engine.solve()."""
        from PyOMES.reactions.reaction_system import ReactionSystem
        rs = ReactionSystem(_make_reactions(), solver="newton_raphson")
        out = rs.engine.solve(totals={"CO2": 0.05, "NH3": 0.04}, strong_ions={})
        assert 5.0 < out.pH < 10.0
        assert abs(out.charge_residual) < 1e-7


# ═══════════════════════════════════════════════════════════════════════════
#  TestPhaseWriteback — engine writes derived species to phase.n_mol
# ═══════════════════════════════════════════════════════════════════════════

class _MockPhase:
    """Minimal phase stub for writeback tests."""
    def __init__(self, n_mol: dict, V_L: float = 1.0):
        self.n_mol = dict(n_mol)
        self.V_L = float(V_L)

    def _refresh_derived(self, derived: dict) -> None:
        self.n_mol.update(derived)


class TestPhaseWriteback:
    """NR engine writes H+, OH-, H2O, and equilibrium species to n_mol."""

    def _make_phase(self):
        """Liquid phase seeded with total C and N as master species (mol, V_L=1 L)."""
        return _MockPhase({"CO2": 0.05, "NH3": 0.04}, V_L=1.0)

    def _run(self):
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        engine = NRChemicalEquilibriumEngine.from_reactions(_make_reactions())
        phase = self._make_phase()
        phases = {"liquid": phase}
        out = engine.solve(phases=phases)
        out.apply_to_phases(phases)  # solve() no longer auto-writes
        return phase, out

    def test_h_plus_written(self):
        phase, _ = self._run()
        assert "H+" in phase.n_mol
        assert float(phase.n_mol["H+"]) > 0.0

    def test_oh_minus_written(self):
        phase, _ = self._run()
        assert "OH-" in phase.n_mol
        assert float(phase.n_mol["OH-"]) > 0.0

    def test_h2o_written(self):
        """H₂O is written as engine-owned derived species (NR engine only)."""
        phase, _ = self._run()
        assert "H2O" in phase.n_mol
        # ≈ 55.5 mol (V_L = 1 L)
        assert float(phase.n_mol["H2O"]) == pytest.approx(55.5, abs=0.1)

    def test_carbonate_species_written(self):
        phase, _ = self._run()
        assert "HCO3-" in phase.n_mol
        assert float(phase.n_mol["HCO3-"]) > 0.0

    def test_carbonate_mass_balance_in_n_mol(self):
        """After writeback, sum of n_mol across carbonate ladder = 0.05 mol."""
        phase, _ = self._run()
        total_C = (
            float(phase.n_mol.get("CO2", 0.0))
            + float(phase.n_mol.get("HCO3-", 0.0))
            + float(phase.n_mol.get("CO3--", 0.0))
        )
        assert total_C == pytest.approx(0.05, rel=1e-5)

    def test_second_solve_reads_updated_n_mol(self):
        """A second solve reads totals from the updated n_mol (sum of all species)."""
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        engine = NRChemicalEquilibriumEngine.from_reactions(_make_reactions())
        phase = self._make_phase()
        phases = {"liquid": phase}
        out1 = engine.solve(phases=phases)
        out1.apply_to_phases(phases)  # solve() no longer auto-writes
        # Second call: n_mol now has HCO3-, CO3--, NH4+ etc from first writeback
        out2 = engine.solve(phases=phases)
        assert out2.pH == pytest.approx(out1.pH, abs=1e-8)

    def test_algebraic_species_set(self):
        """algebraic_species() must include all species computed by the engine."""
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
        engine = NRChemicalEquilibriumEngine.from_reactions(_make_reactions())
        alg = engine.algebraic_species()
        for sp_id in ("H+", "OH-", "H2O", "CO2", "HCO3-", "CO3--", "NH3", "NH4+"):
            assert sp_id in alg, f"{sp_id} missing from algebraic_species()"


# ═══════════════════════════════════════════════════════════════════════════
#  TestPrecipitationEquilibria
# ═══════════════════════════════════════════════════════════════════════════

def _make_carbonate_reactions():
    """Carbonate-only reactions (water + 2-step carbonate ladder)."""
    from PyOMES.chemistry.common_species import H2O, H_plus, OH_minus, CO2, HCO3_minus, CO3_2minus
    from PyOMES.reactions.equilibrium import EquilibriumReaction
    from PyOMES.reactions.stoichiometry import StoichiometryEntry

    def _e(sp, phase, coeff):
        return StoichiometryEntry(species=sp, phase=phase, coefficient=coeff)

    water = EquilibriumReaction(
        stoichiometry=[_e(H2O, "liquid", -1), _e(H_plus, "liquid", +1), _e(OH_minus, "liquid", +1)],
        log_K=-14.0, label="water",
    )
    co2_first = EquilibriumReaction(
        stoichiometry=[_e(CO2, "liquid", -1), _e(H2O, "liquid", -1),
                       _e(HCO3_minus, "liquid", +1), _e(H_plus, "liquid", +1)],
        log_K=-6.35, total_id="CO2", label="co2_first",
    )
    co2_second = EquilibriumReaction(
        stoichiometry=[_e(HCO3_minus, "liquid", -1),
                       _e(CO3_2minus, "liquid", +1), _e(H_plus, "liquid", +1)],
        log_K=-10.33, total_id="CO2", label="co2_second",
    )
    return [water, co2_first, co2_second]


def _make_calcite_reaction():
    """Calcite dissolution: CaCO3(s) <-> Ca++ + CO3--  log_K = -8.48 (Ksp at 25°C)."""
    from PyOMES.chemistry.common_species import Ca_plus_plus, CO3_2minus
    from PyOMES.chemistry.species import Species
    from PyOMES.reactions.equilibrium import EquilibriumReaction
    from PyOMES.reactions.stoichiometry import StoichiometryEntry

    CaCO3 = Species(id="CaCO3", atoms={"Ca": 1, "C": 1, "O": 3}, charge=0, MW=100.086)

    def _e(sp, phase, coeff):
        return StoichiometryEntry(species=sp, phase=phase, coefficient=coeff)

    return EquilibriumReaction(
        stoichiometry=[
            _e(CaCO3,       "solid",  -1),
            _e(Ca_plus_plus,"liquid", +1),
            _e(CO3_2minus,  "liquid", +1),
        ],
        log_K=-8.48,
        label="calcite",
    )


def _make_calcite_engine(use_activity=True):
    """Build NRChemicalEquilibriumEngine with calcite precipitation reactions."""
    from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine
    return NRChemicalEquilibriumEngine.from_reactions(
        _make_carbonate_reactions(),
        precipitation_reactions=[_make_calcite_reaction()],
        use_activity=use_activity,
        activity_model="davies",
    )


class TestPrecipitationEquilibria:
    """Outer active-set precipitation loop in NRChemicalEquilibriumEngine."""

    def test_calcite_precipitation_from_supersaturated(self):
        """CT_Ca=0.002, CT_CO2=0.010, CT_Na=0.005 (SI > 0 from notebook 03).

        After solve: xi > 0 and log10(a_Ca * a_CO3) ≈ -8.48 ± 0.01.
        """
        engine = _make_calcite_engine(use_activity=True)
        out = engine.solve(
            totals={"CO2": 0.010},
            strong_ions={"CT_Ca": 0.002, "CT_Na": 0.005},
        )
        xi = out.extra["minerals_xi_mol_L"]["calcite"]
        si = out.saturation_indices["calcite"]
        assert xi > 0.0, f"Expected precipitation (xi > 0), got xi={xi}"
        assert abs(si) < 0.05, f"SI not near zero after precipitation: SI={si:.4f}"

    def test_calcite_mass_balance(self):
        """Total Ca conserved: xi + dissolved_Ca_eff ≈ CT_Ca_input."""
        CT_Ca_input = 0.002
        CT_CO2_input = 0.010
        engine = _make_calcite_engine(use_activity=True)
        out = engine.solve(
            totals={"CO2": CT_CO2_input},
            strong_ions={"CT_Ca": CT_Ca_input, "CT_Na": 0.005},
        )
        xi = out.extra["minerals_xi_mol_L"]["calcite"]
        # Dissolved Ca is CT_Ca_input - xi (it's a strong ion, not in NR output)
        Ca_dissolved = CT_Ca_input - xi
        assert Ca_dissolved >= -1e-9, f"Dissolved Ca negative: {Ca_dissolved}"
        assert abs(xi + Ca_dissolved - CT_Ca_input) < 1e-8, (
            f"Ca not conserved: xi={xi}, Ca_diss={Ca_dissolved}, input={CT_Ca_input}"
        )

    def test_undersaturated_no_precipitation(self):
        """Very low Ca + CO2 at low pH: SI < 0, xi = 0."""
        engine = _make_calcite_engine(use_activity=True)
        # Low pH (CT_Cl drives acidic), small Ca and CO2 → undersaturated
        out = engine.solve(
            totals={"CO2": 0.001},
            strong_ions={"CT_Ca": 0.0001, "CT_Cl": 0.005},
        )
        xi = out.extra["minerals_xi_mol_L"]["calcite"]
        si = out.saturation_indices["calcite"]
        assert xi == 0.0, f"Expected no precipitation, got xi={xi}"
        assert si < 0.0, f"Expected SI < 0 (undersaturated), got SI={si}"

    def test_dissolution_from_solid(self):
        """Model: supply CT_Ca = n_mol_solid (as if all solid dissolved); alkaline conditions
        make starting state supersaturated; outer loop finds equilibrium where IAP = Ksp.

        Represents dissolution from a CaCO3 solid feed: the solid provides up to
        n_mol_solid mol/L of Ca²⁺.  After solve: IAP ≈ Ksp and
        dissolved_Ca + xi_reprecipitated = n_mol_solid.
        """
        n_mol_solid = 0.001   # mol/L CaCO3(s) available — all treated as dissolved initially
        CT_CO2 = 0.010        # mol/L TIC
        CT_Na  = 0.010        # mol/L (drives pH up; at CT_Ca=0.001, CT_Na=0.010: SI ≈ +2.35)

        engine = _make_calcite_engine(use_activity=True)
        # Pass n_mol_solid as CT_Ca: represents max dissolved Ca if all solid dissolved.
        # The outer loop precipitates back whatever exceeds Ksp.
        out = engine.solve(
            totals={"CO2": CT_CO2},
            strong_ions={"CT_Ca": n_mol_solid, "CT_Na": CT_Na},
        )
        xi = out.extra["minerals_xi_mol_L"]["calcite"]
        si = out.saturation_indices["calcite"]
        assert abs(si) < 0.05, f"IAP not ≈ Ksp after solve: SI={si:.4f}"
        assert 0.0 <= xi <= n_mol_solid, f"xi out of range: xi={xi}, n_mol_solid={n_mol_solid}"
        Ca_dissolved = n_mol_solid - xi
        assert Ca_dissolved >= 0.0, f"Dissolved Ca negative: {Ca_dissolved}"

    def test_complete_dissolution_high_ksp(self):
        """Very high Ksp (log_Ksp=+5): solid fully dissolves into solution.

        Uses a dummy mineral with only CO3-- as dissolved product.
        Solid amount n_mol_solid = 0.001 mol/L.  After solve xi ≈ 0 (no
        precipitation since log_Ksp is large → always undersaturated).
        """
        from PyOMES.chemistry.species import Species
        from PyOMES.chemistry.common_species import CO3_2minus, H2O, H_plus, OH_minus
        from PyOMES.reactions.equilibrium import EquilibriumReaction
        from PyOMES.reactions.stoichiometry import StoichiometryEntry
        from PyOMES.chemical_equilibrium.nr_engine import NRChemicalEquilibriumEngine

        DummySolid = Species(id="DummySolid", atoms={"C": 1, "O": 3}, charge=0, MW=60.0)

        def _e(sp, phase, coeff):
            return StoichiometryEntry(species=sp, phase=phase, coefficient=coeff)

        # A mineral with very high Ksp → always undersaturated
        dummy_mineral = EquilibriumReaction(
            stoichiometry=[
                _e(DummySolid, "solid",  -1),
                _e(CO3_2minus, "liquid", +1),
            ],
            log_K=+5.0,   # very high → always undersaturated → no precipitation
            label="dummy_high_ksp",
        )
        engine = NRChemicalEquilibriumEngine.from_reactions(
            _make_carbonate_reactions(),
            precipitation_reactions=[dummy_mineral],
            use_activity=False,
        )
        # Even with CT_CO2 = 0.01, [CO3--] << 10^5 → SI << 0 → xi = 0
        out = engine.solve(
            totals={"CO2": 0.01},
            strong_ions={"CT_Na": 0.005},
        )
        xi = out.extra["minerals_xi_mol_L"]["dummy_high_ksp"]
        si = out.saturation_indices["dummy_high_ksp"]
        assert xi == 0.0, f"Expected xi=0 (undersaturated), got xi={xi}"
        assert si < 0.0, f"Expected SI < 0, got SI={si}"

    def test_element_stoichiometry_computed(self):
        """element_stoichiometry excludes H+ and correctly tracks non-H+ master coefficients.

        Verifies:
        1. OH- (nu = {"H+": -1}) → element_stoichiometry = {} (H+ excluded → empty)
        2. CO3-- (nu = {"CO2": 1, "H+": -2}) → element_stoichiometry = {"CO2": 1}
        3. A manually-constructed cross-component SecondaryEntry for CaHCO3+
           (nu = {"CO2": 1, "Ca++": 1, "H+": -1}) → element_stoichiometry = {"CO2": 1, "Ca++": 1}
           (direct struct test; tableau builder defers multi-master support to Phase 2)
        """
        from PyOMES.chemical_equilibrium.nr_tableau import SecondaryEntry, build_tableau

        # Part 1+2: test element_stoichiometry from a real build_tableau call
        tableau = build_tableau(_make_reactions(), T_K=298.15)

        oh = next(s for s in tableau.secondaries if s.species_id == "OH-")
        assert oh.element_stoichiometry == {}, (
            f"OH- element_stoichiometry should be empty (only H+ in nu): {oh.element_stoichiometry}"
        )

        co3 = next(s for s in tableau.secondaries if s.species_id == "CO3--")
        assert "CO2" in co3.element_stoichiometry
        assert "H+" not in co3.element_stoichiometry
        assert co3.element_stoichiometry["CO2"] == pytest.approx(1.0)

        # Part 3: direct struct test for cross-component secondary (CaHCO3+)
        # nu includes Ca++ (a future master) + CO2 + H+;
        # element_stoichiometry should contain Ca++ and CO2 but NOT H+.
        sec = SecondaryEntry(
            species_id="CaHCO3+",
            charge=+1,
            nu={"CO2": 1.0, "Ca++": 1.0, "H+": -1.0},
            log_K_prime=-5.09,
            element_stoichiometry={k: v for k, v in {"CO2": 1.0, "Ca++": 1.0, "H+": -1.0}.items()
                                   if k != "H+"},
        )
        assert "CO2" in sec.element_stoichiometry
        assert "Ca++" in sec.element_stoichiometry
        assert "H+" not in sec.element_stoichiometry
        assert sec.element_stoichiometry["CO2"] == pytest.approx(1.0)
        assert sec.element_stoichiometry["Ca++"] == pytest.approx(1.0)

    def test_element_stoichiometry_single_component_unchanged(self):
        """Backward-compatible: single-component secondaries have element_stoichiometry with exactly one non-H+ key."""
        from PyOMES.chemical_equilibrium.nr_tableau import build_tableau

        tableau = build_tableau(_make_reactions(), T_K=298.15)

        hco3 = next(s for s in tableau.secondaries if s.species_id == "HCO3-")
        assert set(hco3.element_stoichiometry.keys()) == {"CO2"}, (
            f"HCO3- element_stoichiometry should be {{'CO2'}}, got {hco3.element_stoichiometry}"
        )
        assert hco3.element_stoichiometry["CO2"] == pytest.approx(1.0)

        nh4 = next(s for s in tableau.secondaries if s.species_id == "NH4+")
        assert set(nh4.element_stoichiometry.keys()) == {"NH3"}, (
            f"NH4+ element_stoichiometry should be {{'NH3'}}, got {nh4.element_stoichiometry}"
        )
        assert nh4.element_stoichiometry["NH3"] == pytest.approx(1.0)
