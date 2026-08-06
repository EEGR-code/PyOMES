# -*- coding: utf-8 -*-
"""
Unit test: Level 2 complexation is totals-consistent.

What it checks:
1) In a CaCl2 + NaHCO3 solution, Level 2 forms carbonate complexes (CaCO3(aq), CaHCO3+).
2) If the totals-consistent coupling patch is applied, the solver reports:
   - CT_TIC_bound_complexes > 0
   - CT_TIC_free_pool_used < CT_TIC_total_input

This test will FAIL with a helpful message if you haven't applied the patch.

Run with:
- pytest -q
or simply run this file in Spyder / Python.
"""

import os
import sys

# Ensure local package imports work if running from different working dirs
HERE = os.path.dirname(os.path.abspath(__file__))
PKG_ROOT = os.path.join(HERE, "v8_split_16")  # adjust if your package folder differs
if os.path.isdir(PKG_ROOT) and PKG_ROOT not in sys.path:
    sys.path.insert(0, PKG_ROOT)

from scripts.chem_recipe import recipe_to_totals
from speciation.engine import BisectionChemicalEquilibriumEngine


def _solve_case(recipe_mol_L, T_C=25.0, use_activity=True):
    totals = recipe_to_totals(recipe_mol_L, acid_pKas={})
    eng = BisectionChemicalEquilibriumEngine(T_C=float(T_C))
    sol = eng.solve(
        acid_totals=totals.get("acid_totals", {}) or {},
        acid_pKas=totals.get("acid_pKas", {}) or {},
        CT_TIC=float(totals.get("CT_TIC", 0.0)),
        CT_P=float(totals.get("CT_P", 0.0)),
        CT_NH_T=float(totals.get("CT_NH_T", 0.0)),
        use_activity=bool(use_activity),
        **(totals.get("strong_ions", {}) or {}),
    )
    return sol


def test_level2_carbonate_complexes_form():
    """
    Ca2+ + (bi)carbonate should produce CaHCO3+ and CaCO3(aq) at equilibrium.
    """
    recipe = {
        "NaHCO3": 0.050,  # 50 mM bicarbonate
        "CaCl2":  0.010,  # 10 mM calcium
        "NaCl":   0.100,  # background ionic strength
    }
    sol = _solve_case(recipe, T_C=25.0, use_activity=False)

    # Check complexes exist and are non-trivial
    ca_hco3 = float(sol.get("CaHCO3+", 0.0))
    ca_co3  = float(sol.get("CaCO3(aq)", 0.0))

    assert (ca_hco3 > 1e-8) or (ca_co3 > 1e-8), (
        "Expected carbonate complexes (CaHCO3+ and/or CaCO3(aq)) to form, "
        "but both were ~0. This may mean LEVEL2_ASSOCIATIONS not applied, "
        "or the Level 2 model isn't being used."
    )


def test_level2_complexation_is_totals_consistent():
    """
    This test specifically checks the totals-consistent patch:
    CT_TIC_total_input, CT_TIC_bound_complexes, CT_TIC_free_pool_used.
    """
    recipe = {
        "NaHCO3": 0.050,  # 50 mM bicarbonate
        "CaCl2":  0.010,  # 10 mM calcium
        "NaCl":   0.100,
    }
    sol = _solve_case(recipe, T_C=25.0, use_activity=False)

    required_keys = ["CT_TIC_total_input", "CT_TIC_bound_complexes", "CT_TIC_free_pool_used"]
    missing = [k for k in required_keys if k not in sol]

    assert not missing, (
        "Totals-consistent complexation diagnostics not found in solution output: "
        f"missing {missing}. \n\n"
        "This usually means you have NOT applied the totals-consistent coupling patch.\n"
        "Expected patch behaviour:\n"
        "- associations_level2.py: add bound_ligands_from_complexes(...)\n"
        "- chemistry_level2.py: wrap apply_level2_associations in an iteration that\n"
        "  subtracts bound ligands from CT_TIC/CT_P/acid_totals before recomputing ladders.\n"
    )

    CT_tot = float(sol["CT_TIC_total_input"])
    CT_bnd = float(sol["CT_TIC_bound_complexes"])
    CT_free = float(sol["CT_TIC_free_pool_used"])

    assert CT_tot > 0.0, "CT_TIC_total_input should be > 0 for this recipe."
    assert CT_bnd > 0.0, (
        "Expected some TIC to be bound into Ca/Mg carbonate complexes, but CT_TIC_bound_complexes is 0. "
        "Either complexation is not active or ligand accounting is not working."
    )
    assert CT_free < CT_tot, (
        "Expected free TIC pool to be smaller than total TIC after complexation. "
        f"Got CT_TIC_free_pool_used={CT_free} and CT_TIC_total_input={CT_tot}."
    )

    # Optional: tight mass-balance sanity (free + bound approx total)
    assert abs((CT_free + CT_bnd) - CT_tot) < 1e-6, (
        "TIC mass balance check failed: free + bound should approximately equal total. "
        f"free={CT_free}, bound={CT_bnd}, total={CT_tot}"
    )


def test_no_calcium_means_no_tic_bound():
    """
    Control: without Ca2+/Mg2+, carbonate complexes should not bind measurable TIC.
    """
    recipe = {
        "NaHCO3": 0.050,
        "NaCl":   0.100,
    }
    sol = _solve_case(recipe, T_C=25.0, use_activity=False)

    # If patch is applied, bound should be ~0
    if "CT_TIC_bound_complexes" in sol:
        CT_bnd = float(sol["CT_TIC_bound_complexes"])
        assert CT_bnd < 1e-10, (
            f"Expected ~0 TIC bound without Ca/Mg, but got CT_TIC_bound_complexes={CT_bnd}."
        )


# -----------------
# Allow running in Spyder without pytest
# -----------------
if __name__ == "__main__":
    print("Running Level 2 complexation unit tests...\n")

    try:
        test_level2_carbonate_complexes_form()
        print("PASS: carbonate complexes form")

        test_level2_complexation_is_totals_consistent()
        print("PASS: complexation is totals-consistent (TIC)")

        test_no_calcium_means_no_tic_bound()
        print("PASS: control case (no Ca) behaves")

        print("\nAll tests passed ✅")

    except AssertionError as e:
        print("\nTEST FAILED ❌")
        print(str(e))
        raise
