# -*- coding: utf-8 -*-
"""
Unit test: Level 2 phosphate complexation is totals-consistent.

What it checks:
1) In a Ca2+ + phosphate solution, Level 2 forms phosphate complexes (e.g., CaHPO4(aq), CaPO4-).
2) Totals-consistent patch produces diagnostics:
   - CT_P_total_input
   - CT_P_bound_complexes
   - CT_P_free_pool_used
   and ensures free + bound ≈ total.

This test supports two modes:
A) Recipe-based phosphate (preferred): uses a phosphate salt that sets CT_P via CHEM_DB (e.g., Na2HPO4, KH2PO4, NaH2PO4).
B) Fallback: if no phosphate salts are available, it bypasses recipe_to_totals for phosphate and sets CT_P directly.

Run with:
- pytest -q
or run this file directly in Spyder / Python.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG_ROOT = os.path.join(HERE, "v8_split_17")  # adjust if your package folder differs
if os.path.isdir(PKG_ROOT) and PKG_ROOT not in sys.path:
    sys.path.insert(0, PKG_ROOT)

from scripts.chem_recipe import recipe_to_totals
from speciation.engine import BisectionChemicalEquilibriumEngine


def _solve_from_recipe(recipe_mol_L, T_C=25.0, use_activity=True):
    totals = recipe_to_totals(recipe_mol_L, acid_pKas={})
    eng = BisectionChemicalEquilibriumEngine(T_C=float(T_C))
    sol = eng.solve(
        acid_totals=totals.get("acid_totals", {}) or {},
        acid_pKas=totals.get("acid_pKas", {}) or {},
        CT_TIC=float(totals.get("CT_TIC", 0.0)),
        CT_P=float(totals.get("CT_P", 0.0)),
        CT_NH_T=float(totals.get("CT_NH_T", 0.0)),
        use_activity=bool(use_activity),
        activity_model="davies",
        **(totals.get("strong_ions", {}) or {}),
    )
    return sol, totals


def _solve_with_manual_CT_P(*, CT_P, strong_ions, T_C=25.0, use_activity=True):
    """
    Fallback solver if you cannot introduce phosphate through CHEM_DB recipes.
    This still validates the *mass-balance coupling logic* for phosphate.
    """
    eng = BisectionChemicalEquilibriumEngine(T_C=float(T_C))
    sol = eng.solve(
        acid_totals={},
        acid_pKas={},
        CT_TIC=0.0,
        CT_P=float(CT_P),
        CT_NH_T=0.0,
        use_activity=bool(use_activity),
        activity_model="davies",
        **(strong_ions or {}),
    )
    return sol


def _has_any_key(sol, keys):
    return any((k in sol) and (float(sol.get(k, 0.0)) > 1e-12) for k in keys)


def test_level2_phosphate_complexes_form_and_mass_balance():
    """
    Primary test:
    - uses CaCl2 + phosphate (via recipe if possible)
    - asserts phosphate complexes form
    - asserts CT_P mass balance diagnostics exist and are consistent
    """
    # Preferred: use phosphate salts from CHEM_DB
    # Try a few likely keys depending on your CHEM_DB setup
    phosphate_candidates = [
        {"Na2HPO4": 0.010},       # 10 mM dibasic sodium phosphate
        {"KH2PO4": 0.010},        # 10 mM monobasic potassium phosphate
        {"NaH2PO4": 0.010},       # if you added this to CHEM_DB
        {"Na2HPO4·2H2O": 0.010},  # if you added this to CHEM_DB
    ]

    recipe_base = {
        "CaCl2": 0.010,  # 10 mM Ca2+
        "NaCl":  0.100,  # background ionic strength
    }

    sol = None
    totals = None
    used_recipe = None

    for phos in phosphate_candidates:
        recipe = dict(recipe_base)
        recipe.update(phos)
        try:
            sol, totals = _solve_from_recipe(recipe, T_C=25.0, use_activity=True)
            used_recipe = recipe
            # Require that CT_P is actually nonzero from recipe_to_totals
            if float(totals.get("CT_P", 0.0)) > 1e-9:
                break
            sol = None
        except KeyError:
            sol = None

    # Fallback: manual CT_P if no phosphate salt recipe works
    if sol is None:
        # 10 mM phosphate total, with CaCl2 and NaCl as strong ions
        strong_ions = {"CT_Ca": 0.010, "CT_Cl": 0.020 + 0.100, "CT_Na": 0.100}  # approx: CaCl2 adds 2 Cl-
        sol = _solve_with_manual_CT_P(CT_P=0.010, strong_ions=strong_ions, T_C=25.0, use_activity=True)

    # 1) Complexes should form (depending on your association list / naming)
    phosphate_complex_keys = [
        "CaH2PO4+",
        "CaHPO4(aq)",
        "CaPO4-",
        "MgHPO4(aq)",
        "MgPO4-",
    ]

    assert _has_any_key(sol, phosphate_complex_keys), (
        "Expected phosphate complexes (e.g., CaHPO4(aq), CaPO4-) to form, but none were found. "
        "This may mean:\n"
        "- phosphate associations are not in LEVEL2_ASSOCIATIONS, or\n"
        "- complex names differ from the expected keys in this test, or\n"
        "- phosphate is not actually entering CT_P.\n"
        f"Recipe used (if any): {used_recipe}"
    )

    # 2) Totals-consistent diagnostics must exist
    required = ["CT_P_total_input", "CT_P_bound_complexes", "CT_P_free_pool_used"]
    missing = [k for k in required if k not in sol]
    assert not missing, (
        "Totals-consistent phosphate diagnostics missing: "
        f"{missing}. This usually means bound ligand accounting/coupling for phosphate "
        "is not wired (bound_ligands_from_complexes doesn't count phosphate ligands, "
        "or chemistry_level2 coupling loop isn't updating CT_P)."
    )

    CT_tot = float(sol["CT_P_total_input"])
    CT_bnd = float(sol["CT_P_bound_complexes"])
    CT_free = float(sol["CT_P_free_pool_used"])

    assert CT_tot > 0.0, "CT_P_total_input should be > 0 for this test."
    assert CT_bnd > 0.0, (
        "Expected some phosphate to be bound into Ca/Mg phosphate complexes, but CT_P_bound_complexes is 0. "
        "Either complexation isn't active, or phosphate ligands aren't being counted."
    )
    assert CT_free < CT_tot, (
        "Expected free phosphate pool to be smaller than total phosphate after complexation. "
        f"Got CT_P_free_pool_used={CT_free}, CT_P_total_input={CT_tot}."
    )

    # Tight mass balance check (free + bound approx total)
    assert abs((CT_free + CT_bnd) - CT_tot) < 1e-6, (
        "Phosphate mass balance check failed: free + bound should ~ equal total. "
        f"free={CT_free}, bound={CT_bnd}, total={CT_tot}"
    )


def test_no_calcium_means_no_phosphate_bound():
    """
    Control: without Ca2+/Mg2+, phosphate should not bind into metal-phosphate complexes.
    """
    # Try to introduce phosphate via recipe; fallback to manual.
    recipe = {"Na2HPO4": 0.010, "NaCl": 0.100}

    try:
        sol, totals = _solve_from_recipe(recipe, T_C=25.0, use_activity=True)
        if float(totals.get("CT_P", 0.0)) < 1e-9:
            raise KeyError("CT_P not set from recipe")
    except Exception:
        # manual: 10 mM phosphate, NaCl background
        strong_ions = {"CT_Na": 0.100, "CT_Cl": 0.100}
        sol = _solve_with_manual_CT_P(CT_P=0.010, strong_ions=strong_ions, T_C=25.0, use_activity=True)

    if "CT_P_bound_complexes" in sol:
        assert float(sol["CT_P_bound_complexes"]) < 1e-10, (
            f"Expected ~0 phosphate bound without Ca/Mg, but got CT_P_bound_complexes={sol['CT_P_bound_complexes']}."
        )


# --------------
# Allow running in Spyder without pytest
# --------------
if __name__ == "__main__":
    print("Running Level 2 phosphate complexation unit tests...\n")

    try:
        test_level2_phosphate_complexes_form_and_mass_balance()
        print("PASS: phosphate complexes form + mass balance holds")

        test_no_calcium_means_no_phosphate_bound()
        print("PASS: control case (no Ca) behaves")

        print("\nAll tests passed ✅")

    except AssertionError as e:
        print("\nTEST FAILED ❌")
        print(str(e))
        raise
