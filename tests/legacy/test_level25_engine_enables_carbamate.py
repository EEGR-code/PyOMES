# -*- coding: utf-8 -*-
"""
Unit tests for engine-level Level-2.5 injection.

This file verifies:
  (A) level=2.5 (or enable_level25=True) injects enable_carbamate=True and carbamate_eq
  (B) level=2 alone does NOT automatically enable carbamate
"""

import os
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(THIS_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from speciation.engine import BisectionChemicalEquilibriumEngine


def _minimal_totals():
    return dict(
        CT_TIC=0.05,
        CT_NH_T=0.05,
        CT_P=0.0,
        acid_totals={},
        CT_Na=0.0, CT_K=0.0, CT_Cl=0.0, CT_NO3=0.0,
        CT_Ca=0.0, CT_Mg=0.0, CT_Zn=0.0, CT_Mn=0.0, CT_Co=0.0, CT_Mo7O24=0.0,
    )


def test_engine_level25_enables_carbamate():
    eng = BisectionChemicalEquilibriumEngine(
        level=2.5,                 # <-- MUST be 2.5 for this test
        use_activity=True,
        activity_model="davies",
        T_C=25.0,
        carbamate_eq={"logK_25": 0.0, "K_basis": "activity", "source": "test"},
    )

    sol = eng.solve(**_minimal_totals())

    dbg = sol.get("carbamate_debug", {})
    assert dbg.get("enabled") is True, (
        "Expected carbamate enabled when BisectionChemicalEquilibriumEngine(level=2.5). "
        f"Got carbamate_debug={dbg}"
    )
    assert "NH2COO-" in sol, "Expected NH2COO- key present when carbamate is enabled."


def test_engine_level2_does_not_enable_carbamate_by_default():
    eng = BisectionChemicalEquilibriumEngine(
        level=2,                   # <-- Level 2 only
        use_activity=True,
        activity_model="davies",
        T_C=25.0,
        carbamate_eq={"logK_25": 0.0, "K_basis": "activity", "source": "test"},
    )

    sol = eng.solve(**_minimal_totals())

    dbg = sol.get("carbamate_debug", {})
    # In pure level=2 runs, carbamate should not be enabled unless caller explicitly requests it.
    assert dbg.get("enabled") in (False, None), (
        "Expected carbamate disabled by default when BisectionChemicalEquilibriumEngine(level=2). "
        f"Got carbamate_debug={dbg}"
    )


if __name__ == "__main__":
    print("Running Level 2.5 engine injection tests...")
    test_engine_level25_enables_carbamate()
    print("PASS: level=2.5 enables carbamate")
    test_engine_level2_does_not_enable_carbamate_by_default()
    print("PASS: level=2 does not enable carbamate by default")
    print("All tests passed ✅")
