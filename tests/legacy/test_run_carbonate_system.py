# -*- coding: utf-8 -*-
"""
test_run_carbonate_system.py

Example carbonate-system validation tests against published pH values.

Cases:
1) 0.1 M NaHCO3 (freshly prepared) -> pH ~ 8.3 (25C)
2) NIST carbonate buffer: 0.025 M NaHCO3 + 0.025 M Na2CO3 -> pH ~ 10.02 (25C)

Run: python test_run_carbonate_system.py
"""

import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

from speciation.engine import BisectionChemicalEquilibriumEngine


def run_case(
    name: str,
    *,
    CT_TIC: float,
    CT_Na: float,
    expected_pH: float,
    tol: float = 0.15,
    level: int = 2,
    use_activity: bool = True,
    activity_model: str = "davies",
    T_C: float = 25.0,
):
    engine = BisectionChemicalEquilibriumEngine(
        level=level,
        use_activity=use_activity,
        activity_model=activity_model,
        T_C=T_C,
    )

    out = engine.solve(
        acid_totals={},
        acid_pKas={},
        CT_TIC=CT_TIC,
        CT_P=0.0,
        CT_NH_T=0.0,
        CT_Na=CT_Na,
        CT_Cl=0.0,
        CT_K=0.0,
        CT_NO3=0.0,
        CT_SO4=0.0,
        CT_Mg=0.0,
        CT_Ca=0.0,
        CT_Zn=0.0,
        CT_Mn=0.0,
        CT_Co=0.0,
        CT_Mo7O24=0.0,
        pH_min=0.0,
        pH_max=14.0,
        n_scan=600,
        tol=1e-12,
    )

    pH = float(out.get("pH", np.nan))
    I = float(out.get("IonicStrength", np.nan))

    ok = abs(pH - expected_pH) <= tol

    print(f"\n{name}")
    print(f"  Inputs: CT_TIC={CT_TIC:.6g} M, CT_Na={CT_Na:.6g} M")
    print(f"  Pred:   pH={pH:.4f}, I={I:.4g} M")
    print(f"  Target: pH≈{expected_pH:.2f} (±{tol:.2f})  ->  {'PASS' if ok else 'FAIL'}")

    # Useful debug species
    for k in ["CO2aq", "HCO3-", "CO3--", "H+", "OH-", "Na+", "aH"]:
        if k in out:
            print(f"    {k:6s}: {out[k]}")

    return ok, out


if __name__ == "__main__":
    # Case 1: 0.1 M NaHCO3
    # CT_TIC = 0.1 M total inorganic carbon (as bicarbonate salt)
    # CT_Na  = 0.1 M sodium counterion
    run_case(
        "Case 1: 0.1 M NaHCO3 (fresh)",
        CT_TIC=0.100,
        CT_Na=0.100,
        expected_pH=8.30,   # PubChem “freshly prepared 0.1 M solution: pH 8.3”
        tol=0.20,
    )

    # Case 2: NIST carbonate buffer standard:
    # 0.025 M NaHCO3 + 0.025 M Na2CO3
    # TIC total = 0.025 + 0.025 = 0.050 M
    # Na total  = 0.025*1 + 0.025*2 = 0.075 M
    run_case(
        "Case 2: NIST buffer 0.025 M NaHCO3 + 0.025 M Na2CO3",
        CT_TIC=0.050,
        CT_Na=0.075,
        expected_pH=10.02,  # NIST reference at 25°C
        tol=0.15,
    )
