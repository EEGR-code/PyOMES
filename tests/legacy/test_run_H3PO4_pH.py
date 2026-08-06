# -*- coding: utf-8 -*-
"""
h3po4_pH_vs_wtpercent.py

Compute predicted pH of H3PO4 solutions from WT_MIN to WT_MAX (wt%)
using your BisectionChemicalEquilibriumEngine (no fermentation_unit import).

Notes:
- Uses a simple density(wt%) approximation via a small lookup table + linear interpolation.
- At higher ionic strength, Davies can be unreliable (engine will warn accordingly).
"""

import os
import numpy as np

# --- make relative imports work when running from IDE ---
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

from speciation.engine import BisectionChemicalEquilibriumEngine

# -----------------------------
# USER SETTINGS
# -----------------------------
LEVEL = 2
USE_ACTIVITY = True
ACTIVITY_MODEL = "davies"
T_C = 25.0

WT_MIN = 0.1
WT_MAX = 50.0
N_POINTS = 20  # number of wt% points

# -----------------------------
# Density model (kg/L) for H3PO4(aq) at ~20–25°C
# Approximate values; linear interpolation is used.
# Replace with your own correlation/data if you have it.
# -----------------------------
_DENSITY_TABLE = np.array([
    # wt%, density kg/L
    (0.0,  0.998),
    (5.0,  1.030),
    (10.0, 1.060),
    (15.0, 1.090),
    (20.0, 1.120),
    (25.0, 1.150),
    (30.0, 1.180),
    (35.0, 1.210),
    (40.0, 1.240),
    (45.0, 1.270),
    (50.0, 1.300),
], dtype=float)


def density_h3po4_kg_L(wt_percent: float) -> float:
    """Linear interpolation over _DENSITY_TABLE; clamps outside range."""
    wt = float(wt_percent)
    x = _DENSITY_TABLE[:, 0]
    y = _DENSITY_TABLE[:, 1]
    if wt <= x.min():
        return float(y[0])
    if wt >= x.max():
        return float(y[-1])
    return float(np.interp(wt, x, y))


def recipe_h3po4_from_wt_percent(wt_percent: float, *, density_kg_L: float) -> dict:
    """
    Return recipe in g/L: {"H3PO4": g_per_L}.
    wt% is mass fraction of H3PO4 in total solution.
    """
    wt_percent = float(wt_percent)
    if wt_percent < 0 or wt_percent > 100:
        raise ValueError("wt_percent must be between 0 and 100.")
    mass_solution_g_L = float(density_kg_L) * 1000.0
    mass_h3po4_g_L = mass_solution_g_L * (wt_percent / 100.0)
    return {"H3PO4": mass_h3po4_g_L}


def gL_to_molL_h3po4(g_per_L: float) -> float:
    MW_H3PO4 = 97.994  # g/mol
    return float(g_per_L) / MW_H3PO4


import math

def _safe(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return float(default)

def diagnose_phosphate_solution(out: dict, CT_P: float, label: str = "") -> None:
    """
    Prints:
      - phosphate mass balance
      - charge balance residual
      - apparent Ka1/Ka2/Ka3 (from concentrations)
      - a quick weak-acid estimate for pH (using Ka1 only)
    Assumes no other acids/strong ions unless present in out.
    """

    H  = _safe(out.get("H+"))
    OH = _safe(out.get("OH-"))

    H3PO4  = _safe(out.get("H3PO4"))
    H2PO4  = _safe(out.get("H2PO4-"))
    HPO4   = _safe(out.get("HPO4--"))
    PO4    = _safe(out.get("PO4---"))

    # Strong ions (if present)
    K   = _safe(out.get("K+"))
    Na  = _safe(out.get("Na+"))
    Mg  = _safe(out.get("Mg++"))
    Ca  = _safe(out.get("Ca++"))
    Zn  = _safe(out.get("Zn++"))
    Mn  = _safe(out.get("Mn++"))
    Co  = _safe(out.get("Co++"))

    Cl  = _safe(out.get("Cl-"))
    NO3 = _safe(out.get("NO3-"))

    # Other systems if present
    HCO3 = _safe(out.get("HCO3-"))
    CO3  = _safe(out.get("CO3--"))
    NH4  = _safe(out.get("NH4+"))
    HSO4 = _safe(out.get("HSO4-"))
    SO4  = _safe(out.get("SO4--"))
    Mo7  = _safe(out.get("Mo7O24------"))

    # 1) Phosphate mass balance
    P_sum = H3PO4 + H2PO4 + HPO4 + PO4
    P_err = P_sum - float(CT_P)

    # 2) Charge balance (cation - anion), should be ~0
    cations = H + NH4 + K + Na + 2*Mg + 2*Ca + 2*Zn + 2*Mn + 2*Co
    anions  = (
        OH + Cl + NO3
        + HCO3 + 2*CO3
        + H2PO4 + 2*HPO4 + 3*PO4
        + HSO4 + 2*SO4
        + 6*Mo7
    )
    charge_resid = cations - anions

    # 3) Apparent Ka values from concentrations (ideal-form)
    def Ka_from(num, den):
        if den <= 0 or num <= 0:
            return float("nan")
        return num / den

    Ka1_app = Ka_from(H * H2PO4, H3PO4)
    Ka2_app = Ka_from(H * HPO4,  H2PO4)
    Ka3_app = Ka_from(H * PO4,   HPO4)

    def pKa(x):
        if not (x > 0.0) or not math.isfinite(x):
            return float("nan")
        return -math.log10(x)

    # 4) Weak-acid pH estimate using Ka1 only: H ~ sqrt(Ka1*C)
    # (This explains the “flattening” you saw at high CT_P.)
    if CT_P > 0 and math.isfinite(Ka1_app) and Ka1_app > 0:
        H_est = math.sqrt(Ka1_app * CT_P)
        pH_est = -math.log10(H_est) if H_est > 0 else float("nan")
    else:
        pH_est = float("nan")

    print("\n=== DIAGNOSTICS {} ===".format(label))
    print("pH           :", _safe(out.get("pH")))
    print("[H+] (M)     :", H)
    print("CT_P (M)     :", CT_P)
    print("P_sum (M)    :", P_sum, "   (error:", P_err, ")")
    print("Charge resid :", charge_resid, "   (should be ~0)")
    print("Ka1_app, pKa1:", Ka1_app, pKa(Ka1_app))
    print("Ka2_app, pKa2:", Ka2_app, pKa(Ka2_app))
    print("Ka3_app, pKa3:", Ka3_app, pKa(Ka3_app))
    print("Weak-acid pH est (Ka1 only):", pH_est)


# -----------------------------
# Run sweep
# -----------------------------
engine = BisectionChemicalEquilibriumEngine(
    level=LEVEL,
    use_activity=USE_ACTIVITY,
    activity_model=ACTIVITY_MODEL,
    T_C=T_C,
)

wts = np.linspace(WT_MIN, WT_MAX, N_POINTS)
rows = []

for wt in wts:
    rho = density_h3po4_kg_L(wt)                 # kg/L
    rec_gL = recipe_h3po4_from_wt_percent(wt, density_kg_L=rho)
    h3po4_gL = rec_gL["H3PO4"]
    ct_p_molL = gL_to_molL_h3po4(h3po4_gL)        # total phosphate (mol/L as P)

    out = engine.solve(
        acid_totals={},
        acid_pKas={},

        CT_TIC=0.0,
        CT_P=ct_p_molL,   # <- total phosphate drives H3PO4 system
        CT_NH_T=0.0,

        CT_Na=0.0,
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
        n_scan=500,
        tol=1e-12,
    )
    
    # Run diagnostics on a few representative points
    diagnose_phosphate_solution(out, CT_P=ct_p_molL, label=f"wt%={wt:.3f}")


    # Debug-style printout similar to your other scripts:
    print(
        out.get("pH"),
        out.get("H+"),
        out.get("H3PO4"),
        out.get("H2PO4-"),
        out.get("HPO4--"),
        out.get("PO4---"),
        out.get("aH"),
    )

    rows.append({
        "wt_percent": float(wt),
        "density_kg_L": float(rho),
        "H3PO4_g_L": float(h3po4_gL),
        "CT_P_mol_L": float(ct_p_molL),
        "pH": float(out.get("pH", np.nan)),
        "IonicStrength": float(out.get("IonicStrength", np.nan)),
    })

# -----------------------------
# Print a copy/paste friendly table
# -----------------------------
print("\nH3PO4 pH vs wt%  (T = {:.1f} °C, level={}, activity={}, model={})".format(
    T_C, LEVEL, USE_ACTIVITY, ACTIVITY_MODEL
))
print("{:>8s} {:>10s} {:>12s} {:>12s} {:>8s} {:>12s}".format(
    "wt%", "rho(kg/L)", "H3PO4(g/L)", "CT_P(mol/L)", "pH", "I (mol/L)"
))
print("-" * 74)

for r in rows:
    print("{:8.3f} {:10.3f} {:12.3f} {:12.4f} {:8.3f} {:12.4f}".format(
        r["wt_percent"],
        r["density_kg_L"],
        r["H3PO4_g_L"],
        r["CT_P_mol_L"],
        r["pH"],
        r["IonicStrength"],
    ))

# -----------------------------
# Plot with Ionic Strength on secondary y-axis
# -----------------------------
try:
    import matplotlib.pyplot as plt

    wt_vals = [r["wt_percent"] for r in rows]
    pH_vals = [r["pH"] for r in rows]
    I_vals  = [r["IonicStrength"] for r in rows]

    fig, ax1 = plt.subplots()

    ax1.plot(wt_vals, pH_vals)
    ax1.set_xlabel("H3PO4 (wt%)")
    ax1.set_ylabel("Predicted pH")
    ax1.set_title("H3PO4 pH and Ionic Strength vs wt%")

    ax2 = ax1.twinx()
    ax2.plot(wt_vals, I_vals, linestyle="--")
    ax2.set_ylabel("Ionic Strength (mol/L)")

    plt.tight_layout()
    plt.show()

except Exception as e:
    print("\n(Plot skipped:", e, ")")
