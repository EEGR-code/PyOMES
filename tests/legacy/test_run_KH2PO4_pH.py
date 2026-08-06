# -*- coding: utf-8 -*-
"""
test_run_K_phosphate_pH.py

Sweep wt% potassium phosphate and compute pH + ionic strength using BisectionChemicalEquilibriumEngine.

Select which potassium phosphate salt via SALT_FORM:
  - "K3PO4"  : K3PO4 (tribasic potassium phosphate)        -> CT_K = 3*CT_P
  - "K2HPO4" : K2HPO4 (dipotassium hydrogen phosphate)     -> CT_K = 2*CT_P
  - "KH2PO4" : KH2PO4 (monopotassium dihydrogen phosphate) -> CT_K = 1*CT_P

Notes:
- Density correlation here is a rough approximation. Replace the density table/correlation
  with your preferred data if you need quantitative accuracy at high wt%.
- For basic salts (especially K3PO4), pH can approach the upper end of the scan.
  Increase pH_max or n_scan if you see clipping.
"""

import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

from speciation.engine import BisectionChemicalEquilibriumEngine

# -----------------------------
# USER SETTINGS
# -----------------------------
SALT_FORM = "K3PO4"  # "K3PO4", "K2HPO4", or "KH2PO4"

LEVEL = 2
USE_ACTIVITY = True
ACTIVITY_MODEL = "davies"
T_C = 25.0

WT_MIN = 0.1
WT_MAX = 30.0      # many potassium phosphates have limited solubility; adjust if needed
N_POINTS = 80

# pH scan bounds (salts can be basic)
PH_MIN = 0.0
PH_MAX = 16.0
N_SCAN = 700

# -----------------------------
# Molecular weights (g/mol) & stoichiometry
# -----------------------------
SALT_DATA = {
    "K3PO4":  {"MW": 212.27, "K_PER_P": 3},
    "K2HPO4": {"MW": 174.18, "K_PER_P": 2},
    "KH2PO4": {"MW": 136.09, "K_PER_P": 1},
}

if SALT_FORM not in SALT_DATA:
    raise ValueError(f"Unknown SALT_FORM={SALT_FORM}. Choose one of: {list(SALT_DATA.keys())}")

MW_SALT = SALT_DATA[SALT_FORM]["MW"]
K_PER_P = SALT_DATA[SALT_FORM]["K_PER_P"]

# -----------------------------
# Density model (kg/L) - rough generic for "salt solutions"
# Replace with better data if you have it.
# -----------------------------
_DENSITY_TABLE = np.array([
    # wt%, density kg/L (approx at ~20–25 °C)
    (0.0,  0.998),
    (5.0,  1.035),
    (10.0, 1.075),
    (15.0, 1.115),
    (20.0, 1.155),
    (25.0, 1.200),
    (30.0, 1.245),
], dtype=float)


def density_kg_L(wt_percent: float) -> float:
    """Linear interpolation over _DENSITY_TABLE; clamps outside range."""
    wt = float(wt_percent)
    x = _DENSITY_TABLE[:, 0]
    y = _DENSITY_TABLE[:, 1]
    if wt <= x.min():
        return float(y[0])
    if wt >= x.max():
        return float(y[-1])
    return float(np.interp(wt, x, y))


def salt_gL_from_wt_percent(wt_percent: float, rho_kg_L: float) -> float:
    """Convert wt% of salt in solution to g/L using density rho_kg_L."""
    mass_solution_g_L = float(rho_kg_L) * 1000.0
    return mass_solution_g_L * (float(wt_percent) / 100.0)


def gL_to_molL_salt(g_per_L: float) -> float:
    return float(g_per_L) / MW_SALT


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
    rho = density_kg_L(wt)
    salt_gL = salt_gL_from_wt_percent(wt, rho_kg_L=rho)
    salt_molL = gL_to_molL_salt(salt_gL)

    # Interpret "salt_molL" as total phosphate added (1 P per formula unit),
    # with K counter-ions according to salt stoichiometry.
    CT_P = salt_molL
    CT_K = K_PER_P * salt_molL

    out = engine.solve(
        acid_totals={},
        acid_pKas={},

        CT_TIC=0.0,
        CT_P=CT_P,
        CT_NH_T=0.0,

        CT_Na=0.0,
        CT_Cl=0.0,
        CT_K=CT_K,
        CT_NO3=0.0,
        CT_SO4=0.0,

        CT_Mg=0.0,
        CT_Ca=0.0,
        CT_Zn=0.0,
        CT_Mn=0.0,
        CT_Co=0.0,
        CT_Mo7O24=0.0,

        pH_min=PH_MIN,
        pH_max=PH_MAX,
        n_scan=N_SCAN,
        tol=1e-12,
    )

    # Debug print (similar spirit to your others)
    print(
        out.get("pH"),
        out.get("H+"),
        out.get("OH-"),
        out.get("H3PO4"),
        out.get("H2PO4-"),
        out.get("HPO4--"),
        out.get("PO4---"),
        out.get("IonicStrength"),
    )

    rows.append({
        "wt_percent": float(wt),
        "density_kg_L": float(rho),
        "salt_g_L": float(salt_gL),
        "salt_mol_L": float(salt_molL),
        "CT_P_mol_L": float(CT_P),
        "CT_K_mol_L": float(CT_K),
        "pH": float(out.get("pH", np.nan)),
        "IonicStrength": float(out.get("IonicStrength", np.nan)),
    })

# -----------------------------
# Print table
# -----------------------------
print(f"\n{SALT_FORM} pH vs wt%  (T = {T_C:.1f} °C, level={LEVEL}, activity={USE_ACTIVITY}, model={ACTIVITY_MODEL})")
print("{:>8s} {:>10s} {:>12s} {:>12s} {:>10s} {:>8s} {:>12s}".format(
    "wt%", "rho(kg/L)", "salt(g/L)", "salt(mol/L)", "CT_K(M)", "pH", "I (mol/L)"
))
print("-" * 82)

for r in rows:
    print("{:8.3f} {:10.3f} {:12.3f} {:12.4f} {:10.4f} {:8.3f} {:12.4f}".format(
        r["wt_percent"],
        r["density_kg_L"],
        r["salt_g_L"],
        r["salt_mol_L"],
        r["CT_K_mol_L"],
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
    ax1.set_xlabel(f"{SALT_FORM} (wt%)")
    ax1.set_ylabel("Predicted pH")
    ax1.set_title(f"{SALT_FORM}: pH and Ionic Strength vs wt%")

    ax2 = ax1.twinx()
    ax2.plot(wt_vals, I_vals, linestyle="--")
    ax2.set_ylabel("Ionic Strength (mol/L)")

    plt.tight_layout()
    plt.show()

except Exception as e:
    print("\n(Plot skipped:", e, ")")
