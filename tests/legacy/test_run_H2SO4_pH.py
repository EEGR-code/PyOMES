# -*- coding: utf-8 -*-
"""
h2so4_pH_vs_wtpercent.py

Compute predicted pH of H2SO4 solutions from WT_MIN to WT_MAX (wt%)
using your BisectionChemicalEquilibriumEngine (no fermenter_unit import).

Notes:
- Uses a simple density(wt%) approximation via lookup table + linear interpolation.
- For concentrated H2SO4 (high ionic strength), Davies/ideal activity models become questionable.
  Treat results above ~1–2 M as "model-extrapolated" unless you have a better activity model.
- H2SO4 is treated as adding total sulfate (S(VI)) with no counter-cations; electroneutrality drives [H+].
  The bisulfate equilibrium (HSO4- <-> H+ + SO4--) is handled by Level 2 via pKa_HSO4.
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
N_POINTS = 80  # number of wt% points

# -----------------------------
# Density model (kg/L) for H2SO4(aq) at ~20–25°C
# (Approximate values; linear interpolation is used.)
# Replace with your own correlation or data if needed.
# -----------------------------
_DENSITY_TABLE = np.array([
    # wt%, density kg/L
    (0.0,   0.998),
    (5.0,   1.030),
    (10.0,  1.066),
    (15.0,  1.102),
    (20.0,  1.139),
    (25.0,  1.179),
    (30.0,  1.219),
    (35.0,  1.261),
    (40.0,  1.303),
    (45.0,  1.350),
    (50.0,  1.398),
], dtype=float)


def density_h2so4_kg_L(wt_percent: float) -> float:
    """Linear interpolation over _DENSITY_TABLE; clamps outside range."""
    wt = float(wt_percent)
    x = _DENSITY_TABLE[:, 0]
    y = _DENSITY_TABLE[:, 1]
    if wt <= x.min():
        return float(y[0])
    if wt >= x.max():
        return float(y[-1])
    return float(np.interp(wt, x, y))


def recipe_h2so4_from_wt_percent(wt_percent: float, *, density_kg_L: float) -> dict:
    """
    Return recipe in g/L: {"H2SO4": g_per_L}.
    wt% is mass fraction of H2SO4 in total solution.
    """
    wt_percent = float(wt_percent)
    if wt_percent < 0 or wt_percent > 100:
        raise ValueError("wt_percent must be between 0 and 100.")
    mass_solution_g_L = float(density_kg_L) * 1000.0
    mass_h2so4_g_L = mass_solution_g_L * (wt_percent / 100.0)
    return {"H2SO4": mass_h2so4_g_L}


def gL_to_molL_h2so4(g_per_L: float) -> float:
    MW_H2SO4 = 98.079  # g/mol
    return float(g_per_L) / MW_H2SO4


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
    rho = density_h2so4_kg_L(wt)                 # kg/L
    rec_gL = recipe_h2so4_from_wt_percent(wt, density_kg_L=rho)
    h2so4_gL = rec_gL["H2SO4"]
    h2so4_molL = gL_to_molL_h2so4(h2so4_gL)

    # Treat H2SO4 addition as total sulfate CT_SO4 (mol/L).
    # Level 2 includes bisulfate equilibrium via pKa_HSO4; [H+] emerges from electroneutrality.
    out = engine.solve(
        acid_totals={},
        acid_pKas={},
        CT_TIC=0.0,
        CT_P=0.0,
        CT_NH_T=0.0,

        CT_Na=0.0,
        CT_Cl=0.0,
        CT_K=0.0,
        CT_NO3=0.0,

        CT_SO4=h2so4_molL,  # total S(VI) from H2SO4
        CT_Mg=0.0,
        CT_Ca=0.0,
        CT_Zn=0.0,
        CT_Mn=0.0,
        CT_Co=0.0,
        CT_Mo7O24=0.0,

        # For strong acids, allow negative pH search space
        pH_min=-3.0,
        pH_max=14.0,
        n_scan=600,
        tol=1e-12,
    )

    print(out.get("pH"), out.get("H+"), out.get("HSO4-"), out.get("SO4--"), out.get("aH"))

    rows.append({
        "wt_percent": float(wt),
        "density_kg_L": float(rho),
        "H2SO4_g_L": float(h2so4_gL),
        "H2SO4_mol_L": float(h2so4_molL),
        "pH": float(out.get("pH", np.nan)),
        "IonicStrength": float(out.get("IonicStrength", np.nan)),
    })

# -----------------------------
# Print a copy/paste friendly table
# -----------------------------
print("\nH2SO4 pH vs wt%  (T = {:.1f} °C, level={}, activity={}, model={})".format(
    T_C, LEVEL, USE_ACTIVITY, ACTIVITY_MODEL
))
print("{:>8s} {:>10s} {:>12s} {:>12s} {:>8s} {:>12s}".format(
    "wt%", "rho(kg/L)", "H2SO4(g/L)", "H2SO4(mol/L)", "pH", "I (mol/L)"
))
print("-" * 74)

for r in rows:
    print("{:8.3f} {:10.3f} {:12.3f} {:12.4f} {:8.3f} {:12.4f}".format(
        r["wt_percent"],
        r["density_kg_L"],
        r["H2SO4_g_L"],
        r["H2SO4_mol_L"],
        r["pH"],
        r["IonicStrength"],
    ))

# -----------------------------
# Optional plot (matplotlib)
# -----------------------------
try:
    import matplotlib.pyplot as plt

    wt_vals = [r["wt_percent"] for r in rows]
    pH_vals = [r["pH"] for r in rows]
    I_vals  = [r["IonicStrength"] for r in rows]

    fig, ax1 = plt.subplots()

    # Primary y-axis: pH
    ax1.plot(wt_vals, pH_vals)
    ax1.set_xlabel("H2SO4 (wt%)")
    ax1.set_ylabel("Predicted pH")
    ax1.set_title("H2SO4 pH and Ionic Strength vs wt%")

    # Secondary y-axis: Ionic Strength
    ax2 = ax1.twinx()
    ax2.plot(wt_vals, I_vals, linestyle="--")
    ax2.set_ylabel("Ionic Strength (mol/L)")

    plt.tight_layout()
    plt.show()

except Exception as e:
    print("\n(Plot skipped:", e, ")")

