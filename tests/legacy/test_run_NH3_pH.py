# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 18:01:49 2026

@author: k2473520
"""

# -*- coding: utf-8 -*-
"""
test_run_NH3_pH.py

Compute predicted pH of aqueous ammonia (NH3) solutions from WT_MIN to WT_MAX (wt%)
using your BisectionChemicalEquilibriumEngine (no fermenter import).

Notes:
- This uses a simple density(wt%) approximation via a small lookup table + linear interpolation.
  Replace the density table with your preferred correlation / datasheet values if you have them.
- For concentrated NH3(aq) (high ionic strength), Davies/ideal activity models can be questionable.
  Treat results at high wt% as "model-extrapolated" unless you have a better activity model.

Run in Spyder: press Run.
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
WT_MAX = 30.0      # typical commercial ammonia solutions go up to ~28–30 wt%
N_POINTS = 80      # number of wt% points

# -----------------------------
# Density model (kg/L) for NH3(aq) at ~20–25°C
# (Approximate values; linear interpolation is used.)
# Replace with a datasheet or correlation for your temperature & concentration range.
# -----------------------------
_DENSITY_TABLE = np.array([
    # wt%, density kg/L
    (0.0,  0.998),
    (1.0,  0.995),
    (5.0,  0.985),
    (10.0, 0.970),
    (15.0, 0.955),
    (20.0, 0.940),
    (25.0, 0.925),
    (28.0, 0.917),
    (30.0, 0.912),
], dtype=float)


def density_nh3_kg_L(wt_percent: float) -> float:
    """Linear interpolation over _DENSITY_TABLE; clamps outside range."""
    wt = float(wt_percent)
    x = _DENSITY_TABLE[:, 0]
    y = _DENSITY_TABLE[:, 1]
    if wt <= x.min():
        return float(y[0])
    if wt >= x.max():
        return float(y[-1])
    return float(np.interp(wt, x, y))


def recipe_nh3_from_wt_percent(wt_percent: float, *, density_kg_L: float) -> dict:
    """
    Return recipe in g/L: {"NH3": g_per_L}.
    wt% is mass fraction of NH3 in total solution.
    """
    wt_percent = float(wt_percent)
    if wt_percent < 0 or wt_percent > 100:
        raise ValueError("wt_percent must be between 0 and 100.")
    mass_solution_g_L = float(density_kg_L) * 1000.0
    mass_nh3_g_L = mass_solution_g_L * (wt_percent / 100.0)
    return {"NH3": mass_nh3_g_L}


def gL_to_molL_nh3(g_per_L: float) -> float:
    MW_NH3 = 17.031  # g/mol
    return float(g_per_L) / MW_NH3


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
    rho = density_nh3_kg_L(wt)  # kg/L
    rec_gL = recipe_nh3_from_wt_percent(wt, density_kg_L=rho)
    nh3_gL = rec_gL["NH3"]
    nh3_molL = gL_to_molL_nh3(nh3_gL)

    # NH3 is neutral; feed it as total ammonia CT_NH_T.
    # The solver will distribute CT_NH_T between NH3(aq) and NH4+ based on pH,
    # and enforce electroneutrality (OH- emerges accordingly).
    out = engine.solve(
        acid_totals={},
        acid_pKas={},
        CT_TIC=0.0,
        CT_P=0.0,
        CT_NH_T=nh3_molL,  # <-- ammonia total (NH3 + NH4+)
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
        n_scan=400,
        tol=1e-12,
    )

    # Helpful debug line (comment out if noisy)
    print(out.get("pH"), out.get("H+"), out.get("OH-"), out.get("NH4+"), out.get("NH3"), out.get("aH"))

    rows.append({
        "wt_percent": float(wt),
        "density_kg_L": float(rho),
        "NH3_g_L": float(nh3_gL),
        "NH3_mol_L": float(nh3_molL),
        "pH": float(out.get("pH", np.nan)),
        "IonicStrength": float(out.get("IonicStrength", np.nan)),
        "NH4+_mol_L": float(out.get("NH4+", np.nan)),
        "NH3aq_mol_L": float(out.get("NH3", np.nan)),
    })

# -----------------------------
# Print a copy/paste friendly table
# -----------------------------
print("\nNH3(aq) pH vs wt%  (T = {:.1f} °C, level={}, activity={}, model={})".format(
    T_C, LEVEL, USE_ACTIVITY, ACTIVITY_MODEL
))
print("{:>8s} {:>10s} {:>12s} {:>12s} {:>8s} {:>12s} {:>12s} {:>12s}".format(
    "wt%", "rho(kg/L)", "NH3(g/L)", "NH3(mol/L)", "pH", "I (mol/L)", "NH4+(M)", "NH3(aq)(M)"
))
print("-" * 110)

for r in rows:
    print("{:8.3f} {:10.3f} {:12.3f} {:12.4f} {:8.3f} {:12.4f} {:12.4f} {:12.4f}".format(
        r["wt_percent"],
        r["density_kg_L"],
        r["NH3_g_L"],
        r["NH3_mol_L"],
        r["pH"],
        r["IonicStrength"],
        r["NH4+_mol_L"],
        r["NH3aq_mol_L"],
    ))

# -----------------------------
# Optional plot (matplotlib)
# -----------------------------
try:
    import matplotlib.pyplot as plt

    wt = [r["wt_percent"] for r in rows]
    pH = [r["pH"] for r in rows]
    I  = [r["IonicStrength"] for r in rows]

    fig, ax1 = plt.subplots()

    # Primary axis: pH
    ax1.plot(wt, pH, linestyle="-", label='pH')
    ax1.set_xlabel("NH3 (wt%)")
    ax1.set_ylabel("Predicted pH")
    ax1.tick_params(axis="y")

    # Secondary axis: Ionic Strength
    ax2 = ax1.twinx()
    ax2.plot(wt, I, linestyle="--", label='I')
    ax2.set_ylabel("Ionic Strength (mol/L)")
    ax2.tick_params(axis="y")

    plt.title("NH3(aq) pH and Ionic Strength vs wt%")
    plt.legend()
    fig.tight_layout()
    plt.show()
    

except Exception as e:
    print("\n(Plot skipped:", e, ")")
