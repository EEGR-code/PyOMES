# -*- coding: utf-8 -*-
"""
naoh_pH_vs_wtpercent.py

Compute predicted pH of NaOH solutions from 0.1 wt% to 50 wt%
using your BisectionChemicalEquilibriumEngine (no fermentation_unit import).

Notes:
- This uses a simple density(wt%) approximation via a small lookup table + linear interpolation.
- For concentrated NaOH (high ionic strength), Davies/ideal activity models become questionable.
  Treat results above ~1–2 M as "model-extrapolated" unless you have a better activity model.
- Works in Spyder: press Run.
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
# Density model (kg/L) for NaOH(aq) at ~20–25°C
# (Approximate values; linear interpolation is used.)
# You can replace this with your own correlation or CRM data.
# -----------------------------
_DENSITY_TABLE = np.array([
    # wt%, density kg/L
    (0.0,  0.998),
    (1.0,  1.010),
    (5.0,  1.053),
    (10.0, 1.108),
    (15.0, 1.166),
    (20.0, 1.219),
    (25.0, 1.274),
    (30.0, 1.328),
    (35.0, 1.382),
    (40.0, 1.437),
    (45.0, 1.493),
    (50.0, 1.538),
], dtype=float)


def density_naoh_kg_L(wt_percent: float) -> float:
    """Linear interpolation over _DENSITY_TABLE; clamps outside range."""
    wt = float(wt_percent)
    x = _DENSITY_TABLE[:, 0]
    y = _DENSITY_TABLE[:, 1]
    if wt <= x.min():
        return float(y[0])
    if wt >= x.max():
        return float(y[-1])
    return float(np.interp(wt, x, y))


def recipe_naoh_from_wt_percent(wt_percent: float, *, density_kg_L: float) -> dict:
    """
    Return recipe in g/L: {"NaOH": g_per_L}.
    wt% is mass fraction of NaOH in total solution.
    """
    wt_percent = float(wt_percent)
    if wt_percent < 0 or wt_percent > 100:
        raise ValueError("wt_percent must be between 0 and 100.")
    mass_solution_g_L = float(density_kg_L) * 1000.0
    mass_naoh_g_L = mass_solution_g_L * (wt_percent / 100.0)
    return {"NaOH": mass_naoh_g_L}


def gL_to_molL_naoh(g_per_L: float) -> float:
    MW_NaOH = 40.00  # g/mol
    return float(g_per_L) / MW_NaOH


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
    rho = density_naoh_kg_L(wt)                 # kg/L
    rec_gL = recipe_naoh_from_wt_percent(wt, density_kg_L=rho)
    naoh_gL = rec_gL["NaOH"]
    naoh_molL = gL_to_molL_naoh(naoh_gL)

    # Your CHEM_DB mapping for NaOH should be: NaOH -> CT_Na += 1 (no anion)
    # Here we call the solver directly with totals (no need for chem_recipe).
    out = engine.solve(
        acid_totals={},
        acid_pKas={},
        CT_TIC=0.0,
        CT_P=0.0,
        CT_NH_T=0.0,
        CT_Na=naoh_molL,   # NaOH contributes Na+; OH- emerges from electroneutrality
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
    
    print(out.get("pH"), out.get("H+"), out.get("OH-"), out.get("Na+"), out.get("aH"))

    rows.append({
        "wt_percent": float(wt),
        "density_kg_L": float(rho),
        "NaOH_g_L": float(naoh_gL),
        "NaOH_mol_L": float(naoh_molL),
        "pH": float(out.get("pH", np.nan)),
        "IonicStrength": float(out.get("IonicStrength", np.nan)),
    })

# -----------------------------
# Print a copy/paste friendly table
# -----------------------------
print("\nNaOH pH vs wt%  (T = {:.1f} °C, level={}, activity={}, model={})".format(
    T_C, LEVEL, USE_ACTIVITY, ACTIVITY_MODEL
))
print("{:>8s} {:>10s} {:>12s} {:>12s} {:>8s} {:>12s}".format(
    "wt%", "rho(kg/L)", "NaOH(g/L)", "NaOH(mol/L)", "pH", "I (mol/L)"
))
print("-" * 70)

for r in rows:
    print("{:8.3f} {:10.3f} {:12.3f} {:12.4f} {:8.3f} {:12.4f}".format(
        r["wt_percent"],
        r["density_kg_L"],
        r["NaOH_g_L"],
        r["NaOH_mol_L"],
        r["pH"],
        r["IonicStrength"],
    ))

# -----------------------------
# Optional plot (matplotlib)
# -----------------------------
try:
    import matplotlib.pyplot as plt

    plt.figure()
    plt.plot([r["wt_percent"] for r in rows], [r["pH"] for r in rows])
    plt.xlabel("NaOH (wt%)")
    plt.ylabel("Predicted pH")
    plt.title("NaOH pH vs wt%")
    plt.show()
except Exception as e:
    print("\n(Plot skipped:", e, ")")
