# -*- coding: utf-8 -*-
"""
pH sensitivity + uncertainty propagation (Spyder-friendly)

Fixes:
- recipe_to_totals expects Dict[str, float] (NOT nested dicts)
- ensures imports come from THIS folder (avoids accidentally importing v8_split_15)
- optional g/L bounds converted using existing recipe_g_L_to_mol_L

Outputs:
- nominal pH (midpoint recipe)
- mean/std pH
- 95% CI
- sensitivity ranking (standardized linear coefficients)
"""

import os
import sys
import json
import numpy as np
import warnings

# ------------------------------------------------------------------
# FORCE LOCAL IMPORTS (prevents accidentally importing v8_split_15)
# Put this script inside v8_split_16/ and run it from there.
# ------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
os.chdir(HERE)

# ------------------------------------------------------------------
# SPYDER DEFAULTS – EDIT THESE
# ------------------------------------------------------------------

USE_BOUNDS_FILE = False           # True → read bounds from JSON
BOUNDS_FILE = "bounds.json"       # Only used if USE_BOUNDS_FILE=True

# Bounds for MATERIALS ADDED (same names as used in CHEM_DB / recipe scripts)
# Provide as (lower, upper) in either mol/L or g/L depending on UNITS.
BOUNDS_DICT = {
    "NaHCO3": (0.045, 0.055),
    "Na2CO3": (0.018, 0.022),
    "NaCl":   (0.095, 0.105),
}

recipe_g_L = {"KH2PO4":   20.0,
               "NH4Cl":     5.5,
               "K2SO4":     0.4,
               "MgSO4":     0.4 * 120/246.5, # HEPTAHYDRATE!
               "K2HPO4":    5.0,
               }  # 1 g/L each chemical

pct = 0.05
BOUNDS_DICT = {key: (val*(1-pct), val*(1+pct)) for key, val in recipe_g_L.items()}

UNITS = "mol_L"                   # "mol_L" or "g_L"
N_SAMPLES = 5000
RANDOM_SEED = 1

# Only needed for acid systems that require pKas (same as NIST recipe script pattern)
ACID_PKAS = {
    # "AceticAcid": 4.76,
}

# ------------------------------------------------------------------
# IMPORT EXISTING CODE (keeps pH consistency with codebase)
# ------------------------------------------------------------------
from scripts.chem_recipe import recipe_to_totals, recipe_g_L_to_mol_L
from speciation.engine import BisectionChemicalEquilibriumEngine


def load_bounds():
    if USE_BOUNDS_FILE:
        with open(BOUNDS_FILE, "r") as f:
            raw = json.load(f)
        return {k: tuple(v) for k, v in raw.items()}
    return BOUNDS_DICT


def recipe_from_sample(sample_dict):
    """
    Return recipe_mol_L: Dict[str, float] as expected by recipe_to_totals().

    If UNITS == "g_L", convert using existing helper recipe_g_L_to_mol_L().
    """
    if UNITS == "mol_L":
        return {k: float(v) for k, v in sample_dict.items()}
    elif UNITS == "g_L":
        recipe_gL = {k: float(v) for k, v in sample_dict.items()}
        return recipe_g_L_to_mol_L(recipe_gL)
    else:
        raise ValueError("UNITS must be 'mol_L' or 'g_L'")


def compute_pH_from_recipe_mol_L(recipe_mol_L):
    """
    Uses the same flow as run_test_NIST_buffer_standards_with_recipe.py:
      recipe -> recipe_to_totals -> engine.solve(**inputs) -> pH
    """
    totals = recipe_to_totals(recipe_mol_L, acid_pKas=ACID_PKAS)

    strong = totals["strong_ions"]
    engine = BisectionChemicalEquilibriumEngine()  # defaults should match your tests unless you changed them

    solve_kwargs = dict(
        acid_totals=totals["acid_totals"],
        acid_pKas=totals["acid_pKas"],
        CT_TIC=totals["CT_TIC"],
        CT_P=totals["CT_P"],
        CT_NH_T=totals["CT_NH_T"],
    )

    # Pass any strong ions that are present (CT_Na, CT_K, CT_Cl, CT_SO4, ...)
    for k, v in strong.items():
        solve_kwargs[k] = float(v)

    sol = engine.solve(**solve_kwargs)
    return float(sol["pH"])


def run_sensitivity_analysis():
    np.random.seed(RANDOM_SEED)

    bounds = load_bounds()
    materials = list(bounds.keys())

    # midpoint nominal
    midpoint = {m: 0.5 * (bounds[m][0] + bounds[m][1]) for m in materials}
    nominal_recipe_mol_L = recipe_from_sample(midpoint)
    pH_nominal = compute_pH_from_recipe_mol_L(nominal_recipe_mol_L)

    # Monte Carlo
    pH_samples = np.empty(N_SAMPLES, dtype=float)
    X_samples = np.empty((N_SAMPLES, len(materials)), dtype=float)

    for i in range(N_SAMPLES):
        sample = {}
        for j, m in enumerate(materials):
            lo, hi = bounds[m]
            val = np.random.uniform(lo, hi)
            sample[m] = val
            X_samples[i, j] = val  # store in original units (mol/L or g/L) for sensitivity reporting

        recipe_mol_L = recipe_from_sample(sample)

        try:
            pH_samples[i] = compute_pH_from_recipe_mol_L(recipe_mol_L)
        except Exception:
            pH_samples[i] = np.nan

    mask = np.isfinite(pH_samples)
    pH_samples = pH_samples[mask]
    X_samples = X_samples[mask, :]

    if pH_samples.size == 0:
        raise RuntimeError(
            "All samples failed to solve. Check material names, bounds, and ACID_PKAS."
        )

    # Uncertainty stats
    pH_mean = float(np.mean(pH_samples))
    pH_std = float(np.std(pH_samples))
    pH_ci = tuple(np.percentile(pH_samples, [2.5, 97.5]).tolist())

    # Sensitivity ranking: standardized linear coefficients
    Xz = (X_samples - X_samples.mean(axis=0)) / X_samples.std(axis=0)
    yz = (pH_samples - pH_samples.mean()) / pH_samples.std()
    coeffs = np.linalg.lstsq(Xz, yz, rcond=None)[0]
    sens = dict(zip(materials, coeffs))

    # Print
    print("\n================ pH SENSITIVITY ANALYSIS ================\n")
    print(f"Working directory: {os.getcwd()}")
    print(f"Recipe units (bounds): {UNITS}")
    print(f"Valid solves: {pH_samples.size} / {N_SAMPLES}\n")

    print(f"Nominal pH (midpoint recipe): {pH_nominal:8.4f}\n")

    print("Uncertainty (Monte Carlo):")
    print(f"  Mean pH : {pH_mean:8.4f}")
    print(f"  Std  pH : {pH_std:8.4f}")
    print(f"  95% CI  : [{pH_ci[0]:.4f}, {pH_ci[1]:.4f}]\n")
    print(f"  Min.Max.: [{min(pH_samples):.4f}, {max(pH_samples):.4f}]\n")

    print("Sensitivity ranking (|standardized coefficient|):")
    for m, c in sorted(sens.items(), key=lambda x: abs(x[1]), reverse=True):
        print(f"  {m:15s}: {c:+.4f}")

    print("\n========================================================\n")

    return {
        "pH_nominal": pH_nominal,
        "pH_mean": pH_mean,
        "pH_std": pH_std,
        "pH_CI": pH_ci,
        "sensitivity": sens,
        "pH_samples": pH_samples,
    }


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        results = run_sensitivity_analysis()
