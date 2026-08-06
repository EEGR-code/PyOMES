# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 13:36:06 2026

@author: k2473520
"""

# optimize_profile_v8_split_10_with_plot.py
# -------------------------------------------------------------------------
# Single-file demonstration:
# 1) Builds the SAME BioSTEAM flowsheet as test_run_v8_split_9.py (replicated here)
# 2) Generates synthetic "training data" using the CURRENT kinetic parameterization
# 3) Fits parameters using scipy.optimize.differential_evolution by matching the
#    ENTIRE biomass growth profile
# 4) Plots biomass profile for TRUE vs FITTED parameters
#
# Notes:
# - This script tries to extract (t, X) from your fermenter robustly.
#   If it can't find the profile, it raises a clear error telling you what to change.
# -------------------------------------------------------------------------

import os
import numpy as np
import biosteam as bst
from dataclasses import dataclass
from typing import Dict, Mapping, Tuple, Any
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt

# ---- Make sure we import from the same folder this script lives in ----
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)

from v8_split_11.fermenter_unit import CUFermentationSpeciation as Fermenter
from v8_split_11.kinetics.core import ParameterSpec, Environment


# =============================================================================
# Plug-in kinetic model matching your legacy split-kinetics RHS
# =============================================================================
@dataclass(frozen=True)
class LegacySplitKineticsModel:
    name: str = "legacy_split_kinetics_model"
    state_ids: Tuple[str, ...] = ("X", "Sac", "Spro", "Sbut")

    param_schema: Tuple[ParameterSpec, ...] = (
        ParameterSpec("mu_m_ac",  default=0.5,  lower=0.0, units="1/h"),
        ParameterSpec("Y_ac",     default=0.9 * 0.400, lower=1e-12, units=""),
        ParameterSpec("Ks_ac",    default=5e-3, lower=0.0, units=""),

        ParameterSpec("mu_m_pro", default=0.29, lower=0.0, units="1/h"),
        ParameterSpec("Y_pro",    default=1.11 * 0.486, lower=1e-12, units=""),
        ParameterSpec("Ks_pro",   default=5e-3, lower=0.0, units=""),

        ParameterSpec("mu_m_but", default=0.5,  lower=0.0, units="1/h"),
        ParameterSpec("Y_but",    default=1.05 * 0.545, lower=1e-12, units=""),
        ParameterSpec("Ks_but",   default=5e-3, lower=0.0, units=""),
    )

    @staticmethod
    def _monod(mu_m: float, Ks: float, S: float) -> float:
        S = max(0.0, float(S))
        Ks = max(0.0, float(Ks))
        denom = Ks + S
        if denom <= 0.0:
            return 0.0
        return float(mu_m) * (S / denom)

    def rhs(self, t_h: float, y: np.ndarray, env: Environment, params: Mapping[str, float]):
        X_i, Sac_i, Spro_i, Sbut_i = [float(v) for v in y]
        X_i = max(0.0, X_i)
        Sac_i = max(0.0, Sac_i)
        Spro_i = max(0.0, Spro_i)
        Sbut_i = max(0.0, Sbut_i)

        mu_ac  = self._monod(params["mu_m_ac"],  params["Ks_ac"],  Sac_i)
        mu_pro = self._monod(params["mu_m_pro"], params["Ks_pro"], Spro_i)
        mu_but = self._monod(params["mu_m_but"], params["Ks_but"], Sbut_i)

        mu_tot = mu_ac + mu_pro + mu_but

        dXdt     = mu_tot * X_i
        dSacdt   = -(1.0 / max(params["Y_ac"],  1e-30)) * mu_ac  * X_i
        dSprodt  = -(1.0 / max(params["Y_pro"], 1e-30)) * mu_pro * X_i
        dSbutdt  = -(1.0 / max(params["Y_but"], 1e-30)) * mu_but * X_i

        dydt = np.array([dXdt, dSacdt, dSprodt, dSbutdt], dtype=float)
        outputs: Dict[str, float] = {"mu_tot": mu_tot, "mu_ac": mu_ac, "mu_pro": mu_pro, "mu_but": mu_but}
        return dydt, outputs


# =============================================================================
# Helpers copied from your test run
# =============================================================================
def solution_stream(ID_name: str, water_kgph: float, solutes_g_per_Lwater: dict):
    s = bst.Stream(ID_name, Water=water_kgph, units="kg/hr")
    Lwater_per_hr = water_kgph  # assume 1 kg/L
    for chem_id, g_per_L in solutes_g_per_Lwater.items():
        s.imass[chem_id] = (g_per_L * Lwater_per_hr) / 1000.0
    return s


# =============================================================================
# Build the same flowsheet
# =============================================================================
def build_fermenter():
    bst.nbtutorial()

    AmmoniumMolybdate = bst.Chemical("AmmoniumMolybdate", phase="s", search_db=True)
    AmmoniumMolybdate.V.add_model(1e-4)

    bst.settings.set_thermo([
        bst.Chemical("Water"),
        bst.Chemical("AceticAcid", CAS="64-19-7", phase="l"),
        bst.Chemical("PropionicAcid", CAS="79-09-4", phase="l"),
        bst.Chemical("ButyricAcid", CAS="107-92-6", phase="l"),
        bst.Chemical("O2", phase="g"),
        bst.Chemical("CO2", phase="g"),
        bst.Chemical("N2", phase="g"),

        bst.Chemical("AmmoniumSulfate", CAS="7783-20-2", phase="s"),
        bst.Chemical("KH2PO4", CAS="7778-77-0", phase="s"),
        bst.Chemical("MgSO4", CAS="10034-99-8", phase="s"),
        bst.Chemical("ZnSO4", CAS="7446-20-0", phase="s"),
        bst.Chemical("CaCl2", CAS="10043-52-4", phase="s"),
        bst.Chemical("MnCl2", CAS="7773-01-5", phase="s"),
        bst.Chemical("CoCl2", CAS="7646-79-9", phase="s"),
        bst.Chemical("CaSO4", CAS="7778-18-9", phase="s"),
        bst.Chemical("H3PO4", CAS="7664-38-2", phase="l"),
        AmmoniumMolybdate,

        bst.Chemical("Yeast", phase="s", search_db=False,
                     formula="CH1.61O0.56N0.16", rho=1540, Cp=1.5, default=True),
    ])

    composition_A = {"AceticAcid": 1.0, "AmmoniumSulfate": 2.0, "Yeast": 0.1}
    composition_B = {"KH2PO4": 20.0, "MgSO4": 10.0}

    MW_tetra = 1235.86
    MW_anhyd = AmmoniumMolybdate.MW
    composition_C = {
        "ZnSO4": 0.22,
        "CaCl2": 0.55,
        "MnCl2": 0.5,
        "CaSO4": 0.1,
        "CoCl2": 0.1,
        "AmmoniumMolybdate": 0.01 * (MW_anhyd / MW_tetra),
    }

    A = solution_stream("A", water_kgph=950.0, solutes_g_per_Lwater=composition_A)
    B = solution_stream("B", water_kgph=50.0, solutes_g_per_Lwater=composition_B)
    C = solution_stream("C", water_kgph=5.0, solutes_g_per_Lwater=composition_C)

    M1 = bst.units.Mixer("M1", ins=(A, B, C), outs=("substrate",))
    M1.simulate()
    substrate = M1.outs[0]

    H3PO4_makeup = bst.Stream("H3PO4_makeup", O2=0, CO2=0, Water=0, units="kg/hr")
    M2 = bst.units.Mixer("M2", ins=(substrate, H3PO4_makeup), outs=("feed",))
    M2.simulate()
    feed = M2.outs[0]

    gas_supply = bst.Stream("gas_supply", O2=0, CO2=0, Water=0, units="kg/hr")
    H3PO4_dose = bst.Stream("H3PO4_dose", H3PO4=0, Water=0, units="kg/hr")

    F = Fermenter(
        "F",
        ins=(feed, gas_supply, H3PO4_dose),
        outs=("vent", "broth"),
        tau=9.75,
        iskinetic=True,
        N=4,

        organism_id="Yeast",
        balance_basis="CHNO",
        n_source_id="NH3",

        yO2_init=1.0,
        yCO2_init=0.0,
        O2_vvm_min=2.0,

        pH_control=True,
        pH_setpoint=3.5,
        pH_Kp=0.6,
        pH_Ki=0.0,
        pH_max_add_molL_hr=0.3,
        pH_CT_P_max=1.0,

        speciation_level=1,
        use_activity=True,
        activity_model="davies",

        n_steps=300,
        use_split_solver=True,
    )
    return F


# =============================================================================
# Extract biomass growth profile from the fermenter (robust search)
# =============================================================================
def _try_get(obj: Any, names):
    for n in names:
        if isinstance(obj, dict) and n in obj:
            return obj[n]
        if hasattr(obj, n):
            return getattr(obj, n)
    return None


def extract_biomass_profile(F) -> Tuple[np.ndarray, np.ndarray]:
    containers = [
        _try_get(F, ["profile", "profiles", "history", "hist", "spec_profile", "spec", "last_spec", "speciation_profile"]),
        _try_get(F, ["_profile", "_profiles", "_history", "_spec_profile", "_spec"]),
    ]
    containers = [c for c in containers if c is not None]

    TIME_KEYS = ["t", "t_h", "time", "time_h", "times", "times_h", "t_eval", "time_grid", "ts"]
    X_KEYS = ["X_g_L", "X_g_per_L", "biomass_g_L", "biomass_g_per_L", "X", "biomass", "Yeast_g_L", "Yeast"]

    for c in containers:
        t = _try_get(c, TIME_KEYS)
        x = _try_get(c, X_KEYS)
        if t is not None and x is not None:
            t = np.asarray(t, dtype=float).ravel()
            x = np.asarray(x, dtype=float).ravel()
            if t.size == x.size and t.size >= 3:
                return t, x

    t = _try_get(F, TIME_KEYS) or _try_get(F, ["_t", "_t_h", "_time"])
    x = _try_get(F, X_KEYS) or _try_get(F, ["_X", "_biomass", "_Yeast"])
    if t is not None and x is not None:
        t = np.asarray(t, dtype=float).ravel()
        x = np.asarray(x, dtype=float).ravel()
        if t.size == x.size and t.size >= 3:
            return t, x

    sol = _try_get(F, ["sol", "solution", "_sol", "ode_sol", "_ode_sol"])
    if sol is not None:
        t2 = _try_get(sol, ["t"])
        y2 = _try_get(sol, ["y"])
        if t2 is not None and y2 is not None:
            t2 = np.asarray(t2, dtype=float).ravel()
            y2 = np.asarray(y2, dtype=float)
            if y2.ndim == 2 and y2.shape[1] == t2.size and t2.size >= 3:
                X_state = y2[0, :]
                return t2, np.asarray(X_state, dtype=float).ravel()

    raise RuntimeError(
        "Could not extract a biomass growth profile from the fermenter after simulation.\n"
        "Edit TIME_KEYS / X_KEYS in extract_biomass_profile() to match your internal names,\n"
        "or store solve_ivp outputs onto the fermenter (self.sol, self.profile).\n"
    )
    

def extract_acetate_profile(F) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns (t_hours, Sac_profile).

    Tries, in order:
      1) Containers on F that might store profiles as dict-like objects.
      2) Direct attributes on F.
      3) solve_ivp-like solution objects on F (sol, _sol, ode_sol, ...), using y[1,:]
         assuming Sac is the 2nd state (matches LegacySplitKineticsModel state_ids).
    """
    # Same container strategy as extract_biomass_profile:
    containers = [
        _try_get(F, ["profile", "profiles", "history", "hist", "spec_profile", "spec", "last_spec", "speciation_profile"]),
        _try_get(F, ["_profile", "_profiles", "_history", "_spec_profile", "_spec"]),
    ]
    containers = [c for c in containers if c is not None]

    TIME_KEYS = ["t", "t_h", "time", "time_h", "times", "times_h", "t_eval", "time_grid", "ts"]

    # Expand acetate key candidates (your code may label it in different ways)
    AC_KEYS = [
        "Sac", "S_ac", "S_Ac", "Acetate", "acetate",
        "AceticAcid", "AceticAcid_mol_L", "AceticAcid_mol",
        "Acetate_mol_L", "Acetate_mol",
        "VFA_Ac", "VFA_Ac_mol_L",
    ]

    # 1) Search dict-like containers
    for c in containers:
        t = _try_get(c, TIME_KEYS)
        ac = _try_get(c, AC_KEYS)
        if t is not None and ac is not None:
            t = np.asarray(t, dtype=float).ravel()
            ac = np.asarray(ac, dtype=float).ravel()
            if t.size == ac.size and t.size >= 3:
                return t, ac

    # 2) Search direct attributes on the fermenter
    t = _try_get(F, TIME_KEYS) or _try_get(F, ["_t", "_t_h", "_time"])
    ac = _try_get(F, AC_KEYS) or _try_get(F, ["_Sac", "_acetate"])
    if t is not None and ac is not None:
        t = np.asarray(t, dtype=float).ravel()
        ac = np.asarray(ac, dtype=float).ravel()
        if t.size == ac.size and t.size >= 3:
            return t, ac

    # 3) Search for solution objects (solve_ivp-style)
    sol = _try_get(F, ["sol", "solution", "_sol", "ode_sol", "_ode_sol", "sol_kin", "_sol_kin"])
    if sol is not None:
        t2 = _try_get(sol, ["t"])
        y2 = _try_get(sol, ["y"])
        if t2 is not None and y2 is not None:
            t2 = np.asarray(t2, dtype=float).ravel()
            y2 = np.asarray(y2, dtype=float)
            if y2.ndim == 2 and y2.shape[1] == t2.size and t2.size >= 3:
                # Assumption: acetate is 2nd state in y (index 1), consistent with your kinetic model
                Sac = y2[1, :]
                return t2, np.asarray(Sac, dtype=float).ravel()

    raise RuntimeError(
        "Could not extract acetate profile from PyOMES after simulation.\n"
        "This means acetate time-series is not being stored in F.profile/history/spec or as a solve_ivp solution.\n\n"
        "Two options:\n"
        "  (A) Add the correct key name to AC_KEYS in extract_acetate_profile().\n"
        "  (B) (Recommended) store the kinetics solution in fermenter_unit.py, e.g. self.sol = sol or self.profile['Sac']=...\n"
    )


# =============================================================================
# Loss on entire profile
# =============================================================================
def combined_profile_loss(
    t_ref,
    X_ref,
    S_ref,
    t_pred,
    X_pred,
    S_pred,
    w_X=1.0,
    w_S=1.0,
):
    # Interpolate predictions onto reference grid
    X_hat = np.interp(t_ref, t_pred, X_pred)
    S_hat = np.interp(t_ref, t_pred, S_pred)

    rX = X_hat - X_ref
    rS = S_hat - S_ref

    return float(w_X * np.sum(rX**2) + w_S * np.sum(rS**2))



# =============================================================================
# Plotting helper
# =============================================================================
# def plot_profiles(t_true, x_true, t_fit, x_fit, title="Biomass growth profile"):
#     plt.figure()
#     plt.scatter(t_true, x_true, label="True parameters")
#     plt.plot(t_fit, x_fit, "k", label="Fitted parameters")
#     plt.xlabel("Time (h)")
#     plt.ylabel("Biomass (as stored)")
#     plt.title(title)
#     plt.legend()
#     plt.tight_layout()
#     plt.show()
    
#     plt.figure()
#     plt.plot(x_true, x_fit, label="Parity")
#     plt.xlabel("True Biomass (h)")
#     plt.ylabel("Prdicted Biomass (as stored)")
#     plt.title(title)
#     plt.legend()
#     plt.tight_layout()
#     plt.show()
    
def plot_profiles(t_true, X_true, S_true,
                t_fit,  X_fit,  S_fit,
                ):
    fig, ax = plt.subplots(2, 1, sharex=True)

    ax[0].plot(t_true, X_true, label="True")
    ax[0].plot(t_fit,  X_fit,  "--", label="Fitted")
    ax[0].set_ylabel("Biomass")
    ax[0].legend()

    ax[1].plot(t_true, S_true, label="True")
    ax[1].plot(t_fit,  S_fit,  "--", label="Fitted")
    ax[1].set_ylabel("Acetate")
    ax[1].set_xlabel("Time (h)")
    ax[1].legend()

    plt.suptitle("Biomass and acetate: true vs fitted")
    plt.tight_layout()
    plt.show()



#%%

# =============================================================================
# Main
# =============================================================================
def main():
    model = LegacySplitKineticsModel()

    true_params = {
        "mu_m_ac": 0.5,
        "Y_ac": 0.9 * 0.400,
        "Ks_ac": 5e-3,
        "mu_m_pro": 0.29,
        "Y_pro": 1.11 * 0.486,
        "Ks_pro": 5e-3,
        "mu_m_but": 0.5,
        "Y_but": 1.05 * 0.545,
        "Ks_but": 5e-3,
    }

    # --- Generate synthetic training profile ---
    F_true = build_fermenter()
    F_true.set_kinetic_model(model=model, params=true_params)
    F_true.simulate()
    t_train, x_train = extract_biomass_profile(F_true)
    tS_train, S_train = extract_acetate_profile(F_true)

    print(f"\nTraining profile extracted: {t_train.size} points")
    print(f"  t range: {t_train.min():.3g} to {t_train.max():.3g} h")

    # --- Fit a subset of parameters (demo) ---
    fit_names = ["mu_m_ac", "Y_ac", "Ks_ac"]
    bounds = [
        (0.05, 1.5),    # mu_m_ac
        (0.05, 1.0),    # Y_ac
        (1e-5, 1e-1),   # Ks_ac
    ]

    def theta_to_params(theta):
        p = dict(true_params)
        for k, v in zip(fit_names, theta):
            p[k] = float(v)
        return p

    def objective(theta):
        params = theta_to_params(theta)
        F = build_fermenter()
        F.set_kinetic_model(model=model, params=params)
    
        try:
            F.simulate()
    
            tX, X_pred = extract_biomass_profile(F)
            tS, S_pred = extract_acetate_profile(F)
    
            return combined_profile_loss(
                t_ref=t_train,
                X_ref=x_train,    # <-- fixed name
                S_ref=S_train,
                t_pred=tX,
                X_pred=X_pred,
                S_pred=S_pred,
                w_X=1.0,
                w_S=1.0,
            )
        except Exception:
            return 1e30



    print("\nRunning differential evolution on full biomass profile...")
    result = differential_evolution(
        objective,
        bounds=bounds,
        maxiter=12,
        popsize=10,
        polish=True,
        disp=True,
        updating="deferred",
        workers=1,
        seed=1,
    )

    best_params = theta_to_params(result.x)
    print("\nOptimization finished.")
    print("Best objective:", result.fun)
    print("Best fitted parameters:")
    for k in fit_names:
        print(f"  {k:10s} = {best_params[k]:.6g}")

    # --- Simulate with best parameters and plot ---
    F_best = build_fermenter()
    F_best.set_kinetic_model(model=model, params=best_params)
    F_best.simulate()
    
    t_best, x_best = extract_biomass_profile(F_best)
    tS_best, S_best = extract_acetate_profile(F_best)
    
    plot_profiles(
        t_true=t_train, X_true=x_train, S_true=S_train,
        t_fit=t_best,   X_fit=x_best,   S_fit=S_best,
    )



if __name__ == "__main__":
    main()
