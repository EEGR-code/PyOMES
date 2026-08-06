#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Demo: 0-D continuous stirred-tank reactor (CSTR).

Migrated to the new ``Simulation`` orchestrator (simulation-class C12;
this demo completed post-ship as one of the three remaining demos).

Tutorial walkthrough of a chemostat: a perfectly mixed reactor with a
constant liquid feed of fresh substrate and a matching liquid drain,
so the working volume is held constant. At steady state, biomass and
substrate concentrations balance growth against washout.

Walks through the same pattern as ``batch_fermenter.py``, plus the
boundary plumbing for feed/drain and a DO controller that adjusts
kLa to maintain a dissolved-O2 setpoint.

Run from the repo root::

    python demos/builder/cstr_fermenter.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import _bootstrap  # noqa: F401, E402

import numpy as np

from vlmodels.fermenter.config import FermenterBuilder
from PyOMES.core import Simulation
from PyOMES.core.boundaries import PressureReliefVent, LiquidFeed, LiquidDrain
from PyOMES.control.cv_loops import PHController, DOAgitationController


# ── Configuration ─────────────────────────────────────────────────────

USE_PH_CONTROL = True
PH_SETPOINT = 6.0

USE_DO_CONTROL = True
DO_SETPOINT_MOL_L = 1.5e-4   # target dissolved O2 (mol/L)


# ─────────────────────────────────────────────────────────────────────
# Step 1 — Sizing and operating-point parameters
# ─────────────────────────────────────────────────────────────────────
# Standard chemostat sizing: dilution rate D = Q / V_liq sets the
# residence time tau = 1/D. We size everything off the dilution rate
# so changing one number rescales the whole problem consistently.

V_total_L = 2000.0
headspace_frac = 0.20
V_liq = V_total_L * (1.0 - headspace_frac)
T_K = 305.15

MW_AcOH = 60.052
MW_yeast = 26.868

C_AcOH_feed_gL = 5.0
C_AcOH_feed = C_AcOH_feed_gL / MW_AcOH

D_per_h = 0.1                            # dilution rate (1/h)
Q_L_per_h = D_per_h * V_liq              # volumetric flow (L/h)
tau_residence = 1.0 / D_per_h            # mean residence time (h)

# Run for five residence times so transients die out.
#
# n_steps is set by stability, not just accuracy: at steady state the
# Monod-limited substrate sits at S_ss = Ks*D/(mu_max-D) -- with
# Ks=5e-3 g/L that's a near-zero concentration for the whole run, and
# the reaction sub-step is explicit Euler. At dt_h=0.025 (tau_sim/2000)
# each step overshoots the tiny available substrate, gets clamped to
# avoid negative inventory, but the paired biomass yield isn't rescaled
# to match -- biomass then compounds without bound (reaches tens of
# thousands of g/L, an AccuracyWarning fires but doesn't flag *this* as
# wrong). dt_h=0.0025 (n_steps=20000) is small enough to stay stable and
# reproduce the analytical Monod S_ss/X_ss below to ~1-2%.
tau_sim = 5.0 * tau_residence
n_steps = 20000


# ─────────────────────────────────────────────────────────────────────
# Step 2 — Build the fermenter ControlVolume
# ─────────────────────────────────────────────────────────────────────

cv = (
    FermenterBuilder()
    .vessel(V_total_L=V_total_L, headspace_frac=headspace_frac, T_K=T_K)
    .gas_feed(vvm_min=1.0, composition={"O2": 0.21, "N2": 0.79})
    .transfer_kinetic(kLa_O2=90.0, kLa_CO2_ratio=1.0)
    .transfer_species("N2", mode="kinetic", kLa_per_h=90.0)
    .chemistry()
    .organism("Yeast", balance_basis="CHO")
    .substrate("AceticAcid", mu_max=0.5, Ks=5e-3, yield_gX_gS=0.36)
    .label("cstr")
    .build()
)


# ─────────────────────────────────────────────────────────────────────
# Step 3 — Boundaries: pressure relief + feed + drain
# ─────────────────────────────────────────────────────────────────────
# Three boundaries. The pressure relief vent keeps headspace pressure
# bounded. LiquidFeed adds substrate at C_AcOH_feed at rate Q. The
# matching LiquidDrain removes broth at the same Q so V_liq stays
# constant — this is the canonical chemostat setup.

cv.boundaries.append(PressureReliefVent(P_set_atm=1.10, mode="instant"))
cv.boundaries.append(LiquidFeed(
    Q_L_per_h=Q_L_per_h,
    feed_conc_mol_L={"AceticAcid": C_AcOH_feed},
    label="substrate_feed",
))
cv.boundaries.append(LiquidDrain(Q_L_per_h=Q_L_per_h, label="broth_drain"))


# ─────────────────────────────────────────────────────────────────────
# Step 4 — Configure controllers (pH PI + DO agitation cascade)
# ─────────────────────────────────────────────────────────────────────
# DOAgitationController adjusts kLa(O2) on the gas-liquid link to
# maintain a target dissolved-O2 concentration. PHController doses
# acid or base to track a pH setpoint.
#
# Pattern 1 collapse (simulation-class decision 2): each controller's
# compute(snapshot, dt_h) returns a ControlAction directly; the
# orchestrator dispatches its flux_applied / params_changed via the
# C9 path resolver. No step_callback indirection.

controllers = []
if USE_PH_CONTROL:
    controllers.append(PHController(
        setpoint=PH_SETPOINT, Kp=0.5, Ki=0.2,
        chemical_id="H3PO4", base_chemical_id="NaOH",
        max_add_molL_hr=0.05,
    ))
if USE_DO_CONTROL:
    controllers.append(DOAgitationController(
        setpoint_mol_L=DO_SETPOINT_MOL_L,
        Kp=1e5, Ki=2e4,
    ))


# ─────────────────────────────────────────────────────────────────────
# Step 5 — Initial conditions
# ─────────────────────────────────────────────────────────────────────
# Seed the reactor at feed concentration and inoculate with a small
# biomass charge.

cv.phases["liquid"].n_mol["AceticAcid"] = C_AcOH_feed * V_liq
cv.phases["liquid"].n_mol["Yeast"] = (0.01 / MW_yeast) * V_liq


# ─────────────────────────────────────────────────────────────────────
# Step 6 — Wrap in a Simulation and run
# ─────────────────────────────────────────────────────────────────────

sim = Simulation(cvs={"main": cv}, controllers=controllers, label="cstr_demo")
result = sim.run(tau_h=tau_sim, n_steps=n_steps)


# ─────────────────────────────────────────────────────────────────────
# Step 7 — Report
# ─────────────────────────────────────────────────────────────────────
# Per-CV nested BatchResult: result.liquid_mol["main"][species].

print("=" * 60)
print("  CSTR — Continuous stirred-tank reactor (new framework)")
print("=" * 60)
print(f"  V = {V_liq:.0f} L, Q = {Q_L_per_h:.1f} L/h, D = {D_per_h:.3f} /h")
print(f"  Feed substrate: {C_AcOH_feed_gL:.1f} g/L acetic acid")
print(f"  Simulation: {tau_sim:.0f} h ({n_steps} steps)")
print(f"  Runtime: {result.runtime_s:.2f} s")
print()

cv_key = "main"
n_AcOH = result.liquid_mol[cv_key].get("AceticAcid", np.zeros(1))
n_Yeast = result.liquid_mol[cv_key].get("Yeast", np.zeros(1))
print(f"  Substrate: {n_AcOH[0]*MW_AcOH/V_liq:.3f} -> {n_AcOH[-1]*MW_AcOH/V_liq:.4f} g/L")
print(f"  Biomass:   {n_Yeast[0]*MW_yeast/V_liq:.3f} -> {n_Yeast[-1]*MW_yeast/V_liq:.4f} g/L")

pH_arr = result.pH[cv_key]
valid_pH = pH_arr[np.isfinite(pH_arr)]
if len(valid_pH) > 0:
    print(f"  pH:        {valid_pH[0]:.2f} -> {valid_pH[-1]:.2f}")

n_O2_liq = result.liquid_mol[cv_key].get("O2", np.zeros(1))
print(f"  DO:        {n_O2_liq[0]/V_liq:.2e} -> {n_O2_liq[-1]/V_liq:.2e} mol/L")

# Analytical Monod steady-state for comparison.
mu_max, Ks = 0.5, 5e-3
if mu_max > D_per_h:
    S_ss = Ks * D_per_h / (mu_max - D_per_h) * MW_AcOH
    print(f"\n  Monod S_ss (analytical): {S_ss:.4f} g/L")

if USE_PH_CONTROL:
    print(f"\n  pH control: setpoint = {PH_SETPOINT}")
if USE_DO_CONTROL:
    from PyOMES.core.transfer_models import KineticTransferModel
    _tm = cv.transfer_models.get("O2")
    kLa_final = _tm.k_transfer if isinstance(_tm, KineticTransferModel) else 0.0
    print(f"\n  DO control: setpoint = {DO_SETPOINT_MOL_L:.2e} mol/L")
    print(f"    Final kLa(O2): {kLa_final:.1f} /h")

print("\n  Done.")
