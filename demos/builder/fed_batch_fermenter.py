#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Demo: 0-D fed-batch fermenter, aerobic growth on acetic acid.

Migrated to the new ``Simulation`` orchestrator (simulation-class C12;
this demo completed post-ship).

Tutorial walkthrough of fed-batch operation: continuous substrate
addition with no liquid removal, so the working volume grows over
time as feed accumulates. Useful when you want to keep substrate
available without diluting biomass or pushing past inhibitory
substrate concentrations.

This demo intentionally tracks more carefully than batch:

  - Both pH and DO controllers are on, so you can see the closed-loop
    behaviour as biomass and oxygen demand climb.
  - The closing report includes a substrate balance (initial + fed
    in − final = consumed), which is the classic fed-batch sanity
    check.

Run from the repo root::

    python demos/builder/fed_batch_fermenter.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import _bootstrap  # noqa: F401, E402

import numpy as np

from vlmodels.fermenter.config import FermenterBuilder
from PyOMES.core import Simulation
from PyOMES.core.boundaries import PressureReliefVent, LiquidFeed
from PyOMES.control.cv_loops import PHController, DOAgitationController


# ── Configuration ─────────────────────────────────────────────────────

USE_PH_CONTROL = True
PH_SETPOINT = 5.5

USE_DO_CONTROL = True
DO_SETPOINT_MOL_L = 1e-4


# ─────────────────────────────────────────────────────────────────────
# Step 1 — Sizing and feed parameters
# ─────────────────────────────────────────────────────────────────────
# Concentrated substrate (50 g/L) fed at 5 L/h into a 1600 L vessel —
# a classic high-cell-density setup where the feed is the only
# substrate source after the inoculum is consumed.

V_total_L = 2000.0
headspace_frac = 0.20
V_liq = V_total_L * (1.0 - headspace_frac)
T_K = 305.15

MW_AcOH = 60.052
MW_yeast = 26.868

C_AcOH_feed_gL = 50.0
C_AcOH_feed = C_AcOH_feed_gL / MW_AcOH
Q_feed = 5.0  # L/h

tau_h = 10.0
n_steps = 1000


# ─────────────────────────────────────────────────────────────────────
# Step 2 — Build the fermenter ControlVolume
# ─────────────────────────────────────────────────────────────────────

cv = (
    FermenterBuilder()
    .vessel(V_total_L=V_total_L, headspace_frac=headspace_frac, T_K=T_K)
    .gas_feed(vvm_min=1.0, composition={"O2": 0.21, "N2": 0.79})
    .transfer_kinetic(kLa_O2=150.0)
    .chemistry()
    .organism("Yeast", balance_basis="CHO")
    .substrate("AceticAcid", mu_max=0.5, Ks=5e-3, yield_gX_gS=0.36)
    .label("fed_batch")
    .build()
)


# ─────────────────────────────────────────────────────────────────────
# Step 3 — Boundaries: pressure relief + substrate feed (no drain)
# ─────────────────────────────────────────────────────────────────────
# Two boundaries: pressure relief (constant) and a liquid feed at
# Q_feed L/h. Crucially there is NO matching drain — the vessel
# fills up over the simulation, and ``cv.phases["liquid"].V_L``
# climbs as moles accumulate. Run too long and you'll overflow the
# vessel; this is normal fed-batch trade-off territory.

cv.boundaries.append(PressureReliefVent(P_set_atm=1.10, mode="instant"))
cv.boundaries.append(LiquidFeed(
    Q_L_per_h=Q_feed,
    feed_conc_mol_L={"AceticAcid": C_AcOH_feed},
    label="substrate_feed",
))


# ─────────────────────────────────────────────────────────────────────
# Step 4 — Configure controllers (pH PI + DO agitation cascade)
# ─────────────────────────────────────────────────────────────────────

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

cv.phases["liquid"].n_mol["AceticAcid"] = (1.0 / MW_AcOH) * V_liq   # 1 g/L
cv.phases["liquid"].n_mol["Yeast"]      = (0.1 / MW_yeast) * V_liq  # 0.1 g/L


# ─────────────────────────────────────────────────────────────────────
# Step 6 — Wrap in a Simulation and run
# ─────────────────────────────────────────────────────────────────────
# state-unification C4 removed chem_env / chem_env_fn entirely. The
# engine reads totals (acid totals + strong ions) straight off
# phases["liquid"].n_mol at each solve, and pH derives from
# n_mol["H+"] populated by the engine writeback.

sim = Simulation(
    cvs={"main": cv}, controllers=controllers, label="fed_batch_demo",
)
result = sim.run(tau_h=tau_h, n_steps=n_steps)


# ─────────────────────────────────────────────────────────────────────
# Step 7 — Report
# ─────────────────────────────────────────────────────────────────────

print("=" * 60)
print("  FED-BATCH — Continuous substrate feed, no drain (new framework)")
print("=" * 60)
print(f"  V = {V_liq:.0f} L (initial), Feed: {C_AcOH_feed_gL:.0f} g/L "
      f"at {Q_feed:.1f} L/h")
print(f"  Total fed: {C_AcOH_feed_gL * Q_feed * tau_h / 1000:.2f} kg "
      f"over {tau_h:.0f} h")
print(f"  Runtime: {result.runtime_s:.2f} s")
print()

cv_key = "main"
n_AcOH = result.liquid_mol[cv_key].get("AceticAcid", np.zeros(1))
n_Yeast = result.liquid_mol[cv_key].get("Yeast", np.zeros(1))
V_liq_final = cv.phases["liquid"].V_L
print(f"  Substrate: {n_AcOH[0]*MW_AcOH/V_liq:.3f} -> "
      f"{n_AcOH[-1]*MW_AcOH/V_liq_final:.3f} g/L")
print(f"  Biomass:   {n_Yeast[0]*MW_yeast/V_liq:.3f} -> "
      f"{n_Yeast[-1]*MW_yeast/V_liq_final:.3f} g/L")

pH_arr = result.pH[cv_key]
valid_pH = pH_arr[np.isfinite(pH_arr)]
if len(valid_pH) > 0:
    print(f"  pH:        {valid_pH[0]:.2f} -> {valid_pH[-1]:.2f}")

n_O2_liq = result.liquid_mol[cv_key].get("O2", np.zeros(1))
print(f"  DO:        {n_O2_liq[0]/V_liq:.2e} -> "
      f"{n_O2_liq[-1]/V_liq_final:.2e} mol/L")

# Substrate balance — the canonical fed-batch sanity check.
n_fed = C_AcOH_feed * Q_feed * tau_h
n_consumed = n_AcOH[0] + n_fed - n_AcOH[-1]
print(f"\n  Substrate balance:")
print(f"    Initial:  {n_AcOH[0] * MW_AcOH / 1000:.3f} kg")
print(f"    Fed:      {n_fed * MW_AcOH / 1000:.3f} kg")
print(f"    Final:    {n_AcOH[-1] * MW_AcOH / 1000:.3f} kg")
print(f"    Consumed: {n_consumed * MW_AcOH / 1000:.3f} kg")

if USE_PH_CONTROL:
    print(f"\n  pH control: setpoint = {PH_SETPOINT}")
if USE_DO_CONTROL:
    from PyOMES.core.transfer_models import KineticTransferModel
    _tm = cv.transfer_models.get("O2")
    kLa_final = _tm.k_transfer if isinstance(_tm, KineticTransferModel) else 0.0
    print(f"\n  DO control: setpoint = {DO_SETPOINT_MOL_L:.2e} mol/L")
    print(f"    Final kLa(O2): {kLa_final:.1f} /h")
