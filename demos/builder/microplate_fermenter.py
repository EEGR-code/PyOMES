#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Demo: 0-D microplate well with a gas-permeable membrane.

Migrated to the new ``Simulation`` orchestrator (simulation-class C12;
this demo completed post-ship).

Tutorial walkthrough of a single 96-well microplate well sealed with
a gas-permeable membrane. An organism (here E. coli) grows on
glucose, consuming O2 and producing CO2. There's no sparging — gas
exchange happens entirely through the membrane.

This demo is the smallest topology in the set:

  - Microlitre-scale volumes (200 µL liquid, 100 µL headspace)
  - Henry-equilibrium internal transfer (thin film, no kLa limitation)
  - Membrane boundary as the only gas exchange path
  - No acid-base equilibria seeded on the speciation engine (glucose
    doesn't dissociate); the engine falls back to a pure-water
    charge balance so pH stays at ~7 throughout

It's a good demo to skim if you want to see what falls out of the
framework when chemistry is *minimal*.

Run from the repo root::

    python demos/builder/microplate_fermenter.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import _bootstrap  # noqa: F401, E402

import numpy as np

from vlmodels.fermenter.config import FermenterBuilder
from PyOMES.core import Simulation
from PyOMES.core.boundaries import MembraneGasBoundary
from PyOMES.core.phases import R_L_ATM_MOL_K
from PyOMES.control.cv_loops import PHController


# ── Configuration ─────────────────────────────────────────────────────

USE_PH_CONTROL = False       # microplate cultures are typically uncontrolled
PH_SETPOINT = 7.0            # only consulted if control is enabled


# ─────────────────────────────────────────────────────────────────────
# Step 1 — Well geometry
# ─────────────────────────────────────────────────────────────────────
# A standard 96-well plate well is ~6.35 mm in diameter, ~10 mm tall.
# We use 200 µL of liquid + 100 µL of headspace under the membrane.

well_diameter_mm = 6.35
well_area_m2 = np.pi * (well_diameter_mm / 2 / 1000) ** 2

V_liquid_uL = 200.0
V_headspace_uL = 100.0
V_liquid_L = V_liquid_uL / 1e6
V_headspace_L = V_headspace_uL / 1e6
V_total_L = V_liquid_L + V_headspace_L
headspace_frac = V_headspace_L / V_total_L
T_K = 310.15  # 37C, mammalian culture-style temperature


# ─────────────────────────────────────────────────────────────────────
# Step 2 — Membrane boundary
# ─────────────────────────────────────────────────────────────────────
# The membrane lets O2/CO2/N2 cross at species-specific permeabilities
# (mol/atm/h per m² per atm of partial-pressure difference). The
# external atmosphere is set to incubator-typical 5 % CO2.

membrane = MembraneGasBoundary(
    permeability={"O2": 0.05, "CO2": 0.25, "N2": 0.01},
    area_m2=well_area_m2,
    external_atmosphere={"O2": 0.1995, "CO2": 0.05, "N2": 0.7505},
    label="well_seal",
)


# ─────────────────────────────────────────────────────────────────────
# Step 3 — Organism + substrate parameters
# ─────────────────────────────────────────────────────────────────────
# Pass atoms / MW explicitly because E. coli isn't in the default
# chemical registry. Glucose is in the registry but we declare it
# explicitly here for symmetry.

MW_glucose = 180.156
MW_ecoli = 23.7


# ─────────────────────────────────────────────────────────────────────
# Step 4 — Build the fermenter ControlVolume
# ─────────────────────────────────────────────────────────────────────
# Note ``transfer_equilibrium()`` rather than ``transfer_kinetic()``:
# in a thin liquid film the gas-liquid kinetics are fast compared to
# everything else, so we treat dissolved-gas / headspace partitioning
# as instantaneous. ``no_gas_feed()`` skips the sparger.

cv = (
    FermenterBuilder()
    .vessel(V_total_L=V_total_L, headspace_frac=headspace_frac, T_K=T_K,
            yO2_init=0.1995, yCO2_init=0.05)
    .no_gas_feed()
    .transfer_equilibrium()
    .chemistry()
    .organism("E_coli", atoms={"C": 1, "H": 1.77, "O": 0.49, "N": 0.24},
              MW=MW_ecoli, balance_basis="CHO")
    .substrate("Glucose", atoms={"C": 6, "H": 12, "O": 6}, MW=MW_glucose,
               mu_max=0.8, Ks=0.02, yield_gX_gS=0.5)
    .label("microplate_well")
    .build()
)

cv.boundaries.append(membrane)


# ─────────────────────────────────────────────────────────────────────
# Step 5 — Configure controllers (optional pH PI)
# ─────────────────────────────────────────────────────────────────────
# pH control is off by default — microplate cultures usually are. If
# you flip USE_PH_CONTROL to True, you'll also want to install the
# carbonate equilibrium reactions on the speciation engine (see
# ``batch_fermenter.py`` for the recipe), otherwise the controller has
# no pH signal to act on.

controllers = []
if USE_PH_CONTROL:
    controllers.append(PHController(
        setpoint=PH_SETPOINT, Kp=0.1, Ki=0.05,
        chemical_id="H3PO4", base_chemical_id="NaOH",
        max_add_molL_hr=0.01,
    ))


# ─────────────────────────────────────────────────────────────────────
# Step 6 — Initial conditions
# ─────────────────────────────────────────────────────────────────────

C_glc_0 = 10.0 / MW_glucose       # 10 g/L glucose
C_ecoli_0 = 0.01 / MW_ecoli       # 0.01 g/L inoculum
cv.phases["liquid"].n_mol["Glucose"] = C_glc_0 * V_liquid_L
cv.phases["liquid"].n_mol["E_coli"]  = C_ecoli_0 * V_liquid_L


# ─────────────────────────────────────────────────────────────────────
# Step 7 — Wrap in a Simulation and run
# ─────────────────────────────────────────────────────────────────────

tau_h = 12.0
n_steps = 2000

sim = Simulation(
    cvs={"main": cv}, controllers=controllers, label="microplate_demo",
)
result = sim.run(tau_h=tau_h, n_steps=n_steps)


# ─────────────────────────────────────────────────────────────────────
# Report
# ─────────────────────────────────────────────────────────────────────

print("=" * 60)
print("  MICROPLATE WELL — E. coli on glucose, membrane aeration")
print("=" * 60)
print(f"  Well: {V_liquid_uL:.0f} uL liquid, {V_headspace_uL:.0f} uL headspace")
print(f"  Membrane: {well_area_m2*1e6:.1f} mm2, T = {T_K-273.15:.1f} C")
print(f"  Simulation: {tau_h:.0f} h ({n_steps} steps)")
print(f"  Runtime: {result.runtime_s:.3f} s")
print()

cv_key = "main"
n_glc = result.liquid_mol[cv_key].get("Glucose", np.zeros(1))
n_ecoli = result.liquid_mol[cv_key].get("E_coli", np.zeros(1))
C_glc_gL = n_glc * MW_glucose / V_liquid_L
C_ecoli_gL = n_ecoli * MW_ecoli / V_liquid_L

print(f"  Glucose:  {C_glc_gL[0]:.3f} -> {C_glc_gL[-1]:.3f} g/L")
print(f"  Biomass:  {C_ecoli_gL[0]:.4f} -> {C_ecoli_gL[-1]:.4f} g/L")

pH_arr = result.pH[cv_key]
valid_pH = pH_arr[np.isfinite(pH_arr)]
if len(valid_pH) > 0:
    print(f"  pH:       {valid_pH[0]:.2f} -> {valid_pH[-1]:.2f}")

# Headspace O2 partial pressure.
n_O2_gas = result.gas_mol[cv_key].get("O2", np.zeros(1))
V_hs = cv.phases["gas"].V_L
if V_hs > 0:
    pO2_0 = n_O2_gas[0] * R_L_ATM_MOL_K * T_K / V_hs
    pO2_f = n_O2_gas[-1] * R_L_ATM_MOL_K * T_K / V_hs
    print(f"  p(O2):    {pO2_0:.4f} -> {pO2_f:.4f} atm")

n_O2_liq = result.liquid_mol[cv_key].get("O2", np.zeros(1))
print(f"  DO:       {n_O2_liq[0]/V_liquid_L:.2e} -> "
      f"{n_O2_liq[-1]/V_liquid_L:.2e} mol/L")

if USE_PH_CONTROL:
    print(f"\n  pH control: setpoint = {PH_SETPOINT}")
else:
    print("\n  pH control: OFF (uncontrolled microplate culture)")

print("\n  Done.")
