#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Demo: 0-D batch fermenter, aerobic growth on acetic acid (new framework).

Migrated to the new ``Simulation`` orchestrator (simulation-class C12).
Walks through the typical fermenter pattern under the new API:

    1. Build the fermenter ControlVolume + Simulation via
       ``FermenterBuilder.build_simulation_and_run()``.
    2. Wire pH control via the CV-native ``PHController`` from
       ``PyOMES.control.cv_loops``.
    3. Read the per-CV time series from the returned ``BatchResult``.

Run from the repo root after ``pip install -e .``::

    python demos/builder/batch_fermenter.py

Set ``USE_PH_CONTROL = False`` near the top to see what happens
without intervention.
"""

# Bootstrap so the demo runs when src/ isn't on sys.path yet.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import _bootstrap  # noqa: F401, E402

import numpy as np

from vlmodels.fermenter.config import FermenterBuilder
from PyOMES.control.cv_loops import PHController


# ── Configuration ─────────────────────────────────────────────────────

USE_PH_CONTROL = True       # set False to see the unbuffered drift
PH_SETPOINT = 5.0           # PI controller target when control is on
TAU_H = 5.0                 # simulation duration (hours)
N_STEPS = 1000              # output grid resolution


# ── Build the Simulation ──────────────────────────────────────────────

def build() -> FermenterBuilder:
    """Build a FermenterBuilder configured for the demo.

    Vessel: 2 L, 305.15 K. Sparged with air at 1 vvm. Kinetic
    gas-liquid transfer with kLa(O2)=150/h. Yeast growing on
    acetic acid; pH-active equilibria built from the chemistry
    preset.
    """
    builder = (
        FermenterBuilder()
        .vessel(V_total_L=2000, T_K=305.15)
        .gas_feed(vvm_min=1.0, composition={"O2": 0.21, "N2": 0.79})
        .transfer_kinetic(kLa_O2=150.0)
        .chemistry()
        .organism("Yeast")
        .substrate("AceticAcid", mu_max=0.5, Ks=5e-3, yield_gX_gS=0.36)
        .label("batch_demo")
    )
    if USE_PH_CONTROL:
        builder = builder.controller(
            PHController(setpoint=PH_SETPOINT, Kp=0.5, Ki=0.0)
        )
    return builder


# ── Run and report ────────────────────────────────────────────────────

def main() -> None:
    result = build().build_simulation_and_run(tau_h=TAU_H, n_steps=N_STEPS)
    print(
        f"Demo: batch fermenter, t = 0 -> {TAU_H} h, "
        f"{N_STEPS} steps, pH control = {USE_PH_CONTROL}"
    )
    # Per-CV nested BatchResult: result.gas_mol["main"][species]
    cv_key = "main"
    print("\nFinal liquid inventory (mol):")
    for sp in sorted(result.liquid_mol[cv_key]):
        print(f"  {sp:>16}: {result.liquid_mol[cv_key][sp][-1]:.4f}")
    pH_final = result.pH[cv_key][-1]
    pH_str = f"{pH_final:.3f}" if np.isfinite(pH_final) else "n/a"
    print(f"\nFinal pH: {pH_str}")
    print(f"Runtime: {result.runtime_s:.3f} s wall-clock")


if __name__ == "__main__":
    main()
