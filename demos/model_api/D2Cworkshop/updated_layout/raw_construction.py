#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Demo: build a fermenter Simulation from PyOMES primitives — no builder.

Mirrors [batch_fermenter.py](../../builder/batch_fermenter.py) topology
(0-D sparged batch fermenter, aerobic growth on acetic acid, PI pH
control), but constructs every object explicitly using the underlying
:mod:`PyOMES.core` types instead of :class:`FermenterBuilder`.

The chemistry is **imported** rather than redeclared — the
:class:`ReactionSystem` factories live in the sibling
[reaction_system.py](../chemistry/reaction_system.py) demo, and this
file simply wires them into a CV. That import demonstrates the
intended pattern: chemistry definitions are reusable artifacts,
written once and consumed by any number of simulation topologies.

What this demo shows that the builder hides:

- Explicit :class:`GasPhase` and :class:`LiquidPhase` construction
  with seeded ``n_mol`` dicts.
- Explicit :class:`KineticGasLiquidLink` with Henry constants and
  per-species kLa values.
- Boundary attachment (:class:`GasFeed`, :class:`PressureReliefVent`)
  by direct list append.
- :class:`ControlVolume` assembly from the parts above.
- :class:`Simulation` construction with a CV-native
  :class:`PHController`.

What this demo intentionally omits:

- Temperature-corrected Henry constants. Values are pinned to 305.15 K;
  for a real model use ``HenryPartition(H_ref=..., dlnH=...)`` with the
  appropriate van 't Hoff coefficient, or rely on the builder which reads
  from the chemistry database automatically.

Run from the repo root after ``pip install -e .``::

    python demos/model_api/D2Cworkshop/raw_construction.py
"""

# Bootstrap so the demo runs when src/ isn't on sys.path yet.
import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_repo_root / "demos"))
import _bootstrap  # noqa: F401, E402

# Sibling chemistry/ folder isn't a package; put it on sys.path so the
# bare `reaction_system` import resolves. Two demos sharing chemistry
# without a parent package — see [reaction_system.py](../chemistry/reaction_system.py)
# for the source of truth on what gets defined.
_chem_dir = Path(__file__).resolve().parents[1] / "chemistry"
sys.path.insert(0, str(_chem_dir))
import reaction_system as chem  # noqa: E402

import numpy as np  # noqa: E402

from PyOMES.core import (  # noqa: E402
    ControlVolume,
    EquilibriumTransferModel,
    GasPhase,
    KineticTransferModel,
    LiquidPhase,
    Simulation,
)
from PyOMES.core.boundaries import GasFeed, PressureReliefVent  # noqa: E402
from PyOMES.core.phases import R_L_ATM_MOL_K  # noqa: E402
from PyOMES.control.cv_loops import PHController  # noqa: E402
from PyOMES.reactions import ReactionSystem  # noqa: E402
from PyOMES.chemistry import HenryPartition  # noqa: E402


# ── Configuration ─────────────────────────────────────────────────────

T_K = 305.15
V_TOTAL_L = 2.0
HEADSPACE_FRAC = 0.20
V_GAS = V_TOTAL_L * HEADSPACE_FRAC
V_LIQ = V_TOTAL_L * (1.0 - HEADSPACE_FRAC)

PH_SETPOINT = 5.0
TAU_H = 5.0
N_STEPS = 1000

# Henry constants at 305.15 K, pinned so this demo stays self-contained.
# Real models should use HenryPartition with dlnH for temperature correction.
HENRY_MOL_L_ATM: dict = {
    "O2": 1.07e-3,
    "CO2": 2.94e-2,
    "N2": 5.27e-4,
}
KLA_PER_H = {
    "O2": 150.0,
    "CO2": 135.0,  # 0.9 × kLa(O2) — diffusivity ratio convention.
}


# ── Step 1: build the phases ──────────────────────────────────────────


def build_gas_phase() -> GasPhase:
    """Build the headspace from ideal-gas law + air composition.

    Air at 1 atm: 21% O₂, 0.04% CO₂, 79% N₂. Trace CO₂ matters
    because the carbonate ladder pulls liquid pH downward when
    sparging starts.
    """
    n_total = (1.0 * V_GAS) / (R_L_ATM_MOL_K * T_K)
    return GasPhase(
        n_mol={
            "O2":  n_total * 0.2095,
            "CO2": n_total * 0.0004,
            "N2":  n_total * 0.7901,
        },
        V_L=V_GAS,
        T_K=T_K,
    )


def build_liquid_phase(gas_phase: GasPhase) -> LiquidPhase:
    """Build the liquid with substrate, biomass, and Henry-equilibrated gases.

    Gas-phase partial pressures + Henry constants → initial dissolved
    O₂, CO₂, N₂. Substrate (acetic acid) at 1.2 g/L; biomass at
    0.1 g/L. A trace H⁺ seed gives the speciation engine something to
    overwrite on first solve.
    """
    p_atm = gas_phase.p_atm

    # Acetic acid: 1.2 g/L total acetate (MW = 60.052 g/mol).
    n_acetate_total = (1.2 / float(chem.ACETIC_ACID.MW)) * V_LIQ
    # Yeast inoculum: 0.1 g/L (MW = 24.626 g/mol).
    n_yeast = (0.1 / float(chem.YEAST.MW)) * V_LIQ

    n_mol = {
        # Substrate + organism.
        chem.ACETIC_ACID.id: n_acetate_total,
        chem.ACETATE_MINUS.id: 0.0,  # engine fills this from speciation
        chem.YEAST.id: n_yeast,
        # Dissolved gases at Henry equilibrium with the headspace.
        "O2":  HENRY_MOL_L_ATM["O2"]  * p_atm.get("O2", 0.0)  * V_LIQ,
        "CO2": HENRY_MOL_L_ATM["CO2"] * p_atm.get("CO2", 0.0) * V_LIQ,
        "N2":  HENRY_MOL_L_ATM["N2"]  * p_atm.get("N2", 0.0)  * V_LIQ,
        # n_mol["CO2"] is the dissolved molecular CO2 form (phase-agnostic
        # Species ID convention from chemistry-unification-3b). The engine
        # populates it on the first speciation solve.
        "CO2": 0.0,
        # H+ placeholder. The speciation engine writes this on every
        # solve; the seed value never persists past the first step.
        "H+": 1.0e-7 * V_LIQ,
    }
    return LiquidPhase(n_mol=n_mol, V_L=V_LIQ, T_K=T_K)


# ── Step 2: build the transfer_models dict ────────────────────────────


def build_transfer_models() -> dict:
    """Build the per-species transfer model dict.

    O₂ and CO₂ are kinetic (kLa-limited); N₂ is at instantaneous
    equilibrium.
    """
    henry = {
        sp: HenryPartition(H_ref=kH * 1000.0 / 101325.0, dlnH=0.0)
        for sp, kH in HENRY_MOL_L_ATM.items()
    }
    return {
        "O2":  KineticTransferModel(henry["O2"],  k_transfer=KLA_PER_H["O2"]),
        "CO2": KineticTransferModel(henry["CO2"], k_transfer=KLA_PER_H["CO2"]),
        "N2":  EquilibriumTransferModel(henry["N2"]),
    }


# ── Step 3: assemble the ControlVolume and Simulation ─────────────────


def build() -> Simulation:
    """Wire phases + transfer_models + chemistry + boundaries + controller into a Simulation."""
    gas = build_gas_phase()
    liquid = build_liquid_phase(gas)

    # ReactionSystem from the sibling chemistry module. The factories
    # produce one kinetic, one single-phase equilibrium, and one
    # cross-phase equilibrium reaction; the system pre-buckets them
    # at construction and the CV routes each bucket appropriately.
    rxn_system = ReactionSystem(
        [
            chem.make_aerobic_growth_on_acetate(
                mu_max_per_h=0.5, Ks_g_per_L=5e-3, yield_gX_gS=0.36,
            ),
            chem.make_acetate_dissociation(pKa=4.756),
            chem.make_co2_partition(),
        ],
        label="raw_construction_chemistry",
    )

    cv = ControlVolume(
        phases={"gas": gas, "liquid": liquid},
        transfer_models=build_transfer_models(),
        reaction_system=rxn_system,
        label="raw_construction",
    )

    # Boundaries: append after construction (mirrors the builder demos).
    cv.boundaries.append(GasFeed(
        vvm_min=1.0,
        y={"O2": 0.21, "N2": 0.79},
        P_inlet_atm=1.0,
        phase_key="gas",
        liquid_phase_key="liquid",
        label="air_sparge",
    ))
    cv.boundaries.append(PressureReliefVent(P_set_atm=1.10, mode="instant"))

    # PI pH controller — CV-native (Pattern 1). H3PO4 / NaOH dosing.
    controllers = [
        PHController(
            setpoint=PH_SETPOINT, Kp=0.5, Ki=0.0,
            chemical_id="H3PO4", base_chemical_id="NaOH",
            max_add_molL_hr=0.05,
        ),
    ]

    return Simulation(
        cvs={"main": cv},
        controllers=controllers,
        label="raw_construction_demo",
    )


# ── Step 4: run and report ────────────────────────────────────────────


def main() -> None:
    sim = build()
    result = sim.run(tau_h=TAU_H, n_steps=N_STEPS)

    cv_key = "main"
    print(
        f"Demo: raw construction (no FermenterBuilder), "
        f"t = 0 -> {TAU_H} h, {N_STEPS} steps."
    )
    print(f"  V_gas = {V_GAS:.3f} L, V_liq = {V_LIQ:.3f} L, T = {T_K:.2f} K")
    print()

    print("Initial -> final liquid inventory (mol):")
    for sp in sorted(result.liquid_mol[cv_key]):
        n0 = result.liquid_mol[cv_key][sp][0]
        nf = result.liquid_mol[cv_key][sp][-1]
        print(f"  {sp:>14}: {n0:>10.4e}  ->  {nf:>10.4e}")
    print()

    pH = result.pH[cv_key]
    pH_valid = pH[np.isfinite(pH)]
    if len(pH_valid) > 0:
        print(f"pH:  {pH_valid[0]:.3f}  ->  {pH_valid[-1]:.3f}")
    print(f"Runtime: {result.runtime_s:.3f} s wall-clock")


if __name__ == "__main__":
    main()
