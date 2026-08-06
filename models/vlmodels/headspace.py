"""Gas headspace calculations.

Provides utility functions for ideal-gas headspace thermodynamics:
total moles, mole fractions, and headspace pressure from the gas
inventory.  All calculations use (atm, L, K, mol) units to match
typical bioprocess conventions.

Also provides :func:`clamp`, a general-purpose numeric clamp used
throughout the control system.
"""
from __future__ import annotations

import math
from typing import Dict

# NOTE: Keep headspace calculations in (atm, L, K, mol) units.
# This avoids repeated conversions and matches typical bioprocess inputs.
from PyOMES.units import R_L_ATM_PER_MOL_K

def total_moles(n_gas: Dict[str, float]) -> float:
    return float(sum(max(0.0, float(v)) for v in n_gas.values()))

def mole_fractions(n_gas: Dict[str, float]) -> Dict[str, float]:
    nT = total_moles(n_gas)
    if nT <= 0.0:
        return {k: 0.0 for k in n_gas.keys()}
    return {k: max(0.0, float(v)) / nT for k, v in n_gas.items()}

def headspace_pressure_atm(n_gas: Dict[str, float], T_K: float, V_hs_L: float) -> float:
    nT = total_moles(n_gas)
    V = max(float(V_hs_L), 1e-30)
    T = max(float(T_K), 1e-9)
    return float(nT * R_L_ATM_PER_MOL_K * T / V)

def clamp(x: float, lo: float, hi: float) -> float:
    return float(lo if x < lo else hi if x > hi else x)
