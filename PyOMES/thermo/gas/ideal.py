# -*- coding: utf-8 -*-
"""IdealGasEOS: the ideal-gas equation of state (Z = 1)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

from PyOMES.units import R_L_ATM_PER_MOL_K as R

from .protocols import GasEOS


@dataclass(frozen=True)
class IdealGasEOS(GasEOS):
    """Ideal gas EOS (Z=1)."""

    def pressure_atm(self, n_tot_mol: float, *, T_K: float, V_L: float) -> float:
        V_L = max(float(V_L), 1e-30)
        return float(n_tot_mol) * float(R) * float(T_K) / V_L

    def partial_pressures_atm(self, n_gas_mol: Dict[str, float], *, T_K: float, V_L: float) -> Dict[str, float]:
        V_L = max(float(V_L), 1e-30)
        return {
            str(k): float(v) * float(R) * float(T_K) / V_L
            for k, v in (n_gas_mol or {}).items()
        }
