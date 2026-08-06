"""Gas-phase EOS and Henry-law VLE helpers.

This module isolates the gas-phase model from the coupled equilibrator so that
future non-ideal EOS options (e.g., Redlich-Kwong) can be introduced without
rewriting equilibrium coupling logic.

Current implementation: ideal gas only.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict

from PyOMES.units import R_L_ATM_PER_MOL_K


class GasEOS:
    """Abstract gas EOS interface."""

    def pressure_atm(self, n_tot_mol: float, *, T_K: float, V_L: float) -> float:
        raise NotImplementedError

    def partial_pressures_atm(self, n_gas_mol: Dict[str, float], *, T_K: float, V_L: float) -> Dict[str, float]:
        raise NotImplementedError


@dataclass(frozen=True)
class IdealGasEOS(GasEOS):
    """Ideal gas EOS (Z=1)."""

    def pressure_atm(self, n_tot_mol: float, *, T_K: float, V_L: float) -> float:
        V_L = max(float(V_L), 1e-30)
        return float(n_tot_mol) * float(R_L_ATM_PER_MOL_K) * float(T_K) / V_L

    def partial_pressures_atm(self, n_gas_mol: Dict[str, float], *, T_K: float, V_L: float) -> Dict[str, float]:
        V_L = max(float(V_L), 1e-30)
        return {
            str(k): float(v) * float(R_L_ATM_PER_MOL_K) * float(T_K) / V_L
            for k, v in (n_gas_mol or {}).items()
        }


@dataclass(frozen=True)
class HenryIdealVLE:
    """Henry-law + ideal-gas VLE relationships for sparingly soluble gases."""

    eos: GasEOS
    henry_mol_per_L_atm: Dict[str, float]  # species -> kH

    def analytic_partition_n_gas(self, total_mol: float, *, species: str, T_K: float, V_hs_L: float, V_liq_L: float) -> float:
        """Analytic split for species where C* = kH * p and p = nRT/V."""
        kH = float(self.henry_mol_per_L_atm[species])
        V_hs_L = max(float(V_hs_L), 1e-30)
        coeff = 1.0 + (kH * float(R_L_ATM_PER_MOL_K) * float(T_K) * (float(V_liq_L) / V_hs_L))
        if coeff <= 0.0:
            return float(max(0.0, total_mol))
        return float(max(0.0, total_mol / coeff))

    def p_total_atm(self, n_gas_mol: Dict[str, float], *, T_K: float, V_hs_L: float) -> float:
        n_tot = sum(float(v) for v in (n_gas_mol or {}).values())
        return self.eos.pressure_atm(n_tot, T_K=T_K, V_L=V_hs_L)

    def p_i_atm(self, n_i_mol: float, *, T_K: float, V_hs_L: float) -> float:
        return self.eos.pressure_atm(float(n_i_mol), T_K=T_K, V_L=V_hs_L)

    def Cstar_mol_L(self, *, species: str, p_i_atm: float) -> float:
        return float(self.henry_mol_per_L_atm[species]) * float(p_i_atm)
