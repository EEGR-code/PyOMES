# -*- coding: utf-8 -*-
"""GasEOS: the interface for gas-phase equations of state.

``GasEOS`` is the abstract interface both implementations satisfy:
``pressure_atm`` (total pressure from n_total) and
``partial_pressures_atm`` (per species). The two return different
quantities: ``IdealGasEOS.partial_pressures_atm`` returns partial
pressures (y_i × P); ``PengRobinsonEOS.partial_pressures_atm`` returns
fugacities (f_i = y_i × φ_i × P), the correct thermodynamic driving
force for Henry-law VLE. For ideal gas conditions (low pressure),
φ_i → 1 and the Peng-Robinson result reduces to the ideal one.
"""

from __future__ import annotations

from typing import Dict


class GasEOS:
    """Abstract gas EOS interface."""

    def pressure_atm(self, n_tot_mol: float, *, T_K: float, V_L: float) -> float:
        raise NotImplementedError

    def partial_pressures_atm(self, n_gas_mol: Dict[str, float], *, T_K: float, V_L: float) -> Dict[str, float]:
        raise NotImplementedError
