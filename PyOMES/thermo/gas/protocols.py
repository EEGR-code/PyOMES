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

from typing import Dict, Protocol, runtime_checkable


@runtime_checkable
class GasEOS(Protocol):
    """Gas-phase equation of state: total and per-species pressure from amounts.

    Implementations satisfy this structurally; they do not subclass it.
    ``partial_pressures_atm`` does not return the same quantity in every
    implementation: ``IdealGasEOS`` returns partial pressures (y_i × P) and
    ``PengRobinsonEOS`` returns fugacities (y_i × φ_i × P).
    """

    def pressure_atm(self, n_tot_mol: float, *, T_K: float, V_L: float) -> float:
        """Total pressure (atm) of ``n_tot_mol`` in ``V_L`` litres at ``T_K``."""
        ...

    def partial_pressures_atm(self, n_gas_mol: Dict[str, float], *, T_K: float, V_L: float) -> Dict[str, float]:
        """Per-species pressure (atm): partial pressure or fugacity, by implementation."""
        ...
