"""Mass-action residuals for the standalone ideal-activity solver."""

from __future__ import annotations

import numpy as np

from .network import ComponentReactionNetwork
from .tableau import MultiComponentTableau


def mass_action_residuals(
    network: ComponentReactionNetwork,
    log10_concentrations: np.ndarray,
) -> np.ndarray:
    """Return ``log10(Q) - log10(K)`` for each declared equilibrium.

    Concentrations are represented in log space, so the reaction quotient is a
    matrix multiplication and never underflows for dilute species.
    """
    values = np.asarray(log10_concentrations, dtype=float)
    if values.shape != (len(network.species),):
        raise ValueError("log10_concentrations has incompatible shape")
    log_q = network.stoichiometric_matrix.T @ values
    log_k = np.asarray([reaction.log_k for reaction in network.reactions])
    return log_q - log_k


def component_balance_residuals(
    tableau: MultiComponentTableau,
    log10_concentrations: np.ndarray,
    component_totals_mol_L: dict[str, float],
) -> np.ndarray:
    """Return calculated-minus-specified component totals in mol/L."""
    values = np.asarray(log10_concentrations, dtype=float)
    if values.shape != (len(tableau.species),):
        raise ValueError("log10_concentrations has incompatible shape")
    concentrations = np.power(10.0, values)
    matrix = np.asarray(tableau.component_matrix, dtype=float)
    calculated = matrix.T @ concentrations
    expected = np.asarray([component_totals_mol_L.get(key, 0.0) for key in tableau.basis.ids])
    return calculated - expected


def charge_balance_residual(
    tableau: MultiComponentTableau,
    log10_concentrations: np.ndarray,
    strong_charge_mol_L: float = 0.0,
) -> float:
    """Return net aqueous charge in mol/L equivalents."""
    values = np.asarray(log10_concentrations, dtype=float)
    charges = np.asarray([species.charge for species in tableau.species], dtype=float)
    return float(charges @ np.power(10.0, values) + strong_charge_mol_L)
