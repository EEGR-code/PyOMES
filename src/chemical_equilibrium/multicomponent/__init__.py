"""Standalone multi-component aqueous-equilibrium pathway.

This package is intentionally opt-in.  It does not modify the existing
``BisectionChemicalEquilibriumEngine`` or Newton-Raphson pathways while its component-basis,
complexation, and precipitation capabilities are developed and validated.
"""

from .components import ComponentBasis, ComponentSpecies
from .reactions import ComponentEquilibriumReaction
from .tableau import MultiComponentTableau
from .network import ComponentReactionNetwork
from .inventory import component_totals
from .solver import IdealMassActionSolver, MultiComponentProblem
from .residuals import mass_action_residuals, component_balance_residuals, charge_balance_residual

__all__ = [
    "ComponentBasis", "ComponentSpecies", "ComponentEquilibriumReaction",
    "MultiComponentTableau",
    "ComponentReactionNetwork",
    "component_totals",
    "IdealMassActionSolver", "MultiComponentProblem",
    "mass_action_residuals",
    "component_balance_residuals", "charge_balance_residual",
]
