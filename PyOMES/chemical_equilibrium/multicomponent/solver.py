"""Standalone multi-component equilibrium-solver public API.

No existing solver imports or constructs this class; it is intentionally opt-in.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from .components import ComponentBasis, ComponentSpecies
from .network import ComponentReactionNetwork
from .tableau import MultiComponentTableau
from .residuals import mass_action_residuals, component_balance_residuals, charge_balance_residual


@dataclass(frozen=True)
class MultiComponentProblem:
    """Inputs for a standalone aqueous equilibrium calculation."""

    basis: ComponentBasis
    tableau: MultiComponentTableau
    network: ComponentReactionNetwork
    component_totals_mol_L: Mapping[str, float]
    strong_charge_mol_L: float = 0.0
    temperature_K: float = 298.15

    def validate(self) -> None:
        unknown = set(self.component_totals_mol_L) - set(self.basis.ids)
        if unknown:
            raise ValueError(f"Unknown component totals: {sorted(unknown)}")
        if self.temperature_K <= 0:
            raise ValueError("temperature_K must be positive")
        if any(value < 0 for value in self.component_totals_mol_L.values()):
            raise ValueError("Component totals must be non-negative")


class IdealMassActionSolver:
    """Reserved standalone ideal-activity mass-action solver.

    The residual/Jacobian implementation is added next; this API is deliberately
    separate from ``NRChemicalEquilibriumEngine``.
    """

    def solve(self, problem: MultiComponentProblem) -> Mapping[str, float]:
        problem.validate()
        from scipy.optimize import least_squares
        import numpy as np

        species_count = len(problem.network.species)
        if problem.tableau.species != problem.network.species:
            raise ValueError("Problem tableau and reaction network species must match")
        scale = max(1e-12, *problem.component_totals_mol_L.values(), abs(problem.strong_charge_mol_L))
        guess = np.full(species_count, np.log10(scale / max(1, species_count)))

        def residual(log_c: np.ndarray) -> np.ndarray:
            mass_action = mass_action_residuals(problem.network, log_c)
            component = component_balance_residuals(
                problem.tableau, log_c, dict(problem.component_totals_mol_L),
            ) / scale
            charge = charge_balance_residual(
                problem.tableau, log_c, problem.strong_charge_mol_L) / scale
            return np.concatenate((mass_action, component, [charge]))

        result = least_squares(residual, guess, max_nfev=2000, xtol=1e-12, ftol=1e-12, gtol=1e-12)
        if not result.success or np.max(np.abs(residual(result.x))) > 1e-8:
            raise RuntimeError(f"Ideal mass-action solve failed: {result.message}")
        return {species.species_id: float(10.0 ** value) for species, value in zip(problem.network.species, result.x)}
