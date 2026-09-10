"""Equilibrium-reaction declarations for the standalone component pathway."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from .components import ComponentBasis, ComponentSpecies


@dataclass(frozen=True)
class ComponentEquilibriumReaction:
    """A mass-action equilibrium reaction using :class:`ComponentSpecies`.

    ``stoichiometry`` maps species to conventional signed coefficients:
    reactants are negative and products are positive.
    """

    label: str
    stoichiometry: Mapping[ComponentSpecies, float]
    log_k: float

    def component_residual(self, basis: ComponentBasis) -> tuple[float, ...]:
        """Return the conserved-component residual of the reaction."""
        residual = [0.0] * len(basis.ids)
        for species, coefficient in self.stoichiometry.items():
            for index, amount in enumerate(basis.vector(species)):
                residual[index] += float(coefficient) * amount
        return tuple(residual)

    def validate_component_balance(self, basis: ComponentBasis, *, tol: float = 1e-12) -> None:
        """Raise when the declared reaction loses or creates a component."""
        residual = self.component_residual(basis)
        if any(abs(value) > tol for value in residual):
            detail = dict(zip(basis.ids, residual))
            raise ValueError(
                f"Reaction {self.label!r} is not component-balanced: {detail}"
            )
