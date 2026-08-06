"""Component-total inventory accounting for the standalone pathway."""

from __future__ import annotations

from typing import Mapping

from .components import ComponentBasis, ComponentSpecies


def component_totals(
    basis: ComponentBasis,
    amounts_mol: Mapping[ComponentSpecies, float],
) -> dict[str, float]:
    """Sum a species inventory into conserved-component totals in moles."""
    totals = {component_id: 0.0 for component_id in basis.ids}
    for species, amount in amounts_mol.items():
        if amount < 0:
            raise ValueError(f"Negative inventory for {species.species_id!r}")
        for component_id, coefficient in zip(basis.ids, basis.vector(species)):
            totals[component_id] += float(amount) * coefficient
    return totals
