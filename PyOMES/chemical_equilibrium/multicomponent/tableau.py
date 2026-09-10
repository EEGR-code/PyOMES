"""Validated component/species tableau for the standalone pathway."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from .components import ComponentBasis, ComponentSpecies
from .reactions import ComponentEquilibriumReaction


@dataclass(frozen=True)
class MultiComponentTableau:
    """Species component matrix plus validated equilibrium declarations."""

    basis: ComponentBasis
    species: tuple[ComponentSpecies, ...]
    reactions: tuple[ComponentEquilibriumReaction, ...]
    component_matrix: tuple[tuple[float, ...], ...]

    @classmethod
    def build(
        cls,
        basis: ComponentBasis,
        species: Iterable[ComponentSpecies],
        reactions: Iterable[ComponentEquilibriumReaction],
    ) -> "MultiComponentTableau":
        species_tuple = tuple(species)
        if len({species_item.species_id for species_item in species_tuple}) != len(species_tuple):
            raise ValueError("MultiComponentTableau species IDs must be unique")
        reactions_tuple = tuple(reactions)
        for reaction in reactions_tuple:
            reaction.validate_component_balance(basis)
        return cls(
            basis=basis,
            species=species_tuple,
            reactions=reactions_tuple,
            component_matrix=tuple(basis.vector(item) for item in species_tuple),
        )
