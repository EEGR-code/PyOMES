"""Reaction-matrix validation for the standalone multi-component pathway."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .components import ComponentSpecies
from .reactions import ComponentEquilibriumReaction


@dataclass(frozen=True)
class ComponentReactionNetwork:
    """Species/reaction stoichiometric matrix with rank validation."""

    species: tuple[ComponentSpecies, ...]
    reactions: tuple[ComponentEquilibriumReaction, ...]
    stoichiometric_matrix: np.ndarray

    @classmethod
    def build(
        cls,
        species: Iterable[ComponentSpecies],
        reactions: Iterable[ComponentEquilibriumReaction],
    ) -> "ComponentReactionNetwork":
        species_tuple = tuple(species)
        reactions_tuple = tuple(reactions)
        index = {item: position for position, item in enumerate(species_tuple)}
        if len(index) != len(species_tuple):
            raise ValueError("Network species declarations must be unique")
        matrix = np.zeros((len(species_tuple), len(reactions_tuple)))
        for column, reaction in enumerate(reactions_tuple):
            for item, coefficient in reaction.stoichiometry.items():
                if item not in index:
                    raise ValueError(
                        f"Reaction {reaction.label!r} references undeclared species "
                        f"{item.species_id!r}"
                    )
                matrix[index[item], column] = coefficient
        if reactions_tuple and np.linalg.matrix_rank(matrix) != len(reactions_tuple):
            raise ValueError("Network reactions are linearly dependent")
        return cls(species_tuple, reactions_tuple, matrix)
