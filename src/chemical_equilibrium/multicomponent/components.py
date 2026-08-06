"""Conserved-component bookkeeping for the standalone equilibrium pathway.

Unlike the legacy NR tableau, this module does not infer one conserved quantity
from a connected reaction graph.  Each species explicitly records how many
moles of every independently conserved component it contains.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping


@dataclass(frozen=True)
class ComponentSpecies:
    """A charged aqueous or solid species expressed in conserved components."""

    species_id: str
    charge: int
    components: Mapping[str, float]

    def __hash__(self) -> int:
        """Allow species declarations to be reaction-stoichiometry keys.

        Species identity is its stable ID; tableau construction separately
        rejects duplicate IDs, preventing ambiguous declarations.
        """
        return hash(self.species_id)

    def coefficient(self, component_id: str) -> float:
        """Return this species' stoichiometric coefficient for a component."""
        return float(self.components.get(component_id, 0.0))


class ComponentBasis:
    """Ordered independent conserved-component basis.

    The basis excludes proton and water, which are handled respectively by
    electroneutrality and the solvent-activity convention.
    """

    def __init__(self, component_ids: tuple[str, ...] | list[str]):
        ids = tuple(str(component_id) for component_id in component_ids)
        if not ids:
            raise ValueError("ComponentBasis requires at least one component")
        if len(set(ids)) != len(ids):
            raise ValueError("ComponentBasis component IDs must be unique")
        self._ids = ids
        self._index = {component_id: index for index, component_id in enumerate(ids)}

    @property
    def ids(self) -> tuple[str, ...]:
        return self._ids

    def vector(self, species: ComponentSpecies) -> tuple[float, ...]:
        """Return a species' component stoichiometry in basis order.

        Unknown component names are rejected.  This prevents silently losing a
        component when a new complex or mineral is declared.
        """
        unknown = set(species.components) - set(self._ids)
        if unknown:
            raise ValueError(
                f"{species.species_id!r} references components not in the basis: "
                f"{sorted(unknown)}"
            )
        return tuple(species.coefficient(component_id) for component_id in self._ids)
