# -*- coding: utf-8 -*-
"""ChemistryDatabase: frozen bundle of species, reactions, and ThermoFramework.

A :class:`ChemistryDatabase` is the user-facing chemistry configuration
object that a :class:`~PyOMES.core.ControlVolume` accepts.  It bundles:

- a :class:`~PyOMES.thermo.ThermoFramework` (thermodynamic conventions),
- a species dict (``{species_id: Species}``), and
- a :class:`~PyOMES.reactions.ReactionSet` of equilibrium declarations.

Databases are composed by extension, not mutation::

    from PyOMES.chemistry.databases.aqueous import AQUEOUS_DEFAULT

    MY_DB = AQUEOUS_DEFAULT.extend(
        species={"ButyricAcid": Species(id="ButyricAcid",
                                         atoms={"C":4,"H":8,"O":2},
                                         charge=0, MW=88.106)},
    )
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from .species import Species
    from .partition import PartitionModel
    from ..reactions.reaction_system import ReactionSystem
    from ..thermo.framework import ThermoFramework


@dataclass(frozen=True)
class ChemistryDatabase:
    """Frozen bundle of thermodynamic framework, species, reactions, and partition models.

    Parameters
    ----------
    thermo : ThermoFramework
        Activity model and standard-state conventions.
    species : dict
        ``{species_id: Species}`` — species declared in this database.
    reactions : ReactionSet
        Equilibrium reactions (acid-base, partitioning).
    partition_models : dict
        ``{species_id: PartitionModel}`` — phase-partition relationships
        (Henry's law, Langmuir, etc.) for gas-transferable species.
    """

    thermo: "ThermoFramework"
    species: Dict[str, "Species"] = field(default_factory=dict)
    reactions: Optional["ReactionSystem"] = None
    partition_models: Dict[str, "PartitionModel"] = field(default_factory=dict)

    def extend(
        self,
        *,
        species: Optional[Dict[str, "Species"]] = None,
        reactions: Optional[Iterable] = None,  # type: ignore[assignment]
        thermo: Optional["ThermoFramework"] = None,
        partition_models: Optional[Dict[str, "PartitionModel"]] = None,
    ) -> "ChemistryDatabase":
        """Return a new database extending this one.

        Parameters
        ----------
        species : dict, optional
            Additional ``{species_id: Species}`` entries.  Merged with
            existing species; new entries override on key collision.
        reactions : iterable, optional
            Additional reactions.  Appended to existing reaction set.
        thermo : ThermoFramework, optional
            Override the thermodynamic framework.  If omitted, inherits
            from this database.
        partition_models : dict, optional
            Additional ``{species_id: PartitionModel}`` entries.  Merged
            with existing partition_models; child entries override parent
            entries for the same species key.

        Returns
        -------
        ChemistryDatabase
            New frozen database; this one is unchanged.
        """
        from ..reactions.reaction_system import ReactionSystem

        new_thermo = thermo if thermo is not None else self.thermo

        merged_species = dict(self.species)
        if species:
            merged_species.update(species)

        existing_rxns = list(self.reactions) if self.reactions is not None else []
        new_rxns = list(reactions) if reactions is not None else []
        merged_rxns = ReactionSystem(existing_rxns + new_rxns) if (existing_rxns or new_rxns) else None

        merged_pm = dict(self.partition_models)
        if partition_models:
            merged_pm.update(partition_models)

        return ChemistryDatabase(
            thermo=new_thermo,
            species=merged_species,
            reactions=merged_rxns,
            partition_models=merged_pm,
        )
