"""Chemistry-facing data structures.

This subpackage provides species declarations, a compound database
(``ChemicalRegistry``), phase-partition models, and acid-base
equilibrium sets used to construct feeds and initial conditions.
"""

from .partition import (
    PartitionModel,
    HenryEquilibrium, RaoultEquilibrium, KspEquilibrium,
    MultispeciesPartitionModel, MultispeciesVLEPartition,
)
from .species import Species, SpeciesConflictError
from .species_check import check_species_consistency
from . import common_species

__all__ = [
    "PartitionModel",
    "HenryEquilibrium",
    "RaoultEquilibrium",
    "KspEquilibrium",
    "MultispeciesPartitionModel",
    "MultispeciesVLEPartition",
    "Species",
    "SpeciesConflictError",
    "check_species_consistency",
    "common_species",
]
