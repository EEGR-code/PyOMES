"""Chemistry-facing data structures.

This subpackage provides species declarations, phase-partition models,
and acid-base equilibrium sets used to construct feeds and initial
conditions. The standalone compound database (``ChemicalRegistry``)
lives at :mod:`PyOMES.compounds` — it has no dependency on ``Species``
or anything else here.
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
