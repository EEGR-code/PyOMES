"""Chemistry-facing data structures.

This subpackage provides species declarations and phase-partition
models used to construct feeds and initial conditions. The standalone compound database (``ChemicalRegistry``)
lives at :mod:`PyOMES.compounds` — it has no dependency on ``Species``
or anything else here.
"""

from .partition import (
    PartitionModel,
    MultispeciesPartitionModel, MultispeciesVLEPartition,
)
from .species import Species, SpeciesConflictError
from .species_check import check_species_consistency, merge_species
from . import common_species

__all__ = [
    "PartitionModel",
    "MultispeciesPartitionModel",
    "MultispeciesVLEPartition",
    "Species",
    "SpeciesConflictError",
    "check_species_consistency",
    "merge_species",
    "common_species",
]
