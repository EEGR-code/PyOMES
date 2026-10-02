"""Chemistry-facing data structures.

This subpackage provides species declarations and phase-partition
models used to construct feeds and initial conditions.
"""

from .partition import (
    PartitionModel,
    MultispeciesPartitionModel, MultispeciesVLEPartition,
)
from .species import Species, SpeciesConflictError
from .species_check import check_species_consistency, merge_species
from .dose import as_dose, check_dose

__all__ = [
    "PartitionModel",
    "MultispeciesPartitionModel",
    "MultispeciesVLEPartition",
    "Species",
    "SpeciesConflictError",
    "check_species_consistency",
    "merge_species",
    "as_dose",
    "check_dose",
]
