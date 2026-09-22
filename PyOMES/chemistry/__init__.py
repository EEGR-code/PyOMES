"""Chemistry-facing data structures.

This subpackage provides species declarations, phase-partition models,
acid-base equilibrium sets, and the compound registry used to construct
feeds and initial conditions.
"""

from .registry import COMPOUND_DB, resolve_compound
from .compounds import Chemical, ChemicalRegistry
from .equilibria import EquilibriumSet, EquilibriumDef
from .partition import (
    PartitionModel,
    HenryEquilibrium, RaoultEquilibrium, KspEquilibrium,
    MultispeciesPartitionModel, MultispeciesVLEPartition,
)
from .species import Species, SpeciesConflictError
from .species_check import check_species_consistency
from .database import ChemistryDatabase
from . import common_species

__all__ = [
    "COMPOUND_DB",
    "resolve_compound",
    "Chemical",
    "ChemicalRegistry",
    "EquilibriumSet",
    "EquilibriumDef",
    "PartitionModel",
    "HenryEquilibrium",
    "RaoultEquilibrium",
    "KspEquilibrium",
    "MultispeciesPartitionModel",
    "MultispeciesVLEPartition",
    "Species",
    "SpeciesConflictError",
    "check_species_consistency",
    "ChemistryDatabase",
    "common_species",
]
