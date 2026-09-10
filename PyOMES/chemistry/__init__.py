"""Chemistry-facing data structures.

This subpackage provides *typed* inputs/outputs for chemistry engines that can be
used both inside the fermenter simulation and as standalone tools.
"""

from .types import AqueousTotals, AqueousTotalsUser, AqueousEquilibrium
from .registry import (
    SALT_DISSOCIATION_MAP,
    ION_TO_ENGINE_KEY,
    normalize_ion_label,
    ion_to_engine_key,
    map_user_ions_to_engine,
    COMPOUND_DB,
    resolve_compound,
)

from .recipe import SolutionRecipe, g, mg, kg, L, mL
from .compounds import Chemical, ChemicalRegistry
from .thermo_params import (
    ThermodynamicConfig, AcidDefinition, WaterDefinition,
    validate_thermodynamics, collect_thermo_params, ThermoSnapshot,
)
from .equilibria import EquilibriumSet, EquilibriumDef
from .partition import (
    PartitionModel, HenryPartition, RaoultPartition,
    HenryEquilibrium, RaoultEquilibrium, KspEquilibrium,
    MultispeciesPartitionModel, MultispeciesVLEPartition,
)
from .species import Species, SpeciesConflictError
from .species_check import check_species_consistency
from .database import ChemistryDatabase
from . import common_species

__all__ = [
    "AqueousTotals",
    "AqueousTotalsUser",
    "AqueousEquilibrium",
    "SALT_DISSOCIATION_MAP",
    "ION_TO_ENGINE_KEY",
    "normalize_ion_label",
    "ion_to_engine_key",
    "map_user_ions_to_engine",
    "COMPOUND_DB",
    "resolve_compound",
    "SolutionRecipe",
    "g",
    "mg",
    "kg",
    "L",
    "mL",
    "Chemical",
    "ChemicalRegistry",
    "ThermodynamicConfig",
    "validate_thermodynamics",
    "collect_thermo_params",
    "ThermoSnapshot",
    "EquilibriumSet",
    "EquilibriumDef",
    "PartitionModel",
    "HenryPartition",
    "RaoultPartition",
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
