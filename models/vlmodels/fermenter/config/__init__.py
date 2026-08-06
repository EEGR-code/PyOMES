# -*- coding: utf-8 -*-
"""Configuration dataclasses for fermenter construction.

Each config owns one concern (vessel, gas feed, chemistry, transfer,
organism, substrate, simulation) and can be validated independently.
"""

from .configs import (
    TransferMode,
    VesselConfig,
    GasFeedConfig,
    SpeciesTransferConfig,
    TransferConfig,
    ChemistryConfig,
    OrganismConfig,
    SubstrateConfig,
    SimulationConfig,
)
from .factory import FermenterFactory
from vlmodels.fermenter.profiles import (
    TimeProfile, ConstantProfile, StepProfile, RampProfile,
    TableProfile, PeriodicProfile, CompositeProfile,
    ProfileTarget, ProfileSet,
    TemperatureProfileTarget, PressureSetpointTarget,
)
from .builder import FermenterBuilder
from .kinetics import (
    GrowthKinetics,
    Monod,
    Contois,
    Andrews,
    ContoisAndrews,
    Tessier,
    Moser,
    Blackman,
    DualSubstrateMonod,
)

__all__ = [
    "TransferMode",
    "VesselConfig",
    "GasFeedConfig",
    "SpeciesTransferConfig",
    "TransferConfig",
    "ChemistryConfig",
    "OrganismConfig",
    "SubstrateConfig",
    "SimulationConfig",
    "FermenterFactory",
    "FermenterBuilder",
    "GrowthKinetics",
    "Monod",
    "Contois",
    "Andrews",
    "ContoisAndrews",
    "Tessier",
    "Moser",
    "Blackman",
    "DualSubstrateMonod",
    # Stage 17: Time profiles
    "TimeProfile",
    "ConstantProfile",
    "StepProfile",
    "RampProfile",
    "TableProfile",
    "PeriodicProfile",
    "CompositeProfile",
    "ProfileTarget",
    "ProfileSet",
    # Stage 17c: Environmental profile targets
    "TemperatureProfileTarget",
    "PressureSetpointTarget",
]
