# -*- coding: utf-8 -*-
"""Stirred-tank template — ControlVolume-based gas+liquid vessel construction.

:class:`StirredTankBuilder` (fluent API) and :class:`StirredTankFactory`
(config-dataclass API) both construct a :class:`~PyOMES.core.ControlVolume`
configured as a single well-mixed vessel. Currently supports gas+liquid
vessels only — a third phase (e.g. :class:`~PyOMES.core.phases.SolidPhase`)
is not yet wired through either class.
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
from .factory import StirredTankFactory
from .profiles import (
    TimeProfile, ConstantProfile, StepProfile, RampProfile,
    TableProfile, PeriodicProfile, CompositeProfile,
    ProfileTarget, ProfileSet,
    TemperatureProfileTarget, PressureSetpointTarget,
)
from .builder import StirredTankBuilder
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
    "StirredTankFactory",
    "StirredTankBuilder",
    "GrowthKinetics",
    "Monod",
    "Contois",
    "Andrews",
    "ContoisAndrews",
    "Tessier",
    "Moser",
    "Blackman",
    "DualSubstrateMonod",
    "TimeProfile",
    "ConstantProfile",
    "StepProfile",
    "RampProfile",
    "TableProfile",
    "PeriodicProfile",
    "CompositeProfile",
    "ProfileTarget",
    "ProfileSet",
    "TemperatureProfileTarget",
    "PressureSetpointTarget",
]
