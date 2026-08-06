"""Fermenter framework — building blocks for bioprocess simulation."""

from __future__ import annotations

__version__ = "0.12.5"

from .stream_adapter import FeedState
from .chemistry.compounds import ChemicalRegistry, Chemical

from .control import (
    PressureReliefController,
    InstantPressureReliefController,
    SmoothPressureReliefController,
    PHController,
    DOAgitationController,
    DOController,
    DOCascadeController,
)
from . import units
from .config import config, Config, WarningConfig
from .thermo import ThermoFramework
from .monitoring import (
    AccuracyMonitor,
    AccuracyWarning,
    print_accuracy_summary,
    reset_accuracy_summary,
    ConservationMonitor,
    ConservationWarning,
    print_conservation_summary,
    reset_conservation_summary,
)

from .properties import (
    BrothState,
    ViscosityResult,
    ViscosityModel,
    WaterViscosity,
    ArrheniusBiomassViscosity,
    PowerLawBiomassViscosity,
)

__all__ = [
    'FeedState',
    'ChemicalRegistry',
    'Chemical',
    'PressureReliefController',
    'InstantPressureReliefController',
    'SmoothPressureReliefController',
    'PHController',
    'DOAgitationController',
    'DOController',
    'DOCascadeController',
    'BrothState',
    'ViscosityResult',
    'ViscosityModel',
    'WaterViscosity',
    'ArrheniusBiomassViscosity',
    'PowerLawBiomassViscosity',
    'units',
    'config',
    'Config',
    'WarningConfig',
    'AccuracyMonitor',
    'AccuracyWarning',
    'print_accuracy_summary',
    'reset_accuracy_summary',
    'ConservationMonitor',
    'ConservationWarning',
    'print_conservation_summary',
    'reset_conservation_summary',
]
