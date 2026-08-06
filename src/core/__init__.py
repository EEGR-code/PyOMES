"""Core abstractions: phases, interfaces, control volumes, and boundaries.

This package defines the foundational types for multi-phase,
multi-zone reactor modelling.
"""

from .phases import (
    Phase,
    GasPhase,
    LiquidPhase,
    SolidPhase,
)
from .interfaces import (
    PhaseInterface,
    TransferDiagnostics,
    AdvanceResult,
)
from .control_volume import ControlVolume, OrchestrationWarning
from .transfer_models import KineticTransferModel, EquilibriumTransferModel
from .property_calculator import PropertyCalculator
from .links import CVLink, LinkFlowRecord, AdvectiveLink, DiffusiveLink
from .gas_liquid_link import KineticGasLiquidLink
from .simulation import Simulation, RunContext
from .snapshot import (
    CVSnapshot,
    SimulationSnapshot,
    build_cv_snapshot,
    build_simulation_snapshot,
)
from .recorder import (
    Recorder, BatchRecorder, BatchResult,
    StreamingFileRecorder, load_run,
    SparseRecorder,
    SummaryRecorder, SummaryResult,
)
from .solvers import (
    SimultaneousEulerSolver,
    SimultaneousAdaptiveSolver,
    SequentialAdvanceSolver,
    StepSolver,
)
from .system_env import SystemEnv
from .system_solver import (
    SystemSolver,
    EventScheduler,
    ExplicitEulerSystemSolver,
    StrangSplittingSystemSolver,
    MultirateSystemSolver,
    ImplicitTransportSystemSolver,
    MonolithicODESolver,
    _pack_state,
    _unpack_state,
)
from .boundaries import (
    ExternalBoundary,
    ExternalFluxRecord,
    apply_boundary,
    GasFeed,
    PressureReliefVent,
    MembraneGasBoundary,
    LiquidFeed,
    LiquidDrain,
    ProportionalGasOutlet,
)

__all__ = [
    "Phase",
    "GasPhase",
    "LiquidPhase",
    "SolidPhase",
    "PhaseInterface",
    "TransferDiagnostics",
    "AdvanceResult",
    "ControlVolume",
    "OrchestrationWarning",
    "KineticTransferModel",
    "EquilibriumTransferModel",
    "PropertyCalculator",
    "CVLink",
    "LinkFlowRecord",
    "AdvectiveLink",
    "DiffusiveLink",
    "KineticGasLiquidLink",
    "Simulation",
    "RunContext",
    "CVSnapshot",
    "SimulationSnapshot",
    "build_cv_snapshot",
    "build_simulation_snapshot",
    "Recorder",
    "BatchRecorder",
    "BatchResult",
    "StreamingFileRecorder",
    "load_run",
    "SparseRecorder",
    "SummaryRecorder",
    "SummaryResult",
    "SimultaneousEulerSolver",
    "SimultaneousAdaptiveSolver",
    "SequentialAdvanceSolver",
    "StepSolver",
    "SystemEnv",
    "SystemSolver",
    "EventScheduler",
    "ExplicitEulerSystemSolver",
    "StrangSplittingSystemSolver",
    "MultirateSystemSolver",
    "ImplicitTransportSystemSolver",
    "MonolithicODESolver",
    "ExternalBoundary",
    "ExternalFluxRecord",
    "apply_boundary",
    "GasFeed",
    "PressureReliefVent",
    "MembraneGasBoundary",
    "LiquidFeed",
    "LiquidDrain",
    "ProportionalGasOutlet",
]
