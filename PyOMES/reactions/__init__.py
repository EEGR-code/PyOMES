# -*- coding: utf-8 -*-
"""Reaction framework for stoichiometrically validated chemistry models.

This package provides:

- :class:`KineticReaction` — a single integrated reaction with
  construction-time elemental balance validation and a callable rate
  law.
- :class:`EquilibriumReaction` — a single algebraic equilibrium
  constraint with an equilibrium constant and optional Van 't Hoff
  temperature parameters, routed to
  :class:`~PyOMES.chemical_equilibrium.engine.BisectionChemicalEquilibriumEngine` (single-phase) or
  :class:`~PyOMES.core.gas_liquid_link.KineticGasLiquidLink`
  (cross-phase partition declarations).
- :class:`ReactionSystem` — single attach point for all reactions on a
  :class:`~PyOMES.core.control_volume.ControlVolume`. Pre-buckets its
  reactions by type at construction; immutable post-attach.
- :class:`ReactionBuilder` — convenience factories for common reaction
  types (aerobic growth, etc.) that derive stoichiometric coefficients
  from substrate/biomass formulas and yields.
- :class:`BlackBoxReactionModel` — adapter for opaque external
  simulators with runtime mass-balance checking.
- :class:`ReactionEnvironment` — immutable snapshot of local conditions
  passed to rate laws.

The three reaction declaration classes are independent — no shared
base class. Shared validation logic lives privately in
``_shared.py``. Concrete classes implementing the
:class:`ReactionModel` protocol (i.e. carrying a ``compute_rates``
method) can be attached to a :class:`~PyOMES.core.ControlVolume`.
"""

# Core types
from .stoichiometry import StoichiometryEntry, StoichiometryError, validate_balance
from .environment import ReactionEnvironment

# Protocol
from .protocols import ReactionModel

# Concrete reaction declarations (independent classes, no shared base)
from .kinetic import KineticReaction
from .equilibrium import EquilibriumReaction
from .blackbox import BlackBoxReactionModel, FluxEntry, MassBalanceWarning, MassBalanceError

# Container and helpers
from .reaction_system import ReactionSystem
from .builder import ReactionBuilder

__all__ = [
    # Stoichiometry
    "StoichiometryEntry",
    "StoichiometryError",
    "validate_balance",
    # Environment
    "ReactionEnvironment",
    # Protocol
    "ReactionModel",
    # Reaction declarations
    "KineticReaction",
    "EquilibriumReaction",
    "BlackBoxReactionModel",
    # Container and helpers
    "ReactionSystem",
    "ReactionBuilder",
    # Black-box auxiliaries
    "FluxEntry",
    "MassBalanceWarning",
    "MassBalanceError",
]
