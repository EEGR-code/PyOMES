# -*- coding: utf-8 -*-
"""Protocol definition for reaction models.

Any object that implements :class:`ReactionModel` can be attached to a
ControlVolume.  The CV calls ``compute_rates(env)`` and applies the
returned source terms to its phases.

Implementations carrying ``compute_rates`` in this package:

- :class:`~PyOMES.reactions.kinetic.KineticReaction` — single
  stoichiometric reaction with construction-time elemental balance
  validation and a callable rate law.
- :class:`~PyOMES.reactions.reaction_system.ReactionSystem` —
  pre-bucketed container of reaction declarations attached to a
  ControlVolume; ``compute_rates`` sums kinetic + black-box
  constituents (equilibrium constituents are routed elsewhere via
  the system's internal pre-bucketing, not through this protocol).
- :class:`~PyOMES.reactions.blackbox.BlackBoxReactionModel` — adapter
  for opaque external simulators with runtime balance checking.

:class:`~PyOMES.reactions.equilibrium.EquilibriumReaction` deliberately
does **not** implement this protocol — equilibrium reactions are
algebraic constraints, not rate-producers.
"""

from __future__ import annotations

from typing import Dict, Protocol, runtime_checkable

from .environment import ReactionEnvironment


@runtime_checkable
class ReactionModel(Protocol):
    """Protocol for anything that computes reaction source terms.

    Implementations must provide ``compute_rates(env)`` which returns
    a nested dict ``{phase_key: {species_id: mol_per_h}}``.

    Positive values = production, negative values = consumption.
    """

    def compute_rates(
        self, env: ReactionEnvironment
    ) -> Dict[str, Dict[str, float]]:
        """Evaluate reaction rates at the given conditions.

        Parameters
        ----------
        env : ReactionEnvironment
            Current local conditions (temperature, concentrations, pH, etc.).

        Returns
        -------
        dict of dict
            ``{phase_key: {species_id: rate_mol_per_h}}``.
            Phase keys must match ControlVolume phase keys (typically
            ``"liquid"`` and/or ``"gas"``).
        """
        ...
