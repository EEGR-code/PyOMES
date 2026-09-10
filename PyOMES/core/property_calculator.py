# -*- coding: utf-8 -*-
"""PropertyCalculator protocol — scalar derived properties on a phase.

A :class:`PropertyCalculator` computes one scalar derived property
(viscosity, density, heat capacity, compressibility, …) from the
current phase state. Results land on ``phase.properties`` keyed by
the calculator's ``key``; kinetic rate laws can then read them via
``env.prop("viscosity")``.

Introduced in state-unification C5 to replace the deleted
``PropertySolver`` protocol with a narrower, focused surface.

Distinctions:

- **Speciation is NOT a PropertyCalculator.** Acid-base equilibria
  are state-completion handled by
  :class:`~PyOMES.chemical_equilibrium.engine.BisectionChemicalEquilibriumEngine` via
  ``cv.reaction_system.engine.solve(phases=...)`` followed by an explicit
  ``result.apply_to_phases(phases)`` commit. The commit writes back to
  ``phase.n_mol``, not to ``phase.properties``.
- **PropertyCalculator is for everything else.** Scalar physical
  properties whose values depend on phase composition / temperature
  / pressure but don't feed back into mass balance.

Lifecycle:

- Attach via ``cv.property_calculators.append(my_calc)`` (or pass
  ``property_calculators=[...]`` to the ControlVolume constructor).
- Invoked **once per ``advance()`` step**, before kinetic reactions,
  so the latest values are visible to rate laws.
- No ``update_after`` opt-in in v1. If a concrete need for post-step
  re-evaluation emerges, the protocol can be extended additively.

Example
-------
>>> class ViscosityCalculator:
...     key = "viscosity"
...     phase_key = "liquid"
...     def compute(self, phase, T_K, P_atm):
...         X_gL = phase.concentrations_g_L.get("Biomass", 0.0)
...         return 1e-3 * (1.0 + 0.05 * X_gL)  # Pa·s
>>> cv.property_calculators.append(ViscosityCalculator())
>>> cv.advance(dt_h=0.01, t_h=0.0)
>>> cv.phases["liquid"].properties["viscosity"]  # set by the call
0.00125
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from .phases import Phase


@runtime_checkable
class PropertyCalculator(Protocol):
    """Computes one scalar derived property for one phase.

    Implementations supply:

    - ``key``: the string identifier used to store the result on
      ``phase.properties[key]`` and read it via ``env.prop(key)``.
      Conventionally lowercase snake_case (``"viscosity"``,
      ``"density_g_L"``, ``"heat_capacity_J_K_kg"``).
    - ``phase_key``: which phase this calculator targets
      (``"liquid"``, ``"gas"``, ``"solid"``, or any custom key).
      The runner skips the calculator silently if that phase is absent.
    - ``compute(phase, T_K, P_atm) -> float``: pure function of the
      phase state. Must not mutate ``phase``.
    """

    key: str
    phase_key: str

    def compute(self, phase: Phase, T_K: float, P_atm: float) -> float:
        """Compute the scalar property value.

        Parameters
        ----------
        phase : Phase
            The phase to read state from. Read-only access — do not
            mutate.
        T_K : float
            Phase temperature (Kelvin).
        P_atm : float
            System pressure (atm). For liquid-phase calculators this
            is typically the headspace pressure; for gas-phase
            calculators it is the partial pressure in the phase.

        Returns
        -------
        float
            The computed property value, stored on
            ``phase.properties[self.key]`` after the call.
        """
        ...
