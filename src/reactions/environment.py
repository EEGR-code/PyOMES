# -*- coding: utf-8 -*-
"""Reaction environment — the local conditions visible to a rate law.

A :class:`ReactionEnvironment` is built by the ControlVolume from its
current phase state and the last property-solver results.  Rate law
functions receive it as their sole argument and use it to compute
reaction rates.

The environment is intentionally read-only: rate laws observe conditions
but do not mutate state.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass(frozen=True)
class ReactionEnvironment:
    """Immutable snapshot of local conditions for rate-law evaluation.

    Parameters
    ----------
    T_K : float
        Temperature (Kelvin).
    V_L : float
        Liquid volume (litres).  Used to convert between extensive (mol)
        and intensive (mol/L) quantities.
    pH : float or None
        Current pH (from the speciation property solver).
    concentrations : dict
        Species concentrations keyed by species ID.  The expected unit
        convention is **mol/L** for dissolved species and **g/L** for
        biomass / solids, matching the fermenter's historical convention.
        Rate laws should document which convention they expect.
    properties : dict
        Additional derived quantities (ionic strength, dissolved O₂,
        viscosity, etc.) keyed by descriptive strings.
    t_h : float or None
        Current simulation time (hours).  Available for time-dependent
        rate laws (e.g. lag phases).
    """

    T_K: float = 298.15
    V_L: float = 1.0
    pH: Optional[float] = None
    concentrations: Dict[str, float] = field(default_factory=dict)
    properties: Dict[str, float] = field(default_factory=dict)
    t_h: Optional[float] = None

    @property
    def has_pH(self) -> bool:
        """``True`` iff a speciation pH is available for this step.

        Rate functions that depend on pH should branch on this rather
        than relying on a sentinel default — a CV without a speciation
        solver will surface ``pH=None`` here.
        """
        return self.pH is not None

    # ── Convenience accessors ──────────────────────────────────────────

    def S(self, species_id: str, default: float = 0.0) -> float:
        """Return the concentration of a substrate/metabolite.

        Shorthand for ``self.concentrations.get(species_id, default)``.
        Intended for readable rate-law expressions::

            rate = mu_max * env.S("AceticAcid") / (Ks + env.S("AceticAcid"))
        """
        return float(self.concentrations.get(species_id, default))

    def X(self, organism_id: str = "Biomass", default: float = 0.0) -> float:
        """Return the concentration of an organism (biomass).

        Alias for :meth:`S` with a default key of ``"Biomass"``.
        """
        return float(self.concentrations.get(organism_id, default))

    def prop(self, key: str, default: float = 0.0) -> float:
        """Return a derived property value.

        Shorthand for ``self.properties.get(key, default)``.
        """
        return float(self.properties.get(key, default))
