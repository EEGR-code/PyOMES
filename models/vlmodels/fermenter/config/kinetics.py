# -*- coding: utf-8 -*-
"""Pluggable growth kinetics models for substrate-biomass systems.

Each kinetics class encapsulates a specific growth rate law (μ as a
function of substrate and biomass concentrations).  The factory and
builder use these to construct rate functions for
:meth:`ReactionBuilder.aerobic_growth`.

The kinetics object computes only the specific growth rate μ (1/h).
The factory handles the conversion to substrate consumption rate
(mol/h), yield, and molecular weight bookkeeping.

Usage with the builder
----------------------
>>> from PyOMES.config.kinetics import Monod, Contois, Andrews
>>>
>>> FermenterBuilder()
...     .substrate("Glucose", mu_max=0.8, Ks=0.02, yield_gX_gS=0.5,
...                kinetics=Monod())                          # default
...     .substrate("Glucose", yield_gX_gS=0.5,
...                kinetics=Contois(mu_max=0.113, Ks=0.085))  # all params on kinetics
...     .substrate("Glucose", yield_gX_gS=0.5,
...                kinetics=Andrews(mu_max=0.8, Ks=0.02, Ki=50.0))

Usage standalone
----------------
>>> kin = Contois(mu_max=0.113, Ks=0.085)
>>> rate_fn = kin.make_rate_fn(
...     organism_id="E_coli", substrate_id="Glucose",
...     MW_organism=23.7, MW_substrate=180.156, yield_gX_gS=0.5)
>>> rxn = ReactionBuilder.aerobic_growth(..., rate_fn=rate_fn)
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Callable, Dict, Optional, Protocol, runtime_checkable


# ════════════════════════════════════════════════════════════════════════
#  Protocol
# ════════════════════════════════════════════════════════════════════════

@runtime_checkable
class GrowthKinetics(Protocol):
    """Protocol for a growth rate law.

    Implementations must provide:
    - ``mu(S_gL, X_gL) -> float``: specific growth rate (1/h)
    - ``make_rate_fn(...)`` : build a rate function for ReactionBuilder
    - ``label``: human-readable name
    """

    @property
    def label(self) -> str: ...

    def mu(self, S_gL: float, X_gL: float) -> float:
        """Compute specific growth rate.

        Parameters
        ----------
        S_gL : float
            Substrate concentration (g/L).
        X_gL : float
            Biomass concentration (g/L).

        Returns
        -------
        float
            Specific growth rate μ (1/h).
        """
        ...

    def make_rate_fn(
        self,
        organism_id: str,
        substrate_id: str,
        MW_organism: float,
        MW_substrate: float,
        yield_gX_gS: float,
    ) -> Callable:
        """Build a rate function compatible with ReactionBuilder.

        Returns a callable ``rate_fn(env) -> float`` that returns
        the extensive substrate consumption rate in mol/h.
        """
        ...


# ════════════════════════════════════════════════════════════════════════
#  Base class with shared rate_fn construction
# ════════════════════════════════════════════════════════════════════════

class _KineticsBase:
    """Shared logic for converting μ(S,X) to mol substrate consumed per hour."""

    def make_rate_fn(
        self,
        organism_id: str,
        substrate_id: str,
        MW_organism: float,
        MW_substrate: float,
        yield_gX_gS: float,
    ) -> Callable:
        """Build a rate function for ReactionBuilder.aerobic_growth.

        The rate function:
          1. Reads concentrations from the ReactionEnvironment
          2. Converts mol/L to g/L
          3. Calls self.mu(S_gL, X_gL) for the specific growth rate
          4. Converts to mol substrate consumed per hour:
             rate = (μ / Y) × X_gL / MW_substrate × V_L

        Parameters
        ----------
        organism_id : str
            Species ID for the organism in the environment.
        substrate_id : str
            Species ID for the substrate in the environment.
        MW_organism : float
            Molecular weight of the organism (g/mol).
        MW_substrate : float
            Molecular weight of the substrate (g/mol).
        yield_gX_gS : float
            Biomass yield (g_biomass / g_substrate).

        Returns
        -------
        callable
            ``rate_fn(env) -> float`` returning mol substrate / h.
        """
        # Capture all parameters in the closure
        kin = self
        mw_x = float(MW_organism)
        mw_s = float(MW_substrate)
        Y = float(yield_gX_gS)
        org_id = str(organism_id)
        sub_id = str(substrate_id)

        def rate_fn(env):
            C_S = env.concentrations.get(sub_id, 0.0)  # mol/L
            C_X = env.concentrations.get(org_id, 0.0)   # mol/L

            S_gL = C_S * mw_s   # g/L
            X_gL = C_X * mw_x   # g/L

            if X_gL <= 1e-30 or S_gL <= 0.0:
                return 0.0

            mu_val = kin.mu(S_gL, X_gL)

            if mu_val <= 0.0:
                return 0.0

            # dS/dt = -(μ/Y) × X_gL  (g_substrate/L/h)
            # in mol/h: × V_L / MW_substrate
            return (mu_val / Y) * X_gL / mw_s * env.V_L

        return rate_fn


# ════════════════════════════════════════════════════════════════════════
#  Monod
# ════════════════════════════════════════════════════════════════════════

@dataclass
class Monod(_KineticsBase):
    """Classic Monod kinetics.

    μ = μ_max × S / (Ks + S)

    Parameters
    ----------
    mu_max : float
        Maximum specific growth rate (1/h).
    Ks : float
        Half-saturation constant (g_substrate/L).
    """

    mu_max: float = 0.5
    Ks: float = 5e-3

    @property
    def label(self) -> str:
        return "Monod"

    def mu(self, S_gL: float, X_gL: float) -> float:
        denom = self.Ks + S_gL
        if denom <= 0.0:
            return 0.0
        return self.mu_max * S_gL / denom


# ════════════════════════════════════════════════════════════════════════
#  Contois
# ════════════════════════════════════════════════════════════════════════

@dataclass
class Contois(_KineticsBase):
    """Contois (density-dependent) kinetics.

    μ = μ_max × (S/X) / (Ks + S/X)

    where S and X are in g/L, and Ks is in g_substrate/g_biomass.

    Parameters
    ----------
    mu_max : float
        Maximum specific growth rate (1/h).
    Ks : float
        Half-saturation constant (g_substrate / g_biomass).
    """

    mu_max: float = 0.5
    Ks: float = 0.1

    @property
    def label(self) -> str:
        return "Contois"

    def mu(self, S_gL: float, X_gL: float) -> float:
        if X_gL <= 1e-30:
            return 0.0
        ratio = S_gL / X_gL
        denom = self.Ks + ratio
        if denom <= 0.0:
            return 0.0
        return self.mu_max * ratio / denom


# ════════════════════════════════════════════════════════════════════════
#  Andrews (Haldane) — substrate inhibition
# ════════════════════════════════════════════════════════════════════════

@dataclass
class Andrews(_KineticsBase):
    """Andrews (Haldane) kinetics — substrate inhibition.

    μ = μ_max × S / (Ks + S + S²/Ki)

    Growth increases with substrate at low S (Monod-like) but decreases
    at high S due to substrate inhibition.

    Parameters
    ----------
    mu_max : float
        Maximum specific growth rate (1/h).  Note: the actual peak
        growth rate is μ_max × √(Ki/Ks) / (1 + √(Ks/Ki))², which is
        less than μ_max when Ki is finite.
    Ks : float
        Half-saturation constant (g_substrate/L).
    Ki : float
        Substrate inhibition constant (g_substrate/L).  Higher Ki means
        less inhibition.
    """

    mu_max: float = 0.5
    Ks: float = 5e-3
    Ki: float = 50.0

    @property
    def label(self) -> str:
        return "Andrews"

    def mu(self, S_gL: float, X_gL: float) -> float:
        denom = self.Ks + S_gL + S_gL * S_gL / self.Ki
        if denom <= 0.0:
            return 0.0
        return self.mu_max * S_gL / denom


# ════════════════════════════════════════════════════════════════════════
#  ContoisAndrews — density-dependent + substrate inhibition
# ════════════════════════════════════════════════════════════════════════

@dataclass
class ContoisAndrews(_KineticsBase):
    """Contois kinetics with Andrews-type substrate inhibition.

    μ = μ_max × r / (Ks + r + r²/Ki)

    where r = S/X (g_substrate/g_biomass).

    Parameters
    ----------
    mu_max : float
        Maximum specific growth rate (1/h).
    Ks : float
        Half-saturation constant (g_substrate / g_biomass).
    Ki : float
        Inhibition constant (g_substrate / g_biomass).
    """

    mu_max: float = 0.5
    Ks: float = 0.1
    Ki: float = 50.0

    @property
    def label(self) -> str:
        return "ContoisAndrews"

    def mu(self, S_gL: float, X_gL: float) -> float:
        if X_gL <= 1e-30:
            return 0.0
        ratio = S_gL / X_gL
        denom = self.Ks + ratio + ratio * ratio / self.Ki
        if denom <= 0.0:
            return 0.0
        return self.mu_max * ratio / denom


# ════════════════════════════════════════════════════════════════════════
#  Tessier — exponential saturation
# ════════════════════════════════════════════════════════════════════════

@dataclass
class Tessier(_KineticsBase):
    """Tessier kinetics — exponential saturation.

    μ = μ_max × (1 − exp(−S/Ks))

    Approaches μ_max exponentially rather than hyperbolically.

    Parameters
    ----------
    mu_max : float
        Maximum specific growth rate (1/h).
    Ks : float
        Characteristic substrate concentration (g/L).
    """

    mu_max: float = 0.5
    Ks: float = 5e-3

    @property
    def label(self) -> str:
        return "Tessier"

    def mu(self, S_gL: float, X_gL: float) -> float:
        import math
        if self.Ks <= 0.0:
            return self.mu_max
        return self.mu_max * (1.0 - math.exp(-S_gL / self.Ks))


# ════════════════════════════════════════════════════════════════════════
#  Moser — sigmoidal (generalised Monod with Hill exponent)
# ════════════════════════════════════════════════════════════════════════

@dataclass
class Moser(_KineticsBase):
    """Moser kinetics — sigmoidal response.

    μ = μ_max × S^n / (Ks + S^n)

    Reduces to Monod when n=1.  For n>1, the response is sigmoidal
    (cooperative substrate binding).

    Parameters
    ----------
    mu_max : float
        Maximum specific growth rate (1/h).
    Ks : float
        Half-saturation constant (g^n / L^n when n ≠ 1).
    n : float
        Hill coefficient (dimensionless).  n=1 is Monod.
    """

    mu_max: float = 0.5
    Ks: float = 5e-3
    n: float = 1.0

    @property
    def label(self) -> str:
        return "Moser"

    def mu(self, S_gL: float, X_gL: float) -> float:
        Sn = S_gL ** self.n
        denom = self.Ks + Sn
        if denom <= 0.0:
            return 0.0
        return self.mu_max * Sn / denom


# ════════════════════════════════════════════════════════════════════════
#  Blackman — linear up to saturation
# ════════════════════════════════════════════════════════════════════════

@dataclass
class Blackman(_KineticsBase):
    """Blackman kinetics — piecewise linear.

    μ = μ_max × min(1, S / (2×Ks))

    Linear growth rate at low S, capped at μ_max.  The factor of 2
    ensures that μ = μ_max/2 at S = Ks (same half-saturation
    definition as Monod).

    Parameters
    ----------
    mu_max : float
        Maximum specific growth rate (1/h).
    Ks : float
        Half-saturation constant (g/L).
    """

    mu_max: float = 0.5
    Ks: float = 5e-3

    @property
    def label(self) -> str:
        return "Blackman"

    def mu(self, S_gL: float, X_gL: float) -> float:
        if self.Ks <= 0.0:
            return self.mu_max
        return self.mu_max * min(1.0, S_gL / (2.0 * self.Ks))


# ════════════════════════════════════════════════════════════════════════
#  DualSubstrate — multiplicative limitation by two substrates
# ════════════════════════════════════════════════════════════════════════

@dataclass
class DualSubstrateMonod(_KineticsBase):
    """Double Monod — growth limited by two substrates multiplicatively.

    μ = μ_max × S/(Ks+S) × DO/(Ko+DO)

    Commonly used for oxygen-limited growth where both substrate and
    dissolved O₂ must be available.

    Parameters
    ----------
    mu_max : float
        Maximum specific growth rate (1/h).
    Ks : float
        Substrate half-saturation constant (g_substrate/L).
    secondary_id : str
        Species ID of the second limiting substrate (e.g. ``"O2"``).
    Ko : float
        Second substrate half-saturation constant (mol/L for dissolved
        gases, g/L for other substrates).
    secondary_in_mol_L : bool
        If True, the secondary substrate concentration is read in mol/L
        directly (appropriate for dissolved gases like O₂).  If False,
        converted from mol/L to g/L using MW.
    secondary_MW : float
        Molecular weight of secondary substrate (only used if
        secondary_in_mol_L is False).
    """

    mu_max: float = 0.5
    Ks: float = 5e-3
    secondary_id: str = "O2"
    Ko: float = 1e-5           # mol/L for dissolved O₂
    secondary_in_mol_L: bool = True
    secondary_MW: float = 32.0

    @property
    def label(self) -> str:
        return f"DualSubstrateMonod(+{self.secondary_id})"

    def mu(self, S_gL: float, X_gL: float) -> float:
        # This is overridden — see make_rate_fn below
        denom = self.Ks + S_gL
        if denom <= 0.0:
            return 0.0
        return self.mu_max * S_gL / denom

    def make_rate_fn(
        self,
        organism_id: str,
        substrate_id: str,
        MW_organism: float,
        MW_substrate: float,
        yield_gX_gS: float,
    ) -> Callable:
        """Override to read the secondary substrate from the environment."""
        kin = self
        mw_x = float(MW_organism)
        mw_s = float(MW_substrate)
        Y = float(yield_gX_gS)
        org_id = str(organism_id)
        sub_id = str(substrate_id)
        sec_id = str(self.secondary_id)
        Ko = float(self.Ko)
        sec_mol = bool(self.secondary_in_mol_L)
        sec_mw = float(self.secondary_MW)

        def rate_fn(env):
            C_S = env.concentrations.get(sub_id, 0.0)
            C_X = env.concentrations.get(org_id, 0.0)
            C_sec = env.concentrations.get(sec_id, 0.0)  # mol/L

            S_gL = C_S * mw_s
            X_gL = C_X * mw_x

            if X_gL <= 1e-30 or S_gL <= 0.0:
                return 0.0

            # Primary substrate (Monod)
            mu_val = kin.mu_max * S_gL / (kin.Ks + S_gL)

            # Secondary limitation
            if sec_mol:
                sec_val = C_sec  # already mol/L
            else:
                sec_val = C_sec * sec_mw  # convert to g/L

            if Ko + sec_val > 0:
                mu_val *= sec_val / (Ko + sec_val)
            else:
                mu_val = 0.0

            if mu_val <= 0.0:
                return 0.0

            return (mu_val / Y) * X_gL / mw_s * env.V_L

        return rate_fn
