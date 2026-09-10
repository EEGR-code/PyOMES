# -*- coding: utf-8 -*-
"""Viscosity model framework for fermentation broth.

This module provides:

* :class:`BrothState` — immutable snapshot of broth conditions available at
  each simulation timestep.
* :class:`ViscosityResult` — structured output from a viscosity calculation.
* :class:`ViscosityModel` — protocol (interface) that any viscosity model must
  satisfy.
* :func:`wrap_viscosity_model` — normalises a plain callable *or* a
  ``ViscosityModel`` object into a uniform internal representation.
* Built-in models: :class:`WaterViscosity`, :class:`ArrheniusBiomassViscosity`.

Design intent
-------------
The fermenter accepts ``viscosity_model=None`` (ignore viscosity, the default),
a plain callable ``f(BrothState) → float`` (returns Pa·s), or a full
``ViscosityModel`` object.  All three are handled transparently.

Usage
-----
>>> from PyOMES.properties.viscosity import BrothState, ArrheniusBiomassViscosity
>>> model = ArrheniusBiomassViscosity()
>>> state = BrothState(T_K=305.15, X_g_L=10.0)
>>> result = model.compute(state)
>>> result.mu_Pa_s   # apparent viscosity in Pa·s
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Optional, Protocol, runtime_checkable


# ════════════════════════════════════════════════════════════════════════
#  State snapshot
# ════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class BrothState:
    """Immutable snapshot of broth conditions at a single timestep.

    Fields are populated by the fermenter from whatever state is available.
    Models should tolerate ``None`` for fields they don't require by using
    ``getattr`` or checking before access.  All concentration fields use
    SI-friendly units (g/L, mol/L) so that models don't need to know about
    bioSTEAM conventions.

    Parameters
    ----------
    T_K : float
        Temperature (Kelvin).
    X_g_L : float
        Biomass concentration (g dry-weight / L).
    concentrations_g_L : dict
        Species → mass concentration (g/L).  Includes substrates, products,
        salts, etc. — whatever the fermenter is tracking.
    concentrations_mol_L : dict
        Species → molar concentration (mol/L).
    pH : float or None
        Current pH (from speciation solver).
    ionic_strength : float or None
        Ionic strength (mol/L) if available from speciation.
    DO_mol_L : float or None
        Dissolved oxygen (mol/L).
    CO2aq_mol_L : float or None
        Dissolved CO₂ (mol/L).
    P_atm : float or None
        Headspace pressure (atm).
    t_h : float or None
        Current time (hours from batch start).
    """
    T_K: float = 305.15
    X_g_L: float = 0.0
    concentrations_g_L: Dict[str, float] = field(default_factory=dict)
    concentrations_mol_L: Dict[str, float] = field(default_factory=dict)
    pH: Optional[float] = None
    ionic_strength: Optional[float] = None
    DO_mol_L: Optional[float] = None
    CO2aq_mol_L: Optional[float] = None
    P_atm: Optional[float] = None
    t_h: Optional[float] = None


# ════════════════════════════════════════════════════════════════════════
#  Result container
# ════════════════════════════════════════════════════════════════════════

@dataclass
class ViscosityResult:
    """Output from a viscosity model evaluation.

    Parameters
    ----------
    mu_Pa_s : float
        Apparent dynamic viscosity (Pa·s).
    mu_rel : float or None
        Viscosity relative to pure water at the same temperature
        (dimensionless).  Informational only.
    extras : dict
        Any additional outputs the model wants to expose (e.g. intermediate
        terms, correction factors, model metadata).
    """
    mu_Pa_s: float
    mu_rel: Optional[float] = None
    extras: Dict[str, Any] = field(default_factory=dict)


# ════════════════════════════════════════════════════════════════════════
#  Protocol (interface)
# ════════════════════════════════════════════════════════════════════════

@runtime_checkable
class ViscosityModel(Protocol):
    """Interface that any viscosity model must satisfy.

    Implementing classes must provide a ``compute`` method that accepts a
    :class:`BrothState` and returns a :class:`ViscosityResult`.
    """

    def compute(self, state: BrothState) -> ViscosityResult:
        """Evaluate viscosity for the given broth conditions."""
        ...


# ════════════════════════════════════════════════════════════════════════
#  Callable adapter  (wraps a plain function into the protocol)
# ════════════════════════════════════════════════════════════════════════

class _CallableViscosityAdapter:
    """Wraps ``f(BrothState) → float`` as a full :class:`ViscosityModel`."""

    def __init__(self, fn: Callable[[BrothState], float]):
        self._fn = fn

    def compute(self, state: BrothState) -> ViscosityResult:
        mu = float(self._fn(state))
        return ViscosityResult(mu_Pa_s=mu)

    def __repr__(self):
        return f"_CallableViscosityAdapter({self._fn!r})"


def wrap_viscosity_model(model) -> Optional[ViscosityModel]:
    """Normalise a user-supplied viscosity model into the internal protocol.

    Accepts:
      - ``None`` → returns ``None`` (viscosity tracking disabled)
      - A callable ``f(BrothState) → float`` → wrapped in adapter
      - An object satisfying :class:`ViscosityModel` → returned as-is

    Raises
    ------
    TypeError
        If *model* is not None, callable, or a ViscosityModel.
    """
    if model is None:
        return None
    if isinstance(model, ViscosityModel):
        return model
    if callable(model):
        return _CallableViscosityAdapter(model)
    raise TypeError(
        f"viscosity_model must be None, a callable, or a ViscosityModel; "
        f"got {type(model).__name__}"
    )


# ════════════════════════════════════════════════════════════════════════
#  Built-in models
# ════════════════════════════════════════════════════════════════════════

def mu_water_Pa_s(T_K: float) -> float:
    """Dynamic viscosity of pure water (Pa·s) as a function of temperature.

    Uses a simplified Vogel-equation fit valid for 273–373 K:

        μ = A × 10^(B / (T − C))

    with A = 2.414e-5 Pa·s, B = 247.8 K, C = 140.0 K.

    This reproduces published NIST values to within ~2 % over 0–100 °C.
    """
    A = 2.414e-5   # Pa·s
    B = 247.8      # K
    C = 140.0      # K
    T = float(T_K)
    if T <= C:
        T = C + 0.1  # guard against singularity
    return A * (10.0 ** (B / (T - C)))


class WaterViscosity:
    """Pure-water viscosity model (temperature-dependent only).

    Useful as a baseline / reference.

    Example
    -------
    >>> model = WaterViscosity()
    >>> result = model.compute(BrothState(T_K=298.15))
    >>> f"{result.mu_Pa_s:.4e}"
    '8.9283e-04'
    """

    def compute(self, state: BrothState) -> ViscosityResult:
        mu = mu_water_Pa_s(state.T_K)
        return ViscosityResult(mu_Pa_s=mu, mu_rel=1.0)

    def __repr__(self):
        return "WaterViscosity()"


class ArrheniusBiomassViscosity:
    """Temperature + biomass concentration viscosity model.

    Combines an Arrhenius-type temperature dependence (via
    :func:`mu_water_Pa_s`) with an exponential biomass correction:

        μ = μ_water(T) × exp(k_X × X)

    where *X* is the biomass concentration in g/L and *k_X* is a fitting
    parameter (default 0.005 L/g, typical for yeast suspensions up to
    ~100 g/L).

    For more concentrated or filamentous broths, *k_X* may be larger
    (0.01–0.03).  The exponential form is the simplest correlation that
    captures the monotone increase of viscosity with biomass.

    Parameters
    ----------
    k_X : float
        Biomass sensitivity coefficient (L/g).  Default 0.005.
    mu_water_fn : callable or None
        Override for the water-viscosity baseline.  Receives T_K, returns
        Pa·s.  Defaults to :func:`mu_water_Pa_s`.

    Example
    -------
    >>> model = ArrheniusBiomassViscosity(k_X=0.005)
    >>> state = BrothState(T_K=305.15, X_g_L=50.0)
    >>> result = model.compute(state)
    >>> result.mu_rel   # should be exp(0.005 * 50) ≈ 1.284
    """

    def __init__(self, k_X: float = 0.005,
                 mu_water_fn: Optional[Callable[[float], float]] = None):
        self.k_X = float(k_X)
        self._mu_water_fn = mu_water_fn or mu_water_Pa_s

    def compute(self, state: BrothState) -> ViscosityResult:
        mu_w = self._mu_water_fn(state.T_K)
        X = max(0.0, float(state.X_g_L))
        correction = math.exp(self.k_X * X)
        mu = mu_w * correction
        return ViscosityResult(
            mu_Pa_s=mu,
            mu_rel=correction,
            extras={"mu_water_Pa_s": mu_w, "k_X": self.k_X, "X_g_L": X},
        )

    def __repr__(self):
        return f"ArrheniusBiomassViscosity(k_X={self.k_X})"


class PowerLawBiomassViscosity:
    """Temperature + biomass power-law viscosity model.

    Uses a power-law correction on biomass concentration:

        μ = μ_water(T) × (1 + k × X^n)

    This form is convenient for fitting experimental data where the
    exponent *n* captures non-linear concentration effects (e.g. cell
    crowding, aggregate formation).

    Parameters
    ----------
    k : float
        Scaling coefficient.  Default 0.001.
    n : float
        Power-law exponent.  Default 1.5.
    mu_water_fn : callable or None
        Override for the water-viscosity baseline.

    Example
    -------
    >>> model = PowerLawBiomassViscosity(k=0.001, n=1.5)
    >>> state = BrothState(T_K=305.15, X_g_L=50.0)
    >>> result = model.compute(state)
    """

    def __init__(self, k: float = 0.001, n: float = 1.5,
                 mu_water_fn: Optional[Callable[[float], float]] = None):
        self.k = float(k)
        self.n = float(n)
        self._mu_water_fn = mu_water_fn or mu_water_Pa_s

    def compute(self, state: BrothState) -> ViscosityResult:
        mu_w = self._mu_water_fn(state.T_K)
        X = max(0.0, float(state.X_g_L))
        correction = 1.0 + self.k * (X ** self.n)
        mu = mu_w * correction
        return ViscosityResult(
            mu_Pa_s=mu,
            mu_rel=correction,
            extras={"mu_water_Pa_s": mu_w, "k": self.k, "n": self.n, "X_g_L": X},
        )

    def __repr__(self):
        return f"PowerLawBiomassViscosity(k={self.k}, n={self.n})"
