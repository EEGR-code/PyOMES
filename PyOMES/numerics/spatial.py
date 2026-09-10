# -*- coding: utf-8 -*-
"""Spatial discretisation schemes for 1D finite-volume transport models.

Provides pluggable advection and dispersion schemes for any model that
discretises a 1D domain into N finite-volume cells.  Each scheme is a
factory function that takes grid parameters at setup time and returns
a vectorised numpy callable for use inside an ODE right-hand side.

The schemes are registered by name and selected via string keywords,
so model classes (e.g. ``HPLCColumn``) can expose the scheme choice
as a simple constructor parameter without coupling to any specific
numerical implementation.

Advection schemes
-----------------
- ``"upwind"`` — first-order upwind differencing.  Simple, stable,
  but introduces numerical diffusion D_num = u × Δz / 2.
- ``"tvd_vanleer"`` — second-order TVD with van Leer flux limiter.
  Negligible numerical diffusion; oscillation-free near sharp fronts.

Dispersion schemes
------------------
- ``"central_2nd"`` — standard 3-point central difference, O(Δz²).
- ``"central_4th"`` — compact 5-point central difference, O(Δz⁴).

Usage
-----
From a model's ODE builder::

    from PyOMES.numerics.spatial import get_advection_scheme, get_dispersion_scheme

    advect = get_advection_scheme("tvd_vanleer", u=3.14, dz=0.15, N=200, n_sp=3)
    disperse = get_dispersion_scheme("central_4th", D_ax=D_ax_array, dz=0.15, N=200, n_sp=3)

    def f(t, y):
        ...
        adv = advect(C, C_in)
        disp = disperse(C, C_in)
        dCdt = adv + disp - phase_ratio * dqdt
        ...

Registering a custom scheme::

    from PyOMES.numerics.spatial import register_advection_scheme

    @register_advection_scheme("my_scheme")
    def build_my_scheme(u, dz, N, n_sp):
        def advect(C, C_in):
            ...
            return dCdt_adv
        return advect

See Also
--------
:class:`~fermenter.models.hplc_column.HPLCColumn` — primary consumer.
:class:`~fermenter.core.links.AdvectiveLink` — link-based transport
    (for arbitrary topologies; schemes assume a regular 1D grid).
"""

from __future__ import annotations

import numpy as np


# ════════════════════════════════════════════════════════════════════════
#  Registries
# ════════════════════════════════════════════════════════════════════════

_ADVECTION_SCHEMES = {}   # name → factory(u, dz, N, n_sp) → callable
_DISPERSION_SCHEMES = {}  # name → factory(D_ax, dz, N, n_sp) → callable


# ════════════════════════════════════════════════════════════════════════
#  Registration decorators
# ════════════════════════════════════════════════════════════════════════

def register_advection_scheme(name):
    """Decorator to register an advection scheme factory.

    The decorated function must have signature::

        def factory(u, dz, N, n_sp) -> callable

    and return a callable with signature::

        def advect(C, C_in) -> dCdt_adv

    where ``C`` is (N, n_sp), ``C_in`` is (1, n_sp), and the
    return value is (N, n_sp).
    """
    def decorator(fn):
        _ADVECTION_SCHEMES[name] = fn
        return fn
    return decorator


def register_dispersion_scheme(name):
    """Decorator to register a dispersion scheme factory.

    The decorated function must have signature::

        def factory(D_ax, dz, N, n_sp) -> callable

    where ``D_ax`` is an ndarray of shape (n_sp,) giving per-species
    axial dispersion coefficients.

    The returned callable must have signature::

        def disperse(C, C_in) -> dCdt_disp

    where ``C`` is (N, n_sp), ``C_in`` is (1, n_sp), and the
    return value is (N, n_sp).
    """
    def decorator(fn):
        _DISPERSION_SCHEMES[name] = fn
        return fn
    return decorator


# ════════════════════════════════════════════════════════════════════════
#  Public query / access API
# ════════════════════════════════════════════════════════════════════════

def available_advection_schemes():
    """Return sorted list of registered advection scheme names."""
    return sorted(_ADVECTION_SCHEMES.keys())


def available_dispersion_schemes():
    """Return sorted list of registered dispersion scheme names."""
    return sorted(_DISPERSION_SCHEMES.keys())


def get_advection_scheme(name, u, dz, N, n_sp):
    """Look up an advection scheme by name and build it.

    Parameters
    ----------
    name : str
        Registered scheme name (e.g. ``"tvd_vanleer"``).
    u : float
        Interstitial velocity (length/time).
    dz : float
        Cell length (same length unit as u).
    N : int
        Number of cells.
    n_sp : int
        Number of species.

    Returns
    -------
    callable
        ``advect(C, C_in) -> dCdt_adv`` where C is (N, n_sp).

    Raises
    ------
    ValueError
        If the scheme name is not registered.
    """
    name = name.lower().strip()
    if name not in _ADVECTION_SCHEMES:
        raise ValueError(
            f"Unknown advection scheme {name!r}. "
            f"Available: {available_advection_schemes()}"
        )
    return _ADVECTION_SCHEMES[name](u, dz, N, n_sp)


def get_dispersion_scheme(name, D_ax, dz, N, n_sp):
    """Look up a dispersion scheme by name and build it.

    Parameters
    ----------
    name : str
        Registered scheme name (e.g. ``"central_4th"``).
    D_ax : ndarray, shape (n_sp,)
        Per-species axial dispersion coefficients (length²/time).
    dz : float
        Cell length (same length unit as D_ax).
    N : int
        Number of cells.
    n_sp : int
        Number of species.

    Returns
    -------
    callable
        ``disperse(C, C_in) -> dCdt_disp`` where C is (N, n_sp).

    Raises
    ------
    ValueError
        If the scheme name is not registered.
    """
    name = name.lower().strip()
    if name not in _DISPERSION_SCHEMES:
        raise ValueError(
            f"Unknown dispersion scheme {name!r}. "
            f"Available: {available_dispersion_schemes()}"
        )
    return _DISPERSION_SCHEMES[name](D_ax, dz, N, n_sp)


# ════════════════════════════════════════════════════════════════════════
#  Built-in advection schemes
# ════════════════════════════════════════════════════════════════════════

@register_advection_scheme("upwind")
def build_upwind(u, dz, N, n_sp):
    """First-order upwind advection.

    Numerical diffusion: D_num = u × Δz / 2.
    Effective plates ≈ N_cells when D_num >> D_ax.

    Simple and unconditionally stable, but introduces substantial
    artificial peak broadening.
    """
    adv_coeff = u / dz

    def advect(C, C_in):
        C_left = np.vstack([C_in, C[:-1, :]])
        return adv_coeff * (C_left - C)

    return advect


@register_advection_scheme("tvd_vanleer")
def build_tvd_vanleer(u, dz, N, n_sp):
    """Second-order TVD advection with van Leer flux limiter.

    Numerical diffusion: O(Δz²), negligible for practical cell counts.
    Effective plates determined by physical D_ax, not grid spacing.

    The scheme reconstructs a piecewise-linear concentration profile
    within each cell, computes face fluxes from the reconstructed
    values, and applies the van Leer limiter to prevent spurious
    oscillations near sharp fronts.

    Van Leer limiter: φ(r) = (r + |r|) / (1 + |r|)
      - r < 0 (local extremum): φ = 0 → falls back to 1st-order upwind
      - r = 1 (smooth):          φ = 1 → full 2nd-order correction
      - r → ∞ (step):            φ → 2 → maximum correction (limited)
    """
    u_over_dz = u / dz
    _eps = 1e-30

    def advect(C, C_in):
        C_pad = np.vstack([C_in, C_in, C, C[-1:, :]])
        delta = np.diff(C_pad, axis=0)
        delta_up = delta[:N + 1, :]
        delta_down = delta[1:N + 2, :]

        safe = np.abs(delta_up) > _eps
        r = np.where(safe, delta_down / np.where(safe, delta_up, 1.0), 0.0)
        phi = (r + np.abs(r)) / (1.0 + np.abs(r))

        C_upwind = C_pad[1:N + 2, :]
        C_face = C_upwind + 0.5 * phi * delta_up

        flux_in = C_face[:N, :]
        flux_out = C_face[1:N + 1, :]
        return u_over_dz * (flux_in - flux_out)

    return advect


# ════════════════════════════════════════════════════════════════════════
#  Built-in dispersion schemes
# ════════════════════════════════════════════════════════════════════════

@register_dispersion_scheme("central_2nd")
def build_dispersion_2nd(D_ax, dz, N, n_sp):
    """Standard 2nd-order central difference for axial dispersion.

    Stencil: d²C/dz² ≈ (C_{i-1} - 2C_i + C_{i+1}) / Δz²

    Truncation error: O(Δz²).
    Boundary conditions:
      - Inlet (left):  C_ghost = C_in  (Danckwerts)
      - Outlet (right): C_ghost = C_N  (zero-gradient)
    """
    coeff = D_ax / (dz * dz)

    def disperse(C, C_in):
        C_ext = np.vstack([C_in, C, C[-1:, :]])
        return coeff[None, :] * (
            C_ext[:-2, :] - 2 * C_ext[1:-1, :] + C_ext[2:, :]
        )

    return disperse


@register_dispersion_scheme("central_4th")
def build_dispersion_4th(D_ax, dz, N, n_sp):
    """4th-order central difference for axial dispersion.

    Stencil:
        d²C/dz² ≈ (-C_{i-2} + 16C_{i-1} - 30C_i + 16C_{i+1} - C_{i+2}) / (12 Δz²)

    Truncation error: O(Δz⁴).
    Boundary conditions (2 ghost cells per side):
      - Inlet:  C_{-2} = C_{-1} = C_in  (constant extrapolation)
      - Outlet: C_{N} = C_{N+1} = C_{N-1}  (zero-gradient)
    """
    coeff = D_ax / (12.0 * dz * dz)

    def disperse(C, C_in):
        C_ext = np.vstack([
            C_in, C_in,
            C,
            C[-1:, :], C[-1:, :]
        ])
        return coeff[None, :] * (
            -C_ext[:-4, :]
            + 16.0 * C_ext[1:-3, :]
            - 30.0 * C_ext[2:-2, :]
            + 16.0 * C_ext[3:-1, :]
            - C_ext[4:, :]
        )

    return disperse
