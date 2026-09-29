# -*- coding: utf-8 -*-
"""Correct a quantity known at a reference temperature to another temperature.

Three relations share one form,

    X(T) = X(T_ref) · exp(−(E/R) · (1/T − 1/T_ref))

and differ only in what X and E are:

- **van 't Hoff** — X an equilibrium constant (K, Ka, Kw, Ksp, a Henry
  constant), E the standard reaction enthalpy ΔH°;
- **Clausius–Clapeyron** — X a saturation vapour pressure, E the enthalpy of
  vaporisation;
- **Arrhenius** — X a rate constant, E the activation energy
  (:func:`PyOMES.reactions.kinetic.rate_laws.arrhenius_factor`, which calls
  :func:`ln_correction`).

:func:`ln_correction` is the one place the relation is written. It takes E/R in
kelvin, the quantity every caller either has (Sander's ``dlnH``, ``ΔH_vap/R``,
``Ea/R``) or computes first (``ΔH / R`` for an enthalpy in J/mol). The wrappers
below apply it to a particular quantity, and each keeps its own documented edge
rules.
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np

from PyOMES.units import R_J_PER_MOL_K as _R_J_MOL_K

_LOG10_E = np.log10(np.e)        # 1/ln(10), used to convert ln K to log10 K


def ln_correction(E_over_R_K, T_K, T_ref_K):
    """Return ``ln X(T) − ln X(T_ref)`` for an energy term ``E_over_R_K = E/R`` (K).

    ``−(E/R) · (1/T − 1/T_ref)``. Works on floats and on numpy arrays of ``T_K``.
    """
    return -E_over_R_K * (1.0 / T_K - 1.0 / T_ref_K)


def vant_hoff_K(
    K_ref: float,
    dH_J_per_mol: float,
    T_K: float,
    T_ref_K: float = 298.15,
) -> float:
    """Temperature-correct an equilibrium constant using the van 't Hoff relation.

    Notes:
      - ΔH° should be the *standard enthalpy change of the equilibrium reaction*.
      - If dH_J_per_mol is 0, this reduces to K(T) = K_ref (no temperature effect).
      - ``K_ref`` is also returned unchanged when it is non-finite or <= 0, when
        ``dH_J_per_mol`` is non-finite, or when ``T_K`` is non-finite or <= 0.
    """
    K_ref = float(K_ref)
    dH_J_per_mol = float(dH_J_per_mol)
    T_K = float(T_K)
    T_ref_K = float(T_ref_K)
    if not np.isfinite(K_ref) or K_ref <= 0.0:
        return float(K_ref)
    if not np.isfinite(dH_J_per_mol) or abs(dH_J_per_mol) < 1e-30:
        return float(K_ref)
    if not np.isfinite(T_K) or T_K <= 0.0:
        return float(K_ref)
    return float(K_ref * np.exp(ln_correction(dH_J_per_mol / _R_J_MOL_K, T_K, T_ref_K)))


def vant_hoff_log_K(
    log_K_ref: float,
    dH_J_per_mol: Optional[float],
    T_K: float,
    T_ref_K: float,
) -> float:
    """Apply the van 't Hoff correction to log10(K).

    Returns ``log_K_ref`` unchanged when ``dH_J_per_mol`` is None or ~0, or when
    ``T_K`` is at ``T_ref_K``.
    """
    if dH_J_per_mol is None or abs(dH_J_per_mol) < 1e-30:
        return float(log_K_ref)
    if abs(T_K - T_ref_K) < 1e-10:
        return float(log_K_ref)
    ln_shift = ln_correction(float(dH_J_per_mol) / _R_J_MOL_K, float(T_K), float(T_ref_K))
    return float(log_K_ref) + ln_shift * _LOG10_E


def henry_constant(kH_ref: float, dlnH_K: float, T_K: float, T_ref_K: float) -> float:
    """Temperature-correct a Henry solubility constant (van 't Hoff, Sander's form).

    ``kH(T) = kH_ref · exp(dlnH · (1/T − 1/T_ref))``, where ``dlnH_K`` is Sander's
    ``d ln kH / d(1/T)`` in kelvin (``= −ΔH_sol/R``). ``kH_ref`` can be in any units;
    the result has the same units.
    """
    return kH_ref * math.exp(ln_correction(-dlnH_K, T_K, T_ref_K))


def clausius_clapeyron(
    P_ref: float, dH_vap_over_R_K: float, T_K: float, T_ref_K: float,
) -> float:
    """Saturation vapour pressure at ``T_K`` from its value ``P_ref`` at ``T_ref_K``.

    ``P(T) = P_ref · exp(−(ΔH_vap/R) · (1/T − 1/T_ref))``, with ``ΔH_vap/R`` in
    kelvin. ``P_ref`` can be in any pressure unit; the result has the same unit.
    """
    return P_ref * math.exp(ln_correction(dH_vap_over_R_K, T_K, T_ref_K))
