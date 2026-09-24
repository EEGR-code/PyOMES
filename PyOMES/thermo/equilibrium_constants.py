# -*- coding: utf-8 -*-
"""Van 't Hoff temperature correction of equilibrium constants.

    ln K(T) = ln K(T_ref) - (ΔH°/R) (1/T - 1/T_ref)

The correction itself is written once, in :func:`vant_hoff_delta_ln_K`. Two
thin wrappers apply it to K (:func:`vant_hoff_K`) and to log10(K)
(:func:`vant_hoff_log_K`). They keep their own edge-case rules on purpose:
the K-space wrapper also returns ``K_ref`` unchanged for a non-positive or
non-finite ``K_ref`` or ``T_K``, whereas the log-space wrapper accepts
``dH_J_per_mol=None`` and skips the correction at ``T_K == T_ref_K``.

Note: ``PyOMES.reactions.equilibrium.constraint.vant_hoff_log_K(constraint, T_K)`` is a
separate function of the same stem that takes an ``EquilibriumConstraint``;
it uses a log10(e) constant that differs in the last bit, so it is not
folded into this module (see ``docs/dev/implementation/OPEN_WORK.md``).
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from ..units import R_J_PER_MOL_K as _R_J_MOL_K

_LOG10_E = np.log10(np.e)        # 1/ln(10), used to convert ln K to log10 K


def vant_hoff_delta_ln_K(dH_J_per_mol: float, T_K: float, T_ref_K: float) -> float:
    """Return ``ln K(T) - ln K(T_ref)`` for a reaction enthalpy ``dH_J_per_mol``."""
    return -(float(dH_J_per_mol) / _R_J_MOL_K) * (1.0 / float(T_K) - 1.0 / float(T_ref_K))


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
    return float(K_ref * np.exp(vant_hoff_delta_ln_K(dH_J_per_mol, T_K, T_ref_K)))


def vant_hoff_log_K(
    log_K_ref: float,
    dH_J_per_mol: Optional[float],
    T_K: float,
    T_ref_K: float,
) -> float:
    """Apply Van't Hoff correction to log10(K).

    Returns ``log_K_ref`` unchanged when ``dH_J_per_mol`` is None or ~0, or when
    ``T_K`` is at ``T_ref_K``.
    """
    if dH_J_per_mol is None or abs(dH_J_per_mol) < 1e-30:
        return float(log_K_ref)
    if abs(T_K - T_ref_K) < 1e-10:
        return float(log_K_ref)
    return float(log_K_ref) + vant_hoff_delta_ln_K(dH_J_per_mol, T_K, T_ref_K) * _LOG10_E
