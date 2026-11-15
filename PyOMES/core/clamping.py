# -*- coding: utf-8 -*-
"""Swappable non-negativity clamping for discrete-step StepSolvers.

A ``clamp_fn`` restores a physical invariant (non-negative inventory)
that the raw discrete-step numerics don't guarantee on their own — a
large ``dt_h`` combined with an otherwise well-posed rate law can
overdraw a species below zero over one step. This is a distinct
concern from *composition* (how sub-systems combine) and *integration
method* (how the combined system is advanced); see
docs/dev/implementation/shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md item 5 for
the three-dimension terminology this module is one leg of.

Every function here has the signature
``clamp_fn(deltas, current_mol, dt_h) -> deltas`` and must operate
correctly on an **arbitrary subset** of species — never assume
whole-dict ownership. This is what lets a caller compose a bespoke
per-species convention (e.g. matching an external reference model's
clamping rule) out of these functions as ingredients, rather than
being stuck with an all-or-nothing choice::

    def bsm2_style_clamp(deltas, current_mol, dt_h):
        result = dict(deltas)
        for sp in NEVER_CLAMP:
            result[sp] = deltas[sp]
        remaining = {k: v for k, v in deltas.items() if k not in NEVER_CLAMP}
        result.update(proportional_clamp(remaining, current_mol, dt_h))
        return result

    SimultaneousEulerSolver(clamp_fn=bsm2_style_clamp)

No solver clamps by default: a ``clamp_fn`` is applied only when it is
passed to a step solver. The functions here scale removals species by
species, so they do not conserve mass.
"""

from __future__ import annotations

from typing import Dict


def proportional_clamp(
    deltas: Dict[str, float],
    current_mol: Dict[str, float],
    dt_h: float,
) -> Dict[str, float]:
    """Scale each overdrawn species' removal so its amount lands at zero.

    For each species where the net delta would overdraw the available
    inventory, that species' removal rate is scaled down so that the
    result is exactly zero (never negative). Each species is scaled on its
    own: the other species of the reaction or flux that overdrew it are
    applied in full, so the step's proportions are broken and material is
    created (for a reaction that consumes the scaled species, its products
    are made as if it had been consumed in full). The result is the same as
    :func:`floor_clamp` with ``eps=0``.

    Parameters
    ----------
    deltas : dict
        ``{species: flux_mol_per_h}``.
    current_mol : dict
        ``{species: current_moles}``.
    dt_h : float
        Timestep duration (hours).

    Returns
    -------
    dict
        Clamped fluxes (same format and species subset as input).
    """
    clamped = dict(deltas)
    for species, rate in deltas.items():
        n_current = current_mol.get(species, 0.0)
        delta_mol = rate * dt_h
        n_after = n_current + delta_mol
        if n_after < 0.0 and delta_mol < 0.0:
            max_removal = n_current
            if abs(delta_mol) > 1e-30:
                scale = max_removal / abs(delta_mol)
                scale = max(0.0, min(1.0, scale))
                clamped[species] = rate * scale
    return clamped


def floor_clamp(
    deltas: Dict[str, float],
    current_mol: Dict[str, float],
    dt_h: float,
    *,
    eps: float = 0.0,
) -> Dict[str, float]:
    """Per-species floor: scale each removal rate that would cross
    ``eps`` independently, ignoring the effect on any other species.

    Simpler and cheaper than :func:`proportional_clamp`, but does not
    preserve relative stoichiometry between simultaneously-changing
    species — matches the blunt per-species floor
    :meth:`~PyOMES.core.phases.Phase.apply_flux` has always applied by
    default. ``eps`` supports reference models that floor at a small
    epsilon rather than exact zero.
    """
    clamped = dict(deltas)
    for species, rate in deltas.items():
        n_current = current_mol.get(species, 0.0)
        delta_mol = rate * dt_h
        n_after = n_current + delta_mol
        if n_after < eps and delta_mol < 0.0:
            max_removal = n_current - eps
            if abs(delta_mol) > 1e-30:
                scale = max_removal / abs(delta_mol)
                scale = max(0.0, min(1.0, scale))
                clamped[species] = rate * scale
    return clamped
