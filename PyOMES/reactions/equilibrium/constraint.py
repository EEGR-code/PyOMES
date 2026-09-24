# -*- coding: utf-8 -*-
"""The equilibrium-constraint protocol and the helpers that work on it.

:class:`EquilibriumConstraint` is the structural contract shared by every
mass-action equilibrium: a stoichiometry, a reference ``log_K``, an optional
Van 't Hoff enthalpy and the reference temperature. Both
:class:`~PyOMES.reactions.equilibrium.reaction.EquilibriumReaction` and the
named interphase constraints in :mod:`PyOMES.reactions.equilibrium.interphase`
satisfy it without sharing a base class.

- :func:`vant_hoff_log_K` returns ``log_K`` corrected to a temperature.
- :func:`classify_equilibrium_constraint` sorts a constraint into
  ``"acid_base"``, ``"gas_liquid"`` or ``"solid_liquid"`` from the phase tags
  of its stoichiometry.

:class:`~PyOMES.reactions.reaction_system.ReactionSystem` and the
chemical-equilibrium engines accept and sort equilibria through these names,
not by concrete type. One exception: the Bisection engine accepts only
``EquilibriumReaction`` for single-phase (acid-base) items, because only it
carries acid-base pKa semantics.
"""

from __future__ import annotations

import math
from typing import Literal, Optional, Protocol, Sequence, runtime_checkable

from PyOMES.units import R_J_PER_MOL_K as _R_J_MOL_K
from PyOMES.reactions.stoichiometry import StoichiometryEntry
from PyOMES.reactions._shared import phases_from_entries

_LOG10_E = 1.0 / math.log(10.0)


@runtime_checkable
class EquilibriumConstraint(Protocol):
    """Structural contract for the mass-action equilibrium family.

    Any type exposing these four attributes — a stoichiometry, a
    reference-temperature ``log_K``, an optional Van 't Hoff
    ``dH_J_per_mol``, and the reference temperature they were measured
    at — can be routed through :func:`vant_hoff_log_K` and the
    equilibrium-classification/tableau-building machinery in
    ``PyOMES.chemical_equilibrium``, regardless of whether it also satisfies
    :class:`~PyOMES.chemistry.partition.PartitionModel` (as
    :class:`~PyOMES.reactions.phase_equilibria.HenryEquilibrium`,
    :class:`~PyOMES.reactions.phase_equilibria.KspEquilibrium`, and
    :class:`~PyOMES.reactions.phase_equilibria.RaoultEquilibrium` all do).

    ``log_K``/``dH_J_per_mol`` are plain attributes carrying the
    *reference* mass-action constant (at ``T_ref_K``) — not methods.
    Composition-dependent corrections (``γ_i``, ``φ_i``) are not part
    of this protocol; temperature correction is handled separately by
    :func:`vant_hoff_log_K`.
    """

    stoichiometry: Sequence[StoichiometryEntry]
    log_K: float
    dH_J_per_mol: Optional[float]
    T_ref_K: float


def vant_hoff_log_K(constraint: EquilibriumConstraint, T_K: float) -> float:
    """Van 't Hoff temperature-corrected log10(K) for any EquilibriumConstraint.

    Returns ``constraint.log_K`` unchanged when ``dH_J_per_mol`` is
    ``None``/~0 or when ``T_K`` is at the reference temperature.
    """
    log_K_ref = float(constraint.log_K)
    dH = constraint.dH_J_per_mol
    T_ref_K = float(constraint.T_ref_K)
    if dH is None or abs(dH) < 1e-30:
        return log_K_ref
    if abs(T_K - T_ref_K) < 1e-10:
        return log_K_ref
    # ln K(T) = ln K(T_ref) − (ΔH/R) (1/T − 1/T_ref)
    delta_ln_K = -(float(dH) / _R_J_MOL_K) * (1.0 / float(T_K) - 1.0 / T_ref_K)
    return log_K_ref + delta_ln_K * _LOG10_E


def classify_equilibrium_constraint(
    item: EquilibriumConstraint,
) -> Literal["acid_base", "gas_liquid", "solid_liquid"]:
    """Classify an :class:`EquilibriumConstraint` from its stoichiometry's phase tags.

    Not an ``isinstance(item, EquilibriumReaction)`` check —
    :class:`~PyOMES.reactions.phase_equilibria.HenryEquilibrium`,
    :class:`~PyOMES.reactions.phase_equilibria.KspEquilibrium`, and
    :class:`~PyOMES.reactions.phase_equilibria.RaoultEquilibrium` are siblings,
    not subclasses, of :class:`~PyOMES.reactions.equilibrium.reaction.EquilibriumReaction`. Any
    ``EquilibriumConstraint``-conforming item is classified purely from
    the distinct phases present across its ``stoichiometry``:

    - any entry with ``phase == "solid"`` → ``"solid_liquid"``
      (takes priority — a solid entry always means a precipitation
      equilibrium, regardless of what else is present).
    - more than one distinct phase (and no solid) → ``"gas_liquid"``.
    - a single phase → ``"acid_base"``.

    Raises
    ------
    ValueError
        If ``item.stoichiometry`` is empty (e.g. a
        ``PartitionModel``-only ``HenryEquilibrium``/``RaoultEquilibrium``
        constructed without ``gas_species``/``liquid_species`` — such an
        instance has no reaction row to classify).
    """
    entries = item.stoichiometry
    if not entries:
        raise ValueError(
            f"classify_equilibrium_constraint: {item!r} has an empty "
            "stoichiometry — nothing to classify (e.g. a PartitionModel-"
            "only HenryEquilibrium/RaoultEquilibrium with gas_species/"
            "liquid_species unset)."
        )
    phases = phases_from_entries(entries)
    if "solid" in phases:
        return "solid_liquid"
    if len(phases) > 1:
        return "gas_liquid"
    return "acid_base"
