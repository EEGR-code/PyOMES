# -*- coding: utf-8 -*-
"""Shared validation helpers for the independent reaction classes.

Private module — not part of the public ``PyOMES.reactions`` surface.
Carries the free-function utilities common to
:class:`~PyOMES.reactions.kinetic.KineticReaction` and
:class:`~PyOMES.reactions.equilibrium.EquilibriumReaction` so the two
classes can remain independent (no ABC, no shared base) without
duplicating logic.

The helpers cover:

- :func:`coerce_and_validate` — list-coerce a stoichiometry sequence
  and run elemental balance validation at construction time.
- :func:`species_ids_from_entries` / :func:`phases_from_entries` —
  derive the sorted unique ID / phase keys from a stoichiometry.
- :func:`is_cross_phase_from_entries` — derived flag for whether the
  stoichiometry spans more than one phase. Today's behaviour
  preserved verbatim from the deleted ``Reaction.is_cross_phase``
  property.
"""

from __future__ import annotations

from typing import List, Optional, Sequence

from .stoichiometry import StoichiometryEntry, validate_balance


def _infer_elements(entries: Sequence[StoichiometryEntry]) -> List[str]:
    """Return the sorted union of all element keys across *entries*."""
    return sorted({elem for e in entries for elem in e.species.atoms})


def coerce_and_validate(
    stoichiometry: Sequence[StoichiometryEntry],
    balance_elements: Optional[Sequence[str]],
    balance_atol: float,
) -> List[StoichiometryEntry]:
    """Coerce a stoichiometry to a list and validate elemental closure.

    Parameters
    ----------
    stoichiometry : sequence of StoichiometryEntry
        Reaction participants.
    balance_elements : sequence of str or None
        Elements to check via :func:`validate_balance`.  ``None`` means
        infer from the atoms present across all entries — every element
        that appears in any participant is checked.
    balance_atol : float
        Absolute tolerance for the residual (mol).

    Returns
    -------
    list of StoichiometryEntry
        The materialised, validated stoichiometry list.

    Raises
    ------
    StoichiometryError
        If any element's residual exceeds ``balance_atol``.
    """
    entries = list(stoichiometry)
    elements = _infer_elements(entries) if balance_elements is None else balance_elements
    validate_balance(entries, elements, balance_atol)
    return entries


def species_ids_from_entries(
    entries: Sequence[StoichiometryEntry],
) -> List[str]:
    """Return the sorted, unique species IDs present in *entries*."""
    return sorted({e.species.id for e in entries})


def phases_from_entries(
    entries: Sequence[StoichiometryEntry],
) -> List[str]:
    """Return the sorted, unique phase keys present in *entries*."""
    return sorted({e.phase for e in entries})


def show_balance_from_entries(
    entries: List[StoichiometryEntry],
    label: str,
    elements: Optional[Sequence[str]],
) -> None:
    """Print a per-element residual table for *entries*.

    Shared implementation backing :meth:`EquilibriumReaction.show_balance`
    and :meth:`KineticReaction.show_balance`.

    Parameters
    ----------
    entries : list of StoichiometryEntry
        Stoichiometry to inspect.
    label : str
        Reaction label, printed as a header.
    elements : sequence of str or None
        Elements to check.  ``None`` means discover from the atoms present
        across all entries.
    """
    if elements is None:
        _elements: Sequence[str] = sorted(
            {elem for e in entries for elem in e.species.atoms}
        )
    else:
        _elements = elements

    print(f"{label}:")
    for elem in _elements:
        residual = sum(
            e.coefficient * float(e.species.atoms.get(elem, 0))
            for e in entries
        )
        parts = [
            f"{e.species.id}: {e.coefficient:+.3g} × {float(e.species.atoms[elem]):.4g}"
            for e in entries
            if elem in e.species.atoms
        ]
        status = "OK" if abs(residual) < 1e-10 else f"FAIL (residual={residual:+.2e})"
        print(
            f"  {elem:2s}  residual={residual:+.2e}  [{status}]"
            + (f"  {',  '.join(parts)}" if parts else "")
        )


def fmt_stoichiometry_string(
    entries: Sequence[StoichiometryEntry],
    arrow: str = "->",
) -> str:
    """Return a human-readable stoichiometry string for *entries*.

    Shared implementation backing :meth:`KineticReaction.show`,
    :meth:`EquilibriumReaction.show`, and
    :meth:`ReactionSystem.show_reactions`.

    Parameters
    ----------
    entries : sequence of StoichiometryEntry
        Stoichiometry to format.
    arrow : str
        Arrow string separating reactants from products (``"->"`` for
        kinetic reactions, ``"<->"`` for equilibrium reactions).
    """
    _PHASE_ALIAS = {"liquid": "aq", "gas": "g", "solid": "s"}

    def _term(e: StoichiometryEntry) -> str:
        alias = _PHASE_ALIAS.get(e.phase, e.phase)
        c = abs(e.coefficient)
        return (
            f"{e.species.id},{alias}"
            if abs(c - 1.0) < 1e-9
            else f"{c:.3g} {e.species.id},{alias}"
        )

    reactants = [e for e in entries if e.coefficient < 0]
    products  = [e for e in entries if e.coefficient > 0]
    lhs = " + ".join(_term(e) for e in reactants)
    rhs = " + ".join(_term(e) for e in products)
    return f"{lhs} {arrow} {rhs}"


def is_cross_phase_from_entries(
    entries: Sequence[StoichiometryEntry],
) -> bool:
    """Return ``True`` iff *entries* span more than one phase.

    Used by both reaction classes to expose ``is_cross_phase`` as a
    derived property. Cross-phase equilibrium reactions are
    partition declarations consumed by
    :class:`~PyOMES.core.gas_liquid_link.KineticGasLiquidLink`, not by
    :class:`~PyOMES.chemical_equilibrium.engine.BisectionChemicalEquilibriumEngine`.
    """
    return len({e.phase for e in entries}) > 1
