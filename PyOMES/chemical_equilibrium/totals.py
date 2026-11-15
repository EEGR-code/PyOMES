# -*- coding: utf-8 -*-
"""The engines' response to a negative component total.

An engine distributes each component's total among the component's species. A
total below zero cannot be distributed. Both engines read their totals from a
phase's amounts through :func:`checked_total`, which treats a negative total in
one of two ways:

- within roundoff of zero (its magnitude at most :data:`ROUNDOFF_RTOL` of the
  sum of the magnitudes of the amounts it was summed from), it is read as zero,
  which changes it by no more than roundoff;
- beyond that, :class:`NegativeTotalError` is raised, naming the component and
  the amounts. Nothing is changed to make the total fit.

A total of a single species (a strong ion) has no roundoff margin of its own, so
any negative amount of it raises.
"""
from __future__ import annotations

from typing import Iterable

#: A negative total within this fraction of the sum of its amounts' magnitudes
#: is roundoff and reads as zero.
ROUNDOFF_RTOL = 1e-12


class NegativeTotalError(RuntimeError):
    """A component total an engine was asked to distribute is below zero.

    A ``RuntimeError`` rather than a ``ValueError``: code that treats a
    ``ValueError`` from a solve as a failure to recover from (the simulation's
    first solve warns and continues) does not swallow it.
    """


def checked_total(component: str, amounts: Iterable[float], *, engine: str) -> float:
    """The sum of *amounts* (mol), as an engine may use it.

    Returns the sum when it is not negative, ``0.0`` when it is negative within
    roundoff, and raises :class:`NegativeTotalError` otherwise.
    """
    values = [float(a) for a in amounts]
    total = sum(values)
    if total >= 0.0:
        return total
    gross = sum(abs(v) for v in values)
    if -total <= ROUNDOFF_RTOL * gross:
        return 0.0
    raise NegativeTotalError(
        f"{engine} engine: the total of {component!r} is {total:.6e} mol, below "
        f"zero by more than roundoff (summed from {len(values)} amount(s) whose "
        f"magnitudes add up to {gross:.6e} mol). An engine cannot distribute a "
        f"negative amount among species, and nothing is changed to make it fit. "
        f"A smaller step, or a rate law that limits itself, keeps the amounts "
        f"non-negative."
    )
