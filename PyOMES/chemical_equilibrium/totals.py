# -*- coding: utf-8 -*-
"""The engines' response to a negative component total.

An engine distributes each component's total among the component's species. A
total below zero cannot be distributed. Both engines read their totals from a
phase's amounts through :func:`checked_total`, which treats a negative total in
one of two ways:

- within tolerance of zero, it is read as zero. The tolerance is the larger of
  an absolute one, ``PyOMES.config.negative_total_tolerance_mol`` (default
  1e-9 mol, the adaptive solver's default ``atol``), and a relative one,
  :data:`ROUNDOFF_RTOL` of the sum of the magnitudes of the amounts the total
  was summed from. The absolute tolerance covers an integrator that takes a
  species to zero and leaves trial totals within its own tolerance of zero on
  either side; the relative one covers roundoff when large amounts cancel.
  Reading such a total as zero changes it by no more than the tolerance;
- beyond that, :class:`NegativeTotalError` is raised, naming the component and
  the amounts. Nothing is changed to make the total fit.
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


def _absolute_tolerance() -> float:
    # Read the live singleton each call (the PyOMES package re-exports the
    # config instance under the submodule's name).
    from PyOMES.config import config
    return float(config.negative_total_tolerance_mol)


def checked_total(component: str, amounts: Iterable[float], *, engine: str) -> float:
    """The sum of *amounts* (mol), as an engine may use it.

    Returns the sum when it is not negative, ``0.0`` when it is negative within
    tolerance, and raises :class:`NegativeTotalError` otherwise.
    """
    values = [float(a) for a in amounts]
    total = sum(values)
    if total >= 0.0:
        return total
    gross = sum(abs(v) for v in values)
    absolute = _absolute_tolerance()
    if -total <= max(absolute, ROUNDOFF_RTOL * gross):
        return 0.0
    raise NegativeTotalError(
        f"{engine} engine: the total of {component!r} is {total:.6e} mol, below "
        f"zero by more than the tolerance (summed from {len(values)} amount(s) "
        f"whose magnitudes add up to {gross:.6e} mol; tolerance "
        f"{absolute:.1e} mol, PyOMES.config.negative_total_tolerance_mol). An "
        f"engine cannot distribute a negative amount among species, and nothing "
        f"is changed to make it fit. A smaller step, a tighter solver tolerance, "
        f"or a rate law that limits itself, keeps the amounts non-negative."
    )
