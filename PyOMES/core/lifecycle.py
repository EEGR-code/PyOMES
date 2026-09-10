# -*- coding: utf-8 -*-
"""Lifecycle primitives for Simulation runtime gating (C6).

Defines:

* :class:`RunContext` — the run-bound state object shared across all
  objects owned by a single :class:`~PyOMES.core.simulation.Simulation`.
  Created once per Simulation; ``is_running`` flips at run-entry and
  back at run-exit.
* :class:`_LockableList` — a list subclass that consults its
  :class:`RunContext` on every mutating method, raising
  :class:`RuntimeError` if the simulation is running.
* :func:`raise_if_running` — free helper invoked from gated mutators
  on individual objects (``Phase``, ``KineticGasLiquidLink``, etc.).

This module is the *one* home for the gating mechanism. Mutator-side
guards consult ``self._context.is_running`` (set by the Simulation at
run boundaries). The orchestrator-mediated mutation path (Pattern B,
unchecked underscore-prefixed methods) bypasses the gate intentionally
— see SIMULATION_CLASS.md decision 2 / 9.

The module exists separately from
:mod:`~PyOMES.core.simulation` to break what would otherwise be a
circular import chain (``ControlVolume`` and ``Phase`` need
``_LockableList`` / ``RunContext`` references; ``Simulation`` needs
them too).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass
class RunContext:
    """Run-bound state shared across all objects owned by a single Simulation.

    See ``SIMULATION_CLASS.md`` decision 9 for the design rationale.

    Attributes
    ----------
    is_running : bool
        ``True`` while a ``Simulation.run()`` call is in progress.
    label : str
        Human-readable label propagated from the owning Simulation
        for diagnostic error messages.

    Notes
    -----
    Future fields (deferred until consumers appear): ``current_t_h``,
    ``current_step_index``, ``recorder`` — landing when a real
    consumer wants them, not anticipated in C1.
    """

    is_running: bool = False
    label: str = ""


def raise_if_running(obj, attr_label: str) -> None:
    """Free helper: raise ``RuntimeError`` if ``obj`` is in a running Simulation.

    Each gated mutator on a lockable object calls this at entry,
    typically before delegating to a ``_<name>_unchecked`` sibling
    (Pattern B). ``obj`` is duck-typed: must expose a ``_context``
    attribute pointing at a :class:`RunContext` (or ``None``).

    The error message includes the owning Simulation's label so
    diagnostics across ensemble runs are unambiguous.
    """
    ctx = getattr(obj, "_context", None)
    if ctx is not None and ctx.is_running:
        raise RuntimeError(
            f"Cannot mutate {obj.__class__.__name__}.{attr_label} while "
            f"Simulation {ctx.label!r} is running. Use a Controller or "
            f"Profile to modify state mid-run."
        )


class _LockableList(list):
    """A list subclass that refuses mutation while its RunContext is running.

    Used to wrap :attr:`ControlVolume.boundaries`,
    :attr:`ControlVolume.property_calculators`,
    :attr:`Simulation.controllers`, and :attr:`Simulation.profiles`.
    Construction takes the initial iterable and a ``label`` for
    error messages; the ``_context`` is wired by the Simulation at
    its ``__init__``.

    Iterating and reading (``__getitem__``, ``__iter__``, ``len``)
    are unaffected — only mutating methods consult the gate.
    """

    def __init__(self, iterable=(), *, label: str = "<unknown>") -> None:
        super().__init__(iterable)
        self._context: Optional[RunContext] = None
        self._label: str = str(label)

    def _raise_if_running(self, op: str) -> None:
        ctx = self._context
        if ctx is not None and ctx.is_running:
            raise RuntimeError(
                f"Cannot {op} on {self._label} while Simulation "
                f"{ctx.label!r} is running. Use a Controller or "
                f"Profile to modify state mid-run."
            )

    # Mutating methods — all gated.
    def append(self, item) -> None:
        self._raise_if_running("append")
        super().append(item)

    def pop(self, index: int = -1):
        self._raise_if_running("pop")
        return super().pop(index)

    def clear(self) -> None:
        self._raise_if_running("clear")
        super().clear()

    def __setitem__(self, idx, value) -> None:
        self._raise_if_running("__setitem__")
        super().__setitem__(idx, value)

    def __delitem__(self, idx) -> None:
        self._raise_if_running("__delitem__")
        super().__delitem__(idx)

    def extend(self, iterable) -> None:
        self._raise_if_running("extend")
        super().extend(iterable)

    def insert(self, idx: int, item) -> None:
        self._raise_if_running("insert")
        super().insert(idx, item)

    def remove(self, value) -> None:
        self._raise_if_running("remove")
        super().remove(value)

    def __iadd__(self, other):
        self._raise_if_running("__iadd__")
        return super().__iadd__(other)

    def __imul__(self, count):
        self._raise_if_running("__imul__")
        return super().__imul__(count)
