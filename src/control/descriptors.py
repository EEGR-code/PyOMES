# -*- coding: utf-8 -*-
"""Descriptor-declared mutable attributes with lifecycle gating.

Provides two descriptors for classes whose attributes are mutable
mid-run only by the orchestrator (Pattern B):

* :class:`MutableScalar` — a single typed value (float, int, …).
* :class:`MutableDict` — a ``{key: value}`` dict updated per-key or
  wholesale.

Each descriptor carries the storage attribute, lifecycle gate, and
(in C2) the class-level ``ParamPath`` builder in one declaration,
replacing the ``@property / @x.setter / _set_x_unchecked`` triplet.

Class-level access returns the descriptor itself (matching SQLAlchemy's
``User.name → Column`` pattern); C2 promotes the descriptor directly as
a ``ParamPath`` handle.  For ``MutableDict``, ``Class.attr["key"]``
returns a :class:`_DictItemRef` stub that C2 promotes to a ``ParamPath``
leaf.

Storage convention: the backing attribute is ``f"_{name}"`` — e.g.
``T_K = MutableScalar(float, positive=True)`` stores in ``self._T_K``.
This mirrors the existing ``self._T_K`` naming used by the ``@property``
triplets being replaced, ensuring any direct ``obj._T_K`` reads continue
to work unchanged.
"""

from __future__ import annotations

from typing import Any, Optional, Type


_UNSET = object()


# ══════════════════════════════════════════════════════════════════════
#  MutableScalar
# ══════════════════════════════════════════════════════════════════════

class MutableScalar:
    """Descriptor for a scalar mutable attribute with lifecycle gating.

    Parameters
    ----------
    type_ : callable
        Type coercion applied on every write (e.g. ``float``, ``int``).
    positive : bool
        If ``True``, raises :exc:`ValueError` on values ``<= 0`` after
        coercion.
    non_negative : bool
        If ``True``, raises :exc:`ValueError` on values ``< 0`` after
        coercion.  Ignored when ``positive=True``.

    Usage
    -----
    ::

        class GasPhase:
            T_K = MutableScalar(float, positive=True)
            V_L = MutableScalar(float, positive=True)

        # Instance-level: returns the stored value.
        phase.T_K          # -> 308.15

        # Class-level: returns the descriptor itself.
        GasPhase.T_K       # -> MutableScalar instance (the descriptor)
    """

    def __init__(
        self,
        type_: Type = float,
        *,
        positive: bool = False,
        non_negative: bool = False,
    ) -> None:
        self.type_ = type_
        self.positive = positive
        self.non_negative = non_negative
        self.name: Optional[str] = None
        self.storage_attr: Optional[str] = None

    def __set_name__(self, owner: type, name: str) -> None:
        self.name = name
        self.storage_attr = f"_{name}"
        self._owner = owner

    def __get__(self, obj: Any, objtype: Optional[type] = None) -> Any:
        if obj is None:
            # Class-level access: return the descriptor itself (C2 uses it
            # directly as a ParamPath handle, matching the SQLAlchemy pattern).
            return self
        val = obj.__dict__.get(self.storage_attr, _UNSET)
        if val is _UNSET:
            raise AttributeError(
                f"{type(obj).__name__}.{self.name} has not been initialised."
            )
        return val

    def __set__(self, obj: Any, value: Any) -> None:
        if isinstance(value, MutableScalar):
            # Dataclass passing the descriptor itself as a default — skip.
            return
        from PyOMES.core.lifecycle import raise_if_running
        raise_if_running(obj, self.name)
        self._set_unchecked(obj, value)

    def _set_unchecked(self, obj: Any, value: Any) -> None:
        """Orchestrator-mediated unchecked path (Pattern B)."""
        obj.__dict__[self.storage_attr] = self._validate(value)

    def _validate(self, value: Any) -> Any:
        value = self.type_(value)
        if self.positive and value <= 0:
            raise ValueError(
                f"{self.name} must be positive, got {value!r}."
            )
        if self.non_negative and not self.positive and value < 0:
            raise ValueError(
                f"{self.name} must be non-negative, got {value!r}."
            )
        return value


# ══════════════════════════════════════════════════════════════════════
#  MutableDict
# ══════════════════════════════════════════════════════════════════════

class MutableDict:
    """Descriptor for a dict mutable attribute with lifecycle gating.

    Supports whole-dict replacement (gated ``__set__``) and single-key
    update via :meth:`_set_item_unchecked` (orchestrator path).

    Parameters
    ----------
    value_type : callable
        Type coercion applied to each value on write.  Default ``float``.

    Usage
    -----
    ::

        class KineticGasLiquidLink:
            kLa: Dict[str, float] = MutableDict(value_type=float)

        # When used as a dataclass default value, @dataclass generates:
        #   self.kLa = kLa_arg
        # which goes through MutableDict.__set__, storing in self._kLa.

        # Instance-level: returns the underlying dict.
        link.kLa           # -> {"O2": 150.0, "CO2": 135.0}

        # Class-level: returns the descriptor itself.
        KineticGasLiquidLink.kLa       # -> MutableDict instance (the descriptor)
        KineticGasLiquidLink.kLa["O2"] # -> _DictItemRef(...)
    """

    def __init__(self, *, value_type: Type = float) -> None:
        self.value_type = value_type
        self.name: Optional[str] = None
        self.storage_attr: Optional[str] = None

    def __set_name__(self, owner: type, name: str) -> None:
        self.name = name
        self.storage_attr = f"_{name}"
        self._owner: Optional[type] = owner

    def __get__(self, obj: Any, objtype: Optional[type] = None) -> Any:
        if obj is None:
            # Class-level access: return the descriptor itself so that
            # @dataclass can read it back via getattr(cls, name) and
            # isinstance(value, MutableDict) fires correctly in __set__.
            return self
        return obj.__dict__.get(self.storage_attr, {})

    def __getitem__(self, key: str) -> "_DictItemRef":
        """``Class.attr["key"]`` → typed path stub for a specific dict key."""
        return _DictItemRef(self, self._owner, key)

    def __set__(self, obj: Any, value: Any) -> None:
        if isinstance(value, MutableDict):
            # Dataclass passing the descriptor itself as the default —
            # initialise to empty dict without gating.
            obj.__dict__[self.storage_attr] = {}
            return
        from PyOMES.core.lifecycle import raise_if_running
        raise_if_running(obj, self.name)
        self._set_unchecked(obj, value)

    def _set_unchecked(self, obj: Any, value: Any) -> None:
        """Replace the whole dict. Orchestrator-mediated unchecked path."""
        obj.__dict__[self.storage_attr] = dict(value)

    def _set_item_unchecked(self, obj: Any, key: str, value: Any) -> None:
        """Update a single key. Orchestrator-mediated unchecked path."""
        d = obj.__dict__.get(self.storage_attr)
        if d is None:
            d = {}
            obj.__dict__[self.storage_attr] = d
        d[key] = self.value_type(value)


# ══════════════════════════════════════════════════════════════════════
#  Class-level reference stub (C2 promotes to ParamPath)
# ══════════════════════════════════════════════════════════════════════
#
# Class-level access to a MutableScalar or MutableDict returns the
# descriptor itself (matching SQLAlchemy's User.name → Column pattern).
# The only additional stub needed is _DictItemRef for a descriptor +
# specific key, returned by  MyClass.attr["key"].

class _DictItemRef:
    """A ``_DictRef`` with a specific key selected.

    Returned by ``MyClass.attr["key"]``.  C2 promotes to a full
    ``ParamPath`` leaf.
    """

    def __init__(
        self,
        descriptor: MutableDict,
        owner: Optional[type],
        key: str,
    ) -> None:
        self._descriptor = descriptor
        self._owner = owner
        self._key = key

    def __repr__(self) -> str:
        owner_name = self._owner.__name__ if self._owner else "?"
        return (
            f"_DictItemRef({owner_name}.{self._descriptor.name}"
            f"[{self._key!r}])"
        )
