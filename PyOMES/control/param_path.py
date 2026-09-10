# -*- coding: utf-8 -*-
"""Typed path representation for the reflective param-change dispatcher.

Provides:

* :exc:`ParamPathError` — raised by the walker (C3) on unknown or
  malformed paths (Q3: fail-loud, no silent ignore).
* :class:`AttrSegment`, :class:`ListSelector`, :class:`IndexSelector`
  — frozen, hashable segment types that make up a :class:`ParamPath`.
* :class:`ParamPath` — an immutable, hashable sequence of segments
  parsed from a string path or built from descriptor references.

String path syntax (Q1 — full reflective form, no alias table)::

    boundaries[GasFeed].vvm_min
    internal_interfaces[KineticGasLiquidLink].kLa.O2
    phases.liquid.T_K

List-selector syntax (Q2)::

    [CamelCaseIdent]      — first element where isinstance(elem, Class) is True
    [CamelCaseIdent:N]    — Nth (0-indexed) isinstance match
    [N]                   — pure numeric index (escape hatch)

The walker (C3) lives in :mod:`PyOMES.control.walker`; this module is
purely a data-structure + parser layer.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional, Tuple, Union


# ══════════════════════════════════════════════════════════════════════
#  Exception
# ══════════════════════════════════════════════════════════════════════

class ParamPathError(Exception):
    """Raised when a param-change path cannot be resolved on a CV.

    Includes the raw path string and the CV label in the message so
    mistyped paths are caught immediately rather than silently mutating
    nothing (Q3: error, not warn-once or ignore).
    """


# ══════════════════════════════════════════════════════════════════════
#  Segment types
# ══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class AttrSegment:
    """Traverse one step by attribute name (or dict key as fallback).

    At walk time the resolver tries ``getattr(target, name)`` first; if
    that raises ``AttributeError`` and *target* is a ``dict``, it falls
    back to ``target[name]``.  This lets ``phases.liquid.T_K`` work
    even though ``cv.phases`` is a ``dict``.
    """

    name: str

    def resolve(self, target: object) -> object:
        """Return the next target by attribute name (dict key as fallback).

        Raises :exc:`ParamPathError` if neither ``getattr`` nor dict lookup
        succeeds.
        """
        try:
            return getattr(target, self.name)
        except AttributeError:
            pass
        if isinstance(target, dict) and self.name in target:
            return target[self.name]
        raise ParamPathError(
            f"Cannot resolve {self.name!r} on "
            f"{type(target).__name__!r}."
        )

    def __repr__(self) -> str:
        return f".{self.name}"


@dataclass(frozen=True)
class ListSelector:
    """Select an element from a list by type name (and optional ordinal).

    Parameters
    ----------
    type_name : str
        CamelCase name looked up in the class registry at walk time
        (Q2: v1 registry scans ``PyOMES.core.__all__``).
    index : int or None
        ``None`` — first matching element.
        ``N`` — Nth (0-indexed) matching element.
    """

    type_name: str
    index: Optional[int] = None

    def resolve(self, target: object, class_registry: dict) -> object:
        """Return the first (or Nth) element of *target* that is an instance
        of the registered class named :attr:`type_name`.

        Raises :exc:`ParamPathError` if the type name is unknown, *target*
        is not a list, or the Nth match does not exist.
        """
        if not isinstance(target, (list, tuple)):
            raise ParamPathError(
                f"ListSelector[{self.type_name!r}] applied to "
                f"{type(target).__name__!r}, not a list."
            )
        cls = class_registry.get(self.type_name)
        if cls is None:
            raise ParamPathError(
                f"Class {self.type_name!r} not found in class registry. "
                f"Available: {sorted(class_registry.keys())}"
            )
        matches = [elem for elem in target if isinstance(elem, cls)]
        idx = 0 if self.index is None else self.index
        if idx >= len(matches):
            raise ParamPathError(
                f"ListSelector[{self.type_name!r}:{idx}]: found "
                f"{len(matches)} match(es) in list of {len(target)}; "
                f"index {idx} out of range."
            )
        return matches[idx]

    def __repr__(self) -> str:
        if self.index is None:
            return f"[{self.type_name}]"
        return f"[{self.type_name}:{self.index}]"


@dataclass(frozen=True)
class IndexSelector:
    """Select a list element by pure numeric index (escape hatch, Q2)."""

    index: int

    def resolve(self, target: object) -> object:
        """Return ``target[index]``, raising :exc:`ParamPathError` on failure."""
        try:
            return target[self.index]  # type: ignore[index]
        except (IndexError, TypeError) as exc:
            raise ParamPathError(
                f"IndexSelector[{self.index}] on {type(target).__name__!r}: "
                f"{exc}"
            ) from exc

    def __repr__(self) -> str:
        return f"[{self.index}]"


Segment = Union[AttrSegment, ListSelector, IndexSelector]


# ══════════════════════════════════════════════════════════════════════
#  Class registry (for ListSelector resolution at walk time)
# ══════════════════════════════════════════════════════════════════════

_CLASS_REGISTRY: Optional[dict] = None


def _get_class_registry() -> dict:
    """Lazy registry: CamelCase name → class, scanned from PyOMES.core.__all__.

    Q2: v1 closed to plugins; a ``@register_param_class`` decorator can widen
    this if a plugin ecosystem emerges.
    """
    global _CLASS_REGISTRY
    if _CLASS_REGISTRY is None:
        _CLASS_REGISTRY = {}
        try:
            import PyOMES.core as _core
            for _name in getattr(_core, "__all__", []):
                _obj = getattr(_core, _name, None)
                if isinstance(_obj, type):
                    _CLASS_REGISTRY[_name] = _obj
        except ImportError:
            pass
    return _CLASS_REGISTRY


# ══════════════════════════════════════════════════════════════════════
#  Parser internals
# ══════════════════════════════════════════════════════════════════════

# Matches a single dot-separated token of the form:
#   identifier            → AttrSegment
#   identifier[selector]  → AttrSegment + List/IndexSelector
_TOKEN_RE = re.compile(
    r'([a-zA-Z_]\w*)'        # group 1: attribute name
    r'(?:\[([^\]]+)\])?'     # group 2: optional selector (contents only)
)

# Matches a pure-numeric selector: [0], [12], …
_INDEX_SEL_RE = re.compile(r'^\d+$')

# Matches a CamelCase type-name selector with optional :N ordinal:
#   GasFeed  /  KineticGasLiquidLink  /  GasFeed:1
_TYPE_SEL_RE = re.compile(r'^([A-Z]\w*)(?::(\d+))?$')


def _parse_token(part: str, path_raw: str) -> Tuple[Segment, ...]:
    """Parse one dot-separated token into one or two segments."""
    m = _TOKEN_RE.fullmatch(part)
    if not m:
        raise ParamPathError(
            f"Invalid path segment {part!r} in {path_raw!r}. "
            f"Expected: identifier or identifier[selector]."
        )

    attr_name: str = m.group(1)
    selector_text: Optional[str] = m.group(2)

    result: list = [AttrSegment(attr_name)]

    if selector_text is not None:
        if _INDEX_SEL_RE.fullmatch(selector_text):
            result.append(IndexSelector(int(selector_text)))
        else:
            m2 = _TYPE_SEL_RE.fullmatch(selector_text)
            if not m2:
                raise ParamPathError(
                    f"Invalid list selector [{selector_text!r}] in "
                    f"{path_raw!r}. Expected: digits, CamelCaseName, "
                    f"or CamelCaseName:N (0-indexed ordinal)."
                )
            type_name: str = m2.group(1)
            idx: Optional[int] = (
                int(m2.group(2)) if m2.group(2) is not None else None
            )
            result.append(ListSelector(type_name, idx))

    return tuple(result)


# ══════════════════════════════════════════════════════════════════════
#  ParamPath
# ══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class ParamPath:
    """An immutable, hashable sequence of traversal segments.

    Constructed from a string path via :meth:`parse` or assembled
    directly from segment objects.  Both forms feed the reflective
    walker (C3) identically.

    Parameters
    ----------
    segments : tuple of Segment
        Non-empty ordered sequence of :class:`AttrSegment`,
        :class:`ListSelector`, and :class:`IndexSelector` objects.
    raw : str
        Original string (if parsed from one).  Used for ``repr`` and
        error messages; empty string when built programmatically.

    Examples
    --------
    >>> p = ParamPath.parse("internal_interfaces[KineticGasLiquidLink].kLa.O2")
    >>> p.segments
    (AttrSegment('internal_interfaces'), ListSelector('KineticGasLiquidLink'),
     AttrSegment('kLa'), AttrSegment('O2'))
    """

    segments: Tuple[Segment, ...]
    raw: str = ""

    # ── Construction ──────────────────────────────────────────────────

    @classmethod
    def parse(cls, s: str) -> "ParamPath":
        """Parse a full reflective path string into a :class:`ParamPath`.

        Parameters
        ----------
        s : str
            Full reflective path, e.g.
            ``"internal_interfaces[KineticGasLiquidLink].kLa.O2"``.

        Raises
        ------
        ParamPathError
            If the string is empty, contains empty segments, or has
            invalid selector syntax.
        """
        s = str(s)
        if not s.strip():
            raise ParamPathError("Path string must not be empty.")

        dot_parts = s.split(".")
        segments: list = []
        for part in dot_parts:
            if not part:
                raise ParamPathError(
                    f"Empty segment (consecutive dots?) in path {s!r}."
                )
            segments.extend(_parse_token(part, s))

        if not segments:
            raise ParamPathError(f"Path {s!r} produced no segments.")

        return cls(segments=tuple(segments), raw=s)

    # ── Accessors ─────────────────────────────────────────────────────

    @property
    def leaf(self) -> Segment:
        """The last segment — the one the walker writes to."""
        return self.segments[-1]

    @property
    def head(self) -> "Tuple[Segment, ...]":
        """All segments except the leaf — the traversal prefix."""
        return self.segments[:-1]

    # ── Representation ────────────────────────────────────────────────

    def __str__(self) -> str:
        if self.raw:
            return self.raw
        return "".join(repr(s) for s in self.segments).lstrip(".")

    def __repr__(self) -> str:
        return f"ParamPath({str(self)!r})"
