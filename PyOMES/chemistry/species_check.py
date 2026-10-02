# -*- coding: utf-8 -*-
"""Cross-reaction species consistency check.

When reaction sets are composed from multiple model files, two failure
modes can occur silently:

- **Hard conflict:** the same ``id`` is declared with different
  ``atoms`` or ``charge`` in different files. The composition is
  internally inconsistent — a ``KineticReaction`` may produce ``CO2``
  with ``{"C":1,"O":2}`` while the speciation engine consumes the same
  ``"CO2"`` key as ``{"C":1,"O":1}``. Always an error.
- **Soft conflict:** the same id is declared with the same ``atoms``
  and ``charge`` but as distinct ``Species`` objects (because the
  user redeclared instead of importing). Functionally correct (because
  ``Species`` has value equality), but signals that future divergence
  is one careless edit away. Warning by default.

:func:`check_species_consistency` walks every ``StoichiometryEntry`` in
a list of reactions, groups by id, and applies the rules above.
:func:`merge_species` combines species given directly (e.g. a model's own
``species=`` and a database's ``species``) into one dict, raising on the
same id with different data.
:class:`~PyOMES.core.control_volume.ControlVolume` runs this check
automatically in ``__init__`` after a reaction model is attached.

Usage
-----
>>> from PyOMES.chemistry import check_species_consistency
>>> check_species_consistency(rxn_set.reactions)
"""

from __future__ import annotations

import warnings
from typing import Dict, Iterable, List, Mapping, Union

from .species import Species, SpeciesConflictError


_VALID_SOFT = ("warn", "raise", "ignore")


def check_species_consistency(
    reactions: Iterable,
    *,
    soft_conflicts: str = "warn",
) -> None:
    """Validate cross-reaction ``Species`` consistency.

    Parameters
    ----------
    reactions : iterable
        Anything yielding objects with a ``stoichiometry`` attribute,
        whose entries each carry a ``species: Species`` field.
        ``ReactionSystem`` and a list of reaction declarations
        (``KineticReaction`` / ``EquilibriumReaction``) both work.
    soft_conflicts : {"warn", "raise", "ignore"}
        Policy when distinct ``Species`` objects with equal data share
        an id. ``"warn"`` (default) emits a ``UserWarning``,
        ``"raise"`` promotes to ``SpeciesConflictError``, ``"ignore"``
        suppresses the diagnostic.

    Raises
    ------
    SpeciesConflictError
        On any hard conflict (same id, different atoms/charge), or on
        any soft conflict when ``soft_conflicts="raise"``.
    """
    if soft_conflicts not in _VALID_SOFT:
        raise ValueError(
            f"soft_conflicts must be one of {_VALID_SOFT}, "
            f"got {soft_conflicts!r}"
        )

    # Group every Species object by id.
    by_id: dict = {}
    for rxn in reactions:
        stoich = getattr(rxn, "stoichiometry", None)
        if stoich is None:
            continue
        for entry in stoich:
            sp = getattr(entry, "species", None)
            if sp is None:
                continue
            by_id.setdefault(sp.id, []).append(sp)

    for sp_id, instances in by_id.items():
        if len(instances) < 2:
            continue
        # Reduce by object identity first — a single object referenced
        # 100 times is the happy path.
        unique_objs = []
        seen_ids = set()
        for sp in instances:
            if id(sp) not in seen_ids:
                seen_ids.add(id(sp))
                unique_objs.append(sp)
        if len(unique_objs) == 1:
            continue

        # Multiple distinct objects share this id. Check whether they
        # are value-equal (soft) or differ (hard).
        ref = unique_objs[0]
        hard = []
        for other in unique_objs[1:]:
            if (dict(other.atoms) != dict(ref.atoms)
                    or other.charge != ref.charge
                    or other.MW != ref.MW):
                hard.append(other)

        if hard:
            details = _format_hard_conflict(sp_id, ref, hard)
            raise SpeciesConflictError(
                f"Species id {sp_id!r} declared with conflicting data "
                f"across reactions:\n{details}",
                species_id=sp_id,
                details=details,
            )

        # All soft — same data, distinct objects.
        if soft_conflicts == "ignore":
            continue
        msg = (
            f"Species id {sp_id!r} declared in {len(unique_objs)} places "
            f"with equal data but as distinct Species objects. Import "
            f"from a shared module instead of redeclaring."
        )
        if soft_conflicts == "raise":
            raise SpeciesConflictError(msg, species_id=sp_id, details=msg)
        warnings.warn(msg, UserWarning, stacklevel=2)


def merge_species(
    *sources: Union[Mapping[str, Species], Iterable[Species]],
) -> Dict[str, Species]:
    """Merge species from several sources into one ``{id: Species}`` dict.

    Each source is a ``{id: Species}`` mapping (such as a database's
    ``species``) or an iterable of ``Species``. The same id may appear in
    several sources only with equal data; the first object seen is kept.

    Raises
    ------
    SpeciesConflictError
        If one id arrives with different atoms, charge or MW.
    ValueError
        If a mapping holds a ``Species`` under a key other than its id.
    """
    merged: Dict[str, Species] = {}
    for source in sources:
        if isinstance(source, Mapping):
            for key, sp in source.items():
                if key != sp.id:
                    raise ValueError(
                        f"Species {sp.id!r} is listed under the key {key!r}; "
                        "a species mapping must key each Species by its id."
                    )
            items = source.values()
        else:
            items = source
        for sp in items:
            ref = merged.get(sp.id)
            if ref is None:
                merged[sp.id] = sp
            elif ref is not sp and (
                dict(ref.atoms) != dict(sp.atoms)
                or ref.charge != sp.charge
                or ref.MW != sp.MW
            ):
                details = _format_hard_conflict(sp.id, ref, [sp])
                raise SpeciesConflictError(
                    f"Species id {sp.id!r} is given twice with different "
                    f"data:\n{details}",
                    species_id=sp.id,
                    details=details,
                )
    return merged


def _format_hard_conflict(sp_id: str, ref: Species, others: List[Species]) -> str:
    """Render a per-instance breakdown of a hard conflict."""
    lines = [f"  ref:    atoms={dict(ref.atoms)}, charge={ref.charge}, MW={ref.MW}"]
    for sp in others:
        lines.append(
            f"  other:  atoms={dict(sp.atoms)}, charge={sp.charge}, MW={sp.MW}"
        )
    return "\n".join(lines)
