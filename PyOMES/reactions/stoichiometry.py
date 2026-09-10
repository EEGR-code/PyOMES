# -*- coding: utf-8 -*-
"""Stoichiometric entry definitions and elemental balance validation.

This module provides the foundational data structures for defining reaction
stoichiometry and validating that reactions close elementally.

A :class:`StoichiometryEntry` represents one participant in a reaction
(substrate consumed, product formed, etc.).  Each entry carries a
reference to a :class:`~PyOMES.chemistry.Species` object — atoms,
charge, and molecular weight live on the species, not on the entry —
plus the entry's phase and stoichiometric coefficient.

:func:`validate_balance` checks that the stoichiometry is mass-consistent
for a specified set of elements (and optionally for charge),
raising :class:`StoichiometryError` with detailed diagnostics if it
does not close.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

from ..chemistry.species import Species


class StoichiometryError(ValueError):
    """Raised when a reaction stoichiometry fails elemental balance validation.

    Attributes
    ----------
    element : str
        The element that does not balance, or ``"charge"`` for the
        charge-conservation pass.
    residual : float
        Net excess (positive) or deficit (negative) in mol of the element.
    details : str
        Per-species breakdown showing each participant's contribution.
    """

    def __init__(self, message: str, *, element: str = "", residual: float = 0.0,
                 details: str = ""):
        super().__init__(message)
        self.element = element
        self.residual = residual
        self.details = details


@dataclass(frozen=True)
class StoichiometryEntry:
    """One participant in a reaction.

    Parameters
    ----------
    species : Species
        The chemical species. Atoms, charge, and molecular weight live
        on the species — no per-entry duplication.
    phase : str
        Target phase key: ``"liquid"``, ``"gas"``, or any CV phase key.
        For biological state that is transported with the liquid but not
        equilibrated (e.g. biomass), use ``"liquid"`` — the equilibrium
        solver will ignore species for which it has no Henry constant.
    coefficient : float
        Stoichiometric coefficient per unit extent of reaction.
        **Positive = produced, negative = consumed.**
    """

    species: Species
    phase: str
    coefficient: float


def validate_balance(
    entries: Sequence[StoichiometryEntry],
    elements: Sequence[str] = ("C", "H", "O"),
    atol: float = 1e-10,
    *,
    check_charge: bool = False,
) -> None:
    """Validate that a set of stoichiometric entries closes elementally.

    For each element in *elements*, computes the net balance:

        Σ_i  (coefficient_i × atoms_i[element])

    If this sum exceeds *atol* for any element, a :class:`StoichiometryError`
    is raised with a detailed per-species breakdown.

    When ``check_charge=True``, the same iteration pattern runs against
    ``entry.species.charge`` as a peer conservation pass.  One mechanism,
    two laws — element conservation and charge conservation.

    Parameters
    ----------
    entries : sequence of StoichiometryEntry
        All participants in the reaction.
    elements : sequence of str
        Elements to check (default: C, H, O).
    atol : float
        Absolute tolerance for the residual (mol).
    check_charge : bool, optional
        If True, additionally validate that the net charge is zero
        (within ``atol``).  Default False — many reactions are written
        in net-neutral form where the H⁺ that balances charge is
        absorbed into the property channel rather than the
        stoichiometry, so charge balance is opt-in.

    Raises
    ------
    StoichiometryError
        If any element's residual (or charge, when checked) exceeds
        ``atol``.
    """
    for elem in elements:
        residual = 0.0
        parts = []
        for e in entries:
            atom_count = float(e.species.atoms.get(elem, 0.0))
            contribution = e.coefficient * atom_count
            residual += contribution
            if atom_count != 0.0:
                parts.append(
                    f"  {e.species.id}: coeff={e.coefficient:+.6g}, "
                    f"{elem}={atom_count:.4g}, contribution={contribution:+.6e}"
                )

        if abs(residual) > atol:
            details = "\n".join(parts)
            raise StoichiometryError(
                f"Element '{elem}' does not balance: residual = {residual:+.6e} mol.\n"
                f"Participants:\n{details}",
                element=elem,
                residual=residual,
                details=details,
            )

    if check_charge:
        residual = 0.0
        parts = []
        for e in entries:
            q = float(e.species.charge)
            contribution = e.coefficient * q
            residual += contribution
            if q != 0.0:
                parts.append(
                    f"  {e.species.id}: coeff={e.coefficient:+.6g}, "
                    f"charge={q:+g}, contribution={contribution:+.6e}"
                )
        if abs(residual) > atol:
            details = "\n".join(parts)
            raise StoichiometryError(
                f"Charge does not balance: residual = {residual:+.6e}.\n"
                f"Participants:\n{details}",
                element="charge",
                residual=residual,
                details=details,
            )


# ── String stoichiometry parser ───────────────────────────────────────────────

_PHASE_ALIASES: Dict[str, str] = {
    "aq": "liquid",
    "l":  "liquid",
    "g":  "gas",
    "s":  "solid",
}

_COEFF_RE = re.compile(r"^(\d+(?:\.\d*)?|\.\d+)\s+")


def _get_common_species() -> Dict[str, Species]:
    """Return all Species objects from chemistry.common_species keyed by id."""
    from ..chemistry import common_species as _cs_mod
    return {
        v.id: v
        for v in vars(_cs_mod).values()
        if isinstance(v, Species)
    }


def _parse_stoichiometry(
    s: str,
    species: Optional[Dict[str, Species]],
    *,
    reaction_type: Optional[str] = None,
) -> List[StoichiometryEntry]:
    """Parse a human-readable stoichiometry string into StoichiometryEntry objects.

    Format::

        [coeff] species_id,phase [+ [coeff] species_id,phase ...] <-> ...

    Reactant terms (left of arrow) receive negated coefficients; product
    terms (right of arrow) receive positive coefficients.

    Parameters
    ----------
    s : str
        Stoichiometry string, e.g.
        ``"CO2,aq + H2O,aq <-> HCO3-,aq + H+,aq"``.
    species : dict[str, Species] or None
        Caller-supplied species for locally declared IDs.  Looked up after
        ``common_species``; caller entries override common ones.
    reaction_type : str or None
        ``"equilibrium"`` or ``"kinetic"``.  When provided, the arrow
        direction is validated against the calling reaction class.
    """
    s = s.strip()

    # ── detect arrow ─────────────────────────────────────────────────────────
    # Replace <-> with placeholder so any remaining -> is unambiguously kinetic.
    s_tmp = s.replace("<->", "\x00")
    has_equil = "\x00" in s_tmp
    has_kinetic = "->" in s_tmp

    if has_equil and has_kinetic:
        raise ValueError(
            f"Stoichiometry string contains both '<->' and '->'; "
            f"use exactly one arrow: {s!r}"
        )
    if not has_equil and not has_kinetic:
        raise ValueError(
            f"No arrow token found in stoichiometry string: {s!r}. "
            "Use '<->' for EquilibriumReaction or '->' for KineticReaction."
        )

    if has_equil:
        arrow = "<->"
        parts = s.split("<->")
    else:
        arrow = "->"
        parts = s.split("->")

    if len(parts) != 2:
        raise ValueError(
            f"Expected exactly one '{arrow}' in stoichiometry string, "
            f"found {len(parts) - 1}: {s!r}"
        )

    # ── validate arrow direction vs reaction class ────────────────────────────
    if reaction_type == "equilibrium" and arrow == "->":
        raise ValueError(
            "Used '->' in an EquilibriumReaction stoichiometry string. "
            "Use '<->' for equilibrium reactions, "
            "or switch to KineticReaction for unidirectional reactions."
        )
    if reaction_type == "kinetic" and arrow == "<->":
        raise ValueError(
            "Used '<->' in a KineticReaction stoichiometry string. "
            "Use '->' for kinetic (unidirectional) reactions, "
            "or switch to EquilibriumReaction."
        )

    lhs_str, rhs_str = parts[0].strip(), parts[1].strip()

    # ── build species lookup ──────────────────────────────────────────────────
    common = _get_common_species()
    lookup: Dict[str, Species] = dict(common)
    if species:
        lookup.update(species)

    # ── parse one side of the arrow ───────────────────────────────────────────
    def _parse_side(side_str: str, sign: float) -> List[StoichiometryEntry]:
        # Split on whitespace-surrounded '+' to avoid splitting species like H+.
        raw_terms = re.split(r"\s+\+\s+", side_str.strip())
        result: List[StoichiometryEntry] = []
        for raw in raw_terms:
            raw = raw.strip()
            if not raw:
                continue

            last_comma = raw.rfind(",")
            if last_comma == -1:
                raise ValueError(
                    f"Missing phase suffix in term {raw!r}. "
                    "Every species must have a phase suffix "
                    "(e.g. ',aq', ',l', ',g', ',s')."
                )

            prefix = raw[:last_comma].strip()
            phase_alias = raw[last_comma + 1:].strip()

            if phase_alias not in _PHASE_ALIASES:
                raise ValueError(
                    f"Unknown phase suffix ',{phase_alias}' in term {raw!r}. "
                    f"Valid options: {sorted(_PHASE_ALIASES)}."
                )
            phase = _PHASE_ALIASES[phase_alias]

            m = _COEFF_RE.match(prefix)
            if m:
                coeff = float(m.group(1))
                if coeff <= 0:
                    raise ValueError(
                        f"Stoichiometric coefficient must be positive, "
                        f"got {coeff} in term {raw!r}."
                    )
                sp_id = prefix[m.end():].strip()
            else:
                if re.match(r"^\d", prefix):
                    raise ValueError(
                        f"Malformed coefficient in term {raw!r}: "
                        "coefficient must be separated from the species ID by a space."
                    )
                coeff = 1.0
                sp_id = prefix

            if not sp_id:
                raise ValueError(f"Empty species ID in term {raw!r}.")

            if sp_id not in lookup:
                extra = (
                    " Pass a 'species' dict for locally-declared species."
                    if not species
                    else f" Caller-supplied IDs: {sorted(species)}."
                )
                raise ValueError(
                    f"Species ID {sp_id!r} not found.{extra} "
                    f"Common species available: {sorted(common)}."
                )

            result.append(StoichiometryEntry(
                species=lookup[sp_id],
                phase=phase,
                coefficient=sign * coeff,
            ))
        return result

    entries: List[StoichiometryEntry] = []
    entries.extend(_parse_side(lhs_str, -1.0))
    entries.extend(_parse_side(rhs_str, +1.0))
    return entries
