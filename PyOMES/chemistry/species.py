# -*- coding: utf-8 -*-
"""Intrinsic chemical species identity.

A :class:`Species` carries the identity of a chemical: its string id,
elemental composition, charge, and molecular weight. It is the
single source of truth for these intrinsic properties across the
declaration surface: reaction stoichiometry, the speciation engine,
and any gas-liquid partitioning logic all reference the same
``Species`` object rather than redeclaring its atoms or charge.

`Species` is frozen and equality-comparable on its fields, so two
``Species(id="CO2", atoms={"C":1,"O":2}, charge=0)`` constructions
in different files compare equal even when they are distinct
objects. The :func:`~PyOMES.chemistry.species_check.check_species_consistency`
utility detects when this happens (a *soft conflict*: redeclared
instead of imported) and warns.

Usage
-----
>>> from PyOMES.chemistry import Species
>>> CO2 = Species(id="CO2", atoms={"C": 1, "O": 2}, charge=0)
>>> CO2.MW   # auto-computed: 12.011 + 2*15.999
44.009
>>> CO2.atoms["C"]
1
>>> CO2 == Species(id="CO2", atoms={"C": 1, "O": 2}, charge=0)
True

For universal inorganic aqueous species, prefer importing from
:mod:`~PyOMES.chemistry.common_species`:

>>> from PyOMES.chemistry.common_species import CO2, H2O, NH3
"""

from __future__ import annotations

import types
from dataclasses import dataclass, field
from typing import Dict, Mapping, Optional

# IUPAC 2021 standard atomic weights (g/mol).
# Used to compute MW automatically when not supplied explicitly.
_ATOMIC_WEIGHTS: Dict[str, float] = {
    "H":  1.008,
    "C":  12.011,
    "N":  14.007,
    "O":  15.999,
    "P":  30.974,
    "S":  32.065,
    "Na": 22.990,
    "Mg": 24.305,
    "Cl": 35.45,
    "K":  39.098,
    "Ca": 40.078,
    "Mn": 54.938,
    "Fe": 55.845,
    "Co": 58.933,
    "Ni": 58.693,
    "Cu": 63.546,
    "Zn": 65.38,
    "Mo": 95.96,
}


class SpeciesConflictError(ValueError):
    """Raised when two ``Species`` declarations share an id but differ
    in atoms, charge, or MW.

    Attributes
    ----------
    species_id : str
        The conflicting id.
    details : str
        Per-field breakdown of the disagreement.
    """

    def __init__(self, message: str, *, species_id: str = "", details: str = ""):
        super().__init__(message)
        self.species_id = species_id
        self.details = details


@dataclass(frozen=True, eq=True)
class Species:
    """Frozen identity record for a chemical species.

    Parameters
    ----------
    id : str
        Unique identifier (e.g. ``"CO2"``, ``"H+"``, ``"S_ac"``).
        Used as the dict key in phase ``n_mol`` mappings and as the
        engine output key, so it must match the convention the
        speciation engine and downstream consumers expect.
    atoms : Mapping[str, int]
        Elemental composition (e.g. ``{"C": 1, "O": 2}`` for CO₂).
        Stored as an immutable mapping so the frozen dataclass
        remains hashable.
    charge : int
        Formal charge (e.g. ``+1`` for H⁺, ``-2`` for CO₃²⁻).
        Defaults to ``0`` for neutral species.
    MW : float or None
        Molecular weight (g/mol).  When omitted (the default), MW is
        computed automatically from ``atoms`` using the IUPAC 2021
        standard atomic weights in :data:`_ATOMIC_WEIGHTS`.  Provide
        an explicit value to override — useful when the formula unit
        does not match the true stoichiometry (e.g. empirical biomass
        formulas) or for virtual charge-carrier species with no atoms.
        Raises :class:`ValueError` if ``atoms`` contains an element
        not present in :data:`_ATOMIC_WEIGHTS` and no explicit ``MW``
        is given.
    """

    id: str
    atoms: Mapping[str, int] = field(default_factory=dict)
    charge: int = 0
    MW: Optional[float] = None

    def __post_init__(self):
        # Coerce atoms to an immutable MappingProxyType so the
        # dataclass stays hashable. Use object.__setattr__ because the
        # dataclass is frozen.
        if not isinstance(self.atoms, types.MappingProxyType):
            object.__setattr__(
                self, "atoms",
                types.MappingProxyType(dict(self.atoms)),
            )
        # Auto-compute MW from atoms when not supplied.
        if self.MW is None:
            unknown = set(self.atoms) - _ATOMIC_WEIGHTS.keys()
            if unknown:
                raise ValueError(
                    f"Species {self.id!r}: cannot compute MW — no atomic weight "
                    f"for element(s) {unknown}. Supply MW explicitly."
                )
            computed = sum(
                _ATOMIC_WEIGHTS[el] * float(n) for el, n in self.atoms.items()
            )
            object.__setattr__(self, "MW", computed)
        else:
            object.__setattr__(self, "MW", float(self.MW))

    def __hash__(self):
        # Hash on id + sorted atoms + charge + MW. MappingProxyType
        # isn't directly hashable, so build a tuple of items.
        return hash((
            self.id,
            tuple(sorted(self.atoms.items())),
            self.charge,
            self.MW,
        ))

    def __repr__(self):
        atoms_str = "".join(
            f"{el}{n}" if n != 1 else el
            for el, n in sorted(self.atoms.items())
        )
        if self.charge > 0:
            charge_str = f"^{self.charge}+" if self.charge != 1 else "^+"
        elif self.charge < 0:
            charge_str = f"^{abs(self.charge)}-" if self.charge != -1 else "^-"
        else:
            charge_str = ""
        return f"Species({self.id!r}: {atoms_str}{charge_str})"
