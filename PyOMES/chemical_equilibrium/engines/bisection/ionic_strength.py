# -*- coding: utf-8 -*-
"""Ionic strength of a Bisection-engine speciation dict.

The Bisection engine returns concentrations keyed by species id, and some of
those ids are synthesised rather than declared (``{name}_HA``, ``{name}_A-``,
``{name}_H2A-``). The charge of each key is therefore read from its trailing
``+`` / ``-`` tokens, with a small override table for ids that carry no such
token. The convention holds for the ids that engine emits; an id that writes its
charge as a number (``Fe2+`` for Fe²⁺) would be read as +1. The NR engine takes
charges from declared ``Species.charge`` instead, and PHREEQC reports its own
ionic strength.

Activity-coefficient models live in :mod:`PyOMES.thermo`.
"""

from __future__ import annotations

import warnings
from typing import Any, Dict


# Ions whose ids do not carry a trailing +/- charge token and so
# cannot be parsed by :func:`_charge_from_suffix`. Add entries here
# for ions whose canonical id breaks the suffix convention. BSM2's
# ``Cation(inert)`` / ``Anion(inert)`` are the only ones in the
# current codebase.
_CHARGE_OVERRIDES: Dict[str, int] = {
    "Cation(inert)": +1,
    "Anion(inert)": -1,
}


def _charge_from_suffix(key: str) -> int:
    """Infer ion charge from trailing +/- tokens on a species id.

    Returns positive count of ``'+'`` characters if ``key`` ends in
    ``'+'``, negative count of ``'-'`` characters if it ends in
    ``'-'``, and 0 otherwise. Designed to handle both canonical
    (``"HCO3-"``, ``"CO3--"``, ``"Mg++"``, ``"PO4---"``,
    ``"Mo7O24------"``) and generic
    (``"S_ac_A-"``, ``"AceticAcid_A-"``) species names with one
    rule.

    Examples
    --------
    >>> _charge_from_suffix("H+")
    1
    >>> _charge_from_suffix("CO3--")
    -2
    >>> _charge_from_suffix("PO4---")
    -3
    >>> _charge_from_suffix("S_ac_A-")
    -1
    >>> _charge_from_suffix("CO2aq")
    0
    """
    if key.endswith("+"):
        return len(key) - len(key.rstrip("+"))
    if key.endswith("-"):
        return -(len(key) - len(key.rstrip("-")))
    return 0


def ionic_strength_from_speciation(sp: Dict[str, Any]) -> float:
    """Compute ionic strength (mol/L) from a speciation dict.

    Charge for each species is inferred from the trailing ``+``/``-``
    tokens on its key (``HCO3-`` → -1, ``CO3--`` → -2, ``Mg++`` →
    +2, ``PO4---`` → -3). Ions whose ids don't follow the convention
    are listed in :data:`_CHARGE_OVERRIDES`. Neutral species (keys
    not ending in a charge token and not in the override table)
    contribute 0.

    Charge inference follows a single rule: adding a new ion that
    follows the suffix convention needs no edit here — its charge is
    inferred automatically. Adding an ion whose id doesn't follow the
    convention requires one entry in :data:`_CHARGE_OVERRIDES`.
    """
    I_sum = 0.0
    for key, v in sp.items():
        if not isinstance(key, str):
            continue
        z = _CHARGE_OVERRIDES.get(key)
        if z is None:
            z = _charge_from_suffix(key)
        if z == 0:
            continue
        try:
            I_sum += float(v) * z * z
        except (TypeError, ValueError) as e:
            warnings.warn(
                f"Ionic strength: could not process {key}={v!r}: {e}. "
                f"Skipping — ionic strength may be underestimated.",
                RuntimeWarning, stacklevel=2,
            )
    return 0.5 * float(I_sum)
