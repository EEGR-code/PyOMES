# -*- coding: utf-8 -*-
"""Activity-related utilities (ionic strength + warnings).

Activity-coefficient models live in :mod:`PyOMES.thermo`. Current fermenter behavior
uses ideal-solution speciation; the warning function helps users spot when activity
effects may become important.
"""

from __future__ import annotations

import numpy as np
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


def warn_if_high_ionic_strength(
    I_molL: float,
    *,
    threshold_molL: float = 0.1,
    t_hr: float | None = None,
    pH: float | None = None,
) -> bool:
    """Emit a one-time warning if ionic strength is high enough that activity effects may matter."""
    try:
        I = float(I_molL)
    except (TypeError, ValueError, AttributeError, KeyError):
        return False
    if not np.isfinite(I):
        return False
    if I < float(threshold_molL):
        return False

    ctx = []
    if t_hr is not None and np.isfinite(t_hr):
        ctx.append(f"t={float(t_hr):.3g} h")
    if pH is not None and np.isfinite(pH):
        ctx.append(f"pH={float(pH):.3g}")
    ctx_s = (", " + ", ".join(ctx)) if ctx else ""

    warnings.warn(
        f"[activity] Ionic strength is {I:.3g} mol/L{ctx_s}. "
        f"At this ionic strength, solution activity effects (activity coefficients / pKa shifts) "
        f"may become important; results from ideal-solution speciation may be biased.",
        RuntimeWarning,
        stacklevel=2,
    )
    return True
