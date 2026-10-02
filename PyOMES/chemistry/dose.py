# -*- coding: utf-8 -*-
"""Doses: what one mole of an added reagent puts into a phase.

A dose is a composition ``{species_id: mol per mol of reagent}``, e.g.
``{"Na+": 1, "OH-": 1}`` for sodium hydroxide; a plain id is one mole of
that species (``"H3PO4"`` is ``{"H3PO4": 1}``). What is added is exactly
what the composition says. :func:`check_dose` checks it against the model's
species and warns when it is not charge-neutral, since a real reagent is.
"""

from __future__ import annotations

import warnings
from typing import Dict, Mapping, Union

from .species import Species

Dose = Union[str, Mapping[str, float]]


def as_dose(dose: Dose) -> Dict[str, float]:
    """Return *dose* as ``{species_id: mol per mol of reagent}``."""
    if isinstance(dose, str):
        if not dose:
            raise ValueError("A dose id must not be empty.")
        return {dose: 1.0}
    if not isinstance(dose, Mapping) or not dose:
        raise TypeError(
            f"A dose is a species id or a non-empty {{species_id: mol}} "
            f"mapping, got {dose!r}"
        )
    out = {str(sp): float(n) for sp, n in dose.items()}
    bad = {sp: n for sp, n in out.items() if n <= 0.0}
    if bad:
        raise ValueError(f"Dose amounts must be > 0, got {bad}")
    return out


def check_dose(
    dose: Dose,
    species: Mapping[str, Species],
    *,
    label: str = "dose",
) -> Dict[str, float]:
    """Return *dose* as a composition, checked against *species*.

    Raises
    ------
    ValueError
        If an id in the dose is not among *species*.

    Warns
    -----
    UserWarning
        If the dose's net charge per mole of reagent is not zero (e.g.
        ``{"Na+": 1}`` without its ``OH-``).
    """
    composition = as_dose(dose)
    missing = sorted(sp for sp in composition if sp not in species)
    if missing:
        available = ", ".join(sorted(species)) or "no species were passed"
        raise ValueError(
            f"{label} {composition} names {missing}, not among the species "
            f"passed to this model (available: {available}). Add a Species "
            f"for each, e.g. species=[Species(id={missing[0]!r}, ...)]."
        )
    net_charge = sum(n * species[sp].charge for sp, n in composition.items())
    if abs(net_charge) > 1e-12:
        warnings.warn(
            f"{label} {composition} adds a net charge of {net_charge:+g} per "
            f"mole of reagent. A real reagent is neutral: include its "
            f"counter-ion (e.g. {{'Na+': 1, 'OH-': 1}} for NaOH).",
            UserWarning,
            stacklevel=3,
        )
    return composition
