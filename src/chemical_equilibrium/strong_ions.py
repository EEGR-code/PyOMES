# -*- coding: utf-8 -*-
"""Strong-ion inference utilities.

Supports both **BioSTEAM streams** and the standalone :class:`FeedState`.

BioSTEAM streams typically contain neutral salts (e.g., NaCl) rather than explicit ionic
species (Na+, Cl-). This module provides helpers to infer effective strong-ion totals in
mol/L assuming full dissociation into principal ions.
"""

from __future__ import annotations

import numpy as np
from typing import Dict

from PyOMES.chemistry.registry import SALT_DISSOCIATION_MAP


def _mol_L_from_object(feed, ID: str) -> float:
    """Return mol/L of *ID* from either a BioSTEAM stream or a FeedState.

    Handles the duck-typing required to support both object types.
    """
    # FeedState path (has .mol_L method)
    if hasattr(feed, "mol_L") and callable(getattr(feed, "mol_L", None)):
        try:
            return float(feed.mol_L(ID))
        except (TypeError, ValueError, AttributeError, KeyError):
            return 0.0

    # Plain-dict path (concentrations_mol_L)
    if isinstance(feed, dict):
        return float(feed.get(ID, 0.0))

    # BioSTEAM Stream path (imol in kmol/hr, F_vol in m3/hr -> kmol/m3 = mol/L)
    try:
        return float(feed.imol[ID] / feed.F_vol)
    except (TypeError, ValueError, AttributeError, KeyError):
        return 0.0


def strong_ions_from_feed_molL(feed) -> Dict[str, float]:
    """Infer strong-ion totals (mol/L) from a feed object.

    Accepts:
      - A :class:`~fermenter.stream_adapter.FeedState`
      - A BioSTEAM/thermosteam ``Stream``
      - A plain ``dict`` mapping species-ID -> mol/L

    Returns a dict with keys matching speciation expectations:
    CT_K, CT_Na, CT_Cl, CT_NO3, CT_SO4, CT_Mg, CT_Ca, CT_Zn, CT_Mn, CT_Co, CT_Mo7O24

    Notes
    -----
    - Assumes full dissociation of salts into principal ions.
    - If explicit ionic species exist in the feed, they are included too.
    - Extend ``SALT_DISSOCIATION_MAP`` in the registry to cover your naming scheme.
    """
    def _mol_L(ID: str) -> float:
        return _mol_L_from_object(feed, ID)

    # Start with any explicit ions if present
    ions = {
        "CT_K": _mol_L("K+") + _mol_L("K"),
        "CT_Na": _mol_L("Na+") + _mol_L("Na"),
        "CT_Cl": _mol_L("Cl-") + _mol_L("Cl"),
        "CT_NO3": _mol_L("NO3-") + _mol_L("NO3"),
        "CT_SO4": _mol_L("SO4--") + _mol_L("SO4"),
        "CT_Mg": _mol_L("Mg++") + _mol_L("Mg"),
        "CT_Ca": _mol_L("Ca++") + _mol_L("Ca"),
        "CT_Zn": _mol_L("Zn++") + _mol_L("Zn"),
        "CT_Mn": _mol_L("Mn++") + _mol_L("Mn"),
        "CT_Co": _mol_L("Co++") + _mol_L("Co"),
        "CT_Mo7O24": _mol_L("Mo7O24------") + _mol_L("Mo7O24"),
    }

    # Add contributions from salts if present in feed
    for salt_id, contribs in SALT_DISSOCIATION_MAP.items():
        c = _mol_L(salt_id)
        if c and np.isfinite(c):
            for key, nu in contribs:
                ions[key] = float(ions.get(key, 0.0)) + float(nu) * float(c)

    return ions
