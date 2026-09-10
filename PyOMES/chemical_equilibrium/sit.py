# -*- coding: utf-8 -*-
"""Specific Ion Interaction Theory (SIT) activity coefficient model.

The canonical implementation is now ``SITLiquidModel`` in PyOMES/thermo/.
This module re-exports it under the legacy ``SITActivityModel`` name for
backward compatibility with any callers that import from here.

    SITActivityModel  ≡  SITLiquidModel

See PyOMES/thermo/sit_liquid_model.py for the full implementation and docstring.
"""

from __future__ import annotations

# Re-export everything callers may have imported from this module.
from PyOMES.thermo.sit_liquid_model import (  # noqa: F401
    SIT_EPSILON,
    ION_CHARGES,
    _get_epsilon,
    SITLiquidModel,
)

# Legacy alias.
SITActivityModel = SITLiquidModel

__all__ = [
    "SIT_EPSILON",
    "ION_CHARGES",
    "_get_epsilon",
    "SITLiquidModel",
    "SITActivityModel",
]
