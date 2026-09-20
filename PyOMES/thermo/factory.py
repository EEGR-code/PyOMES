# -*- coding: utf-8 -*-
"""String-to-model factory for per-ion activity models."""

from __future__ import annotations

from .liquid_phase_model import ActivityModel, DaviesLiquidModel, IdealLiquidModel
from .sit_liquid_model import SITLiquidModel


def make_activity_model(use_activity: bool, activity_model: str) -> ActivityModel:
    """Factory for activity models used by speciation engines."""
    if not use_activity:
        return IdealLiquidModel()

    name = (activity_model or "").strip().lower()
    if name in ("davies", "daviesactivitymodel", "daviesliquidmodel"):
        return DaviesLiquidModel()

    if name in ("sit", "sitactivitymodel", "sitliquidmodel"):
        return SITLiquidModel()

    raise ValueError(
        f"Unknown activity_model={activity_model!r}. "
        "Supported: 'davies', 'sit' (or use_activity=False for ideal)."
    )
