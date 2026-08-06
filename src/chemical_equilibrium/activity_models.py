# -*- coding: utf-8 -*-
"""
Activity model implementations.

Design:
- ActivityModel is a strategy object used by chemistry solvers.
- IdealActivityModel / IdealLiquidModel: gamma = 1 for all species.
- DaviesActivityModel / DaviesLiquidModel: Davies equation vs ionic strength.

The canonical implementations now live in src/thermo/ (LiquidPhaseModel
protocol); the names here are re-exported aliases for backward compatibility.

    DaviesActivityModel  ≡  DaviesLiquidModel   (satisfies both protocols)
    IdealActivityModel   ≡  IdealLiquidModel     (satisfies both protocols)

Callers within src/speciation/ should continue to import from here;
external callers are encouraged to use PyOMES.thermo directly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Protocol

# Water property helpers — moved to src/thermo/; re-exported for backward compat.
from PyOMES.thermo.water_properties import (  # noqa: F401
    water_dielectric_constant,
    water_density_kg_per_m3,
    debye_huckel_A,
    ionic_strength_molal_from_molar,
)

# LiquidPhaseModel implementations that also satisfy ActivityModel.
from PyOMES.thermo.liquid_phase_model import (
    IdealLiquidModel,
    DaviesLiquidModel,
)

# ── ActivityModel protocol (per-ion interface used by NRChemicalEquilibriumEngine) ──

class ActivityModel(Protocol):
    """Protocol for per-ion activity coefficient models (Davies, SIT, ideal)."""
    name: str

    def gamma(self, z: float, I_molL: float, *, T_K: float) -> float:
        """Return activity coefficient gamma for an ion with charge z."""


# Backward-compatible aliases — these names existed before the LiquidPhaseModel
# refactor and are still used throughout src/speciation/.
IdealActivityModel = IdealLiquidModel
DaviesActivityModel = DaviesLiquidModel


def make_activity_model(use_activity: bool, activity_model: str) -> ActivityModel:
    """Factory for activity models used by speciation engines."""
    if not use_activity:
        return IdealLiquidModel()

    name = (activity_model or "").strip().lower()
    if name in ("davies", "daviesactivitymodel", "daviesliquidmodel"):
        return DaviesLiquidModel()

    if name in ("sit", "sitactivitymodel", "sitliquidmodel"):
        from PyOMES.thermo.sit_liquid_model import SITLiquidModel
        return SITLiquidModel()

    raise ValueError(
        f"Unknown activity_model={activity_model!r}. "
        "Supported: 'davies', 'sit' (or use_activity=False for ideal)."
    )
