# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 17:11:11 2026

@author: k2473520
"""

# v8_split_8/kinetics/core.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Mapping, Optional, Sequence, Tuple, Any, Protocol
import numpy as np


@dataclass(frozen=True)
class ParameterSpec:
    """Schema entry for a kinetic model parameter (name, default, bounds, units)."""
    name: str
    default: float
    lower: Optional[float] = None
    upper: Optional[float] = None
    units: str = ""


class ParameterSet(dict):
    """
    Validated parameter dict. Keeps functionality simple: it's just a dict,
    but validates against a schema once at construction.
    """
    @classmethod
    def from_schema(cls, schema: Sequence[ParameterSpec], overrides: Optional[Mapping[str, float]] = None) -> "ParameterSet":
        data = {p.name: float(p.default) for p in schema}
        if overrides:
            for k, v in overrides.items():
                if k not in data:
                    raise KeyError(f"Unknown parameter '{k}'. Allowed: {sorted(data)}")
                data[k] = float(v)

        # bounds check (optional)
        for p in schema:
            val = data[p.name]
            if p.lower is not None and val < p.lower:
                raise ValueError(f"Parameter '{p.name}'={val} < lower bound {p.lower}")
            if p.upper is not None and val > p.upper:
                raise ValueError(f"Parameter '{p.name}'={val} > upper bound {p.upper}")

        return cls(data)


@dataclass(frozen=True)
class Environment:
    """
    Immutable snapshot of reactor conditions passed into RHS.
    Extend freely; models can ignore fields they don't need.
    """
    T_K: float
    V_L: float
    pH: float
    I_molal: float
    DO_mol_L: float
    CO2aq_mol_L: float
    P_atm: float

    # Optional: if your kinetics need these
    kla_O2_per_h: Optional[float] = None
    kla_CO2_per_h: Optional[float] = None
    D_h: Optional[float] = None  # dilution rate etc.
    DO_sat_mol_L: Optional[float] = None
    pO2_atm: Optional[float] = None


class KineticModel(Protocol):
    """
    Contract for any kinetic model. The fermenter only relies on this.
    """
    name: str
    state_ids: Tuple[str, ...]
    param_schema: Tuple[ParameterSpec, ...]

    def rhs(
        self,
        t_h: float,
        y: np.ndarray,
        env: Environment,
        params: Mapping[str, float],
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Returns:
            dydt: np.ndarray (same shape as y)
            outputs: dict of diagnostics (can be empty)
        """
        ...
