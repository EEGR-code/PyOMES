"""Index-based state mapping for plug-in kinetic models.

:class:`StateMapping` maps between a kinetic model's internal ODE
state vector (arbitrary species names and ordering) and the canonical
reactor quantities used by the fermenter (e.g. ``X_g_L``,
``SAc_g_L``).  This allows kinetic models to define their own state
layout without coupling to the fermenter's naming conventions.

The :meth:`StateMapping.auto` classmethod creates a mapping by
matching state IDs to a set of canonical names.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Sequence, Tuple
import numpy as np

DEFAULT_CANONICAL_IDS: Tuple[str, ...] = ("X_g_L", "SAc_g_L", "SPr_g_L", "SBu_g_L")

@dataclass(frozen=True)
class StateMapping:
    """Index-based mapping between a kinetic model's ODE state vector and canonical reactor quantities."""

    state_ids: Tuple[str, ...]
    canonical_to_state_id: Dict[str, str]

    @classmethod
    def auto(cls, state_ids: Sequence[str], canonical_ids: Sequence[str] = DEFAULT_CANONICAL_IDS) -> "StateMapping":
        state_ids = tuple(state_ids)
        state_set = set(state_ids)
        mapping: Dict[str, str] = {}
        for canon in canonical_ids:
            if canon in state_set:
                mapping[canon] = canon
        return cls(state_ids=state_ids, canonical_to_state_id=mapping)

    def index(self) -> Dict[str, int]:
        return {sid: i for i, sid in enumerate(self.state_ids)}

    def canonical_index(self) -> Dict[str, int]:
        idx = self.index()
        out: Dict[str, int] = {}
        for canon, sid in self.canonical_to_state_id.items():
            if sid in idx:
                out[canon] = idx[sid]
        return out

    def make_y0(self, canonical_init: Dict[str, float], fill: float = 0.0) -> np.ndarray:
        y0 = np.full(len(self.state_ids), float(fill), dtype=float)
        cidx = self.canonical_index()
        for canon, val in canonical_init.items():
            j = cidx.get(canon, None)
            if j is not None:
                y0[j] = float(val)
        return y0

    def get_value(self, y: np.ndarray, canonical_id: str, default: float = 0.0) -> float:
        cidx = self.canonical_index()
        j = cidx.get(canonical_id, None)
        if j is None:
            return float(default)
        return float(y[j])

    def get_series(self, y_mat: np.ndarray, canonical_id: str, default: float = 0.0) -> np.ndarray:
        cidx = self.canonical_index()
        j = cidx.get(canonical_id, None)
        if j is None:
            return np.full(y_mat.shape[1], float(default), dtype=float)
        return np.asarray(y_mat[j, :], dtype=float)
