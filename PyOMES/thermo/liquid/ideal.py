# -*- coding: utf-8 -*-
"""IdealLiquidModel: ideal-solution liquid, γ_i = 1 for every species."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict


@dataclass(frozen=True)
class IdealLiquidModel:
    """Ideal-solution liquid model: γ_i = 1.0 for all species at all conditions.

    Returns an empty dict — the absent-means-one convention means callers
    obtain γ_i = 1.0 via ``result.get(species_id, 1.0)``.  Correct for
    dilute non-ionic systems (CH₄, H₂, O₂) at bioprocess concentrations.

    Also satisfies ``ActivityModel`` (per-ion protocol):
        gamma(z, I_molL, *, T_K) → 1.0
    """

    name: str = "ideal"

    def gamma_all(
        self,
        x_mol: Dict[str, float],
        T_K: float,
        *,
        charge: Dict[str, int],
    ) -> Dict[str, float]:
        return {}

    def gamma(self, z: float, I_molL: float, *, T_K: float) -> float:
        """ActivityModel compatibility — always returns 1.0."""
        return 1.0
