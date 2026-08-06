"""Factory for creating speciation engines.

This factory is a thin convenience layer that selects a numerical *policy*:

- ``policy='process'``: warm-start-friendly defaults intended for ODE coupling.
- ``policy='equilibrium'``: stricter defaults intended for standalone use.
"""

from __future__ import annotations

from dataclasses import dataclass

from .engine import BisectionChemicalEquilibriumEngine
from .api import SpeciationEngineAdapter


@dataclass(frozen=True)
class SpeciationFactory:
    """Create a configured speciation engine adapter."""

    @staticmethod
    def create(
        *,
        policy: str = "process",
        use_activity: bool = True,
        activity_model: str = "davies",
        T_C: float = 25.0,
    ) -> SpeciationEngineAdapter:
        engine = BisectionChemicalEquilibriumEngine(
            use_activity=bool(use_activity),
            activity_model=str(activity_model) if activity_model is not None else None,
            T_C=float(T_C),
        )

        if policy not in ("process", "equilibrium"):
            raise ValueError("policy must be 'process' or 'equilibrium'")

        if policy == "process":
            return SpeciationEngineAdapter(
                engine=engine,
                policy=policy,
                default_tol=1e-10,
                default_n_scan=600,
            )

        return SpeciationEngineAdapter(
            engine=engine,
            policy=policy,
            default_tol=1e-12,
            default_n_scan=900,
        )
