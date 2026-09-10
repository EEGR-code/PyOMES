"""Public API wrapper for the speciation engine.

The existing :class:`PyOMES.chemical_equilibrium.engine.BisectionChemicalEquilibriumEngine` is optimized
for being called inside ODE loops via a keyword-heavy ``solve(...)`` method.

For the orchestrator refactor we want a *stable*, typed interface that can be
used:

1) inside fermenter simulations (process-focused, warm-started)
2) as a standalone aqueous equilibrium tool (equilibrium-focused)

This module provides an adapter that wraps the existing engine without changing
its numerical behavior.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

from ..chemistry.types import AqueousTotals, AqueousTotalsUser, AqueousEquilibrium
from .engine import BisectionChemicalEquilibriumEngine


@dataclass
class SpeciationEngineAdapter:
    """Adapter exposing a typed ``equilibrate`` API."""

    engine: BisectionChemicalEquilibriumEngine
    policy: str = "process"  # "process" or "equilibrium"

    # Default numerical settings by policy. These are intentionally conservative;
    # callers may override via `options=`.
    default_tol: float = 1e-12
    default_n_scan: int = 900
    default_pH_min: float = 0.0
    default_pH_max: float = 14.0

    def equilibrate(
        self,
        totals: AqueousTotals | AqueousTotalsUser,
        *,
        initial_guess: Optional[Dict[str, Any]] = None,
        options: Optional[Dict[str, Any]] = None,
    ) -> AqueousEquilibrium:
        """Compute aqueous equilibrium for the given totals.

        Notes
        -----
        - ``initial_guess`` is currently limited to ``{"logH_guess": ...}``.
        - For process mode, callers typically supply ``totals.logH_guess`` from
          the previous timestep.
        """
        # Accept either engine-style totals (AqueousTotals) or user-facing totals
        # (AqueousTotalsUser). The latter is translated without changing chemistry.
        if isinstance(totals, AqueousTotalsUser):
            totals_eng = totals.to_engine()
        else:
            totals_eng = totals

        opts = dict(options or {})
        logH_guess = totals_eng.logH_guess
        if initial_guess and ("logH_guess" in initial_guess):
            logH_guess = float(initial_guess["logH_guess"])

        # Flatten strong ions to match current engine kwargs
        strong_kwargs = {k: float(v) for k, v in (totals_eng.strong_ions or {}).items()}

        out = self.engine.solve(
            acid_totals=dict(totals_eng.acid_totals or {}),
            acid_pKas=dict(totals_eng.acid_pKas or {}),
            CT_TIC=float(totals_eng.CT_TIC),
            CT_P=float(totals_eng.CT_P),
            CT_NH_T=float(totals_eng.CT_NH_T),
            logH_guess=logH_guess,
            pH_min=float(opts.get("pH_min", self.default_pH_min)),
            pH_max=float(opts.get("pH_max", self.default_pH_max)),
            n_scan=int(opts.get("n_scan", self.default_n_scan)),
            tol=float(opts.get("tol", self.default_tol)),
            **strong_kwargs,
        )

        # "CO2" is in the canonical writeback tuple but the legacy kwargs
        # path (used here) emits "CO2aq" instead, which isn't — so it may
        # land in `extra` rather than `species_mol_L`; check both.
        CO2_val = out.species_mol_L.get("CO2", out.extra.get("CO2"))
        if CO2_val is None:
            CO2_val = out.species_mol_L.get("CO2aq", out.extra.get("CO2aq"))

        return AqueousEquilibrium(
            totals=totals_eng,
            pH=float(out.pH),
            CO2=float(CO2_val) if CO2_val is not None else None,
            ionic_strength_molL=float(out.ionic_strength) if out.ionic_strength is not None else None,
            raw=out.to_dict(),
            converged=True,
            n_iter=int(out.n_iter) if out.n_iter is not None else None,
            message="",
        )

    def get_CO2aq_from_totals(self, totals: AqueousTotals | AqueousTotalsUser, *, options: Optional[Dict[str, Any]] = None) -> float:
        """Convenience wrapper for Henry-law coupling.

        Mirrors :meth:`PyOMES.chemical_equilibrium.engine.BisectionChemicalEquilibriumEngine.get_CO2aq_from_totals`.
        """
        if isinstance(totals, AqueousTotalsUser):
            totals_eng = totals.to_engine()
        else:
            totals_eng = totals

        opts = dict(options or {})
        logH_guess = totals_eng.logH_guess
        if "logH_guess" in opts and opts["logH_guess"] is not None:
            logH_guess = float(opts["logH_guess"])
        strong_kwargs = {k: float(v) for k, v in (totals_eng.strong_ions or {}).items()}
        return float(
            self.engine.get_CO2aq_from_totals(
                acid_totals=dict(totals_eng.acid_totals or {}),
                acid_pKas=dict(totals_eng.acid_pKas or {}),
                CT_TIC=float(totals_eng.CT_TIC),
                CT_P=float(totals_eng.CT_P),
                CT_NH_T=float(totals_eng.CT_NH_T),
                logH_guess=logH_guess,
                pH_min=float(opts.get("pH_min", self.default_pH_min)),
                pH_max=float(opts.get("pH_max", self.default_pH_max)),
                n_scan=int(opts.get("n_scan", self.default_n_scan)),
                tol=float(opts.get("tol", self.default_tol)),
                **strong_kwargs,
            )
        )
