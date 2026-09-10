"""Typed inputs/outputs for aqueous equilibrium calculations.

Why this exists
---------------
Historically the codebase passed many chemistry-related values around as ad-hoc
keyword arguments (``CT_TIC``, ``CT_P``, ``acid_totals``...). That is convenient
inside an ODE loop, but makes it harder to:

1) reuse the speciation solver as a standalone tool
2) switch between process-focused vs equilibrium-focused solving policies

The dataclasses below provide a stable schema that matches the *current* engine
inputs, while still allowing future evolution toward a more "totals/alkalinity"
interface (e.g., TIC + TA + strong ions).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .registry import map_user_ions_to_engine


def _T_K_to_C(T_K: float) -> float:
    return float(T_K) - 273.15


@dataclass(frozen=True)
class AqueousTotals:
    """Conserved totals and weak-acid systems for aqueous speciation.

    This schema mirrors the current :class:`PyOMES.chemical_equilibrium.engine.BisectionChemicalEquilibriumEngine`
    call signature, so it can be adopted without changing chemistry logic.

    Parameters
    ----------
    T_C:
        Temperature in degrees C.
    CT_TIC:
        Total inorganic carbon in the *liquid phase* (mol/L). (In coupled VLE,
        this is the liquid share of a system inventory.)
    CT_P:
        Total phosphate in the liquid (mol/L as total P).
    CT_NH_T:
        Total ammonia in the liquid (mol/L as total N).
    acid_totals:
        Weak-acid "systems" (e.g. acetate/citrate) represented as total acid
        concentrations (mol/L). Keys must match the names understood by the
        speciation engine.
    acid_pKas:
        pKas for each acid system (list of floats).
    strong_ions:
        Strong-ion totals in mol/L, using the same naming convention as the
        engine kwargs (e.g., ``CT_Na``, ``CT_Cl``, ``CT_K``...).
    """

    # Thermodynamic conditions
    T_C: float = 25.0

    # Core conserved totals (liquid phase)
    CT_TIC: float = 0.0
    CT_P: float = 0.0
    CT_NH_T: float = 0.0

    # Weak acid systems (optional)
    acid_totals: Dict[str, float] = field(default_factory=dict)
    acid_pKas: Dict[str, List[float]] = field(default_factory=dict)

    # Strong ions (engine-style kwargs)
    strong_ions: Dict[str, float] = field(default_factory=dict)

    # Optional solver hints
    logH_guess: Optional[float] = None


@dataclass(frozen=True)
class AqueousTotalsUser:
    """User-facing totals schema for standalone aqueous equilibrium.

    This is a translation layer that makes it easier to use speciation without
    knowing the internal engine's keyword conventions.

    Notes
    -----
    - This class does **not** change chemistry logic. It maps to the existing
      :class:`AqueousTotals` (engine-style) inputs via :meth:`to_engine`.
    - At this stage, alkalinity (TA) is not part of the core engine schema.
      You may still use this class for systems whose charge balance is governed
      by the provided strong ions.
    """

    # Thermodynamic conditions
    T_K: float = 298.15

    # Conserved totals (mol/L)
    TIC_mol_L: float = 0.0
    P_tot_mol_L: float = 0.0
    TAN_mol_L: float = 0.0

    # Optional alkalinity (equivalents/L).
    # Not consumed by the current engine-style schema, but included now so the
    # user-facing API can grow toward a more conventional (TIC + TA + ions)
    # interface without breaking downstream code.
    TA_eq_L: Optional[float] = None

    # Strong ions using chemistry-style labels (e.g., 'Na+', 'Cl-')
    strong_ions_mol_L: Dict[str, float] = field(default_factory=dict)

    # Optional weak-acid systems (names must match what the current engine understands)
    organic_acids_mol_L: Dict[str, float] = field(default_factory=dict)
    organic_acid_pKas: Dict[str, List[float]] = field(default_factory=dict)

    # Optional solver hint
    logH_guess: Optional[float] = None

    def to_engine(self) -> AqueousTotals:
        """Translate to the engine-style :class:`AqueousTotals`."""
        strong = map_user_ions_to_engine(self.strong_ions_mol_L or {})
        return AqueousTotals(
            T_C=_T_K_to_C(self.T_K),
            CT_TIC=float(self.TIC_mol_L),
            CT_P=float(self.P_tot_mol_L),
            CT_NH_T=float(self.TAN_mol_L),
            acid_totals=dict(self.organic_acids_mol_L or {}),
            acid_pKas=dict(self.organic_acid_pKas or {}),
            strong_ions=strong,
            logH_guess=self.logH_guess,
        )


@dataclass(frozen=True)
class AqueousEquilibrium:
    """Results from an aqueous speciation equilibrium calculation."""

    totals: AqueousTotals

    # Primary outputs
    pH: float
    CO2: Optional[float] = None

    # Optional derived values
    ionic_strength_molL: Optional[float] = None

    # Raw engine output for advanced users / debugging
    raw: Dict[str, float] = field(default_factory=dict)

    # Diagnostics
    converged: bool = True
    n_iter: Optional[int] = None
    message: str = ""
