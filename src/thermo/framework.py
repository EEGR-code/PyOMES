# -*- coding: utf-8 -*-
"""ThermoFramework: frozen thermodynamic convention holder for a CV.

Holds the liquid-phase activity model and gas EOS as protocol objects,
plus temperature-correction helpers.  Frozen by design — a thermodynamic
framework does not change mid-simulation.  Runtime overrides use
``dataclasses.replace()``:

    from dataclasses import replace
    from PyOMES.thermo import ThermoFramework, DaviesLiquidModel
    alt = replace(default_thermo, liquid_activity=DaviesLiquidModel())
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

from ..units import R_J_PER_MOL_K as _R_J
from .liquid_phase_model import LiquidPhaseModel, IdealLiquidModel, DaviesLiquidModel

if TYPE_CHECKING:
    from PyOMES.equilibria.vle import GasEOS


@dataclass(frozen=True)
class ThermoFramework:
    """Frozen holder for thermodynamic conventions within a CV.

    Parameters
    ----------
    liquid_activity : LiquidPhaseModel
        Liquid-phase non-ideality model.  Default: ``IdealLiquidModel()``
        (γ_i = 1 for all species).  Use ``DaviesLiquidModel()`` for
        Davies equation or ``SITLiquidModel()`` for SIT corrections.
    gas_eos : GasEOS or None
        Gas-phase equation of state.  Default ``None`` (no gas phase or
        uses ideal-gas via the KineticGasLiquidLink default).
    standard_T_K : float
        Standard-state temperature (K).  Default 298.15 K (25 °C).
    standard_P_atm : float
        Standard-state pressure (atm).  Default 1.0.

    Backward-compatible properties
    --------------------------------
    ``use_activity`` and ``activity_model`` are read-only properties for
    callers that still reference the old string-based API.  They are
    derived from ``liquid_activity`` and cannot be set.
    """

    liquid_activity: LiquidPhaseModel = field(default_factory=IdealLiquidModel)
    gas_eos: Optional[object] = None  # GasEOS at runtime; annotated via TYPE_CHECKING
    standard_T_K: float = 298.15
    standard_P_atm: float = 1.0

    # ── Backward-compatible read-only properties ───────────────────────

    @property
    def use_activity(self) -> bool:
        """True when ``liquid_activity`` applies non-trivial corrections."""
        return not isinstance(self.liquid_activity, IdealLiquidModel)

    @property
    def activity_model(self) -> str:
        """Name of the active liquid model (e.g. ``"davies"``, ``"sit"``)."""
        return self.liquid_activity.name

    # ── Temperature correction helpers ─────────────────────────────────

    def pKa_at_T(
        self,
        pKa_ref: float,
        dH_J_per_mol: float,
        T_K: float,
        T_ref_K: float = 298.15,
    ) -> float:
        """Compute pKa at temperature T_K via Van 't Hoff.

        Parameters
        ----------
        pKa_ref : float
            Reference pKa at ``T_ref_K``.
        dH_J_per_mol : float
            Enthalpy of dissociation (J/mol).  Zero → no correction.
        T_K : float
            Target temperature (K).
        T_ref_K : float
            Reference temperature (K).  Default 298.15 K.
        """
        if abs(dH_J_per_mol) < 1e-12 or abs(T_K - T_ref_K) < 0.01:
            return pKa_ref
        Ka_ref = 10.0 ** (-pKa_ref)
        Ka_T = Ka_ref * math.exp(
            -dH_J_per_mol / _R_J * (1.0 / T_K - 1.0 / T_ref_K)
        )
        return -math.log10(max(Ka_T, 1e-30))

    def Kw_at_T(self, T_K: float, T_ref_K: float = 298.15) -> float:
        """Compute water autoionisation constant at T_K.

        Uses the BSM2-canonical Van 't Hoff enthalpy (55 900 J/mol).
        """
        _DH_W = 55900.0
        Kw_ref = 1e-14
        if abs(T_K - T_ref_K) < 0.01:
            return Kw_ref
        return Kw_ref * math.exp(
            -_DH_W / _R_J * (1.0 / T_K - 1.0 / T_ref_K)
        )


# Pre-built instances — import these rather than constructing ThermoFramework
# with string kwargs.

# Ideal: no activity corrections, no gas EOS.
THERMO_IDEAL = ThermoFramework()

# Davies: ionic-strength-based activity corrections; valid to I ≈ 0.5 mol/kg.
THERMO_DAVIES = ThermoFramework(liquid_activity=DaviesLiquidModel())
