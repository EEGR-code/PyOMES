# -*- coding: utf-8 -*-
"""SystemEnv — read-only system state view for controller state_rates().

Phase B (CONTROLLER_STATE_PROTOCOL) declares this type. Phase C
(SYSTEM_SOLVER_PROTOCOL) imports it and adds helper methods such as
``liquid_pH()``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, Optional

if TYPE_CHECKING:
    from .control_volume import ControlVolume


@dataclass(frozen=True)
class SystemEnv:
    """Immutable snapshot of system state passed to ``controller.state_rates()``.

    Parameters
    ----------
    cvs : dict
        ``{cv_key: ControlVolume}`` — read-only view of all CVs at the
        current ODE sub-step. Do not mutate.
    t_h : float
        Current simulation time (hours) at the sub-step evaluation point.

    Notes
    -----
    Phase C (``SYSTEM_SOLVER_PROTOCOL``) extends this class with
    convenience accessors such as ``liquid_pH(cv_key)`` that read pH
    from a CV's liquid phase without boilerplate. Any such helpers are
    added to this file when Phase C lands; the dataclass fields above
    are the stable contract shared with Phase B.
    """

    cvs: Dict[str, Any]  # Dict[str, ControlVolume] at runtime
    t_h: float

    def liquid_pH(self, cv_key: str) -> Optional[float]:
        """Return the liquid-phase pH for *cv_key*, or ``None`` if unavailable.

        Returns ``None`` when the CV is absent, has no liquid phase, or the
        liquid phase has no ``pH`` attribute (i.e. the speciation engine has
        not yet run for this phase).
        """
        cv = self.cvs.get(cv_key)
        if cv is None:
            return None
        liq = getattr(cv, "phases", {}).get("liquid")
        if liq is None:
            return None
        pH = getattr(liq, "pH", None)
        return float(pH) if pH is not None else None
