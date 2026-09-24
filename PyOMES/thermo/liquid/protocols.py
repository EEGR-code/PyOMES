# -*- coding: utf-8 -*-
"""LiquidPhaseModel protocol and liquid-phase non-ideality implementations.

Liquid-side EOS, symmetric with GasEOS in PyOMES/thermo/gas_eos.py.

Unit convention
---------------
``x_mol`` arguments contain concentrations in **mol/L** (PyOMES's native unit).
Implementations are responsible for any internal unit conversion (e.g.
DaviesLiquidModel converts to molality before applying Debye-Hückel).
Callers always pass mol/L and do not convert.

Species absent from a ``gamma_all`` return dict are treated as γ_i = 1.0
by callers:

    gamma = model.gamma_all(x_mol, T_K, charge=charges).get(species_id, 1.0)

Dual-protocol implementations
------------------------------
``DaviesLiquidModel`` satisfies both:
- ``LiquidPhaseModel`` — ``gamma_all(x_mol, T_K, *, charge)``
- ``ActivityModel`` (this module) — ``gamma(z, I_molL, *, T_K)``

This allows ThermoFramework to hold a single ``liquid_activity`` object that
works for both the phase-level LiquidPhaseModel and the per-ion ActivityModel
interfaces used internally by NRChemicalEquilibriumEngine.
"""
from __future__ import annotations

import numpy as np
from typing import Dict, Protocol, runtime_checkable


@runtime_checkable
class DifferentiableLiquidModel(Protocol):
    """Extension of :class:`LiquidPhaseModel` exposing an analytic Jacobian.

    Added by CP2 of ``LAYER1_GAP_CLOSURE`` per §8.4 of
    ``THERMODYNAMIC_MODEL_ARCHITECTURE.md``: the ``∂γ_i/∂C_j`` contribution
    is needed once gas-liquid constraints are folded into the NR tableau
    (this phase). Parallel to ``SplitJacobianCapable`` on the engine side.
    Not every :class:`LiquidPhaseModel` need satisfy this — only those with
    a cheap analytic derivative (Davies, SIT). NRTL would need the full
    ``∂γ_i/∂x_j`` matrix, deferred.

    This is standalone groundwork, not yet wired into ``engines/nr/solver.py``'s
    inner Newton loop: that loop already achieves correct convergence via
    the existing outer (ionic-strength fixed-point) / inner (NR) split,
    which treats γ as frozen within each inner solve rather than
    differentiating through it. This Jacobian is exposed for future
    white-box/DAE consumers (§10.4 of ``MASS_EXCHANGE_ARCHITECTURE.md``)
    that need ``∂g/∂z`` including activity sensitivity.
    """

    def jacobian_dgamma_dx(
        self,
        x_mol: Dict[str, float],
        T_K: float,
        *,
        charge: Dict[str, int],
    ) -> np.ndarray:
        """Return ``∂γ_i/∂C_j`` as a dense matrix.

        Parameters
        ----------
        x_mol, T_K, charge : see :meth:`LiquidPhaseModel.gamma_all`.

        Returns
        -------
        ndarray, shape (n, n)
            Row/column order is ``sorted(x_mol)`` (alphabetical by species
            id), matching the convention used by
            :class:`~PyOMES.chemical_equilibrium.protocols.SpeciationJacobian`
            elsewhere in this codebase. Row *i*, column *j* is
            ``∂γ_i/∂C_j`` evaluated at the given composition.
        """
        ...


@runtime_checkable
class LiquidPhaseModel(Protocol):
    """Liquid-phase non-ideality model (liquid-side EOS).

    Symmetric with ``GasEOS``: takes full phase composition in mol/L and
    returns per-species activity coefficients.  Both protocols are
    non-ideality corrections on opposite sides of the fugacity equality:

        φ_i(T,P,y) × y_i × P  =  γ_i(T,x) × x_i × f_i^ref
              ↑ GasEOS                  ↑ LiquidPhaseModel

    Implementations
    ---------------
    - ``IdealLiquidModel``   — γ_i = 1 for all species (default)
    - ``DaviesLiquidModel``  — Davies equation; compresses x_mol to ionic
                               strength internally (in CP2)
    - ``SITLiquidModel``     — Specific Ion Interaction; also ionic-strength
                               based (in CP2)
    - ``NRTLLiquidModel``    — Non-Random Two Liquid; uses full x_mol (future)
    """

    name: str

    def gamma_all(
        self,
        x_mol: Dict[str, float],
        T_K: float,
        *,
        charge: Dict[str, int],
    ) -> Dict[str, float]:
        """Activity coefficients γ_i for all species at the given composition.

        Parameters
        ----------
        x_mol : dict
            Species concentrations in mol/L (PyOMES native unit).
        T_K : float
            Temperature in Kelvin.
        charge : dict
            Charge number z_i per species ID.  Non-ionic models may ignore
            this kwarg without raising.

        Returns
        -------
        dict
            Mapping ``species_id → γ_i``.  Species absent from the returned
            dict are treated as γ_i = 1.0 by callers.
        """
        ...


class ActivityModel(Protocol):
    """Protocol for per-ion activity coefficient models (Davies, SIT, ideal)."""
    name: str

    def gamma(self, z: float, I_molL: float, *, T_K: float) -> float:
        """Return activity coefficient gamma for an ion with charge z."""
