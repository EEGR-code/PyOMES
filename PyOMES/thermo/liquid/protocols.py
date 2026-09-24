# -*- coding: utf-8 -*-
"""Protocols for liquid-phase non-ideality: LiquidPhaseModel, ActivityModel and
DifferentiableLiquidModel.

``LiquidPhaseModel`` is the liquid-side counterpart of
:class:`~PyOMES.thermo.gas.protocols.GasEOS`. The models that satisfy these
protocols live next to this module: ``ideal.py``, ``davies.py`` and ``sit.py``.

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
``IdealLiquidModel``, ``DaviesLiquidModel`` and ``SITLiquidModel`` each satisfy both:
- ``LiquidPhaseModel`` — ``gamma_all(x_mol, T_K, *, charge)``
- ``ActivityModel`` — ``gamma(z, I_molL, *, T_K)``

This allows ThermoFramework to hold a single ``liquid_activity`` object that
works for both the phase-level LiquidPhaseModel interface (used by
``HenryEquilibrium``) and the per-ion ActivityModel interface used by the NR and
Bisection chemical-equilibrium engines.
"""
from __future__ import annotations

import numpy as np
from typing import Dict, Protocol, runtime_checkable


@runtime_checkable
class DifferentiableLiquidModel(Protocol):
    """Extension of :class:`LiquidPhaseModel` exposing an analytic Jacobian.

    ``∂γ_i/∂C_j`` is the activity contribution to the Jacobian of an
    equilibrium residual; this protocol is the liquid-side parallel of
    :class:`~PyOMES.chemical_equilibrium.protocols.SplitJacobianCapable`.
    Not every :class:`LiquidPhaseModel` need satisfy it — only those with a
    cheap analytic derivative (Davies, SIT). A composition-based model such as
    NRTL would need the full ``∂γ_i/∂x_j`` matrix.

    Nothing in the package calls it today. The NR solver
    (``chemical_equilibrium/engines/nr/solver.py``) updates γ in an outer
    fixed-point loop on ionic strength and holds it fixed within each inner
    Newton solve, so it converges without differentiating through γ. The
    Jacobian is for callers that need ``∂g/∂z`` including activity
    sensitivity; the tests check it against finite differences.
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
                               strength internally
    - ``SITLiquidModel``     — Specific Ion Interaction; also ionic-strength
                               based

    A model that uses the full x_mol (for example NRTL) would satisfy the
    same protocol; none is implemented.
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
