# -*- coding: utf-8 -*-
"""LiquidPhaseModel protocol and liquid-phase non-ideality implementations.

Liquid-side EOS, symmetric with GasEOS in src/equilibria/vle.py.

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
- ``ActivityModel`` (src/speciation/) — ``gamma(z, I_molL, *, T_K)``

This allows ThermoFramework to hold a single ``liquid_activity`` object that
works for both the phase-level LiquidPhaseModel and the per-ion ActivityModel
interfaces used internally by NRChemicalEquilibriumEngine.
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Dict, Optional, Protocol, runtime_checkable

from .water_properties import (
    debye_huckel_A,
    ionic_strength_molal_from_molar,
    water_density_kg_per_m3,
)


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

    This is standalone groundwork, not yet wired into ``nr_solver.py``'s
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


def _kg_per_L(T_K: float) -> float:
    """Water density in kg/L at T_K — the I_molL -> I_molal conversion factor."""
    rho_kg_m3 = water_density_kg_per_m3(float(T_K))
    if not np.isfinite(rho_kg_m3) or rho_kg_m3 <= 0.0:
        return 1.0
    return float(rho_kg_m3) / 1000.0


_LN10 = float(np.log(10.0))


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


@dataclass(frozen=True)
class DaviesLiquidModel:
    """Davies equation activity model satisfying both LiquidPhaseModel and ActivityModel.

    Implements the Davies equation (base-10):
        log₁₀(γ_i) = −A(T) z² ( √I/(1+√I) − 0.3 I )

    where I is molal ionic strength and A(T) is the temperature-dependent
    Debye–Hückel constant computed from water density and dielectric constant.

    Valid to I ≈ 0.5 mol/kg.  Use SITLiquidModel for higher ionic strength.

    Dual-protocol
    -------------
    - ``gamma_all(x_mol, T_K, *, charge)`` — computes I from composition,
      returns per-species γ for all ionic species (z ≠ 0 per ``charge`` dict).
      Non-ionic species are absent (caller gets γ = 1.0 via default).
    - ``gamma(z, I_molL, *, T_K)`` — per-ion interface (ActivityModel protocol);
      takes pre-computed ionic strength in mol/L, converts to molality internally.
    """

    name: str = "davies"

    def gamma_all(
        self,
        x_mol: Dict[str, float],
        T_K: float,
        *,
        charge: Dict[str, int],
    ) -> Dict[str, float]:
        # Compute ionic strength in mol/L from composition
        I_molL = 0.0
        for s, c in x_mol.items():
            z = charge.get(s, 0)
            if z != 0:
                I_molL += 0.5 * z * z * max(0.0, float(c))

        result = {}
        for s in x_mol:
            z = charge.get(s, 0)
            if z != 0:
                result[s] = self.gamma(float(z), I_molL, T_K=T_K)
        return result

    def gamma(self, z: float, I_molL: float, *, T_K: float) -> float:
        """Activity coefficient for a single ion (ActivityModel protocol).

        Parameters
        ----------
        z : float
            Ion charge number.
        I_molL : float
            Ionic strength in mol/L; converted to molality internally.
        T_K : float
            Temperature in Kelvin.
        """
        z = float(z)
        I_molL = float(I_molL)
        if z == 0.0 or not np.isfinite(I_molL) or I_molL <= 0.0:
            return 1.0

        I = ionic_strength_molal_from_molar(I_molL, T_K=float(T_K))
        if not np.isfinite(I) or I <= 0.0:
            return 1.0

        A = debye_huckel_A(float(T_K))
        sqrtI = float(np.sqrt(I))
        log10_gamma = -A * (z**2) * (sqrtI / (1.0 + sqrtI) - 0.3 * I)
        return float(10.0**log10_gamma)

    def jacobian_dgamma_dx(
        self,
        x_mol: Dict[str, float],
        T_K: float,
        *,
        charge: Dict[str, int],
    ) -> np.ndarray:
        """``∂γ_i/∂C_j = (∂γ_i/∂I)(z_j²/2)`` — §8.4 of THERMODYNAMIC_MODEL_ARCHITECTURE.md.

        Differentiates the Davies equation's own I-dependence analytically
        (``I = 0.5 Σ z² C`` is itself linear in composition, so
        ``∂I/∂C_j = z_j²/2`` exactly). See :class:`DifferentiableLiquidModel`
        for why this is standalone groundwork, not yet consumed by the
        inner NR loop.
        """
        species_ids = sorted(x_mol)
        n = len(species_ids)
        J = np.zeros((n, n))

        I_molL = 0.0
        for s, c in x_mol.items():
            z = charge.get(s, 0)
            if z != 0:
                I_molL += 0.5 * z * z * max(0.0, float(c))
        if I_molL <= 0.0:
            return J   # dilute limit: gamma=1 everywhere, all partials zero

        I_m = ionic_strength_molal_from_molar(I_molL, T_K=float(T_K))
        if not np.isfinite(I_m) or I_m <= 0.0:
            return J

        A = debye_huckel_A(float(T_K))
        sqrtIm = float(np.sqrt(I_m))
        dIm_dImolL = 1.0 / _kg_per_L(T_K)
        # d(log10 gamma)/dI_m, from log10(gamma) = -A z^2 (sqrt(Im)/(1+sqrt(Im)) - 0.3*Im)
        d_log10_dIm = 1.0 / (2.0 * sqrtIm * (1.0 + sqrtIm) ** 2) - 0.3

        for i, sp_i in enumerate(species_ids):
            z_i = charge.get(sp_i, 0)
            if z_i == 0:
                continue
            gamma_i = self.gamma(float(z_i), I_molL, T_K=T_K)
            dgamma_i_dImolL = gamma_i * _LN10 * (-A * (z_i ** 2) * d_log10_dIm) * dIm_dImolL
            for j, sp_j in enumerate(species_ids):
                z_j = charge.get(sp_j, 0)
                if z_j == 0:
                    continue
                J[i, j] = dgamma_i_dImolL * (z_j ** 2) / 2.0

        return J
