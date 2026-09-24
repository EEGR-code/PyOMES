# -*- coding: utf-8 -*-
"""DaviesLiquidModel: Davies-equation activity coefficients for the liquid phase."""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Dict

from .water_properties import (
    debye_huckel_A,
    ionic_strength_molal_from_molar,
    water_kg_per_L,
)


_LN10 = float(np.log(10.0))


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
        dIm_dImolL = 1.0 / water_kg_per_L(T_K)
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
