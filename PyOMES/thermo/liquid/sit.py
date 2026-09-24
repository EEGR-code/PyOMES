# -*- coding: utf-8 -*-
"""SITLiquidModel — Specific Ion Interaction Theory activity model.

Implements both LiquidPhaseModel and ActivityModel protocols.

SIT equation (Brønsted–Guggenheim–Scatchard):
    log₁₀(γ_j) = −z_j² × A(T) × √I / (1 + 1.5√I) + Σ_k ε(j,k) × m_k

where k sums over counter-ions and ε(j,k) are ion-pair interaction
coefficients (kg/mol) from NEA-TDB compilations.

Valid to I ≈ 3 mol/kg; use DaviesLiquidModel for I < 0.5 mol/kg.

References
----------
- Grenthe et al., "Chemical Thermodynamics of Uranium", NEA-TDB, OECD (1992)
- Ciavatta, Ann. Chim. (Rome) 70, 551-567 (1980)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Tuple

import numpy as np

from .water_properties import (
    debye_huckel_A,
    ionic_strength_molal_from_molar,
    water_density_kg_per_m3,
)

_LN10 = math.log(10.0)


def _kg_per_L(T_K: float) -> float:
    """Water density in kg/L at T_K — the I_molL -> I_molal conversion factor."""
    rho_kg_m3 = water_density_kg_per_m3(float(T_K))
    if not np.isfinite(rho_kg_m3) or rho_kg_m3 <= 0.0:
        return 1.0
    return float(rho_kg_m3) / 1000.0


# ── SIT interaction coefficient database (ε values) ────────────────────────
# Units: kg/mol.  Format: {("cation", "anion"): ε}.  Symmetric: stored once.
# Temperature dependence of ε neglected (NEA-TDB standard practice;
# Δε/ΔT < 0.001 kg/mol/K).

SIT_EPSILON: Dict[Tuple[str, str], float] = {
    # Na⁺ pairs (most important for NaOH/NaHCO₃ dosed AD)
    ("Na+", "Cl-"):          0.03,
    ("Na+", "OH-"):          0.04,
    ("Na+", "HCO3-"):        0.00,
    ("Na+", "CO3--"):       -0.08,
    ("Na+", "SO4--"):       -0.12,
    ("Na+", "HS-"):          0.03,
    ("Na+", "H2PO4-"):      -0.08,
    ("Na+", "HPO4--"):      -0.15,
    ("Na+", "Acetate-"):     0.08,
    ("Na+", "Propionate-"):  0.10,
    ("Na+", "Butyrate-"):    0.12,
    ("Na+", "Valerate-"):    0.13,
    ("Na+", "NO3-"):        -0.04,
    # K⁺ pairs
    ("K+", "Cl-"):           0.00,
    ("K+", "OH-"):           0.09,
    ("K+", "HCO3-"):        -0.01,
    ("K+", "CO3--"):         0.02,
    ("K+", "SO4--"):        -0.06,
    # NH₄⁺ pairs (important for high-TAN AD)
    ("NH4+", "Cl-"):        -0.01,
    ("NH4+", "SO4--"):      -0.12,
    ("NH4+", "HCO3-"):       0.00,
    ("NH4+", "OH-"):         0.00,
    ("NH4+", "Acetate-"):    0.00,
    # H⁺ pairs
    ("H+", "Cl-"):           0.12,
    ("H+", "SO4--"):         0.02,
    ("H+", "NO3-"):          0.07,
    ("H+", "Acetate-"):      0.00,
    # Ca²⁺ pairs
    ("Ca++", "Cl-"):         0.14,
    ("Ca++", "OH-"):        -0.15,
    ("Ca++", "HCO3-"):       0.00,
    ("Ca++", "SO4--"):      -0.08,
    # Mg²⁺ pairs
    ("Mg++", "Cl-"):         0.19,
    ("Mg++", "OH-"):        -0.10,
    ("Mg++", "SO4--"):      -0.10,
}

# Charge lookup for backward-compat compute_gammas() method.
ION_CHARGES: Dict[str, int] = {
    "H+": +1, "Na+": +1, "K+": +1, "NH4+": +1,
    "Ca++": +2, "Mg++": +2, "Zn++": +2, "Mn++": +2, "Co++": +2,
    "OH-": -1, "Cl-": -1, "NO3-": -1, "HS-": -1,
    "HCO3-": -1, "Acetate-": -1, "Propionate-": -1,
    "Butyrate-": -1, "Valerate-": -1, "H2PO4-": -1, "HSO4-": -1,
    "CO3--": -2, "SO4--": -2, "HPO4--": -2,
    "PO4---": -3,
}


def _get_epsilon(
    eps_dict: Dict[Tuple[str, str], float],
    ion1: str,
    ion2: str,
) -> float:
    """Look up ε(ion1, ion2).  Symmetric: tries both orderings."""
    return eps_dict.get((ion1, ion2), eps_dict.get((ion2, ion1), 0.0))


@dataclass
class SITLiquidModel:
    """SIT activity model satisfying both LiquidPhaseModel and ActivityModel.

    Dual-protocol
    -------------
    - ``gamma_all(x_mol, T_K, *, charge)`` — full SIT with composition;
      uses the ``charge`` dict to identify cations/anions and computes the
      ion-pair interaction sum per ionic species.  Non-ionic species absent.
    - ``gamma(z, I_molL, *, T_K)`` — extended Debye-Hückel term only
      (charge-based; no ion-pair specificity).  ActivityModel compatibility.
    - ``compute_gammas(composition_molL, I_molL, T_K)`` — full SIT over
      ``ION_CHARGES`` ions; backward-compatible with NRChemicalEquilibriumEngine.

    Parameters
    ----------
    epsilon : dict, optional
        Ion-pair interaction coefficients ``{("cation", "anion"): ε}``.
        Default: :data:`SIT_EPSILON` (NEA-TDB values for AD-relevant ions).
    """

    name: str = "sit"
    epsilon: Dict[Tuple[str, str], float] = field(
        default_factory=lambda: dict(SIT_EPSILON)
    )

    def gamma_all(
        self,
        x_mol: Dict[str, float],
        T_K: float,
        *,
        charge: Dict[str, int],
    ) -> Dict[str, float]:
        """Activity coefficients for all ionic species using full SIT.

        Uses the caller-supplied ``charge`` dict (from the chemistry layer)
        rather than the hard-coded ``ION_CHARGES`` table, so custom species
        are handled automatically provided their charges are supplied.
        """
        # Compute ionic strength in mol/L
        I_molL = 0.0
        for s, c in x_mol.items():
            z = charge.get(s, 0)
            if z != 0:
                I_molL += 0.5 * z * z * max(0.0, float(c))

        I_m = ionic_strength_molal_from_molar(I_molL, T_K=T_K)
        if not (I_m > 0):
            return {}

        A = debye_huckel_A(T_K)
        sqrtI = math.sqrt(I_m)
        dh_factor = -A * sqrtI / (1.0 + 1.5 * sqrtI)

        rho_kg_L = water_density_kg_per_m3(T_K) / 1000.0
        comp_m = {
            s: max(0.0, float(c)) / rho_kg_L
            for s, c in x_mol.items()
            if float(c) > 0
        }

        cations = {s: m for s, m in comp_m.items() if charge.get(s, 0) > 0}
        anions  = {s: m for s, m in comp_m.items() if charge.get(s, 0) < 0}

        result = {}
        for s in x_mol:
            z = charge.get(s, 0)
            if z == 0:
                continue
            log10_gamma = dh_factor * z * z
            if z > 0:
                for anion, m_k in anions.items():
                    log10_gamma += _get_epsilon(self.epsilon, s, anion) * m_k
            else:
                for cation, m_k in cations.items():
                    log10_gamma += _get_epsilon(self.epsilon, cation, s) * m_k
            result[s] = 10.0**log10_gamma
        return result

    def gamma(self, z: float, I_molL: float, *, T_K: float) -> float:
        """Extended Debye-Hückel term only (ActivityModel protocol).

        Charge-based — no ion-pair specificity.  Use gamma_all() or
        compute_gammas() for full SIT corrections.
        """
        z = float(z)
        I_molL = float(I_molL)
        if z == 0.0 or not np.isfinite(I_molL) or I_molL <= 0.0:
            return 1.0

        I_m = ionic_strength_molal_from_molar(I_molL, T_K=T_K)
        if not np.isfinite(I_m) or I_m <= 0.0:
            return 1.0

        A = debye_huckel_A(T_K)
        sqrtI = math.sqrt(I_m)
        return 10.0 ** (-A * z * z * sqrtI / (1.0 + 1.5 * sqrtI))

    def jacobian_dgamma_dx(
        self,
        x_mol: Dict[str, float],
        T_K: float,
        *,
        charge: Dict[str, int],
    ) -> np.ndarray:
        """``∂γ_i/∂C_j = (∂γ_i/∂I)(z_j²/2)`` — §8.4 of THERMODYNAMIC_MODEL_ARCHITECTURE.md.

        Differentiates only the extended Debye-Hückel term's I-dependence
        analytically (matching :meth:`gamma`'s DH-only treatment); the
        ion-pair ``ε`` cross-terms used by :meth:`gamma_all` are not
        differentiated (each ``ε(j,k)·m_k`` term is itself linear in a
        *different* species' composition, not I — a full treatment would
        need the per-pair partials, which this "cheap analytically"
        extension (per the design doc) does not attempt). See
        :class:`~PyOMES.thermo.liquid_phase_model.DifferentiableLiquidModel`
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
        sqrtIm = math.sqrt(I_m)
        dIm_dImolL = 1.0 / _kg_per_L(T_K)
        # d(log10 gamma_DH)/dI_m, from log10(gamma_DH) = -A z^2 sqrt(Im)/(1+1.5 sqrt(Im))
        d_log10_dIm = 1.0 / (2.0 * sqrtIm * (1.0 + 1.5 * sqrtIm) ** 2)

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

    def compute_gammas(
        self,
        composition_molL: Dict[str, float],
        I_molL: float,
        T_K: float,
    ) -> Dict[str, float]:
        """Full SIT for the fixed ``ION_CHARGES`` ion set (backward compat).

        Called by NRChemicalEquilibriumEngine when an SIT model is active.  Returns
        gammas for every ion in ``ION_CHARGES``; ions absent from
        ``composition_molL`` still receive the DH-only term.
        """
        T_K = float(T_K)
        I_m = ionic_strength_molal_from_molar(float(I_molL), T_K=T_K)
        if not np.isfinite(I_m) or I_m <= 0.0:
            return {ion: 1.0 for ion in ION_CHARGES}

        A = debye_huckel_A(T_K)
        sqrtI = math.sqrt(I_m)
        dh_factor = -A * sqrtI / (1.0 + 1.5 * sqrtI)

        rho_kg_L = water_density_kg_per_m3(T_K) / 1000.0
        comp_m = {
            ion: max(0.0, float(c)) / rho_kg_L
            for ion, c in composition_molL.items()
            if float(c) > 0
        }
        cations = {
            ion: m for ion, m in comp_m.items()
            if ION_CHARGES.get(ion, 0) > 0
        }
        anions = {
            ion: m for ion, m in comp_m.items()
            if ION_CHARGES.get(ion, 0) < 0
        }

        result = {}
        for ion, z in ION_CHARGES.items():
            log10_gamma = dh_factor * z * z
            if z > 0:
                for anion, m_k in anions.items():
                    log10_gamma += _get_epsilon(self.epsilon, ion, anion) * m_k
            elif z < 0:
                for cation, m_k in cations.items():
                    log10_gamma += _get_epsilon(self.epsilon, cation, ion) * m_k
            result[ion] = 10.0**log10_gamma
        return result
