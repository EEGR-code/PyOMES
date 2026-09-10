# -*- coding: utf-8 -*-
"""Peng-Robinson (1976) cubic equation of state for gas mixtures.

Stage 13: Non-ideal gas EoS for pressurised AD, biogas compression,
and biogas upgrading applications.  At atmospheric pressure (1 atm),
the ideal gas law gives errors <1.2% — PR is only needed above ~5 atm.

Implements the :class:`~fermenter.equilibria.vle.GasEOS` protocol so it
plugs into the existing VLE and equilibria infrastructure without
modifying any other code.

Usage
-----
>>> from PyOMES.equilibria.peng_robinson import PengRobinsonEOS, BIOGAS_SPECIES
>>> eos = PengRobinsonEOS(BIOGAS_SPECIES)
>>> P = eos.pressure_atm(0.1, T_K=308.15, V_L=0.4)
>>> p = eos.partial_pressures_atm({"CH4": 0.06, "CO2": 0.04}, T_K=308.15, V_L=0.4)

The ``partial_pressures_atm`` method returns fugacities (f_i = y_i × φ_i × P),
which is the correct thermodynamic driving force for Henry-law VLE.
For ideal gas conditions (low pressure), φ_i → 1 and this reduces to
y_i × P (identical to IdealGasEOS).

References
----------
- Peng & Robinson, Ind. Eng. Chem. Fundam. 15(1), 59-64 (1976)
- Reid, Prausnitz & Poling, The Properties of Gases and Liquids, 5th ed.
- Michelsen & Mollerup, Thermodynamic Models (2007), Ch. 3-4
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

# Gas constant in L·atm/(mol·K)
R = 0.0820574

# PR constants
OMEGA_A = 0.45724
OMEGA_B = 0.07780


# ════════════════════════════════════════════════════════════════════════
#  Species critical properties
# ════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class CriticalProperties:
    """Pure-component critical properties for Peng-Robinson.

    Parameters
    ----------
    Tc : float
        Critical temperature (K).
    Pc : float
        Critical pressure (atm).
    omega : float
        Pitzer acentric factor (dimensionless).
    """
    Tc: float
    Pc: float
    omega: float


# Built-in critical properties for biogas species
BIOGAS_SPECIES: Dict[str, CriticalProperties] = {
    "CH4":         CriticalProperties(Tc=190.56, Pc=45.99,  omega=0.012),
    "CO2":         CriticalProperties(Tc=304.21, Pc=72.83,  omega=0.224),
    "H2":          CriticalProperties(Tc=33.19,  Pc=12.97,  omega=-0.216),
    "N2":          CriticalProperties(Tc=126.20, Pc=33.54,  omega=0.037),
    "O2":          CriticalProperties(Tc=154.58, Pc=49.77,  omega=0.022),
    "H2S":         CriticalProperties(Tc=373.53, Pc=88.20,  omega=0.094),
    "NH3":         CriticalProperties(Tc=405.50, Pc=111.30, omega=0.253),
    "H2O":         CriticalProperties(Tc=647.14, Pc=217.75, omega=0.345),
    "AceticAcid":  CriticalProperties(Tc=592.71, Pc=57.07,  omega=0.467),
    "Propionate":  CriticalProperties(Tc=600.81, Pc=46.17,  omega=0.536),
    "Butyrate":    CriticalProperties(Tc=615.70, Pc=39.34,  omega=0.600),
    "Valerate":    CriticalProperties(Tc=639.60, Pc=33.50,  omega=0.650),
    "Ethanol":     CriticalProperties(Tc=513.92, Pc=61.48,  omega=0.645),
}

# Published binary interaction parameters (kij) for biogas mixtures
# Sources: Knapp et al. (1982), Carroll & Mather (1995),
#          Søreide & Whitson (1992), Guillevic et al. (1985)
BIOGAS_KIJ: Dict[Tuple[str, str], float] = {
    ("CH4", "CO2"):  0.092,
    ("CH4", "H2S"):  0.080,
    ("CH4", "H2O"):  0.485,
    ("CH4", "N2"):   0.031,
    ("CH4", "H2"):   0.001,
    ("CO2", "H2S"):  0.096,
    ("CO2", "H2O"):  0.190,
    ("CO2", "N2"):  -0.017,
    ("CO2", "H2"):  -0.160,
    ("H2S", "H2O"):  0.040,
    ("NH3", "H2O"): -0.256,
    ("N2",  "O2"):   0.000,
    # Ethanol pairs (Reid, Prausnitz & Poling; DECHEMA VLE Data Collection)
    ("CH4",  "Ethanol"):  0.035,
    ("CO2",  "Ethanol"):  0.095,
    ("H2O",  "Ethanol"): -0.070,
    ("N2",   "Ethanol"):  0.080,
}


def _get_kij(kij_dict: Dict[Tuple[str, str], float], sp1: str, sp2: str) -> float:
    """Look up kij (symmetric: kij = kji)."""
    if sp1 == sp2:
        return 0.0
    return kij_dict.get((sp1, sp2), kij_dict.get((sp2, sp1), 0.0))


# ════════════════════════════════════════════════════════════════════════
#  Pure-component PR parameters
# ════════════════════════════════════════════════════════════════════════

def _kappa(omega: float) -> float:
    """PR κ parameter from acentric factor."""
    return 0.37464 + 1.54226 * omega - 0.26992 * omega * omega


def _alpha(T_K: float, Tc: float, omega: float) -> float:
    """PR α(T) function."""
    k = _kappa(omega)
    Tr_sqrt = math.sqrt(max(T_K / Tc, 1e-30))
    return (1.0 + k * (1.0 - Tr_sqrt)) ** 2


def _a_pure(Tc: float, Pc: float) -> float:
    """PR a_i at critical point (without α)."""
    return OMEGA_A * R * R * Tc * Tc / Pc


def _b_pure(Tc: float, Pc: float) -> float:
    """PR b_i."""
    return OMEGA_B * R * Tc / Pc


# ════════════════════════════════════════════════════════════════════════
#  Cubic root solver
# ════════════════════════════════════════════════════════════════════════

def _solve_cubic_Z(A: float, B: float) -> float:
    """Solve the PR cubic in Z and return the vapour (largest real) root.

    The PR equation in terms of Z = PV/(nRT) is:
        Z³ - (1-B)Z² + (A - 3B² - 2B)Z - (AB - B² - B³) = 0

    Parameters
    ----------
    A : float
        a_mix × P / (R²T²)
    B : float
        b_mix × P / (RT)

    Returns
    -------
    float
        Largest real root (vapour root).
    """
    # Coefficients: Z³ + c2·Z² + c1·Z + c0 = 0
    c2 = -(1.0 - B)
    c1 = A - 3.0 * B * B - 2.0 * B
    c0 = -(A * B - B * B - B * B * B)

    # Cardano's method
    p = c1 - c2 * c2 / 3.0
    q = c0 - c2 * c1 / 3.0 + 2.0 * c2 * c2 * c2 / 27.0

    discriminant = q * q / 4.0 + p * p * p / 27.0

    if discriminant > 0:
        # One real root
        sqrt_disc = math.sqrt(discriminant)
        u = -q / 2.0 + sqrt_disc
        v = -q / 2.0 - sqrt_disc
        u_sign = 1.0 if u >= 0 else -1.0
        v_sign = 1.0 if v >= 0 else -1.0
        Z = (u_sign * abs(u) ** (1.0 / 3.0) +
             v_sign * abs(v) ** (1.0 / 3.0) - c2 / 3.0)
        return max(Z, B + 1e-10)  # Z must be > B
    else:
        # Three real roots — take the largest (vapour root)
        if abs(p) < 1e-30:
            Z = -c2 / 3.0
            return max(Z, B + 1e-10)
        r = math.sqrt(max(-p * p * p / 27.0, 0.0))
        if r < 1e-30:
            Z = -c2 / 3.0
            return max(Z, B + 1e-10)
        theta = math.acos(max(-1.0, min(1.0, -q / (2.0 * r))))
        r_cbrt = r ** (1.0 / 3.0)
        roots = []
        for k in range(3):
            Zk = 2.0 * r_cbrt * math.cos((theta + 2.0 * math.pi * k) / 3.0) - c2 / 3.0
            if Zk > B:
                roots.append(Zk)
        if roots:
            return max(roots)
        return max(-c2 / 3.0, B + 1e-10)


# ════════════════════════════════════════════════════════════════════════
#  PengRobinsonEOS class
# ════════════════════════════════════════════════════════════════════════

@dataclass
class PengRobinsonEOS:
    """Peng-Robinson (1976) cubic equation of state for gas mixtures.

    Implements the :class:`~fermenter.equilibria.vle.GasEOS` protocol.

    Parameters
    ----------
    species_params : dict
        ``{species_id: CriticalProperties}`` for each species that may
        appear in the gas phase.  Use :data:`BIOGAS_SPECIES` for the
        built-in biogas parameter set.
    kij : dict, optional
        Binary interaction parameters ``{(sp1, sp2): kij}``.
        Symmetric: only one ordering needed.  Default: :data:`BIOGAS_KIJ`.

    Examples
    --------
    >>> eos = PengRobinsonEOS(BIOGAS_SPECIES)
    >>> P = eos.pressure_atm(0.1, T_K=308.15, V_L=0.4)
    >>> pp = eos.partial_pressures_atm({"CH4": 0.06, "CO2": 0.04},
    ...                                 T_K=308.15, V_L=0.4)
    """
    species_params: Dict[str, CriticalProperties]
    kij: Dict[Tuple[str, str], float] = field(default_factory=lambda: dict(BIOGAS_KIJ))

    def _mixture_params(
        self,
        n_gas_mol: Dict[str, float],
        T_K: float,
    ) -> Tuple[float, float, Dict[str, float], Dict[str, float]]:
        """Compute mixture a_mix, b_mix, and per-species a_i*alpha_i, b_i.

        Returns (a_mix, b_mix, a_alpha_i dict, b_i dict).
        """
        species_list = [sp for sp, n in n_gas_mol.items() if n > 0]
        n_total = sum(max(0.0, n_gas_mol.get(sp, 0.0)) for sp in species_list)
        if n_total <= 0:
            return 0.0, 0.0, {}, {}

        y = {sp: max(0.0, n_gas_mol.get(sp, 0.0)) / n_total for sp in species_list}

        # Per-species parameters
        a_alpha = {}
        b_i = {}
        for sp in species_list:
            cp = self.species_params.get(sp)
            if cp is None:
                # Unknown species — treat as ideal (a=0, b=0)
                a_alpha[sp] = 0.0
                b_i[sp] = 0.0
            else:
                a_alpha[sp] = _a_pure(cp.Tc, cp.Pc) * _alpha(T_K, cp.Tc, cp.omega)
                b_i[sp] = _b_pure(cp.Tc, cp.Pc)

        # van der Waals one-fluid mixing rules with kij
        a_mix = 0.0
        for i in species_list:
            for j in species_list:
                kij_val = _get_kij(self.kij, i, j)
                a_ij = math.sqrt(max(a_alpha[i] * a_alpha[j], 0.0)) * (1.0 - kij_val)
                a_mix += y[i] * y[j] * a_ij

        b_mix = sum(y[sp] * b_i[sp] for sp in species_list)

        return a_mix, b_mix, a_alpha, b_i

    def pressure_atm(
        self,
        n_tot_mol: float,
        *,
        T_K: float,
        V_L: float,
    ) -> float:
        """Total pressure from PR equation (atm).

        For a single-component or when composition is unknown, this uses
        the ideal gas law as a fallback (consistent with the protocol
        requirement that total pressure is computed from n_total).
        """
        V_L = max(float(V_L), 1e-30)
        n = max(float(n_tot_mol), 0.0)
        if n <= 0:
            return 0.0
        # Without composition, fall back to ideal gas
        # (pressure_atm takes n_total only — no composition info)
        return n * R * float(T_K) / V_L

    def pressure_mixture_atm(
        self,
        n_gas_mol: Dict[str, float],
        *,
        T_K: float,
        V_L: float,
    ) -> float:
        """Total pressure from PR equation using full composition (atm).

        This is the composition-aware pressure calculation.

        P = nRT/(V - n*b_mix) - n²*a_mix / (V² + 2*n*b_mix*V - n²*b_mix²)
        """
        V = max(float(V_L), 1e-30)
        T = float(T_K)
        n_total = sum(max(0.0, float(v)) for v in n_gas_mol.values())
        if n_total <= 0:
            return 0.0

        a_mix, b_mix, _, _ = self._mixture_params(n_gas_mol, T)
        if b_mix <= 0:
            return n_total * R * T / V

        nb = n_total * b_mix
        V_nb = V - nb
        if V_nb <= 0:
            V_nb = 1e-10  # avoid singularity

        P = (n_total * R * T / V_nb
             - n_total * n_total * a_mix / (V * V + 2.0 * nb * V - nb * nb))
        return max(P, 0.0)

    def partial_pressures_atm(
        self,
        n_gas_mol: Dict[str, float],
        *,
        T_K: float,
        V_L: float,
    ) -> Dict[str, float]:
        """Fugacity-corrected partial pressures (atm).

        Returns ``{species: f_i}`` where ``f_i = y_i × φ_i × P`` is the
        fugacity of species *i*.  For ideal gas (φ=1), this reduces to
        ``y_i × P``.  For non-ideal gas, the fugacity coefficient φ_i
        captures the deviation from ideality.

        These fugacity values are the correct thermodynamic driving
        force for Henry-law VLE calculations.
        """
        T = float(T_K)
        V = max(float(V_L), 1e-30)
        species_list = [sp for sp, n in n_gas_mol.items() if float(n) > 0]
        n_total = sum(max(0.0, float(n_gas_mol.get(sp, 0.0))) for sp in species_list)

        if n_total <= 0:
            return {sp: 0.0 for sp in n_gas_mol}

        y = {sp: max(0.0, float(n_gas_mol.get(sp, 0.0))) / n_total for sp in species_list}

        a_mix, b_mix, a_alpha, b_i = self._mixture_params(n_gas_mol, T)

        # Total pressure (use composition-aware version)
        P = self.pressure_mixture_atm(n_gas_mol, T_K=T, V_L=V)
        if P <= 0:
            return {sp: 0.0 for sp in n_gas_mol}

        # At very low pressure, fugacity coefficients → 1
        if P < 0.01 or b_mix <= 0:
            return {sp: y.get(sp, 0.0) * P for sp in n_gas_mol}

        # Compressibility factor
        A = a_mix * P / (R * R * T * T)
        B = b_mix * P / (R * T)

        Z = _solve_cubic_Z(A, B)

        # Precompute sum_j(y_j * a_ij) for each i
        sum_ya = {}
        for i in species_list:
            s = 0.0
            for j in species_list:
                kij_val = _get_kij(self.kij, i, j)
                a_ij = math.sqrt(max(a_alpha.get(i, 0.0) * a_alpha.get(j, 0.0), 0.0)) * (1.0 - kij_val)
                s += y[j] * a_ij
            sum_ya[i] = s

        # Fugacity coefficients (Michelsen & Mollerup, eq. 3.37)
        sqrt2 = math.sqrt(2.0)
        ZpB1 = Z + (1.0 + sqrt2) * B
        ZpB2 = Z + (1.0 - sqrt2) * B
        if ZpB1 <= 0 or ZpB2 <= 0 or (Z - B) <= 0:
            # Fallback to ideal
            return {sp: y.get(sp, 0.0) * P for sp in n_gas_mol}

        log_ratio = math.log(ZpB1 / ZpB2)

        result = {}
        for sp in n_gas_mol:
            if sp not in species_list:
                result[sp] = 0.0
                continue

            bi = b_i.get(sp, 0.0)
            bi_over_bm = bi / b_mix if b_mix > 0 else 0.0
            sum_ya_i = sum_ya.get(sp, 0.0)

            ln_phi = (bi_over_bm * (Z - 1.0)
                      - math.log(Z - B)
                      - A / (2.0 * sqrt2 * B) * (2.0 * sum_ya_i / a_mix - bi_over_bm) * log_ratio)

            phi = math.exp(min(ln_phi, 50.0))  # clamp to avoid overflow
            result[sp] = y[sp] * phi * P

        return result

    def fugacity_coefficients(
        self,
        n_gas_mol: Dict[str, float],
        *,
        T_K: float,
        V_L: float,
    ) -> Dict[str, float]:
        """Compute fugacity coefficients φ_i for each species.

        Returns
        -------
        dict
            ``{species: phi_i}`` where φ_i = 1.0 for ideal gas.
        """
        T = float(T_K)
        V = max(float(V_L), 1e-30)
        species_list = [sp for sp, n in n_gas_mol.items() if float(n) > 0]
        n_total = sum(max(0.0, float(n_gas_mol.get(sp, 0.0))) for sp in species_list)

        if n_total <= 0:
            return {sp: 1.0 for sp in n_gas_mol}

        P = self.pressure_mixture_atm(n_gas_mol, T_K=T, V_L=V)
        if P < 0.01:
            return {sp: 1.0 for sp in n_gas_mol}

        pp = self.partial_pressures_atm(n_gas_mol, T_K=T, V_L=V)
        y = {sp: max(0.0, float(n_gas_mol.get(sp, 0.0))) / n_total for sp in species_list}

        result = {}
        for sp in n_gas_mol:
            yi = y.get(sp, 0.0)
            if yi > 0 and P > 0:
                result[sp] = pp.get(sp, 0.0) / (yi * P)
            else:
                result[sp] = 1.0
        return result

    def compressibility_factors(
        self,
        n_gas_mol: Dict[str, float],
        *,
        T_K: float,
        V_L: float,
    ) -> Tuple[float, Dict[str, float]]:
        """Compute mixture Z-factor and per-species fugacity coefficients.

        Returns
        -------
        tuple
            ``(Z_mix, {species: phi_i})``
        """
        T = float(T_K)
        V = max(float(V_L), 1e-30)
        n_total = sum(max(0.0, float(v)) for v in n_gas_mol.values())
        if n_total <= 0:
            return 1.0, {sp: 1.0 for sp in n_gas_mol}

        a_mix, b_mix, _, _ = self._mixture_params(n_gas_mol, T)
        P = self.pressure_mixture_atm(n_gas_mol, T_K=T, V_L=V)
        if P <= 0 or b_mix <= 0:
            return 1.0, {sp: 1.0 for sp in n_gas_mol}

        A = a_mix * P / (R * R * T * T)
        B = b_mix * P / (R * T)
        Z = _solve_cubic_Z(A, B)

        phis = self.fugacity_coefficients(n_gas_mol, T_K=T, V_L=V)
        return Z, phis
