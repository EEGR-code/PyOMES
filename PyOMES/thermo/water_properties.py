# -*- coding: utf-8 -*-
"""Water property correlations for temperature-dependent Debye–Hückel constants.

Lightweight polynomial correlations valid over 0–100 °C (bioprocess range).
Kept in PyOMES/thermo/ so activity models in this package can use them without
importing from PyOMES/speciation/ (which would create a circular dependency).
"""
from __future__ import annotations

import numpy as np


def water_dielectric_constant(T_K: float) -> float:
    """Relative dielectric constant of pure water (dimensionless).

    Common polynomial correlation in T(°C), accurate over ~0–100 °C.
    """
    T_C = float(T_K) - 273.15
    return float(87.740 - 0.40008 * T_C + 9.398e-4 * T_C**2 - 1.410e-6 * T_C**3)


def water_density_kg_per_m3(T_K: float) -> float:
    """Density of pure water (kg/m³) vs temperature.

    Kell-type correlation in T(°C), accurate over ~0–100 °C.
    """
    T_C = float(T_K) - 273.15
    return float(
        1000.0
        * (
            1.0
            - ((T_C + 288.9414) / (508929.2 * (T_C + 68.12963)))
            * (T_C - 3.9863) ** 2
        )
    )


def debye_huckel_A(T_K: float) -> float:
    """Debye–Hückel constant A (base-10) for water at T_K.

    Returns A for molality-based ionic strength (mol/kg water):
        log10(γ_i) = −A z² f(I)

    The standard expression uses ρ in g/cm³; water_density_kg_per_m3()
    returns kg/m³, so we divide by 1000 before use.  A ≈ 0.509 at 25 °C.
    """
    T = float(T_K)
    rho_kg_m3 = water_density_kg_per_m3(T)
    eps_r = water_dielectric_constant(T)

    if (
        not np.isfinite(rho_kg_m3)
        or not np.isfinite(eps_r)
        or rho_kg_m3 <= 0.0
        or eps_r <= 0.0
        or T <= 0.0
    ):
        return 0.509  # safe fallback: ~25 °C value

    rho_g_cm3 = rho_kg_m3 / 1000.0
    return float(1.82483e6 * np.sqrt(rho_g_cm3) / (eps_r * T) ** 1.5)


def ionic_strength_molal_from_molar(I_molL: float, *, T_K: float) -> float:
    """Convert ionic strength from mol/L to mol/kg-water using water density."""
    I = float(I_molL)
    if not np.isfinite(I) or I <= 0.0:
        return 0.0
    rho = water_density_kg_per_m3(float(T_K))  # kg/m³
    kg_per_L = rho / 1000.0
    if not np.isfinite(kg_per_L) or kg_per_L <= 0.0:
        return I
    return float(I / kg_per_L)
