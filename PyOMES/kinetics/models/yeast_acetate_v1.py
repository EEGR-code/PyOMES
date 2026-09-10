# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 17:12:07 2026

@author: k2473520
"""

# v8_split_8/kinetics/models/yeast_acetate_v1.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Mapping, Tuple
import numpy as np

from ..core import KineticModel, ParameterSpec, Environment


@dataclass(frozen=True)
class YeastAcetateV1:
    """
    Wraps the existing RHS you currently compute inside fermenter_unit.py.
    You will copy the math into rhs(...) so the fermenter can call it.
    """
    name: str = "yeast_acetate_v1"

    # Choose state IDs that match your current y-vector meaning.
    # Example (change to your real ones):
    state_ids: Tuple[str, ...] = (
        "X_gDW",          # biomass amount (or mass)
        "Acetate_mol",    # total acetate in reactor liquid (moles)
        "O2_mol",         # dissolved O2 in liquid (moles) or total O2 in liquid
        "CO2_mol",        # dissolved CO2 in liquid (moles)
    )

    # Define the parameters your current RHS uses:
    param_schema: Tuple[ParameterSpec, ...] = (
        ParameterSpec("mu_max_per_h", default=0.30, lower=0.0),
        ParameterSpec("Ks_acetate_mol_L", default=0.01, lower=0.0),
        ParameterSpec("Ko2_mol_L", default=1e-4, lower=0.0),
        ParameterSpec("Yxs_gDW_per_molAc", default=10.0, lower=0.0),
        ParameterSpec("m_maint_molAc_per_gDW_h", default=0.0, lower=0.0),
    )

    def rhs(
        self,
        t_h: float,
        y: np.ndarray,
        env: Environment,
        params: Mapping[str, float],
    ):
        # Unpack states
        X = float(y[0])
        Ac_mol = float(y[1])
        O2_mol = float(y[2])
        CO2_mol = float(y[3])

        # Convert to concentrations if you need Monod terms
        V = max(env.V_L, 1e-12)
        Ac = max(Ac_mol / V, 0.0)
        O2 = max(env.DO_mol_L, 0.0)  # you may use env DO instead of O2_mol/V

        # Unpack parameters
        mu_max = float(params["mu_max_per_h"])
        Ks = float(params["Ks_acetate_mol_L"])
        Ko2 = float(params["Ko2_mol_L"])
        Yxs = float(params["Yxs_gDW_per_molAc"])
        m = float(params["m_maint_molAc_per_gDW_h"])

        # --- COPY YOUR EXISTING KINETIC EXPRESSIONS HERE ---
        # Example only:
        monod_S = Ac / (Ks + Ac) if (Ks + Ac) > 0 else 0.0
        monod_O2 = O2 / (Ko2 + O2) if (Ko2 + O2) > 0 else 0.0
        mu = mu_max * monod_S * monod_O2

        dX = mu * X
        # substrate uptake supporting growth + maintenance (example)
        qS = (mu / max(Yxs, 1e-12)) + m
        dAc = -(qS * X) * V  # mol/h (since X in gDW, qS in mol/gDW/h, times V?? adjust to your units)

        # O2/CO2 changes may be handled elsewhere (gas transfer + speciation),
        # but you can include net biological consumption/production terms here if you already do.
        dO2 = 0.0
        dCO2 = 0.0

        dydt = np.array([dX, dAc, dO2, dCO2], dtype=float)

        outputs: Dict[str, float] = {
            "mu_per_h": mu,
            "qS_mol_per_gDW_h": qS,
        }
        return dydt, outputs
