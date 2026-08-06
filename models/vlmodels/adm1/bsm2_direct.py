# -*- coding: utf-8 -*-
"""BSM2-native ADM1 model with explicit carbon/nitrogen balance.

This module provides :class:`BSM2DirectModel`, which implements the
BSM2 ADM1 Petersen matrix in its native COD-fraction units and
computes inorganic carbon (S_IC) and nitrogen (S_IN) production from
explicit carbon and nitrogen balances — exactly as BSM2 does.

This eliminates the TIC drift that arises from the fermenter's
default elemental closure approach (``_petersen_to_molar``), which
computes CO₂ stoichiometric coefficients by simultaneously balancing
C, H, O, and N, causing small residuals to accumulate in the CO₂
term.

The model satisfies the :class:`~fermenter.reactions.protocols.ReactionModel`
protocol and can be attached directly to a fermenter ``ControlVolume``.

Usage
-----
>>> from PyOMES.models.adm1_bsm2_direct import build_bsm2_direct_model
>>> model = build_bsm2_direct_model()
>>> cv = build_bsm2_cv_with_model(model, ...)

Reference
---------
Rosen, C. & Jeppsson, U. (2006). "Aspects on ADM1 Implementation
within the BSM2 Framework."  Technical Report, Dept. of Industrial
Electrical Engineering and Automation, Lund University.
"""

from __future__ import annotations

import math
from typing import Any, Dict, Optional

from vlmodels.adm1.bsm2 import (
    SPECIES_BSM2, BSM2_BIO_ATOMS, BSM2_BIO_THOD, BSM2_ORGS,
    BSM2_STOICH, BSM2_KINETICS,
    _thod, _atoms, _nC, _nN,
)
from PyOMES.reactions.environment import ReactionEnvironment


# ════════════════════════════════════════════════════════════════════════
#  Carbon and nitrogen content tables (mol_element / gCOD)
# ════════════════════════════════════════════════════════════════════════
# These are used to compute the S_IC and S_IN rows of the Petersen
# matrix via explicit element balances.
#
# IMPORTANT: BSM2 defines carbon and nitrogen content coefficients as
# independent parameters (Rosen & Jeppsson 2006, Table 6.1), NOT
# derived from molecular formulae.  Several species have C/N content
# that differs significantly from their molecular formula because:
#
#   - C_xc is chosen to close the carbon balance for disintegration
#     (C_xc = Σ f_i × C_i for products), NOT from X_xc's formula.
#   - N_xc is similarly a weighted average of product N contents.
#   - S_I, X_I use N_I = 0.06/14, not the biomass formula N content.
#   - BSM2 uses rounded empirical values for several species (C_aa,
#     C_pr, C_li, C_sI, C_xI) that don't exactly match their ThOD.
#
# Using molecular-formula-derived values causes systematic drift in
# S_IC (~2+ mM/d) and S_IN (~1.8+ mM/d) that compounds to ~10%
# error over 10 days.  This was the root cause of the BSM2 benchmark
# divergence (investigated across 13 hypotheses before discovery).

def _build_C_content_bsm2() -> Dict[str, float]:
    """Carbon content: mol_C per gCOD for each species.

    Uses BSM2's canonical values (Rosen & Jeppsson 2006, Table 6.1),
    NOT molecular formula derivations.
    """
    # BSM2 Table 6.1 carbon content coefficients (kmol_C / kgCOD = mol_C / gCOD)
    cc = {
        "S_su":  0.0313,   # C_su
        "S_aa":  0.03,     # C_aa
        "S_fa":  0.0217,   # C_fa
        "S_va":  0.024,    # C_va
        "S_bu":  0.025,    # C_bu
        "S_pro": 0.0268,   # C_pro
        "S_ac":  0.0313,   # C_ac
        "S_h2":  0.0,      # no carbon
        "S_ch4": 0.0156,   # C_ch4
        "S_I":   0.03,     # C_sI
        "X_xc":  0.02786,  # C_xc (weighted average of product C contents)
        "X_ch":  0.0313,   # C_ch
        "X_pr":  0.03,     # C_pr
        "X_li":  0.022,    # C_li
        "X_I":   0.03,     # C_xI
    }
    # All biomass organisms use C_bac
    for org in BSM2_ORGS:
        cc[org] = 0.0313   # C_bac
    return cc


def _build_N_content_bsm2() -> Dict[str, float]:
    """Nitrogen content: mol_N per gCOD for each species.

    Uses BSM2's canonical values (Rosen & Jeppsson 2006, Table 6.1).
    """
    # BSM2 nitrogen content coefficients (kmol_N / kgCOD = mol_N / gCOD)
    nc = {
        "S_su":  0.0,
        "S_aa":  0.007,     # N_aa
        "S_fa":  0.0,
        "S_va":  0.0,
        "S_bu":  0.0,
        "S_pro": 0.0,
        "S_ac":  0.0,
        "S_h2":  0.0,
        "S_ch4": 0.0,
        "S_I":   0.06 / 14,   # N_I = 0.004286
        "X_xc":  0.0376 / 14, # N_xc = 0.002686 (weighted average)
        "X_ch":  0.0,
        "X_pr":  0.007,     # N_aa (proteins have same N as amino acids)
        "X_li":  0.0,
        "X_I":   0.06 / 14, # N_I = 0.004286
    }
    # All biomass organisms use N_bac
    for org in BSM2_ORGS:
        nc[org] = 0.08 / 14  # N_bac = 0.005714
    return nc


def _build_H_content() -> Dict[str, float]:
    """Hydrogen content: mol_H per gCOD for each species.

    Derived from molecular formulae (BSM2 does not define explicit H
    content coefficients — H₂O is a secondary balance closure).
    """
    hc = {}
    for sp, (at, mw, thod) in SPECIES_BSM2.items():
        if thod > 0:
            hc[sp] = at.get("H", 0) / thod
    for org in BSM2_ORGS:
        hc[org] = BSM2_BIO_ATOMS["H"] / BSM2_BIO_THOD
    return hc


_C_CONTENT = _build_C_content_bsm2()
_N_CONTENT = _build_N_content_bsm2()
_H_CONTENT = _build_H_content()


# ════════════════════════════════════════════════════════════════════════
#  Petersen matrix builder
# ════════════════════════════════════════════════════════════════════════

def _build_petersen_rows(p):
    """Build the Petersen matrix as a list of (label, row_dict) tuples.

    Each row_dict maps ``{species: coefficient}`` in COD fractions.
    Coefficients are dimensionless: gCOD_product / gCOD_substrate
    consumed (or gCOD_consumed / gCOD_consumed = -1 for substrate).

    The process rate ρ_j [gCOD/(L·d)] multiplied by ν_ij gives the
    rate of change dS_i/dt [gCOD/(L·d)] for species i.
    """
    rows = []

    # ── R1: Disintegration ────────────────────────────────────────
    # X_xc → f_ch X_ch + f_pr X_pr + f_li X_li + f_xI X_I + f_sI S_I
    rows.append(("disintegration", {
        "X_xc": -1.0,
        "X_ch": p["f_ch_xc"],
        "X_pr": p["f_pr_xc"],
        "X_li": p["f_li_xc"],
        "X_I":  p["f_xI_xc"],
        "S_I":  p["f_sI_xc"],
    }))

    # ── R2: Hydrolysis carbohydrates ──────────────────────────────
    # X_ch → S_su (1:1 in COD)
    rows.append(("hyd_ch", {"X_ch": -1.0, "S_su": 1.0}))

    # ── R3: Hydrolysis proteins ───────────────────────────────────
    # X_pr → S_aa (1:1 in COD)
    rows.append(("hyd_pr", {"X_pr": -1.0, "S_aa": 1.0}))

    # ── R4: Hydrolysis lipids ─────────────────────────────────────
    # X_li → (1-f_fa_li) S_su + f_fa_li S_fa
    rows.append(("hyd_li", {
        "X_li": -1.0,
        "S_su": 1.0 - p["f_fa_li"],
        "S_fa": p["f_fa_li"],
    }))

    # ── R5: Sugar uptake ──────────────────────────────────────────
    rows.append(("uptake_su", {
        "S_su": -1.0,
        "S_bu":  (1 - p["Y_su"]) * p["f_bu_su"],
        "S_pro": (1 - p["Y_su"]) * p["f_pro_su"],
        "S_ac":  (1 - p["Y_su"]) * p["f_ac_su"],
        "S_h2":  (1 - p["Y_su"]) * p["f_h2_su"],
        "X_su":  p["Y_su"],
    }))

    # ── R6: Amino acid uptake ─────────────────────────────────────
    rows.append(("uptake_aa", {
        "S_aa": -1.0,
        "S_va":  (1 - p["Y_aa"]) * p["f_va_aa"],
        "S_bu":  (1 - p["Y_aa"]) * p["f_bu_aa"],
        "S_pro": (1 - p["Y_aa"]) * p["f_pro_aa"],
        "S_ac":  (1 - p["Y_aa"]) * p["f_ac_aa"],
        "S_h2":  (1 - p["Y_aa"]) * p["f_h2_aa"],
        "X_aa":  p["Y_aa"],
    }))

    # ── R7: LCFA uptake ───────────────────────────────────────────
    f_h2_fa = 1.0 - p["f_ac_fa"]
    rows.append(("uptake_fa", {
        "S_fa": -1.0,
        "S_ac": (1 - p["Y_fa"]) * p["f_ac_fa"],
        "S_h2": (1 - p["Y_fa"]) * f_h2_fa,
        "X_fa": p["Y_fa"],
    }))

    # ── R8: Valerate uptake (competitive C4) ──────────────────────
    f_h2_va = 1.0 - p["f_pro_va"] - p["f_ac_va"]
    rows.append(("uptake_va", {
        "S_va":  -1.0,
        "S_pro": (1 - p["Y_c4"]) * p["f_pro_va"],
        "S_ac":  (1 - p["Y_c4"]) * p["f_ac_va"],
        "S_h2":  (1 - p["Y_c4"]) * f_h2_va,
        "X_c4":  p["Y_c4"],
    }))

    # ── R9: Butyrate uptake (competitive C4) ──────────────────────
    f_h2_bu = 1.0 - p["f_ac_bu"]
    rows.append(("uptake_bu", {
        "S_bu": -1.0,
        "S_ac": (1 - p["Y_c4"]) * p["f_ac_bu"],
        "S_h2": (1 - p["Y_c4"]) * f_h2_bu,
        "X_c4": p["Y_c4"],
    }))

    # ── R10: Propionate uptake ────────────────────────────────────
    f_h2_pro = 1.0 - p["f_ac_pro"]
    rows.append(("uptake_pro", {
        "S_pro": -1.0,
        "S_ac":  (1 - p["Y_pro"]) * p["f_ac_pro"],
        "S_h2":  (1 - p["Y_pro"]) * f_h2_pro,
        "X_pro": p["Y_pro"],
    }))

    # ── R11: Acetate uptake (methanogenesis) ──────────────────────
    rows.append(("uptake_ac", {
        "S_ac":  -1.0,
        "S_ch4": (1 - p["Y_ac"]),
        "X_ac":  p["Y_ac"],
    }))

    # ── R12: H₂ uptake (methanogenesis) ──────────────────────────
    rows.append(("uptake_h2", {
        "S_h2":  -1.0,
        "S_ch4": (1 - p["Y_h2"]),
        "X_h2":  p["Y_h2"],
    }))

    # ── R13–R19: Decay (7 organisms → X_xc) ──────────────────────
    for org in BSM2_ORGS:
        rows.append((f"decay_{org}", {org: -1.0, "X_xc": 1.0}))

    return rows


# ════════════════════════════════════════════════════════════════════════
#  Rate functions (BSM2-native, gCOD/(L·d))
# ════════════════════════════════════════════════════════════════════════

_R_J = 8.31446


def _pKa_NH4(T_K):
    Ka_ref = 10.0 ** (-9.25)
    Ka_T = Ka_ref * math.exp(-(51965.0 / _R_J) * (1.0 / T_K - 1.0 / 298.15))
    return -math.log10(max(Ka_T, 1e-30))


def _build_process_rates(kp, petersen_rows):
    """Build a list of (label, rate_fn) corresponding to petersen_rows.

    Each rate_fn(env) returns ρ_j in gCOD/(L·d) — the BSM2 process
    rate for that process.
    """
    # ── Inhibition closures ───────────────────────────────────────
    K_pH_aa = 10 ** (-(kp["pH_UL_aa"] + kp["pH_LL_aa"]) / 2.0)
    n_aa = 3.0 / (kp["pH_UL_aa"] - kp["pH_LL_aa"])
    K_pH_ac = 10 ** (-(kp["pH_UL_ac"] + kp["pH_LL_ac"]) / 2.0)
    n_ac = 3.0 / (kp["pH_UL_ac"] - kp["pH_LL_ac"])
    K_pH_h2 = 10 ** (-(kp["pH_UL_h2"] + kp["pH_LL_h2"]) / 2.0)
    n_h2 = 3.0 / (kp["pH_UL_h2"] - kp["pH_LL_h2"])

    K_S_NH3 = kp["K_S_NH3"]
    K_I_nh3 = kp["K_I_nh3"]

    def _H(env):
        # H+ lives on the speciation property channel after
        # chemistry-unification-1, not on phase.n_mol. (B2.)
        if env.has_pH:
            return 10.0 ** (-float(env.pH))
        return 1e-7

    def _I_pH_aa(env):
        H = _H(env)
        return K_pH_aa ** n_aa / (H ** n_aa + K_pH_aa ** n_aa)

    def _I_pH_ac(env):
        H = _H(env)
        return K_pH_ac ** n_ac / (H ** n_ac + K_pH_ac ** n_ac)

    def _I_pH_h2(env):
        H = _H(env)
        return K_pH_h2 ** n_h2 / (H ** n_h2 + K_pH_h2 ** n_h2)

    def _I_IN(env):
        S_IN = env.concentrations.get("NH3", 0)
        return 1.0 / (1.0 + K_S_NH3 / max(S_IN, 1e-12))

    def _I_nh3(env):
        S_IN = env.concentrations.get("NH3", 0)
        pH = float(env.pH) if env.has_pH else 7.0
        pKa = _pKa_NH4(env.T_K)
        f_nh3 = 1.0 / (1.0 + 10.0 ** (pKa - pH))
        S_nh3 = S_IN * f_nh3
        return 1.0 / (1.0 + S_nh3 / K_I_nh3)

    def _I_h2(env, K_I):
        S_h2_cod = env.concentrations.get("S_h2", 0) * _thod("S_h2")
        return 1.0 / (1.0 + S_h2_cod / K_I)

    # Combined inhibition functions (BSM2 naming)
    def I_5(env):  return _I_pH_aa(env) * _I_IN(env)
    def I_7(env):  return _I_pH_aa(env) * _I_IN(env) * _I_h2(env, kp["K_I_h2_fa"])
    def I_8(env):  return _I_pH_aa(env) * _I_IN(env) * _I_h2(env, kp["K_I_h2_c4"])
    def I_10(env): return _I_pH_aa(env) * _I_IN(env) * _I_h2(env, kp["K_I_h2_pro"])
    def I_11(env): return _I_pH_ac(env) * _I_IN(env) * _I_nh3(env)
    def I_12(env): return _I_pH_h2(env) * _I_IN(env)

    # ── Rate builders ─────────────────────────────────────────────

    def _cod(env, sp):
        """Species concentration in gCOD/L."""
        return env.concentrations.get(sp, 0.0) * _thod(sp)

    def _cod_bio(env, org):
        """Biomass concentration in gCOD/L."""
        return env.concentrations.get(org, 0.0) * BSM2_BIO_THOD

    rates = []

    # R1: Disintegration (first-order on X_xc)
    k_dis = kp["k_dis"]
    rates.append(lambda env, _k=k_dis:
                 _k * _cod_bio(env, "X_xc"))

    # R2: Hydrolysis CH
    k_hyd_ch = kp["k_hyd_ch"]
    rates.append(lambda env, _k=k_hyd_ch:
                 _k * _cod(env, "X_ch"))

    # R3: Hydrolysis PR
    k_hyd_pr = kp["k_hyd_pr"]
    rates.append(lambda env, _k=k_hyd_pr:
                 _k * _cod(env, "X_pr"))

    # R4: Hydrolysis LI
    k_hyd_li = kp["k_hyd_li"]
    rates.append(lambda env, _k=k_hyd_li:
                 _k * _cod(env, "X_li"))

    # R5: Sugar uptake
    rates.append(lambda env:
        kp["k_m_su"] * _cod(env, "S_su") / (kp["K_S_su"] + _cod(env, "S_su"))
        * _cod_bio(env, "X_su") * I_5(env)
        if _cod(env, "S_su") > 0 and _cod_bio(env, "X_su") > 0 else 0.0)

    # R6: Amino acid uptake
    rates.append(lambda env:
        kp["k_m_aa"] * _cod(env, "S_aa") / (kp["K_S_aa"] + _cod(env, "S_aa"))
        * _cod_bio(env, "X_aa") * I_5(env)
        if _cod(env, "S_aa") > 0 and _cod_bio(env, "X_aa") > 0 else 0.0)

    # R7: LCFA uptake
    rates.append(lambda env:
        kp["k_m_fa"] * _cod(env, "S_fa") / (kp["K_S_fa"] + _cod(env, "S_fa"))
        * _cod_bio(env, "X_fa") * I_7(env)
        if _cod(env, "S_fa") > 0 and _cod_bio(env, "X_fa") > 0 else 0.0)

    # R8: Valerate uptake (competitive C4)
    def _rate_va(env):
        S_va_cod = _cod(env, "S_va")
        S_bu_cod = _cod(env, "S_bu")
        X_cod = _cod_bio(env, "X_c4")
        if S_va_cod <= 0 or X_cod <= 0:
            return 0.0
        alloc = S_va_cod / (S_va_cod + S_bu_cod + 1e-12)
        return (kp["k_m_c4"] * S_va_cod / (kp["K_S_c4"] + S_va_cod)
                * X_cod * alloc * I_8(env))
    rates.append(_rate_va)

    # R9: Butyrate uptake (competitive C4)
    def _rate_bu(env):
        S_bu_cod = _cod(env, "S_bu")
        S_va_cod = _cod(env, "S_va")
        X_cod = _cod_bio(env, "X_c4")
        if S_bu_cod <= 0 or X_cod <= 0:
            return 0.0
        alloc = S_bu_cod / (S_bu_cod + S_va_cod + 1e-12)
        return (kp["k_m_c4"] * S_bu_cod / (kp["K_S_c4"] + S_bu_cod)
                * X_cod * alloc * I_8(env))
    rates.append(_rate_bu)

    # R10: Propionate uptake
    rates.append(lambda env:
        kp["k_m_pro"] * _cod(env, "S_pro") / (kp["K_S_pro"] + _cod(env, "S_pro"))
        * _cod_bio(env, "X_pro") * I_10(env)
        if _cod(env, "S_pro") > 0 and _cod_bio(env, "X_pro") > 0 else 0.0)

    # R11: Acetate uptake
    rates.append(lambda env:
        kp["k_m_ac"] * _cod(env, "S_ac") / (kp["K_S_ac"] + _cod(env, "S_ac"))
        * _cod_bio(env, "X_ac") * I_11(env)
        if _cod(env, "S_ac") > 0 and _cod_bio(env, "X_ac") > 0 else 0.0)

    # R12: H₂ uptake
    rates.append(lambda env:
        kp["k_m_h2"] * _cod(env, "S_h2") / (kp["K_S_h2"] + _cod(env, "S_h2"))
        * _cod_bio(env, "X_h2") * I_12(env)
        if _cod(env, "S_h2") > 0 and _cod_bio(env, "X_h2") > 0 else 0.0)

    # R13–R19: Decay
    k_dec = kp["k_dec"]
    for org in BSM2_ORGS:
        rates.append(lambda env, _org=org, _k=k_dec:
                     _k * _cod_bio(env, _org))

    return rates


# ════════════════════════════════════════════════════════════════════════
#  BSM2DirectModel
# ════════════════════════════════════════════════════════════════════════

class BSM2DirectModel:
    """BSM2-native ADM1 Petersen matrix model.

    Computes reaction source terms in BSM2's native COD units and
    derives S_IC and S_IN from explicit carbon and nitrogen balances
    (matching BSM2 exactly), avoiding the elemental closure approach
    that causes TIC drift.

    Satisfies the ``ReactionModel`` protocol: attach directly to a
    fermenter ``ControlVolume`` via its ``reaction_system`` field.

    Parameters
    ----------
    stoich : dict, optional
        Override BSM2 stoichiometric parameters.
    kinetics : dict, optional
        Override BSM2 kinetic parameters.
    """

    def __init__(self, stoich: Dict = None, kinetics: Dict = None):
        self._p = dict(BSM2_STOICH)
        if stoich:
            self._p.update(stoich)
        self._kp = dict(BSM2_KINETICS)
        if kinetics:
            self._kp.update(kinetics)

        # Build Petersen matrix: list of (label, {species: coeff})
        self._petersen = _build_petersen_rows(self._p)

        # Build rate functions: list of callables, same order as _petersen
        self._rates = _build_process_rates(self._kp, self._petersen)

        # Pre-compute S_IC and S_IN coefficients for each process
        # ν_IC,j = -Σ_i (C_i × ν_ij)
        # ν_IN,j = -Σ_i (N_i × ν_ij)
        self._nu_IC = []
        self._nu_IN = []
        self._nu_H = []  # for H₂O balance
        for label, row in self._petersen:
            sum_C = sum(_C_CONTENT.get(sp, 0.0) * coeff
                        for sp, coeff in row.items())
            sum_N = sum(_N_CONTENT.get(sp, 0.0) * coeff
                        for sp, coeff in row.items())
            sum_H = sum(_H_CONTENT.get(sp, 0.0) * coeff
                        for sp, coeff in row.items())
            self._nu_IC.append(-sum_C)
            self._nu_IN.append(-sum_N)
            self._nu_H.append(sum_H)

    def compute_rates(
        self, env: ReactionEnvironment
    ) -> Dict[str, Dict[str, float]]:
        """Evaluate all BSM2 process rates and return source terms.

        Parameters
        ----------
        env : ReactionEnvironment
            Current conditions (concentrations in mol/L, pH, T_K, V_L).

        Returns
        -------
        dict
            ``{"liquid": {species_id: rate_mol_per_h}}``.
        """
        V_L = env.V_L

        # 1. Evaluate all 19 process rates ρ_j [gCOD/(L·d)]
        rho = [rate_fn(env) for rate_fn in self._rates]

        # 2. Accumulate dS_i/dt for COD species [gCOD/(L·d)]
        dS_cod: Dict[str, float] = {}
        for j, (label, row) in enumerate(self._petersen):
            rho_j = rho[j]
            if abs(rho_j) < 1e-30:
                continue
            for sp, nu_ij in row.items():
                dS_cod[sp] = dS_cod.get(sp, 0.0) + nu_ij * rho_j

        # 3. S_IC from carbon balance [mol/(L·d)]
        dS_IC = sum(self._nu_IC[j] * rho[j] for j in range(len(rho)))

        # 4. S_IN from nitrogen balance [mol/(L·d)]
        dS_IN = sum(self._nu_IN[j] * rho[j] for j in range(len(rho)))

        # 5. H₂O from hydrogen balance [mol/(L·d)]
        # After CO₂ (0 H) and NH₃ (3 H per mol), remaining H → H₂O:
        #   H_organic = Σ_j (ρ_j × Σ_i(H_i × ν_ij))  [mol_H/(L·d)]
        #   dn_H2O/dt = -(H_organic + 3 × dS_IN) / 2
        H_organic = sum(self._nu_H[j] * rho[j] for j in range(len(rho)))
        dS_H2O = -(H_organic + 3.0 * dS_IN) / 2.0

        # 6. Convert to mol/h (extensive)
        # COD species: rate_mol_h = dS_i [gCOD/(L·d)] / ThOD [gCOD/mol] × V_L / 24
        sources: Dict[str, float] = {}
        factor = V_L / 24.0

        for sp, dS in dS_cod.items():
            if abs(dS) < 1e-30:
                continue
            thod = _thod(sp) if sp in SPECIES_BSM2 else BSM2_BIO_THOD
            if thod <= 0:
                continue
            sources[sp] = dS / thod * factor

        # Inorganics: already in mol/(L·d), just convert to mol/h
        sources["CO2"] = dS_IC * factor
        sources["NH3"] = dS_IN * factor

        # H₂O
        if abs(dS_H2O) > 1e-30:
            sources["H2O"] = dS_H2O * factor

        return {"liquid": sources}

    def process_rates(self, env: ReactionEnvironment) -> Dict[str, float]:
        """Return individual process rates for diagnostics.

        Returns
        -------
        dict
            ``{process_label: rho_gCOD_per_L_per_d}``.
        """
        rho = [rate_fn(env) for rate_fn in self._rates]
        return {label: rho[j] for j, (label, _) in enumerate(self._petersen)}

    def __repr__(self):
        return f"BSM2DirectModel({len(self._petersen)} processes)"


# ════════════════════════════════════════════════════════════════════════
#  Factory function
# ════════════════════════════════════════════════════════════════════════

def build_bsm2_direct_model(
    *,
    stoich: Dict = None,
    kinetics: Dict = None,
) -> BSM2DirectModel:
    """Construct a BSM2-native Petersen matrix model.

    Parameters
    ----------
    stoich : dict, optional
        Override BSM2 stoichiometric parameters.
    kinetics : dict, optional
        Override BSM2 kinetic parameters.

    Returns
    -------
    BSM2DirectModel
        A model satisfying the ``ReactionModel`` protocol.
    """
    return BSM2DirectModel(stoich=stoich, kinetics=kinetics)
