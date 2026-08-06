# -*- coding: utf-8 -*-
"""BSM2-canonical ADM1 implementation in molar units.

Provides an ADM1 model that is structurally consistent with the
BSM2 Matlab/Simulink reference (Rosén & Jeppsson, 2006) while using
the molar-based framework of the fermenter package.

Key differences from the custom ``adm1.py`` model:

1. **7 organism populations** (BSM2 lumping): X_su, X_aa, X_fa, X_c4,
   X_pro, X_ac, X_h2 (not 12 separate populations)
2. **Separate hydrolysis**: 3 first-order hydrolysis reactions + 1
   disintegration, distinct from the 8 Monod uptake reactions
3. **COD-fraction product distributions**: f_bu_su, f_pro_su, etc.
   converted to molar stoichiometric coefficients
4. **Competitive C4 uptake**: X_c4 degrades both butyrate and valerate
   with competitive allocation
5. **kLa-based gas transfer**: kinetic (not equilibrium) for H₂, CH₄, CO₂
6. **Inorganic nitrogen limitation**: I_IN_lim on all uptake reactions
7. **BSM2 parameter values**: k_m, K_S, K_I in correct converted units

Usage
-----
>>> from PyOMES.models.adm1_bsm2 import (
...     build_bsm2_reactions, build_bsm2_cv, seed_bsm2_strong_ions,
... )

Note
----
``make_bsm2_callback`` was removed in chemistry-unification-1: the
stateless snapshot model no longer writes pH back to ``phase.n_mol``.
Equilibrium-output species (H⁺, OH⁻, HCO₃⁻, CO₃²⁻, NH₄⁺) live in
``PropertyResult`` only. Rate functions read pH from the
``ReactionEnvironment`` (which is built from the current step's
property results); users between steps read pH from the most recent
``AdvanceResult``.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple

from PyOMES.reactions.kinetic import KineticReaction
from PyOMES.reactions.equilibrium import EquilibriumReaction
from PyOMES.reactions.reaction_system import ReactionSystem
from PyOMES.reactions.stoichiometry import StoichiometryEntry
from PyOMES.chemistry.species import Species
from PyOMES.chemistry.common_species import (
    H_plus, OH_minus, H2O as H2O_sp,
    CO2 as CO2_sp, HCO3_minus,
    NH3 as NH3_sp, NH4_plus,
)

logger = logging.getLogger(__name__)


# ════════════════════════════════════════════════════════════════════════
#  Species definitions (BSM2-consistent)
# ════════════════════════════════════════════════════════════════════════

# BSM2 biomass: C₅H₇O₂N (conventional activated sludge formula)
BSM2_BIO_ATOMS = {"C": 5, "H": 7, "O": 2, "N": 1}
BSM2_BIO_MW = 113.0
BSM2_BIO_THOD = (4*5 + 7 - 2*2 - 3*1) * 8  # = 160 gCOD/mol

SPECIES_BSM2 = {
    # Soluble substrates
    "S_su":  ({"C": 6, "H": 12, "O": 6},          180.16,  192),   # sugars (glucose)
    "S_aa":  ({"C": 5, "H": 9, "O": 2, "N": 1},   115.13,  176),   # amino acids
    "S_fa":  ({"C": 16, "H": 32, "O": 2},          256.42,  736),   # LCFA (palmitate)
    "S_va":  ({"C": 5, "H": 10, "O": 2},           102.13,  208),   # valerate
    "S_bu":  ({"C": 4, "H": 8, "O": 2},             88.11,  160),   # butyrate
    "S_pro": ({"C": 3, "H": 6, "O": 2},             74.08,  112),   # propionate
    "S_ac":  ({"C": 2, "H": 4, "O": 2},             60.05,   64),   # acetate
    "S_h2":  ({"H": 2},                               2.016,  16),   # hydrogen
    "S_ch4": ({"C": 1, "H": 4},                      16.04,   64),   # methane
    # chemistry-unification-3b: phase-agnostic Species ID convention.
    # "CO2" is molecular CO2(aq) in the liquid phase — same id as gas CO2;
    # phase is tracked by the n_mol dict it lives in, not the key name.
    # The carbonate ladder is (CO2, HCO3-, CO3--); total TIC is implicit.
    # Same for "NH3" — molecular NH3, one species in the ammonia ladder
    # (NH3, NH4+); total ammoniacal N is implicit.
    "CO2": ({"C": 1, "O": 2},                       44.01,    0),   # molecular CO2(aq)
    "NH3":  ({"N": 1, "H": 3},                      17.03,    0),   # molecular NH3
    "S_I":   ({"C": 5, "H": 7, "O": 2, "N": 1},    113.0,  160),   # soluble inert (≈ biomass formula)
    "H2O":   ({"H": 2, "O": 1},                      18.015,   0),
    # Particulates
    "X_xc":  ({"C": 5, "H": 7, "O": 2, "N": 1},   113.0,  160),   # composites
    "X_ch":  ({"C": 6, "H": 10, "O": 5},           162.14,  192),   # carbohydrates
    "X_pr":  ({"C": 5, "H": 9, "O": 2, "N": 1},   115.13,  176),   # proteins
    "X_li":  ({"C": 51, "H": 98, "O": 6},          807.32, 2432),   # lipids (tripalmitin)
    "X_I":   ({"C": 5, "H": 7, "O": 2, "N": 1},   113.0,  160),   # particulate inert
}

# ── Species objects (chemistry-unification-1) ───────────────────────
# Derived from SPECIES_BSM2 so the model-specific ThOD table and the
# universal Species type stay in sync. Inorganic CO₂, NH₃, H₂O alias
# the universal common_species objects so cross-model composition does
# not trip the soft-conflict check.
SPECIES: Dict[str, Species] = {}
for _sp_id, (_atoms_dict, _mw, _thod_val) in SPECIES_BSM2.items():
    if _sp_id == "CO2":
        # chemistry-unification-3b: phase-agnostic id. Use the shared
        # common_species CO2 object directly so no soft-conflict warning.
        SPECIES["CO2"] = CO2_sp
    elif _sp_id == "NH3":
        SPECIES["NH3"] = NH3_sp
    elif _sp_id == "H2O":
        SPECIES["H2O"] = H2O_sp
    else:
        SPECIES[_sp_id] = Species(
            id=_sp_id, atoms=dict(_atoms_dict), charge=0, MW=float(_mw),
        )
# Each of the 7 BSM2 organism populations: same atoms/MW.
_BSM2_BIO_SPECIES = {
    org: Species(id=org, atoms=dict(BSM2_BIO_ATOMS), charge=0, MW=BSM2_BIO_MW)
    for org in ("X_su", "X_aa", "X_fa", "X_c4", "X_pro", "X_ac", "X_h2")
}
SPECIES.update(_BSM2_BIO_SPECIES)
del _sp_id, _atoms_dict, _mw, _thod_val

def _mw(sp):
    return SPECIES_BSM2[sp][1] if sp in SPECIES_BSM2 else BSM2_BIO_MW

def _thod(sp):
    return SPECIES_BSM2[sp][2] if sp in SPECIES_BSM2 else BSM2_BIO_THOD

def _atoms(sp):
    return SPECIES_BSM2[sp][0] if sp in SPECIES_BSM2 else BSM2_BIO_ATOMS

# Number of C and N atoms per mol
def _nC(sp):
    return _atoms(sp).get("C", 0)

def _nN(sp):
    return _atoms(sp).get("N", 0)


# ════════════════════════════════════════════════════════════════════════
#  BSM2 organism names (7 populations)
# ════════════════════════════════════════════════════════════════════════

ORG_BSM2 = {
    "X_su": "X_su", "X_aa": "X_aa", "X_fa": "X_fa",
    "X_c4": "X_c4", "X_pro": "X_pro", "X_ac": "X_ac", "X_h2": "X_h2",
}
BSM2_ORGS = list(ORG_BSM2.values())


# ════════════════════════════════════════════════════════════════════════
#  BSM2 default parameters (Rosén & Jeppsson 2006)
# ════════════════════════════════════════════════════════════════════════

# --- Stoichiometric (COD-fraction product distributions) ---
BSM2_STOICH = dict(
    # Disintegration fractions
    f_sI_xc=0.10, f_xI_xc=0.20, f_ch_xc=0.20, f_pr_xc=0.20, f_li_xc=0.30,
    # Lipid hydrolysis
    f_fa_li=0.95,
    # Sugar uptake product fractions (COD basis)
    f_h2_su=0.19, f_bu_su=0.13, f_pro_su=0.27, f_ac_su=0.41,
    # Amino acid uptake product fractions
    f_h2_aa=0.06, f_va_aa=0.23, f_bu_aa=0.26, f_pro_aa=0.05, f_ac_aa=0.40,
    # LCFA oxidation
    f_ac_fa=0.70,   # rest is H₂ (implicitly 0.30)
    # Valerate oxidation
    f_pro_va=0.54, f_ac_va=0.31,  # rest is H₂ (0.15)
    # Butyrate oxidation
    f_ac_bu=0.80,   # rest is H₂ (0.20)
    # Propionate oxidation
    f_ac_pro=0.57,  # rest is H₂ (0.43)
    # Yields (kgCOD_X / kgCOD_S)
    Y_su=0.10, Y_aa=0.08, Y_fa=0.06, Y_c4=0.06, Y_pro=0.04, Y_ac=0.05, Y_h2=0.06,
)

# --- Kinetic parameters (BSM2 defaults) ---
BSM2_KINETICS = dict(
    # Disintegration + hydrolysis (d⁻¹)
    k_dis=0.5, k_hyd_ch=10.0, k_hyd_pr=10.0, k_hyd_li=10.0,
    # Uptake rates (d⁻¹, kgCOD_S per kgCOD_X per day)
    k_m_su=30.0, k_m_aa=50.0, k_m_fa=6.0, k_m_c4=20.0,
    k_m_pro=13.0, k_m_ac=8.0, k_m_h2=35.0,
    # Half-saturation (kgCOD/m³ = gCOD/L)
    K_S_su=0.5, K_S_aa=0.3, K_S_fa=0.4, K_S_c4=0.2,
    K_S_pro=0.1, K_S_ac=0.15, K_S_h2=7e-6,
    # Inhibition
    K_I_h2_fa=5e-6, K_I_h2_c4=1e-5, K_I_h2_pro=3.5e-6,  # kgCOD/m³
    K_I_nh3=0.0018,   # mol/L (M)
    K_S_NH3=1e-4,       # mol/L (M) — inorganic nitrogen limitation
    # Decay (d⁻¹)
    k_dec=0.02,
    # pH inhibition
    pH_UL_aa=5.5, pH_LL_aa=4.0,
    pH_UL_ac=7.0, pH_LL_ac=6.0,
    pH_UL_h2=6.0, pH_LL_h2=5.0,
    # Gas transfer
    k_L_a=200.0,  # d⁻¹
)


# ════════════════════════════════════════════════════════════════════════
#  Stoichiometry: COD-fraction → molar conversion
# ════════════════════════════════════════════════════════════════════════

def _petersen_to_molar(substrate, products_cod, Y, organism):
    """Convert BSM2 COD-fraction product distribution to molar stoichiometry.

    Parameters
    ----------
    substrate : str
        Substrate species name (e.g., "S_su").
    products_cod : dict
        ``{species: cod_fraction}`` — COD fraction of (1-Y) going to each product.
    Y : float
        Biomass yield (kgCOD_X / kgCOD_S).
    organism : str
        Organism species name.

    Returns
    -------
    tuple
        (coeffs, elements) where coeffs is list of (species, coeff) tuples.
    """
    ThOD_sub = _thod(substrate)
    ThOD_X = BSM2_BIO_THOD

    coeffs = [(substrate, -1.0)]

    # Products from catabolism
    C_out = 0.0
    N_out = 0.0
    H_out = 0.0
    O_out = 0.0

    for sp, f_cod in products_cod.items():
        ThOD_p = _thod(sp)
        nu = (1 - Y) * f_cod * ThOD_sub / ThOD_p
        if abs(nu) > 1e-12:
            coeffs.append((sp, +nu))
            C_out += nu * _nC(sp)
            N_out += nu * _nN(sp)
            at = _atoms(sp)
            H_out += nu * at.get("H", 0)
            O_out += nu * at.get("O", 0)

    # Biomass
    nu_X = Y * ThOD_sub / ThOD_X
    coeffs.append((organism, +nu_X))
    C_out += nu_X * _nC(organism)
    N_out += nu_X * _nN(organism)
    H_out += nu_X * _atoms(organism).get("H", 0)
    O_out += nu_X * _atoms(organism).get("O", 0)

    # Carbon balance → molecular CO2(aq) (state-unification C3:
    # canonical molecular form; engine redistributes between CO2 /
    # HCO3- / CO3-- each speciation solve, preserving total).
    C_in = _nC(substrate)
    nu_CO2 = C_in - C_out
    if abs(nu_CO2) > 1e-10:
        coeffs.append(("CO2", +nu_CO2))
        O_out += nu_CO2 * 2
        # H from CO2(aq) = 0

    # Nitrogen balance → NH₃
    N_in = _nN(substrate)
    nu_NH3 = N_in - N_out
    if abs(nu_NH3) > 1e-10:
        coeffs.append(("NH3", +nu_NH3))
        H_out += nu_NH3 * 3

    # Hydrogen balance → H₂O
    H_in = _atoms(substrate).get("H", 0)
    nu_H2O = (H_in - H_out) / 2.0
    if abs(nu_H2O) > 1e-10:
        if nu_H2O > 0:
            coeffs.append(("H2O", +nu_H2O))
        else:
            coeffs.append(("H2O", nu_H2O))

    return tuple(coeffs)


def _stoich_disintegration(p):
    """X_xc → products (first-order, BSM2 fractions)."""
    ThOD_xc = _thod("X_xc")
    nu_ch = p["f_ch_xc"] * ThOD_xc / _thod("X_ch")
    nu_pr = p["f_pr_xc"] * ThOD_xc / _thod("X_pr")
    nu_li = p["f_li_xc"] * ThOD_xc / _thod("X_li")
    nu_xI = p["f_xI_xc"] * ThOD_xc / _thod("X_I")
    nu_sI = p["f_sI_xc"] * ThOD_xc / _thod("S_I")
    return (("X_xc", -1), ("X_ch", +nu_ch), ("X_pr", +nu_pr),
            ("X_li", +nu_li), ("X_I", +nu_xI), ("S_I", +nu_sI))


def _stoich_hyd_ch():
    """X_ch → S_su (first-order hydrolysis)."""
    nu = _thod("X_ch") / _thod("S_su")
    return (("X_ch", -1), ("S_su", +nu))


def _stoich_hyd_pr():
    """X_pr → S_aa (first-order hydrolysis)."""
    nu = _thod("X_pr") / _thod("S_aa")
    return (("X_pr", -1), ("S_aa", +nu))


def _stoich_hyd_li(p):
    """X_li → S_su + S_fa (first-order hydrolysis)."""
    ThOD_li = _thod("X_li")
    nu_su = (1 - p["f_fa_li"]) * ThOD_li / _thod("S_su")
    nu_fa = p["f_fa_li"] * ThOD_li / _thod("S_fa")
    return (("X_li", -1), ("S_su", +nu_su), ("S_fa", +nu_fa))


def _stoich_decay(org):
    """Biomass → X_xc (first-order decay)."""
    nu = BSM2_BIO_THOD / _thod("X_xc")
    return ((org, -1), ("X_xc", +nu))


# ════════════════════════════════════════════════════════════════════════
#  Rate function factories
# ════════════════════════════════════════════════════════════════════════

_R_J = 8.31446

def _pKa_NH4(T_K):
    Ka_ref = 10.0 ** (-9.25)
    Ka_T = Ka_ref * math.exp(-(51965.0 / _R_J) * (1.0 / T_K - 1.0 / 298.15))
    return -math.log10(max(Ka_T, 1e-30))


def _build_rate_functions(kp):
    """Build all BSM2-canonical rate function closures."""

    # pH inhibition (BSM2 formulation)
    K_pH_aa = 10 ** (-(kp["pH_UL_aa"] + kp["pH_LL_aa"]) / 2.0)
    n_aa = 3.0 / (kp["pH_UL_aa"] - kp["pH_LL_aa"])
    K_pH_ac = 10 ** (-(kp["pH_UL_ac"] + kp["pH_LL_ac"]) / 2.0)
    n_ac = 3.0 / (kp["pH_UL_ac"] - kp["pH_LL_ac"])
    K_pH_h2 = 10 ** (-(kp["pH_UL_h2"] + kp["pH_LL_h2"]) / 2.0)
    n_h2 = 3.0 / (kp["pH_UL_h2"] - kp["pH_LL_h2"])

    K_S_NH3 = kp["K_S_NH3"]
    K_I_nh3 = kp["K_I_nh3"]

    def _get_pH(env):
        return float(env.pH) if env.has_pH else 7.0

    def _get_H(env):
        # H+ is a speciation-output species, not a phase species; read
        # it via the pH channel rather than from concentrations. (B2.)
        if env.has_pH:
            return 10.0 ** (-float(env.pH))
        return 1e-7

    def _I_pH_aa(env):
        H = _get_H(env)
        return K_pH_aa ** n_aa / (H ** n_aa + K_pH_aa ** n_aa)

    def _I_pH_ac(env):
        H = _get_H(env)
        return K_pH_ac ** n_ac / (H ** n_ac + K_pH_ac ** n_ac)

    def _I_pH_h2(env):
        H = _get_H(env)
        return K_pH_h2 ** n_h2 / (H ** n_h2 + K_pH_h2 ** n_h2)

    def _I_IN(env):
        # S_IN stored as mol NH₃-equivalent; 1 N per molecule → mol/L = mol_N/L
        S_IN = env.concentrations.get("NH3", 0)
        return 1.0 / (1.0 + K_S_NH3 / max(S_IN, 1e-12))

    def _I_nh3(env):
        # S_IN in mol/L = mol_N/L
        S_IN_molL = env.concentrations.get("NH3", 0)
        pH = _get_pH(env)
        pKa = _pKa_NH4(env.T_K)
        f_nh3 = 1.0 / (1.0 + 10.0 ** (pKa - pH))
        S_nh3 = S_IN_molL * f_nh3
        return 1.0 / (1.0 + S_nh3 / K_I_nh3)

    def _I_h2(env, K_I):
        S_h2_cod = env.concentrations.get("S_h2", 0) * _thod("S_h2")
        return 1.0 / (1.0 + S_h2_cod / K_I)

    # ── First-order rate factories ────────────────────────────────────

    def rate_first_order(sp, k_d):
        """Rate = k × C_sp × V [mol/h]."""
        k_h = k_d / 24.0
        def r(env):
            C = env.concentrations.get(sp, 0)
            return k_h * C * env.V_L if C > 0 else 0.0
        return r

    # ── Monod uptake rate factories (BSM2 formulation) ────────────────

    def rate_uptake(org, sub, k_m_d, K_S_cod, I_fn):
        """BSM2 Monod: k_m × S_cod/(K_S + S_cod) × X_cod × I.

        Returns mol_substrate / h.
        """
        ThOD_sub = _thod(sub)
        ThOD_X = BSM2_BIO_THOD
        k_m_h = k_m_d / 24.0
        def r(env):
            S_cod = env.concentrations.get(sub, 0) * ThOD_sub
            X_cod = env.concentrations.get(org, 0) * ThOD_X
            if S_cod <= 0 or X_cod <= 0:
                return 0.0
            I = I_fn(env)
            Rho = k_m_h * S_cod / (K_S_cod + S_cod) * X_cod * I
            return Rho / ThOD_sub * env.V_L
        return r

    def rate_c4_competitive(org, sub, other_sub, k_m_d, K_S_cod, I_fn):
        """BSM2 competitive C4 uptake: adds S_i/(S_bu+S_va) factor."""
        ThOD_sub = _thod(sub)
        ThOD_other = _thod(other_sub)
        ThOD_X = BSM2_BIO_THOD
        k_m_h = k_m_d / 24.0
        def r(env):
            S_cod = env.concentrations.get(sub, 0) * ThOD_sub
            S_other_cod = env.concentrations.get(other_sub, 0) * ThOD_other
            X_cod = env.concentrations.get(org, 0) * ThOD_X
            if S_cod <= 0 or X_cod <= 0:
                return 0.0
            I = I_fn(env)
            alloc = S_cod / (S_cod + S_other_cod + 1e-12)
            Rho = k_m_h * S_cod / (K_S_cod + S_cod) * X_cod * alloc * I
            return Rho / ThOD_sub * env.V_L
        return r

    # ── Inhibition combinations (BSM2) ────────────────────────────────

    def I_5(env): return _I_pH_aa(env) * _I_IN(env)                      # su, aa
    def I_7(env): return _I_pH_aa(env) * _I_IN(env) * _I_h2(env, kp["K_I_h2_fa"])  # fa
    def I_8(env): return _I_pH_aa(env) * _I_IN(env) * _I_h2(env, kp["K_I_h2_c4"])  # c4
    def I_10(env): return _I_pH_aa(env) * _I_IN(env) * _I_h2(env, kp["K_I_h2_pro"]) # pro
    def I_11(env): return _I_pH_ac(env) * _I_IN(env) * _I_nh3(env)       # ac
    def I_12(env): return _I_pH_h2(env) * _I_IN(env)                      # h2

    return {
        "first_order": rate_first_order,
        "uptake": rate_uptake,
        "c4_competitive": rate_c4_competitive,
        "I_5": I_5, "I_7": I_7, "I_8": I_8,
        "I_10": I_10, "I_11": I_11, "I_12": I_12,
    }


# ════════════════════════════════════════════════════════════════════════
#  Entry builder
# ════════════════════════════════════════════════════════════════════════

def _make_entries(coeff_tuples):
    return [
        StoichiometryEntry(species=SPECIES[sp], phase="liquid", coefficient=coeff)
        for sp, coeff in coeff_tuples
    ]


# ════════════════════════════════════════════════════════════════════════
#  Equilibrium reactions (Phase 1 — chemistry-unification-1)
# ════════════════════════════════════════════════════════════════════════

# pKa / Van 't Hoff dH values per Rosen & Jeppsson (2006), matching
# EquilibriumSet.bsm2_default(). The pKa references are at 25°C (298.15 K);
# the engine applies Van 't Hoff to reach the operating T at solve time.
_BSM2_PKW   = 14.0
_BSM2_DH_W  = 55900.0
_BSM2_PKA_CO2_1 = 6.35
_BSM2_DH_CO2_1  = 7646.0
_BSM2_PKA_NH4   = 9.25
_BSM2_DH_NH4    = 51965.0
_BSM2_PKA_VFA = {"S_ac": 4.76, "S_pro": 4.88, "S_bu": 4.82, "S_va": 4.86}


def _build_bsm2_equilibrium_reactions() -> List[EquilibriumReaction]:
    """Equilibrium reactions for BSM2 chemistry: water, carbonate (1st step),
    ammonium, and four VFAs. Returned in insertion order; this order
    seeds the EquilibriumSet that ``BisectionChemicalEquilibriumEngine.from_reactions``
    builds.

    The VFA reactions use a synthetic ``<sp>_anion`` Species for the
    conjugate base — it appears in the stoichiometry for elemental
    bookkeeping only, and is **not** tracked in any phase's
    ``n_mol``. BSM2 keeps the total VFA under the protonated species
    id (``"S_ac"``) and the engine computes the dissociation split
    internally from the pKa.
    """
    rxns: List[EquilibriumReaction] = []

    # Water: H2O ⇌ H+ + OH-
    rxns.append(EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=SPECIES["H2O"], phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=H_plus,        phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=OH_minus,      phase="liquid", coefficient=+1.0),
        ],
        log_K=-_BSM2_PKW, dH_J_per_mol=_BSM2_DH_W,
        balance_elements=("H", "O"),
        label="eq_water",
    ))

    # CO2: CO2 + H2O ⇌ HCO3- + H+ (first dissociation, n_active=1 in
    # BSM2). total_id="CO2" labels the carbonate ladder so the engine
    # resolves Species IDs (CO2, HCO3-, CO3--) for the
    # speciation solve and writeback.
    rxns.append(EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=SPECIES["CO2"], phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=SPECIES["H2O"],   phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=HCO3_minus,       phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=H_plus,           phase="liquid", coefficient=+1.0),
        ],
        log_K=-_BSM2_PKA_CO2_1, dH_J_per_mol=_BSM2_DH_CO2_1,
        total_id="CO2",
        balance_elements=("C", "H", "O"),
        label="eq_CO2",
    ))

    # NH4+ ⇌ NH3 + H+ — total tracker is n_mol["NH3"] (BSM2 convention).
    rxns.append(EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=NH4_plus,       phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=SPECIES["NH3"], phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=H_plus,         phase="liquid", coefficient=+1.0),
        ],
        log_K=-_BSM2_PKA_NH4, dH_J_per_mol=_BSM2_DH_NH4,
        total_id="NH3",
        balance_elements=("N", "H"),
        label="eq_NH4",
    ))

    # VFAs: each as a monoprotic acid HA ⇌ A- + H+, with the total
    # tracked under the protonated species id (BSM2 convention).
    for vfa_id, pKa in _BSM2_PKA_VFA.items():
        ha = SPECIES[vfa_id]
        a_minus_atoms = dict(ha.atoms)
        # Remove one H to form the conjugate base.
        a_minus_atoms["H"] = int(a_minus_atoms.get("H", 0)) - 1
        a_minus = Species(
            id=f"{vfa_id}-",
            atoms=a_minus_atoms,
            charge=-1,
            MW=float(ha.MW) - 1.008,
        )
        rxns.append(EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=ha,      phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=a_minus, phase="liquid", coefficient=+1.0),
                StoichiometryEntry(species=H_plus,  phase="liquid", coefficient=+1.0),
            ],
            log_K=-pKa,
            balance_elements=("C", "H", "O"),
            label=f"eq_{vfa_id}",
        ))

    # Cross-phase partition: CO2 gas ⇌ CO2 liquid. No log_K — Henry's
    # law constant lives on the link. This declaration auto-wires
    # speciation_keys["CO2"] = "CO2" and builds the carbonate ladder
    # from the single-phase acid-base reactions above (chemistry-
    # unification-3b, C7).
    rxns.append(EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=CO2_sp,  phase="gas",    coefficient=-1.0),
            StoichiometryEntry(species=CO2_sp,  phase="liquid", coefficient=+1.0),
        ],
        balance_elements=("C", "O"),
        label="partition_CO2",
    ))

    return rxns


# ════════════════════════════════════════════════════════════════════════
#  Public API
# ════════════════════════════════════════════════════════════════════════

def build_bsm2_reactions(
    *,
    stoich: Dict = None,
    kinetics: Dict = None,
    verbose: bool = True,
) -> ReactionSystem:
    """Build a BSM2-canonical ADM1 reaction system in molar units.

    Returns 19 biochemical + 7 decay = 26 reactions total:
    - 1 disintegration (X_xc → X_ch + X_pr + X_li + X_I + S_I)
    - 3 hydrolysis (X_ch→S_su, X_pr→S_aa, X_li→S_su+S_fa)
    - 8 uptake (Monod kinetics with BSM2 inhibition combinations)
    - 7 decay (first-order back to X_xc)

    Parameters
    ----------
    stoich : dict, optional
        Override BSM2 stoichiometric parameters (merged with defaults).
    kinetics : dict, optional
        Override BSM2 kinetic parameters (merged with defaults).
    verbose : bool
        Print reaction summary.
    """
    p = dict(BSM2_STOICH)
    if stoich:
        p.update(stoich)

    kp = dict(BSM2_KINETICS)
    if kinetics:
        kp.update(kinetics)

    rf = _build_rate_functions(kp)

    reactions = []

    # ── Disintegration ────────────────────────────────────────────────
    coeffs = _stoich_disintegration(p)
    entries = _make_entries(coeffs)
    rxn = KineticReaction(stoichiometry=entries,
                   rate_fn=rf["first_order"]("X_xc", kp["k_dis"]),
                   balance_elements=(), label="R1_disintegration")
    reactions.append(rxn)

    # ── Hydrolysis ────────────────────────────────────────────────────
    for label, coeffs, sp, k in [
        ("R2_hyd_ch", _stoich_hyd_ch(), "X_ch", kp["k_hyd_ch"]),
        ("R3_hyd_pr", _stoich_hyd_pr(), "X_pr", kp["k_hyd_pr"]),
        ("R4_hyd_li", _stoich_hyd_li(p), "X_li", kp["k_hyd_li"]),
    ]:
        entries = _make_entries(coeffs)
        rxn = KineticReaction(stoichiometry=entries,
                       rate_fn=rf["first_order"](sp, k),
                       balance_elements=(), label=label)
        reactions.append(rxn)

    # ── Uptake reactions ──────────────────────────────────────────────

    # R5: Sugar uptake
    coeffs = _petersen_to_molar("S_su", {
        "S_bu": p["f_bu_su"], "S_pro": p["f_pro_su"],
        "S_ac": p["f_ac_su"], "S_h2": p["f_h2_su"],
    }, p["Y_su"], "X_su")
    entries = _make_entries(coeffs)
    reactions.append(KineticReaction(stoichiometry=entries,
        rate_fn=rf["uptake"]("X_su", "S_su", kp["k_m_su"], kp["K_S_su"], rf["I_5"]),
        balance_elements=(), label="R5_uptake_su"))

    # R6: Amino acid uptake
    coeffs = _petersen_to_molar("S_aa", {
        "S_va": p["f_va_aa"], "S_bu": p["f_bu_aa"], "S_pro": p["f_pro_aa"],
        "S_ac": p["f_ac_aa"], "S_h2": p["f_h2_aa"],
    }, p["Y_aa"], "X_aa")
    entries = _make_entries(coeffs)
    reactions.append(KineticReaction(stoichiometry=entries,
        rate_fn=rf["uptake"]("X_aa", "S_aa", kp["k_m_aa"], kp["K_S_aa"], rf["I_5"]),
        balance_elements=(), label="R6_uptake_aa"))

    # R7: LCFA uptake
    f_h2_fa = 1.0 - p["f_ac_fa"]
    coeffs = _petersen_to_molar("S_fa", {
        "S_ac": p["f_ac_fa"], "S_h2": f_h2_fa,
    }, p["Y_fa"], "X_fa")
    entries = _make_entries(coeffs)
    reactions.append(KineticReaction(stoichiometry=entries,
        rate_fn=rf["uptake"]("X_fa", "S_fa", kp["k_m_fa"], kp["K_S_fa"], rf["I_7"]),
        balance_elements=(), label="R7_uptake_fa"))

    # R8: Valerate uptake (competitive C4)
    f_h2_va = 1.0 - p["f_pro_va"] - p["f_ac_va"]
    coeffs = _petersen_to_molar("S_va", {
        "S_pro": p["f_pro_va"], "S_ac": p["f_ac_va"], "S_h2": f_h2_va,
    }, p["Y_c4"], "X_c4")
    entries = _make_entries(coeffs)
    reactions.append(KineticReaction(stoichiometry=entries,
        rate_fn=rf["c4_competitive"]("X_c4", "S_va", "S_bu",
                                      kp["k_m_c4"], kp["K_S_c4"], rf["I_8"]),
        balance_elements=(), label="R8_uptake_va"))

    # R9: Butyrate uptake (competitive C4)
    f_h2_bu = 1.0 - p["f_ac_bu"]
    coeffs = _petersen_to_molar("S_bu", {
        "S_ac": p["f_ac_bu"], "S_h2": f_h2_bu,
    }, p["Y_c4"], "X_c4")
    entries = _make_entries(coeffs)
    reactions.append(KineticReaction(stoichiometry=entries,
        rate_fn=rf["c4_competitive"]("X_c4", "S_bu", "S_va",
                                      kp["k_m_c4"], kp["K_S_c4"], rf["I_8"]),
        balance_elements=(), label="R9_uptake_bu"))

    # R10: Propionate uptake
    f_h2_pro = 1.0 - p["f_ac_pro"]
    coeffs = _petersen_to_molar("S_pro", {
        "S_ac": p["f_ac_pro"], "S_h2": f_h2_pro,
    }, p["Y_pro"], "X_pro")
    entries = _make_entries(coeffs)
    reactions.append(KineticReaction(stoichiometry=entries,
        rate_fn=rf["uptake"]("X_pro", "S_pro", kp["k_m_pro"], kp["K_S_pro"], rf["I_10"]),
        balance_elements=(), label="R10_uptake_pro"))

    # R11: Acetate uptake (methanogenesis)
    coeffs = _petersen_to_molar("S_ac", {
        "S_ch4": 1.0,
    }, p["Y_ac"], "X_ac")
    entries = _make_entries(coeffs)
    reactions.append(KineticReaction(stoichiometry=entries,
        rate_fn=rf["uptake"]("X_ac", "S_ac", kp["k_m_ac"], kp["K_S_ac"], rf["I_11"]),
        balance_elements=(), label="R11_uptake_ac"))

    # R12: H₂ uptake (methanogenesis)
    coeffs = _petersen_to_molar("S_h2", {
        "S_ch4": 1.0,
    }, p["Y_h2"], "X_h2")
    entries = _make_entries(coeffs)
    reactions.append(KineticReaction(stoichiometry=entries,
        rate_fn=rf["uptake"]("X_h2", "S_h2", kp["k_m_h2"], kp["K_S_h2"], rf["I_12"]),
        balance_elements=(), label="R12_uptake_h2"))

    # ── Decay ─────────────────────────────────────────────────────────
    k_dec_h = kp["k_dec"] / 24.0
    for org in BSM2_ORGS:
        coeffs = _stoich_decay(org)
        entries = _make_entries(coeffs)
        rxn = KineticReaction(stoichiometry=entries,
                       rate_fn=rf["first_order"](org, kp["k_dec"]),
                       balance_elements=(), label=f"R_decay_{org}")
        reactions.append(rxn)

    # ── Equilibrium reactions (chemistry-unification-1) ──────────────
    # BSM2 pKa values at 25°C reference; Van 't Hoff dH carries the
    # 35°C correction at solve time. Matches EquilibriumSet.bsm2_default
    # bit-for-bit. CO2 is monoprotic (n_active=1 in legacy form) — the
    # second dissociation (pKa2≈10.33) is not part of the charge
    # balance. The molecular CO2 fraction used for gas transfer is read
    # from PropertyResult.alphas["CO2aq"] (chemistry-unification-2).
    reactions.extend(_build_bsm2_equilibrium_reactions())

    _log = logger.info if verbose else logger.debug
    n_bio = 12  # 1 dis + 3 hyd + 8 uptake
    _log("BSM2-canonical ADM1: %d reactions (%d biochemical + %d decay + %d equilibrium)",
         len(reactions), n_bio, len(BSM2_ORGS), 7)
    for rxn in reactions:
        coeffs_str = []
        for e in rxn.stoichiometry:
            if e.coefficient < 0:
                coeffs_str.append(f"{abs(e.coefficient):.3f} {e.species.id}")
        products_str = []
        for e in rxn.stoichiometry:
            if e.coefficient > 0:
                products_str.append(f"{e.coefficient:.3f} {e.species.id}")
        _log("  %s: %s → %s", rxn.label,
             " + ".join(coeffs_str), " + ".join(products_str))

    return ReactionSystem(reactions, label="ADM1_BSM2")


def build_bsm2_cv(
    reaction_system: ReactionSystem,
    *,
    V_total_L: float = 3700.0,
    headspace_frac: float = 300.0 / 3700.0,
    T_K: float = 308.15,
    k_L_a_per_d: float = 200.0,
    activity_model: str = "davies",
    thermo: "ThermodynamicConfig" = None,
):
    """Build a BSM2-canonical fermenter ControlVolume with gas transfer.

    Uses kLa-based kinetic transfer (not equilibrium) with BSM2
    Henry constants for H₂, CH₄, and CO₂.

    Parameters
    ----------
    reaction_system : ReactionSystem
        BSM2 reaction model.
    V_total_L : float
        Total vessel volume (L).
    headspace_frac : float
        Gas volume fraction of total.
    T_K : float
        Operating temperature (K).
    k_L_a_per_d : float
        Volumetric mass transfer coefficient (1/d).
    activity_model : str
        Activity model for speciation (``"ideal"`` or ``"davies"``).
    thermo : ThermodynamicConfig, optional
        Single source of truth for all thermodynamic parameters (pKa,
        dH, n_active).  When provided, speciation corrections on the
        transfer link are populated exclusively from this config via
        ``apply_to_cv()``.  **Recommended for all new usage.**

        When not provided, a deprecation warning is emitted and a
        default ``ThermodynamicConfig`` is created automatically.  This
        fallback exists for backward compatibility but will be removed
        in a future version.

    Returns
    -------
    ControlVolume
    """
    from vlmodels.fermenter.config.builder import FermenterBuilder
    from vlmodels.fermenter.config.configs import TransferConfig
    from PyOMES.core.boundaries import PressureReliefVent
    from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine

    use_activity = activity_model.lower() != "ideal"
    kLa_h = k_L_a_per_d / 24.0

    # BSM2 Henry constants at 35°C converted to mol/L/atm
    H_co2 = 0.02751
    H_ch4 = 0.001177
    H_h2 = 0.000748

    # Start with empty species (no default O2/N2/CO2 transfer).
    # speciation_keys is no longer passed — the CO2 cross-phase partition
    # reaction in build_bsm2_reactions() drives derive_speciation_keys()
    # at CV construction (chemistry-unification-3b C7).
    tc = TransferConfig(species={})

    # The ReactionSystem has pre-bucketed kinetic vs equilibrium
    # reactions at construction. The whole system is attached to the
    # CV (kinetic + black-box drive compute_rates; equilibria are
    # silently skipped by ReactionSystem.compute_rates). Equilibria
    # also seed the speciation engine below via from_reactions().
    equilibrium_rxns = (
        list(reaction_system.single_phase_equilibria)
        + list(reaction_system.cross_phase_equilibria)
    )

    cv = (FermenterBuilder()
          .vessel(V_total_L=V_total_L, headspace_frac=headspace_frac,
                  T_K=T_K, yO2_init=0.0, yCO2_init=0.0)
          .no_gas_feed()
          .transfer(tc)
          .transfer_species("S_ch4", mode="kinetic", kLa_per_h=kLa_h,
                            henry_mol_L_atm=H_ch4)
          .transfer_species("S_h2", mode="kinetic", kLa_per_h=kLa_h,
                            henry_mol_L_atm=H_h2)
          .transfer_species("CO2", mode="equilibrium",
                            henry_mol_L_atm=H_co2)
          .chemistry(use_activity=use_activity,
                     activity_model=activity_model)
          .reaction_system(reaction_system)
          .label("ADM1_BSM2")
          .build())

    # Attach the BSM2-configured speciation engine to the
    # ReactionSystem. state-unification C4d: cv.reaction_system.engine
    # is the canonical engine attachment; the legacy
    # SpeciationPropertySolver wrapper is gone.
    engine = BisectionChemicalEquilibriumEngine.from_reactions(
        equilibrium_rxns,
        activity_model=activity_model,
        use_activity=use_activity,
        T_K=T_K,
    )
    cv.reaction_system.attach_engine(engine)

    cv.boundaries.append(PressureReliefVent(P_set_atm=1.013/1.01325, mode="instant"))

    # Attach ThermodynamicConfig as a back-reference for chem_env_fn's
    # temperature-tracking accessors (T_K_live).  Post-Phase-2 the
    # gas-liquid link no longer reads pKa data from this config; the
    # molecular fraction for CO2 transfer is read from
    # PropertyResult.alphas["CO2aq"], populated by the speciation
    # engine using the equilibrium reactions declared in
    # build_bsm2_reactions(). Phase 3 will absorb ThermodynamicConfig
    # into ChemistryDatabase and the back-reference can move there.
    if thermo is not None:
        thermo.apply_to_cv(cv)
    else:
        from PyOMES.chemistry.thermo_params import ThermodynamicConfig
        ThermodynamicConfig.bsm2_default().apply_to_cv(cv)

    return cv


def seed_bsm2_strong_ions(
    cv,
    *,
    CT_cation: float = 0.0,
    CT_anion: float = 0.00521,
) -> None:
    """Seed BSM2 strong-ion totals into the CV's liquid n_mol.

    state-unification C4: strong ions live as Species in
    ``phase.n_mol`` (PHREEQC convention; see Option C of the C4
    sub-decisions). Replaces the deleted
    ``make_bsm2_chem_env_fn`` which returned a chem_env dict
    threading ``strong_kwargs`` through ``cv.advance``.

    Parameters
    ----------
    cv : ControlVolume
        The BSM2 CV to seed.
    CT_cation, CT_anion : float
        Monovalent-equivalent inert strong-ion concentrations
        (mol/L). Default 0.0 / 0.00521 (BSM2 values). Stored under
        ``cv.phases["liquid"].n_mol["S_cat"]`` /
        ``["S_an"]`` (charge ±1 lump species, atoms={}).
    """
    liq = cv.phases.get("liquid")
    if liq is None:
        return
    V_L = float(getattr(liq, "V_L", 1.0))
    if CT_cation:
        liq.n_mol["S_cat"] = float(CT_cation) * V_L
    if CT_anion:
        liq.n_mol["S_an"] = float(CT_anion) * V_L
