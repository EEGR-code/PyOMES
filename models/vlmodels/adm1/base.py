# -*- coding: utf-8 -*-
"""ADM1 (Anaerobic Digestion Model No. 1) reaction model.

Provides a complete ADM1 implementation as a reusable module with
optional extensions that can be enabled via keyword arguments.

Quick start (new framework, simulation-class C13)
-------------------------------------------------
>>> from vlmodels.adm1 import build_adm1_reactions, build_adm1_cv
>>> from vlmodels.adm1 import seed_adm1_strong_ions
>>> from PyOMES.core import Simulation
>>>
>>> rxn_set = build_adm1_reactions()
>>> cv = build_adm1_cv(rxn_set, V_total_L=2.0, T_K=308.15)
>>> seed_adm1_strong_ions(cv, CT_Na=0.100, CT_Cl=0.030)
>>> sim = Simulation(cvs={"main": cv})
>>> result = sim.run(tau_h=1440, n_steps=14400)

Notes
-----
Post-STATE_UNIFICATION C4: ``chem_env`` and ``chem_env_fn`` are
gone; strong ions live in ``cv.phases["liquid"].n_mol`` as Species
entries (seeded via ``seed_adm1_strong_ions``). The legacy
``make_adm1_callback`` / ``run_batch`` flow is no longer required —
``Simulation`` owns the time loop and routes ``cv.advance(dt_h,
t_h)`` calls directly.

Extensions
----------
- ``sulfate_reduction=True``: ADM1-S sulfate-reducing bacteria (Stage 11)
- ``soluble_inerts=True``: S_I production from uptake reactions (Stage 10)
- ``h2s_inhibition=True``: H₂S inhibition on methanogens and SRB
  (auto-enabled when ``sulfate_reduction=True``)
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
from PyOMES.chemistry.common_species import CO2 as _CO2_sp, NH3 as _NH3_sp, H2O as _H2O_sp

logger = logging.getLogger(__name__)


# ════════════════════════════════════════════════════════════════════════
#  Species definitions
# ════════════════════════════════════════════════════════════════════════

# {name: (atoms_dict, MW)}
SPECIES = {
    "Starch":      ({"C":6, "H":10, "O":5},                162.14),
    "Glucose":     ({"C":6, "H":12, "O":6},                180.156),
    "Protein":     ({"C":5, "H":9,  "O":2, "N":1},         115.13),
    "AminoAcid":   ({"C":5, "H":11, "O":3, "N":1},         133.15),
    "Lipid":       ({"C":51,"H":98, "O":6},                 807.32),
    "Glycerol":    ({"C":3, "H":8,  "O":3},                 92.094),
    "LCFA":        ({"C":16,"H":32, "O":2},                 256.424),
    "AceticAcid":  ({"C":2, "H":4,  "O":2},                 60.052),
    "Propionate":  ({"C":3, "H":6,  "O":2},                 74.079),
    "Butyrate":    ({"C":4, "H":8,  "O":2},                 88.106),
    "Valerate":    ({"C":5, "H":10, "O":2},                 102.132),
    "CO2":         ({"C":1, "O":2},                          44.009),
    "CH4":         ({"C":1, "H":4},                          16.043),
    "H2":          ({"H":2},                                  2.016),
    "H2O":         ({"H":2, "O":1},                          18.015),
    "NH3":         ({"N":1, "H":3},                          17.031),
    # Decay/disintegration pools (same formula as biomass)
    "X_composite": ({"C":1, "H":1.8, "O":0.5, "N":0.2},    24.6),
    "X_I":         ({"C":1, "H":1.8, "O":0.5, "N":0.2},    24.6),
    "S_I":         ({"C":1, "H":1.8, "O":0.5, "N":0.2},    24.6),
    # Sulfur species (Stage 11)
    "SO4":         ({"S":1, "O":4},                          96.06),
    "H2S":         ({"S":1, "H":2},                          34.08),
    # Ethanol extension
    "Ethanol":     ({"C":2, "H":6, "O":1},                  46.068),
}

BIO_MW = 24.6
BIO_CHO = {"C":1, "H":1.8, "O":0.5}           # CHO organisms
BIO_CHON = {"C":1, "H":1.8, "O":0.5, "N":0.2} # CHON organisms

def _mw(sp):
    if sp in SPECIES:
        return SPECIES[sp][1]
    return BIO_MW

def _at(sp):
    if sp in SPECIES:
        return SPECIES[sp][0]
    return BIO_CHO


# ════════════════════════════════════════════════════════════════════════
#  Organism definitions
# ════════════════════════════════════════════════════════════════════════

# Organism short-keys → species names
ORG = {
    "hydro":    "Hydrolyzer",
    "proteo":   "Proteolyzer",
    "lipo":     "Lipolyzer",
    "acid_c":   "Acidogen_C",
    "acid_p":   "Acidogen_P",
    "gly_f":    "GlyFermenter",
    "syn_lcfa": "Syn_LCFA",
    "syn_prop": "Syn_Prop",
    "syn_but":  "Syn_But",
    "syn_val":  "Syn_Val",
    "meth_a":   "Meth_Ace",
    "meth_h2":  "Meth_H2",
    # Extension organisms (Stage 11)
    "srb_ace":  "SRB_Ace",
    "srb_h2":   "SRB_H2",
    # Extension organisms (ethanol)
    "etoh_f":   "EtOH_Ferm",    # ethanol fermenter (glucose → ethanol)
    "syn_etoh": "Syn_EtOH",     # syntrophic ethanol oxidiser (ethanol → acetate)
}

# Core ADM1 organisms (always present)
CORE_ORGS = [ORG[k] for k in [
    "hydro", "proteo", "lipo",
    "acid_c", "acid_p", "gly_f", "syn_lcfa",
    "syn_prop", "syn_but", "syn_val",
    "meth_a", "meth_h2",
]]

# SRB organisms (only when sulfate_reduction=True)
SRB_ORGS = [ORG["srb_ace"], ORG["srb_h2"]]

# Ethanol organisms (only when ethanol=True)
ETOH_ORGS = [ORG["etoh_f"], ORG["syn_etoh"]]

# CHON organisms (produce nitrogen-containing biomass)
CHON_ORGS = {ORG["proteo"], ORG["acid_p"]}


# ════════════════════════════════════════════════════════════════════════
#  Default kinetic parameters (at T_ref = 35°C)
# ════════════════════════════════════════════════════════════════════════

DEFAULT_T_REF_K = 308.15

DEFAULT_KINETICS = {
    ORG["hydro"]:    {"mu": 0.10, "Ks": 2.0,  "Y": 0.10, "sub": "Starch",
                      "Ea_R": 5000.0},
    ORG["proteo"]:   {"mu": 0.08, "Ks": 1.5,  "Y": 0.05, "sub": "Protein",
                      "Ea_R": 5000.0},
    ORG["lipo"]:     {"mu": 0.06, "Ks": 3.0,  "Y": 0.05, "sub": "Lipid",
                      "Ea_R": 5000.0},
    ORG["acid_c"]:   {"mu": 0.30, "Ks": 0.5,  "Y": 0.10, "sub": "Glucose",
                      "Ki_ace": 15.0, "Ea_R": 7000.0},
    ORG["acid_p"]:   {"mu": 0.20, "Ks": 0.8,  "Y": 0.10, "sub": "AminoAcid",
                      "Ki_ace": 12.0, "Ea_R": 7000.0},
    ORG["gly_f"]:    {"mu": 0.25, "Ks": 0.3,  "Y": 0.10, "sub": "Glycerol",
                      "Ki_ace": 15.0, "Ea_R": 7000.0},
    ORG["syn_lcfa"]: {"mu": 0.010, "Ks": 0.5, "Y": 0.02, "sub": "LCFA",
                      "Ki_H2": 0.5, "Ea_R": 8000.0},
    ORG["syn_prop"]: {"mu": 0.015, "Ks": 0.2, "Y": 0.02, "sub": "Propionate",
                      "Ki_H2": 0.5, "Ea_R": 8000.0},
    ORG["syn_but"]:  {"mu": 0.018, "Ks": 0.3, "Y": 0.02, "sub": "Butyrate",
                      "Ki_H2": 0.5, "Ea_R": 8000.0},
    ORG["syn_val"]:  {"mu": 0.012, "Ks": 0.4, "Y": 0.02, "sub": "Valerate",
                      "Ki_H2": 0.5, "Ea_R": 8000.0},
    ORG["meth_a"]:   {"mu": 0.020, "Ks": 0.3, "Y": 0.04, "sub": "AceticAcid",
                      "Ea_R": 11000.0},
    ORG["meth_h2"]:  {"mu": 0.035, "Ks": 0.01,"Y": 0.01, "sub": "H2",
                      "Y_basis": "CO2", "Ea_R": 11000.0},
    # Stage 11: SRB
    ORG["srb_ace"]:  {"mu": 0.018, "Ks_ace": 0.2, "Ks_SO4": 0.2,
                      "Y": 0.05, "sub": "AceticAcid", "Ea_R": 9000.0},
    ORG["srb_h2"]:   {"mu": 0.025, "Ks_H2": 0.01, "Ks_SO4": 0.2,
                      "Y": 0.01, "sub": "H2", "Y_basis": "SO4", "Ea_R": 9000.0},
    # Ethanol extension
    # Ethanol fermenter: fast glucose fermentation to ethanol + CO₂
    # Competes with VFA acidogenesis for glucose; dominates at low pH
    ORG["etoh_f"]:   {"mu": 0.50, "Ks": 0.5, "Y": 0.10, "sub": "Glucose",
                      "Ea_R": 7000.0},
    # Syntrophic ethanol oxidiser: ethanol → acetate + H₂
    # Thermodynamically similar to propionate oxidation, H₂-inhibited
    ORG["syn_etoh"]: {"mu": 0.020, "Ks": 0.3, "Y": 0.02, "sub": "Ethanol",
                      "Ki_H2": 0.5, "Ea_R": 8000.0},
}

# pH inhibition parameters {organism: (pH_LL, pH_UL) or None}
DEFAULT_PH_INHIB = {
    ORG["hydro"]:    None,  ORG["proteo"]:   None,  ORG["lipo"]:     None,
    ORG["acid_c"]:   (4.0, 5.5), ORG["acid_p"]:   (4.0, 5.5),
    ORG["gly_f"]:    (4.0, 5.5), ORG["syn_lcfa"]: (4.0, 5.5),
    ORG["syn_prop"]: (4.0, 5.5), ORG["syn_but"]:  (4.0, 5.5),
    ORG["syn_val"]:  (4.0, 5.5),
    ORG["meth_a"]:   (6.0, 7.0), ORG["meth_h2"]:  (5.0, 6.0),
    ORG["srb_ace"]:  (5.5, 6.5), ORG["srb_h2"]:   (5.0, 6.0),
    # Ethanol: fermenter tolerant (yeast-like), oxidiser like acetogens
    ORG["etoh_f"]:   (3.5, 5.0),   # very acid-tolerant
    ORG["syn_etoh"]: (4.0, 5.5),   # same as VFA oxidisers
}

# Decay and disintegration defaults
DEFAULT_K_DEC_PER_H = 0.02 / 24.0
DEFAULT_EA_R_DEC = 5000.0
DEFAULT_K_DIS_PER_H = 0.5 / 24.0
DEFAULT_EA_R_DIS = 5000.0
DEFAULT_F_XI_DECAY = 0.10
DEFAULT_DIS_FRACTIONS = dict(f_ch=0.20, f_pr=0.20, f_li=0.25, f_xI=0.25, f_sI=0.10)

# Inhibition constants
DEFAULT_KI_NH3 = 0.0018       # mol/L free NH₃
DEFAULT_KI_H2S_METH = 0.04    # g/L free H₂S (methanogens)
DEFAULT_KI_H2S_SRB = 0.20     # g/L free H₂S (SRB)

# Soluble inert fractions (Stage 10)
DEFAULT_F_SI_PROTEOLYSIS = 0.05
DEFAULT_F_SI_AA_ACID = 0.10

# Gas constant
_R_J = 8.31446


# ════════════════════════════════════════════════════════════════════════
#  Stoichiometry functions
# ════════════════════════════════════════════════════════════════════════

def _stoich_starch_hydrolysis(Y):
    nu_X = Y * _mw("Starch") / BIO_MW
    nu_W = (2 - 0.6*nu_X) / 2
    nu_G = (5 + nu_W - 0.9*nu_X) / 6
    nu_C = 6 - nu_X - 6*nu_G
    return [(("Starch",-1),("H2O",-nu_W),
             (ORG["hydro"],+nu_X),("Glucose",+nu_G),("CO2",+nu_C)),
            ("C","H","O")]

def _stoich_proteolysis(Y, f_sI=0.0):
    nu_X = Y * _mw("Protein") / BIO_MW
    nu_sI = f_sI * nu_X
    nu_X_active = nu_X - nu_sI
    nu_AA = 1.0 - 0.2*nu_X
    nu_C = 5.0 - nu_X - 5.0*nu_AA
    nu_W = 0.5*nu_X + 3.0*nu_AA + 2.0*nu_C - 2.0
    nu_H = (9.0 + 2.0*nu_W - 1.8*nu_X - 11.0*nu_AA) / 2.0
    coeffs = [("Protein",-1),("H2O",-nu_W),
             (ORG["proteo"],+nu_X_active),("AminoAcid",+nu_AA),
             ("CO2",+nu_C),("H2",+nu_H)]
    if nu_sI > 0:
        coeffs.append(("S_I",+nu_sI))
    return [tuple(coeffs), ("C","H","O","N")]

def _stoich_lipid_hydrolysis(Y):
    nu_X = Y * _mw("Lipid") / BIO_MW
    nu_G = 1.0; nu_L = (51 - nu_X - 3*nu_G) / 16
    nu_W = 0.5*nu_X + 3*nu_G + 2*nu_L - 6
    nu_H = (98 + 2*nu_W - 1.8*nu_X - 8*nu_G - 32*nu_L) / 2
    return [(("Lipid",-1),("H2O",-nu_W),
             (ORG["lipo"],+nu_X),("Glycerol",+nu_G),
             ("LCFA",+nu_L),("H2",+nu_H)),
            ("C","H","O")]

def _stoich_glucose_acidogenesis(Y, f_ac=0.41, f_pro=0.27, f_bu=0.13, f_h2=0.19):
    nu_X = Y * _mw("Glucose") / BIO_MW
    r_pro = (f_pro / 112.0) / (f_ac / 64.0)
    r_bu = (f_bu / 160.0) / (f_ac / 64.0)
    c_per_A = 2.0 + 3.0*r_pro + 4.0*r_bu
    o_per_A = 2.0 + 2.0*r_pro + 2.0*r_bu
    h_per_A = 4.0 + 6.0*r_pro + 8.0*r_bu
    nu_A = (-6 + 1.5*nu_X) / (o_per_A - 2*c_per_A)
    nu_P = r_pro * nu_A; nu_B = r_bu * nu_A
    nu_C = 6 - nu_X - c_per_A * nu_A
    nu_H = (12 - 1.8*nu_X - h_per_A*nu_A) / 2
    return [(("Glucose",-1),
             (ORG["acid_c"],+nu_X),("AceticAcid",+nu_A),
             ("Propionate",+nu_P),("Butyrate",+nu_B),
             ("CO2",+nu_C),("H2",+nu_H)),
            ("C","H","O")]

def _stoich_aa_acidogenesis(Y, r_P=0.30, r_B=0.15, r_V=0.10, f_sI=0.0):
    nu_X = Y * _mw("AminoAcid") / BIO_MW
    nu_sI = f_sI * nu_X; nu_X_active = nu_X - nu_sI
    nu_NH3 = 1.0-0.2*nu_X
    nu_A=1.0; nu_P=r_P; nu_B=r_B; nu_V=r_V
    c_used = nu_X+2*nu_A+3*nu_P+4*nu_B+5*nu_V
    nu_C = 5-c_used
    o_rhs = 0.5*nu_X+2*nu_A+2*nu_P+2*nu_B+2*nu_V+2*nu_C
    nu_W = o_rhs-3
    h_rhs = 1.8*nu_X+4*nu_A+6*nu_P+8*nu_B+10*nu_V+3*nu_NH3
    nu_H = (11+2*nu_W-h_rhs)/2
    coeffs = [("AminoAcid",-1),("H2O",-nu_W),
             (ORG["acid_p"],+nu_X_active),("AceticAcid",+nu_A),
             ("Propionate",+nu_P),("Butyrate",+nu_B),
             ("Valerate",+nu_V),("NH3",+nu_NH3),
             ("CO2",+nu_C),("H2",+nu_H)]
    if nu_sI > 0:
        coeffs.append(("S_I",+nu_sI))
    return [tuple(coeffs), ("C","H","O","N")]

def _stoich_glycerol_fermentation(Y):
    nu_X = Y * _mw("Glycerol") / BIO_MW
    nu_A = (3-nu_X)/2; nu_W = 0.5*nu_X
    nu_H = (8-1.8*nu_X-4*nu_A-2*nu_W)/2
    return [(("Glycerol",-1),
             (ORG["gly_f"],+nu_X),("AceticAcid",+nu_A),
             ("H2",+nu_H),("H2O",+nu_W)),
            ("C","H","O")]

def _stoich_lcfa_oxidation(Y):
    nu_X = Y * _mw("LCFA") / BIO_MW
    nu_A = (16-nu_X)/2; nu_W = 0.5*nu_X+2*nu_A-2
    nu_H = (32+2*nu_W-1.8*nu_X-4*nu_A)/2
    return [(("LCFA",-1),("H2O",-nu_W),
             (ORG["syn_lcfa"],+nu_X),("AceticAcid",+nu_A),("H2",+nu_H)),
            ("C","H","O")]

def _stoich_propionate_oxidation(Y):
    nu_X = Y * _mw("Propionate") / BIO_MW
    nu_A=1.0; nu_C=3.0-nu_X-2*nu_A
    a = 0.5*nu_X+2*nu_A+2*nu_C-2
    nu_H = (6+2*a-1.8*nu_X-4*nu_A)/2
    return [(("Propionate",-1),("H2O",-a),
             (ORG["syn_prop"],+nu_X),("AceticAcid",+nu_A),
             ("CO2",+nu_C),("H2",+nu_H)),
            ("C","H","O")]

def _stoich_butyrate_oxidation(Y):
    nu_X = Y * _mw("Butyrate") / BIO_MW
    nu_A = (4-nu_X)/2; a = 0.5*nu_X+2*nu_A-2
    nu_H = (8+2*a-1.8*nu_X-4*nu_A)/2
    return [(("Butyrate",-1),("H2O",-a),
             (ORG["syn_but"],+nu_X),("AceticAcid",+nu_A),("H2",+nu_H)),
            ("C","H","O")]

def _stoich_valerate_oxidation(Y):
    nu_X = Y * _mw("Valerate") / BIO_MW
    nu_P=1.0; nu_A=(5-nu_X-3*nu_P)/2
    a = 0.5*nu_X+2*nu_A+2*nu_P-2
    nu_H = (10+2*a-1.8*nu_X-4*nu_A-6*nu_P)/2
    return [(("Valerate",-1),("H2O",-a),
             (ORG["syn_val"],+nu_X),("AceticAcid",+nu_A),
             ("Propionate",+nu_P),("H2",+nu_H)),
            ("C","H","O")]

def _stoich_aceticlastic_meth(Y):
    nu_X = Y * _mw("AceticAcid") / BIO_MW
    C_r=2-nu_X; H_r=4-1.8*nu_X; O_r=2-0.5*nu_X
    nu_CO2 = (4*C_r+2*O_r-H_r)/8; nu_CH4=C_r-nu_CO2; nu_W=O_r-2*nu_CO2
    return [(("AceticAcid",-1),
             (ORG["meth_a"],+nu_X),("CH4",+nu_CH4),
             ("CO2",+nu_CO2),("H2O",+nu_W)),
            ("C","H","O")]

def _stoich_h2_meth(Y_per_CO2):
    nu_X = Y_per_CO2 * _mw("CO2") / BIO_MW
    nu_CH4=1-nu_X; nu_W=2-0.5*nu_X
    nu_H2 = (1.8*nu_X+4*nu_CH4+2*nu_W)/2
    return [(("H2",-nu_H2),("CO2",-1),
             (ORG["meth_h2"],+nu_X),("CH4",+nu_CH4),("H2O",+nu_W)),
            ("C","H","O")]

# Stage 11: SRB stoichiometry
def _stoich_ace_sulfate_reduction(Y):
    """Acetotrophic sulfate reduction (CHOS balanced).

    AceticAcid + SO₄ → Biomass + CO₂ + H₂S + H₂O
    Solving: nu_SO4 = (8 - 4.8×nu_X) / 10
    """
    nu_X = Y * _mw("AceticAcid") / BIO_MW
    nu_CO2 = 2.0 - nu_X
    nu_SO4 = (8.0 - 4.8 * nu_X) / 10.0
    nu_H2S = nu_SO4
    nu_H2O = 4.0 * nu_SO4 - 2.0 + 1.5 * nu_X
    return [(("AceticAcid", -1), ("SO4", -nu_SO4),
             (ORG["srb_ace"], +nu_X), ("CO2", +nu_CO2),
             ("H2S", +nu_H2S), ("H2O", +nu_H2O)),
            ("C", "H", "O", "S")]

def _stoich_h2_sulfate_reduction(Y_per_SO4):
    """Hydrogenotrophic sulfate reduction (CHOS balanced).

    H₂ + SO₄ + CO₂ → Biomass + H₂S + H₂O
    Per mol SO₄: nu_H2O = 4 + 1.5×nu_X, nu_H2 = 2.4×nu_X + 5
    """
    nu_X = Y_per_SO4 * _mw("SO4") / BIO_MW
    nu_CO2 = nu_X
    nu_H2S = 1.0
    nu_H2O = 4.0 + 1.5 * nu_X
    nu_H2 = 2.4 * nu_X + 5.0
    return [(("H2", -nu_H2), ("SO4", -1), ("CO2", -nu_CO2),
             (ORG["srb_h2"], +nu_X), ("H2S", +nu_H2S), ("H2O", +nu_H2O)),
            ("C", "H", "O", "S")]

# Ethanol extension stoichiometry
def _stoich_ethanol_fermentation(Y):
    """Glucose fermentation to ethanol + CO₂ (CHO balanced).

    C₆H₁₂O₆ → Biomass + Ethanol + CO₂ + H₂O

    Per mol Glucose:
        nu_X = Y × MW_glc / bio_MW
        C: 6 = nu_X + 2×nu_E + nu_CO2
        H: 12 = 1.8×nu_X + 6×nu_E + 2×nu_W
        O: 6 = 0.5×nu_X + nu_E + 2×nu_CO2 + nu_W

    Solving:
        nu_E = 2 - 0.4×nu_X
        nu_CO2 = 2 - 0.2×nu_X
        nu_W = 0.3×nu_X  (produced, sign is positive)
    """
    nu_X = Y * _mw("Glucose") / BIO_MW
    nu_E = 2.0 - 0.4 * nu_X
    nu_CO2 = 2.0 - 0.2 * nu_X
    nu_W = 0.3 * nu_X
    return [(("Glucose", -1),
             (ORG["etoh_f"], +nu_X), ("Ethanol", +nu_E),
             ("CO2", +nu_CO2), ("H2O", +nu_W)),
            ("C", "H", "O")]

def _stoich_ethanol_oxidation(Y):
    """Syntrophic ethanol oxidation to acetate + H₂ (CHO balanced).

    C₂H₅OH + H₂O → Biomass + CH₃COOH + H₂

    Per mol Ethanol:
        nu_X = Y × MW_etoh / bio_MW
        C: 2 = nu_X + 2×nu_A
        H: 6 + 2×nu_W = 1.8×nu_X + 4×nu_A + 2×nu_H2
        O: 1 + nu_W = 0.5×nu_X + 2×nu_A

    Solving:
        nu_A = (2 - nu_X) / 2 = 1 - 0.5×nu_X
        nu_W = 0.5×nu_X + 2×nu_A - 1 = 1 - 0.5×nu_X
        nu_H2 = (6 + 2×nu_W - 1.8×nu_X - 4×nu_A) / 2 = 2 - 0.4×nu_X
    """
    nu_X = Y * _mw("Ethanol") / BIO_MW
    nu_A = 1.0 - 0.5 * nu_X
    nu_W = 1.0 - 0.5 * nu_X
    nu_H2 = 2.0 - 0.4 * nu_X
    return [(("Ethanol", -1), ("H2O", -nu_W),
             (ORG["syn_etoh"], +nu_X), ("AceticAcid", +nu_A), ("H2", +nu_H2)),
            ("C", "H", "O")]

# Decay and disintegration stoichiometry
def _stoich_decay(org_name, f_XI=DEFAULT_F_XI_DECAY):
    return ((org_name, -1), ("X_composite", +(1.0-f_XI)), ("X_I", +f_XI))

def _stoich_disintegration(f_ch=0.20, f_pr=0.20, f_li=0.25, f_xI=0.25, f_sI=0.10):
    nu_S = f_ch * _mw("X_composite") / _mw("Starch")
    nu_P = f_pr * _mw("X_composite") / _mw("Protein")
    nu_L = f_li * _mw("X_composite") / _mw("Lipid")
    nu_xI = f_xI * _mw("X_composite") / _mw("X_I")
    nu_sI = f_sI * _mw("X_composite") / _mw("S_I")
    return (("X_composite",-1),("Starch",+nu_S),("Protein",+nu_P),
            ("Lipid",+nu_L),("X_I",+nu_xI),("S_I",+nu_sI))


# ════════════════════════════════════════════════════════════════════════
#  Rate function helpers
# ════════════════════════════════════════════════════════════════════════

def _ph_inhibition(pH, pH_UL, pH_LL):
    if pH >= pH_UL: return 1.0
    if pH <= 0.0: return 0.0
    return math.exp(-3.0 * ((pH - pH_UL) / (pH_UL - pH_LL)) ** 2)

def _get_pH_from_env(env):
    # H+ is a speciation-output species; read pH off env.pH rather
    # than env.concentrations (post chemistry-unification-1, B2).
    return float(env.pH) if env.has_pH else 7.0

def _arrhenius(T_K, Ea_R, T_ref_K):
    if Ea_R <= 0.0 or T_K <= 0.0: return 1.0
    return math.exp(Ea_R * (1.0 / T_ref_K - 1.0 / T_K))

def _pKa_NH4(T_K):
    Ka_ref = 10.0 ** (-9.25)
    Ka_T = Ka_ref * math.exp(-(51965.0 / _R_J) * (1.0 / T_K - 1.0 / 298.15))  # BSM2-corrected
    return -math.log10(max(Ka_T, 1e-30))

def _pKa_H2S(T_K):
    Ka_ref = 10.0 ** (-7.0)
    Ka_T = Ka_ref * math.exp(-(20000.0 / _R_J) * (1.0 / T_K - 1.0 / 298.15))
    return -math.log10(max(Ka_T, 1e-30))


# ════════════════════════════════════════════════════════════════════════
#  Rate function factories
# ════════════════════════════════════════════════════════════════════════

def _make_rate_functions(
    kin: Dict,
    ph_inhib: Dict,
    T_ref_K: float,
    Ki_NH3: float,
    Ki_H2S_meth: float,
    Ki_H2S_srb: float,
    h2s_inhibition: bool,
):
    """Build all rate function closures from kinetic parameters."""

    H2S_KI = {
        ORG["meth_a"]: Ki_H2S_meth, ORG["meth_h2"]: Ki_H2S_meth,
        ORG["srb_ace"]: Ki_H2S_srb, ORG["srb_h2"]: Ki_H2S_srb,
    }

    def _apply_ph(org, pH):
        params = ph_inhib.get(org)
        if params is None: return 1.0
        return _ph_inhibition(pH, params[1], params[0])

    def _free_nh3(env, pH):
        C_TAN = env.concentrations.get("NH3", 0.0) * _mw("NH3") / 14.007
        if C_TAN <= 0: return 0.0
        pKa = _pKa_NH4(env.T_K)
        return C_TAN / (1.0 + 10.0 ** (pKa - pH))

    def _free_h2s_gL(env, pH):
        C_total = env.concentrations.get("H2S", 0.0) * _mw("H2S")
        if C_total <= 0: return 0.0
        pKa = _pKa_H2S(env.T_K)
        return C_total / (1.0 + 10.0 ** (pH - pKa))

    def _h2s_inhib(org, env, pH):
        if not h2s_inhibition: return 1.0
        Ki = H2S_KI.get(org, 0.0)
        if Ki <= 0: return 1.0
        return Ki / (Ki + _free_h2s_gL(env, pH))

    # ── Individual rate factories ─────────────────────────────────────

    def rate_monod(org, sub, sub_mw, Y, mu, Ks, Ea_R):
        def r(env):
            S = env.concentrations.get(sub, 0) * sub_mw
            X = env.concentrations.get(org, 0) * BIO_MW
            if X <= 1e-30 or S <= 0: return 0.0
            pH = _get_pH_from_env(env)
            f_T = _arrhenius(env.T_K, Ea_R, T_ref_K)
            return (mu * f_T * S / (Ks + S) / Y) * X / sub_mw * env.V_L * _apply_ph(org, pH)
        return r

    def rate_inhibited(org, sub, sub_mw, Y, mu, Ks, inh, inh_mw, Ki, Ea_R):
        def r(env):
            S = env.concentrations.get(sub, 0) * sub_mw
            P = env.concentrations.get(inh, 0) * inh_mw
            X = env.concentrations.get(org, 0) * BIO_MW
            if X <= 1e-30 or S <= 0: return 0.0
            pH = _get_pH_from_env(env)
            f_T = _arrhenius(env.T_K, Ea_R, T_ref_K)
            return (mu * f_T * S / (Ks + S) * Ki / (Ki + P) / Y) * X / sub_mw * env.V_L * _apply_ph(org, pH)
        return r

    def rate_syntroph(org, sub, sub_mw, Y, mu, Ks, Ki_H2, Ea_R):
        def r(env):
            S = env.concentrations.get(sub, 0) * sub_mw
            H2g = env.concentrations.get("H2", 0) * _mw("H2")
            X = env.concentrations.get(org, 0) * BIO_MW
            if X <= 1e-30 or S <= 0: return 0.0
            pH = _get_pH_from_env(env)
            f_T = _arrhenius(env.T_K, Ea_R, T_ref_K)
            return (mu * f_T * S / (Ks + S) * Ki_H2 / (Ki_H2 + H2g) / Y) * X / sub_mw * env.V_L * _apply_ph(org, pH)
        return r

    def rate_ace_meth(org, sub, sub_mw, Y, mu, Ks, Ea_R):
        def r(env):
            S = env.concentrations.get(sub, 0) * sub_mw
            X = env.concentrations.get(org, 0) * BIO_MW
            if X <= 1e-30 or S <= 0: return 0.0
            pH = _get_pH_from_env(env)
            I_pH = _apply_ph(org, pH)
            C_NH3 = _free_nh3(env, pH)
            I_NH3 = Ki_NH3 / (Ki_NH3 + C_NH3)
            I_H2S = _h2s_inhib(org, env, pH)
            f_T = _arrhenius(env.T_K, Ea_R, T_ref_K)
            return (mu * f_T * S / (Ks + S) / Y) * X / sub_mw * env.V_L * I_pH * I_NH3 * I_H2S
        return r

    def rate_h2_meth(org, Y_co2, mu, Ks_h2, Ea_R):
        mw_h2 = _mw("H2"); mw_co2 = _mw("CO2")
        def r(env):
            H2g = env.concentrations.get("H2", 0) * mw_h2
            X = env.concentrations.get(org, 0) * BIO_MW
            CO2 = env.concentrations.get("CO2", 0)
            if X <= 1e-30 or H2g <= 0 or CO2 <= 0: return 0.0
            pH = _get_pH_from_env(env)
            I_pH = _apply_ph(org, pH)
            I_H2S = _h2s_inhib(org, env, pH)
            f_T = _arrhenius(env.T_K, Ea_R, T_ref_K)
            return (mu * f_T * H2g / (Ks_h2 + H2g) / Y_co2) * X / mw_co2 * env.V_L * I_pH * I_H2S
        return r

    def rate_ace_srb(org, Y, mu, Ks_ace, Ks_SO4, Ea_R):
        mw_ace = _mw("AceticAcid"); mw_so4 = _mw("SO4")
        def r(env):
            S_ace = env.concentrations.get("AceticAcid", 0) * mw_ace
            S_SO4 = env.concentrations.get("SO4", 0) * mw_so4
            X = env.concentrations.get(org, 0) * BIO_MW
            if X <= 1e-30 or S_ace <= 0 or S_SO4 <= 0: return 0.0
            pH = _get_pH_from_env(env)
            I_pH = _apply_ph(org, pH)
            I_H2S = _h2s_inhib(org, env, pH)
            f_T = _arrhenius(env.T_K, Ea_R, T_ref_K)
            monod = (S_ace / (Ks_ace + S_ace)) * (S_SO4 / (Ks_SO4 + S_SO4))
            return (mu * f_T * monod / Y) * X / mw_ace * env.V_L * I_pH * I_H2S
        return r

    def rate_h2_srb(org, Y_SO4, mu, Ks_H2, Ks_SO4, Ea_R):
        mw_h2 = _mw("H2"); mw_so4 = _mw("SO4")
        def r(env):
            S_H2 = env.concentrations.get("H2", 0) * mw_h2
            S_SO4 = env.concentrations.get("SO4", 0) * mw_so4
            X = env.concentrations.get(org, 0) * BIO_MW
            if X <= 1e-30 or S_H2 <= 0 or S_SO4 <= 0: return 0.0
            pH = _get_pH_from_env(env)
            I_pH = _apply_ph(org, pH)
            I_H2S = _h2s_inhib(org, env, pH)
            f_T = _arrhenius(env.T_K, Ea_R, T_ref_K)
            monod = (S_H2 / (Ks_H2 + S_H2)) * (S_SO4 / (Ks_SO4 + S_SO4))
            return (mu * f_T * monod / Y_SO4) * X / mw_so4 * env.V_L * I_pH * I_H2S
        return r

    def rate_decay(org_name, k_dec, Ea_R):
        def r(env):
            X = env.concentrations.get(org_name, 0.0)
            if X <= 0: return 0.0
            f_T = _arrhenius(env.T_K, Ea_R, T_ref_K)
            return k_dec * f_T * X * env.V_L
        return r

    def rate_first_order(sp_name, k_rate, sp_mw, Ea_R):
        def r(env):
            C = env.concentrations.get(sp_name, 0.0)
            if C <= 0: return 0.0
            f_T = _arrhenius(env.T_K, Ea_R, T_ref_K)
            return k_rate * f_T * C * env.V_L
        return r

    return {
        "monod": rate_monod, "inhibited": rate_inhibited,
        "syntroph": rate_syntroph, "ace_meth": rate_ace_meth,
        "h2_meth": rate_h2_meth, "ace_srb": rate_ace_srb,
        "h2_srb": rate_h2_srb, "decay": rate_decay,
        "first_order": rate_first_order,
    }


# ════════════════════════════════════════════════════════════════════════
#  Stoichiometry entry builder
# ════════════════════════════════════════════════════════════════════════

def _make_entries(coeff_tuples):
    """Convert (species, coeff) tuples to StoichiometryEntry list."""
    entries = []
    for sp, coeff in coeff_tuples:
        sp_obj = _get_species(sp)
        entries.append(
            StoichiometryEntry(species=sp_obj, phase="liquid", coefficient=coeff)
        )
    return entries


_ADM1_SPECIES_CACHE: Dict[str, Species] = {}


def _get_species(sp_id: str) -> Species:
    """Return (or build and cache) a ``Species`` object for an ADM1 species id."""
    cached = _ADM1_SPECIES_CACHE.get(sp_id)
    if cached is not None:
        return cached
    # Universal inorganics: alias the common_species objects so
    # cross-model composition does not trip the soft-conflict check.
    if sp_id == "CO2":
        sp_obj = _CO2_sp
    elif sp_id == "NH3":
        sp_obj = _NH3_sp
    elif sp_id == "H2O":
        sp_obj = _H2O_sp
    elif sp_id in SPECIES:
        atoms, mw = SPECIES[sp_id]
        sp_obj = Species(id=sp_id, atoms=dict(atoms), charge=0, MW=float(mw))
    elif sp_id in CHON_ORGS:
        sp_obj = Species(id=sp_id, atoms=dict(BIO_CHON), charge=0, MW=float(BIO_MW))
    else:
        sp_obj = Species(id=sp_id, atoms=dict(BIO_CHO), charge=0, MW=float(BIO_MW))
    _ADM1_SPECIES_CACHE[sp_id] = sp_obj
    return sp_obj


# ════════════════════════════════════════════════════════════════════════
#  Public API: build_adm1_reactions()
# ════════════════════════════════════════════════════════════════════════

def build_adm1_reactions(
    *,
    sulfate_reduction: bool = True,
    ethanol: bool = False,
    soluble_inerts: bool = True,
    h2s_inhibition: bool = None,
    kinetics: Dict = None,
    ph_inhib: Dict = None,
    T_ref_K: float = DEFAULT_T_REF_K,
    Ki_NH3: float = DEFAULT_KI_NH3,
    Ki_H2S_meth: float = DEFAULT_KI_H2S_METH,
    Ki_H2S_srb: float = DEFAULT_KI_H2S_SRB,
    k_dec_per_h: float = DEFAULT_K_DEC_PER_H,
    Ea_R_dec: float = DEFAULT_EA_R_DEC,
    k_dis_per_h: float = DEFAULT_K_DIS_PER_H,
    Ea_R_dis: float = DEFAULT_EA_R_DIS,
    f_sI_proteolysis: float = DEFAULT_F_SI_PROTEOLYSIS,
    f_sI_aa_acid: float = DEFAULT_F_SI_AA_ACID,
    dis_fractions: Dict = None,
    verbose: bool = True,
) -> ReactionSystem:
    """Build an ADM1 reaction system with optional extensions.

    Parameters
    ----------
    sulfate_reduction : bool
        Enable ADM1-S sulfate reduction (2 SRB reactions + 2 decay).
    ethanol : bool
        Enable ethanol extension (2 reactions + 2 decay):
        glucose → ethanol fermentation and syntrophic ethanol → acetate
        oxidation.  Adds EtOH_Ferm and Syn_EtOH organisms.
    soluble_inerts : bool
        Enable S_I production from proteolysis and amino acid fermentation.
    h2s_inhibition : bool or None
        Enable H₂S inhibition on methanogens and SRB.
        If None, auto-enabled when ``sulfate_reduction=True``.
    kinetics : dict, optional
        Override default kinetic parameters. Merged with defaults.
    ph_inhib : dict, optional
        Override pH inhibition parameters.
    verbose : bool
        Print reaction summary.

    Returns
    -------
    ReactionSystem
    """
    if h2s_inhibition is None:
        h2s_inhibition = sulfate_reduction

    kin = dict(DEFAULT_KINETICS)
    if kinetics:
        kin.update(kinetics)

    phi = dict(DEFAULT_PH_INHIB)
    if ph_inhib:
        phi.update(ph_inhib)

    dis_frac = dict(DEFAULT_DIS_FRACTIONS)
    if dis_fractions:
        dis_frac.update(dis_fractions)

    f_sI_pro = f_sI_proteolysis if soluble_inerts else 0.0
    f_sI_aa = f_sI_aa_acid if soluble_inerts else 0.0

    rf = _make_rate_functions(kin, phi, T_ref_K, Ki_NH3,
                              Ki_H2S_meth, Ki_H2S_srb, h2s_inhibition)
    k = kin

    # ── Core biochemical reactions ────────────────────────────────────
    rxn_defs = [
        ("0a_starch_hyd", _stoich_starch_hydrolysis, k[ORG["hydro"]]["Y"],
         rf["monod"](ORG["hydro"],"Starch",_mw("Starch"),
                     k[ORG["hydro"]]["Y"],k[ORG["hydro"]]["mu"],
                     k[ORG["hydro"]]["Ks"],k[ORG["hydro"]]["Ea_R"])),

        ("0b_proteolysis", lambda Y: _stoich_proteolysis(Y, f_sI=f_sI_pro),
         k[ORG["proteo"]]["Y"],
         rf["monod"](ORG["proteo"],"Protein",_mw("Protein"),
                     k[ORG["proteo"]]["Y"],k[ORG["proteo"]]["mu"],
                     k[ORG["proteo"]]["Ks"],k[ORG["proteo"]]["Ea_R"])),

        ("0c_lipid_hyd", _stoich_lipid_hydrolysis, k[ORG["lipo"]]["Y"],
         rf["monod"](ORG["lipo"],"Lipid",_mw("Lipid"),
                     k[ORG["lipo"]]["Y"],k[ORG["lipo"]]["mu"],
                     k[ORG["lipo"]]["Ks"],k[ORG["lipo"]]["Ea_R"])),

        ("1a_glc_acid", _stoich_glucose_acidogenesis, k[ORG["acid_c"]]["Y"],
         rf["inhibited"](ORG["acid_c"],"Glucose",_mw("Glucose"),
                         k[ORG["acid_c"]]["Y"],k[ORG["acid_c"]]["mu"],
                         k[ORG["acid_c"]]["Ks"],"AceticAcid",_mw("AceticAcid"),
                         k[ORG["acid_c"]]["Ki_ace"],k[ORG["acid_c"]]["Ea_R"])),

        ("1b_aa_acid", lambda Y: _stoich_aa_acidogenesis(Y, f_sI=f_sI_aa),
         k[ORG["acid_p"]]["Y"],
         rf["inhibited"](ORG["acid_p"],"AminoAcid",_mw("AminoAcid"),
                         k[ORG["acid_p"]]["Y"],k[ORG["acid_p"]]["mu"],
                         k[ORG["acid_p"]]["Ks"],"AceticAcid",_mw("AceticAcid"),
                         k[ORG["acid_p"]]["Ki_ace"],k[ORG["acid_p"]]["Ea_R"])),

        ("1c_gly_ferm", _stoich_glycerol_fermentation, k[ORG["gly_f"]]["Y"],
         rf["inhibited"](ORG["gly_f"],"Glycerol",_mw("Glycerol"),
                         k[ORG["gly_f"]]["Y"],k[ORG["gly_f"]]["mu"],
                         k[ORG["gly_f"]]["Ks"],"AceticAcid",_mw("AceticAcid"),
                         k[ORG["gly_f"]]["Ki_ace"],k[ORG["gly_f"]]["Ea_R"])),

        ("1d_lcfa_ox", _stoich_lcfa_oxidation, k[ORG["syn_lcfa"]]["Y"],
         rf["syntroph"](ORG["syn_lcfa"],"LCFA",_mw("LCFA"),
                        k[ORG["syn_lcfa"]]["Y"],k[ORG["syn_lcfa"]]["mu"],
                        k[ORG["syn_lcfa"]]["Ks"],k[ORG["syn_lcfa"]]["Ki_H2"],
                        k[ORG["syn_lcfa"]]["Ea_R"])),

        ("2a_prop_ox", _stoich_propionate_oxidation, k[ORG["syn_prop"]]["Y"],
         rf["syntroph"](ORG["syn_prop"],"Propionate",_mw("Propionate"),
                        k[ORG["syn_prop"]]["Y"],k[ORG["syn_prop"]]["mu"],
                        k[ORG["syn_prop"]]["Ks"],k[ORG["syn_prop"]]["Ki_H2"],
                        k[ORG["syn_prop"]]["Ea_R"])),

        ("2b_but_ox", _stoich_butyrate_oxidation, k[ORG["syn_but"]]["Y"],
         rf["syntroph"](ORG["syn_but"],"Butyrate",_mw("Butyrate"),
                        k[ORG["syn_but"]]["Y"],k[ORG["syn_but"]]["mu"],
                        k[ORG["syn_but"]]["Ks"],k[ORG["syn_but"]]["Ki_H2"],
                        k[ORG["syn_but"]]["Ea_R"])),

        ("2c_val_ox", _stoich_valerate_oxidation, k[ORG["syn_val"]]["Y"],
         rf["syntroph"](ORG["syn_val"],"Valerate",_mw("Valerate"),
                        k[ORG["syn_val"]]["Y"],k[ORG["syn_val"]]["mu"],
                        k[ORG["syn_val"]]["Ks"],k[ORG["syn_val"]]["Ki_H2"],
                        k[ORG["syn_val"]]["Ea_R"])),

        ("3a_ace_meth", _stoich_aceticlastic_meth, k[ORG["meth_a"]]["Y"],
         rf["ace_meth"](ORG["meth_a"],"AceticAcid",_mw("AceticAcid"),
                        k[ORG["meth_a"]]["Y"],k[ORG["meth_a"]]["mu"],
                        k[ORG["meth_a"]]["Ks"],k[ORG["meth_a"]]["Ea_R"])),

        ("3b_h2_meth", _stoich_h2_meth, k[ORG["meth_h2"]]["Y"],
         rf["h2_meth"](ORG["meth_h2"],k[ORG["meth_h2"]]["Y"],
                       k[ORG["meth_h2"]]["mu"],k[ORG["meth_h2"]]["Ks"],
                       k[ORG["meth_h2"]]["Ea_R"])),
    ]

    # Stage 11: SRB reactions
    if sulfate_reduction:
        rxn_defs.extend([
            ("6a_ace_srb", _stoich_ace_sulfate_reduction, k[ORG["srb_ace"]]["Y"],
             rf["ace_srb"](ORG["srb_ace"],k[ORG["srb_ace"]]["Y"],
                           k[ORG["srb_ace"]]["mu"],k[ORG["srb_ace"]]["Ks_ace"],
                           k[ORG["srb_ace"]]["Ks_SO4"],k[ORG["srb_ace"]]["Ea_R"])),
            ("6b_h2_srb", _stoich_h2_sulfate_reduction, k[ORG["srb_h2"]]["Y"],
             rf["h2_srb"](ORG["srb_h2"],k[ORG["srb_h2"]]["Y"],
                          k[ORG["srb_h2"]]["mu"],k[ORG["srb_h2"]]["Ks_H2"],
                          k[ORG["srb_h2"]]["Ks_SO4"],k[ORG["srb_h2"]]["Ea_R"])),
        ])

    # Ethanol extension: glucose → ethanol fermentation + ethanol → acetate oxidation
    if ethanol:
        rxn_defs.extend([
            ("4a_etoh_ferm", _stoich_ethanol_fermentation, k[ORG["etoh_f"]]["Y"],
             rf["monod"](ORG["etoh_f"], "Glucose", _mw("Glucose"),
                         k[ORG["etoh_f"]]["Y"], k[ORG["etoh_f"]]["mu"],
                         k[ORG["etoh_f"]]["Ks"], k[ORG["etoh_f"]]["Ea_R"])),
            ("4b_etoh_ox", _stoich_ethanol_oxidation, k[ORG["syn_etoh"]]["Y"],
             rf["syntroph"](ORG["syn_etoh"], "Ethanol", _mw("Ethanol"),
                            k[ORG["syn_etoh"]]["Y"], k[ORG["syn_etoh"]]["mu"],
                            k[ORG["syn_etoh"]]["Ks"], k[ORG["syn_etoh"]]["Ki_H2"],
                            k[ORG["syn_etoh"]]["Ea_R"])),
        ])

    # Build biochemical reactions
    reactions = []
    _log = logger.info if verbose else logger.debug
    _log("ADM1 biochemical reactions:")
    for label, stoich_fn, Y, rate_fn in rxn_defs:
        coeffs, bal = stoich_fn(Y)
        entries = _make_entries(coeffs)
        rxn = KineticReaction(stoichiometry=entries, rate_fn=rate_fn,
                       balance_elements=bal, label=label)
        reactions.append(rxn)
        consumed = " + ".join(f"{abs(c):.3f} {s}" for s,c in coeffs if c<0)
        produced = " + ".join(f"{c:.3f} {s}" for s,c in coeffs if c>0)
        _log("  %s: %s → %s", label, consumed, produced)

    # ── Decay reactions ───────────────────────────────────────────────
    all_orgs = list(CORE_ORGS)
    if sulfate_reduction:
        all_orgs.extend(SRB_ORGS)
    if ethanol:
        all_orgs.extend(ETOH_ORGS)

    _log("Decay reactions (%d):", len(all_orgs))
    for org_name in all_orgs:
        label = f"4_decay_{org_name}"
        coeffs = _stoich_decay(org_name)
        entries = _make_entries(coeffs)
        rfn = rf["decay"](org_name, k_dec_per_h, Ea_R_dec)
        rxn = KineticReaction(stoichiometry=entries, rate_fn=rfn,
                       balance_elements=(), label=label)
        reactions.append(rxn)
        consumed = " + ".join(f"{abs(c):.3f} {s}" for s,c in coeffs if c<0)
        produced = " + ".join(f"{c:.3f} {s}" for s,c in coeffs if c>0)
        _log("  %s: %s → %s", label, consumed, produced)

    # ── Disintegration ────────────────────────────────────────────────
    dis_coeffs = _stoich_disintegration(**dis_frac)
    dis_entries = _make_entries(dis_coeffs)
    dis_rfn = rf["first_order"]("X_composite", k_dis_per_h,
                                _mw("X_composite"), Ea_R_dis)
    dis_rxn = KineticReaction(stoichiometry=dis_entries, rate_fn=dis_rfn,
                       balance_elements=(), label="5_disintegration")
    reactions.append(dis_rxn)
    consumed = " + ".join(f"{abs(c):.3f} {s}" for s,c in dis_coeffs if c<0)
    produced = " + ".join(f"{c:.3f} {s}" for s,c in dis_coeffs if c>0)
    _log("Disintegration: %s → %s", consumed, produced)
    n_bio = len(rxn_defs)
    _log("Total: %d reactions (%d biochemical + %d decay + 1 disintegration)",
         len(reactions), n_bio, len(all_orgs))

    # Cross-phase partition declarations (chemistry-unification-3b C7).
    # These auto-wire speciation_keys and build speciation_ladders at CV
    # construction via derive_speciation_keys. No log_K — Henry's law
    # constant lives on the gas-liquid link.
    reactions.append(EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=_CO2_sp, phase="gas",    coefficient=-1.0),
            StoichiometryEntry(species=_CO2_sp, phase="liquid", coefficient=+1.0),
        ],
        balance_elements=("C", "O"),
        label="partition_CO2",
    ))
    if sulfate_reduction:
        h2s_sp = _get_species("H2S")
        reactions.append(EquilibriumReaction(
            stoichiometry=[
                StoichiometryEntry(species=h2s_sp, phase="gas",    coefficient=-1.0),
                StoichiometryEntry(species=h2s_sp, phase="liquid", coefficient=+1.0),
            ],
            balance_elements=("S", "H"),
            label="partition_H2S",
        ))

    return ReactionSystem(reactions, label="ADM1")


# ════════════════════════════════════════════════════════════════════════
#  Public API: build_adm1_cv()
# ════════════════════════════════════════════════════════════════════════

def build_adm1_cv(
    reaction_system: ReactionSystem,
    *,
    V_total_L: float = 2.0,
    headspace_frac: float = 0.20,
    T_K: float = 308.15,
    sulfate_reduction: bool = True,
    ethanol: bool = False,
    activity_model: str = "davies",
):
    """Build a fermenter ControlVolume configured for ADM1 simulation.

    Parameters
    ----------
    reaction_system : ReactionSystem
        From :func:`build_adm1_reactions`.
    V_total_L : float
        Total vessel volume (L).
    headspace_frac : float
        Fraction of vessel that is headspace.
    T_K : float
        Temperature (K).
    sulfate_reduction : bool
        Include H₂S gas-liquid transfer and speciation correction.
    ethanol : bool
        Include ethanol gas-liquid transfer (no speciation correction
        needed — ethanol does not dissociate).
    activity_model : str
        Activity model: ``"davies"``, ``"sit"``, or ``"ideal"``.

    Returns
    -------
    ControlVolume
    """
    from vlmodels.fermenter.config.builder import FermenterBuilder
    from PyOMES.core.boundaries import PressureReliefVent
    from PyOMES.chemical_equilibrium.engine import BisectionChemicalEquilibriumEngine

    use_activity = activity_model.lower() != "ideal"

    # Equilibria are pre-bucketed by the ReactionSystem; the engine
    # consumes the single-phase + cross-phase lists directly (the
    # engine silently skips cross-phase reactions internally).
    equilibrium_rxns = (
        list(reaction_system.single_phase_equilibria)
        + list(reaction_system.cross_phase_equilibria)
    )

    # chemistry-unification-2 migration: .speciation_correction_pka()
    # (which carried pKas on the transfer link) was removed.  The link
    # now reads alphas from PropertyResult.alphas via speciation_keys;
    # the pKa values themselves live on declared equilibrium reactions
    # in build_adm1_reactions() and are consumed by the speciation
    # engine via BisectionChemicalEquilibriumEngine.from_reactions(). For VFAs/H2S to
    # see a non-trivial alpha here, build_adm1_reactions() must declare
    # equilibrium reactions for them — Phase 3 work (ChemistryDatabase
    # rollout). Until then, the link will fall back to alpha=1.0 for
    # any species lacking a declared equilibrium reaction.
    # chemistry-unification-3b C7: speciation_correction() calls removed.
    # Cross-phase partition reactions declared in build_adm1_reactions()
    # drive derive_speciation_keys() at CV construction; no manual
    # speciation_keys wiring needed.
    b = (FermenterBuilder()
         .vessel(V_total_L=V_total_L, headspace_frac=headspace_frac,
                 T_K=T_K, yO2_init=0.0, yCO2_init=0.0)
         .no_gas_feed()
         .transfer_equilibrium()
         .transfer_species("CH4")
         .transfer_species("H2"))

    if sulfate_reduction:
        b = b.transfer_species("H2S")

    if ethanol:
        b = b.transfer_species("Ethanol")  # no speciation correction — no dissociation

    cv = (b.chemistry(use_activity=use_activity,
                      activity_model=activity_model)
            .reaction_system(reaction_system)
            .label("ADM1")
            .build())

    # Install equilibrium reactions on the speciation engine. If the
    # builder did not declare any equilibrium reactions, leave the
    # default lazy-build path in place. state-unification C4d:
    # cv.reaction_system.attach_engine is the canonical attachment.
    if len(equilibrium_rxns) > 0:
        engine = BisectionChemicalEquilibriumEngine.from_reactions(
            equilibrium_rxns,
            activity_model=activity_model,
            use_activity=use_activity,
            T_K=T_K,
        )
        cv.reaction_system.attach_engine(engine)

    cv.boundaries.append(PressureReliefVent(P_set_atm=1.05, mode="instant"))

    # H2O gas-liquid transfer (CP4 of LAYER1_GAP_CLOSURE) — replaces
    # WaterVapourBoundary, retired per §6.2 of MASS_EXCHANGE_ARCHITECTURE.md.
    # WaterVapourBoundary "conjured" gas-phase water from nowhere without
    # decrementing the liquid phase (never mass-conserving; VentWaterLoss
    # was a partial patch this model didn't even use). RaoultEquilibrium's
    # PartitionModel role moves mass between the *actual* tracked liquid
    # and gas H2O pools instead. That requires a real liquid H2O pool to
    # draw from — this model didn't track one before this phase — seeded
    # here at the standard pure-water concentration (matches
    # nr_engine.py's _C_WATER_MOL_L).
    liq = cv.phases["liquid"]
    liq.n_mol["H2O"] = 55.51 * float(liq.V_L)
    gas = cv.phases["gas"]
    gas.n_mol.setdefault("H2O", 0.0)

    import warnings as _warnings
    from PyOMES.chemistry import RaoultEquilibrium
    from PyOMES.core.gas_liquid_link import KineticGasLiquidLink

    with _warnings.catch_warnings():
        # KineticGasLiquidLink is deprecated in favour of the
        # transfer_models kwarg — but transfer_models is only consumable
        # at ControlVolume construction time, and cv is already built by
        # FermenterBuilder above (which doesn't support a generic
        # PartitionModel, only Henry constants, in its transfer_species()
        # API). Direct construction here is the deliberate, internal
        # escape hatch, matching how ControlVolume's own
        # _build_transfer_link does the same suppression for its
        # transfer_models-driven link.
        _warnings.simplefilter("ignore", DeprecationWarning)
        water_link = KineticGasLiquidLink(
            gas_cv_key="ADM1", gas_phase_key="gas",
            liquid_cv_key="ADM1", liquid_phase_key="liquid",
            partition_models={"H2O": RaoultEquilibrium()},
            equilibrium_species={"H2O"},
            # Water does dissociate (H2O <-> H+ + OH-), but Kw=1e-14 means
            # the dissociated fraction (~1.8e-9 of the ~55.5 M pool) is
            # negligible — unlike CO2, where a significant fraction
            # converts to HCO3-/CO3-- at typical pH. No alpha correction
            # (speciation_keys entry) is needed for the Raoult transfer.
            speciation_keys={},
            _label="ADM1_water_vapour",
        )
    # Appended to both lists: internal_interfaces is what advance()/
    # step_internal_transfer() actually iterate; _explicit_interfaces is
    # what cv.snapshot() uses to reconstruct a copy (see ControlVolume.
    # snapshot()) — omitting it would silently drop water transfer after
    # any snapshot/restore cycle (e.g. inside an adaptive ODE solver).
    cv._explicit_interfaces.append(water_link)
    cv.internal_interfaces.append(water_link)

    return cv


# ════════════════════════════════════════════════════════════════════════
#  Public API: seed_adm1_strong_ions()
# ════════════════════════════════════════════════════════════════════════

def seed_adm1_strong_ions(
    cv,
    *,
    CT_Na: float = 0.100,
    CT_Cl: float = 0.030,
) -> None:
    """Seed ADM1 strong-ion totals into the CV's liquid n_mol.

    state-unification C4: strong ions live as Species in
    ``phase.n_mol`` (PHREEQC convention; see Option C of the C4
    sub-decisions). Replaces the deleted
    ``make_adm1_chem_env_fn`` which returned a chem_env dict
    threading ``strong_kwargs`` through ``cv.advance``.

    Parameters
    ----------
    cv : ControlVolume
        The ADM1 CV to seed.
    CT_Na, CT_Cl : float
        Strong-ion totals (mol/L). Stored under
        ``cv.phases["liquid"].n_mol["Na+"]`` / ``["Cl-"]``.
    """
    liq = cv.phases.get("liquid")
    if liq is None:
        return
    V_L = float(getattr(liq, "V_L", 1.0))
    if CT_Na:
        liq.n_mol["Na+"] = float(CT_Na) * V_L
    if CT_Cl:
        liq.n_mol["Cl-"] = float(CT_Cl) * V_L


# ════════════════════════════════════════════════════════════════════════
#  Public API: get_organism_list()
# ════════════════════════════════════════════════════════════════════════

def get_organism_list(sulfate_reduction: bool = True, ethanol: bool = False) -> List[str]:
    """Return the list of organism species names for the given configuration."""
    orgs = list(CORE_ORGS)
    if sulfate_reduction:
        orgs.extend(SRB_ORGS)
    if ethanol:
        orgs.extend(ETOH_ORGS)
    return orgs
