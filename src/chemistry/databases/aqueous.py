# -*- coding: utf-8 -*-
"""AQUEOUS_DEFAULT: core inorganic aqueous chemistry database.

Covers the universal acid-base systems present in essentially every
aqueous model: water dissociation, carbonate (monoprotic, BSM2-style),
and ammonium/ammonia.

Usage::

    from PyOMES.chemistry.databases.aqueous import AQUEOUS_DEFAULT
    db = AQUEOUS_DEFAULT.extend(species={"MyAcid": ...}, reactions=[...])
"""
from __future__ import annotations

from ..database import ChemistryDatabase
from ..common_species import (
    H_plus, OH_minus, H2O,
    CO2, HCO3_minus, CO3_2minus,
    NH3, NH4_plus,
)
from ..species import Species
from ...reactions.equilibrium import EquilibriumReaction
from ...reactions.reaction_system import ReactionSystem
from ...reactions.stoichiometry import StoichiometryEntry
from ...thermo.framework import ThermoFramework

# BSM2-canonical pKa and dH values (Rosen & Jeppsson 2006)
_PKW = 14.0
_DH_W = 55900.0        # J/mol
_PKA_CO2_1 = 6.35
_DH_CO2_1 = 7646.0     # J/mol
_PKA_NH4 = 9.25
_DH_NH4 = 51965.0      # J/mol

_T_REF_K = 298.15

_THERMO = ThermoFramework()  # ideal by default; override via ChemistryDatabase.extend(thermo=...)

_SPECIES: dict = {
    "H+":    H_plus,
    "OH-":   OH_minus,
    "H2O":   H2O,
    "CO2":   CO2,
    "HCO3-": HCO3_minus,
    "CO3--": CO3_2minus,
    "NH3":   NH3,
    "NH4+":  NH4_plus,
}

_REACTIONS = ReactionSystem([
    # H2O ⇌ H+ + OH-
    EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=H2O,       phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=H_plus,     phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=OH_minus,   phase="liquid", coefficient=+1.0),
        ],
        log_K=-_PKW, dH_J_per_mol=_DH_W, T_ref_K=_T_REF_K,
        balance_elements=("H", "O"),
        label="eq_water",
    ),
    # CO2 + H2O ⇌ HCO3- + H+  (1st dissociation, monoprotic treatment)
    EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=CO2,        phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=H2O,        phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=HCO3_minus, phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=H_plus,     phase="liquid", coefficient=+1.0),
        ],
        log_K=-_PKA_CO2_1, dH_J_per_mol=_DH_CO2_1, T_ref_K=_T_REF_K,
        total_id="CO2",
        balance_elements=("C", "H", "O"),
        label="eq_CO2",
    ),
    # NH4+ ⇌ NH3 + H+
    EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=NH4_plus, phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=NH3,      phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=H_plus,   phase="liquid", coefficient=+1.0),
        ],
        log_K=-_PKA_NH4, dH_J_per_mol=_DH_NH4, T_ref_K=_T_REF_K,
        total_id="NH3",
        balance_elements=("N", "H"),
        label="eq_NH4",
    ),
])

AQUEOUS_DEFAULT = ChemistryDatabase(
    thermo=_THERMO,
    species=_SPECIES,
    reactions=_REACTIONS,
)