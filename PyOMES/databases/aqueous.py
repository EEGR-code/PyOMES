# -*- coding: utf-8 -*-
"""AQUEOUS_DEFAULT: core inorganic aqueous chemistry database.

Covers the universal acid-base systems present in essentially every
aqueous model: water dissociation, carbonate (monoprotic, BSM2-style),
and ammonium/ammonia.

Usage::

    from PyOMES.databases.aqueous import AQUEOUS_DEFAULT
    db = AQUEOUS_DEFAULT.extend(species={"MyAcid": ...}, reactions=[...])
"""
from __future__ import annotations

from .database import ChemistryDatabase
from PyOMES.chemistry.species import Species
from PyOMES.reactions.equilibrium.reaction import EquilibriumReaction
from PyOMES.reactions.reaction_system import ReactionSystem
from PyOMES.reactions.stoichiometry import StoichiometryEntry
from PyOMES.thermo.framework import ThermoFramework

# Water, the carbonate system and the ammonium system. MW is computed from
# the atoms. Other databases and models import these objects from here, so
# one definition is shared by identity.
H_plus   = Species(id="H+",    atoms={"H": 1},                  charge=+1)
OH_minus = Species(id="OH-",   atoms={"O": 1, "H": 1},          charge=-1)
H2O      = Species(id="H2O",   atoms={"H": 2, "O": 1},          charge=0)

CO2        = Species(id="CO2",    atoms={"C": 1, "O": 2},              charge=0)
HCO3_minus = Species(id="HCO3-", atoms={"H": 1, "C": 1, "O": 3},      charge=-1)
CO3_2minus = Species(id="CO3--", atoms={"C": 1, "O": 3},               charge=-2)

NH3      = Species(id="NH3",   atoms={"N": 1, "H": 3},          charge=0)
NH4_plus = Species(id="NH4+",  atoms={"N": 1, "H": 4},          charge=+1)

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