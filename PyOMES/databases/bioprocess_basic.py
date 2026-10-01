# -*- coding: utf-8 -*-
"""BIOPROCESS_BASIC: aqueous chemistry extended for aerobic bioprocesses.

Extends :data:`AQUEOUS_DEFAULT` with:

- Phosphate and bisulfate equilibrium ladders.
- Species for aerobic growth models: the organic acids used as substrates
  (acetic, propionic, butyric and citric acid), yeast biomass (:data:`Yeast`,
  CHNO, and :data:`Yeast_CHO`), and the gases O2 and N2.
- :data:`AIR`, a gas composition (``{species id: mole fraction}``) for a
  headspace or feed that starts with, or is fed, air.

The lump salt species :data:`NH4Cl`, :data:`KH2PO4`, and :data:`NaOH` are
defined as :class:`~PyOMES.chemistry.Species` objects so that their molar
masses are available for medium-recipe calculations, but they carry no
dissolution reactions in this database.  Initial medium compositions should
be specified in ionic form directly (``NH4+`` + ``Cl-``,
``K+`` + ``H2PO4-``, etc.).  NaOH pH correction is handled by
:meth:`~PyOMES.core.ControlVolume.equilibrate_to_pH` via the built-in
strong-corrector map, which adds ``Na+`` to the charge balance.

Usage::

    from PyOMES.databases.bioprocess_basic import BIOPROCESS_BASIC
    from PyOMES.databases.bioprocess_basic import NH4Cl, KH2PO4, NaOH
"""
from __future__ import annotations

from PyOMES.chemistry.common_species import (
    H_plus,
    H3PO4, H2PO4_minus, HPO4_2minus, PO4_3minus,
    HSO4_minus, SO4_2minus,
    K_plus, Cl_minus, Na_plus,
)
from PyOMES.chemistry.species import Species
from PyOMES.reactions.equilibrium.reaction import EquilibriumReaction
from PyOMES.reactions.equilibrium.interphase import HenryEquilibrium
from PyOMES.reactions.reaction_system import ReactionSystem
from PyOMES.reactions.stoichiometry import StoichiometryEntry
from .aqueous import AQUEOUS_DEFAULT

_T_REF_K = 298.15

# ---------------------------------------------------------------------------
# Undissociated salt / corrector species
# ---------------------------------------------------------------------------
NH4Cl  = Species(id="NH4Cl",  atoms={"N": 1, "H": 4, "Cl": 1},       charge=0)
KH2PO4 = Species(id="KH2PO4", atoms={"K": 1, "H": 2, "P": 1, "O": 4}, charge=0)
NaOH   = Species(id="NaOH",   atoms={"Na": 1, "O": 1, "H": 1},        charge=0)

# ---------------------------------------------------------------------------
# Substrates, biomass and gases for aerobic growth
# ---------------------------------------------------------------------------
AceticAcid    = Species(id="AceticAcid",    atoms={"C": 2, "H": 4, "O": 2})
PropionicAcid = Species(id="PropionicAcid", atoms={"C": 3, "H": 6, "O": 2})
ButyricAcid   = Species(id="ButyricAcid",   atoms={"C": 4, "H": 8, "O": 2})
CitricAcid    = Species(id="CitricAcid",    atoms={"C": 6, "H": 8, "O": 7})

# Yeast biomass, one C-mol; MW is the formula weight of the atoms.
Yeast     = Species(id="Yeast",     atoms={"C": 1, "H": 1.61, "O": 0.56, "N": 0.16})
Yeast_CHO = Species(id="Yeast_CHO", atoms={"C": 1, "H": 1.61, "O": 0.56})

O2 = Species(id="O2", atoms={"O": 2})
N2 = Species(id="N2", atoms={"N": 2})

# Dry air as a gas composition ({species id: mole fraction}): O2 and CO2 at
# their atmospheric fractions, N2 the balance. For a vessel or feed that a
# model says starts with, or is fed, air.
AIR = {"O2": 0.2095, "CO2": 0.0004, "N2": 1.0 - 0.2095 - 0.0004}

_EXTRA_SPECIES = {
    "H3PO4":  H3PO4,
    "H2PO4-": H2PO4_minus,
    "HPO4--": HPO4_2minus,
    "PO4---": PO4_3minus,
    "HSO4-":  HSO4_minus,
    "SO4--":  SO4_2minus,
    "K+":     K_plus,
    "Cl-":    Cl_minus,
    "Na+":    Na_plus,
    "AceticAcid":    AceticAcid,
    "PropionicAcid": PropionicAcid,
    "ButyricAcid":   ButyricAcid,
    "CitricAcid":    CitricAcid,
    "Yeast":         Yeast,
    "Yeast_CHO":     Yeast_CHO,
    "O2":            O2,
    "N2":            N2,
}

_EXTRA_REACTIONS = ReactionSystem([
    # H3PO4 ⇌ H2PO4- + H+
    EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=H3PO4,      phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=H2PO4_minus, phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=H_plus,      phase="liquid", coefficient=+1.0),
        ],
        log_K=-2.15, T_ref_K=_T_REF_K,
        balance_elements=("P", "H", "O"),
        total_id="phosphate",
        label="eq_phosphate_1",
    ),
    # H2PO4- ⇌ HPO4-- + H+
    EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=H2PO4_minus, phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=HPO4_2minus,  phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=H_plus,       phase="liquid", coefficient=+1.0),
        ],
        log_K=-7.20, T_ref_K=_T_REF_K,
        balance_elements=("P", "H", "O"),
        total_id="phosphate",
        label="eq_phosphate_2",
    ),
    # HPO4-- ⇌ PO4--- + H+
    EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=HPO4_2minus, phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=PO4_3minus,  phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=H_plus,      phase="liquid", coefficient=+1.0),
        ],
        log_K=-12.35, T_ref_K=_T_REF_K,
        balance_elements=("P", "H", "O"),
        total_id="phosphate",
        label="eq_phosphate_3",
    ),
    # HSO4- ⇌ SO4-- + H+
    EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=HSO4_minus, phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=SO4_2minus, phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=H_plus,     phase="liquid", coefficient=+1.0),
        ],
        log_K=-1.99, T_ref_K=_T_REF_K,
        balance_elements=("S", "H", "O"),
        label="eq_bisulfate",
    ),
])

_PARTITION_MODELS = {
    "O2": HenryEquilibrium(H_ref=1.3e-5, dlnH=1500.0),
    "N2": HenryEquilibrium(H_ref=6.4e-6, dlnH=1300.0),
}

BIOPROCESS_BASIC = AQUEOUS_DEFAULT.extend(
    species=_EXTRA_SPECIES,
    reactions=_EXTRA_REACTIONS,
    partition_models=_PARTITION_MODELS,
)
