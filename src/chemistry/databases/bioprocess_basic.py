# -*- coding: utf-8 -*-
"""BIOPROCESS_BASIC: aqueous chemistry extended for aerobic bioprocesses.

Extends :data:`AQUEOUS_DEFAULT` with:

- Phosphate and bisulfate equilibrium ladders.

The lump salt species :data:`NH4Cl`, :data:`KH2PO4`, and :data:`NaOH` are
defined as :class:`~PyOMES.chemistry.Species` objects so that their molar
masses are available for medium-recipe calculations, but they carry no
dissolution reactions in this database.  Initial medium compositions should
be specified in ionic form directly (``NH4+`` + ``Cl-``,
``K+`` + ``H2PO4-``, etc.).  NaOH pH correction is handled by
:meth:`~PyOMES.core.ControlVolume.equilibrate_to_pH` via the built-in
strong-corrector map, which adds ``Na+`` to the charge balance.

Usage::

    from PyOMES.chemistry.databases.bioprocess_basic import BIOPROCESS_BASIC
    from PyOMES.chemistry.databases.bioprocess_basic import NH4Cl, KH2PO4, NaOH
"""
from __future__ import annotations

from ..common_species import (
    H_plus,
    H3PO4, H2PO4_minus, HPO4_2minus, PO4_3minus,
    HSO4_minus, SO4_2minus,
    K_plus, Cl_minus, Na_plus,
)
from ..species import Species
from ..partition import HenryEquilibrium
from ...reactions.equilibrium import EquilibriumReaction
from ...reactions.reaction_system import ReactionSystem
from ...reactions.stoichiometry import StoichiometryEntry
from .aqueous import AQUEOUS_DEFAULT

_T_REF_K = 298.15

# ---------------------------------------------------------------------------
# Undissociated salt / corrector species
# ---------------------------------------------------------------------------
NH4Cl  = Species(id="NH4Cl",  atoms={"N": 1, "H": 4, "Cl": 1},       charge=0)
KH2PO4 = Species(id="KH2PO4", atoms={"K": 1, "H": 2, "P": 1, "O": 4}, charge=0)
NaOH   = Species(id="NaOH",   atoms={"Na": 1, "O": 1, "H": 1},        charge=0)

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
