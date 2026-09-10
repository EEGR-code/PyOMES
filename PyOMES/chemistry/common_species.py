# -*- coding: utf-8 -*-
"""Universal inorganic aqueous species declarations.

This module provides module-level :class:`~PyOMES.chemistry.species.Species`
objects for chemicals that appear in essentially every aqueous-chemistry
model: water, hydronium, hydroxide, the carbonate system, and the
ammonium system. Importing from here (rather than redeclaring in each
model file) means cross-file ``ReactionSystem`` compositions share the
same ``Species`` objects by identity, which the
:func:`~PyOMES.chemistry.species_check.check_species_consistency`
utility then has nothing to warn about.

Species IDs match the convention the speciation engine and downstream
consumers already use (``"H+"``, ``"OH-"``, ``"HCO3-"``, ``"CO3--"``,
``"NH4+"``) so engine output keys do not change.

MW values are computed automatically from ``atoms`` via the IUPAC 2021
atomic weight table in :mod:`~PyOMES.chemistry.species`.

Usage
-----
>>> from PyOMES.chemistry.common_species import H2O, CO2, NH3
>>> from PyOMES.chemistry import StoichiometryEntry
>>> entry = StoichiometryEntry(species=CO2, phase="liquid", coefficient=-1.0)

Organic / model-specific species (acetate, propionate, biomass
populations) live in their model's file — they are not "universal".
"""

from __future__ import annotations

from .species import Species


H_plus   = Species(id="H+",    atoms={"H": 1},                  charge=+1)
OH_minus = Species(id="OH-",   atoms={"O": 1, "H": 1},          charge=-1)
H2O      = Species(id="H2O",   atoms={"H": 2, "O": 1},          charge=0)

CO2        = Species(id="CO2",    atoms={"C": 1, "O": 2},              charge=0)
HCO3_minus = Species(id="HCO3-", atoms={"H": 1, "C": 1, "O": 3},      charge=-1)
CO3_2minus = Species(id="CO3--", atoms={"C": 1, "O": 3},               charge=-2)

NH3      = Species(id="NH3",   atoms={"N": 1, "H": 3},          charge=0)
NH4_plus = Species(id="NH4+",  atoms={"N": 1, "H": 4},          charge=+1)

H3PO4       = Species(id="H3PO4",  atoms={"H": 3, "P": 1, "O": 4}, charge=0)
H2PO4_minus = Species(id="H2PO4-", atoms={"H": 2, "P": 1, "O": 4}, charge=-1)
HPO4_2minus = Species(id="HPO4--", atoms={"H": 1, "P": 1, "O": 4}, charge=-2)
PO4_3minus  = Species(id="PO4---", atoms={"P": 1, "O": 4},          charge=-3)

HSO4_minus = Species(id="HSO4-", atoms={"H": 1, "S": 1, "O": 4}, charge=-1)
SO4_2minus = Species(id="SO4--", atoms={"S": 1, "O": 4},          charge=-2)

H2S      = Species(id="H2S", atoms={"H": 2, "S": 1}, charge=0)
HS_minus = Species(id="HS-", atoms={"H": 1, "S": 1}, charge=-1)

K_plus   = Species(id="K+",  atoms={"K":  1},         charge=+1)
Cl_minus = Species(id="Cl-", atoms={"Cl": 1},         charge=-1)
Na_plus  = Species(id="Na+", atoms={"Na": 1},         charge=+1)

Ca_plus_plus = Species(id="Ca++", atoms={"Ca": 1}, charge=+2)
Mg_plus_plus = Species(id="Mg++", atoms={"Mg": 1}, charge=+2)
Zn_plus_plus = Species(id="Zn++", atoms={"Zn": 1}, charge=+2)
Mn_plus_plus = Species(id="Mn++", atoms={"Mn": 1}, charge=+2)
Cu_plus_plus = Species(id="Cu++", atoms={"Cu": 1}, charge=+2)
Co_plus_plus = Species(id="Co++", atoms={"Co": 1}, charge=+2)
MoO4_2minus = Species(id="MoO4--", atoms={"Mo": 1, "O": 4}, charge=-2)
