# -*- coding: utf-8 -*-
"""AD_BASIC: aqueous chemistry for anaerobic digestion gas-liquid systems.

Extends :data:`~PyOMES.chemistry.databases.bioprocess_basic.BIOPROCESS_BASIC`
with a CO₂ gas ⇌ liquid Henry declaration and the H₂S ⇌ HS⁻ acid-base
equilibrium, enabling correct alpha computation in
:class:`~PyOMES.core.gas_liquid_link.KineticGasLiquidLink` for dissolved
sulfide.

The CO₂ partition is declared as a single
:class:`~PyOMES.chemistry.partition.HenryEquilibrium` instance
(``_CO2_HENRY``), used both in ``partition_models`` (feeding
``KineticGasLiquidLink``/``transfer_models=``) and in the reaction list
(feeding ``ChemicalEquilibriumEngine``/``NRChemicalEquilibriumEngine`` as a gas-liquid
``EquilibriumConstraint``) — one declaration, two roles
(EQUILIBRIUM_CONSTRAINT_UNIFICATION CP3). Previously these were two
independently-parameterized objects (a ``HenryPartition`` here, a
separate placeholder ``EquilibriumReaction`` with no ``log_K`` in the
reaction list) that nothing validated agreed with each other.

VFA acid-base equilibria (acetate, propionate, butyrate, valerate) are
not included — VFA gas stripping is negligible at all realistic pH
(Henry constants ~160 000× larger than CO₂) so the alpha correction
has no material effect on model predictions.

Usage::

    from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC
"""
from __future__ import annotations

from ..common_species import H_plus, H2S, HS_minus
from ..partition import HenryEquilibrium
from ...reactions.equilibrium import EquilibriumReaction
from ...reactions.reaction_system import ReactionSystem
from ...reactions.stoichiometry import StoichiometryEntry
from .bioprocess_basic import BIOPROCESS_BASIC

_T_REF_K = 298.15

_EXTRA_SPECIES = {
    "H2S": H2S,
    "HS-": HS_minus,
}

# One instance, two roles: fed into both _PARTITION_MODELS (below) and
# _EXTRA_REACTIONS' flat list — the same H_ref/dlnH back both the
# KineticGasLiquidLink partition_ratio() and the EquilibriumConstraint
# log_K/dH_J_per_mol the speciation engines classify as gas-liquid.
_CO2_HENRY = HenryEquilibrium(
    H_ref=3.4e-4, dlnH=2400.0, gas_species="CO2", liquid_species="CO2",
    label="partition_CO2",
)

_PARTITION_MODELS = {
    "CO2": _CO2_HENRY,
    "CH4": HenryEquilibrium(H_ref=1.4e-5, dlnH=1600.0),
    "H2":  HenryEquilibrium(H_ref=7.7e-6, dlnH=500.0),
    "NH3": HenryEquilibrium(H_ref=5.9e-1, dlnH=4200.0),
    "H2S": HenryEquilibrium(H_ref=1.0e-3, dlnH=2100.0),
}

_EXTRA_REACTIONS = ReactionSystem([
    _CO2_HENRY,  # gas ⇌ liquid — same object as _PARTITION_MODELS["CO2"]
    # H2S ⇌ HS⁻ + H⁺  (pKa = 7.0, van 't Hoff dH = 20 kJ/mol)
    EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=H2S,      phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=HS_minus,  phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=H_plus,    phase="liquid", coefficient=+1.0),
        ],
        log_K=-7.0,
        dH_J_per_mol=20000.0,
        T_ref_K=_T_REF_K,
        balance_elements=("H", "S"),
        label="eq_H2S",
    ),
])

AD_BASIC = BIOPROCESS_BASIC.extend(
    species=_EXTRA_SPECIES,
    reactions=_EXTRA_REACTIONS,
    partition_models=_PARTITION_MODELS,
)
