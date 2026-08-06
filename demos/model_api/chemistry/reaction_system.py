#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Demo: declaring a ReactionSystem — kinetic, equilibrium, cross-phase.

This is the canonical "here are the building blocks of chemistry" demo
for the direct-API layer. It demonstrates each reaction kind the
framework knows about and how a :class:`ReactionSystem` routes them:

- :class:`KineticReaction` — a single integrated reaction with a
  callable rate law. Built here via
  :meth:`ReactionBuilder.aerobic_growth`, which derives the
  stoichiometric coefficients from substrate / biomass elemental
  formulas and a yield.
- :class:`EquilibriumReaction` (single-phase) — an algebraic
  constraint with ``log_K``. Acetic acid dissociation
  (``AceticAcid ⇌ Acetate⁻ + H⁺``) is the example here.
- :class:`EquilibriumReaction` (cross-phase) — a partition
  declaration that spans two phase keys. Gives the gas-liquid link
  its molecular-form mapping (``{"CO2": "CO2"}``) without the
  user wiring it by hand. ``log_K`` is omitted — the partition
  constant lives on the link (Henry's law), not the reaction.

The module exposes three factory functions so other demos
(:mod:`demos.model_api.D2Cworkshop.raw_construction`) can
import the chemistry and plumb it into a Simulation without
redeclaring stoichiometry. Running the module directly prints the
post-bucketing inventory of a :class:`ReactionSystem` constructed
from all three.

Run from the repo root after ``pip install -e .``::

    python demos/model_api/chemistry/reaction_system.py
"""

# Bootstrap so the demo runs when src/ isn't on sys.path yet.
import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_repo_root / "demos"))
import _bootstrap  # noqa: F401, E402

from PyOMES.chemistry import Species  # noqa: E402
from PyOMES.chemistry.common_species import H_plus  # noqa: E402
from PyOMES.reactions import (  # noqa: E402
    EquilibriumReaction,
    KineticReaction,
    ReactionBuilder,
    ReactionSystem,
    StoichiometryEntry,
)


# ── Reaction-level declarations (importable factories) ────────────────


# Acetic acid: HA (neutral) and its conjugate base. Acetate⁻ isn't in
# common_species (which only holds universal inorganics), so it's
# declared here alongside the dissociation it participates in.
ACETIC_ACID = Species(
    id="AceticAcid", atoms={"C": 2, "H": 4, "O": 2}, charge=0, MW=60.052,
)
ACETATE_MINUS = Species(
    id="Acetate-", atoms={"C": 2, "H": 3, "O": 2}, charge=-1, MW=59.044,
)
# CO2 declared locally (rather than imported from common_species) so
# its MW matches the value ReactionBuilder.aerobic_growth uses when
# it constructs its own internal CO2 — without the match, the CV's
# species-consistency check raises on the tiny MW disagreement
# (44.009 vs common_species' 44.01).
# chemistry-unification-3b: phase-agnostic Species IDs — "CO2" is used
# for both gas-phase CO2 and dissolved CO2(aq). The carbonate ladder is
# (CO2, HCO3-, CO3--).
CO2 = Species(id="CO2", atoms={"C": 1, "O": 2}, charge=0, MW=44.009)

# Biomass: a CHO pseudo-molecule "Yeast" with the same elemental
# composition the FermenterBuilder uses by default.
YEAST = Species(
    id="Yeast", atoms={"C": 1, "H": 1.61, "O": 0.56}, charge=0, MW=24.626,
)


def make_aerobic_growth_on_acetate(
    mu_max_per_h: float = 0.5,
    Ks_g_per_L: float = 5e-3,
    yield_gX_gS: float = 0.36,
) -> KineticReaction:
    """Aerobic growth of yeast on acetic acid (Monod kinetics).

    Stoichiometry is derived from elemental balance — the user supplies
    only the substrate / biomass formulas and the yield; CO₂, O₂ and
    H₂O coefficients fall out of the balance.

    The rate function evaluates Monod growth on the substrate mass
    concentration and returns extensive substrate consumption (mol/h),
    matching the convention :meth:`ReactionBuilder.aerobic_growth`
    expects.

    Parameters
    ----------
    mu_max_per_h : float
        Maximum specific growth rate (1/h).
    Ks_g_per_L : float
        Half-saturation constant on substrate mass concentration (g/L).
    yield_gX_gS : float
        Biomass yield (g biomass / g substrate).

    Returns
    -------
    KineticReaction
    """
    MW_S = float(ACETIC_ACID.MW)
    MW_X = float(YEAST.MW)
    Y = float(yield_gX_gS)

    def rate_fn(env):
        # env.concentrations are mol/L; convert substrate to g/L for
        # Monod, biomass to g/L for the specific-rate-to-extensive-rate
        # conversion. Matches FermenterFactory._make_rate_fn closure.
        C_S = env.concentrations.get("AceticAcid", 0.0)
        C_X = env.concentrations.get("Yeast", 0.0)
        S_gL = C_S * MW_S
        X_gL = C_X * MW_X
        if X_gL <= 1e-30 or S_gL <= 0.0:
            return 0.0
        mu = mu_max_per_h * S_gL / (Ks_g_per_L + S_gL)
        # Extensive substrate consumption (mol/h): (mu/Y) * X * V / MW_S
        return (mu / Y) * X_gL / MW_S * env.V_L

    return ReactionBuilder.aerobic_growth(
        substrate_id=ACETIC_ACID.id,
        substrate_atoms=dict(ACETIC_ACID.atoms),
        MW_substrate=MW_S,
        biomass_id=YEAST.id,
        biomass_atoms=dict(YEAST.atoms),
        MW_biomass=MW_X,
        yield_gX_gS=Y,
        rate_fn=rate_fn,
        balance="CHO",
        label="growth_on_AceticAcid",
    )


def make_acetate_dissociation(pKa: float = 4.756) -> EquilibriumReaction:
    """Acetic acid dissociation: HA ⇌ A⁻ + H⁺.

    ``log_K = -pKa`` by the products/reactants convention. The
    speciation engine consumes this as a single-phase algebraic
    constraint and writes back ``n_mol["H+"]``, ``n_mol["Acetate-"]``,
    and ``n_mol["AceticAcid"]`` (the molecular form) after each solve.

    Parameters
    ----------
    pKa : float
        Acid dissociation pKa (default 4.756, the canonical aqueous
        acetic acid value at 25 °C).

    Returns
    -------
    EquilibriumReaction
    """
    return EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(
                species=ACETIC_ACID, phase="liquid", coefficient=-1.0
            ),
            StoichiometryEntry(
                species=ACETATE_MINUS, phase="liquid", coefficient=+1.0
            ),
            StoichiometryEntry(
                species=H_plus, phase="liquid", coefficient=+1.0
            ),
        ],
        log_K=-float(pKa),
        balance_elements=("C", "H", "O"),
        label="eq_AceticAcid",
    )


def make_co2_partition() -> EquilibriumReaction:
    """Cross-phase partition: CO₂(gas) ⇌ CO₂aq(liquid).

    A cross-phase :class:`EquilibriumReaction` is a *partition
    declaration*, not a thermodynamic constraint — it tells a
    :class:`KineticGasLiquidLink` which liquid-phase species id
    corresponds to the molecular form of the gas-phase species.
    ``log_K`` is omitted because the partition constant lives on
    the link (Henry's law), not on the reaction.

    After this reaction is attached to a :class:`ReactionSystem`
    and the system is attached to a :class:`ControlVolume` that
    carries a gas-liquid link, the link's
    ``derive_speciation_keys()`` call (run in
    :meth:`ControlVolume.__init__`) populates
    ``link.speciation_keys["CO2"] = "CO2"`` from this declaration
    (phase-agnostic id: gas and dissolved CO2 share the same Species).

    Returns
    -------
    EquilibriumReaction
    """
    return EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=CO2, phase="gas", coefficient=-1.0),
            StoichiometryEntry(species=CO2, phase="liquid", coefficient=+1.0),
        ],
        balance_elements=("C", "O"),
        label="partition_CO2",
    )


def build_reaction_system() -> ReactionSystem:
    """Assemble all three reactions into a single :class:`ReactionSystem`.

    The system pre-buckets reactions by type at construction. After
    this call, three properties expose the type-specific projections:
    :attr:`~ReactionSystem.kinetic_reactions`,
    :attr:`~ReactionSystem.single_phase_equilibria`,
    :attr:`~ReactionSystem.cross_phase_equilibria`, and
    :attr:`~ReactionSystem.blackbox_models` (empty here — see the
    ``fba/`` sibling demos for that bucket).

    Returns
    -------
    ReactionSystem
    """
    return ReactionSystem(
        [
            make_aerobic_growth_on_acetate(),
            make_acetate_dissociation(),
            make_co2_partition(),
        ],
        label="reaction_system_demo",
    )


# ── Run and report ────────────────────────────────────────────────────


def _format_stoichiometry(rxn) -> str:
    """Render a reaction's stoichiometry as ``-1.0 A (liquid) + ..."""
    parts = []
    for entry in rxn.stoichiometry:
        parts.append(
            f"{entry.coefficient:+.3g} {entry.species.id} ({entry.phase})"
        )
    return "  " + " ".join(parts)


def main() -> None:
    system = build_reaction_system()

    print(
        "Demo: declaring a ReactionSystem from scratch."
    )
    print(
        "  3 reactions declared; the system pre-buckets them by type "
        "at construction."
    )
    print()
    print(f"  {system!r}")
    print()

    print("Kinetic reactions ({}):".format(len(system.kinetic_reactions)))
    for rxn in system.kinetic_reactions:
        print(f"  [{rxn.label}]")
        print(_format_stoichiometry(rxn))
    print()

    print(
        "Single-phase equilibria ({}):".format(
            len(system.single_phase_equilibria)
        )
    )
    for rxn in system.single_phase_equilibria:
        print(f"  [{rxn.label}]  log_K = {rxn.log_K:+.4g}")
        print(_format_stoichiometry(rxn))
    print()

    print(
        "Cross-phase equilibria ({}):".format(
            len(system.cross_phase_equilibria)
        )
    )
    for rxn in system.cross_phase_equilibria:
        log_K_str = (
            f"{rxn.log_K:+.4g}" if rxn.log_K is not None
            else "n/a (partition declaration; constant lives on the link)"
        )
        print(f"  [{rxn.label}]  log_K = {log_K_str}")
        print(_format_stoichiometry(rxn))
    print()

    print("Black-box models ({}):".format(len(system.blackbox_models)))
    if not system.blackbox_models:
        print(
            "  (none — see demos/model_api/chemistry/fba/ for "
            "BlackBoxReactionModel examples)"
        )
    print()

    print(
        "Aggregated species ({}): {}".format(
            len(system.species_ids), ", ".join(system.species_ids)
        )
    )
    print(
        "Aggregated phases ({}): {}".format(
            len(system.phases), ", ".join(system.phases)
        )
    )


if __name__ == "__main__":
    main()
