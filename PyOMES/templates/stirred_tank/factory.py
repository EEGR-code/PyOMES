# -*- coding: utf-8 -*-
"""Factory for constructing stirred tanks from config dataclasses.

:class:`StirredTankFactory` consumes the config hierarchy
(:class:`~PyOMES.templates.stirred_tank.VesselConfig`,
:class:`~PyOMES.templates.stirred_tank.TransferConfig`, etc.) and
produces a ready-to-use :class:`~PyOMES.core.ControlVolume` whose
``phases`` dict is ``{"gas": GasPhase, "liquid": LiquidPhase}`` with
gas-liquid transfer configured via the ``transfer_models`` kwarg.

The model's chemistry is whatever it is given: a ``chemistry_db``, its own
``species``, or both. There is no default database. The headspace starts
with the gases in the vessel's ``gas_composition`` (none by default), plus
every transfer species at zero; only the species in ``transfer`` transfer.
Growth reactions take their oxygen, carbon dioxide and water by the
organism's ``o2_id`` / ``co2_id`` / ``h2o_id`` (``"O2"``, ``"CO2"``,
``"H2O"`` unless set); the species passed must define them. The organism
and substrates are named by an id among those species, or defined by an id
with ``atoms`` or a ``Species``, which joins them; a definition that
differs from a species of the same id raises unless ``overwrite=True``.

Example
-------
>>> from PyOMES.templates.stirred_tank import *
>>> from PyOMES.databases.anaerobic_digestion import AD_BASIC
>>> from PyOMES.databases.bioprocess_basic import AIR
>>> from PyOMES.core import Simulation
>>>
>>> cv = StirredTankFactory.create_volume(
...     vessel=VesselConfig(V_total_L=2000, T_K=305.15, gas_composition=AIR),
...     gas_feed=GasFeedConfig(vvm_min=1.0, composition={"O2": 0.21, "N2": 0.79}),
...     transfer=TransferConfig.kinetic({"O2": 150.0, "CO2": 135.0}, equilibrium=["N2"]),
...     organism=OrganismConfig("Yeast"),
...     substrates=[SubstrateConfig("AceticAcid", yield_gX_gS=0.36)],
...     chemistry_db=AD_BASIC,
... )
>>> result = Simulation(cvs={"fermenter": cv}).run(tau_h=5.0, n_steps=1000)
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Sequence

from .configs import (
    VesselConfig,
    GasFeedConfig,
    TransferConfig,
    TransferMode,
    ChemistryConfig,
    OrganismConfig,
    SubstrateConfig,
    SimulationConfig,
)
from PyOMES.core.phases import GasPhase, LiquidPhase
from PyOMES.units import R_L_ATM_PER_MOL_K
from PyOMES.core.control_volume import ControlVolume
from PyOMES.core.transfer_models import KineticTransferModel, EquilibriumTransferModel
from PyOMES.core.boundaries import GasFeed
from PyOMES.chemistry.partition import PartitionModel
from PyOMES.chemistry.species import Species, SpeciesConflictError
from PyOMES.chemistry.species_check import merge_species
from PyOMES.reactions.equilibrium.interphase import HenryEquilibrium
from PyOMES.databases.database import ChemistryDatabase


# ════════════════════════════════════════════════════════════════════════
#  StirredTankFactory
# ════════════════════════════════════════════════════════════════════════

class StirredTankFactory:
    """Factory for constructing a stirred-tank ControlVolume from config dataclasses.

    Currently supports gas+liquid vessels only — always produces a CV
    with ``phases={"gas": ..., "liquid": ...}``. The unqualified name
    doesn't imply a third phase (e.g.
    :class:`~PyOMES.core.phases.SolidPhase`) is wired through today;
    it's chosen so this class doesn't need a second rename if one is
    added later.
    """

    @staticmethod
    def create_volume(
        vessel: VesselConfig,
        transfer: TransferConfig,
        chemistry: Optional[ChemistryConfig] = None,
        organism: Optional[OrganismConfig] = None,
        substrates: Optional[Sequence[SubstrateConfig]] = None,
        gas_feed: Optional[GasFeedConfig] = None,
        reaction_system: Optional[Any] = None,
        controllers: Optional[list] = None,
        chemistry_db: Optional[ChemistryDatabase] = None,
        label: str = "fermenter",
        species: Optional[Any] = None,
    ) -> ControlVolume:
        """Create a configured fermenter ControlVolume from config dataclasses.

        Parameters
        ----------
        vessel : VesselConfig
            Vessel geometry, temperature, and initial gas composition.
        transfer : TransferConfig
            Per-species gas-liquid transfer configuration.
        chemistry : ChemistryConfig, optional
            Speciation engine settings.  If ``None``, uses defaults.
        organism : OrganismConfig, optional
            Organism identity and composition.  Required if building
            reactions from substrates.
        substrates : list of SubstrateConfig, optional
            Substrate definitions with kinetic parameters.  If provided
            (along with ``organism``), a ``ReactionSystem`` is built via
            ``ReactionBuilder.aerobic_growth()``.
        gas_feed : GasFeedConfig, optional
            Sparging parameters.  If ``None`` or ``vvm_min=0``, no gas
            feed boundary is created.
        reaction_system : object, optional
            Pre-built reaction system.  If provided, ``organism``/
            ``substrates`` are ignored for reaction building.
        controllers : list, optional
            Reserved for future use (controller integration).
        chemistry_db : ChemistryDatabase, optional
            The model's chemistry database, if it uses one. It supplies
            :class:`~PyOMES.chemistry.partition.PartitionModel` objects for
            transfer species whose ``henry_mol_L_atm`` is not set in
            ``TransferConfig``, and its ``species`` join the model's species.
            There is no default: a transfer species with neither a Henry
            constant nor a partition model here raises.
        label : str
            Human-readable label.
        species : mapping or iterable of Species, optional
            Species the model defines itself, merged with
            ``chemistry_db.species`` (the same id with different data
            raises :class:`~PyOMES.chemistry.SpeciesConflictError`). Growth
            reactions resolve the organism, substrates, O2, CO2, H2O and
            nitrogen source against these species.

        Returns
        -------
        ControlVolume
            A CV whose ``phases`` dict is ``{"gas": ..., "liquid": ...}``
            with a :class:`KineticGasLiquidLink` as an internal interface.
        """
        chemistry = chemistry or ChemistryConfig()
        T_K = vessel.T_K
        model_species = merge_species(
            chemistry_db.species if chemistry_db is not None else {},
            species if species is not None else {},
        )
        db_partition_models = (
            chemistry_db.partition_models if chemistry_db is not None else {}
        )

        # ── 1. Partition models ────────────────────────────────────────
        partition_models_dict: Dict[str, PartitionModel] = {}
        for sp, sp_cfg in transfer.species.items():
            if sp_cfg.mode == TransferMode.NONE:
                continue
            if sp_cfg.henry_mol_L_atm is not None:
                kH_val = float(sp_cfg.henry_mol_L_atm)
                partition_models_dict[sp] = HenryEquilibrium(
                    H_ref=kH_val * 1000.0 / 101325.0, dlnH=0.0
                )
            elif sp in db_partition_models:
                partition_models_dict[sp] = db_partition_models[sp]
            else:
                raise ValueError(
                    f"No partition model for {sp!r}. Provide henry_mol_L_atm "
                    f"in TransferConfig or add a PartitionModel to chemistry_db."
                    + ("" if chemistry_db is not None
                       else " (no chemistry_db was passed)")
                )

        # ── 2. kLa and equilibrium sets ────────────────────────────────
        kLa_dict: Dict[str, float] = {}
        eq_species: set = set()
        for sp, sp_cfg in transfer.species.items():
            if sp_cfg.mode == TransferMode.EQUILIBRIUM:
                eq_species.add(sp)
            elif sp_cfg.mode == TransferMode.KINETIC:
                kLa_dict[sp] = float(sp_cfg.kLa_per_h)

        # ── 3. Phases ──────────────────────────────────────────────────
        V_gas = vessel.V_headspace_L
        V_liq = vessel.V_liquid_L

        # Initial gas moles from the ideal gas law, split by the vessel's
        # gas composition (normalised by its sum).
        n_total_gas = (vessel.P_init_atm * V_gas) / (R_L_ATM_PER_MOL_K * T_K)
        y_sum = sum(vessel.gas_composition.values())
        gas_n_mol: Dict[str, float] = {
            sp: (n_total_gas * (y / y_sum) if y_sum > 0 else 0.0)
            for sp, y in vessel.gas_composition.items()
        }

        # Ensure all transfer species exist in the gas phase (at zero
        # if not already present).  This allows the gas-liquid link to
        # handle CH₄, H₂, NH₃, H₂S, etc. from the first timestep.
        for sp in partition_models_dict:
            if sp not in gas_n_mol:
                gas_n_mol[sp] = 0.0

        gas_phase = GasPhase(
            n_mol=gas_n_mol,
            V_L=V_gas,
            T_K=T_K,
        )

        # Initial liquid: dissolved gases at Henry equilibrium
        liq_n_mol: Dict[str, float] = {}
        for sp, model in partition_models_dict.items():
            p_i = gas_phase.p_atm.get(sp, 0.0)
            liq_n_mol[sp] = model._kH_mol_L_atm(T_K) * p_i * V_liq
        liquid_phase = LiquidPhase(n_mol=liq_n_mol, V_L=V_liq, T_K=T_K)

        # ── 4. Reaction system ────────────────────────────────────────
        # state-unification C4d: speciation engine attaches to
        # cv.reaction_system; the property_solvers list path is
        # gone. Builders that need an explicit engine config
        # (e.g. BSM2 picks level=1) call reaction_system.attach_engine
        # post-build with their preferred BisectionChemicalEquilibriumEngine
        # construction.
        rxn_system = reaction_system
        if rxn_system is None and organism is not None and substrates:
            rxn_system = StirredTankFactory._build_reaction_system(
                organism, substrates, T_K, model_species,
            )

        # Pre-configure the lazy-engine defaults from the chemistry
        # config when a ReactionSystem is present.
        if rxn_system is not None and hasattr(rxn_system, "configure_engine"):
            rxn_system.configure_engine(activity_model=chemistry.activity_model)

        # ── 6. Boundaries ─────────────────────────────────────────────
        boundaries = []
        if gas_feed is not None and gas_feed.vvm_min > 0:
            boundaries.append(GasFeed(
                vvm_min=gas_feed.vvm_min,
                y=dict(gas_feed.composition),
                P_inlet_atm=gas_feed.P_inlet_atm,
                phase_key="gas",
                liquid_phase_key="liquid",
                label="gas_feed",
            ))

        # ── 7. Assemble ControlVolume ─────────────────────────────────
        transfer_models: Dict[str, Any] = {}
        for sp, pm in partition_models_dict.items():
            if sp in eq_species:
                transfer_models[sp] = EquilibriumTransferModel(pm)
            else:
                transfer_models[sp] = KineticTransferModel(
                    partition_model=pm,
                    k_transfer=float(kLa_dict.get(sp, 0.0)),
                )
        return ControlVolume(
            phases={"gas": gas_phase, "liquid": liquid_phase},
            transfer_models=transfer_models,
            boundaries=list(boundaries),
            reaction_system=rxn_system,
            label=label,
        )

    @staticmethod
    def _build_reaction_system(
        organism: OrganismConfig,
        substrates: Sequence[SubstrateConfig],
        T_K: float,
        model_species: Dict[str, Any],
    ) -> Any:
        """Build a ReactionSystem from organism + substrate configs.

        Uses :meth:`ReactionBuilder.aerobic_growth` for each substrate
        with Monod kinetics as the rate law. The organism and substrates
        are resolved against *model_species* by
        :meth:`_resolve_definition`, which adds new definitions to it. O2,
        CO2, H2O and the nitrogen source are the ``Species`` of the
        organism's ``o2_id`` / ``co2_id`` / ``h2o_id`` / ``n_source_id`` in
        *model_species*.

        Returns
        -------
        ReactionSystem or KineticReaction
            A single KineticReaction if one substrate, or a
            ReactionSystem if multiple.
        """
        from PyOMES.reactions import ReactionSystem, ReactionBuilder, Monod

        org = organism
        biomass = StirredTankFactory._resolve_definition(
            "organism", org.organism, org.atoms, org.MW, org.overwrite,
            model_species,
        )
        resolved_substrates = [
            (sub, StirredTankFactory._resolve_definition(
                "substrate", sub.substrate, sub.atoms, sub.MW, sub.overwrite,
                model_species,
            ))
            for sub in substrates
        ]

        gases = {}
        for role, sp_id in (("o2", organism.o2_id), ("co2", organism.co2_id),
                            ("h2o", organism.h2o_id)):
            if sp_id not in model_species:
                raise ValueError(
                    StirredTankFactory._not_among_species(sp_id, model_species)
                    + f" Aerobic growth needs its "
                    f"{role.upper()} ({role}_id={sp_id!r}): pass "
                    f"Species(id={sp_id!r}, ...) in species=, a chemistry_db "
                    f"that defines it, or another {role}_id."
                )
            gases[role] = model_species[sp_id]

        # N source (for CHNO mode)
        n_source = None
        if org.balance_basis == "CHNO":
            n_source = model_species.get(org.n_source_id)
            if n_source is None:
                raise ValueError(
                    StirredTankFactory._not_among_species(
                        org.n_source_id, model_species
                    )
                    + f" CHNO growth needs its nitrogen source "
                    f"(n_source_id={org.n_source_id!r}): pass "
                    f"Species(id={org.n_source_id!r}, ...) in species=, a "
                    f"chemistry_db that defines it, or another n_source_id."
                )

        org_MW = float(biomass.MW)
        organism_id = biomass.id

        reactions: list = []
        for sub, substrate in resolved_substrates:
            sub_MW = float(substrate.MW)
            substrate_id = substrate.id

            # Build rate function from kinetics object or default Monod
            if sub.kinetics is not None:
                # Pluggable kinetics: the kinetics object builds the rate_fn
                rate_fn = sub.kinetics.make_rate_fn(
                    organism_id=organism_id,
                    substrate_id=substrate_id,
                    MW_organism=org_MW,
                    MW_substrate=sub_MW,
                    yield_gX_gS=float(sub.yield_gX_gS),
                )
            else:
                # Default Monod kinetics (backward compatible)
                rate_fn = Monod(mu_max=float(sub.mu_max), Ks=float(sub.Ks)).make_rate_fn(
                    organism_id=organism_id,
                    substrate_id=substrate_id,
                    MW_organism=org_MW,
                    MW_substrate=sub_MW,
                    yield_gX_gS=float(sub.yield_gX_gS),
                )

            rxn =ReactionBuilder.aerobic_growth(
                substrate, biomass,
                o2=gases["o2"], co2=gases["co2"], h2o=gases["h2o"],
                yield_gX_gS=float(sub.yield_gX_gS),
                rate_fn=rate_fn,
                balance=org.balance_basis,
                n_source=n_source,
                label=f"growth_on_{substrate_id}",
            )
            reactions.append(rxn)

        if len(reactions) == 1:
            return reactions[0]
        return ReactionSystem(reactions, label="aerobic_growth")

    @staticmethod
    def _not_among_species(sp_id: str, model_species: Dict[str, Any]) -> str:
        available = ", ".join(sorted(model_species)) or "no species were passed"
        return (f"{sp_id!r} is not among the species passed to this model "
                f"(available: {available}).")

    @staticmethod
    def _resolve_definition(
        name: str,
        given: Any,
        atoms: Optional[Dict[str, float]],
        MW: Optional[float],
        overwrite: bool,
        model_species: Dict[str, Any],
    ) -> Species:
        """Resolve an organism or substrate against the model's species.

        An id alone must be in *model_species*. A definition (an id with
        *atoms*, MW computed from them unless given, or a ``Species``) is
        added to *model_species*; if a species of that id is already there
        with the same atoms, charge and MW, that species is used; if its
        data differ, the definition replaces it when *overwrite* is true and
        raises :class:`~PyOMES.chemistry.SpeciesConflictError` otherwise.
        """
        if isinstance(given, Species):
            candidate = given
        elif atoms is not None:
            candidate = Species(id=given, atoms=dict(atoms), MW=MW)
        else:
            existing = model_species.get(given)
            if existing is None:
                raise ValueError(
                    StirredTankFactory._not_among_species(given, model_species)
                    + f" Define the {name} with atoms= (e.g. .{name}({given!r}, "
                    f"atoms={{...}})), pass Species(id={given!r}, ...) in "
                    f"species=, or a chemistry_db that defines it."
                )
            return existing

        existing = model_species.get(candidate.id)
        if existing is None or existing is candidate or overwrite:
            model_species[candidate.id] = candidate
            return candidate
        if (dict(existing.atoms) == dict(candidate.atoms)
                and existing.charge == candidate.charge
                and existing.MW == candidate.MW):
            return existing
        details = (
            f"  model's: atoms={dict(existing.atoms)}, charge={existing.charge}, "
            f"MW={existing.MW}\n"
            f"  {name}'s: atoms={dict(candidate.atoms)}, "
            f"charge={candidate.charge}, MW={candidate.MW}"
        )
        raise SpeciesConflictError(
            f"The {name} {candidate.id!r} differs from the model's species of "
            f"that id:\n{details}\nPass the id alone to use the model's "
            f"definition, give the {name} another id, or pass overwrite=True "
            f"to replace the model's.",
            species_id=candidate.id,
            details=details,
        )

