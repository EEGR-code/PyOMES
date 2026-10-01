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
Growth reactions use ``"O2"``, ``"CO2"`` and ``"H2O"``; the species passed
must define them.

Example
-------
>>> from PyOMES.templates.stirred_tank import *
>>> from PyOMES.databases.anaerobic_digestion import AD_BASIC
>>> from PyOMES.databases.bioprocess_basic import AIR
>>> from PyOMES.core import Simulation
>>>
>>> cv = StirredTankFactory.create_volume(
...     vessel=VesselConfig(V_total_L=2000, T_K=305.15, gas_composition=AIR),
...     gas_feed=GasFeedConfig(vvm_min=1.0),
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
            reactions take O2, CO2 and H2O from these species by id.

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
        with Monod kinetics as the rate law. O2, CO2 and H2O are the
        ``Species`` of those ids in *model_species*.

        Returns
        -------
        ReactionSystem or KineticReaction
            A single KineticReaction if one substrate, or a
            ReactionSystem if multiple.
        """
        from PyOMES.chemistry.species import Species
        from PyOMES.reactions import ReactionSystem, ReactionBuilder, Monod
        from PyOMES.compounds import ChemicalRegistry

        registry = ChemicalRegistry.default()

        gases = {}
        for sp_id in ("O2", "CO2", "H2O"):
            if sp_id not in model_species:
                available = (", ".join(sorted(model_species))
                             or "no species were passed")
                raise ValueError(
                    f"{sp_id!r} is not among the species passed to this model "
                    f"(available: {available}). Aerobic growth needs O2, CO2 "
                    f"and H2O: pass Species(id={sp_id!r}, ...) in species=, or "
                    "a chemistry_db that defines it."
                )
            gases[sp_id] = model_species[sp_id]

        org = organism.resolve(registry)
        org_atoms = dict(org.atoms)
        org_MW = float(org.MW)
        biomass = Species(id=org.organism_id, atoms=org_atoms, MW=org_MW)

        reactions: list = []
        for sub_cfg in substrates:
            sub = sub_cfg.resolve(registry)
            sub_atoms = dict(sub.atoms)
            sub_MW = float(sub.MW)
            organism_id = org.organism_id
            substrate_id = sub.substrate_id
            substrate = Species(id=substrate_id, atoms=sub_atoms, MW=sub_MW)

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

            # N source (for CHNO mode)
            n_source = None
            if org.balance_basis == "CHNO":
                n_chem = registry[org.n_source_id]
                n_source = Species(id=org.n_source_id, atoms=dict(n_chem.atoms))

            rxn = ReactionBuilder.aerobic_growth(
                substrate, biomass,
                o2=gases["O2"], co2=gases["CO2"], h2o=gases["H2O"],
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

