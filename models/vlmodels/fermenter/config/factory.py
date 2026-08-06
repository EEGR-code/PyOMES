# -*- coding: utf-8 -*-
"""Factory for constructing fermenters from config dataclasses.

:class:`FermenterFactory` consumes the config hierarchy
(:class:`~fermenter.config.VesselConfig`,
:class:`~fermenter.config.TransferConfig`, etc.) and produces a
ready-to-use :class:`~PyOMES.core.ControlVolume` whose ``phases`` dict
is ``{"gas": GasPhase, "liquid": LiquidPhase}`` with gas-liquid
transfer configured via the ``transfer_models`` kwarg.

Example
-------
>>> from PyOMES.config import *
>>> from PyOMES.core import Simulation
>>>
>>> cv = FermenterFactory.create_volume(
...     vessel=VesselConfig(V_total_L=2000, T_K=305.15),
...     gas_feed=GasFeedConfig(vvm_min=1.0),
...     transfer=TransferConfig.default_kinetic(kLa_O2=150.0),
...     chemistry=ChemistryConfig(speciation_level=1),
...     organism=OrganismConfig("Yeast"),
...     substrates=[SubstrateConfig("AceticAcid", yield_gX_gS=0.36)],
... )
>>> result = Simulation(cvs={"fermenter": cv}).run(tau_h=5.0, n_steps=1000)
"""

from __future__ import annotations

import warnings
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
from PyOMES.core.phases import GasPhase, LiquidPhase, R_L_ATM_MOL_K
from PyOMES.core.control_volume import ControlVolume
from PyOMES.core.transfer_models import KineticTransferModel, EquilibriumTransferModel
from PyOMES.core.boundaries import GasFeed
from PyOMES.chemistry.partition import HenryEquilibrium, PartitionModel
from PyOMES.chemistry.database import ChemistryDatabase
from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC


# ════════════════════════════════════════════════════════════════════════
#  FermenterFactory
# ════════════════════════════════════════════════════════════════════════

class FermenterFactory:
    """Factory for constructing a fermenter ControlVolume from config dataclasses."""

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
            Database supplying :class:`~PyOMES.chemistry.partition.PartitionModel`
            objects for species whose ``henry_mol_L_atm`` is not set in
            ``TransferConfig``.  Defaults to
            :data:`~PyOMES.chemistry.databases.anaerobic_digestion.AD_BASIC`.
        label : str
            Human-readable label.

        Returns
        -------
        ControlVolume
            A CV whose ``phases`` dict is ``{"gas": ..., "liquid": ...}``
            with a :class:`KineticGasLiquidLink` as an internal interface.

        Notes
        -----
        The ``create_volume`` name is a Phase 7 holdover — it returned a
        ``GasLiquidVolume`` before that class was deleted, and now
        returns a ``ControlVolume``.  A cosmetic rename (e.g.
        ``create_cv`` / ``create_fermenter``) is deferred because the
        caller surface is broad (every system script, every test).
        Bundle it with a future builder/factory naming sweep rather
        than touching it standalone.
        """
        chemistry = chemistry or ChemistryConfig()
        chemistry_db = chemistry_db or AD_BASIC
        T_K = vessel.T_K

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
            elif sp in chemistry_db.partition_models:
                partition_models_dict[sp] = chemistry_db.partition_models[sp]
            else:
                raise ValueError(
                    f"No partition model for {sp!r}. Provide henry_mol_L_atm "
                    f"in TransferConfig or add a PartitionModel to chemistry_db."
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

        # Initial gas moles from ideal gas law
        n_total_gas = (vessel.P_init_atm * V_gas) / (R_L_ATM_MOL_K * T_K)
        y_sum = vessel.yO2_init + vessel.yCO2_init + vessel.yN2_init
        if y_sum > 0:
            yO2 = vessel.yO2_init / y_sum
            yCO2 = vessel.yCO2_init / y_sum
            yN2 = vessel.yN2_init / y_sum
        else:
            yO2, yCO2, yN2 = 0.0, 0.0, 0.0

        gas_n_mol: Dict[str, float] = {
            "O2": n_total_gas * yO2,
            "CO2": n_total_gas * yCO2,
            "N2": n_total_gas * yN2,
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
            rxn_system = FermenterFactory._build_reaction_system(
                organism, substrates, T_K,
            )

        # Pre-configure the lazy-engine defaults from the chemistry
        # config when a ReactionSystem is present.
        if rxn_system is not None and hasattr(rxn_system, "configure_engine"):
            rxn_system.configure_engine(
                use_activity=chemistry.use_activity,
                activity_model=chemistry.activity_model,
            )

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
    ) -> Any:
        """Build a ReactionSystem from organism + substrate configs.

        Uses :meth:`ReactionBuilder.aerobic_growth` for each substrate
        with Monod kinetics as the rate law.

        Returns
        -------
        ReactionSystem or KineticReaction
            A single KineticReaction if one substrate, or a
            ReactionSystem if multiple.
        """
        from PyOMES.reactions import ReactionSystem, ReactionBuilder
        from PyOMES.chemistry.compounds import ChemicalRegistry

        registry = ChemicalRegistry.default()

        org = organism.resolve(registry)
        org_atoms = dict(org.atoms)
        org_MW = float(org.MW)

        reactions: list = []
        for sub_cfg in substrates:
            sub = sub_cfg.resolve(registry)
            sub_atoms = dict(sub.atoms)
            sub_MW = float(sub.MW)
            organism_id = org.organism_id
            substrate_id = sub.substrate_id

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
                mu_max = float(sub.mu_max)
                Ks = float(sub.Ks)
                MW_s = sub_MW

                def _make_rate_fn(mu_m, ks, mw_s, mw_x, org_id, sub_id, Y):
                    """Create a Monod rate closure with captured parameters."""
                    def rate_fn(env):
                        C_S = env.concentrations.get(sub_id, 0.0)  # mol/L
                        C_X = env.concentrations.get(org_id, 0.0)  # mol/L
                        S_gL = C_S * mw_s  # g/L
                        X_gL = C_X * mw_x  # g/L
                        if X_gL <= 1e-30 or S_gL <= 0.0:
                            return 0.0
                        mu = mu_m * S_gL / (ks + S_gL) if (ks + S_gL) > 0 else 0.0
                        return (mu / Y) * X_gL / mw_s * env.V_L
                    return rate_fn

                rate_fn = _make_rate_fn(
                    mu_max, Ks, sub_MW, org_MW,
                    organism_id, substrate_id, float(sub.yield_gX_gS),
                )

            # N source atoms (for CHNO mode)
            n_source_atoms = None
            if org.balance_basis == "CHNO":
                try:
                    n_chem = registry[org.n_source_id]
                    n_source_atoms = dict(getattr(n_chem, "atoms", {}) or {})
                except (KeyError, AttributeError) as e:
                    warnings.warn(
                        f"Could not resolve N source '{org.n_source_id}' "
                        f"for {organism_id}: {e}. CHNO stoichiometry may "
                        f"be incomplete.",
                        RuntimeWarning, stacklevel=2,
                    )

            rxn = ReactionBuilder.aerobic_growth(
                substrate_id=substrate_id,
                substrate_atoms=sub_atoms,
                MW_substrate=sub_MW,
                biomass_id=organism_id,
                biomass_atoms=org_atoms,
                MW_biomass=org_MW,
                yield_gX_gS=float(sub.yield_gX_gS),
                rate_fn=rate_fn,
                balance=org.balance_basis,
                n_source_id=org.n_source_id,
                n_source_atoms=n_source_atoms,
                label=f"growth_on_{substrate_id}",
            )
            reactions.append(rxn)

        if len(reactions) == 1:
            return reactions[0]
        return ReactionSystem(reactions, label="aerobic_growth")

