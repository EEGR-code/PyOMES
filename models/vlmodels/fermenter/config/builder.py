# -*- coding: utf-8 -*-
"""Fluent builder for fermenter construction.

:class:`FermenterBuilder` provides a chainable API that collects
configuration incrementally and produces a
:class:`~PyOMES.core.ControlVolume` (via ``build()``) or a
:class:`~PyOMES.core.simulation.Simulation` (via ``build_simulation()``).

Example
-------
>>> from PyOMES.config.builder import FermenterBuilder
>>> result = (
...     FermenterBuilder()
...     .vessel(V_total_L=2000, T_K=305.15)
...     .gas_feed(vvm_min=1.0, composition={"O2": 0.21, "N2": 0.79})
...     .transfer_kinetic(kLa_O2=150.0)
...     .chemistry(speciation_level=1)
...     .organism("Yeast")
...     .substrate("AceticAcid", mu_max=0.5, Ks=5e-3, yield_gX_gS=0.36)
...     .build_simulation_and_run(tau_h=5.0, n_steps=1000)
... )
>>> result.pH[-1]
"""

from __future__ import annotations

import warnings
from typing import Any, Callable, Dict, List, Optional

from .configs import (
    VesselConfig,
    GasFeedConfig,
    TransferConfig,
    TransferMode,
    SpeciesTransferConfig,
    ChemistryConfig,
    OrganismConfig,
    SubstrateConfig,
    SimulationConfig,
)
from .factory import FermenterFactory
from PyOMES.core.control_volume import ControlVolume


class FermenterBuilder:
    """Fluent builder for constructing a fermenter from incremental calls.

    Each method stores parameters and returns ``self`` for chaining.
    ``build()`` creates the :class:`ControlVolume`.
    ``build_and_run()`` creates the volume and runs a batch.

    Any method can be called in any order; ``build()`` validates
    completeness.
    """

    def __init__(self):
        self._vessel_kw: Dict[str, Any] = {}
        self._gas_feed_kw: Optional[Dict[str, Any]] = None
        self._transfer_cfg: Optional[TransferConfig] = None
        self._chemistry_kw: Dict[str, Any] = {}
        self._organism_kw: Optional[Dict[str, Any]] = None
        self._substrates: List[Dict[str, Any]] = []
        self._controllers: List[Any] = []
        # simulation-class C11: new framework fluent state.
        # _controllers is reused for both the legacy build() path
        # and the new build_simulation() path (Simulation.controllers).
        self._profiles: List[Any] = []
        self._recorder: Optional[Any] = None
        self._reaction_system: Optional[Any] = None
        self._label: str = "fermenter"
        self._solver: Optional[Any] = None

    # ── Solver ─────────────────────────────────────────────────────────

    def solver(
        self,
        solver_type: str = "euler",
        **kwargs,
    ) -> "FermenterBuilder":
        """Set the time-stepping solver.

        Parameters
        ----------
        solver_type : str
            ``"euler"`` — explicit Euler with snapshot + clamping (default).
            ``"dop853"`` — adaptive DOP853 via scipy (high-order explicit).
            ``"radau"`` — implicit Radau solver (for stiff systems like AD).
            ``"bdf"`` — implicit BDF solver (for stiff systems).
        **kwargs
            Forwarded to the solver constructor (e.g. ``rtol``, ``atol``,
            ``max_step`` for adaptive solvers).

        Returns
        -------
        FermenterBuilder
        """
        from PyOMES.core.solvers import SimultaneousEulerSolver
        if solver_type == "euler":
            self._solver = SimultaneousEulerSolver(**kwargs)
        elif solver_type in ("dop853", "radau", "bdf"):
            from PyOMES.core.solvers import SimultaneousAdaptiveSolver
            method_map = {"dop853": "DOP853", "radau": "Radau", "bdf": "BDF"}
            self._solver = SimultaneousAdaptiveSolver(method=method_map[solver_type], **kwargs)
        else:
            raise ValueError(
                f"Unknown solver_type={solver_type!r}; "
                f"expected 'euler', 'dop853', 'radau', or 'bdf'"
            )
        return self

    # ── Vessel ────────────────────────────────────────────────────────

    def vessel(
        self,
        V_total_L: float = 2.0,
        headspace_frac: float = 0.20,
        T_K: float = 305.15,
        P_init_atm: float = 1.0,
        yO2_init: float = 0.2095,
        yCO2_init: float = 0.0004,
        yN2_init: Optional[float] = None,
    ) -> "FermenterBuilder":
        """Set vessel geometry, temperature, and initial gas composition."""
        self._vessel_kw = {
            "V_total_L": V_total_L,
            "headspace_frac": headspace_frac,
            "T_K": T_K,
            "P_init_atm": P_init_atm,
            "yO2_init": yO2_init,
            "yCO2_init": yCO2_init,
            "yN2_init": yN2_init,
        }
        return self

    # ── Gas feed ──────────────────────────────────────────────────────

    def gas_feed(
        self,
        vvm_min: float = 1.0,
        composition: Optional[Dict[str, float]] = None,
        P_inlet_atm: float = 1.0,
    ) -> "FermenterBuilder":
        """Set continuous gas feed (sparging) parameters.

        Use ``vvm_min=0`` for no sparging (e.g. well plate).
        """
        self._gas_feed_kw = {
            "vvm_min": vvm_min,
            "P_inlet_atm": P_inlet_atm,
        }
        if composition is not None:
            self._gas_feed_kw["composition"] = dict(composition)
        return self

    def no_gas_feed(self) -> "FermenterBuilder":
        """Explicitly disable gas feed (well plate, sealed vessel)."""
        self._gas_feed_kw = None
        return self

    # ── Gas-liquid transfer ───────────────────────────────────────────

    def transfer_kinetic(
        self,
        kLa_O2: float = 150.0,
        kLa_CO2_ratio: float = 0.9,
    ) -> "FermenterBuilder":
        """Set kinetic O₂/CO₂ transfer with equilibrium N₂."""
        self._transfer_cfg = TransferConfig.default_kinetic(
            kLa_O2=kLa_O2, kLa_CO2_ratio=kLa_CO2_ratio,
        )
        return self

    def transfer_equilibrium(self) -> "FermenterBuilder":
        """Set all species to instantaneous Henry equilibrium."""
        self._transfer_cfg = TransferConfig.default_equilibrium()
        return self

    def transfer(self, config: TransferConfig) -> "FermenterBuilder":
        """Set a custom TransferConfig directly."""
        self._transfer_cfg = config
        return self

    def transfer_species(
        self,
        species_id: str,
        mode: str = "equilibrium",
        kLa_per_h: float = 0.0,
        henry_mol_L_atm: Optional[float] = None,
    ) -> "FermenterBuilder":
        """Add or update a species in the gas-liquid transfer configuration.

        Call after ``.transfer_equilibrium()`` or ``.transfer_kinetic()``
        to register additional gas species.  For anaerobic digestion,
        this is how to add CH₄ and H₂ to the transfer::

            .transfer_equilibrium()
            .transfer_species("CH4")   # uses Henry constant from lookup
            .transfer_species("H2")

        Parameters
        ----------
        species_id : str
            Gas species identifier (must have a Henry constant in
            ``_HENRY_PARAMS`` unless ``henry_mol_L_atm`` is provided).
        mode : str
            ``"equilibrium"`` or ``"kinetic"`` (default ``"equilibrium"``).
        kLa_per_h : float
            kLa for kinetic mode (1/h).  Ignored for equilibrium mode.
        henry_mol_L_atm : float or None
            Override Henry constant.  If None, looked up from the
            built-in ``_HENRY_PARAMS`` table at build time.

        Returns
        -------
        FermenterBuilder
            self (for chaining).
        """
        if self._transfer_cfg is None:
            self._transfer_cfg = TransferConfig.default_equilibrium()

        mode_lower = mode.strip().lower()
        if mode_lower in ("equilibrium", "eq"):
            tm = TransferMode.EQUILIBRIUM
        elif mode_lower in ("kinetic", "kin", "kla"):
            tm = TransferMode.KINETIC
        else:
            raise ValueError(f"mode must be 'equilibrium' or 'kinetic', got {mode!r}")

        sp_cfg = SpeciesTransferConfig(
            mode=tm,
            kLa_per_h=kLa_per_h if tm == TransferMode.KINETIC else 0.0,
            henry_mol_L_atm=henry_mol_L_atm,
        )
        self._transfer_cfg.species[species_id] = sp_cfg
        return self

    def speciation_correction(
        self,
        species_id: str,
        molecular_key: str,
    ) -> "FermenterBuilder":
        """Deprecated — no-op stub kept for call-site migration.

        chemistry-unification-3b C7: speciation_keys are now auto-derived
        from declared cross-phase equilibrium reactions at CV construction
        via ``derive_speciation_keys()``.  This method is a no-op and will
        be removed in a future cleanup pass.  Callers should declare
        cross-phase reactions in their reaction builder instead.
        """
        return self

    # ── Chemistry ─────────────────────────────────────────────────────

    def chemistry(
        self,
        use_activity: bool = False,
        activity_model: str = "davies",
    ) -> "FermenterBuilder":
        """Set speciation and aqueous chemistry settings.

        pKa values are no longer threaded through ``chemistry()`` —
        they live on declared equilibrium reactions consumed by
        :meth:`BisectionChemicalEquilibriumEngine.from_reactions`. Model builders that
        need equilibria install them on the engine after ``.build()``.
        """
        self._chemistry_kw = {
            "use_activity": use_activity,
            "activity_model": activity_model,
        }
        return self

    # ── Organism ──────────────────────────────────────────────────────

    def organism(
        self,
        organism_id: str = "Yeast",
        atoms: Optional[Dict[str, float]] = None,
        MW: Optional[float] = None,
        balance_basis: str = "CHO",
        n_source_id: str = "NH3",
    ) -> "FermenterBuilder":
        """Set the organism for reaction building."""
        self._organism_kw = {
            "organism_id": organism_id,
            "balance_basis": balance_basis,
            "n_source_id": n_source_id,
        }
        if atoms is not None:
            self._organism_kw["atoms"] = dict(atoms)
        if MW is not None:
            self._organism_kw["MW"] = MW
        return self

    # ── Substrates ────────────────────────────────────────────────────

    def substrate(
        self,
        substrate_id: str = "AceticAcid",
        atoms: Optional[Dict[str, float]] = None,
        MW: Optional[float] = None,
        mu_max: float = 0.5,
        Ks: float = 5e-3,
        yield_gX_gS: float = 0.36,
        kinetics: Optional[Any] = None,
    ) -> "FermenterBuilder":
        """Add a substrate with kinetic parameters.

        Can be called multiple times for multiple substrates.

        Parameters
        ----------
        substrate_id : str
            Chemical identifier (e.g. ``"Glucose"``).
        atoms : dict or None
            Elemental composition.  Looked up from registry if None.
        MW : float or None
            Molecular weight (g/mol).  Looked up from registry if None.
        mu_max : float
            Maximum specific growth rate (1/h).  Used only if
            ``kinetics`` is None (default Monod).
        Ks : float
            Half-saturation constant (g/L for Monod).  Used only if
            ``kinetics`` is None.
        yield_gX_gS : float
            Biomass yield (g_biomass / g_substrate).
        kinetics : GrowthKinetics or None
            Pluggable kinetics model.  If provided, ``mu_max`` and
            ``Ks`` are ignored (the kinetics object owns those
            parameters).  If None, default Monod is used.

            Available models::

                from PyOMES.config.kinetics import (
                    Monod, Contois, Andrews, ContoisAndrews,
                    Tessier, Moser, Blackman, DualSubstrateMonod,
                )

        Returns
        -------
        FermenterBuilder
            self (for chaining).
        """
        kw: Dict[str, Any] = {
            "substrate_id": substrate_id,
            "mu_max": mu_max,
            "Ks": Ks,
            "yield_gX_gS": yield_gX_gS,
        }
        if atoms is not None:
            kw["atoms"] = dict(atoms)
        if MW is not None:
            kw["MW"] = MW
        if kinetics is not None:
            kw["kinetics"] = kinetics
        self._substrates.append(kw)
        return self

    # ── Controllers ───────────────────────────────────────────────────

    def controller(self, ctrl: Any) -> "FermenterBuilder":
        """Add a controller (pressure relief, pH, DO, etc.).

        Can be called multiple times.
        """
        self._controllers.append(ctrl)
        return self

    # ── Profiles (simulation-class C11) ───────────────────────────────

    def profile(self, prof: Any) -> "FermenterBuilder":
        """Add an open-loop profile (e.g. TemperatureRamp, VVMSchedule).

        Consumed by :meth:`build_simulation`. The legacy :meth:`build`
        path doesn't carry profiles — those are a new-framework
        concept introduced in C10.
        """
        self._profiles.append(prof)
        return self

    # ── Recorder (simulation-class C11) ───────────────────────────────

    def recorder(self, rec: Any) -> "FermenterBuilder":
        """Set a custom recorder (default: ``BatchRecorder`` constructed
        per :meth:`Simulation.run` call).

        Consumed only by :meth:`build_simulation`; ignored by the
        legacy :meth:`build` / :meth:`build_and_run` path.
        """
        self._recorder = rec
        return self

    # ── Custom reaction system ────────────────────────────────────────

    def reaction_system(self, system: Any) -> "FermenterBuilder":
        """Set a pre-built reaction system (overrides organism/substrates)."""
        self._reaction_system = system
        return self

    # ── Label ─────────────────────────────────────────────────────────

    def label(self, name: str) -> "FermenterBuilder":
        """Set a human-readable label for the fermenter."""
        self._label = str(name)
        return self

    # ── Build ─────────────────────────────────────────────────────────

    def _build_configs(self):
        """Create config dataclasses from stored parameters.

        Returns (vessel, gas_feed_or_None, transfer, chemistry,
                 organism_or_None, substrates_list).
        """
        # Vessel (required)
        vessel = VesselConfig(**self._vessel_kw) if self._vessel_kw else VesselConfig()

        # Gas feed (optional)
        gas_feed = GasFeedConfig(**self._gas_feed_kw) if self._gas_feed_kw else None

        # Transfer (required, default to equilibrium)
        transfer = self._transfer_cfg or TransferConfig.default_equilibrium()

        # Chemistry
        chem = ChemistryConfig(**self._chemistry_kw) if self._chemistry_kw else ChemistryConfig()

        # Organism
        org = OrganismConfig(**self._organism_kw) if self._organism_kw else None

        # Substrates
        subs = [SubstrateConfig(**kw) for kw in self._substrates] if self._substrates else None

        # Validate: substrates without organism
        if subs and org is None:
            warnings.warn(
                "Substrates defined but no organism set. Reactions will not "
                "be built. Call .organism() before .build().",
                UserWarning,
                stacklevel=3,
            )
            subs = None

        return vessel, gas_feed, transfer, chem, org, subs

    def build(self) -> ControlVolume:
        """Build and return a configured :class:`ControlVolume`.

        The CV's ``phases`` dict is ``{"gas": ..., "liquid": ...}`` with
        a :class:`KineticGasLiquidLink` as an internal interface.

        Raises
        ------
        ValueError
            If the configuration is incomplete or inconsistent.

        Returns
        -------
        ControlVolume

        Notes
        -----
        The class and method name (``FermenterBuilder.build``) are
        Phase 7 holdovers — they predate the deletion of
        ``GasLiquidVolume``.  A cosmetic rename is deferred; bundle
        it with the parallel rename of
        ``FermenterFactory.create_volume`` rather than churning the
        caller surface twice.
        """
        vessel, gas_feed, transfer, chem, org, subs = self._build_configs()

        return FermenterFactory.create_volume(
            vessel=vessel,
            transfer=transfer,
            chemistry=chem,
            organism=org,
            substrates=subs,
            gas_feed=gas_feed,
            reaction_system=self._reaction_system,
            controllers=self._controllers if self._controllers else None,
            label=self._label,
        )

    # ── New-framework build ───────────────────────────────────────────

    def build_simulation(self, *, label: Optional[str] = None) -> Any:
        """Build the configured CV and wrap it in a :class:`Simulation`.

        Parameters
        ----------
        label : str, optional
            Override the simulation label (defaults to the
            ``FermenterBuilder.label`` value).

        Returns
        -------
        PyOMES.core.simulation.Simulation
            A fresh Simulation owning the constructed CV. Controllers,
            profiles, solver, and recorder set via the fluent API are
            attached. Strong-ion seeding (post-STATE_UNIFICATION C4)
            and chem_env (deleted in STATE_UNIFICATION C4e) are no
            longer builder concerns — strong ions live in
            ``cv.phases["liquid"].n_mol`` as Species entries from
            CV construction.

        Notes
        -----
        Use :meth:`build_simulation_and_run` to build and immediately run.
        """
        # Lazy import — PyOMES.core.simulation imports from
        # PyOMES.core.lifecycle, which transitively touches this module
        # via the package init at import time. Lazy resolves the cycle.
        from PyOMES.core.simulation import Simulation

        cv = self.build()
        return Simulation(
            cvs={"main": cv},
            controllers=list(self._controllers) if self._controllers else None,
            profiles=list(self._profiles) if self._profiles else None,
            solver=self._solver,
            recorder=self._recorder,
            label=str(label) if label is not None else self._label,
        )

    def build_simulation_and_run(
        self,
        tau_h: float = 5.0,
        n_steps: int = 1000,
        *,
        label: Optional[str] = None,
    ) -> Any:
        """Build the Simulation and immediately call ``.run(tau_h, n_steps)``.

        Convenience wrapper for the new framework (mirrors
        :meth:`build_and_run` on the legacy path). Returns whatever
        the recorder's ``finalize()`` returns; for the default
        :class:`~PyOMES.core.recorder.BatchRecorder` this is a
        :class:`~PyOMES.core.recorder.BatchResult`.
        """
        sim = self.build_simulation(label=label)
        return sim.run(tau_h=tau_h, n_steps=n_steps)

    # ── Inspection ────────────────────────────────────────────────────

    def get_configs(self):
        """Return the config dataclasses without building.

        Useful for inspection, serialisation, or modification before
        passing to ``FermenterFactory.create_volume()`` directly.

        Returns
        -------
        dict
            ``{vessel, gas_feed, transfer, chemistry, organism, substrates}``.
        """
        vessel, gas_feed, transfer, chem, org, subs = self._build_configs()
        return {
            "vessel": vessel,
            "gas_feed": gas_feed,
            "transfer": transfer,
            "chemistry": chem,
            "organism": org,
            "substrates": subs,
        }

    def __repr__(self):
        parts = []
        if self._vessel_kw:
            parts.append(f"vessel({self._vessel_kw.get('V_total_L', '?')}L)")
        if self._gas_feed_kw:
            parts.append(f"gas_feed(vvm={self._gas_feed_kw.get('vvm_min', '?')})")
        if self._transfer_cfg:
            modes = [f"{sp}:{c.mode.value}" for sp, c in self._transfer_cfg.species.items()]
            parts.append(f"transfer([{', '.join(modes)}])")
        if self._organism_kw:
            parts.append(f"organism({self._organism_kw.get('organism_id', '?')})")
        if self._substrates:
            ids = [s.get("substrate_id", "?") for s in self._substrates]
            parts.append(f"substrates([{', '.join(ids)}])")
        if self._controllers:
            parts.append(f"controllers({len(self._controllers)})")
        desc = ", ".join(parts) if parts else "empty"
        return f"FermenterBuilder({desc})"
