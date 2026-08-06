# -*- coding: utf-8 -*-
"""Gas-liquid transfer link with per-species kinetic/equilibrium modes.

:class:`KineticGasLiquidLink` connects a gas-phase CV to a liquid-phase CV
and computes mass transfer fluxes for each dissolved gas species.

.. note::
    Direct construction is deprecated — use the ``transfer_models=`` kwarg
    on :class:`~PyOMES.core.control_volume.ControlVolume` instead (see
    ``docs/design/TRANSFER_MODEL.md``); it builds this link internally.
    This module's mechanics and ``partition_models=`` contract are
    unchanged and remain the thing ``transfer_models=`` drives underneath.

Two transfer modes are supported per species:

**Kinetic mode** — kLa-limited transfer, where the flux is driven by the
difference between the Henry equilibrium saturation concentration and the
actual dissolved concentration::

    flux_i = kLa_i × V_liq × (kH_i × p_i  −  C_i)

**Equilibrium mode** — instantaneous Henry partition, where the flux drives
the system to exact gas-liquid equilibrium within the timestep.

For species with acid-base equilibria in the liquid phase (CO₂, NH₃, H₂S,
VFAs), only the molecular (un-ionised) form participates in gas-liquid
transfer.  The effective Henry constant is corrected by the molecular
fraction::

    kH_eff = kH / alpha

where ``alpha = n_mol[mol_key] / sum(n_mol[s] for s in ladder)`` is
computed inline from the liquid phase's ``n_mol`` dict — populated by the
speciation engine's writeback during each ``advance()`` step.  The
``speciation_keys`` mapping (``{gas_species: liquid_mol_key}``, e.g.
``{"CO2": "CO2"}``) identifies which ``n_mol`` key holds the molecular
form.  Ladder membership is derived at CV construction from the declared
single-phase equilibrium reactions via
:meth:`~PyOMES.core.gas_liquid_link.KineticGasLiquidLink.derive_speciation_keys`
and stored in :attr:`speciation_ladders`.

Sign convention: positive flux = gas → liquid (absorption);
negative flux = liquid → gas (stripping).

The class implements two call protocols:

* **``compute_flow(cvs, dt_h)``** — primary path.  Accepts a ``{cv_key:
  ControlVolume}`` dict; used by :class:`~PyOMES.core.simulation.Simulation`
  and external callers.
* **``compute_flux(state_a, state_b, dt_h)``** — PhaseInterface adapter.
  Lets the link sit as an internal interface inside a single
  :class:`~PyOMES.core.control_volume.ControlVolume` that holds both gas and
  liquid phases.  Builds a synthetic ``cvs`` dict and delegates to
  ``compute_flow``.

Example
-------
One ``HenryEquilibrium`` instance per species serves both roles:
``partition_models=`` here, and (for species whose gas-liquid partition
also needs to be visible to the speciation engine, e.g. because a
downstream acid-base ladder depends on the molecular fraction) the same
instance goes into the CV's ``ReactionSystem`` reaction list too —
``EQUILIBRIUM_CONSTRAINT_UNIFICATION`` CP3. There is no separate,
independently-parameterized declaration to keep in sync; species without
that need (O2, N2, CH4, H2, NH3 below) simply omit ``gas_species``/
``liquid_species`` and stay ``PartitionModel``-only.

>>> from PyOMES.chemistry import HenryEquilibrium
>>> from PyOMES.core.gas_liquid_link import KineticGasLiquidLink
>>> co2_henry = HenryEquilibrium(
...     H_ref=3.4e-4, dlnH=2400.0, gas_species="CO2", liquid_species="CO2",
... )  # also passed into the CV's ReactionSystem reaction list
>>> link = KineticGasLiquidLink(
...     gas_cv_key="headspace", gas_phase_key="gas",
...     liquid_cv_key="liquid",  liquid_phase_key="liquid",
...     partition_models={
...         "O2": HenryEquilibrium(H_ref=1.3e-5, dlnH=1500.0),
...         "CO2": co2_henry,
...         "N2": HenryEquilibrium(H_ref=6.4e-6, dlnH=1300.0),
...         "CH4": HenryEquilibrium(H_ref=1.4e-5, dlnH=1600.0),
...         "H2":  HenryEquilibrium(H_ref=7.7e-6, dlnH=500.0),
...         "NH3": HenryEquilibrium(H_ref=5.9e-1, dlnH=4200.0),
...     },
...     kLa={"O2": 150.0, "CO2": 135.0},
...     equilibrium_species={"N2", "CH4", "H2", "NH3"},
...     speciation_keys={"CO2": "CO2", "NH3": "NH3"},
... )
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Sequence, Union

from .phases import R_L_ATM_MOL_K
from ..control.descriptors import MutableDict as _MutableDict
from ..chemistry.partition import PartitionModel



@dataclass
class KineticGasLiquidLink:
    """Gas-liquid transfer with per-species kinetic or equilibrium mode.

    Parameters
    ----------
    gas_cv_key : str
        Key of the gas-phase ControlVolume in the system.
    gas_phase_key : str
        Phase key within the gas CV (typically ``"gas"``).
    liquid_cv_key : str
        Key of the liquid-phase ControlVolume.
    liquid_phase_key : str
        Phase key within the liquid CV (typically ``"liquid"``).
    partition_models : dict
        ``{species_id: PartitionModel}`` for each transferable species.
        Species not listed are not transferred.  Use
        :class:`~PyOMES.chemistry.HenryEquilibrium` for Henry-law systems.
    kLa : dict
        Volumetric mass transfer coefficients ``{species_id: 1/h}`` for
        kinetic-mode species.  Species listed in ``equilibrium_species``
        ignore their kLa value.
    equilibrium_species : set or None
        Species to be partitioned at instantaneous Henry equilibrium
        (no kLa limitation).  Default: empty (all species are kinetic).
    speciation_keys : dict
        Mapping from gas species to the speciation result key for
        the molecular (volatile) form.  For species with acid-base
        equilibria, only the undissociated/molecular form participates
        in gas-liquid transfer.  The effective Henry constant is::

            kH_eff = kH / α₀

        where α₀ = C_molecular / C_total, read from the speciation
        result using the provided key.

        Default: ``{"CO2": "CO2"}``.

        Common entries for AD modelling::

            {
                "CO2": "CO2",            # dissolved CO2 / TIC
                "NH3": "NH3",            # molecular NH₃ / total NH
                "H2S": "H2S_HA",         # molecular H₂S / total sulfide
                "AceticAcid": "AceticAcid_HA",
                "Propionate": "Propionate_HA",
                "Butyrate": "Butyrate_HA",
                "Valerate": "Valerate_HA",
            }

    _label : str
        Human-readable label.

    molecular_driving_force : set
        Species for which the kinetic transfer driving force is computed
        from the **molecular** (un-ionised) concentration rather than
        the total dissolved concentration.  This matches the BSM2
        formulation where kLa is applied to dissolved CO₂(aq) only::

            flux = kLa × V × (kH × p_i  −  C_molecular)
                 = kLa × V × (kH × p_i  −  f_mol × C_total)

        rather than the default (total-IC) formulation::

            flux = kLa × V × (kH_eff × p_i  −  C_total)
                 = kLa × V × (kH/f_mol × p_i  −  C_total)

        Both reach the same equilibrium, but the molecular formulation
        has a relaxation rate ~f_mol× slower (matching BSM2 kinetics).
        Default: empty set.  Add ``"CO2"`` for BSM2-style CO₂ transfer.

    Notes
    -----
    The speciation correction uses operator splitting: the molecular
    fraction α₀ is read from the **previous** timestep's speciation
    result.  This is accurate when the timestep is small relative to
    pH changes — the same approach used for CO₂ since the original
    implementation.
    """

    gas_cv_key: str
    gas_phase_key: str
    liquid_cv_key: str
    liquid_phase_key: str
    partition_models: Dict[str, PartitionModel] = field(default_factory=dict)
    kLa: Dict[str, float] = _MutableDict(value_type=float)
    equilibrium_species: Set[str] = field(default_factory=set)
    speciation_keys: Dict[str, str] = field(
        default_factory=lambda: {"CO2": "CO2"},
    )
    # Populated by derive_speciation_keys: gas_id → ordered list of
    # all liquid ladder Species ids for alpha denominator computation.
    speciation_ladders: Dict[str, List[str]] = field(default_factory=dict)
    molecular_driving_force: Set[str] = field(default_factory=set)
    _label: str = ""
    _validated: bool = field(default=False, repr=False)
    # C6 lifecycle gating: set by Simulation.__init__; None when
    # unowned. Lockable mutators consult this via raise_if_running.
    # Replaced the legacy _simulation_running bool + _enter_simulation
    # / _exit_simulation / _warn_if_simulating trio (deleted at C6).
    _context: Any = field(default=None, repr=False)

    def __post_init__(self) -> None:
        import warnings
        warnings.warn(
            "KineticGasLiquidLink is deprecated and will be removed in a future release. "
            "Use the transfer_models kwarg on ControlVolume instead: "
            "ControlVolume(phases=..., transfer_models={'CO2': KineticTransferModel(...), ...}). "
            "See docs/design/TRANSFER_MODEL.md.",
            DeprecationWarning,
            stacklevel=2,
        )

    # ── CVLink protocol properties ────────────────────────────────────

    @property
    def source_cv_key(self) -> str:
        return self.gas_cv_key

    @property
    def source_phase_key(self) -> str:
        return self.gas_phase_key

    @property
    def sink_cv_key(self) -> str:
        return self.liquid_cv_key

    @property
    def sink_phase_key(self) -> str:
        return self.liquid_phase_key

    @property
    def label(self) -> str:
        if self._label:
            return self._label
        return (f"gl_transfer:{self.gas_cv_key}.{self.gas_phase_key}"
                f"->{self.liquid_cv_key}.{self.liquid_phase_key}")

    # ── PhaseInterface protocol properties ────────────────────────────
    #
    # Lets a single ControlVolume hold both gas and liquid phases and use
    # this link as an internal interface — see CV_UPDATE.md §"Phase 3".
    # The keys identify which entries of ``cv.phases`` to pass as
    # ``state_a`` / ``state_b`` to ``compute_flux``.

    @property
    def phase_a_key(self) -> str:
        return self.gas_phase_key

    @property
    def phase_b_key(self) -> str:
        return self.liquid_phase_key

    # ── Chemistry derivation ──────────────────────────────────────────

    def derive_speciation_keys(self, reaction_system) -> None:
        """Populate :attr:`speciation_keys` from declared cross-phase
        equilibrium reactions.

        Reads ``reaction_system.cross_phase_equilibria`` directly —
        the :class:`~PyOMES.reactions.reaction_system.ReactionSystem`
        has already pre-bucketed cross-phase reactions at
        construction. For each cross-phase reaction whose
        stoichiometry spans this link's gas and liquid phase keys,
        records the implied gas → liquid molecular-form mapping. A
        cross-phase reaction like::

            EquilibriumReaction(
                stoichiometry=[
                    StoichiometryEntry(CO2_gas, "gas", -1),
                    StoichiometryEntry(CO2aq,   "liquid", +1),
                ],
            )

        yields ``self.speciation_keys["CO2"] = "CO2aq"`` (the gas
        entry's species id keys the liquid entry's species id).

        The method is additive — it leaves any pre-existing entries
        in ``speciation_keys`` alone unless overridden by a declared
        reaction. Reactions whose phase keys do not match this
        link's ``gas_phase_key`` / ``liquid_phase_key`` are skipped
        (defensive — the same reaction system may serve multiple
        links in a multi-CV system, each with its own phase wiring).

        Called automatically by
        :meth:`~PyOMES.core.control_volume.ControlVolume.__init__`
        after the reaction system is attached.

        Parameters
        ----------
        reaction_system : ReactionSystem
            The CV's reaction system. Objects that don't expose
            ``cross_phase_equilibria`` (e.g. a bare ``KineticReaction``
            attached directly for test/fixture convenience) are
            silently ignored — they have no partition declarations
            to contribute.
        """
        cross_phase = getattr(reaction_system, "cross_phase_equilibria", None)
        if cross_phase is None:
            # Not a ReactionSystem (e.g. a bare KineticReaction
            # attached directly via test fixtures). No partition
            # declarations possible — return without touching
            # ``speciation_keys``.
            return

        for rxn in cross_phase:
            gas_entries = [
                e for e in rxn.stoichiometry
                if e.phase == self.gas_phase_key
            ]
            liq_entries = [
                e for e in rxn.stoichiometry
                if e.phase == self.liquid_phase_key
            ]
            if len(gas_entries) != 1 or len(liq_entries) != 1:
                continue
            self.speciation_keys[gas_entries[0].species.id] = (
                liq_entries[0].species.id
            )

        # Option A — seed an identity mapping for any partition_models species
        # that has no entry yet (neither pre-set nor from a cross-phase reaction
        # above).  Under the phase-agnostic species-ID convention the gas-phase
        # and liquid-phase IDs are always the same string, so a cross-phase
        # declaration for Henry-partitioned species carries no information
        # beyond "sp → sp". Seeding here removes the boilerplate requirement
        # for callers to declare a routing-only reaction just to enable the
        # BFS in the ladder-building block below.
        #
        # Option B (the unified GasLiquidPartition type bundling the
        # thermodynamic constant and the routing in one object, eliminating
        # the conceptual split between partition_models and speciation_keys)
        # shipped as HenryEquilibrium/RaoultEquilibrium
        # (EQUILIBRIUM_CONSTRAINT_UNIFICATION CP1/CP3): the same instance in
        # partition_models can also be placed directly in the reaction list,
        # in which case it's picked up by the cross_phase loop above instead
        # of falling through to this identity-seeding path. This seeding
        # remains for species that only need the PartitionModel role (no
        # gas_species/liquid_species set) and for legacy cross-phase
        # EquilibriumReaction declarations.
        for sp in self.partition_models:
            if sp not in self.speciation_keys:
                self.speciation_keys[sp] = sp

        # Build speciation_ladders from single-phase equilibria.
        # For each gas → liq_id pair, walk the equilibrium reaction
        # graph to find all Species in the same acid-base ladder as
        # liq_id.  Species connected through shared reactions are
        # BFS-expanded; H+, OH-, H2O and solvents are excluded.
        # Ladder is ordered by charge descending (most protonated first).
        single_phase = getattr(
            reaction_system, "single_phase_equilibria", []
        )
        _SKIP = {"H+", "OH-", "H2O"}
        # Build a fast lookup: species_id → list of reactions it appears in
        sp_to_rxns: Dict[str, list] = {}
        for rxn in single_phase:
            for e in rxn.stoichiometry:
                if e.species.id not in _SKIP:
                    sp_to_rxns.setdefault(e.species.id, []).append(rxn)

        for gas_id, liq_id in list(self.speciation_keys.items()):
            if liq_id in self.speciation_ladders:
                continue  # already derived
            # BFS from liq_id through connected equilibria
            ladder_sp: Dict[str, Any] = {}  # id → Species
            queue = [liq_id]
            visited_ids: set = set()
            while queue:
                current_id = queue.pop()
                if current_id in visited_ids:
                    continue
                visited_ids.add(current_id)
                for rxn in sp_to_rxns.get(current_id, []):
                    for e in rxn.stoichiometry:
                        sid = e.species.id
                        if sid not in _SKIP and sid not in visited_ids:
                            ladder_sp[sid] = e.species
                            queue.append(sid)
                # Ensure liq_id itself is in the ladder even if it has no
                # single-phase reactions (graceful degradation).
                if current_id == liq_id and liq_id not in ladder_sp:
                    # find the Species object from any reaction
                    for rxn in sp_to_rxns.get(liq_id, []):
                        for e in rxn.stoichiometry:
                            if e.species.id == liq_id:
                                ladder_sp[liq_id] = e.species
                                break
                        if liq_id in ladder_sp:
                            break

            if not ladder_sp:
                # No single-phase equilibria connect to liq_id — one-
                # element ladder so alpha = 1.0 gracefully.
                self.speciation_ladders[gas_id] = [liq_id]
            else:
                # Sort by charge descending (most protonated first).
                sorted_sps = sorted(
                    ladder_sp.values(),
                    key=lambda sp: -int(sp.charge),
                )
                self.speciation_ladders[gas_id] = [
                    sp.id for sp in sorted_sps
                ]

        self._validated = False

    # ── Configuration validation ────────────────────────────────────────

    def validate(self) -> list:
        """Check for inconsistent transfer link configuration.

        Detects three classes of misconfiguration:

        1. **Equilibrium/kinetic overlap**: species present in both
           ``equilibrium_species`` and ``kLa``.  Equilibrium silently
           takes precedence — the kLa value is ignored.
        2. **Orphan kLa entries**: species with a kLa value but no
           matching Henry constant.  These are never used because
           ``compute_flow`` iterates over ``henry`` keys.
        3. **Orphan molecular_driving_force entries**: species in
           ``molecular_driving_force`` but not in
           ``speciation_keys``.  The molecular driving force requires
           an alpha lookup via ``PropertyResult.alphas[mol_key]``.

        Returns
        -------
        list of str
            Warning messages (empty if consistent).
        """
        import warnings
        issues = []

        # 1. Equilibrium/kinetic overlap
        overlap = self.equilibrium_species & set(self.kLa.keys())
        if overlap:
            msg = (f"Gas transfer: species {overlap} appear in both "
                   f"equilibrium_species and kLa. Equilibrium mode takes "
                   f"precedence — kLa values are silently ignored.")
            issues.append(msg)
            warnings.warn(msg, UserWarning, stacklevel=2)

        # 2. Orphan kLa entries (no matching partition model)
        orphan_kLa = set(self.kLa.keys()) - set(self.partition_models.keys())
        if orphan_kLa:
            msg = (f"Gas transfer: kLa entries {orphan_kLa} have no "
                   f"matching partition model in link.partition_models "
                   f"and will be ignored. Check species names.")
            issues.append(msg)
            warnings.warn(msg, UserWarning, stacklevel=2)

        # 3. Orphan molecular_driving_force entries
        orphan_mdf = (self.molecular_driving_force
                      - set(self.speciation_keys.keys()))
        if orphan_mdf:
            msg = (f"Gas transfer: molecular_driving_force species "
                   f"{orphan_mdf} have no speciation_keys entry. "
                   f"f_mol cannot be looked up in PropertyResult.alphas — "
                   f"these species will use f_mol=1.0 (no correction).")
            issues.append(msg)
            warnings.warn(msg, UserWarning, stacklevel=2)

        self._validated = True
        return issues

    # ── Simulation lifecycle (C6: unified RunContext) ─────────────────
    # The legacy _enter_simulation / _exit_simulation / _warn_if_simulating
    # methods and the _simulation_running bool were deleted at C6. The
    # gate is now consulted via raise_if_running(self, ...) on every
    # locked mutator below, reading self._context. Set by
    # Simulation.__init__ when this link's owning CV joins a Simulation.

    # ── Setter methods (recommended API) ──────────────────────────────

    def set_transfer_mode(
        self,
        species: str,
        *,
        mode: str,
        kLa: float = None,
        molecular_driving_force: bool = False,
    ) -> None:
        """Set the transfer mode for a species with immediate validation.

        This is the recommended way to configure transfer modes.  It
        handles all bookkeeping atomically and validates immediately.

        Parameters
        ----------
        species : str
            Species identifier (must match a key in ``self.henry``).
        mode : str
            ``"kinetic"`` or ``"equilibrium"``.
        kLa : float, optional
            Volumetric mass transfer coefficient (1/h).  Required for
            kinetic mode; ignored for equilibrium mode.
        molecular_driving_force : bool
            If True (kinetic mode only), the driving force is computed
            from the molecular concentration (BSM2-style) rather than
            total dissolved concentration.  Requires a speciation
            correction to be registered for this species.

        Raises
        ------
        ValueError
            If ``mode`` is not ``"kinetic"`` or ``"equilibrium"``.
            If ``kLa`` is not provided for kinetic mode.
        """
        import warnings

        from .lifecycle import raise_if_running
        raise_if_running(self, f"set_transfer_mode({species!r})")

        mode = mode.lower().strip()
        if mode not in ("kinetic", "equilibrium"):
            raise ValueError(
                f"mode must be 'kinetic' or 'equilibrium', got {mode!r}"
            )

        # Warn if species has no partition model
        if species not in self.partition_models:
            warnings.warn(
                f"set_transfer_mode('{species}'): species has no partition "
                f"model in link.partition_models. Available species: "
                f"{sorted(self.partition_models.keys())}. Transfer will have "
                f"no effect.",
                UserWarning, stacklevel=2,
            )

        if mode == "kinetic":
            if kLa is None:
                raise ValueError(
                    f"kLa is required for kinetic mode (species='{species}')"
                )
            # Remove from equilibrium, add kLa
            self.equilibrium_species.discard(species)
            self.kLa[species] = float(kLa)

            # Handle molecular driving force
            if molecular_driving_force:
                self.molecular_driving_force.add(species)
                if species not in self.speciation_keys:
                    warnings.warn(
                        f"set_transfer_mode('{species}'): "
                        f"molecular_driving_force=True but no "
                        f"speciation_keys entry for this species. f_mol "
                        f"will default to 1.0 (no correction). Add "
                        f"speciation_keys[{species!r}] = '<mol_form_key>' "
                        f"(e.g. 'CO2aq') so the link can look up the "
                        f"alpha in PropertyResult.alphas.",
                        UserWarning, stacklevel=2,
                    )
            else:
                self.molecular_driving_force.discard(species)

        elif mode == "equilibrium":
            # Add to equilibrium, clean up kinetic settings
            self.equilibrium_species.add(species)
            self.kLa.pop(species, None)
            self.molecular_driving_force.discard(species)

        # Invalidate cached validation
        self._validated = False

    def set_henry(self, species: str, partition_model: PartitionModel) -> None:
        """Set the partition model for a species.

        Parameters
        ----------
        species : str
            Species identifier.
        partition_model : PartitionModel
            Partition model (e.g. a :class:`~PyOMES.chemistry.HenryEquilibrium`).
        """
        from .lifecycle import raise_if_running
        raise_if_running(self, f"set_henry({species!r})")
        self.partition_models[species] = partition_model
        self._validated = False

    # ── PhaseInterface entry point (Phase 3) ──────────────────────────

    def compute_flux(
        self,
        state_a: Any,
        state_b: Any,
        dt_h: float,
        *,
        instantaneous: bool = False,
    ) -> Dict[str, float]:
        """PhaseInterface adapter for ``compute_flow``.

        Lets a single :class:`ControlVolume` use this link as an internal
        interface (see CV_UPDATE.md §"Phase 3"). Reads alpha values
        inline from ``state_b.n_mol`` (the liquid phase's canonical
        species), populated by the speciation engine's writeback
        during the CV's ``advance()`` step. No
        :class:`PropertyResult` threading — state-unification C4
        removed that channel.

        Parameters
        ----------
        state_a : GasPhase
            Gas phase state.
        state_b : LiquidPhase
            Liquid phase state.
        dt_h : float
            Timestep duration (hours).
        instantaneous : bool, keyword-only
            Mirrors :meth:`compute_flow`'s ``instantaneous`` flag.  When
            True, returns the instantaneous transfer rate
            ``kLa × V × (C* − C)`` rather than the analytical
            exponential step solution.  Used by ODE-based solvers that
            manage their own internal time-stepping (e.g.
            :class:`~PyOMES.core.solvers.SimultaneousAdaptiveSolver`).

        Returns
        -------
        dict
            ``{species_id: flux_mol_per_h}`` — positive = gas → liquid.
        """
        from types import SimpleNamespace

        # Build a single phases dict that holds both phases.  When
        # gas_cv_key == liquid_cv_key (single-CV case, the canonical
        # PhaseInterface usage), the shared SimpleNamespace keeps both
        # phases visible to compute_flow's two separate lookups; when
        # the keys differ (multi-CV case), the same object answers
        # both cv-key lookups, again with both phases visible.
        shared = SimpleNamespace(phases={
            self.gas_phase_key: state_a,
            self.liquid_phase_key: state_b,
        })
        cvs = {
            self.gas_cv_key: shared,
            self.liquid_cv_key: shared,
        }
        return self.compute_flow(
            cvs, dt_h,
            instantaneous=instantaneous,
        )

    # ── Main computation ──────────────────────────────────────────────

    def compute_flow(
        self,
        cvs: Dict[str, Any],
        dt_h: float,
        instantaneous: bool = False,
    ) -> Dict[str, float]:
        """Compute gas-liquid transfer fluxes for all species.

        Alpha values (molecular fraction of each acid-base ladder)
        are computed inline from the liquid phase's ``n_mol`` —
        canonical species populated by the speciation engine's
        writeback during the CV's ``advance()`` step. The legacy
        ``speciation_override`` PropertyResult channel was removed
        in state-unification C4.

        Parameters
        ----------
        cvs : dict
            ``{cv_key: ControlVolume}`` — the system's CVs.
        dt_h : float
            Timestep duration (hours).  Used by the analytical step
            solution (default mode).  Ignored when ``instantaneous=True``.
        instantaneous : bool
            If True, return the **instantaneous transfer rate** at the
            current state — the simple ``kLa × V × (C* − C)`` driving
            force without the analytical exponential step solution.

            Use ``instantaneous=True`` when calling from within an ODE
            derivative function (e.g. ``SimultaneousAdaptiveSolver``), where the
            integrator manages time-stepping itself.

            Use ``instantaneous=False`` (default) when calling from the
            Euler snapshot solver, where the analytical step solution
            gives a better approximation over the finite timestep.

        Returns
        -------
        dict
            ``{species_id: flux_mol_per_h}`` — positive = gas→liquid.
        """
        gas_cv = cvs.get(self.gas_cv_key)
        liq_cv = cvs.get(self.liquid_cv_key)
        if gas_cv is None or liq_cv is None:
            return {}

        # Run configuration validation once on first call
        if not self._validated:
            self.validate()

        gas_phase = gas_cv.phases.get(self.gas_phase_key)
        liq_phase = liq_cv.phases.get(self.liquid_phase_key)
        if gas_phase is None or liq_phase is None:
            return {}

        V_liq = float(liq_phase.V_L)
        V_gas = float(gas_phase.V_L)
        T_K = float(gas_phase.T_K)
        dt_h = max(float(dt_h), 1e-30)

        flow: Dict[str, float] = {}

        for species, model in self.partition_models.items():
            n_gas_i = float(gas_phase.n_mol.get(species, 0.0))
            n_liq_i = self._liquid_total(species, liq_phase)

            if instantaneous:
                flux = self._instantaneous_rate(
                    species, model, n_gas_i, n_liq_i, V_gas, V_liq, T_K,
                    liq_phase,
                )
            elif species in self.equilibrium_species:
                flux = self._equilibrium_flux(
                    species, model, n_gas_i, n_liq_i, V_gas, V_liq, T_K,
                    dt_h, liq_phase,
                )
            else:
                flux = self._kinetic_flux(
                    species, model, n_gas_i, n_liq_i, V_gas, V_liq, T_K,
                    dt_h, liq_phase,
                )

            # Safety clamp: don't transfer more than available
            # (skip for instantaneous — ODE solver handles limits)
            if not instantaneous:
                flux = self._clamp_flux(flux, n_gas_i, n_liq_i, dt_h)

            if abs(flux) > 0.0:
                flow[species] = flux

        return flow

    # ── Kinetic mode ──────────────────────────────────────────────────

    def _kinetic_flux(
        self,
        species: str,
        model: PartitionModel,
        n_gas: float,
        n_liq: float,
        V_gas: float,
        V_liq: float,
        T_K: float,
        dt_h: float,
        liq_phase: Any,
    ) -> float:
        """Compute kLa-driven flux using the coupled analytical solution.

        Models gas and liquid as coupled reservoirs with first-order
        kinetic transfer.  The analytical solution is an exponential
        approach to the equilibrium partition::

            n_liq(dt) = n_liq_eq + (n_liq_0 − n_liq_eq) × exp(−B × dt)

        Two driving force formulations are supported:

        **Default (total-IC basis):**

        Uses model.partition_ratio(V_liq, V_gas, T_K, alpha) which incorporates the
        speciation correction (kH_eff = kH / alpha) internally::

            β = model.partition_ratio(V_liq, V_gas, T_K, alpha)
            n_liq_eq = β × n_total / (1 + β)
            B = kLa × (1 + β)

        **Molecular driving force** (``species in self.molecular_driving_force``):

        Uses model.partition_ratio with alpha=1 to get r_raw (intrinsic Henry
        constant, no amplification), and f_mol = alpha to represent the
        molecular fraction.  The BSM2 formulation::

            β_raw = model.partition_ratio(V_liq, V_gas, T_K, alpha=1)
            n_liq_eq = β_raw × n_total / (β_raw + f_mol)
            B = kLa × (β_raw + f_mol)

        If model.partition_ratio() returns None (non-linear model), falls back to
        _instantaneous_rate().

        Positive = gas→liquid (absorption).
        """
        kla = float(self.kLa.get(species, 0.0))
        if kla <= 0.0:
            return 0.0

        n_total = n_gas + n_liq
        if n_total <= 0.0:
            return 0.0

        use_molecular = species in self.molecular_driving_force
        alpha = self._alpha_for(species, liq_phase)

        if use_molecular:
            f_mol = max(1e-12, float(alpha)) if alpha is not None else 1.0
            beta_raw = model.partition_ratio(V_liq, V_gas, T_K, alpha=1.0)
            if beta_raw is None:
                return self._instantaneous_rate(
                    species, model, n_gas, n_liq, V_gas, V_liq, T_K, liq_phase
                )
            n_liq_eq = beta_raw * n_total / (beta_raw + f_mol)
            B = kla * (beta_raw + f_mol)
        else:
            alpha_val = float(alpha) if alpha is not None else 1.0
            beta = model.partition_ratio(V_liq, V_gas, T_K, alpha=alpha_val)
            if beta is None:
                return self._instantaneous_rate(
                    species, model, n_gas, n_liq, V_gas, V_liq, T_K, liq_phase
                )
            n_liq_eq = beta * n_total / (1.0 + beta)
            B = kla * (1.0 + beta)

        delta = n_liq_eq - n_liq
        if abs(delta) < 1e-30:
            return 0.0

        exponent = B * dt_h
        if exponent > 50.0:
            frac = 1.0
        else:
            frac = 1.0 - math.exp(-exponent)

        return delta * frac / dt_h

    # ── Instantaneous rate (for ODE derivatives) ──────────────────

    _KLA_EQUILIBRIUM: float = 1e6  # effective kLa for equilibrium species in ODE mode

    def _instantaneous_rate(
        self,
        species: str,
        model: PartitionModel,
        n_gas: float,
        n_liq: float,
        V_gas: float,
        V_liq: float,
        T_K: float,
        liq_phase: Any,
    ) -> float:
        """Instantaneous transfer rate for use in ODE derivative functions.

        Returns the first-order driving force::

            flux = kLa × V_liq × (C_star − C)

        where ``C_star = model.partition_ratio(...) × n_gas / V_liq`` is the Henry
        equilibrium saturation at the current gas-phase partial pressure
        and ``C = n_liq / V_liq`` is the dissolved concentration.

        For equilibrium-mode species, a very high effective kLa
        (``_KLA_EQUILIBRIUM = 1e6 /h``) is used so the ODE integrator
        drives the system to near-equilibrium.

        Positive = gas→liquid (absorption).
        """
        if species in self.equilibrium_species:
            kla = self._KLA_EQUILIBRIUM
        else:
            kla = float(self.kLa.get(species, 0.0))
        if kla <= 0.0:
            return 0.0

        use_molecular = species in self.molecular_driving_force
        alpha = self._alpha_for(species, liq_phase)

        if use_molecular:
            f_mol = max(1e-12, float(alpha)) if alpha is not None else 1.0
            beta_raw = model.partition_ratio(V_liq, V_gas, T_K, alpha=1.0)
            if beta_raw is not None:
                C_star = beta_raw * n_gas / V_liq if V_liq > 0 else 0.0
            else:
                n_total = n_gas + n_liq
                n_liq_eq = model.equilibrium_a_moles(n_total, V_liq, V_gas, T_K, alpha=1.0)
                C_star = n_liq_eq / V_liq if V_liq > 0 else 0.0
                f_mol = 1.0
            C_mol = f_mol * n_liq / V_liq if V_liq > 0 else 0.0
            return kla * V_liq * (C_star - C_mol)
        else:
            alpha_val = float(alpha) if alpha is not None else 1.0
            beta = model.partition_ratio(V_liq, V_gas, T_K, alpha=alpha_val)
            if beta is not None:
                C_star = beta * n_gas / V_liq if V_liq > 0 else 0.0
            else:
                n_total = n_gas + n_liq
                n_liq_eq = model.equilibrium_a_moles(n_total, V_liq, V_gas, T_K, alpha=alpha_val)
                C_star = n_liq_eq / V_liq if V_liq > 0 else 0.0
            C = n_liq / V_liq if V_liq > 0 else 0.0
            return kla * V_liq * (C_star - C)

    # ── Equilibrium mode ──────────────────────────────────────────────

    def _equilibrium_flux(
        self,
        species: str,
        model: PartitionModel,
        n_gas: float,
        n_liq: float,
        V_gas: float,
        V_liq: float,
        T_K: float,
        dt_h: float,
        liq_phase: Any,
    ) -> float:
        """Compute flux to reach Henry equilibrium in one timestep.

        flux = (n_liq_eq − n_liq_current) / dt_h

        where n_liq_eq = model.equilibrium_a_moles(...) incorporating
        the alpha correction for speciating species.
        """
        n_total = max(0.0, n_gas + n_liq)
        if n_total <= 0.0:
            return 0.0

        alpha = self._alpha_for(species, liq_phase)
        alpha_val = float(alpha) if alpha is not None else 1.0
        n_liq_eq = model.equilibrium_a_moles(n_total, V_liq, V_gas, T_K, alpha=alpha_val)
        return (n_liq_eq - n_liq) / dt_h

    def _liquid_total(self, species: str, liq_phase: Any) -> float:
        """Sum n_mol over the equilibrium ladder for a gas species.

        For gas species that participate in an acid-base ladder
        (e.g. ``CO2``), the liquid-side TOTAL inventory is the sum of
        all related molecular species in n_mol (e.g. ``CO2aq`` +
        ``HCO3-`` + ``CO3--``). The ladder is looked up by gas-side
        key via :data:`~PyOMES.chemical_equilibrium.acid_base._CANONICAL_NAMES`;
        species without a recognised ladder fall back to a single
        n_mol entry under the same key (e.g. ``O2``, ``N2``).

        Introduced in state-unification C3 to replace the previous
        single-key read (``liq_phase.n_mol[species]``), which relied
        on the BSM2-conflated convention where the gas-side key
        also stored the conserved total. After C3, that convention
        is gone — totals are implicit (summed) and the link sums
        explicitly here.
        """
        n_mol = getattr(liq_phase, "n_mol", {})
        ladder = self.speciation_ladders.get(species)
        if ladder is None:
            return float(n_mol.get(species, 0.0))
        ladder_sum = float(sum(n_mol.get(sp, 0.0) for sp in ladder))
        if ladder_sum > 0.0:
            return ladder_sum
        # Pre-speciation state (ladder species not yet populated by an
        # engine solve) — fall back to the bare-name total.
        return float(n_mol.get(species, 0.0))

    def _alpha_for(self, species: str, liq_phase: Any) -> Optional[float]:
        """Compute the molecular fraction for a gas species inline
        from the liquid phase's ``n_mol``.

        alpha = ``n_mol[mol_key] / sum(n_mol[s] for s in ladder)``,
        where ``mol_key = speciation_keys[species]`` (the molecular
        form's species id) and ``ladder`` is the canonical-name set
        for the gas species, looked up via
        :data:`~PyOMES.chemical_equilibrium.acid_base._CANONICAL_NAMES`.

        Returns ``None`` when:

        - The species is not registered in ``speciation_keys``
          (callers treat ``None`` as "no correction" — single-species
          transfers like O2, N2 hit this branch).
        - The species has no canonical ladder in
          :data:`_CANONICAL_NAMES`.
        - The ladder species are not yet populated (pre-speciation
          state, no engine has run yet on this liquid phase).
        """
        mol_key = self.speciation_keys.get(species)
        if mol_key is None:
            return None
        ladder = self.speciation_ladders.get(species)
        if ladder is None:
            return None
        n_mol = getattr(liq_phase, "n_mol", {})
        total = sum(float(n_mol.get(sp, 0.0)) for sp in ladder)
        if total <= 0.0:
            return None
        mol_n = float(n_mol.get(mol_key, 0.0))
        return min(1.0, max(0.0, mol_n / total))

    # ── Safety clamping ───────────────────────────────────────────────

    def _clamp_flux(
        self,
        flux: float,
        n_gas: float,
        n_liq: float,
        dt_h: float,
    ) -> float:
        """Clamp flux so it doesn't transfer more than available.

        Positive flux = gas→liquid: limited by gas inventory.
        Negative flux = liquid→gas: limited by liquid inventory.
        """
        if flux > 0.0:
            # Absorption: limited by gas moles
            max_flux = float(n_gas) / dt_h
            return min(flux, max_flux)
        elif flux < 0.0:
            # Stripping: limited by liquid moles
            max_flux = float(n_liq) / dt_h
            return max(flux, -max_flux)
        return 0.0

    # ── Mutators for controller integration ────────────────────────────

    def set_kLa(self, species: str, value: float) -> None:
        """Update the kLa for a specific species.

        Gated by the RunContext lifecycle (C6): raises ``RuntimeError``
        if called mid-``Simulation.run``.  The orchestrator's
        controller-action apply path uses the reflective walker
        (PARAM_PATH_DISPATCHER) to bypass the gate.
        """
        from .lifecycle import raise_if_running
        raise_if_running(self, f"set_kLa({species!r})")
        self.kLa[species] = float(value)

    def set_kLa_with_co2_ratio(
        self, kLa_O2: float, co2_ratio: float = 0.9
    ) -> None:
        """Set kLa for O₂ and CO₂ simultaneously via a fixed ratio.

        Gated by the RunContext (C6).
        """
        from .lifecycle import raise_if_running
        raise_if_running(self, "set_kLa_with_co2_ratio")
        self.kLa["O2"] = float(kLa_O2)
        self.kLa["CO2"] = float(kLa_O2) * float(co2_ratio)

    # ── Representation ────────────────────────────────────────────────

    def __repr__(self):
        species = sorted(self.partition_models.keys())
        modes = []
        for sp in species:
            m = "eq" if sp in self.equilibrium_species else f"kLa={self.kLa.get(sp, 0):.1f}"
            modes.append(f"{sp}:{m}")
        mode_str = ", ".join(modes)
        return (f"KineticGasLiquidLink("
                f"{self.gas_cv_key}.{self.gas_phase_key} -> "
                f"{self.liquid_cv_key}.{self.liquid_phase_key}, "
                f"[{mode_str}])")
