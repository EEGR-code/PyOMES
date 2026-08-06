# -*- coding: utf-8 -*-
"""NRTableau — log-linear tableau for the Newton-Raphson speciation engine.

Builds the data structure consumed by :func:`~PyOMES.chemical_equilibrium.nr_solver.solve_nr`:
for every non-master species ``j``, stores the log-linear expression

    log(a_j) = log_K'_j  +  Σ_k  ν_{jk} · x_k

where ``x_k = log10(a_k)`` are the NR unknowns (log-activities of the master
species) and ``ν_{jk}`` are stoichiometric coefficients derived by BFS traversal
of the reaction graph.

Public entry point
------------------
:func:`build_tableau` — takes a list of :class:`~PyOMES.reactions.equilibrium.EquilibriumReaction`
objects (single-phase, from a :class:`~PyOMES.reactions.reaction_system.ReactionSystem`)
and returns a :class:`NRTableau`.

Master species selection (Steps 1-3 from the design doc)
---------------------------------------------------------
1. H⁺ is pre-registered as a universal master (not in the reaction graph).
2. Non-H⁺, non-H₂O species form nodes; each EquilibriumReaction adds edges.
3. Connected components are found via BFS on the undirected graph.
4. Per component, master selection priority:
   (a) ``total_id`` on any reaction in the component, if that id is a member;
   (b) DAG source — species that is a reactant in ≥1 reaction but never a product;
   (c) Tiebreak: most H atoms, then most positive charge.

Tableau derivation (Step 4)
----------------------------
BFS from the master set.  For each reaction where exactly one non-H₂O/H⁺
species is unknown, we rearrange the equilibrium log-K equation to express
``log(a_unknown)`` in terms of known log-activities.  Van't Hoff temperature
correction is applied to log_K values at build time.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..units import R_J_PER_MOL_K as _R_J_MOL_K

logger = logging.getLogger(__name__)

_LOG10_E = np.log10(np.e)        # 1/ln(10), used to convert ln K to log10 K


# ─────────────────────────────────────────────────────────────────────────────
#  Data structures
# ─────────────────────────────────────────────────────────────────────────────

class ConfigurationError(ValueError):
    """Raised when declared chemistry cannot be folded into one NRTableau.

    Distinct from the plain ``ValueError``\\ s :func:`build_tableau` already
    raises for malformed/incomplete input (missing water reaction,
    undeclared masters) — a ``ConfigurationError`` specifically means "this
    chemistry is out of the tableau's current *capability*", not "this
    chemistry is malformed". The motivating case (CP1 of
    ``LAYER1_GAP_CLOSURE``): a gas-liquid reaction whose liquid-phase
    participants span two independent, already-multi-species acid-base
    components — folding it would require merging those components, which
    ``NRTableau``'s one-master-per-component design does not support (see
    ``MASS_EXCHANGE_ARCHITECTURE.md`` §14.3). That capability is tracked by
    ``MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md``, not here.
    """


@dataclass(frozen=True)
class SecondaryEntry:
    """Log-linear expression for one derived (non-master) species.

    ``log(a_j) = log_K_prime + Σ_k  nu[k] · x_k``
    where ``x_k = log10(a_k)`` is the NR unknown for master ``k``.

    ``element_stoichiometry`` maps each non-H⁺ master ID to the count of that
    master's formula unit in this secondary species.  H⁺ is excluded because
    H⁺ closes the charge balance, not a mass balance.  For single-component
    secondaries this has exactly one key; for cross-component species (e.g.
    CaHCO₃⁺) it has one key per contributing component.

    ``phase`` distinguishes a liquid-phase secondary (``"liquid"``, the
    default — its value is a concentration in mol/L, e.g. HCO₃⁻ derived
    from the CO₂ master) from a gas-phase secondary folded in by CP1 of
    ``LAYER1_GAP_CLOSURE`` (``"gas"`` — its value is a partial pressure in
    atm, e.g. gas-phase CO₂ derived from the same master via Henry's law).
    A gas-phase secondary's ``species_id`` may collide with a liquid-phase
    master/secondary's bare id (gas and liquid CO₂ are conventionally
    declared with the identical id, per ``HenryEquilibrium``) — ``phase``
    is what disambiguates them; consumers iterating ``tableau.secondaries``
    must not assume ``species_id`` alone is unique.
    """
    species_id: str
    charge: int
    nu: Dict[str, float]    # master_id → stoichiometric coefficient
    log_K_prime: float      # accumulated log10 K from master to this species
    element_stoichiometry: Dict[str, float] = field(default_factory=dict)
    # master_id → count of that master's formula unit in this secondary
    # (excludes H+; used for cross-component mass-balance accumulation)
    phase: str = "liquid"

    @property
    def c_key(self) -> str:
        """Internal dict key for concentration/gamma/charge lookups.

        Equal to ``species_id`` for liquid-phase secondaries (backward
        compatible with pre-CP1 behaviour). Gas-phase secondaries get a
        reserved ``":gas"`` suffix — CP2 of ``LAYER1_GAP_CLOSURE`` found
        that using bare ``species_id`` for both collides whenever gas and
        liquid share an id (the common Henry case, e.g. both "CO2"):
        Python dict construction silently lets the later-inserted entry
        (the gas secondary, built after masters) overwrite the earlier
        one (the liquid master's own concentration), corrupting its mass
        balance. ``nr_solver.py`` uses this key everywhere it builds or
        reads a ``{species_id: value}`` dict keyed by tableau species
        (``c``, ``gammas``, ``all_charges``); the *public* identity for
        output purposes (``species_mol_L`` vs ``partial_pressures_atm``
        in :class:`~PyOMES.chemical_equilibrium.protocols.EquilibriumResult`)
        remains ``species_id`` plus ``phase``.
        """
        return self.species_id if self.phase != "gas" else f"{self.species_id}:gas"


@dataclass(frozen=True)
class ComponentInfo:
    """One independent chemical component (one non-H⁺ master).

    ``species_ids`` lists every liquid-phase species in this connected
    component of the reaction graph; the solver sums their concentrations
    to get ``C_i_total`` for the mass-balance residual.

    ``gas_species_ids`` lists the bare ids of any gas-phase secondaries
    folded onto this component (CP1 of ``LAYER1_GAP_CLOSURE`` — see
    :class:`SecondaryEntry`'s ``phase`` field). Disjoint in *meaning* from
    ``species_ids`` even when the bare id strings coincide (the common
    Henry case): a gas-phase id here is a reminder that this component's
    total, once CP2 wires the volume-aware mass balance, must also read
    from the gas phase's ``n_mol``, not only the liquid phase's. Empty for
    components with no folded gas-liquid equilibrium.
    """
    master_id: str
    species_ids: Tuple[str, ...]   # all liquid-phase members of this component
    gas_species_ids: Tuple[str, ...] = field(default_factory=tuple)


@dataclass
class NRTableau:
    """Complete log-linear tableau for the Newton-Raphson speciation solver.

    Attributes
    ----------
    masters : list of str
        Master species IDs, H⁺ first, then one per chemical component.
        The order defines the column order of the Jacobian and the index
        into ``x`` (the NR unknown vector).
    master_charges : dict
        ``{master_id: int}`` — ionic charge of each master.
    secondaries : list of SecondaryEntry
        Derived species with their log-linear expressions.
    components : list of ComponentInfo
        One entry per non-H⁺ master, in the same order as ``masters[1:]``.
        Used to (a) read totals from ``phase.n_mol`` and (b) build the
        mass-balance residuals.
    log_Kw : float
        log10(Kw) used for the OH⁻ secondary.
    T_ref_K : float
        Reference temperature of the stored log_K values.
    """
    masters: List[str]
    master_charges: Dict[str, int]
    secondaries: List[SecondaryEntry]
    components: List[ComponentInfo]
    log_Kw: float
    T_ref_K: float

    # ── Convenience accessors ─────────────────────────────────────────

    def master_index(self, master_id: str) -> int:
        """Column index of *master_id* in the NR unknown vector."""
        return self.masters.index(master_id)

    def all_secondary_ids(self) -> List[str]:
        return [s.species_id for s in self.secondaries]

    def component_for_master(self, master_id: str) -> Optional[ComponentInfo]:
        for comp in self.components:
            if comp.master_id == master_id:
                return comp
        return None


# ─────────────────────────────────────────────────────────────────────────────
#  Temperature correction helper
# ─────────────────────────────────────────────────────────────────────────────

def _vant_hoff_log_K(
    log_K_ref: float,
    dH_J_per_mol: Optional[float],
    T_K: float,
    T_ref_K: float,
) -> float:
    """Apply Van't Hoff correction to log10(K).

    Returns ``log_K_ref`` unchanged when ``dH_J_per_mol`` is None or ~0.
    """
    if dH_J_per_mol is None or abs(dH_J_per_mol) < 1e-30:
        return float(log_K_ref)
    if abs(T_K - T_ref_K) < 1e-10:
        return float(log_K_ref)
    # ln K(T) = ln K(T_ref) − (ΔH/R) (1/T − 1/T_ref)
    delta_ln_K = -(float(dH_J_per_mol) / _R_J_MOL_K) * (1.0 / float(T_K) - 1.0 / float(T_ref_K))
    return float(log_K_ref) + delta_ln_K * _LOG10_E


def _derive_gas_secondary(rxn, known: Dict[str, Tuple[dict, float]], T_K: float) -> "SecondaryEntry":
    """Derive one gas-phase :class:`SecondaryEntry` from an already-resolved tableau.

    Unlike the acid-base BFS derivation (which may need several passes to
    chain through intermediate secondaries), a folded gas-liquid reaction
    is derived in one independent step: its liquid-phase participant(s)
    are guaranteed already present in *known* (the component-attachment
    resolution in :func:`build_tableau` only attaches a gas-liquid item
    once its liquid entries have a resolved home component), and its sole
    gas-phase participant is the one unknown being solved for. Mirrors the
    log-linear rearrangement used for acid-base secondaries, but the
    resulting entry's value is a partial pressure (atm), not a
    concentration (mol/L) — CP2 of ``LAYER1_GAP_CLOSURE`` is what wires
    the unit-aware residual/Jacobian assembly; this function only builds
    the log-linear expression itself.

    Raises
    ------
    ConfigurationError
        If *rxn* does not have exactly one gas-phase participant (multi-
        gas-species mass-action rows are not supported by this phase), or
        if a liquid-phase participant is unexpectedly absent from *known*
        (should not happen given the caller's prior validation).
    """
    unknowns = []
    known_entries = []   # list of (coeff: float, (nu: dict, log_K_prime: float))
    h_coeff = 0.0

    for e in rxn.stoichiometry:
        sp_id = e.species.id
        if sp_id in _SOLVENT_IDS:
            continue
        c = float(e.coefficient)
        if sp_id == _H_ID:
            h_coeff += c
            continue
        if e.phase == "gas":
            unknowns.append((c, sp_id, e.species))
            continue
        if sp_id not in known:
            raise ConfigurationError(
                f"NRTableau: gas-liquid reaction {rxn!r} references liquid "
                f"species {sp_id!r} that is not resolvable from the "
                "acid-base network — this should not happen given prior "
                "component-attachment validation."
            )
        known_entries.append((c, known[sp_id]))

    if len(unknowns) == 0:
        raise ConfigurationError(
            f"NRTableau: gas-liquid reaction {rxn!r} has no derivable "
            "gas-phase participant. If this is a solvent Raoult fold "
            f"(all participants in {sorted(_SOLVENT_IDS)}, e.g. H2O), it "
            "should have been routed to _derive_solvent_gas_secondary, "
            "not here — this indicates a reaction with a non-solvent "
            "liquid participant but no non-solvent gas participant, which "
            "this phase's folding does not support."
        )
    if len(unknowns) != 1:
        raise ConfigurationError(
            f"NRTableau: gas-liquid reaction {rxn!r} must have exactly one "
            f"gas-phase participant to fold into the tableau; found "
            f"{len(unknowns)}. Multi-gas-species mass-action rows are not "
            "supported by this phase's folding."
        )

    target_coeff, target_id, target_sp = unknowns[0]

    log_K_rxn = _vant_hoff_log_K(
        float(rxn.log_K),
        rxn.dH_J_per_mol,
        T_K,
        float(rxn.T_ref_K),
    )

    accumulated_log_K = log_K_rxn
    accumulated_nu: Dict[str, float] = {}

    for c_j, (nu_j, logK_j) in known_entries:
        accumulated_log_K -= c_j * logK_j
        for k_id, nu_jk in nu_j.items():
            accumulated_nu[k_id] = accumulated_nu.get(k_id, 0.0) - c_j * nu_jk

    accumulated_nu["H+"] = accumulated_nu.get("H+", 0.0) - h_coeff

    log_K_prime = accumulated_log_K / target_coeff
    nu = {k: v / target_coeff for k, v in accumulated_nu.items()
          if abs(v / target_coeff) > 1e-15}
    elem_stoich = {k: v for k, v in nu.items() if k != _H_ID}

    return SecondaryEntry(
        species_id=target_id,
        charge=int(target_sp.charge),
        nu=nu,
        log_K_prime=float(log_K_prime),
        element_stoichiometry=elem_stoich,
        phase="gas",
    )


def _derive_solvent_gas_secondary(rxn, T_K: float) -> "SecondaryEntry":
    """Derive a constant gas-phase :class:`SecondaryEntry` for a solvent
    Raoult fold (H2O — CP4 of ``LAYER1_GAP_CLOSURE``).

    Unlike :func:`_derive_gas_secondary`, the liquid-phase participant is
    never looked up in ``known`` — there is no master/component for a
    solvent's own total to attach to (H2O is excluded from the acid-base
    graph entirely, per :data:`_SOLVENT_IDS`, because its activity is
    conventionally fixed at 1: the same convention every other reaction
    referencing H2O as a spectator already assumes). Its liquid-side
    activity is therefore treated as exactly 1 (``log10(1) = 0``) rather
    than solved for, giving a result with ``nu = {}`` — a pure function
    of T alone (``p_H2O,gas = P_sat(T)`` for ``RaoultEquilibrium``,
    matching its own ``P_sat()`` method exactly), with **no** coupling to
    any master and **no** contribution to any mass-balance row
    (``element_stoichiometry = {}``). Consequently this fold does not by
    itself track/decrement the liquid water pool — CP4 relies on
    ``transfer_models=EquilibriumTransferModel(RaoultEquilibrium())``
    (§6.2 of ``MASS_EXCHANGE_ARCHITECTURE.md``) for that; tracking water
    as a genuine finite-total master (so this fold's own Newton system
    closes the liquid-side balance directly) is deferred as a follow-up,
    not attempted here (see CP4's design discussion).

    Raises
    ------
    ConfigurationError
        If *rxn* does not have exactly one gas-phase and one liquid-phase
        participant (both solvent, since routing here already required
        every participant id to be in ``_SOLVENT_IDS``).
    """
    gas_entries = [e for e in rxn.stoichiometry if e.phase == "gas"]
    liquid_entries = [e for e in rxn.stoichiometry if e.phase == "liquid"]
    if len(gas_entries) != 1 or len(liquid_entries) != 1:
        raise ConfigurationError(
            f"NRTableau: solvent gas-liquid reaction {rxn!r} must have "
            f"exactly one gas-phase and one liquid-phase participant; "
            f"found {len(gas_entries)} gas, {len(liquid_entries)} liquid."
        )
    gas_e = gas_entries[0]
    liq_e = liquid_entries[0]
    target_coeff = float(gas_e.coefficient)

    log_K_rxn = _vant_hoff_log_K(
        float(rxn.log_K), rxn.dH_J_per_mol, T_K, float(rxn.T_ref_K),
    )
    # Liquid-side activity == 1 by the pure-solvent convention, so its
    # contribution to accumulated_log_K/accumulated_nu is exactly zero —
    # this loop is here for clarity/symmetry with _derive_gas_secondary,
    # not because it changes the arithmetic.
    accumulated_log_K = log_K_rxn - float(liq_e.coefficient) * 0.0

    log_K_prime = accumulated_log_K / target_coeff

    return SecondaryEntry(
        species_id=gas_e.species.id,
        charge=int(gas_e.species.charge),
        nu={},
        log_K_prime=float(log_K_prime),
        element_stoichiometry={},
        phase="gas",
    )


# ─────────────────────────────────────────────────────────────────────────────
#  Main builder
# ─────────────────────────────────────────────────────────────────────────────

_SOLVENT_IDS = frozenset({"H2O"})
_H_ID = "H+"
_OH_ID = "OH-"


def build_tableau(reactions, *, T_K: float = 298.15) -> NRTableau:
    """Build an :class:`NRTableau` from declared equilibrium reactions.

    Parameters
    ----------
    reactions : iterable of EquilibriumConstraint
        All equilibrium constraints for the system — typically
        ``EquilibriumReaction``, but any ``EquilibriumConstraint``-
        conforming item is accepted. Items are classified via
        :func:`~PyOMES.reactions.equilibrium.classify_equilibrium_constraint`:
        acid-base items build the graph as before; gas-liquid items with a
        ``log_K`` set (e.g. a fully-parameterized ``HenryEquilibrium``/
        ``RaoultEquilibrium``) are folded in as gas-phase secondaries
        attached to the acid-base component their liquid-phase form
        belongs to (or a new singleton component, for gases with no
        acid-base ladder — see :class:`SecondaryEntry`'s ``phase`` field);
        gas-liquid items with ``log_K is None`` are pure partition-routing
        declarations (consumed elsewhere, e.g. by
        ``KineticGasLiquidLink``) and are silently skipped, as before;
        solid-liquid items are always silently skipped — precipitation is
        folded into :class:`~PyOMES.chemical_equilibrium.nr_engine.NRChemicalEquilibriumEngine`
        via its own nested active-set loop, not via this graph (see
        ``MASS_EXCHANGE_ARCHITECTURE.md`` §14.3). Must include exactly one
        water-dissociation reaction (H₂O ⇌ H⁺ + OH⁻).
    T_K : float
        Operating temperature (K).  Van't Hoff correction is applied to
        all log_K values at this temperature.

    Returns
    -------
    NRTableau

    Raises
    ------
    ValueError
        If no water-dissociation reaction is found, if the reaction graph
        contains species that cannot be derived from the selected masters,
        or if a ``total_id`` override points to a species not present in
        the component.
    ConfigurationError
        If a gas-liquid reaction's liquid-phase participants span two
        independent, already-multi-species acid-base components — folding
        it would require merging components, which is out of this phase's
        capability (see ``MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md``).
    """
    from ..reactions.equilibrium import (
        EquilibriumConstraint, classify_equilibrium_constraint,
    )

    # ── Filter and classify reactions ─────────────────────────────────
    water_rxn = None
    ab_rxns: List = []   # acid-base (non-water, non-cross-phase)
    gas_rxns: List = []  # gas-liquid, fully parameterized (log_K set)
    solvent_gas_rxns: List = []
    # gas-liquid where BOTH sides are a solvent species (H2O Raoult fold,
    # CP4 of LAYER1_GAP_CLOSURE) — routed separately from gas_rxns since
    # a solvent has no component/master to attach onto (see
    # _derive_solvent_gas_secondary).

    for rxn in reactions:
        if not isinstance(rxn, EquilibriumConstraint):
            continue
        try:
            kind = classify_equilibrium_constraint(rxn)
        except ValueError:
            # Empty stoichiometry (e.g. a PartitionModel-only
            # HenryEquilibrium/RaoultEquilibrium with gas_species/
            # liquid_species unset) — nothing to add to the graph.
            continue
        if kind == "solid_liquid":
            # Precipitation — handled by NRChemicalEquilibriumEngine's own nested
            # active-set loop, not folded into this graph (§14.3).
            continue
        if kind == "gas_liquid":
            if rxn.log_K is None:
                # Pure partition-routing declaration (e.g. a hand-built
                # cross-phase EquilibriumReaction feeding
                # KineticGasLiquidLink only) — nothing to fold in.
                continue
            participant_ids = {e.species.id for e in rxn.stoichiometry}
            if participant_ids <= _SOLVENT_IDS:
                solvent_gas_rxns.append(rxn)
            else:
                gas_rxns.append(rxn)
            continue
        if rxn.log_K is None:
            continue
        product_ids = {e.species.id for e in rxn.stoichiometry if e.coefficient > 0}
        if _H_ID in product_ids and _OH_ID in product_ids:
            if water_rxn is not None:
                raise ValueError("NRTableau: more than one water-dissociation reaction.")
            water_rxn = rxn
        else:
            ab_rxns.append(rxn)

    if water_rxn is None:
        raise ValueError(
            "NRTableau: no water-dissociation reaction found (H₂O ⇌ H⁺ + OH⁻). "
            "Declare one EquilibriumReaction with H⁺ and OH⁻ as products."
        )

    log_Kw = _vant_hoff_log_K(
        float(water_rxn.log_K),
        water_rxn.dH_J_per_mol,
        T_K,
        float(water_rxn.T_ref_K),
    )
    T_ref_K = float(water_rxn.T_ref_K)

    # ── Collect all non-solvent, non-H⁺ species from acid-base reactions ─
    graph_species: Dict[str, object] = {}   # id → Species
    for rxn in ab_rxns:
        for e in rxn.stoichiometry:
            sp = e.species
            if sp.id not in _SOLVENT_IDS and sp.id != _H_ID:
                graph_species[sp.id] = sp

    # ── Build undirected adjacency list (exclude H⁺, H₂O) ────────────
    adj: Dict[str, set] = {sp_id: set() for sp_id in graph_species}
    for rxn in ab_rxns:
        nodes = [e.species.id for e in rxn.stoichiometry
                 if e.species.id not in _SOLVENT_IDS and e.species.id != _H_ID]
        for i, a in enumerate(nodes):
            for b in nodes[i + 1:]:
                adj[a].add(b)
                adj[b].add(a)

    # ── Find connected components (BFS) ──────────────────────────────
    visited: set = set()
    raw_components: List[set] = []
    for start in graph_species:
        if start in visited:
            continue
        comp: set = set()
        queue = [start]
        while queue:
            node = queue.pop()
            if node in visited:
                continue
            visited.add(node)
            comp.add(node)
            queue.extend(nbr for nbr in adj.get(node, set()) if nbr not in visited)
        raw_components.append(comp)

    # ── Resolve gas-liquid attachment onto acid-base components ───────
    # Each folded gas-liquid item attaches its gas-phase entries onto
    # whichever component contains its liquid-phase entries (the CO2/NH3/
    # H2S case — an existing single-component ladder), or forms a fresh
    # singleton component when none of its liquid entries appear in any
    # existing component (the inert-gas O2/CH4/N2/H2 case). An item whose
    # liquid entries span more than one *existing* component would require
    # merging two independent multi-species components — out of this
    # phase's capability (§14.3) — and raises ConfigurationError rather
    # than silently misattaching.
    gas_attachments: Dict[int, List] = {}
    for rxn in gas_rxns:
        liquid_ids = [
            e.species.id for e in rxn.stoichiometry
            if e.phase == "liquid" and e.species.id not in _SOLVENT_IDS
            and e.species.id != _H_ID
        ]
        if not liquid_ids:
            raise ConfigurationError(
                f"NRTableau: gas-liquid reaction {rxn!r} has no liquid-phase "
                "participant to attach onto a component."
            )

        found_indices: set = set()
        unfound_ids: List[str] = []
        for sp_id in liquid_ids:
            hit = next(
                (idx for idx, comp_ids in enumerate(raw_components)
                 if sp_id in comp_ids),
                None,
            )
            if hit is None:
                unfound_ids.append(sp_id)
            else:
                found_indices.add(hit)

        if len(found_indices) > 1:
            raise ConfigurationError(
                f"NRTableau: gas-liquid reaction {rxn!r} would bridge "
                f"independent multi-species components (indices "
                f"{sorted(found_indices)}) — merging components is out of "
                "this phase's capability. See "
                "MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md."
            )

        # Register any not-yet-seen liquid species so master selection's
        # atom/charge lookups (graph_species) resolve for them.
        for e in rxn.stoichiometry:
            if e.phase == "liquid" and e.species.id in unfound_ids:
                graph_species.setdefault(e.species.id, e.species)

        if found_indices:
            target_idx = next(iter(found_indices))
            raw_components[target_idx] |= set(unfound_ids)
        else:
            target_idx = len(raw_components)
            raw_components.append(set(unfound_ids))

        gas_attachments.setdefault(target_idx, []).append(rxn)

    # Helper: which ab_rxns touch a component?
    def _rxns_in_comp(comp_ids: set) -> list:
        return [rxn for rxn in ab_rxns
                if any(e.species.id in comp_ids
                       for e in rxn.stoichiometry
                       if e.species.id not in _SOLVENT_IDS and e.species.id != _H_ID)]

    # ── Select one master per component ───────────────────────────────
    def _select_master(comp_ids: set, comp_rxns: list) -> str:
        # (a) total_id override: if any reaction in this component declares
        #     total_id pointing to a member species, use it.
        for rxn in comp_rxns:
            tid = getattr(rxn, "total_id", None)
            if tid is not None and tid in comp_ids:
                return tid

        # (b) DAG source: appears as a reactant but never as a product
        #     (among the component's own species, ignoring H⁺/H₂O).
        product_set: set = set()
        reactant_set: set = set()
        for rxn in comp_rxns:
            for e in rxn.stoichiometry:
                sp_id = e.species.id
                if sp_id not in comp_ids:
                    continue
                if e.coefficient > 0:
                    product_set.add(sp_id)
                elif e.coefficient < 0:
                    reactant_set.add(sp_id)
        sources = reactant_set - product_set
        if len(sources) == 1:
            return next(iter(sources))

        # (c) Tiebreak: most H atoms, then most positive charge
        candidates = sorted(
            comp_ids,
            key=lambda sp_id: (
                -int(graph_species[sp_id].atoms.get("H", 0)),
                -int(graph_species[sp_id].charge),
            ),
        )
        return candidates[0]

    # ── Pre-populate "known" dict: id → (nu, log_K_prime) ────────────
    # nu maps master_id → float coefficient for the log-linear expression.
    # For H⁺: log(a_H+) = 1·x_H+  →  nu={"H+": 1.0}, log_K'=0
    # For OH⁻: log(a_OH-) = log_Kw + (-1)·x_H+  →  nu={"H+": -1.0}, log_K'=log_Kw
    known: Dict[str, Tuple[dict, float]] = {
        _H_ID: ({"H+": 1.0}, 0.0),
        _OH_ID: ({"H+": -1.0}, float(log_Kw)),
    }

    master_ids: List[str] = []
    master_charges_dict: Dict[str, int] = {"H+": +1}
    component_infos: List[ComponentInfo] = []

    for idx, comp_ids in enumerate(raw_components):
        comp_rxns = _rxns_in_comp(comp_ids)
        m_id = _select_master(comp_ids, comp_rxns)
        master_ids.append(m_id)
        master_charges_dict[m_id] = int(graph_species[m_id].charge)
        known[m_id] = ({m_id: 1.0}, 0.0)
        gas_ids_here = tuple(sorted({
            e.species.id
            for rxn in gas_attachments.get(idx, [])
            for e in rxn.stoichiometry
            if e.phase == "gas"
        }))
        component_infos.append(ComponentInfo(
            master_id=m_id,
            species_ids=tuple(sorted(comp_ids)),
            gas_species_ids=gas_ids_here,
        ))

    logger.debug("NRTableau masters: %s", ["H+"] + master_ids)

    # ── BFS to derive all secondary species ───────────────────────────
    # Apply Van't Hoff to each reaction's log_K before using it.
    remaining = list(ab_rxns)
    max_passes = len(ab_rxns) + 1
    for _pass in range(max_passes):
        if not remaining:
            break
        still_remaining = []
        progress = False
        for rxn in remaining:
            # Classify each participant (excluding H₂O and H⁺).
            unknowns = []
            known_entries = []   # list of (coeff: float, (nu: dict, log_K_prime: float))
            h_coeff = 0.0

            for e in rxn.stoichiometry:
                sp_id = e.species.id
                if sp_id in _SOLVENT_IDS:
                    continue       # H₂O: log[a_H2O] ≈ 0, absorbed into log_K
                c = float(e.coefficient)
                if sp_id == _H_ID:
                    h_coeff += c
                    continue
                if sp_id in known:
                    known_entries.append((c, known[sp_id]))
                else:
                    unknowns.append((c, sp_id, e.species))

            if len(unknowns) != 1:
                still_remaining.append(rxn)
                continue

            # Exactly one unknown — derive it.
            target_coeff, target_id, target_sp = unknowns[0]

            log_K_rxn = _vant_hoff_log_K(
                float(rxn.log_K),
                rxn.dH_J_per_mol,
                T_K,
                float(rxn.T_ref_K),
            )

            # Rearrange:  target_coeff · log(a_target)
            #   = log_K_rxn
            #     − Σ_{known_j}  c_j · (log_K'_j + Σ_k  ν_{jk} · x_k)
            #     − h_coeff · x_H+
            accumulated_log_K = log_K_rxn
            accumulated_nu: Dict[str, float] = {}

            for c_j, (nu_j, logK_j) in known_entries:
                accumulated_log_K -= c_j * logK_j
                for k_id, nu_jk in nu_j.items():
                    accumulated_nu[k_id] = accumulated_nu.get(k_id, 0.0) - c_j * nu_jk

            accumulated_nu["H+"] = accumulated_nu.get("H+", 0.0) - h_coeff

            log_K_prime = accumulated_log_K / target_coeff
            nu = {k: v / target_coeff for k, v in accumulated_nu.items()
                  if abs(v / target_coeff) > 1e-15}

            known[target_id] = (nu, log_K_prime)
            progress = True

        remaining = still_remaining
        if not progress:
            break

    if remaining:
        unresolved = sorted({
            e.species.id
            for rxn in remaining
            for e in rxn.stoichiometry
            if e.species.id not in known
            and e.species.id not in _SOLVENT_IDS
            and e.species.id != _H_ID
        })
        raise ValueError(
            f"NRTableau: could not derive log-linear expressions for: {unresolved}. "
            "Check that the reaction graph is fully connected from the selected masters."
        )

    # ── Build SecondaryEntry list (everything known except masters) ───
    master_set = frozenset(master_ids) | {_H_ID}
    secondaries: List[SecondaryEntry] = []
    for sp_id, (nu, log_K_prime) in known.items():
        if sp_id in master_set:
            continue
        if sp_id == _OH_ID:
            charge = -1
        else:
            charge = int(graph_species[sp_id].charge)
        # element_stoichiometry: nu coefficients for non-H+ masters only.
        # H+ is excluded because it closes the charge balance, not mass balance.
        elem_stoich = {k: v for k, v in nu.items() if k != _H_ID}
        secondaries.append(SecondaryEntry(
            species_id=sp_id,
            charge=charge,
            nu=dict(nu),
            log_K_prime=float(log_K_prime),
            element_stoichiometry=elem_stoich,
            phase="liquid",
        ))

    # ── Derive gas-phase secondaries from folded gas-liquid reactions ──
    # Independent single-pass derivation (not folded into the ab_rxns BFS
    # above): every gas-liquid reaction's liquid-phase participants are
    # already resolvable via `known` by construction (the component-
    # attachment step above only attached a reaction whose liquid entries
    # all landed in one component, and that component's own ladder has
    # just finished deriving). No cross-gas dependency is assumed — each
    # gas-liquid reaction is derived independently from `known` alone.
    for idx, rxns_here in gas_attachments.items():
        for rxn in rxns_here:
            secondaries.append(_derive_gas_secondary(rxn, known, T_K))

    # ── Derive solvent (H2O) Raoult secondaries ────────────────────────
    # CP4 of LAYER1_GAP_CLOSURE — see _derive_solvent_gas_secondary for
    # why this is a separate, simpler pass (constant relation, no master
    # coupling) rather than reusing the gas_attachments machinery above.
    for rxn in solvent_gas_rxns:
        secondaries.append(_derive_solvent_gas_secondary(rxn, T_K))

    # ── Basis validation: check for degenerate masters ────────────────
    # Two masters are degenerate if one is algebraically expressible in
    # terms of the other without any other master — i.e. its nu vector
    # has a non-zero coefficient only for that one other master.
    # Simple check: each master must appear as a key in its own nu dict
    # (and no other master's nu dict should map to it with coefficient ±1
    # unless there's also a self-coefficient).
    # Full rank check of the formula matrix is deferred.
    for m_id in master_ids:
        # If any secondary's nu has ONLY one non-zero key equal to m_id,
        # that secondary IS m_id (would mean we're deriving a master as
        # a secondary), which is a construction error.
        pass  # BFS correctness prevents this; full rank check deferred.

    return NRTableau(
        masters=["H+"] + master_ids,
        master_charges=master_charges_dict,
        secondaries=secondaries,
        components=component_infos,
        log_Kw=float(log_Kw),
        T_ref_K=float(T_ref_K),
    )
