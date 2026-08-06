# -*- coding: utf-8 -*-
"""Per-entry activity dispatch — phase-aware, evaluate-at-a-point helper.

:func:`activity_for_entry` gives one :class:`~PyOMES.reactions.stoichiometry.StoichiometryEntry`
its activity/fugacity correction for whatever phase it lives in:

- ``phase == "liquid"`` → ``thermo.liquid_activity.gamma_all(...)`` (activity
  coefficient γ_i).
- ``phase == "gas"`` → ``thermo.gas_eos.partial_pressures_atm(...)`` (partial
  pressure p_i, atm).
- ``phase == "solid"`` → ``1.0`` (pure-solid convention: solid activity is 1).

This is groundwork for Phase 2 of Layer 1 gap closure
(``LAYER1_GAP_CLOSURE.md``), which wires this dispatch into
``build_tableau()``'s Newton residual so gas-liquid/solid-liquid rows can
be solved simultaneously with acid-base rows. This module only proves the
dispatch is correct in isolation — it is not called from
``build_tableau()``/``solve_nr()`` yet.

Why ``charge`` is a separate keyword, not derived internally
--------------------------------------------------------------
``LiquidPhaseModel.gamma_all`` needs a ``{species_id: charge}`` map
covering every ionic species in the liquid composition to compute ionic
strength (Davies/SIT) — but ``Phase`` objects store only ``n_mol`` floats,
no ``Species``/charge metadata. That map has to come from wherever the
reaction declarations are (each ``StoichiometryEntry.species.charge``),
which for standalone use here means the caller supplies it explicitly.
Phase 2's wiring passes the full charge map already available on the
built ``NRTableau`` (``tableau.master_charges``); without a caller-
supplied map, this function falls back to the single ``entry`` species'
own charge, matching every other species in ``x_mol`` to ``0`` — correct
mass-action behaviour for that species alone, but underestimates ionic
strength if other charged species are present. Pass ``charge=`` explicitly
whenever more than one ionic species is in play.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Dict, Optional

if TYPE_CHECKING:
    from ..core.phases import Phase
    from ..reactions.stoichiometry import StoichiometryEntry
    from ..thermo import ThermoFramework

_VALID_PHASES = frozenset({"liquid", "gas", "solid"})


def activity_for_entry(
    entry: "StoichiometryEntry",
    phases: Dict[str, "Phase"],
    thermo: "ThermoFramework",
    T_K: float,
    *,
    charge: Optional[Dict[str, int]] = None,
) -> float:
    """Dispatch the activity/fugacity correction for one stoichiometry entry.

    Parameters
    ----------
    entry : StoichiometryEntry
        The participant to evaluate. ``entry.phase`` selects the dispatch
        branch; ``entry.species.id`` selects which value to return from
        the per-species result dict.
    phases : dict
        ``{phase_key: Phase}``, keyed the same way as ``entry.phase``
        (e.g. ``{"gas": GasPhase(...), "liquid": LiquidPhase(...)}``).
        Unused for solid entries.
    thermo : ThermoFramework
        Supplies ``liquid_activity`` (for liquid entries) and ``gas_eos``
        (for gas entries; falls back to :class:`~PyOMES.equilibria.vle.IdealGasEOS`
        when ``thermo.gas_eos`` is ``None`` per the framework's own
        documented convention).
    T_K : float
        Temperature (K).
    charge : dict, optional
        ``{species_id: charge}`` covering every ionic species present in
        the liquid phase's ``n_mol`` — required for a correct ionic-
        strength computation whenever more than one ionic species is
        present. Defaults to just ``{entry.species.id: entry.species.charge}``
        when omitted (see module docstring).

    Returns
    -------
    float
        ``γ_i`` for a liquid entry, ``p_i`` (atm) for a gas entry, or
        ``1.0`` for a solid entry.

    Raises
    ------
    ValueError
        If ``entry.phase`` is not one of ``"liquid"``, ``"gas"``, ``"solid"``.
    KeyError
        If ``entry.phase`` is ``"liquid"``/``"gas"`` and that key is
        missing from ``phases``.
    """
    if entry.phase == "solid":
        return 1.0

    if entry.phase not in _VALID_PHASES:
        raise ValueError(
            f"activity_for_entry: unrecognised phase {entry.phase!r} for "
            f"entry {entry!r}. Expected one of {sorted(_VALID_PHASES)}."
        )

    phase_obj = phases[entry.phase]
    species_id = entry.species.id

    if entry.phase == "liquid":
        charge_map = (
            charge if charge is not None
            else {species_id: int(entry.species.charge)}
        )
        gammas = thermo.liquid_activity.gamma_all(
            phase_obj.n_mol, T_K, charge=charge_map,
        )
        return gammas.get(species_id, 1.0)

    # entry.phase == "gas"
    gas_eos = thermo.gas_eos
    if gas_eos is None:
        from ..equilibria.vle import IdealGasEOS
        gas_eos = IdealGasEOS()
    pressures = gas_eos.partial_pressures_atm(
        phase_obj.n_mol, T_K=T_K, V_L=phase_obj.V_L,
    )
    return pressures.get(species_id, 0.0)
