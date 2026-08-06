# -*- coding: utf-8 -*-
"""Structured per-step records produced by controllers and profiles.

This module hosts the frozen dataclasses used by the per-step
recorder channel (``BatchResult.controller_actions`` and
``profile_actions``) and emitted by the controller / profile
protocols at C7-C10.

Stubs at C3; the full :class:`ControlAction` schema lands at C7
when the controller-protocol port begins. The :class:`ProfileRecord`
schema lands at C10. The schemas were chosen at design time
(SIMULATION_CLASS.md decision 2 + the ControlAction record section);
the fields enumerated here match those decisions.

Note: this module is **not** ``src/sim/control.py`` — that module
is part of the legacy ``CUFermentationSpeciation`` island
(``src/sim/`` is out of scope for this phase). The new framework's
:class:`ControlAction` lives here so it stays close to the
controller port in ``src/control/``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict


@dataclass(frozen=True)
class ControlAction:
    """Auditable record of one controller's per-step action.

    A frozen record — the fields are populated by the controller's
    :meth:`compute` method (or by the orchestrator's apply path
    for orchestrator-mediated additions). The dataclass itself is
    immutable at the attribute level; the contained dicts are not
    deeply frozen (Python convention: trust the producer).

    Attributes
    ----------
    controller_label : str
        Human-readable identifier of the producing controller.
    target_cv_key : str
        CV key the action targets. Empty for multi-CV controllers
        that act on more than one CV.
    t_h : float
        Wall-clock time at which the action was computed.
    dt_h : float
        Step duration this action covers.
    flux_applied : dict
        ``{phase_key: {species_id: mol_per_h}}`` — external fluxes
        the orchestrator will apply via
        :meth:`ControlVolume.apply_external_flux` after this action
        is recorded.
    params_changed : dict
        ``{param_path: new_value}`` — parameter mutations the
        orchestrator will route through the unchecked setter path
        (C9; Pattern B). Param paths follow
        ``"cv_key.attr.subattr"`` syntax resolved by
        ``_resolve_param_path`` at C9.
    vented_mol : dict
        ``{species_id: total_mol_vented_over_dt}`` — net moles
        leaving the system through this action (informational;
        useful for pressure-relief and outlet diagnostics).
    dosed_mol : dict
        ``{compound_id: total_mol_dosed_over_dt}`` — net moles
        added through this action (acid/base dosing audit).

    Notes
    -----
    This is the stub form at C3 — the dataclass is defined so the
    :class:`BatchResult` schema and the recorder can reference it.
    Controllers do not produce real :class:`ControlAction` instances
    until C8 lands; until then the recorder always gets empty
    lists.
    """

    controller_label: str = ""
    target_cv_key: str = ""
    t_h: float = 0.0
    dt_h: float = 0.0
    flux_applied: Dict[str, Dict[str, float]] = field(default_factory=dict)
    params_changed: Dict[str, Any] = field(default_factory=dict)
    vented_mol: Dict[str, float] = field(default_factory=dict)
    dosed_mol: Dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class ProfileRecord:
    """Auditable record of one profile's per-step action.

    A frozen record produced by a :class:`Profile`'s ``apply()``
    method. Reports what was mutated and to what value.

    Attributes
    ----------
    profile_label : str
        Human-readable identifier of the producing profile.
    t_h : float
        Wall-clock time at which the profile applied.
    targets : dict
        ``{target_path: new_value}`` — what was changed and to
        what. Exact key conventions settled when concrete profiles
        port at C10.

    Notes
    -----
    Stub at C3; full schema lands at C10.
    """

    profile_label: str = ""
    t_h: float = 0.0
    targets: Dict[str, Any] = field(default_factory=dict)
