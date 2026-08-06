# -*- coding: utf-8 -*-
"""CV-native open-loop profiles (C10).

A :class:`Profile` is an open-loop time-varying mutator: each step
the orchestrator invokes ``profile.apply(t_h, sim)`` *before*
advancing the CVs. The profile reads the current time, computes
its target value, applies it through the **Pattern B unchecked
setter** so the lifecycle gate (set by :class:`Simulation.run`)
doesn't fire, and returns a :class:`ProfileRecord` for audit.

Concrete profiles shipped at C10:

* :class:`TemperatureRamp` — linear ramp of ``Phase.T_K`` for one
  or more phases of a target CV.
* :class:`VVMSchedule` — sets ``GasFeed.vvm_min`` on the unique
  :class:`GasFeed` boundary of a target CV.
* :class:`SetpointTrajectory` — sets a named attribute on a
  controller (e.g., ``ph_ctrl.setpoint``) from a waypoint
  schedule.

Note (decision 11 / STATE_UNIFICATION): no ``ChemEnvProfile`` —
the ``chem_env`` dict was deleted system-wide. Time-varying
strong-ion totals (if anyone ever needs them) are an additive
extension via a dedicated ``strong_ion_profile`` hook, not
shipped here.

The Profile Protocol consumes a :class:`Simulation` reference
duck-typed (no formal type import) so this module doesn't
introduce a circular dependency with
:mod:`PyOMES.core.simulation`.
"""
from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from typing import Any, List, Optional, Protocol, Sequence, Tuple, runtime_checkable

from PyOMES.control.actions import ProfileRecord


# ════════════════════════════════════════════════════════════════════════
#  Profile Protocol
# ════════════════════════════════════════════════════════════════════════


@runtime_checkable
class Profile(Protocol):
    """Open-loop time-varying mutator.

    Implementations:

    * Carry their waypoint / schedule on ``self``.
    * Implement ``apply(t_h, sim) -> ProfileRecord`` — orchestrator
      calls this once per step *before* CV advancement so the
      mutated state is what the integration sees.
    * Mutate state through Pattern B unchecked setters
      (``phase._set_T_K_unchecked``, etc.) so the lifecycle gate
      isn't tripped during the run.

    The formal Protocol declares only ``apply``; ``label`` is
    expected (duck-typed) but not enforced.
    """

    def apply(self, t_h: float, sim: Any) -> ProfileRecord: ...


# ════════════════════════════════════════════════════════════════════════
#  Helpers
# ════════════════════════════════════════════════════════════════════════


def _interp_waypoints(
    t_h: float, waypoints: Sequence[Tuple[float, float]],
) -> float:
    """Linear interpolation between sorted (t_h, value) waypoints.

    Before the first waypoint: returns the first value.
    After the last waypoint: returns the last value (hold).
    Empty waypoints: returns 0.0.
    """
    if not waypoints:
        return 0.0
    pts = sorted(waypoints, key=lambda p: p[0])
    if t_h <= pts[0][0]:
        return float(pts[0][1])
    if t_h >= pts[-1][0]:
        return float(pts[-1][1])
    # Binary search for the bracket.
    times = [p[0] for p in pts]
    idx = bisect.bisect_right(times, float(t_h))
    t0, v0 = pts[idx - 1]
    t1, v1 = pts[idx]
    if t1 == t0:
        return float(v0)
    frac = (float(t_h) - t0) / (t1 - t0)
    return float(v0) + float(frac) * (float(v1) - float(v0))


# ════════════════════════════════════════════════════════════════════════
#  TemperatureRamp
# ════════════════════════════════════════════════════════════════════════


@dataclass
class TemperatureRamp:
    """Ramp one or more phases' ``T_K`` linearly between waypoints.

    Parameters
    ----------
    target_cv_key : str
        Which CV in ``sim.cvs`` to mutate.
    waypoints : sequence of (t_h, T_K) tuples
        Times and target temperatures. Sorted internally.
    phases : tuple of str
        Phase keys on the target CV to mutate. Default
        ``("liquid", "gas")`` covers the common fermenter shape.
    label : str
        Human-readable identifier propagated to
        :class:`ProfileRecord`.

    Notes
    -----
    Mutation uses the ``MutableScalar._set_unchecked`` descriptor path
    (Pattern B). Phases missing from the CV are silently skipped.
    """

    target_cv_key: str
    waypoints: List[Tuple[float, float]]
    phases: Tuple[str, ...] = ("liquid", "gas")
    label: str = "temperature_ramp"

    def apply(self, t_h: float, sim: Any) -> ProfileRecord:
        from ..control.descriptors import MutableScalar
        T_K = _interp_waypoints(t_h, self.waypoints)
        cv = sim.cvs.get(self.target_cv_key)
        targets: dict = {}
        if cv is not None:
            for pkey in self.phases:
                phase = cv.phases.get(pkey)
                if phase is None:
                    continue
                # Pattern B: bypass the lifecycle gate via descriptor.
                _d = type(phase).__dict__.get("T_K")
                if isinstance(_d, MutableScalar):
                    _d._set_unchecked(phase, float(T_K))
                targets[f"cv.{self.target_cv_key}.phases.{pkey}.T_K"] = float(T_K)
        return ProfileRecord(
            profile_label=self.label,
            t_h=float(t_h),
            targets=targets,
        )


# ════════════════════════════════════════════════════════════════════════
#  VVMSchedule
# ════════════════════════════════════════════════════════════════════════


@dataclass
class VVMSchedule:
    """Schedule the unique :class:`GasFeed`'s ``vvm_min`` on a CV.

    Parameters
    ----------
    target_cv_key : str
        Which CV in ``sim.cvs`` to mutate.
    waypoints : sequence of (t_h, vvm_min) tuples
        Times and target vvm values. Sorted internally.
    label : str
        Human-readable identifier propagated to
        :class:`ProfileRecord`.

    Notes
    -----
    Finds the *first* :class:`~PyOMES.core.boundaries.GasFeed` in
    ``cv.boundaries`` and writes ``vvm_min`` directly. ``GasFeed``
    is not currently lifecycle-gated, so a plain attribute write
    suffices.
    """

    target_cv_key: str
    waypoints: List[Tuple[float, float]]
    label: str = "vvm_schedule"

    def apply(self, t_h: float, sim: Any) -> ProfileRecord:
        from PyOMES.core.boundaries import GasFeed
        vvm = _interp_waypoints(t_h, self.waypoints)
        cv = sim.cvs.get(self.target_cv_key)
        targets: dict = {}
        if cv is not None:
            for boundary in cv.boundaries:
                if isinstance(boundary, GasFeed):
                    from ..control.descriptors import MutableScalar
                    _d = type(boundary).__dict__.get("vvm_min")
                    if isinstance(_d, MutableScalar):
                        _d._set_unchecked(boundary, float(vvm))
                    targets[
                        f"cv.{self.target_cv_key}.gas_feed.vvm_min"
                    ] = float(vvm)
                    break
        return ProfileRecord(
            profile_label=self.label,
            t_h=float(t_h),
            targets=targets,
        )


# ════════════════════════════════════════════════════════════════════════
#  SetpointTrajectory
# ════════════════════════════════════════════════════════════════════════


@dataclass
class SetpointTrajectory:
    """Vary a named attribute on a controller over time.

    The classic use case is a pH-setpoint schedule that varies
    the setpoint of a :class:`PHController` during fermentation
    (e.g., ramp from 6.5 to 5.0 over the first 24 hours).

    Parameters
    ----------
    controller : object
        The controller whose attribute will be set. Direct
        attribute write — controllers are not lifecycle-gated.
    attribute : str
        Name of the attribute to set (e.g., ``"setpoint"``,
        ``"setpoint_mol_L"``).
    waypoints : sequence of (t_h, value) tuples
        Times and target values for the named attribute.
    label : str
        Human-readable identifier propagated to
        :class:`ProfileRecord`.
    """

    controller: Any
    attribute: str
    waypoints: List[Tuple[float, float]]
    label: str = "setpoint_trajectory"

    def apply(self, t_h: float, sim: Any) -> ProfileRecord:
        value = _interp_waypoints(t_h, self.waypoints)
        setattr(self.controller, self.attribute, float(value))
        ctrl_label = getattr(self.controller, "label", "?")
        return ProfileRecord(
            profile_label=self.label,
            t_h=float(t_h),
            targets={
                f"controller.{ctrl_label}.{self.attribute}": float(value),
            },
        )
