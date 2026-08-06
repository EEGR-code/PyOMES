# -*- coding: utf-8 -*-
"""Time-profile infrastructure for dynamic operating conditions.

Stage 17a: Composable time profiles that modify simulation parameters
(feed rates, temperature, pressure setpoints, etc.) as functions of time.

Profile types
-------------
- :class:`ConstantProfile` — fixed value (baseline)
- :class:`StepProfile` — step changes at specified times
- :class:`RampProfile` — linear interpolation between waypoints
- :class:`TableProfile` — arbitrary time series (numpy-based)
- :class:`PeriodicProfile` — repeating pattern (diurnal, intermittent)
- :class:`CompositeProfile` — sum or product of other profiles

Dispatcher
----------
- :class:`ProfileSet` — collects named profiles and applies them to
  target objects at each timestep via :meth:`apply`.

Usage
-----
>>> from PyOMES.profiles import ProfileSet, RampProfile, StepProfile
>>> ps = ProfileSet()
>>> ps.add("liquid_feed.Q_L_per_h", RampProfile([(0, 0), (48, 100)]))
>>> ps.add("T_K", StepProfile([(0, 308.15), (720, 328.15)]))
>>> # In run_batch, ps.apply(t_h, cv) is called before each advance()
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np


# ════════════════════════════════════════════════════════════════════════
#  Profile protocol and implementations
# ════════════════════════════════════════════════════════════════════════

class TimeProfile:
    """Abstract base for time profiles.

    All profiles are callables: ``value = profile(t_h)``.
    """

    def __call__(self, t_h: float) -> float:
        raise NotImplementedError


@dataclass(frozen=True)
class ConstantProfile(TimeProfile):
    """Fixed value at all times.

    Parameters
    ----------
    value : float
        The constant value returned for any *t_h*.
    """
    value: float

    def __call__(self, t_h: float) -> float:
        return self.value

    def __repr__(self):
        return f"ConstantProfile({self.value})"


@dataclass(frozen=True)
class StepProfile(TimeProfile):
    """Step changes at specified times.

    Returns the value of the most recent step at or before *t_h*.
    Before the first step time, returns the first step's value.

    Parameters
    ----------
    steps : sequence of (t_h, value) tuples
        Must be sorted by time (ascending).  At least one step required.

    Example
    -------
    >>> p = StepProfile([(0, 100), (24, 200), (48, 150)])
    >>> p(12)   # 100 (between t=0 and t=24)
    >>> p(30)   # 200 (between t=24 and t=48)
    >>> p(100)  # 150 (after t=48)
    """
    steps: Tuple[Tuple[float, float], ...]

    def __init__(self, steps: Sequence[Tuple[float, float]]):
        sorted_steps = tuple(sorted(steps, key=lambda s: s[0]))
        if len(sorted_steps) == 0:
            raise ValueError("StepProfile requires at least one step")
        object.__setattr__(self, "steps", sorted_steps)
        object.__setattr__(self, "_times", [s[0] for s in sorted_steps])

    def __call__(self, t_h: float) -> float:
        idx = bisect.bisect_right(self._times, float(t_h)) - 1
        if idx < 0:
            return self.steps[0][1]
        return self.steps[idx][1]

    def __repr__(self):
        return f"StepProfile({len(self.steps)} steps)"


@dataclass(frozen=True)
class RampProfile(TimeProfile):
    """Linear interpolation between waypoints.

    Between waypoints, values are linearly interpolated.
    Before the first waypoint, returns the first value.
    After the last waypoint, returns the last value (hold).

    Parameters
    ----------
    points : sequence of (t_h, value) tuples
        Must contain at least 2 points.  Sorted internally.

    Example
    -------
    >>> p = RampProfile([(0, 0), (24, 100), (48, 100), (72, 0)])
    >>> p(12)   # 50.0 (linear ramp from 0 to 100)
    >>> p(36)   # 100.0 (flat hold)
    >>> p(60)   # 50.0 (linear ramp from 100 to 0)
    """
    points: Tuple[Tuple[float, float], ...]

    def __init__(self, points: Sequence[Tuple[float, float]]):
        sorted_pts = tuple(sorted(points, key=lambda p: p[0]))
        if len(sorted_pts) < 2:
            raise ValueError("RampProfile requires at least 2 points")
        object.__setattr__(self, "points", sorted_pts)
        object.__setattr__(self, "_times", np.array([p[0] for p in sorted_pts]))
        object.__setattr__(self, "_values", np.array([p[1] for p in sorted_pts]))

    def __call__(self, t_h: float) -> float:
        return float(np.interp(float(t_h), self._times, self._values))

    def __repr__(self):
        return f"RampProfile({len(self.points)} points)"


@dataclass(frozen=True)
class TableProfile(TimeProfile):
    """Arbitrary time series with numpy interpolation.

    Parameters
    ----------
    t_h : array-like
        Time points (hours).  Must be monotonically increasing.
    values : array-like
        Values at each time point.  Same length as *t_h*.
    """
    t_h_arr: np.ndarray
    values_arr: np.ndarray

    def __init__(self, t_h: Sequence[float], values: Sequence[float]):
        t_arr = np.asarray(t_h, dtype=float)
        v_arr = np.asarray(values, dtype=float)
        if len(t_arr) != len(v_arr):
            raise ValueError("t_h and values must have same length")
        if len(t_arr) < 1:
            raise ValueError("TableProfile requires at least 1 point")
        object.__setattr__(self, "t_h_arr", t_arr)
        object.__setattr__(self, "values_arr", v_arr)

    def __call__(self, t_h: float) -> float:
        return float(np.interp(float(t_h), self.t_h_arr, self.values_arr))

    def __repr__(self):
        return f"TableProfile({len(self.t_h_arr)} points)"


@dataclass(frozen=True)
class PeriodicProfile(TimeProfile):
    """Repeating pattern with a fixed period.

    Wraps a base profile so that ``t_h`` is taken modulo ``period_h``
    before evaluation.

    Parameters
    ----------
    period_h : float
        Period (hours).
    base_profile : TimeProfile
        Profile to evaluate on the wrapped time.

    Example
    -------
    >>> daily = RampProfile([(0, 0), (12, 100), (24, 0)])
    >>> p = PeriodicProfile(24.0, daily)
    >>> p(36)   # same as daily(12) = 100
    """
    period_h: float
    base_profile: TimeProfile

    def __call__(self, t_h: float) -> float:
        return self.base_profile(float(t_h) % self.period_h)

    def __repr__(self):
        return f"PeriodicProfile(period={self.period_h}h, {self.base_profile!r})"


@dataclass(frozen=True)
class CompositeProfile(TimeProfile):
    """Combine multiple profiles by summation or multiplication.

    Parameters
    ----------
    profiles : sequence of TimeProfile
        Profiles to combine.
    mode : str
        ``"sum"`` (default) or ``"product"``.

    Example
    -------
    >>> baseline = ConstantProfile(100)
    >>> perturbation = RampProfile([(0, 0), (24, 20), (48, 0)])
    >>> p = CompositeProfile([baseline, perturbation], mode="sum")
    >>> p(12)   # 110.0
    """
    profiles: Tuple[TimeProfile, ...]
    mode: str = "sum"

    def __init__(self, profiles: Sequence[TimeProfile], mode: str = "sum"):
        if mode not in ("sum", "product"):
            raise ValueError(f"mode must be 'sum' or 'product', got {mode!r}")
        object.__setattr__(self, "profiles", tuple(profiles))
        object.__setattr__(self, "mode", mode)

    def __call__(self, t_h: float) -> float:
        values = [p(float(t_h)) for p in self.profiles]
        if self.mode == "sum":
            return sum(values)
        else:
            result = 1.0
            for v in values:
                result *= v
            return result

    def __repr__(self):
        return f"CompositeProfile({len(self.profiles)} profiles, mode={self.mode!r})"


# ════════════════════════════════════════════════════════════════════════
#  ProfileTarget — describes what a profile controls
# ════════════════════════════════════════════════════════════════════════

@dataclass
class ProfileTarget:
    """Describes how to apply a profile value to a simulation object.

    Parameters
    ----------
    obj : object
        The mutable target (e.g., a LiquidFeed boundary, a GasPhase).
    attr : str
        Simple attribute name on *obj* (e.g., ``"Q_L_per_h"``).
    key : str or None
        If the attribute is a dict, the key within it
        (e.g., ``"AceticAcid"`` for ``feed_conc_mol_L["AceticAcid"]``).
    """
    obj: Any
    attr: str
    key: Optional[str] = None

    def apply(self, value: float) -> None:
        """Set the target to *value*."""
        if self.key is not None:
            getattr(self.obj, self.attr)[self.key] = value
        else:
            setattr(self.obj, self.attr, value)

    def current(self) -> float:
        """Read the current value."""
        if self.key is not None:
            return getattr(self.obj, self.attr).get(self.key, 0.0)
        return getattr(self.obj, self.attr)


# ════════════════════════════════════════════════════════════════════════
#  ProfileSet — collects and applies profiles
# ════════════════════════════════════════════════════════════════════════

class ProfileSet:
    """Collection of named time profiles with targets.

    Each profile maps a name to a ``(ProfileTarget, TimeProfile)`` pair.
    Calling :meth:`apply` evaluates every profile at the current time
    and sets the corresponding target attribute.

    state-unification C4e: the chem_env-override path
    (``add_chem_env`` / ``apply_chem_env``) was removed alongside
    the chem_env dict itself. Time-varying chemistry state (strong
    ions, totals) should be threaded via apply_flux on phase.n_mol
    entries from a boundary or step_callback instead.

    Usage
    -----
    >>> ps = ProfileSet()
    >>> ps.add_target("feed_Q", target, RampProfile([(0, 0), (48, 100)]))
    >>> ps.apply(t_h=24.0)  # sets target attributes
    """

    def __init__(self):
        self._profiles: Dict[str, Tuple[ProfileTarget, TimeProfile]] = {}

    def add_target(
        self,
        name: str,
        target: ProfileTarget,
        profile: TimeProfile,
    ) -> "ProfileSet":
        """Register a named profile with an explicit target.

        Parameters
        ----------
        name : str
            Unique identifier (for diagnostics and removal).
        target : ProfileTarget
            What to set when the profile is evaluated.
        profile : TimeProfile
            The time-varying value.

        Returns
        -------
        ProfileSet
            self (for chaining).
        """
        self._profiles[name] = (target, profile)
        return self

    def add(
        self,
        name: str,
        obj: Any,
        attr: str,
        profile: TimeProfile,
        key: Optional[str] = None,
    ) -> "ProfileSet":
        """Convenience: register a profile by object + attribute name.

        Parameters
        ----------
        name : str
            Unique identifier.
        obj : object
            The mutable target object.
        attr : str
            Attribute on *obj* to set.
        profile : TimeProfile
            The time-varying value.
        key : str or None
            Dict key if *attr* is a dict.

        Returns
        -------
        ProfileSet
            self (for chaining).
        """
        return self.add_target(name, ProfileTarget(obj, attr, key), profile)

    def apply(self, t_h: float) -> Dict[str, float]:
        """Evaluate all attribute profiles at *t_h* and apply to targets.

        Returns
        -------
        dict
            ``{name: value}`` — the values that were applied.
        """
        applied = {}
        for name, (target, profile) in self._profiles.items():
            value = profile(float(t_h))
            target.apply(value)
            applied[name] = value
        return applied

    def evaluate(self, t_h: float) -> Dict[str, float]:
        """Evaluate all profiles without applying (read-only).

        Returns
        -------
        dict
            ``{name: value}`` — what would be applied (attributes + context).
        """
        return {
            name: profile(float(t_h))
            for name, (_, profile) in self._profiles.items()
        }

    def remove(self, name: str) -> None:
        """Remove a named profile."""
        self._profiles.pop(name, None)

    @property
    def names(self) -> List[str]:
        """List of registered profile names."""
        return list(self._profiles.keys())

    def __len__(self) -> int:
        return len(self._profiles)

    def __repr__(self):
        return f"ProfileSet({len(self._profiles)} attribute profiles)"


# ════════════════════════════════════════════════════════════════════════
#  Stage 17c: Specialised targets for environmental profiles
# ════════════════════════════════════════════════════════════════════════

class TemperatureProfileTarget(ProfileTarget):
    """Sets temperature (K) on both gas/liquid phases and speciation engine.

    When temperature changes, this target updates:
    - ``gas_phase.T_K`` — affects ideal gas law, partial pressures
    - ``liquid_phase.T_K`` — affects reaction environment T_K (Arrhenius)
    - Speciation engine ``T_C`` — affects pKa temperature corrections
    - Henry constants are recomputed per-step by the link (Stage 8)

    Parameters
    ----------
    cv : ControlVolume
        The fermenter CV to control.
    """

    def __init__(self, cv: Any):
        self._cv = cv

    def apply(self, value: float) -> None:
        T_K = float(value)
        self._cv.phases["gas"].T_K = T_K
        self._cv.phases["liquid"].T_K = T_K
        # Update speciation engine temperature (for pKa corrections).
        # state-unification C4d: engine attaches via
        # cv.reaction_system.attach_engine; no more property_solvers
        # list to iterate.
        T_C = T_K - 273.15
        rxn_system = getattr(self._cv, "reaction_system", None)
        if rxn_system is not None:
            engine = getattr(rxn_system, "engine", None)
            if engine is not None:
                engine.T_C = T_C
                model = getattr(engine, "model", None)
                if model is not None:
                    model.T_C = T_C

    def current(self) -> float:
        return self._cv.phases["gas"].T_K


class PressureSetpointTarget(ProfileTarget):
    """Sets pressure setpoint on a PressureReliefVent boundary.

    Parameters
    ----------
    vent : PressureReliefVent
        The pressure relief vent boundary to control.
    """

    def __init__(self, vent: Any):
        self._vent = vent

    def apply(self, value: float) -> None:
        self._vent.P_set_atm = float(value)

    def current(self) -> float:
        return self._vent.P_set_atm

