"""Package-level configuration surface.

Holds the global ``config`` singleton with accuracy-warning
thresholds and throttle behaviour. Initialised from the
``VLSIM_WARNINGS`` environment variable at import time (flat
presets only — ``silent`` / ``verbose`` / ``production``).

Fine-grained tuning is via attribute assignment:

    import PyOMES
    PyOMES.config.warnings.pH_change_threshold = 0.5
    PyOMES.config.warnings.throttle = "silent"

Or via a preset constructor:

    PyOMES.config.warnings = PyOMES.WarningConfig.production()

CI / batch usage:

    $ VLSIM_WARNINGS=silent python run_simulation.py
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field
from typing import Literal


@dataclass
class WarningConfig:
    """Global accuracy-warning configuration.

    Thresholds and throttle behaviour for ``AccuracyWarning``
    emissions. Modify ``PyOMES.config.warnings`` before running
    simulations.
    """

    # Thresholds — per check
    pH_change_threshold: float = 0.3            # pH units per step
    newton_iters_threshold: int = 15            # iters per speciation solve
    charge_residual_threshold: float = 1e-8     # |sum z_i C_i|
    dt_over_tau_min_threshold: float = 0.1      # dt / tau_min
    scipy_rejection_threshold: float = 0.3      # rejected / accepted

    # Ionic-strength regime (migrated from speciation/engine.py)
    ionic_strength_ideal_threshold: float = 0.10   # mol/L
    ionic_strength_davies_threshold: float = 0.50  # mol/L

    # Clamping (STEP_SOLVER_INTERFACE_REFINEMENT.md item 8)
    negative_mole_tolerance: float = 1e-9    # mol; only checked when clamp_fn=None
    clamp_invoked_tolerance: float = 1e-9    # mol/h; rate change from clamp_fn

    # Conservation (state-unification C6) — relative-plus-max-clamp
    # factors applied to ``max(total_element_mol, 1)`` so dilute
    # systems do not spuriously trip on noise-level absolute
    # residuals. Per-step ratio guards against single-step
    # stoichiometry breakage; cumulative ratio guards against
    # accumulating float roundoff over a long run.
    conservation_per_step_factor: float = 1e-8
    conservation_cumulative_factor: float = 1e-6

    # Throttle — controls how often each distinct category is emitted
    throttle: Literal["once", "first_N", "always", "silent"] = "once"
    first_N: int = 10

    @classmethod
    def silent(cls) -> "WarningConfig":
        return cls(throttle="silent")

    @classmethod
    def verbose(cls) -> "WarningConfig":
        return cls(throttle="always")

    @classmethod
    def production(cls) -> "WarningConfig":
        return cls(throttle="first_N", first_N=3)


@dataclass
class Config:
    warnings: WarningConfig = field(default_factory=WarningConfig)


def _warning_config_from_env() -> WarningConfig:
    """Construct a ``WarningConfig`` from the ``VLSIM_WARNINGS``
    environment variable, if set.

    Recognised values (case-insensitive, whitespace-stripped):
    ``silent``, ``verbose``, ``production``. Empty string is
    treated as unset. Unknown values emit a single ``UserWarning``
    and fall back to the default config.
    """
    raw = os.environ.get("VLSIM_WARNINGS")
    if raw is None:
        return WarningConfig()
    preset = raw.strip().lower()
    if not preset:
        return WarningConfig()
    if preset == "silent":
        return WarningConfig.silent()
    if preset == "verbose":
        return WarningConfig.verbose()
    if preset == "production":
        return WarningConfig.production()
    warnings.warn(
        f"VLSIM_WARNINGS={raw!r}: unknown preset. "
        f"Expected one of 'silent', 'verbose', 'production'. "
        f"Falling back to default WarningConfig.",
        UserWarning,
        stacklevel=2,
    )
    return WarningConfig()


config = Config(warnings=_warning_config_from_env())
