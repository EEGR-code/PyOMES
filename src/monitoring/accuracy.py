"""``AccuracyMonitor`` and ``AccuracyWarning``.

Each ``ControlVolume`` owns one ``AccuracyMonitor`` (constructed
in ``ControlVolume.__init__``) that runs cheap per-step checks:
pH change, Newton iterations, charge residual, dt vs tau_min,
scipy step rejections, and ionic-strength regime.

Thresholds and throttle mode are read off
``PyOMES.config.warnings`` at every emission — configuration
changes propagate without re-instantiation.

A module-level ``_summary_counter`` aggregates emissions across
CVs so multi-zone runs get one meaningful end-of-simulation
summary line.
"""

from __future__ import annotations

import math
import warnings
from collections import Counter
from typing import Optional, Set


class AccuracyWarning(UserWarning):
    """Emitted when a numerical accuracy check flags a potential
    issue.

    Use Python's standard ``warnings`` machinery to filter or
    promote::

        import warnings
        from PyOMES import AccuracyWarning

        warnings.simplefilter("ignore", AccuracyWarning)  # silence
        warnings.simplefilter("error",  AccuracyWarning)  # promote
    """


# Module-level emission counter, keyed by category. Persists across
# CVs so a MultiCVSystem run gets one meaningful summary line.
# Pytest sessions that need a clean slate between tests should call
# ``reset_accuracy_summary()`` from a fixture.
_summary_counter: Counter = Counter()


def _get_warning_config():
    """Read the current ``WarningConfig`` from the package singleton.

    Imports lazily because the ``PyOMES`` package re-exports the
    ``config`` instance as a sibling attribute, which shadows the
    submodule of the same name; the direct ``from PyOMES.config
    import config`` form sidesteps that and reads the live
    singleton each call so users can mutate
    ``PyOMES.config.warnings`` between simulations without
    re-instantiating monitors.
    """
    from PyOMES.config import config as _config_singleton

    return _config_singleton.warnings


def print_accuracy_summary() -> None:
    """Emit the end-of-run summary line for the ``"once"`` throttle.

    Prints one line listing per-category emission counts; no-op if
    nothing was emitted. Counts come from the module-level
    ``_summary_counter`` so multi-CV runs aggregate naturally.
    """
    if not _summary_counter:
        return
    parts = [f"{cat}={count}" for cat, count in sorted(_summary_counter.items())]
    total = sum(_summary_counter.values())
    print(
        f"AccuracyWarning summary: {total} emission(s) across "
        f"[{', '.join(parts)}]"
    )


def reset_accuracy_summary() -> None:
    """Clear the module-level summary counter.

    Useful between pytest tests, or between distinct simulation
    runs in the same process when the counter would otherwise
    confuse the ``"first_N"`` throttle.
    """
    _summary_counter.clear()


class AccuracyMonitor:
    """Per-CV accuracy check runner.

    Holds per-CV step history (``last_pH``, ``step_count``,
    ``_warned_categories``). Thresholds and throttle mode are
    *not* cached here — they're read off
    ``PyOMES.config.warnings`` at every ``_emit`` call so
    configuration changes between simulations take effect
    without re-instantiation.

    Per-check methods (``check_pH_jump``, ``check_newton_iters``,
    etc.) land in checkpoint 3. This skeleton ships ``_emit`` and
    ``reset`` only — enough for checkpoint 2's tests to exercise
    the throttle plumbing end-to-end.
    """

    def __init__(self) -> None:
        self.last_pH: Optional[float] = None
        self.step_count: int = 0
        self._warned_categories: Set[str] = set()

    def reset(self) -> None:
        """Clear per-monitor state.

        Does NOT touch the module-level ``_summary_counter`` —
        call ``reset_accuracy_summary()`` for that.
        """
        self.last_pH = None
        self.step_count = 0
        self._warned_categories.clear()

    def _emit(self, category: str, message: str) -> None:
        """Apply the configured throttle then dispatch
        ``warnings.warn``.

        Throttle semantics:

        - ``"silent"``: no emissions
        - ``"once"``: each category emits at most once per monitor
          (per-CV semantics — multiple CVs each emit once)
        - ``"first_N"``: the module-level counter caps total
          emissions across all monitors at ``first_N``
        - ``"always"``: every call emits

        The summary counter is incremented on every actual
        emission (post-throttle), so ``print_accuracy_summary``
        reflects what users saw.
        """
        cfg = _get_warning_config()
        throttle = cfg.throttle
        if throttle == "silent":
            return
        if throttle == "once" and category in self._warned_categories:
            return
        if throttle == "first_N":
            # Global cap across all monitors — preserves "you'll
            # see the first N occurrences total" semantics for
            # multi-CV runs.
            if _summary_counter[category] >= cfg.first_N:
                return
        self._warned_categories.add(category)
        _summary_counter[category] += 1
        warnings.warn(message, AccuracyWarning, stacklevel=3)

    # ── Per-step checks ──────────────────────────────────────────

    def check_pH_jump(self, current_pH: Optional[float]) -> None:
        """Flag large pH shifts between consecutive steps.

        A jump larger than ``WarningConfig.pH_change_threshold``
        suggests that operator splitting is losing information,
        or that the kinetic step is too coarse for the pH-shifting
        dynamics. No-op on the first call (no prior pH to compare),
        on ``None`` (the speciation solver didn't produce a pH),
        and on NaN.
        """
        self.step_count += 1
        if current_pH is None:
            return
        try:
            current = float(current_pH)
        except (TypeError, ValueError):
            return
        if math.isnan(current):
            return
        prior = self.last_pH
        self.last_pH = current
        if prior is None:
            return
        delta = abs(current - prior)
        threshold = _get_warning_config().pH_change_threshold
        if delta > threshold:
            self._emit(
                "pH_change",
                f"pH changed by {delta:.3f} between consecutive steps "
                f"(threshold {threshold:.3f}). Operator splitting may "
                f"be losing information; consider reducing dt or "
                f"switching to SimultaneousAdaptiveSolver if currently using "
                f"SimultaneousEulerSolver.",
            )

    def check_newton_iters(self, iters: Optional[int]) -> None:
        """Flag when the underlying Newton-style solver took too
        many iterations.

        High iteration counts indicate warm-start drift or a sharp
        composition change. No-op when the solver doesn't expose
        an iteration count (``None`` or negative — by design;
        scipy ``brentq`` doesn't surface it cleanly).
        """
        if iters is None:
            return
        try:
            n = int(iters)
        except (TypeError, ValueError):
            return
        if n < 0:
            return
        threshold = _get_warning_config().newton_iters_threshold
        if n > threshold:
            self._emit(
                "newton_iters",
                f"Speciation solver took {n} iterations "
                f"(threshold {threshold}). Warm-start may be far from "
                f"the solution; expect sharper composition change or "
                f"stiffness in this region.",
            )

    def check_charge_residual(self, residual: Optional[float]) -> None:
        """Flag a charge-balance residual that exceeds tolerance.

        Measures direct DAE-constraint drift after a kinetic step,
        before the end-of-step re-solve. No-op on ``None`` / NaN.
        """
        if residual is None:
            return
        try:
            r = float(residual)
        except (TypeError, ValueError):
            return
        if math.isnan(r):
            return
        threshold = _get_warning_config().charge_residual_threshold
        if abs(r) > threshold:
            self._emit(
                "charge_residual",
                f"Charge-balance residual {abs(r):.3e} mol/L exceeds "
                f"threshold {threshold:.3e}. The kinetic step drifted "
                f"off the algebraic constraint; the end-of-step "
                f"re-solve will correct it but the splitting error "
                f"may be larger than expected.",
            )

    def check_dt_vs_tau_min(
        self,
        dt_h: float,
        tau_min_h: Optional[float],
    ) -> None:
        """Flag when the integration step is large relative to the
        fastest kinetic timescale.

        Computed once at CV construction from a coarse
        ``1 / max_rate`` Jacobian estimate. No-op when the
        kinetic envelope can't be estimated (``tau_min_h=None``)
        — typical for property-solver-only CVs.
        """
        if tau_min_h is None or tau_min_h <= 0.0:
            return
        try:
            dt = float(dt_h)
            tau = float(tau_min_h)
        except (TypeError, ValueError):
            return
        if dt <= 0.0 or math.isnan(dt) or math.isnan(tau):
            return
        ratio = dt / tau
        threshold = _get_warning_config().dt_over_tau_min_threshold
        if ratio > threshold:
            self._emit(
                "dt_over_tau_min",
                f"Integration step dt={dt:.3g} h is {ratio:.2f}x the "
                f"fastest kinetic timescale tau_min={tau:.3g} h "
                f"(threshold ratio {threshold}). Consider reducing dt "
                f"or switching to a stiff-aware solver "
                f"(SimultaneousAdaptiveSolver method=Radau or BDF).",
            )

    def check_scipy_rejections(
        self,
        n_accepted: int,
        n_attempted: int,
    ) -> None:
        """Flag a high scipy ``solve_ivp`` step-rejection rate.

        Uses the ratio ``(n_attempted - n_accepted) / n_accepted``
        as a proxy for adaptive-step difficulty. ``n_attempted``
        is approximated by ``solve_ivp.nfev`` upstream; this is a
        rough indicator that spikes when the integrator is fighting
        stiffness or a tolerance boundary. No-op when the run was
        trivially short (``n_attempted <= 1``).
        """
        try:
            accepted = int(n_accepted)
            attempted = int(n_attempted)
        except (TypeError, ValueError):
            return
        if attempted <= 1 or accepted <= 0:
            return
        rejection_ratio = (attempted - accepted) / accepted
        threshold = _get_warning_config().scipy_rejection_threshold
        if rejection_ratio > threshold:
            self._emit(
                "scipy_rejections",
                f"SimultaneousAdaptiveSolver attempted {attempted} function "
                f"evaluations for {accepted} accepted steps "
                f"(ratio {rejection_ratio:.2f} above threshold "
                f"{threshold:.2f}). Suggests stiffness or a tight "
                f"tolerance boundary; consider switching to Radau/BDF "
                f"or relaxing rtol/atol.",
            )

    def check_ionic_strength(
        self,
        I_molL: Optional[float],
        activity_model: str,
        use_activity: bool,
    ) -> None:
        """Flag when ionic strength exceeds the regime of the
        configured activity model.

        Replaces the per-instance ``_warned_high_I`` dedup that
        used to live on ``BisectionChemicalEquilibriumEngine``; the throttle takes
        over.

        - ``use_activity=False`` is treated as the ideal-solution
          assumption (``ionic_strength_ideal_threshold``,
          default 0.10 mol/L).
        - ``davies`` activity model uses
          ``ionic_strength_davies_threshold`` (default 0.50 mol/L).
        - Unknown activity models silently no-op.
        """
        if I_molL is None:
            return
        try:
            I = float(I_molL)
        except (TypeError, ValueError):
            return
        if math.isnan(I):
            return
        cfg = _get_warning_config()
        if not use_activity:
            model_key = "ideal"
            threshold = cfg.ionic_strength_ideal_threshold
        else:
            model_key = (activity_model or "").strip().lower()
            if model_key == "davies":
                threshold = cfg.ionic_strength_davies_threshold
            else:
                return
        if I <= threshold:
            return
        if model_key == "ideal":
            msg = (
                f"Ionic strength {I:.2f} mol/L exceeds the threshold "
                f"~{threshold:.2f} mol/L where the ideal-solution "
                f"assumption is typically reasonable. Results may be "
                f"biased; consider enabling an activity model (e.g. "
                f"Davies)."
            )
        else:
            msg = (
                f"Ionic strength {I:.2f} mol/L exceeds the typical "
                f"recommended range for the {model_key.capitalize()} "
                f"activity model (threshold ~{threshold:.2f} mol/L). "
                f"Predictions may be unreliable at this ionic strength."
            )
        self._emit(f"ionic_strength_{model_key}", msg)

    def check_negative_mole(self, n_mol: "dict[str, float]") -> None:
        """Flag any species below zero beyond a small tolerance.

        Only meaningful when ``clamp_fn=None`` was configured on the
        active discrete-step ``StepSolver`` — with clamping enabled
        this should never fire. Disabling clamping is a real way to
        silently feed unphysical (negative) state into downstream
        chemistry with nothing watching for it; this check is that
        watcher.
        """
        tol = _get_warning_config().negative_mole_tolerance
        offenders = {sp: n for sp, n in n_mol.items() if n < -tol}
        if not offenders:
            return
        worst_sp = min(offenders, key=lambda sp: offenders[sp])
        self._emit(
            "negative_mole",
            f"Species {worst_sp!r} went negative ({offenders[worst_sp]:.3e} "
            f"mol, tolerance {tol:.1e}) after a step with clamping "
            f"disabled (clamp_fn=None). {len(offenders)} species "
            f"affected: {sorted(offenders)}.",
        )

    def check_clamp_invoked(
        self,
        deltas_before: "dict[str, float]",
        deltas_after: "dict[str, float]",
        dt_h: float,
    ) -> None:
        """Flag species whose delta ``clamp_fn`` changed by more than
        tolerance.

        Useful beyond visibility — it signals that the chosen
        ``dt_h``/kinetics combination is aggressive enough that
        clamping is doing real work, an early hint that a smaller
        ``dt_h`` might be needed rather than relying on the
        correction.
        """
        tol = _get_warning_config().clamp_invoked_tolerance
        changed = {
            sp: (before, deltas_after.get(sp, before))
            for sp, before in deltas_before.items()
            if abs(deltas_after.get(sp, before) - before) > tol
        }
        if not changed:
            return
        parts = ", ".join(
            f"{sp}: {before:.3g}->{after:.3g} mol/h"
            for sp, (before, after) in sorted(changed.items())
        )
        self._emit(
            "clamp_invoked",
            f"clamp_fn scaled removal rate(s) to prevent negative "
            f"inventory over dt_h={dt_h:.3g}: {parts}. Consider a "
            f"smaller dt_h if this fires often, rather than relying on "
            f"the correction.",
        )
