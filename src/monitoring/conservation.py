# -*- coding: utf-8 -*-
"""``ConservationMonitor`` and ``ConservationWarning``.

Introduced in state-unification C6. Tracks element and charge
balances per :meth:`ControlVolume.advance` step and emits an
:class:`ConservationWarning` when drift exceeds the threshold
configured on ``PyOMES.config.warnings``. Mirrors the throttle /
summary plumbing of :class:`~PyOMES.monitoring.accuracy.AccuracyMonitor`.

Three failure modes the monitor catches:

- **Feeds without matching counterions.** A boundary that adds
  ``Na+`` to the liquid without adding ``Cl-`` (or vice-versa)
  breaks charge balance; the engine's writeback will rebalance pH
  to compensate but the underlying state becomes physically wrong.
- **Kinetic rates that violate stoichiometric balance.** A
  ``KineticReaction`` whose ``compute_rates`` returns source terms
  that don't elementally close. The construction-time
  :func:`validate_balance` check guards against bad stoichiometry
  for declared reactions, but black-box reaction models and rate
  laws that overload species can still produce non-conservative
  source terms.
- **Accumulated float roundoff.** Long simulations can drift below
  the threshold per step but cross the cumulative threshold over
  many steps — a calibration mismatch in the integrator's
  tolerances or a subtle conservation bug.

Defaults — relative-plus-max-clamp:

- **Per-step**: ``1e-8 × max(total_element_mol, 1)`` /
  ``1e-8 × max(|sum z·n|_initial, 1)`` for charge.
- **Cumulative**: ``1e-6 × max(total_element_mol, 1)``.

The ``max(..., 1)`` clamp prevents dilute systems from tripping on
noise-level absolute residuals (which would otherwise fire when
``total_element_mol < 1`` and the relative threshold becomes
sub-floating-point).
"""

from __future__ import annotations

import math
import warnings
from collections import Counter
from typing import Dict, List, Optional, Set

from ..chemistry.species import Species


class ConservationWarning(UserWarning):
    """Emitted when element or charge balance drifts beyond the
    configured threshold for a CV's advance() step.

    Filter or promote via standard ``warnings`` machinery::

        import warnings
        from PyOMES import ConservationWarning

        warnings.simplefilter("ignore", ConservationWarning)  # silence
        warnings.simplefilter("error",  ConservationWarning)  # promote
    """


# Module-level emission counter (parallel to AccuracyMonitor's).
_summary_counter: Counter = Counter()


def _get_warning_config():
    """Read the current ``WarningConfig`` from the package singleton."""
    from PyOMES.config import config as _config_singleton

    return _config_singleton.warnings


def print_conservation_summary() -> None:
    """Emit the end-of-run summary line for the ``"once"`` throttle."""
    if not _summary_counter:
        return
    parts = [f"{cat}={count}" for cat, count in sorted(_summary_counter.items())]
    total = sum(_summary_counter.values())
    print(
        f"ConservationWarning summary: {total} emission(s) across "
        f"[{', '.join(parts)}]"
    )


def reset_conservation_summary() -> None:
    """Clear the module-level summary counter."""
    _summary_counter.clear()


class ConservationMonitor:
    """Per-CV conservation check runner.

    Tracks element totals and charge balance across phases. On
    each ``check_step`` call (driven by ``ControlVolume.advance``),
    compares current totals to baseline / previous step and emits
    a :class:`ConservationWarning` when drift exceeds the
    configured threshold.

    Attach via ``cv.reaction_system.attach_conservation_monitor(monitor)``
    (mirrors :meth:`attach_monitor` for the AccuracyMonitor).

    Species composition comes from the
    :attr:`_species_registry` dict, populated automatically by
    :meth:`ControlVolume.__init__` by walking the reaction
    stoichiometries. Species in ``phase.n_mol`` without a registry
    entry are skipped — typical for unnamed strong-ion lumps like
    ``S_cat``/``S_an`` (``atoms={}``, only charge contributes;
    those are looked up from the engine's strong-ion mapping in
    a separate code path if needed).
    """

    def __init__(self) -> None:
        self.step_count: int = 0
        self._baseline_element_totals: Dict[str, float] = {}
        self._prior_element_totals: Dict[str, float] = {}
        self._baseline_set: bool = False
        self._species_registry: Dict[str, Species] = {}
        self._warned_categories: Set[str] = set()

    def reset(self) -> None:
        """Clear per-monitor state (baseline, history, warning dedup).

        Does NOT touch the module-level summary counter — call
        :func:`reset_conservation_summary` for that.
        """
        self.step_count = 0
        self._baseline_element_totals.clear()
        self._prior_element_totals.clear()
        self._baseline_set = False
        self._warned_categories.clear()

    def set_species_registry(
        self, registry: Dict[str, Species]
    ) -> None:
        """Provide / update the ``species_id → Species`` map used
        for element + charge accounting. Called by
        :meth:`ControlVolume.__init__` after reactions are attached.
        """
        self._species_registry = dict(registry)

    # ── Internals ────────────────────────────────────────────────

    def _emit(self, category: str, message: str) -> None:
        cfg = _get_warning_config()
        throttle = cfg.throttle
        if throttle == "silent":
            return
        if throttle == "once" and category in self._warned_categories:
            return
        if throttle == "first_N":
            if _summary_counter[category] >= cfg.first_N:
                return
        self._warned_categories.add(category)
        _summary_counter[category] += 1
        warnings.warn(message, ConservationWarning, stacklevel=3)

    def _element_totals(self, phases) -> Dict[str, float]:
        """Sum element moles across all phases.

        For each species in each phase's n_mol, multiply mole count
        by the atom count for each element on the species, and
        accumulate. Returns ``{element: total_mol}``.
        """
        totals: Dict[str, float] = {}
        for phase in phases.values():
            n_mol = getattr(phase, "n_mol", {})
            for sp_id, n in n_mol.items():
                species = self._species_registry.get(sp_id)
                if species is None:
                    continue
                for elem, count in species.atoms.items():
                    if count == 0:
                        continue
                    totals[elem] = totals.get(elem, 0.0) + float(n) * float(count)
        return totals

    def _charge_residual(self, phases) -> float:
        """Sum ``z × n`` across all species in all phases.

        Includes species whose ``Species`` record carries a
        non-zero charge (including the unnamed-lump strong-ion
        species ``S_cat`` / ``S_an`` if they are in the registry).
        """
        residual = 0.0
        for phase in phases.values():
            n_mol = getattr(phase, "n_mol", {})
            for sp_id, n in n_mol.items():
                species = self._species_registry.get(sp_id)
                if species is None:
                    continue
                charge = int(getattr(species, "charge", 0))
                if charge == 0:
                    continue
                residual += charge * float(n)
        return residual

    # ── Per-step check ───────────────────────────────────────────

    def check_step(self, phases) -> None:
        """Run the per-step element + charge conservation checks.

        Called once per :meth:`ControlVolume.advance` step at the
        end of the sequential body. Establishes the baseline on
        first call.

        Element drift:

        - Per-step: warn if ``|current - prior| > ratio_step ×
          max(prior, 1)`` for any element.
        - Cumulative: warn if ``|current - baseline| > ratio_cum ×
          max(baseline, 1)`` for any element.

        Charge drift:

        - Per-step: warn if ``|current - prior| > ratio_step ×
          max(|baseline|, 1)``.
        - Cumulative: warn if ``|current| > ratio_cum ×
          max(|baseline|, 1)``.

        No-op if the species registry is empty (the CV has no
        reactions / no Species records to drive accounting).
        """
        self.step_count += 1
        if not self._species_registry:
            return

        cfg = _get_warning_config()
        ratio_step = float(cfg.conservation_per_step_factor)
        ratio_cum = float(cfg.conservation_cumulative_factor)

        elem_totals = self._element_totals(phases)
        charge_total = self._charge_residual(phases)

        if not self._baseline_set:
            self._baseline_element_totals = dict(elem_totals)
            self._prior_element_totals = dict(elem_totals)
            self._baseline_charge = charge_total
            self._baseline_set = True
            return

        # Element checks
        for elem, current in elem_totals.items():
            prior = self._prior_element_totals.get(elem, 0.0)
            baseline = self._baseline_element_totals.get(elem, 0.0)
            scale = max(abs(baseline), 1.0)
            delta_step = abs(current - prior)
            delta_cum = abs(current - baseline)
            if delta_step > ratio_step * scale:
                self._emit(
                    f"element_step_{elem}",
                    f"Element '{elem}' total drifted by "
                    f"{delta_step:.3e} mol between consecutive steps "
                    f"(threshold {ratio_step * scale:.3e}, "
                    f"baseline {baseline:.3e}). Likely cause: a "
                    f"kinetic reaction or boundary that does not "
                    f"close on '{elem}'.",
                )
            if delta_cum > ratio_cum * scale:
                self._emit(
                    f"element_cum_{elem}",
                    f"Element '{elem}' cumulative drift "
                    f"{delta_cum:.3e} mol since simulation start "
                    f"(threshold {ratio_cum * scale:.3e}, "
                    f"baseline {baseline:.3e}). Likely cause: "
                    f"accumulated float roundoff or a slow "
                    f"persistent {elem} leak.",
                )

        # Charge checks
        baseline_charge = getattr(self, "_baseline_charge", 0.0)
        prior_charge = getattr(self, "_prior_charge", baseline_charge)
        charge_scale = max(abs(baseline_charge), 1.0)
        delta_charge_step = abs(charge_total - prior_charge)
        if delta_charge_step > ratio_step * charge_scale:
            self._emit(
                "charge_step",
                f"Charge balance drifted by "
                f"{delta_charge_step:.3e} mol between consecutive "
                f"steps (threshold {ratio_step * charge_scale:.3e}). "
                f"Likely cause: a feed without a matching counterion.",
            )
        if abs(charge_total) > ratio_cum * charge_scale:
            self._emit(
                "charge_cum",
                f"Charge residual {abs(charge_total):.3e} mol "
                f"exceeds cumulative threshold "
                f"{ratio_cum * charge_scale:.3e}. Likely cause: "
                f"accumulated charge imbalance from unmatched "
                f"feeds or stoichiometric breakage.",
            )

        self._prior_element_totals = dict(elem_totals)
        self._prior_charge = charge_total
