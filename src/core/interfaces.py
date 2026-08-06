# -*- coding: utf-8 -*-
"""Phase interface protocol for inter-phase and inter-CV transfer.

A :class:`PhaseInterface` sits between two phases (which may be in the
same ControlVolume or in different ones) and computes molar fluxes based
on the current state of both sides.

Examples of internal interfaces (within a single CV):
  - Gas-liquid: Henry equilibrium or kLa kinetic transfer
  - Liquid-solid: dissolution, precipitation, adsorption
  - Gas-solid: heterogeneous reaction (e.g. direct biogas–solid contact)

Examples of external interfaces (between a CV and the environment):
  - Gas feed: fixed-rate or pressure-driven gas supply
  - Vent: pressure-relief valve releasing gas to atmosphere
  - Liquid feed: substrate addition
  - Liquid drain: product removal

The distinction between internal and external is semantic, not structural —
both satisfy the same :class:`PhaseInterface` protocol.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Protocol, runtime_checkable

from .phases import Phase
from .links import LinkFlowRecord
from .boundaries import ExternalFluxRecord


# ════════════════════════════════════════════════════════════════════════
#  Protocol
# ════════════════════════════════════════════════════════════════════════

@runtime_checkable
class PhaseInterface(Protocol):
    """Contract for any interface between two phases.

    Implementations must define:
      - ``phase_a_key``, ``phase_b_key``: string keys identifying the
        phases in the owning ControlVolume's ``phases`` dict.
      - ``compute_flux(state_a, state_b, dt_h)``: returns species molar
        fluxes (mol/h) **from phase_a to phase_b**.  Positive = transfer
        from a→b.  Negative = transfer from b→a.

    Convention
    ----------
    The flux dict is applied as:
      - ``phase_a.apply_flux(-flux, dt)``   (loses material)
      - ``phase_b.apply_flux(+flux, dt)``   (gains material)

    This ensures mass conservation: what leaves one phase enters the other.
    """

    @property
    def phase_a_key(self) -> str: ...

    @property
    def phase_b_key(self) -> str: ...

    def compute_flux(
        self,
        state_a: Phase,
        state_b: Phase,
        dt_h: float,
        *,
        instantaneous: bool = False,
    ) -> Dict[str, float]:
        """Compute species molar fluxes from phase_a to phase_b.

        Parameters
        ----------
        state_a, state_b : Phase
            Current state of each phase. Kinetic interfaces that
            need pH or speciation read them inline from
            ``state_b.n_mol`` (canonical species populated by the
            speciation engine writeback during the CV's
            ``advance()`` step).
        dt_h : float
            Timestep duration (hours).  Needed by rate-based (kinetic)
            interfaces.  Equilibrium interfaces may ignore it.
        instantaneous : bool
            If True, return the instantaneous rate at the current
            state rather than a step-averaged rate — needed by
            continuous-RHS solvers (e.g.
            :class:`~PyOMES.core.solvers.SimultaneousAdaptiveSolver`),
            which need ``dn/dt``, not a discrete-step-averaged flux.
            Mirrors :meth:`~PyOMES.core.boundaries.ExternalBoundary.compute_flux`'s
            identically-named parameter. Default False (discrete-step
            solvers' usage).

        Returns
        -------
        dict
            ``{species_id: flux_mol_per_h}`` — positive means a→b.
        """
        ...


# ════════════════════════════════════════════════════════════════════════
#  Transfer result (returned by ControlVolume.step_internal_transfer)
# ════════════════════════════════════════════════════════════════════════

@dataclass
class TransferDiagnostics:
    """Mass-balance diagnostics for a single ``step_internal_transfer`` call.

    Attributes
    ----------
    species_before : dict
        Total moles of each species across all phases *before* the step.
    species_after : dict
        Total moles of each species across all phases *after* the step.
    residual : dict
        Per-species residual: ``after - before``.  Should be ~0 for all
        species that are only transferred (not created or destroyed).
    fluxes : list of dict
        The flux dict returned by each internal interface, in order.
    total_mol_before : float
        Sum of all species moles before.
    total_mol_after : float
        Sum of all species moles after.
    """
    species_before: Dict[str, float] = field(default_factory=dict)
    species_after: Dict[str, float] = field(default_factory=dict)
    residual: Dict[str, float] = field(default_factory=dict)
    fluxes: list = field(default_factory=list)
    total_mol_before: float = 0.0
    total_mol_after: float = 0.0

    @property
    def total_residual(self) -> float:
        """Absolute total residual (mol).  Should be ~0."""
        return self.total_mol_after - self.total_mol_before

    @property
    def max_abs_residual(self) -> float:
        """Largest absolute residual across all species (mol)."""
        if not self.residual:
            return 0.0
        return max(abs(v) for v in self.residual.values())

    def is_conserved(self, tol_mol: float = 1e-12) -> bool:
        """Check whether all species are conserved within tolerance."""
        return self.max_abs_residual <= tol_mol

    def summary(self) -> str:
        """Human-readable summary string."""
        lines = ["Transfer Diagnostics:"]
        lines.append(f"  Total mol before: {self.total_mol_before:.8e}")
        lines.append(f"  Total mol after:  {self.total_mol_after:.8e}")
        lines.append(f"  Total residual:   {self.total_residual:.4e}")
        lines.append(f"  Max |residual|:   {self.max_abs_residual:.4e}")
        lines.append(f"  Conserved:        {self.is_conserved()}")
        if self.residual:
            lines.append("  Per-species residual:")
            for sp, r in sorted(self.residual.items(), key=lambda x: -abs(x[1])):
                if abs(r) > 0.0:
                    lines.append(f"    {sp:20s}  {r:+.6e} mol")
        return "\n".join(lines)


# ════════════════════════════════════════════════════════════════════════
#  Advance result (returned by ControlVolume.advance)
# ════════════════════════════════════════════════════════════════════════

@dataclass
class AdvanceResult:
    """Result of a single :meth:`ControlVolume.advance` call.

    Aggregates the outputs of each sub-step (reactions, external source
    terms, internal equilibrium transfer) so the orchestrator can
    inspect what happened without reaching into the CV.

    Attributes
    ----------
    transfer : TransferDiagnostics or None
        Mass-balance diagnostics from ``step_internal_transfer``.
        ``None`` if the CV has no internal interfaces.  This is the
        per-step mass-balance audit used by the sequential-body path;
        it is **not** the same as ``transfer_record`` below.
    reaction_sources : dict or None
        Net source terms applied from the reaction model, as
        ``{phase_key: {species_id: mol_per_h}}``.  ``None`` if the CV
        has no reaction model.
    transfer_record : LinkFlowRecord or None
        Diagnostics from the gas-liquid (or other kinetic) transfer
        link populated by snapshot-based solvers
        (e.g. :class:`SimultaneousEulerSolver`).  Despite the similar name,
        this differs from ``transfer`` above:
        ``transfer`` is a :class:`TransferDiagnostics` (mass-balance
        residuals across all internal interfaces, used by the
        sequential body), while ``transfer_record`` is a
        :class:`LinkFlowRecord` (per-link species flow totals, used
        by the snapshot solvers).  Default ``None``.
    boundary_records : list of ExternalFluxRecord
        Per-boundary application diagnostics populated by
        snapshot-based solvers.  The default sequential body does not
        emit per-boundary records (it applies boundary fluxes through
        :meth:`Phase.apply_flux` and leaves diagnostics empty).

    Notes
    -----
    The ``properties`` field was removed in state-unification C4. pH
    and other derived chemistry-state lives on the CV's phases
    directly: read ``cv.phases["liquid"].pH`` (raises if no
    speciation engine populated ``n_mol["H+"]``) and
    ``cv.phases["liquid"].n_mol[species_id]`` for individual species
    concentrations.
    """

    transfer: Optional[TransferDiagnostics] = None
    reaction_sources: Optional[Dict[str, Dict[str, float]]] = None
    transfer_record: Optional[LinkFlowRecord] = None
    boundary_records: List[ExternalFluxRecord] = field(default_factory=list)
