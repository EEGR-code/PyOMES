# -*- coding: utf-8 -*-
"""Speciation engine protocol hierarchy.

Defines the capability tiers that speciation engines may satisfy, plus the
data containers they return.  All protocols are @runtime_checkable so that
callers can use ``isinstance(engine, GrayBoxEngineProtocol)`` without
importing concrete engine types.

Tier summary
------------
ChemicalEquilibriumEngineProtocol   — black box: solve() only
GrayBoxEngineProtocol      — adds jacobian_dz_dy() (total sensitivity)
WhiteBoxEngineProtocol     — adds residual(), jacobian_dg_dz(), jacobian_dg_dy()

See docs/dev/ideas/CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md for the full design.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, FrozenSet, Optional, Tuple

import numpy as np
from scipy.sparse import spmatrix
from typing import Protocol, runtime_checkable


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SpeciationJacobian:
    """Dense total-sensitivity matrix: dz_dy[i, j] = ∂c_i / ∂T_j.

    Rows are indexed by algebraic_ids (species concentrations);
    columns are indexed by component_ids (component totals).
    """
    dz_dy: np.ndarray
    algebraic_ids: Tuple[str, ...]
    component_ids: Tuple[str, ...]


@dataclass(frozen=True)
class SparseJacobian:
    """Sparse partial Jacobian returned by split white-box methods.

    Wraps a scipy.sparse CSR matrix with named row and column indices
    for alignment with the solver state vector.
    """
    matrix: spmatrix
    row_ids: Tuple[str, ...]
    col_ids: Tuple[str, ...]


# ---------------------------------------------------------------------------
# EquilibriumResult — immutable solve() output.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EquilibriumResult:
    """Immutable snapshot of a converged equilibrium solve.

    Returned by ``solve()`` on :class:`BisectionChemicalEquilibriumEngine`,
    :class:`NRChemicalEquilibriumEngine`, and :class:`PHREEQCChemicalEquilibriumEngine`. ``solve()``
    makes no side-effecting writes to phases; state is committed only by an
    explicit call to :meth:`apply_to_phases`.

    Field provenance
    ----------------
    ``pH``, ``pH_conc``, ``logH``, ``aH``, ``gamma_H``, ``gamma_OH``,
    ``ionic_strength``, ``charge_residual`` are the common meta keys emitted
    by both the legacy acid-base solver (``acid_base.py``) and the NR solver
    (``nr_solver.py``). ``n_iter`` is carried for forward-compatibility with
    callers that already guard on it (e.g. ``speciation/api.py``); it is not
    currently populated by any engine.

    ``species_mol_L`` holds equilibrated species concentrations (mol/L),
    keyed by species ID — this is what :meth:`apply_to_phases` writes back.
    Its exact membership is engine-specific: :class:`BisectionChemicalEquilibriumEngine`
    restricts it to the canonical species tuple it has always written back;
    :class:`NRChemicalEquilibriumEngine` includes every tableau master + secondary
    (+ H2O); :class:`PHREEQCChemicalEquilibriumEngine` includes every PHREEQC species.

    ``saturation_indices`` holds mineral SI values (from precipitation
    equilibria, when declared). ``extra`` is a catch-all for engine-specific
    diagnostics that don't warrant a named field (e.g. mineral
    ``xi_mol_L``, or generic polyprotic-ladder keys that were never part of
    the phase writeback).
    """

    pH: float
    pH_conc: Optional[float] = None
    logH: Optional[float] = None
    aH: Optional[float] = None
    gamma_H: Optional[float] = None
    gamma_OH: Optional[float] = None
    ionic_strength: Optional[float] = None
    charge_residual: Optional[float] = None
    n_iter: Optional[int] = None
    species_mol_L: Dict[str, float] = field(default_factory=dict)
    partial_pressures_atm: Dict[str, float] = field(default_factory=dict)
    saturation_indices: Dict[str, float] = field(default_factory=dict)
    alphas: Dict[str, float] = field(default_factory=dict)
    extra: Dict[str, Any] = field(default_factory=dict)

    def apply_to_phases(self, phases: Any, liquid_key: str = "liquid") -> None:
        """Commit this result to phase state.

        Called explicitly by the owning solver/CV when a step is accepted —
        not on every speciation solve, so that rejected trial steps inside
        adaptive/implicit ODE solvers do not corrupt phase state.

        Writes ``species_mol_L`` (mol/L) to ``phases[liquid_key].n_mol``
        (mol) via ``_refresh_derived``, scaled by ``phases[liquid_key].V_L``.
        No-op if the liquid phase is absent, has no volume, or doesn't
        support derived-species writeback.

        .. note::
            Does **not** write ``partial_pressures_atm`` back to a gas
            phase — CP2 of ``LAYER1_GAP_CLOSURE`` (which first populates
            that field, for folded gas-liquid secondaries) scoped gas-
            phase write-back out of this phase entirely. CP3 added the
            related ``step_internal_transfer()`` engine-owned-species
            scope-filter fix but deliberately did not add gas write-back
            alongside it (a separate, still-open follow-up). Until it's
            implemented, a folded ``NRChemicalEquilibriumEngine`` solve is
            a correct diagnostic (its returned split is right) but has no
            effect on gas-phase ``n_mol`` when committed via this method.
        """
        liq = phases.get(liquid_key) if hasattr(phases, "get") else None
        if liq is None or not hasattr(liq, "_refresh_derived"):
            return
        V_L = float(getattr(liq, "V_L", 0.0))
        if V_L <= 0.0:
            return
        derived_mol = {sp: conc * V_L for sp, conc in self.species_mol_L.items()}
        if derived_mol:
            liq._refresh_derived(derived_mol)

    def to_dict(self) -> Dict[str, Any]:
        """Return a plain-dict snapshot of every field.

        A narrow, explicit concession for callers that need a dict-shaped
        result (e.g. :class:`~PyOMES.chemical_equilibrium.api.SpeciationEngineAdapter`'s
        ``raw=`` field) — this is not a ``Mapping`` implementation:
        ``EquilibriumResult`` does not support ``result["pH"]`` or
        ``result.get(...)``, only attribute access and this explicit method.
        """
        return asdict(self)


# ---------------------------------------------------------------------------
# Capability mixin protocols
# ---------------------------------------------------------------------------

@runtime_checkable
class ResidualCapable(Protocol):
    """White-box: evaluate the equilibrium residual g(y, z) without re-solving."""

    def residual(self, **kwargs: Any) -> np.ndarray:
        ...


@runtime_checkable
class SolutionJacobianCapable(Protocol):
    """Gray-box: return the total sensitivity ∂z/∂y at the current operating point."""

    def jacobian_dz_dy(self, **kwargs: Any) -> SpeciationJacobian:
        ...


@runtime_checkable
class SplitJacobianCapable(Protocol):
    """White-box: return the split partial Jacobians ∂g/∂z and ∂g/∂y."""

    def jacobian_dg_dz(self, **kwargs: Any) -> SparseJacobian:
        ...

    def jacobian_dg_dy(self, **kwargs: Any) -> SparseJacobian:
        ...


# ---------------------------------------------------------------------------
# Named tier protocols
# ---------------------------------------------------------------------------

@runtime_checkable
class ChemicalEquilibriumEngineProtocol(Protocol):
    """Black-box speciation engine: solve() interface only.

    Both BisectionChemicalEquilibriumEngine and NRChemicalEquilibriumEngine satisfy this protocol
    structurally without any changes (verified 2026-06-25).
    """

    n_solve_calls: int

    def solve(self, **kwargs: Any) -> "EquilibriumResult":
        ...

    def algebraic_species(self) -> FrozenSet[str]:
        ...

    def reset_cache(self) -> None:
        ...

    def reset_counters(self) -> None:
        ...


@runtime_checkable
class GrayBoxEngineProtocol(ChemicalEquilibriumEngineProtocol, SolutionJacobianCapable, Protocol):
    """Gray-box engine: black-box interface plus total sensitivity jacobian_dz_dy()."""


@runtime_checkable
class WhiteBoxEngineProtocol(
    GrayBoxEngineProtocol,
    ResidualCapable,
    SplitJacobianCapable,
    Protocol,
):
    """White-box engine: full residual + split Jacobians for DAE coupling."""
