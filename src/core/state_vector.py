# -*- coding: utf-8 -*-
"""Unified state-vector packing/unpacking for ODE-style solvers.

Replaces the duplicated logic behind three near-identical
implementations that had started to drift apart:
``_pack_state``/``_unpack_state``/``_state_index_map``
(``system_solver.py``, multi-CV, no exclusion set),
``_pack_extended_state``/``_unpack_extended_state``
(``system_solver.py``, multi-CV + controller differential-state
extension, Phase E), and ``_StateVector`` (``solvers.py``, single-CV,
phase-scoped, with an ``algebraic_species`` exclusion set — used by
``SimultaneousAdaptiveSolver``, now deleted). The first two pairs
became thin wrapper functions delegating to :class:`StateVector` below
(checkpoint 7), keeping their original signatures so no existing
caller (``ImplicitTransportSystemSolver``, ``MonolithicODESolver``,
``_build_rhs``) needed to change. ``_StateVector`` was migrated
separately (checkpoint 7b), alongside generalizing
``SimultaneousAdaptiveSolver`` off the gas/liquid assumption — its
species lists were read as raw attributes (``sv.gas_species``,
``sv.n_gas``, …) throughout that solver's ODE hot loop and Jacobian
machinery, not through clean function calls, so unifying it could not
be separated from that generalization. See
docs/phases-upcoming/STEP_SOLVER_REFINEMENT_CHECKLIST.md.

Species-ordering convention throughout: CV (``cvs`` dict insertion
order) → phase (alphabetical) → species (alphabetical).
"""

from __future__ import annotations

from typing import Any, Dict, FrozenSet, List, Optional, Tuple

import numpy as np


class StateVector:
    """Packs/unpacks CV phase state (+ optional controller differential
    state) into a flat numpy array.

    The species layout is computed once at construction time from the
    CVs' current ``phase.n_mol`` keys — callers that run many
    pack/unpack cycles over an unchanging species set (e.g. one
    ``solve_ivp`` integration) should construct one instance and reuse
    it, rather than reconstructing per call.

    Parameters
    ----------
    cvs : dict[str, ControlVolume]
        CVs in scope, keyed by ``cv_key``, in the order they should
        appear in the packed vector (dict insertion order).
    exclude_species : frozenset[str], optional
        Species keys excluded from every phase's packing (e.g.
        ``{"H+"}`` for an algebraically-determined species). Default:
        no exclusion.
    ctrl_list : list of (controller, sorted_keys), optional
        Controller differential-state extension (Phase E), appended
        after the CV species block in ``sim.controllers`` order
        (already filtered/sorted by the caller, e.g. via
        ``_ctrl_list(sim)``). Default: no controller extension.
    """

    def __init__(
        self,
        cvs: Dict[str, Any],
        *,
        exclude_species: FrozenSet[str] = frozenset(),
        ctrl_list: Optional[List[tuple]] = None,
    ) -> None:
        self.cvs = cvs
        self.exclude_species = exclude_species
        self.ctrl_list: List[tuple] = list(ctrl_list) if ctrl_list else []

        layout: List[Tuple[str, str, str]] = []
        for cv_key, cv in cvs.items():
            for pk in sorted(cv.phases):
                n_mol = cv.phases[pk].n_mol
                for sk in sorted(n_mol):
                    if sk in self.exclude_species:
                        continue
                    layout.append((cv_key, pk, sk))
        self._layout: List[Tuple[str, str, str]] = layout

        self.n_cv: int = len(layout)
        n_ctrl = sum(len(keys) for _, keys in self.ctrl_list)
        self.n_total: int = self.n_cv + n_ctrl

    def pack(self) -> np.ndarray:
        """Pack current CV (+ controller) state into a flat array."""
        y = np.zeros(self.n_total, dtype=np.float64)
        for i, (cv_key, pk, sk) in enumerate(self._layout):
            y[i] = float(self.cvs[cv_key].phases[pk].n_mol.get(sk, 0.0))
        offset = self.n_cv
        for ctrl, keys in self.ctrl_list:
            state = ctrl.differential_state()
            for k in keys:
                y[offset] = float(state[k])
                offset += 1
        return y

    def unpack(self, y: np.ndarray, *, floor: bool = False) -> None:
        """Write a packed array back into CV phases (+ controllers), in place.

        Inverse of :meth:`pack`. ``y`` must have the same length and
        ordering as an array produced by :meth:`pack` on this same
        instance (i.e. the species set must not have changed since
        construction).

        Parameters
        ----------
        floor : bool
            If True, floor each unpacked value at 0.0 before writing
            (guards against adaptive-integrator overshoot on the final
            accepted solution). Default False, matching
            ``_unpack_state``'s original behavior — callers that need
            the floor (e.g. writing a final ``solve_ivp`` result back
            to the CV) opt in explicitly.
        """
        if floor:
            from .clamping import floor_nonnegative
            y = floor_nonnegative(y)
        for i, (cv_key, pk, sk) in enumerate(self._layout):
            self.cvs[cv_key].phases[pk].n_mol[sk] = float(y[i])
        offset = self.n_cv
        for ctrl, keys in self.ctrl_list:
            ctrl.set_state({k: float(y[offset + i]) for i, k in enumerate(keys)})
            offset += len(keys)

    @property
    def index_map(self) -> Dict[Tuple[str, str, str], int]:
        """``{(cv_key, phase_key, species_key): index}`` for this layout."""
        return {key: i for i, key in enumerate(self._layout)}

    def phase_species(self, cv_key: str) -> Dict[str, List[str]]:
        """Return ``{phase_key: [species_key, ...]}`` for one CV in this
        instance's scope, in the same order used by :meth:`pack`/
        :meth:`unpack`.

        For hot loops that need direct per-phase index arithmetic
        instead of going through the full ``(cv_key, phase_key,
        species_key)`` tuple layout on every call (e.g.
        :class:`~PyOMES.core.solvers.SimultaneousAdaptiveSolver`'s ODE
        derivative closure) — combine with :meth:`phase_offset`.
        """
        result: Dict[str, List[str]] = {}
        for k, pk, sk in self._layout:
            if k == cv_key:
                result.setdefault(pk, []).append(sk)
        return result

    def phase_offset(self, cv_key: str) -> Dict[str, int]:
        """Return ``{phase_key: starting_index}`` for one CV — the flat
        index each phase's species block begins at, matching
        :meth:`phase_species`'s ordering."""
        offsets: Dict[str, int] = {}
        offset = 0
        for pk, species in self.phase_species(cv_key).items():
            offsets[pk] = offset
            offset += len(species)
        return offsets
