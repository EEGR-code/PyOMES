# -*- coding: utf-8 -*-
"""Axis 2 system-integrator extension point and concrete implementations.

Phase C (SYSTEM_SOLVER_PROTOCOL) delivers:

* :class:`SystemSolver` — protocol that ``Simulation._step()`` dispatches to.
* :class:`EventScheduler` — priority-queue manager for periodic controller
  update events (ZOH semantics).
* :func:`_pack_state` / :func:`_unpack_state` — state-vector packing
  convention shared with Phases D and E.
* :class:`ExplicitEulerSystemSolver` — transparent wrapper around the
  default ``Simulation._step_default()`` body; validates protocol wiring.
* :class:`StrangSplittingSystemSolver` — second-order Strang split of
  link-flux from per-CV chemistry.
* :class:`MultirateSystemSolver` — adaptive subcycling of link flux at
  CFL-safe sub-steps while advancing CVs once per macro step.

See ``docs/dev/implementation/shipped/SYSTEM_SOLVER_PROTOCOL.md`` for the design.
"""

from __future__ import annotations

import heapq
import math
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Protocol, runtime_checkable

import numpy as np

if TYPE_CHECKING:
    from .simulation import Simulation
    from .interfaces import AdvanceResult


# ════════════════════════════════════════════════════════════════════════
#  SystemSolver protocol
# ════════════════════════════════════════════════════════════════════════

@runtime_checkable
class SystemSolver(Protocol):
    """Axis 2 extension point: advances the full multi-CV system one step.

    ``advance_system`` is fully responsible for advancing all CVs, applying
    all links, firing controller events, and updating controller states for
    one macro timestep.  It returns the same 4-tuple as
    ``Simulation._step_default()``:

        (results, link_records, controller_actions, profile_actions)

    where ``results`` is ``Dict[str, AdvanceResult]``.  Implementors may
    delegate to ``sim._invoke_profiles``, ``sim._apply_links``,
    ``sim._invoke_controllers``, and ``sim._apply_controller_action``
    to reuse all existing Simulation bookkeeping.
    """

    def advance_system(
        self,
        sim: "Simulation",
        dt_h: float,
        t_h: float,
    ) -> Any: ...  # returns (results, link_records, ctrl_actions, prof_actions)


# ════════════════════════════════════════════════════════════════════════
#  EventScheduler
# ════════════════════════════════════════════════════════════════════════

class EventScheduler:
    """Priority-queue manager for periodic controller update events.

    Only controllers that expose a positive ``update_period_h`` (Phase B
    protocol) are placed in the queue.  Controllers without
    ``update_period_h`` (or with ``update_period_h = None``) fire every
    macro step via the existing ``Simulation._should_fire`` path.

    ZOH cache: ``action_cache[id(ctrl)]`` holds the last
    :class:`~PyOMES.control.actions.ControlAction` emitted by each
    periodic controller.  Between T_c events the cached action is
    re-applied unchanged by the system solver.

    Parameters
    ----------
    controllers : iterable
        The controller list from ``Simulation.controllers``.
    t_start_h : float
        Wall-clock time at run entry.  First events scheduled at
        ``t_start_h + update_period_h``.
    """

    def __init__(self, controllers, t_start_h: float) -> None:
        self._heap: List = []          # (next_fire_t_h, idx, ctrl)
        self._period: Dict[int, float] = {}  # id(ctrl) → period_h
        self._counter: int = 0         # monotone tie-breaker
        self.action_cache: Dict[int, Any] = {}  # id(ctrl) → ControlAction

        for ctrl in controllers:
            period = getattr(ctrl, "update_period_h", None)
            if period is not None:
                period = float(period)
                if period > 0.0:
                    next_t = t_start_h + period
                    heapq.heappush(
                        self._heap, (next_t, self._counter, ctrl)
                    )
                    self._counter += 1
                    self._period[id(ctrl)] = period

    def next_event_h(self) -> Optional[float]:
        """Smallest scheduled fire time, or ``None`` if queue is empty."""
        return self._heap[0][0] if self._heap else None

    def controllers_due(self, t_h: float, tol: float = 1e-12) -> List:
        """Pop and return every controller whose next fire time ≤ ``t_h``."""
        due = []
        while self._heap and self._heap[0][0] <= t_h + tol:
            _, _, ctrl = heapq.heappop(self._heap)
            due.append(ctrl)
        return due

    def advance_clocks(self, fired: List, t_h: float) -> None:
        """Re-schedule each fired controller to ``t_h + T_c``."""
        for ctrl in fired:
            period = self._period.get(id(ctrl))
            if period is not None:
                next_t = t_h + period
                heapq.heappush(
                    self._heap, (next_t, self._counter, ctrl)
                )
                self._counter += 1


# ════════════════════════════════════════════════════════════════════════
#  State-vector packing convention (shared with Phases D and E)
# ════════════════════════════════════════════════════════════════════════

def _pack_state(sim: "Simulation") -> np.ndarray:
    """Pack all CV species into a flat state vector.

    Ordering convention (deterministic across calls):

    * CVs in ``sim.cvs`` insertion order.
    * Phases within each CV in alphabetical order of ``phase_key``.
    * Species within each phase in alphabetical order of ``species_key``.

    Returns a 1-D ``float64`` array of mole values. Thin wrapper over
    :class:`~PyOMES.core.state_vector.StateVector` (STEP_SOLVER_INTERFACE_REFINEMENT.md
    item 4) — kept as a standalone function so every existing caller's
    signature is unaffected.
    """
    from .state_vector import StateVector
    return StateVector(sim.cvs).pack()


def _unpack_state(vec: np.ndarray, sim: "Simulation") -> None:
    """Write a packed state vector back into the simulation's CV phases.

    Inverse of :func:`_pack_state`.  Mutates ``phase.n_mol`` in-place.
    ``vec`` must have the same length and ordering as a vector produced
    by ``_pack_state(sim)``.
    """
    from .state_vector import StateVector
    StateVector(sim.cvs).unpack(vec)


def _state_index_map(sim: "Simulation") -> Dict[tuple, int]:
    """Return ``{(cv_key, phase_key, species_key): index}`` for a packed vector.

    Convenience for tests and Phase D matrix assembly.
    """
    from .state_vector import StateVector
    return StateVector(sim.cvs).index_map


# ════════════════════════════════════════════════════════════════════════
#  ExplicitEulerSystemSolver
# ════════════════════════════════════════════════════════════════════════

class ExplicitEulerSystemSolver:
    """Transparent wrapper around ``Simulation._step_default()``.

    Passing this solver explicitly produces byte-for-byte identical
    results to the default (``system_solver=None``) path.  Its sole
    purpose is to validate that the :class:`SystemSolver` dispatch
    wiring is correct.
    """

    def advance_system(
        self,
        sim: "Simulation",
        dt_h: float,
        t_h: float,
    ) -> Any:
        return sim._step_default(dt_h, t_h)


# ════════════════════════════════════════════════════════════════════════
#  StrangSplittingSystemSolver
# ════════════════════════════════════════════════════════════════════════

class StrangSplittingSystemSolver:
    """Second-order Strang splitting of link-flux from per-CV chemistry.

    Sequence per macro step::

        profiles(t) → links(dt/2) → cv.advance(dt) → links(dt/2) → controllers

    Cost over explicit Euler: one extra ``_apply_links`` call per step —
    essentially free compared to speciation/kinetic solves.

    Periodic controllers (``update_period_h`` set) fire once at the macro-step
    endpoint.  Sub-step controller events (T_c < dt) are Phase E territory.
    """

    def advance_system(
        self,
        sim: "Simulation",
        dt_h: float,
        t_h: float,
    ) -> Any:
        # 1. Apply profiles at step start (same as default Euler).
        profile_actions = sim._invoke_profiles(t_h=t_h)

        # 2. First half-step of link flux.
        link_records = sim._apply_links(dt_h / 2.0)

        # 3. Full CV advance (chemistry).
        results: Dict[str, "AdvanceResult"] = {}
        for cv_key, cv in sim.cvs.items():
            solver = sim._solver_for(cv_key)
            # Trusted orchestrator path - bypasses the OrchestrationWarning
            # ownership guard on the public cv.advance().
            results[cv_key] = cv._advance_unchecked(dt_h, t_h, solver=solver)

        # 4. Second half-step of link flux.
        link_records += sim._apply_links(dt_h / 2.0)

        # 5. Controllers observe fully-resolved post-step state.
        controller_actions = [
            a for a in sim._invoke_controllers(t_h=t_h, dt_h=dt_h, results=results)
            if a is not None
        ]
        for action in controller_actions:
            sim._apply_controller_action(action, dt_h=dt_h)

        return results, link_records, controller_actions, profile_actions


# ════════════════════════════════════════════════════════════════════════
#  MultirateSystemSolver  (Phase C)
# ════════════════════════════════════════════════════════════════════════

class MultirateSystemSolver:
    """Adaptive subcycling: M link evaluations at ``dt_h / M`` per macro step.

    M is chosen each macro step from the CFL number so that each
    transport sub-step is CFL-stable::

        M = ceil(dt_h / min_τ_CFL)

    where ``τ_CFL = V_source_phase / Q`` for :class:`~PyOMES.core.links.AdvectiveLink`
    and ``τ_CFL = V_source_phase / (kLa_max × V_eff)`` for
    :class:`~PyOMES.core.links.DiffusiveLink`.  Falls back to ``M = 1``
    when no links define a CFL timescale (same cost as explicit Euler).

    Sequence per macro step::

        profiles(t) → cv.advance(dt) → [links(dt/M)] × M → controllers

    Per-CV chemistry is advanced once for the full ``dt_h``; the CFL
    instability lives entirely in the inter-CV transport, which is
    subcycled at the sub-step scale.
    """

    def __init__(self, max_subcycles: int = 1000) -> None:
        self._max_subcycles = max_subcycles

    def _compute_M(self, sim: "Simulation", dt_h: float) -> int:
        """Compute the number of link sub-steps required for CFL stability."""
        min_tau = math.inf
        for link in sim.links:
            src_cv = sim.cvs.get(link.source_cv_key)
            if src_cv is None:
                continue
            src_phase = src_cv.phases.get(link.source_phase_key)
            if src_phase is None:
                continue
            V_src = float(src_phase.V_L)

            # AdvectiveLink: τ = V_src / Q
            Q = getattr(link, "Q_L_per_h", None)
            if Q is not None and float(Q) > 0.0:
                tau = V_src / float(Q)
                if tau > 0.0:
                    min_tau = min(min_tau, tau)
                continue

            # DiffusiveLink: τ = V_src / (kLa_max * V_eff)
            kLa_dict = getattr(link, "kLa", None)
            if kLa_dict:
                kLa_max = max(float(v) for v in kLa_dict.values() if float(v) > 0.0)
                if kLa_max > 0.0:
                    V_eff_attr = getattr(link, "V_eff_L", None)
                    if V_eff_attr is not None:
                        V_eff = max(float(V_eff_attr), 1e-30)
                    else:
                        snk_cv = sim.cvs.get(link.sink_cv_key)
                        if snk_cv is not None:
                            snk_phase = snk_cv.phases.get(link.sink_phase_key)
                            V_snk = float(snk_phase.V_L) if snk_phase else V_src
                        else:
                            V_snk = V_src
                        V_eff = 2.0 * V_src * V_snk / (V_src + V_snk)
                    tau = V_src / (kLa_max * V_eff)
                    if tau > 0.0:
                        min_tau = min(min_tau, tau)

        if math.isinf(min_tau) or min_tau <= 0.0:
            return 1
        M = math.ceil(dt_h / min_tau)
        return max(1, min(M, self._max_subcycles))

    def advance_system(
        self,
        sim: "Simulation",
        dt_h: float,
        t_h: float,
    ) -> Any:
        # 1. Apply profiles at step start.
        profile_actions = sim._invoke_profiles(t_h=t_h)

        # 2. Advance CVs once for the full macro step (chemistry is not CFL-stiff).
        results: Dict[str, "AdvanceResult"] = {}
        for cv_key, cv in sim.cvs.items():
            solver = sim._solver_for(cv_key)
            # Trusted orchestrator path - bypasses the OrchestrationWarning
            # ownership guard on the public cv.advance().
            results[cv_key] = cv._advance_unchecked(dt_h, t_h, solver=solver)

        # 3. Subcycle link flux at CFL-safe sub-steps.
        M = self._compute_M(sim, dt_h)
        dt_sub = dt_h / M
        link_records = []
        for _ in range(M):
            link_records += sim._apply_links(dt_sub)

        # 4. Controllers observe fully-resolved post-step state.
        controller_actions = [
            a for a in sim._invoke_controllers(t_h=t_h, dt_h=dt_h, results=results)
            if a is not None
        ]
        for action in controller_actions:
            sim._apply_controller_action(action, dt_h=dt_h)

        return results, link_records, controller_actions, profile_actions


# ════════════════════════════════════════════════════════════════════════
#  ImplicitTransportSystemSolver  (Phase D)
# ════════════════════════════════════════════════════════════════════════

class ImplicitTransportSystemSolver:
    """IMEX Axis 2 solver: inter-CV transport implicit, per-CV chemistry explicit.

    Lie–Trotter IMEX split per macro step::

        1. profiles(t)
        2. cv.advance(dt) for each CV   [chemistry, Axis 1]
        3. assemble sparse A from link parameters
        4. solve (I + dt·A)·n* = n_after_chem
        5. write n* back to CV phases
        6. controllers at macro-step boundary

    Transport is unconditionally stable — no CFL constraint — at the cost of
    one LU factorisation per unique (topology, dt_h) pair.  For a fixed-topology
    run the factorisation is cached; subsequent steps pay only the O(N²)
    backsolve.

    Matrix convention — positive diagonal means net *loss* from that entry:

    **AdvectiveLink** (Q L h⁻¹, CV i → CV j):
      A[i,s, i,s] += Q / V_i   (source drains)
      A[j,s, i,s] -= Q / V_i   (sink gains)

    **DiffusiveLink** (kLa per species, between CV i and CV j):
      E_s = kLa_s × V_eff
      A[i,s, i,s] += E_s/V_i    A[i,s, j,s] -= E_s/V_j
      A[j,s, j,s] += E_s/V_j    A[j,s, i,s] -= E_s/V_i

    Note: the off-diagonal subscripts in IMPLICIT_TRANSPORT.md lines 77–84 are
    transposed; the formulas above are correct (derived from
    dn_i/dt = −(E/V_i)·n_i + (E/V_j)·n_j).

    Controller events
    -----------------
    Controllers fire at macro-step boundaries via ``sim._invoke_controllers``
    — identical to StrangSplittingSystemSolver and MultirateSystemSolver.
    ``EventScheduler`` (Phase E) is not used; T_c < dt_h is not supported.
    """

    def __init__(self) -> None:
        self._cached_topo_hash: Optional[int] = None
        self._cached_lu: Optional[Any] = None

    # ── Internal helpers ───────────────────────────────────────────────

    def _topo_hash(self, sim: "Simulation", dt_h: float) -> int:
        """Hash over link topology, CV volumes, and dt_h for cache keying."""
        link_tuples = tuple(
            (
                type(lnk).__name__,
                lnk.source_cv_key,
                lnk.source_phase_key,
                lnk.sink_cv_key,
                lnk.sink_phase_key,
                getattr(lnk, "Q_L_per_h", None),
                tuple(sorted(getattr(lnk, "kLa", {}).items())),
                getattr(lnk, "V_eff_L", None),
            )
            for lnk in sim.links
        )
        vol_tuples = tuple(
            (ck, pk, float(ph.V_L))
            for ck, cv in sim.cvs.items()
            for pk, ph in sorted(cv.phases.items())
        )
        return hash((link_tuples, vol_tuples, dt_h))

    def _assemble_transport_matrix(self, sim: "Simulation") -> Any:
        """Build the sparse transport rate matrix A (CSR format).

        Reads link attributes directly; does not call ``compute_flow()``.
        ``compute_flow`` returns clamped fluxes, not the rate coefficients
        needed for matrix assembly.
        """
        import scipy.sparse as sp
        from .links import AdvectiveLink, DiffusiveLink

        idx = _state_index_map(sim)
        N = len(idx)

        rows: List[int] = []
        cols: List[int] = []
        data_vals: List[float] = []

        def _add(r: int, c: int, v: float) -> None:
            rows.append(r)
            cols.append(c)
            data_vals.append(v)

        for lnk in sim.links:
            src_cv = sim.cvs.get(lnk.source_cv_key)
            snk_cv = sim.cvs.get(lnk.sink_cv_key)
            if src_cv is None or snk_cv is None:
                continue
            src_phase = src_cv.phases.get(lnk.source_phase_key)
            snk_phase = snk_cv.phases.get(lnk.sink_phase_key)
            if src_phase is None or snk_phase is None:
                continue

            if isinstance(lnk, AdvectiveLink):
                Q = float(lnk.Q_L_per_h)
                if Q <= 0.0:
                    continue
                V_src = float(src_phase.V_L)
                if V_src <= 0.0:
                    continue
                rate = Q / V_src
                species = (
                    lnk.species_filter
                    if lnk.species_filter is not None
                    else list(src_phase.n_mol.keys())
                )
                for sk in species:
                    key_i = (lnk.source_cv_key, lnk.source_phase_key, sk)
                    key_j = (lnk.sink_cv_key, lnk.sink_phase_key, sk)
                    if key_i not in idx:
                        continue
                    i = idx[key_i]
                    _add(i, i, rate)
                    if key_j in idx:
                        _add(idx[key_j], i, -rate)

            elif isinstance(lnk, DiffusiveLink):
                V_src = float(src_phase.V_L)
                V_snk = float(snk_phase.V_L)
                if V_src <= 0.0 or V_snk <= 0.0:
                    continue
                V_eff = (
                    max(float(lnk.V_eff_L), 1e-30)
                    if lnk.V_eff_L is not None
                    else 2.0 * V_src * V_snk / (V_src + V_snk)
                )
                for sk, kla in lnk.kLa.items():
                    kla_f = float(kla)
                    if kla_f <= 0.0:
                        continue
                    E = kla_f * V_eff
                    key_i = (lnk.source_cv_key, lnk.source_phase_key, sk)
                    key_j = (lnk.sink_cv_key, lnk.sink_phase_key, sk)
                    if key_i not in idx or key_j not in idx:
                        continue
                    i = idx[key_i]
                    j = idx[key_j]
                    _add(i, i,  E / V_src)
                    _add(i, j, -E / V_snk)
                    _add(j, j,  E / V_snk)
                    _add(j, i, -E / V_src)

        if not data_vals:
            return sp.csr_matrix((N, N), dtype=np.float64)
        A = sp.coo_matrix(
            (data_vals, (rows, cols)), shape=(N, N), dtype=np.float64
        )
        return A.tocsr()

    def _get_lu(self, sim: "Simulation", dt_h: float) -> Any:
        """Return cached SuperLU for (I + dt·A); recompute on topology/dt change."""
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla

        h = self._topo_hash(sim, dt_h)
        if h != self._cached_topo_hash:
            A = self._assemble_transport_matrix(sim)
            N = A.shape[0]
            M = (
                sp.eye(N, format="csc", dtype=np.float64)
                + dt_h * A.tocsc()
            )
            self._cached_lu = spla.splu(M)
            self._cached_topo_hash = h
        return self._cached_lu

    # ── SystemSolver protocol ──────────────────────────────────────────

    def advance_system(
        self,
        sim: "Simulation",
        dt_h: float,
        t_h: float,
    ) -> Any:
        # 1. Profiles at step start.
        profile_actions = sim._invoke_profiles(t_h=t_h)

        # 2. Per-CV chemistry (Axis 1).
        results: Dict[str, "AdvanceResult"] = {}
        for cv_key, cv in sim.cvs.items():
            solver = sim._solver_for(cv_key)
            # Trusted orchestrator path - bypasses the OrchestrationWarning
            # ownership guard on the public cv.advance().
            results[cv_key] = cv._advance_unchecked(dt_h, t_h, solver=solver)

        # 3–5. Implicit transport solve (skip when no links present).
        if sim.links:
            rhs = _pack_state(sim)
            n_star = self._get_lu(sim, dt_h).solve(rhs)
            _unpack_state(n_star, sim)

        # 6. Controllers at macro-step boundary.
        controller_actions = [
            a
            for a in sim._invoke_controllers(t_h=t_h, dt_h=dt_h, results=results)
            if a is not None
        ]
        for action in controller_actions:
            sim._apply_controller_action(action, dt_h=dt_h)

        return results, [], controller_actions, profile_actions


# ════════════════════════════════════════════════════════════════════════
#  Extended state packing (Phase E — shared with MonolithicODESolver)
# ════════════════════════════════════════════════════════════════════════

def _ctrl_list(sim: "Simulation") -> List[tuple]:
    """Return [(ctrl, sorted_keys)] for every controller with differential state.

    Only controllers where ``differential_state()`` returns a non-empty dict
    are included.  Iteration order matches ``sim.controllers``; keys within
    each controller are sorted alphabetically so the packing convention is
    deterministic.
    """
    result: List[tuple] = []
    for ctrl in sim.controllers:
        ds_fn = getattr(ctrl, "differential_state", None)
        if ds_fn is None:
            continue
        state = ds_fn()
        if state:
            result.append((ctrl, sorted(state.keys())))
    return result


def _pack_extended_state(sim: "Simulation", ctrl_list: List[tuple]) -> np.ndarray:
    """Pack CV species + controller differential states into one flat vector.

    Layout::

        y = [ CV species block (Phase C order) | ctrl_0 states | ctrl_1 states | … ]

    The CV block is identical to :func:`_pack_state`.  Controller states are
    appended in ``sim.controllers`` order (filtered to those with differential
    state), with keys alphabetical within each controller. Thin wrapper over
    :class:`~PyOMES.core.state_vector.StateVector`.
    """
    from .state_vector import StateVector
    return StateVector(sim.cvs, ctrl_list=ctrl_list).pack()


def _unpack_extended_state(
    y: np.ndarray,
    sim: "Simulation",
    ctrl_list: List[tuple],
    n_cv: int,
) -> None:
    """Write a packed extended state vector back into CV phases and controllers.

    Inverse of :func:`_pack_extended_state`. ``n_cv`` (the length of the
    CV species block) is accepted for signature compatibility with
    existing callers that precompute it, but is no longer needed
    internally — :class:`~PyOMES.core.state_vector.StateVector` computes
    its own offset from ``sim.cvs``.
    """
    from .state_vector import StateVector
    StateVector(sim.cvs, ctrl_list=ctrl_list).unpack(y)


# ════════════════════════════════════════════════════════════════════════
#  RHS assembler (Phase E)
# ════════════════════════════════════════════════════════════════════════

def _build_rhs(sim: "Simulation", ctrl_list: List[tuple], n_cv: int):
    """Return the ODE RHS closure for :class:`MonolithicODESolver`.

    The returned callable ``rhs(t_h, y) -> dydt`` assembles three
    contributions:

    1. Per-CV kinetic chemistry rates from ``cv.compute_rhs()``.
    2. Inter-CV transport rates read directly from link parameters
       (bypasses ``compute_flow()`` to avoid the discrete-step safety cap).
    3. Controller differential-state rates from ``ctrl.state_rates()``.

    ``n_cv`` must equal ``len(_pack_state(sim))`` and is precomputed by
    the caller so it is not recomputed on every RHS evaluation.
    """
    from .links import AdvectiveLink, DiffusiveLink
    from .system_env import SystemEnv

    idx = _state_index_map(sim)  # (cv_key, phase_key, species_key) → int
    n_ctrl = sum(len(keys) for _, keys in ctrl_list)
    n_total = n_cv + n_ctrl

    def rhs(t_h: float, y: np.ndarray) -> np.ndarray:
        _unpack_extended_state(y, sim, ctrl_list, n_cv)

        dydt = np.zeros(n_total, dtype=np.float64)

        # 1. Per-CV chemistry rates.
        for cv_key, cv in sim.cvs.items():
            rates = cv.compute_rhs(float(t_h))
            for phase_key, species_rates in rates.items():
                for species_key, rate in species_rates.items():
                    i = idx.get((cv_key, phase_key, species_key))
                    if i is not None:
                        dydt[i] += rate

        # 2. Inter-CV transport rates — analytical from link parameters.
        for lnk in sim.links:
            src_cv = sim.cvs.get(lnk.source_cv_key)
            snk_cv = sim.cvs.get(lnk.sink_cv_key)
            if src_cv is None or snk_cv is None:
                continue
            src_phase = src_cv.phases.get(lnk.source_phase_key)
            snk_phase = snk_cv.phases.get(lnk.sink_phase_key)
            if src_phase is None or snk_phase is None:
                continue

            if isinstance(lnk, AdvectiveLink):
                Q = float(lnk.Q_L_per_h)
                if Q <= 0.0:
                    continue
                V_src = float(src_phase.V_L)
                if V_src <= 0.0:
                    continue
                species = (
                    lnk.species_filter
                    if lnk.species_filter is not None
                    else list(src_phase.n_mol.keys())
                )
                for sp in species:
                    n = float(src_phase.n_mol.get(sp, 0.0))
                    if n <= 0.0:
                        continue
                    flux = (n / V_src) * Q  # mol/h
                    i_src = idx.get((lnk.source_cv_key, lnk.source_phase_key, sp))
                    i_snk = idx.get((lnk.sink_cv_key, lnk.sink_phase_key, sp))
                    if i_src is not None:
                        dydt[i_src] -= flux
                    if i_snk is not None:
                        dydt[i_snk] += flux

            elif isinstance(lnk, DiffusiveLink):
                V_src = float(src_phase.V_L)
                V_snk = float(snk_phase.V_L)
                if V_src <= 0.0 or V_snk <= 0.0:
                    continue
                V_eff = (
                    max(float(lnk.V_eff_L), 1e-30)
                    if lnk.V_eff_L is not None
                    else 2.0 * V_src * V_snk / (V_src + V_snk)
                )
                for sp, kla in lnk.kLa.items():
                    kla_f = float(kla)
                    if kla_f <= 0.0:
                        continue
                    C_src = float(src_phase.n_mol.get(sp, 0.0)) / V_src
                    C_snk = float(snk_phase.n_mol.get(sp, 0.0)) / V_snk
                    flux = kla_f * V_eff * (C_src - C_snk)  # mol/h
                    i_src = idx.get((lnk.source_cv_key, lnk.source_phase_key, sp))
                    i_snk = idx.get((lnk.sink_cv_key, lnk.sink_phase_key, sp))
                    if i_src is not None:
                        dydt[i_src] -= flux
                    if i_snk is not None:
                        dydt[i_snk] += flux

        # 3. Controller differential-state rates.
        if ctrl_list:
            env = SystemEnv(cvs=dict(sim.cvs), t_h=float(t_h))
            offset = n_cv
            for ctrl, keys in ctrl_list:
                rates = ctrl.state_rates(env, float(t_h))
                for i, k in enumerate(keys):
                    dydt[offset + i] = float(rates.get(k, 0.0))
                offset += len(keys)

        return dydt

    return rhs


# ════════════════════════════════════════════════════════════════════════
#  MonolithicODESolver  (Phase E)
# ════════════════════════════════════════════════════════════════════════

_MONO_TOL: float = 1e-10  # floating-point tolerance for event boundary checks


class MonolithicODESolver:
    """Co-integrating Axis 2 solver for tight control loops and multi-rate events.

    Integrates all CV species concentrations, controller integral states,
    and inter-CV transport in a single ``scipy.solve_ivp`` call, firing
    periodic controllers at exact T_c boundaries via an outer loop over
    integration segments.

    Use this solver when:

    * A controller has ``update_period_h`` < ``dt_h`` — macro-step ZOH is
      too coarse and integral-state accuracy matters.
    * Per-CV kinetics are stiff and benefit from adaptive sub-stepping
      at the system level (select ``method="LSODA"``).

    For slow regulatory control (T_c ≫ dt_h), the existing default
    ``ExplicitEulerSystemSolver`` is sufficient and much cheaper.

    Parameters
    ----------
    method : str
        ODE method passed to ``scipy.solve_ivp``. ``"RK45"`` (default)
        for smooth problems; ``"LSODA"`` auto-selects stiff/non-stiff and
        is appropriate when per-CV kinetics are stiff.
    default_period_h : float or None
        Fallback update period (hours) for controllers that do not declare
        ``update_period_h``.  If ``None`` (default) and any controller
        lacks a positive ``update_period_h``, :meth:`advance_system`
        raises :exc:`ValueError` at entry — declare the period on the
        controller or pass this fallback.

    Notes
    -----
    Transport appears as explicit derivatives in the ODE RHS, computed
    analytically from link parameters rather than via ``compute_flow()``.
    This is appropriate when the solver can take sub-steps shorter than
    the CFL timescale.  For unconditionally stable transport without
    tight controllers, prefer :class:`ImplicitTransportSystemSolver`.

    ``advance_system`` returns ``results = {}``; per-step
    ``AdvanceResult`` kinetics are not populated (``cv.advance()`` is not
    called — the ODE integrator owns the full step).
    """

    def __init__(
        self,
        method: str = "RK45",
        default_period_h: Optional[float] = None,
    ) -> None:
        self.method = method
        self.default_period_h = default_period_h

    # ── Internal helpers ───────────────────────────────────────────────

    def _effective_period(self, ctrl) -> float:
        """Resolve update period for *ctrl*; raise if none can be found."""
        period = getattr(ctrl, "update_period_h", None)
        if period is None:
            period = self.default_period_h
        if period is None:
            raise ValueError(
                f"Controller {ctrl!r} has no update_period_h and "
                f"MonolithicODESolver.default_period_h is None. "
                f"Declare update_period_h on the controller or pass "
                f"default_period_h to MonolithicODESolver()."
            )
        return float(period)

    def _fire_controllers(
        self,
        sim: "Simulation",
        due_ctrls: List,
        t_h: float,
        dt_h: float,
        zoh_cache: Dict,
    ) -> List:
        """Build snapshot, call compute() for each due controller, cache result."""
        from .snapshot import build_simulation_snapshot

        if not due_ctrls:
            return []
        sim_snap = build_simulation_snapshot(sim, t_h=t_h, results={})
        actions = []
        for ctrl in due_ctrls:
            tcv = getattr(ctrl, "target_cv_key", None)
            if tcv is not None and tcv in sim_snap.cvs:
                state = sim_snap.cvs[tcv]
            elif tcv is None and len(sim_snap.cvs) == 1:
                state = next(iter(sim_snap.cvs.values()))
            else:
                state = sim_snap
            action = ctrl.compute(state, dt_h)
            zoh_cache[id(ctrl)] = action
            actions.append(action)
        return actions

    # ── SystemSolver protocol ──────────────────────────────────────────

    def advance_system(
        self,
        sim: "Simulation",
        dt_h: float,
        t_h: float,
    ) -> Any:
        from scipy.integrate import solve_ivp

        # Unlike every other SystemSolver, MonolithicODESolver never calls
        # cv.advance()/_advance_unchecked() — it integrates every CV
        # directly via cv.compute_rhs() inside one shared solve_ivp call
        # (see _build_rhs() above). A per-CV sim.solver=... is therefore
        # never consulted; catch it here rather than let it be silently
        # moot. Checked at entry (not construction/setter time) so it is
        # correct regardless of which of sim.solver / sim.system_solver
        # is assigned first or reassigned later.
        if sim.solver is not None:
            raise ValueError(
                "MonolithicODESolver integrates every CV directly via "
                "cv.compute_rhs() inside one shared solve_ivp() call; it "
                "never calls cv.advance(), so a per-CV Simulation(solver=...) "
                "is never consulted. Remove solver=... or use a different "
                "system_solver."
            )

        # 1. Profiles at step start.
        profile_actions = sim._invoke_profiles(t_h=t_h)

        # 2. Validate and resolve all controller periods early so we fail fast.
        ctrl_periods: Dict[int, float] = {}
        for ctrl in sim.controllers:
            ctrl_periods[id(ctrl)] = self._effective_period(ctrl)

        # 3. Precompute state layout and build RHS closure.
        ctrl_list = _ctrl_list(sim)
        n_cv = len(_pack_state(sim))
        rhs = _build_rhs(sim, ctrl_list, n_cv)

        # 4. Build per-controller schedule: {id(ctrl): {"ctrl", "period", "next_fire"}}.
        t_end_macro = t_h + dt_h
        ctrl_schedule: List[Dict] = [
            {
                "ctrl": ctrl,
                "period": ctrl_periods[id(ctrl)],
                "next_fire": t_h + ctrl_periods[id(ctrl)],
            }
            for ctrl in sim.controllers
        ]

        # 5. Outer integration loop.
        y = _pack_extended_state(sim, ctrl_list)
        t_now = t_h
        zoh_cache: Dict[int, Any] = {}
        controller_actions: List = []

        while t_now < t_end_macro - _MONO_TOL:
            # Nearest upcoming controller event within this macro step.
            if ctrl_schedule:
                t_next_event = min(
                    min(s["next_fire"] for s in ctrl_schedule),
                    t_end_macro,
                )
            else:
                t_next_event = t_end_macro
            t_seg_end = min(t_end_macro, t_next_event)

            sol = solve_ivp(
                rhs,
                [t_now, t_seg_end],
                y,
                method=self.method,
                dense_output=False,
            )
            if not sol.success:
                raise RuntimeError(
                    f"MonolithicODESolver: solve_ivp failed at "
                    f"t={t_now:.6g}: {sol.message}"
                )

            y = sol.y[:, -1]
            t_now = sol.t[-1]

            # Write integrated state back into the simulation.
            _unpack_extended_state(y, sim, ctrl_list, n_cv)

            # Fire controllers whose next_fire ≤ t_now.
            due_entries = [
                s for s in ctrl_schedule
                if s["next_fire"] <= t_now + _MONO_TOL
            ]
            if due_entries:
                due_ctrls = [s["ctrl"] for s in due_entries]
                acts = self._fire_controllers(sim, due_ctrls, t_now, dt_h, zoh_cache)
                controller_actions.extend(acts)
                for action in acts:
                    sim._apply_controller_action(action, dt_h=dt_h)
                # Re-schedule fired controllers.
                for s in due_entries:
                    s["next_fire"] = t_now + s["period"]
                # Re-pack to capture any post-action state changes.
                y = _pack_extended_state(sim, ctrl_list)

        return {}, [], controller_actions, profile_actions
