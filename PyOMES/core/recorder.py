# -*- coding: utf-8 -*-
"""Recorder protocol and default :class:`BatchRecorder`.

The :class:`Recorder` protocol is the pluggable per-step time-series
capture surface for :class:`~PyOMES.core.simulation.Simulation`. The
default :class:`BatchRecorder` pre-allocates per-CV NumPy arrays at
``record_init`` time and writes per-step into them; ``finalize``
returns a :class:`BatchResult`.

This module lands C3 of SIMULATION_CLASS. ``Simulation.run()``
wires the recorder at C4. Controller / profile / accuracy /
conservation channels exist on the result but are populated only
once their respective checkpoints land (C9/C10).

Note: the legacy ``BatchResult`` in
[`models/vlmodels/fermenter/config/factory.py`](../../models/vlmodels/fermenter/config/factory.py)
is **not** moved here. The two classes coexist on the
simulation-class branch until C14 deletes ``run_batch`` and its
``BatchResult``; this avoids a backwards-compat shim per the
hard-breaks discipline. The two have slightly different schemas
(the new one is nested per-CV; the legacy one is flat single-CV).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Protocol, runtime_checkable

import numpy as np

from .interfaces import AdvanceResult
from .links import LinkFlowRecord
from .boundaries import ExternalFluxRecord
from ..control.actions import ControlAction, ProfileRecord


# ════════════════════════════════════════════════════════════════════════
#  Recorder protocol
# ════════════════════════════════════════════════════════════════════════


@runtime_checkable
class Recorder(Protocol):
    """Per-step recorder interface consumed by Simulation.run().

    Three methods, called in this order over a single ``.run()``:

    1. ``record_init(sim, n_steps)`` — once at run entry. Allows the
       recorder to inspect the topology, pre-allocate storage, and
       record the t=0 state.
    2. ``record_step(step_index, t_h, dt_h, advance_results,
       controller_actions, profile_actions)`` — once per step,
       called after ``cv.advance(...)`` has returned for every CV
       and after controllers/profiles have applied for the step.
    3. ``finalize()`` — once at run exit (in a ``try/finally``).
       Returns whatever result object the recorder produces.

    Future recorders (``SparseRecorder``, ``StreamingRecorder``,
    ``SummaryRecorder``, etc.) implement the same three methods.
    """

    def record_init(self, sim: Any, n_steps: int, start_t_h: float = 0.0) -> None: ...

    def record_step(
        self,
        step_index: int,
        t_h: float,
        dt_h: float,
        advance_results: Dict[str, AdvanceResult],
        link_records: List[LinkFlowRecord],
        controller_actions: List[ControlAction],
        profile_actions: List[ProfileRecord],
    ) -> None: ...

    def finalize(self) -> Any: ...


# ════════════════════════════════════════════════════════════════════════
#  BatchResult dataclass
# ════════════════════════════════════════════════════════════════════════


@dataclass
class BatchResult:
    """Time-series results returned by :meth:`BatchRecorder.finalize`.

    Nested per-CV by construction (single-CV simulations populate
    every nested dict under their sole CV key — typically
    ``"main"``). This is intentionally different from the legacy
    :class:`models.vlmodels.fermenter.config.factory.BatchResult`
    (flat single-CV) — the two coexist until C14 deletes the
    legacy.

    Attributes
    ----------
    t_h : np.ndarray
        Time grid of shape ``(n_steps + 1,)`` including ``t=0``.
    phase_mol : dict
        ``{cv_key: {phase_key: {species_id: ndarray}}}`` — primary
        per-CV, per-phase mole storage. Covers every phase the CV
        exposes, including non-standard names (``"mobile"``,
        ``"solid"``, etc.).
    gas_mol : dict
        ``{cv_key: {species_id: ndarray}}`` — computed alias for
        ``phase_mol[cv_key].get("gas", {})``. Read-only; do not
        assign into the returned dict.
    liquid_mol : dict
        ``{cv_key: {species_id: ndarray}}`` — computed alias for
        ``phase_mol[cv_key].get("liquid", {})``. Read-only.
    pH : dict
        ``{cv_key: ndarray}`` per-CV pH series. NaN where
        speciation was not computed.
    ionic_strength : dict
        ``{cv_key: ndarray}`` per-CV ionic strength series.
    P_atm : dict
        ``{cv_key: ndarray}`` per-CV total gas pressure series.
    transfer_records : list
        Per-step list of :class:`LinkFlowRecord`. Each step's list
        aggregates every link flow active that step (internal
        gas-liquid + inter-CV).
    boundary_records : list
        Per-step ``{cv_key: List[ExternalFluxRecord]}``.
    advance_results : list
        Per-step ``{cv_key: AdvanceResult}`` for inspection.
    controller_actions : list
        Per-step ``List[ControlAction]`` from every firing
        controller. Empty until C9.
    profile_actions : list
        Per-step ``List[ProfileRecord]`` from every firing
        profile. Empty until C10.
    accuracy_records : list
        Per-step accuracy-monitor outputs. Empty until wired in
        C4 / C9.
    conservation_records : list
        Per-step conservation-monitor outputs (element + charge
        accounting). Empty until wired in C4 / C9.
    runtime_s : float
        Wall-clock time from ``record_init`` to ``finalize``.
    """

    t_h: np.ndarray = field(default_factory=lambda: np.array([]))
    phase_mol: Dict[str, Dict[str, Dict[str, np.ndarray]]] = field(default_factory=dict)
    pH: Dict[str, np.ndarray] = field(default_factory=dict)
    ionic_strength: Dict[str, np.ndarray] = field(default_factory=dict)
    P_atm: Dict[str, np.ndarray] = field(default_factory=dict)
    transfer_records: List[List[LinkFlowRecord]] = field(default_factory=list)
    boundary_records: List[Dict[str, List[ExternalFluxRecord]]] = field(default_factory=list)
    advance_results: List[Dict[str, AdvanceResult]] = field(default_factory=list)
    controller_actions: List[List[ControlAction]] = field(default_factory=list)
    profile_actions: List[List[ProfileRecord]] = field(default_factory=list)
    accuracy_records: List[List[Any]] = field(default_factory=list)
    conservation_records: List[List[Any]] = field(default_factory=list)
    runtime_s: float = 0.0

    # ── Convenience aliases for gas / liquid (computed on access) ──────

    @property
    def gas_mol(self) -> Dict[str, Dict[str, np.ndarray]]:
        return {cv_key: phases.get("gas", {})
                for cv_key, phases in self.phase_mol.items()}

    @property
    def liquid_mol(self) -> Dict[str, Dict[str, np.ndarray]]:
        return {cv_key: phases.get("liquid", {})
                for cv_key, phases in self.phase_mol.items()}

    # ── Export helpers ─────────────────────────────────────────────────

    def to_dataframe(self):
        """Return a long-form pandas DataFrame with one row per
        (t_h, cv_key, channel_kind, phase_key, species, value).

        Covers the five NumPy array channels:
          - ``phase_mol``   gas and liquid moles per species
          - ``pH``          per-CV pH series
          - ``ionic_strength``  per-CV ionic strength
          - ``P_atm``       per-CV total gas pressure

        ``phase_key`` is ``"gas"`` or ``"liquid"`` for mol rows and
        ``None`` for scalar channels.  ``species`` is ``None`` for
        scalar channels.

        Requires ``pandas``; raises ``ImportError`` with an install
        hint if it is absent.  The underlying NumPy arrays remain
        accessible without pandas via the dataclass fields.
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "pandas is required for BatchResult.to_dataframe(). "
                "Install it with: pip install pandas"
            )

        rows = []

        # phase_mol: all phases (gas, liquid, and any custom phase keys)
        for cv_key, phase_dict in self.phase_mol.items():
            for phase_key, species_dict in phase_dict.items():
                for species, arr in species_dict.items():
                    for i, t in enumerate(self.t_h):
                        rows.append((t, cv_key, "phase_mol", phase_key, species, arr[i]))

        # Scalar channels: pH, ionic_strength, P_atm
        for cv_key, arr in self.pH.items():
            for i, t in enumerate(self.t_h):
                rows.append((t, cv_key, "pH", None, None, arr[i]))

        for cv_key, arr in self.ionic_strength.items():
            for i, t in enumerate(self.t_h):
                rows.append((t, cv_key, "ionic_strength", None, None, arr[i]))

        for cv_key, arr in self.P_atm.items():
            for i, t in enumerate(self.t_h):
                rows.append((t, cv_key, "P_atm", None, None, arr[i]))

        return pd.DataFrame(
            rows,
            columns=["t_h", "cv_key", "channel_kind", "phase_key", "species", "value"],
        )

    def to_wide(self, cv_key: str, phase_key: str):
        """Return a wide-form DataFrame for one CV and phase.

        Index is ``t_h``; columns are species names.  Useful for
        plotting a single phase's time series directly.

        ``phase_key`` can be any phase key present in the CV (e.g.
        ``"gas"``, ``"liquid"``, ``"mobile"``, ``"solid"``).
        """
        df = self.to_dataframe()
        mask = (
            (df["cv_key"] == cv_key)
            & (df["channel_kind"] == "phase_mol")
            & (df["phase_key"] == phase_key)
        )
        return (
            df.loc[mask]
            .pivot(index="t_h", columns="species", values="value")
            .rename_axis(None, axis="columns")
        )

    @classmethod
    def concat(cls, *results: "BatchResult") -> "BatchResult":
        """Concatenate two or more ``BatchResult`` objects into one trajectory.

        ``t_h`` and all nested NumPy arrays are concatenated; the first time
        point of each subsequent result is dropped because it duplicates the
        last point of the preceding result in a typical resume workflow. List
        channels (``transfer_records``, ``advance_results``, etc.) are
        concatenated without adjustment (one entry per step, no overlap).

        Species or CV keys present in one result but absent from another are
        filled with zeros for the missing time range (``NaN`` for ``pH`` and
        ``ionic_strength``). ``runtime_s`` is the sum of all inputs.

        Parameters
        ----------
        *results : BatchResult
            Two or more results to merge. At least one is required.

        Returns
        -------
        BatchResult
            A new ``BatchResult`` covering the combined time range.
        """
        if not results:
            raise ValueError("BatchResult.concat requires at least one result.")
        if len(results) == 1:
            return results[0]

        out = cls()

        # t_h: drop index 0 from every result after the first
        out.t_h = np.concatenate(
            [results[0].t_h] + [r.t_h[1:] for r in results[1:]]
        )

        # Collect all CV keys across every result
        all_cv_keys: set = set()
        for r in results:
            all_cv_keys.update(r.phase_mol.keys())
            all_cv_keys.update(r.pH.keys())
            all_cv_keys.update(r.ionic_strength.keys())
            all_cv_keys.update(r.P_atm.keys())

        for cv_key in all_cv_keys:
            # phase_mol: collect all (phase_key, species) pairs across results
            phase_species: Dict[str, set] = {}
            for r in results:
                for pk, sp_dict in r.phase_mol.get(cv_key, {}).items():
                    phase_species.setdefault(pk, set()).update(sp_dict)

            out.phase_mol[cv_key] = {}
            for phase_key, species_set in phase_species.items():
                out.phase_mol[cv_key][phase_key] = {}
                for species in species_set:
                    parts = []
                    for idx, r in enumerate(results):
                        arr = r.phase_mol.get(cv_key, {}).get(phase_key, {}).get(species)
                        if arr is None:
                            arr = np.zeros(len(r.t_h))
                        parts.append(arr if idx == 0 else arr[1:])
                    out.phase_mol[cv_key][phase_key][species] = np.concatenate(parts)

            # Scalar per-CV channels: pH / ionic_strength (fill NaN) / P_atm (fill 0)
            for attr, fill in (("pH", np.nan), ("ionic_strength", np.nan), ("P_atm", 0.0)):
                parts = []
                for idx, r in enumerate(results):
                    arr = getattr(r, attr).get(cv_key)
                    if arr is None:
                        arr = np.full(len(r.t_h), fill)
                    parts.append(arr if idx == 0 else arr[1:])
                getattr(out, attr)[cv_key] = np.concatenate(parts)

        # List channels: plain concatenation (no index-0 overlap)
        for attr in (
            "transfer_records", "boundary_records", "advance_results",
            "controller_actions", "profile_actions",
            "accuracy_records", "conservation_records",
        ):
            combined = []
            for r in results:
                combined.extend(getattr(r, attr))
            setattr(out, attr, combined)

        out.runtime_s = sum(r.runtime_s for r in results)
        return out

    def to_csv(self, path, **kwargs) -> None:
        """Write results to a CSV file.

        Wraps :meth:`to_dataframe` and writes with
        ``float_format="%.10e"`` by default to preserve numerical
        precision.  Pass ``float_format=...`` in ``kwargs`` to
        override.

        Requires ``pandas``; raises ``ImportError`` with an install
        hint if it is absent.
        """
        kwargs.setdefault("float_format", "%.10e")
        kwargs.setdefault("index", False)
        self.to_dataframe().to_csv(path, **kwargs)

    def to_parquet(self, path, **kwargs) -> None:
        """Write results to a Parquet file.

        Wraps :meth:`to_dataframe` and writes via the ``pyarrow``
        engine.  Parquet is binary and columnar: smaller on disk than
        CSV, lossless for floats, and fast to read back with
        ``pandas.read_parquet()``.

        Requires both ``pandas`` and ``pyarrow``; raises
        ``ImportError`` with an install hint if either is absent.
        """
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            raise ImportError(
                "pyarrow is required for BatchResult.to_parquet(). "
                "Install it with: pip install pandas pyarrow"
            )
        self.to_dataframe().to_parquet(path, **kwargs)


# ════════════════════════════════════════════════════════════════════════
#  BatchRecorder
# ════════════════════════════════════════════════════════════════════════


class BatchRecorder:
    """Default in-memory recorder.

    Pre-allocates ``(n_steps + 1,)``-shaped NumPy arrays per
    (CV key, species) tuple at :meth:`record_init`, writes them
    in place during :meth:`record_step`, and returns the
    :class:`BatchResult` from :meth:`finalize`.

    New species appearing during the run are accommodated: arrays
    are allocated on first observation, with values prior to the
    first sighting left at the array's default (``0.0`` for
    species inventories, ``NaN`` for pH/IS).

    The recorder is **single-use**: ``record_init`` may be called
    at most once. To run a Simulation twice, construct a fresh
    recorder (or pass ``recorder=None`` to :meth:`Simulation.run`,
    which will default-construct one at run entry from C4 onward).
    """

    def __init__(self) -> None:
        self._initialised: bool = False
        self._finalized: bool = False
        self._sim: Any = None
        self._n_steps: int = 0
        self._start_time: float = 0.0
        self._result: BatchResult = BatchResult()

    # ── Lifecycle ──────────────────────────────────────────────────────

    def record_init(self, sim: Any, n_steps: int, start_t_h: float = 0.0) -> None:
        if self._initialised:
            raise RuntimeError(
                "BatchRecorder.record_init called twice on the same recorder. "
                "Construct a fresh BatchRecorder per Simulation.run() call."
            )
        if n_steps < 1:
            raise ValueError(
                f"BatchRecorder.record_init requires n_steps >= 1; got {n_steps!r}."
            )
        self._sim = sim
        self._n_steps = int(n_steps)
        n = self._n_steps + 1

        # Time grid: t_h[0] is the start of the run (start_t_h, default 0).
        # Indices 1..n are written by record_step.
        self._result.t_h = np.zeros(n)
        self._result.t_h[0] = float(start_t_h)

        # Per-CV pre-allocation: scan every phase on every CV.
        for cv_key, cv in sim.cvs.items():
            self._result.phase_mol[cv_key] = {}
            self._result.pH[cv_key] = np.full(n, np.nan)
            self._result.ionic_strength[cv_key] = np.full(n, np.nan)
            self._result.P_atm[cv_key] = np.zeros(n)

            for phase_key, phase in cv.phases.items():
                self._result.phase_mol[cv_key][phase_key] = {
                    sp: np.zeros(n) for sp in phase.n_mol
                }

        # Initial state at index 0 (t=0). Any chemistry the caller
        # solved before record_init is captured; pH[0] / I[0] stay
        # NaN if no speciation engine has populated H+ yet.
        self._record_cv_state(0)

        self._start_time = time.perf_counter()
        self._initialised = True

    def record_step(
        self,
        step_index: int,
        t_h: float,
        dt_h: float,
        advance_results: Dict[str, AdvanceResult],
        link_records: Optional[List[LinkFlowRecord]] = None,
        controller_actions: Optional[List[ControlAction]] = None,
        profile_actions: Optional[List[ProfileRecord]] = None,
    ) -> None:
        if not self._initialised:
            raise RuntimeError(
                "BatchRecorder.record_step called before record_init."
            )
        if self._finalized:
            raise RuntimeError(
                "BatchRecorder.record_step called after finalize."
            )
        if not (1 <= step_index <= self._n_steps):
            raise IndexError(
                f"BatchRecorder.record_step: step_index={step_index!r} "
                f"outside range 1..{self._n_steps}."
            )

        link_records = list(link_records or [])
        controller_actions = list(controller_actions or [])
        profile_actions = list(profile_actions or [])

        self._result.t_h[step_index] = float(t_h)
        self._record_cv_state(step_index)

        # Aggregate per-step parallel channels.
        self._result.advance_results.append(dict(advance_results))
        self._result.controller_actions.append(controller_actions)
        self._result.profile_actions.append(profile_actions)

        # Transfer records: inter-CV link records (from Simulation._step)
        # come first; then per-CV internal gas-liquid link records
        # (from each AdvanceResult.transfer_record).
        step_transfer: List[LinkFlowRecord] = list(link_records)
        step_boundary: Dict[str, List[ExternalFluxRecord]] = {}
        for cv_key, ar in advance_results.items():
            if ar.transfer_record is not None:
                step_transfer.append(ar.transfer_record)
            step_boundary[cv_key] = list(ar.boundary_records)
        self._result.transfer_records.append(step_transfer)
        self._result.boundary_records.append(step_boundary)

    def finalize(self) -> BatchResult:
        if not self._initialised:
            raise RuntimeError(
                "BatchRecorder.finalize called before record_init."
            )
        if self._finalized:
            raise RuntimeError(
                "BatchRecorder.finalize called twice on the same recorder."
            )
        self._result.runtime_s = time.perf_counter() - self._start_time
        self._finalized = True
        return self._result

    # ── Internal: CV state capture at one timestep index ──────────────

    def _record_cv_state(self, idx: int) -> None:
        """Read every CV's current state into the pre-allocated arrays."""
        n_alloc = self._n_steps + 1
        for cv_key, cv in self._sim.cvs.items():
            # All phases into phase_mol (on-demand alloc for new species).
            for phase_key, phase in cv.phases.items():
                phase_dict = self._result.phase_mol[cv_key].setdefault(phase_key, {})
                for sp, n in phase.n_mol.items():
                    arr = phase_dict.get(sp)
                    if arr is None:
                        arr = np.zeros(n_alloc)
                        phase_dict[sp] = arr
                    arr[idx] = float(n)

            # Gas-specific scalar: total pressure.
            gas = cv.phases.get("gas")
            if gas is not None:
                self._result.P_atm[cv_key][idx] = float(gas.P_atm)

            # Liquid-specific scalars: pH and ionic strength.
            liquid = cv.phases.get("liquid")
            if liquid is not None:
                if "H+" in liquid.n_mol:
                    try:
                        self._result.pH[cv_key][idx] = float(liquid.pH)
                    except ValueError:
                        pass
                # Ionic strength: PropertyCalculator output preferred,
                # legacy speciation key as fallback.
                is_val = liquid.properties.get("ionic_strength")
                if is_val is None:
                    is_val = liquid.speciation.get("IonicStrength")
                if is_val is not None:
                    self._result.ionic_strength[cv_key][idx] = float(is_val)


# ════════════════════════════════════════════════════════════════════════
#  Helpers
# ════════════════════════════════════════════════════════════════════════


def _collect_run_metadata() -> Dict[str, str]:
    import os
    import platform
    import sys

    try:
        from importlib.metadata import version as _pkg_ver
        vlsim_version = _pkg_ver("PyOMES")
    except Exception:
        vlsim_version = "unknown"
    return {
        "vlsim_version": vlsim_version,
        "python_version": sys.version,
        "platform": platform.platform(),
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS", "1"),
    }


# ════════════════════════════════════════════════════════════════════════
#  StreamingFileRecorder
# ════════════════════════════════════════════════════════════════════════


class StreamingFileRecorder:
    """Recorder that flushes Parquet chunks to disk every N steps.

    Buffers per-step state rows in memory and writes one Parquet file per
    flush to *path*.  A JSON manifest tracks chunk boundaries.  On crash,
    all flushed chunks are recoverable; the in-flight buffer is lost.

    Column layout: one row per step (wide format).  Columns are
    ``step_index``, ``t_h``, one ``{cv}|{phase}|{species}`` column per
    mole quantity, and ``{cv}|pH``, ``{cv}|ionic_strength``,
    ``{cv}|P_atm`` per CV.

    Use :func:`load_run` to reassemble chunks into a :class:`BatchResult`.

    Parameters
    ----------
    path : str or pathlib.Path
        Directory to write chunks and manifest into.  Created if absent.
    flush_every_n_steps : int
        Number of buffered rows before a flush.  Default 1000.
    compression : {"snappy", "zstd", "none"}
        Parquet compression codec.  ``"snappy"`` is the default.

    Notes
    -----
    ``pyarrow`` is required.  A clear :exc:`ImportError` is raised at
    construction time if it is absent.

    List channels (``advance_results``, ``controller_actions``, etc.) are
    not written to Parquet; the :class:`BatchResult` returned by
    :meth:`finalize` will have empty lists for those channels.
    """

    def __init__(
        self,
        path: Any,
        flush_every_n_steps: int = 1000,
        compression: str = "snappy",
    ) -> None:
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            raise ImportError(
                "pyarrow is required for StreamingFileRecorder. "
                "Install it with: pip install pyarrow"
            )
        from pathlib import Path as _Path
        self._path: Any = _Path(path)
        self._flush_every: int = int(flush_every_n_steps)
        self._compression: str = compression
        self._buffer: List[Dict[str, Any]] = []
        self._chunks: List[Dict[str, Any]] = []
        self._chunk_index: int = 0
        self._col_names: List[str] = []
        self._cv_phase_species: Dict[str, Dict[str, List[str]]] = {}
        self._initialised: bool = False
        self._finalized: bool = False
        self._start_time: float = 0.0
        self._sim: Any = None

    # ── Protocol ──────────────────────────────────────────────────────

    def record_init(self, sim: Any, n_steps: int, start_t_h: float = 0.0) -> None:
        if self._initialised:
            raise RuntimeError(
                "StreamingFileRecorder.record_init called twice. "
                "Construct a fresh recorder per Simulation.run() call."
            )
        self._sim = sim
        self._path.mkdir(parents=True, exist_ok=True)

        # Build column schema from initial CV/phase/species topology.
        self._col_names = ["step_index", "t_h"]
        for cv_key, cv in sim.cvs.items():
            self._cv_phase_species[cv_key] = {}
            for phase_key, phase in cv.phases.items():
                species_list = list(phase.n_mol.keys())
                self._cv_phase_species[cv_key][phase_key] = species_list
                for sp in species_list:
                    self._col_names.append(f"{cv_key}|{phase_key}|{sp}")
            for ch in ("pH", "ionic_strength", "P_atm"):
                self._col_names.append(f"{cv_key}|{ch}")

        self._write_manifest({
            "schema_version": 1,
            "cv_phase_species": {
                cv_key: {pk: list(sps) for pk, sps in phases.items()}
                for cv_key, phases in self._cv_phase_species.items()
            },
            "columns": self._col_names,
            "chunks": [],
            "n_steps_total": None,
            "metadata": _collect_run_metadata(),
        })

        self._buffer.append(self._capture_row(0, float(start_t_h)))
        self._start_time = time.perf_counter()
        self._initialised = True

    def record_step(
        self,
        step_index: int,
        t_h: float,
        dt_h: float,
        advance_results: Dict[str, AdvanceResult],
        link_records: Optional[List[LinkFlowRecord]] = None,
        controller_actions: Optional[List[ControlAction]] = None,
        profile_actions: Optional[List[ProfileRecord]] = None,
    ) -> None:
        if not self._initialised:
            raise RuntimeError(
                "StreamingFileRecorder.record_step called before record_init."
            )
        if self._finalized:
            raise RuntimeError(
                "StreamingFileRecorder.record_step called after finalize."
            )
        self._buffer.append(self._capture_row(step_index, float(t_h)))
        if len(self._buffer) >= self._flush_every:
            self._flush()

    def finalize(self) -> BatchResult:
        if not self._initialised:
            raise RuntimeError(
                "StreamingFileRecorder.finalize called before record_init."
            )
        if self._finalized:
            raise RuntimeError(
                "StreamingFileRecorder.finalize called twice."
            )
        if self._buffer:
            self._flush()

        manifest = self._read_manifest()
        manifest["chunks"] = self._chunks
        manifest["n_steps_total"] = (
            self._chunks[-1]["end_step"] if self._chunks else 0
        )
        manifest["runtime_s"] = time.perf_counter() - self._start_time
        self._write_manifest(manifest)

        self._finalized = True
        return load_run(self._path)

    # ── Pickle support ────────────────────────────────────────────────

    def __getstate__(self) -> Dict[str, Any]:
        state = self.__dict__.copy()
        # In-flight buffer is transient; flushed chunks on disk are durable.
        state["_buffer"] = []
        return state

    # ── Internal ──────────────────────────────────────────────────────

    def _capture_row(self, step_index: int, t_h: float) -> Dict[str, Any]:
        row: Dict[str, Any] = {"step_index": step_index, "t_h": t_h}
        for cv_key, cv in self._sim.cvs.items():
            for phase_key, phase in cv.phases.items():
                for sp in self._cv_phase_species.get(cv_key, {}).get(phase_key, []):
                    row[f"{cv_key}|{phase_key}|{sp}"] = float(phase.n_mol.get(sp, 0.0))

            gas = cv.phases.get("gas")
            row[f"{cv_key}|P_atm"] = float(gas.P_atm) if gas is not None else 0.0

            liquid = cv.phases.get("liquid")
            if liquid is not None:
                ph_val = float("nan")
                if "H+" in liquid.n_mol:
                    try:
                        ph_val = float(liquid.pH)
                    except (ValueError, AttributeError):
                        pass
                row[f"{cv_key}|pH"] = ph_val
                is_val = liquid.properties.get("ionic_strength")
                if is_val is None:
                    is_val = liquid.speciation.get("IonicStrength")
                row[f"{cv_key}|ionic_strength"] = (
                    float(is_val) if is_val is not None else float("nan")
                )
            else:
                row[f"{cv_key}|pH"] = float("nan")
                row[f"{cv_key}|ionic_strength"] = float("nan")
        return row

    def _flush(self) -> None:
        import pyarrow as pa
        import pyarrow.parquet as pq

        if not self._buffer:
            return

        col_lists: Dict[str, list] = {col: [] for col in self._col_names}
        for row in self._buffer:
            for col in self._col_names:
                col_lists[col].append(row.get(col, 0.0))

        arrays = []
        fields = []
        for col in self._col_names:
            if col == "step_index":
                arrays.append(pa.array(col_lists[col], type=pa.int64()))
                fields.append(pa.field(col, pa.int64()))
            else:
                arrays.append(pa.array(col_lists[col], type=pa.float64()))
                fields.append(pa.field(col, pa.float64()))

        table = pa.table(
            dict(zip(self._col_names, arrays)),
            schema=pa.schema(fields),
        )

        chunk_file = f"chunk_{self._chunk_index:05d}.parquet"
        pq.write_table(table, self._path / chunk_file, compression=self._compression)

        self._chunks.append({
            "index": self._chunk_index,
            "file": chunk_file,
            "start_step": int(self._buffer[0]["step_index"]),
            "end_step": int(self._buffer[-1]["step_index"]),
            "start_t_h": float(self._buffer[0]["t_h"]),
            "end_t_h": float(self._buffer[-1]["t_h"]),
            "n_rows": len(self._buffer),
        })
        self._chunk_index += 1
        self._buffer = []

    def _write_manifest(self, manifest: Dict[str, Any]) -> None:
        import json
        (self._path / "manifest.json").write_text(json.dumps(manifest, indent=2))

    def _read_manifest(self) -> Dict[str, Any]:
        import json
        return json.loads((self._path / "manifest.json").read_text())


# ════════════════════════════════════════════════════════════════════════
#  load_run
# ════════════════════════════════════════════════════════════════════════


def load_run(path: Any) -> BatchResult:
    """Reassemble :class:`StreamingFileRecorder` chunks into a :class:`BatchResult`.

    Reads the manifest and all Parquet chunks written by a
    :class:`StreamingFileRecorder`, then reconstructs the same nested
    structure as :class:`BatchRecorder` produces (``phase_mol``, ``pH``,
    ``ionic_strength``, ``P_atm``, ``t_h``).  List channels
    (``advance_results``, ``controller_actions``, etc.) are empty; only
    the numeric time-series channels are recovered from Parquet.

    Parameters
    ----------
    path : str or pathlib.Path
        Directory written by :class:`StreamingFileRecorder`.

    Returns
    -------
    BatchResult

    Raises
    ------
    FileNotFoundError
        If *path* does not contain a ``manifest.json``.
    ImportError
        If ``pyarrow`` is not installed.
    """
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError:
        raise ImportError(
            "pyarrow is required for load_run. "
            "Install it with: pip install pyarrow"
        )
    import json
    from pathlib import Path as _Path

    src = _Path(path)
    manifest_path = src / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"load_run: manifest.json not found in {str(src)!r}. "
            "Is this a StreamingFileRecorder output directory?"
        )

    manifest = json.loads(manifest_path.read_text())
    chunks = manifest.get("chunks", [])

    result = BatchResult()
    if not chunks:
        return result

    tables = [pq.read_table(src / c["file"]) for c in chunks]
    full_table = pa.concat_tables(tables)
    data = full_table.to_pydict()

    result.t_h = np.array(data["t_h"])
    n = len(result.t_h)

    cv_phase_species: Dict[str, Dict[str, List[str]]] = manifest.get(
        "cv_phase_species", {}
    )
    for cv_key, phases in cv_phase_species.items():
        result.phase_mol[cv_key] = {}
        result.pH[cv_key] = np.full(n, np.nan)
        result.ionic_strength[cv_key] = np.full(n, np.nan)
        result.P_atm[cv_key] = np.zeros(n)

        for phase_key, species_list in phases.items():
            result.phase_mol[cv_key][phase_key] = {}
            for sp in species_list:
                col = f"{cv_key}|{phase_key}|{sp}"
                vals = data.get(col, [0.0] * n)
                result.phase_mol[cv_key][phase_key][sp] = np.array(vals, dtype=float)

        for ch, attr in (
            ("pH", "pH"),
            ("ionic_strength", "ionic_strength"),
            ("P_atm", "P_atm"),
        ):
            col = f"{cv_key}|{ch}"
            vals = data.get(col)
            if vals is not None:
                getattr(result, attr)[cv_key] = np.array(vals, dtype=float)

    result.runtime_s = float(manifest.get("runtime_s", 0.0))
    return result


# ════════════════════════════════════════════════════════════════════════
#  SparseRecorder
# ════════════════════════════════════════════════════════════════════════


class SparseRecorder:
    """Recorder that captures every M-th step plus the initial and final states.

    Reduces in-memory footprint and post-run analysis cost for runs where
    intermediate detail is not needed.  The returned :class:`BatchResult`
    contains only the recorded time points; ``t_h`` and all per-species
    arrays have length equal to the number of recorded steps.

    Parameters
    ----------
    every_m_steps : int
        Recording interval.  Steps where ``step_index % every_m_steps == 0``
        are recorded, plus always the initial (step 0) and final states.
        Default 10.
    """

    def __init__(self, every_m_steps: int = 10) -> None:
        self._m: int = int(every_m_steps)
        self._n_steps: int = 0
        self._rows: List[Dict[str, Any]] = []
        self._cv_phase_species: Dict[str, Dict[str, List[str]]] = {}
        self._initialised: bool = False
        self._finalized: bool = False
        self._start_time: float = 0.0
        self._sim: Any = None

    # ── Protocol ──────────────────────────────────────────────────────

    def record_init(self, sim: Any, n_steps: int, start_t_h: float = 0.0) -> None:
        if self._initialised:
            raise RuntimeError(
                "SparseRecorder.record_init called twice. "
                "Construct a fresh recorder per Simulation.run() call."
            )
        self._sim = sim
        self._n_steps = int(n_steps)

        for cv_key, cv in sim.cvs.items():
            self._cv_phase_species[cv_key] = {
                pk: list(ph.n_mol.keys()) for pk, ph in cv.phases.items()
            }

        self._rows.append(self._capture_row(0, float(start_t_h)))
        self._start_time = time.perf_counter()
        self._initialised = True

    def record_step(
        self,
        step_index: int,
        t_h: float,
        dt_h: float,
        advance_results: Dict[str, AdvanceResult],
        link_records: Optional[List[LinkFlowRecord]] = None,
        controller_actions: Optional[List[ControlAction]] = None,
        profile_actions: Optional[List[ProfileRecord]] = None,
    ) -> None:
        if not self._initialised:
            raise RuntimeError(
                "SparseRecorder.record_step called before record_init."
            )
        if self._finalized:
            raise RuntimeError(
                "SparseRecorder.record_step called after finalize."
            )
        if step_index % self._m == 0 or step_index == self._n_steps:
            self._rows.append(self._capture_row(step_index, float(t_h)))

    def finalize(self) -> BatchResult:
        if not self._initialised:
            raise RuntimeError(
                "SparseRecorder.finalize called before record_init."
            )
        if self._finalized:
            raise RuntimeError(
                "SparseRecorder.finalize called twice."
            )
        self._finalized = True
        return self._build_result()

    # ── Internal ──────────────────────────────────────────────────────

    def _capture_row(self, step_index: int, t_h: float) -> Dict[str, Any]:
        row: Dict[str, Any] = {"step_index": step_index, "t_h": t_h}
        for cv_key, cv in self._sim.cvs.items():
            for phase_key, phase in cv.phases.items():
                for sp in self._cv_phase_species.get(cv_key, {}).get(phase_key, []):
                    row[f"{cv_key}|{phase_key}|{sp}"] = float(phase.n_mol.get(sp, 0.0))

            gas = cv.phases.get("gas")
            row[f"{cv_key}|P_atm"] = float(gas.P_atm) if gas is not None else 0.0

            liquid = cv.phases.get("liquid")
            if liquid is not None:
                ph_val = float("nan")
                if "H+" in liquid.n_mol:
                    try:
                        ph_val = float(liquid.pH)
                    except (ValueError, AttributeError):
                        pass
                row[f"{cv_key}|pH"] = ph_val
                is_val = liquid.properties.get("ionic_strength")
                if is_val is None:
                    is_val = liquid.speciation.get("IonicStrength")
                row[f"{cv_key}|ionic_strength"] = (
                    float(is_val) if is_val is not None else float("nan")
                )
            else:
                row[f"{cv_key}|pH"] = float("nan")
                row[f"{cv_key}|ionic_strength"] = float("nan")
        return row

    def _build_result(self) -> BatchResult:
        n = len(self._rows)
        result = BatchResult()
        result.t_h = np.array([r["t_h"] for r in self._rows])
        result.runtime_s = time.perf_counter() - self._start_time

        for cv_key, phases in self._cv_phase_species.items():
            result.phase_mol[cv_key] = {}
            result.pH[cv_key] = np.array(
                [r[f"{cv_key}|pH"] for r in self._rows]
            )
            result.ionic_strength[cv_key] = np.array(
                [r[f"{cv_key}|ionic_strength"] for r in self._rows]
            )
            result.P_atm[cv_key] = np.array(
                [r[f"{cv_key}|P_atm"] for r in self._rows]
            )
            for phase_key, species_list in phases.items():
                result.phase_mol[cv_key][phase_key] = {}
                for sp in species_list:
                    col = f"{cv_key}|{phase_key}|{sp}"
                    result.phase_mol[cv_key][phase_key][sp] = np.array(
                        [r.get(col, 0.0) for r in self._rows]
                    )
        return result


# ════════════════════════════════════════════════════════════════════════
#  SummaryResult + SummaryRecorder
# ════════════════════════════════════════════════════════════════════════


@dataclass
class SummaryResult:
    """Result produced by :class:`SummaryRecorder`.

    Holds the initial and final CV state plus cumulative boundary flux totals
    for the run.  Suitable for parameter-sweep / Monte-Carlo studies where
    only endpoints matter and storing thousands of full trajectories is
    wasteful.

    Attributes
    ----------
    start_t_h : float
        Simulation time at the start of the run.
    end_t_h : float
        Simulation time at the end of the run.
    n_steps : int
        Number of steps executed.
    runtime_s : float
        Wall-clock duration.
    initial_phase_mol : dict
        ``{cv_key: {phase_key: {species: float}}}`` — moles at t=0.
    final_phase_mol : dict
        ``{cv_key: {phase_key: {species: float}}}`` — moles at t=end.
    boundary_totals : dict
        ``{cv_key: {boundary_label: {species: float}}}`` — cumulative
        ``mol_applied`` summed over all steps from every
        :class:`~PyOMES.core.boundaries.ExternalFluxRecord`.
    """

    start_t_h: float = 0.0
    end_t_h: float = 0.0
    n_steps: int = 0
    runtime_s: float = 0.0
    initial_phase_mol: Dict[str, Dict[str, Dict[str, float]]] = field(
        default_factory=dict
    )
    final_phase_mol: Dict[str, Dict[str, Dict[str, float]]] = field(
        default_factory=dict
    )
    boundary_totals: Dict[str, Dict[str, Dict[str, float]]] = field(
        default_factory=dict
    )


class SummaryRecorder:
    """Recorder that captures only initial/final state plus boundary totals.

    Stores the CV state at ``t=0`` and at the final step, and accumulates
    cumulative ``mol_applied`` from every
    :class:`~PyOMES.core.boundaries.ExternalFluxRecord` across all steps.
    Memory usage is O(1) in the number of steps.

    Returns a :class:`SummaryResult` from :meth:`finalize` (not a
    :class:`BatchResult`).  Use this recorder when trajectory detail is
    not needed — e.g. parameter sweeps where only endpoint moles matter.
    """

    def __init__(self) -> None:
        self._initialised: bool = False
        self._finalized: bool = False
        self._start_time: float = 0.0
        self._n_steps: int = 0
        self._sim: Any = None
        self._result: SummaryResult = SummaryResult()

    # ── Protocol ──────────────────────────────────────────────────────

    def record_init(self, sim: Any, n_steps: int, start_t_h: float = 0.0) -> None:
        if self._initialised:
            raise RuntimeError(
                "SummaryRecorder.record_init called twice. "
                "Construct a fresh recorder per Simulation.run() call."
            )
        self._sim = sim
        self._n_steps = int(n_steps)
        self._result.start_t_h = float(start_t_h)
        self._result.initial_phase_mol = self._snapshot_phase_mol(sim)
        self._start_time = time.perf_counter()
        self._initialised = True

    def record_step(
        self,
        step_index: int,
        t_h: float,
        dt_h: float,
        advance_results: Dict[str, AdvanceResult],
        link_records: Optional[List[LinkFlowRecord]] = None,
        controller_actions: Optional[List[ControlAction]] = None,
        profile_actions: Optional[List[ProfileRecord]] = None,
    ) -> None:
        if not self._initialised:
            raise RuntimeError(
                "SummaryRecorder.record_step called before record_init."
            )
        if self._finalized:
            raise RuntimeError(
                "SummaryRecorder.record_step called after finalize."
            )
        # Accumulate boundary flux totals from this step.
        for cv_key, ar in advance_results.items():
            cv_totals = self._result.boundary_totals.setdefault(cv_key, {})
            for rec in ar.boundary_records:
                label_totals = cv_totals.setdefault(rec.boundary_label, {})
                for sp, mol in rec.mol_applied.items():
                    label_totals[sp] = label_totals.get(sp, 0.0) + float(mol)

        # Always overwrite with the latest CV state — finalize will keep the last.
        self._result.end_t_h = float(t_h)
        self._result.n_steps = step_index

    def finalize(self) -> SummaryResult:
        if not self._initialised:
            raise RuntimeError(
                "SummaryRecorder.finalize called before record_init."
            )
        if self._finalized:
            raise RuntimeError(
                "SummaryRecorder.finalize called twice."
            )
        self._result.final_phase_mol = self._snapshot_phase_mol(self._sim)
        self._result.runtime_s = time.perf_counter() - self._start_time
        self._finalized = True
        return self._result

    # ── Internal ──────────────────────────────────────────────────────

    @staticmethod
    def _snapshot_phase_mol(sim: Any) -> Dict[str, Dict[str, Dict[str, float]]]:
        snap: Dict[str, Dict[str, Dict[str, float]]] = {}
        for cv_key, cv in sim.cvs.items():
            snap[cv_key] = {}
            for phase_key, phase in cv.phases.items():
                snap[cv_key][phase_key] = {
                    sp: float(n) for sp, n in phase.n_mol.items()
                }
        return snap
