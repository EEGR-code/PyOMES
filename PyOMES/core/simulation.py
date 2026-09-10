# -*- coding: utf-8 -*-
"""Simulation orchestrator.

Single orchestration pathway for all CV-based models. Subsumes the
legacy :class:`~PyOMES.core.multi_cv.MultiCVSystem` (deleted at C14)
and the legacy ``run_batch`` time-loop body. Controllers and profiles
are first-class members of the :class:`Simulation`; lifecycle gating
is implemented through a shared :class:`RunContext` that every
lockable owned object holds a reference to.

See [`docs/phases-upcoming/SIMULATION_CLASS.md`](../../docs/phases-upcoming/SIMULATION_CLASS.md)
for the design rationale and
[`docs/phases-upcoming/SIMULATION_CLASS_CHECKLIST.md`](../../docs/phases-upcoming/SIMULATION_CLASS_CHECKLIST.md)
for the implementation checkpoint plan.

This module currently lands C1 of the phase: skeleton + RunContext
wiring. ``.run()`` and ``_step()`` land at C4.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from .control_volume import ControlVolume
from .interfaces import AdvanceResult
from .lifecycle import RunContext, _LockableList, raise_if_running
from .links import CVLink, LinkFlowRecord
from .recorder import BatchRecorder
from .snapshot import CVSnapshot, build_cv_snapshot, build_simulation_snapshot
from ..control.actions import ControlAction, ProfileRecord


class Simulation:
    """Single orchestration class for all CV-based models.

    Parameters
    ----------
    cvs : dict
        ``{key: ControlVolume}`` — the named CVs in this simulation.
        Each CV must not already be owned by another Simulation
        (use ``cv.snapshot()`` to fork a CV across simulations).
    links : list of CVLink, optional
        Inter-CV transport links. Endpoints (``source_cv_key``,
        ``sink_cv_key``, ``source_phase_key``, ``sink_phase_key``)
        are validated at construction.
    controllers : list, optional
        Controller protocol implementations. Consumed at C9. Today
        stored as-is; no run-time effect until ``.run()`` lands at C4.
    profiles : list, optional
        Profile protocol implementations. Consumed at C10. Today
        stored as-is.
    solver : StepSolver, dict, or future MultiCVStepSolver, optional
        Step-solver configuration. ``None`` defaults to the
        sequential ``cv.advance()`` body. A dict ``{cv_key: solver}``
        routes per-CV. A ``MultiCVStepSolver`` (placeholder type)
        targets a future system-level integrator — not implemented
        in this phase; passing one raises ``NotImplementedError``
        at run-time.
    recorder : Recorder, optional
        Pluggable recorder for time-series capture. Default
        :class:`~PyOMES.core.recorder.BatchRecorder` (lands at C3)
        is constructed at ``.run()`` entry if ``None`` is passed.
    label : str
        Human-readable label.

    Notes
    -----
    Construction wires every CV to a freshly-created
    :class:`RunContext`. The wiring is one-way (Sim → CV
    via ``cv._context``); CVs never hold a back-reference to
    the Simulation.
    """

    def __init__(
        self,
        cvs: Dict[str, ControlVolume],
        links: Optional[List[CVLink]] = None,
        controllers: Optional[List[Any]] = None,
        profiles: Optional[List[Any]] = None,
        solver: Optional[Any] = None,
        recorder: Optional[Any] = None,
        label: str = "",
        system_solver: Optional[Any] = None,
    ):
        self.cvs: Dict[str, ControlVolume] = dict(cvs)
        self.links: List[CVLink] = list(links or [])
        # Lockable wrappers for the list-shaped configuration surfaces
        # (decision 9: controllers / profiles cannot be appended /
        # popped mid-run). Storage is exposed via gated properties
        # below so reassignment is also blocked.
        self._controllers = _LockableList(
            controllers or [], label="Simulation.controllers"
        )
        self._profiles = _LockableList(
            profiles or [], label="Simulation.profiles"
        )
        self._solver: Any = solver
        self._recorder: Any = recorder
        self._system_solver: Any = system_solver
        self.label: str = str(label)

        # Wall-clock time accumulator. Simulation owns this; passed
        # directly to cv.advance(dt_h, t_h) each step. Initialised
        # to 0; .run() resets to 0 at entry (decision 10, always
        # destructive). Updated by _step at C4.
        self._t_h: float = 0.0

        # Pre-flight validation: link endpoints reference existing
        # CV keys and phase keys (raise before any state mutation).
        for link in self.links:
            self._validate_link_endpoints(link)

        # chemistry-unification-3b C10: warn when a link connects CVs
        # with different ThermoFramework instances.  Expected for CVs
        # representing physically distinct phases (e.g. aqueous + organic);
        # potentially a bug if they represent zones of the same phase.
        self._warn_thermo_mismatch()

        # Pre-flight validation: no CV is already owned by another
        # Simulation. Raise before wiring so a half-owned partial
        # Simulation never escapes.
        for key, cv in self.cvs.items():
            if cv._context is not None:
                raise RuntimeError(
                    f"ControlVolume {key!r} is already owned by another "
                    f"Simulation (context label={cv._context.label!r}). "
                    f"Each CV can belong to only one Simulation at a time; "
                    f"use cv.snapshot() to fork."
                )

        # Create the shared RunContext and wire every owned object.
        self._context: RunContext = RunContext(
            is_running=False,
            label=self.label,
        )
        # Lockable lists owned by this Simulation directly.
        self._controllers._context = self._context
        self._profiles._context = self._context
        # CVs, their phases, internal interfaces, and lockable lists.
        for cv in self.cvs.values():
            self._wire_cv_context(cv)

    def _wire_cv_context(self, cv: ControlVolume) -> None:
        """Propagate ``self._context`` to a CV and every lockable
        object it owns.

        Called for every CV at ``Simulation.__init__``. Handles:

        * ``cv._context`` itself (set previously in C1; reassigned
          here for clarity)
        * every Phase in ``cv.phases``
        * every PhaseInterface in ``cv.internal_interfaces`` that
          carries a ``_context`` attribute (e.g.
          :class:`KineticGasLiquidLink`)
        * ``cv.boundaries`` (a ``_LockableList`` after C6)
        * ``cv.property_calculators`` (a ``_LockableList`` after C6)
        """
        cv._context = self._context
        for phase in cv.phases.values():
            phase._context = self._context
        for iface in cv.internal_interfaces:
            if hasattr(iface, "_context"):
                iface._context = self._context
        if hasattr(cv.boundaries, "_context"):
            cv.boundaries._context = self._context
        for boundary in cv.boundaries:
            if hasattr(boundary, "_context"):
                boundary._context = self._context
        if hasattr(cv.property_calculators, "_context"):
            cv.property_calculators._context = self._context

    def _validate_link_endpoints(self, link: CVLink) -> None:
        """Carry the link-endpoint validation from MultiCVSystem verbatim."""
        src_key = link.source_cv_key
        snk_key = link.sink_cv_key
        if src_key not in self.cvs:
            raise KeyError(
                f"Link {link.label!r} references source_cv_key={src_key!r}, "
                f"but Simulation only contains: {sorted(self.cvs.keys())}"
            )
        if snk_key not in self.cvs:
            raise KeyError(
                f"Link {link.label!r} references sink_cv_key={snk_key!r}, "
                f"but Simulation only contains: {sorted(self.cvs.keys())}"
            )
        src_cv = self.cvs[src_key]
        snk_cv = self.cvs[snk_key]
        if link.source_phase_key not in src_cv.phases:
            raise KeyError(
                f"Link {link.label!r} references source phase "
                f"{link.source_phase_key!r} in CV {src_key!r}, "
                f"but that CV only has phases: {src_cv.phase_keys}"
            )
        if link.sink_phase_key not in snk_cv.phases:
            raise KeyError(
                f"Link {link.label!r} references sink phase "
                f"{link.sink_phase_key!r} in CV {snk_key!r}, "
                f"but that CV only has phases: {snk_cv.phase_keys}"
            )

    def _warn_thermo_mismatch(self) -> None:
        """Warn if any link connects CVs with different ThermoFramework instances."""
        import warnings
        for link in self.links:
            src_cv = self.cvs.get(link.source_cv_key)
            snk_cv = self.cvs.get(link.sink_cv_key)
            if src_cv is None or snk_cv is None:
                continue
            src_db = getattr(src_cv, "chemistry_db", None)
            snk_db = getattr(snk_cv, "chemistry_db", None)
            if src_db is None or snk_db is None:
                continue
            if src_db.thermo is not snk_db.thermo and src_db.thermo != snk_db.thermo:
                warnings.warn(
                    f"Link {getattr(link, 'label', repr(link))!r} connects "
                    f"CVs with different ThermoFrameworks. "
                    f"Expected if the CVs represent physically distinct phases "
                    f"(e.g. aqueous + organic); potentially a bug if they "
                    f"represent zones of the same phase.",
                    UserWarning,
                    stacklevel=3,
                )

    # ── Accessors ──────────────────────────────────────────────────────

    def __getitem__(self, key: str) -> ControlVolume:
        return self.cvs[key]

    def __contains__(self, key: str) -> bool:
        return key in self.cvs

    @property
    def cv_keys(self) -> List[str]:
        return list(self.cvs.keys())

    # ── Gated configuration surfaces (C6) ──────────────────────────────

    @property
    def controllers(self) -> _LockableList:
        return self._controllers

    @controllers.setter
    def controllers(self, value) -> None:
        raise_if_running(self, "controllers")
        self._set_controllers_unchecked(value)

    def _set_controllers_unchecked(self, value) -> None:
        new_list = _LockableList(
            value or [], label="Simulation.controllers"
        )
        new_list._context = self._context
        self._controllers = new_list

    @property
    def profiles(self) -> _LockableList:
        return self._profiles

    @profiles.setter
    def profiles(self, value) -> None:
        raise_if_running(self, "profiles")
        self._set_profiles_unchecked(value)

    def _set_profiles_unchecked(self, value) -> None:
        new_list = _LockableList(value or [], label="Simulation.profiles")
        new_list._context = self._context
        self._profiles = new_list

    @property
    def solver(self) -> Any:
        return self._solver

    @solver.setter
    def solver(self, value: Any) -> None:
        raise_if_running(self, "solver")
        self._set_solver_unchecked(value)

    def _set_solver_unchecked(self, value: Any) -> None:
        self._solver = value

    @property
    def system_solver(self) -> Any:
        return self._system_solver

    @system_solver.setter
    def system_solver(self, value: Any) -> None:
        raise_if_running(self, "system_solver")
        self._set_system_solver_unchecked(value)

    def _set_system_solver_unchecked(self, value: Any) -> None:
        self._system_solver = value

    @property
    def recorder(self) -> Any:
        return self._recorder

    @recorder.setter
    def recorder(self, value: Any) -> None:
        raise_if_running(self, "recorder")
        self._set_recorder_unchecked(value)

    def _set_recorder_unchecked(self, value: Any) -> None:
        self._recorder = value

    # ── System-wide inventory ──────────────────────────────────────────

    def total_mol(self) -> Dict[str, float]:
        """Sum moles of each species across every CV and every phase.

        Returns ``{species_id: total_mol}`` aggregated over the entire
        system. Carries the semantic of
        :meth:`~PyOMES.core.multi_cv.MultiCVSystem.total_mol`.
        """
        totals: Dict[str, float] = {}
        for cv in self.cvs.values():
            for species, n in cv.total_mol().items():
                totals[species] = totals.get(species, 0.0) + float(n)
        return totals

    # ── Snapshot ───────────────────────────────────────────────────────

    def snapshot(self) -> "Simulation":
        """Return an independent simulation with deep-copied CVs.

        Each CV is deep-copied via ``cv.snapshot()`` (which gives the
        copy ``_context = None``). The new Simulation creates its own
        :class:`RunContext` in ``__init__`` and wires the snapshotted
        CVs. Links, controllers, profiles, solver, and recorder are
        *shared* (not copied) — they are configuration/behaviour
        objects, not state.
        """
        return Simulation(
            cvs={k: cv.snapshot() for k, cv in self.cvs.items()},
            links=list(self.links),
            controllers=list(self.controllers),
            profiles=list(self.profiles),
            solver=self.solver,
            recorder=self._recorder,
            label=self.label,
            system_solver=self._system_solver,
        )

    # ── Checkpoint ────────────────────────────────────────────────────

    @property
    def checkpoint_t_h(self) -> Optional[float]:
        """Wall-clock time (hours) recorded by the most recent checkpoint.

        ``None`` on a freshly constructed :class:`Simulation`. Set to
        the manifest's ``t_h`` value by :meth:`load_checkpoint`. Use
        directly as ``start_t_h`` when resuming::

            sim = Simulation.load_checkpoint("run1/")
            r2 = sim.run(tau_h=5.0, n_steps=500,
                         start_t_h=sim.checkpoint_t_h)
        """
        return getattr(self, "_checkpoint_t_h", None)

    def save_checkpoint(
        self,
        path: Any,
        *,
        mode: str = "data+env",
        compression: str = "gzip",
    ) -> None:
        """Save the simulation state to a checkpoint directory.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination directory. Created (with parents) if absent.
        mode : {"data", "data+env", "source", "full"}
            Fidelity level. ``"data"`` pickles the full object graph
            and writes a manifest with PyOMES version, timestamp, and
            BLAS info. ``"data+env"`` additionally captures
            ``pip freeze``, Python version, platform, and
            OMP_NUM_THREADS in an ``environment/`` subdirectory.
            ``"source"`` and ``"full"`` raise
            :exc:`NotImplementedError` in this iteration.
        compression : {"none", "gzip", "zstd"}
            Compression applied to ``simulation.pkl``. ``"gzip"``
            uses the standard library; ``"zstd"`` requires
            ``zstandard`` (``pip install zstandard``); ``"none"``
            writes raw pickle bytes.

        Notes
        -----
        Call between :meth:`run` invocations, not mid-step. The
        checkpoint captures the complete object graph including
        controllers, links, profiles, solver, and recorder. Objects
        with external state (open file handles, database connections)
        need ``__getstate__``/``__setstate__``; typical pure-Python
        PyOMES controllers never require this.
        """
        import datetime
        import gzip
        import json
        import pickle
        from pathlib import Path

        _VALID_MODES = ("data", "data+env", "source", "full")
        _VALID_COMPRESSIONS = ("none", "gzip", "zstd")

        if mode not in _VALID_MODES:
            raise ValueError(
                f"save_checkpoint: mode must be one of {_VALID_MODES!r}; "
                f"got {mode!r}."
            )
        if mode in ("source", "full"):
            raise NotImplementedError(
                f"save_checkpoint mode={mode!r} is not yet implemented. "
                f"Only 'data' and 'data+env' are available in this iteration."
            )
        if compression not in _VALID_COMPRESSIONS:
            raise ValueError(
                f"save_checkpoint: compression must be one of "
                f"{_VALID_COMPRESSIONS!r}; got {compression!r}."
            )

        dest = Path(path)
        dest.mkdir(parents=True, exist_ok=True)

        # Pickle the full object graph.
        raw = pickle.dumps(self, protocol=pickle.HIGHEST_PROTOCOL)
        if compression == "gzip":
            payload = gzip.compress(raw)
        elif compression == "zstd":
            try:
                import zstandard as zstd  # type: ignore[import]
            except ImportError:
                raise ImportError(
                    "zstandard is required for compression='zstd'. "
                    "Install with: pip install zstandard"
                )
            payload = zstd.ZstdCompressor().compress(raw)
        else:
            payload = raw

        (dest / "simulation.pkl").write_bytes(payload)

        # PyOMES version
        try:
            from importlib.metadata import version as _pkg_ver
            vlsim_version = _pkg_ver("PyOMES")
        except Exception:
            vlsim_version = "unknown"

        # BLAS info
        import numpy as np_mod
        try:
            blas_info = str(np_mod.__config__.blas_opt_info)
        except AttributeError:
            blas_info = "unknown"

        manifest = {
            "schema_version": 1,
            "mode": mode,
            "vlsim_version": vlsim_version,
            "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "t_h": float(self._t_h),
            "compression": compression,
            "blas_info": blas_info,
        }

        # If the recorder writes streaming chunks, record the chunk directory.
        # The buffer is excluded from the pickle via StreamingFileRecorder.__getstate__.
        from .recorder import StreamingFileRecorder as _SFR
        if isinstance(self._recorder, _SFR):
            manifest["streaming_recorder_path"] = str(self._recorder._path)

        (dest / "manifest.json").write_text(json.dumps(manifest, indent=2))

        # mode "data+env": capture runtime environment
        if mode == "data+env":
            import os
            import platform
            import subprocess
            import sys

            env_dir = dest / "environment"
            env_dir.mkdir(exist_ok=True)

            (env_dir / "python_version.txt").write_text(sys.version)
            (env_dir / "platform.txt").write_text(platform.platform())
            (env_dir / "omp_num_threads.txt").write_text(
                os.environ.get("OMP_NUM_THREADS", "1")
            )
            try:
                result = subprocess.run(
                    [sys.executable, "-m", "pip", "freeze"],
                    capture_output=True, text=True, timeout=30,
                )
                pip_freeze = result.stdout
            except Exception:
                pip_freeze = ""
            (env_dir / "pip_freeze.txt").write_text(pip_freeze)

    @classmethod
    def load_checkpoint(
        cls,
        path: Any,
        *,
        verify: str = "strict",
    ) -> "Simulation":
        """Load a :class:`Simulation` from a checkpoint directory.

        Parameters
        ----------
        path : str or pathlib.Path
            Directory written by :meth:`save_checkpoint`.
        verify : {"strict", "warn", "skip"}
            Drift-validation policy on load.

            ``"strict"`` (default) raises :exc:`RuntimeError` on
            PyOMES version mismatch. ``"warn"`` issues a
            :class:`UserWarning` instead. ``"skip"`` bypasses all
            checks.

            Python-version, third-party-package, platform,
            OMP-thread-count, and BLAS checks require mode
            ``"data+env"`` and are silently skipped for ``"data"``
            mode checkpoints.

        Returns
        -------
        Simulation
            The loaded simulation with :attr:`checkpoint_t_h` set
            to the wall-clock time recorded at save.
        """
        import gzip
        import json
        import pickle
        import warnings
        from pathlib import Path

        _VALID_VERIFY = ("strict", "warn", "skip")
        if verify not in _VALID_VERIFY:
            raise ValueError(
                f"load_checkpoint: verify must be one of {_VALID_VERIFY!r}; "
                f"got {verify!r}."
            )

        src = Path(path)
        manifest_path = src / "manifest.json"
        pkl_path = src / "simulation.pkl"

        if not manifest_path.exists():
            raise FileNotFoundError(
                f"load_checkpoint: manifest.json not found in {str(src)!r}."
            )
        if not pkl_path.exists():
            raise FileNotFoundError(
                f"load_checkpoint: simulation.pkl not found in {str(src)!r}."
            )

        manifest = json.loads(manifest_path.read_text())

        # PyOMES version check (available for all modes from the manifest).
        if verify != "skip":
            try:
                from importlib.metadata import version as _pkg_ver
                current_version = _pkg_ver("PyOMES")
            except Exception:
                current_version = "unknown"
            saved_version = manifest.get("vlsim_version", "unknown")
            if saved_version != current_version:
                msg = (
                    f"PyOMES version mismatch: checkpoint was saved with "
                    f"PyOMES {saved_version!r}, current is {current_version!r}. "
                    f"The pickle may fail or produce incorrect results."
                )
                if verify == "strict":
                    raise RuntimeError(msg)
                warnings.warn(msg, UserWarning, stacklevel=2)

        # mode "data+env": check runtime environment drift
        if manifest.get("mode") == "data+env" and verify != "skip":
            import os
            import platform
            import sys

            env_dir = src / "environment"

            # Python version — strict raises, warn warns
            py_path = env_dir / "python_version.txt"
            if py_path.exists():
                saved_py = py_path.read_text().strip()
                current_py = sys.version.strip()
                if saved_py != current_py:
                    msg = (
                        f"Python version mismatch: checkpoint used "
                        f"{saved_py!r}, current is {current_py!r}. "
                        f"The pickle may fail or produce incorrect results."
                    )
                    if verify == "strict":
                        raise RuntimeError(msg)
                    warnings.warn(msg, UserWarning, stacklevel=2)

            # Platform — warn only (different OS ≠ wrong results usually)
            plat_path = env_dir / "platform.txt"
            if plat_path.exists():
                saved_plat = plat_path.read_text().strip()
                current_plat = platform.platform().strip()
                if saved_plat != current_plat:
                    warnings.warn(
                        f"Platform mismatch: checkpoint was created on "
                        f"{saved_plat!r}, current is {current_plat!r}.",
                        UserWarning, stacklevel=2,
                    )

            # OMP_NUM_THREADS — warn only (affects BLAS summation order)
            omp_path = env_dir / "omp_num_threads.txt"
            if omp_path.exists():
                saved_omp = omp_path.read_text().strip()
                current_omp = os.environ.get("OMP_NUM_THREADS", "1").strip()
                if saved_omp != current_omp:
                    warnings.warn(
                        f"OMP_NUM_THREADS mismatch: checkpoint used "
                        f"{saved_omp!r}, current is {current_omp!r}. "
                        f"Threaded BLAS summation order may differ.",
                        UserWarning, stacklevel=2,
                    )

            # Third-party packages — warn only (set comparison avoids
            # false positives from pip freeze ordering changes)
            pip_path = env_dir / "pip_freeze.txt"
            if pip_path.exists():
                saved_pkgs = set(pip_path.read_text().strip().splitlines())
                try:
                    import subprocess
                    result = subprocess.run(
                        [sys.executable, "-m", "pip", "freeze"],
                        capture_output=True, text=True, timeout=30,
                    )
                    current_pkgs = set(result.stdout.strip().splitlines())
                except Exception:
                    current_pkgs = set()
                if saved_pkgs != current_pkgs:
                    warnings.warn(
                        "Third-party package versions differ between "
                        "checkpoint and current environment. Numerical "
                        "results may not be bit-for-bit reproducible.",
                        UserWarning, stacklevel=2,
                    )

            # BLAS — warn only
            saved_blas = manifest.get("blas_info", "unknown")
            import numpy as _np_ck
            try:
                current_blas = str(_np_ck.__config__.blas_opt_info)
            except AttributeError:
                current_blas = "unknown"
            if saved_blas != current_blas:
                warnings.warn(
                    "BLAS configuration mismatch: summation order may "
                    "differ. Use mode='full' for binary-hash-level "
                    "fidelity (not yet implemented).",
                    UserWarning, stacklevel=2,
                )

        # Load pickle with decompression matching the manifest.
        compression = manifest.get("compression", "none")
        raw_bytes = pkl_path.read_bytes()
        if compression == "gzip":
            raw = gzip.decompress(raw_bytes)
        elif compression == "zstd":
            try:
                import zstandard as zstd  # type: ignore[import]
            except ImportError:
                raise ImportError(
                    "zstandard is required to load this checkpoint. "
                    "Install with: pip install zstandard"
                )
            raw = zstd.ZstdDecompressor().decompress(raw_bytes)
        else:
            raw = raw_bytes

        sim = pickle.loads(raw)
        sim._checkpoint_t_h = float(manifest["t_h"])
        return sim

    # ── Run loop (C4) ──────────────────────────────────────────────────

    def run(self, tau_h: float, n_steps: int, *, start_t_h: float = 0.0) -> Any:
        """Advance every CV for ``tau_h`` hours in ``n_steps`` steps.

        Single-CV path at C4; multi-CV link handling lands at C5.
        Controller and profile invocation land at C9 and C10.

        Parameters
        ----------
        tau_h : float
            Total simulation duration (hours). Must be > 0.
        n_steps : int
            Number of timesteps. Must be >= 1. For
            :class:`SimultaneousEulerSolver` this is the number of
            explicit steps; for adaptive solvers (e.g.
            ``SimultaneousAdaptiveSolver``) this is the output grid resolution.
        start_t_h : float, optional
            Start of the output time grid (hours). Default ``0.0``
            preserves all existing behaviour (grid runs
            ``[0, tau_h]``). Pass a non-zero value to shift the grid
            to ``[start_t_h, start_t_h + tau_h]`` — used for
            checkpoint resume::

                sim2 = Simulation.load_checkpoint("run1/")
                r2 = sim2.run(tau_h=5.0, n_steps=500,
                              start_t_h=sim2.checkpoint_t_h)

        Returns
        -------
        Any
            Whatever the recorder's :meth:`finalize` returns. For
            the default :class:`BatchRecorder` this is a
            :class:`~PyOMES.core.recorder.BatchResult`.

        Notes
        -----
        Decision 10 (always destructive): the wall-clock accumulator
        ``self._t_h`` resets to ``start_t_h`` at run entry (default
        ``0.0``). This is intentional and not a gap for HPC checkpoint
        resume — the time grid is shifted externally via ``start_t_h``
        rather than by accumulating internal state across calls.

        Decision 11: ``Simulation`` owns ``t_h`` and threads it to
        ``cv.advance(dt_h, t_h)`` directly; no ``chem_env`` /
        ``chem_env_fn`` plumbing.

        The lifecycle flag ``self._context.is_running`` is flipped
        ``True`` at run entry and cleared in a ``try/finally`` at
        exit. C6 wires the lockable-mutator consumers; until then
        the flag is set but not consulted.
        """
        if int(n_steps) < 1:
            raise ValueError(
                f"Simulation.run: n_steps must be >= 1; got {n_steps!r}."
            )
        if float(tau_h) <= 0.0:
            raise ValueError(
                f"Simulation.run: tau_h must be > 0; got {tau_h!r}."
            )

        # Resolve recorder. If user passed one, use it. Otherwise
        # construct fresh per call — BatchRecorder is single-use.
        recorder = (
            self._recorder if self._recorder is not None else BatchRecorder()
        )

        # Initial speciation solve at t=0 so the recorder's
        # record_init captures pH[0] and ionic_strength[0].
        # Carries the equivalent from
        # models/vlmodels/fermenter/config/factory.py:run_batch
        # lines 542-557.
        for cv in self.cvs.values():
            self._initial_solve(cv)

        recorder.record_init(self, n_steps, start_t_h=float(start_t_h))

        # C9: reset every controller (decision 10: .run() always
        # destructive). Controllers without a reset() are tolerated.
        for ctrl in self._controllers:
            if hasattr(ctrl, "reset"):
                ctrl.reset()
        # Sample-period bookkeeping per controller (ZOH state).
        # Keyed by id(ctrl); reset every .run().
        self._controller_sample_state: Dict[int, Dict[str, Any]] = {}

        # Build the time grid (linspace keeps t[-1] exact).
        # start_t_h shifts the grid to [start_t_h, start_t_h + tau_h].
        _start = float(start_t_h)
        t = np.linspace(_start, _start + float(tau_h), int(n_steps) + 1)
        dt_h = float(tau_h) / int(n_steps)

        # Reset the wall-clock accumulator to start_t_h (decision 10:
        # always destructive). Flip the lifecycle flag.
        self._t_h = _start
        self._context.is_running = True
        try:
            for i in range(1, int(n_steps) + 1):
                t_step_end = float(t[i])
                self._t_h = t_step_end
                (
                    results, link_records,
                    controller_actions, profile_actions,
                ) = self._step(dt_h, t_step_end)
                recorder.record_step(
                    step_index=i,
                    t_h=t_step_end,
                    dt_h=dt_h,
                    advance_results=results,
                    link_records=link_records,
                    controller_actions=controller_actions,
                    profile_actions=profile_actions,
                )
        finally:
            self._context.is_running = False

        return recorder.finalize()

    def _step(
        self, dt_h: float, t_h: float,
    ) -> "tuple[Dict[str, AdvanceResult], List[LinkFlowRecord], List[ControlAction], List[ProfileRecord]]":
        """Advance the system by one timestep.

        Dispatches to the active :class:`~PyOMES.core.system_solver.SystemSolver`
        when ``system_solver`` is set; otherwise falls through to
        :meth:`_step_default` (the sequential explicit-Euler body).
        """
        if self._system_solver is not None:
            return self._system_solver.advance_system(self, dt_h, t_h)
        return self._step_default(dt_h, t_h)

    def _step_default(
        self, dt_h: float, t_h: float,
    ) -> "tuple[Dict[str, AdvanceResult], List[LinkFlowRecord], List[ControlAction], List[ProfileRecord]]":
        """Sequential explicit-Euler step (default Axis 2 behaviour).

        Order per step:

        1. Apply every profile (open-loop time-varying mutators).
           Profiles fire *before* advancement so the integration
           sees the updated T_K, vvm, setpoints, etc. Mutation is
           through Pattern B unchecked setters so the lifecycle
           gate isn't tripped.
        2. Apply every inter-CV link's flow symmetrically (source
           loses, sink gains) — yields one :class:`LinkFlowRecord`
           per link for diagnostics.
        3. Advance each CV individually with its dispatched solver.
        4. Build a SimulationSnapshot from the post-advance state.
        5. Invoke each controller (respecting its sample period),
           collecting a :class:`ControlAction` per firing.
        6. Apply each ControlAction's flux_applied (via
           cv.apply_external_flux) and params_changed (via
           the Pattern B unchecked-setter dispatch).

        Returns
        -------
        results : dict
            ``{cv_key: AdvanceResult}`` from each CV's advance call.
        link_records : list of LinkFlowRecord
            Inter-CV link flow diagnostics for this step (empty
            if no links).
        controller_actions : list of ControlAction
            Per-controller actions emitted this step.
        profile_actions : list of ProfileRecord
            Per-profile records emitted this step.
        """
        # 1. Apply profiles before advance.
        profile_actions = self._invoke_profiles(t_h=t_h)

        # 2. Apply inter-CV links before advancing CVs.
        link_records = self._apply_links(dt_h)

        # 3. Advance each CV with its dispatched solver. Trusted
        #    orchestrator path — bypasses the OrchestrationWarning
        #    ownership guard on the public cv.advance().
        results: Dict[str, AdvanceResult] = {}
        for cv_key, cv in self.cvs.items():
            solver = self._solver_for(cv_key)
            results[cv_key] = cv._advance_unchecked(dt_h, t_h, solver=solver)

        # 4-6. Controller invocation + action application.
        # Filter None (never-fired placeholder) before apply and recorder.
        controller_actions = [
            a for a in self._invoke_controllers(
                t_h=t_h, dt_h=dt_h, results=results,
            )
            if a is not None
        ]
        for action in controller_actions:
            self._apply_controller_action(action, dt_h=dt_h)

        return results, link_records, controller_actions, profile_actions

    def _invoke_profiles(self, t_h: float) -> List[ProfileRecord]:
        """Apply every profile and collect the resulting records.

        Profiles mutate state through Pattern B unchecked setters
        (so the lifecycle gate isn't tripped) and return a
        :class:`ProfileRecord` summarising the targets they
        touched. Empty profile list → empty record list.
        """
        if not self._profiles:
            return []
        records: List[ProfileRecord] = []
        for profile in self._profiles:
            records.append(profile.apply(float(t_h), self))
        return records

    # ── Controller invocation + dispatch (C9) ──────────────────────────

    def _invoke_controllers(
        self,
        t_h: float,
        dt_h: float,
        results: Dict[str, AdvanceResult],
    ) -> List[Optional[ControlAction]]:
        """Build a post-advance :class:`SimulationSnapshot`; invoke
        each controller (respecting sample-period gating); return
        the collected :class:`ControlAction` list.

        Sample-period filtering (zero-order hold): controllers with
        ``sample_period_h`` or ``sample_period_s`` only fire at
        sampling instants. Between samples, the orchestrator reissues
        the controller's last action (ZOH).

        Never-fired slots return ``None`` so recorder consumers can
        distinguish "ZOH re-issue" from "has never fired". The caller
        must filter ``None`` entries before applying actions or passing
        to the recorder.
        """
        if not self._controllers:
            return []
        sim_snap = build_simulation_snapshot(self, t_h=t_h, results=results)
        actions: List[Optional[ControlAction]] = []
        for ctrl in self._controllers:
            sample_state = self._controller_sample_state.setdefault(
                id(ctrl), {"last_fire_t_h": None, "last_action": None},
            )
            should_fire = self._should_fire(ctrl, t_h, sample_state)
            if should_fire:
                # Per-CV controllers get a CVSnapshot; multi-CV get
                # the SimulationSnapshot (controllers select via
                # target_cv_key internally).
                tcv = getattr(ctrl, "target_cv_key", None)
                if tcv is not None and tcv in sim_snap.cvs:
                    state = sim_snap.cvs[tcv]
                elif tcv is None and len(sim_snap.cvs) == 1:
                    # Single-CV simulation: hand the lone CVSnapshot.
                    state = next(iter(sim_snap.cvs.values()))
                else:
                    state = sim_snap
                action = ctrl.compute(state, dt_h)
                sample_state["last_fire_t_h"] = float(t_h)
                sample_state["last_action"] = action
            else:
                # Zero-order hold: re-issue the last action.  None means
                # the controller has never fired; callers skip None entries.
                action = sample_state["last_action"]
            actions.append(action)
        return actions

    @staticmethod
    def _should_fire(
        ctrl: Any, t_h: float, sample_state: Dict[str, Any],
    ) -> bool:
        """Decide whether ``ctrl`` should fire this step based on
        its (optional) ``sample_period_h`` / ``sample_period_s``."""
        period_h = getattr(ctrl, "sample_period_h", None)
        period_s = getattr(ctrl, "sample_period_s", None)
        if period_h is None and period_s is not None:
            period_h = float(period_s) / 3600.0
        if period_h is None or float(period_h) <= 0.0:
            return True  # No sampling period → fire every step.
        last = sample_state["last_fire_t_h"]
        if last is None:
            return True  # First call always fires.
        # Fire when at least one full period has elapsed (with a
        # tiny epsilon to absorb floating-point grid noise).
        return (float(t_h) - float(last)) >= (float(period_h) - 1e-12)

    def _apply_controller_action(
        self, action: ControlAction, dt_h: float,
    ) -> None:
        """Apply a single controller action to the target CV.

        - ``flux_applied`` → :meth:`ControlVolume.apply_external_flux`
          per phase (already orchestrator-public; not gated).
        - ``params_changed`` → :meth:`_apply_param_change` dispatch
          via the path-resolver (Pattern B; unchecked setters).
        """
        if not action.target_cv_key:
            return
        cv = self.cvs.get(action.target_cv_key)
        if cv is None:
            return  # Stale target_cv_key (should never happen in v1).

        # 1. flux_applied: apply per phase.
        for phase_key, species_flux in (action.flux_applied or {}).items():
            if not species_flux:
                continue
            if phase_key in cv.phases:
                cv.apply_external_flux(phase_key, dict(species_flux), dt_h)

        # 2. params_changed: dispatch via path resolver.
        for path, value in (action.params_changed or {}).items():
            self._apply_param_change(cv, path, value)

    def _apply_param_change(
        self, cv: ControlVolume, path: Any, value: Any,
    ) -> None:
        """Dispatch a ``params_changed`` entry onto a CV via the reflective
        walker (PARAM_PATH_DISPATCHER phase, C3).

        Accepts three path forms (Q1, Q5):

        * **string** — full reflective path, e.g.
          ``"internal_interfaces[KineticGasLiquidLink].kLa.O2"``.
          Parsed to a :class:`~PyOMES.control.param_path.ParamPath` and walked.
        * :class:`~PyOMES.control.param_path.ParamPath` — pre-parsed; walked
          directly.
        * **typed descriptor handle** — a
          :class:`~PyOMES.control.descriptors._DictItemRef`,
          :class:`~PyOMES.control.descriptors.MutableScalar`, or
          :class:`~PyOMES.control.descriptors.MutableDict` obtained at
          class level (e.g. ``KineticGasLiquidLink.kLa["O2"]``).

        Unknown or malformed paths raise
        :exc:`~PyOMES.control.param_path.ParamPathError` (Q3).
        """
        from ..control.param_path import ParamPath, ParamPathError
        from ..control.descriptors import MutableScalar, MutableDict, _DictItemRef

        if isinstance(path, str):
            self._walk_and_write(cv, ParamPath.parse(path), value)
        elif isinstance(path, ParamPath):
            self._walk_and_write(cv, path, value)
        elif isinstance(path, _DictItemRef):
            owner_cls = path._owner
            target = self._find_owner_in_cv(cv, owner_cls)
            if target is None:
                raise ParamPathError(
                    f"params_changed path {path!r}: no "
                    f"{owner_cls.__name__} found in CV {cv.label!r}."
                )
            path._descriptor._set_item_unchecked(target, path._key, value)
        elif isinstance(path, MutableScalar):
            owner_cls = getattr(path, "_owner", None)
            if owner_cls is None:
                raise ParamPathError(
                    f"MutableScalar path has no _owner; only class-level "
                    f"descriptor references are valid as path keys."
                )
            target = self._find_owner_in_cv(cv, owner_cls)
            if target is None:
                raise ParamPathError(
                    f"params_changed path for {path.name!r}: no "
                    f"{owner_cls.__name__} found in CV {cv.label!r}."
                )
            path._set_unchecked(target, value)
        elif isinstance(path, MutableDict):
            owner_cls = getattr(path, "_owner", None)
            if owner_cls is None:
                raise ParamPathError(
                    f"MutableDict path has no _owner."
                )
            target = self._find_owner_in_cv(cv, owner_cls)
            if target is None:
                raise ParamPathError(
                    f"params_changed path for {path.name!r}: no "
                    f"{owner_cls.__name__} found in CV {cv.label!r}."
                )
            path._set_unchecked(target, value)
        else:
            raise ParamPathError(
                f"params_changed key must be a string, ParamPath, or "
                f"typed descriptor handle; got {type(path).__name__!r}."
            )

    def _walk_and_write(
        self, cv: ControlVolume, path: Any, value: Any,
    ) -> None:
        """Walk *path* from *cv* and write *value* at the leaf."""
        from ..control.param_path import (
            ParamPath, ParamPathError,
            AttrSegment, ListSelector, IndexSelector,
            _get_class_registry,
        )

        registry = _get_class_registry()
        target: Any = cv

        for segment in path.head:
            if isinstance(segment, AttrSegment):
                target = segment.resolve(target)
            elif isinstance(segment, ListSelector):
                target = segment.resolve(target, registry)
            else:  # IndexSelector
                target = segment.resolve(target)

        leaf = path.leaf
        if not isinstance(leaf, AttrSegment):
            raise ParamPathError(
                f"Path {path!r} has a non-attribute leaf "
                f"({type(leaf).__name__}); the last segment must be "
                f"an identifier to write to."
            )
        self._write_leaf(target, leaf.name, value, path)

    def _write_leaf(
        self, target: Any, attr_name: str, value: Any, path: Any,
    ) -> None:
        """Apply the leaf write on *target*.

        Priority:
        1. Target is a dict → direct key write (e.g. ``link.kLa["O2"]``).
        2. Attribute is a :class:`~PyOMES.control.descriptors.MutableScalar`
           or :class:`~PyOMES.control.descriptors.MutableDict` → call
           ``_set_unchecked``.
        3. Pre-migration fallback: ``_set_{attr_name}_unchecked(value)``
           exists on *target* → call it.
        4. No writable path found → :exc:`~PyOMES.control.param_path.ParamPathError`.
        """
        from ..control.param_path import ParamPathError
        from ..control.descriptors import MutableScalar, MutableDict

        if isinstance(target, dict):
            target[attr_name] = value
            return

        descriptor = None
        for cls in type(target).__mro__:
            if attr_name in cls.__dict__:
                descriptor = cls.__dict__[attr_name]
                break

        if isinstance(descriptor, MutableScalar):
            descriptor._set_unchecked(target, value)
            return

        if isinstance(descriptor, MutableDict):
            descriptor._set_unchecked(target, value)
            return

        unchecked_fn = getattr(target, f"_set_{attr_name}_unchecked", None)
        if unchecked_fn is not None:
            unchecked_fn(float(value))
            return

        raise ParamPathError(
            f"No writable path for {attr_name!r} on "
            f"{type(target).__name__!r} in path {path!r}. "
            f"Declare {attr_name!r} with MutableScalar or MutableDict, "
            f"or provide _set_{attr_name}_unchecked."
        )

    def _find_owner_in_cv(
        self, cv: ControlVolume, owner_cls: type,
    ) -> Optional[Any]:
        """Search standard CV collections for the first instance of *owner_cls*.

        Searches ``cv.internal_interfaces``, ``cv.boundaries``,
        and ``cv.phases.values()`` in that order.
        """
        for iface in cv.internal_interfaces:
            if isinstance(iface, owner_cls):
                return iface
        for boundary in cv.boundaries:
            if isinstance(boundary, owner_cls):
                return boundary
        for phase in cv.phases.values():
            if isinstance(phase, owner_cls):
                return phase
        return None

    def _apply_links(self, dt_h: float) -> List[LinkFlowRecord]:
        """Compute and apply all inter-CV link flows.

        Carried from
        :meth:`~PyOMES.core.multi_cv.MultiCVSystem._apply_links`
        verbatim. Returns a per-link :class:`LinkFlowRecord` list
        for the recorder; an empty list when ``self.links`` is empty
        (the common single-CV path).
        """
        records: List[LinkFlowRecord] = []
        for link in self.links:
            flow = link.compute_flow(self.cvs, dt_h)
            if flow:
                src_cv = self.cvs[link.source_cv_key]
                snk_cv = self.cvs[link.sink_cv_key]
                neg_flow = {sp: -rate for sp, rate in flow.items()}
                src_cv.phases[link.source_phase_key].apply_flux(neg_flow, dt_h)
                snk_cv.phases[link.sink_phase_key].apply_flux(flow, dt_h)
            total_transferred = (
                sum(abs(v) * float(dt_h) for v in flow.values()) if flow else 0.0
            )
            records.append(
                LinkFlowRecord(
                    link_label=link.label,
                    source=f"{link.source_cv_key}.{link.source_phase_key}",
                    sink=f"{link.sink_cv_key}.{link.sink_phase_key}",
                    flow_mol_per_h=dict(flow),
                    total_mol_transferred=total_transferred,
                )
            )
        return records

    def _solver_for(self, cv_key: str) -> Optional[Any]:
        """Resolve the solver for a specific CV (decision 8 dispatch).

        - ``self.solver = None`` → all CVs use the sequential body.
        - ``self.solver`` is a single :class:`StepSolver` → all
          CVs use that solver.
        - ``self.solver`` is a ``Dict[cv_key, StepSolver]`` → per-CV
          dispatch. CVs absent from the dict fall back to the
          sequential body (``None``).
        - ``self.solver`` is a future ``MultiCVStepSolver`` (third
          arm of the union annotation; placeholder for a system-
          level integrator that sees all CVs at once) → raise
          ``NotImplementedError``. No implementation in this phase
          per decision 8; the placeholder prevents a future
          backwards-compatible break for ``isinstance(solver, (...))``
          dispatch in user code.
        """
        solver = self.solver
        if solver is None:
            return None
        if isinstance(solver, dict):
            return solver.get(cv_key)  # absent → None → sequential body
        from .system_solver import SystemSolver
        if isinstance(solver, SystemSolver):
            # Axis 2 (SystemSolver / MultiCVStepSolver) objects must not be
            # passed as the Axis 1 `solver` param — use `system_solver=...`.
            raise NotImplementedError(
                "A MultiCVStepSolver / SystemSolver (Axis 2) was passed as "
                "the `solver` (Axis 1) parameter. Use `system_solver=...` "
                "for Axis 2 solvers. Pass a StepSolver or "
                "Dict[cv_key, StepSolver] for `solver`."
            )
        return solver

    def _initial_solve(self, cv: ControlVolume) -> None:
        """Run an initial speciation solve on a CV so pH[0] / I[0]
        are populated before the recorder's ``record_init`` reads them.

        No-op if the CV lacks a reaction_system or speciation engine.
        Mirrors the legacy run_batch initial-solve loop (lines 542-557).
        """
        if cv.reaction_system is None:
            return
        engine = getattr(cv.reaction_system, "engine", None)
        if engine is None:
            return
        liq = cv.phases.get("liquid")
        if liq is None:
            return
        try:
            result = engine.solve(phases=cv.phases, T_K=float(liq.T_K))
            result.apply_to_phases(cv.phases)
        except (TypeError, ValueError, AttributeError, KeyError) as exc:
            import warnings
            from ..monitoring.accuracy import AccuracyWarning
            warnings.warn(
                f"Initial speciation solve failed for CV {cv.label!r}: {exc}. "
                f"pH[0] and I[0] will be NaN; the run will continue.",
                AccuracyWarning,
                stacklevel=2,
            )

    # ── Repr ───────────────────────────────────────────────────────────

    def __repr__(self) -> str:
        n_links = len(self.links)
        n_ctrl = len(self.controllers)
        n_prof = len(self.profiles)
        keys = ", ".join(sorted(self.cvs.keys()))
        return (
            f"Simulation(label={self.label!r}, "
            f"cvs=[{keys}], links={n_links}, "
            f"controllers={n_ctrl}, profiles={n_prof})"
        )
