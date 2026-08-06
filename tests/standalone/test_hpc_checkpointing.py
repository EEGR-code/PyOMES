# -*- coding: utf-8 -*-
"""HPC_CHECKPOINTING C3/C4 tests.

Covers save_checkpoint / load_checkpoint (modes 'data' and 'data+env'),
checkpoint_t_h property, verify policy, and the full resume integration
test (save → load → run(start_t_h) → BatchResult.concat).
"""

import json

import numpy as np
import pytest


# ════════════════════════════════════════════════════════════════════════
#  Helpers
# ════════════════════════════════════════════════════════════════════════

def _make_sim(n_mol=None, label="main"):
    """Simple liquid-only Simulation ready to run."""
    from PyOMES.core import ControlVolume, LiquidPhase, Simulation
    n_mol = n_mol or {"S": 1.0}
    liq = LiquidPhase(n_mol=dict(n_mol), V_L=1.0, T_K=298.15)
    cv = ControlVolume(phases={"liquid": liq}, label=label)
    return Simulation(cvs={label: cv})


# ════════════════════════════════════════════════════════════════════════
#  C3 — checkpoint_t_h property
# ════════════════════════════════════════════════════════════════════════

class TestCheckpointTHProperty:

    def test_fresh_simulation_returns_none(self):
        """A freshly constructed Simulation has no checkpoint; property is None."""
        sim = _make_sim()
        assert sim.checkpoint_t_h is None

    def test_property_is_read_only(self):
        """checkpoint_t_h cannot be set directly — it is a read-only property."""
        sim = _make_sim()
        with pytest.raises(AttributeError):
            sim.checkpoint_t_h = 5.0  # type: ignore[misc]

    def test_set_internally_after_load(self, tmp_path):
        """After load_checkpoint, checkpoint_t_h returns the manifest's t_h."""
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=3.0, n_steps=6)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data")
        loaded = Simulation.load_checkpoint(tmp_path / "ckpt", verify="skip")
        assert loaded.checkpoint_t_h == pytest.approx(3.0)


# ════════════════════════════════════════════════════════════════════════
#  C3 — save_checkpoint, mode='data'
# ════════════════════════════════════════════════════════════════════════

class TestSaveCheckpointData:

    def test_creates_directory_and_required_files(self, tmp_path):
        """save_checkpoint creates the target directory and both required files."""
        sim = _make_sim()
        sim.run(tau_h=1.0, n_steps=4)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data")
        assert (tmp_path / "ckpt" / "simulation.pkl").exists()
        assert (tmp_path / "ckpt" / "manifest.json").exists()

    def test_creates_parent_directories(self, tmp_path):
        """save_checkpoint creates nested parent directories as needed."""
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        sim.save_checkpoint(tmp_path / "a" / "b" / "ckpt", mode="data")
        assert (tmp_path / "a" / "b" / "ckpt" / "simulation.pkl").exists()

    def test_manifest_has_all_required_fields(self, tmp_path):
        """manifest.json contains schema_version, mode, vlsim_version,
        timestamp, t_h, compression, and blas_info."""
        sim = _make_sim()
        sim.run(tau_h=1.0, n_steps=4)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data")
        manifest = json.loads((tmp_path / "ckpt" / "manifest.json").read_text())
        for field in ("schema_version", "mode", "vlsim_version",
                      "timestamp", "t_h", "compression", "blas_info"):
            assert field in manifest, f"manifest missing field {field!r}"

    def test_manifest_schema_version_is_1(self, tmp_path):
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data")
        manifest = json.loads((tmp_path / "ckpt" / "manifest.json").read_text())
        assert manifest["schema_version"] == 1

    def test_manifest_mode_is_data(self, tmp_path):
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data")
        manifest = json.loads((tmp_path / "ckpt" / "manifest.json").read_text())
        assert manifest["mode"] == "data"

    def test_manifest_t_h_matches_sim_accumulator(self, tmp_path):
        """t_h in the manifest equals sim._t_h at the moment of save."""
        sim = _make_sim()
        sim.run(tau_h=2.5, n_steps=10)
        t_h_at_save = sim._t_h
        sim.save_checkpoint(tmp_path / "ckpt", mode="data")
        manifest = json.loads((tmp_path / "ckpt" / "manifest.json").read_text())
        assert manifest["t_h"] == pytest.approx(t_h_at_save)

    def test_manifest_t_h_zero_before_run(self, tmp_path):
        """Saving before any run records t_h=0.0 in the manifest."""
        sim = _make_sim()
        sim.save_checkpoint(tmp_path / "ckpt", mode="data")
        manifest = json.loads((tmp_path / "ckpt" / "manifest.json").read_text())
        assert manifest["t_h"] == pytest.approx(0.0)

    def test_manifest_compression_field_matches_argument(self, tmp_path):
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        for comp in ("none", "gzip"):
            sim.save_checkpoint(tmp_path / f"ckpt_{comp}", mode="data", compression=comp)
            m = json.loads((tmp_path / f"ckpt_{comp}" / "manifest.json").read_text())
            assert m["compression"] == comp

    def test_gzip_produces_valid_compressed_bytes(self, tmp_path):
        """simulation.pkl written with gzip is decompressible."""
        import gzip
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data", compression="gzip")
        raw = (tmp_path / "ckpt" / "simulation.pkl").read_bytes()
        decompressed = gzip.decompress(raw)
        assert len(decompressed) > 0

    def test_none_compression_produces_raw_pickle(self, tmp_path):
        """simulation.pkl written with compression='none' is a valid pickle."""
        import pickle
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data", compression="none")
        raw = (tmp_path / "ckpt" / "simulation.pkl").read_bytes()
        loaded = pickle.loads(raw)
        from PyOMES.core import Simulation
        assert isinstance(loaded, Simulation)

    def test_not_implemented_source(self, tmp_path):
        sim = _make_sim()
        with pytest.raises(NotImplementedError):
            sim.save_checkpoint(tmp_path / "ckpt", mode="source")

    def test_not_implemented_full(self, tmp_path):
        sim = _make_sim()
        with pytest.raises(NotImplementedError):
            sim.save_checkpoint(tmp_path / "ckpt", mode="full")

    def test_invalid_mode_raises_value_error(self, tmp_path):
        sim = _make_sim()
        with pytest.raises(ValueError, match="mode"):
            sim.save_checkpoint(tmp_path / "ckpt", mode="invalid")

    def test_invalid_compression_raises_value_error(self, tmp_path):
        sim = _make_sim()
        with pytest.raises(ValueError, match="compression"):
            sim.save_checkpoint(tmp_path / "ckpt", mode="data", compression="bz2")


# ════════════════════════════════════════════════════════════════════════
#  C3 — load_checkpoint, mode='data'
# ════════════════════════════════════════════════════════════════════════

class TestLoadCheckpointData:

    def _save(self, tmp_path, sim, *, compression="gzip"):
        ckpt = tmp_path / "ckpt"
        sim.save_checkpoint(ckpt, mode="data", compression=compression)
        return ckpt

    def test_returns_simulation_instance(self, tmp_path):
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=1.0, n_steps=4)
        ckpt = self._save(tmp_path, sim)
        loaded = Simulation.load_checkpoint(ckpt, verify="skip")
        assert isinstance(loaded, Simulation)

    def test_checkpoint_t_h_set_from_manifest(self, tmp_path):
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=2.5, n_steps=10)
        ckpt = self._save(tmp_path, sim)
        loaded = Simulation.load_checkpoint(ckpt, verify="skip")
        assert loaded.checkpoint_t_h == pytest.approx(2.5)

    def test_cv_state_preserved_through_round_trip(self, tmp_path):
        """CV n_mol values survive save → pickle → load."""
        from PyOMES.core import Simulation
        sim = _make_sim({"S": 0.75, "X": 0.2})
        sim.run(tau_h=0.1, n_steps=2)
        ckpt = self._save(tmp_path, sim)
        loaded = Simulation.load_checkpoint(ckpt, verify="skip")
        orig_liq = sim["main"].phases["liquid"].n_mol
        load_liq = loaded["main"].phases["liquid"].n_mol
        assert load_liq["S"] == pytest.approx(orig_liq["S"])
        assert load_liq["X"] == pytest.approx(orig_liq["X"])

    def test_round_trip_compression_none(self, tmp_path):
        """compression='none' round-trip preserves state and checkpoint_t_h."""
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.5, n_steps=5)
        ckpt = self._save(tmp_path, sim, compression="none")
        loaded = Simulation.load_checkpoint(ckpt, verify="skip")
        assert isinstance(loaded, Simulation)
        assert loaded.checkpoint_t_h == pytest.approx(sim._t_h)

    def test_verify_skip_ignores_version_mismatch(self, tmp_path):
        """verify='skip' loads successfully even when vlsim_version is wrong."""
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        ckpt = self._save(tmp_path, sim)
        manifest = json.loads((ckpt / "manifest.json").read_text())
        manifest["vlsim_version"] = "0.0.0-bogus"
        (ckpt / "manifest.json").write_text(json.dumps(manifest))
        loaded = Simulation.load_checkpoint(ckpt, verify="skip")
        assert loaded.checkpoint_t_h == pytest.approx(sim._t_h)

    def test_verify_strict_raises_on_version_mismatch(self, tmp_path):
        """verify='strict' raises RuntimeError when vlsim_version differs."""
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        ckpt = self._save(tmp_path, sim)
        manifest = json.loads((ckpt / "manifest.json").read_text())
        manifest["vlsim_version"] = "0.0.0-bogus"
        (ckpt / "manifest.json").write_text(json.dumps(manifest))
        with pytest.raises(RuntimeError, match="version mismatch"):
            Simulation.load_checkpoint(ckpt, verify="strict")

    def test_verify_warn_warns_on_version_mismatch(self, tmp_path):
        """verify='warn' issues a UserWarning and still loads."""
        import warnings
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        ckpt = self._save(tmp_path, sim)
        manifest = json.loads((ckpt / "manifest.json").read_text())
        manifest["vlsim_version"] = "0.0.0-bogus"
        (ckpt / "manifest.json").write_text(json.dumps(manifest))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            loaded = Simulation.load_checkpoint(ckpt, verify="warn")
        user_warns = [w for w in caught if issubclass(w.category, UserWarning)]
        assert len(user_warns) == 1
        assert "version mismatch" in str(user_warns[0].message).lower()
        assert loaded.checkpoint_t_h == pytest.approx(sim._t_h)

    def test_verify_strict_does_not_raise_when_versions_match(self, tmp_path):
        """verify='strict' is silent when vlsim_version in manifest equals current."""
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        ckpt = self._save(tmp_path, sim)
        # Default save records the actual version → no mismatch
        loaded = Simulation.load_checkpoint(ckpt, verify="strict")
        assert loaded.checkpoint_t_h == pytest.approx(sim._t_h)

    def test_missing_manifest_raises_file_not_found(self, tmp_path):
        """FileNotFoundError when manifest.json is absent."""
        from PyOMES.core import Simulation
        ckpt = tmp_path / "ckpt"
        ckpt.mkdir()
        (ckpt / "simulation.pkl").write_bytes(b"x")
        with pytest.raises(FileNotFoundError, match="manifest.json"):
            Simulation.load_checkpoint(ckpt, verify="skip")

    def test_missing_pkl_raises_file_not_found(self, tmp_path):
        """FileNotFoundError when simulation.pkl is absent."""
        from PyOMES.core import Simulation
        ckpt = tmp_path / "ckpt"
        ckpt.mkdir()
        manifest = {
            "schema_version": 1, "mode": "data", "t_h": 0.0,
            "compression": "none", "vlsim_version": "unknown",
            "timestamp": "2026-01-01T00:00:00Z", "blas_info": "x",
        }
        (ckpt / "manifest.json").write_text(json.dumps(manifest))
        with pytest.raises(FileNotFoundError, match="simulation.pkl"):
            Simulation.load_checkpoint(ckpt, verify="skip")

    def test_invalid_verify_raises_value_error(self, tmp_path):
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        ckpt = self._save(tmp_path, sim)
        with pytest.raises(ValueError, match="verify"):
            Simulation.load_checkpoint(ckpt, verify="maybe")

    def test_loaded_sim_can_run(self, tmp_path):
        """The loaded simulation can immediately be passed to run()."""
        from PyOMES.core import BatchResult, Simulation
        sim = _make_sim()
        sim.run(tau_h=0.5, n_steps=5)
        ckpt = self._save(tmp_path, sim)
        loaded = Simulation.load_checkpoint(ckpt, verify="skip")
        result = loaded.run(tau_h=0.5, n_steps=5)
        assert isinstance(result, BatchResult)


# ════════════════════════════════════════════════════════════════════════
#  C3 — Integration: full HPC resume round-trip
# ════════════════════════════════════════════════════════════════════════

class TestCheckpointResumeIntegration:

    def test_save_load_run_concat_produces_clean_trajectory(self, tmp_path):
        """Full round-trip: run → save → load → run(start_t_h) → concat.

        Verifies that the concatenated trajectory has:
        - No duplicate boundary time point
        - Correct t_h[0] and t_h[-1]
        - Species arrays of the expected length
        """
        from PyOMES.core import BatchResult, Simulation

        sim1 = _make_sim({"S": 1.0})
        r1 = sim1.run(tau_h=1.0, n_steps=4)   # t_h: [0, 0.25, 0.5, 0.75, 1.0]
        sim1.save_checkpoint(tmp_path / "ckpt", mode="data")

        sim2 = Simulation.load_checkpoint(tmp_path / "ckpt", verify="skip")
        r2 = sim2.run(tau_h=1.0, n_steps=4,   # t_h: [1.0, 1.25, 1.5, 1.75, 2.0]
                      start_t_h=sim2.checkpoint_t_h)

        full = BatchResult.concat(r1, r2)

        # 5 + 4 = 9 time points (boundary at t=1.0 deduplicated)
        assert full.t_h.shape == (9,)
        assert full.t_h[0] == pytest.approx(0.0)
        assert full.t_h[4] == pytest.approx(1.0)
        assert full.t_h[-1] == pytest.approx(2.0)
        assert np.sum(np.abs(full.t_h - 1.0) < 1e-10) == 1  # exactly one boundary point

        # Species arrays match t_h length
        for arr in full.liquid_mol["main"].values():
            assert arr.shape == (9,)

    def test_checkpoint_t_h_equals_first_result_end(self, tmp_path):
        """sim.checkpoint_t_h after load equals r1.t_h[-1]."""
        from PyOMES.core import Simulation

        sim1 = _make_sim()
        r1 = sim1.run(tau_h=2.0, n_steps=8)
        sim1.save_checkpoint(tmp_path / "ckpt", mode="data")

        sim2 = Simulation.load_checkpoint(tmp_path / "ckpt", verify="skip")
        assert sim2.checkpoint_t_h == pytest.approx(r1.t_h[-1])

    def test_multiple_resume_segments(self, tmp_path):
        """Three sequential save→load→run segments concatenate cleanly."""
        from PyOMES.core import BatchResult, Simulation

        sim = _make_sim({"S": 1.0})
        r1 = sim.run(tau_h=1.0, n_steps=4)
        sim.save_checkpoint(tmp_path / "ckpt1", mode="data")

        sim = Simulation.load_checkpoint(tmp_path / "ckpt1", verify="skip")
        r2 = sim.run(tau_h=1.0, n_steps=4, start_t_h=sim.checkpoint_t_h)
        sim.save_checkpoint(tmp_path / "ckpt2", mode="data")

        sim = Simulation.load_checkpoint(tmp_path / "ckpt2", verify="skip")
        r3 = sim.run(tau_h=1.0, n_steps=4, start_t_h=sim.checkpoint_t_h)

        full = BatchResult.concat(r1, r2, r3)
        # 5 + 4 + 4 = 13 points covering [0, 3.0]
        assert full.t_h.shape == (13,)
        assert full.t_h[0] == pytest.approx(0.0)
        assert full.t_h[-1] == pytest.approx(3.0)


# ════════════════════════════════════════════════════════════════════════
#  C4 — save_checkpoint, mode='data+env'
# ════════════════════════════════════════════════════════════════════════

class TestSaveCheckpointDataPlusEnv:

    def test_creates_environment_subdirectory(self, tmp_path):
        """mode='data+env' creates an environment/ subdirectory."""
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data+env")
        assert (tmp_path / "ckpt" / "environment").is_dir()

    def test_all_four_env_files_exist(self, tmp_path):
        """All required environment files are written."""
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data+env")
        env = tmp_path / "ckpt" / "environment"
        for fname in ("python_version.txt", "platform.txt",
                      "omp_num_threads.txt", "pip_freeze.txt"):
            assert (env / fname).exists(), f"missing {fname}"

    def test_python_version_matches_current(self, tmp_path):
        """environment/python_version.txt contains sys.version."""
        import sys
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data+env")
        saved = (tmp_path / "ckpt" / "environment" / "python_version.txt").read_text()
        assert saved.strip() == sys.version.strip()

    def test_platform_matches_current(self, tmp_path):
        """environment/platform.txt contains platform.platform()."""
        import platform
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data+env")
        saved = (tmp_path / "ckpt" / "environment" / "platform.txt").read_text()
        assert saved.strip() == platform.platform().strip()

    def test_omp_num_threads_recorded(self, tmp_path):
        """environment/omp_num_threads.txt contains OMP_NUM_THREADS or '1'."""
        import os
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data+env")
        saved = (tmp_path / "ckpt" / "environment" / "omp_num_threads.txt").read_text()
        assert saved.strip() == os.environ.get("OMP_NUM_THREADS", "1")

    def test_manifest_mode_is_data_plus_env(self, tmp_path):
        """manifest.json records mode='data+env'."""
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data+env")
        manifest = json.loads((tmp_path / "ckpt" / "manifest.json").read_text())
        assert manifest["mode"] == "data+env"

    def test_simulation_pkl_and_manifest_also_present(self, tmp_path):
        """mode='data+env' still writes simulation.pkl and manifest.json."""
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        sim.save_checkpoint(tmp_path / "ckpt", mode="data+env")
        assert (tmp_path / "ckpt" / "simulation.pkl").exists()
        assert (tmp_path / "ckpt" / "manifest.json").exists()


# ════════════════════════════════════════════════════════════════════════
#  C4 — load_checkpoint, mode='data+env', drift checks
# ════════════════════════════════════════════════════════════════════════

class TestLoadCheckpointDataPlusEnv:

    def _save_env(self, tmp_path, sim):
        ckpt = tmp_path / "ckpt"
        sim.save_checkpoint(ckpt, mode="data+env")
        return ckpt

    def _corrupt_env_file(self, ckpt, filename, content):
        (ckpt / "environment" / filename).write_text(content)

    def test_clean_round_trip_no_warnings(self, tmp_path):
        """save → load with identical env produces no warnings or raises."""
        import warnings
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.5, n_steps=5)
        ckpt = self._save_env(tmp_path, sim)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            loaded = Simulation.load_checkpoint(ckpt, verify="strict")
        user_warns = [w for w in caught if issubclass(w.category, UserWarning)]
        assert user_warns == []
        assert loaded.checkpoint_t_h == pytest.approx(sim._t_h)

    def test_verify_strict_raises_on_python_version_mismatch(self, tmp_path):
        """Python version mismatch raises RuntimeError with verify='strict'."""
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        ckpt = self._save_env(tmp_path, sim)
        self._corrupt_env_file(ckpt, "python_version.txt",
                               "3.0.0 (fake, fake, fake)")
        with pytest.raises(RuntimeError, match="Python version mismatch"):
            Simulation.load_checkpoint(ckpt, verify="strict")

    def test_verify_warn_warns_on_python_version_mismatch(self, tmp_path):
        """Python version mismatch issues UserWarning with verify='warn'."""
        import warnings
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        ckpt = self._save_env(tmp_path, sim)
        self._corrupt_env_file(ckpt, "python_version.txt",
                               "3.0.0 (fake, fake, fake)")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            loaded = Simulation.load_checkpoint(ckpt, verify="warn")
        user_warns = [w for w in caught if issubclass(w.category, UserWarning)]
        assert any("Python version mismatch" in str(w.message) for w in user_warns)
        assert loaded.checkpoint_t_h == pytest.approx(sim._t_h)

    def test_verify_skip_ignores_python_version_mismatch(self, tmp_path):
        """verify='skip' loads without checking any env files."""
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        ckpt = self._save_env(tmp_path, sim)
        self._corrupt_env_file(ckpt, "python_version.txt", "1.0.0")
        loaded = Simulation.load_checkpoint(ckpt, verify="skip")
        assert loaded.checkpoint_t_h == pytest.approx(sim._t_h)

    def test_verify_strict_warns_not_raises_on_platform_mismatch(self, tmp_path):
        """Platform mismatch is a warning even with verify='strict'."""
        import warnings
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        ckpt = self._save_env(tmp_path, sim)
        self._corrupt_env_file(ckpt, "platform.txt",
                               "Linux-99.0-fake-x86_64-with-fake")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            loaded = Simulation.load_checkpoint(ckpt, verify="strict")
        user_warns = [w for w in caught if issubclass(w.category, UserWarning)]
        assert any("Platform mismatch" in str(w.message) for w in user_warns)
        assert loaded.checkpoint_t_h == pytest.approx(sim._t_h)

    def test_verify_strict_warns_on_omp_mismatch(self, tmp_path):
        """OMP_NUM_THREADS mismatch is a warning even with verify='strict'."""
        import warnings
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        ckpt = self._save_env(tmp_path, sim)
        self._corrupt_env_file(ckpt, "omp_num_threads.txt", "64")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            loaded = Simulation.load_checkpoint(ckpt, verify="strict")
        user_warns = [w for w in caught if issubclass(w.category, UserWarning)]
        assert any("OMP_NUM_THREADS" in str(w.message) for w in user_warns)
        assert loaded.checkpoint_t_h == pytest.approx(sim._t_h)

    def test_verify_strict_warns_on_pip_freeze_mismatch(self, tmp_path):
        """Different pip freeze output is a warning even with verify='strict'."""
        import warnings
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        ckpt = self._save_env(tmp_path, sim)
        # Add a fake package that wouldn't be in the current env
        orig = (ckpt / "environment" / "pip_freeze.txt").read_text()
        (ckpt / "environment" / "pip_freeze.txt").write_text(
            orig + "\nfake-nonexistent-pkg==999.0.0\n"
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            loaded = Simulation.load_checkpoint(ckpt, verify="strict")
        user_warns = [w for w in caught if issubclass(w.category, UserWarning)]
        assert any("package" in str(w.message).lower() for w in user_warns)
        assert loaded.checkpoint_t_h == pytest.approx(sim._t_h)

    def test_verify_strict_warns_on_blas_mismatch(self, tmp_path):
        """BLAS info mismatch is a warning even with verify='strict'."""
        import warnings
        from PyOMES.core import Simulation
        sim = _make_sim()
        sim.run(tau_h=0.1, n_steps=2)
        ckpt = self._save_env(tmp_path, sim)
        manifest = json.loads((ckpt / "manifest.json").read_text())
        manifest["blas_info"] = "FakeBLAS-99"
        (ckpt / "manifest.json").write_text(json.dumps(manifest))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            loaded = Simulation.load_checkpoint(ckpt, verify="strict")
        user_warns = [w for w in caught if issubclass(w.category, UserWarning)]
        assert any("BLAS" in str(w.message) for w in user_warns)
        assert loaded.checkpoint_t_h == pytest.approx(sim._t_h)

    def test_integration_data_plus_env_resume(self, tmp_path):
        """Full data+env round-trip: run → save → load → run → concat."""
        from PyOMES.core import BatchResult, Simulation
        sim1 = _make_sim({"S": 1.0})
        r1 = sim1.run(tau_h=1.0, n_steps=4)
        sim1.save_checkpoint(tmp_path / "ckpt", mode="data+env")

        sim2 = Simulation.load_checkpoint(tmp_path / "ckpt", verify="strict")
        r2 = sim2.run(tau_h=1.0, n_steps=4, start_t_h=sim2.checkpoint_t_h)

        full = BatchResult.concat(r1, r2)
        assert full.t_h.shape == (9,)
        assert full.t_h[0] == pytest.approx(0.0)
        assert full.t_h[-1] == pytest.approx(2.0)
