#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Demo: export BatchResult to DataFrame, CSV, and Parquet.

Shows the three export helpers added in the RESULT_EXPORT phase and
the common analytical patterns built on top of them.

Run from the repo root after ``pip install -e . pandas pyarrow``::

    python demos/model_api/export_results.py
"""

import sys
import tempfile
import pathlib
from pathlib import Path

import numpy as np

# Bootstrap: put models/ on sys.path so vlmodels is importable.
_repo_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_repo_root / "demos"))
import _bootstrap  # noqa: F401, E402

# ── 1. Build a small single-CV simulation ────────────────────────────
#
# Sparged batch fermenter, aerobic growth on acetic acid.  Same
# topology as demos/builder/batch_fermenter.py but constructed here
# from the factory helpers for brevity.

from vlmodels.fermenter.config import (
    FermenterFactory,
    VesselConfig,
    TransferConfig,
    GasFeedConfig,
)
from PyOMES.core import Simulation, SimultaneousEulerSolver

cv = FermenterFactory.create_volume(
    vessel=VesselConfig(V_total_L=10.0, T_K=305.15),
    transfer=TransferConfig.default_kinetic(kLa_O2=120.0),
    gas_feed=GasFeedConfig(vvm_min=0.5, composition={"O2": 0.21, "N2": 0.79}),
)

sim = Simulation(cvs={"main": cv}, solver=SimultaneousEulerSolver())

# ── 2. Run for 2 hours, 200 steps ────────────────────────────────────

result = sim.run(tau_h=2.0, n_steps=200)

print(f"Run complete: {len(result.t_h)} time points, "
      f"runtime {result.runtime_s*1000:.1f} ms")
print(f"Gas species:    {sorted(result.gas_mol['main'])}")
print(f"Liquid species: {sorted(result.liquid_mol['main'])}")
print()

# ── 3. Long-form DataFrame ────────────────────────────────────────────

df = result.to_dataframe()

print("Long-form DataFrame — first 5 rows:")
print(df.head())
print(f"\nShape: {df.shape}  (rows × columns)")
print(f"channel_kind values: {sorted(df['channel_kind'].unique())}")
print()

# ── 4. Groupby aggregation across species ────────────────────────────
#
# Mean value of each channel_kind across the whole run.

summary = (
    df.groupby(["cv_key", "channel_kind", "species"], dropna=False)["value"]
    .mean()
    .reset_index()
    .rename(columns={"value": "mean_value"})
)
print("Per-channel mean values (all CVs):")
print(summary.to_string(index=False))
print()

# ── 5. Wide-form for plotting ─────────────────────────────────────────
#
# to_wide() pivots the liquid phase_mol rows into a t_h × species
# DataFrame — the shape plotting libraries expect.

liquid_wide = result.to_wide("main", "liquid")
print("Wide-form liquid mol (first 3 rows):")
print(liquid_wide.head(3))
print()

gas_wide = result.to_wide("main", "gas")
print("Wide-form gas mol (first 3 rows):")
print(gas_wide.head(3))
print()

# Demonstrate: reconstruct a wide table manually from the long form
# (the pattern users would apply for multi-CV or cross-phase queries)
import pandas as pd

manual_wide = (
    df.query("cv_key == 'main' and channel_kind == 'phase_mol' and phase_key == 'liquid'")
    .pivot(index="t_h", columns="species", values="value")
    .rename_axis(None, axis="columns")
)
pd.testing.assert_frame_equal(liquid_wide, manual_wide)
print("Manual pivot matches to_wide() — OK")
print()

# ── 6. CSV round-trip ─────────────────────────────────────────────────

with tempfile.TemporaryDirectory() as tmpdir:
    csv_path = pathlib.Path(tmpdir) / "run.csv"
    result.to_csv(str(csv_path))

    df_back = pd.read_csv(str(csv_path))
    print(f"CSV written: {csv_path.name}  ({csv_path.stat().st_size // 1024} KB)")
    print(f"CSV re-read: {df_back.shape[0]} rows")

    # Verify numerical round-trip within %.10e precision
    orig_vals = df["value"].dropna().to_numpy()
    back_vals = df_back["value"].dropna().to_numpy()
    # %.10e preserves ~10 significant decimal digits; relative error ≲ 5e-11.
    # Use absolute tolerance of 1e-8 to cover values up to ~200 (P_atm).
    max_err = np.max(np.abs(orig_vals - back_vals))
    print(f"CSV max |value| round-trip error: {max_err:.2e}")
    assert max_err < 1e-8, f"CSV precision loss: {max_err}"
    print()

# ── 7. Parquet round-trip ─────────────────────────────────────────────

with tempfile.TemporaryDirectory() as tmpdir:
    pq_path = pathlib.Path(tmpdir) / "run.parquet"
    result.to_parquet(str(pq_path))

    df_pq = pd.read_parquet(str(pq_path))
    print(f"Parquet written: {pq_path.name}  ({pq_path.stat().st_size // 1024} KB)")
    print(f"Parquet re-read: {df_pq.shape[0]} rows")

    # Parquet is lossless — exact float64 equality
    orig_vals = df["value"].dropna().to_numpy()
    pq_vals = df_pq["value"].dropna().to_numpy()
    np.testing.assert_array_equal(orig_vals, pq_vals)
    print("Parquet round-trip: exact float64 equality — OK")
    print()

# ── 8. Selective Parquet column read ─────────────────────────────────
#
# Parquet's columnar layout lets pandas load only the columns you
# need, which matters for wide tables or very long runs.

with tempfile.TemporaryDirectory() as tmpdir:
    pq_path = pathlib.Path(tmpdir) / "run.parquet"
    result.to_parquet(str(pq_path))
    df_slim = pd.read_parquet(str(pq_path), columns=["t_h", "species", "value"])
    print(f"Selective Parquet read (3 columns): shape {df_slim.shape}")

print("\nAll export checks passed.")
