# -*- coding: utf-8 -*-
"""Tests for BatchResult export helpers (RESULT_EXPORT phase).

Covers:
- to_dataframe(): column inventory, row counts, value correctness,
  NaN handling for absent pH / ionic_strength
- to_wide(): pivot shape and index/columns
- to_csv(): round-trip value equality, float_format override
- to_parquet(): round-trip value equality
- Missing-dep behaviour: helpful ImportError messages
"""

import math
import sys
import types

import numpy as np
import pytest


# ═══════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════

def _make_batch_result(
    cv_keys=("main",),
    gas_species=("O2", "N2"),
    liq_species=("S", "X"),
    n_steps=4,
    seed=0,
):
    """Build a minimal BatchResult with deterministic NumPy arrays.

    Produces *n_steps + 1* time points.  pH and ionic_strength are
    populated; P_atm is set to a constant per CV.  All values are
    strictly positive to simplify round-trip assertions.
    """
    from PyOMES.core import BatchResult

    rng = np.random.default_rng(seed)
    n = n_steps + 1
    t_h = np.linspace(0.0, 1.0, n)

    phase_mol = {}
    pH = {}
    ionic_strength = {}
    P_atm = {}

    for cv in cv_keys:
        phase_mol[cv] = {
            "gas": {sp: rng.random(n) + 0.1 for sp in gas_species},
            "liquid": {sp: rng.random(n) + 0.1 for sp in liq_species},
        }
        pH[cv] = np.full(n, 7.0 + rng.random())
        ionic_strength[cv] = np.full(n, 0.05 + rng.random() * 0.01)
        P_atm[cv] = np.full(n, 1.0 + rng.random() * 0.1)

    return BatchResult(
        t_h=t_h,
        phase_mol=phase_mol,
        pH=pH,
        ionic_strength=ionic_strength,
        P_atm=P_atm,
    )


def _make_batch_result_no_scalars(n_steps=3):
    """BatchResult with pH / ionic_strength left as NaN arrays
    (as produced when no speciation engine runs)."""
    from PyOMES.core import BatchResult

    n = n_steps + 1
    t_h = np.linspace(0.0, 0.5, n)
    return BatchResult(
        t_h=t_h,
        phase_mol={"main": {
            "gas": {"O2": np.ones(n)},
            "liquid": {"S": np.ones(n) * 2.0},
        }},
        pH={"main": np.full(n, np.nan)},
        ionic_strength={"main": np.full(n, np.nan)},
        P_atm={"main": np.ones(n) * 1.05},
    )


# ═══════════════════════════════════════════════════════════════════════
#  to_dataframe — column inventory and shape
# ═══════════════════════════════════════════════════════════════════════

class TestToDataframeShape:

    def test_column_names(self):
        pytest.importorskip("pandas")
        result = _make_batch_result()
        df = result.to_dataframe()
        assert list(df.columns) == [
            "t_h", "cv_key", "channel_kind", "phase_key", "species", "value"
        ]

    def test_row_count_single_cv(self):
        """Total rows = sum of rows from all channels.

        Single CV, 2 gas species, 2 liquid species, n_steps=4 → n=5 time points.
          phase_mol gas:    2 species × 5 = 10
          phase_mol liquid: 2 species × 5 = 10
          pH:               1 × 5 = 5
          ionic_strength:   1 × 5 = 5
          P_atm:            1 × 5 = 5
          total = 35
        """
        pytest.importorskip("pandas")
        n_steps = 4
        result = _make_batch_result(n_steps=n_steps)
        df = result.to_dataframe()
        n = n_steps + 1
        expected = (2 + 2 + 1 + 1 + 1) * n  # 6 channels × 5 time points
        assert len(df) == expected

    def test_row_count_multi_cv(self):
        pytest.importorskip("pandas")
        n_steps = 2
        result = _make_batch_result(cv_keys=("A", "B"), n_steps=n_steps)
        df = result.to_dataframe()
        n = n_steps + 1
        # Per CV: (2 gas + 2 liq + pH + IS + P_atm) × n = 7 channels
        expected = 2 * 7 * n
        assert len(df) == expected

    def test_channel_kind_values(self):
        pytest.importorskip("pandas")
        result = _make_batch_result()
        df = result.to_dataframe()
        kinds = set(df["channel_kind"].unique())
        assert kinds == {"phase_mol", "pH", "ionic_strength", "P_atm"}

    def test_phase_key_for_mol_rows(self):
        pytest.importorskip("pandas")
        result = _make_batch_result()
        df = result.to_dataframe()
        mol_rows = df[df["channel_kind"] == "phase_mol"]
        phase_keys = set(mol_rows["phase_key"].unique())
        assert phase_keys == {"gas", "liquid"}

    def test_phase_key_none_for_scalar_channels(self):
        pytest.importorskip("pandas")
        result = _make_batch_result()
        df = result.to_dataframe()
        scalar_rows = df[df["channel_kind"].isin(["pH", "ionic_strength", "P_atm"])]
        assert scalar_rows["phase_key"].isna().all()

    def test_species_none_for_scalar_channels(self):
        pytest.importorskip("pandas")
        result = _make_batch_result()
        df = result.to_dataframe()
        scalar_rows = df[df["channel_kind"].isin(["pH", "ionic_strength", "P_atm"])]
        assert scalar_rows["species"].isna().all()

    def test_species_populated_for_mol_rows(self):
        pytest.importorskip("pandas")
        result = _make_batch_result(gas_species=("O2",), liq_species=("S",))
        df = result.to_dataframe()
        mol_rows = df[df["channel_kind"] == "phase_mol"]
        assert set(mol_rows["species"].unique()) == {"O2", "S"}

    def test_t_h_range(self):
        pytest.importorskip("pandas")
        result = _make_batch_result(n_steps=4)
        df = result.to_dataframe()
        assert df["t_h"].min() == pytest.approx(0.0)
        assert df["t_h"].max() == pytest.approx(1.0)


# ═══════════════════════════════════════════════════════════════════════
#  to_dataframe — value correctness
# ═══════════════════════════════════════════════════════════════════════

class TestToDataframeValues:

    def test_gas_mol_values_match_arrays(self):
        pytest.importorskip("pandas")
        result = _make_batch_result(gas_species=("O2",), liq_species=(), n_steps=3)
        df = result.to_dataframe()
        rows = df[
            (df["channel_kind"] == "phase_mol")
            & (df["phase_key"] == "gas")
            & (df["species"] == "O2")
        ].sort_values("t_h")
        np.testing.assert_allclose(
            rows["value"].to_numpy(),
            result.gas_mol["main"]["O2"],
        )

    def test_liquid_mol_values_match_arrays(self):
        pytest.importorskip("pandas")
        result = _make_batch_result(gas_species=(), liq_species=("S",), n_steps=3)
        df = result.to_dataframe()
        rows = df[
            (df["channel_kind"] == "phase_mol")
            & (df["phase_key"] == "liquid")
            & (df["species"] == "S")
        ].sort_values("t_h")
        np.testing.assert_allclose(
            rows["value"].to_numpy(),
            result.liquid_mol["main"]["S"],
        )

    def test_pH_values_match_array(self):
        pytest.importorskip("pandas")
        result = _make_batch_result(n_steps=3)
        df = result.to_dataframe()
        rows = df[df["channel_kind"] == "pH"].sort_values("t_h")
        np.testing.assert_allclose(
            rows["value"].to_numpy(),
            result.pH["main"],
        )

    def test_nan_pH_preserved(self):
        pytest.importorskip("pandas")
        result = _make_batch_result_no_scalars()
        df = result.to_dataframe()
        ph_rows = df[df["channel_kind"] == "pH"]
        assert ph_rows["value"].isna().all()

    def test_P_atm_values_match_array(self):
        pytest.importorskip("pandas")
        result = _make_batch_result(n_steps=2)
        df = result.to_dataframe()
        rows = df[df["channel_kind"] == "P_atm"].sort_values("t_h")
        np.testing.assert_allclose(
            rows["value"].to_numpy(),
            result.P_atm["main"],
        )

    def test_empty_gas_mol_produces_no_phase_mol_gas_rows(self):
        pytest.importorskip("pandas")
        result = _make_batch_result(gas_species=(), liq_species=("S",), n_steps=2)
        df = result.to_dataframe()
        gas_rows = df[
            (df["channel_kind"] == "phase_mol") & (df["phase_key"] == "gas")
        ]
        assert len(gas_rows) == 0


# ═══════════════════════════════════════════════════════════════════════
#  to_wide
# ═══════════════════════════════════════════════════════════════════════

class TestToWide:

    def test_index_is_t_h(self):
        pytest.importorskip("pandas")
        result = _make_batch_result(n_steps=4)
        wide = result.to_wide("main", "liquid")
        assert wide.index.name == "t_h"
        np.testing.assert_allclose(wide.index.to_numpy(), result.t_h)

    def test_columns_are_species(self):
        pytest.importorskip("pandas")
        result = _make_batch_result(liq_species=("S", "X"), n_steps=2)
        wide = result.to_wide("main", "liquid")
        assert set(wide.columns) == {"S", "X"}

    def test_values_match_liquid_mol(self):
        pytest.importorskip("pandas")
        result = _make_batch_result(liq_species=("S",), n_steps=3)
        wide = result.to_wide("main", "liquid")
        np.testing.assert_allclose(wide["S"].to_numpy(), result.liquid_mol["main"]["S"])

    def test_gas_phase_key(self):
        pytest.importorskip("pandas")
        result = _make_batch_result(gas_species=("O2", "N2"), n_steps=2)
        wide = result.to_wide("main", "gas")
        assert set(wide.columns) == {"O2", "N2"}
        np.testing.assert_allclose(
            wide["O2"].to_numpy(), result.gas_mol["main"]["O2"]
        )

    def test_shape(self):
        pytest.importorskip("pandas")
        n_steps = 5
        result = _make_batch_result(
            liq_species=("S", "X", "P"), n_steps=n_steps
        )
        wide = result.to_wide("main", "liquid")
        assert wide.shape == (n_steps + 1, 3)


# ═══════════════════════════════════════════════════════════════════════
#  to_csv round-trip
# ═══════════════════════════════════════════════════════════════════════

class TestToCsvRoundTrip:

    def test_round_trip_value_equality(self, tmp_path):
        pd = pytest.importorskip("pandas")
        result = _make_batch_result(n_steps=3)
        path = tmp_path / "out.csv"
        result.to_csv(str(path))
        df_back = pd.read_csv(str(path))
        df_orig = result.to_dataframe()
        # Same shape
        assert df_back.shape == df_orig.shape
        # Numeric values preserved within float64 precision
        np.testing.assert_allclose(
            df_back["value"].dropna().to_numpy(),
            df_orig["value"].dropna().to_numpy(),
            rtol=1e-9,
        )

    def test_default_float_format_preserves_precision(self, tmp_path):
        """%.10e gives 10 significant decimal digits — check that a
        known value survives the round-trip within that tolerance."""
        pd = pytest.importorskip("pandas")
        from PyOMES.core import BatchResult

        sentinel = 1.23456789012345678
        result = BatchResult(
            t_h=np.array([0.0]),
            phase_mol={"main": {"gas": {"O2": np.array([sentinel])}, "liquid": {}}},
            pH={"main": np.array([np.nan])},
            ionic_strength={"main": np.array([np.nan])},
            P_atm={"main": np.array([1.0])},
        )
        path = tmp_path / "precision.csv"
        result.to_csv(str(path))
        df_back = pd.read_csv(str(path))
        recovered = df_back.loc[
            (df_back["channel_kind"] == "phase_mol")
            & (df_back["species"] == "O2"),
            "value",
        ].iloc[0]
        assert abs(recovered - sentinel) < 1e-10

    def test_float_format_override(self, tmp_path):
        pd = pytest.importorskip("pandas")
        result = _make_batch_result(n_steps=2)
        path = tmp_path / "short.csv"
        result.to_csv(str(path), float_format="%.3f")
        # File exists and is non-empty
        assert path.stat().st_size > 0
        df = pd.read_csv(str(path))
        assert len(df) > 0

    def test_csv_has_header_row(self, tmp_path):
        pytest.importorskip("pandas")
        result = _make_batch_result(n_steps=1)
        path = tmp_path / "hdr.csv"
        result.to_csv(str(path))
        first_line = path.read_text().splitlines()[0]
        assert "t_h" in first_line
        assert "cv_key" in first_line


# ═══════════════════════════════════════════════════════════════════════
#  to_parquet round-trip
# ═══════════════════════════════════════════════════════════════════════

class TestToParquetRoundTrip:

    def test_round_trip_value_equality(self, tmp_path):
        pd = pytest.importorskip("pandas")
        pytest.importorskip("pyarrow")
        result = _make_batch_result(n_steps=3)
        path = tmp_path / "out.parquet"
        result.to_parquet(str(path))
        df_back = pd.read_parquet(str(path))
        df_orig = result.to_dataframe()
        assert df_back.shape == df_orig.shape
        # Parquet is lossless for float64 — exact equality
        np.testing.assert_array_equal(
            df_back["value"].dropna().to_numpy(),
            df_orig["value"].dropna().to_numpy(),
        )

    def test_parquet_file_smaller_than_csv(self, tmp_path):
        pytest.importorskip("pandas")
        pytest.importorskip("pyarrow")
        result = _make_batch_result(n_steps=20)
        csv_path = tmp_path / "out.csv"
        pq_path = tmp_path / "out.parquet"
        result.to_csv(str(csv_path))
        result.to_parquet(str(pq_path))
        assert pq_path.stat().st_size < csv_path.stat().st_size

    def test_column_names_preserved(self, tmp_path):
        pd = pytest.importorskip("pandas")
        pytest.importorskip("pyarrow")
        result = _make_batch_result(n_steps=2)
        path = tmp_path / "cols.parquet"
        result.to_parquet(str(path))
        df = pd.read_parquet(str(path))
        assert list(df.columns) == [
            "t_h", "cv_key", "channel_kind", "phase_key", "species", "value"
        ]


# ═══════════════════════════════════════════════════════════════════════
#  Missing-dependency behaviour
# ═══════════════════════════════════════════════════════════════════════

class TestMissingDeps:

    def _block_import(self, monkeypatch, module_name):
        """Make ``import <module_name>`` raise ImportError."""
        sentinel = types.ModuleType(module_name)
        monkeypatch.setitem(sys.modules, module_name, None)

    def test_to_dataframe_raises_without_pandas(self, monkeypatch):
        self._block_import(monkeypatch, "pandas")
        result = _make_batch_result()
        with pytest.raises(ImportError, match="pip install pandas"):
            result.to_dataframe()

    def test_to_csv_raises_without_pandas(self, monkeypatch, tmp_path):
        self._block_import(monkeypatch, "pandas")
        result = _make_batch_result()
        with pytest.raises(ImportError, match="pip install pandas"):
            result.to_csv(str(tmp_path / "out.csv"))

    def test_to_parquet_raises_without_pyarrow(self, monkeypatch, tmp_path):
        self._block_import(monkeypatch, "pyarrow")
        result = _make_batch_result()
        with pytest.raises(ImportError, match="pip install pandas pyarrow"):
            result.to_parquet(str(tmp_path / "out.parquet"))

    def test_error_message_mentions_install_command(self, monkeypatch):
        self._block_import(monkeypatch, "pandas")
        result = _make_batch_result()
        with pytest.raises(ImportError) as exc_info:
            result.to_dataframe()
        assert "pip install" in str(exc_info.value)

    def test_numpy_arrays_accessible_without_pandas(self, monkeypatch):
        """The underlying NumPy arrays must remain usable when pandas
        is absent — to_dataframe() failing must not corrupt state."""
        self._block_import(monkeypatch, "pandas")
        result = _make_batch_result(n_steps=2)
        # Trigger the ImportError
        with pytest.raises(ImportError):
            result.to_dataframe()
        # Arrays are still intact
        assert result.liquid_mol["main"]["S"].shape == (3,)
        assert result.t_h.shape == (3,)
