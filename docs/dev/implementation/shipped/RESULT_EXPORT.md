# Result Export — Planning Note

> **Status: SHIPPED 2026-06-01.** Tag: `result-export-shipped`
> (merge `b9f44c4`). `to_dataframe` / `to_csv` / `to_parquet` /
> `to_wide` added to `BatchResult`; 5 NumPy channels; pandas +
> pyarrow optional; 32 tests; demo at
> `demos/model_api/export_results.py`. 1009/0 tests post-merge.

> ~~**Status:** Trigger-gated, low-cost. Surfaced 2026-05-28
> during the post-simulation-class exploration session. The
> scope is small (~80 LOC + 1 demo) and can be picked up at
> any time; the only caveat is that
> [RUN_HISTORY.md](RUN_HISTORY.md) will eventually reshape
> `BatchResult`, at which point these helpers refactor with~~
> it.

## Context

`BatchResult` ([`src/core/recorder.py`](../../src/core/recorder.py))
is a plain dataclass holding NumPy arrays nested by CV key
and phase key. Users wanting their data in pandas / CSV /
Parquet / Excel write the conversion themselves. The
conversion is mechanical but writing it once and shipping it
as part of VLsim saves every user the same boilerplate, and
the conversion lives in one place where its shape can evolve
with `BatchResult`.

## Scope

### Source surface

Add to [`src/core/recorder.py`](../../src/core/recorder.py)
on `BatchResult`:

- `to_dataframe(self) -> pandas.DataFrame` — long-form
  conversion. Each row is a `(t_h, cv_key, channel_kind,
  phase_key, species, value)` tuple. `channel_kind` ∈
  `{"phase_mol", "phase_property", "pH", "ionic_strength",
  "P_atm", "boundary_flux", "transfer_flow",
  "controller_action", "profile"}` (matching the existing
  `BatchResult` fields). Long form rather than wide because
  the species inventory varies across CVs; a long-form table
  joins cleanly across CVs/phases.
- `to_csv(self, path, **kwargs) -> None` — thin wrapper
  around `to_dataframe().to_csv(path, **kwargs)`.
- `to_parquet(self, path, **kwargs) -> None` — thin wrapper
  around `to_dataframe().to_parquet(path, **kwargs)`.
- `pandas` and `pyarrow` as optional dependencies. If absent,
  `to_dataframe()` raises `ImportError` with a clear pointer
  (`pip install pandas pyarrow`); the underlying NumPy arrays
  remain accessible via the existing fields, so users can
  still write their own export.

### Demo

`demos/model_api/export_results.py` (~120 LOC), or a sibling
subtree (e.g. `demos/model_api/results/`) once result-handling
accumulates more than this single file. Original scope used
`demos/api/`; renamed to `demos/model_api/` by the
[DEMO_RESTRUCTURE.md](DEMO_RESTRUCTURE.md) ship on 2026-05-29:

1. Build a small CV simulation (single CV, batch fermenter
   shape, ~10 species).
2. Run it for a short horizon.
3. Show `result.to_dataframe()`, `.head()`,
   `.groupby("species")` aggregation.
4. Show `to_csv()` and `to_parquet()` round-trips.
5. Show how to recover a wide-form per-species time series
   from the long-form DataFrame for plotting.

### Tests

`tests/standalone/test_result_export.py`:

- Round-trip: write → read → check value-level equality.
- Long-form shape: assert the column inventory matches the
  spec above.
- Missing-dep behaviour: monkeypatch `import pandas` to raise
  and assert the error message points the user at the
  install command.

## The refactor-risk banner

This work ships **with awareness** that the long-form schema
is keyed on today's `BatchResult` shape:

- `BatchResult.gas_mol`, `liquid_mol`, `pH`, `ionic_strength`,
  `P_atm`, `boundary_records`, `transfer_records`,
  `controller_records`, `profile_records`.

[SNAPSHOT_PHASE_AGNOSTIC.md](SNAPSHOT_PHASE_AGNOSTIC.md)
proposes replacing `gas_mol` / `liquid_mol` with a general
`phase_mol[cv_key][phase_key][species]` channel (keeping
the named aliases as views). When that ships, the
long-form schema needs one extra column (`phase_key`) and
the alias-channel rows fold into the general channel.

[RUN_HISTORY.md](RUN_HISTORY.md)'s richer-recorder work
may extend the channel inventory further (StreamingFileRecorder
introduces a chunk manifest; SummaryRecorder produces a
final-state-only view that doesn't have a time index).

**The shape chosen here — long-form with a `channel_kind`
column — is the migration-friendly target.** When SNAPSHOT
ships, `phase_mol` rows just gain a `phase_key` value that
isn't `None`; existing `gas_mol` / `liquid_mol` rows can
either stay as legacy aliases or be remapped at read time.
The user-facing API surface (`to_dataframe`, `to_csv`,
`to_parquet`) does not change.

## Open design questions

1. **Long form vs wide form by default?** Long form
   (proposed) is migration-friendly and joins across CVs;
   wide form (one column per species per phase per CV) is
   what most plotting code expects. Probably ship long-form
   as canonical with a `result.to_wide(cv_key, phase_key) →
   DataFrame` helper for the common single-CV plot.
2. **Should `to_dataframe()` materialise the full table at
   once, or yield row chunks for very long runs?** Long
   runs aren't the natural shape for `BatchResult` (which
   pre-allocates everything in memory); for those, the
   `StreamingFileRecorder` from RUN_HISTORY is the right
   answer, and its Parquet chunks are already pandas-ready.
   So `to_dataframe()` materialises by default; chunked
   reading is the streaming recorder's job.
3. **CSV float format.** Pandas default is `%g` which loses
   precision. Should `to_csv()` set `float_format="%.10e"`
   by default? Lean yes — users who want short numbers can
   override.
4. **Excel?** `to_excel()` would need `openpyxl` and gives
   one sheet per channel-kind for usability. Skip for v1
   (Excel is the worst format for time-series data); add
   only if a real user asks.

## Trigger conditions

Any one of:

1. **User-facing pull.** Anyone using VLsim asks "how do I
   get this into pandas / CSV / a notebook?". Ship today's
   scope as a one-day shipment.
2. **Bundled with `demos/model_api/` writes.** The
   `export_results.py` demo lives in `demos/model_api/`; if that
   subtree is being populated anyway (see
   [DEMO_RESTRUCTURE.md](DEMO_RESTRUCTURE.md)), shipping the
   helper alongside is incremental.
3. **Folded into RUN_HISTORY.** When RUN_HISTORY ships,
   these helpers are part of the recorder-suite redesign.
   Shipping them now is decoupled but acceptable; shipping
   them with RUN_HISTORY saves one refactor.

## Relationship to other phases

- **Independent of** every shipped or trigger-gated phase
  except as noted below.
- **Refactor-coupled with** [RUN_HISTORY.md](RUN_HISTORY.md)
  and [SNAPSHOT_PHASE_AGNOSTIC.md](SNAPSHOT_PHASE_AGNOSTIC.md).
  Both will reshape `BatchResult`; the long-form export
  schema is designed to migrate cleanly.
- **Demo lives in** the `demos/model_api/` subtree shipped via
  [DEMO_RESTRUCTURE.md](DEMO_RESTRUCTURE.md) on 2026-05-29
  (originally scoped as `demos/api/`, renamed during the
  restructure).
- **Follows** `simulation-class` (shipped 2026-05-27), which
  installed `BatchResult` in its current shape.
