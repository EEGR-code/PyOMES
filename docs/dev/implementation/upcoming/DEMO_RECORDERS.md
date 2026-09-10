# Recorder Demo — Planning Note

> **Status:** Planned. Not yet implemented.

A short demo script illustrating all four recorder variants side-by-side
on the same simulation.  Intended audience: users who want to choose the
right recorder for a given use-case.

## What this is

Three new recorder implementations shipped in RUN_HISTORY (2026-06-08):
`StreamingFileRecorder`, `SparseRecorder`, and `SummaryRecorder`.  None of
the existing demos in `demos/` exercise them — the builder demos all use the
default `BatchRecorder` implicitly.  A dedicated demo closes that gap.

## Proposed location

```
demos/model_api/recorder_comparison.py
```

Sits alongside `export_results.py` in `demos/model_api/`.

## Demo outline

Run the same simple fed-batch simulation four times (or once with four
independent `Simulation` snapshots), each with a different recorder:

1. **`BatchRecorder` (default)** — baseline; show `result.t_h` shape and
   `result.phase_mol`.
2. **`StreamingFileRecorder`** — point at a `tmp/` directory; show chunk
   files written to disk; call `load_run(path)` to reconstruct.  Print
   chunk count and manifest summary.
3. **`SparseRecorder(every_m_steps=50)`** — show that `result.t_h` is
   much shorter than the batch result; plot or print selected time points.
4. **`SummaryRecorder`** — show `SummaryResult.initial_phase_mol`,
   `final_phase_mol`, and `boundary_totals` (if a feed boundary is present).

## Key points the demo should make legible

- `StreamingFileRecorder` is crash-resilient: the Parquet chunks on disk
  are readable even if `finalize()` is never called.  Show this by manually
  calling `_flush()` and then reading the chunk with `pandas.read_parquet`.
- `SparseRecorder` trades resolution for memory: a 10 000-step run with
  `every_m_steps=100` keeps 102 rows instead of 10 001.
- `SummaryRecorder` returns `SummaryResult`, not `BatchResult` — the return
  type is different by design.  Useful for parameter sweeps.
- All four satisfy `Simulation(recorder=<recorder>)` — the API surface is
  identical at the call site.

## Dependency

Requires `pyarrow` for `StreamingFileRecorder` / `load_run`.  Demo should
call `pytest.importorskip`-equivalent at the top (or guard with a try/except)
so it fails gracefully if `pyarrow` is absent.

## Trigger

Write this when the recorder infrastructure has been exercised in a real
model (or when a user asks for an example).  Not blocking anything.
