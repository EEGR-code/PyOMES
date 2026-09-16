# Results

Post-processing a `Simulation`'s output — the concern is orthogonal to
every other tutorial folder (chemistry, topology, solver mechanics): once
you have a `BatchResult`, what do you do with it?

| Notebook | Covers |
|---|---|
| [01_exporting_results.ipynb](01_exporting_results.ipynb) | The three export helpers `BatchResult` carries: a tidy long-form `DataFrame` (`.to_dataframe()`), a pivoted wide-form table for plotting (`.to_wide()`), and file export (`.to_csv()`/`.to_parquet()`) with round-trip precision guarantees. |

Needs `pandas` and `pyarrow`: `pip install -e ".[export]"` (the project's
existing extras group for this) — the only tutorial content with a
non-core dependency besides the optional PHREEQC notebooks.

## Running

```bash
jupyter lab docs/tutorials/results/
```
