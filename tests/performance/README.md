# Performance

Runtime/timing comparisons (e.g. engine-vs-engine or
activity-model-vs-activity-model speed benchmarks). Distinct from
[`tests/validation/`](../validation/), which checks numerical *accuracy*
against analytical/reference solutions — this directory is for *speed*
claims only.

- `bench_adaptive_solver.py`: `SimultaneousAdaptiveSolver` on iron oxidation, BSM2
  and a weak-acid model, with controllers, a profile and a feed and drain. It gates
  on engine solves and error against a stored reference, or on a bit-level
  fingerprint with `--exact`, and reports CPU time without gating on it. Stored
  numbers are in `baselines/adaptive_solver.json`. The fingerprints depend on the
  machine and library versions; the references and solve counts are the portable part.

Measurement scripts for the `conservation-ledger` phase
(`docs/dev/implementation/upcoming/CONSERVATION_LEDGER_CHECKLIST.md`). They
measure results rather than speed, and live here because they are run by hand
like the benchmark:

- `ledger_fingerprint.py`: runs the phase's fixed cases and writes every
  amount, pH and warning to a JSON file (`capture`), and compares two such
  files bit for bit (`compare`). Output files are not committed.

`demos/usecases/04_compare_runtime_by_usecase.ipynb`, which used to hold
engine-vs-engine comparisons, was deleted rather than migrated (see
`docs/dev/implementation/upcoming/TUTORIALS_REORG_CHECKLIST.md`).
