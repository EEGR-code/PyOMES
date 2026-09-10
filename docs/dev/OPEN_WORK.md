# Open Work — 2026-06-10

> **Stale as of 2026-07-09.** This snapshot predates `partition-model`,
> the full solver-architecture A–E track (all shipped 2026-06-11/12,
> contradicting the "trigger-gated" status still written below),
> `chemistry-unification-3b`, `THERMODYNAMIC_MODEL_ARCHITECTURE`,
> `CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE`, `LAYER1_GAP_CLOSURE`, and
> `step-solver-interface-refinement`. **Current source of truth:**
> [`phases-upcoming/README.md`](phases-upcoming/README.md)'s "Recently
> shipped" list and the individual design docs in
> [`design/`](design/) and [`phases-upcoming/`](phases-upcoming/).
> Left in place as a historical snapshot rather than deleted or
> rewritten line-by-line; verify anything below against `git log`/code
> before trusting it, per the pattern that already caught the
> HPC_CHECKPOINTING entry being wrong (see `phases-upcoming/README.md`'s
> corrected entry).

Summary of trigger-gated and open work items after
`run-history` shipped.  Current suite:
**1226/0** standalone tests.

---

## Recently shipped

| Phase | Tag | Date |
|---|---|---|
| `cufermenter-sunset` | `cufermenter-sunset-shipped` | 2026-06-01 |
| `framework-polish` | `framework-polish-shipped` | 2026-06-01 |
| `result-export` | `result-export-shipped` | 2026-06-01 |
| `speciation-level-retirement` | `speciation-level-retirement-shipped` | 2026-06-03 |
| `chemistry-unification-3b` | `chemistry-unification-3b-shipped` | 2026-06-03 |
| `property-snapshot-phase-agnostic` | `property-snapshot-phase-agnostic-shipped` | 2026-06-04 |
| `param-path-dispatcher` | `param-path-dispatcher-shipped` | 2026-06-08 |
| `run-history` | `run-history-shipped` | 2026-06-08 |
| `partition-model` | `partition-model-shipped` | 2026-06-16 |

Design docs archived in [`../phases-shipped/`](../phases-shipped/).

---

## Ready to implement

### HPC_CHECKPOINTING
**Doc:** [`HPC_CHECKPOINTING.md`](../phases-shipped/HPC_CHECKPOINTING.md)

**Trigger:** A batch job where mid-run checkpoint/restart is needed
(e.g. HPC cluster with walltime limits, or a long optimisation sweep
that may be interrupted).

**Scope:** Parquet-based checkpoint writes at configurable intervals;
`Simulation.resume_from(path)` factory; compatible with the existing
`StreamingFileRecorder` / `SparseRecorder` infrastructure added in
RESULT_EXPORT.

---

## Solver architecture track

Six sequenced phases forming the two-axis integration architecture.
See [`SOLVER_ARCHITECTURE.md`](SOLVER_ARCHITECTURE.md) for full context.

Phases A and B are independent of each other; C requires both.
D and E both require C but are independent of each other.
F and G are SUNDIALS opt-in; no work scheduled pending concrete trigger.

| Phase | Doc | Axis | Prerequisites | Trigger |
|---|---|---|---|---|
| A — CV_COMPUTE_INTERFACE | [`CV_COMPUTE_INTERFACE.md`](CV_COMPUTE_INTERFACE.md) | 1+2 foundation | — | Start of any solver work |
| B — CONTROLLER_STATE_PROTOCOL | [`CONTROLLER_STATE_PROTOCOL.md`](CONTROLLER_STATE_PROTOCOL.md) | 2 | — | Controller with integral state + tight loop |
| C — SYSTEM_SOLVER_PROTOCOL | [`SYSTEM_SOLVER_PROTOCOL.md`](SYSTEM_SOLVER_PROTOCOL.md) | 2 | A, B | Multi-CV model needing better than Euler |
| D — IMPLICIT_TRANSPORT | [`IMPLICIT_TRANSPORT.md`](IMPLICIT_TRANSPORT.md) | 2 | C | CFL subcycle M ≳ 10 (HPLC, fast-circulation fermenter) |
| E — MONOLITHIC_ODE | [`MONOLITHIC_ODE.md`](MONOLITHIC_ODE.md) | 2 | C | T_c < dt_h (tight pH/DO control, multi-rate digital controller) |
| F — PER_CV_DAE *(opt-in)* | see SOLVER_ARCHITECTURE.md | 1 | A | Operator-splitting error demonstrably non-negligible |
| G — SYSTEM_DAE *(opt-in)* | see SOLVER_ARCHITECTURE.md | 1+2 | E, F | Phase F in use + coupled multi-CV speciation bottleneck |

---

## Open (no trigger gate)

### `build_bsm2_cv` dual-solver wiring
No dedicated doc.

`build_bsm2_cv` (in `models/vlmodels/adm1/bsm2.py`) previously wired:
- a `from_reactions` `BisectionChemicalEquilibriumEngine` (correct, production path), and
- a legacy `Level 1` solver via `FermenterBuilder.chemistry(
  speciation_level=1)`.

Both were attached to the same CV — a leftover from the incremental
chemistry-unification migration.  This was targeted for clean-up in
SPECIATION_LEVEL_RETIREMENT (shipped 2026-06-03).  Verify that the
legacy wiring was actually removed in that phase; if not, it is now
an unblocked standalone clean-up.

---

## Chemistry unification — complete

All six branches shipped:

| Branch | Tag | Date |
|---|---|---|
| `chemistry-unification-1` | `chemistry-unification-1-shipped` | 2026-05-13 |
| `chemistry-unification-2` | `chemistry-unification-2-shipped` | 2026-05-15 |
| `chemistry-unification-3` | `chemistry-unification-3-shipped` | 2026-05-15 |
| `chemistry-unification-4` | `chemistry-unification-4-shipped` | 2026-05-16 |
| `chemistry-unification-3b` | `chemistry-unification-3b-shipped` | 2026-06-03 |

Design doc and plan archived in
[`../phases-shipped/`](../phases-shipped/).
