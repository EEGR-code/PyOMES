# Solver Promotion — Phase 6 design note

> **Status: Shipped 2026-05-01.** The promotion described below was
> implemented as Phase 6 of the [CV refactor](CV_UPDATE.md); see
> [PHASE6_CHECKLIST.md](PHASE6_CHECKLIST.md) for the implementation
> log.  This file is the original design note, kept as a historical
> record.

## What this is

A short planning note capturing the idea of promoting
`EulerSnapshotSolver` and `ScipyODESolver` from being
`GasLiquidVolume`-attached strategies to being something a
`ControlVolume` can use directly via
`cv.advance(solver=...)`.

This was the smallest concrete piece of the unbundling described in
[CONTAINER_LAYERING.md](../upcoming/CONTAINER_LAYERING.md) and
the natural follow-on to the Phase 5 work captured in
[CV_UPDATE.md](CV_UPDATE.md).  Phase 6 used this note as the seed of
its checklist (similar in shape to
[PHASE5_CHECKLIST.md](PHASE5_CHECKLIST.md)).

The original suggestion was flagged in
[ORDERING_CRITIQUE.md](ORDERING_CRITIQUE.md) Review 3 ("Consider
making the strategy pluggable").

---

## What changes

Today the `StepSolver` protocol takes a `GasLiquidVolume` and
reaches through that wrapper for state. After this promotion the
solver protocol takes a `ControlVolume` directly. The CV grows a
`solver=` parameter on `advance()` that defaults to the existing
sequential operator-split behaviour but can be set to
`EulerSnapshotSolver` or `ScipyODESolver`.

Rough endpoint:

```python
# Default (sequential operator splitting — current cv.advance() behaviour):
result = cv.advance(dt_h, chem_env)

# Snapshot-Euler with proportional clamping:
result = cv.advance(dt_h, chem_env, solver=EulerSnapshotSolver())

# Adaptive ODE with implicit method:
result = cv.advance(dt_h, chem_env, solver=ScipyODESolver(method="Radau"))
```

`GasLiquidVolume` shrinks further: it stops owning solver dispatch
and becomes essentially a boundaries holder + result cache + thin
forwarding `advance()`. The boundaries layer gets handled inside
whichever solver the caller picks (the snapshot solver already
snapshots boundaries; the sequential default would either keep
boundaries unaware, or grow a boundaries pass — design choice).

---

## The design question — CV-level or Unit-level?

Where the solver dispatches:

| Option | Solver attaches to | Scope solver sees |
|---|---|---|
| A | `ControlVolume` | One CV (its phases + internal interfaces + boundaries) |
| B | A merged `Vessel`/`Unit` class | The unit's full topology (1+ CVs + links + boundaries) |
| C | The topology object directly (CV for single-CV, multi-CV graph for N-CV, 1D field for HPLC) | Whatever the topology is |

The principled answer is that **a solver operates at whatever scope
owns the topology being integrated**. A Unit-level attachment is
fine *if* the Unit owns the topology, but in that case the Unit is
just forwarding `advance()` to its topology's solver — the solver
itself is still operating on the topology, not on the unit.

Options A and C converge for the single-CV case (the CV *is* the
topology). Option B is the right *packaging* once a Unit class
exists that owns the topology, but it isn't orthogonal to A — Unit
just dispatches through to whichever scope holds the topology.

---

## Recommendation

**Do option A first** — promote the solvers to operate on a
`ControlVolume` via `cv.advance(solver=...)`. Reasons:

- **Smallest contained piece of work.** The current solvers are
  already single-CV after Phase 5; the promotion is essentially
  "stop pretending `GasLiquidVolume` is special; let any CV use
  these solvers."
- **Immediately useful.** Any single-CV user — including
  post-Phase-5 GLV, single-phase liquid CVs with feeds and
  reactions, multi-phase CVs without GLV's gas-liquid mould — can
  pick a solver without a wrapper class.
- **A foundation, not a commitment.** When (if) a Unit class merges
  `GasLiquidVolume` and `MultiCVSystem`, it will hold a topology and
  forward `advance()` through whichever solver the topology
  supports. The CV-level protocol is exactly what that forward call
  invokes.
- **Doesn't pre-commit to the rename.** This work can be done before
  any decision on `GasLiquidVolume` → `Vessel`/`Unit` naming.

**Then, if and when a use case appears, generalise to multi-CV
topologies.** This is genuinely new capability — "snapshot Euler
across a multi-zone system" doesn't exist anywhere today. Options:

- A separate solver tier that operates on a multi-CV topology
  object (snapshots all CVs, evaluates boundaries + internal
  interfaces + inter-CV links + reactions from frozen state, applies
  combined deltas).
- A single `StepSolver` protocol that takes either a CV or a
  multi-CV topology, branching internally.

Either is workable; pick when the use case is concrete enough to
inform the choice.

---

## What this covers

After option A, the four common shapes are all handled cleanly:

| Shape | API after option A |
|---|---|
| Single CV, single phase | `cv.advance(dt_h, solver=...)` |
| Single CV, mixed phase (post-Phase-5 GLV pattern) | `cv.advance(dt_h, solver=...)` |
| Multi CV, single phase per CV | `MultiCVSystem.advance_all(dt_h, solvers={...})` — per-CV solver choice |
| Multi CV, mixed phase per CV | Same as above |

What option A does **not** cover:

- Snapshot integration across multiple CVs as one Euler step. That
  requires the multi-CV-aware solver tier described above. Gated on
  the trigger conditions in
  [CONTAINER_LAYERING.md](../upcoming/CONTAINER_LAYERING.md).

---

## Triggers — when to actually do this

The CV-level promotion (option A) is worth doing when one of these
appears:

- A model that wants snapshot or adaptive integration but doesn't
  fit the GLV gas+liquid mould — e.g. a single-phase liquid CV with
  feeds and reactions wanting clamped Euler, or a custom-phase CV
  wanting `solve_ivp`.
- A wish to slim `GasLiquidVolume` further (it shrinks to little
  more than a boundaries layer + result cache + thin forwarder).
- A planned change to `MultiCVSystem` that wants per-CV solver
  choice (different zones using different integration strategies).
- Cleanup motivation — the GLV-shaped tangling of solvers is the
  last big thing the Phase 5 work didn't fix.

The multi-CV-aware solver tier is gated separately on a multi-zone
model that needs snapshot semantics across zones — see
[CONTAINER_LAYERING.md](../upcoming/CONTAINER_LAYERING.md).

---

## Scope estimate

Roughly Phase-5-sized:

- Extend the `StepSolver` protocol so `solve_step` takes a
  `ControlVolume` (or accepts both for a transition period — though
  the established pattern in this codebase is no shims, so probably
  switch outright).
- `ControlVolume.advance(dt_h, chem_env, ..., solver=None)` — when
  `solver` is not `None`, dispatch to it; otherwise run the existing
  sequential body. The sequential body could itself be packaged as
  a `SequentialAdvanceSolver` for symmetry.
- Boundaries: decide whether they remain a GLV-level concept (so
  `cv.advance(solver=...)` doesn't see them) or become a CV-level
  concept too. The latter is cleaner long-term but a larger change.
- `GasLiquidVolume.advance` becomes a thin forwarder:
  ```python
  def advance(self, dt_h, chem_env=None, ...):
      result = self._cv.advance(dt_h, chem_env, ..., solver=self._solver)
      ...
      return wrap_in_GasLiquidAdvanceResult(result)
  ```
- Test updates: ~10s of edits, similar surface to Phase 5 changes
  in `tests/standalone/`.
- Test for the new pattern: a non-GLV single-CV using
  `EulerSnapshotSolver`. Currently impossible; this proves the
  abstraction holds.

One PR is plausible, with a checkpoint structure modelled on
[PHASE5_CHECKLIST.md](PHASE5_CHECKLIST.md).

---

## Open design questions to resolve before starting

When this becomes active work, these need decisions:

1. **Boundaries on CV or on GLV?** Today they're GLV-only. Moving
   them to CV is cleaner but is a separate piece of scope.
   Realistic answer: leave them on GLV for now, design the CV-level
   solver protocol to accept an optional boundaries list as a
   parameter, and let GLV pass them through.
2. **Should the existing sequential `cv.advance()` body become a
   `SequentialAdvanceSolver` strategy, or stay inline as the default
   when `solver=None`?** Symmetry argues for the former; pragmatism
   argues for the latter (less churn).
3. **Does the new protocol take `dt_h` or a time interval?** Current
   `StepSolver.solve_step` takes `dt_h`; CV's `advance` takes
   `dt_h`. Aligned. No change.
4. **What does `compute_reaction_rates` (the read-only sibling of
   `advance`) do under a snapshot solver?** Currently it builds a
   `ReactionEnvironment` directly. Probably untouched — it's not
   an integration path.

---

## Relationship to other phases

- [CV_UPDATE.md](CV_UPDATE.md) — Phases 1–5 shipped 2026-04-29. Made
  GLV a thin single-CV wrapper, which makes this promotion natural.
- [PHASE5_CHECKLIST.md](PHASE5_CHECKLIST.md) — implementation log
  for Phase 5, the pattern this future phase would follow.
- [CONTAINER_LAYERING.md](../upcoming/CONTAINER_LAYERING.md) — the larger
  unbundling. This promotion is its first slice; merging GLV +
  MultiCVSystem into a `Unit` class is a separate, larger slice
  gated on different triggers.
- [docs/solvers.md](../solvers.md) — characterises the three
  integration paths whose home this work is changing. Will need a
  small update once promotion lands (the "Future direction" section
  describes this very work).
- [ORDERING_CRITIQUE.md](ORDERING_CRITIQUE.md) Review 3 — the
  original "splitting strategy pluggable" suggestion this note is
  acting on.
