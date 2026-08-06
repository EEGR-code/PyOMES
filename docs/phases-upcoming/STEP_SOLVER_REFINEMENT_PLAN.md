# Step-Solver Interface Refinement — Implementation Plan

## Status

> **Stage 0 and Stage 1 shipped 2026-07-08.** See the Sequencing table
> below for tags; [STEP_SOLVER_REFINEMENT_CHECKLIST.md](STEP_SOLVER_REFINEMENT_CHECKLIST.md)
> for Stage 1's full implementation log. Stages 2–3 remain as scoped
> below — Stage 2 (Q7 design decision) has not been resolved, so Stage
> 4 (SIA) is not yet ready to checklist.

Implementation plan for the work described in
[../phases-shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md](../phases-shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md)
(shipped; moved from `phases-upcoming/` when Stage 1 shipped) and the
"interleaving" thread in
[`../design/MASS_EXCHANGE_ARCHITECTURE.md`](../design/MASS_EXCHANGE_ARCHITECTURE.md)
§10–13 and [`../design/SOLVER_ARCHITECTURE.md`](../design/SOLVER_ARCHITECTURE.md)'s
"Identified extensions". Those two docs are the source of truth for
*why*; this plan pins the **slicing**, the **branch/tag layout**, and
what's ready to checklist now versus what still needs a design
decision first.

This plan was produced by a scoping review (2026-07-06) that
re-verified every claim in the two source docs against the current
code (not just trusted the prose) — see the *Verification notes*
section below for what was confirmed and the one drift found.

## Sequencing

| # | Name | Branch / tag | Content | Readiness |
|---|---|---|---|---|
| 0 | Merge rename | *(no branch — finished `rename-simultaneous-step-solvers`)* | `EulerSnapshotSolver`→`SimultaneousEulerSolver`, `ScipyODESolver`→`SimultaneousAdaptiveSolver`. | **Shipped 2026-07-06** — tag `rename-simultaneous-step-solvers-shipped`. Turned out to be entirely uncommitted working-tree state (not actually merged as the design note's banner implied) — split into a rename commit + one unrelated doc-fix commit before merging; see checklist checkpoint 0. |
| 1 | Step-solver interface refinement | `step-solver-interface-refinement` / `step-solver-interface-refinement-shipped` | Items 1, 2, 3, 4, 6, 7, 8, 9 from `STEP_SOLVER_INTERFACE_REFINEMENT.md` (ownership guard, Monolithic validation, `SequentialAdvanceSolver` reification, state-vector unification, gas/liquid generalization, shared clamp module, clamp diagnostics, demo notebook). | **Shipped 2026-07-08** — tag `step-solver-interface-refinement-shipped`. 10 checkpoints (0–9, including 7b split out from 4); full log in [STEP_SOLVER_REFINEMENT_CHECKLIST.md](STEP_SOLVER_REFINEMENT_CHECKLIST.md). 1950 → 2003 tests. |
| 2 | Z-staleness decision (Q7) | *(no branch — design discussion only)* | Resolve `MASS_EXCHANGE_ARCHITECTURE.md` §12 Q7: does `SequentialIterativeSystemSolver` re-solve speciation (z) each iteration, or reuse the first pass? | **Blocks Stage 4.** Not code — a design conversation to have before opening the SIA branch. |
| 3 | Reactive D_eff transport | `reactive-deff-transport` / `reactive-deff-transport-shipped` | `DispersiveFlow.compute_flow()` gains an optional engine reference; when the engine satisfies `GrayBoxEngineProtocol`, compute `D_eff,CT = Σ(∂z_i/∂CT)·D_i` via the already-shipped `jacobian_dz_dy()` before assembling flux. | **Not blocked** — `MonolithicODESolver` (its stated prerequisite) already shipped 2026-06-12. Scoped in `SOLVER_ARCHITECTURE.md` "Identified extensions"; not yet checklisted at checkpoint granularity. |
| 4 | SIA (`SequentialIterativeSystemSolver`) | `sequential-iterative-system-solver` / `...-shipped` | New `SystemSolver` iterating `(transport → cv.advance())` to convergence. | **Blocked on Stage 2.** Scoped in `SOLVER_ARCHITECTURE.md`; not checklisted — the checkpoint list depends on the Q7 answer (whether a re-solve hook needs adding to `cv.advance()`/`compute_rhs()` as part of the same phase). |

Item 5 from `STEP_SOLVER_INTERFACE_REFINEMENT.md` (composition /
integration-method / clamp_fn as three independently swappable
dimensions) is **explicitly deferred**, per the design note's own
recommendation — it's a bigger lift that should only be built when a
concrete model needs a (composition, method) pair the two shipped
`StepSolver`s can't express, not built reflexively. Not scheduled in
this plan.

Recorder demos (`DEMO_RECORDERS.md`), precipitation CV integration
(`NR_PRECIPITATION_CV_INTEGRATION.md`), and the multi-component
complexation engine (`MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md`)
are unrelated tracks — no file overlap with Stages 1/3/4 except
`control_volume.py`, and there only at different methods (`__init__`
vs `advance()`/`step_internal_transfer()`). No sequencing dependency
either direction; mentioned here only so a future reader doesn't
wonder why they're absent from this plan.

## Why this order

1. **Stage 0 first, unconditionally.** It's already done and tested;
   leaving it unmerged risks a conflicting fork the moment Stage 1
   opens a branch touching the same files (`solvers.py`,
   `system_solver.py`, `control_volume.py`).
2. **Stage 1 before Stages 3/4.** Items 1 and 2 in the bundle are
   live, confirmed-in-code correctness gaps today (`cv.advance()` has
   no ownership guard; `MonolithicODESolver` silently discards
   `solver=`) — not hypothetical future needs. Everything else in the
   bundle is independent, low-risk, and explicitly recommended to ship
   together by the source design note.
3. **Stage 2 is a prerequisite for Stage 4, not a phase in its own
   right.** SIA's whole value proposition is O(dt²) convergence from
   iterating `(transport → cv.advance())`; if z is reused stale across
   iterations (one candidate answer to Q7), iterating only converges
   *y*, not the true simultaneous solution the note promises — so the
   answer to Q7 changes what SIA's checkpoint list even needs to
   contain (possibly a new re-solve hook on `cv.advance()`). Opening
   the SIA branch before this is resolved risks discovering the
   re-solve requirement mid-implementation.
4. **Stage 3 does not need to wait for Stage 2.** The docs present
   D_eff and SIA side-by-side as "the interleaving thread," which
   reads as if they're coupled — they aren't. D_eff's stated
   prerequisite is "MonolithicODE **or** SIA iteration to be
   consistent," and MonolithicODE already shipped (`monolithic-ode-shipped`,
   2026-06-12). D_eff can be built and used today against the
   Monolithic path without SIA existing at all. It's the
   better-scoped of the two "interleaving" items if that thread is
   picked up before Q7 is settled.

## Verification notes (2026-07-06 scoping review)

Confirmed directly against code before this plan was written (not
assumed from the docs):

- `cv._context` exists ([`control_volume.py:312`](../../src/core/control_volume.py#L312))
  but `advance()` ([`control_volume.py:516`](../../src/core/control_volume.py#L516))
  has no check against it — item 1's gap is real.
- `MonolithicODESolver.advance_system()` ([`system_solver.py:858`](../../src/core/system_solver.py#L858))
  never calls `sim._solver_for`/reads `sim.solver` — item 2's gap is
  real. Contrast with the four other `SystemSolver`s, all of which do
  ([`system_solver.py:243`](../../src/core/system_solver.py#L243) inside
  `StrangSplittingSystemSolver`, [`:346`](../../src/core/system_solver.py#L346)
  inside `MultirateSystemSolver`, [`:557`](../../src/core/system_solver.py#L557)
  inside `ImplicitTransportSystemSolver`, and
  [`simulation.py:859`](../../src/core/simulation.py#L859) inside
  `Simulation._step_default` — the path `ExplicitEulerSystemSolver`
  transparently forwards to).
- `_build_rhs()` ([`system_solver.py:714-736`](../../src/core/system_solver.py#L714-L736))
  moves each species independently at its own `kLa`/`Q`, no gray-box
  Jacobian consulted — the D_eff gap (§10.3) is real and unchanged
  since `MASS_EXCHANGE_ARCHITECTURE.md` §10.4 was last corrected
  (2026-07-06).
- No `SequentialIterativeSystemSolver` exists anywhere in
  `system_solver.py` — SIA is genuinely unbuilt, matching the docs.
- **One drift found:** `MASS_EXCHANGE_ARCHITECTURE.md` §13.3/§13.4
  draft an `EquilibriumScope` dataclass for the engine-scope
  declaration question, but §12 Q8 records that type as superseded
  2026-06-30 — `algebraic_species()` is the sole ownership mechanism
  now. The underlying question in §13.3/§13.4 (should scope-overlap
  between engine-owned species and registered `PhaseInterface`s be
  validated at construction time?) is still open and unaffected by
  this plan — it's a `MASS_EXCHANGE_ARCHITECTURE.md` §11 priority-2
  item, not part of the step-solver track, and not scheduled here —
  but any future implementation of it should key off
  `algebraic_species()`, not the dropped `EquilibriumScope` type.
- Additionally discovered while tracing item 1's mechanism: **all
  four `SystemSolver`s other than `MonolithicODESolver` call the
  public `cv.advance(dt_h, t_h, solver=solver)`**, not a trusted
  internal entry point. Once the ownership guard lands, all four call
  sites listed above must switch to calling `cv._advance_unchecked(...)`
  directly, or every orchestrated `Simulation.run()` will spuriously
  emit `OrchestrationWarning` on every step. This is folded into
  Stage 1's checklist as its own checkpoint — it isn't optional
  polish, it's required for the guard to work at all.

## What ships where — detail

- **Stage 1**: see [STEP_SOLVER_REFINEMENT_CHECKLIST.md](STEP_SOLVER_REFINEMENT_CHECKLIST.md)
  for the full checkpoint-by-checkpoint breakdown.
- **Stage 3 (D_eff)**: scope transcribed from `SOLVER_ARCHITECTURE.md`'s
  "Identified extensions" — extend `DispersiveFlow.compute_flow()` to
  accept an optional engine reference; when the engine satisfies
  `GrayBoxEngineProtocol`, compute `D_eff` per conserved total before
  assembling the flux dict; default behaviour (per-species `kLa`,
  no `D_eff`) unchanged unless the flag/engine-reference is present.
  Not checklisted at checkpoint granularity yet — write
  `REACTIVE_DEFF_TRANSPORT_CHECKLIST.md` when this stage opens,
  modelled on this plan's Stage 1 checklist.
- **Stage 4 (SIA)**: scope transcribed from `SOLVER_ARCHITECTURE.md`
  — new `SystemSolver` wrapping any inner `SystemSolver`, convergence
  loop over `(transport → cv.advance())`, tolerance/max-iteration
  parameters, compatible with any per-CV `StepSolver`. Checklist
  cannot be written with full confidence until Stage 2 resolves
  whether an explicit re-solve hook is needed on `cv.advance()`/
  `compute_rhs()` for iteration N>1 to be meaningfully different from
  N=1.

## How Stage 1 was started and shipped (historical — both stages done)

Followed [README.md](README.md)'s branching convention:

1. `git checkout main && git merge --no-ff rename-simultaneous-step-solvers -m "Merge rename-simultaneous-step-solvers: SimultaneousEulerSolver/SimultaneousAdaptiveSolver rename"`
2. `git tag rename-simultaneous-step-solvers-shipped <commit-hash>`
3. `git push && git push --tags`
4. `git branch -d rename-simultaneous-step-solvers && git push origin --delete rename-simultaneous-step-solvers`
5. `git checkout -b step-solver-interface-refinement` off the now-updated `main`
6. Worked through [STEP_SOLVER_REFINEMENT_CHECKLIST.md](STEP_SOLVER_REFINEMENT_CHECKLIST.md)
   in order (10 checkpoints, one commit each), then shipped the same way:
   merge `--no-ff` to `main`, tag `step-solver-interface-refinement-shipped`,
   push, delete branch.

## How to start Stage 2 (next, when picked up)

Stage 2 is a design discussion, not a branch — resolve
`MASS_EXCHANGE_ARCHITECTURE.md` §12 Q7 first (see the Sequencing table
above), then either open `reactive-deff-transport` (Stage 3, unblocked
already) or `sequential-iterative-system-solver` (Stage 4, once Q7 is
resolved), following the same checklist-first convention.
