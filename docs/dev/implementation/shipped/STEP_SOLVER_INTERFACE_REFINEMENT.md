# Step-Solver Interface Refinement — Planning Note

> **Status: Shipped 2026-07-06.** Item 10 (rename) shipped separately as
> Stage 0 (`rename-simultaneous-step-solvers-shipped`). Items 1, 2, 3, 4, 6,
> 7, 8, 9 shipped as Stage 1 of
> [`STEP_SOLVER_REFINEMENT_PLAN.md`](../upcoming/STEP_SOLVER_REFINEMENT_PLAN.md)
> on branch `step-solver-interface-refinement`
> (tag `step-solver-interface-refinement-shipped`) — see
> [`STEP_SOLVER_REFINEMENT_CHECKLIST.md`](../upcoming/STEP_SOLVER_REFINEMENT_CHECKLIST.md)
> for the full implementation log, including several scope corrections found
> during implementation (notably: item 6's generalization split into two
> checkpoints, 4 and 7b, because `SimultaneousAdaptiveSolver`'s ODE hot loop
> turned out to be far more deeply hardcoded to gas/liquid than "mechanical"
> implied). **Item 5** (composition/integrator/clamp_fn three-dimension
> factoring) remains explicitly deferred, per this note's own recommendation.
> The deferred interleaving/multi-CV-aware solver tier discussion below is
> Stages 2–4 of the parent plan, not part of this shipped stage.

## What this is

A consolidated record of a multi-round design review of VLsim's per-CV /
per-system time-stepping interfaces — `StepSolver` (Axis 1) and `SystemSolver`
(Axis 2) — triggered by auditing exactly how customizable `cv.advance()` /
`Simulation` solver selection really is today. This extends the shipped
two-axis design (`SOLVER_ARCHITECTURE.md`, Phases A–E) with a set of concrete
follow-ups surfaced by walking through real customization scenarios: bespoke
operator-splitting order, reproducing an external reference model's numerics,
per-species clamping control, and an unguarded footgun in `cv.advance()`.

This is **not** the SUNDIALS/DAE track (Phase F/G placeholders in
`SOLVER_ARCHITECTURE.md`) — everything here is about refining the existing
pure-scipy Axis 1/2 implementations and their extension surface, independent
of whether DAE solvers ever get built.

**Next planned session:** continue exploring the interleaving / multi-CV-aware
solver tier (see [Deferred: interleaving](#deferred-the-interleaving--multi-cv-aware-solver-tier)
below) — flagged explicitly by the user as the next thread to pick up.

---

## Background — current architecture (brief recap)

Two orthogonal `typing.Protocol`s:

- **`StepSolver`** (`src/core/solvers.py`) — Axis 1, per-CV physics. One method,
  `solve_step(cv, dt_h, t_h, external_source_terms) -> AdvanceResult`. Three
  points on this axis today: the inline sequential operator-split default,
  `SimultaneousEulerSolver` (simultaneous snapshot + proportional clamp), and
  `SimultaneousAdaptiveSolver` (continuous RHS via `solve_ivp`, adaptive/implicit methods).
  Selected via `cv.advance(dt_h, t_h, solver=...)`.
- **`SystemSolver`** (`src/core/system_solver.py`) — Axis 2, whole-system
  integration. Five shipped implementations from `ExplicitEulerSystemSolver`
  (transparent wrapper) up to `MonolithicODESolver` (full co-integration of CV
  species + controller state + event scheduling). Selected via
  `Simulation(system_solver=...)`.
- Per-CV Axis-1 choice composes with Axis-2 choice via
  `Simulation(solver=...)` + `Simulation._solver_for(cv_key)` — a single
  `StepSolver`, a `Dict[cv_key, StepSolver]`, or `None`.

Full detail of the original investigation (three shipped `StepSolver` rungs,
five shipped `SystemSolver` rungs, the DAE Phase F/G placeholders) is in the
artifact produced at the start of this review; the items below are what came
out of pressure-testing that architecture in conversation.

---

## Resolved findings and recommended changes

### 1. `cv.advance()` is an unguarded footgun for CVs owned by a `Simulation`

**Finding (confirmed in code, not assumed):** `ControlVolume.advance()` has no
`is_running` or `_context` gate. Nothing stops calling `cv.advance(dt_h, t_h)`
directly on a CV that is already wired into a `Simulation` with inter-CV
links/controllers/profiles. It executes, returns a well-formed
`AdvanceResult`, and silently treats the CV as isolated — no inter-CV mass
transfer, no controller action, no profile application that step. This is not
duplicated logic (every `SystemSolver` except `MonolithicODESolver` correctly
delegates to `cv.advance()` as a subroutine) — it's a duplicated *entry point*
that produces a complete-looking but physically wrong answer if misused.

**Resolution:** guard the public `advance()` method, but **warn, don't raise**.
A hard raise would block a legitimate use case (deliberately stepping one CV
in isolation as a diagnostic, ignoring its coupling, e.g. "what would this
CV's kinetics alone look like"). A warning is consistent with the precedent
this codebase already has for exactly this class of problem —
`AccuracyMonitor`/`ConservationMonitor` emit `UserWarning` subclasses by
default, throttled via `VLsim.config.warnings`/`VLSIM_WARNINGS`, and
promotable to a hard stop via `warnings.simplefilter("error", VLsim.SomeWarning)`.
Follow the same shape (new category, tentatively `OrchestrationWarning`).

**Mechanism — public/internal split**, mirroring the "Pattern B unchecked
setter" convention already used elsewhere in this codebase (public setter
gated, paired `_set_unchecked` for trusted internal callers):

```python
class ControlVolume:
    def advance(self, dt_h, t_h=0.0, *, external_source_terms=None, solver=None):
        if self._context is not None:
            warnings.warn(
                "cv.advance() called directly on a CV owned by a Simulation; "
                "inter-CV links/controllers/profiles will NOT be applied this "
                "step. Use sim.run(...) unless this is intentional.",
                OrchestrationWarning, stacklevel=2,
            )
        return self._advance_unchecked(dt_h, t_h, external_source_terms=external_source_terms, solver=solver)

    def _advance_unchecked(self, dt_h, t_h=0.0, *, external_source_terms=None, solver=None):
        solver = solver if solver is not None else SequentialAdvanceSolver()   # see item 3
        return solver.solve_step(self, dt_h, t_h, external_source_terms=external_source_terms)
```

`SystemSolver.advance_system()` implementations call `_advance_unchecked()`
directly (they are the trusted orchestrator) — the warning must **not** fire
on that internal path, only on direct user calls to the public `advance()`.

**Why not just delete `advance()` and force everything through `Simulation`?**
This was explicitly considered and rejected. `compute_rhs()` /
`compute_differential_rhs()` (the "black box" RHS surface used by
`MonolithicODESolver`) and `StepSolver.solve_step()` (used by
`SimultaneousEulerSolver`/the sequential default) are **two different integration
paradigms, not two spellings of the same thing**:

- `compute_rhs()` returns a continuous derivative (`dn/dt`) for a generic ODE
  integrator to consume.
- `SimultaneousEulerSolver`'s proportional clamping is inherently a **discrete**
  operation — it cannot be expressed as a continuous derivative without
  reintroducing exactly the "fake infinite rates" pathology
  `SOLVER_ARCHITECTURE.md` explicitly rejects as worse than operator splitting
  ("reintroduces the stiffness the algebraic solve avoids").

So collapsing everything to `compute_rhs`-only would be a real regression for
the discrete-step family, not a simplification. `advance()`/`StepSolver` stays
as the discrete-step interface; `compute_rhs`/`compute_differential_rhs` stays
as the continuous-RHS interface for Axis 2's `MonolithicODESolver`. The
resolution is the narrower ownership guard above, not removal.

**Practical idiom to adopt going forward:** `cv.advance()` remains available
and legitimate for **unowned** CVs — unit tests, notebooks exploring a bare
CV's behavior, and as the trusted internal building block Axis-2 solvers call.
For **owned** CVs, the recommended/default way to run anything — even a
single CV with zero links — is always through `Simulation(cvs={"only": cv}).run(...)`,
because that's where `t_h` bookkeeping, the recorder, and `RunContext`
lifecycle safety actually live.

### 2. `MonolithicODESolver` silently ignores per-CV `solver=` configuration

**Finding:** confirmed by reading `MonolithicODESolver.advance_system()` in
full — unlike the other four `SystemSolver`s, it never calls `cv.advance()`
and never calls `sim._solver_for(cv_key)`. It integrates every CV directly via
`cv.compute_rhs()` inside one shared `solve_ivp` call. So if a user configures
`Simulation(solver={"reactor": SimultaneousAdaptiveSolver(...)}, system_solver=MonolithicODESolver())`,
the `solver={...}` dict is read by nothing — **no error, no warning, it is
just silently moot.** This is the one place the "two axes are orthogonal"
framing in `SOLVER_ARCHITECTURE.md` doesn't hold, and it currently fails
silently rather than loudly.

**Resolution:** `MonolithicODESolver.advance_system()` (or `Simulation` at
`system_solver` assignment time) should warn — or raise, this is more
defensible as a hard error than the ownership-guard case above, since there's
no legitimate reason to pass a per-CV `solver=` that will provably never be
consulted — if `sim.solver` is non-`None` when `MonolithicODESolver` is the
active `system_solver`. Small, self-contained fix; no architectural change
needed, just a validation check at entry to `advance_system()`.

### 3. Reify the sequential default as `SequentialAdvanceSolver`

**Finding:** `solver=None` today runs ~70 lines of inline logic inside
`ControlVolume.advance()` that is not reachable as an object — it can't be
passed explicitly, can't be wrapped/decorated (e.g. a generic timing/logging
wrapper `timed(inner_solver)`), and can't be placed into a config-driven
`{"sequential": ???, "snapshot": SimultaneousEulerSolver()}` dispatch dict.

**Resolution:** promote it to a real `StepSolver`-conforming class,
`SequentialAdvanceSolver`, with `solver=None` becoming pure sugar for
`SequentialAdvanceSolver()` rather than a separate code path. This is not a
novel idea for this codebase — Axis 2 already solved the identical problem:
`ExplicitEulerSystemSolver.advance_system()` is a one-line transparent forward
to `sim._step_default()`, existing purely so the default is reachable through
the `SystemSolver` protocol (its own docstring: "validates that the
`SystemSolver` dispatch wiring is correct"). Axis 1 never got the equivalent —
do the same thing for symmetry.

```python
class SequentialAdvanceSolver:
    """The operator-split default, promoted to a StepSolver instance."""
    def solve_step(self, cv, dt_h, t_h=0.0, external_source_terms=None):
        ...  # today's inline body, moved here verbatim, unchanged behavior
        return AdvanceResult(...)
```

After this: `cv.advance(dt_h, t_h)` and
`cv.advance(dt_h, t_h, solver=SequentialAdvanceSolver())` are byte-identical
calls.

### 4. Unify state-vector packing/unpacking (three near-duplicate implementations)

**Finding:** the same species-ordering convention (CV → phase → species,
alphabetical) is implemented three times with slightly different scope:
`_StateVector` in `solvers.py` (single-CV, phase-scoped, has an
`algebraic_species` exclusion set), `_pack_state`/`_unpack_state`/`_state_index_map`
in `system_solver.py` (multi-CV, no exclusion set), and
`_pack_extended_state`/`_unpack_extended_state` (multi-CV + controller
differential-state extension, Phase E). Not urgent today, but risks drifting
further apart every time a new solver is added.

**Resolution (user confirmed this is worth doing):** one parameterized utility
— scope (single-CV or all-CVs), optional exclusion set, optional
controller-state extension — replacing all three call sites. Target module:
`src/core/state_vector.py`.

### 5. Three composable dimensions inside a `StepSolver` — currently bundled, not orthogonal

**Correct terminology** (standard numerical-methods vocabulary, introduced
during this review to disambiguate what was being conflated):

| Dimension | Term | Answers |
|---|---|---|
| How sub-processes combine into one step | **Operator-splitting scheme** (Lie-Trotter/sequential, Strang/symmetric, unsplit-simultaneous, fully-coupled DAE) | "In what order/combination do instantaneous reactions, kinetics, boundaries, transfer interact?" |
| How a (possibly split) system is numerically advanced | **Time-integration method** (explicit Euler, RK45/DOP853, implicit BDF/Radau, adaptive step control) | "Given the composition, what numerical scheme moves it forward by dt?" |
| Restoring a physical invariant the raw numerics don't guarantee | **Positivity-preserving projection** ("clamping") | "After integrating, does the result respect non-negativity/conservation?" |

**Finding:** Axis 2 already achieves clean separation between *composition*
(which `SystemSolver`) and *per-CV integration method* (which `StepSolver`,
and within `SimultaneousAdaptiveSolver`, its `method=` parameter) — confirmed in code,
`StrangSplittingSystemSolver`/`MultirateSystemSolver`/`ImplicitTransportSystemSolver`
all still call `sim._solver_for(cv_key)` and honor whatever per-CV `StepSolver`
is configured. But **within Axis 1 itself, composition and integration method
are not separated** — `SimultaneousEulerSolver` welds {simultaneous-snapshot
composition} to {explicit Euler}; `SimultaneousAdaptiveSolver` welds {continuous-RHS/unsplit
composition} to {your choice of DOP853/Radau/BDF/LSODA}, but only within that
one composition style. Only two (composition, method) *pairs* exist today, not
two independently selectable knobs — there's no way to ask for "simultaneous
snapshot composition + implicit BDF" without writing a wholly new class.

**Resolution (target, not immediate):** if this becomes a real need, factor
`StepSolver` into three independently swappable pieces:
`StepSolver(composition=SimultaneousSnapshot(), integrator=ImplicitEuler(),
clamp_fn=proportional_clamp)`. This is a bigger lift than anything else in
this list — it's inventing a small composition-strategy sub-framework, not
just reifying/generalizing existing classes — so it should only be pursued
if/when a concrete model needs a (composition, method) combination the two
shipped classes can't provide, not built reflexively.

### 6. Generalize `SimultaneousEulerSolver` / `SimultaneousAdaptiveSolver` off the gas/liquid assumption

**Finding (from the original investigation, carried forward):** both solvers
`raise ValueError` unless `cv.phases` contains exactly `"gas"` and `"liquid"`
— a leftover from their pre-promotion origin attached to the now-deleted
`GasLiquidVolume` wrapper (`SOLVER_PROMOTION.md`). The boundaries-handling loop
inside these solvers is *already* phase-agnostic
(`for boundary in cv.boundaries: ... boundary.phase_key`); only the snapshot
construction (`cv.phases["gas"]`/`cv.phases["liquid"]` direct lookups) and the
transfer-link search (assumes exactly one `KineticGasLiquidLink`) are
hardcoded.

**Resolution:** mechanical, not a redesign — snapshot `cv.phases` generically
(`{k: p.snapshot() for k, p in cv.phases.items()}`), iterate
`cv.internal_interfaces` generically instead of searching for one named link
type. Once done, any CV phase shape (single-phase liquid with feeds/reactions,
solids-bearing CVs, etc.) gets the same snapshot-Euler / adaptive-ODE options
a gas-liquid CV has today, with zero change to the caller-facing
`cv.advance(solver=...)` surface.

### 7. Shared clamping module — built to not foreclose per-species customization

**Finding:** non-negativity handling is inconsistent across all three current
solvers, and it is **not** solely a kinetics/ODE concern — checked which
mechanisms can actually overdraw inventory over a discrete `dt_h`: kinetic
reactions, external boundary fluxes (feed/vent/membrane), and inter-phase/
inter-CV transfer fluxes **all** can. `SimultaneousEulerSolver` already gets this
right structurally (its `_clamp_deltas()` operates on the combined delta dict
from boundaries + transfer + reactions + external source terms together) —
but the *technique* differs pointlessly across solvers: `SimultaneousEulerSolver`
does proportional scale-down (preserves relative stoichiometry);
the sequential default instead relies on `Phase.apply_flux`'s own blunt
per-species `max(0.0, …)` floor; `SimultaneousAdaptiveSolver` has no batch-delta clamping
at all, only a floor on unpack (a structurally different, narrower concern —
guarding against adaptive-integrator overshoot, not against a large discrete
`dt_h` overdrawing inventory for an otherwise well-posed rate law).

**Resolution:** pull `clamp_fn(deltas, current_mol, dt_h) -> deltas` out as a
shared, explicit, swappable constructor parameter for the two discrete-step
solvers (`SequentialAdvanceSolver`, `SimultaneousEulerSolver`) — default
`proportional_clamp`, with `floor_clamp`/`None` (disabled) as alternatives.
Keep `SimultaneousAdaptiveSolver`'s floor conceptually separate (different failure mode)
but equally named/swappable rather than a bare hardcoded `max(0.0, …)`.

**Critical design requirement raised by the user, and confirmed necessary:**
the shared module's building blocks must operate correctly on an **arbitrary
subset** of species, not assume whole-dict ownership — otherwise a shared
module *would* limit customizability, which is the opposite of the goal. This
matters concretely for reproducing an external reference model's (e.g.
BSM2/ADM1) specific per-species clamping convention, which this codebase needs
for cross-implementation comparability:

```python
def bsm2_style_clamp(deltas, current_mol, dt_h):
    """Match a specific reference model's per-species convention."""
    result = dict(deltas)
    for sp in NEVER_CLAMP:               # species the reference model lets go negative
        result[sp] = deltas[sp]
    for sp in FLOOR_AT_EPSILON:          # reference model uses eps, not exact zero
        result[sp] = _floor(deltas[sp], current_mol[sp], dt_h, eps=1e-9)
    remaining = {k: v for k, v in deltas.items() if k not in NEVER_CLAMP | FLOOR_AT_EPSILON}
    result.update(proportional_clamp(remaining, current_mol, dt_h))  # shared module's default, used as an ingredient
    return result

SimultaneousEulerSolver(clamp_fn=bsm2_style_clamp)
```

Because `clamp_fn` stays a **plain callable** (not a closed enum or a Protocol
class — deliberately, per the class-vs-function line-drawing discussed below)
a user can always write this kind of composite regardless of what's in the
shared module; the module just needs its own defaults to be usable as
ingredients, not monolithic take-it-or-leave-it functions.

### 8. Negative-mole check (disabled clamping) + "clamping was invoked" diagnostic

**Finding:** disabling clamping (`clamp_fn=None`) is a real way to silently
get unphysical (negative) state feeding into downstream chemistry, with
nothing currently watching for it.

**Resolution (both confirmed easy, additive, fit the existing `AccuracyMonitor`
`UserWarning`-subclass pattern — same table `docs/solvers.md` already
documents: pH jump, Newton iterations, charge-balance residual, scipy
rejection rate, ionic-strength regime):**

(a) If `clamp_fn=None`, check post-step whether any species went negative
(beyond a small tolerance) and warn.

(b) Independently: compare the delta dict before/after `clamp_fn` is applied;
if changed beyond tolerance, emit a diagnostic naming which species were
clamped and by how much (e.g. "removal rate for `S_lcfa` scaled from −12.3 to
−8.1 mol/h to prevent negative inventory over dt_h=0.5"). This is useful
beyond visibility — it signals that the chosen `dt_h`/kinetics combination is
aggressive enough that clamping is doing real work, an early hint that a
smaller `dt_h` might be needed rather than relying on the correction.

### 9. Demo notebook: authoring a custom `StepSolver` / `SystemSolver`

**Finding:** `demos/model_api/solver_comparison.py` only compares the five
*shipped* `SystemSolver`s against each other. No demo currently walks through
*writing your own* — a real gap for a codebase whose whole extension model is
"the protocol is thin enough that you write your own."

**Resolution (user wants this as a feature-demonstration notebook):**

- A custom `StepSolver` with a deliberately different intra-CV ordering than
  any shipped solver (e.g. algebraic instantaneous-equilibrium solve, then a
  hand-written integrator for kinetics only, then boundaries, then transfer).
- A custom `SystemSolver` with bespoke inter-CV interleaving beyond
  `StrangSplittingSystemSolver`'s symmetric half-step (e.g. an asymmetric
  N-stage split).

```python
class MyOrderedStepSolver:
    """Solve instantaneous equilibria, THEN feed, THEN kinetics with a
    custom integrator, THEN transfer — a bespoke ordering matching
    [some reference model / paper convention]."""

    def solve_step(self, cv, dt_h, t_h, external_source_terms=None):
        engine = cv.reaction_system.engine
        result = engine.solve(phases=cv.phases, T_K=cv.phases["liquid"].T_K)
        result.apply_to_phases(cv.phases)                      # 1. instantaneous
        # ... apply external_source_terms directly here ...    # 2. feed
        # ... your own kinetic integration, any method ...      # 3. kinetics
        # ... internal transfer last ...                       # 4. transfer
        return AdvanceResult(...)
```

### 10. Rename `EulerSnapshotSolver` → `SimultaneousEulerSolver`, `ScipyODESolver` → `SimultaneousAdaptiveSolver` — **SHIPPED 2026-07-06**

**Finding:** both names leaked implementation history rather than describing
current behavior. `ScipyODESolver` names the library it happens to be built
on (`scipy.integrate.solve_ivp`), not what it conceptually does — if the ODE
backend ever changed, the name would be actively wrong. `EulerSnapshotSolver`
bakes its (currently fixed) integration method into the name, which would
become misleading the moment it grows an `integrator=` parameter (item 5).
Comparing the two also surfaced something not previously obvious: they share
the *same* composition scheme (unsplit/simultaneous — every sub-system reads
one common state) and differ only in integration method (single discrete
Euler batch vs. continuous adaptive stepping) — a relationship the old names
(`Snapshot` vs. `Scipy`) didn't reveal at all.

**Shipped:** renamed to `SimultaneousEulerSolver` and
`SimultaneousAdaptiveSolver` — consistent `<CompositionScheme><IntegrationMethod>Solver`
convention, shared "Simultaneous" root now visible. Clean rename, no
backwards-compat aliases (per this repo's no-shims convention). Also fixed
the stale "this is the default solver" docstring on
`SimultaneousEulerSolver` while touching it (see the now-resolved bullet
below).

Scope of the rename: `src/core/solvers.py` (definitions), all call sites in
`src/core/` (`control_volume.py`, `simulation.py`, `gas_liquid_link.py`,
`interfaces.py`, `__init__.py`) and `src/monitoring/accuracy.py`, one model
builder (`models/vlmodels/fermenter/config/builder.py`), tests
(`test_simulation.py`, `test_cv_advance.py`, and
`test_scipy_ode_solver_jac.py` → renamed to
`test_simultaneous_adaptive_solver_jac.py`), demos (`export_results.py`,
two `.ipynb` notebooks, validated as well-formed JSON after edit), and live
docs (`docs/solvers.md`, `SOLVER_ARCHITECTURE.md`,
`CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md`, both class-diagram docs,
`architecture.md`, `README.md`, this note). `docs/shipped/` and
`docs/obsolete/` deliberately left untouched — frozen historical
record of the names as they were at each shipping point. Full test suite
run after: 1950 passed, 27 skipped, 0 failed — including a runtime-generated
`AccuracyWarning` message confirmed to correctly reflect the new class name.
Branch: `rename-simultaneous-step-solvers`, not yet merged (awaiting
confirmation).

---

## Minor documentation/consistency findings (low priority, from the original investigation)

These don't need design discussion, just cleanup whenever convenient:

- ~~`SimultaneousEulerSolver`'s class docstring in `solvers.py` says "This is
  the default solver"~~ — **fixed as part of item 10.** Was stale, leftover
  from when it was `GasLiquidVolume`'s default before that class was deleted
  in Phase 7. The actual default (`solver=None`) is the sequential body (see
  item 3, still unshipped). `docs/solvers.md` already had this right; the
  in-code docstring now does too.
- `docs/solvers.md` only documents Axis 1 (`StepSolver`) — no mention of
  `SystemSolver`/`system_solver.py`/the five Axis-2 solvers anywhere in it.
  `SOLVER_ARCHITECTURE.md` is the only doc tying both axes together; someone
  reading only `docs/solvers.md` wouldn't discover Axis 2 exists.
- Two independent controller-event-firing mechanisms exist:
  `EventScheduler` (priority-queue, used by the Strang/Multirate/
  ImplicitTransport solvers) and `MonolithicODESolver`'s own inline
  `ctrl_schedule` list — not a bug, but a duplicated concept that could drift.
- `StepSolver` is a plain `Protocol`; `SystemSolver` is `@runtime_checkable`.
  Inconsistent — the one place this matters
  (`Simulation._solver_for` telling Axis-1 and Axis-2 objects apart) works
  around it with `hasattr(solver, "advance_system")` duck-typing instead of
  `isinstance(solver, SystemSolver)`, which would work fine given the
  decorator already present.
- `ControlVolume.compute_jacobian`'s docstring says it's "used by
  `MonolithicODESolver` (Phase E)" — but `MonolithicODESolver.advance_system`
  never calls it; it doesn't pass a `jac` argument to `solve_ivp` at all. The
  actual Jacobian consumer today is Axis-1's `SimultaneousAdaptiveSolver`
  (`use_engine_jacobian=True`, via `engine.jacobian_dz_dy()`).

---

## Open questions / known limitations (not resolved into "will fix" — worth being aware of)

Surfaced while discussing whether to keep `cv.advance()`'s signature as-is;
neither has an agreed resolution yet:

- **Adaptive macro-stepping isn't expressible.** Every solver — including
  `SimultaneousAdaptiveSolver`, which does adaptive sub-stepping *internally* — is handed
  a fixed `dt_h` by the outer loop and must fully consume it. There's no
  channel for a solver to tell the orchestrator "I only trust this state for
  `dt_h/3`, call me again with a smaller step." If genuine error-controlled
  macro-stepping is ever needed, `AdvanceResult` would need an "actual dt
  consumed" field and both `Simulation._step_default()` and every
  `SystemSolver` would need to read it. No trigger for this yet.
- **Solvers are coupled to the concrete `ControlVolume` class, not an abstract
  interface.** `StepSolver.solve_step` is typed to take a `ControlVolume`, but
  implementations reach into `cv.phases`, `cv.internal_interfaces`,
  `cv.boundaries`, `cv.reaction_system`, `cv._accuracy_monitor` as concrete
  attribute access, not through a narrow Protocol. A different CV-shaped
  container (not subclassing `ControlVolume`) couldn't reuse
  `SimultaneousEulerSolver`/`SimultaneousAdaptiveSolver` without duck-typing the whole class.
  No trigger for this yet either — flagging because it's expensive to
  retrofit later if it does come up.

---

## Deferred: the interleaving / multi-CV-aware solver tier

**Flagged by the user as the next thread to continue exploring.** Summary of
where this stands:

**What's missing:** true *sub-CV-step* interleaving — e.g. applying half an
inter-CV flux mid-way through integrating one CV's own kinetic ODE, rather
than only before/after a whole `cv.advance()` call. `StepSolver.solve_step()`
only ever sees one CV plus a flat `external_source_terms` rate dict; it has
zero visibility into other CVs or the links between them, by design (Axis 1
is scoped to exactly one CV's topology). This is `SOLVER_PROMOTION.md`'s
deferred "Option C" — "a separate solver tier that operates on a multi-CV
topology object" — never built.

**What is already possible today, for context:**
- Ordering *within* one CV (instantaneous reactions, kinetics, boundaries,
  transfer, in any order/method) — fully open via a custom `StepSolver`.
- Ordering of inter-CV flux *relative to* a whole CV step (before/after/split
  around) — fully open via a custom `SystemSolver`, following the precedent
  `StrangSplittingSystemSolver` already sets (`links(dt/2) → cv.advance(dt) →
  links(dt/2)`); an "extreme"/bespoke version (asymmetric splits, more
  stages) just needs a bespoke `advance_system()` in the same style.

**Why this matters — two distinct motivations, worth keeping separate when
scoping:**
1. **Physical necessity**, when transport and kinetic timescales are
   comparable rather than well-separated: compartmental fermenter models
   where circulation timescale ≈ reaction timescale (mixing-time/
   reaction-time ratio driving observed heterogeneity — a real bioprocess
   scale-down research topic), and HPLC-style column cells where fast
   adsorption/desorption is comparable to inter-cell diffusion timescale.
   Both already named as transport-stiff examples in `SOLVER_ARCHITECTURE.md`.
2. **Reproducibility against an external reference implementation or paper**
   that itself uses a specific discrete splitting order — independent of
   whether the timescales in your own model are actually close.

**`MonolithicODESolver` is not a substitute for this, in either direction —**
worth restating precisely since it came up directly: Monolithic *eliminates*
splitting bias (co-integrates everything in one continuous ODE, converging
toward the unsplit solution as tolerance tightens). A specific interleaving
scheme *deliberately reproduces a chosen bias* to match a reference model that
itself splits (BSM2, ADM1, PHREEQC-style — literally the reference
implementations this codebase's own default is built to be comparable with).
**Decision rule:** if the reference model you're matching is itself a genuine
unsplit/DAE solve, `MonolithicODESolver` is the right tool. If the reference
model uses its own discrete operator-splitting convention, Monolithic will
actively *not* reproduce it — you need to control the interleaving directly.
Secondary practical differences: Monolithic can't easily express literal
discrete jumps at reference-model-specific points (only controller-firing
segment boundaries today), and co-integrating everything is exactly the
costlier, harder-to-converge regime the design doc cites as the reason
speciation stays algebraic rather than folded into one big ODE.

---

## Relationship to other docs

- [`../design/SOLVER_ARCHITECTURE.md`](../design/SOLVER_ARCHITECTURE.md) — the
  shipped two-axis anchor doc (Phases A–E) this note extends. Not superseding
  it; this is a refinement pass over the existing pure-scipy track, separate
  from the SUNDIALS/DAE Phase F/G placeholders it defines. (This link was
  bare/broken — `SOLVER_ARCHITECTURE.md` lives in `docs/design/`, not this
  file's directory, in either its pre- or post-move location — fixed while
  moving this file to `shipped/`.)
- [`SOLVER_PROMOTION.md`](SOLVER_PROMOTION.md) — origin of the gas/liquid
  coupling (item 6) and the deferred multi-CV-topology solver tier
  ("Option C", see the interleaving section above). Now same-directory
  (both live in `shipped/`).
- [`../solvers.md`](../solvers.md) — Axis-1-only summary doc when this note
  was written; now has an Axis-2 section too, and items 3/6/7 are documented
  (see `STEP_SOLVER_REFINEMENT_CHECKLIST.md` checkpoint 9).
- [`../design/MASS_EXCHANGE_ARCHITECTURE.md`](../design/MASS_EXCHANGE_ARCHITECTURE.md) —
  already identifies a related, independent extension (`SequentialIterativeSystemSolver`
  / SIA) as a candidate next Axis-2 solver; not covered by this note but worth
  cross-checking scope against if a solver-architecture phase is opened, since
  both would touch `system_solver.py`.

## Suggested next steps

1. Continue the interleaving/multi-CV-aware solver tier discussion (per the
   user's stated intent) before finalizing scope — it may reshape how items
   1, 2, and 5 are sequenced if a genuinely new solver tier gets scoped.
2. If/when this is picked up as an actual phase: items 1–4 and 6–9 are
   independent, low-risk, and could ship together as one phase following this
   repo's usual branch-per-phase convention. Item 5 (composition/integrator
   factoring) is a bigger, separate lift and shouldn't block the others.
3. Follow the standard phase workflow: write a checklist file modelled on
   existing `*_CHECKLIST.md` files, branch off `main`, merge with `--no-ff`,
   tag `<phase-name>-shipped` on completion.
