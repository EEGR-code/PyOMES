# Speciation Engine ↔ ReactionModel Boundary — Design Reference

> **Status:** Architecture decision record, 2026-05-18 (updated 2026-06-11).
> **No work scheduled.** The current architecture is intentional. This note
> records why the boundary exists and why it is shaped the way it is.

## Context

After chemistry-unification-1..4 (shipped through 2026-05-16) one
architectural question kept surfacing: why does `BisectionChemicalEquilibriumEngine` remain a
distinct type, when acid-base equilibria are now `EquilibriumReaction` objects
sharing declaration machinery with kinetic reactions?

## What the chemistry-unification track actually unified

| Layer | Pre-unification | Post-unification |
|---|---|---|
| Declaration type | `Reaction` (kinetic only) + bespoke pKa registries | `KineticReaction` + `EquilibriumReaction` — two typed classes (STATE_UNIFICATION, 2026-05-22) |
| Validation | Atoms balance, kinetic only | Atoms + (opt) charge, both kinds |
| Stoichiometry | `StoichiometryEntry(species_id, ...)` | `StoichiometryEntry(species: Species, ...)` — shared |
| Container | `ReactionSet` (kinetic) + `EquilibriumSet` (separate) | `ReactionSystem` (pre-bucketed; no partition step) |
| Runtime home | `cv.reaction_model` only | `cv.reaction_system` (kinetic) + `cv.reaction_system.engine` (equilibrium) |
| Solver | scipy ODE | single Euler / scipy ODE (kinetic) + brentq + Davies fixed-point (equilibrium) |

**Unified**: declaration type, validation infrastructure, stoichiometry,
container, `Species`.

**Deliberately not unified**: runtime homes, solvers, output shapes.

This was scoped from the start.
[`CHEMISTRY_UNIFICATION.md:7-14`](../implementation/shipped/CHEMISTRY_UNIFICATION.md#L7-L14):

> "This document proposes unifying the declaration surface ... while keeping
> the existing algebraic solver intact. This is not a solver rewrite."

## Why the solvers stay split

Speciation operates 10+ orders of magnitude faster than kinetic bioreactions.
The runtime split is justified by stiffness (DAE-vs-ODE timescale separation),
not historical artefact. See
[`SOLVER_ARCHITECTURE.md`](SOLVER_ARCHITECTURE.md)
§ "Why operator splitting is the right default for speciation" for the full
argument and the stiffness taxonomy.

## Structural shapes side by side

[`ReactionModel.compute_rates(env)`](../../src/reactions/protocols.py#L43):
`(state) → {phase: {species: mol/h}}`. Produces *derivatives*. The integrator
applies them to `phase.n_mol`.

[`BisectionChemicalEquilibriumEngine.solve(**kwargs)`](../../src/chemical_equilibrium/engine.py#L193):
`(state) → {pH, IonicStrength, species, alphas, ...}`. Produces
*constraint-satisfying derived properties*. Equilibrium-output species (H⁺,
OH⁻, HCO₃⁻, …) are **written to `phase.n_mol`** via the engine's
`_refresh_derived` path (STATE_UNIFICATION, 2026-05-22; supersedes B2
resolution in [`CHEMISTRY_UNIFICATION_PLAN.md:271-293`](../implementation/shipped/CHEMISTRY_UNIFICATION_PLAN.md#L271-L293)).

Different output shapes, different side-effects, different solvers — these are
not the same kind of object.

## Why `PropertyCalculator`, not `ReactionModel`

Two-part justification:

- **Consumer order.** Kinetic rates and gas-liquid links read pH / alphas as
  *inputs*. The speciation engine runs first (step 1 of `advance()`), then
  property calculators (step 1b), before any reaction integration
  ([`control_volume.py:449-463`](../../src/core/control_volume.py#L449-L463)).
- **Output contract.** `PropertyCalculator.compute(phase, T_K, P_atm) → float`
  produces a single derived property value that lands on
  `phase.properties[calc.key]`. `ReactionModel.compute_rates` returns rates
  and does not.

`PropertyCalculator` covers viscosity, density, and other genuinely
property-flavoured computations — it survives any future change to the
equilibrium-solver story.

## Runtime routing: where equilibrium reactions actually live

`EquilibriumReaction` deliberately does not implement the `ReactionModel`
protocol — it has no `compute_rates`. So:

- An `EquilibriumReaction` is **declared** in the unified framework
  but **never executed** through the `ReactionModel` protocol.
- [`ReactionSystem`](../../src/reactions/reaction_system.py) pre-buckets at
  construction time, routing `EquilibriumReaction` instances to
  [`BisectionChemicalEquilibriumEngine.from_reactions(equilibrium_subset)`](../../src/chemical_equilibrium/engine.py#L95).
  `from_reactions` reads `log_K`, `dH_J_per_mol`, classifies the reaction
  into water / acid / cation_acid, and constructs `EquilibriumDef` entries
  on an `EquilibriumSet`. The original `EquilibriumReaction` objects are
  discarded post-construction.
- `cv.reaction_system` drives **kinetic-only** reactions through
  `compute_rates`. Equilibrium reactions live inside `engine._equilibrium_set`,
  accessed via `cv.reaction_system.engine`.

## Naming asymmetry: `ReactionModel` vs `BisectionChemicalEquilibriumEngine`

Three layered reasons the names look different:

1. **Type shape.** `ReactionModel` is a **protocol** with multiple implementers
   (`KineticReaction`, `ReactionSystem`, `BlackBoxReactionModel`). `BisectionChemicalEquilibriumEngine`
   is a **concrete class** with two solve paths: the declared-reaction path (via
   `EquilibriumSet`; primary) and a legacy kwargs path. "Model" reads as a
   protocol noun; "Engine" reads as a single named workhorse.
2. **Output contracts differ.** Rates vs properties. They cannot satisfy each
   other's protocols. Renaming `BisectionChemicalEquilibriumEngine` to `EquilibriumReactionModel`
   would suggest a shared protocol that doesn't exist.
3. **History + design intent.** `BisectionChemicalEquilibriumEngine` predates the unification
   track. The track did not rename it — consistent with the design doc's
   explicit "keep the solver split" stance. The asymmetric names quietly
   encode the asymmetric runtime contracts.

A lighter-touch fix that captures most of the discoverability win: sharpen
`ReactionModel`'s docstring to "kinetic reactions only; equilibria are routed
via `BisectionChemicalEquilibriumEngine.from_reactions` and never flow through this protocol."
No API churn.

## Resolved positions

- The chemistry-unification track did exactly what it announced: declarations
  unified, runtime split preserved. There is no missing Phase 5.
- The runtime split is justified by stiffness (DAE-vs-ODE timescale
  separation), not historical artefact or engineering convenience.
- `PropertyCalculator` (which replaced `PropertySolver` in STATE_UNIFICATION,
  2026-05-22) covers viscosity, density, and other property computations — it
  survives any future change to the equilibrium-solver story.
- Type split of `Reaction` into `KineticReaction` / `EquilibriumReaction` was
  done in STATE_UNIFICATION (2026-05-22); runtime errors became type errors.
- Naming asymmetry of `ReactionModel` vs `BisectionChemicalEquilibriumEngine` is intentional
  encoding of runtime contract asymmetry; renaming would mislead.

## Cross-references

- [`SOLVER_ARCHITECTURE.md`](SOLVER_ARCHITECTURE.md) —
  two-axis integration architecture; stiffness taxonomy; why operator splitting
  is the right default; Phase F (per-CV DAE opt-in) trigger conditions
- [`../shipped/CHEMISTRY_UNIFICATION.md`](../implementation/shipped/CHEMISTRY_UNIFICATION.md) — rationale doc
- [`../shipped/CHEMISTRY_UNIFICATION_PLAN.md`](../implementation/shipped/CHEMISTRY_UNIFICATION_PLAN.md) —
  branch slicing, resolved knobs (B1/B2/B3)
- [`../shipped/CHEMISTRY_UNIFICATION_1_CHECKLIST.md`](../implementation/shipped/CHEMISTRY_UNIFICATION_1_CHECKLIST.md) —
  declaration-unification implementation log
- [`../shipped/CHEMISTRY_UNIFICATION_3_CHECKLIST.md`](../implementation/shipped/CHEMISTRY_UNIFICATION_3_CHECKLIST.md) —
  cross-phase equilibrium reaction infrastructure
