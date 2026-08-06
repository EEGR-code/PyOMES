# Speciation Engine Level Retirement — Design Note

> **Status:** Shipped 2026-06-03. Tag `speciation-level-retirement-shipped`.
> Branch merged into `main` with `--no-ff`.

## Context

`SpeciationEngine.__init__` carries a `level` parameter (1, 2,
2.5, or 3) that dispatches to one of three internal solver
implementations:

- **Level 1** → `ChemistryLevel1Model`
  ([`src/speciation/chemistry_level1.py`](../../src/speciation/chemistry_level1.py))
  — water + carbonate + ammonium + VFAs + strong ions.
- **Level 2** → `ChemistryLevel2Model`
  ([`src/speciation/chemistry_level2.py`](../../src/speciation/chemistry_level2.py))
  — Level 1 + phosphate + bisulfate, broader strong-ion table,
  slightly different activity-coefficient handling.
- **Level 3** → `LegacySpeciationModel`
  ([`src/speciation/legacy_adapter.py`](../../src/speciation/legacy_adapter.py))
  — fallback to the original direct-solve adapter.

After the chemistry-unification track (2026-05-13 through
2026-05-16), the set of supported equilibria became a function of
the declared `EquilibriumReaction` instances on the
`ReactionSystem`. The `level` knob survives as a runtime
implementation choice, but **the chemistry it expresses is now
redundant** — every "level 1 chemistry" is just a particular
set of equilibrium reactions, and the same holds for level 2.

## The reframing

The three solver implementations are not three *engines*. They
are three **chemistry presets** masquerading as engines:

- "Level 1" is a preset for water + carbonate + ammonium + VFAs.
- "Level 2" is that preset *plus* phosphate, bisulfate, and a
  wider strong-ion table.
- "Level 3" is "anything goes via the legacy adapter".

The numerical math that solves a given set of equilibria —
Newton iteration on charge balance, Davies / SIT activity
coefficients, Van 't Hoff temperature correction — is shared
across the three. What differs is *which chemistry is hardcoded
into the dispatcher*. That hardcoding is the whole point of the
level concept, and the whole point of removing it.

Once chemistry is fully expressed through
`EquilibriumReaction` declarations (state-unification C3
completed this for the engine input side), there is a single
solver. No dispatch. The "level" parameter goes away entirely.

The level-N hardcoded chemistry tables become **`ChemistryDatabase`
presets** — packaged bundles of `EquilibriumReaction` objects
the user can import to seed a `ReactionSystem` without
re-declaring every acid-base pair. For example:

```python
from VLsim.chemistry.presets import bsm2_chemistry, full_aqueous_chemistry

# What "level 1" really was — the BSM2-style preset.
system = ReactionSystem([
    *bsm2_chemistry,          # 5 EquilibriumReactions: water + CO2 + NH4 + 4 VFAs
    KineticReaction(...),      # user-supplied kinetics
])

# What "level 2" really was — extended preset.
system = ReactionSystem([
    *full_aqueous_chemistry,  # adds phosphate, bisulfate, broader strong ions
    KineticReaction(...),
])
```

This is the same `ChemistryDatabase` rollout described in the
trigger-gated chemistry-unification-3b phase
([`CHEMISTRY_UNIFICATION_PLAN.md`](CHEMISTRY_UNIFICATION_PLAN.md)).
Retiring the `level` parameter is downstream of that rollout, but
the work overlaps substantially.

## Why this is debt rather than a bug

It works. Tests pass. The three implementations dispatch
correctly. The state-unification phase's canonical-naming refactor
(C3) operates above the level dispatch.

The cost shows up when:

1. **A new chemistry feature has to land in all three
   implementations.** Most recent additions have touched level 1
   + level 2 and left level 3 lagging.
2. **Numerical-accuracy regressions are level-specific.** The
   three implementations have slightly different convergence
   characteristics — pH at ionic strength I = 0.1 can differ by
   1e-3 units between levels.
3. **Users have to choose `level` without knowing what it means.**
   The `state-unification` C4
   `cv.reaction_system.configure_engine` API hides the knob
   behind a default (level=2), but the duplication is still in
   the codebase.

## Concrete reproducer (added 2026-05-29)

The Q1 demo
[`demos/model_api/model_construction/raw_construction.py`](../../demos/model_api/model_construction/raw_construction.py)
provides a ready-made reproducer for the central claim above —
that `level` materially changes solver behaviour even when the
declared `EquilibriumReaction` instances are the same.

The demo declares one acetate dissociation
(`AceticAcid ⇌ Acetate⁻ + H⁺`, `log_K=-4.756`) plus one
cross-phase CO₂ partition. With the seeded initial state
(~20 mM total acetate, no strong ions, `n_mol["H+"]` seeded at
`1e-7 × V_L`):

| `configure_engine(level=...)` | Initial pH | Final pH | `n_mol["H+"]` updated each step? |
|---|---|---|---|
| `1` | 3.234 (consistent with the equilibrium) | 7.000 (substrate exhausted) | yes |
| `2` (default) | 7.000 (matches seed) | 7.000 | **no** — pinned at seed |

Same `EquilibriumReaction` declaration, same CV topology, same
controller, same boundaries. Only `level` differs. At `level=2`
the PI pH controller dosed 0.4 mol H₃PO₄ over the run with no
effect on pH, because H⁺ was never updated by the engine.

The demo currently pins `level=1` with a comment pointing at
this planning note. When the retirement ships, **removing the
`configure_engine(level=1)` line should leave demo behaviour
unchanged** — Q1 is a built-in regression target.

This is stronger than the "duplication" framing in the
[Why this is debt](#why-this-is-debt-rather-than-a-bug) section.
At `level=2`, declared equilibria appear to be either suppressed
by, or silently no-op'd because of, the auto-added BSM2-style
default ladders that the level-2 dispatcher pre-loads (the
demo's `n_mol` doesn't carry the expected total trackers like
`CT_TIC` / `CT_NH_T`). Worth investigating the root cause as
part of step 2 of the [Approach](#approach-when-triggered) —
"pick one solver implementation" — so the unified solver makes
declared reactions authoritative regardless of which auto-
chemistry would otherwise have been loaded.

## Trigger conditions to revisit

Pick this up when:

- **`ChemistryDatabase` (chemistry-unification-3b) starts.** The
  natural place to package the level-N chemistry tables as
  importable bundles. Retiring the `level` parameter becomes a
  side effect of completing that rollout.
- **A chemistry feature needs to land in all three solvers.** If
  the next chemistry-feature pass touches more than one
  implementation, fold the retirement into the same branch.
- **The solver-implementation divergence causes a real bug.** If
  a user-reported issue traces to "this works at level 1 but
  not level 2", the cost of maintaining the divergence has
  crossed the threshold.

## Approach when triggered

The work is **retirement**, not unification:

1. **Extract the chemistry tables.** Each level's hardcoded
   chemistry (which equilibria, which pKa values, which strong
   ions) becomes a `tuple[EquilibriumReaction, ...]` bundle in
   `src/chemistry/presets/`. The numerical math that consumes
   those equilibria stays in *one* solver module.
2. **Pick one solver implementation.** The simplest path is to
   take Level 1's lean implementation and extend it to cover the
   Level 2 chemistries — most of which can be expressed by
   declaring more `EquilibriumReaction` instances. The Level 2
   activity-coefficient nuances that genuinely differ get folded
   into the unified solver as configurable knobs (`use_activity`,
   `activity_model` — already exposed through
   `configure_engine`).
3. **Migrate Level 3 (legacy adapter) callers.** Audit which
   chemistries the legacy adapter is currently the only one
   handling. Add `EquilibriumReaction` declarations for them.
   Delete `legacy_adapter.py`.
4. **Delete `ChemistryLevel1Model`, `ChemistryLevel2Model`,
   `LegacySpeciationModel`.** Delete the level dispatch in
   `SpeciationEngine.__init__`. The engine just *solves*; it
   doesn't choose its chemistry.
5. **Drop the `level` constructor parameter.** Remove from the
   public API. `configure_engine(use_activity, activity_model)`
   stays.
6. **Update the AccuracyMonitor.** The level-specific
   tolerances become single-solver tolerances. Re-baseline any
   tests that drift.

After the retirement, the engine has **one** solver, the chemistry
is **always** defined by declared `EquilibriumReaction` instances,
and the user's choice of "which chemistry" is the choice of
"which preset to import from `ChemistryDatabase`."

## Out of scope for `state-unification`

`state-unification` C4 introduces
`cv.reaction_system.configure_engine(use_activity, activity_model)`
as the user-facing engine-config API. The `level` parameter stays
behind the default (level=2 since that covers more chemistry) and
is not exposed in the new API. That hides the knob but does not
retire it — the three solver implementations remain on disk and
remain dispatched-to by the engine's internal level dispatch.

This note records the deferred retirement. No implementation is
scheduled until `ChemistryDatabase` (chemistry-unification-3b) is
activated.

## Cross-references

- [`CHEMISTRY_UNIFICATION_PLAN.md`](CHEMISTRY_UNIFICATION_PLAN.md) —
  Phase 3b carries the `ChemistryDatabase` rollout that this
  retirement folds into.
- [`STATE_UNIFICATION.md`](STATE_UNIFICATION.md) — phase that
  introduced `configure_engine` and hid `level` behind defaults.
- [`../design/SPECIATION_REACTIONMODEL_BOUNDARY.md`](../design/SPECIATION_REACTIONMODEL_BOUNDARY.md) —
  design reference on the runtime split between kinetic and
  equilibrium reactions. Orthogonal to the level concept; the
  split survives the retirement.
