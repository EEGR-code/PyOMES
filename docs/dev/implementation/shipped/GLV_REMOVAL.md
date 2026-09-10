# GLV Removal — Note for a Future Phase

> **Status: Shipped 2026-05-12** — Phase 7 deleted `GasLiquidVolume`
> entirely.  Implementation tracked in
> [PHASE7_CHECKLIST.md](PHASE7_CHECKLIST.md); this note is kept as
> the historical design record.

## What this is

A planning note for the eventual removal of `GasLiquidVolume` as a
class. After Phase 6 (see [PHASE6_CHECKLIST.md](PHASE6_CHECKLIST.md))
ships, GLV is reduced to a near-transparent shim — solvers and
boundaries live on `ControlVolume`, and `GasLiquidAdvanceResult` is
gone. The remaining responsibilities of GLV are small enough that a
future phase can delete the class entirely with mostly-mechanical
caller migration.

This is **not** active work. It depends on Phase 6 having shipped
first; that's the only hard precondition.

## Post-Phase-6 GLV — what's left

After Phase 6, `GasLiquidVolume` holds:

1. A reference to its inner `ControlVolume` (`self._cv`).
2. A reference to the configured `StepSolver` (`self._solver`),
   which is forwarded into `cv.advance(solver=...)`.
3. A `_last_result` cache populated each `advance()` call.
4. The `pH` and `ionic_strength` accessors, both of which read
   `self._last_result.properties.get("speciation")`.
5. A `transfer_link` accessor (returning `self._cv.internal_interfaces[0]`
   or similar — convenience for controllers that adjust kLa).
6. The `gas_phase` / `liquid_phase` accessors (returning
   `self._cv.phases["gas"]` / `["liquid"]`).
7. The `phases` / `phase_keys` / `__getitem__` / `__contains__`
   forwarders.
8. The `apply_external_flux` forwarder.
9. The `total_mol` forwarder.
10. The `snapshot()` method.
11. The `.create(henry=..., kLa=..., ...)` classmethod factory that
    builds the inner CV with a configured `KineticGasLiquidLink`.

Everything except (1)–(4) and (11) is a one-line forwarder to
`self._cv`. (4) is the only behaviour that doesn't trivially live
on a CV today.

## What removal entails

Three pieces of structural work, plus caller migration:

### 1. Re-home the `pH` / `ionic_strength` accessors

This is the only non-trivial design question. Options:

**A. Move to `ControlVolume`.** `cv.pH` and `cv.ionic_strength` read
from a CV-level `_last_result` cache populated by `cv.advance()`.

- Pros: matches Phase 6's "solvers on CV" framing — the result
  cache lives where the integration happens.
- Cons: `pH` and `ionic_strength` are speciation-flavoured; making
  them universal accessors on every CV (including pure gas CVs,
  HPLC cells, etc.) is conceptually noisy. Likely best to make
  these accessors return `None` when no `"speciation"` result is
  present — same as today's GLV behaviour.

**B. Move to a `RunResult` object the orchestrator carries.**
The orchestrator collects per-step `AdvanceResult` objects into a
time series and exposes `.pH(t)` / `.ionic_strength(t)` on the
result.

- Pros: cleanly separates per-step state from the cumulative
  history. Doesn't add fields to `ControlVolume`.
- Cons: callers that today read `glv.pH` mid-simulation (between
  `advance()` calls — controllers, especially) need a different
  access pattern. Significant API change.
- **Deferred** as a follow-on phase rather than addressed here;
  see [RUN_HISTORY.md](RUN_HISTORY.md) for the standalone design
  note.

**C. Drop the convenience accessors entirely.** Callers read
`result.properties.get("speciation").pH` directly off the
`AdvanceResult`.

- Pros: simplest. No new abstraction.
- Cons: every existing caller using `glv.pH` updates. Significant
  caller churn.

**Recommendation pending a real design decision.** A is least
disruptive and most natural after Phase 6's promotion. B is
architecturally cleanest but a bigger change. Pick when this phase
becomes active.

### 2. Replace `.create(henry=..., kLa=...)` with a free-function factory

Today: `GasLiquidVolume.create(...)` is a classmethod that builds a
configured GLV (which itself wraps a CV with a `KineticGasLiquidLink`).

Replacement:
```python
def make_gas_liquid_cv(
    gas_phase, liquid_phase,
    henry, kLa=None, equilibrium_species=None,
    speciation_keys=None, speciation_corrections=None,
    henry_params=None,
    property_solvers=None, reaction_model=None,
    boundaries=None,
    label="",
) -> ControlVolume:
    ...
```

A free function (likely in
`src/core/gas_liquid_link.py` or a new
`src/core/gas_liquid_factories.py`) that builds the
`KineticGasLiquidLink` with the right phase keys, then returns a
`ControlVolume` wired up with the link, the supplied phases, and
the supplied other components. No wrapper class involved.

### 3. Migrate every caller

Sites using `GasLiquidVolume`:

- `models/vlmodels/fermenter/config/factory.py` — the main
  `FermenterFactory` builds a GLV via `GasLiquidVolume.create(...)`.
  Replace with the new free function returning a CV.
- `models/vlmodels/fermenter/config/builder.py` — the
  `FermenterBuilder.build()` returns a `GasLiquidVolume`. Change
  return type to `ControlVolume`.
- `models/vlmodels/fermenter/profiles.py` —
  `TemperatureSetpointTarget` and similar already use
  `getattr(glv, "property_solvers", [])` (post-Phase-5); these
  patterns work on any CV.
- `models/vlmodels/adm1/base.py`, `bsm2.py`, `bsm2_direct.py` — all
  construct GLVs via the public ctor or `.create`. Replace with
  the new free function.
- `systems/*.py` — every example/system script constructs a GLV
  via `FermenterBuilder().build()`. After (2), they receive a CV
  instead. Variable names like `glv` should be renamed (e.g. to
  `cv`), but that's stylistic.
- `tests/standalone/test_*.py` — many tests construct GLVs.
  Migrate construction; result-reading patterns should already
  work post-Phase-6.
- `src/chemistry/thermo_params.py` — already accesses GLV via
  `getattr(glv, "property_solvers", [])` etc.; works on any CV.
- `src/control/*.py` — controllers that adjust kLa via
  `glv.transfer_link.set_kLa(...)`; the `transfer_link` access
  pattern needs an equivalent on CV (e.g.
  `cv.internal_interfaces[0]` with a type check).

Estimated migration: ~30–50 files touched; mostly mechanical
once the design questions in (1) are settled.

## Triggers — when to actually do this

Don't do this speculatively. Phase 6's reduction of GLV makes the
removal cheap; doing it before Phase 6 ships would be a much bigger
disruption.

The trigger is one of:

- **Cleanup motivation.** GLV is a 30-line shim that does almost
  nothing; deleting it removes a layer of indirection from every
  fermenter / ADM1 code path.
- **A model that wants the unified `Vessel` / `Unit` shape** flagged
  in [CONTAINER_LAYERING.md](CONTAINER_LAYERING.md). GLV removal is
  on the path to that work.
- **Pressure to simplify the package's public API.** If users find
  the GLV / CV split confusing, removal collapses two concepts to
  one.

## Scope estimate

Larger than Phase 6 but smaller than the full
[CONTAINER_LAYERING.md](CONTAINER_LAYERING.md) unbundling:

- Two structural changes (re-home accessors, replace factory).
- Substantial caller migration (~30–50 files).
- Test updates (mostly mechanical: variable rename `glv` → `cv` and
  `glv.method()` patterns that were already CV-shaped after
  Phase 6).
- Documentation refresh: `docs/architecture.md`,
  `docs/class_diagrams.md` lose the `GasLiquidVolume` box.

Could be done as one PR with a checkpoint structure modelled on
[PHASE5_CHECKLIST.md](../shipped/PHASE5_CHECKLIST.md), or
split into two: (a) CV-level accessor design + free-function
factory, (b) caller migration sweep.

## Open design questions to resolve before starting

1. **Where do `pH` / `ionic_strength` live?** See section 1 above.
   This decision shapes much of the rest of the work.
2. **Naming of the free-function factory.** `make_gas_liquid_cv`?
   `build_gas_liquid_cv`? `GasLiquidCV` (a callable with PascalCase
   to mimic constructor syntax)? Convention to choose.
3. **What about `_last_result`?** Whether to keep it on the CV (as a
   private cache) or eliminate (callers must keep their own
   reference). Linked to question 1.
4. **Should the post-removal world rename `KineticGasLiquidLink` to
   something like `GasLiquidPhaseInterface`?** No — the name
   describes the physics (kinetic, with kLa), and renaming
   wouldn't change anything functional. Keep as-is.
5. **Should `GasLiquidAdvanceResult` come back?** It was dropped in
   Phase 6. After GLV removal, `cv.advance()` already returns
   `AdvanceResult`. No reason to bring it back.

## Relationship to other phases

- [PHASE7_CHECKLIST.md](PHASE7_CHECKLIST.md) — the implementation
  plan derived from this design note; turn here when work begins.
- [../shipped/PHASE6_CHECKLIST.md](../shipped/PHASE6_CHECKLIST.md)
  — shipped 2026-05-01.  Its reduction of GLV to a shim is what
  makes this work cheap.
- [../shipped/SOLVER_PROMOTION.md](../shipped/SOLVER_PROMOTION.md)
  — the design note that motivated Phase 6; Phase 7 is the
  natural follow-on.
- [RUN_HISTORY.md](RUN_HISTORY.md) — the deferred unified-history
  idea (option B in §1) that Phase 7 declined in favour of option
  C.  Picked up as a future phase if the trigger appears.
- [CHEMISTRY_UNIFICATION.md](CHEMISTRY_UNIFICATION.md) — installs
  `phase.pH` / `cv.current_pH()` with staleness guardrails.
  Sequenced after Phase 7; will replace the
  `result.properties["speciation"].pH` access pattern that Phase
  7 leaves callers in.
- [CONTAINER_LAYERING.md](CONTAINER_LAYERING.md) — the larger
  unbundling.  After GLV is removed, the next step in that
  direction is generalising `EulerSnapshotSolver` to arbitrary
  phase combos and providing a multi-CV-aware solver tier.
- [../shipped/CV_UPDATE.md](../shipped/CV_UPDATE.md) —
  the original five-phase plan that started the cleanup.
