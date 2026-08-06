# ControlVolume Advance Order Refactor

> **Status: Shipped 2026-04-29.** All five phases landed; the standalone
> test suite is at 740 passing (`python -m pytest tests/standalone -q`).
> This document is preserved as the historical record of the plan.
> See [PHASE5_CHECKLIST.md](PHASE5_CHECKLIST.md) for the implementation
> log and [CONTAINER_LAYERING.md](../phases-upcoming/CONTAINER_LAYERING.md) for the open
> follow-on layering question.

## Summary

The current `ControlVolume.advance()` runs property solvers (speciation, pH) **last**,
after reactions and internal transfer have already run. This forces the codebase to cache
results in `_last_properties` so the next step can read them — creating a one-step lag,
and requiring explicit workarounds: seeding initial pH in ADM1 models, injecting
speciation into snapshot CVs, and the gas-liquid link reaching directly into CV internals.

The fix is to run property solvers **first** in `advance()`, passing results directly
through the call chain. This eliminates `_last_properties` entirely and — critically —
unblocks making `KineticGasLiquidLink` a proper `PhaseInterface`, which means gas-liquid
transfer can live inside a single `ControlVolume` rather than requiring the
`GasLiquidVolume` wrapper.

### Why `_last_properties` is not needed

`_last_properties` exists only because speciation runs last. If it runs first:

- The result is live in the same call stack — passed directly to reactions and interfaces
- The next call to `advance()` recomputes fresh speciation from whatever the moles are then
- No cache, no lag, no seeding workarounds
- `AdvanceResult.properties` reflects the state that **drove** the step (pre-mutation),
  which is physically correct for an explicit Euler scheme

---

## Naming Conventions

This refactor introduces a second downstream channel in `advance()` — the property
solver's output. To keep the two channels legible and avoid overloading, both
channels are given explicit names:

| Channel | Name | Direction | Type | Contents |
|---|---|---|---|---|
| Chemistry inputs | `chem_env` | orchestrator → solvers | `Dict[str, Any]` | `acid_totals`, `acid_pKas`, `CT_P`, `CT_NH_T`, `strong_kwargs`, `t_h`, `logH_guess` |
| Solver outputs | `property_results` | solvers → reactions / interfaces | `Dict[str, PropertyResult]` | `{"speciation": PropertyResult(pH=..., ionic_strength=..., species={...})}` |

### Why `chem_env` (rename from `context`)

The existing parameter is called `context` — a widely recognised code smell for
"bag of miscellaneous stuff." What it actually carries is the chemistry framework
parameters the CV doesn't own: totals, pKas, strong ions, time. `chem_env` names
the concept — *the CV's chemical environment* — without cataloguing the contents.

- Shorter than `chem_params`, more descriptive than `chemistry` (which could read
  as a reference to the `src/chemistry/` module).
- Pairs conceptually with the existing `ReactionEnvironment` dataclass: `chem_env`
  is an input to `ReactionEnvironment` construction, not the same thing. The
  shared word signals a real relationship worth a reader noticing — both describe
  the chemical environment, at different stages of the pipeline.

### Why `property_results` (new parameter)

A dict of `PropertyResult` objects — naming it `property_results` makes the type
obvious at every call site. `properties` alone was considered but rejected as
too generic in general software (could mean object attributes, config, etc.).
`property_results` aligns one-to-one with the existing `PropertyResult` /
`PropertySolver` types.

### Resulting signatures

```python
def advance(self, dt_h, chem_env=None):
    # chem_env flows into property solvers
    # property_results flows out to reactions and interfaces

def step_internal_transfer(self, dt_h, property_results=None):
    ...

def compute_flux(self, state_a, state_b, dt_h, property_results=None):
    ...
```

### Scope of the rename

The existing `context` parameter on `advance()`, `PropertySolver.solve()`, and
everywhere else it appears in the codebase is renamed to `chem_env` as part of
this refactor. This is mechanical find-replace across:

- `src/core/control_volume.py`, `src/core/property_solver.py`
- `src/core/speciation_solver.py`, `src/speciation/*`
- `models/vlmodels/**` orchestrator call sites
- All tests that construct `context` dicts

---

## New `advance()` Sequence

```
advance(dt_h, chem_env):
  1. Run property solvers on current phase state  →  property_results
     (solvers receive chem_env as their chemistry framework)
  2. Apply external source terms                  (feed arrives before reactions)
  3. Integrate reactions                          (receives property_results directly — no cache read;
                                                   reactions see post-feed moles)
  4. step_internal_transfer(dt_h, property_results=property_results)
       └─ iface.compute_flux(phase_a, phase_b, dt_h, property_results=property_results)
  5. Return AdvanceResult(properties=property_results, ...)
```

This sequence follows `ORDERING.md`: feed runs before reactions because, in
a continuously fed CV, material that arrives in the timestep is immediately
available for reaction.  Note (per `ORDERING_CRITIQUE.md`) that this
sequential operator splitting is **not** equivalent to the simultaneous
explicit-Euler step that `EulerSnapshotSolver` performs — both are valid
O(dt) Euler variants, and they differ by O(dt²) cross-terms.  A future
iteration may make the splitting strategy pluggable on the CV (e.g. a
`splitting="sequential" | "snapshot" | "strang"` parameter); deferred
because no current model needs the alternatives.

---

## `compute_reaction_rates` caveat

`ControlVolume.compute_reaction_rates` is the read-only mirror of the
reactions step inside `advance`.  Its new signature is:

```python
def compute_reaction_rates(self, chem_env=None, property_results=None) -> ...
```

The method does **not** run property solvers internally — it builds a
`ReactionEnvironment` directly from the supplied `property_results`.
When `property_results=None`:

- `pH` on the resulting environment is `None`
- No species concentrations from speciation are forwarded

Reaction models that depend on `pH` will then either error out or fall
back to a hard-coded default — both are footguns.  Callers using this
method as a snapshot-style sibling of `advance` should run the property
solvers first and pass the dict through, e.g.

```python
prop = {s.key: s.solve(cv.phases, chem_env, dt_h) for s in cv.property_solvers}
rates = cv.compute_reaction_rates(chem_env=chem_env, property_results=prop)
```

The contract is intentional: the method is read-only by design and
leaves the choice of *when* to run speciation to the caller.

---

## Implementation Plan

### Phase 1 — Fix the advance order and eliminate `_last_properties`

#### `src/core/interfaces.py`
- Add optional `property_results: dict = None` to `PhaseInterface.compute_flux()`
  signature

#### `src/core/property_solver.py`
- Rename `context` parameter on `PropertySolver.solve()` to `chem_env`
- Update docstring accordingly

#### `src/core/control_volume.py`
- Remove `self._last_properties` from `__init__`
- Rename the `context` parameter on `advance()` to `chem_env`
- Reorder `advance()`:
  1. Run property solvers (passing `chem_env`) → `property_results`
  2. Pass `property_results` into `_build_reaction_environment()` directly
     (remove the `self._last_properties.get(...)` read at line 355)
  3. Integrate reactions
  4. Apply external source terms
  5. Call `step_internal_transfer(dt_h, property_results=property_results)`
  6. Return `AdvanceResult(properties=property_results, ...)`
- Update `step_internal_transfer()`:
  - Add `property_results=None` parameter
  - Pass to each `iface.compute_flux(..., property_results=property_results)`
- Update `_build_reaction_environment()`:
  - Accept `property_results` and `chem_env` as direct parameters instead of
    reading from `self._last_properties` and an ambiguous `context` dict

---

### Phase 2 — Remove all `_last_properties` workarounds

#### `models/vlmodels/adm1/base.py` and `models/vlmodels/adm1/bsm2.py`
- Remove the `_SeedSpec` / `_last_properties` seeding blocks
- No longer needed: the first `advance()` call computes speciation from initial moles
  naturally, the same as every subsequent step

#### `models/vlmodels/fermenter/config/factory.py`
- Remove construction-time `_last_properties` seeding (lines 551–552)

#### `src/core/gas_liquid_volume.py`
- Update `.pH` and `.ionic_strength` properties — instead of reading from
  `_cv_liquid._last_properties`, read from a stored `self._last_result` set after
  each `advance()` call

#### `src/core/solvers.py` (`EulerSnapshotSolver`, `ScipyODESolver`)
- Remove explicit writes to `_cv_liquid._last_properties`
- Store the final `GasLiquidAdvanceResult` on the `GasLiquidVolume` instance as
  `self._last_result` instead

#### `src/chemistry/thermo_params.py`
- Update the `getattr(cv_liq, "_last_properties", {})` read (line 837) to use the
  new result store

---

### Phase 3 — Make `KineticGasLiquidLink` implement `PhaseInterface`

#### `src/core/gas_liquid_link.py`
- Add `phase_a_key = "gas"` and `phase_b_key = "liquid"` properties
- Add `compute_flux(self, state_a, state_b, dt_h, property_results=None)` method:
  - Wraps the existing `compute_flow` logic
  - Reads speciation from `property_results` instead of reaching into
    `_last_properties`
- The class now satisfies both `CVLink` (existing) and `PhaseInterface` (new)
- No existing code breaks — dual-protocol, purely additive

This enables:
```python
cv = ControlVolume(
    phases={"gas": GasPhase(...), "liquid": LiquidPhase(...)},
    internal_interfaces=[KineticGasLiquidLink(...)],
    property_solvers=[SpeciationPropertySolver(...)],
)
result = cv.advance(dt_h=0.01, chem_env=chem_env_dict)
```

---

### Phase 4 — Update tests

#### `tests/standalone/test_cv_advance.py`
- Remove tests that set or assert on `_last_properties` directly
- Add tests confirming `AdvanceResult.properties` is populated from the pre-step state
- Add end-to-end test: single CV with gas + liquid phases and `KineticGasLiquidLink`
  as an internal interface

#### `tests/standalone/test_gas_liquid_link.py`
- Remove manual `_last_properties` injection (lines 150, 283)
- Replace with passing speciation via the `property_results` parameter to
  `compute_flux()`

#### `tests/standalone/test_property_solvers.py`
- Remove `test_cv_last_properties_initially_empty` — concept no longer exists

---

### Phase 5 — Simplify `GasLiquidVolume` *(follow-on, not blocking)*

Once Phase 3 is complete, `GasLiquidVolume` can be reduced to a thin factory that
constructs a single `ControlVolume` with both phases and `KineticGasLiquidLink` as
an internal interface. The following can then be removed:

- The two-CV split (`_cv_gas`, `_cv_liquid`)
- `MultiCVSystem` coordination wrapper
- `EulerSnapshotSolver`'s bespoke two-CV coordination logic
- `ScipyODESolver`'s manual CV key references

This is the largest change and is best done as a separate PR once Phases 1–4 are
validated by the test suite.

---

## Risk Assessment

| Phase | Risk | Reason |
|-------|------|--------|
| 1 | Low | No existing CV has both interfaces and property solvers — reordering has no effect on current behaviour. The `context` → `chem_env` rename is mechanical find-replace. |
| 2 | Low | Removing workarounds that only exist to patch the ordering problem |
| 3 | Low | Purely additive (dual-protocol) — no existing call sites change |
| 4 | Low | Test updates follow directly from the code changes |
| 5 | Medium | Largest structural change; requires redirecting all `GasLiquidVolume` users |

Phases 1–4 can be merged together. Phase 5 is best done separately once the test
suite confirms the new `ControlVolume` path is correct.
