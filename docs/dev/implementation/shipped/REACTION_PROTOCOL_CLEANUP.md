# Reaction Protocol Cleanup — Plan

> **Status: Shipped 2026-05-12** — `integration_mode` rate/delta
> dispatch deleted. 708 standalone tests passing post-cleanup (down
> from 712: 4 dead-path tests removed). Tag:
> `reaction-protocol-cleanup-shipped`.

## Role

Precursor cleanup branch that shipped **before** `chemistry-unification-1`
began. **Not** strictly part of chemistry unification — fixed a leaky
abstraction in the `Reaction` protocol that happened to share files
(`reaction.py`, `reaction_set.py`, `protocols.py`, `control_volume.py`,
`solvers.py`) with the equilibrium-reaction work, so landing it first
de-risked Phase 1 of chemistry unification by removing one moving part
from those dispatch sites.

## Motivation

Today `Reaction`, `ReactionSet`, and `BlackBoxReactionModel` carry an
`integration_mode: str` field with values `"rate"` (rate_fn returns
mol/h) or `"delta"` (rate_fn returns mol per timestep). The two modes
are dispatched on at three sites in the orchestrator:

- [`control_volume.py:346–354`](../../src/core/control_volume.py#L346) —
  delta-mode applies sources atomically via `_apply_deltas`, rate-mode
  integrates them as derivatives.
- [`solvers.py:241–256`](../../src/core/solvers.py#L241) — `EulerSnapshotSolver`
  treats delta as fake-rate by dividing by `dt_h`, then integrates as
  a derivative.
- [`solvers.py:483–484`](../../src/core/solvers.py#L483) — `ScipyODESolver`
  caches `is_delta = (rxn_mode == "delta")` and branches on it inside
  the hot loop.

The two solver paths handle delta-mode **inconsistently** —
`control_volume.advance()` honours delta-as-atomic-application, while
`solvers.py` treats delta as fake-rate. Identical reaction declarations
produce different numerical behaviour depending on which solver is
attached. Real bug-magnet.

`"delta"` mode was originally added to support external models that
internally integrate over `dt` and emit net mole changes. But in the
current codebase **nothing actually uses it**:

- BSM2 ([`bsm2_direct.py:493`](../../models/vlmodels/adm1/bsm2_direct.py#L493))
  declares `integration_mode = "rate"` and returns `mol/h`.
- The only `integration_mode="delta"` references in the codebase are
  three dead-code path tests:
  [`test_delta_mode_applies_directly`](../../tests/standalone/test_cv_advance.py#L277),
  [`test_integration_mode_delta`](../../tests/standalone/test_reactions.py#L329),
  and the `integration_mode = "delta"` declaration at
  [`test_cv_advance.py:285`](../../tests/standalone/test_cv_advance.py#L285)
  inside a synthetic `DeltaModel` test class.

The flag is a feature for zero current users.

## Approach

Delete delta mode entirely. The reaction protocol becomes single-mode:
`compute_rates(env)` returns `mol/h`, full stop. No
`integration_mode` field on any class. No dispatch branches in CV or
solvers. No `_apply_deltas` method.

Future delta-emitting external models, if they ever appear, write a
five-line wrapper that divides their output by their step size and
return rates — that wrapper lives with the model, not in the generic
adapter.

This is consistent with the project rule from prior phases: "don't
add features for hypothetical future requirements." Adding back a
delta path would be one branch when a real consumer appears.

## Scope

1. **Remove `integration_mode` field** from:
   - `Reaction` ([`reaction.py:61`](../../src/reactions/reaction.py#L61))
   - `ReactionSet` ([`reaction_set.py:35`](../../src/reactions/reaction_set.py#L35))
   - `BlackBoxReactionModel` constructor signature, body, and docstring
     ([`blackbox.py:84, 95–98, 113, 121`](../../src/reactions/blackbox.py#L84))
   - The protocol method at
     [`protocols.py:55`](../../src/reactions/protocols.py#L55).

2. **Remove dispatch branches**:
   - [`control_volume.py:346–354`](../../src/core/control_volume.py#L346)
     collapses to the rate-only path (the `else` branch becomes
     unconditional).
   - [`solvers.py:241–256`](../../src/core/solvers.py#L241) drops the
     `rxn_mode == "delta"` division-by-dt branch.
   - [`solvers.py:483–484`](../../src/core/solvers.py#L483) drops the
     `is_delta` cache and any subsequent branching on it (verify hot
     loop downstream).

3. **Remove `_apply_deltas`** from
   [`control_volume.py:589–604`](../../src/core/control_volume.py#L589) —
   confirmed sole caller is the delta-mode dispatch site (grep over
   the codebase).

4. **Remove BSM2's dead `integration_mode` property** at
   [`bsm2_direct.py:492–494`](../../models/vlmodels/adm1/bsm2_direct.py#L492).
   Was always `"rate"` and is now ignored entirely.

5. **Delete dead-path tests**:
   - [`test_integration_mode_default`](../../tests/standalone/test_reactions.py#L134) —
     tested that default field value was `"rate"`. No field anymore.
   - [`test_integration_mode`](../../tests/standalone/test_reactions.py#L277) —
     same for ReactionSet.
   - [`test_integration_mode_delta`](../../tests/standalone/test_reactions.py#L329) —
     tested that the field could be set to `"delta"`. Path deleted.
   - [`TestAdvanceDeltaMode.test_delta_mode_applies_directly`](../../tests/standalone/test_cv_advance.py#L277)
     — entire class can be removed (one method).
   - The `integration_mode = "rate"` line in
     [`test_factory.py:231`](../../tests/standalone/test_factory.py#L231)
     and [`test_builder.py:261`](../../tests/standalone/test_builder.py#L261)
     — removed if they're setting the attribute on a custom test class;
     no harm left if the test class no longer needs it.

6. **Update class diagrams** at
   [`docs/class_diagrams.md:363, 371, 379, 389`](../class_diagrams.md#L363) —
   drop `+integration_mode: str` from each `Reaction`-flavoured class
   (Reaction, ReactionSet, BlackBoxReactionModel, and one other —
   confirm at edit time).

## Out of scope

- Anything related to chemistry unification proper (`kind`/`log_K`/
  charge balance/equilibrium reactions/staleness machinery/database
  packaging). All of that lives in the `chemistry-unification-N`
  branches.
- Adding `dt_h` to `ReactionEnvironment`. No consumer needs it after
  delta mode is gone. Revisit only when a real consumer appears.
- Renaming any other field on `ReactionEnvironment` or `Reaction`.

## Branch and tag

- Branch: `reaction-protocol-cleanup` (off `main`)
- Tag at ship: `reaction-protocol-cleanup-shipped`

Per the convention in [README.md](README.md), this is a substantive
phase and gets its own branch + `--no-ff` merge + tag.

## Risks

- **No production numerical regression.** Confirmed BSM2 declares
  `"rate"` and was never on the delta path. The synthetic `DeltaModel`
  test class is the only "user."
- **Hidden delta consumers in `models/` directory beyond BSM2.** Mitigation:
  grep for `integration_mode` in `models/` before declaring done.
  Already grepped the whole tree; no further consumers found.
- **`solvers.py:483` `is_delta` cache is referenced downstream.** Need
  to verify the full extent of branching before deletion. The grep
  shows two references within solvers.py; both should collapse to the
  rate path cleanly. Verify at edit time.

## Estimated size

~100–200 LoC, almost all deletions. One focused branch; expected to
ship in a single session after the design is locked.
