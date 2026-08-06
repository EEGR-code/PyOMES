# Reaction Protocol Cleanup — Checklist

> **Status: Shipped 2026-05-12** — all 7 checkpoints landed. 708
> standalone tests passing post-cleanup (down from 712: 4 dead-path
> tests removed). Tag: `reaction-protocol-cleanup-shipped`.

Working checklist for the precursor cleanup that shipped before
`chemistry-unification-1`. Conceptual framing is in
[REACTION_PROTOCOL_CLEANUP.md](REACTION_PROTOCOL_CLEANUP.md). This
file is the implementation plan.

## Goal

Delete `integration_mode` from the `Reaction` protocol and every
class that carries it. After this branch ships:

- `Reaction.compute_rates(env)` returns `mol/h` only — no field, no
  flag, no per-call mode.
- The orchestrator (`ControlVolume.advance`) and every solver
  (`EulerSnapshotSolver`, `ScipyODESolver`) integrate sources as
  derivatives unconditionally.
- The `_apply_deltas` path in `ControlVolume` is gone.
- The synthetic delta-mode tests are gone (zero production
  consumers).
- BSM2's dead `integration_mode = "rate"` property is gone.

## Out of scope

- Anything related to chemistry unification (`kind` field, `log_K`,
  charge balance, equilibrium reactions, staleness machinery,
  `ChemistryDatabase`). All of that lives in the
  `chemistry-unification-N` branches.
- Adding `dt_h` to `ReactionEnvironment`. No consumer needs it.
- Renaming any other field on `ReactionEnvironment` or `Reaction`.
- Adding BSM2 regression test coverage (separate concern; BSM2 was
  never on the delta path so this branch introduces no BSM2 risk).

## Resolved decisions

- **Option A: pure deletion** of delta mode, no replacement
  mechanism. Confirmed by user during planning. No `returns_deltas`
  flag, no subclass, no output converter — just delete. Consistent
  with prior phases' "no features for hypothetical users" rule.
- **No `dt_h` on `ReactionEnvironment`.** Without delta mode, no
  reaction adapter needs to know `dt_h`. Add it back if/when a real
  consumer appears.
- **Synthetic delta-mode tests deleted**, not migrated. They tested
  the existence and behaviour of a path that no longer exists; there
  is nothing to preserve.

## Checkpoints

### 1. Remove `integration_mode` from `Reaction`-flavoured classes

- [ ] [`src/reactions/reaction.py:61`](../../src/reactions/reaction.py#L61)
      — delete the class attribute `integration_mode: str = "rate"`.
- [ ] [`src/reactions/reaction_set.py:35`](../../src/reactions/reaction_set.py#L35)
      — delete the same class attribute.
- [ ] [`src/reactions/blackbox.py`](../../src/reactions/blackbox.py)
      — delete the constructor parameter (line 113), the assignment
      (line 121), and the docstring entry for it (lines 95–98). Also
      clean the docstring fragment at line 84 ("(or mol if
      integration_mode='delta')").
- [ ] [`src/reactions/protocols.py:55`](../../src/reactions/protocols.py#L55)
      — delete the `integration_mode` property declaration on the
      `ReactionModel` protocol.

Sanity check: `Grep "integration_mode" src/` returns zero matches
after this checkpoint.

### 2. Remove dispatch branches from `ControlVolume`

- [ ] [`control_volume.py:343–354`](../../src/core/control_volume.py#L343)
      — collapse the `if mode == "delta": ... else: ...` block to the
      rate-only path. The `else` body becomes unconditional; the
      `getattr(...)` mode lookup goes away.
- [ ] [`control_volume.py:589–604`](../../src/core/control_volume.py#L589)
      — delete `_apply_deltas` method entirely.

### 3. Remove dispatch branches from solvers

- [ ] [`solvers.py:240–256`](../../src/core/solvers.py#L240) —
      drop the `rxn_mode == "delta"` branch in `EulerSnapshotSolver`.
      The `getattr(cv.reaction_model, "integration_mode", "rate")`
      lookup goes away. Sources are added to deltas unconditionally
      as rates.
- [ ] [`solvers.py:483–484`](../../src/core/solvers.py#L483) —
      drop the `rxn_mode` lookup and `is_delta` cache. Trace any
      subsequent `is_delta` references in the hot loop and collapse
      them to the rate path. Verify no other branches downstream.

Sanity check: `Grep "is_delta|integration_mode|rxn_mode" src/core/`
returns zero matches after this checkpoint.

### 4. Remove dead `integration_mode` property from BSM2

- [ ] [`models/vlmodels/adm1/bsm2_direct.py:492–494`](../../models/vlmodels/adm1/bsm2_direct.py#L492)
      — delete the `@property def integration_mode(self) -> str:
      return "rate"` block. Now that the protocol no longer carries
      this field, the property is dead code.

Sanity check: `Grep "integration_mode" models/` returns zero matches
after this checkpoint.

### 5. Delete dead-path tests

- [ ] [`tests/standalone/test_reactions.py:134–135`](../../tests/standalone/test_reactions.py#L134)
      — delete `test_integration_mode_default`.
- [ ] [`tests/standalone/test_reactions.py:277–278`](../../tests/standalone/test_reactions.py#L277)
      — delete `test_integration_mode`.
- [ ] [`tests/standalone/test_reactions.py:329–334`](../../tests/standalone/test_reactions.py#L329)
      — delete `test_integration_mode_delta`.
- [ ] [`tests/standalone/test_cv_advance.py:271–298`](../../tests/standalone/test_cv_advance.py#L271)
      — delete the `TestAdvanceDeltaMode` class (one method,
      `test_delta_mode_applies_directly`). Drop the section header
      comment block at lines 271–273.
- [ ] [`tests/standalone/test_factory.py:231`](../../tests/standalone/test_factory.py#L231)
      — read the surrounding context. If the line sets
      `integration_mode = "rate"` on a custom test class, just
      delete the line. If it asserts on a Reaction's integration_mode
      attribute, delete the assertion.
- [ ] [`tests/standalone/test_builder.py:261`](../../tests/standalone/test_builder.py#L261)
      — same treatment as `test_factory.py:231`.
- [ ] Drop the `- advance() with delta-mode (BlackBoxReactionModel)`
      bullet from the module docstring at
      [`tests/standalone/test_cv_advance.py:9`](../../tests/standalone/test_cv_advance.py#L9).

Sanity check: `Grep "integration_mode|delta_mode|delta-mode" tests/`
returns zero matches after this checkpoint.

### 6. Update class diagrams

- [ ] [`docs/class_diagrams.md:363, 371, 379, 389`](../class_diagrams.md#L363)
      — drop the `+integration_mode: str` line from each Reaction-
      flavoured class block. There are four mentions; one per class
      (Reaction, ReactionSet, BlackBoxReactionModel, and one other
      inheritor — confirm which at edit time).

Sanity check: `Grep "integration_mode" docs/` returns zero matches
after this checkpoint.

### 7. Final test sweep + ship

- [ ] Run the full standalone test suite:
      `python -m pytest tests/standalone -q`. Expect 712 minus the
      number of deleted tests (4 individual tests + 1 test class
      method = 5 tests deleted, leaving ~707 passing).
- [ ] If any other test fails, investigate — it should not.
- [ ] Update [README.md](README.md) priority list: remove the
      `REACTION_PROTOCOL_CLEANUP.md` entry; add a one-line note
      under "Currently in flight" indicating chemistry-unification-1
      is next.
- [ ] Move [REACTION_PROTOCOL_CLEANUP.md](REACTION_PROTOCOL_CLEANUP.md)
      and this checklist to [`../phases-shipped/`](../phases-shipped/).
      Add "Status: Shipped" banners to each.
- [ ] Ship via the convention in [README.md](README.md):
      `git checkout main && git merge --no-ff reaction-protocol-cleanup
      -m "Merge reaction-protocol-cleanup: delete integration_mode rate/delta dispatch"`.
- [ ] Tag: `git tag reaction-protocol-cleanup-shipped <commit-hash>`.
- [ ] Push: `git push && git push --tags`.
- [ ] Delete branch:
      `git branch -d reaction-protocol-cleanup` and
      `git push origin --delete reaction-protocol-cleanup`.

## Final test count expectation

712 → ~707 (5 tests deleted: 3 in `test_reactions.py`, 1 in
`test_cv_advance.py`, plus possibly 0–2 from
`test_factory.py`/`test_builder.py` depending on what those lines
do). Final number confirmed empirically at checkpoint 7.
