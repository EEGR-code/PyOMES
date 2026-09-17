# PHController Corrector Chemistry Validation — Implementation Plan

> Status: design note, not yet started. No branch, no checklist, no code
> yet. Written on `main` (independently of the `tutorials-followups` phase
> that surfaced the motivating bug) so it doesn't get bundled into that
> phase's diff — see "How to start one" below when picked up.

## Motivation

While fixing `docs/tutorials/D2C_workshop/raw_construction.py` (tutorials-followups
checkpoint 3), a silent-failure class of bug surfaced: `PHController` was
configured with `chemical_id="H3PO4"` (acid corrector) and
`base_chemical_id="NaOH"` (base corrector). `NaOH` is a recognised
strong-corrector alias — `ControlVolume.apply_external_flux`
(`PyOMES/core/control_volume.py:495-527`) resolves it to `Na+` via the
`_STRONG_CORRECTOR_ION` map (`control_volume.py:35-38`), which shifts pH
purely through charge balance, no declared chemistry needed. `H3PO4` is
**not** in that map — it only shifts pH by actually dissociating, which
requires the phosphate ladder (`H3PO4 ⇌ H2PO4- ⇌ HPO4-- ⇌ PO4---`) to be
declared as `EquilibriumReaction`s in the CV's `reaction_system`. The
script's manually-built `ReactionSystem` never declared it, so every dose
of `H3PO4` accumulated as inert neutral acid (confirmed: 0.356 mol by
t=5h) with zero effect on pH — the acid half of the PI loop was a silent
no-op, and pH ran away to 12.089 with nothing able to correct it back
down.

**The interesting part: this exact check already exists in the codebase,
just not on the path that caught the bug.**
`ControlVolume.equilibrate_to_pH` (`control_volume.py:1190-1206`) validates
its own `corrector_id` argument before running: if it's not a
strong-corrector alias, it must appear in some declared
`EquilibriumReaction`'s stoichiometry, or a `ValueError` is raised. That
guard is only exercised by `equilibrate_to_pH`'s standalone one-shot
probing utility — `PHController`'s actual per-step dosing path
(`compute()` → `Simulation._apply_controller_action` →
`cv.apply_external_flux`) never runs it. Porting the same check onto
`PHController` would have caught this bug immediately instead of letting
it silently produce a 12.089 pH.

## Design

- **Scope: `PHController`-specific, not a general dosing check.** The new
  validation only ever inspects `PHController.chemical_id` and
  `PHController.base_chemical_id` — the two fields the class already owns
  for this exact purpose. It has no visibility into, and no reason to
  touch, anything dosed by other controllers or boundaries (a nutrient
  `LiquidFeed`, a future non-pH controller reusing `flux_applied`, etc.).
  Confirmed `PHController` is the only CV-native controller in
  `cv_loops.py` with a "dose a raw chemical id" pattern —
  `DOAgitationController` doses RPM/kLa, not a chemical.

- **New method:** `PHController.validate_against_cv(cv: ControlVolume) -> None`
  in `PyOMES/control/cv_loops.py`. For each configured id
  (`chemical_id`, and `base_chemical_id` if set), checks:
  1. Is it a strong-corrector alias (`_STRONG_CORRECTOR_ION`)? If yes, OK.
  2. Otherwise, does it appear in `cv.reaction_system.single_phase_equilibria`'s
     stoichiometry? (Same set `equilibrate_to_pH` already computes at
     `control_volume.py:1196-1200`.) If yes, OK.
  3. Otherwise, warn — see Open Question 1 for warn-vs-raise.

- **Trigger point:** `Simulation.run()`'s existing per-run controller setup
  loop:
  ```python
  # simulation.py:759-761 (current)
  for ctrl in self._controllers:
      if hasattr(ctrl, "reset"):
          ctrl.reset()
  ```
  Extend this loop to also call `ctrl.validate_against_cv(cv)` when the
  controller exposes that method (duck-typed via `hasattr`, mirroring the
  `reset()` check immediately next to it) and its target CV can be
  resolved. Runs once per `Simulation.run()` call, not once per step —
  matches `reset()`'s own once-per-run cadence and avoids re-emitting the
  same warning on every one of (e.g.) 1000 steps.

- **CV resolution:** reuse the exact `target_cv_key` resolution pattern
  already used to dispatch `ctrl.compute()`
  (`simulation.py:923-930`): `target_cv_key` set → that CV; `target_cv_key`
  is `None` and `self.cvs` has exactly one entry → that CV; `None` with
  multiple CVs → ambiguous, skip validation silently (matches how
  `PHController._resolve_cv` itself only raises when a real snapshot is
  actually dispatched to it, not proactively).

- **Import-cycle check needed:** `_STRONG_CORRECTOR_ION` lives in
  `PyOMES/core/control_volume.py`; `PHController` lives in
  `PyOMES/control/cv_loops.py`, which already imports from
  `PyOMES.core.snapshot` — so `core` → `control` isn't a new direction.
  Still need to confirm `control_volume.py` never imports anything from
  `PyOMES.control` before adding the reverse import (should be clean per
  Pattern 1's design, but verify before implementing).

## Open questions (need a decision before implementing)

1. **Warn vs. raise.** `equilibrate_to_pH`'s existing precedent is a hard
   `ValueError`. Three options:
   - **(a) Always warn** (`UserWarning`, non-fatal) — matches the
     "printed once, ignorable" ethos of `AccuracyWarning`/`ConservationWarning`.
   - **(b) Always raise** — fails fast, matches the `equilibrate_to_pH`
     precedent exactly, but is a behavior-breaking risk if any
     currently-working setup relies on a corrector that isn't declared in
     equilibria for a reason this check can't see.
   - **(c) Warn by default, opt-in `strict=True` to raise** — safest
     migration path, more surface area.

   **Recommendation: (a)** for a first cut — smallest, safest, immediately
   catches the bug class without any risk of breaking an existing working
   setup. (c) is a natural follow-up if warnings prove too easy to miss in
   practice.

2. **New warning category, or plain `UserWarning`?** `AccuracyWarning` and
   `ConservationWarning` are the two existing `UserWarning` subclasses
   wired into `PyOMES.config.warnings`'s throttle machinery
   (`silent`/`once`/`first_N`). A new category (e.g.
   `PHControllerConfigWarning`) would need the same throttle wiring for
   consistency.

   **Recommendation:** plain `UserWarning` for a first cut — the check
   fires at most once per controller per `Simulation.run()` call (not
   per-step), so the throttle machinery that exists specifically to tame
   per-step spam isn't really needed here.

3. Confirmed no other `cv_loops.py` controller needs this (see Scope
   above) — flagging in case that's wrong.

## Implementation steps (once the questions above are resolved)

1. Add `PHController.validate_against_cv(cv)` to `cv_loops.py` — the
   three-step check described above; warning message should name the
   missing corrector and what's needed to fix it, mirroring
   `equilibrate_to_pH`'s existing message wording
   (`control_volume.py:1202-1206`) for consistency.
2. Wire the call into `Simulation.run()`'s per-run controller loop
   (`simulation.py:757-761`), duck-typed via `hasattr`, using the CV
   resolution described above.
3. Unit tests:
   - `PHController(chemical_id="H3PO4")` on a CV with no phosphate
     equilibria → warning fires.
   - Same, but with the phosphate ladder declared → no warning.
   - `chemical_id="NaOH"` / `"KOH"` → no warning even with no declared
     equilibria (strong-corrector path).
   - `target_cv_key=None` with multiple CVs → no crash, validation
     silently skipped.
4. Regression check: full suite, plus the tutorials that exercise
   `PHController` (`raw_construction.py` post-fix,
   `templates/batch_fermenter.py`) — confirm zero new spurious warnings on
   already-correct chemistry.
5. One-line cross-reference added to `PHController`'s docstring pointing
   at the new check.

## How to start one

Per this folder's usual convention
([README.md](README.md#how-to-start-one)): write a checklist file
(`PHCONTROLLER_CORRECTOR_VALIDATION_CHECKLIST.md`), cut a branch off
`main` (suggested name: `phcontroller-corrector-validation`), work the
checkpoints there, and add a "Currently in flight" pointer to this
folder's `README.md`. Small and self-contained — no dependency on any
other in-flight phase, ready to pick up any time.
