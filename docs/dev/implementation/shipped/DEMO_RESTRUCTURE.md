# Demo Restructure — Planning Note

> **⚠ Reconciliation banner (2026-05-29).** **Shipped, with
> deviations from option β below.** The restructure executed
> during the 2026-05-29 session split what this doc called
> `demos/api/` into two thematic subtrees inside a renamed
> `demos/model_api/` parent. The four file moves and the Q1 /
> Q2 demo writes did happen as scoped; the directory layout
> diverged.
>
> **What shipped:**
>
> ```
> demos/
>   README.md
>   _bootstrap.py
>   builder/                       (4 fermenter demos, README)
>     batch_fermenter.py
>     cstr_fermenter.py
>     fed_batch_fermenter.py
>     microplate_fermenter.py
>   model_api/                     (renamed from option β's api/)
>     README.md
>     chemistry/                   (reaction defs, FBA nested here)
>       reaction_system.py         ← Q2 (importable factories
>                                       + runnable inspector)
>       fba/
>         fba_toy.py
>         fba_ecoli_core.py
>         ecoli_core.json
>     model_construction/          (manual Simulation assembly)
>       raw_construction.py        ← Q1 (imports from
>                                       ../chemistry/reaction_system)
> ```
>
> **How it deviated from option β:**
>
> - `api/` was renamed `model_api/` and stopped being a flat
>   container for the new demos. Instead it became an umbrella
>   over two peer subtrees, `chemistry/` and `model_construction/`,
>   that name **what's being defined** rather than **what layer**.
> - The Q6 FBA demos moved from `api/fba/` to
>   `chemistry/fba/` — FBA is now categorised as chemistry
>   (`BlackBoxReactionModel` definitions), not a parallel demo
>   class.
> - Three READMEs shipped (top-level `demos/`, `demos/builder/`,
>   `demos/model_api/`), not the two scoped here. The leaf
>   subfolders (`chemistry/`, `model_construction/`) inherit
>   context from `model_api/README.md`.
> - The Q1 demo (`raw_construction.py`) **imports** chemistry
>   from the Q2 module rather than redeclaring it — a pattern
>   that wasn't in the original scope but landed naturally given
>   the chemistry/construction split.
>
> **What's still valid:**
>
> - The principle that **top-level subtrees name which framework
>   layer the demo teaches** holds, with `model_api/` being the
>   layer rather than the layer's split.
> - The CI ask was a no-op — nothing in the repo (`.github/`,
>   `Makefile`, `scripts/`, `tox.ini`, `conftest.py`,
>   `pyproject.toml`) automatically walks `demos/`. Verified
>   2026-05-29.
> - Open design question 1 (no aliases, one-shot move) was
>   honoured. Open design question 2 (subdivision of `api/`)
>   was answered yes via the chemistry/construction split.
>   Open design question 3 (builder README) was honoured.
> - The forward-looking demos still scoped to land here —
>   `export_results.py` from
>   [RESULT_EXPORT.md](RESULT_EXPORT.md), `chemistry_database.py`
>   from [CHEMISTRY_UNIFICATION_PLAN.md](CHEMISTRY_UNIFICATION_PLAN.md)
>   Phase 3b — should target `model_api/` subtrees, not the now-
>   defunct `api/` path. Update those notes if they haven't been
>   already.
>
> The body below is the original 2026-05-28 plan, retained for
> historical context.

---

> **Status (original):** Trigger-gated. Surfaced 2026-05-28 during the
> post-simulation-class exploration session. The new
> `demos/api/` subtree is being populated as new demos are
> written; the migration of the four existing fermenter demos
> to `demos/builder/` is the remaining work. Pick up when one
> of the trigger conditions below fires.

## Context

The `demos/` directory is currently flat:

```
demos/
  batch_fermenter.py        (builder-based)
  cstr_fermenter.py         (builder-based)
  fed_batch_fermenter.py    (builder-based)
  microplate_fermenter.py   (builder-based)
```

All four use `FermenterBuilder` / `FermenterFactory` as the
entry point. None demonstrate the lower layers
(`ControlVolume` + `Phase` + `ReactionSystem` +
`KineticGasLiquidLink` + `Simulation` constructed directly).
A new user reading the demos cannot tell from the directory
listing which layer of the framework they're learning.

The exploration session identified three new direct-API demos
that would land naturally in a dedicated subtree:

- A "build a model without the builder" walkthrough (the Q1
  demo from the exploration session — `raw_construction.py`).
- A `ReactionSystem` declaration walkthrough (Q2 —
  `reaction_system.py`).
- FBA / `BlackBoxReactionModel` demos under
  `demos/api/fba/` (Q6 — `fba_toy.py`,
  `fba_ecoli_core.py`). These have been written and live
  there now.

A later demo (`chemistry_database.py`) lands when
`chemistry-unification-3b` ships; see
[CHEMISTRY_UNIFICATION_PLAN.md](CHEMISTRY_UNIFICATION_PLAN.md)
Phase 3b mandate.

## The restructure (option β from the exploration discussion)

Target layout:

```
demos/
  README.md                     ← which layer to start with, when
  builder/
    batch_fermenter.py
    cstr_fermenter.py
    fed_batch_fermenter.py
    microplate_fermenter.py
  api/
    raw_construction.py         ← Q1 (to write)
    reaction_system.py          ← Q2 (to write)
    export_results.py           ← from RESULT_EXPORT.md
    fba/
      fba_toy.py                ← Q6 (shipped 2026-05-28)
      fba_ecoli_core.py         ← Q6 (shipped 2026-05-28)
      ecoli_core.json           ← data for the above
    chemistry_database.py       ← Q7, lands with chem-unification-3b
```

The principle: top-level subfolders name **which framework
layer** the demo teaches.

- `builder/` — quick-start path; user touches
  `FermenterBuilder` and gets a configured `Simulation` back.
  The current four demos are the canonical examples.
- `api/` — direct framework usage; user constructs `Phase`s,
  `ReactionSystem`s, `ControlVolume`s, `Simulation`s
  themselves. Pedagogical "what's underneath" demos plus
  demos that don't fit the fermenter pattern (FBA, future
  custom-CV examples).

A top-level `demos/README.md` (currently absent) explains
which subtree a new user should start with.

## Scope when picked up

### Already in flight (shipped 2026-05-28)
- `demos/api/fba/fba_toy.py`
- `demos/api/fba/fba_ecoli_core.py`
- `demos/api/fba/ecoli_core.json`

### File moves
- `demos/batch_fermenter.py` → `demos/builder/batch_fermenter.py`
- `demos/cstr_fermenter.py` → `demos/builder/cstr_fermenter.py`
- `demos/fed_batch_fermenter.py` → `demos/builder/fed_batch_fermenter.py`
- `demos/microplate_fermenter.py` → `demos/builder/microplate_fermenter.py`

### New writes (pure pedagogy, no framework changes)
- `demos/api/raw_construction.py` (~150 LOC) — direct
  construction walkthrough mirroring `batch_fermenter.py` but
  bypassing the builder.
- `demos/api/reaction_system.py` (~120 LOC) — declare
  aerobic-growth `KineticReaction`, an acid-base
  `EquilibriumReaction`, and a cross-phase CO₂ equilibrium;
  show pre-bucketing inspection
  (`system.kinetic_reactions`, `system.single_phase_equilibria`,
  `system.cross_phase_equilibria`,
  `system.blackbox_models`).
- `demos/README.md` (~50 LOC) — short orientation: when to
  start with `builder/`, when to read `api/`.

### Demo-CI updates
- Whatever script runs `demos/*.py` in CI needs to walk both
  `demos/builder/*.py` and `demos/api/**/*.py`.

## Open design questions

1. **Should the existing fermenter demos stay top-level as
   aliases until docs catch up?** Probably not — a one-shot
   move with a short note in CHANGELOG is cleaner than a
   long-lived symlink/duplicate pattern. The repo doesn't
   have CHANGELOG conventions today, so a note in
   `README.md` at the top level suffices.
2. **Does `demos/api/` warrant further subdivision?** Today's
   list:
   - `fba/` exists for thematic grouping (FBA-specific demos
     + data).
   - `raw_construction.py`, `reaction_system.py`,
     `export_results.py`, `chemistry_database.py` are flat.
   When `api/` accumulates more than ~5-6 flat files,
   sub-grouping (e.g. `api/construction/`, `api/recording/`,
   `api/chemistry/`) may pay off. Defer the question until
   then.
3. **Should `demos/builder/` get its own README?** Probably
   yes once the move happens — a short list of the four
   builder demos and what each one covers (batch vs CSTR vs
   fed-batch vs microplate).

## Trigger conditions

Any one of:

1. **A natural pause where the demo subtree can be
   reorganised without disrupting in-flight work.**
   Touching the demo files mid-phase risks merge conflicts;
   between phases is the right window.
2. **Q1 / Q2 demos being written.** Writing
   `raw_construction.py` and `reaction_system.py` forces the
   `demos/api/` subtree to exist (it already does via the
   FBA demos), at which point the asymmetry with the
   top-level builder demos becomes the visible friction.
3. **A new contributor needs orientation.** The point of
   the restructure is mostly external clarity; if someone
   joins the project and asks "which demo should I start
   with?", the answer is awkward today.

## Relationship to other phases

- **Independent of** every shipped or trigger-gated phase. The
  restructure is purely organisational; no source code
  changes.
- **Provides homes for** the new pedagogical demos surfaced in
  the 2026-05-28 exploration session (Q1, Q2, Q6, Q7) and
  the `export_results.py` demo accompanying
  [RESULT_EXPORT.md](RESULT_EXPORT.md).
