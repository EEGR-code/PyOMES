# Framework Polish Bundle — Design Note

> **Status:** **Shipped 2026-06-01.** Tag `framework-polish-shipped`
> (merge `7643669`). P1–P9 all implemented; 1021/0 tests.
> Notable: P8 also fixed a silent bug — `gas_feed.y.<species>`
> param path was discarding writes because `GasFeed.y` returned
> a copy.

## Context

The 2026-05-27 audit surfaced several small surgical issues
that share three traits: each is independently testable, none
requires a design decision, and the whole set is too small to
warrant its own phase per item but too useful to leave
scattered. This note collects them so the bundle has a
persistent home and an unambiguous scope.

Items requiring design discussion (param-path dispatcher,
property-calculator phase-agnostic runner, snapshot/recorder
phase-agnostic surfaces) are **not** in this bundle — they have
their own planning notes:

- [`PARAM_PATH_DISPATCHER.md`](PARAM_PATH_DISPATCHER.md)
- [`PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md`](PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md)
- [`SNAPSHOT_PHASE_AGNOSTIC.md`](SNAPSHOT_PHASE_AGNOSTIC.md)

## The bundle (P1–P9)

Approximate total surface: **~200 LOC across ~8 files**. No API
changes for external callers. No BSM2 sentinel re-baseline
expected (these are framework-quality, not numerical, edits).

### P1 — Fail-loud initial speciation solve  *(severity: HIGH)*

[`src/core/simulation.py:708-712`](../../src/core/simulation.py#L708-L712)
— `Simulation._initial_solve` catches `(TypeError, ValueError,
AttributeError, KeyError)` and `pass`es. Any failure of the
initial solve leaves `pH[0]`/`I[0]` as NaN and the run
continues silently.

Fix: replace the `pass` with a single warning emit through
`VLsim.config.warnings.emit("speciation_initial_solve_failed", ...)`.
Use the existing `AccuracyWarning` machinery so users can
silence via the standard env-var preset if desired.

Estimated: ~10 LOC + 1 test asserting the warning fires when
the engine raises.

### P2 — `_build_reaction_environment` gas-only CV defaults  *(severity: MEDIUM)*

[`src/core/control_volume.py:608-623`](../../src/core/control_volume.py#L608-L623)
— when no `"liquid"` phase is present, `V_L` defaults to `1.0`
and `T_K` to `298.15`. Rate laws reading `env.V_L` /
`env.T_K` get dimensionally meaningless values silently.

Fix: prefer the gas phase's `V_L` / `T_K` when no liquid is
present; raise if neither phase is present (the CV has no
state for the rate law to read at all).

Estimated: ~15 LOC + 2 tests (gas-only-CV concentration
sanity, both-absent raise).

**Interaction with
[PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md](PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md).**
P2 and that note share the same "framework handles CVs
without a liquid phase" thread at different layers (reaction
env vs property runner). Same conceptual root, same testing
infrastructure (a gas-only CV fixture). If both ship, do them
on the same branch.

### P3 — Cache species-vector map per-`ReactionSystem`  *(severity: LOW)*

[`src/core/control_volume.py:718-756`](../../src/core/control_volume.py#L718-L756)
— `_build_species_vector_map` is called every `advance()` step
and runs a full dummy `compute_rates` plus a sort. Species
membership is invariant for the life of a `ReactionSystem`
(lifecycle gating prevents mid-run mutation).

Fix: cache the result on `self` keyed by
`id(self._reaction_system)`; invalidate when
`_set_reaction_system_unchecked` is called.

Estimated: ~20 LOC + 1 test (assert single-call per run with a
spy).

### P4 — Rewrite `gas_liquid_link.py` module docstring  *(severity: MEDIUM)*

[`src/core/gas_liquid_link.py:22-47`](../../src/core/gas_liquid_link.py#L22-L47)
— the module docstring references deleted machinery:
`PropertyResult.alphas`, `property_results` kwarg,
`speciation_override` kwarg, `MultiCVSystem`. The actual code
reads alphas inline from `phase.n_mol` via `_CANONICAL_NAMES`.

Fix: rewrite the docstring around the current code path.
Documentation-only change.

Estimated: ~25 LOC (mostly deletion + rewrite of the example
block).

### P5 — Fix `cv.advance()` docstring  *(severity: LOW)*

[`src/core/control_volume.py:400-401`](../../src/core/control_volume.py#L400-L401)
— the docstring says `Owner is the orchestrator (MultiCVSystem
until SIMULATION_CLASS ships)`. `MultiCVSystem` was deleted at
SIMULATION_CLASS C14 and that phase has shipped.

Fix: replace with `Owner is the Simulation orchestrator`.
Also sweep for any other "until SIMULATION_CLASS ships"
stragglers in `src/core/`.

Estimated: ~5 LOC.

### P6 — Apply architecture.md follow-up edits  *(severity: HIGH)*

The three edits from the audit follow-up
([architecture.md](../architecture.md):
"Equilibrium pathways" rewrite, tense fix on SIMULATION_CLASS
line, "Directory layout" replacement) — **already applied
2026-05-27 on `main`** as part of the audit follow-up commit
that includes this note. Listed here so the bundle scope is
complete on paper.

Estimated: ~80 LOC (already landed).

### P7 — Refresh `CUFERMENTER_SUNSET.md` demos paragraph  *(severity: LOW)*

[`CUFERMENTER_SUNSET.md:43-47`](CUFERMENTER_SUNSET.md#L43-L47)
— stated the 3 fermenter demos were broken at module load.
Commit `5fd8581` migrated them; all four demos now run.

**Already applied 2026-05-27** as part of the audit
follow-up. Listed for completeness.

### P8 — Pattern B completion on boundary classes  *(severity: MEDIUM)*

[`src/core/boundaries.py`](../../src/core/boundaries.py)
— `GasFeed.y` and `GasFeed.vvm_min` lack `_set_y_unchecked` /
`_set_vvm_min_unchecked` siblings. The orchestrator's
`_apply_param_change` writes them directly
([simulation.py:613, 622](../../src/core/simulation.py#L613)),
inconsistent with the gated pattern used elsewhere
(`_set_T_K_unchecked`, `_set_kLa_unchecked`).

Fix: add gated public setters + unchecked siblings on
`GasFeed`, and sweep other boundary classes for the same gap.
The orchestrator can switch to the unchecked path once the
siblings exist.

Estimated: ~30 LOC + 3 tests (gate fires, unchecked bypasses,
orchestrator uses unchecked).

**Interaction with the param-path dispatcher rewrite:** the
[descriptor-based design in PARAM_PATH_DISPATCHER.md](PARAM_PATH_DISPATCHER.md)
replaces the `@property` + `_set_<name>_unchecked` pair on each
gated attribute with a single descriptor declaration. If the
dispatcher rewrite happens first, P8 dissolves entirely (no
need to add unchecked siblings that the descriptor would
immediately replace). If P8 ships in this polish bundle first,
the descriptor migration consumes the new unchecked setters
during its own migration step. Either ordering is fine; do not
do P8 in parallel with the dispatcher rewrite.

### P9 — ZOH-hold placeholder disambiguation  *(severity: LOW)*

[`src/core/simulation.py:516-525`](../../src/core/simulation.py#L516-L525)
— when a controller hasn't fired yet (sample-period gating
hasn't elapsed and `last_action is None`), the orchestrator
appends a default-constructed `ControlAction` to the per-step
record. Recorder consumers can't distinguish "never fired" from
"fired this step with no effect".

Fix: return `None` from `_invoke_controllers` for never-fired
slots; have the recorder skip `None` entries. Alternative:
keep the placeholder but flag it with a
`is_placeholder: bool = False` field on `ControlAction`.

Estimated: ~15 LOC + 1 test.

**Interaction with
[RUN_HISTORY.md](RUN_HISTORY.md).** RUN_HISTORY's "richer
default recorder" reframing would necessarily revisit this
data-model question. If P9 ships in this polish bundle first,
RUN_HISTORY's recorder design consumes the answer; if
RUN_HISTORY goes first, P9 dissolves into whatever shape that
phase picks. Do not ship P9 in parallel with a RUN_HISTORY
checkpoint that's also touching the per-step controller
record.

## Items deliberately excluded

- **Param-path dispatcher rewrite** — needs Q1–Q5 design
  discussion; see PARAM_PATH_DISPATCHER.md.
- **PropertyCalculator gas-phase runner** — needs the
  per-phase question pinned; see
  [`PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md`](PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md).
- **Snapshot/Recorder phase-key generalisation** — needs
  alias-as-property vs eager-mirror decision; see
  [`SNAPSHOT_PHASE_AGNOSTIC.md`](SNAPSHOT_PHASE_AGNOSTIC.md).
- **`peng_robinson.py` unreachable, `coupled.py` orphan
  function** — folds into `CUFERMENTER_SUNSET`.
- **`_CANONICAL_NAMES` / `_TOTAL_KEY_OVERRIDES` hardcoding** —
  folds into `chemistry-unification-3b`.

## Trigger conditions

Any one of:

1. **Immediate.** P1 (fail-loud initial solve) catches a
   real bug. The other items bundle alongside.
2. **Convenience.** A free afternoon — the bundle is small
   enough to ship in one branch without disrupting other work.
3. **Bundled with another phase that touches the same
   files.** P3 and P5 sit in `control_volume.py`; if that
   file is being edited anyway, fold them in.

## Branching convention

Per the standard pattern
([upcoming/README.md](README.md) §Branching and tagging
convention): one branch `framework-polish` off `main`, all P1–P9
implementation on it, `--no-ff` merge with tag
`framework-polish-shipped`. No checklist file needed — this
note IS the checklist; tick items inline as they land.

## Relationship to other phases

- **Independent of** all five trigger-gated phases
  (`CHEMISTRY_UNIFICATION_PLAN`, `CONTAINER_LAYERING`,
  `RUN_HISTORY`, `SPECIATION_LEVEL_RETIREMENT`,
  `CUFERMENTER_SUNSET`).
- **Prerequisite for** the param-path dispatcher rewrite via
  P8 (Pattern B completion on boundary classes).
- **Naturally bundles with** the
  PROPERTY_CALCULATORS_PHASE_AGNOSTIC and
  SNAPSHOT_PHASE_AGNOSTIC work if a single "generalisation +
  polish" pass is taken.
- **Follows** `simulation-class` (shipped 2026-05-27).
