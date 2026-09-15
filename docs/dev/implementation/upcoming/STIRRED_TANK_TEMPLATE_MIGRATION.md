# Stirred-tank template migration — design discussion

> Status: pre-phase design discussion, 2026-09-15. No branch, no checklist,
> no code yet. Written up from a planning conversation (chat, not a design
> session against code) that audited `models/vlmodels/fermenter/` and its
> naming. When this is picked up: follow `README.md`'s "How to start one" —
> write a checklist file from `PHASE_KICKOFF_TEMPLATE.md`, branch, implement,
> ship.

## Commit discipline for this phase

**The repo owner runs every `git add`, `git commit`, and `git push` for
this phase personally — an assistant executing this plan should not run
those commands.** At each checkpoint below: make the file edits, report
exactly what changed (and run the checkpoint's sanity check), then stop
and wait for the owner to review, stage, and commit before moving to the
next checkpoint. This is so progress can be verified step by step before
it affects the repo, not committed in an automated batch. `git mv` counts
as a staging operation for this purpose — prefer a plain filesystem move
plus a separate `git add` the owner runs, over `git mv` run on the owner's
behalf, unless they say otherwise.

## Goal

Move `models/vlmodels/fermenter/` into `PyOMES` core as
`PyOMES/templates/stirred_tank/`, renaming `FermenterBuilder` →
`StirredTankBuilder` and `FermenterFactory` → `StirredTankFactory` along
the way.

## Why

- **It isn't fermentation-specific.** `FermenterBuilder`/`FermenterFactory`
  build a general reacting, gas-liquid, well-mixed vessel — batch, CSTR,
  and fed-batch are already reachable from the *same* builder purely by
  which boundaries get attached after `.build()` (confirmed: `factory.py`
  has zero topology-specific branching). The name overclaims specificity
  it doesn't have and undersells the generality it does have.
- **It's already flagged as legacy in the code itself.** `builder.py`'s
  `build()` docstring says outright: *"The class and method name
  (`FermenterBuilder.build`) are Phase 7 holdovers — they predate the
  deletion of `GasLiquidVolume`. A cosmetic rename is deferred; bundle it
  with the parallel rename of `FermenterFactory.create_volume`..."* This
  phase is that bundle.
- **It resolves a test-bank ambiguity for free.** Several `test_simulation.py`
  classes (`TestKineticGasLiquidLinkGating`, `TestC9FluxApplyDispatch`,
  `TestC9ParamPathDispatch`, `TestC10TemperatureRamp`, `TestC10VVMSchedule`,
  `TestC10ProfileOrderingBeforeAdvance`, `TestC11FermenterBuilderSimulation`)
  were flagged elsewhere as needing individual triage between "testing core
  via a `vlmodels` fixture" and "testing `vlmodels` itself" — precisely
  because the builder straddled that boundary. Once it's unambiguously
  `PyOMES` core, all of those are unambiguously core tests. Only
  `TestC13ADM1Simulation` (uses `vlmodels.adm1.base`, a genuine external
  model) remains a real triage case, and it's independent of this move.

## Naming decisions already made

- Package: `stirred_tank` (describes the vessel/mixing topology — a single
  well-mixed CV; stays accurate regardless of phase composition).
- Classes: **`StirredTankBuilder`, `StirredTankFactory` — no phase
  qualifier.** Considered `GasLiquidStirredTank...` (today's `factory.py`
  unconditionally constructs both a `GasPhase` and a `LiquidPhase`, no
  branch produces anything else — confirmed by reading the code) but
  rejected it: `PyOMES.core.phases` already has a `SolidPhase`, and
  `ControlVolume` is already phase-key-agnostic, so a third phase getting
  added to *this same class* later (rather than as a sibling) is plausible
  enough that baking "GasLiquid" into the name risks a second rename. The
  asymmetry favors the unqualified name — the cost of being unqualified
  today is a reader checking the docstring; the cost of over-qualifying is
  another rename later. **The current gas+liquid-only constraint must be
  stated explicitly in the class docstring** to compensate for the name no
  longer carrying that information itself.
- `configs.py`'s dataclasses (`VesselConfig`, `GasFeedConfig`,
  `TransferConfig`, `ChemistryConfig`, `OrganismConfig`, `SubstrateConfig`)
  keep their current names — already appropriately generic (or, for
  `Organism`/`SubstrateConfig`, appropriately scoped to what this template
  actually targets). No renames needed there.

## Scope of the move

Everything currently under `models/vlmodels/fermenter/` moves as one unit:

```
models/vlmodels/fermenter/__init__.py            -> PyOMES/templates/stirred_tank/__init__.py (rewritten, not moved verbatim)
models/vlmodels/fermenter/config/__init__.py     -> merged into the above (rewritten, not moved verbatim)
models/vlmodels/fermenter/config/builder.py      -> PyOMES/templates/stirred_tank/builder.py
models/vlmodels/fermenter/config/configs.py      -> PyOMES/templates/stirred_tank/configs.py
models/vlmodels/fermenter/config/factory.py      -> PyOMES/templates/stirred_tank/factory.py
models/vlmodels/fermenter/config/kinetics.py     -> PyOMES/templates/stirred_tank/kinetics.py (see open question below)
models/vlmodels/fermenter/profiles.py            -> PyOMES/templates/stirred_tank/profiles.py
```

`profiles.py` moves along with the rest structurally — it's part of the
same package — but its separate open question (zero test coverage
anywhere, flagged in the demos/models pruning discussion) is untouched by
this move and should be resolved independently.

`kinetics.py` (Monod, Contois, Andrews, Tessier, Moser, Blackman,
DualSubstrateMonod) is a genuine open question, not a settled part of this
move — see "Open questions" below. Default placement above is inside
`stirred_tank/` (least disruptive); the alternative is merging it into the
already-existing `PyOMES/kinetics/` subpackage (which has its own
`KineticModel` protocol and a concrete model, `YeastAcetateV1`, tested by
`test_kinetics.py`) so there's one kinetics home in core instead of two
disconnected ones.

## Full blast radius (confirmed by repo-wide grep, not estimated)

Every file referencing `vlmodels.fermenter`, `FermenterBuilder`, or
`FermenterFactory`, split by risk:

**Mechanical (import path / class name swap only, no behavior change expected):**

- `tests/standalone/test_builder.py`, `test_configs.py` — stay exactly
  where they are; only their import lines change, no file move
- `tests/standalone/test_simulation.py` — the seven classes listed above
- `demos/builder/{batch,cstr,fed_batch,microplate}_fermenter.py`,
  `demos/model_api/export_results.py` — import line, plus each can drop
  its `import _bootstrap` line entirely (none of them touch `vlmodels` for
  anything else once this lands)
- `models/vlmodels/__init__.py` — drop the `fermenter` re-export entirely
  (the subpackage no longer exists in `vlmodels`) and the now-inaccurate
  "built on the fermenter framework" docstring sentence
- `README.md`, `demos/README.md`, `demos/builder/README.md`,
  `demos/_bootstrap.py`'s own docstring — narrative/index updates
- `PyOMES/core/recorder.py` — one stale docstring cross-reference, cosmetic

**Leave alone (frozen per this repo's own convention):**

- `docs/dev/implementation/shipped/EQUILIBRIUM_RESULT.md`,
  `SIMULATION_CLASS_CHECKLIST.md` — shipped docs are frozen except for
  factual corrections; a historical reference to the old name is not an
  error to fix.

**Needs real verification, not just find-and-replace:**

1. `models/vlmodels/adm1/base.py` (line 976) and `adm1/bsm2.py`
   (lines 804-805) **functionally depend** on this code — both lazily
   import `FermenterBuilder` (`bsm2.py` also imports `TransferConfig`)
   inside their `build_*_cv()` functions to construct their control
   volumes. This means the two flagship canonical models are built using
   this template internally, not just demoed alongside it.
2. `builder.py`'s `build_simulation()` has a documented circular-import
   workaround: a lazy `from PyOMES.core.simulation import Simulation`,
   with a comment explaining that a top-level import would collide with
   `PyOMES.core.lifecycle`'s own import chain. That workaround was written
   for a cross-package situation (`vlmodels` importing `PyOMES`). Moving
   the file bodily into `PyOMES` turns it into an intra-package import,
   which changes the shape of the problem — it likely still needs to stay
   lazy, but this needs to be checked in the new location, not assumed to
   carry over unchanged.

## Checkpoints

Matches this repo's `PHASE_KICKOFF_TEMPLATE.md` convention — one commit
per checkpoint, sanity check attached to each, owner commits per the
"Commit discipline" section above.

1. **Decide the `kinetics.py` destination** (open question below) —
   resolve before anything touches its importers.
2. **Create skeletons**: `PyOMES/templates/__init__.py` and
   `PyOMES/templates/stirred_tank/__init__.py`, both empty for now — real
   content drafted in checkpoint 4, once the re-export list is final.
   *Sanity check:* `python -c "import PyOMES.templates.stirred_tank"`.
3. **Move the four (or five, per checkpoint 1) implementation files**
   into `PyOMES/templates/stirred_tank/`, preserving history. Delete the
   two old `__init__.py` files outright (content is being rewritten, not
   moved). *Sanity check:* `models/vlmodels/fermenter/` no longer exists.
4. **Fix intra-package references inside the moved files**, file by file:
   - `builder.py`: `class FermenterBuilder` → `class StirredTankBuilder`;
     `from .factory import FermenterFactory` → `... import StirredTankFactory`;
     the internal call `FermenterFactory.create_volume(...)` (line 470);
     every `-> "FermenterBuilder"` return annotation; the `__repr__`
     literal; remove the now-resolved "Phase 7 holdover, rename deferred"
     docstring note (lines 459–466) instead of leaving it stale; add the
     gas+liquid-only constraint as an explicit docstring statement (per
     the naming decision above); re-examine the lazy `PyOMES.core.simulation`
     import per risk #2, and update its comment regardless since the
     "cross-package" framing no longer applies.
   - `factory.py`: `class FermenterFactory` → `class StirredTankFactory`;
     check for its own matching deferred-rename note (referenced from
     `builder.py`) and remove it.
   - `configs.py`: sweep for any stray `Fermenter`-branded docstring text
     (dataclass names themselves are unchanged, per the naming decision).
   - `kinetics.py`: fix the `>>> from PyOMES.config.kinetics import ...`
     docstring example.
   - `profiles.py`: fix the `>>> from PyOMES.profiles import ...`
     docstring example; strip the "Stage 17a" dev-stage label from the
     module docstring while it's open.
   - Write `stirred_tank/__init__.py`'s real content: re-export
     `StirredTankBuilder`, `StirredTankFactory`, the configs, kinetics
     classes, profile classes — mirroring the old `__all__` list, minus
     its stale "# Stage 17"/"# Stage 17c" comments.
   *Sanity check:* `python -c "from PyOMES.templates.stirred_tank import StirredTankBuilder"`
   succeeds on its own, before any external importer is touched.
5. **Fix the two functionally-dependent call sites**: `adm1/base.py:976`
   and `adm1/bsm2.py:804-805` (plus its `TransferConfig` import).
   *Sanity check:* run `pytest tests/standalone/test_bsm2_reference.py -v`
   **in isolation, right here** — not deferred to a final full-suite run —
   since this is the one check that would catch numerical drift from the
   move.
6. **Fix `vlmodels/__init__.py`**. *Sanity check:* `python -c "import vlmodels"`.
7. **Fix test imports**: `test_builder.py`, `test_configs.py` (import
   lines only, no move), and the seven named classes in `test_simulation.py`.
   *Sanity check:* `pytest tests/standalone/test_builder.py tests/standalone/test_configs.py tests/standalone/test_simulation.py -v`.
8. **Fix demo imports**: the four `demos/builder/*.py` files plus
   `export_results.py` — update the import line, drop `import _bootstrap`
   from all five. *Sanity check:* run each script directly and confirm it
   still prints its report.
9. **Fix remaining doc cross-references**: `PyOMES/core/recorder.py`'s
   stale docstring mention, `README.md`, `demos/README.md`,
   `demos/builder/README.md`.
10. **Full suite**: root `pytest`. Then ship per the branching/tagging
    convention (`--no-ff` merge, tag `stirred-tank-template-shipped`),
    moving this doc + its checklist to `../shipped/` with a "Shipped"
    banner.

## Open questions for whoever picks this up

1. **`kinetics.py` destination** — stay under the new `stirred_tank/`
   package, or merge into the existing `PyOMES/kinetics/` subpackage
   (which already has a `KineticModel` protocol and `YeastAcetateV1`)? The
   latter avoids two disconnected "kinetics" homes in core but is a larger
   change than a pure move.
2. **`hplc/column.py`'s fate** — flagged in the wider demos/models
   discussion as worth the same core-vs-example test this move just
   applied to the fermenter builder (it's generic textbook physics, not
   tied to one published dataset). Not resolved here; independent decision.
3. **`fermenter/profiles.py`'s test coverage** — moves structurally with
   everything else, but the "zero test coverage anywhere" finding from the
   earlier pruning discussion is untouched by this move and needs its own
   decision (test-and-keep vs. park).
4. **Future template ideas** (liquid-only stirred tank, precipitation-aware
   liquid tank, generic anaerobic-digestion stirred tank, a discretized 1D
   flow-reactor template generalizing `hplc/column.py`, a flash/VLE
   template off `PyOMES.equilibria`, a reactor-cascade/network template for
   true multi-zone topologies, and now potentially a three-phase
   `stirred_tank` variant per the naming discussion above) — named and
   deliberately deferred, not scoped for this phase. Recorded here only so
   they aren't re-discovered from scratch later.
5. **Relationship to the broader demos/`vlmodels` pruning discussion** —
   this migration was scoped out of that larger conversation as its own
   self-contained phase. The two are compatible (this move only makes the
   later pruning pass's core-vs-example sort cleaner) but are tracked
   separately; the wider pruning plan has not yet been written up as its
   own doc.
