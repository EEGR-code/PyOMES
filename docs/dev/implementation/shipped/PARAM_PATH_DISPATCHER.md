# Param-Path Dispatcher Rewrite — Design Note

> **Status:** Trigger-gated. Surfaced 2026-05-27 from the
> post-`simulation-class` audit. All five design questions
> (Q1–Q5) are pinned. Pick up when a second application domain
> needs to mutate parameters mid-run, when the closed-vocabulary
> brittleness of the current dispatcher causes a real bug, or
> bundled with another phase already touching the orchestrator.

## Context

[`Simulation._apply_param_change`](../../src/core/simulation.py)
([src/core/simulation.py:575-627](../../src/core/simulation.py#L575-L627))
is the orchestrator's Pattern B dispatch point — it takes a
`(path_string, value)` tuple from a
`ControlAction.params_changed` entry and routes the mutation to
the right object's unchecked setter. Today's implementation is
a hand-written if/elif on string prefixes:

```python
if path.startswith("kLa."):
    species = path[len("kLa."):]
    for iface in cv.internal_interfaces:
        if isinstance(iface, KineticGasLiquidLink):
            iface._set_kLa_unchecked(species, float(value))
            return
    return
if path.startswith("gas_feed.y."):
    ...
# Unknown path — silently ignored.
return
```

Three structural problems:

1. **Closed vocabulary.** Every new mutable parameter (HPLC
   column flow `Q`, dialysis `JV`, membrane `kLa_per_area`,
   ADM1 inhibition constants, …) requires editing this
   dispatcher. The orchestrator is the single point that
   knows about every domain's parameters.
2. **Hardcoded type lookups.** `isinstance(iface, KineticGasLiquidLink)`
   and `isinstance(boundary, GasFeed)` couple the orchestrator
   directly to fermenter classes.
3. **Silent failure.** A typo (`"kla.O2"` vs `"kLa.O2"`) or a
   path that was renamed silently mutates nothing.

For a framework intended to span multiple application domains
this is the largest leak in the orchestration layer. The
controller-input leak (`CVSnapshot`, fixed in
[`property-snapshot-phase-agnostic`](../shipped/SNAPSHOT_PHASE_AGNOSTIC.md)
2026-06-04) constrained *what controllers can read*; this leak
constrains *what controllers can write*.

## Proposed shape

Two layered surfaces, both backed by the same class-side
declarations:

### Surface 1 — descriptor-declared mutable attributes (Q4ʹ)

Each class with orchestrator-mutable attributes declares them
via a descriptor (the Pythonic industry-standard pattern,
mirroring SQLAlchemy / Pydantic / Django models):

```python
from VLsim.control.descriptors import MutableScalar, MutableDict

class KineticGasLiquidLink:
    kLa = MutableDict(species_key=str, value_type=float)
    # When mutated by the orchestrator, the descriptor delegates
    # to self._kLa internally and bypasses the lifecycle gate.
    # The public `link.kLa["O2"] = ...` path still runs through
    # raise_if_running.

class GasFeed:
    vvm_min = MutableScalar(float)
    y = MutableDict(species_key=str, value_type=float)

class GasPhase:
    T_K = MutableScalar(float, positive=True)
    V_L = MutableScalar(float, positive=True)
```

The descriptor encapsulates:
- The storage attribute (`self._kLa`, `self._T_K`, …).
- The validation (positive, range, type).
- The gated public-access path (consults `self._context.is_running`).
- The unchecked orchestrator-access path.
- The class-level metadata used by the path builder (Surface 2).

This replaces the current `@property` + `_set_<name>_unchecked`
pair on each gated attribute. Per-class declaration cost is one
line per mutable attribute; the descriptor machinery is in one
shared module
(`src/control/descriptors.py` or similar).

### Surface 2 — typed path builder via class-level access (Q5)

Class-level access to a descriptor yields a `ParamPath` object;
instance-level access yields the underlying value (SQLAlchemy
pattern: `User.name` at class level produces a column reference;
`user.name` at instance level produces the string value).

```python
from VLsim.core import KineticGasLiquidLink, GasFeed

# A controller emits:
return ControlAction(
    target_cv_key="main",
    params_changed={
        KineticGasLiquidLink.kLa["O2"]: 150.0,
        GasFeed.vvm_min: 0.5,
    },
)
```

`KineticGasLiquidLink.kLa["O2"]` returns a typed `ParamPath`
object the dispatcher consumes directly. Renaming `kLa` →
`mass_transfer_coef` breaks every callsite at edit time (IDE
rename), not at runtime.

A **string path form** is also accepted by the dispatcher for
ergonomics and serialisation (logging, replay, telemetry):

```
boundaries[GasFeed].vvm_min
internal_interfaces[KineticGasLiquidLink].kLa.O2
phases.liquid.T_K
```

The string form parses through the same `ParamPath` object the
typed form constructs directly, then both feed the walker
identically.

### The walker

At dispatch time the walker descends the CV graph guided by
the `ParamPath` (whether built typed or parsed from a string):

```python
def _apply_param_change(self, cv, path: ParamPath, value):
    target = cv
    for segment in path.segments[:-1]:
        target = segment.resolve(target)
    leaf = path.segments[-1]
    # The descriptor knows how to do the unchecked write itself:
    if isinstance(leaf, DescriptorSegment):
        leaf.descriptor.__set_unchecked__(target, leaf.key, value)
    elif isinstance(target, dict):
        target[leaf.name] = value
    else:
        raise KeyError(
            f"params_changed path {path!r} does not resolve to a "
            f"writable target on cv {cv.label!r}."
        )
```

Estimated implementation: ~150 LOC for the descriptors module,
~80 LOC for the walker + `ParamPath`, ~40 LOC for the string
parser, plus migration of every existing `@property` +
`_set_X_unchecked` pair to a descriptor declaration (≈5 classes
× 2 attrs each = 10 sites). ~50 tests.

## Open design questions — all pinned

### Q1 — Path syntax: shorthand vs full reflective form

**Decision:** *Full reflective form only.* No alias table.

The string form is the full reflective path
(`internal_interfaces[KineticGasLiquidLink].kLa.O2`); the typed
form is the equivalent Python expression
(`KineticGasLiquidLink.kLa["O2"]`). There is no shorthand alias
table.

Rationale: aliases re-introduce the closed-vocabulary problem on
a smaller surface, and the full form is self-documenting — a
reader can resolve the path by reading the CV's structure
without consulting an alias glossary.

Migration impact: every existing emit in `cv_loops.py` that
writes `"kLa.O2"`, `"gas_feed.y.O2"`, `"gas_feed.vvm_min"`
becomes the full string form (or the typed form). ~15 callsites;
mechanical change.

### Q2 — List selector default

**Decision:** *Type selector by default; type-plus-index for
disambiguation; pure index as escape hatch.*

Three syntaxes layered (string form shown; typed form uses
Python attribute / item syntax to the same effect):

- `[CamelCaseIdent]` — pick the first element where
  `isinstance(elem, CamelCaseIdent)` is `True`. The common
  case.
- `[CamelCaseIdent:N]` — pick the Nth match (zero-indexed)
  among the type filter. Used when a CV has multiple internal
  interfaces or boundaries of the same type.
- `[N]` (pure digits) — pick list index N regardless of type.
  Escape hatch reserved for tests and exotic ordering-sensitive
  cases; normal controller code prefers the type form.

The walker disambiguates by inspecting the bracket contents:
digits → pure index; identifier (CamelCase, or with a colon)
→ type selector.

Class resolution: the walker needs a registry from the bracket
identifier string to the actual Python class. v1 implementation
scans `VLsim.core.__init__.__all__` at module import (closed to
plugins, simple). A `@register_param_class` decorator widening
this to user-defined classes is a follow-on if a plugin
ecosystem emerges.

### Q3 — Unknown-path behaviour

**Decision:** *Error.*

Unknown paths raise `ParamPathError` (a dedicated exception
subclass) immediately, with the path and the CV label in the
message. No silent ignore; no warn-once.

Rationale: this aligns with how unknown
`ControlAction.target_cv_key` already raises, and matches the
"fail-loud configuration" principle from the audit. Controllers
that genuinely want optional paths (e.g. only emit `kLa.CO2`
when CO₂ is in the gas mix) wrap the emit in a `try/except` or
check the path against a known-paths set themselves.

### Q4ʹ — Mutable attribute declaration: introspective or
descriptor-based?

**Decision:** *Descriptors.*

Each class with orchestrator-mutable attributes declares them
via a `MutableScalar` / `MutableDict` / `MutableSequence`
descriptor. The descriptor IS the source of truth for what's
mutable; the dispatcher reads class-level metadata from it.

Considered alternatives:

- **(rejected) Introspective discovery** — walker probes for
  `_set_<name>_unchecked` methods by name. Conventional in
  Ruby and ActiveRecord but unusual in Python frameworks;
  fragile to typos in the underscore method name, no
  type-checker support, no obvious way to surface the "what's
  mutable?" question to a user.
- **(rejected) Class-side `MUTABLE_PATHS` dict** — explicit but
  duplicates information: the dict says "kLa is mutable",
  meanwhile the public `@property` setter and the
  `_set_kLa_unchecked` method both exist independently.
  Three sites to keep in sync.
- **Descriptors** — single source of truth; the descriptor
  defines storage, validation, gated public access, and
  unchecked orchestrator access in one declaration.
  Industry-standard pattern in Python (SQLAlchemy, Pydantic,
  Django ORM).

Migration cost: replace each existing `@property` +
`@x.setter` + `_set_<name>_unchecked` triplet (10 sites across
~5 classes) with a single descriptor declaration. The
descriptor module is new (~150 LOC).

### Q5 — Typed `ParamPath` builder

**Decision:** *Typed builder, driven by the same descriptors as
Q4ʹ.*

Class-level access to a descriptor yields a `ParamPath`;
instance-level access yields the underlying value. This is the
SQLAlchemy pattern (`User.name` at class level produces a
column reference; `user.name` at instance level produces the
string value).

```python
# At class level — produces a typed path:
KineticGasLiquidLink.kLa["O2"]    # -> ParamPath
GasFeed.vvm_min                    # -> ParamPath

# At instance level — produces the value:
link.kLa["O2"]                     # -> 150.0 (a float)
gas_feed.vvm_min                   # -> 0.5 (a float)
```

A **string path form** is also accepted by the dispatcher for
ergonomics, serialisation, and any callsite that doesn't have
the class in scope:

```
boundaries[GasFeed].vvm_min
internal_interfaces[KineticGasLiquidLink].kLa.O2
phases.liquid.T_K
```

Both forms parse to the same `ParamPath` representation, then
walk identically.

Rationale: typed paths get refactor safety (renaming `kLa`
breaks the callsites at edit time) and IDE autocomplete; the
string form keeps replay / logging / telemetry simple. The
duplication cost is bounded by the fact that both forms feed
the same walker — no two separate dispatchers to maintain.

Revisit conditions (these would suggest the typed form is too
expensive to maintain):
1. Plugin ecosystem outside the tree wants to write
   controllers — typed selectors require the plugin's classes
   to be importable, which may not always be true.
2. Round-tripping paths through JSON / YAML / database becomes
   a frequent operation — string form dominates.

If either fires, demote the typed form to a convenience layer
on top of the canonical string form.

## Out of scope

- **Cross-CV param paths.** Today every `params_changed` entry
  targets the action's `target_cv_key` only. Mutating
  parameters on a *different* CV from the same controller
  would need a `(cv_key, path)` tuple form, not covered by
  this note.
- **Mutating `Simulation`-level state from a controller**
  (e.g. swapping the recorder mid-run). Not required by any
  current use case; not designed for.
- **Type coercion of values.** The descriptor enforces type
  on write (e.g. `MutableScalar(float, positive=True)` rejects
  negative floats); but auto-coercion from `int` → `float`
  or from `numpy.float64` → `float` is an implementation
  detail.

## Trigger conditions

Any one of:

1. **A second application domain needs param-path mutation.**
   HPLC adding a controller that adjusts mobile-phase flow,
   wastewater adding a feed-rate controller, etc. — anything
   that today's dispatcher cannot route without an edit.
2. **A real silent-typo bug in `cv_loops.py`.** The closed
   vocabulary is so small today that this is unlikely, but if
   it happens once, the fail-loud rewrite immediately pays
   for itself.
3. **Bundled with another phase that's already touching the
   orchestrator.** If `RUN_HISTORY` or `CONTAINER_LAYERING`
   lands and is already editing `Simulation._step`, folding
   this in costs less than a separate phase.

## Relationship to other phases

- **Completes the audit trilogy with**
  [`../shipped/SNAPSHOT_PHASE_AGNOSTIC.md`](../shipped/SNAPSHOT_PHASE_AGNOSTIC.md)
  and
  [`../shipped/PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md`](../shipped/PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md)
  — all three addressed the same audit finding ("fermenter-shaped
  framework surfaces dressed up as general"). Those two shipped
  2026-06-04; this note is the last of the three.
- **Interacts with
  [`RUN_HISTORY.md`](RUN_HISTORY.md).** The string form of
  `ParamPath` is the natural serialised shape for the
  `params_changed` channel of a richer recorder. If
  RUN_HISTORY ships first (still recording the old
  closed-vocabulary paths), the dispatcher rewrite later has
  to migrate the recorded format. If this note ships first,
  RUN_HISTORY's recorder picks up the typed/string form for
  free.
- **Subsumes** the
  [framework-polish bundle, P8](FRAMEWORK_POLISH.md#p8-pattern-b-completion-on-boundary-classes):
  the Pattern B completion on boundary classes (adding
  `_set_y_unchecked` etc. to `GasFeed`) is no longer a
  prerequisite because the descriptor model replaces the
  property + unchecked-setter pair entirely. If this rewrite
  goes first, P8 dissolves; if P8 ships first as part of the
  polish bundle, the descriptor migration consumes those
  unchecked setters during the migration step.
- **Independent of** all other trigger-gated phases.
- **Follows** `simulation-class` (shipped 2026-05-27), which
  introduced `_apply_param_change` in its current form at C9.
