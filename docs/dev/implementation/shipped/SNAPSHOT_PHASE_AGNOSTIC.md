# Snapshot / Recorder Phase-Agnostic Surfaces — Design Note

> **Status:** Trigger-gated. Surfaced 2026-05-27 from the
> post-`simulation-class` audit as the second of the two
> generalisation gaps where the liquid/gas-phase assumption
> leaked from fermenter use cases into supposedly general
> framework machinery (the first is
> [PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md](PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md)).
> Pick up when a real non-fermenter CV with phases named other
> than `"liquid"` / `"gas"` is built, or pre-emptively bundled
> with the property-calculator note.

## Context

Two framework surfaces — `CVSnapshot` (the controller's read-only
view of a CV) and `BatchRecorder` (the default time-series
recorder) — both hardcode the assumption that every CV has
exactly the phase keys `"liquid"` and `"gas"`, and they bake in
some fermenter-specific sensor derivations on top of that.

### `CVSnapshot` — controller-input surface

[`src/core/snapshot.py:170-244`](../../src/core/snapshot.py#L170-L244)
in `build_cv_snapshot`:

```python
liquid = cv.phases.get("liquid")
gas = cv.phases.get("gas")
...
sensors: Dict[str, Any] = {"T_C": T_K - 273.15}
if V_liq_L > 0.0 and "O2" in n_liq_mol:
    sensors["DO_mol_L"] = n_liq_mol["O2"] / V_liq_L
```

A CV with phases named `"mobile"` / `"stationary"` (HPLC),
`"permeate"` / `"retentate"` (membrane), or
`"intracellular"` / `"extracellular"` (cell culture) yields a
snapshot whose `V_liq_L`, `V_gas_L`, `n_liq_mol`, `n_gas_mol`,
`pH`, `ionic_strength`, `species_properties` fields are all
zero / empty / `None`. The controller has nothing to read; no
controller for that domain can be written without bypassing
`CVSnapshot`.

The `sensors["DO_mol_L"]` derivation is doubly fermenter-shaped:
it assumes the liquid phase tracks `"O2"` and that the user
cares about dissolved oxygen.

### `BatchRecorder` — time-series storage surface

[`src/core/recorder.py:200-216`](../../src/core/recorder.py#L200-L216)
in `record_init`:

```python
gas = cv.phases.get("gas")
liquid = cv.phases.get("liquid")
if gas is not None:
    for sp in gas.n_mol:
        self._result.gas_mol[cv_key][sp] = np.zeros(n)
if liquid is not None:
    for sp in liquid.n_mol:
        self._result.liquid_mol[cv_key][sp] = np.zeros(n)
```

And in `_record_cv_state` ([line 289-329](../../src/core/recorder.py#L289-L329))
the per-step writeback hits only those same two keys. Solid
phases (`SolidPhase` exists in `phases.py` but is never
recorded), arbitrary phase names, and multi-liquid systems are
silently dropped from the recorded trajectory.

`BatchResult` itself reflects the same shape: `gas_mol`,
`liquid_mol`, `pH`, `ionic_strength`, `P_atm` keys with no
neutral "phase_mol" dict, so even if the recorder were fixed,
downstream consumers would still need to learn new field names.

## Why this matters for the framework

Controllers and the recorder are the two surfaces a downstream
user touches first when building a model. A user whose domain
doesn't have `"liquid"` and `"gas"` phases discovers within
minutes that:

- The controller protocol is unusable — `CVSnapshot` is empty.
- The default recorder loses all per-species trajectory.

Both can be worked around (write your own snapshot builder;
write your own recorder), but the workaround means rebuilding
the whole controller-input contract and the whole time-series
schema rather than configuring a small extension point. The
generalisation cost is large because the assumption is baked
into the canonical types, not into a runner.

For cross-domain reuse — wastewater, HPLC, cell culture,
membrane separations — this is the second-largest
fermenter-shaped leak in the framework (after the param-path
dispatcher; see audit follow-up).

## Proposed shape (when picked up)

The minimal change is to make the snapshot and recorder iterate
over every phase the CV exposes, while retaining the
`liquid` / `gas` shorthand for the common case.

### `CVSnapshot` — generalised

```python
@dataclass(frozen=True)
class CVSnapshot:
    cv_key: str
    t_h: float
    T_K: float

    # General per-phase storage
    phase_keys: List[str]
    n_mol: Dict[str, Dict[str, float]]     # {phase_key: {species: mol}}
    V_L: Dict[str, float]                  # {phase_key: V_L}
    properties: Dict[str, Dict[str, float]]  # {phase_key: properties dict}

    # Convenience aliases for the common fermenter shape — None when
    # the phase isn't present. Computed in build_cv_snapshot from the
    # general fields above; controllers MAY use these instead of
    # walking phase_keys themselves.
    pH: Optional[float] = None             # alias for phase("liquid").pH
    ionic_strength: Optional[float] = None
    V_liq_L: Optional[float] = None
    V_gas_L: Optional[float] = None
    P_gas_atm: Optional[float] = None
    n_liq_mol: Optional[Dict[str, float]] = None
    n_gas_mol: Optional[Dict[str, float]] = None
    y_gas: Optional[Dict[str, float]] = None

    sensors: Dict[str, Any] = field(default_factory=dict)
```

The convenience aliases keep every existing fermenter controller
working without an edit. New domains write controllers against
the general `n_mol[phase_key]` / `V_L[phase_key]` /
`properties[phase_key]` shape.

### `sensors` derivation — declarative extension

Currently the sensors dict is populated by a small hand-written
block at the end of `build_cv_snapshot`. Replace with a list of
registered sensor functions on `Simulation` (or per-CV), each of
which takes the CV and returns a (key, value) pair. Default
registrations cover `T_C` and `DO_mol_L` for the fermenter
shape; users register their own (e.g. `column_pressure_drop` for
HPLC) without editing the snapshot builder.

### `BatchRecorder` — generalised storage

```python
@dataclass
class BatchResult:
    t_h: np.ndarray
    # General per-phase storage: {cv_key: {phase_key: {species: ndarray}}}
    phase_mol: Dict[str, Dict[str, Dict[str, np.ndarray]]]
    # Per-phase scalar property time-series:
    # {cv_key: {phase_key: {property_key: ndarray}}}
    phase_properties: Dict[str, Dict[str, Dict[str, np.ndarray]]]
    # Special-cased aliases (kept for fermenter compatibility):
    pH: Dict[str, np.ndarray]
    ionic_strength: Dict[str, np.ndarray]
    P_atm: Dict[str, np.ndarray]
    gas_mol: Dict[str, Dict[str, np.ndarray]]      # alias for phase_mol[cv]["gas"]
    liquid_mol: Dict[str, Dict[str, np.ndarray]]   # alias for phase_mol[cv]["liquid"]
    # … link/boundary/controller/profile channels unchanged …
```

`gas_mol` and `liquid_mol` become views into the general
`phase_mol` storage rather than primary fields. Existing
analysis code keeps working.

## Open design questions to resolve when picked up

1. **Convenience aliases — `Optional` or `0.0`/empty default?**
   Today's API returns `0.0` for `V_liq_L` when no liquid phase
   is present, which is wrong (should be "absent" not "zero").
   Changing to `Optional[float]` is the honest fix but a
   breaking change for the few controllers that assume
   non-`None`. Probably worth the break.

2. **Sensors registration — global, per-Simulation, or per-CV?**
   Per-CV is the most flexible; per-Simulation is simplest;
   global is brittle. Per-CV with a Simulation-level default
   list is a reasonable middle.

3. **`BatchResult.phase_mol` as primary storage** — should the
   `gas_mol` / `liquid_mol` aliases be properties (computed on
   access) or eagerly-mirrored dicts (consume memory but match
   today's debugger view)? Properties is the right answer;
   eager mirroring duplicates ~half the result's storage.

4. **`SolidPhase` recording** — `SolidPhase` exists and has
   `n_mol`, `total_mol()`, `apply_flux()`. The current recorder
   silently ignores it. Once `phase_mol` is the primary
   storage, solids fall in naturally; worth verifying nothing
   else (pressure derivation, conservation monitor) treats
   solids as second-class.

5. **`y_gas` and gas `P_atm`** — these are gas-phase-derived
   convenience quantities. Should they migrate into the general
   `phase_properties` channel (every phase can publish whatever
   derived scalars it likes via `PropertyCalculator`) or stay
   as named aliases on the snapshot?
   **Hard dependency:** this question cannot be answered without
   [`PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md`](PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md)
   shipping a gas-phase property-calculator runner first. The
   migration path it gestures at (publish via `PropertyCalculator`)
   does not exist today.

## Out of scope

- Pluggable / streaming recorders (`SparseRecorder`,
  `StreamingRecorder`, `SummaryRecorder`) — those are the
  `RUN_HISTORY` phase's territory. This note is only about the
  per-phase-key generalisation of the existing in-memory
  `BatchRecorder` and `CVSnapshot` shapes.
- The fermenter-specific `sensors["DO_mol_L"]` value semantics
  themselves (whether "dissolved oxygen as a concentration"
  belongs in a generic snapshot at all) — accept as legacy
  fermenter convenience; just don't hardcode the derivation
  inside the builder.
- Snapshotting `internal_interfaces` state (link kLa, boundary
  setpoints) — not currently in `CVSnapshot`; out of scope
  unless a controller wants to read its own previous-step
  state.

## Trigger conditions

Any one of:

1. A real non-fermenter CV (HPLC mobile/stationary, membrane
   permeate/retentate, cell-culture intracellular/extracellular)
   needs to be controlled by a `Controller` consuming a
   `CVSnapshot`.
2. A solid-phase trajectory needs to be recorded for analysis.
3. The "framework polish" bundle picks this up pre-emptively
   because the audit flagged it and the cost is low (~120 LOC
   + tests; no `BatchResult` consumer that we own will break
   if the aliases-as-properties pattern is followed).

## Relationship to other phases

- **HPLC trigger bundle (added 2026-05-28).** HPLC
  (mobile/stationary phases) triggers this note, the
  per-phase property runner
  ([`PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md`](PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md)),
  and the unit-physics descriptor work in
  [`CONTAINER_LAYERING.md`](CONTAINER_LAYERING.md) all at once.
  Plan the HPLC pickup as a bundle, not as three independent
  phases.
- **Pairs naturally with**
  [`PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md`](PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md)
  — both are runner-side / type-side fixes for the same leak,
  share testing infrastructure (CVs with non-standard phase
  keys), and would be confusing to ship in different orders.
  Q5 above is a **hard dependency** on that note's gas-phase
  runner, not merely a soft pairing.
- **Independent of**
  [`RUN_HISTORY.md`](RUN_HISTORY.md) — that phase reframes
  the result schema for richer recording (streaming, sparse,
  controller-aware). This phase is the prerequisite cleanup:
  phase-key-agnostic storage in the existing recorder shape.
  RUN_HISTORY can land on top.
- **Independent of**
  [`CONTAINER_LAYERING.md`](CONTAINER_LAYERING.md),
  [`CUFERMENTER_SUNSET.md`](CUFERMENTER_SUNSET.md),
  [`CHEMISTRY_UNIFICATION_PLAN.md`](CHEMISTRY_UNIFICATION_PLAN.md)
  Phase 3b, and `SPECIATION_LEVEL_RETIREMENT.md`.
- **Follows** `simulation-class` (shipped 2026-05-27), which
  introduced both `CVSnapshot` and `BatchRecorder` in their
  current shape.
