# RESERVOIR_TYPE — FlowBoundary unification and the Reservoir type

> **Status:** Design discussion, 2026-07-10, revised 2026-07-10. Not yet
> started — no branch, no checklist, no code written.
>
> **Revision note.** This note's original conclusion (§2 below, preserved
> for the record) kept `PhaseInterface`/`CVLink`/`ExternalBoundary` as three
> separate protocols — the three-way merge looked structurally blocked.
> Follow-up discussion revisited that conclusion once `Reservoir` was fully
> specified: **§5 now adopts `FlowBoundary` as a real unifying protocol.**
> `Reservoir` turned out to be the piece that makes it work, not just a
> fix for `ExternalBoundary` in isolation — §5 explains why each of §2's
> three obstacles is actually addressed, not just newly ignored.
>
> **Depends on nothing shipped.** Extends the open unification question
> already on record in
> [`MASS_EXCHANGE_ARCHITECTURE.md`](../design/MASS_EXCHANGE_ARCHITECTURE.md)
> §7 and §12 Q1.
>
> **Explicitly not resolved here:** how `PartitionModel` relates to
> `FlowBoundary` (and to `EquilibriumPhenomena` — see
> [`PHENOMENA_PROTOCOL.md`](PHENOMENA_PROTOCOL.md)). Under active discussion
> as a separate thread; §5 deliberately treats `PartitionModel` as an
> unchanged implementation detail wherever it appears below.

---

## 1. Background — the three protocols as they exist today

`VLsim` has three transport-domain protocols, all returning a
`{species_id: flux_mol_per_h}` dict, differing only in the topology of their
two endpoints (`MASS_EXCHANGE_ARCHITECTURE.md` §2's Axis 2):

| Protocol | Endpoints | Both sides tracked? | Method | Wired at |
|---|---|---|---|---|
| [`PhaseInterface`](../../src/core/interfaces.py) | two phases, same CV | Yes | `compute_flux(state_a: Phase, state_b: Phase, dt_h)` | inside `cv.advance()`, step 4 |
| [`CVLink`](../../src/core/links.py) | two phases, different CVs | Yes | `compute_flow(cvs: Dict[str, ControlVolume], dt_h)` | `Simulation._apply_links()`, before any CV advances |
| [`ExternalBoundary`](../../src/core/boundaries.py) | one phase, CV ↔ outside the model | **No** — only one side exists | `compute_flux(cv: ControlVolume, dt_h)` | inside `cv.advance()`, step 2b |

`PhaseInterface` and `CVLink` are symmetric and mass-conserving — a source
loses exactly what a sink gains, and this is checked
(`TransferDiagnostics.is_conserved()`, `LinkFlowRecord`). `ExternalBoundary`
is explicitly *not* — its own module docstring states the distinction
outright:

> "The distinction is physically meaningful: internal transfers conserve the
> CV's total inventory, while external fluxes intentionally change it."
> ([`boundaries.py:10-11`](../../src/core/boundaries.py))

The question that opened this discussion: does it make sense to unify all
three under one conceptual type (`FlowBoundary` or similar)?

## 2. Why a naive three-way merge looked blocked (original conclusion)

`MASS_EXCHANGE_ARCHITECTURE.md` §7 already sketches a two-way version of
this (`PhaseInterface` + `CVLink` only, "same CV" as a degenerate case of one
`TransportLink` type taking `(source_cv, source_phase)` /
`(sink_cv, sink_phase)`), and leaves it explicitly open in §12 Q1 — *"Or is
the execution-context difference ... a sufficient reason to keep them
separate?"* Three concrete obstacles surfaced in the first pass of
discussion:

1. **Signature mismatch.** `compute_flux(state_a, state_b, dt_h)` takes two
   `Phase` objects; `compute_flow(cvs, dt_h)` takes the whole CV registry;
   `compute_flux(cv, dt_h)` takes one whole `ControlVolume` (several
   boundary implementations, e.g. a pressure-relief vent, need CV-wide
   context — headspace pressure summed across every gas species — not just
   a single phase). *Resolved — see §5.1.*
2. **Execution-context/timing differs, and a shared base type doesn't
   collapse it.** `CVLink` is invoked at the `Simulation` level, before any
   CV advances, because its two endpoints may not share a CV.
   `PhaseInterface` and `ExternalBoundary` are both invoked inside
   `cv.advance()`, but at different steps. *Partially resolved — see §5.2:
   the calling convention unifies; the timing question does not, and isn't
   claimed to.*
3. **`ExternalBoundary` has no second endpoint to name.** Its asymmetry
   isn't an implementation gap — the design intends for mass to
   representationally leave/enter the modeled system at a boundary. Forcing
   it into a `(source, sink)` shape needs *something* to put on the other
   side. *Resolved — see §5.3, this is exactly what `Reservoir` (§3-4) is
   for.*

## 3. The `Reservoir` proposal

Introduce a `Reservoir` type to serve as `ExternalBoundary`'s "other side,"
refined through discussion into a specific, narrow shape:

### 3.1 Passive accumulator, not a depleting state

`Reservoir` does **not** need to model depletion-driven physics to be worth
having. It holds a running total —
`reservoir.cumulative_mol[species] += -flux` — updated at the exact point a
boundary's flux is already applied
([`solvers.py:171-180`](../../src/core/solvers.py), step 2b). This is
*write-only*: nothing reads the accumulator back into `compute_flux`'s
calculation. The boundary still computes its flux exactly as before (fixed
composition for a feed, pressure-driven for a vent); the reservoir only
watches and totals.

This mirrors an existing pattern in the codebase — `AccuracyMonitor` /
`ConservationMonitor` are already attached to `ControlVolume` as passive
observers that never feed back into the physics loop. A `Reservoir`
accumulator is the same pattern applied to the boundary side.

**What this buys:** a genuine capability gap closes. Per
[`control_volume.py`](../../src/core/control_volume.py)'s own
`apply_external_flux` docstring, external fluxes are "not tracked by
`step_internal_transfer` diagnostics (because external fluxes intentionally
change the CV's total inventory)" — currently nothing accumulates the total
drained/fed past a CV boundary for a whole-model mass-balance audit. A
passive `Reservoir` accumulator closes that.

**What this avoids:** deciding "what happens when a reservoir runs dry."
There is no depletion feedback, so there's no clamping/negative-value
question to resolve for the common case (atmosphere, waste line) — the
accumulator just keeps counting, which is the physically correct picture
for an effectively-infinite external world.

### 3.2 Multi-phase accounting borrowed from `Phase`, not from `ControlVolume`

`Reservoir` should reuse `Phase`'s per-phase primitives
(`phases: Dict[str, Phase]`, `apply_flux`) to get real multi-species /
multi-phase bookkeeping — e.g. one shared `Atmosphere` reservoir backing
both a `GasFeed` (drawing from it) and a vent (returning to it), giving
genuine net-exchange accounting instead of two uncorrelated boundary
objects.

It should **not** adopt `ControlVolume`'s full machinery —
`reaction_system`/`PhenomenaSystem` (see
[`PHENOMENA_PROTOCOL.md`](PHENOMENA_PROTOCOL.md) for the rename),
`property_calculators`, `AccuracyMonitor`/`ConservationMonitor` attachment,
its own `boundaries` list, `ParamPath`/lifecycle-locking. All of that is
irrelevant to something that is a boundary condition by definition, and
duplicating it would make `Reservoir` a second name for
`ControlVolume`'s existing responsibility — the same "don't build two
objects that do structurally the same job" lesson already recorded in
`MASS_EXCHANGE_ARCHITECTURE.md` §14 for the Ksp/Henry folding decision.

### 3.3 The honest remaining judgement call

Once scoped this narrowly, `Reservoir` is functionally close to
`ControlVolume(phases={...}, reaction_system=None)` — which already works
today with zero new code, by using the existing `CVLink` machinery to
connect a full CV to a minimal one. The case for a dedicated `Reservoir`
type was never about unlocking otherwise-impossible capability on its own:

- **Type-level intent.** A `Reservoir` structurally guarantees "this can
  never react, it's always a boundary" — the same kind of decision
  `STATE_UNIFICATION` already made splitting `Reaction` into
  `KineticReaction`/`EquilibriumReaction` specifically so "runtime errors
  became type errors," and the same discipline `EquilibriumReaction` vs.
  `HenryEquilibrium`/`KspEquilibrium`/`RaoultEquilibrium` already
  demonstrates (siblings under `EquilibriumConstraint`, not
  isinstance-checked against each other — see
  [`equilibrium.py:115-119`](../../src/reactions/equilibrium.py)).
- **Avoided setup cost.** `ControlVolume.__init__` always constructs
  monitors, lifecycle locks, and caches even when unused; a lean
  `Reservoir` skips all of it.

**This is now additionally justified by §5** — `Reservoir` isn't only a
clarity/ergonomics call anymore; it's the specific piece that makes the
`FlowBoundary` unification possible at all (§5.3). The case for building it
is stronger post-revision than it was when this section was first written.

### 3.4 Which `ExternalBoundary` implementers actually get a `Reservoir`

Not yet audited. Working hypothesis: sort each of the six existing
implementers (`LiquidFeed`, `GasFeed`, `LiquidDrain`, `MembraneGasBoundary`,
`PressureReliefVent`, `ProportionalGasOutlet`) by "would we ever want to
audit what's on the other side" — feed/drain candidates likely yes,
vent/relief candidates likely no beyond the passive accumulator. Every
`ExternalBoundary` still needs *some* `PhaseCarrier` sink once `FlowBoundary`
exists (§5.1) — the open question is only whether that sink is a real,
named, shared `Reservoir` or an anonymous throwaway one per boundary
instance.

---

## 4. What `Reservoir` alone would have unified (superseded framing)

> The original version of this note stopped here, concluding that
> `Reservoir` fixes `ExternalBoundary`'s asymmetry but that a full
> three-way protocol merge remained separately blocked by §2's obstacles 1
> and 2. That framing is superseded by §5 — obstacle 1 turned out to
> resolve *because* of the same `phases: Dict[str, Phase]` shape noted
> here. Kept for the record since it's the bridge that led to §5, not
> because the conclusion still stands as originally stated.

If `Reservoir` exposes the same `phases: Dict[str, Phase]` shape
`ControlVolume` does, the "sink" side of a transport link doesn't need to
be `ControlVolume` *or* `Reservoir` specifically — it can be typed
structurally as "anything with `.phases[key]`," satisfied by both.

---

## 5. Resolution: `FlowBoundary`, enabled by `Reservoir`

Revisiting §2's three obstacles now that `Reservoir` exists as a designed
(if not yet built) type.

### 5.1 A structural `PhaseCarrier` closes the signature gap

§2's obstacle 1 was three different first-argument shapes
(`Phase`+`Phase`, `Dict[str, ControlVolume]`, `ControlVolume`). A thin
structural protocol, satisfied by both `ControlVolume` and `Reservoir`
without either needing to change, resolves it:

```python
@runtime_checkable
class PhaseCarrier(Protocol):
    phases: Dict[str, Phase]
```

```python
@runtime_checkable
class FlowBoundary(Protocol):
    def compute_flow(
        self,
        source: PhaseCarrier, source_phase_key: str,
        sink: PhaseCarrier, sink_phase_key: str,
        dt_h: float,
    ) -> Dict[str, float]: ...
```

One method shape now covers all three topologies:

- **`PhaseInterface`'s case** — `source is sink` (the same `ControlVolume`),
  different `phase_key`s. Exactly the "same CV as a degenerate case"
  framing `MASS_EXCHANGE_ARCHITECTURE.md` §7 already proposed.
- **`CVLink`'s case** — `source`/`sink` are different `ControlVolume`
  instances.
- **`ExternalBoundary`'s case** — `source` is a `ControlVolume`, `sink` is
  a `Reservoir`. This is the piece that was missing before `Reservoir`
  existed — there was no second `PhaseCarrier` to put there.

Implementations needing CV-wide context beyond a single phase (e.g.
`PressureReliefVent` reading headspace pressure across every gas species)
read it off `source.phases[...]` directly — `PhaseCarrier` exposes the
whole `phases` dict, not just the one keyed phase, so this doesn't need a
signature exception.

### 5.2 Execution-context/timing: the calling convention unifies; *when* it fires does not

§2's obstacle 2 said a shared type doesn't collapse the dispatch problem.
That remains half true. **What changes:** every `FlowBoundary` instance is
now called the same way regardless of topology — today's three different
calling conventions (`iface.compute_flux(phase_a, phase_b, dt_h)` /
`link.compute_flow(cvs, dt_h)` / `boundary.compute_flux(cv, dt_h)`) collapse
to one (`fb.compute_flow(source, source_phase_key, sink, sink_phase_key,
dt_h)`). **What does not change:** *when* in the timestep each instance
needs to fire. A same-`ControlVolume` `FlowBoundary` still needs to run
inside `cv.advance()` (same-step chemistry coupling, per everything
established about `advance()`'s sequential body); a different-CV
`FlowBoundary` still needs to run at the `Simulation` level, before any CV
advances. The orchestrator still sorts instances by topology to know which
loop calls them — but now it does so by comparing `source`/`sink` identity,
not by `isinstance` against three different protocols.

### 5.3 `Reservoir` makes `ExternalBoundary`'s case genuinely symmetric

§2's obstacle 3 — no second endpoint to name — is exactly what `Reservoir`
was designed to fix (§3). With it, `ExternalBoundary`'s case is symmetric
and mass-conserving the same way `PhaseInterface`/`CVLink` already are: the
`ControlVolume` side loses via `apply_flux`, the `Reservoir` side gains via
its passive accumulator. This is the load-bearing piece — §5.1 and §5.2
only work because this one does. `Reservoir`'s design was originally
motivated purely by audit-trail accounting (§3.1); unlocking `FlowBoundary`
turned out to be a second, larger consequence of the same design.

### 5.4 Ownership does not change — two arrows, one type

`FlowBoundary` unifies the *protocol*, not *ownership*. `ControlVolume`
still owns same-`ControlVolume` and `ControlVolume`-to-`Reservoir`
instances — today's separate `internal_interfaces` and `boundaries` lists
could merge into one `flow_boundaries: List[FlowBoundary]`, or stay as two
lists of the same type for semantic clarity (not decided, see Open
Questions). `Simulation` still owns cross-CV instances (`sim.links`,
unchanged in ownership, just typed `List[FlowBoundary]` instead of
`List[CVLink]`).

### 5.5 What is deliberately left open by this resolution

`PartitionModel`'s relationship to `FlowBoundary` is not addressed here —
under active discussion as its own thread (see the status banner). Nothing
in §5.1-5.4 depends on how that question resolves; `PartitionModel` is
treated as an unchanged implementation detail wherever a `FlowBoundary`
implementer happens to wrap one, same as today.

---

## Open questions

1. **Naming.** `FlowBoundary` was the working name — but
   `MASS_EXCHANGE_ARCHITECTURE.md` §4 already proposed renaming `CVLink` to
   `InterzonalFlow` (never shipped). Does adopting `FlowBoundary` retire
   that dangling rename proposal, replace it with a different name
   entirely, or do the two need reconciling some other way? `Reservoir`
   itself also wasn't checked against existing identifiers.
2. **`ControlVolume`'s attribute shape.** Does `internal_interfaces` +
   `boundaries` merge into one `flow_boundaries: List[FlowBoundary]`, or
   stay two lists of the same type for readability (§5.4)? Merging removes
   a distinction that's currently informative at a glance (same-CV vs.
   external) even though the type no longer requires it.
3. **`PartitionModel`'s relationship to `FlowBoundary`** — explicitly
   deferred (§5.5), to be resolved in a separate discussion.
4. **Attachment shape for `Reservoir`.** Is a `Reservoir` optional
   per-`ExternalBoundary` instance, or can multiple boundaries share one
   (the `Atmosphere` / `GasFeed`-plus-vent case explicitly wants sharing,
   §3.4)? If shared, where does the shared instance live — owned by the
   `ControlVolume`, by the `Simulation`, or constructed standalone and
   passed to multiple boundaries?
5. **Energy accounting.** The original proposal included "material and/or
   energy." Whether `VLsim` has a general enthalpy/heat-duty balance layer
   for a `Reservoir` to plug into has **not been verified** — `T_K` is
   tracked on `Phase` and `ThermoFramework` handles activity, but no
   heat-duty/enthalpy-flow computation surfaced in anything read this
   session. Recommendation unchanged: keep mass-accounting unification and
   energy-accounting as separate design threads.
6. **Construction/validation parity with `ControlVolume`.** Does a
   `FlowBoundary`-based `Reservoir` need the same phase-key existence
   checks `ControlVolume.__init__` runs for `internal_interfaces`
   ([`control_volume.py:241-253`](../../src/core/control_volume.py))?
7. **Performance.** Does `source.phases[source_phase_key]` lookup inside
   every `FlowBoundary` implementation (vs. receiving `Phase` objects
   directly, as `PhaseInterface.compute_flux` does today) have a real cost
   worth checking, or is it negligible relative to the surrounding
   speciation/reaction work? Not measured.

## Trigger conditions

Pick this up together — `Reservoir` and `FlowBoundary` are no longer
separable motivations. Concrete triggers: a need to audit mass balance
across an `ExternalBoundary` (closing a whole-model mole balance for a
paper figure or regression sentinel), a model needing a shared feed/vent
reservoir (net atmosphere exchange), or a desire to simplify the
`ControlVolume`/`Simulation` orchestration surface by collapsing three
dispatch loops into one shared-type loop. Not worth starting speculatively
— see §3.3.

## Cross-references

- [`MASS_EXCHANGE_ARCHITECTURE.md`](../design/MASS_EXCHANGE_ARCHITECTURE.md)
  §2 (two-axis taxonomy this note extends), §4 (`CVLink` → `InterzonalFlow`
  rename — proposed, not shipped; now in tension with `FlowBoundary`, see
  Open Question 1), §6.2 (the `WaterVapourBoundary`
  infinite-liquid-reservoir mistake — the precedent for why a genuinely
  stateful reservoir should become a real tracked `Phase`, not a stand-in),
  §7 (`TransportLink` sketch — the direct ancestor of §5's resolution),
  §12 Q1 (the unification question §5 now answers).
- [`PHENOMENA_PROTOCOL.md`](PHENOMENA_PROTOCOL.md) — sibling design
  discussion, same session; `PartitionModel`'s relationship to both notes
  is the explicitly-open shared thread (§5.5 here, open question there).
- [`src/core/interfaces.py`](../../src/core/interfaces.py) — `PhaseInterface`
  protocol.
- [`src/core/links.py`](../../src/core/links.py) — `CVLink` protocol.
- [`src/core/boundaries.py`](../../src/core/boundaries.py) —
  `ExternalBoundary` protocol and its six concrete implementers.
- [`src/core/solvers.py`](../../src/core/solvers.py) —
  `SequentialAdvanceSolver.solve_step`, where `ExternalBoundary.compute_flux`
  is actually invoked (step 2b).
- [`src/core/simulation.py`](../../src/core/simulation.py) — `_apply_links`,
  where `CVLink.compute_flow` is actually invoked.
