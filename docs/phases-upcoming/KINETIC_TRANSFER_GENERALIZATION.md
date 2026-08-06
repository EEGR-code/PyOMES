# KINETIC_TRANSFER_GENERALIZATION — arbitrary rate laws for phase transfer

> **Status:** Design discussion, 2026-07-10. Not yet started — no branch,
> no checklist, no code written. Phase 2 of two.
>
> **Depends on:** [`PHENOMENA_PROTOCOL.md`](PHENOMENA_PROTOCOL.md)
> (Phase 1) — needs `KineticPhenomena`'s protocol shape settled before a
> concrete class can be implemented against it. Phase 1 has **not
> shipped**; nothing here should start before it does.
>
> **Related, not merged:** [`RESERVOIR_TYPE.md`](RESERVOIR_TYPE.md) —
> sibling design discussion, same session, now concluding with a
> `FlowBoundary` protocol unifying `PhaseInterface`/`CVLink`/
> `ExternalBoundary`. `PhaseInterface` below is `FlowBoundary`'s same-CV
> case (`RESERVOIR_TYPE.md` §5.1) — a naming/framing note only, nothing in
> this phase's design changes as a result.

---

## Where this landed

`KineticTransferModel` today ([`transfer_models.py`](../../src/core/transfer_models.py))
is pure declarative config — `partition_model`, `k_transfer`,
`transfer_basis` — with no execution logic of its own. All actual work is
deferred to `KineticGasLiquidLink`, built by a factory function
(`_build_transfer_link`,
[`control_volume.py:46-110`](../../src/core/control_volume.py)) that
unpacks a `transfer_models=` dict into the link's `partition_models`/`kLa`/
`equilibrium_species` fields. `KineticReaction`, by contrast, self-executes
via an arbitrary user-supplied `rate_fn(env) -> float` and directly
satisfies `ReactionModel` itself
([`kinetic.py:48-80`](../../src/reactions/kinetic.py)).

This asymmetry is a design choice, not a physical necessity. Transfer
kinetics today is constrained to first-order linear relaxation
(`flux = kLa × V × (C* − C)`) — a strict special case of what an arbitrary
`rate_fn` could already express (the same physics is writable today as an
ordinary `KineticReaction` with cross-phase stoichiometry; nothing wires
that up as "the" transfer path currently). Resolving the asymmetry is
desirable independent of `PHENOMENA_PROTOCOL.md`'s taxonomy question — it
generalizes `VLsim` to model non-linear transfer physics (saturable
membrane transport, biofouling-limited `kLa`, …) that linear `kLa` cannot
express, and matches the existing "sensible default + power-user `rate_fn`
escape hatch" pattern already shipped in `ReactionBuilder.aerobic_growth`.

This is the higher-risk half of the two-phase split: it touches a shipped,
tested, documented API (`transfer_models=`), likely consumed by reference
models in `models/vlmodels` (BSM2/ADM1 use kLa-style gas-liquid transfer).
Any change here must be additive and non-breaking by default.

---

## Background and motivation

Recap of the current call chain (see also the earlier session discussion
captured informally in this doc's git history, not reproduced here):

```
transfer_models={"O2": KineticTransferModel(HenryEquilibrium(...), k_transfer=150.0)}
        │  (per-species config, not itself a PhaseInterface)
        ▼
_build_transfer_link()  →  one KineticGasLiquidLink instance
        │  (the actual PhaseInterface implementer; direct construction
        │   is deprecated — DeprecationWarning fires in __post_init__,
        │   suppressed internally by the factory)
        ▼
appended to cv.internal_interfaces
```

`KineticGasLiquidLink` itself implements **two** protocols on one object —
`PhaseInterface` (`phase_a_key`/`phase_b_key`/`compute_flux`) and `CVLink`
(`source_cv_key`/`source_phase_key`/`sink_cv_key`/`sink_phase_key`/
`compute_flow`) — with `compute_flux` as a thin adapter that builds a
synthetic `cvs` dict and delegates to `compute_flow`
([`gas_liquid_link.py:44-54`](../../src/core/gas_liquid_link.py)).

---

## Key design decisions

### 1. Rate-function signature for transfer

`KineticReaction.rate_fn(env: ReactionEnvironment)` reads one phase's
state. Transfer physics is inherently two-sided — the driving force
depends on both phases, not one. `ReactionEnvironment` doesn't currently
expose that. Options, not yet decided between:

- A new, purpose-built `TransferEnvironment` bundling both phases'
  concentrations/`T_K`/`V_L` (mirrors `ReactionEnvironment`'s role, but
  two-sided).
- A `rate_fn(state_a: Phase, state_b: Phase, T_K: float, dt_h: float, ...) -> Dict[str, float]`
  shape closer to `PhaseInterface.compute_flux`'s own signature directly.

The second option has the advantage of looking like what `KineticPhenomena`
(Phase 1) would actually require methodwise; the first mirrors existing
naming conventions more closely. Needs deciding before implementation.

### 2. `partition_model` optionality

The default linear-relaxation mode still needs `partition_model` (it's
where `C*`, the equilibrium target, comes from). A fully custom `rate_fn`
might not need it (a pure empirical correlation) or might want to consume
it directly (a non-linear function of `C − C*`, still referencing the same
equilibrium target). Recommended default: keep `partition_model` required,
let an optional `rate_fn` override the built-in linear-relaxation
calculation entirely when supplied. **No `rate_fn` supplied → today's exact
linear behavior, byte-identical.** This keeps the change strictly additive.

### 3. `KineticGasLiquidLink`'s role — two candidate shapes

**(a) Full self-execution.** `KineticTransferModel` itself gains
`phase_a_key`/`phase_b_key`/`compute_flux` and directly satisfies
`PhaseInterface` — `FlowBoundary`'s same-CV case, per `RESERVOIR_TYPE.md`
§5.1. `KineticGasLiquidLink`'s role shrinks to a thin per-CV
aggregator: iterate a set of self-executing per-species `KineticPhenomena`/
`EquilibriumPhenomena` instances and merge their fluxes into one
`PhaseInterface`-shaped return. This is the "real" fix — it's what makes
`KineticTransferModel` actually satisfy `KineticPhenomena` (Phase 1's
protocol) rather than staying a special case documented as an exception.

**(b) Aggregator-only generalization.** `KineticTransferModel` stays inert
config; only `KineticGasLiquidLink`'s internals change, to call each
entry's optional `rate_fn` instead of always using the built-in linear
formula. Lower risk, smaller diff, ships without redefining
`KineticGasLiquidLink`'s role or touching the `PhaseInterface` protocol
question at all — but doesn't actually resolve the structural asymmetry
`PHENOMENA_PROTOCOL.md` identified; `KineticTransferModel` still wouldn't
self-execute or literally satisfy `KineticPhenomena`.

**Recommendation:** evaluate (a) as the target shape, but treat (b) as a
legitimate smaller first increment if (a) turns out to require touching
too much of `KineticGasLiquidLink`/`_build_transfer_link` at once. Not
decided — flag explicitly in a checklist if/when this phase starts, rather
than deciding by default during implementation.

### 4. Backward compatibility

`transfer_models=` is shipped, tested, and (per `MASS_EXCHANGE_ARCHITECTURE.md`
§5.3) the "preferred path" for gas-liquid transfer already. Constraints:

- Existing `KineticTransferModel(partition_model, k_transfer)` calls with
  no `rate_fn` must produce identical numerical results.
- `KineticGasLiquidLink`'s deprecation-on-direct-construction behavior must
  be explicitly preserved or explicitly re-decided — not silently dropped
  as a side effect of restructuring its role.
- BSM2/ADM1 reference-model golden-trajectory sentinels (re-baselined
  multiple times across past phases, per project history) should be
  re-verified once any concrete implementation lands, following the
  project's established "trust but verify" numerics discipline.

### 5. Solver-scoping caveat

The "simultaneous RHS" benefit this phase is partly motivated by
(`PHENOMENA_PROTOCOL.md`'s citation of
[`solvers.py:756-799`](../../src/core/solvers.py)) already works today for
*existing* `KineticTransferModel` instances via `internal_interfaces` —
generalizing the rate function doesn't change *which* solvers see it.
`SequentialAdvanceSolver`'s `step_internal_transfer` path must keep working
identically for the default (non-adaptive) case. Confirm this phase implies
no default-solver behavior change — only new capability, not a new
requirement.

### 6. Does `EquilibriumTransferModel` need the same treatment?

Open. It already has full generality via any `PartitionModel.partition_ratio()`
implementation — arguably it doesn't need a `rate_fn` escape hatch the way
the kinetic side does, since "equilibrium" by definition has no rate to
customize (only the target). But leaving it untouched while generalizing
only the kinetic side reintroduces a milder version of the very asymmetry
this phase exists to fix — worth an explicit decision rather than a
default-by-omission.

---

## Open questions

1. **Rate-fn signature** — finalize (§1).
2. **`partition_model` optionality** — finalize (§2), default recommended
   above.
3. **`KineticGasLiquidLink`'s role** — (a) full self-execution vs. (b)
   aggregator-only generalization (§3); may ship (b) first.
4. **`EquilibriumTransferModel` parity** — does it need anything, or does
   `PartitionModel` already cover it (§6)?
5. **Regression scope** — which existing tests/sentinels need re-running
   once a concrete implementation exists (§4).

## Trigger conditions

A concrete modelling need for non-linear transfer kinetics that linear
`kLa` cannot express (membrane fouling, saturable transport, or similar) —
or naturally once `PHENOMENA_PROTOCOL.md` has shipped and
`KineticPhenomena`'s method shape is settled enough to implement against.
Not worth starting speculatively before Phase 1 lands, since this phase's
whole point is making `KineticTransferModel` satisfy a protocol that
doesn't exist yet.

## Cross-references

- [`PHENOMENA_PROTOCOL.md`](PHENOMENA_PROTOCOL.md) — Phase 1, depended on.
- [`MASS_EXCHANGE_ARCHITECTURE.md`](../design/MASS_EXCHANGE_ARCHITECTURE.md)
  §3 (`PartitionModel`/`TransferModel`/`PhaseInterface` invariant — "a
  `TransferModel` never knows chemistry"), §5.3 (gas-liquid transfer
  implementation table).
- [`src/core/transfer_models.py`](../../src/core/transfer_models.py) —
  `KineticTransferModel`/`EquilibriumTransferModel` as they exist today.
- [`src/core/gas_liquid_link.py`](../../src/core/gas_liquid_link.py) —
  `KineticGasLiquidLink`, the actual executor and dual `PhaseInterface`/
  `CVLink` implementer — now `FlowBoundary`'s two topology cases on one
  object, per `RESERVOIR_TYPE.md` §5.1.
- [`src/reactions/kinetic.py`](../../src/reactions/kinetic.py) —
  `KineticReaction`, the existing arbitrary-`rate_fn` precedent this phase
  generalizes transfer models toward.
- [`src/reactions/builder.py`](../../src/reactions/builder.py) —
  `ReactionBuilder`'s custom-`rate_fn` escape-hatch precedent (defaults +
  power-user override, the same pattern this phase applies to transfer).
- [`src/core/solvers.py`](../../src/core/solvers.py) —
  `SimultaneousAdaptiveSolver`'s combined-RHS closure; whatever this phase
  builds must keep working correctly inside it.
