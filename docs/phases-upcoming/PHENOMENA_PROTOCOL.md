# PHENOMENA_PROTOCOL — a shared taxonomy for reactions and transfer models

> **Status:** Design discussion, 2026-07-10. Not yet started — no branch, no
> checklist, no code written. Phase 1 of two.
>
> **Depends on nothing shipped.** Extends the constraint-family taxonomy
> already on record in
> [`MASS_EXCHANGE_ARCHITECTURE.md`](../design/MASS_EXCHANGE_ARCHITECTURE.md)
> §14.
>
> **Followed by:**
> [`KINETIC_TRANSFER_GENERALIZATION.md`](KINETIC_TRANSFER_GENERALIZATION.md)
> (Phase 2) — the behavior-changing half of this idea. Phase 2 needs this
> doc's protocol shape settled first, but Phase 1 is independently useful
> even if Phase 2 never ships (see §5 below).
>
> **Related, not merged:** [`RESERVOIR_TYPE.md`](RESERVOIR_TYPE.md) — a
> sibling design discussion from the same session, now concluding with a
> real `FlowBoundary` unification for `PhaseInterface`/`CVLink`/
> `ExternalBoundary`. `PartitionModel`'s relationship to *both* notes is
> the explicitly-open shared thread between them — not resolved in either.
>
> **Revision note (2026-07-10):** `ReactionSystem` is renamed
> `PhenomenaSystem` in this note's discussion below, reflecting the
> broadened taxonomy — see §7. Its actual responsibilities (bucketing,
> engine ownership, `compute_rates` execution, monitor attachment) are
> unchanged; this is a rename, not a redesign.

---

## Where this landed

Started from: is there a cleaner formulation grouping reactions (kinetic
and equilibrium) and wider mass-transfer phenomena (`HenryEquilibrium`,
`KspEquilibrium`, `RaoultEquilibrium`, …) under one conceptual type, so
they interface more easily with `Reaction`, `PartitionModel`, and
`Species`?

A single flat `ChemicalPhenomena` spanning everything doesn't work — kinetic
reactions (rate-producing, `dy/dt`) and equilibrium constraints (algebraic,
`0 = g(y,z)`) have genuinely different consumption/solver paths (the
differential/algebraic split in `SOLVER_ARCHITECTURE.md`), and isotherms
(Langmuir/Freundlich — not yet implemented anywhere in the codebase) aren't
mass-action-reformattable at all, per §14.1's own ruling.

Resplitting along **equilibrium vs. kinetic** — rather than by domain
(chemistry vs. transport) — is the right axis, because it aligns with what
actually determines the solver path rather than crossing it:

- **`EquilibriumPhenomena` is already substantially real.**
  `HenryEquilibrium`/`KspEquilibrium`/`RaoultEquilibrium` already
  dual-satisfy `EquilibriumConstraint` and `PartitionModel` today
  ([`partition.py:221-253`](../../src/chemistry/partition.py)), and
  `LAYER1_GAP_CLOSURE` already folds gas-liquid equilibrium rows into the
  *same* simultaneous Newton solve as acid-base equilibria — genuine
  solver-level convergence, not just a naming coincidence.
- **`KineticPhenomena` looked asymmetric at first**, but the asymmetry
  turned out to be a design choice, not a physical necessity —
  `SimultaneousEulerSolver`/`SimultaneousAdaptiveSolver` already evaluate
  kinetic reactions and kinetic `PhaseInterface` transfer from one frozen
  snapshot into one combined RHS
  ([`solvers.py:756-799`](../../src/core/solvers.py)) — the kinetic-side
  analog of `LAYER1_GAP_CLOSURE`.

Conclusion: a three-tier `Phenomena → {EquilibriumPhenomena,
KineticPhenomena}` hierarchy is well-supported. The admission criterion for
either branch isn't domain or the kinetic/equilibrium label itself — it's
§14.2's existing litmus test: **is the type a pure evaluate-at-a-point
function, or does it nest an internal solve?** `MultispeciesVLEPartition`
already fails this test today and stays excluded on purpose.

This document scopes **only the taxonomy** — new protocol definitions, no
behavior change, no rename of anything load-bearing without a separate
decision. Making `KineticTransferModel` actually satisfy `KineticPhenomena`
(giving it a real, self-executing rate function instead of being inert
config) is `KINETIC_TRANSFER_GENERALIZATION.md`'s job.

**How `Phenomena` connects to `ControlVolume` (added 2026-07-10):** through
`PhenomenaSystem` (the renamed `ReactionSystem`), not directly.
`ControlVolume` never gains a `phenomena: List[Phenomena]` field of its own
— `PhenomenaSystem` keeps doing the bucketing and owns the resulting
equilibrium engine, exactly as `ReactionSystem` does today. §7 below
records why a direct `ControlVolume ↔ Phenomena` link was considered and
rejected: the equilibrium engine's construction and ownership has nowhere
else sensible to live without either duplicating `PhenomenaSystem`'s job on
`ControlVolume` or leaving the engine ownerless.

---

## Background and motivation

Today, reaction-like declarations are bucketed by `ReactionSystem` via
`isinstance` at construction
([`reaction_system.py:139-143`](../../src/reactions/reaction_system.py)):
`KineticReaction`, `BlackBoxReactionModel`, and an `EquilibriumConstraint`-
satisfying bucket (`EquilibriumReaction` and its
`HenryEquilibrium`/`KspEquilibrium`/`RaoultEquilibrium` siblings). None of
this is unified under a literal shared type — `equilibrium.py`'s own
docstring says so outright: "There is no shared base class; shared
validation logic lives as free functions in `_shared.py`."

That's not an oversight. Those free functions
(`species_ids_from_entries`, `phases_from_entries`, `coerce_and_validate`,
`is_cross_phase_from_entries`) already operate uniformly across all three
declaration classes today — the *structural* commonality already exists in
practice, just not as a named type. `Phenomena` would mostly formalize
something already true, not invent new machinery.

---

## Key design decisions

### 1. The three-tier shape

```
Phenomena  «Protocol»
    stoichiometry: Sequence[StoichiometryEntry]
    species_ids: List[str]
    label: str

├── EquilibriumPhenomena  «Protocol»
│     + log_K: float
│     + dH_J_per_mol: Optional[float]
│     + T_ref_K: float
│     Implementers (today): EquilibriumReaction, HenryEquilibrium,
│                            KspEquilibrium, RaoultEquilibrium
│
└── KineticPhenomena  «Protocol»
      + a rate-evaluation method (already satisfied today by
        KineticReaction.compute_rates; the transfer-side shape is
        KINETIC_TRANSFER_GENERALIZATION.md's job)
      Implementers (today): KineticReaction
      Implementers (Phase 2): generalized KineticTransferModel
```

`Phenomena` itself is thin on purpose — just enough for
species/stoichiometry-level introspection (validation, balance-checking,
"what touches species X regardless of kind") without claiming anything
about execution. This mirrors the black/gray/white-box tiering already
used for `ChemicalEquilibriumEngineProtocol` — capability tiers layered on
a minimal base, not one flat interface.

### 2. Admission criterion

Restate §14.2's litmus test explicitly as the actual boundary for either
sub-protocol, not domain and not the kinetic/equilibrium label: **pure
evaluate-at-a-point, never an internal iteration or convergence loop.**
This is why `MultispeciesVLEPartition.equilibrium_all_a_moles` (which
solves internally) stays excluded from both branches — nesting a solver
inside a generic dispatch is the same anti-pattern §14.2 already names for
the PHREEQC dual-write and precipitation active-set-loop cases.

### 3. `EquilibriumPhenomena` vs. `EquilibriumConstraint` — rename in place, or parallel protocol?

Open question, leaning toward **rename in place** rather than adding a new
protocol alongside the existing one. `EquilibriumConstraint` already *is*
`EquilibriumPhenomena` in every respect that matters — introducing a
second name for the same job would repeat the exact "two objects doing
structurally the same thing" mistake `MASS_EXCHANGE_ARCHITECTURE.md` §14
already flags for the Ksp/Henry folding decision. Needs confirming there
are no external call sites (demos, docs, `EQUILIBRIUM_CONSTRAINT_
UNIFICATION`'s own shipped surface) that would break from a rename versus
being fine with a compatible extension.

### 4. Where does `BlackBoxReactionModel` fit?

Open. It's a `ReactionModel` implementer today (opaque external simulator
adapter, runtime balance checking) — but does it have a meaningful
`stoichiometry` to expose declaratively the way `KineticReaction` does, or
is that exactly the thing that makes it "black box"? May need to stay
outside `Phenomena` entirely, satisfying only `ReactionModel`.

### 5. Isotherms stay explicitly excluded

Record this rather than let it be silently reconsidered later:
Langmuir/Freundlich adsorption isotherms don't belong in either branch —
not `EquilibriumPhenomena` (§14.1: "no log-K exists for a saturable-
capacity curve"), and not `KineticPhenomena` either (they don't exist in
the codebase yet at all — `MASS_EXCHANGE_ARCHITECTURE.md` §5.3 lists them
"❌ not yet"). If/when they're implemented, they remain a distinct
constraint type on `PartitionModel`/`MultispeciesPartitionModel`, per
§14.1's own ruling.

### 6. Is Phase 1 worth doing even if Phase 2 never ships?

Yes, deliberately scoped that way. A named `Phenomena`/`EquilibriumPhenomena`
protocol pair is useful on its own — for documentation, for UML-style
architecture figures (the motivating use case this session), and for any
future generic species-consistency tooling — independent of whether
`KineticTransferModel` ever gets generalized. Sequencing it first, cheaply,
de-risks Phase 2 rather than blocking on it.

### 7. `PhenomenaSystem`: keep a dedicated System-shaped type

Rename `ReactionSystem` → `PhenomenaSystem`, keeping its actual role
unchanged (bucket `Phenomena` implementers, own the resulting engine,
execute kinetic rates, provide the `ControlVolume` attachment point). Not
absorbed into `ControlVolume` directly —
`ChemicalEquilibriumEngineProtocol.from_reactions()` is already a
standalone, public classmethod (used directly in dozens of demos/tests and
inside `models/vlmodels/adm1/base.py`'s own builder), so a dedicated System
type isn't gatekeeping external construction access; it's just keeping
bucketing, `compute_rates` execution, and monitor attachment off
`ControlVolume`'s already-large constructor.

---

## Open questions

1. **Naming.** `Phenomena` was the working name; not checked against
   existing identifiers or for collision with domain vocabulary.
2. **Rename vs. parallel protocol for `EquilibriumConstraint`** (§3 above).
3. **`BlackBoxReactionModel`'s membership** (§4 above).
4. **Is a formal `Protocol` actually warranted**, or would documenting the
   existing informal `_shared.py` convention be enough without adding a new
   type name? Worth asking honestly rather than assuming the type is
   needed — same discipline applied to `Reservoir` in the sibling note.
5. **Sequencing relative to Phase 2** — should this ship independently
   regardless of Phase 2's fate (§6 above says yes; worth confirming that
   holds once Phase 2's actual scope is better known).

## Trigger conditions

Pick this up before starting `KINETIC_TRANSFER_GENERALIZATION.md`'s actual
implementation (that phase needs `KineticPhenomena`'s method shape settled
first) — or independently, whenever a concrete need for
cross-kind species/phenomenon introspection shows up (a generic
consistency-checking tool, or a future architecture-figure/documentation
pass that would benefit from a named, diagrammable taxonomy rather than
`_shared.py`'s implicit one).

## Cross-references

- [`MASS_EXCHANGE_ARCHITECTURE.md`](../design/MASS_EXCHANGE_ARCHITECTURE.md)
  §14 (constraint-family taxonomy and folding limits — the direct ancestor
  of this note), §14.1 (two constraint families), §14.2 (the litmus test
  this note's admission criterion is restated from).
- [`SPECIATION_REACTIONMODEL_BOUNDARY.md`](../design/SPECIATION_REACTIONMODEL_BOUNDARY.md)
  — the differential/algebraic split rationale for why `KineticReaction`
  and equilibrium constraints were never unified in the first place.
- [`../phases-shipped/EQUILIBRIUM_CONSTRAINT_UNIFICATION.md`](../phases-shipped/EQUILIBRIUM_CONSTRAINT_UNIFICATION.md)
  — shipped precedent for exactly this shape of phase (declaration-API
  generalization, no solve-time behavior change).
- [`../phases-shipped/LAYER1_GAP_CLOSURE.md`](../phases-shipped/LAYER1_GAP_CLOSURE.md)
  — shipped precedent for the equilibrium-side solver convergence this
  note cites as evidence.
- [`KINETIC_TRANSFER_GENERALIZATION.md`](KINETIC_TRANSFER_GENERALIZATION.md)
  — Phase 2, depends on this doc.
- [`RESERVOIR_TYPE.md`](RESERVOIR_TYPE.md) — sibling design discussion,
  same session, same underlying question applied to the transport-topology
  protocols instead.
- [`src/reactions/equilibrium.py`](../../src/reactions/equilibrium.py) —
  `EquilibriumConstraint` protocol.
- [`src/reactions/kinetic.py`](../../src/reactions/kinetic.py) —
  `KineticReaction`, today's only real `KineticPhenomena`-shaped object.
- [`src/reactions/reaction_system.py`](../../src/reactions/reaction_system.py)
  — the `isinstance`-based bucketing this note's taxonomy would sit
  alongside (not necessarily replace).
- [`src/chemistry/partition.py`](../../src/chemistry/partition.py) —
  `HenryEquilibrium`/`KspEquilibrium`/`RaoultEquilibrium`'s dual
  `EquilibriumConstraint`/`PartitionModel` satisfaction.
- [`src/core/solvers.py`](../../src/core/solvers.py) —
  `SimultaneousAdaptiveSolver`'s combined-RHS closure, the kinetic-side
  convergence evidence.
