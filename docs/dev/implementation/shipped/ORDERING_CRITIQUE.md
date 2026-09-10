# ORDERING.md — Expert Panel Critique

Three independent reviews of the claims and recommendations in ORDERING.md.

---

## Panel

1. **Bioprocess specialist** — reactor modelling conventions, ADM1/BSM2, fermentation
2. **Numerical methods specialist** — ODE integration, operator splitting, consistency
3. **Software architect** — API design, maintainability, separation of concerns

---

## Review 1: Bioprocess Specialist

### What the document gets right

The distinction between orchestrator-level sequencing (SBR phases, fed-batch triggering,
control actions) and within-timestep ordering is correct and useful. These are genuinely
different concerns and the document is right to separate them.

The claim that feed-before-reactions better reflects a continuously fed vessel is
reasonable for CSTR-type systems. In a CSTR, substrate that enters the vessel is
immediately available for microbial uptake — there is no physical hold-up between
arrival and reaction.

### What the document gets wrong or omits

**The standard CSTR formulation treats feed and reaction as concurrent, not sequential.**
The textbook mass balance for a CSTR is:

```
dS/dt = D × (S_in − S) + r(S)
```

The dilution term `D × (S_in − S)` and the reaction term `r(S)` are both evaluated at
the same state `S` and summed into a single derivative. Neither is "before" or "after"
the other — they are simultaneous contributions to the rate of change. The ordering
question only arises because operator splitting breaks this single ODE into sequential
sub-steps.

This means the document's framing — "feed arrives, then reactions consume what is now
present" — imposes a causal sequence that does not exist in the continuous model. It is
one valid discretisation, but it is not more "physically correct" than the reverse. Both
introduce O(dt) splitting error relative to the true concurrent formulation.

**Reference implementations may deliberately use the reverse ordering.** BSM2 and
PyADM1 integrate the dilution and reaction terms together inside a single ODE
right-hand side — they do not split them at all. If this codebase intends to reproduce
BSM2 results for validation, the ordering choice must be checked against the reference
implementation rather than decided on physical intuition alone.

### Score: 6/10

The recommendation is defensible but overstated. The assertion "there is no physically
meaningful case for feed after reactions" should be softened to "within a sequential
operator-split scheme, feed-before-reactions is the more natural convention for
substrate-limited kinetics in a CSTR."

---

## Review 2: Numerical Methods Specialist

### What the document gets right

The document correctly identifies that all sub-steps share the same O(dt) operator-
splitting error and that pre-step property evaluation is a clean approach for explicit
Euler.

### What the document gets wrong or omits

**The proposed sequence is NOT aligned with the EulerSnapshotSolver.** This is the
most significant error in the document. The final section claims the refactor "brings
`ControlVolume.advance()` into alignment with the snapshot solver's behaviour." It
does not.

The EulerSnapshotSolver works as follows:

1. Freeze state at time t (snapshot)
2. Compute ALL fluxes from the frozen snapshot — boundaries, reactions, transfer —
   all evaluated at the same unmodified state
3. Sum all flux contributions into a single delta vector
4. Apply the combined delta to the live state in one step

This is a **simultaneous** explicit Euler step. All sub-systems see the same state.
No sub-system's output affects another sub-system's input within the step.

The proposed `advance()` sequence is **sequential**:

1. Run property solvers on state at time t
2. Apply feed → **state is now modified** (t*)
3. Integrate reactions on state at t* → **state modified again** (t**)
4. Apply internal transfer on state at t**

Reactions now see post-feed moles, not the original pre-step moles. Transfer sees
post-feed-post-reaction moles. This is a Lie-Trotter operator splitting, not a
simultaneous Euler step. The two approaches produce different results at O(dt):

```
Snapshot:    S(t+dt) = S(t) + dt × [feed_rate(S(t)) + rxn_rate(S(t)) + transfer_rate(S(t))]
Sequential:  S* = S(t) + dt × feed_rate(S(t))
             S(t+dt) = S* + dt × rxn_rate(S*) + dt × transfer_rate(...)
```

These differ by O(dt²) cross-terms between feed and reaction, which is within the
overall O(dt) truncation error of Euler — so it is acceptable — but the document
should not claim equivalence.

**No discussion of splitting order sensitivity.** For substrate-limited Monod kinetics
where both feed addition and reaction consumption affect the same species, the ordering
of feed and reaction matters more than for independent sub-systems. Feed-before-
reactions gives reactions access to the fed substrate within the same step, which is
mildly more accurate for substrate-limited growth. Reactions-before-feed would
underestimate substrate availability. This is the strongest argument for the proposed
ordering — but the document doesn't make it in these terms.

**No mention of Strang splitting.** If splitting-order sensitivity is a concern, Strang
splitting (half-step feed, full-step react, half-step feed) would give second-order
splitting error. Not necessarily worth implementing, but worth acknowledging as the
principled alternative.

### Score: 5/10

The ordering recommendation is acceptable in practice but the claimed alignment with
the snapshot solver is incorrect, and the numerical justification is incomplete. The
strongest argument — that feed-before-reactions reduces splitting error for substrate-
limited kinetics — is absent.

---

## Review 3: Software Architect

### What the document gets right

The separation between orchestrator-level concerns and CV-level concerns is well drawn.
SBR phase sequencing, control actions, and feed-triggering logic genuinely do not
belong inside `advance()`.

The identification of the current ordering as a likely historical accident is plausible
and well-argued. The code shows signs of incremental development where reactions were
the initial concern and feeds were added later.

The proposed sequence is clean, easy to reason about, and will simplify the downstream
refactoring (eliminating `_last_properties`, enabling `KineticGasLiquidLink` as a
`PhaseInterface`).

### What the document gets wrong or omits

**The document frames the ordering as a physics question when it is equally a
software contract question.** The choice of ordering defines the semantics of
`advance()` — what state reactions see, what state transfer sees. This is an API
contract. The document should state this explicitly: "We are choosing the convention
that reactions operate on post-feed state. This is a modelling choice baked into the
`advance()` contract, and users writing custom `ReactionModel` implementations need
to understand it."

**No discussion of whether `advance()` should support multiple strategies.** The
document assumes a single fixed ordering. But the existence of both `EulerSnapshotSolver`
(simultaneous) and the proposed sequential approach suggests that different models may
benefit from different splitting strategies. A cleaner design might make the splitting
strategy configurable on the CV — similar to how `GasLiquidVolume` already accepts a
`solver` parameter — rather than hardcoding one ordering into `advance()`.

**The claim "historical accident" may be unfair.** If the original author was following
a reference implementation (BSM2, PyADM1) that evaluates feed and reaction
simultaneously in a single ODE derivative, the original code's ordering may not have
mattered at all — it was never intended to be the primary integration path. The
`EulerSnapshotSolver` was the intended primary solver, and `advance()` was kept as a
simple fallback. Calling it an accident dismisses the possibility that it was
intentionally low-priority.

### Score: 7/10

The recommendation is sound from a software perspective and will produce cleaner code.
The document would be stronger if it acknowledged the ordering as a deliberate
convention choice rather than a uniquely correct physical answer, and if it considered
whether the splitting strategy should be pluggable.

---

## Consolidated Scoring

| Claim | BP | NM | SA | Avg |
|-------|----|----|-----|-----|
| Feed-before-reactions is more physically natural | 6 | 7 | — | 6.5 |
| No meaningful case exists for the reverse | 5 | 5 | — | 5.0 |
| Counterexamples (SBR, bolus, fed-batch) correctly dismissed | 8 | — | 8 | 8.0 |
| Proposed sequence is correct and complete | 6 | 5 | 7 | 6.0 |
| Alignment with EulerSnapshotSolver | — | 3 | 5 | 4.0 |
| Overall recommendation | 6 | 5 | 7 | 6.0 |

**Legend:** BP = Bioprocess, NM = Numerical Methods, SA = Software Architect

---

## Key Recommendations

1. **Soften the certainty.** Replace "there is no physically meaningful case" with
   "within a sequential operator-split scheme, feed-before-reactions is the more
   natural convention." The continuous model is concurrent; any sequential ordering
   is a discretisation choice.

2. **Correct the EulerSnapshotSolver claim.** The proposed sequential advance is
   NOT equivalent to the snapshot solver. Both are valid explicit Euler variants
   with O(dt) error, but they produce different results. The document should
   acknowledge this difference rather than claiming alignment.

3. **Add the strongest numerical argument.** For substrate-limited kinetics (Monod),
   feed-before-reactions reduces splitting error because reactions see the correct
   substrate availability. This is the most concrete justification and it is missing.

4. **Acknowledge the ordering as a contract.** Users writing custom ReactionModel
   implementations should know that reactions will see post-feed state. This is a
   modelling convention, not a physical law.

5. **Consider making the strategy pluggable.** If the CV already has a configurable
   integrator for reactions, it could equally accept a configurable splitting strategy.
   This would allow the snapshot-style simultaneous approach as an option alongside the
   sequential approach, preserving backward compatibility and supporting different
   modelling conventions.
