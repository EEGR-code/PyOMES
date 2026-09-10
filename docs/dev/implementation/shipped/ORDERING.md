# Advance Step Ordering — Design Note

## The question

Should external source terms (feeds, dosing) be applied before or after reactions
within a single `ControlVolume.advance()` call?

## Conclusion

**Feed before reactions.** There is no physically meaningful case for the reverse
ordering within a single Euler timestep.

## Reasoning

In any continuous or semi-continuous process, feed arrival and biological/chemical
reaction are concurrent. Within a single timestep, material that enters the CV is
immediately available for reaction. Applying reactions first and feed second has no
physical analogue — it implies that feed cannot be consumed until the following step,
which introduces an artificial lag.

Cases that might appear to support feed-after-reactions all resolve at the orchestration
level, not within a single `advance()` call:

- **Sequential batch reactors (SBR):** FILL and REACT are distinct simulation phases
  with different CV configurations. During REACT there is no feed; during FILL there
  may be no reaction. The ordering within a single step is irrelevant because the two
  never coexist in the same call.

- **Bolus pH correction:** Acid or base is added only after observing pH drift. But
  this is a control action executed *between* steps by the control system — the
  controller reads the result of one `advance()` and doses before the next begins.

- **Fed-batch substrate depletion:** Feed is triggered when substrate falls below a
  threshold. Again, this is an orchestrator-level decision about *when* to feed, not
  about the ordering of feed and reaction within a timestep.

In all of these cases the sequencing logic belongs above `advance()`. The current
feed-after-reactions ordering in `ControlVolume.advance()` is a historical accident,
not a deliberate physical modelling choice.

## Correct advance sequence

```
1. Run property solvers   →  prop_results    (speciation, pH from current moles)
2. Apply external sources                    (feed arrives; state-dependent flux
                                              magnitudes computed from prop_results)
3. Integrate reactions                       (react what is now present after feeding)
4. step_internal_transfer                    (gas-liquid equilibration using prop_results)
5. Return AdvanceResult(properties=prop_results)
```

Property solvers run first so that state-dependent boundaries (pressure vents, pH
dosing) can compute their flux magnitude from current speciation before any moles
change. All subsequent steps — feed, reactions, and transfer — then operate on a
consistent set of pre-step properties, with the same O(dt) operator-splitting error
throughout.

## Relationship to EulerSnapshotSolver

The `EulerSnapshotSolver` in `GasLiquidVolume` already applies this logic correctly:
it computes all fluxes (boundaries, transfer, reactions) from a frozen pre-step
snapshot and applies all deltas simultaneously. The plain `ControlVolume.advance()`
never received the same treatment. The refactor described in `CV_UPDATE.md` brings
`ControlVolume.advance()` into alignment with the snapshot solver's behaviour, and
sets the foundation for eventually removing `GasLiquidVolume` as a separate class.
