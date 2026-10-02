# Growth Stoichiometry — Design Note

> Status: design note, not yet started. No branch, no checklist, no code.
> Written 2026-10-01 while settling how the stirred-tank template names its
> nitrogen source; the template keeps `n_source_id` until this lands.

## What growth balances today

`ReactionBuilder.aerobic_growth` (`reactions/kinetic/builder.py`) solves the
growth stoichiometry per mole of substrate in closed form, one unknown per
balanced element:

| Element | Solved for |
|---|---|
| C | CO2 produced |
| N (only with `balance="CHNO"`) | nitrogen source consumed |
| H | H2O produced |
| O | O2 consumed |

The nitrogen source may carry H and O, and both are counted. Beyond that:

- **The participants are fixed.** Substrate, O2, biomass, CO2, H2O and
  optionally a nitrogen source. A product such as ethanol, or a sulfur or
  phosphorus source, cannot be added.
- **Other elements are ignored silently.** The reaction's balance check
  covers only the balanced elements, so a source such as `(NH4)2SO4` leaves
  its S (and the sulfate's O) unbalanced with no error.
- **Biomass nitrogen under `"CHO"` is not flagged.** The stirred-tank
  tutorials use `Yeast` (N 0.16) with `balance_basis="CHO"`; the biomass
  nitrogen appears from nowhere and the check does not see it.
- **No charge balance.** A charged source (NH4+) leaves the reaction
  charged; there is no H+ term to close it.

## Direction

### A solver for one reaction from its participants

The caller lists what a reaction consumes and produces, a basis species
(coefficient −1) and the yields they know; the solver finds the other
coefficients:

```python
ReactionBuilder.balanced_reaction(
    consumed=[GLUCOSE, O2, NH3],
    produced=[YEAST, CO2, H2O],
    basis=GLUCOSE,
    yields={YEAST: 0.49},          # g per g of the basis
    rate_fn=..., label="oxidative_growth",
)
```

The unknowns are the coefficients other than the basis. The equations are
one balance per element present in the participants, a charge balance, and
one equation per yield (a mass yield is linear in the coefficients once the
molar masses are known). The system is linear, so it is solved directly
(no root finder), and its rank decides the outcome:

| Case | Meaning | Outcome |
|---|---|---|
| rank = unknowns, consistent | exactly specified | solved |
| rank < unknowns | under-specified: *k* degrees of freedom | raises, naming how many constraints are missing and which coefficients are free |
| more equations than the rank, consistent | redundant constraints | solved |
| inconsistent | no coefficients satisfy everything | raises, with the residual per element or yield (usually a missing participant, e.g. no nitrogen source for a biomass containing N) |

It also warns when a coefficient has the wrong sign (a listed product is
consumed, typically an impossible yield) and names participants the
balances cannot tell apart (e.g. NH3 and NH4+ without a charge balance).
An under-specified reaction never falls back to a least-squares guess.

`aerobic_growth` and `monod_aerobic_growth` become presets over this solver
with their current participants; their results are unchanged to roundoff.
The unbalanced biomass nitrogen above becomes the inconsistent case, with
its residual.

### Several reactions, not fixed splits

Where a substrate is used in more than one way, each way is its own
reaction with its own rate law, each exactly specified by its balances and
its yield. For *S. cerevisiae* (Sonnleitner–Käppeli): oxidative growth on
glucose; fermentative growth on glucose to ethanol and CO2; oxidative growth
on ethanol. The split between them follows from the rate laws and changes
with conditions (the respiratory quotient rises when respiration saturates).
A single reaction with a fixed product ratio cannot do that, so fixed ratios
are not part of the first version.

### Template

The template resolves each reaction's participants against the model's
species by the rules of `.organism()` / `.substrate()` (an id among the
species, an id with `atoms=`, or a `Species`). `balance_basis` and
`n_source_id` give way to the participants listed.

## Open questions

1. The template's shape for several reactions on one substrate: a
   `pathways=` argument on `.substrate()`, a `.reaction(...)` /
   `.pathway(...)` call per reaction, or both.
2. Whether a non-balanced element in the biomass raises or warns. Raising
   makes the three stirred-tank tutorials that use `Yeast` with `"CHO"`
   fail until they switch to `Yeast_CHO` or list a nitrogen source.
3. Whether fixed ratios between two participants (e.g. a measured
   respiratory quotient) are added later, for deliberately lumped
   reactions fitted at one operating point and as a consistency check on
   measured data. If so, yields and ratios are both linear constraints and
   could share one `constraints=` argument.
4. Whether the charge balance is always on, or on only when a participant
   is charged (with H+ then required among the participants).

## Relationship to other notes

- [EXPLICIT_SPECIES_RESOLUTION.md](EXPLICIT_SPECIES_RESOLUTION.md): its
  checkpoint 9 resolves `n_source_id` against the model's species and
  leaves the shape of growth participants to this note.
