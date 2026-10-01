# Growth Element Sources — Design Note

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

The nitrogen source may carry H and O, and both are counted. Nothing else
is:

- **Two fixed modes.** `"CHO"` and `"CHNO"`; nitrogen is the only element
  beyond C, H and O that can be balanced. Sulfur, phosphorus and the like
  cannot.
- **Other elements are ignored silently.** The reaction's balance check
  covers only the balanced elements, so a source such as `(NH4)2SO4` leaves
  its S (and the sulfate's O) unbalanced with no error.
- **Biomass nitrogen under `"CHO"` is not flagged.** The stirred-tank
  tutorials use `Yeast` (N 0.16) with `balance_basis="CHO"`; the biomass
  nitrogen appears from nowhere and the check does not see it.
- **No charge balance.** A charged source (NH4+) leaves the reaction
  charged; there is no H+ term to close it.

## Direction

- **One mapping from element to source.** `balance=` and `n_source=` are
  replaced by `element_sources={"N": NH3, "S": H2S}`. The balanced elements
  are C, H, O and the mapping's keys. The template takes the same mapping
  of ids (resolved against the model's species) or `Species`, in place of
  `balance_basis` and `n_source_id`.
- **A linear solve.** Unknowns: O2, CO2, H2O and one source per extra
  element. Equations: one balance per element. The square system is solved
  with `numpy.linalg.solve`; for CHNO it reproduces today's closed form to
  roundoff.
- **Raise instead of ignoring.** An element in the substrate or biomass
  that is neither C, H, O nor a key of the mapping; a source that lacks its
  element; a source carrying an element that is not balanced; a singular
  system.
- **Charged sources.** An optional `h_plus=` `Species`, required when a
  source is charged, adds H+ as an unknown and a charge equation.

## Open questions

1. Whether a non-balanced element in the biomass raises or warns. Raising
   makes the three stirred-tank tutorials that use `Yeast` with `"CHO"`
   fail until they switch to `Yeast_CHO` or balance nitrogen.
2. Whether `h_plus=` is enough for charged sources, or a general
   "charge-balancing species" argument is wanted.
3. Whether `monod_aerobic_growth` and the template take the mapping under
   the same name.

## Relationship to other notes

- [EXPLICIT_SPECIES_RESOLUTION.md](EXPLICIT_SPECIES_RESOLUTION.md): its
  checkpoint 9 resolves `n_source_id` against the model's species and
  leaves the argument's shape to this note.
