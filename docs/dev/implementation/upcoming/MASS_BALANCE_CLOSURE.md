# Mass Balance Closure — Design Note

> Status: design note, not yet started. No branch, no checklist, no code.
> Written 2026-10-01. Follows `explicit-species-resolution`, which gives every
> id in `n_mol` a `Species` with atoms; closure is checked element by element
> over those species.

## The standard

Every model conserves every element and charge to roundoff, with what enters
and leaves through boundaries accounted for, whichever speciation engine it
uses. The conservation monitor reports only genuine imbalances. Acceptance
test: the two strict expected failures for H and O in
`tests/standalone/test_user_defined_model.py` pass, and are joined by the same
check for each engine and for runs with feeds, vents and dosing.

## Where closure fails today

1. **The engines treat water as an unlimited solvent.** An equilibrium that
   consumes or forms water (CO2 + H2O <-> HCO3- + H+; OH- + H+ -> H2O) does not
   take it from or add it to `n_mol`. The Bisection engine keeps H2O fixed
   (`_SOLVENT_IDS`); the NR engine writes its own H2O.
2. **The engines do not take H+ or OH- as input.** All three solve pH from a
   charge balance: H+ and OH- (and, in NR, H2O) are written back from the
   solution, so H+ or OH- that a feed or dose adds is overwritten and its H
   and O leave the books. Dosing `{"Na+": 1, "OH-": 1}` gives the right pH,
   but the OH- and the water it forms are not accounted for.
3. **The monitor does not net out boundary flows.** Gas leaving through a vent
   or solute entering with a feed reads as drift.
4. **The liquid volume never changes.** Feeds and doses add moles but no
   volume, so the water balance and every concentration in a fed-batch run
   are off (`docs/tutorials/templates/fed_batch_fermenter.py`: 1600 L at the
   start and at the end with 5 L/h of feed).

pH and dissolved concentrations are not affected by 1 and 2 at dilute
conditions (10 mmol/L of water formed is about 0.02 % of the solvent); element
totals are.

## Direction

- **Water as a tracked participant.** Equilibria that consume or form H2O
  debit or credit it in `n_mol`, in every engine. Water's activity stays
  about 1 at dilute conditions; its amount is booked.
- **Proton and element totals as engine input.** The engines take the
  liquid's totals of H and O (or a proton balance: total H relative to a
  reference state of the components, as PHREEQC does internally) from
  `n_mol`, including any dosed or fed H+ and OH-, and solve speciation
  subject to those totals and charge balance, instead of treating H+ and OH-
  as outputs only.
- **Boundary-aware monitoring.** Each boundary reports the element and charge
  flux it applied; the monitor compares the change in totals with the net
  boundary flux.
- **Variable liquid volume.** Feeds, drains and doses carry a volume (and its
  water); the phase's `V_L` follows. Shared with
  [DOSING_AGENTS.md](DOSING_AGENTS.md).

## Verification

- Strict closure tests per engine (Bisection, NR, PHREEQC) on a user-defined
  model, with and without boundaries.
- Before/after fingerprints of the stirred-tank template, ADM1 and BSM2: pH
  and dissolved concentrations unchanged within a stated tolerance; water
  amounts and H/O totals move, and the movement is recorded.
- BSM2 sentinels and the buffer-standard validation suites unchanged within
  their tolerances.

## Open questions

1. The engines' formulation: a proton balance (PHREEQC style) or explicit H
   and O totals, and whether the Bisection engine's one-dimensional search
   survives it.
2. Whether water's activity departs from 1 at high ionic strength (it is
   already an activity-model question), or stays 1 with only its amount
   booked.
3. The order of the four parts: engines first (1, 2), or the monitor (3)
   first so the gaps are measured before they are closed.
4. How the variable volume is defined: additive volumes of what enters, or
   from mass and density.

## Relationship to other notes

- [EXPLICIT_SPECIES_RESOLUTION.md](EXPLICIT_SPECIES_RESOLUTION.md): every id
  in `n_mol` resolves to a `Species` with atoms, which element closure needs.
- [DOSING_AGENTS.md](DOSING_AGENTS.md): solution dosing needs the variable
  volume; dose compositions are the inputs closure books.
- [STRONG_ION_INFERENCE_GENERALIZATION.md](STRONG_ION_INFERENCE_GENERALIZATION.md):
  the engines' fixed strong-ion id tables.
