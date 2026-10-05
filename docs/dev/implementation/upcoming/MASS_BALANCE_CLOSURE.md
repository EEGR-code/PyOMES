# Mass Balance Closure — Design Note

> **Status:** design settled, not yet started. No branch, no checklist, no
> code.
>
> - **Written:** 2026-10-01.
> - **Audited against the code and measured:** 2026-10-02, main at `cf1c572`.
> - **Decisions taken with the repo owner:** 2026-10-02 to 2026-10-05.
> - **Follows:** `explicit-species-resolution`, which gives every id in
>   `n_mol` a `Species` with atoms; closure is checked element by element over
>   those species.
> - **Phases:** three, each with its own branch and checklist (see
>   [Phases](#phases)).

## The standard

A model conserves every element and its charge to roundoff across everything
the package does to it, with what enters and leaves through boundaries
accounted for, regardless of ChemicalEquilibriumEngine choice (e.g. Bisection or the NR
engine). The conservation monitor reports only genuine imbalances, and names
their source.

One rule sits behind the decisions below: **the package does not repair a
specification behind the user's back.** It does not add or remove material to
make a liquid neutral, to keep an amount non-negative, or to make a volume
fit. Where a specification is unphysical or a model does not close, it warns,
says by how much, and runs what it was given.

What the standard does not promise:

- **Models that do not close by their own definition.** BSM2's and ADM1's
  kinetics follow the published benchmarks, whose stoichiometry is on a COD
  basis; the stirred-tank template's growth ignores biomass nitrogen under
  `"CHO"` ([GROWTH_STOICHIOMETRY.md](GROWTH_STOICHIOMETRY.md)). They warn and
  run.
- **A model that selects the electroneutral closure** (see
  [Decisions](#decisions), "pH closure"). Hydrogen is then not conserved, by
  the model's own choice, and the monitor reports the amount.
- **The PHREEQC engine in a control volume.** It is a reference calculator
  for direct solves, not a core engine (`OPEN_WORK.md`, "The PHREEQC engine is
  a reference calculator, not a control-volume engine").

Acceptance test: the two strict expected failures for H and O in
`tests/standalone/test_user_defined_model.py` pass, joined by the same check
for each engine and for runs with feeds, vents and dosing, at 1e-12 relative.

## Where closure fails today

Checked against the code on 2026-10-02. Items 1 to 4 are the note's original
list, corrected; 5 to 10 were found during the audit.

1. **The engines never book water.** An equilibrium that consumes or forms
   water does not take it from or add it to `n_mol`. The Bisection engine
   neither reads nor writes `H2O` (it is not in `_written_species`;
   `_SOLVENT_IDS` only skips water when classifying reactants). The NR engine
   overwrites liquid water with 55.509 mol/L × `V_L` on every solve
   (`species_mol_L["H2O"] = _C_WATER_MOL_L` in `engines/nr/engine.py`), which
   also erases water that kinetic reactions made.
2. **The engines force the liquid neutral.** Both find pH from a charge
   balance (`solve_from_equilibrium_set` in `engines/bisection/acid_base.py`;
   the charge row in `engines/nr/solver.py`) and ignore the H+ and OH- they
   are given. A liquid with a net charge of Z mol is made neutral by removing
   (Z > 0) or adding (Z < 0) Z mol of hydrogen. Every declared equilibrium is
   charge-balanced, so no redistribution of species can change a liquid's net
   charge; the engine changes it by changing what the liquid contains.
3. **The monitor does not net out boundary flows, and is weaker than that.**
   It compares totals between steps, so vented gas or fed solute reads as
   drift. It runs only under `SequentialAdvanceSolver` (the single
   `check_step` call site). It is inert when the reaction system is a bare
   `KineticReaction` or absent, because its species registry is filled only
   when the reaction system accepts a conservation monitor; the fed-batch
   tutorial's tank has an empty registry and raises nothing. Its baseline is
   taken at the end of step 1, so the initial solve and first step are never
   checked. Its per-step threshold, 1e-8 × max(total, 1), is 1.1e-6 mol of H
   on the acceptance model, above the engine gap of 1.4e-7 mol per step.
4. **The liquid volume never changes.** Nothing in `PyOMES/` or `models/`
   writes a liquid's `V_L` after construction; `LiquidFeed` returns C × Q and
   no volume.
5. **Most liquids hold no solvent water.** The stirred-tank template builds
   its liquid from dissolved gases only: the 1600 L fed-batch tank ends with
   16 mol of water (the growth reaction's), where about 88,800 mol fills it.
   The D2C script holds 0.006 mol, BSM2 none. ADM1 seeds its water.
6. **Amounts are kept non-negative by default, which creates mass.** The step
   solvers' default clamp scales each overdrawn species on its own, and
   several floors say nothing (`OPEN_WORK.md`, "Amounts are kept non-negative
   by default, without the model asking for it").
7. **The NR engine's folded gas-liquid rows lose mass.** The solve splits the
   component between gas and liquid, but `EquilibriumResult.apply_to_phases`
   writes only the liquid share, and the control volume drops the transfer
   model for that species because the engine owns it.
8. **NR precipitation is not stored.** Mineral amounts are reported in the
   result; the dissolved species are written back at their reduced amounts.
   Read in the code, not measured.
9. **NR conserves component totals only to its residual tolerance** (1e-10
   mol/L), where the Bisection engine conserves them to roundoff.
10. **The Bisection solver can return a non-neutral state without saying so.**
    With no water reaction declared and too little acid or base to absorb the
    charge, no neutral pH exists; `solve_pH` returns the pH with the smallest
    residual.

## Measurements

Scratch scripts, 2026-10-02 to 2026-10-05, not committed. Each case runs in a
few seconds.

**User-defined model (1 h in 20 steps), relative drift.**

| Engine | C | H | O |
|---|---|---|---|
| Bisection | -1.7e-16 | +2.5e-8 | +2.5e-8 |
| NR | +1.6e-8 | -2.1e-5 | -2.1e-5 |

Bisection's gap is exactly 2 H and 1 O per HCO3- formed. NR's is the water
reset; its C drift is the residual tolerance.

**Dosing 1e-3 mol** (both engines): Na+ + OH- moves pH from 3.3885 to 3.8664
and loses 2.0e-3 mol H and 1.0e-3 mol O; OH- alone and H+ alone leave pH
unchanged and lose 1.0e-3 mol H.

**Where each model's change comes from** (a ledger of engine write-back,
boundary flow, dosing and remainder):

| Case | Finding |
|---|---|
| ADM1, 0.2 h in 8 steps | Every `ConservationWarning` is vent flow (N 1.911e-4 mol in one step). After netting, no step is over threshold. The first solve creates 0.112 mol of H and O to balance Na+ and Cl-. |
| BSM2, 100 × 0.01 h | The first solve creates 17.7 mol of H for the `S_an` lump. Kinetics do not close: C +0.25 %, H +0.79 %, O +1.0 %, N -0.17 % per hour, in every step. Charge closes to 7e-15 mol. |
| D2C script, 1 h in 200 steps | C, N, P and Na close to 1e-16 after netting. The engine removes 3.924e-2 mol H and 1.962e-2 mol O, exactly 2 H and 1 O per NaOH dosed, in the 58 dosing steps. |
| Fed-batch configuration, 2 h in 100 steps | Monitor silent. C and H close. N +1.743 mol is the yeast nitrogen. O +1.741 mol is the clamp (0.870 mol O2 in 39 steps). Volume stays 1600 L with 10 L fed. |
| NR with a folded Henry row | Total C falls from 1.0e-2 to 2.0e-4 mol in 5 steps. |

**What a net charge costs.** Five 1 L liquids starting with ±4 to 5 mmol of
net charge, solved once with water booked: each ends neutral, hydrogen changes
by exactly the starting charge with the opposite sign (1.008 mg per mmol),
and every other element is unchanged. Acetic acid 10 mmol with Na+ 5 mmol
gives pH 4.76 and loses 5.04 mg; the same liquid specified as sodium acetate
plus acetic acid, or as acetic acid plus NaOH, gives the same pH with no
change in mass.

**The two closures compared** (prototypes wrapped around the engines):

| | pH against today | Bare H+ dose | Start with 5 mmol unbalanced Na+ |
|---|---|---|---|
| Force neutrality, book water | identical | ignored | pH 4.76, 5 mmol H removed |
| Conserve net charge, book water | 2e-15 (Bisection), 1.3e-9 (NR) | pH 3.39 → 2.94, closed | pH 3.39, charge stays, closed |

With the liquid first made neutral, the conserving closure reproduces today's
results: BSM2 to 1.4e-14 in pH (3.5e-14 in every species), D2C to 4.5e-13,
ADM1 to 4e-8. Seeded as they are today, BSM2 goes from pH 3.30 to 4.85 and
ADM1 from 12.64 to 6.89.

Not measured: the batch template, any notebook, the control volume in
`test_iron_oxidation.py`, Bisection states where `n_active` is below the
number of pKas, and how far results move under a derived liquid volume
(estimated at 2e-4 relative for the template tank).

## Decisions

**pH closure.** The solve conserves the liquid's net charge: the charges after
it sum to what they were before it. This is the default, and it is the same
one-dimensional search the Bisection engine does today with one constant
changed, so that engine needs no rewrite. Explicit hydrogen and oxygen totals
as extra unknowns are not needed: with no redox, the oxygen balance only fixes
the water amount, and the hydrogen balance then follows from charge and the
component totals. After the write-back, liquid water is set so the liquid's
oxygen total is unchanged.

An **electroneutral closure** (the charges sum to zero) is kept as a setting a
model selects explicitly. BSM2 and ADM1 are published with pH defined that
way, over totals and ion lumps, and their builders select it. The monitor
reports the hydrogen it adds or removes. `ReactionSystem`'s existing `solver`
argument uses `"charge_balance"` to mean the Bisection engine, so the new
setting needs a different name.

**Charged specifications.** A liquid specified with a net charge describes
something that cannot exist. The control volume warns when it first advances,
and a feed is checked when the simulation is built; the message states the net
charge and that the fix is to correct the specification or add a counter-ion
of the user's choice, real or a lump. Nothing is corrected automatically and
there is no reconciliation tool; a liquid's or a feed's net charge is readable
by the user. `equilibrate_to_pH` and the totals call (below) remain the ways
to build a neutral state deliberately.

**The direct engine call.** `engine.solve(totals=..., strong_ions=...)` keeps
its meaning: the pH at which those totals are neutral. Totals do not say which
forms the material is in, so there is no input charge to conserve.
`engine.solve(phases=...)` follows the closure above. The validation suites
use the first form and are untouched. Distinct names for the two are logged in
`OPEN_WORK.md`.

**Water.** A liquid that tracks its solvent holds its water in `n_mol`,
whichever way its volume is given, and every stream into or out of it carries
its water explicitly.

**Liquid volume.** Every liquid has a volume model, which gives the volume
its contents occupy. The first model is ideal mixing: water's molar volume
from its density at the phase's temperature, solutes from a table of molar
volumes passed to the model (zero if not listed). A liquid is built in one of
two ways, and they differ only in what `V_L` returns:

- with a **stated volume**: `V_L` is that number and stays constant. This is
  today's form and meaning, and the assumption the benchmark digesters make
  (outflow equals inflow, constant density, reaction water neglected);
- **derived**: `V_L` is the volume model's value for the current contents,
  and is read-only.

The two must not become separate code paths. In both, the amount of water
present is calculated from the underlying volume model: the make-up utility
below inverts it to find the water a given volume contains. So a
stated-volume liquid whose reactions name water holds the water its volume
implies, the engines book water the same way for both, and the only branch is
in `V_L`. A stated-volume liquid still knows what its contents would occupy,
so a warning can state the difference. A liquid in a model with no solvent
species holds no water. To be checked in the volume phase: BSM2 and ADM1
carrying their solvent water under a stated volume, against BSM2's sentinels.
If keeping the two forms on one path proves awkward in the code, that is
raised with the owner, not worked around.

The model classes belong in `thermo/liquid/`, which may import only `units`,
so they take plain data. Streams add moles, never litres, so a different
volume model (apparent molar volumes, a density correlation) can be swapped in
later without changing anything else.

A **make-up utility** turns solutes, a solvent, a target volume and a
temperature into amounts. Solutes may be given as amounts (mol), masses (g),
molar concentrations (mol/L) or mass concentrations (g/L), in any mix. The
same calculation with a flow in place of a volume builds a stream. It takes
species in the model; splitting a salt into its ions stays with the
`OPEN_WORK.md` entry on the recipe layer.

A feed or drain with a net volumetric flow attached to a stated-volume liquid
warns that the volume will not change and names the derived form.

**Gas volume.** Stays as constructed. A vessel type that owns the total
volume, the shrinking headspace and a liquid overflow outlet are logged
(`OPEN_WORK.md`, "Container volume, and liquids that keep a stated volume").

**Water activity.** Taken as 1 in this work, as the engines do today. That is
an approximation which worsens with ionic strength, and it is not consistent
with a non-ideal activity model: ion activity coefficients below 1 imply a
water activity below 1. It is left out because closure is about amounts (where
every atom is), while activity moves where an equilibrium sits, and because
adding it would re-baseline every model that uses a non-ideal activity model,
BSM2 included. Logged in `OPEN_WORK.md`, "Water's activity is taken as 1
whatever the activity model".

**Clamp.** No clamp by default: an amount may go negative, with a warning, and
a clamp is something the user passes to the solver. The silent floors get the
same treatment.

**Tolerances.** Closure tests use 1e-12 relative. The monitor's threshold
stays at 1e-8 of the element total per step until the engine phase ships, then
is set from that phase's measurements, scaled on the inventory excluding
water.

**Models that do not close.** They warn and run. A reaction whose
stoichiometry does not balance its elements warns when the control volume is
built, naming the reaction and the elements, so the run-time imbalance has a
stated cause. BSM2's and ADM1's docstrings say they implement the published
benchmarks.

**Other engine failures.** The NR gas write-back (item 7) and NR
precipitation (item 8) are fixed in the engine phase; precipitation starts
with an audit, since
[NR_PRECIPITATION_CV_INTEGRATION.md](NR_PRECIPITATION_CV_INTEGRATION.md)
predates the current code.

## Phases

Three phases, each on its own branch with its own checklist and
`<name>-shipped` tag. The checkpoints below are provisional; each phase's
checklist is written when it starts, and the first phase's measurements inform
the other two.

### 1. `conservation-ledger`

Leaves a monitor that reports only genuine imbalances, with their source,
under every step solver. Results unchanged except in the last checkpoint.

1. Closure tests per engine (plain, with a feed, a vent, a dose) as strict
   expected failures carrying today's numbers.
2. The control volume records the element and charge flux of every external
   change, after any clamp: boundaries, source terms, controller doses,
   inter-CV links. The engine write-back and the clamp are recorded as named
   sources. The record is kept automatically: the user declares no source or
   sink objects for accounting's sake. Each boundary also keeps a cumulative
   total of what it applied, for auditing (how much left through a vent). This
   is the accounting that [RESERVOIR_TYPE.md](RESERVOIR_TYPE.md) §3.1 asks
   for, without its `Reservoir` type or `FlowBoundary` protocol; the
   per-boundary total is kept as its own small record so a later `Reservoir`
   can adopt it. The existing `ExternalFluxRecord` and `LinkFlowRecord` are
   reused where they fit.
3. The monitor compares the change in totals with the recorded flux; its
   baseline is taken before the first solve.
4. The monitor is active for every control volume and every step solver.
5. A reaction that does not balance its elements warns at definition. How
   BSM2's reactions pass the existing balance check is looked at first.
6. No clamp by default, and the silent floors off or warned. Which runs clamp
   today and what each rate law does with a negative concentration are checked
   first; results are measured before and after.
7. Close-out.

### 2. `variable-liquid-volume`

Leaves liquids that dilute, concentrate and change volume with what enters
and leaves, for the models that select it. Unblocks
[DOSING_AGENTS.md](DOSING_AGENTS.md).

1. A prototype measures how far results move under a derived volume on the
   template tank and its tutorials. The owner sees the numbers before anything
   is committed.
2. The recorder keeps the liquid volume at each step.
3. Every `LiquidPhase` has a volume model (ideal mixing first). A stated
   volume holds `V_L` constant; otherwise `V_L` is the model's value for the
   contents. `V_L` is read-only and the public setter goes.
4. The make-up utility and the stream helper.
5. The template tank and its tutorials carry their solvent water and derive
   their volume; feeds and drains carry water. The comment in
   `docs/tutorials/templates/fed_batch_fermenter.py` and the `OPEN_WORK.md`
   entry "The liquid volume never changes" are updated here.
6. The step and system solvers under a changing volume, with a guard on
   `LiquidDrain` as the volume approaches zero. Checked against an inert
   tracer, C0·V0/(V0 + Q·t).
7. The warning for a net volumetric flow on a stated-volume liquid.
8. Close-out.

### 3. `engine-closure`

Leaves the Bisection and NR engines conserving every element and the charge
in a control volume, with pH unchanged within a stated tolerance for models
that start neutral.

1. NR stops overwriting liquid water.
2. Water is booked from the oxygen balance at write-back.
3. NR's write-back conserves component totals exactly.
4. The conserving closure in both engines; the electroneutral closure as an
   explicit setting, selected by the BSM2 and ADM1 builders. The adaptive step
   solver integrates H+ under the conserving closure.
5. The charged-specification warning and the readable net charge; the repo's
   other models and tests specified neutral.
6. The Bisection solver reports when no neutral pH exists (item 10).
7. NR folded gas rows are written back to the gas phase.
8. NR precipitation: audit, then stored.
9. The monitor's threshold set from this phase's measurements.
10. Close-out: the H and O markers removed from
    `test_user_defined_model.py`, engine docstrings, `OPEN_WORK.md` entries,
    the "engine fit" paragraph of `DOSING_AGENTS.md`.

## Verification

- Strict closure tests per engine (Bisection, NR) on a user-defined model,
  with and without boundaries, at 1e-12 relative.
- A fingerprint before and after every checkpoint that changes results: the
  user-defined model under each engine, the D2C script (equilibria and
  dosing), ADM1, BSM2 and the fed-batch configuration. Each records what
  moved, by how much, and why.
- BSM2's sentinels unchanged: it keeps a stated volume and the electroneutral
  closure, so its solve is the one it has today.
- The buffer-standard validation suites unchanged: they use the totals call.

## Relationship to other notes

- [EXPLICIT_SPECIES_RESOLUTION.md](../shipped/EXPLICIT_SPECIES_RESOLUTION.md): every id
  in `n_mol` resolves to a `Species` with atoms, which element closure needs.
- [DOSING_AGENTS.md](DOSING_AGENTS.md): solution dosing needs the derived
  volume; its second open question (is a solution's water listed or implied)
  is answered here: listed, with the make-up calculation as the convenience.
  Under the conserving closure a dose of H+ or OH- is honoured as given.
- [STRONG_ION_INFERENCE_GENERALIZATION.md](STRONG_ION_INFERENCE_GENERALIZATION.md):
  under the conserving closure a control volume's pH no longer depends on the
  fixed strong-ion tables; the totals call and the electroneutral closure
  still use them.
- [GROWTH_STOICHIOMETRY.md](GROWTH_STOICHIOMETRY.md): the template's biomass
  nitrogen, which the monitor will report.
- [NR_PRECIPITATION_CV_INTEGRATION.md](NR_PRECIPITATION_CV_INTEGRATION.md):
  taken into the engine phase.
- `OPEN_WORK.md`: the default clamp and floors; the two meanings of
  `engine.solve()`; container volume and stated-volume liquids; the PHREEQC
  engine.
