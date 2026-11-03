# Liquid Volume Models — Design Note

> **Status:** design discussion. No branch, no checklist, no code.
>
> - **Written:** 2026-11-03, main at `333a62b`.
> - **Owns:** phase 2 (`variable-liquid-volume`) of
>   [MASS_BALANCE_CLOSURE.md](MASS_BALANCE_CLOSURE.md), whose scope and
>   checkpoints are here ([Scope](#scope-and-checkpoints)). The water and
>   liquid-volume decisions in that note are settled ground and are not
>   reopened; this note says what a volume model is, where it lives, how
>   `V_L` behaves in the solvers, and what the phase does.
> - **Order:** after `conservation-ledger` (phase 1), before
>   `engine-closure` (phase 3); not in parallel with either.
> - **Decisions taken with the repo owner:** 2026-11-03 (see
>   [Decisions](#decisions)).
> - **Measured:** nothing yet.

## The problem

A liquid's `V_L` is set at construction and never written again
(`OPEN_WORK.md`, "The liquid volume never changes"). Feeds return C × Q in
mol/h and no volume, `LiquidPhase.apply_flux` changes moles only, and the
adaptive solver reads `V_L` once per step (`liq_V_L` in
`core/solvers.py`). A fed-batch tank never dilutes, and every concentration
it reports is on its starting volume.

The obvious fix, `dV/dt = Σ Q`, adds stream volumes. That is ideal mixing:
it is wrong whenever the volume a material occupies in the liquid differs
from the volume it had in its stream. Strong acids and bases diluted into
water, ethanol and water, and concentrated salt solutions all contract or
expand on mixing, by up to a few percent. In a fed-batch run that error
accumulates, and every concentration inherits it.

The closure note already rules that out: streams add moles, never litres,
and a liquid's volume comes from a volume model of its contents. This note
fills in the model.

## What a volume model is

Volume is a thermodynamic property of a phase's contents, in the same way
as activity coefficients and gas fugacities:

- `V = (∂G/∂P)` at constant T and amounts.
- Each species' partial molar volume is `v̄ᵢ = (∂μᵢ/∂P)` at constant T and
  amounts, and `V = Σ nᵢ·v̄ᵢ`.
- The mixing (excess) volume is `Vᴱ = V − Σ nᵢ·vᵢ°`, the pressure derivative
  of the excess Gibbs energy `Gᴱ` whose amount derivatives give the activity
  coefficients.

So the activity model and the volume model are two derivatives of one
quantity. The established models:

| Model | Basis | Typical use |
|---|---|---|
| Ideal mixing (Amagat) | `V = Σ nᵢ·vᵢ°`, no excess volume | Default, dilute solutions |
| Apparent molar volume (Masson, Redlich–Meyer) | Solute `v̄` varies with ionic strength, with the Debye–Hückel limiting slope | Electrolytes at moderate concentration |
| Redlich–Kister excess volume | Polynomial `Vᴱ` in mole fraction | Non-electrolyte mixtures |
| Pitzer volumetric | Pressure derivative of the Pitzer `Gᴱ` | Concentrated brines; PHREEQC's solution density uses a variant |
| HKF | Standard partial molar volumes of aqueous ions vs T and P | Elevated T and P |
| Liquid equations of state (Peng–Robinson with volume translation, PC-SAFT) | Full equation of state | Organic liquids, high pressure |
| Kell, IAPWS-95 | Pure water | The solvent's own molar volume |

`thermo/liquid/water_properties.py` already has a Kell-type water density,
which the ideal model uses for water.

### Density is the same information

`ρ = Σ nᵢ·MWᵢ / V`. A density correlation is a volume model written the other
way round, so the two approaches do not clash. What would clash is two
sources of truth: a density property calculator and a volume model that
disagree. The volume model is the single source of truth; a density, where
anything needs one, is derived from it.

Volume does not go through the property-calculator hook. Property calculators
run after speciation and write to `phase.properties` for rate laws, but
volume is needed earlier, and by everything that turns moles into
concentrations: the engines, transfer models, links, drains and the
recorder.

### Pairing with the activity model

Ideally the volume model is the volumetric counterpart of the activity model,
for example Pitzer activity with Pitzer volumetric terms. Davies and SIT have
no volumetric form of their own; the matching dilute-solution volume model is
the apparent molar volume with the Debye–Hückel slope. At the pressures
PyOMES runs at, pairing models from different sources is an accepted
approximation, not an error. The framework holds both models next to each
other so the pairing is visible, and does not enforce it (open question 1).

## Where it lives

A modelling area of its own, parallel to the gas equations of state in
`thermo/gas/` (`GasEOS`, `IdealGasEOS`, `PengRobinsonEOS`):

- **Package:** `thermo/liquid/volume/`. The layering rules in
  `tests/standalone/test_package_layering.py` allow it: the subfolder is
  inside `thermo.liquid`, so it may import `water_properties`, absolutely
  (`from PyOMES.thermo.liquid.water_properties import ...`), since relative
  imports of two dots are not allowed. `thermo/` imports only `units`, so the
  models take plain data: per-species parameter tables keyed by id, and the
  solvent's id, passed in by the caller. No fixed `"H2O"` id inside the
  models (compare [GAS_SPECIES_IN_CORE.md](GAS_SPECIES_IN_CORE.md)).
- **Protocol:**

  ```python
  @runtime_checkable
  class LiquidVolumeModel(Protocol):
      name: str

      def volume_L(self, n_mol: Mapping[str, float], *,
                   T_K: float, P_atm: float) -> float:
          """Volume (L) the given amounts occupy at T and P."""

      def partial_molar_volumes_L_per_mol(
              self, n_mol: Mapping[str, float], *,
              T_K: float, P_atm: float) -> Dict[str, float]:
          """v̄ᵢ for each species in n_mol."""
  ```

  `volume_L` must be homogeneous of degree one in the amounts (doubling every
  amount doubles the volume). Every model in the table above is. The
  make-up utility (closure note) inverts `volume_L` for the solvent amount
  that gives a target volume; the volume rises monotonically with the
  solvent amount, so a one-dimensional root search works for any model, and
  a model may provide a closed form.
- **Implementations, in order:**
  1. **Ideal mixing:** the solvent from water's density at T, solutes from a
     table of molar volumes, zero if not listed. As the closure note decided.
  2. **Apparent molar volume:** for electrolytes; needed for the dosing work
     (concentrated NaOH, acids).
  3. **Redlich–Kister excess volume:** for non-electrolyte mixtures, if a
     model needs one.
  4. **Pitzer volumetric:** with a Pitzer activity model, which does not exist
     yet. Liquid equations of state are out of scope.
- **Held by:** `ThermoFramework`, as `liquid_volume` next to
  `liquid_activity` and `gas_eos`, so a control volume's conventions live in
  one place; a liquid phase gets it when bound to its control volume. The
  per-species parameters come from the `Species` declarations, kept to
  intrinsic data (see [Decisions](#decisions)). `ThermoFramework.gas_eos` is held but read by nothing
  today. The volume model is wired in the same phase that adds it, so it does
  not become a second unused field.

## `V_L` in the solvers

**Owner decision (2026-11-03):** `V_L` stays a stored state on the liquid
phase. It is updated from the volume model as the contents change, not
integrated from flows.

- **Stated volume:** unchanged from today; `V_L` is the given number
  (closure note).
- **Derived volume:** the phase stores `V_L` and refreshes it from the model
  whenever the package changes its amounts or temperature: after
  `apply_flux`, after an engine's write-back (which books water once
  `engine-closure` ships), and when T is set. The public setter goes for this
  form (closure note). Storing the value, rather than computing it on every
  read, keeps the cost to one model call per change and gives each step one
  defined volume.
- **`SequentialAdvanceSolver`:** refreshes after each sub-step that applies
  moles (feeds, boundaries, reactions), so each sub-step sees the volume its
  predecessor left.
- **`SimultaneousEulerSolver`:** computes every delta from the start-of-step
  snapshot, including its volume, and refreshes once after applying them.
- **`SimultaneousAdaptiveSolver`:** the state vector stays amounts only. Each
  right-hand-side evaluation computes the volume from the trial amounts with
  the model and uses it for concentrations, the engine and drain dilution;
  after an accepted step the stored `V_L` is set from the accepted amounts.
  This replaces the `liq_V_L` cached once per step.
- **`MonolithicODESolver`:** has its own right-hand side (reactions and links
  only today), so it computes the volume from the trial amounts in the same
  way.

  *Considered and not chosen:* adding `V_L` to the state vector with
  `dV/dt = Σ v̄ᵢ·dnᵢ/dt`. That is exact in principle, but the integrated volume
  drifts away from the model's value for the same amounts by the integrator's
  tolerance, and it holds the same information twice. Computing it from the
  amounts is consistent by construction.
- **System solver and links:** a link that moves liquid refreshes both
  control volumes.
- **Recorder:** keeps `V_L` at each step (checkpoint 2).

### Streams

- **Feeds.** `Q` is the feed's flow at its own conditions. The stream helper
  (closure note) turns `Q` and the feed's composition into mol/h of every
  species, water included, using the feed's volume model. The tank's volume
  then changes by `V(n + Δn) − V(n)`, which differs from `Q·Δt` when mixing
  is non-ideal. That difference is the effect being modelled.
- **Drains.** `Q` is at tank conditions; a drain removes the fraction
  `Q·Δt/V` of every species. Because `volume_L` is homogeneous of degree one,
  the volume falls by exactly `Q·Δt` under any model, so drains need no
  correction for mixing. The exponential step-averaged form in `LiquidDrain`
  assumes a constant `V`; under a shrinking volume it is checked against the
  linear rate, with the guard near zero volume (checkpoint 6).
- **Doses.** A solution dose carries moles, water included
  ([DOSING_AGENTS.md](DOSING_AGENTS.md)); its volume follows the same way.

### Pressure

The protocol takes `P_atm` so a pressure-dependent model fits later. Liquids
are close to incompressible at the pressures PyOMES runs at, and the ideal
and apparent-molar models ignore it. Where the value comes from (the gas
phase's pressure or the standard pressure) is open question 2.

## Scope and checkpoints

Provisional, as in the closure note; the checklist is written when the phase
starts. Checkpoints 1 to 8 are the closure note's former phase 2 list, in its
order, with additions; 9 and 10 moved here from elsewhere in that note.

1. A prototype measures how far results move under a derived volume on the
   template tank and its tutorials. The owner sees the numbers before anything
   is committed.
2. The recorder keeps the liquid volume at each step, with a volume
   diagnostic per step: `V_L` before and after, the nominal volume each
   stream brought (its `Q·Δt` at its own conditions) and the remainder (the
   change from mixing, temperature and reactions). Volume is not conserved,
   so the remainder is reported, never warned on. One check it allows: under
   ideal mixing at constant temperature with no reaction, the remainder is
   zero to roundoff. Phase 1's ledger is not changed for this.
3. The `thermo/liquid/volume/` package with the protocol and the ideal-mixing
   model; every liquid has a volume model. A stated volume holds `V_L`
   constant; otherwise `V_L` is the model's value for the contents, written
   only by the phase's refresh.
4. The make-up utility and the stream helper. The stream helper uses the
   receiving control volume's volume model at the feed's temperature, unless
   a model is passed to the feed. A liquid with a derived volume also gets an
   explicit method that makes it up to a target volume (`make_up_to(V)` or
   similar): it adds or removes solvent water by inverting the model, and the
   ledger records the change as a named source. It is a method, not a `V_L`
   setter, so no write to `V_L` changes the contents silently.
5. The template tank and its tutorials carry their solvent water and derive
   their volume; feeds and drains carry water. BSM2 and ADM1 carry their
   solvent water under a stated volume, checked against BSM2's sentinels (the
   ledger checklist's checkpoint 2e defers BSM2's `H2O` exemption to here).
   If BSM2's sentinels move, work stops and the owner sees the numbers; the
   sentinels are not adjusted to fit. `ProportionalGasOutlet` removes from the
   liquid the water it calculates leaving as vapour
   (`last_water_loss_mol_per_h`, set and read by nothing today, found by the
   ledger checklist's audit): once liquids hold their water, the missing
   debit would be both an unexplained H and O loss and a missing volume sink.
   The comment in `docs/tutorials/templates/fed_batch_fermenter.py` and the
   `OPEN_WORK.md` entry "The liquid volume never changes" are updated here.
6. Every solver path under a changing volume (`SequentialAdvanceSolver`,
   `SimultaneousEulerSolver`, `SimultaneousAdaptiveSolver`,
   `MonolithicODESolver`, the system solver and links), with a guard on
   `LiquidDrain` as the volume approaches zero. Checked against an inert
   tracer, C0·V0/(V0 + Q·t).
7. The warning for a net volumetric flow on a stated-volume liquid.
8. Close-out.
9. **NR stops overwriting liquid water** (formerly the closure note's phase 3
   checkpoint 1). The NR engine sets liquid water to 55.509 mol/L × `V_L` on
   every solve (`species_mol_L["H2O"] = _C_WATER_MOL_L` in
   `engines/nr/engine.py`). Under a derived volume that feeds back: the
   volume sets the water and the water sets the volume. 55.509 mol/L is water
   at 1000 g/L, more than water's molar concentration at 25 °C (about
   55.34 mol/L), and solutes add their own volume on top, so the volume grows
   on every solve with no feed at all. The reset also erases water a feed
   brought in. Lands before checkpoint 5 derives any volume under NR.
10. **The monitor's threshold scaled on turnover, on trial** (replaces the
    closure note's solvent-excluded threshold). Once liquids hold their water
    (about 88,800 mol in the template tank), a threshold of 1e-8 of the
    element total is about 1.8e-3 mol of H per step, four orders of magnitude
    above the engine gap of about 1.4e-7 mol per step, and the monitor is
    nearly blind to H and O until phase 3. The problem is scaling on the
    inventory: a large inert stock makes the threshold coarse. The trial
    replaces it with

    `threshold = max(k·ε·inventory, rtol·turnover)`

    where `k·ε·inventory` is the floating-point floor (ε is machine
    precision; about 1e-10 mol for the template tank, water included) and
    `turnover` is the gross amount of the element that moved in the step:
    stream flows, reaction extents and the engine write-back. No species is
    treated specially, and the monitor stays sensitive to water itself,
    which phase 3's water booking needs. Lands with or before checkpoint 5.
    The reason for the form is stated in the code where the threshold is
    computed, not only here.

    *Tests.* Realistic cases, each run before and after: a genuine imbalance
    the size of the engine gap is caught in a water-holding tank (the
    template fed-batch tank under each engine); a closed run with large
    inert water raises nothing; a step where almost nothing moves does not
    warn on roundoff; vented and fed material nets out; BSM2's and ADM1's
    warnings are recorded and explained. If the trial does not hold up, the
    solvent-excluded scale is the fallback and the owner decides.

    *Kept open.* Checks or reports that leave the solvent out (an inventory
    excluding water, or the balance of everything but water) may still be
    wanted to find problems later; nothing in this checkpoint rules them out.

## Decisions

- **2026-11-03.** Volume models are a separate thermodynamic modelling area
  alongside the activity models and gas equations of state, not a property
  calculator. The volume model is the single source of truth for volume, and
  density is derived from it.
- **2026-11-03.** `V_L` stays a stored state on the phase, refreshed from the
  volume model after each change to its contents, including within the
  adaptive solver's steps.
- **2026-11-03.** This note owns the `variable-liquid-volume` phase. The
  phases run in order: `conservation-ledger`, then this one, then
  `engine-closure`. The NR water reset and the monitor's threshold move into
  this phase (checkpoints 9 and 10), and the threshold's reason is documented
  in the code.

The rest were settled in the review of this note against the closure note
and the ledger checklist:

- **2026-11-03. The monitor's threshold** is trialled on turnover rather than
  the solvent-excluded inventory, with tests on realistic cases
  (checkpoint 10). Water-excluding checks stay possible for later
  diagnosis.
- **2026-11-03. Who holds the model.** The model choice is held by
  `ThermoFramework`, one per control volume, next to the activity model. Its
  per-species parameters come from the `Species` declarations, gathered by
  core and passed to the model as plain tables. `Species` carries only
  intrinsic, single-species data of wide use (a standard molar volume).
  Model-specific parameter sets (apparent-molar-volume coefficients,
  Redlich–Kister or Pitzer interaction terms, many of them per pair of
  species) are passed to the model, not added to `Species`, so the species
  declaration does not grow into a home for optional thermodynamic models. A
  per-phase override is not added: no control volume has two liquid phases,
  and the activity model makes the same assumption; liquid–liquid systems
  would need both changed together.
- **2026-11-03. A feed's volume model** is the receiving control volume's,
  at the feed's temperature, unless one is passed to the feed (checkpoint 4).
- **2026-11-03. Evaporation.** `ProportionalGasOutlet` debits the liquid for
  the water vapour it calculates leaving (checkpoint 5). Its effect on the gas
  headspace (vapour raising the pressure in a fixed volume, or the gas volume
  in a constant-pressure vessel) is later work, logged in `OPEN_WORK.md`,
  "Container volume, and liquids that keep a stated volume".
- **2026-11-03. BSM2 and ADM1** carry their solvent water under a stated
  volume in this phase; if BSM2's sentinels move, the owner sees the numbers
  (checkpoint 5).
- **2026-11-03. Writing `V_L`.** For a derived volume, the phase's refresh is
  the only writer, which reconciles the closure note's "read-only" with
  stored state. Making a liquid up to a volume is an explicit method that
  changes the water and is recorded by the ledger (checkpoint 4). A stated
  volume stays settable, as today.
- **2026-11-03. Volume in the ledger.** A volume diagnostic is recorded in
  this phase, reported and not checked (checkpoint 2); phase 1's ledger stays
  moles only.
- **2026-11-03. Out of scope.** The gas headspace, which needs a vessel
  concept (a rigid vessel's gas volume is the vessel's less the liquid's; the
  equation of state then gives the pressure) and follows this phase, and
  water activity. Both are logged in `OPEN_WORK.md`.

## Open questions

1. **Pairing warning.** Whether a non-ideal activity model with the ideal
   volume model (or the reverse) warns, or is only documented.
2. **Pressure.** The gas phase's pressure, the standard pressure, or a value
   passed by the model.
3. **Second model.** The apparent molar volume is proposed second because
   dosing concentrated reagents needs it; confirm against the models the
   owner expects to run.
4. **Refresh granularity.** Whether `SequentialAdvanceSolver` refreshes after
   each sub-step (proposed) or once per step; measured in the prototype
   (checkpoint 1).

## Relationship to other notes

- [MASS_BALANCE_CLOSURE.md](MASS_BALANCE_CLOSURE.md): the standard and the
  settled decisions this builds on (water held in `n_mol`, stated and derived
  volumes on one code path, streams in moles, the make-up utility). This note
  owns its phase 2. Phase 1 (`conservation-ledger`, with its checklist) comes
  first; phase 3 (`engine-closure`) follows, books water into the liquids
  this phase gives water, and calls this phase's volume refresh after its
  write-back.
- [DOSING_AGENTS.md](DOSING_AGENTS.md): solution dosing needs this phase for
  the volume a dose adds, and phase 3 for a dose of H+ or OH- to be honoured
  as given.
- `OPEN_WORK.md`, "Container volume, and liquids that keep a stated volume":
  the shrinking headspace and an overflow outlet stay out of scope here.
