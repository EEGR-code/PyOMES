# Liquid Volume Models — Design Note

> **Status:** design discussion. No branch, no checklist, no code.
>
> - **Written:** 2026-11-03, main at `333a62b`.
> - **Refines:** the "Liquid volume" decision and phase 2
>   (`variable-liquid-volume`) of
>   [MASS_BALANCE_CLOSURE.md](MASS_BALANCE_CLOSURE.md). Nothing settled there
>   is reopened; this note says what a volume model is, where it lives, and
>   how `V_L` behaves in the solvers.
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
other so the pairing is visible, and does not enforce it (open question 2).

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
  one place; a liquid phase gets it when bound to its control volume
  (open question 1). `ThermoFramework.gas_eos` is held but read by nothing
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

  *Considered and not chosen:* adding `V_L` to the state vector with
  `dV/dt = Σ v̄ᵢ·dnᵢ/dt`. That is exact in principle, but the integrated volume
  drifts away from the model's value for the same amounts by the integrator's
  tolerance, and it holds the same information twice. Computing it from the
  amounts is consistent by construction.
- **System solver and links:** a link that moves liquid refreshes both
  control volumes.
- **Recorder:** keeps `V_L` at each step (closure note, phase 2 checkpoint 2).

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
  linear rate, with the guard near zero volume from the closure note
  (phase 2 checkpoint 6).
- **Doses.** A solution dose carries moles, water included
  ([DOSING_AGENTS.md](DOSING_AGENTS.md)); its volume follows the same way.

### Pressure

The protocol takes `P_atm` so a pressure-dependent model fits later. Liquids
are close to incompressible at the pressures PyOMES runs at, and the ideal
and apparent-molar models ignore it. Where the value comes from (the gas
phase's pressure or the standard pressure) is open question 3.

## Decisions

- **2026-11-03.** Volume models are a separate thermodynamic modelling area
  alongside the activity models and gas equations of state, not a property
  calculator. The volume model is the single source of truth for volume, and
  density is derived from it.
- **2026-11-03.** `V_L` stays a stored state on the phase, refreshed from the
  volume model after each change to its contents, including within the
  adaptive solver's steps.

## Open questions

1. **Where the model is held.** On `ThermoFramework` (one per control volume,
   recommended) or on each `LiquidPhase` (as the closure note's wording
   "every `LiquidPhase` has a volume model" could also be read), or the
   framework's model as a default that a phase can override.
2. **Pairing warning.** Whether a non-ideal activity model with the ideal
   volume model (or the reverse) warns, or is only documented.
3. **Pressure.** The gas phase's pressure, the standard pressure, or a value
   passed by the model.
4. **Second model.** The apparent molar volume is proposed second because
   dosing concentrated reagents needs it; confirm against the models the
   owner expects to run.
5. **Refresh granularity.** Whether `SequentialAdvanceSolver` refreshes after
   each sub-step (proposed) or once per step; measured in the phase 2
   prototype.

## Relationship to other notes

- [MASS_BALANCE_CLOSURE.md](MASS_BALANCE_CLOSURE.md): the settled decisions
  this builds on (water held in `n_mol`, stated and derived volumes on one
  code path, streams in moles, the make-up utility). This note details phase 2
  checkpoints 3, 4 and 6.
- [DOSING_AGENTS.md](DOSING_AGENTS.md): solution dosing depends on the derived
  volume.
- `OPEN_WORK.md`, "Container volume, and liquids that keep a stated volume":
  the shrinking headspace and an overflow outlet stay out of scope here.
