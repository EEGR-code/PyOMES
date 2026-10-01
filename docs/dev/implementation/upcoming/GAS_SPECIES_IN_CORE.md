# Gas Species in Core — Design Note

> Status: design note, not yet started. No branch, no checklist, no code.
> Written 2026-10-01 during `explicit-species-resolution` checkpoint 8, which
> removed the stirred-tank template's fixed gas ids (its gas phase now comes
> from the model's declarations). The same ids are assumed in core, outside
> the template; this note collects those sites.

## The rule

The same rule as [EXPLICIT_SPECIES_RESOLUTION.md](EXPLICIT_SPECIES_RESOLUTION.md):
**a model only knows the species and reactions it was given.** A model whose
oxygen is `Species(id="O2_aq", ...)` should work everywhere a model whose
oxygen is `"O2"` does. Today the template does; the core sites below do not.

## Sites that assume a gas id

| Site | What it does now | What breaks with other ids |
|---|---|---|
| `core/snapshot.py` (`sensors["DO_mol_L"]`) | Dissolved oxygen is `n_liq_mol["O2"] / V_liq_L`, only if `"O2"` is present | `DOAgitationController` / `DOCascadeController` read a missing sensor; DO control does nothing |
| `core/boundaries.py` `GasFeed.__init__` | `y or {"O2": 0.21, "N2": 0.79}`: a feed with no composition feeds air | A feed created without a composition invents O2 and N2 the model never declared |
| `core/boundaries.py` `MembraneGasBoundary._DEFAULT_ATMOSPHERE` | External atmosphere defaults to O2 / CO2 / N2 at dry-air fractions | Same: gases the model never declared, and no flux for differently named ones |
| `core/gas_liquid_link.py` `KineticGasLiquidLink.speciation_keys` | Defaults to `{"CO2": "CO2"}` (molecular CO2 for the alpha correction) | A model naming CO2 otherwise gets no alpha correction unless it overrides the mapping |
| `core/gas_liquid_link.py` `set_kLa_with_co2_ratio(kLa_O2, co2_ratio)` | Writes `kLa["O2"]` and `kLa["CO2"]` | Writes kLa for ids the model may not have |
| `control/cv_loops.py` `DOAgitationController` / `DOCascadeController` | Adjust `kLa.O2` and couple `kLa.CO2` by `kLa_CO2_ratio`; param paths hard-code `kLa.O2` / `kLa.CO2`; the cascade adjusts the feed's O2 fraction (`yO2_*`) | Control acts on ids the model may not have |
| `control/cv_loops.py` `_DEFAULT_GAMMA` / `_DEFAULT_MW_KG` (`_nozzle_gamma`, `_nozzle_mw`) | Heat-capacity ratio and molar mass by gas name (O2, N2, CO2, H2O, "Air") for the nozzle / vent physics | Another id gets γ = 1.35 and the molar mass from `chemicals` if passed, otherwise air's 28.97 g/mol |

## Direction

- **Defaults go where they invent gases.** `GasFeed` with no composition and
  `MembraneGasBoundary` with no external atmosphere raise, as the template's
  `GasFeedConfig` now does; air is passed explicitly (e.g. `AIR` from the
  bioprocess database).
- **Ids the model chooses are passed in.** The DO sensor and the DO
  controllers take the oxygen species' id (and the CO2 one, for the coupled
  kLa); the link's alpha mapping is derived from the model's declared
  cross-phase reactions only (as `derive_speciation_keys` already does), with
  no `{"CO2": "CO2"}` default; `set_kLa_with_co2_ratio` takes the two ids or
  goes.
- **Property tables keyed by name become species data or explicit input.** The
  molar mass is already on `Species`; the heat-capacity ratio needs a home
  (an argument to the vent physics, or a property on the gas species).

## Open questions

1. Whether the DO sensor stays a built-in sensor parametrised by an id, or
   becomes something the DO controllers compute from the species they are
   given (so there is no oxygen-specific sensor in `snapshot.py`).
2. Whether `set_kLa_with_co2_ratio` survives with two id arguments, or goes in
   favour of setting each kLa through `params_changed`.
3. Where the heat-capacity ratio lives (vent argument, `Species` property, or a
   database table the model passes).
4. Whether the template's growth defaults (`o2_id="O2"` etc. on
   `OrganismConfig`) should also lose their defaults once the core sites are
   explicit.

## Relationship to other notes

- [EXPLICIT_SPECIES_RESOLUTION.md](EXPLICIT_SPECIES_RESOLUTION.md): same rule;
  its checkpoint 8 made the stirred-tank template's gas phase come from the
  model's declarations.
- [PHCONTROLLER_CORRECTOR_VALIDATION.md](PHCONTROLLER_CORRECTOR_VALIDATION.md):
  the other controller-side species-id question.
