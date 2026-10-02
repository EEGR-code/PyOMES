# Dosing Agents — Design Note

> Status: design note, not yet started. No branch, no checklist, no code.
> Written 2026-10-01 while settling how a pH corrector names what it adds.
> Depends on the liquid volume becoming variable (OPEN_WORK, "The liquid
> volume never changes").

## Where dosing stands

A corrector (`ControlVolume.equilibrate_to_pH`, `PHController`'s acid and
base) is a composition: the species, and moles of each, that one mole of
reagent adds to the liquid, e.g. `{"Na+": 1, "OH-": 1}` for NaOH. A plain id
is one mole of that species. Every id must be among the model's species, and
a dose whose net charge is not zero warns. What is dosed is exactly what the
composition says; nothing checks it against a reagent formula, because
hydrates, solutions and deliberate lumps would all disagree with one.

Three things that are not covered:

- **Solutions.** A 1 M NaOH stock is a composition per litre of solution
  (`{"Na+": 1.0, "OH-": 1.0, "H2O": 55.3}` mol/L) plus the volume it adds.
  The composition is easy; the volume is not, because nothing in PyOMES
  changes a liquid phase's `V_L` after construction (feeds and doses add
  moles only).
- **Dissolution.** A reagent dosed as a solid (Ca(OH)2, CaCO3, or NaOH
  pellets) that dissolves through a declared equilibrium, e.g. a Ksp for
  `Ca(OH)2(s) <-> Ca++ + 2 OH-`. The dose then goes to a solid phase, and the
  model's dissolution equilibrium moves it into the liquid.
- **Engine fit.** A declaration the engine cannot honour does something else
  without a word. All three engines solve pH from a charge balance and do not
  take H+ or OH- as input (Bisection and NR overwrite both; PHREEQC sets pH by
  charge), so dosed OH- only acts through its counter-ion: the pH is right,
  but its H and O leave the books ([MASS_BALANCE_CLOSURE.md](MASS_BALANCE_CLOSURE.md)).
  A counter-ion moves pH only if the engine recognises it: through a fixed
  strong-ion id table in Bisection and NR, through the user's
  `component_map` in PHREEQC. The Bisection engine has no dissolution or
  precipitation; NR supports precipitation.

## Direction

- **Solution dosing.** A dosing agent given per litre of solution, with the
  volume it adds; `equilibrate_to_pH` finds a volume of stock rather than
  moles of reagent, and `PHController` acts on a volumetric dose rate. Needs
  the variable liquid volume first.
- **Dissolution route.** The corrector doses into a named phase (a solid
  phase for an undissolved reagent); the model's declared dissolution
  equilibrium does the rest.
- **Engine checks.** When a CV is built (or an engine attached), warn when its
  declarations and the engine do not fit: a dissolution or precipitation
  equilibrium with an engine that has none; dosed H+ / OH- with an engine
  that overwrites them. Each warning names the declaration and the engine.

## Open questions

1. How a dosing agent is written: one object carrying composition, basis
   (per mole of reagent or per litre of solution) and phase, or keyword
   arguments on each entry point.
2. Whether a solution's water is listed in its composition or implied by its
   volume (and how the two are kept consistent).
3. Whether the engine checks warn or raise, and whether they run at
   construction or at the first solve.

## Relationship to other notes

- [EXPLICIT_SPECIES_RESOLUTION.md](../shipped/EXPLICIT_SPECIES_RESOLUTION.md): its
  checkpoint 13 makes correctors compositions of the model's species with a
  charge-neutrality warning.
- [PHCONTROLLER_CORRECTOR_VALIDATION.md](PHCONTROLLER_CORRECTOR_VALIDATION.md):
  checking that a controller's dose can move pH.
- [STRONG_ION_INFERENCE_GENERALIZATION.md](STRONG_ION_INFERENCE_GENERALIZATION.md):
  the Bisection and NR engines recognise strong ions through a fixed id table.
- [MASS_BALANCE_CLOSURE.md](MASS_BALANCE_CLOSURE.md): booking the H and O
  of dosed H+ / OH- and the water they form; the variable liquid volume.
