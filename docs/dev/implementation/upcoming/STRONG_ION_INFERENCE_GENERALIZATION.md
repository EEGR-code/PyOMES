# Strong-Ion Inference Generalization — Design Note

> Status: design note, not yet started. No branch, no checklist, no code
> yet. Surfaced during a conversational review of
> `docs/tutorials/ArXiv_preprint/01_predict_ph_simple_liquid.ipynb`
> (2026-09-18), while explaining why `NRChemicalEquilibriumEngine.solve()`
> requires `strong_ions=` as a separate argument from `totals=`. See
> "How to start one" below when picked up.

## Motivation

`01_predict_ph_simple_liquid.ipynb` calls
`engine.solve(totals={...}, strong_ions={"CT_K": CT_P, "CT_Cl": CT_N})` —
the *direct* call pattern. K⁺ and Cl⁻ have to be passed separately from
`totals` because they never appear in any declared `EquilibriumReaction`
(they don't react), so nothing in the reaction network encodes their
presence or amount.

While explaining this, a natural question came up: the engine also has a
*phase-based* call pattern (`solve(phases={"liquid": liquid_phase})`,
`NRChemicalEquilibriumEngine._read_from_phases`,
`PyOMES/chemical_equilibrium/nr_engine.py:530-582`) that takes one
`n_mol` dict of *everything present* and derives both `totals` (by
summing `n_mol` over each tableau component's species) and `strong_ions`
for you. That's a strictly nicer API — you just state the recipe's
composition once, in one dict, and let the engine sort out which species
are reacting vs. spectating. So: why isn't *that* the only interface, and
why does `strong_ions=` need to exist as a concept the caller manages at
all?

Investigating found that the phase-based path's strong-ion derivation is
*not actually derived from the reaction network* — it's a hardcoded,
closed lookup:

```python
# PyOMES/chemical_equilibrium/nr_engine.py:56-75
_STRONG_ION_SPECIES_TO_KEY: Dict[str, str] = {
    "K+": "CT_K", "Na+": "CT_Na", "Cl-": "CT_Cl", "NO3-": "CT_NO3",
    "Mg++": "CT_Mg", "Ca++": "CT_Ca", "Zn++": "CT_Zn", "Mn++": "CT_Mn",
    "Cu++": "CT_Cu", "Co++": "CT_Co", "Fe2+": "CT_Fe2",
    "Mo7O24------": "CT_Mo7O24", "MoO4--": "CT_MoO4",
    "S_cat": "CT_cation", "S_an": "CT_anion",
}
```

`_read_from_phases` only recognizes a species as a strong ion if its id
is literally one of these ~15 entries (`nr_engine.py:576-580`:
`if sp_id in n_mol` gated by iterating this dict, not by iterating
`n_mol`). A species present in `n_mol` that (a) isn't part of any
declared reaction's component *and* (b) isn't on this list silently
contributes nothing to the charge balance — not an error, not a warning,
just dropped.

Compounding this: the paired charge table (`_STRONG_CHARGES`, needed to
turn each `CT_*` key back into a signed charge for the charge-balance
residual and the ionic-strength sum) is a **second** hardcoded dict that
duplicates information the `Species` object already carries on `.charge`
— and it exists as **three independently-written copies**, currently in
sync but with nothing enforcing that:

- `PyOMES/chemical_equilibrium/nr_engine.py:86-92`
- `PyOMES/chemical_equilibrium/nr_solver.py:321-327` (inside `_ionic_strength`)
- `PyOMES/chemical_equilibrium/nr_solver.py:452-458` (inside `solve_nr`)

Every time a new charged strong ion needs to be supported (the most
recent additions look like `Mg++`/`Ca++`/`Zn++`/`Mn++`/`Cu++`/`Co++`/
`Fe2+`/molybdate, presumably added for iron/precipitation work — see
`tests/validation/speciation/07_iron_oxidation.ipynb`/`08_iron_oxidation_and_precipitation.ipynb`),
it has to be added to `_STRONG_ION_SPECIES_TO_KEY` *and* to all three
copies of `_STRONG_CHARGES`, by hand, in sync — even though every one of
those species already declares its own charge in `common_species.py`.

## What "strong ion" actually means, structurally

A strong ion, in this codebase's own terms, is simply: **a species with
nonzero `n_mol`, on the liquid phase, that is not a member of any
tableau component** (i.e., doesn't appear in any declared
`EquilibriumReaction`'s stoichiometry). That's a structural property —
"is this species in the reaction network" — fully determined by the
`NRTableau` already built in `from_reactions()`, combined with a
property (`.charge`) the `Species` object already carries. Nothing about
"is this a strong ion" requires a separately maintained allowlist — the
allowlist is currently doing the job that `tableau component membership`
+ `Species.charge` could already do generically.

**"Strong ion" is inherently an aqueous-phase concept, and the
generalization must preserve that scoping explicitly.** It already holds
true today, if incidentally: `_read_from_phases` splits `n_mol` (from
`phases["liquid"]`) from `gas_n_mol` (from `phases["gas"]`) at the top of
the method, and the strong-ion loop only ever reads `n_mol` — a
gas-phase entry with the same species id is never considered. That's
correct — a gas molecule has no charge-balance role, and a solid-phase
species (once precipitation/`SolidPhase` writeback lands, per
`NR_PRECIPITATION_CV_INTEGRATION.md`) is accounted for via its own ξ
(extent precipitated), not as a dissolved charge carrier. The
generalized rule should state this as an explicit invariant — "only
`phases['liquid'].n_mol` is ever scanned for strong-ion inference" — not
leave it as an implicit side effect of how the dict happens to be built,
so a future edit to `_read_from_phases` (e.g. adding a third phase type)
can't accidentally start picking up gas/solid entries as if they were
dissolved ions.

## Proposed generalization

Replace `_STRONG_ION_SPECIES_TO_KEY` + the three `_STRONG_CHARGES`
copies with a single derivation: for any species id present in
`phases["liquid"].n_mol` (phase-based path) that is not a member of any
`NRTableau` component, treat it as a strong ion automatically, using
`Species.charge` from its own definition for the charge-balance/ionic-
strength contribution. Gas-phase and (once relevant) solid-phase `n_mol`
entries are never scanned for this, by explicit rule, matching the
scoping above. This would:

- Remove three hand-maintained, independently-drifting tables in favor
  of one structural rule plus data already on `Species`.
- Fix the silent-drop gap: any charged species not currently on the
  allowlist (e.g. a future `Li+`, `Br-`, or any locally-declared ion in
  a model file that isn't in `common_species`) would be picked up
  automatically instead of needing an engine-code change first.
- Keep the *direct* call pattern (`solve(totals=..., strong_ions=...)`)
  unaffected in the near term — this note is scoped to how the
  phase-based path derives `strong_ions` from `n_mol`, not to redesigning
  the direct kwarg API (see Open question 3).

## Open questions

1. **`CT_*` key naming.** Downstream code (`strong_ions` dict keys,
   `_STRONG_CHARGES` lookups, any caller passing `strong_ions=` directly
   by name) is written in terms of `CT_K`/`CT_Cl`/etc., not raw species
   ids. Does anything outside `nr_engine.py`/`nr_solver.py` depend on
   these specific key strings, or could the generalized path key
   `strong_ions` by species id directly (`"K+"` instead of `"CT_K"`)?
   Needs a usage audit before deciding whether `CT_*` naming survives as
   a public convention or becomes purely internal.
2. **`S_cat`/`S_an` generic buckets.** These two existing keys
   (`nr_engine.py:73-74`) don't correspond to one real species — they
   look like a generic "some cation"/"some anion" escape hatch for
   models that don't want to name a specific strong ion. Confirm their
   actual call sites and decide how (or whether) they fit the
   species-id-derived scheme, since they have no matching `Species`
   object to read `.charge` from.
3. **Direct call pattern.** Should `solve(totals=..., strong_ions=...)`
   eventually be superseded by a single `solve(n_mol=...)`-shaped direct
   call (mirroring the phase-based path, but without requiring a
   `LiquidPhase` object), so scripting/testing callers get the same
   auto-derivation benefit? Out of scope for this note's minimal fix,
   but worth deciding since it would change the notebook-authoring
   guidance already given in `docs/tutorials/ArXiv_preprint/
   01_predict_ph_simple_liquid.ipynb`.
4. **Neutral, non-reacting species in `n_mol`.** A species with
   `charge=0` that's present in liquid `n_mol` but not in any reaction
   (e.g. a biological state variable transported with the liquid)
   contributes nothing under the generalized rule — correct silently, or
   should the engine warn that a "known but chemically inert here"
   species was seen, in case it was meant to be declared into a reaction
   and wasn't? Related to `PHCONTROLLER_CORRECTOR_VALIDATION.md`'s theme
   of catching silently-inert configuration.

## Trigger conditions

Pick this up when either happens:

- A model needs a strong ion not on the current allowlist and hits the
  silent-drop gap in practice (a real correctness bug, not just a code
  smell).
- Another charged species is added to `common_species.py` for
  precipitation/complexation work (per
  `MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md`) and someone
  has to manually update three hardcoded tables again — a good forcing
  function to do the generalization instead of adding a fourth entry to
  each.

## How to start one

Per this folder's usual convention
([README.md](README.md#how-to-start-one)): resolve the open questions
above, write a checklist file
(`STRONG_ION_INFERENCE_GENERALIZATION_CHECKLIST.md`), cut a branch off
`main` (suggested name: `strong-ion-inference-generalization`), and add a
"Currently in flight" pointer to this folder's `README.md`. Small and
self-contained — no dependency on any other in-flight phase, though it
touches the same `nr_engine.py`/`nr_solver.py` files as the NR
Precipitation and `MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md`
tracks, so re-check for merge overlap if either is in flight when this is
picked up.
