# Strong-Ion Inference Generalization — Design Note

> Status: design note, not yet started. No branch, no checklist, no code
> yet. Surfaced during a conversational review of
> `docs/tutorials/ArXiv_preprint/01_predict_ph_simple_liquid.ipynb`
> (2026-09-18), while explaining why `NRChemicalEquilibriumEngine.solve()`
> requires `strong_ions=` as a separate argument from `totals=`. Scope
> decisions below (de-duplicate now, generalize `NRChemicalEquilibriumEngine`
> only, leave `BisectionChemicalEquilibriumEngine` alone for now) were made
> in the same conversation, same day. See "How to start one" below when
> picked up.

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

## Scope: `NRChemicalEquilibriumEngine` only, for now

The duplication described above isn't confined to `nr_engine.py`/
`nr_solver.py`. `PyOMES/chemical_equilibrium/engine.py` (the older,
bisection-based `BisectionChemicalEquilibriumEngine`) carries its **own**,
independently-written copy of the species→`CT_*` allowlist
(`engine.py:747-758`), and `acid_base.py` goes further — `CT_K`, `CT_Na`,
`CT_Cl`, `CT_cation`, `CT_anion`, etc. are literal **named function
parameters** on its `solve`-shaped functions, not just dict keys.

**This note deliberately does not touch the bisection engine.** Checked
during scoping: `ReactionSystem`'s default is `solver="charge_balance"`
(`PyOMES/reactions/reaction_system.py:110`), which builds
`BisectionChemicalEquilibriumEngine` — meaning it is currently **the
default production solver** for every `ControlVolume`/`ReactionSystem`
model, including `StirredTankBuilder` (BSM2/ADM1 templates), unless a
model explicitly opts into `solver="newton_raphson"`.
`NRChemicalEquilibriumEngine` is the newer engine, used explicitly today
in the standalone tutorial notebooks (this one included) and the
precipitation/`LAYER1_GAP_CLOSURE` work, but it has not been promoted to
the default.

A future transition to `NRChemicalEquilibriumEngine` as the preferred/
default solver — which would make the bisection engine's own copy of
this same duplication moot rather than needing its own fix — is plausible
and worth keeping in mind, but is its own, separately-scoped phase (default-
solver swap across every CV-based production model, plus BSM2/ADM1
golden-trajectory re-validation, in the style of past `chemistry-unification`/
`layer1-gap-closure` solver-numerics changes). Not decided, not scheduled,
not part of this note.

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

## Why the silent-drop gap isn't already caught

Reasonable to expect some existing check would flag this — the ingredients
are already in the codebase. Checked directly, and it isn't caught, for
two separable reasons.

`EquilibriumResult.charge_residual` (`nr_solver.py:573`, literally
`R_final[-1]`, the last row of the Newton residual at convergence) is
*purely the solver's own self-consistency check* on the charge-balance
equation it was actually given. It's near-zero by construction whenever
Newton converges, regardless of whether a strong ion was silently
dropped beforehand — the solver never sees the dropped ion, so its own
equation is satisfied by definition either way. It's informational only:
`01_predict_ph_simple_liquid.ipynb` just prints it ("should be ~0") for a
human to eyeball. `AccuracyMonitor.check_charge_residual`
(`monitoring/accuracy.py:218`) checks the same solved-for value for
drift after a kinetic step — "did the constraint the solver knows about
stay satisfied," not whether that constraint was complete to begin with.
Neither is a ground-truth check.

`ConservationMonitor._charge_residual` (`monitoring/conservation.py`) *is*
independent, and closer to a real ground-truth check: it sums `z × n`
over every species in `n_mol` using each species' own `.charge`, via a
registry that `ControlVolume._collect_species_registry()`
(`control_volume.py:840-878`) explicitly supplements with any
`common_species`-cataloged species found in `n_mol` — the method's own
docstring says this exists specifically so "Cl⁻, Na⁺, K⁺... which never
appear in reaction stoichiometry" still get included "for charge
conservation accounting." So this check would see the *true* charge sum
including an off-allowlist strong ion, independent of
`_STRONG_ION_SPECIES_TO_KEY` entirely.

It still wouldn't catch this bug, because it's wired as a **drift-over-
time check against a per-run baseline**, not an absolute "does this
equal zero" check. A silently-dropped strong ion present at a roughly
constant amount from the first step onward produces a *wrong but
constant* charge offset — baked into the baseline at step 1, matched by
every later step, so drift stays ~zero forever and nothing ever warns.
It's built to catch a conservation law being violated *during a run* (a
real bug that changes total charge over time), not a structurally
incomplete input that's wrong from the very first step in a stable way.

And separately: none of this monitoring is wired up for the *direct*
`solve(totals=..., strong_ions=...)` call pattern
`01_predict_ph_simple_liquid.ipynb` actually uses.
`ConservationMonitor` only attaches via `ControlVolume.__init__`/
`ReactionSystem.attach_conservation_monitor` — the notebook calls the
bare `NRChemicalEquilibriumEngine.solve()` directly, no `ControlVolume`
involved, so there is zero independent cross-check running in that
context regardless of the drift-vs-absolute distinction above.

**This strengthens the case for Phase 1**, and suggests a cheap,
independently-useful interim option worth considering alongside it: a
one-shot *static* check (not drift-based) comparing the true charge sum
(`Σ z·n` over everything in `phases["liquid"].n_mol`, via each species'
own `.charge`) against the charge the solve actually consumed (masters +
secondaries + recognized `strong_ions`), and warning on any nonzero
difference — this would catch exactly the silent-drop case even on the
very first call, unlike the drift-based `ConservationMonitor` check, and
would work for the *direct* call pattern too, not just CV-integrated
runs.

## Phase 0: de-duplicate immediately (decided, low-risk)

Independent of whether/when the full generalization below happens: pull
the three `_STRONG_CHARGES` copies (`nr_engine.py:86-92`,
`nr_solver.py:321-327`, `nr_solver.py:452-458`) into a single canonical
definition — e.g. keep it once in `nr_engine.py` next to
`_STRONG_ION_SPECIES_TO_KEY` and have both `nr_solver.py` call sites
import it instead of retyping it. Pure DRY refactor, no behavior change,
no open questions, ships on its own ahead of (or instead of, if the
bigger generalization stalls) Phase 1. This alone removes the "three
tables can silently drift apart" risk.

## Phase 1: derive strong-ion status structurally

Replace `_STRONG_ION_SPECIES_TO_KEY` + the now-unified `_STRONG_CHARGES`
with a single derivation: for any species id present in
`phases["liquid"].n_mol` (phase-based path) that is not a member of any
`NRTableau` component, treat it as a strong ion automatically, using
`Species.charge` from its own definition for the charge-balance/ionic-
strength contribution. Gas-phase and (once relevant) solid-phase `n_mol`
entries are never scanned for this, by explicit rule, matching the
scoping above. This would:

- Remove the remaining hand-maintained table in favor of one structural
  rule plus data already on `Species`.
- Fix the silent-drop gap: any charged species not currently on the
  allowlist (e.g. a future `Li+`, `Br-`, or any locally-declared ion in
  a model file that isn't in `common_species`) would be picked up
  automatically instead of needing an engine-code change first.
- Extend the same "just declare what's present, no backend naming
  required" ergonomics to the *direct* call pattern too — decided
  direction, see Open question 3 for the remaining shape question.

## Open questions

1. **`CT_*` key naming.** Downstream code (`strong_ions` dict keys,
   `_STRONG_CHARGES` lookups, any caller passing `strong_ions=` directly
   by name) is written in terms of `CT_K`/`CT_Cl`/etc., not raw species
   ids. Does anything outside `nr_engine.py`/`nr_solver.py` depend on
   these specific key strings, or could the generalized path key
   `strong_ions` by species id directly (`"K+"` instead of `"CT_K"`)?
   Needs a usage audit before deciding whether `CT_*` naming survives as
   a public convention or becomes purely internal.
2. **`S_cat`/`S_an` generic buckets — resolved, not actually open.**
   `nr_engine.py:73-74` maps `S_cat`/`S_an` to `CT_cation`/`CT_anion`, but
   neither is declared anywhere as a concrete `Species`. That's the
   actual gap, not a conceptual one:
   `PyOMES/monitoring/conservation.py`'s own charge-balance accounting
   (`conservation.py:104-111`) already anticipates `S_cat`/`S_an` being
   ordinary registry entries — "looked up from the engine's strong-ion
   mapping... `atoms={}`, only charge contributes" — i.e. exactly a
   normal `Species(id="S_cat", atoms={}, charge=+1)` object, no different
   in kind from `K+`/`Cl-`. **Fix: add concrete `Species` declarations
   for the generic lumps** (or let a model declare its own local
   equivalent under whatever id it wants), and Phase 1's structural rule
   ("not in any reaction → strong ion, read `.charge`") picks them up
   automatically — no special case needed.
   `PyOMES/chemistry/thermo_params.py`'s
   `compute_CT_cation_from_charge_balance` (solving the separate,
   upstream problem of *how much* generic charge carrier is needed to
   hit a target pH) is unaffected by this — its output number feeds into
   `n_mol={"S_cat": ...}` exactly as it feeds `strong_ions=
   {"CT_cation": ...}` today.
3. **Direct call pattern — decided direction, shape still open.** Yes:
   `solve(totals=..., strong_ions=...)` should eventually gain the same
   "just declare what's present" ergonomics as the phase-based path, so a
   scripting/testing/notebook caller never has to know the `CT_*`
   convention exists (mirroring what `_read_from_phases` already does
   today — it takes real species ids like `"K+"` in `n_mol`, translating
   to `CT_K` only internally). Still open: exact call shape (a new
   `solve(n_mol=..., V_liq_L=...)` direct-pattern kwarg? reusing
   `totals=` for this?). Would also update the notebook-authoring
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

## Adjacent, out of scope: the salt-recipe maps (`chem_recipe.py`/`registry.py`/`recipe.py`)

Scoping also turned up `chem_recipe.py`'s ~30-compound `ChemSpec`
registry and `registry.py`'s `SALT_DISSOCIATION_MAP`, both mapping a salt
name (`"NaCl"`, `"KH2PO4"`, ...) to its constituent `CT_*` ions — a third
and fourth place `CT_*` strings appear, alongside `recipe.py`'s
`SolutionRecipe` convenience layer built on top of them. Not part of
Phase 0/Phase 1, but the same "hand-maintained table duplicating
something the framework can already express" pattern shows up here too,
worth recording rather than re-discovering later.

A salt's dissociation into ions doesn't have to live in a bespoke
name→ions table — it's expressible either of two ways already native to
the framework:

- **Implicitly**, by how the user populates the phase's contents in the
  first place — stating `n_mol={"Na+": ..., "Cl-": ...}` directly, never
  naming "NaCl" to the framework at all. This is what Phase 1 already
  assumes.
- **As a declared solid-liquid equilibrium reaction** (`NaCl(s) ⇌ Na⁺ +
  Cl⁻`), using the salt's real, measured solubility product — large but
  not fabricated; these are all genuinely highly soluble salts, so the
  true constant already predicts ~complete dissolution at the mmol/L-to-
  low-mol/L doses these recipes use. Mechanically expressible; **but not
  actually a good reason to prefer it here**, since it drags in the full
  precipitation machinery (a tracked `SolidPhase`, an extra tableau
  unknown, the active-set solve in
  `MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md`) for a salt
  that never meaningfully re-precipitates in this regime, and risks the
  kind of large-`log_K` numerical conditioning issue this notebook's own
  §6c PHREEQC discussion already ran into with extreme `-gamma` values —
  real solver cost for zero behavioral difference from just stating the
  ion is present.

Either way, `chem_recipe.py`/`registry.py`'s hardcoded salt→`CT_*` maps
are a third, parallel way of expressing information the framework can
already represent — not something structurally required. Their actual
value is narrower: skipping the step of writing those ion amounts out by
hand for ~30 common lab compounds, given a mass in grams. Whether that
convenience is worth keeping as a separate hardcoded registry is an open
question for a *separate* note if ever pursued — not resolved here, and
not blocking Phase 0/Phase 1 either way.

**Update 2026-09-20 (`chemical-equilibrium-engines-subfolder`, Part B):**
`chemical_equilibrium/strong_ions.py` — the only consumer of
`SALT_DISSOCIATION_MAP`, via `strong_ions_from_feed_molL(feed)` — has been
removed, along with its test and its package-level export. It had no
production callers (the engines take `strong_ions=` as a plain dict), and it
silently dropped anything outside its fixed 11-key `CT_*` set. As a result
`SALT_DISSOCIATION_MAP` (still defined in `chemistry/registry.py` and exported
from `PyOMES.chemistry`) now has **no consumer in the repo**, so the open
question above now covers it too, alongside `chem_recipe.py`'s `ChemSpec`
registry: keep, merge or remove is still undecided and still a separate note.
If BioSTEAM coupling ever needs neutral-salt expansion again, it belongs in
`PyOMES/stream_adapter.py`, emitting species ids (not `CT_*` keys) and warning
on unmapped species — not in the equilibrium package.

**Update 2026-09-20 (same phase, checkpoint 9b):** `chemical_equilibrium/api.py`
(`SpeciationEngineAdapter`) and `factory.py` (`SpeciationFactory`) were also
removed — Bisection-only, no callers, restorable from commit `2e5554a`. That
leaves a second orphaned cluster in the same recipe layer, deliberately **not**
touched by that phase:

- `chemistry/types.py`: `AqueousEquilibrium` is now unused anywhere;
  `AqueousTotalsUser.to_engine()` has no caller; `AqueousTotals` and
  `AqueousTotalsUser` are reachable only through `SolutionRecipe.to_totals_user()`
  / `to_totals()` in `chemistry/recipe.py`, which nothing outside `recipe.py`
  calls. All three are still exported from `PyOMES.chemistry`.
- Stale wording: `chemistry/recipe.py`'s docstring says its totals are "for use
  with the standalone speciation interface", and `chemistry/registry.py`'s
  docstring lists `AqueousTotalsUser` as a consumer; the interface they refer
  to no longer exists.

Together with `SALT_DISSOCIATION_MAP` and `chem_recipe.py`'s `ChemSpec`
registry, this makes the recipe layer (`chem_recipe.py`, `recipe.py`,
`registry.py`, `types.py`) a candidate for the separate "keep, merge or remove"
note described above: most of it now has no consumer inside the repo.

## Trigger conditions

**Phase 0 needs no trigger** — it's a decided, no-risk de-duplication and
can be picked up any time, independent of Phase 1.

Pick up **Phase 1** when either happens:

- A model needs a strong ion not on the current allowlist and hits the
  silent-drop gap in practice (a real correctness bug, not just a code
  smell).
- Another charged species is added to `common_species.py` for
  precipitation/complexation work (per
  `MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md`) and someone
  has to manually update the (by then de-duplicated, but still
  allowlist-based) strong-ion table again — a good forcing function to
  do the generalization instead of adding a fourth entry.

## How to start one

Per this folder's usual convention
([README.md](README.md#how-to-start-one)): resolve the open questions
above, write a checklist file
(`STRONG_ION_INFERENCE_GENERALIZATION_CHECKLIST.md`), cut a branch off
`main` (suggested name: `strong-ion-inference-generalization`). Small and
self-contained — no dependency on any other in-flight phase, though it
touches the same `nr_engine.py`/`nr_solver.py` files as the NR
Precipitation and `MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md`
tracks, so re-check for merge overlap if either is in flight when this is
picked up.
