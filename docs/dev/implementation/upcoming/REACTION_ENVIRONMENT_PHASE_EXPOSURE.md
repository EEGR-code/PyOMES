# ReactionEnvironment Gas/Solid Phase Exposure — Design Note

> Status: design note, not yet started. No branch, no checklist, no code
> yet. Written on `main` (independently of the `tutorials-followups` phase
> that surfaced the motivating gap) — see "How to start one" below when
> picked up.

## Motivation

While fixing `tests/validation/speciation/07_iron_oxidation.ipynb` /
`08_iron_oxidation_and_precipitation.ipynb` (tutorials-followups checkpoint
4), the Singer & Stumm (1970) Fe²⁺ oxidation rate law needed `p(O2)` in atm
(its literature-calibrated form), but `KineticReaction.rate_fn` only ever
receives a `ReactionEnvironment` whose `concentrations` are built
**exclusively from the liquid phase's `n_mol`**
(`ControlVolume._build_reaction_environment`, `control_volume.py:924-956` —
confirmed by reading both call sites, `:700-748` and `:924-1003`; neither
merges gas-phase state in, even when the CV has a `GasPhase`). The fix
applied in that notebook was a manual Henry's-law conversion of the
literature `k` into an aqueous-[O2]-compatible constant, computed explicitly
in the notebook rather than by any framework feature.

That conversion is correct and auditable, but it's also a tax every user
writing a rate law calibrated against a gas-phase quantity has to pay by
hand, with no framework support and no obvious place to look for one. Rate
laws expressed in terms of a gas partial pressure are common in
environmental/atmospheric kinetics literature (O₂, CO₂, ozone reactions are
frequently reported this way) — this will recur.

## Is this actually a gap, or the intended design?

Worth stating plainly: `KineticReaction` stoichiometry in every example in
this codebase is declared `phase='liquid'` exclusively (via the `_e()`
convention seen in every tutorial/notebook that declares kinetics — Monod
growth, this Singer-Stumm reaction, etc.), and gas-liquid coupling is
handled by a deliberately separate mechanism
(`PartitionModel`/`KineticTransferModel`/`EquilibriumTransferModel`) that
adjusts the liquid's dissolved concentration via Henry's law/kLa. So "a
kinetic reaction's rate depends on what's dissolved in the liquid" is a
coherent, consistently-applied design story, not an oversight — gas phase
already influences kinetics indirectly, through its effect on the liquid
concentration a transfer model maintains.

The proposal below doesn't change that story. It adds a narrow, read-only,
opt-in piece of additional context to `ReactionEnvironment` for the cases
where a rate law is more naturally (or more defensibly, per its literature
source) expressed in gas-phase terms directly — without requiring every
`rate_fn` author to reimplement the Henry's-law conversion inline, and
without changing how liquid-only rate laws work at all.

## Proposed design

- **`ReactionEnvironment.properties['p_gas_atm']`** (or a dedicated
  `p_gas: Dict[str, float]` field, still to be decided — see Open Questions)
  — partial pressures (atm) of every species present in a linked `GasPhase`,
  populated automatically by `ControlVolume._build_reaction_environment`
  whenever the CV has a `"gas"` phase, alongside the existing
  liquid-concentrations-only `concentrations` dict. Read-only; a rate_fn
  that doesn't need it simply never looks at it — no behavior change for
  every existing `phase='liquid'`-only kinetic reaction.
- Computed via ideal gas law from the `GasPhase`'s own `n_mol`/`V_L`/`T_K`
  (`p_i = n_i R T / V_gas`), the same relationship `GasPhase.p_atm` already
  exposes elsewhere in the codebase (confirmed used in
  `docs/tutorials/D2C_workshop/raw_construction.py`'s `build_liquid_phase`)
  — this proposal is largely "expose an already-computed quantity to one
  more consumer," not new physics.
- Analogous to `env.S(species_id)` for liquid concentrations, a convenience
  accessor (e.g. `env.p(species_id)`) on `ReactionEnvironment` would keep
  rate-law code readable: `rate = k * env.S('Fe2+') * env.p('O2') * env.S('OH-')**2`
  directly matching the literature's own notation, no inline Henry's-law
  conversion required.

## Solid-phase exposure — same question, raised alongside this one

The same reasoning extends to `SolidPhase` (already exists in
`PyOMES.core.phases`, per its own docstring: "gas+liquid-only constraint is
now stated explicitly... a third phase could plausibly be added to this
same class later"). A kinetic rate law for a surface-catalyzed or
dissolution/precipitation-coupled reaction might equally want to reference
solid-phase state (e.g., total mineral surface area, or a solid's own
`n_mol`) the same way this proposal exposes gas partial pressure. Flagging
this now so the two aren't designed inconsistently if/when both land —
whoever picks this up should design the `ReactionEnvironment` extension
generally (partition by phase key, not two bespoke `p_gas`/`solid_mol`
fields bolted on separately) rather than solving only the gas case and
re-opening the same question for solids later. No solid-phase use case has
actually surfaced yet (unlike the gas case, which has a concrete motivating
bug) — this is a forward-compatibility note, not a second scoped deliverable.

## Open questions (need a decision before implementing)

1. **Field shape.** A single `properties['p_gas_atm']` dict (minimal API
   surface change, consistent with how `ionic_strength` already rides in
   `properties`) vs. a dedicated typed field on `ReactionEnvironment`
   (`p_gas: Dict[str, float]`, more discoverable, slightly bigger API
   surface). Recommend the dedicated field if solid-phase exposure (above)
   is being designed at the same time — a `Dict[str, Dict[str, float]]`
   keyed by phase (`{"gas": {...}, "solid": {...}}`) would generalize
   cleanly to both; a single `properties` dict bag would not.
2. **What about CVs with no gas phase?** Should `env.p(...)` raise, return
   `0.0`, or return `None`/absent-key, when no `GasPhase` is linked? Given
   `env.S()` already defaults missing liquid species to `0.0`, matching
   that convention (`env.p()` defaults to `0.0`) is probably right for
   consistency, but worth confirming it doesn't silently mask a
   rate-law-vs-CV-topology mismatch the way the original K_SS bug did.
3. **Multi-gas-phase CVs.** Does any current or planned topology attach more
   than one `GasPhase` to a single `ControlVolume`? If not (seems unlikely
   given `ControlVolume.phases` is keyed by a fixed set of phase-role
   strings, not an arbitrary list), a flat `{species_id: p_atm}` dict is
   sufficient; if so, this needs a phase-key dimension too.

## Implementation steps (once the questions above are resolved)

1. Resolve Open Question 1 (field shape) — informs everything else.
2. Add gas-phase partial-pressure computation to both
   `ControlVolume._build_reaction_environment` call sites
   (`control_volume.py:700-748` and `:924-1003`) — reuse `GasPhase.p_atm`
   rather than recomputing the ideal-gas relationship inline.
3. Add the `env.p(species_id)` convenience accessor to
   `PyOMES/reactions/environment.py`, mirroring `env.S()`.
4. Regression check: existing `phase='liquid'`-only kinetic reactions
   (Monod growth, this Singer-Stumm reaction) must see zero behavior
   change — this is purely additive.
5. Update `07_iron_oxidation.ipynb`/`08_iron_oxidation_and_precipitation.ipynb`
   (tutorials-followups checkpoint 4) to use `env.p('O2')` directly against
   the literature `k`, dropping the manual Henry's-law pre-conversion —
   the cleanest available end-to-end demonstration that the feature closes
   the gap it was built for.
6. Unit tests: a `rate_fn` referencing `env.p('O2')` on a CV with a linked
   `GasPhase` returns the expected ideal-gas partial pressure; the same on
   a CV with no gas phase returns the Open-Question-2 default without
   raising.

## How to start one

Per this folder's usual convention
([README.md](README.md#how-to-start-one)): write a checklist file
(`REACTION_ENVIRONMENT_PHASE_EXPOSURE_CHECKLIST.md`), cut a branch off
`main` (suggested name: `reaction-environment-phase-exposure`), and work the
checkpoints there. Small and self-contained — no dependency on any other
in-flight phase.
