# Chemistry Unification 3 — Checklist

> **Status: Shipped 2026-05-15** — all 11 active checkpoints landed
> across 9 commits on the `chemistry-unification-3` branch. 778
> standalone tests passing (758 baseline + 20 net new: cross-phase
> reaction coverage in `test_reactions.py` and
> `test_gas_liquid_link.py`, canonical naming + suffix-parser coverage
> in `test_speciation.py`). BSM2 reference sentinels re-baselined with
> the NH4+/I bug-fix rationale documented in `test_bsm2_reference.py`
> (pH 3.3316 → 3.3402, S_h2 ~3% drift, other species 1e-7 to 1e-5
> relative; `RTOL_SENTINEL=1e-9` retained). Checkpoint 7 (builder
> migrations) deferred to Phase 3b alongside the unified Species-ID
> emission convention. Tag: `chemistry-unification-3-shipped`.

Working checklist for `chemistry-unification-3`. Conceptual framing is
in [../phases-upcoming/CHEMISTRY_UNIFICATION.md](../phases-upcoming/CHEMISTRY_UNIFICATION.md) under "Phase E"
and "Leakage That Won't Go Away Automatically → Leak 2". The branch
slicing and resolved scope-shape are in
[../phases-upcoming/CHEMISTRY_UNIFICATION_PLAN.md](../phases-upcoming/CHEMISTRY_UNIFICATION_PLAN.md) under
"Phase 3 — `chemistry-unification-3`". This file is the implementation
plan: file-level edits, ordered checkpoints, and sanity checks.
Modelled on
[CHEMISTRY_UNIFICATION_2_CHECKLIST.md](CHEMISTRY_UNIFICATION_2_CHECKLIST.md).

## Goal

Two structural pieces that share one BSM2 re-baseline event:

1. **Cross-phase equilibrium reactions (Leak 2 closure).** Allow
   equilibrium reactions to span phases (e.g.
   `CO2(gas) ⇌ CO2aq(liquid)`). `KineticGasLiquidLink` derives its
   `speciation_keys` mapping from these declared reactions instead of
   accepting it as a hand-passed kwarg. The link's
   gas → mol_key wiring becomes a downstream consequence of the
   chemistry declaration, not a parallel string mapping.
2. **Ionic-strength dedup at source (Phase 2 follow-up).** Pick one
   canonical naming convention in
   [`_compute_species_eq`](../../src/speciation/acid_base.py#L1112)
   for species the engine knows by name (CO₂ ladder, phosphate
   ladder, sulfate ladder, NH₃/NH₄⁺ ladder). Drop the parallel
   generic `{name}_HA`/`{name}_A-`/`{name}_BH+`/`{name}_B` emissions
   for those recognised species. Rewrite
   [`ionic_strength_from_speciation`](../../src/speciation/activity.py#L52)
   so canonical names count via the z-dict and the `_A-` fallback
   loop survives ONLY for unrecognised VFA-style acids whose names
   the engine has no canonical entry for. Removes the double-count
   bug class the Phase 2 deviation note tagged for this phase.

After this branch ships:

- [`src/reactions/reaction.py`](../../src/reactions/reaction.py)
  accepts equilibrium reactions whose stoichiometry spans two phases.
  `log_K` is optional for cross-phase reactions (the Henry constant
  governing the partition lives on the link, not on the reaction).
  Single-phase equilibrium reactions still require `log_K`.
- [`src/reactions/reaction_set.py`](../../src/reactions/reaction_set.py)'s
  `partition()` still returns `(kinetic, equilibrium)`. Cross-phase
  equilibrium reactions go in the equilibrium bucket; consumers
  filter further by phase span.
- [`src/speciation/engine.py`](../../src/speciation/engine.py)'s
  `from_reactions` skips cross-phase equilibrium reactions — the
  speciation engine only solves single-phase (liquid) acid-base
  equilibria. The link handles cross-phase ones.
- [`src/core/gas_liquid_link.py`](../../src/core/gas_liquid_link.py)
  gains `derive_speciation_keys(reaction_set)`. Iterates over
  equilibrium reactions, finds those with one gas-phase entry and
  one liquid-phase entry, and populates
  `self.speciation_keys[gas_species_id] = liquid_species_id`. Wired
  into `ControlVolume.__init__` after `reaction_model` attach so
  declared chemistry drives link configuration. The `speciation_keys`
  constructor argument is removed (clean break; no manual mapping).
- [`src/speciation/acid_base.py`](../../src/speciation/acid_base.py)'s
  `_compute_species_eq` emits canonical names ONLY (no parallel
  `_HA`/`_A-`/`_BH+`/`_B` emission) for recognised acid systems
  (`eq_def.name ∈ {"CO2", "phosphate", "bisulfate", "NH3", "NH4"}`).
  Unrecognised acids (VFAs, custom user systems) keep the generic
  emission as before.
- [`src/speciation/activity.py`](../../src/speciation/activity.py)'s
  `ionic_strength_from_speciation` is rewritten around
  charge-suffix parsing: trailing `+`/`-` tokens infer the
  ion's charge, with a small override table for non-conforming
  keys (`Cation(inert)` / `Anion(inert)`). The bipartite z-dict
  + `_A-` fallback loop collapses into one pass over the species
  dict. The "Phase 3 cleanup" warning block on the function is
  deleted.
- [`src/speciation/engine.py`](../../src/speciation/engine.py)'s
  `solve()` alpha dispatcher simplifies: the canonical-or-generalised
  `_lookup_mol("CO2aq", "CO2_HA")` fallback chain collapses to a
  direct canonical-name lookup. The deviation-note workaround from
  Phase 2 evaporates.
- **BSM2/ADM1 builder migration is deferred to Phase 3b**
  (decision 2026-05-15, see *Resolved decisions* below). The
  cross-phase reaction infrastructure ships in Phase 3 but no
  production builder declares cross-phase reactions yet; the
  `derive_speciation_keys` hook is exercised by a dedicated
  unit test instead. BSM2 keeps its
  `TransferConfig(species={}, speciation_keys={"CO2": "CO2aq"})`
  construction; ADM1 keeps its
  `.speciation_correction(species, mol_key)` builder calls.
  Phase 3b migrates these alongside the unified Species-ID
  emission convention.
- BSM2 sentinels re-baseline once (canonical naming is a numerical
  no-op for `bsm2_default`-path BSM2 — generic `_H2A`/`_HA-`/`_A--`
  forms it emits today don't collide with `_A-` — but the
  `from_reactions` path used by `build_bsm2_cv` does emit `CO2_A-`
  today, and dropping that emission alongside the canonical
  `HCO3-` addition is a structural change that must be sentinel-
  verified). RTOL_SENTINEL=1e-9 expected to hold.

## Out of scope

Everything reserved for later chemistry-unification phases or
trigger-gated deferral:

- **Phase 3b (deferred, trigger-gated like
  [CONTAINER_LAYERING.md](CONTAINER_LAYERING.md))**: `ChemistryDatabase`
  dataclass with `.extend()`, `cv.chemistry_db` field, stock database
  modules (`aqueous.py`, `bioprocess_basic.py`,
  `anaerobic_digestion.py`), `ThermoFramework` (engine activity-knob
  absorption + `ThermodynamicConfig` absorption), and the
  `core/links.py` `ThermoFramework` mismatch warning. The trigger is
  a concrete pull for run-time chemistry composition that the
  module-level `Species` import pattern shipped in Phase 1 doesn't
  serve. Decision recorded in *Resolved decisions* below.
- **Phase 4 (`chemistry-unification-4`)**: `AccuracyMonitor`,
  `AccuracyWarning`, `WarningConfig`.
- **Deleting `ThermodynamicConfig`.** Phase 2 trimmed
  `apply_to_cv` to a back-reference; Phase 3 leaves it in place
  because `build_bsm2_cv` / `build_adm1_cv` still call
  `pKa_at_T` / `Kw_at_T` / `context_kwargs_at_current_T` indirectly
  via `make_bsm2_chem_env_fn`. `ThermodynamicConfig` absorbs into
  `ChemistryDatabase` in Phase 3b when that phase triggers.
- **Rewriting `validate_thermodynamics`.** The Phase 2 shrink left
  it sourcing from `_thermo_config` and `chem_env_fn` snapshots
  only. Phase 3 leaves it as-is — it will be rewritten around
  `ChemistryDatabase` when 3b triggers.
- **Henry constants in the cross-phase reactions.** The
  `log_K` on a cross-phase equilibrium reaction is optional in
  Phase 3 and consumed by no one. A future phase may use it to
  drive `link.henry[species]`; this branch leaves Henry on
  `TransferConfig` / `link.henry` exactly as Phase 2 does.

## Resolved decisions

Pinned during the planning walkthrough on 2026-05-15:

- **Bundle the ionic-strength dedup with 3a, do not split into a
  precursor branch.** The plan doc calls 3a the "natural home" —
  cross-phase reactions pin the canonical naming convention because
  the link's auto-derived `speciation_keys` need to point at
  whatever `_compute_species_eq` emits. Doing both together
  enforces one convention end-to-end. The Phase 2 deviation note
  already proved any change to `_compute_species_eq`'s emitted
  keys forces a BSM2 sentinel re-baseline; a standalone precursor
  branch buys no risk reduction (it would re-baseline once, then
  3a would re-baseline again). One branch, one re-baseline event.

- **Defer 3b (database packaging) entirely; ship 3a only.**
  `ChemistryDatabase` + stock databases + `ThermoFramework` is pure
  refactor/packaging on top of pieces that already exist (Phase 1
  promoted `Species` to module-level imports, which handles BSM2's
  and ADM1's static-composition case today). No concrete pull
  exists for run-time chemistry composition. Marking 3b as
  trigger-gated like `CONTAINER_LAYERING.md` keeps the branch
  count honest. Concretely: after `chemistry-unification-3`,
  `chemistry-unification-4` stays Phase F (`AccuracyMonitor`)
  as planned in
  [CHEMISTRY_UNIFICATION_PLAN.md](CHEMISTRY_UNIFICATION_PLAN.md);
  3b becomes a deferred phase that gets pulled in when a concrete
  use case appears (composing kits at runtime, a new ThermoFramework
  consumer, etc.). The Phase 3 plan doc is updated post-ship to
  reflect this.

- **Canonical names replace generic emission for recognised
  species, not add to it.** The first-principles answer is one
  convention per species. The recognised-species table (CO₂,
  phosphate, bisulfate, NH₃/NH₄⁺) gets canonical-only emission
  from `_compute_species_eq`; the generic `_HA`/`_A-`/`_BH+`/`_B`
  emission for these species is dropped. Unrecognised VFA-style
  acids (S_ac, S_pro, S_bu, S_va, custom user systems) continue to
  use generic emission — they have no canonical name to replace
  it with, and the `_A-` loop in `ionic_strength_from_speciation`
  picks them up safely (no double-count because no canonical
  name in the z-dict matches). The recognised-species table lives
  inside `_compute_species_eq` for now; when Phase 3b triggers
  it moves onto the `Species` declarations themselves.

- **`speciation_keys` constructor argument retained for this
  phase.** Phase 3 ships the `derive_speciation_keys`
  infrastructure (additive — existing entries survive when no
  cross-phase reaction overrides them) but does not migrate
  BSM2/ADM1 builders away from the explicit kwarg. Reason: the
  alpha-key naming convention question that surfaced during
  implementation (see "Alpha-key convention deferral" below)
  has a clean answer (engine emits Species IDs directly) that
  belongs in Phase 3b's unified ChemistryDatabase rollout.
  Forcing the migration now requires either a translation
  rule (~10 lines, deletable in 3b) or partial pull-forward
  of 3b's emission rewrite. Both expand scope past what
  Phase 3 was carved for. Phase 3b's mandate now explicitly
  includes this migration; 3a ships the infrastructure clean
  with an exercising unit test.

- **Alpha-key convention deferral (decision 2026-05-15).**
  Today's `PropertyResult.alphas` keys (`"CO2aq"`,
  `"S_ac_HA"`) are *synthesized* by `_compute_species_eq`
  rather than reading Species IDs from the reaction
  declarations. This mismatch means `derive_speciation_keys`
  can't trivially map `gas_id → liq_id` — it would need a
  translation table or category-aware lookup. The clean fix
  is to make the engine emit species using their declared
  Species IDs directly (no `_HA`/`_A-`/`_BH+`/`_B`/`aq`
  decorators); then `derive_speciation_keys` is one line and
  no translation rule exists anywhere. That fix requires
  carrying Species references through the engine
  (`EquilibriumDef` gains a `species_refs` field; the
  legacy `bsm2_default()` factory adopts Species-based
  construction; every test asserting on `out["CO2aq"]` /
  `out["S_ac_HA"]` / etc. migrates). Scope is closer to
  Phase 3b than 3a. Discussion 2026-05-15; the long-term
  answer is recorded in
  [CHEMISTRY_UNIFICATION_PLAN.md](CHEMISTRY_UNIFICATION_PLAN.md)
  under Phase 3b's mandate.

- **Cross-phase equilibrium reactions don't require `log_K`.**
  The Henry constant governing the gas-liquid partition lives on
  `link.henry` (with Phase 2's `henry_params` for T-dependence).
  A cross-phase equilibrium reaction is *declarative wiring* — it
  tells the link which liquid molecular form corresponds to a
  given gas species. It is not a thermodynamic constraint that
  the engine solves. `Reaction.__init__` validation allows
  `log_K=None` iff `kind="equilibrium"` AND the stoichiometry
  spans more than one phase. Single-phase equilibrium reactions
  still require `log_K` (they go to `SpeciationEngine.from_reactions`).

- **`SpeciationEngine.from_reactions` filters out cross-phase
  reactions silently.** The engine's existing
  `equilibrium_reactions` iteration in
  [`engine.py:194-211`](../../src/speciation/engine.py#L194)
  already classifies each reaction via `_classify_equilibrium`.
  Add a pre-check: if the reaction spans more than one phase,
  skip it (no warning — cross-phase reactions are valid input,
  just not consumed by this code path). The link picks them up
  via its own iteration.

- **Suffix-based charge parsing in
  `ionic_strength_from_speciation`.** Replace the bipartite
  z-dict + `_A-` fallback loop with a single charge-from-suffix
  rule: trailing `+`/`++`/`-`/`--`/`---` tokens infer the
  ion's charge magnitude and sign. A small `_CHARGE_OVERRIDES`
  table handles ions whose ids don't end in a charge token
  (BSM2's `Cation(inert)` / `Anion(inert)`). One pass over the
  dict instead of two; the 35-line z-dict collapses into a
  4-line override table. Rationale: one rule covers both
  canonical (`HCO3-`, `CO3--`, `Mg++`, `PO4---`) and generic
  (`S_ac_A-`, `S_pro_A-`, …) emissions uniformly. Adding a new
  recognised acid in the future (e.g. H₂S → sulfide ladder)
  needs no edit to this function — its `HS-` / `S--` keys
  parse correctly by construction. Removes the entire
  double-count bug class flagged by the Phase 2 deviation
  note: the function cannot double-count regardless of input
  because each key contributes exactly one charge based on its
  suffix. Discussion: 2026-05-15, resolved with user.
  Considered alternatives (A: trust the upstream
  contract; B: defensive skip-set encoding the recognised /
  generic split in two places) were rejected — A leaves a
  silent-failure mode for test-isolation dict construction; B
  encodes the recognised-acid set in two places (here and in
  `_compute_species_eq`'s canonical-name table) and the user
  pointed out the generalised / readable / simple priorities
  align with the suffix-parsing answer.

## Checkpoints

The checkpoints are ordered so each leaves the test suite in a
runnable state. The dedup (1–3) lands first because it's confined to
the speciation module and exposes the canonical naming surface that
the link will consume in checkpoint 5. Cross-phase reactions and link
auto-derivation (4–6) follow. Builder migrations (7) trigger the
BSM2 sentinel re-baseline (8). Tests and ship (9–11) close out.

### 1. Canonical naming in `_compute_species_eq`

- [ ] [`src/speciation/acid_base.py:1112-1176`](../../src/speciation/acid_base.py#L1112)
      — introduce a module-level constant
      `_CANONICAL_NAMES: Dict[str, Tuple[str, ...]]` mapping
      `eq_def.name` aliases to ordered canonical species keys
      (most protonated first). Examples:

      ```python
      _CANONICAL_NAMES = {
          "CO2":       ("CO2aq", "HCO3-", "CO3--"),
          "CO2aq":     ("CO2aq", "HCO3-", "CO3--"),
          "phosphate": ("H3PO4", "H2PO4-", "HPO4--", "PO4---"),
          "H3PO4":     ("H3PO4", "H2PO4-", "HPO4--", "PO4---"),
          "bisulfate": ("H2SO4", "HSO4-", "SO4--"),
          # cation acids:
          "NH4":       ("NH4+", "NH3"),
          "NH3":       ("NH4+", "NH3"),
      }
      ```

      The aliasing handles the asymmetry between the
      `bsm2_default()` factory path (uses `name="NH4"`,
      `name="CO2"`) and the `from_reactions` path (uses whatever
      `total_id` or acid `Species.id` the reaction declares,
      typically `"NH3"`/`"CO2aq"` depending on the speciation
      direction).
- [ ] Rewrite the polyprotic acid branch (lines 1139–1162) to
      consult `_CANONICAL_NAMES`. When `eq_def.name` matches a
      registered key, emit `out[canonical[i]] = CT * alphas[i]`
      for `i in range(n + 1)` and **skip** the generic
      `_species_key(...)` emission entirely. When `eq_def.name`
      doesn't match, keep the existing generic emission (VFAs).
      Drop the special-case `if eq_def.name == "CO2" and n >= 2`
      block — it's subsumed by the registry.
- [ ] Rewrite the cation_acid branch (lines 1125–1137) the same
      way: when `eq_def.name` matches `_CANONICAL_NAMES`, emit
      `out["NH4+"] = BH` and `out["NH3"] = B` directly; skip
      `_{name}_BH+` / `_{name}_B` generic emission. When it
      doesn't match (no current callers, but for future custom
      cation acids), keep the generic emission.
- [ ] Verify by grep that the only Phase 2 deviation-note
      workaround in `engine.solve()`'s alpha dispatcher
      ([`engine.py:415-440`](../../src/speciation/engine.py#L415))
      becomes unnecessary: the canonical-or-generalised fallback
      chain `_lookup_mol("CO2aq", "CO2_HA")` and
      `_lookup_mol("NH3", "NH3_B", "NH4_B")` will simplify in
      checkpoint 3. Don't simplify yet — wait until activity.py
      and the engine alpha block are touched together.

Sanity check: BSM2 engine via `bsm2_default()` no longer emits
`CO2_H2A` / `CO2_HA-` / `CO2_A--` keys; the species dict carries
`CO2aq` / `HCO3-` / `CO3--` only. `from_reactions`-built engines
emit canonical names for declared CO2/NH3 ladders. VFA emissions
(`S_ac_HA`, `S_ac_A-`, etc.) are unchanged.

### 2. Rewrite `ionic_strength_from_speciation` around suffix parsing

- [ ] [`src/speciation/activity.py:52-153`](../../src/speciation/activity.py#L52)
      — replace the bipartite z-dict + `_A-` fallback loop with
      a single charge-from-suffix pass. Introduce a helper:

      ```python
      def _charge_from_suffix(key: str) -> int:
          """Infer ion charge from trailing +/- tokens.

          Returns positive count of '+' chars if key ends in '+',
          negative count of '-' chars if it ends in '-', 0 otherwise.
          """
          if key.endswith("+"):
              return len(key) - len(key.rstrip("+"))
          if key.endswith("-"):
              return -(len(key) - len(key.rstrip("-")))
          return 0
      ```

      And a small override table for ions whose ids don't carry
      a suffix charge token:

      ```python
      _CHARGE_OVERRIDES = {
          "Cation(inert)": +1,
          "Anion(inert)": -1,
      }
      ```

      The function body collapses to one pass:

      ```python
      def ionic_strength_from_speciation(sp):
          I_sum = 0.0
          for key, v in sp.items():
              z = _CHARGE_OVERRIDES.get(key)
              if z is None:
                  z = _charge_from_suffix(key)
              if z == 0:
                  continue
              try:
                  I_sum += float(v) * z * z
              except (TypeError, ValueError):
                  warnings.warn(...)  # mirror existing behaviour
          return 0.5 * I_sum
      ```

- [ ] Delete the 35-line `z` dict (lines 89–131 in current code)
      — replaced by the suffix parser + 2-entry override table.
- [ ] Delete the `_A-`-suffix fallback loop (lines 141–151) —
      subsumed by the unified pass.
- [ ] Delete the `.. warning::` block in the docstring (Phase
      2's tag-for-Phase-3). Replace with a short description of
      the new contract: charge is inferred from the trailing
      `+`/`-` tokens on each species key; ions whose ids don't
      follow this convention (`Cation(inert)`, `Anion(inert)`)
      are listed in `_CHARGE_OVERRIDES`. Future ions following
      the suffix convention work automatically.
- [ ] Verify the override table is sufficient by grepping for
      every species id written into a `sp` dict that
      `ionic_strength_from_speciation` might see: scan
      `_compute_species_eq` ([`acid_base.py:1112-1176`](../../src/speciation/acid_base.py#L1112))
      and the strong-ion echo block (lines 1163–1174) for every
      key it might emit. Any key that doesn't end in `+`/`-`
      *and* contributes charge must be in `_CHARGE_OVERRIDES`.

Sanity check (worked example):
`ionic_strength_from_speciation({"H+": 1e-7, "OH-": 1e-7,
"HCO3-": 0.05, "CO3--": 0.005, "S_ac_A-": 0.001,
"Cation(inert)": 0.02, "Mg++": 0.003})` should return
`0.5 * (1e-7 + 1e-7 + 0.05 + 0.005*4 + 0.001 + 0.02 + 0.003*4)
= 0.5 * (~0.1130002) ≈ 0.0565` mol/L. Each key contributes
exactly once, by the suffix rule or the override table; the
two summation passes of the legacy implementation collapse
into one with no overlap by construction.

Sanity check (regression): a dict containing both
`{"HCO3-": 0.05, "CO2_A-": 0.05}` (legacy double-emission
scenario, no longer producible by `_compute_species_eq` after
checkpoint 1) is counted as `0.5 * (0.05 + 0.05) = 0.05` — both
entries contribute since suffix parsing treats them as
independent ions. This is the correct semantics: the function
takes the dict as ground truth. Eliminating double-emission is
checkpoint 1's job (single source of truth at the emitter), not
this function's. The behaviour difference vs the legacy
implementation is irrelevant in practice because checkpoint 1
ensures no caller produces such a dict from real chemistry.

### 3. Simplify the engine's alpha dispatcher

- [ ] [`src/speciation/engine.py:415-451`](../../src/speciation/engine.py#L415)
      — simplify `_lookup_mol(*keys)` to a single-key lookup
      `out.get("CO2aq")`, `out.get("NH3")`, `out.get(mol_key)`
      because checkpoint 1 guarantees canonical names are
      always present when the totals are present. Delete the
      `_lookup_mol` helper if it's now a one-liner.
- [ ] Drop the deviation-note comment block (lines 396–400 in
      current code) that motivated the canonical-or-generalised
      fallback — the underlying double-naming problem it
      describes is resolved.
- [ ] Acid totals iteration (lines 442–449) keeps using
      `f"{name}_HA"` — VFAs continue with generic emission.
      Comment to clarify the asymmetry: recognised acid ladders
      use canonical (CO2aq, NH3); unrecognised acids
      (VFAs, custom) use `{name}_HA`.

Sanity check: after checkpoint 3, BSM2 engine still produces
`out["alphas"]["CO2aq"]` ≈ 0.722 and
`out["alphas"]["NH3"]` ≈ value matching the existing dispatcher's
output (the change is the lookup path, not the math).

### 4. Cross-phase support in `Reaction` and `partition()`

- [ ] [`src/reactions/reaction.py`](../../src/reactions/reaction.py)
      — relax `__init__` validation: when `kind="equilibrium"` AND
      the stoichiometry spans more than one distinct phase,
      `log_K` MAY be `None`. Single-phase equilibrium reactions
      still require `log_K`. Add a helper property
      `Reaction.is_cross_phase` returning
      `len({entry.phase for entry in self.stoichiometry}) > 1`.
      Use it in the validation message for clarity.
- [ ] Update the `Reaction` docstring's "Validation" section to
      describe the cross-phase case explicitly: cross-phase
      equilibrium reactions are *partition declarations* consumed
      by gas-liquid links, not by the speciation engine.
- [ ] [`src/reactions/reaction_set.py`](../../src/reactions/reaction_set.py)
      — `partition()` API unchanged (returns
      `(kinetic, equilibrium)`). Add a docstring note: the
      `equilibrium` bucket may contain cross-phase reactions;
      consumers filter further via `rxn.is_cross_phase`.

Sanity check: constructing a cross-phase
`Reaction(kind="equilibrium", log_K=None, stoichiometry=[...gas..., ...liquid...])`
succeeds. Same with `log_K=None` and a single phase fails.
Single-phase `kind="equilibrium"` with `log_K` set still works.

### 5. `SpeciationEngine.from_reactions` filters cross-phase

- [ ] [`src/speciation/engine.py:190-222`](../../src/speciation/engine.py#L190)
      — in `from_reactions`, add a pre-check inside the
      iteration loop: `if rxn.is_cross_phase: continue`. The
      engine consumes single-phase liquid acid-base reactions
      only. Cross-phase reactions in the same bucket are
      silently skipped; the link picks them up.
- [ ] Update the docstring to mention the filter — one bullet
      under "Stoichiometry contract".

Sanity check: a reaction set mixing single-phase liquid
equilibrium (e.g. acetate/acetic acid) with cross-phase
equilibrium (e.g. CO2(gas) ⇌ CO2aq(liquid)) is consumed
correctly: the engine builds an EquilibriumSet from the liquid
ones only; the cross-phase ones don't error.

### 6. `KineticGasLiquidLink.derive_speciation_keys` (infrastructure-only, additive)

> **Scope adjustment 2026-05-15.** Originally specified as
> "remove the `speciation_keys` constructor argument" — that's
> deferred to Phase 3b. Phase 3 ships the new mechanism
> *additively* so existing callers (BSM2's
> `TransferConfig(speciation_keys={"CO2": "CO2aq"})`, ADM1's
> `.speciation_correction(...)` builder calls) keep working.
> Rationale: the alpha-key naming-convention question that
> would let the migration drop the kwarg cleanly has a
> dependency on engine emission semantics that belongs in
> Phase 3b. See *Resolved decisions* → "Alpha-key convention
> deferral".

- [x] [`src/core/gas_liquid_link.py`](../../src/core/gas_liquid_link.py)
      — added method `derive_speciation_keys(self, reaction_set) -> None`:
      iterates over equilibrium reactions in `reaction_set`,
      filters by `rxn.is_cross_phase`, and for each cross-phase
      reaction with exactly one gas-phase entry and one
      liquid-phase entry, sets
      `self.speciation_keys[gas_entry.species.id] = liquid_entry.species.id`.
      Skips reactions whose phase keys don't match
      `(self.gas_phase_key, self.liquid_phase_key)` (defensive —
      the same reaction set may serve multiple links in a
      multi-CV system).
- [x] Method is **additive**: it leaves any pre-existing
      entries in `speciation_keys` alone unless overridden by a
      declared reaction. The constructor argument and its
      default factory `{"CO2": "CO2aq"}` are **retained** —
      Phase 3b deletes them alongside the unified Species-ID
      emission convention.
- [x] [`src/core/control_volume.py`](../../src/core/control_volume.py)
      — `__init__`, after `reaction_model` is attached, iterates
      `self.internal_interfaces`, finds each
      `KineticGasLiquidLink`, and calls
      `link.derive_speciation_keys(rxns)`. Triggers
      automatically; if no cross-phase reactions are declared,
      this is a no-op.

Sanity check: `Grep "speciation_keys=" src/ models/` still
returns the BSM2 / TransferConfig / factory references — that
is the *intended* state for Phase 3. Phase 3b drops them.

### 7. Builder migrations — DEFERRED to Phase 3b

> **Scope adjustment 2026-05-15.** All builder migrations
> originally scoped here are deferred to Phase 3b. The
> infrastructure for cross-phase reactions is shipped in
> checkpoints 4-6; the *consumers* (BSM2 declaring its CO₂
> gas-liquid partition; ADM1 declaring NH₃/H₂S/VFA partitions)
> migrate when Phase 3b's unified Species-ID emission lands.
> Doing it now requires a translation rule that Phase 3b would
> immediately delete.

Items deferred:

- BSM2 builder: drop `speciation_keys={"CO2": "CO2aq"}` from
  `TransferConfig` construction in
  [`models/vlmodels/adm1/bsm2.py:795`](../../models/vlmodels/adm1/bsm2.py#L795);
  declare a cross-phase CO₂ gas-liquid partition reaction in
  `_build_bsm2_equilibrium_reactions`.
- ADM1 builder: replace the three
  `.speciation_correction(species, mol_key)` calls in
  [`models/vlmodels/adm1/base.py:979-987`](../../models/vlmodels/adm1/base.py#L979)
  with cross-phase reaction declarations.
- `TransferConfig.speciation_keys` field
  ([`models/vlmodels/fermenter/config/configs.py:208`](../../models/vlmodels/fermenter/config/configs.py#L208))
  and its `to_dict`/`from_dict` plumbing — delete.
- `FermenterBuilder.speciation_correction` method
  ([`models/vlmodels/fermenter/config/builder.py:237`](../../models/vlmodels/fermenter/config/builder.py#L237))
  — delete.
- `factory.py`'s `speciation_keys = dict(transfer.speciation_keys)`
  gather and `speciation_keys=` constructor kwarg
  ([`models/vlmodels/fermenter/config/factory.py:255-336`](../../models/vlmodels/fermenter/config/factory.py#L255))
  — delete.

Sanity check: skipped — no migration ships this phase.

### 7. Builder migrations

Each builder that previously passed
`TransferConfig.speciation_keys=...` or
`KineticGasLiquidLink(speciation_keys=...)` is rewritten to
declare the cross-phase equilibrium reaction(s) in its reaction
set instead. The link picks them up via the
`ControlVolume.__init__` hook from checkpoint 6.

- [ ] [`models/vlmodels/adm1/bsm2.py:790-795`](../../models/vlmodels/adm1/bsm2.py#L790)
      — `build_bsm2_cv`: drop the `speciation_keys={"CO2": "CO2aq"}`
      hand-passed kwarg on `TransferConfig`. Add a CO2 gas-liquid
      equilibrium reaction to BSM2's `ReactionSet` (alongside the
      existing CO2(aq)/HCO3⁻ liquid-phase equilibrium):

      ```python
      Reaction(
          label="CO2(g) <-> CO2(aq)",
          kind="equilibrium",
          stoichiometry=[
              StoichiometryEntry(CO2_gas, "gas", -1),
              StoichiometryEntry(CO2, "liquid", +1),
          ],
          # log_K=None: Henry's law lives on the link
      )
      ```

      Add a `CO2_gas: Species` declaration to BSM2's species
      module (or import from `common_species` if appropriate —
      check Phase 1's `common_species.py` for an existing entry).
- [ ] [`models/vlmodels/adm1/base.py`](../../models/vlmodels/adm1/base.py)
      — `build_adm1_cv`: the three `.speciation_correction(...)`
      builder calls (which currently set `speciation_keys`) are
      replaced by `Reaction(...)` declarations in the ADM1
      reaction set for each volatile gas (CO2, NH3, H2S, CH4, H2
      where speciation matters). Where the underlying acid-base
      equilibrium is already declared in
      `build_adm1_reactions()`, only the cross-phase wiring needs
      to be added.
- [ ] [`models/vlmodels/fermenter/config/configs.py:208-272`](../../models/vlmodels/fermenter/config/configs.py#L208)
      — delete the `speciation_keys` field from `TransferConfig`.
      Update `to_dict`/`from_dict` accordingly. Update the
      docstring to point at the cross-phase-reaction declaration
      path.
- [ ] [`models/vlmodels/fermenter/config/builder.py:277-280`](../../models/vlmodels/fermenter/config/builder.py#L277)
      — delete the `.speciation_correction(species, mol_key)`
      builder method (after this phase it's a no-op — the link
      derives `speciation_keys` from the reaction set, not from
      builder calls). Callers should declare the cross-phase
      reaction in their reaction set instead. The shipped tag
      `chemistry-unification-2-shipped` is the last point where
      this method has work to do.
- [ ] [`models/vlmodels/fermenter/config/factory.py:253-336`](../../models/vlmodels/fermenter/config/factory.py#L253)
      — drop the `speciation_keys = dict(transfer.speciation_keys)`
      gather and the `speciation_keys=` kwarg on the
      `KineticGasLiquidLink(...)` constructor call.

Sanity check: `Grep "TransferConfig.*speciation_keys|
speciation_keys=" models/ src/` returns matches only in
test-construction sites (which migrate in checkpoint 9) and
docstrings/comments. No live constructor or builder code path
references the field.

### 8. BSM2 sentinel re-baseline

- [ ] [`tests/standalone/test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py)
      — run the test; expect failures.
      `from_reactions`-built BSM2 (the current `build_bsm2_cv`)
      previously emitted `CO2_HA` / `CO2_A-` keys and read alphas
      via the canonical-or-generalised fallback. After
      checkpoints 1–3, it emits `CO2aq` / `HCO3-` canonical keys
      directly. The numerical species dict is *almost* identical
      (same values, different keys); the ionic-strength sum
      changes because the previously-counted `CO2_A-` charge
      contribution now comes through `HCO3-` via the z-dict
      (algebraically equivalent — these are the same charge).
      **Expected drift:** negligible if checkpoint 1's
      implementation is clean; ≤1e-12 from float reordering.
      Document any non-trivial drift in the test's comment
      block, mirroring the Phase 2 re-baseline note structure.
- [ ] If drift exceeds the existing `RTOL_SENTINEL = 1e-9`
      budget, **stop and investigate before rebaselining.** The
      dedup is intended to be a numerical no-op for BSM2
      (canonical = generalised algebraically; only the keys
      change). Any drift > 1e-9 indicates the canonical-name
      table or the `_compute_species_eq` rewrite has a real bug.
- [ ] Add a comment block to `test_bsm2_reference.py` describing
      the chemistry-unification-3 change, regardless of whether
      sentinels needed updating — mirrors the structure of the
      Phase 1 and Phase 2 notes already in that file.

Sanity check: `RTOL_SENTINEL = 1e-9` continues to hold (or
re-baseline with explanation if not). BSM2 trajectory unchanged
to the documented tolerance.

### 9. Tests — rewrites

- [ ] [`tests/standalone/test_speciation.py`](../../tests/standalone/test_speciation.py)
      — tests that assert on `species["CO2_HA"]` / `species["CO2_A-"]`
      / `species["NH3_B"]` / `species["NH4_B"]` rewrite to assert
      on canonical keys (`species["CO2aq"]`, `species["HCO3-"]`,
      `species["NH3"]`, `species["NH4+"]`). Identify with
      `Grep "CO2_HA\|CO2_A-\|NH3_B\|NH4_B" tests/standalone/`.
      VFA tests asserting `S_ac_HA` / `S_ac_A-` are unchanged.
- [ ] [`tests/standalone/test_gas_liquid_link.py`](../../tests/standalone/test_gas_liquid_link.py)
      — tests that construct a `KineticGasLiquidLink` with
      explicit `speciation_keys={"CO2": "CO2aq"}` rewrite to
      construct a reaction set with the cross-phase CO2
      partition declaration and call
      `link.derive_speciation_keys(reaction_set)`. Identify with
      `Grep "speciation_keys=" tests/standalone/`.
- [ ] [`tests/standalone/test_reactions.py`](../../tests/standalone/test_reactions.py)
      — add coverage for `Reaction.is_cross_phase`, cross-phase
      `log_K=None` acceptance, and single-phase `log_K=None`
      rejection.
- [ ] [`tests/standalone/test_property_solvers.py`](../../tests/standalone/test_property_solvers.py)
      — the Phase 2 alpha-channel coverage tests are unchanged
      in expected output; verify they still pass with the
      simplified dispatcher (canonical-name lookups only).
- [ ] [`tests/standalone/test_ionic_strength.py`](../../tests/standalone/test_ionic_strength.py)
      (if exists; else relevant tests under
      `test_speciation.py`) — add regression coverage for the
      new suffix-parsing function:

      - `_charge_from_suffix("HCO3-") == -1`,
        `_charge_from_suffix("CO3--") == -2`,
        `_charge_from_suffix("PO4---") == -3`,
        `_charge_from_suffix("Mo7O24------") == -6`,
        `_charge_from_suffix("Mg++") == +2`,
        `_charge_from_suffix("H+") == +1`,
        `_charge_from_suffix("H2O") == 0`,
        `_charge_from_suffix("CO2aq") == 0`.
      - `_charge_from_suffix("Cation(inert)") == 0` (no
        suffix); `_CHARGE_OVERRIDES["Cation(inert)"] == +1`
        carries it.
      - `ionic_strength_from_speciation({"H+": 1e-7, "OH-":
        1e-7, "HCO3-": 0.05, "CO3--": 0.005})` returns
        `0.5 * (1e-7 + 1e-7 + 0.05 + 0.005*4) = ~0.035`
        mol/L. Demonstrates canonical-name handling without
        any z-dict lookups.
      - `ionic_strength_from_speciation({"S_ac_A-": 0.01,
        "S_pro_A-": 0.005})` returns
        `0.5 * (0.01 + 0.005) = 0.0075` mol/L. Demonstrates
        VFA-style generic-name handling via the same suffix
        rule.
      - `ionic_strength_from_speciation({"Cation(inert)": 0.02,
        "Anion(inert)": 0.02})` returns `0.5 * (0.02 + 0.02)
        = 0.02` mol/L. Demonstrates the override table.

### 10. Tests — new

- [ ] New tests in `test_reactions.py` (~3 tests):
      - `Reaction.is_cross_phase` returns False for single-phase,
        True for cross-phase stoichiometry.
      - `Reaction(kind="equilibrium", log_K=None, ...)` succeeds
        for cross-phase stoichiometry.
      - `Reaction(kind="equilibrium", log_K=None, ...)` raises
        for single-phase stoichiometry.
- [ ] New tests in `test_gas_liquid_link.py` (~3 tests):
      - `link.derive_speciation_keys(reaction_set)` populates
        `speciation_keys` from a single cross-phase reaction.
      - `derive_speciation_keys` skips reactions whose phase
        keys don't match the link's gas/liquid phase keys
        (multi-link defense).
      - `derive_speciation_keys` skips single-phase equilibrium
        reactions (no cross-phase, no contribution).
- [ ] New test in `test_speciation.py` (~2 tests):
      - `_compute_species_eq` for `bsm2_default()`-style
        n_protons=2 CO2 emits `CO2aq`/`HCO3-`/`CO3--` only (no
        `CO2_H2A`/`CO2_HA-`/`CO2_A--`).
      - `_compute_species_eq` for `from_reactions`-style
        n_protons=1 CO2 emits `CO2aq`/`HCO3-` only (no
        `CO2_HA`/`CO2_A-`).
- [ ] New test in `test_control_volume.py` (~1 test):
      - `ControlVolume.__init__` with a reaction set containing
        a cross-phase reaction and an internal
        `KineticGasLiquidLink` populates the link's
        `speciation_keys` automatically.

Estimate: +9 tests net; possibly minus a few from consolidating
canonical vs generic redundancy.

### 11. Update docs

- [ ] [`docs/phases-upcoming/CHEMISTRY_UNIFICATION.md`](CHEMISTRY_UNIFICATION.md)
      — update Leak 2's status: "Resolved 2026-05-XX in
      `chemistry-unification-3`". Note the auto-derivation
      mechanism (cross-phase equilibrium reactions →
      `link.derive_speciation_keys`).
- [ ] [`docs/phases-upcoming/CHEMISTRY_UNIFICATION_PLAN.md`](CHEMISTRY_UNIFICATION_PLAN.md)
      — Phase 3 section: replace the "Scope re-evaluation flagged"
      block with a "Shipped" note pointing at the moved checklist.
      Add a new "deferred phases" subsection (or extend the
      existing one) listing 3b as trigger-gated like
      `CONTAINER_LAYERING.md`.
- [ ] [`docs/class_diagrams.md`](../class_diagrams.md) —
      `KineticGasLiquidLink` block: drop the
      `speciation_keys: Dict[str, str]` constructor field;
      add the `derive_speciation_keys(reaction_set)` method.
      `Reaction` block: add `is_cross_phase` property.
- [ ] [`docs/architecture.md`](../architecture.md) — search for
      `speciation_keys` mentions; update the gas-liquid coupling
      paragraph to describe the cross-phase reaction declaration
      pattern.

### 12. Final test sweep + ship

- [ ] Run the full standalone test suite:
      `python -m pytest tests/standalone -q`. Baseline post
      Phase 2: 758 tests. Phase 3 estimate: 758 → ~765 (+9
      from checkpoint 10, possibly −1 to −2 from consolidation
      in checkpoint 9).
- [ ] **Critical:** confirm
      [`test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py)
      sentinels pass at `RTOL_SENTINEL = 1e-9` budget (or
      re-baseline with explanation per checkpoint 8).
- [ ] Update [README.md](README.md) priority list: mark
      `chemistry-unification-3` as shipped; flag
      `chemistry-unification-4` (Phase F / AccuracyMonitor) as
      next. Add the deferred 3b entry alongside
      `CONTAINER_LAYERING` / `RUN_HISTORY` in the trigger-gated
      list.
- [ ] Move this checklist to
      [`../phases-shipped/CHEMISTRY_UNIFICATION_3_CHECKLIST.md`](../phases-shipped/)
      with a "Shipped" status banner at the top mirroring
      Phase 2's.
- [ ] Ship via the convention in [README.md](README.md):
      `git checkout main && git merge --no-ff chemistry-unification-3
      -m "Merge chemistry-unification-3: cross-phase equilibrium
      reactions, ionic-strength dedup, gas-liquid leak 2 closed"`.
- [ ] Tag: `git tag chemistry-unification-3-shipped <commit-hash>`.
- [ ] Push: `git push && git push --tags`.
- [ ] Delete branch:
      `git branch -d chemistry-unification-3` and
      `git push origin --delete chemistry-unification-3`.

## Final test count expectation

758 → ~765 (estimate: +3 in `test_reactions` for cross-phase
acceptance, +3 in `test_gas_liquid_link` for derivation, +2 in
`test_speciation` for canonical naming, +1 in
`test_control_volume` for the attach-time wiring; minus 0–2
from consolidating redundant canonical-vs-generic assertions in
existing tests). Final number confirmed empirically at
checkpoint 12.
