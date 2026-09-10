# Chemistry Unification — Implementation Plan

## Status

Implementation plan for the work described in
[CHEMISTRY_UNIFICATION.md](CHEMISTRY_UNIFICATION.md). All design
knobs for Phase 1 are now resolved (see *Resolved decisions* under
the Phase 1 section); this plan pins the **slicing**, the
**branch/tag layout**, and the **file-level scope** for each
phase. The design doc remains the source of truth for *why*; this doc
captures *what ships in which branch*.

## Sequencing

Six branches total, shipping in order. Two precursor cleanups are
already shipped; the remaining four are chemistry unification proper.

| # | Branch | Tag | Conceptual content |
|---|---|---|---|
| 0 | `reaction-protocol-cleanup` | `reaction-protocol-cleanup-shipped` | Collapse `integration_mode` rate/delta into single rate path. Shipped 2026-05-12. See [`shipped/REACTION_PROTOCOL_CLEANUP_CHECKLIST.md`](../shipped/REACTION_PROTOCOL_CLEANUP_CHECKLIST.md). |
| 0.5 | `bsm2-reference-test` | `bsm2-reference-test-shipped` | Golden-trajectory regression test for BSM2 reference model + `SpeciationPropertySolver` forward-note docstring. Shipped 2026-05-12. See [`shipped/BSM2_REFERENCE_TEST.md`](../shipped/BSM2_REFERENCE_TEST.md). |
| 1 | `chemistry-unification-1` | `chemistry-unification-1-shipped` | Design-doc Phases A + B + C absorbed into one branch under the stateless snapshot model: `Reaction.kind`/`log_K` fields; `Species` dataclass with `atoms` + `charge`; `StoichiometryEntry` restructured to reference `Species`; `common_species.py` + `species_check.py`; `ReactionSet.partition()`; `SpeciationEngine.from_reactions()`; engine reads acid totals from phases at solve time; `chem_env` shrinks to `(t_h, T, strong-ion params, solver hints)`; `env.has_pH`. **No** staleness machinery, **no** `phase.pH`, **no** `cv.current_pH()`, **no** `_last_*` cached fields. |
| 2 | `chemistry-unification-2` | `chemistry-unification-2-shipped` | Phase D Leak 1: `PropertyResult.alphas` channel; delete `SpeciationCorrection.f_molecular`; revive or delete the dormant `_check_f_molecular_consistency` cross-check. |
| 3 | `chemistry-unification-3` | `chemistry-unification-3-shipped` | Phase E + Phase D Leak 2: `ChemistryDatabase` dataclass with `.extend()`; ship `aqueous.py` / `bioprocess_basic.py` / `anaerobic_digestion.py`; cross-phase equilibrium reactions for partitioning; `ThermoFramework` link-mismatch warnings. **Scope already trimmed by Phase 1 absorbing the `Species` type** — see this phase's section below. |
| 4 | `chemistry-unification-4` | `chemistry-unification-4-shipped` | Phase F: `AccuracyMonitor` with cheap per-step checks; `AccuracyWarning`; `VLsim.config.warnings` package-level configuration with `WarningConfig` thresholds and `VLSIM_WARNINGS` environment variable. |

> **Phase renumbering note (2026-05-13).** The plan previously listed
> five chemistry-unification phases (1–5). Resolution of design knob B3
> merged the original Phase 1 and Phase 2 into a single branch
> (`chemistry-unification-1`) — see *Resolved decisions* under Phase 1.
> Subsequent phases were renumbered: original Phase 3 → new Phase 2,
> original Phase 4 → new Phase 3, original Phase 5 → new Phase 4. The
> design doc's conceptual labels ("Phase A / B / C / D / E / F") are
> unchanged — they describe ideas, not branches.

Each branch ships through the convention in
[README.md](README.md): off `main`, merge `--no-ff`, tag, push, delete
branch.

## Discrepancies between the design doc and current code

The design doc was written before the current code state was inspected
in detail, and before the Phase 1 design knobs were resolved. Five
corrections affect implementation:

1. **`Reaction` already has an `integration_mode` field** controlling
   rate-vs-delta semantics. The design doc proposed a new `mode`
   field; that would be a second `mode`-shaped field on the same
   object. Phase 0 (`reaction-protocol-cleanup`) removes
   `integration_mode` entirely. Chemistry unification then introduces
   a single 2-valued `kind: str` field with values `"kinetic"`
   (default) and `"equilibrium"` — replacing the design doc's
   proposed `mode`.

2. **The speciation solver uses `scipy.optimize.brentq`, not
   Newton-Raphson.** The design doc claims NR throughout
   ([`acid_base.py:157`](../../src/speciation/acid_base.py#L157) is
   actually `brentq`). Same root-finder family; refactor logic
   unaffected. Doc to be corrected when CHEMISTRY_UNIFICATION.md
   moves to `shipped/`.

3. **The "dormant" cross-check is `_check_f_molecular_consistency`**
   at [`thermo_params.py:802`](../../src/chemistry/thermo_params.py#L802),
   not `_check_within_zone_f_molecular`. Phase 2 (was Phase 3) needs
   to verify whether it's actually dormant (the design doc claimed
   `spec_result` is hard-wired to `None`) before deciding
   revive-vs-delete.

4. **`Species` type lands in Phase 1, not Phase 3.** The design doc
   places `Species` (with `formula`, `charge`) inside the
   `ChemistryDatabase` infrastructure of Phase E (this plan's Phase 3,
   was Phase 4). Resolution of design knob A2 promoted `Species` into
   Phase 1 instead, with `atoms: Dict[str, int]`, `charge: int`, and
   `MW: float` as fields. Phase 3 retains only the `ChemistryDatabase`
   wrapper, cross-phase reactions, and `ThermoFramework` — see the
   Phase 3 scope re-evaluation note below.

5. **No staleness machinery; stateless snapshot model.** The design
   doc proposes `LiquidPhase._equilibrium_species`,
   `LiquidPhase._equilibrium_stale`, `phase.pH`, `cv.current_pH()`,
   `StaleEquilibriumError`, and the post-solve writeback of H⁺ to
   `phase.n_mol`. Resolution of design knob B2 eliminates all of these.
   Equilibrium-output species (H⁺, OH⁻, HCO₃⁻, CO₃²⁻, NH₄⁺) do not
   live in `phase.n_mol` — they are pure property outputs returned via
   `PropertyResult` on the `AdvanceResult`. Rate functions read pH
   from `ReactionEnvironment` (which carries property results from the
   current step's speciation solve); users between steps read pH from
   the most recent `AdvanceResult`. Warm-start state lives inside the
   `SpeciationEngine`, not on the phase. This aligns the equilibrium-
   property read with the rest of the post-Phase-7 architecture, which
   was already stateless (no `_last_properties` cache, no `cv.pH`
   accessor).

## Phase-by-phase scope

### Phase 1 — `chemistry-unification-1`: Reaction declaration unified + delivery cleanup

**Conceptual content (from design doc):** Phases A + B + C together
under the stateless snapshot model. Equilibrium reactions become
first-class declarations (Phases A+B), and `chem_env` shrinks to its
minimal contract in the same branch (Phase C absorbed via B3
resolution). This is the largest single phase; original Phases 1 and
2 of the plan merged into one branch.

**File-level scope (high-level — refine into a checklist when work
begins):**

- **New module:**
  [`src/chemistry/species.py`](../../src/chemistry/species.py) —
  introduces the `Species` dataclass (frozen, `eq=True`) carrying the
  intrinsic identity of a chemical: `id: str`, `atoms: Dict[str, int]`,
  `charge: int = 0`, `MW: float = 0.0`. Atomic composition and charge
  live here, not on each `StoichiometryEntry`. (See *Resolved
  decisions* below for the A2 rationale.)
- **New module:**
  [`src/chemistry/common_species.py`](../../src/chemistry/common_species.py) —
  module-level declarations for universal inorganic aqueous species
  (`H_plus`, `H2O`, `OH`, `CO2`, `HCO3`, `CO3`, `NH3`, `NH4`, …) so
  every model can import the same `Species` objects rather than
  redeclaring them. Organics that are model-specific (acetate,
  propionate, biomass populations, …) live in their model files.
- **New module:**
  [`src/chemistry/species_check.py`](../../src/chemistry/species_check.py) —
  `check_species_consistency(reactions, *, soft_conflicts="warn")`
  utility. Detects two failure modes across composed reaction sets:
  *hard conflict* (same `id`, different `atoms`/`charge`) raises;
  *soft conflict* (distinct `Species` objects with equal data — user
  redeclared instead of imported) warns by default. Wire into
  `ControlVolume.__init__` to run automatically after `reaction_model`
  is attached.
- [`src/reactions/reaction.py`](../../src/reactions/reaction.py) — add
  `kind: str = "kinetic"` field; `log_K: Optional[float] = None`
  (required iff `kind="equilibrium"`); make `rate_fn` optional
  (forbidden iff `kind="equilibrium"`). Validate combinations in
  `__init__`.
- [`src/reactions/stoichiometry.py`](../../src/reactions/stoichiometry.py) —
  `StoichiometryEntry` changes from `(species_id, phase, coefficient,
  atoms, MW)` to `(species: Species, phase, coefficient)` — atoms / MW
  / charge move onto the referenced `Species`. `validate_balance`
  iterates over `entry.species.atoms` for element conservation; gains
  an opt-in `check_charge: bool = False` parameter that iterates over
  `entry.species.charge` as a peer conservation pass. One mechanism,
  two laws — no pseudo-element hack, no parallel validator.
- [`src/reactions/reaction_set.py`](../../src/reactions/reaction_set.py) —
  add `partition()` method that splits a `ReactionSet` into kinetic
  and equilibrium subsets. Used by the builder pattern (see B1
  resolution).
- [`src/speciation/engine.py`](../../src/speciation/engine.py) — add
  `SpeciationEngine.from_reactions(equilibrium_rxns, activity_model,
  T_K, ...)` factory that builds the engine's matrix from declared
  equilibrium reactions. Engine reads species charges off
  `entry.species.charge` for activity / Davies coefficients. The
  engine **also reads acid totals (`CT_TIC`, `CT_NH_T`, `acid_totals`)
  directly from the `phases` argument** at each `solve()` call — no
  longer via `chem_env`. Internal warm-start state (`_last_logH`)
  lives on the engine; not on any phase. Activity model configuration
  unchanged — stays at engine construction; Phase 3 promotes it to
  `ThermoFramework`.
- [`src/core/speciation_solver.py`](../../src/core/speciation_solver.py) —
  `SpeciationPropertySolver` becomes a thin wrapper around the engine.
  Removes the legacy `equilibrium_set: EquilibriumSet` constructor
  argument (replaced by the engine being constructed from declared
  reactions upstream).
- [`src/core/property_solver.py`](../../src/core/property_solver.py) —
  `chem_env`'s contract shrinks. Document the new minimal contract:
  `(t_h, T_K, strong_ion params, optional solver hints)`. The
  `acid_totals` / `CT_P` / `CT_NH_T` / `acid_pKas` keys are removed.
- [`src/core/control_volume.py`](../../src/core/control_volume.py) —
  runs `check_species_consistency` after `reaction_model` attach.
  Property solver ordering unchanged (already runs first).
- [`src/reactions/environment.py`](../../src/reactions/environment.py) —
  `env.has_pH` property added (introspection-based — returns `True`
  iff the CV's property results include a non-`None` speciation pH).
- New module: error type `SpeciesConflictError(ValueError)` for hard
  conflicts detected by `species_check`. **Not added:**
  `PHUndefinedError` and `StaleEquilibriumError` — both are obviated
  by the stateless snapshot model (see B2 resolution).
- **Builder migration (B1 resolution):** model builders explicitly
  call `partition()` on their reaction set and pass the equilibrium
  subset to `SpeciationEngine.from_reactions(...)`. The pattern:
  ```python
  kinetic_rxns, equilibrium_rxns = rxn_set.partition()
  engine = SpeciationEngine.from_reactions(
      equilibrium_rxns, activity_model="davies", T_K=308.15,
  )
  solver = SpeciationPropertySolver(engine=engine)
  cv = ControlVolume(
      phases=...,
      reaction_model=rxn_set,        # full set; kinetics use rate_fn
      property_solvers=[solver],
  )
  ```
  Partition lives at the builder/wiring site, not hidden inside the
  solver or the engine.
- **Existing model migration:** BSM2 and ADM1 builders rewrite their
  `SPECIES_BSM2` / `SPECIES_ADM1` dicts as module-level `Species`
  declarations; the `_make_entries` helper that currently denormalises
  atoms into entries becomes a thin pass-through. BSM2's
  `make_bsm2_chem_env_fn` shrinks to populate only `(t_h, T_K,
  strong-ion params)` — no totals, no `acid_pKas`. BSM2's
  `make_bsm2_callback` is **deleted** — H⁺ writeback is gone since
  H⁺ no longer lives in `phase.n_mol`. The BSM2 reference golden test
  ([`tests/standalone/test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py))
  drops its `callback(...)` line and adjusts `chem_env_fn` kwargs;
  sentinel goldens pass bit-for-bit because the changes are
  structural, not numerical.
- **`apply_flux` is *not* a correctness boundary** under the
  stateless snapshot model. Existing direct-dict-access paths to
  `phase.n_mol` (e.g. [`models/vlmodels/fermenter/unit.py:2178`](../../models/vlmodels/fermenter/unit.py#L2178))
  remain valid — no audit / migration step needed.

**Test impact:** `test_reactions.py` (44 tests, +N for kind/log_K and
Species refactor), `test_speciation.py` (~30 tests, substantially
rewritten — migration to declared equilibrium reactions and dropping
`chem_env={"acid_totals": ...}` constructions), `test_property_solvers.py`
(16 tests, contract update for the new `chem_env` shape),
`test_cv_advance.py` (47 tests, ordering + species-check call site +
chem_env contract), `test_bsm2_reference.py` (6 tests, goldens must
still pass bit-for-bit — pure structural change). New: `test_species.py`
(Species identity, atom/charge access), `test_species_check.py`
(conflict detection). Total touched: ~140–170 tests, plus net
additions.

**Resolved decisions:**

- **A2 — Charge balance representation: charge on `Species`, not on
  entries or as a pseudo-element.** `Species` carries both `atoms` and
  `charge`. `validate_balance` runs element checks via `entry.species.
  atoms` and an opt-in charge check via `entry.species.charge` — one
  iteration pattern, two conservation laws. Rationale: the existing
  per-entry `atoms` dict is already a duplication footgun (the same
  species in 10 reactions has its atomic composition declared 10
  times; the BSM2 code works around this with the ad-hoc `SPECIES_BSM2`
  dict). Promoting both atoms and charge to a `Species` type fixes the
  footgun and lands charge balance cleanly in the same change. The
  pseudo-element `"Q"` option was a hack; the parallel `charge: int`
  field on entries would have doubled the conservation infrastructure.
  This expands Phase 1 scope vs the original "additive only" framing
  but is the self-consistent answer the code is already reaching for.
  Discussion: 2026-05-12, resolved with user.

- **Species declaration/sharing pattern: module-level Python objects
  with an opt-in consistency utility.** Species are declared as
  module-level `Species` objects (frozen dataclass with `eq=True` —
  value equality on fields). Cross-file sharing uses Python's normal
  import system; object identity is preserved automatically when a
  shared species is imported. The `check_species_consistency` utility
  validates composed reaction sets at CV construction time, catching
  both *hard conflicts* (same id, different data → error) and *soft
  conflicts* (different objects, same data → warning, signals
  "redeclared instead of imported"). The alternative — a
  `SpeciesRegistry` with mandatory string-keyed lookup — was rejected
  because it reintroduces stringly-typed references that the
  `Species`-as-object design removed, and forces every reaction
  construction site to thread a registry. Discussion: 2026-05-12,
  resolved with user.

- **B1 — Translator location: builder partitions explicitly; engine
  consumes equilibrium reactions at construction; solver is a thin
  wrapper.** The model builder calls `rxn_set.partition()` and passes
  the equilibrium subset to `SpeciationEngine.from_reactions(...)`.
  `SpeciationPropertySolver` is a thin adapter around the engine; it
  does not partition internally. Rationale: keeps each component
  single-responsibility (engine solves equilibria; solver adapts to
  `PropertySolver` protocol; builder wires components together);
  partition is visible at the wiring site rather than hidden inside
  the solver; engine has a clean `from_reactions` factory usable
  standalone. The alternatives (solver partitions internally, engine
  reads equilibrium reactions from `chem_env`, `PropertySolver`
  gains an `on_cv_attach` lifecycle hook, etc.) were considered and
  rejected — see discussion 2026-05-13.

- **B2 — Staleness scope: eliminated; adopt stateless snapshot
  model.** No `_equilibrium_species` field on `LiquidPhase`. No
  `_equilibrium_stale` flag. No `phase.pH` property. No
  `cv.current_pH()` accessor. No `StaleEquilibriumError`. No H⁺
  writeback to `phase.n_mol` after speciation solves. Equilibrium-
  output species (H⁺, OH⁻, HCO₃⁻, CO₃²⁻, NH₄⁺) do not live in
  `phase.n_mol` — they are pure property outputs returned via
  `PropertyResult` on the `AdvanceResult`. Rate functions read pH and
  other equilibrium-output species through `ReactionEnvironment`
  (which carries the current step's property results); users between
  steps read pH from the most recent `AdvanceResult`. Warm-start state
  lives inside `SpeciationEngine`, not on any phase. Rationale: the
  staleness machinery was solving a problem that disappears under
  strict snapshot semantics, which the rest of the post-Phase-7
  architecture already follows. Aligns the equilibrium-property read
  with every other read in `advance()`. This also requires committing
  that kinetic reactions do not have H⁺ (or other equilibrium-output
  species) directly in their stoichiometry — true for BSM2 and ADM1
  today; any future reaction needing explicit H⁺ stoichiometry would
  need to be reframed in totals (deliberate constraint, not silent
  loss). Discussion: 2026-05-13, resolved with user.

- **B3 — Intermediate-state policy: clean break; merge Phase 1 and
  Phase 2.** Engine reads acid totals directly from `phases` at
  `solve()` time. `chem_env` shrinks immediately to its minimal
  contract. No shim, no parallel data paths, no separate cleanup
  phase. Rationale: B2's resolution shrank Phase 1's scope enough
  that absorbing the original Phase 2's work brings Phase 1 back to
  a reasonable size while eliminating an awkward intermediate state.
  Consistent with the project's "no shims, clean breaks" rule applied
  throughout the seven shipped phases of the CV refactor. The
  bigger-PR cost is acceptable given no legacy users constrain
  migration order. Discussion: 2026-05-13, resolved with user.

### Phase 2 — `chemistry-unification-2`: Gas-liquid leak 1 (alphas channel)

**Conceptual content:** Phase D Leak 1 from design doc (was Phase 3 in
the previous plan numbering). Plus a decision on the dormant
`_check_f_molecular_consistency`.

**File-level scope:**

- [`src/core/property_solver.py`](../../src/core/property_solver.py) —
  add `alphas: Dict[str, float]` field to `PropertyResult`.
- [`src/speciation/engine.py`](../../src/speciation/engine.py) —
  populate `alphas` in `solve()` output (CO2aq fraction, NH3 fraction,
  etc.).
- [`src/core/gas_liquid_link.py`](../../src/core/gas_liquid_link.py) —
  `KineticGasLiquidLink` reads
  `property_results["speciation"].alphas[...]` instead of calling
  `correction.f_molecular(pH, T_K)`. Delete
  `SpeciationCorrection.f_molecular` and probably
  `SpeciationCorrection` itself if no other consumer.
- [`src/chemistry/thermo_params.py`](../../src/chemistry/thermo_params.py) —
  decide: revive `_check_f_molecular_consistency` to read from
  `alphas`, or delete the dead-code function. Verify dormancy first
  (the design doc claim that `spec_result` is hard-wired to `None`
  needs confirmation).

**Test impact:** Modest — `test_gas_liquid_link.py` (27 tests,
several rewritten to read `alphas` instead of computing
`f_molecular`).

### Phase 3 — `chemistry-unification-3`: Database packaging + cross-phase reactions

> **Scope re-evaluation flagged 2026-05-12; refreshed 2026-05-15
> post Phase 2 shipping.** Phase 1's resolution of design knob A2
> (Species dataclass with charge + atoms, module-level declarations,
> `common_species.py`, `species_check.py`) absorbed several items
> originally scoped here: the `Species` type itself, the universal
> species declarations, the conflict-detection mechanism, and most
> call-site migration. Phase 2 (shipped 2026-05-15) added the
> `PropertyResult.alphas` channel and collapsed the gas-liquid
> link's two correction paths into a single alpha read, closing
> Leak 1. What's left for this phase is the *packaging layer* on top
> — `ChemistryDatabase` as a frozen dataclass bundling `(species_set,
> reaction_set, ThermoFramework)` with `.extend()` composition, the
> `chemistry_db` field on `ControlVolume`, cross-phase
> equilibrium-reaction validation, the `ThermoFramework`
> link-mismatch warning, and curated stock-database modules.
>
> **Scope decisions resolved 2026-05-15** (see
> [`CHEMISTRY_UNIFICATION_3_CHECKLIST.md`](CHEMISTRY_UNIFICATION_3_CHECKLIST.md)
> "Resolved decisions"). Phase 3 ships:
>
> - **Ionic-strength dedup at source.** Canonical-name emission
>   in `_compute_species_eq` for recognised acids (CO₂,
>   phosphate, sulfate, NH₃/NH₄⁺); generic emission for VFAs.
>   `ionic_strength_from_speciation` rewritten around
>   charge-suffix parsing. Unmasks a latent NH4+/I bug
>   introduced by Phase 1 (BSM2 was undercounting ionic
>   strength).
> - **Cross-phase equilibrium reaction infrastructure.**
>   `Reaction.is_cross_phase`, `log_K=None` for cross-phase
>   equilibria, `SpeciationEngine.from_reactions` filters
>   cross-phase silently,
>   `KineticGasLiquidLink.derive_speciation_keys` is wired
>   into `ControlVolume.__init__`.
>
> Phase 3 **does not** migrate BSM2/ADM1 builders to declare
> cross-phase reactions. The reason is a design dependency on
> the alpha-key naming convention surfaced during
> implementation 2026-05-15: today's
> `PropertyResult.alphas` keys (`"CO2aq"`, `"S_ac_HA"`) don't
> match the Species IDs declared in reaction stoichiometry
> (`"CO2"`, `"S_ac"`). The clean fix is "engine emits species
> using Species IDs directly" (no synthesized
> `_HA`/`_A-`/`_aq` decorators), which requires
> `EquilibriumDef` to carry Species references and touches the
> legacy `bsm2_default()` factory path — closer in shape to
> Phase 3b's `ChemistryDatabase` rollout than to Phase 3a's
> cross-phase wiring.
>
> **Phase 3b mandate (formalized 2026-05-15):** Phase 3b
> ships:
>
> - `ChemistryDatabase` dataclass with `.extend()` composition.
> - Stock-database modules (`aqueous.py`, `bioprocess_basic.py`,
>   `anaerobic_digestion.py`).
> - `ThermoFramework` absorbing the engine-level activity-knob
>   configuration and `ThermodynamicConfig`'s remaining
>   responsibilities.
> - **Unified Species-ID emission convention.** The engine
>   emits species under their declared Species IDs (no more
>   `_HA`/`_A-`/`_BH+`/`_B`/`aq` decorators); alpha-keys are
>   Species IDs; `KineticGasLiquidLink.derive_speciation_keys`
>   becomes the trivial `speciation_keys[gas_id] = liq_id`
>   with no translation table or charge-class lookup. The
>   small `_CANONICAL_NAMES` table that Phase 3 leaves in
>   `acid_base.py` goes away.
> - **BSM2/ADM1 builder migration** to declare cross-phase
>   reactions, drop explicit `speciation_keys=` kwargs, and
>   drop the `FermenterBuilder.speciation_correction(...)`
>   builder method and `TransferConfig.speciation_keys` field.
> - `core/links.py` `ThermoFramework` mismatch warning.
> - **Demo acceptance criterion (added 2026-05-28, path updated
>   2026-05-29).** Ship `demos/model_api/chemistry/chemistry_database.py`
>   demonstrating the create → extend → import lifecycle with one
>   of the stock-database modules (`aqueous.py` /
>   `bioprocess_basic.py` / `anaerobic_digestion.py`). This is the
>   "Q7 demo" deferred from the 2026-05-28 exploration session — it
>   has no home until 3b's `ChemistryDatabase` type exists, so it
>   lives here as part of the acceptance criteria. (Original scope
>   used `demos/api/`; the [DEMO_RESTRUCTURE.md](DEMO_RESTRUCTURE.md)
>   ship on 2026-05-29 renamed that subtree to `demos/model_api/`
>   and nested chemistry-defs under `chemistry/`.)
>
> **Interaction with
> [PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md](PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md)
> — open question (added 2026-05-28).** The `ThermoFramework`
> introduced here absorbs activity-model configuration. The
> Peng–Robinson fugacity use case named in
> PROPERTY_CALCULATORS_PHASE_AGNOSTIC sits at the seam: does
> high-pressure gas-phase fugacity configuration live on
> `ThermoFramework` (this note) or on a gas-phase
> `PropertyCalculator` (the other note)? Neither pre-commits.
> A high-pressure biogas storage workload would force the
> question; whichever phase is triggered first should resolve
> it explicitly rather than picking one path silently.
>
> 3b remains trigger-gated (no concrete pull for runtime
> chemistry composition today), but its scope is now locked.
> When triggered, the work is the unified-emission refactor +
> the database packaging on top, not a re-litigation of the
> naming-convention question.

**Conceptual content:** Phase E from design doc + Phase D Leak 2
(cross-phase equilibrium reactions for gas-liquid partitioning). Was
Phase 4 in the previous plan numbering.

**File-level scope (refine after Phase 1 lands — see re-sizing note
above):**

- New package: `src/chemistry/databases/` with
  `aqueous.py`, `bioprocess_basic.py`, `anaerobic_digestion.py`.
- New module: `src/chemistry/database.py` with `ChemistryDatabase`
  frozen dataclass + `.extend()` composition. **Wraps the `Species`
  type that arrives in Phase 1** — does not introduce it.
- New module: `src/thermo/framework.py` with `ThermoFramework` frozen
  dataclass holding `activity_model`, reference state, standard
  conditions. Replaces the engine-level `use_activity` / `activity_model`
  string-keyed configuration introduced in Phase 1 (engine reads the
  framework off `cv.chemistry_db.thermo`).
- [`src/reactions/reaction.py`](../../src/reactions/reaction.py) —
  extend stoichiometry validation to permit cross-phase equilibrium
  reactions (e.g. `[(CO2_gas_species, "gas", -1), (CO2_species, "liquid", +1)]`).
- [`src/core/gas_liquid_link.py`](../../src/core/gas_liquid_link.py) —
  derive `speciation_keys` automatically from declared cross-phase
  equilibrium reactions; remove the manual `speciation_keys` argument.
- [`src/core/control_volume.py`](../../src/core/control_volume.py) —
  CV gains `chemistry_db: ChemistryDatabase` field; `thermo`
  accessible via `cv.chemistry_db.thermo`.
- [`src/core/links.py`](../../src/core/links.py) (or wherever
  `AdvectiveLink` / `DiffusiveLink` live) — add construction-time
  warning if linked CVs have different `ThermoFramework` instances.
- Move hardcoded data from
  [`src/speciation/acid_base.py:173–179`](../../src/speciation/acid_base.py#L173)
  (pKa defaults) and
  [`src/chemistry/equilibria.py`](../../src/chemistry/equilibria.py)
  (EquilibriumSet presets) into the new database modules.

**Phase 2 follow-ups (added to scope 2026-05-15):**

- **Deduplicate `ionic_strength_from_speciation`'s two counting paths.**
  [`src/speciation/activity.py`](../../src/speciation/activity.py#L52)
  sums charges via a hardcoded `z`-dict for canonical names
  (`HCO3-`, `NH4+`, …) AND a fallback loop over keys ending in
  `_A-`. When both canonical and generalised names are present in
  `sp`, the species is double-counted. Phase 2 worked around this
  by emitting canonical names *into the alphas channel only* (not
  into the species dict). Phase 3 should pick one naming
  convention (canonical) and have `_compute_species_eq` emit only
  those keys — see the warning block on the function itself. See
  [`docs/shipped/CHEMISTRY_UNIFICATION_2_CHECKLIST.md`](../shipped/CHEMISTRY_UNIFICATION_2_CHECKLIST.md)
  checkpoint 2 deviation note for the BSM2 sentinel drift this
  caused when first investigated (2.1e-5 on S_ac).
- **Auto-derive `speciation_keys` from cross-phase equilibrium
  reactions (Leak 2 closure).** Phase 2 kept the manual
  `link.speciation_keys = {"CO2": "CO2aq", …}` mapping; Phase 3's
  cross-phase reaction support means each gas↔liquid equilibrium
  reaction implies a `speciation_keys` entry, and the manual
  argument can be removed. Already listed under the `gas_liquid_link.py`
  bullet above.

**Test impact:** Significant — every model construction site changes
from passing pKas explicitly to passing a `chemistry_db`. Migration
of BSM2 / ADM1 builders. Substantial.

### Phase 4 — `chemistry-unification-4`: Accuracy monitoring

> **Shipped 2026-05-16.** Tag `chemistry-unification-4-shipped`.
> Implementation log:
> [../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md).
> Net: +43 standalone tests (778 → 821); zero BSM2 sentinel drift
> (Phase 4 is purely additive on numerics). Five of the six cheap
> checks specified in Phase F are wired (pH change, Newton iters,
> charge residual, scipy step rejections, ionic-strength regime);
> the `dt` vs τ_min check has a no-op hook with a deferred
> Jacobian-free estimate. The splitting-error proxy
> (pH-pre-vs-pH-post on the same step) is documented as future
> work — it needs solver-internal state that wasn't justified
> for Phase 4's scope. `cv.run_with_diagnostics(...)` was
> deferred entirely (no stub, no half-implementation); the
> docs/solvers.md "Accuracy monitoring" section carries the
> forward note.

## Sequencing of doc updates

- **When Phase 0 ships:** ✅ Done. `reaction-protocol-cleanup` and
  `bsm2-reference-test` are both in
  [`shipped/`](../shipped/) with status banners.
- **When Phase 1 ships:** update [class_diagrams.md](../class_diagrams.md)
  layer 3 — equilibrium reactions become a `Reaction` `kind`; the
  parallel `PropertySolver ←→ ReactionModel` framing collapses.
  Substantial rewrite, not a search-and-replace (per the design doc's
  "Class-Diagram Impact" section). Note that the deferred-refactor
  block at lines 305–323 of the current `class_diagrams.md` is the
  thing being collapsed — read it as a guide to what changes. Also
  update the `SpeciationPropertySolver` docstring forward note (in
  [src/core/speciation_solver.py:51-66](../../src/core/speciation_solver.py#L51-L66))
  to "absorbed" status since the absorption has now occurred.
- **When Phase 3 ships:** ✅ Done 2026-05-15. No architecture.md
  update was needed — Phase 3 deferred `ChemistryDatabase` to
  Phase 3b. The user-facing chemistry config surface remains
  module-level `Species` declarations + per-model reaction
  builders (BSM2's `build_bsm2_reactions`, ADM1's
  `build_adm1_reactions`) until 3b triggers. See
  [`../shipped/CHEMISTRY_UNIFICATION_3_CHECKLIST.md`](../shipped/CHEMISTRY_UNIFICATION_3_CHECKLIST.md)
  for what actually shipped (canonical-naming dedup +
  cross-phase reaction infrastructure).
- **When Phase 3b ships:** ✅ Done 2026-06-03. Tag
  `chemistry-unification-3b-shipped`. `ChemistryDatabase` and
  `ThermoFramework` are now the user-facing chemistry config surface
  (see [`CHEMISTRY_UNIFICATION_3B_CHECKLIST.md`](CHEMISTRY_UNIFICATION_3B_CHECKLIST.md)).
  `_CANONICAL_NAMES` deleted; unified Species-ID emission
  (`CO2` replaces `CO2aq`); ghost-key CO₂ transfer bug fixed.
  See
  [architecture.md](../architecture.md) for updated descriptions.
- **When Phase 4 ships:** ✅ Done 2026-05-16. Tags
  `chemistry-unification-1-shipped` through
  `chemistry-unification-4-shipped` + `chemistry-unification-3b-shipped`
  (2026-06-03) are the permanent checkpoints. The design doc and this
  plan move to [`shipped/`](../shipped/) now that all
  scheduled phases including 3b are complete.

## What this plan does not commit to

- **Specific file:line edits within each phase.** Those go in
  per-phase checklist files (e.g. `CHEMISTRY_UNIFICATION_1_CHECKLIST.md`)
  added when work on that phase begins.
- **Final scope of Phase 3.** The re-evaluation note in that phase's
  section flags possible compression or split into Phase 3a/3b after
  Phase 1 lands. Decision deferred to when work on Phase 3 begins.
