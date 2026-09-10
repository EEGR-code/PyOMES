# Chemistry Unification 1 — Checklist

> **Status: Shipped 2026-05-13** — all 13 checkpoints landed.
> 750 standalone tests passing (up from 714 baseline: +37 net).
> BSM2 reference sentinels re-baselined to remove a hidden one-step
> pH lag from the deleted `make_bsm2_callback` — see commit
> `chemistry-unification-1: migrate BSM2 + ADM1 builders` for the
> drift breakdown. Tag: `chemistry-unification-1-shipped`.

Working checklist for `chemistry-unification-1`. Conceptual framing is
in [../upcoming/CHEMISTRY_UNIFICATION.md](../upcoming/CHEMISTRY_UNIFICATION.md);
the branch slicing and resolved knobs are in
[../upcoming/CHEMISTRY_UNIFICATION_PLAN.md](../upcoming/CHEMISTRY_UNIFICATION_PLAN.md).
This file is the implementation plan: file-level edits, ordered
checkpoints, and sanity checks. Modelled on
[REACTION_PROTOCOL_CLEANUP_CHECKLIST.md](REACTION_PROTOCOL_CLEANUP_CHECKLIST.md).

## Goal

Collapse the kinetic / equilibrium declaration split. After this
branch ships:

- Equilibrium reactions (water dissociation, carbonate, ammonium, VFAs)
  are first-class [`Reaction`](../../src/reactions/reaction.py) objects
  with `kind="equilibrium"` and `log_K`. Kinetic reactions keep `kind="kinetic"`
  (the default) and a `rate_fn`.
- A new [`Species`](../../src/chemistry/species.py) dataclass owns
  `id`, `atoms`, `charge`, `MW`. `StoichiometryEntry` references a
  `Species` instead of duplicating `species_id`/`atoms`/`MW` per
  entry. `validate_balance` iterates over `entry.species.atoms`; charge
  balance is opt-in via the same mechanism.
- Module-level [`common_species.py`](../../src/chemistry/common_species.py)
  declares the universal inorganic aqueous species (H⁺, H₂O, OH⁻, CO₂,
  HCO₃⁻, CO₃²⁻, NH₃, NH₄⁺, …) so every model imports the same `Species`
  objects. A
  [`species_check.py`](../../src/chemistry/species_check.py) utility
  detects hard conflicts (raise) and soft conflicts (warn).
- [`SpeciationEngine.from_reactions(...)`](../../src/speciation/engine.py)
  builds the engine's matrix from declared equilibrium reactions; the
  engine reads acid totals (`CT_TIC`, `CT_NH_T`, `acid_totals`)
  directly from `phases` at each `solve()` call.
- [`SpeciationPropertySolver`](../../src/core/speciation_solver.py) is
  a thin wrapper around the engine. The legacy `equilibrium_set`
  constructor argument is gone.
- `chem_env`'s contract shrinks to `(t_h, T_K, strong-ion params,
  optional solver hints)`. `acid_totals`, `acid_pKas`, `CT_P`,
  `CT_NH_T`, `CT_TIC` are removed.
- BSM2 + ADM1 builders are migrated. `SPECIES_BSM2` becomes
  module-level `Species` declarations. `make_bsm2_chem_env_fn` shrinks
  to populate only `(t_h, T_K, strong_kwargs)`. `make_bsm2_callback`
  is **deleted** (no H⁺ writeback under the stateless snapshot model).
- [`ReactionEnvironment`](../../src/reactions/environment.py) gains
  `env.has_pH` (introspection — `True` iff `pH is not None`).
- `ControlVolume.__init__` runs `check_species_consistency` after
  `reaction_model` is attached.

## Out of scope

Everything reserved for later chemistry-unification phases:

- **Phase 2:** `PropertyResult.alphas`, deletion of
  `SpeciationCorrection.f_molecular`, decision on the dormant
  `_check_f_molecular_consistency`.
- **Phase 3:** `ChemistryDatabase` dataclass with `.extend()`,
  `cv.chemistry_db` field, stock database modules (`aqueous.py`,
  `bioprocess_basic.py`, `anaerobic_digestion.py`), `ThermoFramework`
  (activity model is still configured at engine construction in this
  branch), cross-phase equilibrium reactions for gas-liquid
  partitioning, link-mismatch warnings.
- **Phase 4:** `AccuracyMonitor`, `AccuracyWarning`, `WarningConfig`.
- Rewriting `_petersen_to_molar` or BSM2's reaction *content* (the
  reactions themselves stay byte-equivalent — only their declaration
  surface changes).
- Migrating `SpeciationCorrection` on the gas-liquid link off
  `acid_pKas` (Phase 2/3 work).
- Touching `SpeciationCorrection` consumers in
  [gas_liquid_link.py](../../src/core/gas_liquid_link.py) — they keep
  reading their existing inputs.
- The HPLC, membrane, multi-CV, and headspace test modules — they do
  not touch speciation or `chem_env`.

## Resolved decisions

Pinned in [CHEMISTRY_UNIFICATION_PLAN.md](CHEMISTRY_UNIFICATION_PLAN.md)
under "Phase 1 — Resolved decisions". Recapped here for editing
context:

- **A2 — charge on `Species`, not on entries.** `Species` carries
  `atoms` *and* `charge`. `validate_balance` runs element checks via
  `entry.species.atoms`, with an opt-in `check_charge: bool = False`
  parameter that adds `entry.species.charge` as a peer conservation
  pass.
- **Sharing pattern — module-level objects + opt-in consistency
  utility.** No `SpeciesRegistry`. Cross-file sharing via Python
  imports. `check_species_consistency(reactions, *,
  soft_conflicts="warn")` runs at CV construction.
- **B1 — builder partitions explicitly.** `rxn_set.partition()` →
  `(kinetic, equilibrium)`. Builder passes the equilibrium subset to
  `SpeciationEngine.from_reactions(...)`. Engine has a clean factory.
  Solver is a thin adapter, no internal partitioning.
- **B2 — stateless snapshot model.** No
  `_equilibrium_species`/`_equilibrium_stale` flags. No `phase.pH`,
  no `cv.current_pH()`, no `StaleEquilibriumError`. No H⁺ writeback to
  `phase.n_mol`. Equilibrium-output species are pure property outputs
  on `PropertyResult`. Warm-start state lives inside
  `SpeciationEngine`. **Constraint:** kinetic reactions must not
  declare H⁺ (or other equilibrium-output species) in their
  stoichiometry. True for BSM2 / ADM1 today; verify during checkpoint
  9.
- **B3 — clean break, no shim.** `chem_env` shrinks immediately.
  Engine reads totals from `phases` directly. No parallel data path,
  no separate cleanup phase.

- **`SpeciationPropertySolver` stays as a stable adapter** (resolved
  2026-05-13). Phase 1 makes it a ~20-line wrapper around
  `SpeciationEngine`; it is **not** scheduled for later removal. The
  engine has a richer surface than the `PropertySolver` protocol
  (warm-start cache, `n_solve_calls`, `solve(**kwargs)` used directly
  by `test_speciation.py`) — welding it to the protocol shape would
  cost flexibility the test suite already uses. `SpeciationPropertySolver`
  becomes the canonical example of a protocol adapter for the
  property-flavoured solvers Phase 7 anticipated (viscosity, density).
  Checkpoint 12 rewrites the forward note in `speciation_solver.py`
  to "stable adapter", not "absorbed".

## Checkpoints

The checkpoints are ordered so that each one leaves the test suite in
a runnable state. New types land first (1–2), then the protocol
changes that consume them (3–6), then the model migrations that
exercise the new surface (7–9), then test rewrites + final sweep.

### 1. Add `Species` + `common_species` + `species_check`

- [ ] New module: [`src/chemistry/species.py`](../../src/chemistry/species.py).
      Define `Species` as a `@dataclass(frozen=True, eq=True)` with
      fields `id: str`, `atoms: Dict[str, int]`,
      `charge: int = 0`, `MW: float = 0.0`. `__post_init__` should
      coerce `atoms` to a frozen mapping (e.g. via
      `object.__setattr__(self, "atoms", types.MappingProxyType(dict(atoms)))`)
      so the dataclass remains hashable. Add a
      `SpeciesConflictError(ValueError)` here too.
- [ ] New module:
      [`src/chemistry/common_species.py`](../../src/chemistry/common_species.py).
      Module-level `Species` declarations for: `H_plus`, `OH_minus`,
      `H2O`, `CO2`, `HCO3_minus`, `CO3_2minus`, `NH3`, `NH4_plus`. IDs
      should be the literal strings the speciation engine already uses
      (`"H+"`, `"OH-"`, `"H2O"`, `"CO2"`, `"HCO3-"`, `"CO3--"`,
      `"NH3"`, `"NH4+"`) so engine output keys do not change.
- [ ] New module:
      [`src/chemistry/species_check.py`](../../src/chemistry/species_check.py).
      Implement
      `check_species_consistency(reactions, *, soft_conflicts="warn")`:
      walk every `entry.species` in every reaction, group by `id`,
      raise `SpeciesConflictError` on hard conflicts (same id,
      different atoms or charge), warn on soft conflicts (distinct
      objects, equal data) when `soft_conflicts == "warn"`. Accept
      `soft_conflicts in ("warn", "raise", "ignore")`.
- [ ] [`src/chemistry/__init__.py`](../../src/chemistry/__init__.py)
      — re-export `Species`, `SpeciesConflictError`, and
      `check_species_consistency`. Re-export the `common_species`
      module as a namespace (so callers can write
      `from VLsim.chemistry import common_species as cs`).

Sanity check: `python -c "from VLsim.chemistry import Species,
common_species, check_species_consistency, SpeciesConflictError"`
imports cleanly.

### 2. Restructure `StoichiometryEntry` to reference `Species`

- [ ] [`src/reactions/stoichiometry.py:43-72`](../../src/reactions/stoichiometry.py#L43)
      — replace the `(species_id, phase, coefficient, atoms, MW)`
      fields with `(species: Species, phase: str, coefficient: float)`.
      Remove `atoms` and `MW` from the dataclass (they move onto
      `Species`).
- [ ] [`src/reactions/stoichiometry.py:75-124`](../../src/reactions/stoichiometry.py#L75)
      — `validate_balance(entries, elements=("C","H","O"), atol=1e-10,
      *, check_charge=False)`: replace `e.atoms.get(elem, 0.0)` with
      `e.species.atoms.get(elem, 0)`. When `check_charge=True`, run a
      second pass that sums `e.coefficient * e.species.charge` and
      raises `StoichiometryError(element="charge", ...)` on imbalance.
      One iteration pattern, two conservation laws.
- [ ] [`src/reactions/reaction.py:96-109`](../../src/reactions/reaction.py#L96)
      — replace `entry.species_id` with `entry.species.id` in
      `compute_rates` and the `species_ids` property.
- [ ] [`src/reactions/reaction.py:116-121`](../../src/reactions/reaction.py#L116)
      — repr formatter uses `e.species.id`.
- [ ] [`src/reactions/__init__.py`](../../src/reactions/__init__.py)
      — no API removals; `StoichiometryEntry` stays importable but
      the constructor signature changes.

Sanity check: `Grep "entry\.atoms|entry\.species_id|\.atoms,? MW="
src/reactions/` returns zero matches after this checkpoint.

### 3. Add `kind` and `log_K` to `Reaction`

- [ ] [`src/reactions/reaction.py:31-76`](../../src/reactions/reaction.py#L31)
      — add `kind: str = "kinetic"` and `log_K: Optional[float] = None`
      keyword-only constructor parameters. Make `rate_fn` `Optional`.
      Store on the instance.
- [ ] **Scope extension (2026-05-13 during implementation):** also
      add `dH_J_per_mol: Optional[float] = None` and
      `T_ref_K: float = 298.15` as equilibrium-only constructor
      parameters. Required so that
      `SpeciationEngine.from_reactions(...)` can carry Van 't Hoff
      temperature correction through to the engine's
      `EquilibriumDef.pKas_at_T(T_K)` path — without this, BSM2's 35°C
      operation would drift the sentinels off their 25°C reference
      pKa values (the existing `chem_env_fn` pre-corrects pKa values
      via `ThermodynamicConfig.context_kwargs_at_current_T`, and that
      path goes away in this branch). Phase 3's `ThermoFramework`
      will absorb these fields; in Phase 1 they live on the reaction.
- [ ] In `Reaction.__init__`, validate combinations:
      - `kind="kinetic"` requires `rate_fn` not `None`, forbids
        `log_K`.
      - `kind="equilibrium"` requires `log_K is not None`, forbids
        `rate_fn`.
      - Any other `kind` value raises `ValueError`.
- [ ] [`src/reactions/reaction.py:78-104`](../../src/reactions/reaction.py#L78)
      — `compute_rates(env)` raises a clear
      `RuntimeError("Reaction(kind='equilibrium') has no rate; route
      via SpeciationEngine.from_reactions().")` for equilibrium
      reactions. Kinetic path unchanged.
- [ ] [`src/reactions/protocols.py`](../../src/reactions/protocols.py)
      — no change. `ReactionModel.compute_rates` continues to mean
      "kinetic source terms"; the `partition()` call is upstream.

Sanity check: construct a `Reaction(kind="equilibrium",
stoichiometry=[...], log_K=-14.0)` succeeds; passing both `rate_fn`
and `log_K` raises `ValueError`.

### 4. `ReactionSet.partition()`

- [ ] [`src/reactions/reaction_set.py:23-93`](../../src/reactions/reaction_set.py#L23)
      — add `partition(self) -> Tuple[ReactionSet, ReactionSet]`
      returning `(kinetic_set, equilibrium_set)`. Each subset
      preserves order and inherits `self.label` with a suffix
      (`f"{label}_kinetic"` / `f"{label}_equilibrium"`).
- [ ] `ReactionSet.compute_rates` already calls each constituent's
      `compute_rates`; this will raise on equilibrium reactions, which
      is correct (the integrator should never see them — the builder
      should have routed them out). Add a one-line `__init__` check:
      if any reaction in the set has `kind="equilibrium"`, emit a
      `UserWarning` reminding the caller to `partition()` before
      attaching to a CV. Do not raise — `ReactionSet` is also used as
      the canonical declaration container before partitioning.

Sanity check: `rxn_set.partition()` round-trips:
`len(kinetic) + len(equilibrium) == len(rxn_set)`; ordering
preserved.

### 5. `SpeciationEngine.from_reactions(...)`

- [ ] [`src/speciation/engine.py`](../../src/speciation/engine.py) —
      add a `@classmethod from_reactions(cls, equilibrium_rxns, *,
      activity_model="davies", T_K=298.15, **engine_kwargs)`. Internal
      logic: build an `EquilibriumSet` from the supplied reactions
      (translate each `Reaction(kind="equilibrium")` into an
      `EquilibriumDef` keyed by the protonated species `id`, with
      `pKas=(-log_K,)` and category inferred from the stoichiometry
      sign of `H+`). Construct the engine carrying that
      `EquilibriumSet` and the activity-model parameters.
- [ ] Engine's `solve(**kwargs)` learns to read totals from a
      `phases` argument when supplied. New signature:
      `solve(self, *, phases=None, **kwargs)`. When `phases` is
      provided, derive `acid_totals`, `CT_TIC`, `CT_NH_T` from
      `phases["liquid"].n_mol` divided by `phases["liquid"].V_L`. The
      mapping `species_id → total_key` lives on the engine's
      `EquilibriumSet` (already present as `EquilibriumDef.total_key`).
- [ ] Keep the legacy keyword-arg path (`acid_totals=...`, `CT_TIC=...`)
      working until checkpoint 6 finishes the
      `SpeciationPropertySolver` rewrite. Internal precedence:
      explicit kwargs win over phase-derived values, so existing
      callers keep working until they migrate.
- [ ] Engine's warm-start cache (`_logH_last`, `_I_last`) stays where
      it is. No `LiquidPhase` cache.
- [ ] [`src/speciation/__init__.py`](../../src/speciation/__init__.py)
      — re-export `SpeciationEngine.from_reactions` is automatic via
      the class export.

Sanity check: build a small equilibrium-only `ReactionSet` (water +
acetate), call `SpeciationEngine.from_reactions(eq_rxns)`, then
`engine.solve(phases={"liquid": liq})`, confirm `pH` matches a
direct `EquilibriumSet` build.

### 6. Shrink `chem_env`; rewrite `SpeciationPropertySolver`

- [ ] [`src/core/speciation_solver.py`](../../src/core/speciation_solver.py)
      — rewrite per the docstring forward note (lines 51–66 say
      "When chemistry-unification-1 lands, rewrite this class rather
      than adding to it"). Remove the `acid_pKas` and
      `equilibrium_set` constructor arguments. New signature:
      `SpeciationPropertySolver(engine, *, liquid_phase_key="liquid")`.
- [ ] `solve(self, phases, chem_env, dt_h)` no longer reads
      `acid_totals`, `acid_pKas`, `CT_TIC`, `CT_NH_T`, `CT_P`,
      `acid_dH_J_per_mol`, `Kw`, `pKa1_TIC`, `pKa2_TIC`, `pKa_NH`,
      `equilibrium_set` from `chem_env`. It passes `phases=phases`
      directly to `engine.solve(...)` and forwards only the surviving
      keys: `T_K`, `logH_guess`, `strong_kwargs`. Remove the per-key
      fall-through loop for `dH_*` and pre-corrected pKa overrides
      (those values now live on the `EquilibriumSet` carried by the
      engine, sourced from the equilibrium reactions).
- [ ] [`src/core/property_solver.py:1-22`](../../src/core/property_solver.py#L1)
      — update the module docstring to reflect the new `chem_env`
      contract: `(t_h, T_K, strong-ion params, optional solver
      hints)`. Drop the bullet list mentioning `acid_totals` /
      `acid_pKas` / `CT_NH_T` etc.
- [ ] [`src/core/property_solver.py:84-109`](../../src/core/property_solver.py#L84)
      — `PropertySolver.solve` docstring updated for the new
      `chem_env` contract (drop `acid_totals`, `acid_pKas`,
      `strong_kwargs` references that no longer apply universally —
      only `t_h`, `T_K`, `logH_guess`, and per-solver hints remain
      contractual).
- [ ] [`src/reactions/environment.py:19-50`](../../src/reactions/environment.py#L19)
      — add `has_pH` property. Implementation: `return self.pH is not
      None`. One-liner. Tests in checkpoint 12.

Sanity check: `Grep "acid_totals|acid_pKas|CT_TIC|CT_NH_T|CT_P"
src/core/` should match only [`property_solver.py`'s migration-note
docstring](../../src/core/property_solver.py) (if any) — every other
src/core hit must be gone.

### 7. `ControlVolume` runs species check

- [ ] [`src/core/control_volume.py:75-112`](../../src/core/control_volume.py#L75)
      — at the end of `__init__`, after the interface validation, if
      `reaction_model is not None`, gather all reactions
      (`reaction_model.reactions` for a `ReactionSet`, or
      `[reaction_model]` if it has a `.stoichiometry` attribute,
      otherwise skip) and call `check_species_consistency(reactions,
      soft_conflicts="warn")`. Skip silently for objects that don't
      expose stoichiometry (e.g. `BlackBoxReactionModel`).
- [ ] No call-site changes needed in `advance()` — property solver
      ordering and `_build_reaction_environment` already match the
      stateless snapshot model.

Sanity check: building a CV with a reaction set that redeclares a
species emits a `UserWarning`; constructing with conflicting `atoms`
raises `SpeciesConflictError`.

### 8. Migrate BSM2 + ADM1 builders

- [ ] [`models/vlmodels/adm1/bsm2.py:53-74`](../../models/vlmodels/adm1/bsm2.py#L53)
      — replace the `SPECIES_BSM2` `(atoms, MW, ThOD)` tuple-dict with
      module-level `Species` declarations. Keep ThOD as a side table
      (`THOD_BSM2: Dict[str, float]`) — it isn't part of `Species`
      because it's a model-specific COD convention, not an intrinsic
      property. Bring shared inorganics
      (`CO2`, `NH3`, `H2O`) in via
      `from VLsim.chemistry.common_species import CO2, NH3, H2O`.
- [ ] [`models/vlmodels/adm1/bsm2.py:76-91`](../../models/vlmodels/adm1/bsm2.py#L76)
      — `_atoms`, `_mw` lookup helpers either retire (callers go
      through `Species` directly) or become thin wrappers reading off
      `Species` objects. `_thod` and `_nC`/`_nN` stay (they're
      ThOD-specific or convenience accessors, not Species data).
- [ ] [`models/vlmodels/adm1/bsm2.py:406-412`](../../models/vlmodels/adm1/bsm2.py#L406)
      — `_make_entries(coeff_tuples)` becomes
      `[StoichiometryEntry(species=SPECIES[sp], phase="liquid",
      coefficient=coeff) for sp, coeff in coeff_tuples]`. `SPECIES`
      is the new module-level dict mapping `id` → `Species` object.
- [ ] [`models/vlmodels/adm1/bsm2.py:690-779`](../../models/vlmodels/adm1/bsm2.py#L690)
      — `make_bsm2_chem_env_fn` shrinks to populate only
      `{"strong_kwargs": {"CT_cation": ..., "CT_anion": ...}}` (plus
      `t_h` if the caller threads it). Drop `acid_totals`, `CT_TIC`,
      `CT_NH_T`, `CT_P`, the `tc.context_kwargs_at_current_T()`
      block, and the `H_prev`-based `logH_guess` warm-start
      (warm-start lives on the engine now). The `ThermodynamicConfig`
      requirement check becomes a check that the CV has a
      speciation solver carrying an `EquilibriumSet` derived from
      equilibrium reactions — a softer requirement, just enough to
      give a useful error.
- [ ] [`models/vlmodels/adm1/bsm2.py:782-801`](../../models/vlmodels/adm1/bsm2.py#L782)
      — **delete** `make_bsm2_callback` entirely. Drop the
      docstring's `make_bsm2_callback` reference at line 27.
- [ ] [`models/vlmodels/adm1/bsm2.py:585-687`](../../models/vlmodels/adm1/bsm2.py#L585)
      — `build_bsm2_cv` adopts the partition pattern: after
      `reaction_model(reaction_set)` is wired, it calls
      `kinetic_rxns, equilibrium_rxns = reaction_set.partition()`
      and constructs a `SpeciationEngine.from_reactions(equilibrium_rxns,
      activity_model=activity_model, T_K=T_K)`. The
      `.chemistry(speciation_level=..., use_activity=..., activity_model=...)`
      call on `FermenterBuilder` either drops `speciation_level` (the
      reactions are now self-describing) or stays as a hint to the
      engine; pick one and document. Adding equilibrium reactions to
      `reaction_set` is the responsibility of `build_bsm2_reactions`
      — see next bullet.
- [ ] [`models/vlmodels/adm1/bsm2.py:419-582`](../../models/vlmodels/adm1/bsm2.py#L419)
      — `build_bsm2_reactions(...)` appends equilibrium reactions
      (water + carbonate first dissociation + ammonium + 4 VFAs) at
      the end of the kinetic list. Use `Species` references; pull pKa
      values from `BSM2_KINETICS`-style constants near the top of the
      file (or a new `BSM2_EQUILIBRIA` dict). Sentinel cross-check:
      after migration, `EquilibriumSet.bsm2_default()` (existing
      preset) and the engine built via `from_reactions` should
      produce the same pKa table at identical T_K.
- [ ] [`models/vlmodels/adm1/bsm2_direct.py:36-40`](../../models/vlmodels/adm1/bsm2_direct.py#L36)
      — same import update: `SPECIES_BSM2` is now a `Species` dict,
      not `(atoms, MW, ThOD)` tuples. Verify nothing in
      `bsm2_direct.py` reaches into `_atoms`/`_mw` via the old tuple
      indexing.
- [ ] [`models/vlmodels/adm1/base.py`](../../models/vlmodels/adm1/base.py)
      — same `StoichiometryEntry` constructor migration. Inspect at
      edit time; the file is small enough to handle as one pass.
- [ ] [`models/vlmodels/fermenter/config/factory.py:298-309`](../../models/vlmodels/fermenter/config/factory.py#L298)
      — `FermenterFactory.create_volume` constructs the engine from
      reactions when a `reaction_model` is supplied:
      `engine = SpeciationEngine.from_reactions(rxn_set.partition()[1],
      activity_model=..., T_K=...)`. When no reaction model exists
      (empty fermenter), fall back to a default empty equilibrium
      set so the engine still functions for tests.
- [ ] [`models/vlmodels/fermenter/config/builder.py:354-369`](../../models/vlmodels/fermenter/config/builder.py#L354)
      — `chemistry(...)` builder method drops `acid_pKas` parameter
      (pKas come from equilibrium reactions). Other args
      (`speciation_level`, `use_activity`, `activity_model`) stay.
- [ ] [`models/vlmodels/fermenter/config/factory.py:480-531`](../../models/vlmodels/fermenter/config/factory.py#L480)
      — `run_batch` default `chem_env` shrinks: drop `acid_totals`,
      `acid_pKas`, `CT_P`, `CT_NH_T`. Keep `strong_kwargs` empty.
- [ ] **systems/ scripts** — eight files in the top-level `systems/`
      directory pass legacy `StoichiometryEntry(species_id, phase,
      coeff, atoms, MW)` calls (see grep results in *Cross-cutting*
      below). They each need the constructor migration. Treat them as
      one pass in a single commit; they are not under test, so
      breakage is bounded.

### 9. Verify B2 constraint: no kinetic reaction touches H⁺

The stateless snapshot model assumes kinetic stoichiometries do **not**
contain H⁺ (or other equilibrium-output species like OH⁻, HCO₃⁻,
CO₃²⁻, NH₄⁺). Verify before integration tests run.

- [ ] `Grep '"H\+"' models/` and `Grep '"H\+"' systems/` and inspect
      every match. If any kinetic reaction's stoichiometry mentions
      H⁺, file a separate ticket — that reaction needs reframing in
      totals (per the B2 resolution constraint, this is "deliberate,
      not silent loss"). Do **not** carry the old behaviour forward
      under the stateless model.
- [ ] Same scan for `"OH-"`, `"HCO3-"`, `"CO3--"`, `"NH4+"`. The
      expected outcome is zero matches in stoichiometry positions
      (these may legitimately appear in `species` dicts on
      `PropertyResult` consumers — that's fine).

### 10. Tests — rewrites

The test rewrites land after the production code is migrated so the
suite runs end-to-end at each prior checkpoint with `xfail` markers
where needed. Counts are estimates (refine after edits land).

- [ ] [`tests/standalone/test_reactions.py`](../../tests/standalone/test_reactions.py)
      (44 tests) — every `StoichiometryEntry(species_id, phase, coeff,
      atoms, ...)` construction migrates to
      `StoichiometryEntry(species=Species(id=..., atoms=..., charge=0),
      phase=..., coefficient=...)` or to `StoichiometryEntry(species=A,
      phase=..., coefficient=...)` with `A = Species(...)` defined at
      the top of the helper. Add tests for `kind="equilibrium"` +
      `log_K`, the new validation paths, and `compute_rates` raising
      on equilibrium reactions. Also add tests for `ReactionSet.partition()`.
- [ ] [`tests/standalone/test_speciation.py`](../../tests/standalone/test_speciation.py)
      (~30 tests) — substantially rewritten. The `_solve_simple`
      helper currently wires `acid_totals`/`acid_pKas` through
      `engine.solve(**kwargs)`. After the rewrite, each test either
      keeps using the legacy kwargs (still supported per checkpoint 5,
      so we keep legacy-shape tests as a regression net) **or** is
      ported to the `from_reactions` API (preferred for new tests).
      Decision per test: if it tests the math of pH/CT_TIC, keep
      legacy form; if it tests the *interface*, port. Document the
      split in the file's module docstring.
- [ ] [`tests/standalone/test_property_solvers.py`](../../tests/standalone/test_property_solvers.py)
      (16 tests) — the helper `_make_phases_and_chem_env`'s
      `chem_env` dict drops `acid_totals`, `acid_pKas`, `CT_P`,
      `CT_NH_T`. The `_make_engine_and_solver` helper builds the
      engine via `from_reactions` (small acetate + carbonate
      equilibrium reaction set) so the parity test at lines 90–138
      still has something to compare against. The
      `acid_pKas={"AceticAcid": 4.76}` constructor argument drops
      from `SpeciationPropertySolver(...)` — the engine carries it.
- [ ] [`tests/standalone/test_cv_advance.py`](../../tests/standalone/test_cv_advance.py)
      (47 tests) — every `chem_env = {...}` dict in `cv.advance(...)`
      calls drops the legacy keys. ~20 such constructions
      (`chem_env = {` lines at 132, 145, 164, 463, 472, 490, 501,
      515, 655 plus a few more — confirm at edit time). The `t_h`
      key stays. Where a test explicitly exercises `acid_totals`
      (for property-solver flow), port to the new pattern: pass the
      acid totals via the underlying liquid phase moles instead.
- [ ] [`tests/standalone/test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py)
      (6 tests) — the fixture at lines 79–134:
      drop the `make_bsm2_callback` import; drop the
      `callback = make_bsm2_callback(cv)` call; drop the
      `callback(cv, result, t_h, DT_H, step)` line; the
      `chem_env_fn` call in line 120 still works (its output is
      smaller now). The sentinel values **must remain bit-equivalent**
      — this is the explicit success criterion (per the plan: "pure
      structural change, no numerical change"). If sentinels drift,
      stop and investigate before adjusting them.
- [ ] [`tests/standalone/test_factory.py`](../../tests/standalone/test_factory.py)
      and
      [`tests/standalone/test_builder.py`](../../tests/standalone/test_builder.py)
      — drop `acid_pKas=` from any
      `FermenterBuilder().chemistry(...)` calls. Drop default
      `chem_env` keys from any `run_batch` call sites that pass an
      explicit dict.
- [ ] [`tests/standalone/test_stoichiometry.py`](../../tests/standalone/test_stoichiometry.py)
      — this file tests an unrelated legacy fermenter path; verify it
      doesn't touch `StoichiometryEntry`. (Spot-check: it doesn't.)

### 11. Tests — new

- [ ] New file `tests/standalone/test_species.py` (~12 tests):
      - `Species` field defaults, `eq=True` semantics (two equal
        objects compare equal; modifying one's `atoms` is impossible —
        frozen).
      - `__hash__` works (so `Species` can live in a `set`).
      - `Species(id="X", atoms={"C":1}, charge=-1)` — charge
        accessible.
      - `validate_balance` with charge: balanced
        (`HA → A⁻ + H⁺` with `check_charge=True`) passes; imbalanced
        raises with `element="charge"`.
      - `common_species`: `H_plus.charge == 1`, `OH_minus.charge == -1`,
        `CO3_2minus.charge == -2`, `H2O.atoms == {"H":2, "O":1}`,
        identity preserved across imports.
- [ ] New file `tests/standalone/test_species_check.py` (~8 tests):
      - No conflicts → no warning, no exception.
      - Hard conflict (same id, different atoms) → `SpeciesConflictError`.
      - Hard conflict (same id, different charge) → `SpeciesConflictError`.
      - Soft conflict (distinct objects, equal data) →
        `UserWarning` when `soft_conflicts="warn"`.
      - Soft conflict + `soft_conflicts="raise"` → `SpeciesConflictError`.
      - Soft conflict + `soft_conflicts="ignore"` → silent.
      - Composed `ReactionSet` of imported `common_species` →
        no warnings.
- [ ] New tests in `test_reactions.py` for `ReactionSet.partition()`,
      `Reaction(kind="equilibrium")`, and the validation matrix.
      Estimate: +6 tests in `test_reactions.py`.
- [ ] New test in `test_property_solvers.py`: building
      `SpeciationEngine.from_reactions(eq_rxns)` and confirming a
      solve produces the same pH as a hand-built `EquilibriumSet`
      with identical pKas.
- [ ] New test in `test_property_solvers.py` or
      `test_cv_advance.py`: `env.has_pH` returns `True` after a
      speciation solve, `False` for an environment built with
      `pH=None`.

### 12. Update class diagrams + speciation_solver forward note

- [ ] [`docs/class_diagrams.md:305-323`](../class_diagrams.md#L305) —
      the "Forward note (transient framing)" block. Rewrite to
      describe the **post-Phase-1** state: equilibrium reactions live
      on the same `Reaction` class as kinetic ones (distinguished by
      `kind`), `chem_env` is a near-empty dict
      (`(t_h, T_K, strong-ion params, optional solver hints)`), and
      `SpeciationPropertySolver` is a thin, deliberately retained
      engine adapter (NOT "absorbed" — drop that word). The note
      now describes what *remains* deferred (Phase 2 alphas, Phase 3
      `ChemistryDatabase`, Phase 4 `AccuracyMonitor`) — not what
      lands in this branch.
- [ ] [`docs/class_diagrams.md:410`](../class_diagrams.md#L410) —
      `StoichiometryEntry` class block: replace `+species_id: str` /
      `+atoms: dict` / `+MW: float` with `+species: Species` /
      `+phase: str` / `+coefficient: float`. Add a new `Species`
      class block with `+id: str` / `+atoms: dict` / `+charge: int`
      / `+MW: float`. Add the relationship arrow:
      `StoichiometryEntry "1" --> "1" Species : species`.
- [ ] [`docs/class_diagrams.md:365-373`](../class_diagrams.md#L365) —
      `Reaction` class block: add `+kind: str` and
      `+log_K: Optional[float]` lines. `+rate_fn: callable` becomes
      `+rate_fn: Optional[callable]`.
- [ ] [`docs/class_diagrams.md:374-379`](../class_diagrams.md#L374) —
      `ReactionSet` class block: add `+partition() (kinetic, equilibrium)`.
- [ ] [`docs/class_diagrams.md:342-348`](../class_diagrams.md#L342) —
      `SpeciationPropertySolver` block: drop `acid_pKas` and any
      `equilibrium_set` references. Add the `+from_reactions()$`
      static-method line on the `SpeciationEngine` block (which
      should also be present somewhere in the layer 3 mermaid block).
- [ ] [`src/core/speciation_solver.py:51-66`](../../src/core/speciation_solver.py#L51)
      — rewrite the "Forward note" docstring. Status: **stable
      adapter** as of `chemistry-unification-1` (NOT "absorbed" — the
      class is deliberately retained, see the resolved decision
      block above). The new docstring should describe the
      `PropertySolver`-protocol-adapter role and point to
      `SpeciationEngine` for direct standalone use.
- [ ] [`docs/architecture.md`](../architecture.md) — search for
      `chem_env` / `acid_totals` / `acid_pKas` mentions and update to
      reflect the new minimal contract. Phase 3 will overhaul the
      chemistry section more substantially; this pass is a
      consistency tweak, not a rewrite.

### 13. Final test sweep + ship

- [ ] Run the full standalone test suite:
      `python -m pytest tests/standalone -q`. Expect ~712 tests
      (current baseline) ± the net of additions in checkpoint 11.
      Estimate: +20 tests (test_species + test_species_check +
      partition tests + has_pH test) ≈ ~732 passing.
- [ ] **Critical:** confirm
      [`test_bsm2_reference.py`](../../tests/standalone/test_bsm2_reference.py)
      sentinels pass bit-for-bit. The `RTOL_SENTINEL = 1e-9` budget
      is the success criterion. Any drift means the migration
      changed numerics, which it should not.
- [ ] Update [README.md](README.md) priority list: mark
      `chemistry-unification-1` as shipped; flag
      `chemistry-unification-2` as next.
- [ ] Move
      [CHEMISTRY_UNIFICATION_PLAN.md](CHEMISTRY_UNIFICATION_PLAN.md)
      to [`../shipped/`](../shipped/) only after
      Phase 4 ships (per the plan's "Sequencing of doc updates").
      For this branch, only this checklist moves.
- [ ] Move this checklist to
      [`../shipped/CHEMISTRY_UNIFICATION_1_CHECKLIST.md`](../shipped/)
      with a "Shipped" status banner at the top.
- [ ] Ship via the convention in
      [README.md](README.md):
      `git checkout main && git merge --no-ff chemistry-unification-1
      -m "Merge chemistry-unification-1: unify reaction declaration,
      Species type, stateless snapshot speciation"`.
- [ ] Tag: `git tag chemistry-unification-1-shipped <commit-hash>`.
- [ ] Push: `git push && git push --tags`.
- [ ] Delete branch:
      `git branch -d chemistry-unification-1` and
      `git push origin --delete chemistry-unification-1`.

## Cross-cutting reference: where `StoichiometryEntry` is constructed

For checkpoint 8's systems/ migration. From `Grep "StoichiometryEntry"`:

- `models/vlmodels/adm1/bsm2.py:411` — central `_make_entries` helper.
- `models/vlmodels/adm1/base.py` — direct constructions.
- `systems/twelve_rxn_AD.py`, `systems/twelve_rxn_AD_stage2_3.py`,
  `systems/three_stage_AD.py`, `systems/nine_rxn_AD.py`,
  `systems/five_stage_AD.py`, `systems/eight_rxn_AD.py`,
  `systems/cstr_ad_stage7_17.py`,
  `systems/anaerobic_digestion.py` — all eight use the legacy
  signature. Most likely a search-and-replace: `StoichiometryEntry(sp,
  ph, coeff, atoms_dict, MW)` →
  `StoichiometryEntry(SPECIES[sp], ph, coeff)` with a local `SPECIES`
  dict near the top.
- `tests/standalone/test_reactions.py`, `tests/standalone/test_cv_advance.py`,
  `tests/standalone/test_membrane.py`,
  `tests/standalone/test_multi_cv.py` — covered in checkpoint 10.

## Final test count expectation

712 → ~732 (estimate: +12 in `test_species`, +8 in
`test_species_check`, +6 in `test_reactions`, +1 in
`test_property_solvers`, −0 deletions). Final number confirmed
empirically at checkpoint 13.
