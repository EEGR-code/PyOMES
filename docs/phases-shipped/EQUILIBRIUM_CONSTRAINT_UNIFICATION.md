# EQUILIBRIUM_CONSTRAINT_UNIFICATION — Phase 1 of 2

> **Status:** Shipped 2026-07-02. Tag `equilibrium-constraint-unification-shipped`
>   (after merge). CP1–CP5 all landed; two unplanned findings folded into CP3
>   (see its section below) — a `ReactionSystem` bucketing blocker and an
>   adjacent precipitation-forwarding gap — plus a J/(mol·K) gas-constant
>   consolidation surfaced while testing CP1 (`VLsim.units.R_J_PER_MOL_K`).
> **Design authority:** [MASS_EXCHANGE_ARCHITECTURE.md](../design/MASS_EXCHANGE_ARCHITECTURE.md)
>   §14 (constraint-family taxonomy, wrapping litmus test, precipitation
>   resolution); [CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md](../design/CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md)
>   §3 (`ChemicalEquilibriumSystem`, two-family `reactions` typing).
> **Companion phase:** [LAYER1_GAP_CLOSURE.md](LAYER1_GAP_CLOSURE.md) (Phase 2
>   of 2) — depends on this phase shipping first. Carries the
>   `SpeciationEngine` → `ChemicalEquilibriumEngine` rename (§18 of the
>   architecture doc) at its close, not here — this phase does not change
>   `solve()`-time behavior, so a rename here would have no functional payoff.

---

## Background and motivation

Design discussion on 2026-07-01 (see chat transcript / `MASS_EXCHANGE_ARCHITECTURE.md`
§14 for the full derivation) surfaced a concrete, currently-shipped
inconsistency in how gas-liquid and solid-liquid equilibria are declared
relative to acid-base:

- `SpeciationEngine.from_reactions()` silently drops any `EquilibriumReaction`
  where `is_cross_phase` is true (`src/speciation/engine.py`), with a comment
  that these are "consumed by the gas-liquid link" — but they aren't actually
  consumed by anything; they're just discarded.
- `KineticGasLiquidLink` takes its actual Henry parameters through a
  completely separate `partition_models={"CO2": HenryPartition(H_ref=...,
  dlnH=...)}` dict (`src/core/gas_liquid_link.py`), keyed by species-id string,
  with no reference back to any declared reaction. A user can declare
  `CO2(gas) ⇌ CO2(aq)` as a reaction (discarded) and separately configure
  `HenryPartition(H_ref=...)` with parameters that must be manually kept
  consistent — nothing validates they agree.
- `NRSpeciationEngine.from_reactions()` requires precipitation reactions
  pre-sorted into a separate `precipitation_reactions=` kwarg
  (`src/speciation/nr_engine.py`), rather than being auto-detected from the
  `StoichiometryEntry(phase="solid")` tag that's already present on the
  reaction's own stoichiometry.

Investigation also established that Henry's law and Ksp are not fundamentally
different from an acid-base log-K reaction — they're the same mass-action
structure with domain-conventional parameters (`H_ref`/`dlnH` map directly onto
`log_K`/van't Hoff slope; the shipped `NR_PRECIPITATION_SPECIATION` already
sets `log_K = Ksp`). This phase unifies the *declaration* surface around that
insight. It deliberately does **not** change solve-time behavior: gas-liquid
transfer continues via the existing SNIA sequential correction
(`EquilibriumTransferModel`/`KineticGasLiquidLink`), and precipitation
continues via the existing nested active-set loop in `NRSpeciationEngine`.
Folding these into one simultaneous Newton system is
[LAYER1_GAP_CLOSURE.md](LAYER1_GAP_CLOSURE.md)'s job, not this phase's.

---

## Key design decisions

**Two constraint families.** Mass-action equilibria (acid-base, Henry, Ksp,
Raoult) are one family; isotherms (Langmuir, Freundlich) are a genuinely
separate one — not in scope for this phase, no concrete use case exists yet.

**Shared protocol, not classmethod constructors (2026-07-01 refinement).**
Rather than one concrete type (`EquilibriumReaction`) with several
`from_henry_ref(...)`-style alternate constructors that all collapse into the
same internal `log_K`/`dH_J_per_mol` fields, the mass-action family is a set
of sibling types, each with `__init__` keywords matching its own standard
domain terminology (`H_ref`/`dlnH` for Henry, `Ksp` for solubility products),
all structurally satisfying one shared Protocol:

```python
class EquilibriumConstraint(Protocol):
    stoichiometry: Tuple[StoichiometryEntry, ...]
    def log_K(self, T_K: float) -> float: ...
    def dH_J_per_mol(self) -> Optional[float]: ...
```

`log_K(T_K)` is the *reference* mass-action constant only (temperature-corrected
via van't Hoff) — composition-dependent corrections (`γ_i`, `φ_i`) are not part
of this protocol; they remain Phase 2's job, evaluated directly from
`ThermoFramework` primitives at NR trial points (see the wrapping litmus test
below). This matters because it means a Henry-shaped type never converts away
from its native fields — `henry_obj.H_ref` and `henry_obj.log_K(298.15)` are
both always available, the latter computed fresh from the former.

**Naming: the `*Equilibrium` family, reusing existing `PartitionModel` types
where they already exist.** `EquilibriumReaction` needs no rename (already
fits the shared root, zero blast radius). `HenryPartition` (shipped) is
renamed to `HenryEquilibrium`, with `HenryPartition` kept as a backward-compatible
alias for one phase (matching the `DaviesActivityModel = DaviesLiquidModel`
pattern in `THERMODYNAMIC_MODEL_ARCHITECTURE.md`). `KspPartition`/`RaoultPartition`
(design-stage only, not yet implemented) are introduced directly as
`KspEquilibrium`/`RaoultEquilibrium` — no rename needed since nothing depends
on the old names yet.

**One object serves both roles — this is what fixes CP3's double-declaration
bug.** `HenryEquilibrium`, `RaoultEquilibrium`, and single-ion `KspEquilibrium`
each satisfy **both** `PartitionModel` (their existing role, feeding
`TransferModel`/kinetic use unchanged) **and** the new `EquilibriumConstraint`
(feeding the reaction list/tableau). The same constructed instance is used in
both places — no separate, independently-parameterized declaration, no
derive-or-validate step needed.

**Wrapping litmus test** (governs which types may dual-satisfy both
protocols): a type may satisfy `EquilibriumConstraint` alongside
`PartitionModel` only if its methods are pure evaluate-at-a-point functions,
never ones that internally solve/converge. `HenryPartition`/`RaoultPartition`/
single-ion `KspPartition` all have `partition_ratio() → float` today — the
codebase's own signal of linearity/closed-form — so this is safe.
`MultispeciesVLEPartition` (internally solves the coupled VLE) does **not**
get this treatment; it stays `PartitionModel`-only, and Phase 2's folding for
non-ideal gas mixtures reads `GasEOS`/`LiquidPhaseModel` primitives directly
instead (per §14.2 of `MASS_EXCHANGE_ARCHITECTURE.md`).

**No rename of the engine classes in this phase.** `SpeciationEngine`/
`NRSpeciationEngine` keep their current names throughout — unrelated to the
`HenryPartition` → `HenryEquilibrium` rename above, which is a declaration-type
rename, not an engine rename. The engine rename accompanies the point where
`solve()` actually resolves gas-liquid/solid-liquid simultaneously (Phase 2's
close).

**Precipitation stays exactly as it works today.** This phase does not touch
`NRSpeciationEngine`'s active-set loop's numerics — only how the underlying
reaction gets declared/classified, formalizing the existing `log_K = Ksp`
convention with a named constructor rather than requiring hand-built
`EquilibriumReaction` objects.

---

## Scope

**In scope:**

- `EquilibriumConstraint` Protocol (`stoichiometry`, `log_K(T_K)`,
  `dH_J_per_mol()`) — the shared structural contract for the mass-action
  family. `EquilibriumReaction` satisfies it as-is (no change needed — it
  already stores `log_K`/`dH_J_per_mol` natively).
- `HenryEquilibrium` — renamed from `HenryPartition` (shipped); gains
  `log_K(T_K)`/`dH_J_per_mol()` derived from its existing `H_ref`/`dlnH`
  fields, plus explicit `gas_species`/`liquid_species` stoichiometry fields
  (today it's keyed externally by a dict key in `KineticGasLiquidLink`, not an
  internal field — this phase makes that explicit). Keeps its existing
  `PartitionModel`-satisfying methods (`equilibrium_a_moles`, `partition_ratio`)
  unchanged. `HenryPartition = HenryEquilibrium` kept as a backward-compatible
  alias for one phase.
- `KspEquilibrium`, `RaoultEquilibrium` — new types (design-stage
  `KspPartition`/`RaoultPartition` were never implemented under those names,
  so no rename/alias needed), each satisfying both `PartitionModel` and
  `EquilibriumConstraint` the same way.
- A single auto-classification step replacing the two current patchwork
  mechanisms: given a declared item's type (`EquilibriumReaction`,
  `HenryEquilibrium`, `KspEquilibrium`, `RaoultEquilibrium` — anything
  satisfying `EquilibriumConstraint`) and its `StoichiometryEntry` phase tags,
  classify it as single-phase acid-base / gas-liquid / solid-liquid
  automatically. Both `SpeciationEngine.from_reactions()`'s silent-filter and
  `NRSpeciationEngine.from_reactions()`'s `precipitation_reactions=` kwarg
  requirement are replaced by this one path — callers pass one flat list of
  any conforming type.
- Fixing the `HenryPartition`/reaction double-declaration: since
  `HenryEquilibrium` satisfies both `PartitionModel` and
  `EquilibriumConstraint`, the *same constructed instance* is passed both to
  `KineticGasLiquidLink`'s `partition_models=` dict and into the reaction list
  feeding `NRSpeciationEngine`/`build_tableau()`. No separate declaration, no
  derive-or-validate step — the bug is fixed by construction.
- A per-entry activity-dispatch helper: given a `StoichiometryEntry` and the
  CV's phases/`ThermoFramework`, return the correct activity evaluation
  (`LiquidPhaseModel.gamma_all()` for liquid, `GasEOS.partial_pressures_atm()`
  for gas, `1.0` for solid). Built and unit-tested in isolation here; Phase 2
  wires it into `build_tableau()`'s Newton residual assembly.
- Full regression safety: every existing SNIA-path behavior (kinetic
  gas-liquid transfer, existing precipitation active-set) unchanged; full test
  suite green throughout, including via the `HenryPartition` alias.

**Out of scope (tracked elsewhere):**

- Actually folding gas-liquid/Ksp/Raoult rows into `g(y,z) = 0` —
  [LAYER1_GAP_CLOSURE.md](LAYER1_GAP_CLOSURE.md), Phase 2.
- `SpeciationEngine` → `ChemicalEquilibriumEngine` rename — Phase 2's closing
  checkpoint.
- Isotherm family (Langmuir/Freundlich) protocol — no concrete use case;
  deferred until one exists.
- Multi-ion Ksp / multi-component basis generalization —
  [MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md](MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md),
  separate standalone pathway.
- `GasEOS`/`LiquidPhaseModel` Jacobian primitives (§8.4, §8.5 of
  `THERMODYNAMIC_MODEL_ARCHITECTURE.md`) — needed only once Phase 2 actually
  solves these rows; this phase declares and dispatches, it doesn't solve.
- Surface complexation / site conservation (§12 Q9 of
  `MASS_EXCHANGE_ARCHITECTURE.md`) — deferred, no change here.

---

## Checkpoints

### CP1 — `EquilibriumConstraint` protocol + `*Equilibrium` sibling types

- [x] Add `EquilibriumConstraint` Protocol (`stoichiometry`, `log_K(T_K)`,
      `dH_J_per_mol()`) to `src/reactions/equilibrium.py` (or a new shared
      module if that creates a circular-import issue with `src/chemistry/`).
      `@runtime_checkable`, matching the codebase's existing protocol style.
- [x] Confirm `EquilibriumReaction` satisfies it with no changes (already
      stores `log_K`/`dH_J_per_mol` natively) — add a structural-typing test
      (`isinstance(rxn, EquilibriumConstraint)`) to lock this in.
- [x] Rename `HenryPartition` → `HenryEquilibrium` in `src/chemistry/partition.py`.
      Add `gas_species`/`liquid_species` constructor fields (today only
      implicit via the external `partition_models={"CO2": ...}` dict key —
      make explicit). Add `log_K(T_K)`/`dH_J_per_mol()` methods derived from
      existing `H_ref`/`dlnH` fields per the §14.1 mapping. Keep
      `equilibrium_a_moles`/`partition_ratio` unchanged. Add
      `HenryPartition = HenryEquilibrium` backward-compatible alias with a
      `DeprecationWarning` on use of the old name (matching the
      `DaviesActivityModel = DaviesLiquidModel` pattern).
- [x] Add `KspEquilibrium` (new — `Ksp`, `dH_J_per_mol`, `stoichiometry`
      constructor fields; solid-phase `StoichiometryEntry` with activity fixed
      at 1 by convention) and `RaoultEquilibrium` (new — `P_sat_ref_atm`,
      `dlnPsat` constructor fields, matching `RaoultPartition`'s design).
      Both satisfy `PartitionModel` (single-species, linear —
      `partition_ratio() → float`) and `EquilibriumConstraint`.
- [x] Unit tests: each type's `log_K(T_K)` matches a manual conversion from its
      native domain parameters; `isinstance(x, EquilibriumConstraint)` and
      `isinstance(x, PartitionModel)` both true for all three `*Equilibrium`
      types; `HenryPartition` alias emits `DeprecationWarning` and constructs
      an identical object to `HenryEquilibrium`.

**Two resolved discrepancies against the checklist text above** (settled
before implementation, not re-litigated): (1) `EquilibriumConstraint` shipped
attribute-based (`stoichiometry`/`log_K`/`dH_J_per_mol`/`T_ref_K` as plain
attributes) with a separate `vant_hoff_log_K()` free function doing the
temperature correction, not `log_K(T_K)`/`dH_J_per_mol()` methods — this
means `EquilibriumReaction` needed zero changes, matching its existing plain
attributes. (2) `RaoultEquilibrium` kept `RaoultPartition`'s real shipped
field names (`P_sat_ref`, `dH_vap`, `T_ref`, `C_water_mol_L`), not the
placeholder `P_sat_ref_atm`/`dlnPsat` above — `RaoultPartition` had already
shipped via `THERMODYNAMIC_MODEL_ARCHITECTURE.md` by the time this doc was
written, so it was renamed like `HenryPartition` rather than introduced fresh.

**One unplanned fix surfaced while writing this checkpoint's tests:** the
J/(mol·K) gas constant was defined with two different precisions across four
files (rounded `8.31446` in `thermo_params.py`/`framework.py` vs. CODATA
`8.31446261815324` in `equilibrium.py`/`nr_tableau.py`), so
`HenryEquilibrium.dH_J_per_mol` didn't round-trip exactly through
`vant_hoff_log_K`. Consolidated into one authority,
`VLsim.units.R_J_PER_MOL_K`; re-baselined 4 BSM2 sentinel values that shifted
at the ~1e-9 relative level as a result (`test_bsm2_reference.py`, dated
comment explains why).

### CP2 — Auto-classification / single ingestion path

- [x] Add a classifier (e.g. `classify_equilibrium_constraint(item) ->
      Literal["acid_base", "gas_liquid", "solid_liquid"]`) that accepts
      anything satisfying `EquilibriumConstraint` and inspects
      `StoichiometryEntry.phase` tags — not an `isinstance(rxn,
      EquilibriumReaction)` check, since `HenryEquilibrium`/`KspEquilibrium`/
      `RaoultEquilibrium` are siblings, not subclasses, of `EquilibriumReaction`.
- [x] `SpeciationEngine.from_reactions()`: stop silently dropping cross-phase
      items; route them through the classifier. Gas-liquid/solid-liquid items
      remain outside this engine's tableau for now (this phase doesn't fold
      them) — but they're now visible, not discarded, and accessible via a
      `cross_phase_constraints` attribute for downstream consumers
      (`KineticGasLiquidLink`, CP3).
- [x] `NRSpeciationEngine.from_reactions()`: accept one flat list of any
      `EquilibriumConstraint`-conforming item; auto-classify solid-liquid
      entries instead of requiring the separate `precipitation_reactions=`
      kwarg. Keep `precipitation_reactions=` as a deprecated alias for one
      phase (raises `DeprecationWarning`, still works).
- [x] Update `build_tableau()`'s reaction filter (`src/speciation/nr_tableau.py:200`,
      currently `isinstance(rxn, EquilibriumReaction)` + `is_cross_phase`) to
      use the shared classifier against `EquilibriumConstraint` instead — still
      filters gas-liquid/solid-liquid out of the graph in this phase (Phase 2
      changes that), but via the shared classification path so it recognizes
      `HenryEquilibrium`/`KspEquilibrium`/`RaoultEquilibrium` too, not just
      `EquilibriumReaction`.
- [x] Tests: one flat list mixing `EquilibriumReaction` + `HenryEquilibrium` +
      `KspEquilibrium` classifies each correctly; existing
      `precipitation_reactions=` call sites still work (deprecated-alias path)
      and emit the warning.

### CP3 — Fix the double-declaration bug (one object, two roles)

- [x] Update `KineticGasLiquidLink`'s `partition_models=` dict to accept
      `HenryEquilibrium` instances directly (already satisfies `PartitionModel`
      — no change needed to the dict's contract, just to what's typically
      passed in).
- [x] Update construction guidance/examples: a single `HenryEquilibrium(H_ref=...,
      dlnH=..., gas_species="CO2", liquid_species="CO2")` instance is
      constructed once, then passed both to `KineticGasLiquidLink`'s
      `partition_models=` dict *and* into the flat reaction list feeding
      `NRSpeciationEngine.from_reactions()`/`SpeciationEngine.from_reactions()`
      — CP2's classifier recognizes it as a gas-liquid constraint from the
      same object, so there is nothing left to keep in sync.
- [x] Update `KineticGasLiquidLink`'s module docstring/example
      (`src/core/gas_liquid_link.py`) to show this single-object pattern
      instead of the current independent-`HenryPartition`-dict example.
- [x] Tests: constructing a CV from one `HenryEquilibrium` instance used in
      both places produces a fully consistent configuration; add a regression
      test that would have caught the old bug (two independently-parameterized
      declarations for the same species silently disagreeing) — confirm the
      new pattern makes that class of bug structurally impossible, not just
      caught at validation time.

**Two unplanned findings surfaced and folded in during implementation**
(discussed with the user before proceeding):

1. **Blocking:** `ReactionSystem` (`src/reactions/reaction_system.py`) — the
   actual single reaction-attach point on a CV, not just the bare
   `from_reactions()` factories CP1/CP2 touched — hard-`isinstance`-checked
   `EquilibriumReaction` and raised `TypeError` on any `EquilibriumConstraint`
   sibling. Without a fix, the "one object, two roles" pattern this
   checkpoint targets couldn't actually work through real CV construction.
   Fixed by switching its bucketing to `classify_equilibrium_constraint()`;
   also fixed `show_reactions()`/`show_balance()`, which assumed
   `EquilibriumReaction`-only attributes, and added a `label` field to
   `HenryEquilibrium`/`RaoultEquilibrium` (missing since CP1, would have
   crashed `show_reactions()`).
2. **Adjacent, pre-existing gap** (not part of the double-declaration bug,
   found while reworking the bucketing above): `ReactionSystem.engine`'s
   `newton_raphson` branch never included `precipitation_equilibria` in the
   list passed to `NRSpeciationEngine.from_reactions()` — a
   `ReactionSystem`'s precipitation reactions never reached the active-set
   solver via that path. CP2's flat-list auto-classification made this a
   one-line fix (concatenate the third bucket); folded in since it's now
   trivial, rather than left as a silent latent bug.

`AD_BASIC` (`src/chemistry/databases/anaerobic_digestion.py`) was also
migrated to the one-object pattern: CO₂'s gas⇌liquid partition was
previously declared twice (a `HenryPartition` in `partition_models`, a
separate log_K-less placeholder `EquilibriumReaction` in the reaction list)
— now one shared `HenryEquilibrium` instance. Discovered along the way:
`KineticGasLiquidLink` is already deprecated for direct construction in
favor of `ControlVolume(transfer_models=...)` (shipped separately, before
this phase) — its docstring now notes this while still documenting the
`partition_models=` contract that `transfer_models=` builds internally.

### CP4 — Per-entry activity-dispatch helper (Phase 2 groundwork)

- [x] Add a helper (e.g. `activity_for_entry(entry: StoichiometryEntry,
      phases: Dict[str, Phase], thermo: ThermoFramework, T_K: float) -> float`)
      that dispatches: `entry.phase == "liquid"` → `thermo.liquid_activity.gamma_all(...)`
      for that species; `entry.phase == "gas"` → `thermo.gas_eos.partial_pressures_atm(...)`
      for that species; `entry.phase == "solid"` → `1.0` (pure-solid convention).
- [x] Built and tested standalone here — **not** wired into `build_tableau()`'s
      Newton residual yet (Phase 2 CP1/CP2 does that). This checkpoint only
      needs to prove the dispatch is correct in isolation against known
      `LiquidPhaseModel`/`GasEOS` outputs.
- [x] Tests: dispatch returns matching values to calling `gamma_all`/
      `partial_pressures_atm` directly for each phase type; solid case returns
      `1.0` unconditionally.

### CP5 — Regression suite + doc closeout

- [x] Full test suite green — no behavior change to any existing SNIA
      gas-liquid transfer or precipitation active-set path.
- [x] Update `MASS_EXCHANGE_ARCHITECTURE.md` §14 with a "shipped" note pointing
      here; no change to §10.4's status table yet (that's Phase 2's job, since
      solve-time behavior hasn't changed).
- [x] Move this file to `docs/phases-shipped/`, tag
      `equilibrium-constraint-unification-shipped`, merge `--no-ff` into
      `main`, delete the feature branch.

---

## Files changed

| File | Nature of change |
|---|---|
| `src/reactions/equilibrium.py` | CP1: `EquilibriumConstraint` Protocol; confirm `EquilibriumReaction` conformance |
| `src/chemistry/partition.py` | CP1: `HenryPartition` → `HenryEquilibrium` (+ alias); new `KspEquilibrium`, `RaoultEquilibrium` |
| `src/speciation/engine.py` | CP2: stop silently dropping cross-phase items; expose them |
| `src/speciation/nr_engine.py` | CP2: flat list + auto-classification; deprecate `precipitation_reactions=` |
| `src/speciation/nr_tableau.py` | CP2: shared classifier in `build_tableau()`'s filter |
| `src/core/gas_liquid_link.py` | CP3: single-object construction pattern; docstring/example update |
| `src/speciation/` (new module, e.g. `activity_dispatch.py`) | CP4: per-entry activity-dispatch helper |
| `tests/` | New tests per checkpoint; deprecation-path coverage |
| `docs/design/MASS_EXCHANGE_ARCHITECTURE.md` | CP5: shipped note |

---

## Trigger conditions

Start whenever you're ready — this phase requires no external trigger; it's
the agreed first half of Layer 1 gap closure, and delivers value (fixing the
declaration double-bug, formalizing Ksp/Henry constructors) independent of
whether Phase 2 is picked up immediately after.

## How to start

Follow the standard convention in [README.md](README.md): branch
`equilibrium-constraint-unification` off `main`, work through CP1→CP5 in
order, one commit per checkpoint, full suite run at each. Tag and merge on
completion.
