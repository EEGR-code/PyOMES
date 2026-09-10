# Notebook Demos

Record of Jupyter notebook demos created, plus any codebase changes made
during their development.  Update this file whenever a notebook is added,
revised, or when a demo-driven API change lands.

---

## Session — 2026-07-08

### Notebooks created: `demos/features/SolverProtocols/`

Two-notebook family, from the `step-solver-interface-refinement` phase's
checkpoint 8, following the `ChemicalEquilibriumProtocol/` folder's
`0_README.ipynb` + numbered-notebook shape. Both generated via
`demos/features/SolverProtocols/_generate_notebooks.py`
(`nb()`/`md()`/`code()` helpers, matching
`ChemicalEquilibriumProtocol/_generate_notebooks.py`); executed via
`jupyter nbconvert --execute` under the `biosteam` kernel, 0 errors.

#### `0_README.ipynb`

Architecture overview (pure markdown, no code) — the two orthogonal solver
axes (`StepSolver`/Axis 1, `SystemSolver`/Axis 2), a comparison table for
each axis's shipped rungs, the clamping (`clamp_fn`/`floor_nonnegative`)
model, and a decision guide. Condensed from `docs/solvers.md`. Added on a
second request the same session, after the walkthrough notebook (below)
already existed as the folder's only file.

#### `01_writing_a_custom_solver.ipynb`

Writing a custom `StepSolver` (Axis 1) and `SystemSolver` (Axis 2).
Originally shipped as a `.py` script (`demos/model_api/custom_solver_demo.py`,
matching the sibling `solver_comparison.py` convention); converted to a
notebook at user request since the checklist item was always titled "demo
notebook," then renamed from `step_and_system_solver_demo.ipynb` to its
current name to fit the two-notebook numbering once `0_README.ipynb` was
added. The `.py` script was deleted (no duplicate demo); cross-references
in `docs/solvers.md` and `demos/README.md` updated to point at the
notebooks.

Two custom solvers, each inverting a real documented ordering decision rather
than an arbitrary toy reordering:
- `SpeciationAfterFeedStepSolver` — reverses `SequentialAdvanceSolver`'s
  pre-feed-speciation contract (see `docs/dev/implementation/shipped/ORDERING.md`).
- `AsymmetricThreeStageSystemSolver` — asymmetric 3-stage link/CV-advance
  interleaving, composing the same `sim._apply_links`/`cv._advance_unchecked`/
  `sim._invoke_controllers` primitives `StrangSplittingSystemSolver` uses.

---

## Session — 2026-06-19

### Notebooks created (updated_layout series)

All three `updated_layout` notebooks follow the BioSTEAM pattern:
pre-code design basis → condensed implementation → explicit validation.
Reading order: `0_README.ipynb` → Example1 → Example2 → Example3.

#### `demos/model_api/D2Cworkshop/updated_layout/0_README.ipynb`

All-markdown overview notebook with reading-order table and description of the
common 5-section structure shared by all notebooks in the folder.

#### `demos/model_api/D2Cworkshop/updated_layout/Example1_mtp_well.ipynb`

| Section | Content |
|---|---|
| 1 | Background: closed headspace as a finite O₂ reservoir; EquilibriumTransferModel rationale (1/kLa = 36 s ≪ 1/μmax = 2 h) |
| 2 | Design basis: elemental balance → ν_O₂ = 1.015; O₂ budget (supply vs stoichiometric demand); safety margin ≈ 2.5×; predictions incl. X_pred = 0.18 g/L |
| 3 | Condensed implementation |
| 4 | Validation: limiting resource check, endpoint biomass vs prediction, observed yield vs Y, plots with design-basis reference lines annotated |

#### `demos/model_api/D2Cworkshop/updated_layout/Example2_batch_fermenter.ipynb`

| Section | Content |
|---|---|
| 1 | Background: MTP well → sparged fermenter; three new design questions (kLa, vvm, pH control) |
| 2 | Design basis: OTR/OUR balance → kLa = 150/h gives 2× safety margin at 50% DO; sparging supply ~50× OUR (low O₂ utilisation); Henderson-Hasselbalch at pH 5.0 / pKa 4.756 → 36% undissociated; four predictions |
| 3 | Condensed implementation; ν_O₂ extracted from stoichiometry for use in validation |
| 4 | Validation: OTR and OUR time series reconstructed from simulation; min DO check; yield check; 2×2 plot (OTR/OUR with shaded surplus, DO%, substrate/biomass with predicted endpoint, pH) |
| 5 | Discussion: kLa sensitivity (halving to 75/h as one-parameter test); pH control feasibility |

#### `demos/model_api/D2Cworkshop/updated_layout/Example3_CSTR.ipynb`

| Section | Content |
|---|---|
| 1 | Background: continuous vs batch comparison table; what LiquidFeed + LiquidDrain add to the model |
| 2 | Design basis: washout curve S*(D) = Ks·D/(μmax−D) and X*(D) = Y(S₀−S*) computed and plotted (no PyOMES imports); D=0.2/h chosen (40% μmax, 99.93% conversion, X*=1.80 g/L); S_feed=5 g/L justified; convergence τ ≈ 1/D → 2 HRTs ≈ 86% |
| 3 | Condensed implementation |
| 4 | Validation: SS accuracy check (S*, X*, conversion); normalised convergence (C−C*)/(C₀−C*) vs exp(−Dt) overlay; 2×2 plot (washout curve X*(D) with simulation endpoint star, normalised convergence, absolute concentrations with SS lines, pH) |
| 5 | Discussion: washout curve as design tool; higher/lower D trade-offs; practical convergence rule of thumb |

---

### Notebooks created

#### `demos/model_api/D2Cworkshop/basic_layout/Example3_CSTR.ipynb`

Extends `Example2_batch_fermenter` by adding `LiquidFeed` and `LiquidDrain`
boundaries to convert the sparged batch vessel into a CSTR at a specified
dilution rate.

Structure (16 sections):

| Section | Content |
|---|---|
| 1 | Path setup + imports (adds `LiquidFeed`, `LiquidDrain`) |
| 2 | Physical parameters — vessel geometry (same as Ex2) + CSTR operating point: `D_PER_H=0.2`, `S_FEED_G_L=5.0`, `Q_L_PER_H = D × V_liq` |
| 3–10 | Species, kinetics, equilibria, ReactionSystem, phases, transfer models, CV — all identical to Example2 |
| 11 | Boundaries: gas feed + pressure vent (unchanged) + **`LiquidFeed`** (substrate inlet) + **`LiquidDrain`** (effluent outlet) + pH controller |
| 12–13 | `Simulation` + run (10 h = 2 HRTs) |
| 14 | Theoretical SS: S* = Ks·D/(μmax−D), X* = Y·(S₀−S*); compared against simulation endpoint |
| 15 | Results summary |
| 16 | 2×2 plots — substrate g/L + SS line, biomass g/L + SS line, dissolved O₂, pH |

**Key design decisions:** Equal Q for feed and drain maintains constant volume.
`LiquidDrain` uses exponential formulation for stability. D = 0.2 /h (0.4× μmax)
gives S* ≈ 3.3 mg/L and X* ≈ 1.8 g/L. No biomass in feed (free-cell operation).

---

#### `demos/model_api/D2Cworkshop/basic_layout/Example1_mtp_well.ipynb`

Closed-system MTP well demo. Same Monod aerobic growth on acetic acid as
`Example2_batch_fermenter`, but in a sealed microtitre plate well with no gas
boundaries. Introduces pulling equilibrium reactions from the package database
rather than declaring them manually.

Structure (14 sections):

| Section | Content |
|---|---|
| 1 | Path setup + imports (`AQUEOUS_DEFAULT`, `AD_BASIC`) |
| 2 | Physical parameters — deep-well plate geometry (200 µL liquid, 500 µL headspace) |
| 3 | Species declarations (`ACETIC_ACID`, `ACETATE_MINUS`, `YEAST`) |
| 4 | Kinetic reaction via `ReactionBuilder.monod_aerobic_growth` |
| 5 | **Database reactions** — inspect `AQUEOUS_DEFAULT.reactions`, extract `eq_water` + `eq_CO2` by label |
| 6 | `EquilibriumReaction`: acetic acid dissociation (domain-specific, declared locally) |
| 7 | `ReactionSystem` assembled from growth + database reactions + dissociation |
| 8 | `GasPhase` + `LiquidPhase` (dry air headspace; dissolved gases at Henry equilibrium) |
| 9 | `EquilibriumTransferModel` for O₂/CO₂/N₂ — no kLa (instantaneous thin-film exchange) |
| 10 | `ControlVolume` — no boundaries (closed system) |
| 11–12 | `Simulation` + run |
| 13 | Summary printout (substrate consumed, biomass grown, O₂ fraction remaining) |
| 14 | 2-row single-column plot: vapour fractions (top) + substrate/biomass/pH (bottom) |

**Key design decisions:**

1. **Database reactions** — `AQUEOUS_DEFAULT.reactions` is iterated and filtered by
   label (`"eq_water"`, `"eq_CO2"`). The NH₄⁺/NH₃ reaction is excluded because
   nitrogen species are not seeded in the liquid phase for this model.

2. **Closed headspace** — no `GasFeed` or `PressureReliefVent`. Total moles of
   each gas species (gas + liquid combined) are conserved; O₂ depletion in the
   headspace is visible in the vapour-fraction plot.

3. **Equilibrium transfer** — `EquilibriumTransferModel` for all three species.
   In a sub-millilitre well the gas-liquid film is thin and well-mixed, so
   instantaneous partitioning is the appropriate limit.

4. **Plot layout** — 2 rows, 1 column (`sharex=True`). Top: O₂/CO₂/N₂ vapour
   mole fractions (%). Bottom: substrate + biomass concentrations (g/L, left
   axis) and pH (right axis via `twinx`).

---

## Session — 2026-06-17

### Notebooks created

#### `demos/builder/batch_fermenter.ipynb`

High-level entry-point demo.  Mirrors `demos/builder/batch_fermenter.py` but
in notebook form, adding a 2×2 matplotlib time-series figure (substrate,
biomass, dissolved O₂, pH).

Structure:

| Section | Content |
|---|---|
| 1 | Path setup + imports (repo-root finder, no hard-coded paths) |
| 2 | Parameters (`USE_PH_CONTROL`, `PH_SETPOINT`, `TAU_H`, `N_STEPS`) |
| 3 | `FermenterBuilder` chain |
| 4 | `build_simulation_and_run` |
| 5 | Final-state summary table |
| 6 | 2×2 time-series plots |

Notable: path setup uses a `_find_root()` helper that walks up from `Path.cwd()`
until `pyproject.toml` is found, so the notebook runs correctly regardless of
which directory the VSCode Jupyter kernel sets as CWD.

---

#### `demos/model_api/D2Cworkshop/basic_layout/Example2_batch_fermenter.ipynb`

Explicit construction demo — the same batch fermenter topology built from
primitives without `FermenterBuilder`.  Each layer the builder hides is shown
in its own markdown+code cell pair.

Final structure (36 cells):

| Section | Content |
|---|---|
| 1 | Path setup + imports |
| 2 | Physical parameters (process values only; Henry constants come from database) |
| 3 | Species declarations — locally declared vs `common_species` imports |
| 4 | Kinetic reaction: `ReactionBuilder.monod_aerobic_growth` (Monod rate + stoichiometry in one call) |
| 5a | `EquilibriumReaction`: water dissociation (H₂O ⇌ H⁺ + OH⁻) |
| 5b | `EquilibriumReaction`: carbonate equilibrium (CO₂ + H₂O ⇌ HCO₃⁻ + H⁺) |
| 5c/5d | `EquilibriumReaction`: acetic acid dissociation + CO₂ cross-phase partition |
| 6 | `ReactionSystem` assembly — shows pre-bucketing by type |
| 7 | `GasPhase` from ideal gas law |
| 8 | `LiquidPhase` with Henry-equilibrated dissolved gases |
| 9 | `KineticGasLiquidLink` with `AD_BASIC.partition_models` |
| 10 | `ControlVolume` assembly |
| 11 | Boundaries (`GasFeed`, `PressureReliefVent`) + `PHController` |
| 12 | `Simulation` |
| 13–15 | Run, results table, 2×2 plots |

**Key design decisions recorded during notebook development:**

1. **Henry constants from `AD_BASIC`** — `partition_models` (`O2`, `CO2`, `N2`)
   are pulled from `AD_BASIC.partition_models` rather than declared as
   hard-coded `mol/L/atm` constants.  `kLa` values remain explicit because
   they are process parameters, not physical properties.  The previous
   `HENRY_MOL_L_ATM` dict and the `* 1000 / 101325` unit conversion are gone.

2. **Water dissociation must be explicit** — `BisectionChemicalEquilibriumEngine` always seeds
   `EquilibriumSet` with a default `pKw = 14.0`, so the simulation runs without
   a water reaction.  However, declaring it explicitly (`rxn_water`) makes the
   assumption visible and allows overriding.

3. **Carbonate equilibrium is required for correct CO₂ speciation** — Without
   `rxn_co2_aq` (CO₂ + H₂O ⇌ HCO₃⁻ + H⁺) in `single_phase_equilibria`, the
   BFS inside `KineticGasLiquidLink.derive_speciation_keys` finds no acid-base
   ladder connected to CO₂.  `_alpha_for("CO2")` then returns `1.0` at every
   step, meaning the link treats all dissolved inorganic carbon as free molecular
   CO₂ and overstates the transfer driving force at any pH where bicarbonate is
   significant.  Adding `rxn_co2_aq` seeds the ladder `[CO₂, HCO₃⁻]` and
   enables the α-correction.

4. **Cross-phase `EquilibriumReaction` vs `HenryPartition` — friction
   resolved** — These previously required two separate declarations for CO₂:
   a `HenryPartition` in `partition_models` (thermodynamic constant) and a
   cross-phase `EquilibriumReaction` (identity routing declaration to seed BFS).
   Under the phase-agnostic species-ID convention the routing reaction is always
   an identity mapping (`"CO2" → "CO2"`), so it carries no real information.
   **Option A implemented** (2026-06-17): `derive_speciation_keys` now seeds
   `speciation_keys[sp] = sp` for every species in `partition_models` that
   lacks an explicit entry.  The cross-phase routing reaction is no longer
   required.  Option B (a unified `GasLiquidPartition` type that bundles both
   the thermodynamic constant and routing) is noted as a future improvement in
   `gas_liquid_link.py`.

5. **`HCO3-` and `OH-` seeded as `0.0`** in the initial `LiquidPhase.n_mol`
   so the speciation engine has writeable slots on the first solve step.

6. **`Species.atoms` passed directly** — `ACETIC_ACID.atoms` and `YEAST.atoms`
   are now passed to `ReactionBuilder.aerobic_growth` without wrapping in
   `dict()`.  See codebase change below.

**Generator script:** `C:\Users\k2473520\AppData\Local\Temp\gen_nb.py` — a
Python script that builds the `.ipynb` JSON programmatically.  Edit this script
and re-run it to regenerate the notebook; do not edit the `.ipynb` directly.

**2026-06-17 update — string stoichiometry format:** Sections 5a, 5b, and 5c
were updated to use the new human-readable string format for `EquilibriumReaction`
stoichiometry (see codebase change below).  The import of `H_plus`, `OH_minus`,
`H2O`, and `HCO3_minus` from `common_species` was removed from the setup cell
because the string parser resolves them automatically.  Locally declared species
(`ACETIC_ACID`, `ACETATE_MINUS`) are still passed explicitly via
`species={...}` to the `rxn_dissoc` reaction.

---

### Codebase changes

#### `src/reactions/builder.py` — `Mapping` annotation for atom dicts

**Change:** `substrate_atoms`, `biomass_atoms`, and `n_source_atoms` parameters
of `ReactionBuilder.aerobic_growth` re-annotated from `Dict[str, float]` to
`Mapping[str, float]`.  `Mapping` added to the `typing` import.  Docstring
parameter types updated to match.

**Why:** `Species.atoms` is stored as `types.MappingProxyType` (set in
`Species.__post_init__` to keep the frozen dataclass hashable).
`MappingProxyType` is a read-only `Mapping` and fully supports `.get()` and
iteration, which is all `aerobic_growth` uses internally.  The previous `Dict`
annotation was too narrow — it implied mutability the function never exercises —
and caused notebook authors to write the unnecessary `dict(sp.atoms)` wrapping.
Changing the annotation to `Mapping` closes the gap without any behaviour
change.

**Files touched:** `src/reactions/builder.py` (import line + three parameter
annotations + two docstring lines).

---

#### `src/reactions/stoichiometry.py` — string stoichiometry parser

**Change:** Added `_parse_stoichiometry(s, species, *, reaction_type=None)`
private function.  Parser accepts strings in the form
`"[coeff] species_id,phase [+ ...] <-> [coeff] species_id,phase [+ ...]"`.
Phase suffixes: `,aq`/`,l` → `"liquid"`, `,g` → `"gas"`, `,s` → `"solid"`.
Standard inorganics are resolved automatically from `common_species`; locally
declared species must be passed via the `species` dict.

Sanity checks raise `ValueError` with descriptive messages for: missing/
ambiguous arrow, arrow direction mismatch with the calling reaction class,
missing phase suffix, unknown phase alias, unresolvable species ID.

Also added `_PHASE_ALIASES`, `_COEFF_RE`, and `_get_common_species()` helper.

**Files touched:** `src/reactions/stoichiometry.py`.

#### `src/reactions/equilibrium.py` and `src/reactions/kinetic.py` — string kwarg

**Change:** Both `EquilibriumReaction.__init__` and `KineticReaction.__init__`
now accept `stoichiometry` as either `list[StoichiometryEntry]` (unchanged) or
`str` (new).  When a string is passed, `_parse_stoichiometry` is called before
any existing validation.  Both constructors also gain an optional keyword-only
`species: dict[str, Species] | None = None` argument.

**Files touched:** `src/reactions/equilibrium.py`, `src/reactions/kinetic.py`.

---

### Packages installed

The following packages were installed into the system Python 3.12 interpreter
(`C:\Users\k2473520\AppData\Local\Programs\Python\Python312\python.exe`) to
support notebook execution in VSCode:

| Package | Version | Purpose |
|---|---|---|
| `ipykernel` | 7.3.0 | VSCode Jupyter kernel |
| `matplotlib` | (latest) | Time-series plots in notebooks |
