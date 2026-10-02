# HPLC column template on the CV graph — Implementation Plan (draft)

> Draft plan, 2026-09-28. No branch, no checklist, no code yet. Written from a
> planning conversation on the `vant-hoff-single-source` branch (at `1291b6f`);
> references to `thermo/temperature_correction.py` assume that phase has
> shipped. Direction set by the repo owner the same day: **the template builds a
> standard CV graph run by `Simulation`, uses the standard backend throughout
> (including the built-in temperature corrections), and supports solid–liquid
> mass transfer, starting with retained species declared as their own
> `Species`.** Items marked *proposed* are this draft's suggestion and still need
> agreement; the set-up questions at the end are open.

## Goal

A `PyOMES/templates/hplc_column/` template that builds an HPLC column as N
cells in series, each cell a `ControlVolume` with a mobile (liquid) and a
stationary (solid) phase, with selectable retention models per solute. The
result is an ordinary `Simulation`: controllers, profiles, recorders,
checkpoints and the conservation monitor all apply. The template replaces the
standalone [`models/vlmodels/hplc/column.py`](../../../../models/vlmodels/hplc/column.py)
once it reproduces that model's results.

## The model

Transport-dispersive column with linear-driving-force uptake, the same physics
as `column.py`. For cell *i* and solute *s*:

- **Mobile phase:** advection from cell *i−1* to *i* at flow Q, axial
  dispersion between neighbours, and uptake to the stationary phase.
- **Stationary phase:** `dq/dt = kf · (q_eq(C, T) − q)`, where `q_eq` comes from
  the solute's retention model and may depend on all solutes' concentrations
  (competitive Langmuir), on temperature, and on pH.
- **Inlet:** mobile phase at flow Q into cell 0, with the sample injected into
  it. **Outlet:** flow Q out of cell N−1; the detector signal is cell N−1's
  mobile-phase concentration.

## Mapping onto the CV graph

| Physical element | Backend object |
|---|---|
| Cell *i* | `ControlVolume(label="cell_000"…)` with phases `"liquid"` and `"solid"` |
| Mobile phase in a cell | `LiquidPhase`, volume ε·A·Δz |
| Stationary phase in a cell | `SolidPhase`, volume (1−ε)·A·Δz, with a sorbent mass ρ_b·A·Δz (needs a mass field, see Prerequisites) |
| Free solute | a `Species` (from a database module, a `ChemistryDatabase`, or user-declared) in `"liquid"` |
| Retained solute | its own `Species`, same atoms, id `<solute>_ads`, in `"solid"` |
| Uptake | one cross-phase `KineticReaction` per retained solute, `"<s>,aq -> <s>_ads,s"`, rate from the solute's retention and uptake objects |
| Retention model | a data object in the partition-model family (see Retention models) |
| Advection | `AdvectiveLink` cell *i* → *i+1*, `Q_L_per_h = Q` |
| Axial dispersion | `DiffusiveLink` between neighbours, `kLa[s] = D_ax,s / Δz²`, `V_eff_L` = the cell's liquid volume |
| Mobile-phase supply | a feed boundary on cell 0 (flow Q, mobile-phase composition) |
| Injection | see Injection below |
| Outlet | `LiquidDrain` on cell N−1, flow Q |
| Column temperature | every phase's `T_K`; an oven programme is a temperature profile over all cells |
| Recycle (later) | an `AdvectiveLink` from cell N−1 to cell 0, switched by a controller |

**Uptake rate, in consistent units.** Each retention model declares its loading
basis:

- **Mass basis** (Langmuir, Freundlich): `q` in mol/kg of sorbent;
  capacity = the cell's sorbent mass.
- **Pore-volume basis** (linear partition `K_D`, size/ion exclusion): `q` in
  mol/L of pore liquid; capacity = ε_p × the cell's solid volume.

The reaction rate is then `kf(T) · capacity · (q_eq − q)` in mol/h, with
`q = n_ads / capacity`. This fixes the unit mix in `column.py` (see Findings).

**Dispersion and cell count.** Advective links are first-order upwind, which
adds numerical dispersion `D_num = u·Δz/2`. Plate count is then
`N ≈ L·u / (2·(D_ax + D_num))`. The factory ports `column.py`'s upwind rule for
choosing the number of cells from a target plate count, so the user can give
`target_plates` or `n_cells`.

**Temperature.** Every temperature-dependent parameter is corrected at the
cell's own `T_K` using the backend's functions:
- `kf` via `arrhenius_factor(Ea/R, T, T_ref)`;
- equilibrium constants (Langmuir `K`, partition `K_D`, solute pKa) via
  `vant_hoff_K` / `vant_hoff_log_K`.

With no activation energy or enthalpy given, the parameter is temperature
independent.

**pH** *(proposed)*. First version: the mobile phase's pH is fixed, taken from
the config or computed once at build time by a speciation engine from the
mobile-phase composition (e.g. 5 mM H₂SO₄). A weak-acid solute's retention then
uses its neutral fraction from its declared acid–base equilibrium, with van 't
Hoff correction, instead of `column.py`'s own `_f_protonated`. Per-cell
speciation (sample acids shifting local pH) is a later extension: it runs an
engine per cell and does not vectorise.

**Injection** *(proposed)*. Profiles are applied once per macro step
(`Simulation._invoke_profiles`), and a boundary's `compute_flux` is not given
the time. A typical injection (20 µL at 0.6 mL/min lasts 2 s) is far shorter
than a sensible macro step. Two options:

1. **Initial plug (first version).** Place the sample in the first cells'
   mobile phase at t = 0, spread over as many cells as the injection volume
   fills. Exact injected mass, no framework change, and a good approximation
   when the plug is short relative to the column.
2. **Time-aware feed (later).** A feed boundary whose composition is a function
   of time, evaluated inside the integrator, which needs boundaries to receive
   `t_h` (a Prerequisites item). It reproduces `column.py`'s
   `injection_inlet`, including extra-column mixing.

## Retention models (backend)

Added to the partition-model family, consistent with `PHENOMENA_PROTOCOL.md` §5
and `MASS_EXCHANGE_ARCHITECTURE.md` §14.1–14.2 (isotherms are a distinct
constraint type on `PartitionModel`; competitive Langmuir is evaluate-at-a-point,
so it can be wrapped safely):

| Model | Form | Basis | Notes |
|---|---|---|---|
| Linear partition | `q_eq = K_D · C` | pore volume | ion exclusion, size exclusion; `partition_ratio()` returns a float |
| Langmuir | `q_eq = q_max·K·C / (1 + K·C)` | mass | `partition_ratio()` returns `None` |
| Competitive Langmuir | `q_eq,i = q_max,i·K_i·C_i / (1 + Σ_j K_j·C_j)` | mass | a `MultispeciesPartitionModel` |
| *Later:* bi-Langmuir, Freundlich, linear solvent strength (gradients) | | | |

Each model:
- holds its parameters plus optional temperature data (reference temperature
  and ΔH for its constants);
- exposes a scalar method and an array method for the equilibrium loading;
- for weak acids, takes the neutral fraction as an input, which the template
  supplies from pH.

The array form is not needed until the performance track, but defining it now
keeps the objects data-declared (see Performance track).

Uptake (`kf`, optional `Ea/R`) is a small data object alongside the retention
model. The first version uses it through the cross-phase `KineticReaction`. The
later liquid–solid transfer interface (`KINETIC_TRANSFER_GENERALIZATION.md`)
reuses the same objects.

## Package layout and API *(proposed)*

```
PyOMES/templates/hplc_column/
    __init__.py
    configs.py     # ColumnConfig, PackingConfig, MobilePhaseConfig, SoluteConfig,
                   # InjectionConfig, DiscretisationConfig
    factory.py     # HPLCColumnFactory.create_simulation(configs) -> Simulation
    builder.py     # HPLCColumnBuilder (chainable), .build_simulation(), .run()
    result.py      # HPLCResult: chromatogram, spatial profiles, peak metrics, plots
```

Configs use chromatography units (cm, mL, mL/min, min, µL, mol/L); the factory
converts to the framework's litres and hours. Sketch:

```python
result = (
    HPLCColumnBuilder()
    .column(length_cm=30.0, diameter_cm=0.78)
    .packing(void_fraction=0.30, particle_porosity=0.50, bulk_density_g_mL=0.40)
    .mobile_phase(flow_mL_min=0.6, T_C=60.0, composition={"H2SO4": 0.005})
    .solute("Glucose", retention=LinearPartition(K_D=0.243), kf_per_min=30)
    .solute("LacticAcid", retention=LinearPartition(K_D=0.674), kf_per_min=40)
    .injection(volume_uL=20.0, sample={"Glucose": 0.025, "LacticAcid": 0.020})
    .discretisation(target_plates=500)
    .run(t_end_min=25.0, n_output=500)
)
result.chromatogram.plot()
result.peaks()   # retention time, plates, resolution, asymmetry per solute
```

`HPLCResult` wraps the `BatchResult` (which already keys `phase_mol` by any
phase name, including `"solid"`):
- the chromatogram is the last cell's liquid concentrations over time;
- spatial profiles are all cells at one time;
- the raw `BatchResult` stays accessible.

Peak metrics:
- **Retention time:** from the first moment of the peak.
- **Plates:** from its variance.
- **Resolution and asymmetry:** standard definitions.

## Prerequisites (framework fixes)

Found in the code on 2026-09-28; each is a general fix, not HPLC-specific.

1. **Rate laws see the solid phase.** `_build_reaction_environment`
   (`control_volume.py:924-1003`) and `SimultaneousAdaptiveSolver`'s
   environment (`solvers.py:784-797`) expose only liquid concentrations; the
   uptake rate needs `q`. Coordinate with, or absorb,
   `REACTION_ENVIRONMENT_PHASE_EXPOSURE.md`, which proposed exactly this
   phase-keyed extension but had no solid-phase use case yet.
2. **`SolidPhase` gets a mass** (or bulk density): today it has only `V_L`
   (`phases.py:478-525`). Design together with
   `NR_PRECIPITATION_CV_INTEGRATION.md`, which also touches `SolidPhase`.
3. **The monolithic right-hand side becomes complete.** `_build_rhs`
   (`system_solver.py:630-736`) sums reactions, links and controllers but not
   boundary fluxes or internal interfaces. `SimultaneousAdaptiveSolver` already
   does both (`solvers.py:756-779`).
4. **No silent mass loss.** A link flux into a slot missing from the state
   layout is skipped (`system_solver.py:693-698`). Raise instead. (The template
   declares every species in every cell, so it won't trigger this, but a silent
   loss must not be possible.)
5. **Jacobian sparsity.** `MonolithicODESolver` passes none to `solve_ivp`, so a
   finite-difference Jacobian costs one right-hand-side evaluation per state
   (600 at 100 cells × 6 species). Links and reactions give the pattern directly.
   Moved ahead of the template because it decides whether the reference runs are
   practical.
6. **A multi-cell temperature profile.** `TemperatureRamp`
   (`control/cv_profiles.py`) targets one CV. Either extend it to a set of CVs
   or have the template add one per cell.
7. *Later, for time-aware injection:* boundaries receive `t_h`.

## Validation

1. **Analytical retention times.** For linear retention, the peak's first moment
   matches `t_R = (ε + (1−ε)·ε_p·K_D)·V_col / Q`, and the equivalent Langmuir
   form in the linear limit.
2. **Against `column.py`.** Same parameters, upwind advection and the same cell
   count in both; chromatograms agree to solver tolerance. The injection has to
   match too: an initial plug in the template against `column.py` started from
   the same initial condition.
3. **Mass balance.** Injected moles equal eluted plus retained, checked by the
   conservation monitor and by integrating the outlet flux.
4. **Temperature.** A `K_D` with ΔH shifts the retention time as `vant_hoff_K`
   predicts; a ramp across cells gives the expected shift.
5. **pH.** A weak acid elutes earlier at higher pH, reproducing `column.py`'s
   pH tests.
6. **Port `test_hplc_column.py`** to the template, then retire
   `models/vlmodels/hplc/` and its tests *(proposed; see question 8)*.

## Performance track (after the template)

The CV graph is slower than the standalone column. A 2026-09-28 scratch
benchmark (100 cells × 3 solutes, liquid→solid reaction per solute, advective
and dispersive links, 600 states) measured:
- the CV graph's right-hand side (`_build_rhs`) at 2.29 ms per evaluation, of
  which 1.10 ms is per-CV reactions;
- `column.py`'s vectorised right-hand side at 0.045 ms per evaluation.

The sandbox caps CPU, so only the ratio (about 50×) is meaningful. The plan to
close the gap is a general system-solver feature, to get its own design note
when the template is working:

- **Compile, don't walk.** The system solver compiles the graph once per run.
  CVs with the same structure (phase set, species, reactions) form a group,
  whose block of the `StateVector` layout (`state_vector.py:73-80`) reshapes to
  cells × slots. Data-declared models (retention, uptake, temperature) are
  evaluated as one array expression per group; parameters may differ per cell.
- **Transport as a matrix or a chain link.** Linear links become one sparse
  matrix (`ImplicitTransportSystemSolver` already assembles it). A chain link
  spanning an ordered set of CVs applies `numerics/spatial.py`, which brings
  back TVD advection.
- **Padding is solver-internal and within a phase set.** Species or reaction
  differences inside a group are padded only in the solver's working arrays:
  - gather from the state vector through an index map, with a zero sentinel for
    missing slots;
  - scatter back real slots only;
  - mask reactions per cell (0/1) where they are absent.

  Padding never appears in a CV's `n_mol`, the integrated state or results.
  Phases are never padded.
- **Reachability pass.** A `Simulation` pre-flight step declares, for every
  (CV, phase), each species that can ever occur there, at 0.0. It starts from
  the initial `n_mol` and repeats until nothing changes:
  - feed compositions add to their target phase;
  - every species in any attached reaction's stoichiometry is present in the
    phase it names (over-inclusion is the safe side);
  - black-box flux mappings add their targets;
  - interfaces carry a species across;
  - links carry a species from source to sink, and diffusive links both ways.

  Declared species are real slots; padding covers only what a cell can never
  hold.
- **Validation.** Against the template's plain-graph results, plus one
  deliberately different case (a microplate or a parameter ensemble) so the
  design does not take on HPLC-only assumptions.
- **Limits.** Opaque `rate_fn` lambdas run per cell; speciation engines solve
  per CV and are not batched.

## Phases

1. **Framework prerequisites** 1–6 above (one phase, or split 1–2 and 3–6).
2. **Retention models** in the backend: linear partition, Langmuir, competitive
   Langmuir, with temperature data and scalar and array forms, plus tests.
   Independent of phase 1.
3. **Template core:**
   - `configs.py`, `factory.py`, `builder.py` for isocratic runs with an
     initial-plug injection;
   - `result.py` with the chromatogram and peak metrics;
   - validation items 1–5.
4. **Close-out:** port the tests, write the tutorial in
   `docs/tutorials/templates/`, retire `models/vlmodels/hplc/`, and update
   `docs/architecture.md`, `models/README.md` and the README test table.
5. **Extensions** (each optional, in any order): time-aware injection with
   extra-column mixing; recycle and peak shaving; per-cell speciation;
   gradient elution.
6. **Performance track:** its own design note and phases.

Phase 3 needs 1 and 2. Phase 6 needs 3.

## Findings

- **The Langmuir units in `column.py`.** The phase ratio multiplies `ρ_b`
  (g/mL) by `q` (documented as mol/g) and treats the product as mol/L
  (`column.py:1210`, `1241-1245`). So `q_max` behaves as mol/kg. Its tests pass
  because the predicted and simulated retention times share the formula.
- **The shipped "HPLC trigger bundle" notes.** `SNAPSHOT_PHASE_AGNOSTIC.md` and
  `PROPERTY_CALCULATORS_PHASE_AGNOSTIC.md` sit in `shipped/` but read partly as
  proposals. The property-calculator runner does dispatch by `phase_key`
  (`control_volume.py:891-920`). The recorder still reads
  `phases.get("gas")`/`.get("liquid")` in three places (around `recorder.py:507`,
  `717`, `978`). Check what is still open before phase 3.
- **The built-in phase transfer is gas–liquid only.** `_build_transfer_link`
  requires a `GasPhase` and a `LiquidPhase` (`control_volume.py:46-110`). This
  is why uptake goes through a cross-phase `KineticReaction` for now.

## Open questions

**Set-up information**

1. **Retention models in the first version.** Linear `K_D`, Langmuir and
   competitive Langmuir as proposed, or others? Which modifiers (pH,
   temperature) from the start?
2. **Isocratic or gradient elution.** Gradients make the mobile-phase
   composition a transported variable feeding the retention model (e.g. linear
   solvent strength, log k = log k_w − Sφ).
3. **Column model detail.** Transport-dispersive only, as proposed, or also
   equilibrium-dispersive and general-rate (pore diffusion)?
4. **Inputs.** Is the config split above (column, packing, mobile phase,
   solutes, injection, discretisation) right? Is a detector response (RI/UV
   factors) wanted?
5. **Parameter sources.** User-supplied per solute, or a solute–column database
   (e.g. an Aminex HPX-87H `K_D` table held in a `ChemistryDatabase`)?
6. **Units.** Chromatography units at the config boundary, converted to hours
   and litres by the factory, as proposed?
7. **Injection.** Initial plug first, time-aware feed later, as proposed?
8. **The standalone model's fate.** Retire `models/vlmodels/hplc/` once the
   template reproduces it, with no shims, as proposed?

**Framework**

9. **Controllers and profiles that add material** (e.g. a pH controller dosing
   `NaOH` adds `Na+`): do they declare what they can add? This is needed by the
   reachability pass; not yet checked.
10. **Reference solver for phase 3.** Monolithic (once prerequisite 3 lands), or
    implicit transport, which splits transport from chemistry at first order in
    time and may smear peaks?

## Cross-references

- [`REACTION_ENVIRONMENT_PHASE_EXPOSURE.md`](REACTION_ENVIRONMENT_PHASE_EXPOSURE.md)
  — prerequisite 1.
- [`NR_PRECIPITATION_CV_INTEGRATION.md`](NR_PRECIPITATION_CV_INTEGRATION.md) —
  also changes `SolidPhase`; design prerequisite 2 with it.
- [`KINETIC_TRANSFER_GENERALIZATION.md`](KINETIC_TRANSFER_GENERALIZATION.md),
  [`PHENOMENA_PROTOCOL.md`](PHENOMENA_PROTOCOL.md) — the later liquid–solid
  transfer interface.
- [`../shipped/IMPLICIT_TRANSPORT.md`](../shipped/IMPLICIT_TRANSPORT.md),
  [`../shipped/MONOLITHIC_ODE.md`](../shipped/MONOLITHIC_ODE.md),
  [`../shipped/SYSTEM_SOLVER_PROTOCOL.md`](../shipped/SYSTEM_SOLVER_PROTOCOL.md)
  — the system solvers the template runs on.
- [`../shipped/STIRRED_TANK_TEMPLATE_MIGRATION.md`](../shipped/STIRRED_TANK_TEMPLATE_MIGRATION.md)
  — the template pattern, and the "1D flow-reactor template" idea this could
  later generalise into.
- [`../../ideas/MASS_EXCHANGE_ARCHITECTURE.md`](../../ideas/MASS_EXCHANGE_ARCHITECTURE.md)
  §5.3, §14 — liquid ↔ solid adsorption, isotherms as a constraint type.
