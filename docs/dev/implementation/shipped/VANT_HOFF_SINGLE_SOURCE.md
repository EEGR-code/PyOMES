# One source for the van 't Hoff relation — Design Note

> **Status: Shipped 2026-09-28** — implemented on branch
> `vant-hoff-single-source`, merged into `main` via `git merge --no-ff` as commit
> `90c5915`, tagged `vant-hoff-single-source-shipped`. The checklist,
> [`VANT_HOFF_SINGLE_SOURCE_CHECKLIST.md`](VANT_HOFF_SINGLE_SOURCE_CHECKLIST.md),
> records the implementation, every measured shift and the deviations settled
> along the way (among them: the 0.01 K skip replaced by the kernel's 1e-10 K rule
> rather than removed outright); where the two differ, the checklist wins.

> Design discussion, 2026-09-28. Written
> from a planning conversation that started from the `OPEN_WORK.md` entry "A
> third van 't Hoff copy in `reactions/equilibrium/constraint.py` differs from
> `thermo` in the last bit". Direction approved by the repo owner the same day:
> **one kernel in `PyOMES/thermo/`, which every copy calls instead of writing the
> formula out; van 't Hoff for all equilibrium constants including Henry
> constants; Clausius–Clapeyron and Arrhenius as their own named functions on the
> same kernel; `chemistry/` allowed to import `thermo/`.** Numbers may move, as
> long as the formulation is correct and every shift is attributable either to
> floating-point arithmetic or to a named, deliberate change of formulation.

## The problem

The same relation,

    X(T) = X(T_ref) · exp(−(E/R) · (1/T − 1/T_ref))

is written out by hand in 15 places (searched 2026-09-28, on `main` at
`314fba9`). It is the van 't Hoff equation when X is an equilibrium constant and
E the reaction enthalpy, Clausius–Clapeyron when X is a vapour pressure and E the
enthalpy of vaporisation, and Arrhenius when X is a rate constant and E the
activation energy. Every copy is mathematically correct, with the right sign
convention for its quantity, but they differ in arithmetic order, in how they
convert between K and log10 K, and in their edge-case rules, so the same
quantity can come out differently by path. For a mineral, the NR engine computes
Ksp(T) with `thermo.vant_hoff_log_K` (`nr/engine.py:633`) while
`KspEquilibrium` computes it with `constraint.vant_hoff_log_K`
(`interphase.py:472`), and the two differ in the last bit. Each copy is also a
place a future change (to R, to an edge case) has to be repeated.

## Inventory

| # | Where | Corrects | Form | Edge rules |
|---|---|---|---|---|
| 1 | `thermo/equilibrium_constants.py` `vant_hoff_delta_ln_K`, `vant_hoff_K`, `vant_hoff_log_K` | K, log10 K | canonical: `-(dH/R)·(1/T − 1/T_ref)`; log10 via `np.log10(np.e)` | K: skip if K ≤ 0 or non-finite, dH non-finite or \|dH\| < 1e-30, T ≤ 0 or non-finite. log K: skip if dH is None or \|dH\| < 1e-30, or \|T − T_ref\| < 1e-10 |
| 2 | `reactions/equilibrium/constraint.py:63` `vant_hoff_log_K(constraint, T)` | log10 K of any `EquilibriumConstraint` | same maths; log10 via `1/math.log(10)` | same as #1 log K |
| 3 | `reactions/equilibrium/plots.py:99` `plot_vant_hoff` | log10 K over a T grid | `dH/(R·ln10)·(1/T_ref − 1/T)` | only when dH is not None |
| 4 | `reactions/equilibrium/plots.py:322` `plot_speciation` | pKa at T | as #3 | as #3 |
| 5 | `chemical_equilibrium/engines/bisection/equilibria.py:141` `EquilibriumDef.pKas_at_T` | pKa (Bisection engine) | via Ka: `Ka·exp(-dH/R·(…))`, then `-log10(max(Ka, 1e-30))` | skip if `correction == "none"`, **\|T − T_ref\| < 0.01 K**, or \|dH\| < 1e-10 |
| 6 | `…/bisection/equilibria.py:173` `WaterDef.Kw_at_T` | Kw (Bisection engine) | as #5 | as #5 |
| 7 | `thermo/framework.py:88` `ThermoFramework.pKa_at_T` | pKa | as #5 | \|dH\| < 1e-12 or **\|T − T_ref\| < 0.01 K** |
| 8 | `thermo/framework.py:102` `ThermoFramework.Kw_at_T` | Kw, ΔH = 55 900 J/mol hard-coded | as #5 | **\|T − T_ref\| < 0.01 K** |
| 9 | `chemistry/partition.py:110` `_kH_mol_L_atm_from_ref` | Henry constant (also used by `HenryEquilibrium`) | `exp(dlnH·(1/T − 1/T_ref))`, Sander's `dlnH = −ΔH_sol/R` (K) | none |
| 10 | `reactions/equilibrium/interphase.py:285` `RaoultEquilibrium.P_sat` | water vapour pressure | Clausius–Clapeyron, `-dH_vap/R·(…)` | none |
| 11 | `core/phases.py:98` `water_vapour_P_sat_atm` | water vapour pressure | Clausius–Clapeyron, `(ΔH_vap/R)·(1/T_ref − 1/T)`, ΔH_vap/R = 5290 K | none |
| 12 | `models/vlmodels/adm1/base.py:475` `_arrhenius` | rate factor | Arrhenius, `Ea_R·(1/T_ref − 1/T)` | returns 1 if Ea_R ≤ 0 or T ≤ 0 |
| 13 | `…/adm1/base.py:479, 484` `_pKa_NH4`, `_pKa_H2S` | pKa | as #5, rounded R | none |
| 14 | `…/adm1/bsm2.py:331` `_pKa_NH4` | pKa | as #13 | none |
| 15 | `…/adm1/bsm2_direct.py:277` `_pKa_NH4` | pKa | as #13 | none |

Callers: #1 is used by both engines (NR tableau and precipitation, Bisection's
`acid_base.py`). #2's only production caller is `KspEquilibrium.equilibrium_a_moles`,
reached only when a gas-liquid link uses a Ksp model as a partition model, which
nothing in the repo does; otherwise tests. #5–6 are the Bisection engine's
per-solve pKa/Kw. #7–8 are used only by tests and `chemistry_database.ipynb`.
#9 is used by `MultispeciesVLEPartition` and `HenryEquilibrium`. #11 is the
headspace water vapour pressure; #10 is `RaoultEquilibrium`.

**Of the two log10(e) constants**, `np.log10(np.e)` (#1) is the correctly rounded
double of 1/ln 10 (error 1.1e-17); `1.0/math.log(10.0)` (#2) is 1 ulp low (error
6.6e-17). Folding #2 into #1 moves towards the correct value.

## Design

- **One kernel in a new `thermo/temperature_correction.py`**, which replaces
  `thermo/equilibrium_constants.py` (its four importers are updated; no shim). It
  is written in terms of E/R in kelvin, the quantity every copy either has or
  computes first:

      ln_correction(E_over_R_K, T_K, T_ref_K) = -E_over_R_K * (1/T_K - 1/T_ref_K)

  Callers that hold an enthalpy in J/mol compute `dH / R` and pass it; callers
  that already hold a value in kelvin (Sander's `dlnH`, `ΔH_vap/R`, `Ea/R`) pass
  it directly, so no R is divided out and multiplied back.
- **Named wrappers on the kernel**, each keeping its own documented edge rules.
  In `thermo/temperature_correction.py`: `vant_hoff_K`, `vant_hoff_log_K` (moved
  from `equilibrium_constants.py`), a Henry-constant function, and
  `clausius_clapeyron`. In `reactions/kinetic/rate_laws.py`, where kinetics lives:
  `arrhenius_factor`, calling the kernel (so rate laws that gain temperature
  dependence later have it to hand). These are the only public ways to apply the
  relation; nothing but the kernel writes `(1/T − 1/T_ref)`.
- **Callers keep only what is theirs**: unit conversion (Sander mol m⁻³ Pa⁻¹ to
  mol L⁻¹ atm⁻¹), pKa ↔ Ka, clamps such as `max(Ka, 1e-30)`, and their own data
  (reference values, enthalpies, the ADM1/BSM2 rounded R).
- `constraint.vant_hoff_log_K(constraint, T)` becomes a one-line wrapper that
  unpacks the constraint and calls `thermo`, renamed so that two functions no
  longer share the name `vant_hoff_log_K` with different signatures.
- `chemistry/` may import `thermo/`: the rule in
  `tests/standalone/test_package_layering.py` changes from `{chemistry, units}` to
  `{chemistry, thermo, units}`. `thermo/` imports only `units`, so the order
  units → thermo → chemistry stays acyclic.
- A guard test, in the style of `test_gas_constant_single_source.py`, fails if
  the pattern `1/T − 1/T_ref` appears in `PyOMES/` or `models/` outside
  `thermo/temperature_correction.py`.

## Expected numerical changes

Checked by reading each copy against the kernel; each is measured before and
after in its checkpoint (maximum absolute and relative shift over a fixed grid of
temperatures and enthalpies, plus the full suite and the BSM2 sentinels).

**Arithmetic only (floating-point rounding):**

- #2 → #1: log10(e) moves by 1 ulp, towards the correct value. The `OPEN_WORK.md`
  entry measured this at up to ~4e-15 in log10 K in about a fifth of 20,000 random
  cases; re-measured in checkpoint 1.
- #3, #4: dividing by `R·ln 10` instead of multiplying by log10(e) after dividing
  by R; last-bit differences.
- #5–8 routed through `vant_hoff_K` keep their K-space path; `(-dH)/R` equals
  `-(dH/R)` exactly in IEEE arithmetic, so the exponent is bit-identical and any
  shift comes from where the conversion happens.
- #9–15: bit-identical by construction where the caller passes a kelvin value
  (`-(-dlnH)·x == dlnH·x` and `(1/T_ref − 1/T) == −(1/T − 1/T_ref)` exactly); the
  measurement confirms it.

**A deliberate change of formulation (not arithmetic):**

- #5–8 skip the correction when T is within **0.01 K** of T_ref, a step
  discontinuity: at ΔH = 50 kJ/mol, T = T_ref + 0.0099 K gives no correction,
  where the exact value is 2.9e-4 in log10 K, while T_ref + 0.0101 K gives the full
  3.0e-4. The kernel skips only at
  |T − T_ref| < 1e-10 K. Adopting it makes the relation exact at every
  temperature, and changes results only for temperatures within 0.01 K of the
  reference (298.15 K for every stock reaction). BSM2 runs at 308.15 K and is
  unaffected by this part.

**Not changed:** the ADM1/BSM2 rounded R (`_R_J = 8.31446`) stays in `models/` as
data the caller passes (`51965.0 / _R_J`); whether the benchmark fixes R is the
separate `OPEN_WORK.md` entry "ADM1 / BSM2 use a rounded gas constant". The two
different water vapour-pressure parameter sets (#10: 0.03169 atm and 44 011 J/mol;
#11: 0.0313 bar and ΔH_vap/R = 5290 K) are data, not formula; they stay as they are
and are covered by the `OPEN_WORK.md` entry "No single source for water's physical
constants".

## Out of scope

- Unifying the reference data (R in `models/`, the water constants).
- The other fundamental constants in "Sweep the package for each fundamental
  constant" (273.15, 298.15, 101325, …); this phase touches only the relation.
- The pKa ↔ Ka ↔ log K conversions themselves, beyond calling the kernel.

## Checkpoints (proposed)

1. **Kernel, and the reactions copies.** `thermo/temperature_correction.py` with
   the kernel and its wrappers, replacing `equilibrium_constants.py` (importers
   updated); `arrhenius_factor` in `reactions/kinetic/rate_laws.py`; #2 as a
   renamed wrapper; #3–4 call it. Tests for each wrapper.
2. **Bisection engine** (#5–6), including the 0.01 K change, measured.
3. **`ThermoFramework`** (#7–8), same change.
4. **Henry and vapour pressure** (#9–11), with the layering-rule change.
5. **Models** (#12–15): ADM1's `_arrhenius` calls `arrhenius_factor`; the pKa
   functions call the kernel, passing their own E/R (with their own rounded R).
6. **Guard test and close-out**: the single-source guard; delete the `OPEN_WORK.md`
   entry this phase resolves and trim the van 't Hoff paragraph of "Sweep the
   package for each fundamental constant"; docs.

## Decisions (settled 2026-09-28)

1. **Module: `thermo/temperature_correction.py`**, replacing
   `equilibrium_constants.py`. Every function in it corrects a value known at a
   reference temperature to another temperature; the old name fits van 't Hoff and
   Henry but not Clausius–Clapeyron. The relations' own names stay in the function
   names and docstrings, so a search for "van 't Hoff" still finds them.
2. **`arrhenius_factor` lives in `reactions/kinetic/rate_laws.py`**, calling the
   kernel: Arrhenius is kinetics, so it sits where kinetics lives and where future
   temperature-dependent rate laws would need it (none has any temperature
   dependence today). The formula itself exists only in the kernel.
   `reactions/kinetic/` may import `thermo/` under the layering rules.
3. **`ThermoFramework.pKa_at_T` / `Kw_at_T` stay, as thin callers.** Only tests and
   one notebook use them, and `Kw_at_T` hard-codes ΔH = 55 900 J/mol and
   Kw = 1e-14; removing them is an API change for another day.
4. **Scope** (from the direction above): van 't Hoff for every equilibrium constant,
   including Henry constants; Clausius–Clapeyron and Arrhenius as their own named
   functions on the kernel; `models/` routed through the kernel with their own data.
5. **Numbers may move** when the formulation stays correct and each shift is
   attributable to floating-point arithmetic or to the one named change of
   formulation (dropping the 0.01 K skip).
6. **`chemistry/` may import `thermo/`.**
