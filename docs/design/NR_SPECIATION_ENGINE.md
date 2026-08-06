# Extended Newton-Raphson Speciation Engine — Design Note

**Date:** 2026-06-23
**Status:** Implemented — shipped 2026-06-23, tag `nr-speciation-engine-shipped`.

---

## Motivation

The current `BisectionChemicalEquilibriumEngine` solves equilibrium chemistry by reducing the
system to a single scalar equation: the charge balance in [H⁺].  This works
because the acid-base ladders in the present chemistry (VFAs, carbonate,
ammonia, sulfide) are mutually independent — they interact only through the
shared proton pool, so each ladder contributes a scalar charge term that is a
function of [H⁺] alone.

This assumption breaks for two classes of chemistry that are relevant to
wastewater and AD modelling:

- **Metal complexation** — Zn²⁺ + CO₃²⁻ ⇌ ZnCO₃ couples the zinc and
  carbonate ladders: carbonate's alpha fractions depend on [Zn²⁺] and vice
  versa.  The system can no longer be reduced to a 1D problem.
- **Precipitation/dissolution** — struvite (MgNH₄PO₄), calcite (CaCO₃),
  vivianite (Fe₃(PO₄)₂), iron sulfide (FeS) — each introduces a new
  phase whose presence or absence is itself an unknown.

The extended NR engine is a generalised speciation backend that handles
arbitrary reaction networks by solving a full Newton-Raphson system in
log-activity space.  It is intended to coexist with the current engine as a
user-selectable alternative, not to replace it.

---

## Mathematical approach

### Unknowns

Choose one **master species** per independent chemical component.  Their
log-activities `log aₖ` are the unknowns.  The number of masters equals the
rank of the formula matrix (elements × species), which equals the number of
independent conservation laws.

### Secondary species

Every non-master species `j` is expressed analytically in terms of the master
log-activities via its formation reaction:

```
log cⱼ = log Kⱼ + Σₖ νⱼₖ · log aₖ
```

where νⱼₖ is the stoichiometric coefficient of master `k` in the formation
reaction for species `j`, and log Kⱼ is the cumulative formation constant
(product of individual step constants along the path from masters to `j`).
These coefficients are read directly from the declared `EquilibriumReaction`
stoichiometry — no new thermodynamic data is required.

### Residual equations

For `m` master species, write `m` residual equations:

```
Rᵢ = Σⱼ sᵢⱼ · cⱼ(a) − Cᵢ,total = 0    for i = 1 … m−1   (mass balances)
Rₘ = Σⱼ zⱼ · cⱼ(a)             = 0    (charge balance — closes H⁺)
```

`sᵢⱼ` is the number of moles of component `i` contributed by one mole of
species `j` (a column of the formula matrix).  `Cᵢ,total` is the known total
concentration of component `i`, summed from `phase.n_mol` at the start of
each solve.

### Jacobian

The Jacobian is analytic:

```
∂Rᵢ / ∂(log aₖ) = Σⱼ sᵢⱼ · νⱼₖ · cⱼ
```

For the expected scale (~10–20 master species), this is a dense matrix solved
with a standard LU factorisation at each Newton step.  Activity corrections
(Davies model) enter through the activity coefficient γⱼ and create an outer
ionic-strength iteration, identical to the current engine's fixed-point loop.

### Solve loop

Iteration runs in `log a` space (not `a` space) for numerical stability.  A
backtracking line search damps the Newton step when the residual norm would
increase.  Warmstarting from the previous timestep's solution is essential for
performance inside `cv.advance()` — consecutive steps are nearly identical and
typically converge in 2–4 Newton iterations.

---

## Master species selection

### Graph-based automated selection

The master species are derived automatically from the declared reactions via a
graph traversal — no user annotation is required for standard acid-base and
metal-complexation chemistry.

**Step 1 — Build the reaction graph.**
Nodes are species IDs; each `EquilibriumReaction` adds undirected edges
between all its participants (H₂O and other solvents are excluded).  H⁺ is
pre-registered as a universal master and is not part of the graph.

**Step 2 — Find connected components (BFS).**
Each component is a set of species that can interconvert through the declared
reactions.  An isolated species (no declared equilibria) becomes its own
singleton component and is its own master — it is treated as a conservative
tracer.

**Step 3 — Select one master per component.**
Within each component, inspect the role of each species across the reactions
in that component:

- A species that appears **only as a reactant** (stoichiometric coefficient
  < 0) and never as a product is a DAG source.  This is the natural master —
  it is the most protonated (or least complexed) form of the component.
- If no unique source exists (e.g. a cyclic reaction network), fall back to:
  highest H-atom count from `Species.atoms`, then highest charge.

This selects the correct master for each common case:
- Acid-base ladders: most protonated form (H₃PO₄, NH₄⁺, H₂S, etc.)
- Metal systems: bare ion (Zn²⁺, Ca²⁺, Fe²⁺) since it appears only as
  a reactant in complexation reactions
- Cross-component ligand species (OH⁻, CO₃²⁻) are intermediates, not
  sources — they correctly remain secondary species

**Step 4 — Derive tableau coefficients.**
For each non-master species, traverse the reaction graph from the relevant
master(s) to the species, accumulating stoichiometric coefficients νⱼₖ and
summing log K values along the path.  Cross-component species (e.g., ZnOH⁺)
acquire coefficients from multiple masters:

```
Zn²⁺ + OH⁻ ⇌ ZnOH⁺   (log K_ZnOH)
OH⁻ = H₂O − H⁺        (via Kw: log[OH⁻] = log Kw − log[H⁺])

→  log[ZnOH⁺] = (log K_ZnOH + log Kw) + 1·log[Zn²⁺] − 1·log[H⁺]
                  └── log K' ──┘            └── ν(Zn²⁺) ──┘   └── ν(H⁺) ──┘
```

All coefficients are derived from the existing `EquilibriumReaction`
stoichiometry fields.

**Step 5 — Validate.**
Construct the submatrix of the formula matrix at the selected master columns
and check its rank equals the number of elements.  If two masters are linearly
dependent (same chemical component expressed differently), raise a descriptive
error identifying the pair.

### Override via `total_id`

Where the heuristic would misfire, the existing `total_id` field on
`EquilibriumReaction` serves as an explicit override — if set, it takes
precedence over the graph-derived selection for that component.  For redox
systems (Fe²⁺ vs Fe³⁺ ambiguity), explicit declaration via `total_id` is the
expected path until a pe-tracking mechanism is added.

---

## Explicit H⁺, OH⁻, and H₂O in `n_mol`

In the current engine, H⁺ and OH⁻ are implicit — pH is the solver output but
the corresponding molar counts are not written to `phase.n_mol`.  This makes
it difficult to inspect the charge balance directly or trace proton accounting
during debugging.

The NR engine writes these as **engine-owned derived species** after every
solve:

```python
n_H_plus   = (a_H_plus / gamma_H)       * V_L   # mol
n_OH_minus = (Kw / a_H_plus / gamma_OH) * V_L   # mol
n_H2O      = (rho_water / M_water)      * V_L   # mol  (~55.5 × V_L)
```

"Engine-owned" means the speciation engine overwrites these entries on every
call, regardless of what kinetics may have written to them.  This is the same
treatment as all other derived species (HCO₃⁻, CO₃²⁻, etc.) under the
existing `_refresh_derived` mechanism.

The immediate benefit: the charge balance is checkable by inspection at any
point in the integration:

```python
charge_residual = sum(sp.charge * phase.n_mol[sp.id]
                      for sp in all_species if sp.id in phase.n_mol)
# should be ~0 after any speciation solve
```

**Why not track H⁺ as a kinetic state variable?**
If H⁺ were a conserved species updated by kinetic rates, every kinetic
reaction that produces or consumes protons would need to declare explicit H⁺
stoichiometry — a significant burden on model authors and a source of
mass-balance errors.  More fundamentally, the total-H mass balance is
numerically dominated by water (~55.5 mol/L), making H⁺ tracking a small
perturbation on a large background — precision loss is unavoidable in
double-precision arithmetic.  The charge balance avoids this by never
constructing a proton mass balance across the solvent.

---

## Precipitation (deferred, designed for)

Precipitation is not in the initial scope but the design anticipates it.
The outer phase-detection loop wraps the NR solve:

1. Solve without any active solids.
2. For each candidate precipitate, evaluate the ion activity product (IAP)
   against its Ksp from the tableau.
3. If IAP > Ksp: add the solid mole count as a new unknown, add the Ksp
   saturation constraint, adjust mass balances to include solid moles,
   re-solve.
4. If any active solid has a negative mole count, remove it and re-solve.
5. Repeat until no phase changes occur.

For the target precipitates in bioprocess chemistry (struvite, calcite,
vivianite, FeS, HAP), the candidate set is small and known — combinatorial
explosion of phase combinations is not a practical concern.

---

## Implementation plan

### New files

| File | Contents |
|---|---|
| `src/chemical_equilibrium/nr_engine.py` | `NRChemicalEquilibriumEngine`: `from_reactions()` factory, `.solve()` |
| `src/chemical_equilibrium/nr_tableau.py` | `NRTableau` dataclass; graph traversal; master selection; tableau derivation; basis validation |
| `src/chemical_equilibrium/nr_solver.py` | `solve_nr()`: Newton-Raphson loop, backtracking line search, warmstart cache |

### Modified files

| File | Change |
|---|---|
| `src/reactions/reaction_system.py` | Add `solver: str = "charge_balance"` parameter; route `engine` property to `NRChemicalEquilibriumEngine` when `solver="newton_raphson"` |

### Unchanged

The existing `BisectionChemicalEquilibriumEngine`, `EquilibriumSet`, `acid_base.py`, and all
activity model code are untouched.  The `ReactionSystem.engine` property
returns an object satisfying the same `.solve(phases=...)` interface
regardless of which backend is active.

### Validation strategy

Run both engines on the same `ReactionSystem` for the standard BSM2 chemistry
(VFAs, carbonate, ammonia, sulfide) and assert that pH, species
concentrations, and charge residual agree to within solver tolerance.  This
serves as the regression suite before the NR engine is used for any chemistry
the 1D engine cannot handle.

---

## Scope and deferrals

**In scope for initial implementation:**
- Full NR solve for arbitrary single-phase acid-base + metal complexation
  networks
- Automated master species selection (graph heuristic + `total_id` override)
- H⁺, OH⁻, H₂O written to `n_mol` as engine-owned derived species
- Warmstarting and activity-coefficient outer iteration
- Basis validation at construction time

**Deferred:**
- Precipitation / phase detection outer loop
- Redox (pe as explicit master); requires explicit user declaration of iron /
  sulfur oxidation state convention
- Gibbs energy minimization backend (distinct design path; preferable if
  precipitation scope grows to many competing solids)
- Thermodynamic database layer (required before adding metal complexation
  constants; a curated set for bioprocess-relevant metals is a prerequisite)
