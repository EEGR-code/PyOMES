# Multi-component complexation and precipitation: implementation plan

> **Status: partial skeleton, stalled 2026-06-29.** Commit `23e0cb0`
> ("feat: MULTICOMPONENT_COMPLEXATION — standalone multi-component
> speciation pathway") added
> [`src/chemical_equilibrium/multicomponent/`](../../src/chemical_equilibrium/multicomponent/)
> — `components.py`, `reactions.py`, `tableau.py`, `network.py`,
> `inventory.py`, `solver.py`, `residuals.py` (353 lines total) —
> **directly to `main`**, bypassing this repo's one-branch-per-phase
> convention (see [`README.md`](README.md)'s branching section). Only
> `components.py` has tests (8, in
> `tests/standalone/test_multicomponent_components.py`); the other six
> modules are untested. No `MultiComponentEquilibriumEngine` facade
> class exists yet (the entry point this plan proposes below). No
> checklist file was ever written for this work. Notebook 08
> (`08_iron_oxidation_and_precipitation.ipynb`, referenced throughout
> this plan as the acceptance target) was not touched by that commit —
> `05_precipitation_equilibrium.ipynb` was updated instead.
>
> Discovered and documented 2026-07-09 during a scoping review of
> upcoming work; this banner records the drift found, not a decision to
> resume or abandon the work. If picked up: retroactively branch
> `multicomponent-complexation` off `main`, write a checklist capturing
> what already exists against this plan's Phase 1 scope, and continue
> from there rather than starting over — the existing skeleton follows
> the plan's proposed module layout (§ "New standalone pathway") and
> the compatibility rule below has held (nothing in
> `BisectionChemicalEquilibriumEngine`/`NRChemicalEquilibriumEngine` was
> touched).

## Purpose

Extend the chemistry capability to model coupled aqueous complexation and
equilibrium precipitation in systems such as the iron/phosphate/citrate growth
medium used by `08_iron_oxidation_and_precipitation.ipynb`.

The target chemistry includes independent conserved components such as Fe(II),
Fe(III), Ca, Mg, phosphate, citrate, ammonium, and sulfate.  Aqueous complexes
and minerals can contain more than one of those components.

## Non-negotiable compatibility rule

The implementation **must be created as a standalone, separate code pathway**.
It must be tested and capable of running correctly before any discussion or
decision that could change the functionality of existing pathways.

In particular:

- Do not alter the behaviour of `BisectionChemicalEquilibriumEngine`, existing `NRTableau`, or
  existing `NRChemicalEquilibriumEngine` callers.
- Do not change existing `ReactionSystem(..., solver="newton_raphson")`
  semantics.
- Do not modify the current single-mineral precipitation path in place.
- New functionality must be opt-in, behind a new explicit API and covered by
  new tests.
- Only after the standalone pathway has passed its acceptance suite should a
  separate compatibility review consider sharing, replacing, or refactoring
  existing code.

## Why the current pathway cannot solve this network

`NRTableau` presently selects one master species for each connected reaction
graph.  This works for acid-base ladders and simple single-component examples.
It fails when reactions connect multiple independently conserved components.

For example, these reactions create one connected graph:

```text
Fe2+ + Cit3-    = FeCit-
Ca2+ + HPO4--   = CaHPO4(aq)
Fe3+ + HPO4--   = FeHPO4+
```

but the system requires separate balances for Fe(II), Fe(III), Ca, phosphate,
and citrate.  Selecting one graph master incorrectly reduces this to one mass
balance.

The current precipitation active-set code also returns a dissolved-equilibrium
result and mineral extent without persistent `SolidPhase` write-back.  Using it
inside a kinetic `ControlVolume` therefore risks apparent mass loss.

## New standalone pathway

Create a new package, for example:

```text
src/chemical_equilibrium/multicomponent/
    __init__.py
    components.py
    tableau.py
    solver.py
    precipitation.py
    cv_adapter.py
```

Suggested public entry point:

```python
from VLsim.chemical_equilibrium.multicomponent import MultiComponentEquilibriumEngine
```

The new engine must not be constructed by `ReactionSystem.engine`.  Notebook 08
and dedicated tests should construct it directly until its behaviour is proven.

## Phase 1: component model and aqueous tableau

### 1. Define independent components

Construct a component matrix rather than inferring one component from graph
connectivity.  At minimum, the engine needs a component vector for each
independent conserved quantity:

```text
Fe(II), Fe(III), Ca, Mg, Zn, Mn, Cu, Co,
phosphate, citrate, ammonium, sulfate, molybdate
```

Hydrogen is not a conserved-component row; it is the charge-balance variable.
Water is treated as a solvent with unit activity unless a future model opts out.

Each species carries component stoichiometry.  Examples:

| Species | Fe(II) | Fe(III) | Ca | phosphate | citrate |
|---|---:|---:|---:|---:|---:|
| `Fe2+` | 1 | 0 | 0 | 0 | 0 |
| `FeCit-` | 1 | 0 | 0 | 0 | 1 |
| `FeHPO4+` | 0 | 1 | 0 | 1 | 0 |
| `CaHPO4` | 0 | 0 | 1 | 1 | 0 |

### 2. Derive species activities from master activities

Use a general stoichiometric linear-algebra formulation:

```text
log(a_species) = log(K') + nu_H * log(a_H+) + sum(nu_i * log(a_component_i))
```

The engine must validate that the declared reaction basis is independent and
that every aqueous species is representable from the selected component basis.
Emit clear errors for rank deficiency or incompatible reaction declarations.

### 3. Solve residuals

Solve:

```text
one component mass balance per independent component
+ electroneutrality
```

Activity coefficients are an outer iteration initially, matching the existing
NR convention.  The first implementation may retain Davies activity correction,
but its validation report must state the medium ionic-strength range where that
choice is acceptable.

### 4. Avoid strong-ion double counting

A metal is either:

1. an inert strong ion, contributing fixed charge; or
2. an aqueous component, contributing via its free and complexed species.

It must never be both.  The new engine should make this a construction-time
validation error.

## Phase 2: coupled precipitation equilibrium

### 1. Mineral declarations

Use `EquilibriumReaction`-like declarations, but keep them local to the new
pathway initially.  The first mineral set for notebook 08 is:

```text
ferrihydrite: Fe(OH)3(s) + 3H+ = Fe3+ + 3H2O
strengite:    FePO4·2H2O(s) = Fe3+ + PO4--- + 2H2O
vivianite:    Fe3(PO4)2·8H2O(s) = 3Fe2+ + 2PO4--- + 8H2O
brushite:     CaHPO4·2H2O(s) = Ca2+ + H+ + PO4--- + 2H2O
struvite:     MgNH4PO4·6H2O(s) = Mg2+ + NH4+ + PO4--- + 6H2O
```

All constants must come from one named thermodynamic database/version.  The
notebook should record both the database and the reaction convention used for
each `log_K`.

### 2. Active-set algorithm

For active minerals, solve the coupled complementarity problem:

```text
xi_j >= 0
SI_j <= 0 for absent minerals
SI_j = 0 for present minerals
```

Use the full mineral-by-mineral Jacobian, not a diagonal approximation, because
ferrihydrite and strengite share Fe(III), while brushite and struvite share
phosphate.

Enforce non-negative dissolved component totals and use damped Newton steps or
a bounded least-squares step where necessary.

## Phase 3: persistent solid-phase CV adapter

Create `cv_adapter.py` to connect the standalone equilibrium engine to a
`ControlVolume` only after Phases 1 and 2 pass independently.

At each operator-split timestep:

1. Read conserved totals from `LiquidPhase` plus relevant `SolidPhase` mineral
   inventory.
2. Apply kinetic reaction sources to the liquid.
3. Run the multi-component aqueous-plus-mineral equilibrium solve.
4. Write dissolved species to `LiquidPhase`.
5. Write final mineral amounts, not merely incremental precipitation, to
   `SolidPhase.n_mol`.
6. Run the existing conservation monitor over both phases.

The adapter must support dissolution of existing solids when the liquid becomes
undersaturated.  A solid is therefore an inventory constraint, not a one-way
sink.

## Test plan and acceptance criteria

### Unit tests

- Citrate protonation mass balance and charge balance.
- Fe(II)-citrate and Fe(III)-citrate complexation mass balances.
- Ca/Mg phosphate complexation with independent component totals.
- Construction-time rejection of double-counted component/strong-ion species.
- Rank-deficient and unresolved reaction-network diagnostics.
- One-mineral ferrihydrite equilibrium.
- Two-mineral ferrihydrite/strengite competition.
- Brushite/struvite shared-phosphate competition.

### CV integration tests

- Precipitation conserves every element and charge across liquid plus solid.
- An initially present solid dissolves under undersaturation.
- A kinetic Fe(II)-oxidation source can create Fe(III), then precipitate or
  complex without a conservation warning.
- Existing single-phase, legacy speciation, and current NR test suites remain
  unchanged and passing.

### Notebook acceptance test

Notebook 08 must run from a fresh kernel and demonstrate:

1. initial equilibrium at the chosen pH;
2. aqueous complex distribution;
3. mineral saturation indices and solid amounts;
4. total Fe, P, Ca, Mg, citrate, and charge conservation;
5. a comparison against the current no-precipitation notebook.

## Integration decision gate

Do not change existing pathways until all of the following are true:

- the standalone engine has its own passing test suite;
- notebook 08 runs reproducibly from a clean kernel;
- conservation tests pass for precipitation and dissolution;
- numerical behaviour is compared against the selected external thermodynamic
  database for representative cases;
- maintainers explicitly approve a compatibility plan.

Only then should the project decide whether the new engine remains opt-in,
shares internal utilities with `NRChemicalEquilibriumEngine`, or becomes a future default.

## Note on convergence with the constraint-row framework (2026-07-01)

Per the `MASS_EXCHANGE_ARCHITECTURE.md` §14 design discussion: the reasons this
pathway must stay standalone for now (the master-species-per-connected-component
limitation in `NRTableau`) are unrelated to the compatibility-rule caveats
above and remain fully in force — nothing here changes the acceptance-gate
sequencing. But when this engine does reach its integration decision, the
recommended target is implementing the same constraint-row protocol used
elsewhere in the equilibrium-engine framework (see `MASS_EXCHANGE_ARCHITECTURE.md`
§14 and `CHEMICAL_EQUILIBRIUM_ENGINE_ARCHITECTURE.md`), so this becomes a
pluggable engine backend rather than a fourth independent precipitation
representation. This is a preference for that future decision point, not a
change to the standalone-first rule above.
