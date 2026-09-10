# Chemistry Unification — Rationale and Direction

> **✅ All chemistry-unification phases complete (2026-06-03).** All
> six branches (1 through 4, plus 3b) shipped.
> Tags: `chemistry-unification-1-shipped` through
> `chemistry-unification-4-shipped` and
> `chemistry-unification-3b-shipped`.

> **⚠ Reconciliation banner (2026-06-01).** Several design
> decisions in this document were resolved differently by
> subsequent phases. Read the body as the original rationale;
> use the notes below to understand what actually shipped.
>
> - **Prerequisite note (line 16–18):** References `chem_env`
>   and `property_results` naming conventions from CV_UPDATE.md.
>   Both were **deleted** by STATE_UNIFICATION (2026-05-22).
>   Storage is now `phase.n_mol` (unified); `PropertyResult` /
>   `PropertySolver` / `chem_env` no longer exist.
>
> - **Staleness machinery (§Equilibrium State / B2 resolution):**
>   The proposed `_equilibrium_stale` flag, `StaleEquilibriumError`,
>   and three-state pH distinction were **rejected**. STATE_UNIFICATION
>   chose stateless snapshot semantics instead: `Phase.pH` is a
>   plain `@property` from `n_mol["H+"]`; there are no staleness
>   guards.
>
> - **`ThermoFramework` object:** Described as a future
>   activity-model configuration surface. **Never implemented.**
>   Activity-model configuration remains at `SpeciationEngine`
>   construction time.
>
> - **Phase A–F scope (§Scope and Phasing):** The actual
>   shipping sequence diverged. Phases A–B merged into
>   CHEMISTRY_UNIFICATION Phase 1 (shipped 2026-05-13). The
>   authoritative execution record is
>   [CHEMISTRY_UNIFICATION_PLAN.md](CHEMISTRY_UNIFICATION_PLAN.md).

## Summary

Speciation (pH, acid-base equilibria) is physically a subset of reactions, but the
codebase treats it as a parallel subsystem with its own API, its own hardcoded
configuration, and its own conceptual framing. This document proposes unifying
the **declaration surface** — making equilibrium reactions a mode of the existing
`Reaction` framework — while **keeping the existing algebraic solver intact**.

This is not a solver rewrite. The current Newton-Raphson-in-log[H⁺] speciation
solver is numerically the same approach PHREEQC uses and is kept unchanged. What
changes is how its inputs are gathered: from declared equilibrium reactions rather
than from hardcoded Level 1/2/2.5 logic and externally-supplied totals.

**Prerequisite:** [CV_UPDATE.md](../shipped/CV_UPDATE.md) must land first. That refactor
fixes the `advance()` ordering, eliminates `_last_properties`, and establishes
the `chem_env` / `property_results` naming conventions this document assumes.

---

## Current State

The codebase has two parallel chemistry subsystems:

- **`src/reactions/`** — kinetic reactions. `Reaction` objects carry stoichiometry
  with mass-balance validation at construction. Integrated via RK4 / scipy ODE
  solvers.

- **`src/speciation/`** and **`src/equilibria/`** — acid-base equilibria. Hardcoded
  Level 1, 2, and 2.5 engines. Takes externally-supplied `acid_totals`, `CT_P`,
  `CT_NH_T`, `strong_kwargs` as input. Solves algebraically via Newton-Raphson
  on log[H⁺] with a charge-balance residual.

- **`KineticGasLiquidLink.SpeciationCorrection`** — a mini-speciation engine
  embedded in the gas-liquid transfer code, reimplementing pKa → molecular
  fraction logic inline.

---

## Why the Split Exists

**Legitimate reason (solver-level):** Speciation reactions operate on timescales
10+ orders of magnitude faster than kinetic bioreactions (microseconds vs hours).
Integrating both as ODEs in a unified system produces a stiff problem that is
slow and numerically fragile. The pragmatic solution — treat the fast equilibria
as instantaneous algebraic constraints — is sound engineering.

**Illegitimate reason (conceptual):** The solver-level split got elevated into a
declaration-level split. The codebase has `src/reactions/` and `src/speciation/`
as sibling subsystems with different APIs, different validation, different
construction patterns, and different contribution mechanisms to the CV state.
Users and developers have to learn two frameworks for what is physically one
concept.

---

## Architectural Costs of the Current Split

1. **Chemical system declared twice.** Acid-base pairs live as separate species
   in `phase.n_mol` *and* as totals in `chem_env["acid_totals"]`. The orchestrator
   is responsible for keeping these in sync.

2. **`chem_env` exists largely because of the split.** The reason speciation
   needs external parameters it can't derive from phase moles is that the
   relationships between species (pKas, acid-base groupings) aren't declared
   anywhere the CV can find them.

3. **Speciation logic leaks into the gas-liquid link.** `SpeciationCorrection`
   reimplements `f_molecular(pH, T, pKas)` inline because the speciation API
   only exposes "solve the whole system," not "give me this one ratio."

4. **Elemental validation doesn't extend to speciation.** `Reaction.__init__`
   validates C, H, O, N balance. Acid-base equilibria aren't expressed in that
   framework, so no equivalent check exists. A pKa entry for a reaction that
   doesn't actually balance has no safety net.

5. **Two parallel construction APIs for one concept.** Adding a new weak acid
   with volatile and bioreactive participation requires touching the chemistry
   registry, speciation engine configuration, `speciation_corrections` in the
   transfer link, `Reaction` definitions, and the orchestrator's `acid_totals`
   population — five surfaces for one physical concept.

6. **Operator splitting at the declaration level.** Reactions and speciation
   are composed via last-step pH handoff rather than co-declaration. This is
   numerically OK for small `dt` but it reflects the architectural split back
   into the integration.

---

## Proposed Approach: Unified Declaration, Same Solver

A `Reaction` gains a `mode` field with values `"kinetic"` (default) and
`"equilibrium"`. Kinetic reactions carry a rate function as they do today.
Equilibrium reactions carry a `log_K` (or equivalent) value instead.

```python
# Kinetic — unchanged
Reaction(
    stoichiometry=[("Glucose", -1), ("Ethanol", +2), ("CO2", +2)],
    rate_fn=lambda env: ...,
)

# Equilibrium — new
Reaction(
    mode="equilibrium",
    stoichiometry=[("AceticAcid", -1), ("Acetate-", +1), ("H+", +1)],
    log_K=-4.76,
)
```

A solver dispatcher partitions reactions by mode:

- `mode="equilibrium"` reactions feed the existing algebraic speciation solver.
  Mass balance constraints are derived from the union of equilibrium stoichiometries;
  the charge balance is derived from declared species charges; pKas come from
  `log_K` fields.

- `mode="kinetic"` reactions feed the existing ODE integrator as today.

The two sets are composed via operator splitting (as today) or optionally via
DAE integration (future, not required).

---

## What's Borrowed from PHREEQC

PHREEQC is the reference implementation for declarative aqueous equilibrium.
Numerically it uses the same solver family as VLsim today. Architecturally it
makes the following choices worth borrowing:

**1. Reactions as the declaration surface.** A single `Reaction` type (with a
mode field) is the one way to declare chemistry. Kinetic and equilibrium
reactions live in the same framework, share stoichiometry validation, and share
the `ReactionEnvironment` consumed by rate / equilibrium evaluation.

**2. Conservation derived from stoichiometry.** Totals like `acid_totals["AceticAcid"]`
or `CT_P` are computed by summing across species linked by declared equilibrium
reactions, not passed in externally. `chem_env`'s contents shrink dramatically —
likely reduced to `t_h` and solver hints.

**3. Same algebraic solver, fed declaratively.** The existing Newton-Raphson
speciation solver is kept. What changes is that its matrix (species, charges,
equilibrium constants, conservation rows) is constructed from the declared
equilibrium reactions at solver-setup time rather than from hardcoded Level 1/2/2.5
configurations.

**5. Declarative chemistry database.** Thermodynamic parameters and species
declarations move from hardcoded engine logic into a `ChemistryDatabase` — a
runtime object bundling `ThermoFramework`, species declarations, and the
`ReactionSet`. Default databases ship as Python modules (see *Database
Packaging* below). Users can swap or extend databases the way PHREEQC users
swap `phreeqc.dat` / `minteq.dat` / `pitzer.dat`, but through Python
composition rather than a parsed text format.

---

## What's NOT Borrowed (and why)

**4. Element-level state tracking.** PHREEQC tracks moles of elements rather
than moles of species, and computes species concentrations by distributing
element totals across all species containing each element.

This provides real benefits for PHREEQC's use case — dynamic species creation
via mineral dissolution, surface complexation, redox distribution across
undeclared species — but none of those scenarios apply to VLsim. Here all
species are pre-declared, all transformations are explicit reactions, and all
phases are discrete.

Under the unified-declaration proposal above, element conservation is already
guaranteed at the declaration level: every state change passes through a
`Reaction` with mass-balance validation, or through an external flux with
declared composition. Adding element-level tracking on top would:

- Duplicate an already-enforced invariant
- Introduce a translation layer (element moles ↔ species moles on every step)
- Expand the API surface for no new functionality

For VLsim, element awareness is architectural weight without a corresponding
feature.

---

## Runtime Packaging: `ThermoFramework`

The PHREEQC-style data-file-driven approach needs a runtime object to represent
the thermodynamic conventions carried from the database into a simulation.
`ThermoFramework` is that object. It holds the activity model, reference state,
and standard-state parameters for a physical phase, and is the single source of
truth for "how thermodynamic quantities are interpreted" within a CV.

### Naming

`ThermoContext` was considered and rejected — it collides with the
[CV_UPDATE.md](../shipped/CV_UPDATE.md) decision to remove `context` as a name. The
industry-standard term `PropertyPackage` (Aspen, HYSYS) was also considered but
rejected because it visually collides with VLsim's existing `PropertySolver`
(a dynamic state-derivation concern, not a thermodynamic-convention concern).

**`ThermoFramework`** — unambiguous, no collisions, reads naturally at call
sites.

### Contents

```python
@dataclass(frozen=True)
class ThermoFramework:
    activity_model: ActivityModel
    reference_state: str = "infinite_dilution"
    standard_T_K: float = 298.15
    standard_P_atm: float = 1.0
    # fugacity_model: Optional[FugacityModel] = None   # later
```

Frozen by design — a thermodynamic framework shouldn't change mid-simulation.
Runtime overrides use `dataclasses.replace()`:

```python
thermo_ideal = dataclasses.replace(thermo_default,
                                    activity_model=IdealActivityModel())
```

### Construction paths

```python
# From a shipped database module (the default path)
from VLsim.chemistry.databases.aqueous import AQUEOUS_DEFAULT
thermo = AQUEOUS_DEFAULT.thermo         # packaged alongside reactions

# From scratch (tests, debugging)
thermo = ThermoFramework(
    activity_model=DaviesActivityModel(),
    reference_state="infinite_dilution",
)

# Shared across CVs by reference
cv_zone1 = ControlVolume(phases=..., thermo=thermo)
cv_zone2 = ControlVolume(phases=..., thermo=thermo)
```

### Access path through the object graph

```
ThermoFramework  (held by CV)
     │
     ▼
ReactionEnvironment  (built per step, forwards thermo)
     │
     ▼
rate_fn(env)         env.c(species)  → concentration
speciation solver    env.a(species)  → activity (uses thermo)
```

`ReactionEnvironment` gains a `thermo` field and an `a(species)` method.
Existing rate functions that use `env.c(...)` or `env.concentrations[...]`
continue to work unchanged; those that want activity-corrected kinetics opt in
via `env.a(...)`. The `SpeciationPropertySolver` drops its own `activity_model`
parameter and reads `env.thermo.activity_model` at solve time.

### Effect on `ControlVolume`

The CV gains one field (`thermo`) and zero new methods. Its responsibilities
are unchanged — state ownership, advance orchestration, delegation to solvers
and reactions. Thermodynamic evaluation stays where it belongs: in solvers,
in `ReactionEnvironment`'s accessors, and in the activity model classes.

This parallels how the CV already holds `property_solvers` and `reaction_model`
as external configuration without taking on their computational responsibilities.

### Per-CV thermo and physical consistency

Allowing each CV to hold its own `ThermoFramework` is intentional, but comes
with a consistency concern depending on what the CVs represent:

- **Case A — CVs are zones of the same physical phase** (e.g. sparger + bulk
  zones of one aqueous broth). All CVs *must* share one `ThermoFramework`. The
  solution is the same solution; different conventions would break charge
  balance across zone boundaries and make inter-zone advective transfer
  thermodynamically inconsistent.

- **Case B — CVs are physically distinct phases** (e.g. aqueous broth +
  organic extractant; liquor + solid in a crystallizer). Per-CV frameworks are
  *required*. Aqueous electrolyte models (Davies, Pitzer) and organic activity
  models (NRTL, UNIQUAC) are different thermodynamic worlds; enforcing shared
  thermo would be physically wrong.

The architecture cannot tell Case A from Case B automatically. VLsim's current
use case is Case A; Case B becomes relevant as the framework extends to
multi-phase bioprocesses (ATPS, in-situ product removal, crystallization).

### Consistency enforcement

**Now (Option 3) — warn on link construction.** When a `MultiCVSystem` or an
`AdvectiveLink` / `DiffusiveLink` connects CVs with different
`ThermoFramework` instances, emit a warning:

```
UserWarning: AdvectiveLink connects CVs with different ThermoFrameworks.
This is unusual — expected if the CVs represent physically distinct phases
(e.g., aqueous + organic); potentially a bug if they represent zones of the
same phase.
```

Cheap to implement, catches the realistic mistake without foreclosing valid
multi-phase use.

**Later (Option 4) — link types declare their requirement.** When a Case-B link
type is added (e.g. a hypothetical `ExtractionLink`), link types declare
whether thermo-matching is required:

```python
class AdvectiveLink:
    requires_shared_thermo = True       # same solution recirculating

class ExtractionLink:                    # hypothetical future
    requires_shared_thermo = False      # immiscible phases by design
```

Option 4 is the principled long-term answer but premature while VLsim has only
one link family. Option 3 is sufficient for now.

---

## Database Packaging: `ChemistryDatabase`

The `ChemistryDatabase` is the runtime bundle that holds a `ThermoFramework`,
a species declaration dict, and a `ReactionSet`. It is the "chemistry library"
that a CV is configured with.

### Format: Python modules, not YAML

Shipped databases are plain Python modules. The Python object graph *is* the
database; there is no separate text format to parse, validate, or keep in sync
with the object model.

```python
# VLsim/chemistry/databases/aqueous.py
from VLsim.chemistry import ChemistryDatabase, Species
from VLsim.reactions import Reaction, ReactionSet
from VLsim.thermo import ThermoFramework, DaviesActivityModel

THERMO = ThermoFramework(activity_model=DaviesActivityModel())

SPECIES = {
    "H+":         Species(charge=+1, formula={"H": 1}),
    "OH-":        Species(charge=-1, formula={"H": 1, "O": 1}),
    "AceticAcid": Species(charge= 0, formula={"C": 2, "H": 4, "O": 2}),
    "Acetate-":   Species(charge=-1, formula={"C": 2, "H": 3, "O": 2}),
    # ...
}

REACTIONS = ReactionSet([
    Reaction(mode="equilibrium",
             stoichiometry={"H2O": -1, "H+": +1, "OH-": +1},
             log_K=-14.0),
    Reaction(mode="equilibrium",
             stoichiometry={"AceticAcid": -1, "Acetate-": +1, "H+": +1},
             log_K=-4.76),
    # ...
])

AQUEOUS_DEFAULT = ChemistryDatabase(thermo=THERMO,
                                     species=SPECIES,
                                     reactions=REACTIONS)
```

Usage at the CV level:

```python
from VLsim.chemistry.databases.aqueous import AQUEOUS_DEFAULT

cv = ControlVolume(
    phases={"liquid": liquid},
    chemistry_db=AQUEOUS_DEFAULT,
    property_solvers=[SpeciationPropertySolver()],
)
```

Extension by composition (the mechanism for user customisation):

```python
from VLsim.chemistry.databases.aqueous import AQUEOUS_DEFAULT

MY_CHEM = AQUEOUS_DEFAULT.extend(
    species={"ButyricAcid": Species(charge=0, formula={"C": 4, "H": 8, "O": 2}),
             "Butyrate-":   Species(charge=-1, formula={"C": 4, "H": 7, "O": 2})},
    reactions=[Reaction(mode="equilibrium",
                         stoichiometry={"ButyricAcid": -1, "Butyrate-": +1, "H+": +1},
                         log_K=-4.82)],
)
```

### Why Python modules, not YAML

The choice of Python over YAML is deliberate:

- **VLsim's users already write Python.** The "YAML makes databases accessible
  to non-programmers" argument doesn't apply — anyone customising VLsim
  chemistry is already in a Python environment.
- **Type safety at authoring time.** IDE autocomplete, mypy, and dataclass
  `__post_init__` catch typos (wrong species name, malformed formula, missing
  charge) before the simulation runs. YAML + Pydantic gives the same at load
  time, but later and with more indirection.
- **No parser to maintain.** YAML needs PyYAML + a schema (Pydantic or
  equivalent) + migration handling across schema versions. Python modules
  have none of these costs.
- **Natural composition.** `AQUEOUS_DEFAULT.extend(...)` is more powerful and
  more legible than YAML merge strategies.
- **No new dependency.** PyYAML and Pydantic are both avoided.
- **No data/code boundary to police.** Shipped databases are trusted package
  code; user-authored databases are authored by the same person who runs the
  simulation. The "arbitrary code execution" concern of loading Python as
  data doesn't apply.
- **Consistency with the rest of VLsim.** Everything else in the codebase is
  Python objects composed together; chemistry databases are no exception.

### When YAML (or another format) would be reconsidered

YAML could be added as a **conversion utility** (producing a
`ChemistryDatabase` from a file) if a concrete need arises — e.g. a user
community with a large corpus of existing text databases, a non-Python
frontend, or a requirement to ship databases as user-editable config. None of
those apply today.

The primary representation remains the Python object graph either way. A
hypothetical future `ChemistryDatabase.from_yaml(path)` would be an additive
loader, not a replacement.

### `ChemistryDatabase` structure

```python
@dataclass(frozen=True)
class ChemistryDatabase:
    thermo: ThermoFramework
    species: Dict[str, Species]
    reactions: ReactionSet

    def extend(self,
               species: Optional[Dict[str, Species]] = None,
               reactions: Optional[Iterable[Reaction]] = None,
               thermo: Optional[ThermoFramework] = None,
               ) -> "ChemistryDatabase":
        ...
```

Frozen dataclass — composed by extension, not mutated in place. No loader, no
parser, no schema file. The Python interpreter is the loader; type hints and
dataclass validation are the schema.

### Shipped layout

```
VLsim/
  chemistry/
    databases/
      __init__.py
      aqueous.py              # AQUEOUS_DEFAULT
      bioprocess_basic.py     # BIOPROCESS_BASIC (extends AQUEOUS_DEFAULT)
      anaerobic_digestion.py  # AD_BASIC (extends BIOPROCESS_BASIC)
```

Each database is a module exposing a single top-level constant
(`AQUEOUS_DEFAULT`, `BIOPROCESS_BASIC`, ...). Extensions compose via `.extend()`.

---

## Leakage That Won't Go Away Automatically

Even after the unified-declaration refactor, two specific leaks in the gas-liquid
link persist unless explicitly addressed.

### Leak 1: `SpeciationCorrection.f_molecular()` inside `KineticGasLiquidLink`

**Resolved 2026-05-15 in `chemistry-unification-2`.**
`PropertyResult` now carries an `alphas: Dict[str, float]` channel
populated by `SpeciationEngine.solve()`, keyed by molecular form
(`"CO2aq"`, `"NH3"`, `"{Acid}_HA"`). `KineticGasLiquidLink._effective_henry`
reads `alphas[mol_key]` via the `speciation_keys` mapping — collapsing
the two pre-existing correction paths (inline pKa/pH math via
`SpeciationCorrection.f_molecular`, and dict lookup in `species` for
the concentration-based path) into a single alpha read. The link and
engine now use the same effective-Ka math; the BSM2 trajectory drifted
by 1e-7 to 1e-5 across species (Leak 1 closing — see
`tests/standalone/test_bsm2_reference.py` rebaseline comment).
`SpeciationCorrection` and its plumbing
(`set_speciation_correction`, builder's `speciation_correction_pka`,
`TransferConfig.speciation_corrections`,
`ThermodynamicConfig._apply_to_transfer_link`) deleted entirely.

### Leak 2: `speciation_keys` mapping in `KineticGasLiquidLink.__init__`

**Infrastructure shipped 2026-05-15 in `chemistry-unification-3`;
builder migration deferred to Phase 3b.**

`Reaction` gained an `is_cross_phase` property and accepts `log_K=None`
for cross-phase equilibrium reactions (partition declarations
consumed by gas-liquid links, not by the speciation engine).
`SpeciationEngine.from_reactions` filters cross-phase reactions out
silently. `KineticGasLiquidLink.derive_speciation_keys(reaction_set)`
iterates declared cross-phase equilibrium reactions and records the
gas → liquid mapping. `ControlVolume.__init__` calls
`derive_speciation_keys` on every internal link after the reaction
model is attached. The cross-phase reaction declaration shape is:

```python
Reaction(
    kind="equilibrium",
    stoichiometry=[
        StoichiometryEntry(species=CO2, phase="gas",    coefficient=-1.0),
        StoichiometryEntry(species=CO2, phase="liquid", coefficient=+1.0),
    ],
    # log_K=None: Henry's law lives on the link, not on the reaction
)
```

**What's still deferred** (Phase 3b — unified ChemistryDatabase
rollout): BSM2/ADM1 builder migration to declare these reactions
and drop their explicit `speciation_keys=` constructor kwargs.
The deferral exists because today's `PropertyResult.alphas` keys
(`"CO2aq"`, `"S_ac_HA"`) are *synthesized* by `_compute_species_eq`
rather than matching the Species IDs declared in reaction
stoichiometry (`"CO2"`, `"S_ac"`). The clean fix is to make the
engine emit species using their declared Species IDs directly —
which requires `EquilibriumDef` to carry Species references,
touches the legacy `bsm2_default()` factory path, and migrates
every test asserting on `out["CO2aq"]` / `out["S_ac_HA"]`. Scope
is closer to Phase 3b than to the cross-phase wiring itself, so
the migration ships alongside the unified Species-ID convention.

Phase 3 keeps the explicit `speciation_keys=` kwargs working
through additive semantics: `derive_speciation_keys` leaves
pre-existing entries alone unless overridden by a declared
cross-phase reaction. BSM2/ADM1 continue to function unchanged.

### Leak 3: dormant f_molecular cross-check in `thermo_params._check_f_molecular_consistency`

**Resolved 2026-05-15 in `chemistry-unification-2` — deleted, not
revived.** With `SpeciationCorrection` gone (Leak 1 closure), the
parameter-drift bug class this check guarded against (link's pKas
vs engine's pKas) can no longer occur: the link no longer holds
pKas. The check's purpose evaporated alongside its dormancy.
`_check_f_molecular_consistency` and its call site in
`validate_thermodynamics` are both removed. The remaining
within-zone / across-zone checks in `validate_thermodynamics` still
catch parameter mismatches between `_thermo_config` and
`chem_env_fn` snapshots at declaration time; Phase 3
(`ChemistryDatabase`) will rewrite the data sourcing.

---

## Class-Diagram Impact

[docs/class_diagrams.md](../class_diagrams.md) currently frames
layer 3 as "two parallel chemistry-flavoured plug-ins on the CV"
(`PropertySolver` + `ReactionModel`).  That framing reflects the
post-Phase-7 state and **will need a rewrite, not just an addition**,
when this phase ships.

Concretely:

- Declaration-level parallelism collapses.  Both kinetic and
  equilibrium chemistry are expressed as `Reaction` objects
  distinguished by a `mode` field; one `ReactionSet` holds both
  flavours.  Layer 3's intro paragraph and the parallel
  `PropertySolver ←→ ReactionModel` framing both go.
- Solver-level parallelism stays, internally.  The diagram should
  surface this as a dispatcher inside the CV that routes
  `mode="equilibrium"` reactions to the algebraic speciation solver
  and `mode="kinetic"` reactions to the ODE integrator.
- `SpeciationPropertySolver` and `SpeciationCorrection` mostly
  disappear (absorbed into the unified `Reaction` framework and the
  new `PropertyResult.alphas` channel).  `PropertySolver` survives
  for genuinely-property-flavoured concerns (viscosity, density
  correlations) that don't fit the stoichiometry model.
- The `KineticGasLiquidLink ..> PropertyResult : reads pH` arrow
  becomes a phase-state read (`phase.alphas` or equivalent) — see
  Leak 1 above.
- `chem_env` shrinks to a near-empty dict (likely just `t_h` and
  solver hints).  Any class with a `chem_env: dict` field can have
  it dropped or renamed.

Update the layer-3 intro paragraph and the `<<protocol>>` /
relationship arrows together — the rewrite is more than a search-
and-replace, and bundling it with the implementation PR keeps the
diagrams and the code in lockstep.

---

## Redox: Deferred, With Forward-Compatibility Constraints

Redox chemistry is out of scope for this refactor but **not foreclosed** by it.
The refactor commits to a small set of design constraints that leave redox as
a clean future extension rather than a breaking change.

### Why defer

- No current VLsim use case requires redox — designing in the abstract risks
  an API that doesn't fit real use cases when they arrive
- The chemistry unification is already substantial; each added scope item
  compounds implementation risk
- Redox has its own concerns (reference electrode convention, pe often badly
  constrained in biological systems, non-physical electron activity) that
  deserve their own design discussion when there's a concrete use case

### Forward-compatibility constraints

The refactor must honour these to keep redox available as an additive
extension:

1. **Open species namespace.** No enumeration of allowed species anywhere in
   the registry, reaction framework, or speciation solver. Declaring `"e-"`
   (or any novel species) must register without complaint.

2. **Iterative charge balance.** The speciation solver's charge-balance
   residual iterates over declared species — `∑(z_i × C_i)` — rather than
   summing a hardcoded set of terms. Extending to electrons or new redox
   species is automatic.

3. **Element loops skip formula-less species.** `Reaction.__init__` validates
   elemental balance only for species that declare a formula. Species with
   no formula (e.g. `e-`) are skipped, not rejected. No canonical element set
   (C/H/O/N) is enforced.

4. **Master variables are enumerable.** The speciation solver's Newton
   iteration is structured as a list of master variables — even if that list
   contains a single entry (`log[H+]`) today. Adding `log[e-]` later is a
   list-item addition, not a re-architect.

5. **`PropertyResult.raw` is the extension hook.** Optional output fields
   (`pe`, `E_h`, redox couples) are added to `raw` without schema changes.

6. **`ThermoFramework` accepts optional fields.** A future
   `reference_electrode: str = "SHE"` field can be added via
   `dataclasses.replace()` without breaking existing construction paths.

7. **`log_K` is general.** On an equilibrium reaction, `log_K` represents
   whatever equilibrium relation the stoichiometry expresses. E° → log K
   conversion utilities are additive, not schema-changing.

### What is explicitly not implemented

- No `e-` species in any default data file
- No `pe` or `E_h` computation in the speciation solver output
- No charge-balance extension for electron balance
- No E° / electrode potential utilities or conversions
- No redox-specific `ThermoFramework` fields
- No validation rules biased toward or against redox reactions

These are all additive extensions when a concrete use case arrives.

---

## Equilibrium State: Storage, Staleness, and Warm-Starting

Under the unified chemistry approach, equilibrium-species concentrations
(H⁺, OH⁻, HCO₃⁻, CO₃²⁻, ...) are determined algebraically by the speciation
solver, not integrated kinetically. Once computed they need to be stored
somewhere for downstream consumers (rate functions, gas-liquid transfer,
controllers) to read, and warm-starting the next Newton-Raphson solve needs a
home.

### Storage: in the phase's `n_mol`

The speciation solver writes computed equilibrium-species concentrations
directly into the liquid phase's `n_mol` dict, alongside kinetic species. Two
consequences:

- Warm-start for the next solve is implicit: `[H⁺]` from the previous solve is
  still there, ready to seed the Newton-Raphson initial guess.
- Rate functions reading `env.c("HCO3-")` get the equilibrium-correct value —
  provided they run after speciation.

No separate hint dict, no hint field on `PropertyResult`. `chem_env` sheds
`logH_guess` entirely and shrinks further toward empty.

### Why in-phase storage needs guarding

Equilibrium species in `n_mol` have a subtly different contract than kinetic
species:

- **Kinetic** (`Glucose`, `Ethanol`, ...): always meaningful. Integrated
  through time. `phase.n_mol["Glucose"]` is a physical quantity at any moment.
- **Equilibrium** (`H⁺`, `HCO₃⁻`, ...): meaningful *only immediately after* a
  speciation solve. Between solves, after mutations to kinetic species, the
  stored values no longer satisfy charge balance — they are *stale*.

Unguarded, this would let a user compute `-log10(phase.n_mol["H+"] / V_L)` and
get a stale pH silently. The following guardrails address the three distinct
failure modes (undefined, stale, and safe access) without over-constraining
reactions-only simulations.

### Three-state distinction

The phase carries two pieces of state, not one boolean:

```python
class LiquidPhase:
    n_mol: Dict[str, float]
    V_L: float
    T_K: float
    _equilibrium_species: FrozenSet[str] = frozenset()   # what's algebraic
    _equilibrium_stale: bool = True                      # meaningful iff above is non-empty
```

| State | `_equilibrium_species` | `_equilibrium_stale` | Meaning |
|---|---|---|---|
| **Undefined** | `frozenset()` | (inert) | No equilibrium chemistry configured |
| **Stale** | non-empty | `True` | Equilibrium configured but phase mutated since last solve |
| **Fresh** | non-empty | `False` | Equilibrium species consistent with composition |

`_equilibrium_species` is populated at CV construction, derived from the
declared equilibrium reactions in the attached `ChemistryDatabase`. A
simulation with only kinetic reactions gets an empty set, and the staleness
machinery stays dormant.

### Accessors

**`phase.pH` (property)** — returns current pH, or raises:

```python
@property
def pH(self) -> float:
    if not self._equilibrium_species:
        raise PHUndefinedError(
            "pH is not defined for this phase: no equilibrium chemistry has "
            "been declared. To enable pH, declare acid-base reactions "
            "(including H2O ⇌ H+ + OH-) in the chemistry database and attach "
            "a SpeciationPropertySolver to the CV."
        )
    if self._equilibrium_stale:
        raise StaleEquilibriumError(
            "Phase has been mutated since the last speciation solve. Call "
            "cv.current_pH() to trigger a fresh solve, or wait until after "
            "the next advance()."
        )
    return -math.log10(self.n_mol["H+"] / self.V_L)
```

**`cv.current_pH()`** — runs speciation on demand if stale, returns pH:

```python
def current_pH(self, chem_env=None) -> float:
    liquid = self.phases["liquid"]
    if not liquid._equilibrium_species:
        raise PHUndefinedError(...)
    if liquid._equilibrium_stale:
        self._resolve_equilibrium(chem_env)
    return liquid.pH
```

**`env.has_pH` (property)** — boolean for rate functions that optionally use
pH:

```python
@property
def has_pH(self) -> bool:
    return bool(self._liquid_phase._equilibrium_species)
```

**`phase.n_mol["H+"]` (direct dict access)** — continues to work; returns
whatever is stored with no staleness check. Respects the dict contract; this
is the correct access pattern for speciation solvers warm-starting from the
previous value. Users reading this directly accept responsibility for
interpreting it correctly.

### Mutation discipline

`LiquidPhase.apply_flux()` (and any other mutation path) sets the stale flag,
but only if equilibrium is active:

```python
def apply_flux(self, flux, dt_h):
    # ... apply the flux
    if self._equilibrium_species:
        self._equilibrium_stale = True
```

Speciation solvers clear the flag after writing their results back:

```python
class SpeciationPropertySolver:
    def solve(self, env) -> PropertyResult:
        # ... Newton-Raphson on log[H+]
        phase = env.phases["liquid"]
        phase.n_mol.update(equilibrium_concentrations)
        phase._equilibrium_stale = False
        return PropertyResult(pH=pH, ...)
```

### Error types

- **`PHUndefinedError`** — phase has no equilibrium chemistry configured.
  Raised by `phase.pH`, `env.equilibrium_c(...)`, and `cv.current_pH()`.
- **`StaleEquilibriumError`** — equilibrium is configured but the phase has
  been mutated since last solve. Raised by `phase.pH` and
  `env.equilibrium_c(...)`. Not raised by `cv.current_pH()`, which resolves
  by running the solver.

Both are correctness-level errors, not warnings — silent stale-pH reads are
exactly what this system is designed to prevent.

### Construction-time warning

At `ControlVolume.__init__`, warn if initial `n_mol` contains equilibrium
species:

```python
if liquid._equilibrium_species and any(
    sp in liquid.n_mol for sp in liquid._equilibrium_species
):
    warnings.warn(
        "Liquid phase contains initial values for equilibrium-determined "
        f"species ({sorted(sp for sp in liquid._equilibrium_species if sp in liquid.n_mol)}) "
        "that will be overwritten on the first advance() call. Specify "
        "initial state in terms of species totals (e.g. AceticAcid, not "
        "Acetate-) — the solver will distribute them.",
        UserWarning,
    )
```

Fires only when equilibrium is configured and the user has written to species
the solver will overwrite. Reactions-only simulations never trigger this
warning.

### Speciation solver contract

Custom speciation solvers must meet a contract, documented on the
`PropertySolver` protocol:

> If the solver computes equilibrium species, it must:
> 1. Write all computed equilibrium concentrations into `phase.n_mol` before
>    returning.
> 2. Set `phase._equilibrium_stale = False` before returning.
>
> Failure to meet either condition leaves the phase in an inconsistent state
> that will cause downstream reads to raise `StaleEquilibriumError` or produce
> incorrect warm-starts.

A development-only test utility validates custom solvers meet the contract:

```python
def validate_speciation_solver(solver, test_phase):
    test_phase._equilibrium_stale = True
    solver.solve(env_from(test_phase))
    assert not test_phase._equilibrium_stale, (
        f"{type(solver).__name__} failed to clear the stale flag."
    )
    for species in test_phase._equilibrium_species:
        assert species in test_phase.n_mol, (
            f"{type(solver).__name__} did not write back {species}."
        )
```

### Behaviour for reactions-only simulations

A simulation with only kinetic reactions (no equilibrium reactions declared,
no speciation solver attached) behaves as follows:

- `_equilibrium_species` is empty at construction
- `_equilibrium_stale` is never checked (dormant)
- `apply_flux()` skips the stale-setting branch — no overhead
- `phase.n_mol["H+"]` is an ordinary species value (zero unless declared) —
  no warnings, no errors
- `phase.pH` raises `PHUndefinedError` with a clear message if the user tries
  to read it — no silent `-log10(0)` or `nan`
- `env.has_pH` returns `False`, so rate functions can branch gracefully
- No construction warnings
- No performance penalty

The equilibrium-state machinery scales down cleanly: inactive when not
configured, active when needed.

---

## Equilibrium-Kinetic Coupling and Accuracy Monitoring

### Coupling strategies: the spectrum

Three strategies span the space from loosest to tightest coupling:

| Strategy | Description | VLsim today |
|---|---|---|
| **Operator splitting** | Equilibrium solved once per step; kinetics integrated with frozen pH | `EulerSnapshotSolver`, and `ScipyODESolver(freeze_speciation=True)` |
| **Pseudo-DAE** | Equilibrium re-solved at each ODE derivative evaluation `f(t, y)` | `ScipyODESolver` (default) |
| **True DAE** | Mass-matrix formulation; algebraic + differential satisfied simultaneously | Not implemented |

Both ends of the current implementation spectrum are already in place. The
question is whether to add true DAE.

### Decision: defer true DAE, architect to permit

True DAE (via CasADi, SUNDIALS IDA, or equivalent) offers rigorous
simultaneous satisfaction of differential and algebraic equations and handles
stiff systems well. It also adds a heavy C-library dependency, requires
mass-matrix setup, and is overkill for VLsim's bioprocess domain where
pseudo-DAE with `DOP853` or stiffness-tolerant `Radau`/`BDF` handles nearly
everything.

**Defer true DAE.** Keep `EulerSnapshotSolver` as the fast path and
`ScipyODESolver` as the accurate path.

**Architect the refactor so true DAE can be added later as an additive
extension.** Four constraints:

1. **`Reaction.mode` declares what a reaction is, not how to solve it.**
   A future DAE solver consumes the same `Reaction(mode="equilibrium", ...)`
   declarations, just treats them as algebraic constraints rather than running
   them through separate Newton-Raphson.

2. **State partitioning is derived, not hardcoded.** The split between
   kinetic (integrated) and equilibrium (algebraic) species is already
   captured by `LiquidPhase._equilibrium_species` from the warm-start design.
   A future DAE solver reads that set to build its mass matrix rows.

3. **Solver selection is per-CV.** A CV instantiates one solver. A future
   `DAESolver` plugs in at the same seam as `EulerSnapshotSolver` /
   `ScipyODESolver`.

4. **No implicit operator-splitting semantics leak elsewhere.** Code reads
   `phase.n_mol["H+"]` on a fresh phase (guarded by `_equilibrium_stale` —
   raises if stale). A DAE solver satisfies the same contract.

These constraints are already satisfied by the rest of the refactor. No new
design work needed.

### Accuracy monitoring

Since splitting and pseudo-DAE both have regimes where they become
inaccurate, the refactor adds automatic monitoring to flag when these regimes
are being entered. Two categories of check:

**Cheap heuristics** (always on, negligible cost):

- **pH change per step** — `|pH_{n} - pH_{n-1}| > threshold` suggests
  operator splitting is losing information
- **Newton iteration count** — `n_iters > threshold` suggests warm-start is
  far from the solution; indicator of sharp composition change or stiffness
- **Charge-balance residual after kinetic step** — `|∑(z_i × C_i)| >
  threshold` measures direct DAE constraint drift before the end-of-step
  re-solve
- **pH shift across end-of-step re-solve** — `|pH_post - pH_pre|` on the
  step-closing speciation measures splitting error directly
- **dt vs fastest reaction timescale** — computed at construction from
  initial-state Jacobian; warns if `dt > 0.1 × τ_min`
- **scipy step rejection rate** — for `ScipyODESolver`, exposing
  `solve_ivp`'s existing internal adaptive-step statistics catches stiffness
  or tolerance-boundary problems

**Rigorous diagnostics** (opt-in, expensive):

- **Step-halving (Richardson extrapolation)** — run one step at `dt`, two
  at `dt/2`, compare. 3× cost. Available via `cv.run_with_diagnostics(...)`.
- **Cross-solver comparison** — run `EulerSnapshotSolver` and
  `ScipyODESolver` on the same problem; compare trajectories. Offline
  calibration tool, not per-step.

Warnings are surfaced via a dedicated category:

```python
class AccuracyWarning(UserWarning):
    """Emitted when a numerical accuracy check flags a potential issue."""
```

Users can filter or elevate it via Python's standard `warnings` machinery:

```python
import warnings
from VLsim import AccuracyWarning

warnings.simplefilter("ignore", AccuracyWarning)   # silence
warnings.simplefilter("error",  AccuracyWarning)   # promote to exception
```

### Package-level configuration

Thresholds and throttle behaviour are set in one place, before runtime:

```python
# VLsim/config.py

from dataclasses import dataclass, field
from typing import Literal
import os

@dataclass
class WarningConfig:
    """Global accuracy-warning configuration.
    Modify ``VLsim.config.warnings`` before running simulations."""

    # Thresholds
    pH_change_threshold:       float = 0.3
    newton_iters_threshold:    int   = 15
    charge_residual_threshold: float = 1e-8
    dt_over_tau_min_threshold: float = 0.1

    # Throttle — controls how many times each distinct warning is emitted
    throttle: Literal["once", "first_N", "always", "silent"] = "once"
    first_N:  int = 10   # used only when throttle == "first_N"

    @classmethod
    def silent(cls)     -> "WarningConfig": return cls(throttle="silent")
    @classmethod
    def verbose(cls)    -> "WarningConfig": return cls(throttle="always")
    @classmethod
    def production(cls) -> "WarningConfig": return cls(throttle="first_N", first_N=3)


@dataclass
class Config:
    warnings: WarningConfig = field(default_factory=WarningConfig)


# Module-level singleton — user-modifiable before simulation runs
config = Config()

# Optional: initialise from environment variable
if "VLSIM_WARNINGS" in os.environ:
    preset = os.environ["VLSIM_WARNINGS"].lower()
    if   preset == "silent":     config.warnings = WarningConfig.silent()
    elif preset == "verbose":    config.warnings = WarningConfig.verbose()
    elif preset == "production": config.warnings = WarningConfig.production()
```

Usage patterns:

```python
import VLsim

# Fine-grained threshold tuning
VLsim.config.warnings.pH_change_threshold = 0.5
VLsim.config.warnings.throttle = "silent"

# Preset
VLsim.config.warnings = VLsim.WarningConfig.production()

# Via environment variable (CI, batch jobs)
# $ VLSIM_WARNINGS=silent python run_simulation.py
```

The `AccuracyMonitor` that emits warnings reads these values at check time,
so configuration changes take effect for all subsequent simulations in the
process.

#### Throttle modes

- `"once"` (default) — each distinct warning type emitted on first
  occurrence only; a summary count is printed at the end of the simulation
- `"first_N"` — first N occurrences of each warning type are emitted, then
  a summary
- `"always"` — every occurrence emitted (noisy but complete)
- `"silent"` — all accuracy warnings suppressed

The default `"once"` strikes the right balance — users aren't flooded, but
nothing goes completely unseen. `"silent"` is appropriate for runs where
accuracy is already verified; `"always"` for diagnosis.

#### Precedence

Python's `warnings.simplefilter("ignore", AccuracyWarning)` overrides the
config object — filters apply after the config decides whether to emit.
This means users who want stronger control (per-module filters, conversion
to errors) can layer Python's standard mechanisms on top without fighting
VLsim's config.

---

## Scope and Phasing (Directional)

This is directional only — file-level detail belongs in a subsequent
`CHEMISTRY_UNIFICATION_PLAN.md` once open questions are resolved.

- **Phase A** — Extend `Reaction` with `mode` and `log_K` fields. Extend
  elemental balance validation to include charge balance. Purely additive; no
  existing code paths change.

- **Phase B** — Rewrite `SpeciationEngine` to consume declared equilibrium
  reactions instead of Level 1/2/2.5 configurations. Keep the Newton-Raphson
  solver; change only its input-gathering layer. Add `_equilibrium_species`
  and `_equilibrium_stale` fields to `LiquidPhase`; derive
  `_equilibrium_species` at CV construction from declared equilibrium
  reactions; wire mutation methods to set the stale flag; update the
  speciation solver to write back computed concentrations and clear the flag.
  Add `PHUndefinedError`, `StaleEquilibriumError`, `phase.pH`,
  `cv.current_pH()`, `env.has_pH`, and the construction-time warning.

- **Phase C** — Derive conservation totals (`acid_totals`, `CT_P`, `CT_NH_T`)
  from phase state + declared equilibrium reactions. Remove orchestrator
  bookkeeping of these totals. `chem_env` shrinks to `t_h` + solver hints.

- **Phase D** — Clean up the two gas-liquid link leaks. Expose `alphas` in
  `PropertyResult`. Optionally declare gas-liquid partitioning as equilibrium
  reactions.

- **Phase E** — Move pKas and Henry constants out of hardcoded engine
  configuration and into shipped `ChemistryDatabase` Python modules
  (`VLsim/chemistry/databases/aqueous.py` etc.). Implement the
  `ChemistryDatabase` dataclass with `.extend()` composition. Allow users to
  swap or extend databases by Python import. Add link-construction warnings
  (Option 3) for thermo mismatch between linked CVs.

- **Phase F** — **Resolved 2026-05-16 in `chemistry-unification-4`.**
  Tag `chemistry-unification-4-shipped`. `AccuracyMonitor` is now
  attached to every `ControlVolume` (`cv._accuracy_monitor`) and
  the property solver hands it the engine's per-solve metrics.
  Five of the six cheap checks ship (pH change, Newton iterations,
  charge residual, scipy step rejections, ionic-strength regime);
  the `dt` vs τ_min check has a no-op hook with a deferred
  Jacobian-free 1/max_rate estimate. The splitting-error proxy
  (pH pre-vs-post on the same step) was dropped — needs
  solver-internal state the rest of Phase 4 doesn't justify.
  `AccuracyWarning(UserWarning)` is the category;
  `VLsim.config.warnings` holds per-check thresholds and a
  throttle mode (`once`/`first_N`/`always`/`silent`);
  `VLSIM_WARNINGS={silent,verbose,production}` initialises the
  throttle from the environment.
  `cv.run_with_diagnostics(...)` is **deferred entirely** — the
  docs/solvers.md "Accuracy monitoring" section carries the
  forward note; the rigorous step-halving + cross-solver work
  justifies its own phase when a real diagnostic workflow
  surfaces. The legacy `_warned_high_I` per-instance dedup on
  `SpeciationEngine` is gone; the engine's ionic-strength
  threshold warning migrated to `AccuracyWarning` via the
  property solver. See
  [../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md](../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md).

Phases A and B together constitute the core conceptual unification. Phases C, D,
E are hygiene / polish on top. Phase F is numerical hygiene — cheap to add,
but only meaningful once A–E land.

---

## Open Questions

These should be resolved before a full implementation plan is written.

1. ~~**Activity models.**~~ **Resolved** — see the *Runtime Packaging:
   `ThermoFramework`* section above. Activity models live on a frozen
   `ThermoFramework` object, loaded by default from the chemistry data file,
   overridable at runtime via `dataclasses.replace()`, held by the CV and
   consumed by the `ReactionEnvironment` / solvers.

2. ~~**Redox.**~~ **Resolved — deferred with forward-compatibility
   constraints.** See the *Redox: Deferred, With Forward-Compatibility
   Constraints* section above. Redox is out of scope for this refactor, but
   the refactor commits to seven design constraints (open species namespace,
   iterative charge balance, formula-less species handling, enumerable master
   variables, extensible `PropertyResult.raw` and `ThermoFramework`, general
   `log_K`) that keep redox as a future additive extension rather than a
   breaking change.

3. ~~**Data file format.**~~ **Resolved — Python modules, no separate file
   format.** See the *Database Packaging: `ChemistryDatabase`* section above.
   Shipped databases live as plain Python modules under
   `VLsim/chemistry/databases/`, each exposing a top-level `ChemistryDatabase`
   constant. User customisation is by Python import + `.extend()` composition.
   No YAML, JSON, or PHREEQC-`.dat` parser is implemented — Python's type
   system, dataclass validation, and IDE tooling replace what a schema layer
   would provide. A YAML converter can be added later as an additive utility
   if a concrete need arises.

4. ~~**Warm-start handling.**~~ **Resolved — pH stored in phase's
   `n_mol["H+"]`, with three-state guardrails.** See the *Equilibrium State:
   Storage, Staleness, and Warm-Starting* section above. Equilibrium-species
   concentrations live in the phase's `n_mol` (Option F); warm-start is
   implicit via the previous [H⁺] value; staleness is tracked via
   `_equilibrium_species` + `_equilibrium_stale` to distinguish undefined,
   stale, and fresh states; `PHUndefinedError` and `StaleEquilibriumError`
   surface misuse with clear messages; reactions-only simulations incur no
   overhead. `chem_env` sheds `logH_guess` entirely.

5. ~~**Equilibrium-kinetic coupling strength.**~~ **Resolved — defer true
   DAE, architect to permit, add accuracy monitoring with package-level
   configuration.** See the *Equilibrium-Kinetic Coupling and Accuracy
   Monitoring* section above. `EulerSnapshotSolver` (operator splitting) and
   `ScipyODESolver` (pseudo-DAE) remain the two solver paths. Four
   architectural constraints (mode declares what, not how; state partitioning
   derived; per-CV solver selection; no splitting-semantic assumptions) keep
   true DAE as a future additive extension. An `AccuracyMonitor` emits
   `AccuracyWarning` for cheap per-step checks (pH change, Newton iterations,
   charge residual, splitting-error proxy, `dt` vs τ_min, scipy step
   rejections); thresholds and throttle modes live on a package-level
   `VLsim.config.warnings` object modifiable before runtime and driven
   optionally by the `VLSIM_WARNINGS` environment variable.

---

## Relationship to [CV_UPDATE.md](../shipped/CV_UPDATE.md)

| | [CV_UPDATE.md](../shipped/CV_UPDATE.md) | This document |
|---|---|---|
| Type | Implementation plan | Rationale and direction |
| Scope | `src/core/` — ordering and data flow | `src/reactions/`, `src/speciation/`, data layer |
| File-level detail | Yes — files, lines, phases | No — directional only |
| Open questions | None — ready to implement | Several (see above) |
| Ordering | Must land first | Depends on above |

This document does not commit to implementation. It locks in the **approach** so
that a subsequent implementation plan can be written with confidence once the
open questions are settled.
