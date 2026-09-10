# NR_PRECIPITATION_CV_INTEGRATION — Phase 2 of 2

> **Status:** Upcoming — not yet started.
> **Depends on:** [NR_PRECIPITATION_SPECIATION.md](../shipped/NR_PRECIPITATION_SPECIATION.md)
>   (Phase 1) shipped 2026-06-23 (`nr-precipitation-speciation-shipped`).
> **Tag to follow:** `nr-precipitation-cv-integration-shipped`.

---

## Background and motivation

Phase 1 (`NR_PRECIPITATION_SPECIATION`) adds precipitation equilibrium to
`NRChemicalEquilibriumEngine.solve()` and returns mineral amounts in the output dict.
The engine is fully functional for scripting and notebook use.

Phase 2 wires the precipitation result into the `ControlVolume` lifecycle so
that mineral amounts persist across timesteps in a running simulation:

- `SolidPhase.n_mol` accumulates precipitated mineral across ODE steps.
- `_read_from_phases` sums liquid + solid contributions when computing
  component totals, conserving mass across the solid/liquid boundary.
- The three phase types (`LiquidPhase`, `GasPhase`, `SolidPhase`) behave
  consistently: a CV that could in principle develop a solid/gas phase always
  carries the relevant empty phase object from construction, avoiding
  mid-simulation lazy creation (which would invalidate ODE state references).

---

## Key design decisions

### 1. Phase-type consistency across LiquidPhase, GasPhase, SolidPhase

**Rule:** if a `ReactionSystem` contains reactions that imply a phase, the
`ControlVolume` must be constructed with that phase attached — even if it
starts empty.

| `ReactionSystem` contains | CV must carry |
|---|---|
| Any `KineticReaction` or `EquilibriumReaction` (liquid) | `LiquidPhase` |
| Any `cross_phase_equilibria` or `KineticGasLiquidLink` | `GasPhase` |
| Any `precipitation_equilibria` | `SolidPhase` |

A missing phase that is implied raises a clear `ConfigurationError` at
`ControlVolume.__init__` time (not at solve time).  This is consistent with
how a missing gas phase already errors today when gas-liquid mass transfer is
declared.

**Empty phase construction:** `SolidPhase(n_mol={})` is valid and zero-cost.
The precipitation solver writes to it; the ODE integrator never touches it
(no ODE state in the solid phase in Phase 2 scope).

### 2. _read_from_phases sums liquid + solid contributions

When `NRChemicalEquilibriumEngine._read_from_phases(phases)` is called and a
`SolidPhase` is present, the total for each component includes the solid
contribution via the dissolution reaction stoichiometry.

For calcite:
```
CT_Ca_total  = n_mol_liq["Ca++"] / V_L
             + n_mol_solid["CaCO3"] / V_L   # stoich coeff 1 for Ca
CT_CO2_total = sum(n_mol_liq[sp] for sp in carbonate_species) / V_L
             + n_mol_solid["CaCO3"] / V_L   # stoich coeff 1 for C
```

The stoichiometric coefficients are read from the `precipitation_equilibria`
reactions (the liquid-phase entries).  The engine iterates over all
`precipitation_equilibria` and for each solid species found in
`SolidPhase.n_mol`, adds `n_mol_solid[solid_id] / V_L × |ν_component|` to the
relevant component total.

### 3. SolidPhase writeback after solve

After `_solve_with_precipitation` converges, `_writeback` is extended to update
`SolidPhase.n_mol`:

```python
for mineral_name, result in out["minerals"].items():
    xi = result["xi_mol_L"]                    # mol/L precipitated
    solid_species_id = ...                     # from dissolution reaction
    phases["solid"].n_mol[solid_species_id] = xi * V_L
```

The dissolved totals in `LiquidPhase.n_mol` are adjusted accordingly through
the existing speciation writeback path (the inner NR solve already sees
adjusted effective totals and writes back the correct dissolved concentrations).

### 4. Strong-ion Ca²⁺ tracking across timesteps

Ca²⁺ currently lives in `strong_ions` (charge balance only, not in
`n_mol`).  For precipitation to work correctly across timesteps the *total*
Ca²⁺ (dissolved + precipitated) must be conserved.

`_read_from_phases` reads the current dissolved Ca²⁺ from
`n_mol["Ca++"]` (if present) or from the strong-ion mapping, and adds the
solid CaCO₃ contribution to get `CT_Ca_total`.  After the precipitation solve,
the dissolved Ca²⁺ is implied by `CT_Ca_total − ξ`; the solver
writes nothing to `n_mol["Ca++"]` explicitly — the kinetic ODE state for Ca²⁺
is not tracked (Ca²⁺ has no kinetic reactions).

If a future phase introduces Ca²⁺ kinetics (e.g. Ca²⁺ uptake by biomass),
Ca²⁺ will need to enter `n_mol` as a proper ODE state.  The architecture
is already consistent with this: `n_mol["Ca++"]` is a valid key; adding it
simply enables the ODE integrator to track it.

---

## Scope

**In scope:**

- `ControlVolume.__init__` validates that a `SolidPhase` is present when
  `ReactionSystem.precipitation_equilibria` is non-empty; raises
  `ConfigurationError` otherwise.
- Consistent validation for all three phase types (GasPhase, LiquidPhase,
  SolidPhase) — the same rule applied uniformly.
- `NRChemicalEquilibriumEngine._read_from_phases` sums liquid + solid contributions
  using dissolution reaction stoichiometry.
- `NRChemicalEquilibriumEngine._writeback` extended to write mineral amounts to
  `SolidPhase.n_mol`.
- `SolidPhase` exposed in `ControlVolume` construction API and in
  `cv.phases` dict under key `"solid"`.
- Integration tests: a `ControlVolume` with calcite precipitation runs
  multiple timesteps; mineral amount at steady state matches analytical Ksp
  solution; mass is conserved.

**Out of scope:**

- ODE-tracked solid state (solid mass in the ODE state vector) — the solid
  phase is algebraic in this phase.
- Kinetic dissolution / precipitation rates (separate future phase if needed).
- Multiple co-precipitating minerals sharing a common component — deferred
  until at least one real use case requires it.
- Precipitation-driven pH control loops.

---

## Files changed

| File | Nature of change |
|---|---|
| `src/core/phases.py` | Ensure `SolidPhase` is a fully initialised peer of `LiquidPhase` and `GasPhase` (review any asymmetries) |
| `src/control_volume.py` | `__init__`: accept `solid_phase` kwarg; validate phase-type consistency; expose in `self.phases["solid"]` |
| `src/chemical_equilibrium/nr_engine.py` | `_read_from_phases`: sum solid contribution; `_writeback`: write mineral amounts to `SolidPhase.n_mol` |
| `tests/` | New integration test: multi-timestep CV with calcite precipitation; mass conservation check |
| `docs/upcoming/NR_PRECIPITATION_SPECIATION.md` | Add "Shipped" banner |
| `docs/shipped/NR_PRECIPITATION_SPECIATION.md` | Move on ship |

---

## Trigger conditions

This phase should be started when:

1. Phase 1 (`nr-precipitation-speciation`) is shipped and tagged.
2. There is a concrete use case requiring precipitation to persist across
   timesteps — e.g. a lime-softening simulation, a struvite reactor, or a
   batch model where CaCO₃ scale builds up over time.

There is no urgency to start Phase 2 immediately after Phase 1: the Phase 1
speciation-only implementation is sufficient for all current demo and
assessment use cases.

---

## How to start

1. Confirm Phase 1 tag `nr-precipitation-speciation-shipped` exists on `main`.
2. Create branch `nr-precipitation-cv-integration` off `main`.
3. Work through the files in the order listed in *Files changed*.
4. Run the full test suite (`pytest`) after each file.
5. When green: merge `--no-ff`, tag `nr-precipitation-cv-integration-shipped`,
   delete branch.
