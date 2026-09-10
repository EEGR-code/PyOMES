`# PropertyCalculator Phase-Agnostic Runner — Design Note

> **Status:** Trigger-gated. Surfaced 2026-05-27 from the
> post-`simulation-class` audit as a generalisation gap: the
> `PropertyCalculator` *protocol* is phase-agnostic, but the
> *runner* in `ControlVolume.advance()` only invokes calculators
> against the liquid phase. Pick up when a real cross-phase
> property calculation is needed (gas-phase viscosity for
> compressible flow, gas-phase density for pressure-drop, solid-
> phase porosity, etc.).

## Context

`PropertyCalculator` was introduced in `state-unification` C5 as
a narrow extension protocol:

```python
class PropertyCalculator(Protocol):
    @property
    def key(self) -> str: ...
    def compute(self, phase, T_K: float, P_atm: float) -> float: ...
```

The protocol takes *a phase* — any phase — and returns one scalar
that lands on `phase.properties[key]`. By construction it is
phase-agnostic.

The runner that invokes registered calculators, however, is
liquid-only:

[`src/core/control_volume.py:555-578`](../../src/core/control_volume.py#L555-L578)

```python
def _run_property_calculators(self) -> None:
    liq = self.phases.get("liquid")
    if liq is None:
        return
    gas = self.phases.get("gas")
    P_atm = float(gas.P_atm) if gas is not None else 1.0
    T_K = float(liq.T_K)
    for calc in self.property_calculators:
        value = calc.compute(liq, T_K, P_atm)
        liq.properties[calc.key] = float(value)
```

Three consequences:

1. Gas-phase property calculators (gas viscosity, gas density,
   gas compressibility, fugacity coefficients) cannot be
   registered — the runner unconditionally hands the liquid
   phase to `.compute()`.
2. Solid-phase property calculators (effective diffusivity in a
   porous solid, particle-scale density, internal porosity
   updates as a function of conversion) likewise have nowhere
   to live.
3. Multi-liquid systems (liquid-liquid extraction, biphasic
   solvent systems) have only one of the two liquids visible
   to the calculator.

For today's fermenter use cases none of this matters — the only
shipped `PropertyCalculator` is liquid viscosity, and biomass /
substrate concentration is naturally a liquid-phase property.
The shape of the runner is fermenter-coincident, not framework-
intentional.

## Why this matters for the framework

The package is intended to generalise across process-simulation
domains. Each of these picks up gas-phase or solid-phase property
calculations as a first-class concern:

- **HPLC / packed columns.** Mobile-phase viscosity drives the
  pressure-drop calculation through the Carman–Kozeny relation.
  Mobile-phase density is needed for any compressibility
  correction. Solid-phase porosity matters for the effective
  diffusivity in the stationary phase.
- **Anaerobic digestion compressed-headspace work.** Gas-phase
  fugacity coefficients (Peng–Robinson is already available at
  [`src/equilibria/peng_robinson.py`](../../src/equilibria/peng_robinson.py)
  but unreachable — see the audit) are needed for high-pressure
  biogas storage modelling.
- **Cell culture / perfusion.** Membrane permeability is a
  function of solid-phase properties that change with fouling.
- **Cryogenic / hyperthermophilic systems.** Gas-phase density
  diverges significantly from ideal-gas at non-ambient
  conditions; a calculator that overrides `gas.P_atm`'s ideal-gas
  derivation would matter.

In each case the *protocol* survives unchanged; what fails today
is the *invocation*.

## Proposed shape (when picked up)

The minimal change is to have each `PropertyCalculator` declare
which phase(s) it targets, and have the runner iterate phases:

```python
class PropertyCalculator(Protocol):
    @property
    def key(self) -> str: ...
    @property
    def phase_key(self) -> str:  # e.g. "liquid", "gas", "solid"
        ...
    def compute(self, phase, T_K: float, P_atm: float) -> float: ...
```

```python
def _run_property_calculators(self) -> None:
    for calc in self.property_calculators:
        phase = self.phases.get(calc.phase_key)
        if phase is None:
            continue
        T_K = float(phase.T_K)
        # Pressure: from the same phase if it exposes P_atm;
        # else from the gas phase if any; else 1.0 atm.
        if hasattr(phase, "P_atm"):
            P_atm = float(phase.P_atm)
        else:
            gas = self.phases.get("gas")
            P_atm = float(gas.P_atm) if gas is not None else 1.0
        value = calc.compute(phase, T_K, P_atm)
        phase.properties[calc.key] = float(value)
```

`phase.properties` already exists on `GasPhase` (it doesn't
today, but adding it is a one-line addition mirroring
`LiquidPhase.properties`); the runner change is the load-bearing
part.

## Open design questions to resolve when picked up

1. **Default phase key** — should `phase_key` default to
   `"liquid"` to preserve the current calling convention, or
   should it be required? Required is more honest but breaks
   the existing single shipped calculator.
2. **Multi-phase calculators** — does the protocol allow a
   single calculator to target multiple phases (e.g. one
   instance computing gas *and* liquid density)? Or always
   one-calculator-per-phase? The latter is simpler; the former
   matches how a single thermodynamic model often produces both.
3. **`GasPhase.properties` dict** — `LiquidPhase` has one;
   `GasPhase` and `SolidPhase` do not. Either add `.properties`
   to the `Phase` protocol or accept the asymmetry. Symmetric
   storage is the right answer if gas-phase calculators are
   first-class.
4. **`AccuracyMonitor` / `ConservationMonitor` interaction** —
   neither monitor reads `phase.properties` today, so this
   change is invisible to monitoring. If gas-phase density
   becomes a monitored quantity, the monitor needs updating
   too.

## Out of scope

- Speciation is **not** a `PropertyCalculator`. Acid-base
  equilibria are state-completion handled by the
  `SpeciationEngine` writing molecular species back to
  `phase.n_mol`. This note does not propose moving any
  speciation work into `PropertyCalculator`.
- The pressure-passing convention (`P_atm` argument) — keep as
  is for now; this note doesn't aim to redesign the calculator
  signature.

## Trigger conditions

Any one of:

1. A real HPLC / packed-column workload needs mobile-phase
   viscosity or gas-phase density.
2. A high-pressure biogas storage model needs gas-phase
   fugacity coefficients (Peng–Robinson is already implemented
   but unreachable from the framework runner).
3. A cell-culture or perfusion model needs solid-phase
   property tracking.
4. The single-day cost of just adding gas-phase invocation
   (estimated: ~30 LOC + tests + a planning note like this one
   becoming a checklist) is judged low enough to do
   pre-emptively as part of a "framework polish" bundle.

## Relationship to other phases

- **HPLC trigger bundle (added 2026-05-28).** The three-way
  bundle with
  [`SNAPSHOT_PHASE_AGNOSTIC.md`](SNAPSHOT_PHASE_AGNOSTIC.md)
  and [`CONTAINER_LAYERING.md`](CONTAINER_LAYERING.md) is
  triggered by the same workload (HPLC mobile/stationary
  phases). Picking up HPLC pulls all three in.
- **Hard dependency from
  [`SNAPSHOT_PHASE_AGNOSTIC.md`](SNAPSHOT_PHASE_AGNOSTIC.md)
  Q5.** That note's open question about whether `y_gas` and
  gas `P_atm` should migrate into the general `phase_properties`
  channel "via `PropertyCalculator`" cannot be answered without
  this note's gas-phase runner shipping first. The two notes
  are not merely "pairs naturally" on Q5; SNAPSHOT Q5 depends on
  PROPERTY_CALCULATORS.
- **Shares the "no-liquid CV" thread with
  [FRAMEWORK_POLISH P2](FRAMEWORK_POLISH.md#p2-_build_reaction_environment-gas-only-cv-defaults--severity-medium).**
  P2 generalises `_build_reaction_environment` to gas-only CVs;
  this note generalises the property-calculator runner to the
  same shape. Both fix the same conceptual root (stop treating
  the liquid phase as mandatory) at different layers and share
  testing infrastructure (a gas-only CV fixture). If both ship,
  do them in the same branch.
- **Interaction with
  [CHEMISTRY_UNIFICATION_PLAN.md](CHEMISTRY_UNIFICATION_PLAN.md)
  Phase 3b — open question.** 3b introduces `ThermoFramework`
  absorbing activity-model configuration. The Peng–Robinson
  fugacity use case named under "Why this matters" sits at the
  seam: does fugacity configuration live on `ThermoFramework`
  or on a gas-phase `PropertyCalculator`? The two notes do not
  pre-commit; a high-pressure biogas storage workload would
  force the question.
- **Independent of** `RUN_HISTORY`,
  `SPECIATION_LEVEL_RETIREMENT`, and `CUFERMENTER_SUNSET` phases.
- **Naturally bundles with** the snapshot/recorder phase-key
  generalisation flagged in the same audit — both are framework
  surfaces where the liquid/gas-phase assumption leaked from
  fermenter use cases into supposedly general machinery.
- **Follows** `state-unification` (shipped 2026-05-22), which
  introduced the `PropertyCalculator` protocol in C5. The
  protocol is intentionally narrow and stable; this is a
  runner-side gap, not a protocol redesign.
