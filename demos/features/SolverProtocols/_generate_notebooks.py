"""Generate the SolverProtocols feature notebooks as valid .ipynb JSON.

This folder covers PyOMES's two orthogonal solver protocols:
`StepSolver` (Axis 1 — per-CV physics) and `SystemSolver` (Axis 2 —
whole-system orchestration across CVs, links, and controllers). Both
protocols are deliberately thin — one method each — so a bespoke
solver is a small, self-contained class. See docs/solvers.md for the
full solver landscape these notebooks are drawn from, and
docs/phases-shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md for the design
review that shipped item 8 (the walkthrough notebook).

Run once from the repo root:
    python demos/features/SolverProtocols/_generate_notebooks.py
"""
import json
from pathlib import Path

HERE = Path(__file__).parent

METADATA = {
    "kernelspec": {
        "display_name": "biosteam",
        "language": "python",
        "name": "python3",
    },
    "language_info": {
        "codemirror_mode": {"name": "ipython", "version": 3},
        "file_extension": ".py",
        "mimetype": "text/x-python",
        "name": "python",
        "nbconvert_exporter": "python",
        "pygments_lexer": "ipython3",
        "version": "3.12.3",
    },
}


def nb(*cells):
    return {"nbformat": 4, "nbformat_minor": 5, "metadata": METADATA, "cells": list(cells)}


def md(id_, src):
    return {"cell_type": "markdown", "id": id_, "metadata": {}, "source": src}


def code(id_, src):
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": id_,
        "metadata": {},
        "outputs": [],
        "source": src,
    }


# ═══════════════════════════════════════════════════════════════════════════
#  0  Architecture overview (pure markdown, no code — matches the
#     ChemicalEquilibriumProtocol/0_README.ipynb convention)
# ═══════════════════════════════════════════════════════════════════════════

overview_nb = nb(
    md("title", """\
# `SolverProtocols` — Architecture Overview

PyOMES has two orthogonal solver axes, each a thin, deliberately
minimal `typing.Protocol`:

- **Axis 1 — `StepSolver`.** Per-CV physics: one method,
  `solve_step(cv, dt_h, t_h, external_source_terms) -> AdvanceResult`.
  Advances a single `ControlVolume` by one macro timestep.
- **Axis 2 — `SystemSolver`.** Whole-system orchestration: one method,
  `advance_system(sim, dt_h, t_h)`. Given a whole `Simulation`
  (multiple CVs, inter-CV links, controllers, profiles), decides how
  one macro step gets composed.

They compose orthogonally — `Simulation(solver=..., system_solver=...)`
— pick any `StepSolver` per CV under any `SystemSolver`. Neither axis
has a privileged implementation: `solver=None` and `system_solver=None`
are pure sugar for one particular rung each (`SequentialAdvanceSolver`,
`ExplicitEulerSystemSolver`), not special-cased code paths.

This notebook is the map: what each protocol contract actually is, how
the shipped rungs on each axis compare, and which one to reach for.\
"""),

    md("notebooks", """\
## Notebooks in this folder

| Notebook | Covers |
|---|---|
| [01_writing_a_custom_solver.ipynb](01_writing_a_custom_solver.ipynb) | Writing your own `StepSolver` and `SystemSolver` — two worked examples, each inverting a real shipped solver's documented design decision |

Launch from the repo root with
`jupyter lab demos/features/SolverProtocols/`.\
"""),

    md("axis1", """\
## Axis 1 — `StepSolver`: the three shipped rungs

| Solver | Semantics | Earns its keep when |
|---|---|---|
| `SequentialAdvanceSolver` (`solver=None` default) | Sequential operator split: speciation → feed → react → transfer; each sub-step sees post-previous state | Default path; small `dt_h`; simple models without strong feed/reaction coupling |
| `SimultaneousEulerSolver` | Simultaneous explicit Euler: all sub-systems read the same frozen snapshot, deltas summed and applied together, with swappable clamping | Larger `dt_h` where operator-splitting bias matters; matches the pre-refactor semantics ADM1 was validated against |
| `SimultaneousAdaptiveSolver` | `scipy.integrate.solve_ivp` adaptive; explicit (DOP853) or implicit (Radau, BDF) methods; optional speciation freeze (DAE-style, BSM2/PyADM1 convention) | Stiff models (BSM2 anaerobic digestion, fast H₂ kinetics); long macro timesteps explicit Euler cannot survive |

All three require only a `"liquid"` phase (speciation and reactions
are liquid-scoped) — a gas phase, solid phase, additional phases, or
none of those beyond liquid are all supported.

**Ownership guard.** If a CV is owned by a `Simulation`, calling
`cv.advance(...)` directly (rather than through `sim.run(...)`) emits
`OrchestrationWarning` — inter-CV links, controllers, and profiles
won't be applied for that step. Every `SystemSolver` calls the trusted
internal `cv._advance_unchecked(...)` instead, which bypasses the
guard.\
"""),

    md("axis1-instantiation", """\
### Instantiating each `StepSolver`

All three are plain classes — construct one and pass it as
`cv.advance(dt_h, t_h, solver=...)` or
`Simulation(solver=...)`/`Simulation(solver={cv_key: ...})`. Every
keyword below has a default, so `SolverClass()` with no arguments is
always valid; the "example" column shows a non-default configuration
for context.

| Solver | Default instantiation | Example with keywords | Keyword arguments |
|---|---|---|---|
| `SequentialAdvanceSolver` | `SequentialAdvanceSolver()` | `SequentialAdvanceSolver(clamp_fn=None)` | `clamp_fn` (default `proportional_clamp`) — swappable non-negativity clamp; `floor_clamp` or a bespoke composite are drop-in alternatives; `None` disables clamping entirely |
| `SimultaneousEulerSolver` | `SimultaneousEulerSolver()` | `SimultaneousEulerSolver(clamp_fn=floor_clamp)` | `clamp_fn` (default `proportional_clamp`) — same contract as above, applied once to the combined multi-source delta dict |
| `SimultaneousAdaptiveSolver` | `SimultaneousAdaptiveSolver()` | `SimultaneousAdaptiveSolver(method="Radau", rtol=1e-6, atol=1e-9)` | `method` (default `"DOP853"`) — `scipy.integrate.solve_ivp` method, `"Radau"`/`"BDF"` for stiff systems; `rtol`/`atol` (default `1e-6`/`1e-9`) — solver tolerances; `max_step` (default `1.0` h) — caps the internal adaptive step; `freeze_speciation` (default `False`) — solve speciation once per macro step instead of at every derivative evaluation (BSM2/PyADM1 DAE convention); `use_engine_jacobian` (default `False`) — use the speciation engine's analytical `jacobian_dz_dy()` instead of finite differences (only with a `GrayBoxEngineProtocol` engine and an implicit `method`) |

`clamp_fn`/`floor_clamp`/`proportional_clamp` live in
[`src/core/clamping.py`](../../../src/core/clamping.py) — see the
Clamping section below.\
"""),

    md("clamping", """\
## Clamping — non-negativity handling

All three `StepSolver`s restore a physical invariant (non-negative
inventory) the raw numerics don't guarantee on their own:

| Function | Used by | Behaviour |
|---|---|---|
| `proportional_clamp` (default) | `SequentialAdvanceSolver`, `SimultaneousEulerSolver` | Scales a removal rate down so the result lands at exactly zero, preserving relative stoichiometry |
| `floor_clamp` | Either discrete-step solver, passed explicitly | Per-species floor (optionally at `eps`) — simpler, doesn't preserve relative stoichiometry |
| `floor_nonnegative` | `SimultaneousAdaptiveSolver` (always, not swappable) | Floors a raw ODE state array — a different failure mode (adaptive-integrator overshoot), not a `(deltas, current_mol, dt_h)` decision |

`clamp_fn` is a plain callable operating on an **arbitrary subset** of
species — this is what lets a bespoke composite match an external
reference model's specific per-species convention, using the shared
functions as ingredients rather than an all-or-nothing choice.
`clamp_fn=None` disables clamping entirely (a genuine way to let a
species go negative for diagnosis — `AccuracyMonitor.check_negative_mole`
flags it); `AccuracyMonitor.check_clamp_invoked` independently flags
when clamping actually changed a rate, naming the before/after values.\
"""),

    md("axis2", """\
## Axis 2 — `SystemSolver`: the five shipped rungs

| Solver | Sequence per macro step | Earns its keep when |
|---|---|---|
| `ExplicitEulerSystemSolver` (default) | `profiles(t) → links(dt) → cv.advance(dt) → controllers` | τ_CFL ≫ dt_h and slow controllers — the common case |
| `StrangSplittingSystemSolver` | `profiles(t) → links(dt/2) → cv.advance(dt) → links(dt/2) → controllers` | Free 2nd-order splitting upgrade over Euler, ~zero extra cost |
| `MultirateSystemSolver` | `profiles(t) → cv.advance(dt) → [links(dt/M)] × M → controllers` | Fast inter-CV circulation (τ_CFL comparable to dt_h) without a global step-size cut |
| `ImplicitTransportSystemSolver` | Sparse implicit link-flux solve, unconditionally stable | Many-CV compartmental transport models (HPLC, packed beds), no tight controllers |
| `MonolithicODESolver` | Single `solve_ivp` co-integrating CV species, inter-CV transport, and controller differential state | Tight feedback loops (T_c < dt_h), multi-rate digital controllers, co-integrated integral states (PI/PID, observers) |

All five call `cv._advance_unchecked(...)` internally — except
`MonolithicODESolver`, which integrates every CV directly via
`cv.compute_rhs()` inside its shared `solve_ivp` call and therefore
**raises** if a per-CV `Simulation(solver=...)` is also configured
(there's no legitimate way for it to be silently correct, since it
never consults that config).\
"""),

    md("choosing", """\
## Choosing a solver

| If your model has… | Start with |
|---|---|
| Mild kinetics, small `dt_h`, no special needs | `cv.advance()` (`SequentialAdvanceSolver`, the default) |
| Aggressive feeds or strong feed/reaction coupling, larger `dt_h` | `SimultaneousEulerSolver` |
| Stiff fast kinetics (BSM2 H₂), pH-coupled transfer, long timesteps | `SimultaneousAdaptiveSolver(method="Radau")` or `method="BDF"` |
| Reproducing BSM2 / PyADM1 validation results | `SimultaneousAdaptiveSolver(freeze_speciation=True, method="Radau")` |
| τ_CFL ≫ dt_h, slow controllers | `ExplicitEulerSystemSolver` (the default) |
| Fast inter-CV circulation | `MultirateSystemSolver` |
| Tight feedback loops / multi-rate digital controllers | `MonolithicODESolver` |
| Reproducing a reference model's own discrete interleaving convention | Write your own `SystemSolver` — see [01_writing_a_custom_solver.ipynb](01_writing_a_custom_solver.ipynb) |

If a model works with the defaults and the results look right, there's
no reason to switch. The non-default solvers exist because real models
exist that the defaults cannot integrate accurately or stably.\
"""),

    md("cross-references", """\
## Cross-references

- [`../../../docs/solvers.md`](../../../docs/solvers.md) — full
  solver landscape (strengths/weaknesses, configuration, accuracy +
  conservation monitoring) this notebook is a condensed map of.
- [`../../../docs/design/SOLVER_ARCHITECTURE.md`](../../../docs/design/SOLVER_ARCHITECTURE.md) —
  full two-axis design rationale, the DAE/SUNDIALS Phase F/G
  placeholders, and identified-but-unbuilt extensions (SIA, reactive
  D_eff).
- [`../../../docs/phases-shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md`](../../../docs/phases-shipped/STEP_SOLVER_INTERFACE_REFINEMENT.md) —
  the design review that shipped the ownership guard, the
  gas/liquid generalization, the shared `clamp_fn` module, and this
  folder's walkthrough notebook.
- [`../../../docs/phases-shipped/ORDERING.md`](../../../docs/phases-shipped/ORDERING.md) —
  `SequentialAdvanceSolver`'s operator-splitting ordering rationale,
  inverted by [01_writing_a_custom_solver.ipynb](01_writing_a_custom_solver.ipynb)'s
  Axis 1 example.\
"""),
)

save_path_overview = HERE / "0_README.ipynb"
with open(save_path_overview, "w", encoding="utf-8") as f:
    json.dump(overview_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path_overview}")


# ═══════════════════════════════════════════════════════════════════════════
#  1  Writing a custom StepSolver and SystemSolver (executed, code-level)
# ═══════════════════════════════════════════════════════════════════════════

SETUP = """\
from PyOMES.core import (
    AdvanceResult,
    AdvectiveLink,
    ControlVolume,
    LiquidPhase,
    SequentialAdvanceSolver,
    Simulation,
)
from PyOMES.core.system_solver import ExplicitEulerSystemSolver
from PyOMES.reactions import EquilibriumReaction, ReactionSystem, StoichiometryEntry
from PyOMES.chemistry.common_species import H_plus, OH_minus, H2O
from typing import Any, Dict, Optional

print("Imports OK")
"""

walkthrough_nb = nb(
    md("title", """\
# Writing your own `StepSolver` and `SystemSolver`

See [0_README.ipynb](0_README.ipynb) for the protocol overview and the
full shipped-solver landscape on both axes. This notebook writes one
custom solver on each axis — not toy reorderings, but each one
inverting a real, documented design decision already made by a shipped
solver, so the comparison teaches an actual trade-off rather than an
arbitrary difference:

- **Axis 1:** `SpeciationAfterFeedStepSolver` reverses
  `SequentialAdvanceSolver`'s documented ordering (see
  [`docs/phases-shipped/ORDERING.md`](../../../docs/phases-shipped/ORDERING.md)) —
  the shipped default pins speciation to the *pre-feed* state within
  one step, so a large feed's pH-shifting effect only becomes visible
  on the *next* step. This solver feeds first, so the same-step pH
  already reflects the feed.
- **Axis 2:** `AsymmetricThreeStageSystemSolver` follows the same
  extension pattern `StrangSplittingSystemSolver` established
  (`sim._apply_links` / `cv._advance_unchecked` /
  `sim._invoke_controllers`), just composed into an asymmetric 3-stage
  split instead of one symmetric half-step.\
"""),

    code("setup", SETUP),

    # ── 1. Axis 1 — custom StepSolver ───────────────────────────────────────

    md("s1-md", """\
## 1  Axis 1 — a custom `StepSolver`

`SequentialAdvanceSolver` (the shipped default — `solver=None` is pure
sugar for it) runs speciation *before* applying that step's feed, so
the feed's effect on pH is deliberately invisible until the *next*
step (the "operator-splitting contract" in `ORDERING.md`). Reversing
that order is a one-class change:

```
1. External source terms + boundary fluxes (feed)   <- moved first
2. Speciation solve                                  <- now post-feed
3. Reactions (post-feed, post-speciation state)
4. Internal transfer
```

Neither ordering is "more correct" in the abstract — which one you
want depends on what you're reproducing: PyOMES's own convention, or an
external reference model that feeds first.\
"""),

    code("s1-class", """\
class SpeciationAfterFeedStepSolver:
    \"\"\"Feed the CV *before* solving speciation, reversing
    ``SequentialAdvanceSolver``'s documented pre-feed-speciation
    contract.

    This is the minimum viable custom ``StepSolver`` -- one method,
    ``solve_step(cv, dt_h, t_h, external_source_terms)``, returning an
    ``AdvanceResult``. Everything else (the shipped rungs, this one,
    and any solver a user writes) satisfies the exact same protocol.
    \"\"\"

    def solve_step(
        self,
        cv: Any,
        dt_h: float,
        t_h: float = 0.0,
        external_source_terms: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> AdvanceResult:
        # 1. Feed first.
        if external_source_terms:
            for phase_key, rates in external_source_terms.items():
                if phase_key in cv.phases:
                    cv.phases[phase_key].apply_flux(rates, dt_h)
        for boundary in cv.boundaries:
            flux = boundary.compute_flux(cv, dt_h)
            if boundary.phase_key in cv.phases:
                cv.phases[boundary.phase_key].apply_flux(flux, dt_h)

        # 2. Speciation on the POST-feed state (the reversed step).
        if cv.reaction_system is not None and hasattr(cv.reaction_system, "engine"):
            engine = cv.reaction_system.engine
            if engine is not None:
                liq = cv.phases.get("liquid")
                if liq is not None:
                    result = engine.solve(phases=cv.phases, T_K=float(liq.T_K))
                    result.apply_to_phases(cv.phases)

        # 3. Reactions, against the now-post-feed speciation.
        rxn_sources = None
        if cv.reaction_system is not None:
            rhs = cv.compute_reaction_rates(t_h)
            rxn_sources = {}
            for pk, sp_rates in rhs.items():
                rxn_sources[pk] = {}
                phase = cv.phases[pk]
                for sp, rate in sp_rates.items():
                    phase.n_mol[sp] = max(0.0, phase.n_mol.get(sp, 0.0) + rate * dt_h)
                    rxn_sources[pk][sp] = rate

        # 4. Internal transfer last.
        transfer_diag = cv.step_internal_transfer(dt_h)

        return AdvanceResult(transfer=transfer_diag, reaction_sources=rxn_sources)


def make_acid_base_cv() -> ControlVolume:
    \"\"\"A liquid CV with a water-dissociation equilibrium -- enough to
    show a pH shift from a strong-base feed within a single step.\"\"\"
    liq = LiquidPhase(n_mol={"Na+": 0.0}, V_L=1.0, T_K=298.15)
    rxn = EquilibriumReaction(
        stoichiometry=[
            StoichiometryEntry(species=H2O, phase="liquid", coefficient=-1.0),
            StoichiometryEntry(species=H_plus, phase="liquid", coefficient=+1.0),
            StoichiometryEntry(species=OH_minus, phase="liquid", coefficient=+1.0),
        ],
        log_K=-14.0,
        balance_elements=("H", "O"),
        label="water",
    )
    system = ReactionSystem([rxn], label="water_equilibrium")
    return ControlVolume(phases={"liquid": liq}, reaction_system=system, label="cv")

print("SpeciationAfterFeedStepSolver defined")
"""),

    md("s1-run-md", """\
### Run it

A strong-base feed (`Na+`/`OH-`) large enough to noticeably shift pH
within `dt_h`, applied under both the shipped default and the custom
ordering, starting from the same initial state.\
"""),

    code("s1-run", """\
feed = {"liquid": {"Na+": 0.02, "OH-": 0.02}}  # mol/h, applied over dt_h

cv_default = make_acid_base_cv()
cv_custom = make_acid_base_cv()

# Seed derived species (H+/OH-) on a fresh CV via one speciation solve,
# so pH is readable before either custom ordering runs.
for cv in (cv_default, cv_custom):
    engine = cv.reaction_system.engine
    result = engine.solve(phases=cv.phases, T_K=float(cv.phases["liquid"].T_K))
    result.apply_to_phases(cv.phases)
pH_before = cv_default.phases["liquid"].pH

cv_default.advance(
    dt_h=1.0, t_h=0.0, external_source_terms=feed,
    solver=SequentialAdvanceSolver(),
)
cv_custom.advance(
    dt_h=1.0, t_h=0.0, external_source_terms=feed,
    solver=SpeciationAfterFeedStepSolver(),
)

print(f"pH before the step:                          {pH_before:.3f}")
print(f"pH after (SequentialAdvanceSolver, shipped):  {cv_default.phases['liquid'].pH:.3f}")
print(f"pH after (SpeciationAfterFeedStepSolver):     {cv_custom.phases['liquid'].pH:.3f}")
"""),

    md("s1-takeaway", """\
Both are valid Lie-Trotter splittings of the same physics — they
differ only in which state the same-step speciation solve sees.
`SequentialAdvanceSolver`'s choice is documented in `ORDERING.md`;
reproducing a reference model that feeds first is a one-class change,
not a fork of the framework.\
"""),

    # ── 2. Axis 2 — custom SystemSolver ─────────────────────────────────────

    md("s2-md", """\
## 2  Axis 2 — a custom `SystemSolver`

`StrangSplittingSystemSolver` upgrades `ExplicitEulerSystemSolver` to
2nd-order accuracy with one extra `_apply_links` call:
`links(dt/2) -> cv.advance(dt) -> links(dt/2)`. A bespoke `SystemSolver`
is no more than composing the same three primitives
(`sim._apply_links` / `cv._advance_unchecked` /
`sim._invoke_controllers`) in whatever order your model — or the
reference model you're reproducing — needs. Here: an asymmetric
3-stage split instead of one symmetric half-step.\
"""),

    code("s2-class", """\
class AsymmetricThreeStageSystemSolver:
    \"\"\"Bespoke 3-stage inter-CV interleaving: an asymmetric link split
    around two half-duration CV advances, instead of
    ``StrangSplittingSystemSolver``'s single symmetric half-step.

    Sequence per macro step::

        profiles(t) -> links(dt*w0) -> cv.advance(dt/2)
                    -> links(dt*w1) -> cv.advance(dt/2)
                    -> links(dt*w2) -> controllers

    with asymmetric weights ``w0 + w1 + w2 == 1`` (default
    ``(0.2, 0.5, 0.3)``).
    \"\"\"

    def __init__(self, weights=(0.2, 0.5, 0.3)):
        if len(weights) != 3 or abs(sum(weights) - 1.0) > 1e-9:
            raise ValueError("weights must be a 3-tuple summing to 1.0")
        self.weights = weights

    def advance_system(self, sim: "Simulation", dt_h: float, t_h: float) -> Any:
        w0, w1, w2 = self.weights
        profile_actions = sim._invoke_profiles(t_h=t_h)

        link_records = sim._apply_links(dt_h * w0)

        results: Dict[str, "AdvanceResult"] = {}
        for cv_key, cv in sim.cvs.items():
            solver = sim._solver_for(cv_key)
            results[cv_key] = cv._advance_unchecked(dt_h / 2.0, t_h, solver=solver)

        link_records += sim._apply_links(dt_h * w1)

        for cv_key, cv in sim.cvs.items():
            solver = sim._solver_for(cv_key)
            # Overwrites the first-half AdvanceResult -- callers that need
            # both would extend AdvanceResult or accumulate a list instead.
            results[cv_key] = cv._advance_unchecked(dt_h / 2.0, t_h + dt_h / 2.0, solver=solver)

        link_records += sim._apply_links(dt_h * w2)

        controller_actions = [
            a for a in sim._invoke_controllers(t_h=t_h, dt_h=dt_h, results=results)
            if a is not None
        ]
        for action in controller_actions:
            sim._apply_controller_action(action, dt_h=dt_h)

        return results, link_records, controller_actions, profile_actions


def build_two_cv_sim(system_solver=None) -> Simulation:
    liq_r = LiquidPhase(n_mol={"S": 10.0}, V_L=1.0, T_K=310.0)
    liq_c = LiquidPhase(n_mol={"S": 0.0}, V_L=1.0, T_K=310.0)
    return Simulation(
        cvs={
            "reactor": ControlVolume(phases={"liquid": liq_r}, label="reactor"),
            "recycle": ControlVolume(phases={"liquid": liq_c}, label="recycle"),
        },
        links=[
            AdvectiveLink(
                _source_cv_key="reactor", _source_phase_key="liquid",
                _sink_cv_key="recycle", _sink_phase_key="liquid",
                Q_L_per_h=2.0, _label="fwd",
            ),
            AdvectiveLink(
                _source_cv_key="recycle", _source_phase_key="liquid",
                _sink_cv_key="reactor", _sink_phase_key="liquid",
                Q_L_per_h=2.0, _label="bwd",
            ),
        ],
        system_solver=system_solver,
        label="custom_system_solver_demo",
    )

print("AsymmetricThreeStageSystemSolver defined")
"""),

    md("s2-run-md", """\
### Run it

Same two-CV recirculating model (`reactor` <-> `recycle`) under
`ExplicitEulerSystemSolver` and the custom asymmetric split.\
"""),

    code("s2-run", """\
dt_h, n_steps = 0.1, 10

sim_euler = build_two_cv_sim(ExplicitEulerSystemSolver())
sim_custom = build_two_cv_sim(AsymmetricThreeStageSystemSolver())

result_euler = sim_euler.run(tau_h=dt_h * n_steps, n_steps=n_steps)
result_custom = sim_custom.run(tau_h=dt_h * n_steps, n_steps=n_steps)

S_euler = result_euler.liquid_mol["reactor"]["S"][-1]
S_custom = result_custom.liquid_mol["reactor"]["S"][-1]

print(f"ExplicitEulerSystemSolver          reactor S = {S_euler:.5f} mol")
print(f"AsymmetricThreeStageSystemSolver    reactor S = {S_custom:.5f} mol")
print(f"Relative difference: {abs(S_euler - S_custom) / S_euler:.2e}")
"""),

    md("s2-takeaway", """\
Both conserve mass and converge to the same steady state; they differ
in splitting error at finite `dt_h` — exactly the same trade-off
`StrangSplittingSystemSolver` demonstrates against Euler, just with a
bespoke asymmetric weighting instead of a symmetric one.\
"""),

    md("summary", """\
## Summary

| Axis | Protocol | Shipped default | Custom solver here | What it changes |
|---|---|---|---|---|
| 1 — per-CV physics | `StepSolver.solve_step(cv, dt_h, t_h, external_source_terms)` | `SequentialAdvanceSolver` (speciation before feed) | `SpeciationAfterFeedStepSolver` | Same-step pH visibility of a feed |
| 2 — whole-system orchestration | `SystemSolver.advance_system(sim, dt_h, t_h)` | `ExplicitEulerSystemSolver` / `StrangSplittingSystemSolver` (symmetric split) | `AsymmetricThreeStageSystemSolver` | Inter-CV transport splitting error at finite `dt_h` |

Both protocols are one method each — a bespoke solver is exactly as
much code as the physics/orchestration it needs, no framework
subclassing required. See [0_README.ipynb](0_README.ipynb) for the
full landscape of shipped solvers on both axes.\
"""),
)

save_path_walkthrough = HERE / "01_writing_a_custom_solver.ipynb"
with open(save_path_walkthrough, "w", encoding="utf-8") as f:
    json.dump(walkthrough_nb, f, indent=1, ensure_ascii=False)
print(f"Written: {save_path_walkthrough}")
