#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Demo: dynamic FBA on a toy network, wired through BlackBoxReactionModel.

.. note::

   **Attribution.** The 7-reaction central-metabolism network in this
   demo is an **original pedagogical construction** — not a copy of
   any published model. No direct published counterpart exists for
   this specific 7-reaction topology: Orth, Thiele & Palsson (2010)
   Box 1 presents the S-matrix formalism and steady-state constraint
   equations, not a concrete toy reaction list. The *dFBA dynamic-
   coupling pattern* (Mahadevan et al 2002) and the *FBA framework*
   (Orth et al 2010) are cited references; the reaction network itself
   is an original teaching abstraction. Use this demo to understand
   the *protocol* (how to wire an FBA solver through
   :class:`~PyOMES.reactions.BlackBoxReactionModel`), not as a source
   of verified biochemical stoichiometry.

The dynamic-FBA coupling pattern — solve the LP each timestep with
uptake bounds derived from the current extracellular concentrations,
integrate the resulting fluxes into a mass balance, repeat — follows:

    Mahadevan, R., Edwards, J. S., & Doyle III, F. J. (2002). Dynamic
    flux balance analysis of diauxic growth in Escherichia coli.
    Biophysical Journal, 83(3), 1331-1340.

The general FBA framework (objective function, stoichiometric
matrix, LP formulation, pseudo-steady-state assumption on
intracellular metabolites) follows:

    Orth, J. D., Thiele, I., & Palsson, B. O. (2010). What is flux
    balance analysis? Nature Biotechnology, 28(3), 245-248.

The intracellular metabolites are pseudo-steady-state inside each LP
(S @ v = 0 every timestep). The extracellular state (glucose, O2,
biomass, acetate, CO2) evolves over time. This is what
:class:`~PyOMES.reactions.BlackBoxReactionModel` is designed for — a
static one-shot FBA solve wouldn't drive a simulation at all.

The toy reactions are::

    v_uptake_glc     : Glucose_ext  -> Glucose_int       [0, vmax_glc(env)]
    v_uptake_o2      : O2_ext       -> O2_int            [0, vmax_o2(env)]
    v_resp           : 1 Glc + 6 O2 -> 6 CO2 + 38 ATP    (oxidative)
    v_ferm           : 1 Glc        -> 2 Ac + 2 CO2 + 2 ATP (fermentative)
    v_growth         : (1/6) Glc + 20 ATP -> 1 Biomass   (synthesis)
    v_excrete_ac     : Acetate_int  -> Acetate_ext
    v_excrete_co2    : CO2_int      -> CO2_ext

Internal steady-state mass balances inside each LP::

    glucose_int : v_uptake_glc - v_resp - v_ferm - (1/6) v_growth = 0
    o2_int      : v_uptake_o2  - 6 v_resp                          = 0
    atp_int     : 38 v_resp + 2 v_ferm - 20 v_growth               = 0
    acetate_int : 2 v_ferm   - v_excrete_ac                        = 0
    co2_int     : 6 v_resp + 2 v_ferm - v_excrete_co2              = 0

Objective: max v_growth.

When O2 is plentiful the LP prefers respiration (high ATP yield, low
acetate excretion). When O2 is limited the LP switches to fermentation
(low ATP yield, acetate accumulates). The simulation below shows the
transition as O2 is consumed.

Run from the repo root after ``pip install -e .``::

    python demos/model_api/chemistry/fba/fba_toy.py
"""

# Bootstrap so the demo runs when src/ isn't on sys.path yet.
import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(_repo_root / "demos"))
import _bootstrap  # noqa: F401, E402

import numpy as np  # noqa: E402
from scipy.optimize import linprog  # noqa: E402

from PyOMES.core import Simulation  # noqa: E402
from PyOMES.core.control_volume import ControlVolume  # noqa: E402
from PyOMES.core.phases import LiquidPhase  # noqa: E402
from PyOMES.reactions import (  # noqa: E402
    BlackBoxReactionModel,
    FluxEntry,
    ReactionSystem,
)


# ── The FBA solver (the "black box") ──────────────────────────────────


class ToyFBASolver:
    """5-reaction-plus-2-exchanges central-metabolism toy.

    Exposes a single method, :meth:`solve`, taking the
    :class:`BlackBoxReactionModel` default ``inputs`` dict (CV species
    concentrations plus the environment scalars) and returning a
    ``{flux_id: mol/h}`` dict.
    """

    REACTIONS = [
        "v_uptake_glc",
        "v_uptake_o2",
        "v_resp",
        "v_ferm",
        "v_growth",
        "v_excrete_ac",
        "v_excrete_co2",
    ]
    INTERNAL_METABOLITES = [
        "glucose_int",
        "o2_int",
        "atp_int",
        "acetate_int",
        "co2_int",
    ]

    # Stoichiometric matrix S[i, j]: coefficient of metabolite i in
    # reaction j. Each column is one reaction (see module docstring).
    #                       upt_glc  upt_o2   resp   ferm   grow   ex_ac  ex_co2
    S = np.array(
        [
            [+1.0,    0.0,   -1.0,  -1.0,  -1 / 6,   0.0,    0.0],   # glc_int
            [ 0.0,   +1.0,   -6.0,   0.0,    0.0,    0.0,    0.0],   # o2_int
            [ 0.0,    0.0,  +38.0,  +2.0,  -20.0,    0.0,    0.0],   # atp_int
            [ 0.0,    0.0,    0.0,  +2.0,    0.0,   -1.0,    0.0],   # ac_int
            [ 0.0,    0.0,   +6.0,  +2.0,    0.0,    0.0,   -1.0],   # co2_int
        ]
    )

    def __init__(
        self,
        Vmax_glc: float = 10.0,
        Ks_glc: float = 0.5,
        Vmax_o2: float = 15.0,
        Ks_o2: float = 0.01,
    ):
        self.Vmax_glc = Vmax_glc
        self.Ks_glc = Ks_glc
        self.Vmax_o2 = Vmax_o2
        self.Ks_o2 = Ks_o2

    def solve(self, inputs: dict) -> dict:
        glc = max(inputs.get("Glucose", 0.0), 0.0)
        o2 = max(inputs.get("O2", 0.0), 0.0)
        V_L = float(inputs.get("_volume_L", 1.0))

        # Michaelis-Menten uptake bounds (mol / L / h). This is the
        # dynamic coupling — bounds depend on the current extracellular
        # concentrations rather than being fixed (Mahadevan et al 2002).
        vmax_glc = self.Vmax_glc * glc / (self.Ks_glc + glc + 1e-12)
        vmax_o2 = self.Vmax_o2 * o2 / (self.Ks_o2 + o2 + 1e-12)

        n_rxn = len(self.REACTIONS)
        c = np.zeros(n_rxn)
        c[self.REACTIONS.index("v_growth")] = -1.0  # maximise v_growth

        A_eq = self.S
        b_eq = np.zeros(self.S.shape[0])

        bounds = [(0.0, None)] * n_rxn
        bounds[self.REACTIONS.index("v_uptake_glc")] = (0.0, vmax_glc)
        bounds[self.REACTIONS.index("v_uptake_o2")] = (0.0, vmax_o2)

        result = linprog(
            c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs"
        )
        if not result.success:
            return {flux_id: 0.0 for flux_id in self.REACTIONS}

        # Fluxes are mol/L/h; convert to mol/h source terms.
        return {
            flux_id: float(result.x[i]) * V_L
            for i, flux_id in enumerate(self.REACTIONS)
        }


# ── Wire the FBA into a BlackBoxReactionModel ─────────────────────────


def make_fba_reaction_model() -> BlackBoxReactionModel:
    """Wrap :class:`ToyFBASolver` behind the BlackBox protocol.

    Each FBA flux is mapped to one CV species via :class:`FluxEntry`.
    ``sign=-1`` for uptake fluxes (positive FBA value → consumption);
    ``sign=+1`` for production / excretion fluxes. Atoms are declared
    for the runtime carbon-balance check.

    Internal-only reactions (``v_resp``, ``v_ferm``) are absent from
    the mapping — they're consumed by the steady-state mass balance
    inside the LP and have no CV-side presence.
    """
    flux_mapping = {
        "v_uptake_glc": FluxEntry(
            "Glucose", "liquid", {"C": 6, "H": 12, "O": 6}, sign=-1
        ),
        "v_uptake_o2": FluxEntry(
            "O2", "liquid", {"O": 2}, sign=-1
        ),
        "v_growth": FluxEntry(
            "Biomass", "liquid", {"C": 1, "H": 1.8, "O": 0.5}, sign=+1
        ),
        "v_excrete_ac": FluxEntry(
            "Acetate", "liquid", {"C": 2, "H": 4, "O": 2}, sign=+1
        ),
        "v_excrete_co2": FluxEntry(
            "CO2", "liquid", {"C": 1, "O": 2}, sign=+1
        ),
    }

    return BlackBoxReactionModel(
        external_model=ToyFBASolver(),
        flux_mapping=flux_mapping,
        balance_elements=("C",),  # carbon-only check (no redox cofactors)
        balance_atol=1e-6,
        on_imbalance="warn",
    )


# ── Build the Simulation ──────────────────────────────────────────────


def build() -> Simulation:
    """Build a single-CV simulation with the FBA reaction model.

    The CV has a single liquid phase containing glucose, dissolved O2,
    biomass, acetate and CO2 (the last two start at zero). The FBA
    solver drives every species change. No gas-liquid coupling: with
    no O2 supply, the run is intentionally O2-limited so the LP
    switches from respiration to fermentation as the run progresses.
    """
    liquid = LiquidPhase(
        n_mol={
            "Glucose": 0.050,    # 50 mmol  (50 mM in 1 L)
            "O2":      0.0002,   # 0.2 mmol (saturated 0.2 mM)
            "Biomass": 0.005,    # 5 mmol C-equivalent — seed inoculum
            "Acetate": 0.0,
            "CO2":     0.0,
        },
        V_L=1.0,
        T_K=310.15,
    )
    cv = ControlVolume(
        phases={"liquid": liquid},
        reaction_system=ReactionSystem([make_fba_reaction_model()]),
        label="fba_toy",
    )
    return Simulation(cvs={"main": cv}, label="fba_toy_demo")


# ── Run and report ────────────────────────────────────────────────────


TAU_H = 6.0
N_STEPS = 600


def main() -> None:
    sim = build()
    result = sim.run(tau_h=TAU_H, n_steps=N_STEPS)

    cv_key = "main"
    print(
        f"Demo: dFBA toy network (original 7-reaction central-metabolism "
        f"construction, Mahadevan et al 2002 dynamic coupling). "
        f"t = 0 -> {TAU_H} h, {N_STEPS} steps."
    )
    print("\nTrajectory (mol):")
    print(
        f"  {'t [h]':>6}  "
        f"{'Glucose':>9}  {'O2':>9}  "
        f"{'Biomass':>9}  {'Acetate':>9}  {'CO2':>9}"
    )
    t_h = result.t_h
    sample_idx = np.linspace(0, len(t_h) - 1, 10).astype(int)
    liq = result.liquid_mol[cv_key]
    for i in sample_idx:
        print(
            f"  {t_h[i]:>6.2f}  "
            f"{liq['Glucose'][i]:>9.5f}  {liq['O2'][i]:>9.6f}  "
            f"{liq['Biomass'][i]:>9.5f}  {liq['Acetate'][i]:>9.5f}  "
            f"{liq['CO2'][i]:>9.5f}"
        )

    print(f"\nRuntime: {result.runtime_s:.3f} s wall-clock")
    print(
        "Note: with no oxygen supply, the LP shifts from "
        "respiration to fermentation as O2 is depleted; "
        "acetate accumulates and biomass yield drops."
    )


if __name__ == "__main__":
    main()
