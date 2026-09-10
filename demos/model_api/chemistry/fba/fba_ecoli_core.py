#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Demo: dynamic FBA loaded from a JSON data file (E. coli core scale).

.. note::

   **Attribution.** The ``ecoli_core.json`` shipped with this demo is
   an **original hand-crafted pedagogical construction** — not a
   faithful copy of any published metabolic reconstruction. It covers
   18 reactions: lumped glycolysis (GLYC), pyruvate dehydrogenase
   (PDH), lumped TCA cycle (TCA), electron-transport chain /
   oxidative phosphorylation (ETC), mixed-acid fermentation routes
   (lactate via LDH, formate + acetyl-CoA via PFL, ethanol via ADH,
   acetate via ACK), biomass synthesis (BIOMASS_SYN), ATP maintenance
   (ATP_MAINT), and 8 exchange reactions.

   **Scope and known limitations** — by design, for pedagogical
   clarity:

   - Glycolysis and TCA are each a single lumped reaction; intermediate
     metabolites (PEP, OAA, isocitrate, etc.) and the pentose-
     phosphate / anaplerotic / glyoxylate pathways are absent.
   - The biomass equation uses simplified stoichiometry (1 pyr +
     0.5 accoa + 5 ATP per biomass unit), not the multi-precursor
     BIOMASS_Ecoli_core_w_GAM equation from a real reconstruction.
   - ATP and NADH are cofactors only; no explicit proton, water, or
     phosphate tracking.
   - 18 reactions vs. 95 in the canonical BiGG ``e_coli_core`` model.

   Use this demo to understand the *loader pattern* (how to drive
   :class:`~PyOMES.reactions.BlackBoxReactionModel` from a structured
   data file), not as a source of verified *E. coli* stoichiometry.
   For a canonical 95-reaction model see the BiGG ``e_coli_core``
   (Orth, Fleming & Palsson 2011).

A bigger sibling to ``fba_toy.py``. The same dynamic-FBA pattern
(Mahadevan, Edwards & Doyle 2002 — solve the LP each timestep with
uptake bounds derived from the current extracellular concentrations)
applied to a larger network loaded from ``ecoli_core.json``.

References for the network the shipped subset is *inspired by*:

    Edwards, J. S., & Palsson, B. O. (2000). The Escherichia coli
    MG1655 in silico metabolic genotype: its definition,
    characteristics, and capabilities. PNAS, 97(10), 5528-5533.

    Orth, J. D., Fleming, R. M. T., & Palsson, B. O. (2011).
    Reconstruction and Use of Microbial Metabolic Networks: the Core
    Escherichia coli Metabolic Model as an Educational Guide. EcoSal
    Plus, 4(1).

Conventions follow COBRA:

- Exchange reactions have ``{met: -1}`` stoichiometry.
- Positive exchange flux = secretion. Negative = uptake.
- The CV-side rate for each exchange is just the LP flux value
  (no sign flip): if ``v_EX_glc_e = -10`` (uptake), the CV's
  glucose changes by -10 mol/h.

Run from the repo root after ``pip install -e .``::

    python demos/model_api/chemistry/fba/fba_ecoli_core.py
"""

# Bootstrap so the demo runs when PyOMES isn't installed yet.
import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(_repo_root / "demos"))
import _bootstrap  # noqa: F401, E402

import json  # noqa: E402

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


# ── The FBA solver, loaded from JSON ──────────────────────────────────


class JSONFBASolver:
    """Loads stoichiometry + bounds from a JSON file and solves an LP.

    The JSON schema is:

    .. code-block:: json

        {
          "metabolites": {met_id: {"name": ..., "atoms": {elem: count}}, ...},
          "reactions":   {rxn_id: {"name": ..., "stoichiometry": {met: coef},
                                   "lower_bound": float, "upper_bound": float}, ...},
          "objective":   {rxn_id: coefficient, ...},
          "cv_species_mapping": {
              exchange_rxn_id: {
                  "cv_species": str, "phase": str,
                  "Vmax": float | absent, "Ks": float | absent},
              ...
          }
        }

    Dynamic bounds are applied at every :meth:`solve` call: if a
    mapping entry has ``Vmax`` and ``Ks``, that exchange's lower
    bound is overwritten with ``-Vmax * c / (Ks + c)`` (Michaelis-
    Menten uptake) using the current CV concentration of the mapped
    species. All other reaction bounds stay at the JSON defaults.
    """

    def __init__(self, json_path: Path):
        data = json.loads(json_path.read_text(encoding="utf-8"))
        self.metabolites = list(data["metabolites"].keys())
        self.met_index = {m: i for i, m in enumerate(self.metabolites)}
        self.met_atoms = {
            m: dict(spec.get("atoms", {}))
            for m, spec in data["metabolites"].items()
        }

        self.reactions = list(data["reactions"].keys())
        self.rxn_index = {r: j for j, r in enumerate(self.reactions)}
        self.rxn_specs = data["reactions"]

        self.objective = dict(data["objective"])
        self.cv_mapping = dict(data["cv_species_mapping"])

        # Pre-build the stoichiometric matrix S (M × N).
        M, N = len(self.metabolites), len(self.reactions)
        self.S = np.zeros((M, N))
        for r, spec in self.rxn_specs.items():
            j = self.rxn_index[r]
            for m, coef in spec["stoichiometry"].items():
                self.S[self.met_index[m], j] = float(coef)

        # Base bounds vector (overwritten per-step for dynamic
        # uptake exchanges).
        self.base_lb = np.array(
            [self.rxn_specs[r]["lower_bound"] for r in self.reactions]
        )
        self.base_ub = np.array(
            [self.rxn_specs[r]["upper_bound"] for r in self.reactions]
        )

        # Objective coefficient vector — linprog minimises, so we
        # negate.
        self.c = np.zeros(N)
        for r, coef in self.objective.items():
            self.c[self.rxn_index[r]] = -float(coef)

    def solve(self, inputs: dict) -> dict:
        V_L = float(inputs.get("_volume_L", 1.0))

        lb = self.base_lb.copy()
        ub = self.base_ub.copy()

        # Dynamic uptake bounds (Michaelis-Menten).
        for rxn, m in self.cv_mapping.items():
            if "Vmax" not in m or "Ks" not in m:
                continue
            conc = max(inputs.get(m["cv_species"], 0.0), 0.0)
            v_uptake_max = m["Vmax"] * conc / (m["Ks"] + conc + 1e-12)
            # Override the uptake (negative) lower bound. Bound is
            # the most negative value the flux can take — closer to
            # zero = less uptake capability.
            lb[self.rxn_index[rxn]] = -v_uptake_max

        bounds = list(zip(lb, ub))

        result = linprog(
            self.c, A_eq=self.S, b_eq=np.zeros(self.S.shape[0]),
            bounds=bounds, method="highs",
        )
        if not result.success:
            return {rxn: 0.0 for rxn in self.reactions}

        # Fluxes are mol/L/h; convert to mol/h source terms.
        return {
            rxn: float(result.x[j]) * V_L
            for j, rxn in enumerate(self.reactions)
        }


# ── Wire the FBA into a BlackBoxReactionModel ─────────────────────────


def make_fba_reaction_model(json_path: Path) -> BlackBoxReactionModel:
    """Wrap :class:`JSONFBASolver` behind the BlackBox protocol.

    All exchange reactions get ``sign=+1`` (COBRA convention: the LP
    value is already in the CV's reference frame — negative for
    uptake, positive for secretion). Atoms come from the JSON
    metabolite declarations.
    """
    solver = JSONFBASolver(json_path)

    flux_mapping = {}
    for rxn, m in solver.cv_mapping.items():
        # Find which metabolite this exchange touches and grab its
        # atoms from the JSON.
        stoich = solver.rxn_specs[rxn]["stoichiometry"]
        if len(stoich) != 1:
            raise ValueError(
                f"Exchange reaction {rxn!r} must have exactly one "
                f"metabolite in its stoichiometry; got {list(stoich)}."
            )
        (met_id,) = stoich.keys()
        flux_mapping[rxn] = FluxEntry(
            species_id=m["cv_species"],
            phase=m["phase"],
            atoms=solver.met_atoms[met_id],
            sign=+1.0,
        )

    return BlackBoxReactionModel(
        external_model=solver,
        flux_mapping=flux_mapping,
        balance_elements=("C",),
        balance_atol=1e-6,
        on_imbalance="warn",
    )


# ── Build the Simulation ──────────────────────────────────────────────


JSON_PATH = Path(__file__).resolve().parent / "ecoli_core.json"


def build() -> Simulation:
    """Build a single-CV simulation around the JSON-loaded FBA model.

    The CV's liquid phase initial inventory must include every
    species the JSON's ``cv_species_mapping`` references. Their
    starting moles are the "extracellular conditions" the dFBA pattern
    evolves through time.
    """
    liquid = LiquidPhase(
        n_mol={
            "Glucose": 0.080,    # 80 mmol substrate
            "O2":      0.0002,   # 0.2 mmol — small DO buffer
            "Biomass": 0.005,    # 5 mmol seed inoculum
            "Acetate": 0.0,
            "Lactate": 0.0,
            "CO2":     0.0,
            "Formate": 0.0,
            "Ethanol": 0.0,
        },
        V_L=1.0,
        T_K=310.15,
    )
    cv = ControlVolume(
        phases={"liquid": liquid},
        reaction_system=ReactionSystem(
            [make_fba_reaction_model(JSON_PATH)]
        ),
        label="fba_ecoli_core",
    )
    return Simulation(cvs={"main": cv}, label="fba_ecoli_core_demo")


# ── Run and report ────────────────────────────────────────────────────


TAU_H = 6.0
N_STEPS = 600


def main() -> None:
    if not JSON_PATH.exists():
        raise FileNotFoundError(
            f"Expected {JSON_PATH} to exist next to this demo. The "
            f"file ships with PyOMES under demos/api/fba/."
        )

    sim = build()
    result = sim.run(tau_h=TAU_H, n_steps=N_STEPS)

    cv_key = "main"
    print(
        f"Demo: dFBA on a JSON-loaded central-metabolism subset "
        f"(E. coli core scale). t = 0 -> {TAU_H} h, {N_STEPS} steps."
    )
    print(f"Model file: {JSON_PATH}")
    print("\nTrajectory (mol):")
    print(
        f"  {'t [h]':>6}  "
        f"{'Glucose':>9}  {'O2':>9}  "
        f"{'Biomass':>9}  {'Acetate':>9}  {'Lactate':>9}  "
        f"{'Formate':>9}  {'Ethanol':>9}  {'CO2':>9}"
    )
    t_h = result.t_h
    sample_idx = np.linspace(0, len(t_h) - 1, 10).astype(int)
    liq = result.liquid_mol[cv_key]
    for i in sample_idx:
        print(
            f"  {t_h[i]:>6.2f}  "
            f"{liq['Glucose'][i]:>9.5f}  {liq['O2'][i]:>9.6f}  "
            f"{liq['Biomass'][i]:>9.5f}  {liq['Acetate'][i]:>9.5f}  "
            f"{liq['Lactate'][i]:>9.5f}  "
            f"{liq['Formate'][i]:>9.5f}  {liq['Ethanol'][i]:>9.5f}  "
            f"{liq['CO2'][i]:>9.5f}"
        )

    print(f"\nRuntime: {result.runtime_s:.3f} s wall-clock")
    print(
        "Note: aerobic growth (high ATP yield via ETC) runs until O2 "
        "is depleted, then the LP shifts to mixed-acid fermentation — "
        "PFL routes pyruvate to formate + acetyl-CoA (NADH-free), and "
        "ADH drains the glycolytic NADH via ethanol. Biomass growth "
        "continues at lower yield in the anaerobic phase."
    )


if __name__ == "__main__":
    main()
