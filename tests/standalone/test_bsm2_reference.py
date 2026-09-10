# -*- coding: utf-8 -*-
"""Golden-trajectory regression tests for the BSM2 reference model.

Locks current numerics of ``models/vlmodels/adm1/bsm2.py`` before
chemistry-unification refactors shift the speciation / reaction
declaration surface. Any drift caught here should be intentional —
update the sentinel values explicitly in the same commit that
introduces the numerical change, with a one-line note on why.

See ``docs/dev/implementation/shipped/BSM2_REFERENCE_TEST.md`` for the rationale
and shipping plan.

Approach
--------
- Build the canonical BSM2 CV via ``build_bsm2_reactions`` +
  ``build_bsm2_cv`` with default args.
- Seed the liquid phase with a documented "starter culture" of
  substrates and biomass (the builder leaves the liquid phase
  essentially empty).
- Run ``N_STEPS`` Euler-snapshot steps at ``dt = DT_H``.
- Assert hardcoded end-of-trajectory values for sentinel species,
  pH, and gas-phase mole counts (``TestBSM2Sentinels``), plus
  invariants over the full trajectory (``TestBSM2Invariants``).
"""
from __future__ import annotations

import math

import pytest


# ═══════════════════════════════════════════════════════════════════════
#  Trajectory configuration
# ═══════════════════════════════════════════════════════════════════════

DT_H = 0.01          # 36 s timestep
N_STEPS = 100        # 1 hour of simulated time
SNAPSHOT_EVERY = 10  # snapshot for invariant checks

# Liquid-phase initial concentrations (mol/L).  Documented "starter
# culture" — values are arbitrary but representative of a primed
# anaerobic digester; the test validates self-consistency over time
# rather than the seed itself.
INITIAL_LIQUID_CONC = {
    # Particulate composites and macromolecules
    "X_xc": 5.0e-3,
    "X_ch": 1.0e-3,
    "X_pr": 1.0e-3,
    "X_li": 5.0e-4,
    "X_I":  5.0e-3,
    # Soluble substrates
    "S_su":  2.0e-3,
    "S_aa":  2.0e-3,
    "S_fa":  5.0e-4,
    "S_va":  1.0e-3,
    "S_bu":  1.0e-3,
    "S_pro": 2.0e-3,
    "S_ac":  5.0e-3,
    "S_h2":  1.0e-7,
    "S_ch4": 1.0e-4,
    # Inorganics — chemistry-unification-3b: phase-agnostic Species IDs.
    # n_mol["CO2"] is molecular CO2(aq); the carbonate ladder
    # (CO2 + HCO3- + CO3--) implicitly tracks total inorganic C.
    # n_mol["NH3"] is molecular NH3; the engine redistributes to
    # NH4+ on the first solve at low pH. Total ammoniacal N is
    # implicit (n_mol["NH3"] + n_mol["NH4+"]).
    "CO2": 1.0e-2,
    "NH3":   5.0e-3,
    # Strong-ion lumps — state-unification C4: replaces the legacy
    # chem_env["strong_kwargs"]["CT_anion"] = 0.00521 path.
    # n_mol["S_an"] is the unnamed-anion lump (charge=-1, atoms={}).
    # The engine reads charge contribution from n_mol species
    # directly via _populate_strong_ions_from_phases.
    "S_an": 5.21e-3,
    # Biomass (one population per uptake pathway)
    "X_su":  1.0e-3,
    "X_aa":  1.0e-3,
    "X_fa":  1.0e-3,
    "X_c4":  1.0e-3,
    "X_pro": 1.0e-3,
    "X_ac":  1.0e-3,
    "X_h2":  1.0e-3,
}


# ═══════════════════════════════════════════════════════════════════════
#  Fixture
# ═══════════════════════════════════════════════════════════════════════

@pytest.fixture(scope="module")
def bsm2_trajectory():
    """Build the BSM2 CV, run a short trajectory, return state history.

    Module-scoped to amortise build + run cost across the assertion
    classes — every test reads from this dict rather than rebuilding.
    """
    from vlmodels.adm1.bsm2 import (
        build_bsm2_reactions, build_bsm2_cv,
    )
    from PyOMES.chemistry.thermo_params import ThermodynamicConfig

    rxn_set = build_bsm2_reactions(verbose=False)
    cv = build_bsm2_cv(rxn_set, thermo=ThermodynamicConfig.bsm2_default())

    # Seed liquid substrates — builder leaves the liquid phase empty.
    # state-unification C4e: INITIAL_LIQUID_CONC now also carries
    # the S_an strong-ion lump (formerly threaded through
    # make_bsm2_chem_env_fn's strong_kwargs["CT_anion"]=0.00521).
    V_liq = cv.phases["liquid"].V_L
    for sp, conc in INITIAL_LIQUID_CONC.items():
        cv.phases["liquid"].n_mol[sp] = conc * V_liq

    history = []

    def snapshot(t_h, pH):
        history.append({
            "t_h": t_h,
            "liquid_n_mol": dict(cv.phases["liquid"].n_mol),
            "gas_n_mol": dict(cv.phases["gas"].n_mol),
            "liquid_V_L": cv.phases["liquid"].V_L,
            "gas_V_L": cv.phases["gas"].V_L,
            "pH": pH,
        })

    snapshot(0.0, None)

    last_pH = None
    for step in range(N_STEPS):
        t_h = step * DT_H
        result = cv.advance(DT_H, t_h)

        # state-unification C4: pH lives on phase.pH (derived from
        # n_mol["H+"]), not on result.properties. The chem_env_fn
        # path is vestigial post-C4 (strong ions live in n_mol).
        try:
            last_pH = cv.phases["liquid"].pH
        except ValueError:
            last_pH = None

        if (step + 1) % SNAPSHOT_EVERY == 0:
            snapshot((step + 1) * DT_H, last_pH)

    return {
        "cv": cv,
        "history": history,
        "final_pH": last_pH,
    }


# ═══════════════════════════════════════════════════════════════════════
#  Golden sentinel values
# ═══════════════════════════════════════════════════════════════════════
#
# These values were captured from the current numerics of bsm2.py on
# the bsm2-reference-test branch.  Update intentionally — every diff
# to a sentinel should accompany a commit that explains why the
# numerics shifted.

# End-of-trajectory liquid concentrations (mol/L).
#
# IMPORTANT — these are NOT pyadm1 reference predictions. They are a
# self-consistency snapshot of PyOMES's own BSM2 numerics, captured
# when the test was added (commit 00df677, 2026-05-12) and re-baselined
# only when chemistry refactors intentionally change the math. See
# docs/dev/implementation/shipped/BSM2_REFERENCE_TEST.md for the original framing
# ("Out of scope: Literature validation. This test detects drift from
# *current numerics*, not correctness against Rosén & Jeppsson 2006").
# Cross-validation against pyadm1 or other reference implementations
# would live in a separate test file (e.g. tests/integration/
# test_bsm2_vs_pyadm1.py) with its own reference data files; it
# does not exist yet in this repo.
#
# Re-baselined in chemistry-unification-2 (commit landing the
# gas-liquid link migration to PropertyResult.alphas). The link
# previously read its own SpeciationCorrection.f_molecular(pH, T_K)
# inline — which paired the engine's activity-based pH with a
# *thermodynamic* Ka (no gamma correction on Ka itself). It now
# reads alphas[mol_key] from the speciation engine, which uses
# effective Ka (Ka_eff = Ka / (gamma_H * gamma_HCO3-)). At BSM2's
# ionic strength (~0.1 mol/L, Davies activity), these formulae
# differ by ~1e-5 to 1e-7 relative for the rates/concentrations
# below. This is Leak 1 closing — link and engine now agree on
# the same molecular fraction.
#
# Re-baselined in chemistry-unification-3 (2026-05-15): the
# canonical-naming dedup in _compute_species_eq unmasks a latent
# NH4+/ionic-strength bug introduced by Phase 1. Phase 1's
# migration to BisectionChemicalEquilibriumEngine.from_reactions sets
# eq_def.name = "NH3" for the NH4+/NH3 cation_acid (via
# total_id="NH3" on the reaction); the pre-Phase-1
# bsm2_default() factory used eq_def.name = "NH4", which hit a
# hardcoded `if eq_def.name == "NH4"` branch in
# _compute_species_eq that emitted canonical NH4+/NH3 names.
# The "NH3"-named path missed that branch and emitted only
# NH3_BH+/NH3_B, so ionic_strength_from_speciation's z-dict
# lookup for "NH4+" silently returned 0. BSM2 was undercounting
# I by the NH4+ contribution (~2.5e-3 mol/L, ~25% of total I).
# Phase 3's _CANONICAL_NAMES table aliases "NH3" → (NH4+, NH3)
# so both paths emit NH4+ uniformly. Net numerical effect:
# higher I → smaller Davies γ → smaller effective Ka → less
# dissociation → higher pH (3.332 → 3.340). Largest species
# drift: S_h2 (~3 %); other species 1e-7 to 1e-5 relative.
#
# Re-baselined 2026-05-20 during integrator-removal: the Integrator
# Protocol + RK4Integrator + EulerIntegrator were deleted; the
# sequential cv.advance() body's reaction sub-step collapsed from
# RK4 (4 RHS evaluations per dt_h) to a single forward Euler
# evaluation. Largest drift: S_h2 (~6e-3 relative, both liquid and
# gas) — H2 sits at ~1e-6 mol/L with fast Monod uptake, so per-step
# rate change vs. concentration is highest there, exactly the case
# INTEGRATOR_REMOVAL.md flagged ("Monod kinetics with substantial
# substrate depletion per step will see the largest drift"). Other
# species drift 1e-7 to 5e-6 relative. pH shifts by 1.6e-7 relative
# (3.340171 → 3.340171). Qualitative behaviour unchanged.
#
# Re-baselined 2026-06-02 during chemistry-unification-3b (C5: CO₂
# rename). Phase-agnostic Species-ID convention: n_mol["CO2"] replaces
# n_mol["CO2aq"] — same numerical value, different key. The carbonate
# ladder is now (CO2, HCO3-, CO3--); _CANONICAL_NAMES deleted;
# speciation engine emits via Species.id. Zero numerical drift expected
# (pure structural key rename). Key "CO2aq" → "CO2" throughout.
#
# Previous baselines for reference:
# Re-baselined 2026-05-21 during state-unification C3 (Phase
# unification, B2 removal). Adopting canonical naming for n_mol:
# molecular species live under their molecular IDs (CO2aq, HCO3-,
# CO3--, NH3, NH4+, H+, OH-) with no separate total tracker; total
# inorganic C / N is implicit (sum over the ladder). BSM2 migrated
# from n_mol["CO2"] (total) to n_mol["CO2aq"] (molecular) for the
# carbonate side; NH3 stayed as the canonical molecular name (the
# legacy BSM2 convention was already aligned with engine output).
# The speciation engine now writes derived species directly to
# phase.n_mol via _refresh_derived; KineticGasLiquidLink.compute_flow
# reads liquid totals by summing the canonical ladder via
# _CANONICAL_NAMES. Phase.pH derives from n_mol["H+"]/V_L (raises
# when H+ absent, replacing the old None-from-speciation-dict
# behaviour). Largest drift: S_h2 (1e-7 → 1.7e-8 mol/L, ~factor
# 57 — both H2 levels are well below the simulation's effective
# resolution, so this represents redistribution within numerical
# noise rather than physical change). Gas CO2 +17 % (3.67 → 4.31
# mol) reflecting the engine writeback now driving alpha lookups
# inline from n_mol ratios rather than via PropertyResult.alphas.
# Other species drift 1e-9 to 1.5e-2 relative. pH shifts by 4.4e-4
# absolute (3.340171 → 3.338685). Qualitative behaviour unchanged.
# Note the sentinel keys themselves migrate: "CO2"/"NH3" (formerly
# totals) split into canonical molecular and ladder-summed totals.
# Re-baselined 2026-05-21 during state-unification C4 (property
# layer collapse). cv.advance signature collapsed from
# (dt_h, chem_env) to (dt_h, t_h); strong-ion totals migrated into
# phase.n_mol as Species (n_mol["S_an"] = 5.21e-3 * V_L lumps
# the unnamed anion, charge=-1, atoms={}); the property_solvers
# iteration is gone from advance() -- speciation runs through
# cv.reaction_system.engine.solve(phases=...) directly. The link's
# alpha computation moved inline (PropertyResult.alphas channel
# removed). Drift from C3 baseline: pH 3.338685 -> 3.303580
# (-0.035 absolute, ~1.1% relative); other species 1e-9 to 1e-3
# relative. The shift comes from the engine's solve path now
# flowing through cv.advance directly rather than via the
# SpeciationPropertySolver wrapper; subtle differences in monitor
# wiring + alpha-lookup timing account for the residual.
# Qualitative behaviour unchanged.
SENTINEL_FINAL_LIQUID_CONC = {
    "S_ac":   0.00500000039477684,
    "S_pro":  0.0019999999947897394,
    "S_h2":   1.7012434526578658e-08,
    "S_ch4":  2.4113464446918378e-05,
    "CO2":    0.00850952869589902,
    "HCO3-":  8.468001007344533e-06,
    "NH3":    1.1190790929182218e-08,
    "NH4+":   0.004999988933569145,
}

# End-of-trajectory gas-phase mole counts (mol).
# Re-baselined alongside the liquid concentrations above.
SENTINEL_FINAL_GAS_MOL = {
    "S_ch4": 0.24305015591967846,
    "S_h2":  0.00026888028640737274,
    "CO2":   3.669917042504009,
}

# End-of-trajectory pH.  Sits well below typical AD operating range
# (6.5–8.0) because the starter inoculum carries unbuffered VFAs and
# BSM2's default ``CT_cation`` is 0.  Realistic operational pH is not
# the point of this test — drift detection is.
#
# Re-baselined 2026-05-13 during chemistry-unification-1: the previous
# values were captured with a hidden one-step pH lag from the deleted
# ``make_bsm2_callback`` (which wrote H⁺ to ``n_mol["H+"]`` at the end
# of each step, so rate functions saw step k-1's pH on step k). The
# stateless snapshot model removes this lag; rate functions now read
# the current step's pH from ``env.pH``. Largest species-level shift
# was S_ch4 (~17 %); qualitative behavior unchanged.
#
# Re-baselined 2026-05-15 during chemistry-unification-3: the
# NH4+/ionic-strength bug fix (see comment above SENTINEL_FINAL_LIQUID_CONC)
# pushes pH from 3.3316 to 3.3402 (Δ +8.6e-3 absolute, +2.6e-3
# relative). Higher I → smaller γ → smaller effective Ka → less
# dissociation → higher pH.
#
# Re-baselined 2026-05-20 during integrator-removal: RK4 → single
# Euler shifts pH by 5.4e-7 absolute (1.6e-7 relative). Tiny — the
# reaction sub-step's per-step trajectory difference washes out
# almost completely in the equilibrium speciation re-solve at
# end-of-step.
#
# Re-baselined 2026-05-21 during state-unification C3: canonical
# n_mol naming + engine writeback shifts pH by 1.5e-3 absolute
# (3.340171 → 3.338685). Drift driven by the engine now consuming
# carbonate-ladder totals by summing canonical species in n_mol
# (rather than reading a single conflated key), plus the link's
# inline n_liq summing across the ladder. The new computation is
# the design-intended PHREEQC-style flow; the old value was an
# artefact of the BSM2 conflation. Qualitative behaviour unchanged.
#
# Re-baselined 2026-06-02 during chemistry-unification-3b C5+C7: phase-agnostic
# CO₂ species ID. Previously, gas-liquid CO₂ transfer applied flux to
# n_liq["CO2"] (ghost key) while the speciation engine read n_liq["CO2aq"].
# The two keys were decoupled — the transfer loop was computing the right
# driving force from n_liq["CO2aq"] via _liquid_total, but writing the
# resulting flux to the ghost "CO2" key. This meant dissolved CO₂ never
# received the gas-transfer contribution in the speciation read path.
# After C5, both keys are unified: "CO2" is the dissolved form, the
# transfer writes to n_liq["CO2"], and speciation reads n_liq["CO2"].
# The coupled feedback loop is now correct. Liquid CO₂ shifts ~15 %
# (9.990e-3 → 8.516e-3 mol/L), gas CO₂ ~15 % (4.308 → 3.673 mol),
# S_h2 and S_ch4 ~1.7 %. pH: 3.3036 → 3.3044 (Δ +8e-4 absolute).
#
# Re-baselined 2026-07-02 during equilibrium-constraint-unification CP1:
# consolidated the J/(mol·K) gas constant into a single authoritative
# PyOMES.units.R_J_PER_MOL_K (CODATA 8.31446261815324). Previously
# src/chemistry/thermo_params.py and src/thermo/framework.py used a
# rounded 8.31446 for Van 't Hoff Ka(T) corrections while
# src/reactions/equilibrium.py and src/speciation/nr_tableau.py already
# used the precise value — the two R's didn't cancel on round-trip,
# producing a ~1e-9-level relative error in every temperature-corrected
# Ka. Species closest to this test's RTOL_SENTINEL=1e-9 boundary shift:
# liquid S_h2 (1.70124345582e-8 → 1.70124345266e-8), liquid HCO3-
# (8.46800127e-6 → 8.46800101e-6), liquid NH3 (1.11907933e-8 →
# 1.11907909e-8), gas S_h2 (2.688802867e-4 → 2.688802864e-4). All other
# sentinels (including pH) unchanged within tolerance. Qualitative
# behaviour unchanged — this is the constants becoming self-consistent,
# not a model change.
SENTINEL_FINAL_PH = 3.3044056113459916  # re-baselined C7 (writeback "CO2aq"→"CO2")

RTOL_SENTINEL = 1e-9


# ═══════════════════════════════════════════════════════════════════════
#  Sentinels — fail loudly when numerics shift
# ═══════════════════════════════════════════════════════════════════════

class TestBSM2Sentinels:
    """Hardcoded end-of-trajectory values; tight tolerance."""

    def test_final_liquid_concentrations(self, bsm2_trajectory):
        final = bsm2_trajectory["history"][-1]
        V_liq = final["liquid_V_L"]
        for sp, expected in SENTINEL_FINAL_LIQUID_CONC.items():
            assert expected is not None, \
                f"sentinel for liquid {sp} not yet captured"
            actual = final["liquid_n_mol"].get(sp, 0.0) / V_liq
            assert math.isclose(actual, expected, rel_tol=RTOL_SENTINEL), (
                f"liquid {sp}: expected {expected!r}, got {actual!r} "
                f"(rel diff {abs(actual - expected) / max(abs(expected), 1e-30):.3e})"
            )

    def test_final_gas_mol(self, bsm2_trajectory):
        final = bsm2_trajectory["history"][-1]
        for sp, expected in SENTINEL_FINAL_GAS_MOL.items():
            assert expected is not None, \
                f"sentinel for gas {sp} not yet captured"
            actual = final["gas_n_mol"].get(sp, 0.0)
            assert math.isclose(actual, expected, rel_tol=RTOL_SENTINEL), (
                f"gas {sp}: expected {expected!r}, got {actual!r}"
            )

    def test_final_pH(self, bsm2_trajectory):
        assert SENTINEL_FINAL_PH is not None, \
            "sentinel for pH not yet captured"
        pH = bsm2_trajectory["final_pH"]
        assert pH is not None, "final pH not populated"
        assert math.isclose(pH, SENTINEL_FINAL_PH, rel_tol=RTOL_SENTINEL), (
            f"pH: expected {SENTINEL_FINAL_PH!r}, got {pH!r}"
        )


# ═══════════════════════════════════════════════════════════════════════
#  Invariants — properties any correct BSM2 trajectory should satisfy
# ═══════════════════════════════════════════════════════════════════════

class TestBSM2Invariants:
    """Properties checked over the full trajectory, not just endpoint."""

    def test_pH_physically_valid(self, bsm2_trajectory):
        # Loose range — the inoculum is not operationally realistic
        # (see SENTINEL_FINAL_PH comment).  This invariant only catches
        # gross failure of the speciation solver.
        for snap in bsm2_trajectory["history"]:
            pH = snap["pH"]
            if pH is None:
                continue
            assert 1.0 < pH < 13.0, (
                f"pH {pH!r} outside physically-valid range at t={snap['t_h']:.3f}h"
            )

    def test_concentrations_non_negative(self, bsm2_trajectory):
        for snap in bsm2_trajectory["history"]:
            for sp, n in snap["liquid_n_mol"].items():
                assert n >= -1e-12, (
                    f"liquid {sp} negative ({n!r}) at t={snap['t_h']:.3f}h"
                )
            for sp, n in snap["gas_n_mol"].items():
                assert n >= -1e-12, (
                    f"gas {sp} negative ({n!r}) at t={snap['t_h']:.3f}h"
                )

    def test_speciation_populated_at_endpoint(self, bsm2_trajectory):
        assert bsm2_trajectory["final_pH"] is not None, (
            "speciation result missing at end of trajectory — "
            "indicates SpeciationPropertySolver did not run"
        )
