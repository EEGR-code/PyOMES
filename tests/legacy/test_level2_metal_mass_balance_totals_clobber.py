# -*- coding: utf-8 -*-
"""
Unit test: Level 2 association coupling must preserve metal totals.

This test is designed to catch a subtle but serious bug:
- apply_level2_associations() writes back FREE metal into out["Ca++"], out["Mg++"], ...
- if the totals-consistent coupling loop calls apply_level2_associations() repeatedly
  without restoring totals, later iterations treat FREE as TOTAL and totals are lost.

We validate that the reconstructed total metal after solve matches the input total.
"""

from __future__ import annotations

from typing import Dict, Any

from speciation.engine import BisectionChemicalEquilibriumEngine
from speciation.associations_level2 import LEVEL2_ASSOCIATIONS


def _reconstruct_metal_total(sol: Dict[str, Any], metal_key: str) -> float:
    # Prefer explicit stored total if present
    tot_key = f"{metal_key}_Total"
    if tot_key in sol:
        return float(sol[tot_key])

    # Otherwise reconstruct from free + complexes
    total = float(sol.get(metal_key, 0.0))
    for r in LEVEL2_ASSOCIATIONS:
        if r.metal_free_key == metal_key:
            total += float(sol.get(r.product, 0.0))
    return float(total)



def test_level2_preserves_metal_totals_across_assoc_iterations():
    # Choose a recipe that forms carbonate complexes
    # Equivalent to: 0.02 M NaHCO3 + 0.01 M CaCl2 (so Ca complexes should form)
    eng = BisectionChemicalEquilibriumEngine(
        level=2,
        use_activity=True,
        activity_model="davies",
        T_C=25.0,
    )

    sol = eng.solve(
        CT_Na=0.020,
        CT_Cl=0.020,       # from CaCl2 (2 Cl per Ca)
        CT_Ca=0.010,
        CT_Mg=0.0,
        CT_Zn=0.0,
        CT_Mn=0.0,
        CT_Co=0.0,
        CT_Mo7O24=0.0,
        CT_TIC=0.020,
        CT_P=0.0,
        acid_totals={},

        # This is key: force >1 association-coupling iteration
        max_iter_assoc=6,
        assoc_tol=1e-12,

        # With PHREEQC constants this should be the default, but keep explicit:
        assoc_K_basis="activity",
        use_activity=True,
        activity_model="davies",
    )

    ca_total_recon = _reconstruct_metal_total(sol, "Ca++")

    # Reconstructed total should match the input within a tight tolerance
    assert abs(ca_total_recon - 0.010) < 1e-9, (
        f"Ca mass balance failed: reconstructed={ca_total_recon:.12g}, expected=0.010"
    )

    # Also sanity-check that at least one Ca-carbonate complex formed
    ca_hco3 = float(sol.get("CaHCO3+", 0.0))
    ca_co3 = float(sol.get("CaCO3(aq)", 0.0))
    assert (ca_hco3 + ca_co3) > 0.0, "Expected Ca carbonate complexes to form, but none were present."


if __name__ == "__main__":
    print("Running Level 2 metal totals mass-balance test...")
    test_level2_preserves_metal_totals_across_assoc_iterations()
    print("PASS ✅")
