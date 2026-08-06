"""chemistry_database.py — Q7 demo: ChemistryDatabase create → extend → import.

Demonstrates the user-facing ChemistryDatabase lifecycle:

1. **Import a stock database** — ``AD_BASIC`` provides CO₂/NH₄⁺ aqueous
   chemistry plus cross-phase CO₂ partition out of the box.
2. **Extend with custom species** — add a model-specific VFA (butyric acid)
   without touching the stock module.
3. **Override the ThermoFramework** — switch to activity-corrected Davies
   model for a higher-accuracy run.
4. **Inspect the result** — iterate reactions and species to confirm the
   composition.

The demo is self-contained and produces no simulation output; it exercises
the chemistry declaration API only.  Run with:

    python demos/model_api/chemistry/chemistry_database.py
"""
from __future__ import annotations

import sys
import os

# ── path setup (dev install fallback) ─────────────────────────────────────────
_SRC = os.path.join(os.path.dirname(__file__), "..", "..", "..", "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

# ── 1. Import a stock database ─────────────────────────────────────────────────
from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC
from PyOMES.chemistry.databases.aqueous import AQUEOUS_DEFAULT

print("=== 1. Stock databases ===")
print(f"AQUEOUS_DEFAULT species : {sorted(AQUEOUS_DEFAULT.species.keys())}")
aq_labels = [r.label for r in AQUEOUS_DEFAULT.reactions]
print(f"AQUEOUS_DEFAULT reactions: {aq_labels}")
print()
ad_labels = [r.label for r in AD_BASIC.reactions]
print(f"AD_BASIC reactions       : {ad_labels}")
print()

# ── 2. Extend with a custom species ───────────────────────────────────────────
from PyOMES.chemistry import Species

ButyricAcid = Species(
    id="ButyricAcid",
    atoms={"C": 4, "H": 8, "O": 2},
    charge=0,
    MW=88.106,
)
Butyrate = Species(
    id="Butyrate-",
    atoms={"C": 4, "H": 7, "O": 2},
    charge=-1,
    MW=87.098,
)

MY_DB = AD_BASIC.extend(
    species={
        "ButyricAcid": ButyricAcid,
        "Butyrate-":   Butyrate,
    },
)

print("=== 2. Extended database ===")
print(f"MY_DB species count: {len(MY_DB.species)}")
assert "ButyricAcid" in MY_DB.species, "ButyricAcid missing from extended DB"
assert "CO2" in MY_DB.species,         "CO2 unexpectedly dropped by extension"
print(f"MY_DB has ButyricAcid: {ButyricAcid.id!r}")
print(f"MY_DB has CO2:         {MY_DB.species['CO2'].id!r}")
print()

# ── 3. Override ThermoFramework ────────────────────────────────────────────────
import dataclasses
from PyOMES.thermo import ThermoFramework

activity_thermo = dataclasses.replace(
    AD_BASIC.thermo,
    activity_model="davies",
    use_activity=True,
)
MY_DB_ACTIVE = MY_DB.extend(thermo=activity_thermo)

print("=== 3. Activity-corrected database ===")
print(f"AD_BASIC.thermo.use_activity   : {AD_BASIC.thermo.use_activity}")
print(f"MY_DB_ACTIVE.thermo.use_activity: {MY_DB_ACTIVE.thermo.use_activity}")
assert MY_DB_ACTIVE.thermo.use_activity is True
assert AD_BASIC.thermo.use_activity is False  # original unchanged
print()

# ── 4. Inspect reaction composition ───────────────────────────────────────────
print("=== 4. Reaction composition of MY_DB_ACTIVE ===")
for rxn in MY_DB_ACTIVE.reactions:
    label = rxn.label or "(unlabelled)"
    phases = {e.phase for e in rxn.stoichiometry}
    sp_ids = [e.species.id for e in rxn.stoichiometry]
    print(f"  {label:20s}  phases={sorted(phases)}  species={sp_ids}")

print()

# ── 5. pKa temperature correction via ThermoFramework ─────────────────────────
print("=== 5. Temperature correction ===")
tf = MY_DB_ACTIVE.thermo
pKa_CO2_25 = 6.35
dH_CO2 = 7646.0  # J/mol (BSM2-canonical)
pKa_CO2_35 = tf.pKa_at_T(pKa_ref=pKa_CO2_25, dH_J_per_mol=dH_CO2, T_K=308.15)
print(f"CO2 pKa1 at 25 C: {pKa_CO2_25:.4f}")
print(f"CO2 pKa1 at 35 C: {pKa_CO2_35:.4f}  (Van 't Hoff correction)")
assert pKa_CO2_35 < pKa_CO2_25, "pKa should decrease with temperature for CO2"
print()

print("All assertions passed. ChemistryDatabase lifecycle demo complete.")
