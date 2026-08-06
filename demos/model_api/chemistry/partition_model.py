#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Demo: PartitionModel — H₂S alpha correction and custom partition models.

This demo shows the three things the :class:`PartitionModel` protocol and
:class:`HenryPartition` dataclass enable:

1. **Inspect stock partition models** — access :data:`AD_BASIC`'s built-in
   H₂S, CO₂, CH₄, H₂, and NH₃ ``HenryPartition`` objects.

2. **Temperature dependence** — ``HenryPartition`` carries a van 't Hoff
   coefficient so kH varies with temperature, unlike the old pinned constants
   in the fermenter builder.

3. **H₂S alpha correction** — at pH ≈ pKa (7.0) half of dissolved sulfide
   is in the non-volatile HS⁻ form.  Without the alpha correction the
   link overestimates the gas-phase H₂S fraction by ~2×.  This demo
   quantifies that error across the AD operating pH range.

4. **Custom partition model** — extend a database with a user-declared
   ``HenryPartition`` for a species not in the stock databases.

Run from the repo root after ``pip install -e .``::

    python demos/model_api/chemistry/partition_model.py
"""

import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_repo_root / "demos"))
import _bootstrap  # noqa: F401, E402

import math  # noqa: E402

from PyOMES.chemistry import HenryPartition, PartitionModel  # noqa: E402
from PyOMES.chemistry.databases.anaerobic_digestion import AD_BASIC  # noqa: E402


# ── 1. Stock partition models in AD_BASIC ──────────────────────────────────────

print("=== 1. Partition models in AD_BASIC ===")
for sp_id, model in sorted(AD_BASIC.partition_models.items()):
    kH_25 = model._kH_mol_L_atm(298.15)
    print(f"  {sp_id:5s}  kH(25 C) = {kH_25:.3e} mol/(L.atm)  dlnH = {model.dlnH:.0f} K")
print()


# ── 2. Temperature dependence ──────────────────────────────────────────────────

print("=== 2. H2S temperature dependence ===")
h2s = AD_BASIC.partition_models["H2S"]
for T_C in (15, 25, 35, 45):
    T_K = T_C + 273.15
    kH = h2s._kH_mol_L_atm(T_K)
    print(f"  T = {T_C:2d} C   kH = {kH:.4f} mol/(L.atm)")
print()


# ── 3. H₂S alpha correction across pH ─────────────────────────────────────────
#
# At pH = pKa (7.0) half of total dissolved sulfide is H₂S (volatile) and
# half is HS⁻ (non-volatile).  alpha = [H₂S] / ([H₂S] + [HS⁻]).
#
# For a monoprotic acid:  alpha = 1 / (1 + 10^(pH - pKa))
#
# The gas-liquid link passes alpha to HenryPartition.beta() which divides kH
# by alpha, making the effective Henry constant larger (more liquid-favoured)
# at higher pH.  For CO₂ (pKa₁ = 6.35) the same correction applies but the
# magnitude is small because kH is already large.

PKA_H2S = 7.0
T_K = 308.15   # 35°C — typical AD operating temperature
V_liq = 1.0    # L (normalised)
V_gas = 0.25   # L (headspace fraction 20% of 1.25 L total)

print("=== 3. H2S alpha correction at 35 C ===")
print(f"  pKa(H2S) = {PKA_H2S},  V_liq = {V_liq} L,  V_gas = {V_gas} L")
print()
print(f"  {'pH':>5}  {'alpha':>7}  {'f_gas (alpha=1)':>16}  {'f_gas (corrected)':>18}  {'error':>8}")
print(f"  {'-'*5}  {'-'*7}  {'-'*16}  {'-'*18}  {'-'*8}")

for pH in (5.0, 6.0, 7.0, 7.5, 8.0, 9.0):
    alpha = 1.0 / (1.0 + 10 ** (pH - PKA_H2S))

    # Without correction: alpha fixed at 1.0
    beta_uncorrected = h2s.beta(V_liq, V_gas, T_K, alpha=1.0)
    f_gas_uncorrected = 1.0 - beta_uncorrected / (1.0 + beta_uncorrected)

    # With alpha correction
    beta_corrected = h2s.beta(V_liq, V_gas, T_K, alpha=alpha)
    f_gas_corrected = 1.0 - beta_corrected / (1.0 + beta_corrected)

    if f_gas_corrected > 0:
        error = (f_gas_uncorrected - f_gas_corrected) / f_gas_corrected
    else:
        error = float("inf")

    print(f"  {pH:5.1f}  {alpha:7.3f}  {f_gas_uncorrected:16.4f}  {f_gas_corrected:18.4f}  {error:+8.1%}")

print()
print("  Interpretation: at pH = pKa (7.0) the uncorrected model overestimates")
print("  the gas-phase H2S fraction by ~2x.  At pH 8 the error exceeds 4x.")
print()


# ── 4. Extend a database with a custom partition model ────────────────────────

print("=== 4. Custom partition model ===")

# Ethanol: moderately volatile, no ionisation (no alpha correction needed).
# Sander (2015) kH₀ ≈ 192 mol/(L·atm) at 25°C, dlnH ≈ 6600 K.
EtOH_PARTITION = HenryPartition(H_ref=1.9e0, dlnH=6600.0)

MY_DB = AD_BASIC.extend(partition_models={"Ethanol": EtOH_PARTITION})
print(f"  AD_BASIC partition_models : {sorted(AD_BASIC.partition_models)}")
print(f"  MY_DB   partition_models  : {sorted(MY_DB.partition_models)}")
print()

kH_etoh_25 = MY_DB.partition_models["Ethanol"]._kH_mol_L_atm(298.15)
kH_etoh_35 = MY_DB.partition_models["Ethanol"]._kH_mol_L_atm(308.15)
print(f"  Ethanol kH(25 C) = {kH_etoh_25:.1f} mol/(L.atm)")
print(f"  Ethanol kH(35 C) = {kH_etoh_35:.1f} mol/(L.atm)")
print()

# Protocol check: HenryPartition satisfies PartitionModel
assert callable(MY_DB.partition_models["Ethanol"].beta)
assert callable(MY_DB.partition_models["Ethanol"].equilibrium_a_moles)
print("  HenryPartition satisfies PartitionModel protocol: OK")
print()

print("All assertions passed. PartitionModel demo complete.")
