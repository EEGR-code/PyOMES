"""fermenter.units

Central place for **units**, **conversion helpers**, and shared physical constants.

Why this exists
--------------
The codebase historically contained repeated definitions of common constants
(e.g., the ideal gas constant in L·atm/(mol·K)) and mixed naming patterns such as
`molL`, `mol_L`, and `mol/L` across modules.

This module provides:

* A single, authoritative definition of shared constants.
* A small set of conversion helpers that are safe to call at *boundaries*
  (I/O, stream mapping, diagnostics), while keeping the inner ODE RHS numeric
  and fast.
* A documented naming convention for new variables.

Naming convention (recommended)
-------------------------------
Use a suffix that encodes the unit, e.g.:

* `T_K` (Kelvin)
* `P_atm` or `P_Pa`
* `V_L` or `V_m3`
* `n_mol` (moles)
* `C_mol_L` (concentration)
* `Q_L_min`, `Q_m3_s`, `rate_mol_h`

Existing public attributes are **not renamed** to preserve backward
compatibility, but new code should prefer the suffix style above.
"""

from __future__ import annotations

from typing import Final


# --- Physical constants (authoritative definitions) ---

# Ideal gas constant in L·atm/(mol·K). Kept in this unit because the headspace
# model uses (atm, L, K, mol) throughout.
R_L_ATM_PER_MOL_K: Final[float] = 0.082057366080960

# Ideal gas constant in J/(mol·K) (CODATA 2018 value). Used by van't Hoff
# temperature corrections (log_K(T)) and Clausius-Clapeyron partition models.
R_J_PER_MOL_K: Final[float] = 8.31446261815324

# Exact definitions for pressure/volume conversions.
PA_PER_ATM: Final[float] = 101_325.0
L_PER_M3: Final[float] = 1_000.0

MIN_PER_HR: Final[float] = 60.0
SEC_PER_MIN: Final[float] = 60.0
SEC_PER_HR: Final[float] = MIN_PER_HR * SEC_PER_MIN


# --- Simple conversion helpers ---

def atm_to_Pa(P_atm: float) -> float:
    """Convert pressure from atm to Pa."""
    return float(P_atm) * PA_PER_ATM


def Pa_to_atm(P_Pa: float) -> float:
    """Convert pressure from Pa to atm."""
    return float(P_Pa) / PA_PER_ATM


def L_to_m3(V_L: float) -> float:
    """Convert volume from L to m³."""
    return float(V_L) / L_PER_M3


def m3_to_L(V_m3: float) -> float:
    """Convert volume from m³ to L."""
    return float(V_m3) * L_PER_M3


def g_L_to_mol_L(c_g_L: float, MW_g_mol: float) -> float:
    """Convert concentration from g/L to mol/L."""
    mw = max(float(MW_g_mol), 1e-30)
    return float(c_g_L) / mw


def mol_L_to_g_L(c_mol_L: float, MW_g_mol: float) -> float:
    """Convert concentration from mol/L to g/L."""
    return float(c_mol_L) * float(MW_g_mol)
