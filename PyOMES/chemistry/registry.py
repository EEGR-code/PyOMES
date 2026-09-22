"""Compound registry: name -> molar mass, ion composition, and other
recipe-relevant metadata.

Scope (intentional)
-------------------
This is **not** a full thermodynamic database. It exists to resolve a
user-provided compound name (e.g. ``"KOH"``, ``"MgSO4·7H2O"``) to its
molar mass and other dosing-relevant properties.

Notes
-----
The user requested that this registry **not** rely on BioSTEAM/thermosteam for
properties like molecular weights. Accordingly, this module contains no runtime
hooks to those libraries.
"""

from __future__ import annotations

from typing import Dict


# -----------------------------------------------------------------------------
# Minimal compound registry (recipe parsing)
# -----------------------------------------------------------------------------

# This registry is intentionally small. It exists to support "recipe-style"
# solution definitions (e.g., grams of Na2SO4, mL of 1 M HCl) without depending
# on BioSTEAM/thermosteam for molecular weights.

# Compound records may include:
# - MW_g_mol: molecular weight (g/mol)
# - ions: mapping of user-facing ion labels -> stoichiometry per mole compound
# - TIC_mol_per_mol: mol TIC added per mol compound (e.g. carbonates)
# - P_tot_mol_per_mol: mol total phosphate per mol compound
# - TAN_mol_per_mol: mol total ammonia nitrogen per mol compound
# - alk_eq_per_mol: alkalinity equivalents per mol compound (+ base, - acid)

COMPOUND_DB: Dict[str, dict] = {
    # Strong electrolytes / common salts
    "NaCl": {"MW_g_mol": 58.44277, "ions": {"Na+": 1, "Cl-": 1}},
    "KCl": {"MW_g_mol": 74.5513, "ions": {"K+": 1, "Cl-": 1}},
    "Na2SO4": {"MW_g_mol": 142.04214, "ions": {"Na+": 2, "SO4--": 1}},
    "MgSO4": {"MW_g_mol": 120.366, "ions": {"Mg++": 1, "SO4--": 1}},
    "MgSO4·7H2O": {"MW_g_mol": 246.469, "ions": {"Mg++": 1, "SO4--": 1}},
    "CaCl2": {"MW_g_mol": 110.984, "ions": {"Ca++": 1, "Cl-": 2}},
    "CaCl2·2H2O": {"MW_g_mol": 147.014, "ions": {"Ca++": 1, "Cl-": 2}},

    # Strong acids/bases (alkalinity bookkeeping)
    "HCl": {"MW_g_mol": 36.46094, "ions": {"Cl-": 1}, "alk_eq_per_mol": -1.0},
    "NaOH": {"MW_g_mol": 39.997, "ions": {"Na+": 1}, "alk_eq_per_mol": +1.0},
    "KOH": {"MW_g_mol": 56.1056, "ions": {"K+": 1}, "alk_eq_per_mol": +1.0},

"H2SO4": {"MW_g_mol": 98.07848, "ions": {"SO4--": 1}, "alk_eq_per_mol": -2.0},
"H3PO4": {"MW_g_mol": 97.994, "P_tot_mol_per_mol": 1.0, "alk_eq_per_mol": 0.0},
# Weak acids (represented via acid_totals + acid_pKas in speciation)
"AceticAcid": {"MW_g_mol": 60.052, "acid_system": "AceticAcid"},
"CitricAcid": {"MW_g_mol": 192.124, "acid_system": "CitricAcid"},


    # Carbonate system salts (contribute TIC and alkalinity)
    "NaHCO3": {"MW_g_mol": 84.0066, "ions": {"Na+": 1}, "TIC_mol_per_mol": 1.0, "alk_eq_per_mol": +1.0},
    "Na2CO3": {"MW_g_mol": 105.9888, "ions": {"Na+": 2}, "TIC_mol_per_mol": 1.0, "alk_eq_per_mol": +2.0},
    "KHCO3": {"MW_g_mol": 100.115, "ions": {"K+": 1}, "TIC_mol_per_mol": 1.0, "alk_eq_per_mol": +1.0},
    "K2CO3": {"MW_g_mol": 138.205, "ions": {"K+": 2}, "TIC_mol_per_mol": 1.0, "alk_eq_per_mol": +2.0},

    # Phosphate system salts (contribute phosphate totals)
    "KH2PO4": {"MW_g_mol": 136.0855, "ions": {"K+": 1}, "P_tot_mol_per_mol": 1.0, "alk_eq_per_mol": 0.0},
    "Na2HPO4": {"MW_g_mol": 141.957, "ions": {"Na+": 2}, "P_tot_mol_per_mol": 1.0, "alk_eq_per_mol": +1.0},
    "NaH2PO4": {"MW_g_mol": 119.977, "ions": {"Na+": 1}, "P_tot_mol_per_mol": 1.0, "alk_eq_per_mol": 0.0},

    # Ammonium salts (TAN contribution)
    "NH4Cl": {"MW_g_mol": 53.491, "ions": {"Cl-": 1}, "TAN_mol_per_mol": 1.0, "alk_eq_per_mol": 0.0},
    "(NH4)2SO4": {"MW_g_mol": 132.139, "ions": {"SO4--": 1}, "TAN_mol_per_mol": 2.0, "alk_eq_per_mol": 0.0},
}

_COMPOUND_ALIASES: Dict[str, str] = {
    # Hydrate spellings
    "MgSO4.7H2O": "MgSO4·7H2O",
    "MgSO4*7H2O": "MgSO4·7H2O",
    "CaCl2.2H2O": "CaCl2·2H2O",
    "CaCl2*2H2O": "CaCl2·2H2O",
    # Common names
    "Sodium chloride": "NaCl",
    "Potassium chloride": "KCl",
    "Sodium sulfate": "Na2SO4",
    "Potassium hydroxide": "KOH",
    "Sodium hydroxide": "NaOH",
    "Hydrochloric acid": "HCl",
"Sulfuric acid": "H2SO4",
"H2SO4": "H2SO4",
"Phosphoric acid": "H3PO4",
"H3PO4": "H3PO4",
"Acetic acid": "AceticAcid",
"acetic acid": "AceticAcid",
"Citric acid": "CitricAcid",
"citric acid": "CitricAcid",
"Sodium bicarbonate": "NaHCO3",
"sodium bicarbonate": "NaHCO3",

}


def resolve_compound(name: str) -> dict:
    """Resolve a user-provided compound name to a COMPOUND_DB record."""
    if not isinstance(name, str):
        raise TypeError("Compound name must be a string")
    key = name.strip()
    key = _COMPOUND_ALIASES.get(key, key)
    if key not in COMPOUND_DB:
        key2 = key.replace(" ", "")
        key2 = _COMPOUND_ALIASES.get(key2, key2)
        if key2 in COMPOUND_DB:
            key = key2
        else:
            raise KeyError(
                f"Unknown compound {name!r}. Add it to PyOMES.chemistry.registry.COMPOUND_DB "
                f"or use one of: {sorted(COMPOUND_DB.keys())}"
            )
    return COMPOUND_DB[key]


def validate_compound_id(name: str, *, context: str = "") -> bool:
    """Check whether a compound ID is known in the chemistry registry.

    Parameters
    ----------
    name : str
        Compound name to check (e.g. ``"KOH"``, ``"H3PO4"``).
    context : str, optional
        Human-readable context for the warning message (e.g. ``"PHController acid"``).

    Returns
    -------
    bool
        True if the compound is recognised, False otherwise.

    Notes
    -----
    When the compound is not found, a ``UserWarning`` is issued listing the
    available compounds. This allows downstream code to decide whether to
    raise an error or continue with degraded functionality.
    """
    try:
        resolve_compound(name)
        return True
    except (KeyError, TypeError):
        import warnings
        ctx = f" ({context})" if context else ""
        warnings.warn(
            f"[Chemistry] Compound {name!r}{ctx} is not in COMPOUND_DB. "
            f"Dosing with this compound will have no effect on speciation. "
            f"Register it in PyOMES.chemistry.registry.COMPOUND_DB or use "
            f"one of: {sorted(COMPOUND_DB.keys())}",
            UserWarning,
            stacklevel=2,
        )
        return False

