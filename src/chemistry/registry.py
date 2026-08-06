"""Unified chemistry registry.

This module centralizes small pieces of "chemistry database" information used
across the codebase, to avoid duplicating mappings in multiple places.

Scope (intentional)
-------------------
This is **not** a full thermodynamic database. It only contains the minimal,
pragmatic mappings needed for:

* strong-ion inference (e.g., NaCl -> Na+ + Cl-)
* user-facing ion label normalization (e.g., "SO4^2-" -> "SO4--")
* mapping user-facing ions ("Na+") to engine keyword conventions ("CT_Na")

The goal is to keep these conventions consistent between:

* :mod:`fermenter.speciation.strong_ions`
* :class:`fermenter.chemistry.types.AqueousTotalsUser`
* future recipe/builders that translate lab recipes into totals

Notes
-----
The user requested that this registry **not** rely on BioSTEAM/thermosteam for
properties like molecular weights. Accordingly, this module contains no runtime
hooks to those libraries.
"""

from __future__ import annotations

from typing import Dict, List, Tuple


# -----------------------------------------------------------------------------
# Ion label normalization and mapping to engine keys
# -----------------------------------------------------------------------------

def normalize_ion_label(label: str) -> str:
    """Normalize common ion-label formats.

    Accepts user inputs like::

        "SO4--", "SO4^2-", "SO4-2", "Ca2+", "Ca++", "cl-".

    Returns a canonical-ish label used by :func:`ion_to_engine_key`.
    """
    s = str(label).strip().replace(" ", "")
    if not s:
        return s
    # Remove caret notation: SO4^2-
    s = s.replace("^", "")
    # Normalize charge formatting
    if s.endswith("-2"):
        s = s[:-2] + "--"
    if s.endswith("+2"):
        s = s[:-2] + "++"
    if s.endswith("-3"):
        s = s[:-2] + "---"
    if s.endswith("+3"):
        s = s[:-2] + "+++"
    # Light normalization of case for common ions
    if len(s) >= 2 and s[0].isalpha():
        s = s[0].upper() + s[1:]
    return s


# Common ions in typical fermentation electrolytes.
ION_TO_ENGINE_KEY: Dict[str, str] = {
    # Monovalent
    "Na+": "CT_Na",
    "K+": "CT_K",
    "Cl-": "CT_Cl",
    "NO3-": "CT_NO3",
    # Divalent
    "Ca++": "CT_Ca",
    "Mg++": "CT_Mg",
    "SO4--": "CT_SO4",
    # Trace
    "Zn++": "CT_Zn",
    "Mn++": "CT_Mn",
    "Co++": "CT_Co",
    "Mo7O24------": "CT_Mo7O24",
}


def ion_to_engine_key(ion: str) -> str:
    """Map a user-facing ion label (e.g., 'Na+', 'Cl-', 'Ca++') to engine kwargs."""
    s = normalize_ion_label(ion)
    synonyms = {
        "Ca2+": "Ca++",
        "Mg2+": "Mg++",
        "SO4^2-": "SO4--",
    }
    s = synonyms.get(s, s)
    if s in ION_TO_ENGINE_KEY:
        return ION_TO_ENGINE_KEY[s]
    # Fallback: strip charge markers and prepend CT_
    core = "".join(ch for ch in s if ch.isalpha() or ch.isdigit())
    if not core:
        raise ValueError(f"Unrecognized ion label: {ion!r}")
    return f"CT_{core}"


def map_user_ions_to_engine(user_ions_mol_L: Dict[str, float]) -> Dict[str, float]:
    """Convert user-facing ion labels to engine keyword conventions."""
    out: Dict[str, float] = {}
    for k, v in (user_ions_mol_L or {}).items():
        key = ion_to_engine_key(k)
        out[key] = float(out.get(key, 0.0)) + float(v)
    return out


# -----------------------------------------------------------------------------
# Salt dissociation map (neutral salts -> strong ion totals)
# -----------------------------------------------------------------------------

# Values are (engine_ion_key, stoich coefficient)
SALT_DISSOCIATION_MAP: Dict[str, List[Tuple[str, float]]] = {
    # Chlorides
    "NaCl": [("CT_Na", 1), ("CT_Cl", 1)],
    "KCl": [("CT_K", 1), ("CT_Cl", 1)],
    "NH4Cl": [("CT_Cl", 1)],
    "CaCl2": [("CT_Ca", 1), ("CT_Cl", 2)],
    "MgCl2": [("CT_Mg", 1), ("CT_Cl", 2)],

    # Nitrates
    "NaNO3": [("CT_Na", 1), ("CT_NO3", 1)],
    "KNO3": [("CT_K", 1), ("CT_NO3", 1)],
    "Ca(NO3)2": [("CT_Ca", 1), ("CT_NO3", 2)],
    "Mg(NO3)2": [("CT_Mg", 1), ("CT_NO3", 2)],

    # Sulfates
    "Na2SO4": [("CT_Na", 2), ("CT_SO4", 1)],
    "K2SO4": [("CT_K", 2), ("CT_SO4", 1)],
    "(NH4)2SO4": [("CT_SO4", 1)],
    "MgSO4": [("CT_Mg", 1), ("CT_SO4", 1)],
    "CaSO4": [("CT_Ca", 1), ("CT_SO4", 1)],

    # Trace metal sulfates
    "ZnSO4": [("CT_Zn", 1), ("CT_SO4", 1)],
    "MnSO4": [("CT_Mn", 1), ("CT_SO4", 1)],
    "CoSO4": [("CT_Co", 1), ("CT_SO4", 1)],

    # Molybdate salts (illustrative)
    "(NH4)6Mo7O24": [("CT_Mo7O24", 1)],

    # ---- BioSTEAM chemical-name aliases ----
    # BioSTEAM registers salts by common name rather than formula.
    # These entries ensure strong_ions_from_feed_molL picks them up.
    "AmmoniumSulfate": [("CT_SO4", 1)],          # = (NH4)2SO4
    "AmmoniumMolybdate": [("CT_Mo7O24", 1)],      # = (NH4)6Mo7O24
    "KH2PO4": [("CT_K", 1)],                      # K+ only; phosphate via CT_P
    "MnCl2": [("CT_Mn", 1), ("CT_Cl", 2)],
    "CoCl2": [("CT_Co", 1), ("CT_Cl", 2)],
}


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
                f"Unknown compound {name!r}. Add it to fermenter.chemistry.registry.COMPOUND_DB "
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
            f"Register it in fermenter.chemistry.registry.COMPOUND_DB or use "
            f"one of: {sorted(COMPOUND_DB.keys())}",
            UserWarning,
            stacklevel=2,
        )
        return False


def validate_compound_ids(names: Dict[str, str]) -> Dict[str, bool]:
    """Validate multiple compound IDs and warn for any that are unrecognised.

    Parameters
    ----------
    names : dict
        Mapping of ``{context_label: compound_id}`` pairs, e.g.
        ``{"pH acid": "H3PO4", "pH base": "KOH"}``.

    Returns
    -------
    dict
        Mapping of ``{context_label: is_valid}`` for each input.
    """
    return {ctx: validate_compound_id(cid, context=ctx) for ctx, cid in names.items()}

