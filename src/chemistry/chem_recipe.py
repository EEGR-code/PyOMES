# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 13:56:11 2026

@author: k2473520
"""

# scripts/chem_recipe.py

from dataclasses import dataclass
from typing import Dict, Tuple, Optional


@dataclass(frozen=True)
class ChemSpec:
    """Speciation recipe for a dosing chemical.

    Maps one mole of a chemical (e.g. NaOH, HCl, Na₂CO₃) to its
    contributions to the speciation system: strong ion totals, acid
    systems, inorganic carbon, phosphorus, and ammonia.
    """
    # One of: strong_ions, acid_system, TIC, P, NH_T
    # strong_ions: contributions to CT_* totals
    strong_ions: Dict[str, float]  # e.g. {"CT_Na": +1, "CT_Cl": +1} per mole of chemical
    acid_system: Optional[Tuple[str, float]] = None  # (acid_name, stoich) e.g. ("AceticAcid", 1)
    TIC: float = 0.0   # mol TIC per mol chemical
    P: float = 0.0     # mol P per mol chemical
    NH_T: float = 0.0  # mol total ammonia per mol chemical
    molar_mass_g_mol: Optional[float] = None  # needed if you want g/L entry


# --- Registry: add what you need over time ---
CHEM_DB: Dict[str, ChemSpec] = {
    # Strong electrolytes
    "NaCl":   ChemSpec(strong_ions={"CT_Na": 1.0, "CT_Cl": 1.0}, molar_mass_g_mol=58.44277),
    "KCl":    ChemSpec(strong_ions={"CT_K":  1.0, "CT_Cl": 1.0}, molar_mass_g_mol=74.5513),
    "CaCl2":  ChemSpec(strong_ions={"CT_Ca": 1.0, "CT_Cl": 2.0}, molar_mass_g_mol=110.98),
    "MgSO4":  ChemSpec(strong_ions={"CT_Mg": 1.0, "CT_SO4": 1.0}, molar_mass_g_mol=120.366),

    # Carbonate system feeding TIC
    "NaHCO3": ChemSpec(strong_ions={"CT_Na": 1.0}, TIC=1.0, molar_mass_g_mol=84.0066),
    "Na2CO3": ChemSpec(strong_ions={"CT_Na": 2.0}, TIC=1.0, molar_mass_g_mol=105.9888),

    # Phosphate feeding CT_P (CT_P is "total P", not charge-specific species)
    "KH2PO4": ChemSpec(strong_ions={"CT_K": 1.0}, P=1.0, molar_mass_g_mol=136.0855),
    "K2HPO4": ChemSpec(strong_ions={"CT_K": 2.0}, P=1.0, molar_mass_g_mol=174.176),

    # Total ammonia
    "NH4Cl":  ChemSpec(strong_ions={"CT_Cl": 1.0}, NH_T=1.0, molar_mass_g_mol=53.491),

    # Weak acid system introduced via its salt (example: sodium acetate)
    "NaAc":   ChemSpec(strong_ions={"CT_Na": 1.0}, acid_system=("AceticAcid", 1.0), molar_mass_g_mol=82.0343),

    # Strong acid / base (represented as unpaired ions)
    # HCl -> contributes Cl- only; electroneutrality forces more H+ => lower pH
    "HCl":    ChemSpec(strong_ions={"CT_Cl": 1.0}, molar_mass_g_mol=36.4609),
    # NaOH -> contributes Na+ only; electroneutrality forces more OH- / less H+ => higher pH
    "NaOH":   ChemSpec(strong_ions={"CT_Na": 1.0}, molar_mass_g_mol=39.997),
    
    # --- Phosphate salts (feed CT_P + K strong ion) ---
    "KH2PO4": ChemSpec(strong_ions={"CT_K": 1.0}, P=1.0, molar_mass_g_mol=136.0855),     # potassium dihydrogen phosphate
    "K2HPO4": ChemSpec(strong_ions={"CT_K": 2.0}, P=1.0, molar_mass_g_mol=174.176),      # dipotassium hydrogen phosphate
    
    # --- Ammonium salt (feeds NH_T and Cl-) ---
    "NH4Cl":  ChemSpec(strong_ions={"CT_Cl": 1.0}, NH_T=1.0, molar_mass_g_mol=53.491),   # ammonium chloride
    
    # --- Sulfate salts ---
    "K2SO4":  ChemSpec(strong_ions={"CT_K": 2.0, "CT_SO4": 1.0}, molar_mass_g_mol=174.259),   # potassium sulphate
    
    # magnesium sulphate heptahydrate (very common in media)
    "MgSO4·7H2O": ChemSpec(strong_ions={"CT_Mg": 1.0, "CT_SO4": 1.0}, molar_mass_g_mol=246.474),
    
    # iron(II) sulphate: choose hydrate you actually add
    "FeSO4": ChemSpec(strong_ions={"CT_SO4": 1.0, "CT_Fe": 1.0}, molar_mass_g_mol=151.908),        # anhydrous
    "FeSO4·7H2O": ChemSpec(strong_ions={"CT_SO4": 1.0, "CT_Fe": 1.0}, molar_mass_g_mol=278.014),   # heptahydrate
    
    # --- Chlorides ---
    "ZnCl2":  ChemSpec(strong_ions={"CT_Zn": 1.0, "CT_Cl": 2.0}, molar_mass_g_mol=136.315),
    
    "MnCl2·7H2O": ChemSpec(strong_ions={"CT_Mn": 1.0, "CT_Cl": 2.0}, molar_mass_g_mol=251.007),
    
    # copper(II) chloride: pick correct form
    "CuCl2": ChemSpec(strong_ions={"CT_Cu": 1.0, "CT_Cl": 2.0}, molar_mass_g_mol=134.452),          # anhydrous
    "CuCl2·2H2O": ChemSpec(strong_ions={"CT_Cu": 1.0, "CT_Cl": 2.0}, molar_mass_g_mol=170.482),     # dihydrate
    
    # cobalt(II) chloride: pick correct form
    "CoCl2": ChemSpec(strong_ions={"CT_Co": 1.0, "CT_Cl": 2.0}, molar_mass_g_mol=129.839),          # anhydrous
    "CoCl2·6H2O": ChemSpec(strong_ions={"CT_Co": 1.0, "CT_Cl": 2.0}, molar_mass_g_mol=237.93),      # hexahydrate
    
    # calcium chloride: pick correct form
    "CaCl2": ChemSpec(strong_ions={"CT_Ca": 1.0, "CT_Cl": 2.0}, molar_mass_g_mol=110.98),           # anhydrous
    "CaCl2·2H2O": ChemSpec(strong_ions={"CT_Ca": 1.0, "CT_Cl": 2.0}, molar_mass_g_mol=147.014),     # dihydrate
    
    "Na2MoO4·2H2O": ChemSpec(strong_ions={"CT_Na": 2.0, "CT_MoO4": 1.0}, molar_mass_g_mol=241.95),
    "Na2MoO4": ChemSpec(strong_ions={"CT_Na": 2.0, "CT_MoO4": 1.0}, molar_mass_g_mol=205.92),

    "Biotin": ChemSpec(strong_ions={}, molar_mass_g_mol=244.31),
    
    "CitricAcid": ChemSpec(strong_ions={},                 # pure acid, no counter-ion
                           acid_system=("CitricAcid", 1.0),
                           molar_mass_g_mol=192.124
                           ),

    # Manganese(II) chloride tetrahydrate
    # MnCl2·4H2O -> Mn2+ + 2 Cl-
    "MnCl2·4H2O": ChemSpec(
        strong_ions={"CT_Mn": 1.0, "CT_Cl": 2.0},
        molar_mass_g_mol=197.91,
    ),
    
    # Copper(II) chloride dihydrate
    # CuCl2·2H2O -> Cu2+ + 2 Cl-
    "CuCl2·2H2O": ChemSpec(
        strong_ions={"CT_Cu": 1.0, "CT_Cl": 2.0},
        molar_mass_g_mol=170.48,
    ),
    
    # Cobalt(II) chloride hexahydrate
    # CoCl2·6H2O -> Co2+ + 2 Cl-
    "CoCl2·6H2O": ChemSpec(
        strong_ions={"CT_Co": 1.0, "CT_Cl": 2.0},
        molar_mass_g_mol=237.93,
    ),
    
    # Sodium molybdate dihydrate
    # Na2MoO4·2H2O -> 2 Na+ + MoO4(2-)
    #
    # IMPORTANT:
    # - If your solver expects CT_Mo7O24 (heptamolybdate), we approximate:
    #     1 mol MoO4(2-) corresponds to 1 mol Mo, while 1 mol Mo7O24 corresponds to 7 mol Mo
    #   so we map: CT_Mo7O24 += (1/7)*moles
    "Na2MoO4·2H2O": ChemSpec(
        strong_ions={"CT_Na": 2.0, "CT_Mo7O24": 1.0/7.0},
        molar_mass_g_mol=241.95,
    ),
    
    # Calcium chloride dihydrate
    # CaCl2·2H2O -> Ca2+ + 2 Cl-
    "CaCl2·2H2O": ChemSpec(
        strong_ions={"CT_Ca": 1.0, "CT_Cl": 2.0},
        molar_mass_g_mol=147.02,
    ),

    # Sulfuric acid (adds sulfate with no counter-cation → drives pH down via electroneutrality)
    # H2SO4 → (effectively) total sulfate pool; the solver will split into HSO4- / SO4-- by KaS.
    "H2SO4": ChemSpec(
        strong_ions={"CT_SO4": 1.0},
        molar_mass_g_mol=98.079,   # g/mol
    ),
    
    # --- Additions for standards/benchmark recipes ---
    # Acetic acid (weak acid; no counter-ion)
    "AceticAcid": ChemSpec(
        strong_ions={},
        acid_system=("AceticAcid", 1.0),
        molar_mass_g_mol=60.052,   # CH3COOH
    ),
    
    # Sodium phosphate dibasic (for McIlvaine table stock: 0.2 M Na2HPO4)
    # Note: hydrate form varies in practice; for mol/L stock solutions, treat as "species source".
    "Na2HPO4": ChemSpec(
        strong_ions={"CT_Na": 2.0},
        P=1.0,
        molar_mass_g_mol=141.96,   # anhydrous Na2HPO4
    ),
    
    # Potassium hydroxide (strong base; represented as K+ only, same convention as NaOH in your DB)
    "KOH": ChemSpec(
        strong_ions={"CT_K": 1.0},
        molar_mass_g_mol=56.1056,
    ),
    
    # --- Aliases for common hydrate naming (ASCII '.' instead of Unicode '·') ---
    # Some users type MgSO4.7H2O instead of MgSO4·7H2O on Windows.
    "MgSO4.7H2O": ChemSpec(strong_ions={"CT_Mg": 1.0, "CT_SO4": 1.0}, molar_mass_g_mol=246.474),
    
    # --- Missing phosphate salts (needed for phosphate-buffer standards) ---
    # Sodium dihydrogen phosphate (monobasic)
    "NaH2PO4": ChemSpec(strong_ions={"CT_Na": 1.0}, P=1.0, molar_mass_g_mol=119.98),
    
    # Optional hydrate form used in many recipes
    "NaH2PO4·H2O": ChemSpec(strong_ions={"CT_Na": 1.0}, P=1.0, molar_mass_g_mol=137.99),
    
    # Sodium hydrogen phosphate dihydrate (dibasic dihydrate) used in many media/buffer recipes
    "Na2HPO4·2H2O": ChemSpec(strong_ions={"CT_Na": 2.0}, P=1.0, molar_mass_g_mol=177.99),
    
    # Optional ASCII alias too
    "Na2HPO4.2H2O": ChemSpec(strong_ions={"CT_Na": 2.0}, P=1.0, molar_mass_g_mol=177.99),

    # --- Total ammonia as neutral aqueous ammonia (NH3(aq)) ---
    # This lets recipes represent "aqueous NH3" (total ammonia) without adding counter-ions.
    "NH3": ChemSpec(
        strong_ions={},     # neutral solute, no strong-ion contribution
        NH_T=1.0,           # 1 mol NH_T per mol NH3(aq)
        molar_mass_g_mol=17.031
    ),
    
    # Optional alias (common lab naming)
    "NH4OH": ChemSpec(
        strong_ions={},
        NH_T=1.0,
        molar_mass_g_mol=35.045  # if you ever enter g/L as NH4OH equivalent
    ),


}


def recipe_to_totals(
    recipe_mol_L: Dict[str, float],
    *,
    acid_pKas: Dict[str, float],
) -> Dict:
    """
    Convert a chemical recipe (mol/L per chemical) into inputs for BisectionChemicalEquilibriumEngine.solve.

    Returns dict with:
      strong_ions (CT_*),
      acid_totals, acid_pKas,
      CT_TIC, CT_P, CT_NH_T
    """
    strong_ions: Dict[str, float] = {}
    acid_totals: Dict[str, float] = {}
    CT_TIC = 0.0
    CT_P = 0.0
    CT_NH_T = 0.0

    for chem, c in recipe_mol_L.items():
        if c == 0:
            continue
        if chem not in CHEM_DB:
            raise KeyError(f"Chemical '{chem}' not in CHEM_DB. Add it to scripts/chem_recipe.py")

        spec = CHEM_DB[chem]

        # strong ions
        for k, coeff in spec.strong_ions.items():
            strong_ions[k] = strong_ions.get(k, 0.0) + coeff * c

        # weak acid systems
        if spec.acid_system is not None:
            acid_name, stoich = spec.acid_system
            acid_totals[acid_name] = acid_totals.get(acid_name, 0.0) + stoich * c
            if acid_name not in acid_pKas:
                raise KeyError(
                    f"Missing pKa for acid system '{acid_name}'. "
                    f"Provide it in acid_pKas."
                )

        # global totals
        CT_TIC += spec.TIC * c
        CT_P += spec.P * c
        CT_NH_T += spec.NH_T * c

    return dict(
        strong_ions=strong_ions,
        acid_totals=acid_totals,
        acid_pKas=acid_pKas,
        CT_TIC=CT_TIC,
        CT_P=CT_P,
        CT_NH_T=CT_NH_T,
    )


def recipe_g_L_to_mol_L(recipe_g_L: Dict[str, float]) -> Dict[str, float]:
    """Optional helper if you prefer entering g/L."""
    out: Dict[str, float] = {}
    for chem, gL in recipe_g_L.items():
        if gL == 0:
            continue
        if chem not in CHEM_DB or CHEM_DB[chem].molar_mass_g_mol is None:
            raise KeyError(f"Need molar_mass_g_mol for '{chem}' in CHEM_DB to use g/L.")
        out[chem] = gL / CHEM_DB[chem].molar_mass_g_mol
    return out
