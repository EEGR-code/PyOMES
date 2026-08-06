# -*- coding: utf-8 -*-
"""
Core acid-base speciation solver (Level 1 foundation).

Systems included:
- Water: H+/OH-
- Monoprotic acids (VFAs): HA/A-
- TIC carbonate ladder: CO2(aq)/HCO3-/CO3--
- TAN: NH4+/NH3(aq)
- Phosphate ladder: H3PO4/H2PO4-/HPO4--/PO4---
- Bisulfate: HSO4-/SO4--

Strong ions are treated as fully dissociated fixed totals (free ions) in this solver.

Activity handling:
- If an ActivityModel other than ideal is provided, we solve with activity corrections
  using an outer fixed-point iteration on ionic strength I (Davies).
- pH reported is:
    - ideal: pH = -log10([H+])
    - activity: pH = -log10(aH) = -log10(gamma_H * [H+])
  and we also return pH_conc in both cases for convenience.

Returned dict contains at least:
- pH, pH_conc, H+, OH-, logH
- CO2aq, HCO3-, CO3--
- NH4+, NH3
- phosphate species
- sulfate species
- strong ions echoed
- IonicStrength
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple, Union, Optional, Sequence
AcidPKa = Union[float, Sequence[float]]

import numpy as np
from scipy.optimize import brentq

from .activity import ionic_strength_from_speciation
from .activity_models import ActivityModel


# Gas constant (J/mol/K) for van 't Hoff temperature corrections
_R_J_PER_MOLK = 8.31446261815324


def _vant_hoff_K(
    K_ref: float,
    dH_J_per_mol: float,
    T_K: float,
    T_ref_K: float = 298.15,
) -> float:
    """Temperature-correct an equilibrium constant using the van 't Hoff relation.

    ln K(T) = ln K(T_ref) - (ΔH°/R) * (1/T - 1/T_ref)

    Notes:
      - ΔH° should be the *standard enthalpy change of the equilibrium reaction*.
      - If dH_J_per_mol is 0, this reduces to K(T) = K_ref (no temperature effect).
      - This function intentionally keeps the wider codebase interface unchanged.
    """
    K_ref = float(K_ref)
    dH_J_per_mol = float(dH_J_per_mol)
    T_K = float(T_K)
    T_ref_K = float(T_ref_K)
    if not np.isfinite(K_ref) or K_ref <= 0.0:
        return float(K_ref)
    if not np.isfinite(dH_J_per_mol) or abs(dH_J_per_mol) < 1e-30:
        return float(K_ref)
    if not np.isfinite(T_K) or T_K <= 0.0:
        return float(K_ref)
    return float(K_ref * np.exp(-(dH_J_per_mol / _R_J_PER_MOLK) * (1.0 / T_K - 1.0 / T_ref_K)))


# Charges for common ions (used for gamma application and a few computed keys)
_Z = {
    "H+": +1,
    "OH-": -1,
    "HCO3-": -1,
    "CO3--": -2,
    "NH4+": +1,
    "H2PO4-": -1,
    "HPO4--": -2,
    "PO4---": -3,
    "HSO4-": -1,
    "SO4--": -2,
    "Cu++": +2,
    "Fe++": +2,
    "Mo7O24------": -6,

}


def _safe_float(x: Any, default: float = 0.0) -> float:
    try:
        return float(x)
    except (TypeError, ValueError, AttributeError, KeyError):
        return float(default)


def _pick_bracket_from_guess(
    *,
    logH_guess: Optional[float],
    pH_min: float,
    pH_max: float,
) -> Tuple[float, float]:
    if logH_guess is None or not np.isfinite(logH_guess):
        return pH_min, pH_max
    pH_guess = float(-logH_guess)
    a = max(pH_min, pH_guess - 3.0)
    b = min(pH_max, pH_guess + 3.0)
    if b <= a:
        return pH_min, pH_max
    return a, b

def _as_pka_list(pka: AcidPKa) -> list[float]:
    if isinstance(pka, (list, tuple)):
        return [float(x) for x in pka]
    return [float(pka)]

def _polyprotic_alpha(H: float, Kas: list[float]) -> list[float]:
    """
    Return alpha fractions for H_nA, H_{n-1}A-, ..., A^{n-}
    where Kas are [Ka1, Ka2, ... Kan] for an n-protic acid.
    """
    n = len(Kas)
    # Denominator: sum_{i=0..n} (prod_{j=1..i} Ka_j) * H^{n-i}
    denom = 0.0
    prods = [1.0]
    for k in Kas:
        prods.append(prods[-1] * k)

    for i in range(n + 1):
        denom += prods[i] * (H ** (n - i))

    if denom <= 0.0:
        # fallback: all fully protonated
        return [1.0] + [0.0] * n

    alphas = []
    for i in range(n + 1):
        alphas.append((prods[i] * (H ** (n - i))) / denom)
    return alphas

def _polyprotic_charge(CT: float, H: float, Kas: list[float]) -> float:
    """
    Return total anionic charge contribution (mol/L equivalent of negative charge),
    i.e. 1*[H_{n-1}A-] + 2*[H_{n-2}A2-] + ... + n*[A^{n-}].
    """
    alphas = _polyprotic_alpha(H, Kas)
    # alpha[0] is neutral H_nA; alpha[i] has charge -i
    return sum(i * CT * alphas[i] for i in range(1, len(alphas)))


def solve_acid_base(
    *,
    acid_totals: Dict[str, float],
    acid_pKas: Dict[str, AcidPKa],    
    CT_TIC: float,
    CT_NH_T: float,
    CT_P: float,

    # Strong ions
    CT_K: float = 0.0, CT_Na: float = 0.0, CT_Cl: float = 0.0, CT_NO3: float = 0.0, CT_SO4: float = 0.0,
    CT_Mg: float = 0.0, CT_Ca: float = 0.0, CT_Zn: float = 0.0, CT_Mn: float = 0.0, CT_Co: float = 0.0,
    CT_Mo7O24: float = 0.0,
    CT_cation: float = 0.0, CT_anion: float = 0.0,

    # Equilibrium constants
    Kw: float = 1e-14,
    pKa1_TIC: float = 6.35,
    pKa2_TIC: float = 10.33,
    pKa_NH: float = 9.25,
    pKa1_P: float = 2.15,
    pKa2_P: float = 7.20,
    pKa3_P: float = 12.35,
    pKa_HSO4: float = 1.99,

    # --- NEW: van 't Hoff temperature corrections (defaults keep legacy behavior) ---
    # Provide ΔH° (J/mol) for each equilibrium to enable temperature-dependent K.
    # If left at 0.0, the constant is treated as temperature-invariant (legacy).
    dH_Kw_J_per_mol: float = 0.0,
    dH_TIC1_J_per_mol: float = 0.0,
    dH_TIC2_J_per_mol: float = 0.0,
    dH_NH_J_per_mol: float = 0.0,
    dH_P1_J_per_mol: float = 0.0,
    dH_P2_J_per_mol: float = 0.0,
    dH_P3_J_per_mol: float = 0.0,
    dH_HSO4_J_per_mol: float = 0.0,
    acid_dH_J_per_mol: Optional[Dict[str, float]] = None,
    T_ref_K: float = 298.15,

    # Root solver
    pH_min: float = -0.5,
    pH_max: float = 20.0,
    n_scan: int = 200,
    tol: float = 1e-12,

    # Activity model + temperature
    activity_model: ActivityModel,
    T_K: float,

    # Activity iteration controls
    max_outer: int = 20,
    I_init: float = 0.01,
    I_tol: float = 1e-8,
    damping: float = 0.6,

    # Warm start
    logH_guess: Optional[float] = None,
) -> Dict[str, Any]:
    
    # Ensure every acid present in totals has pKas provided
    missing = [name for name in acid_totals.keys() if name not in acid_pKas]
    if missing:
        raise KeyError(
            "Missing acid_pKas entries for: " + ", ".join(missing)
            + ". Provide pKas for all acids present in acid_totals."
        )

    

    # Normalize inputs
    CT_TIC = _safe_float(CT_TIC)
    CT_NH_T = _safe_float(CT_NH_T)
    CT_P = _safe_float(CT_P)

    # --- van 't Hoff temperature corrections for equilibrium constants ---
    # We keep the wider interface unchanged: callers can optionally provide ΔH° values
    # (J/mol) for each equilibrium. If ΔH° is left at 0.0, the constant remains
    # temperature-invariant (legacy behavior).
    acid_dH = dict(acid_dH_J_per_mol or {})

    # Kw is used directly in the water equilibrium (OH = Kw/H)
    Kw = _vant_hoff_K(Kw, dH_Kw_J_per_mol, float(T_K), float(T_ref_K))

    # pH window preferences from warm-start
    if logH_guess is not None and np.isfinite(logH_guess):
        pH_guess = -float(logH_guess)
        pH_lo_pref = max(float(pH_min), pH_guess - 2.0)
        pH_hi_pref = min(float(pH_max), pH_guess + 2.0)
    else:
        pH_lo_pref, pH_hi_pref = float(pH_min), float(pH_max)
        
    # --- NEW: model-agnostic gamma wrapper (supports gamma(z, I) and gamma(z, I, T_K)) ---
    # Supports:
    #   gamma(z, I, *, T_K=...)   (your DaviesActivityModel)
    #   gamma(z, I, T_K)          (positional temperature)
    #   gamma(z, I)              (no temperature)
    def _gamma(z: int, I_val: float) -> float:
        I_val = max(0.0, float(I_val))
        z = int(z)
        # 1) Preferred: keyword-only temperature (DaviesActivityModel in your codebase)
        try:
            return float(activity_model.gamma(z, I_val, T_K=float(T_K)))
        except TypeError:
            pass
        # 2) Some models accept temperature positionally
        try:
            return float(activity_model.gamma(z, I_val, float(T_K)))
        except TypeError:
            pass
        # 3) Some models don't use temperature at all
        return float(activity_model.gamma(z, I_val))

    
    def gammas_from_I(I_val: float) -> Dict[str, float]:
        I_val = max(0.0, float(I_val))

        # Stage 16: If the activity model has compute_gammas (e.g. SIT),
        # call it with the full solution composition for ion-pair-specific
        # corrections.  Otherwise fall back to charge-only gamma calls.
        if hasattr(activity_model, 'compute_gammas'):
            # Build composition dict from solver scope variables
            comp = {}
            if CT_Na > 0: comp["Na+"] = CT_Na
            if CT_K > 0: comp["K+"] = CT_K
            if CT_Cl > 0: comp["Cl-"] = CT_Cl
            if CT_NO3 > 0: comp["NO3-"] = CT_NO3
            if CT_Ca > 0: comp["Ca++"] = CT_Ca
            if CT_Mg > 0: comp["Mg++"] = CT_Mg
            if CT_SO4 > 0: comp["SO4--"] = CT_SO4
            if CT_NH_T > 0: comp["NH4+"] = CT_NH_T * 0.95  # approx: most is NH₄⁺
            if CT_TIC > 0: comp["HCO3-"] = CT_TIC * 0.8     # approx: most is HCO₃⁻ at AD pH
            # VFA anions from acid_totals (mostly dissociated at AD pH)
            for acid_name, acid_conc in acid_totals.items():
                if acid_conc > 0:
                    anion = acid_name.replace("Acid", "") + "-"
                    if acid_name == "AceticAcid":
                        anion = "Acetate-"
                    elif acid_name == "H2S":
                        anion = "HS-"
                    else:
                        anion = acid_name + "-"
                    comp[anion] = acid_conc * 0.95  # approximate dissociated fraction
            try:
                return activity_model.compute_gammas(comp, I_val, float(T_K))
            except (TypeError, ValueError, AttributeError, KeyError):
                pass  # fall back to charge-only

        return {
            "H+": _gamma(+1, I_val),
            "OH-": _gamma(-1, I_val),
            "HCO3-": _gamma(-1, I_val),
            "CO3--": _gamma(-2, I_val),
            "NH4+": _gamma(+1, I_val),
            "H2PO4-": _gamma(-1, I_val),
            "HPO4--": _gamma(-2, I_val),
            "PO4---": _gamma(-3, I_val),
            "HSO4-": _gamma(-1, I_val),
            "SO4--": _gamma(-2, I_val),
        }


    def effective_constants(gammas: Dict[str, float]) -> Dict[str, Any]:
        """
        Build activity- and temperature-corrected equilibrium constants.

        CHANGE (Update 3):
        - Supports monoprotic acids (pKa is float) AND polyprotic acids (pKa is list/tuple).
        - Ka_vfa_eff[name] becomes a LIST of Ka steps [Ka1, Ka2, ...].
          (Monoprotic acids will have a list of length 1.)
        - Applies activity correction per step using gamma for the product species charge (-i).
        """

        def _as_pka_list(pka_val):
            if isinstance(pka_val, (list, tuple)):
                return [float(x) for x in pka_val]
            return [float(pka_val)]

        Ka_vfa_eff: Dict[str, list] = {}

        for name in acid_totals.keys():
            # guaranteed present due to the validation block above
            pKa_val = acid_pKas[name]
            pka_list = _as_pka_list(pKa_val)
        
            # Temperature-correct each step
            Ka_steps_ref_T = []
            for p in pka_list:
                Ka_ref = 10.0 ** (-float(p))
                Ka_T = _vant_hoff_K(Ka_ref, acid_dH.get(name, 0.0), float(T_K), float(T_ref_K))
                Ka_steps_ref_T.append(float(Ka_T))
        
            # Activity-correct each step: step i produces charge -i
            Ka_steps_eff = []
            for i, Ka_T in enumerate(Ka_steps_ref_T, start=1):
                gamma_Ai = _gamma(-i, I)  # uses current outer-loop ionic strength
                Ka_steps_eff.append(Ka_T / (gammas["H+"] * gamma_Ai))
        
            Ka_vfa_eff[name] = Ka_steps_eff


        Ka1 = _vant_hoff_K(10.0 ** (-float(pKa1_TIC)), dH_TIC1_J_per_mol, float(T_K), float(T_ref_K))
        Ka2 = _vant_hoff_K(10.0 ** (-float(pKa2_TIC)), dH_TIC2_J_per_mol, float(T_K), float(T_ref_K))
        KaN = _vant_hoff_K(10.0 ** (-float(pKa_NH)),   dH_NH_J_per_mol,   float(T_K), float(T_ref_K))
        KaP1 = _vant_hoff_K(10.0 ** (-float(pKa1_P)),  dH_P1_J_per_mol,   float(T_K), float(T_ref_K))
        KaP2 = _vant_hoff_K(10.0 ** (-float(pKa2_P)),  dH_P2_J_per_mol,   float(T_K), float(T_ref_K))
        KaP3 = _vant_hoff_K(10.0 ** (-float(pKa3_P)),  dH_P3_J_per_mol,   float(T_K), float(T_ref_K))
        KaS  = _vant_hoff_K(10.0 ** (-float(pKa_HSO4)), dH_HSO4_J_per_mol, float(T_K), float(T_ref_K))

        return dict(
            Ka_vfa_eff=Ka_vfa_eff,
            Ka1_eff=Ka1 / (gammas["H+"] * gammas["HCO3-"]),
            Ka2_eff=Ka2 * (gammas["HCO3-"] / (gammas["H+"] * gammas["CO3--"])),
            KaN_eff=KaN * (gammas["NH4+"] / gammas["H+"]),
            KaP1_eff=KaP1 / (gammas["H+"] * gammas["H2PO4-"]),
            KaP2_eff=KaP2 * (gammas["H2PO4-"] / (gammas["H+"] * gammas["HPO4--"])),
            KaP3_eff=KaP3 * (gammas["HPO4--"] / (gammas["H+"] * gammas["PO4---"])),
            KaS_eff=KaS * (gammas["HSO4-"] / (gammas["H+"] * gammas["SO4--"])),
        )


    def make_residual(gammas: Dict[str, float], K_eff: Dict[str, Any]):
        def residual(pH: float) -> float:
            return electroneutrality_residual(
                pH=pH,
                acid_totals=acid_totals, Ka_vfa_eff=K_eff["Ka_vfa_eff"],
                CT_TIC=CT_TIC, Ka1_eff=K_eff["Ka1_eff"], Ka2_eff=K_eff["Ka2_eff"],
                CT_NH_T=CT_NH_T, KaN_eff=K_eff["KaN_eff"],
                CT_P=CT_P, KaP1_eff=K_eff["KaP1_eff"], KaP2_eff=K_eff["KaP2_eff"], KaP3_eff=K_eff["KaP3_eff"],
                CT_SO4=_safe_float(CT_SO4), KaS_eff=K_eff["KaS_eff"],
                Kw=Kw, gamma_H=gammas["H+"], gamma_OH=gammas["OH-"],
                strong_ions=dict(
                    CT_K=_safe_float(CT_K), CT_Na=_safe_float(CT_Na), CT_Cl=_safe_float(CT_Cl), CT_NO3=_safe_float(CT_NO3),
                    CT_Mg=_safe_float(CT_Mg), CT_Ca=_safe_float(CT_Ca), CT_Zn=_safe_float(CT_Zn), CT_Mn=_safe_float(CT_Mn),
                    CT_Co=_safe_float(CT_Co), CT_Mo7O24=_safe_float(CT_Mo7O24), CT_cation=_safe_float(CT_cation), CT_anion=_safe_float(CT_anion),
                ),
            )
        return residual
    
    def electroneutrality_residual(
        *,
        pH: float,
        acid_totals: Dict[str, float],
        Ka_vfa_eff: Dict[str, Any],
    
        CT_TIC: float, Ka1_eff: float, Ka2_eff: float,
        CT_NH_T: float, KaN_eff: float,
        CT_P: float, KaP1_eff: float, KaP2_eff: float, KaP3_eff: float,
        CT_SO4: float, KaS_eff: float,
    
        Kw: float,
        gamma_H: float,
        gamma_OH: float,
    
        strong_ions: Dict[str, float],
    ) -> float:
        """
        Return electroneutrality residual:
            sum(cations) - sum(anions)
    
        CHANGE (Update 4):
        - Handles polyprotic acids if Ka_vfa_eff[name] is a list of Ka steps.
        - Monoprotic acids behave exactly as before.
        """
    
        def _polyprotic_charge(CT: float, H: float, Kas: list[float]) -> float:
            """
            For an n-protic acid system with totals CT and steps Kas=[Ka1..Kan],
            return the total anionic charge contribution:
                1*[H_{n-1}A-] + 2*[H_{n-2}A2-] + ... + n*[A^{n-}]
            """
            n = len(Kas)
            if n == 1:
                Ka = float(Kas[0])
                denom = Ka + H
                return (CT * Ka / denom) if denom > 0 else 0.0
    
            # alpha fractions for H_nA, H_{n-1}A-, ..., A^{n-}
            prods = [1.0]
            for k in Kas:
                prods.append(prods[-1] * float(k))
    
            denom = 0.0
            for i in range(n + 1):
                denom += prods[i] * (H ** (n - i))
    
            if denom <= 0.0:
                return 0.0
    
            charge = 0.0
            for i in range(1, n + 1):
                alpha_i = (prods[i] * (H ** (n - i))) / denom
                charge += i * CT * alpha_i
            return charge
    
        H = 10.0 ** (-float(pH))
    
        # Water (activity corrected if gammas != 1)
        OH = Kw / (gamma_H * gamma_OH * H) if H > 0 else 0.0
    
        # VFAs + other acids: monoprotic or polyprotic
        A_charge_sum = 0.0
        for name, CT in acid_totals.items():
            CT = _safe_float(CT)
            Kas = Ka_vfa_eff.get(name, [])
            if Kas is None:
                continue
    
            # Allow either a float (legacy) or list (new)
            if isinstance(Kas, (float, int)):
                Kas_list = [float(Kas)]
            else:
                Kas_list = [float(x) for x in Kas]  # list/tuple/etc.
    
            if not Kas_list:
                continue
    
            A_charge_sum += _polyprotic_charge(CT, H, Kas_list)
    
        # Carbonate ladder (CO2 / HCO3- / CO3--)
        denomC = H * H + Ka1_eff * H + Ka1_eff * Ka2_eff
        if denomC > 0.0:
            HCO3 = CT_TIC * (Ka1_eff * H) / denomC
            CO3 = CT_TIC * (Ka1_eff * Ka2_eff) / denomC
        else:
            HCO3 = 0.0
            CO3 = 0.0
    
        # TAN (NH4+ / NH3)
        denomN = KaN_eff + H
        NH4 = (CT_NH_T * H / denomN) if denomN > 0 else 0.0
    
        # Phosphate ladder
        denomP = H**3 + KaP1_eff * H**2 + KaP1_eff * KaP2_eff * H + KaP1_eff * KaP2_eff * KaP3_eff
        if denomP > 0.0:
            H2PO4 = CT_P * (KaP1_eff * H**2) / denomP
            HPO4 = CT_P * (KaP1_eff * KaP2_eff * H) / denomP
            PO4 = CT_P * (KaP1_eff * KaP2_eff * KaP3_eff) / denomP
        else:
            H2PO4 = 0.0
            HPO4 = 0.0
            PO4 = 0.0
    
        # Bisulfate (HSO4- / SO4--)
        denomS = KaS_eff + H
        HSO4 = (CT_SO4 * H / denomS) if denomS > 0 else 0.0
        SO4 = CT_SO4 - HSO4
    
        # Strong ions (fully dissociated)
        CT_K = _safe_float(strong_ions.get("CT_K", 0.0))
        CT_Na = _safe_float(strong_ions.get("CT_Na", 0.0))
        CT_Cl = _safe_float(strong_ions.get("CT_Cl", 0.0))
        CT_NO3 = _safe_float(strong_ions.get("CT_NO3", 0.0))
    
        CT_Mg = _safe_float(strong_ions.get("CT_Mg", 0.0))
        CT_Ca = _safe_float(strong_ions.get("CT_Ca", 0.0))
        CT_Zn = _safe_float(strong_ions.get("CT_Zn", 0.0))
        CT_Mn = _safe_float(strong_ions.get("CT_Mn", 0.0))
        CT_Co = _safe_float(strong_ions.get("CT_Co", 0.0))
        CT_Mo7O24 = _safe_float(strong_ions.get("CT_Mo7O24", 0.0))
        CT_cation = _safe_float(strong_ions.get("CT_cation", 0.0))
        CT_anion = _safe_float(strong_ions.get("CT_anion", 0.0))
    
        # Charge balance: cations - anions
        cations = (
            H
            + NH4
            + CT_K
            + CT_Na
            + CT_cation
            + 2.0 * CT_Mg
            + 2.0 * CT_Ca
            + 2.0 * CT_Zn
            + 2.0 * CT_Mn
            + 2.0 * CT_Co
        )
    
        anions = (
            OH
            + CT_Cl
            + CT_NO3
            + CT_anion
            + A_charge_sum
            + HCO3
            + 2.0 * CO3
            + H2PO4
            + 2.0 * HPO4
            + 3.0 * PO4
            + HSO4
            + 2.0 * SO4
            + 6.0 * CT_Mo7O24
        )
    
        return float(cations - anions)



    def solve_pH(residual, pH_lo: float, pH_hi: float) -> float:
        a, b = float(pH_lo), float(pH_hi)
        fa, fb = float(residual(a)), float(residual(b))
        if np.isfinite(fa) and np.isfinite(fb) and np.sign(fa) != np.sign(fb):
            return float(brentq(lambda x: residual(x), a, b, xtol=tol))
        xs = np.linspace(float(pH_min), float(pH_max), int(n_scan))
        fs = np.array([residual(float(x)) for x in xs], dtype=float)
        sgn = np.sign(fs)
        idx = np.where(sgn[:-1] * sgn[1:] < 0)[0]
        if len(idx):
            lo = float(xs[idx[0]])
            hi = float(xs[idx[0] + 1])
            return float(brentq(lambda x: residual(x), lo, hi, xtol=tol))
        i0 = int(np.argmin(np.abs(fs)))
        return float(xs[i0])

    # Ideal fast path
    if getattr(activity_model, "name", "ideal") == "ideal":
        I = 0.0
        gammas = gammas_from_I(I)
        K_eff = effective_constants(gammas)
        residual = make_residual(gammas, K_eff)
        pH_sol = solve_pH(residual, pH_lo=pH_lo_pref, pH_hi=pH_hi_pref)

        sp = compute_species(
            pH=float(pH_sol),
            acid_totals=acid_totals,
            Ka_vfa_eff=K_eff["Ka_vfa_eff"],
            CT_TIC=CT_TIC, Ka1_eff=K_eff["Ka1_eff"], Ka2_eff=K_eff["Ka2_eff"],
            CT_NH_T=CT_NH_T, KaN_eff=K_eff["KaN_eff"],
            CT_P=CT_P, KaP1_eff=K_eff["KaP1_eff"], KaP2_eff=K_eff["KaP2_eff"], KaP3_eff=K_eff["KaP3_eff"],
            CT_SO4=_safe_float(CT_SO4), KaS_eff=K_eff["KaS_eff"],
            Kw=Kw, gamma_H=gammas["H+"], gamma_OH=gammas["OH-"],
            strong_ions=dict(
                CT_K=_safe_float(CT_K), CT_Na=_safe_float(CT_Na), CT_Cl=_safe_float(CT_Cl), CT_NO3=_safe_float(CT_NO3),
                CT_Mg=_safe_float(CT_Mg), CT_Ca=_safe_float(CT_Ca), CT_Zn=_safe_float(CT_Zn), CT_Mn=_safe_float(CT_Mn),
                CT_Co=_safe_float(CT_Co), CT_Mo7O24=_safe_float(CT_Mo7O24), CT_cation=_safe_float(CT_cation), CT_anion=_safe_float(CT_anion),
            ),
        )

        H = float(sp["H+"])
        gamma_H = float(gammas["H+"])
        pH_conc = -np.log10(H)
        sp["pH_conc"] = float(pH_conc)
        sp["aH"] = float(gamma_H * H)
        sp["gamma_H"] = float(gamma_H)
        sp["gamma_OH"] = float(gammas["OH-"])
        sp["IonicStrength"] = float(ionic_strength_from_speciation(sp))
        sp["pH"] = float(pH_conc)
        sp["logH"] = float(np.log10(H))
        return sp

    # Davies/activity outer iteration
    I = float(I_init)
    converged = False
    last_delta_I = None
    sp_last = None
    gammas_last = None
    pH_last = None

    for _outer in range(int(max_outer)):
        gammas = gammas_from_I(I)
        K_eff = effective_constants(gammas)
        residual = make_residual(gammas, K_eff)

        pH_sol = solve_pH(residual, pH_lo=pH_lo_pref, pH_hi=pH_hi_pref)

        sp = compute_species(
            pH=float(pH_sol),
            acid_totals=acid_totals,
            Ka_vfa_eff=K_eff["Ka_vfa_eff"],
            CT_TIC=CT_TIC, Ka1_eff=K_eff["Ka1_eff"], Ka2_eff=K_eff["Ka2_eff"],
            CT_NH_T=CT_NH_T, KaN_eff=K_eff["KaN_eff"],
            CT_P=CT_P, KaP1_eff=K_eff["KaP1_eff"], KaP2_eff=K_eff["KaP2_eff"], KaP3_eff=K_eff["KaP3_eff"],
            CT_SO4=_safe_float(CT_SO4), KaS_eff=K_eff["KaS_eff"],
            Kw=Kw, gamma_H=gammas["H+"], gamma_OH=gammas["OH-"],
            strong_ions=dict(
                CT_K=_safe_float(CT_K), CT_Na=_safe_float(CT_Na), CT_Cl=_safe_float(CT_Cl), CT_NO3=_safe_float(CT_NO3),
                CT_Mg=_safe_float(CT_Mg), CT_Ca=_safe_float(CT_Ca), CT_Zn=_safe_float(CT_Zn), CT_Mn=_safe_float(CT_Mn),
                CT_Co=_safe_float(CT_Co), CT_Mo7O24=_safe_float(CT_Mo7O24), CT_cation=_safe_float(CT_cation), CT_anion=_safe_float(CT_anion),
            ),
        )

        # Save last computed snapshot (used on convergence / fallback)
        sp_last = sp
        gammas_last = gammas
        pH_last = float(pH_sol)

        I_new = float(ionic_strength_from_speciation(sp))
        last_delta_I = abs(I_new - I)

        if last_delta_I < float(I_tol):
            I = I_new
            converged = True
            break

        I = float(damping) * I + (1.0 - float(damping)) * I_new

    if getattr(activity_model, "name", "ideal") != "ideal" and not converged:
        import warnings
        warnings.warn(
            f"Davies activity iteration did not converge after {max_outer} iterations "
            f"(last |ΔI|={last_delta_I:.2e} mol/L, I≈{I:.3g} mol/L).",
            RuntimeWarning,
        )

    # Use the last computed snapshot from the outer loop (already at the best-available I).
    # This avoids an extra full recomputation at the end of the activity iteration.
    if sp_last is None or gammas_last is None or pH_last is None:
        gammas_last = gammas_from_I(I)
        K_eff = effective_constants(gammas_last)
        pH_last = float(pH_sol)
        sp_last = compute_species(
            pH=float(pH_last),
            acid_totals=acid_totals,
            Ka_vfa_eff=K_eff["Ka_vfa_eff"],
            CT_TIC=CT_TIC, Ka1_eff=K_eff["Ka1_eff"], Ka2_eff=K_eff["Ka2_eff"],
            CT_NH_T=CT_NH_T, KaN_eff=K_eff["KaN_eff"],
            CT_P=CT_P, KaP1_eff=K_eff["KaP1_eff"], KaP2_eff=K_eff["KaP2_eff"], KaP3_eff=K_eff["KaP3_eff"],
            CT_SO4=_safe_float(CT_SO4), KaS_eff=K_eff["KaS_eff"],
            Kw=Kw, gamma_H=gammas_last["H+"], gamma_OH=gammas_last["OH-"],
            strong_ions=dict(
                CT_K=_safe_float(CT_K), CT_Na=_safe_float(CT_Na), CT_Cl=_safe_float(CT_Cl), CT_NO3=_safe_float(CT_NO3),
                CT_Mg=_safe_float(CT_Mg), CT_Ca=_safe_float(CT_Ca), CT_Zn=_safe_float(CT_Zn), CT_Mn=_safe_float(CT_Mn),
                CT_Co=_safe_float(CT_Co), CT_Mo7O24=_safe_float(CT_Mo7O24), CT_cation=_safe_float(CT_cation), CT_anion=_safe_float(CT_anion),
            ),
        )

    sp = sp_last
    gammas = gammas_last

    H = float(sp["H+"])
    gamma_H = float(gammas["H+"])
    pH_conc = -np.log10(H)
    pH_act = -np.log10(gamma_H * H)

    sp["pH_conc"] = float(pH_conc)
    sp["aH"] = float(gamma_H * H)
    sp["gamma_H"] = float(gamma_H)
    sp["gamma_OH"] = float(gammas["OH-"])
    sp["IonicStrength"] = float(ionic_strength_from_speciation(sp))

    if getattr(activity_model, "name", "ideal") == "ideal":
        sp["pH"] = float(pH_conc)
    else:
        sp["pH"] = float(pH_act)

    sp["logH"] = float(np.log10(H))
    return sp


def compute_species(
    *,
    pH: float,
    acid_totals: Dict[str, float],
    Ka_vfa_eff: Dict[str, Any],

    CT_TIC: float, Ka1_eff: float, Ka2_eff: float,
    CT_NH_T: float, KaN_eff: float,
    CT_P: float, KaP1_eff: float, KaP2_eff: float, KaP3_eff: float,
    CT_SO4: float, KaS_eff: float,

    Kw: float, gamma_H: float, gamma_OH: float,

    strong_ions: Dict[str, float],
) -> Dict[str, Any]:
    """
    CHANGE (Update 5):
    - Supports monoprotic acids (Ka list length 1) with legacy keys:
        {name}_HA and {name}_A-
    - Supports polyprotic acids (Ka list length > 1) and outputs:
        {name}_H{n}A, {name}_H{n-1}A-, ..., {name}_A---...
      Example for citric (n=3):
        CitricAcid_H3A
        CitricAcid_H2A-
        CitricAcid_HA--
        CitricAcid_A---
    """
    def _polyprotic_alphas(H: float, Kas: list[float]) -> list[float]:
        n = len(Kas)
        prods = [1.0]
        for k in Kas:
            prods.append(prods[-1] * float(k))

        denom = 0.0
        for i in range(n + 1):
            denom += prods[i] * (H ** (n - i))

        if denom <= 0.0:
            return [1.0] + [0.0] * n

        alphas = []
        for i in range(n + 1):
            alphas.append((prods[i] * (H ** (n - i))) / denom)
        return alphas

    def _poly_key(name: str, n: int, i: int) -> str:
        """
        i = 0..n, where charge is -i and remaining H count is (n-i)
        """
        remaining_H = n - i
        if remaining_H >= 2:
            base = f"H{remaining_H}A"
        elif remaining_H == 1:
            base = "HA"
        else:
            base = "A"
        suffix = "-" * i
        return f"{name}_{base}{suffix}"

    H = 10.0 ** (-float(pH))
    # Water (activity corrected if gammas != 1)
    OH = Kw / (gamma_H * gamma_OH * H)

    out: Dict[str, Any] = {
        "H+": float(H),
        "OH-": float(OH),
    }

    # Acids (monoprotic or polyprotic)
    for name, CT in acid_totals.items():
        CT = _safe_float(CT)
        Kas = Ka_vfa_eff.get(name, None)
        
        if Kas is None:
            raise KeyError(f"Ka_vfa_eff missing for acid '{name}'. Check acid_pKas/effective_constants.")

        # Accept either float (legacy) or list
        if isinstance(Kas, (float, int)):
            Kas_list = [float(Kas)]
        else:
            Kas_list = [float(x) for x in Kas]

        if not Kas_list:
            # no Ka info => skip
            continue

        if len(Kas_list) == 1:
            Ka_eff = float(Kas_list[0])
            A_minus = CT * Ka_eff / (Ka_eff + H) if (Ka_eff + H) else 0.0
            HA = CT - A_minus
            out[f"{name}_HA"] = float(HA)
            out[f"{name}_A-"] = float(A_minus)
        else:
            alphas = _polyprotic_alphas(H, Kas_list)
            n = len(Kas_list)
            for i in range(n + 1):
                out[_poly_key(name, n, i)] = float(CT * alphas[i])

    # Carbonate ladder
    denom = H * H + Ka1_eff * H + Ka1_eff * Ka2_eff
    if denom:
        CO2 = CT_TIC * (H * H) / denom
        HCO3 = CT_TIC * (Ka1_eff * H) / denom
        CO3 = CT_TIC * (Ka1_eff * Ka2_eff) / denom
    else:
        CO2 = HCO3 = CO3 = 0.0

    out["CO2aq"] = float(CO2)
    out["HCO3-"] = float(HCO3)
    out["CO3--"] = float(CO3)

    # TAN
    denomN = KaN_eff + H
    NH4 = CT_NH_T * H / denomN if denomN else 0.0
    NH3 = CT_NH_T - NH4
    out["NH4+"] = float(NH4)
    out["NH3"] = float(NH3)

    # Phosphate ladder
    denomP = H**3 + KaP1_eff * H**2 + KaP1_eff * KaP2_eff * H + KaP1_eff * KaP2_eff * KaP3_eff
    if denomP:
        H3PO4 = CT_P * (H**3) / denomP
        H2PO4 = CT_P * (KaP1_eff * H**2) / denomP
        HPO4 = CT_P * (KaP1_eff * KaP2_eff * H) / denomP
        PO4 = CT_P * (KaP1_eff * KaP2_eff * KaP3_eff) / denomP
    else:
        H3PO4 = H2PO4 = HPO4 = PO4 = 0.0

    out["H3PO4"] = float(H3PO4)
    out["H2PO4-"] = float(H2PO4)
    out["HPO4--"] = float(HPO4)
    out["PO4---"] = float(PO4)

    # Bisulfate
    denomS = KaS_eff + H
    HSO4 = CT_SO4 * H / denomS if denomS else 0.0
    SO4 = CT_SO4 - HSO4
    out["HSO4-"] = float(HSO4)
    out["SO4--"] = float(SO4)

    # Echo strong ions (fully dissociated assumption)
    out["K+"] = float(strong_ions.get("CT_K", 0.0))
    out["Na+"] = float(strong_ions.get("CT_Na", 0.0))
    out["Cl-"] = float(strong_ions.get("CT_Cl", 0.0))
    out["NO3-"] = float(strong_ions.get("CT_NO3", 0.0))

    out["Mg++"] = float(strong_ions.get("CT_Mg", 0.0))
    out["Ca++"] = float(strong_ions.get("CT_Ca", 0.0))
    out["Zn++"] = float(strong_ions.get("CT_Zn", 0.0))
    out["Mn++"] = float(strong_ions.get("CT_Mn", 0.0))
    out["Co++"] = float(strong_ions.get("CT_Co", 0.0))

    out["Mo7O24------"] = float(strong_ions.get("CT_Mo7O24", 0.0))

    return out


# ════════════════════════════════════════════════════════════════════════
#  Generalised solver using EquilibriumSet
# ════════════════════════════════════════════════════════════════════════

def solve_from_equilibrium_set(
    *,
    equilibrium_set,
    concentrations: Dict[str, float],
    strong_ions: Dict[str, float] = None,
    T_K: float = 308.15,
    activity_model: "ActivityModel" = None,
    pH_min: float = -0.5,
    pH_max: float = 20.0,
    n_scan: int = 200,
    tol: float = 1e-12,
    max_outer: int = 20,
    I_init: float = 0.01,
    I_tol: float = 1e-8,
    damping: float = 0.6,
    logH_guess: Optional[float] = None,
) -> Dict[str, Any]:
    """Solve the acid-base charge balance using an EquilibriumSet.

    This is the generalised solver that replaces the hardcoded
    equilibrium systems in :func:`solve_acid_base`.  All equilibria
    are read from the provided :class:`EquilibriumSet`, using the
    generalised polyprotic formula.

    Parameters
    ----------
    equilibrium_set : EquilibriumSet
        Defines which equilibria participate in the charge balance,
        their pKa values (already temperature-corrected), categories,
        and ``n_active`` parameters.
    concentrations : dict
        ``{total_key: mol/L}`` — total concentrations for each
        equilibrium, keyed by the ``total_key`` from each
        :class:`EquilibriumDef`.  Missing keys are treated as zero.
    strong_ions : dict, optional
        ``{name: mol/L}`` for fully dissociated ions.  Recognised
        keys: ``CT_Na``, ``CT_K``, ``CT_Cl``, ``CT_NO3``, ``CT_Mg``,
        ``CT_Ca``, ``CT_Zn``, ``CT_Mn``, ``CT_Co``, ``CT_Mo7O24``.
    T_K : float
        Operating temperature (K).  Used for activity coefficients.
        pKa values in the EquilibriumSet should already be corrected
        to this temperature (the solver does NOT apply van 't Hoff).
    activity_model : ActivityModel
        Activity coefficient model.  Use ``make_activity_model(False, "ideal")``
        for ideal solution.
    pH_min, pH_max, n_scan, tol : float/int
        Root-finding parameters.
    logH_guess : float, optional
        Warm-start (log10[H+]).

    Returns
    -------
    dict
        Same format as :func:`solve_acid_base`: ``pH``, ``H+``, ``OH-``,
        species concentrations, ``IonicStrength``, ``logH``, etc.
    """
    from .activity import ionic_strength_from_speciation
    from .activity_models import make_activity_model as _make_am

    if activity_model is None:
        activity_model = _make_am(False, "ideal")
    if strong_ions is None:
        strong_ions = {}

    eq_set = equilibrium_set
    Kw = eq_set.water.Kw_at_T(T_K)

    # Build list of (EquilibriumDef, concentration, Kas_at_T) tuples
    eq_data = []
    for eq_def in eq_set:
        CT = _safe_float(concentrations.get(eq_def.total_key, 0.0))
        # pKas at operating T — use all pKas for species computation,
        # but only n_active for the charge balance.
        all_Kas = [10.0 ** (-pk) for pk in eq_def.pKas_at_T(T_K)]
        active_Kas = all_Kas[:eq_def.n_active]
        eq_data.append((eq_def, CT, all_Kas, active_Kas))

    # Parse strong ions
    si = {k: _safe_float(v) for k, v in strong_ions.items()}

    # ── Warm-start bracket ────────────────────────────────────────
    if logH_guess is not None and np.isfinite(logH_guess):
        pH_guess = -float(logH_guess)
        pH_lo = max(float(pH_min), pH_guess - 2.0)
        pH_hi = min(float(pH_max), pH_guess + 2.0)
    else:
        pH_lo, pH_hi = float(pH_min), float(pH_max)

    # ── Activity gamma wrapper ────────────────────────────────────
    def _gamma(z: int, I_val: float) -> float:
        I_val = max(0.0, float(I_val))
        try:
            return float(activity_model.gamma(z, I_val, T_K=float(T_K)))
        except TypeError:
            pass
        try:
            return float(activity_model.gamma(z, I_val, float(T_K)))
        except TypeError:
            pass
        return float(activity_model.gamma(z, I_val))

    # ── Charge balance residual ───────────────────────────────────
    def residual(pH: float, gamma_H: float = 1.0, gamma_OH: float = 1.0,
                 Ka_scale_fn=None) -> float:
        """Electroneutrality residual: cations − anions.

        Ka_scale_fn: if provided, called as Ka_scale_fn(charge_i) to
        get the activity correction factor for each dissociation step.
        """
        H = 10.0 ** (-float(pH))
        OH = Kw / (gamma_H * gamma_OH * H) if H > 0 else 0.0

        cations = H
        anions = OH

        for eq_def, CT, all_Kas, active_Kas in eq_data:
            if CT <= 0:
                continue

            if eq_def.category == "cation_acid":
                # BH⁺ → B + H⁺ : protonated form is a cation
                # Use only active Kas (normally 1 for NH4)
                if active_Kas:
                    Ka = active_Kas[0]
                    denom = Ka + H
                    if denom > 0:
                        cations += CT * H / denom
            elif eq_def.category in ("acid", "inorganic_acid"):
                # Polyprotic acid: use n_active dissociation steps
                # Both denominator and charge sum use only active Kas
                if active_Kas:
                    cations_contrib, anions_contrib = _generalised_acid_charge(
                        H, CT, active_Kas)
                    anions += anions_contrib
            # strong_ion category: handled separately below

        # Strong ions
        cations += _safe_float(si.get("CT_K", 0.0))
        cations += _safe_float(si.get("CT_Na", 0.0))
        cations += _safe_float(si.get("CT_cation", 0.0))
        cations += 2.0 * _safe_float(si.get("CT_Mg", 0.0))
        cations += 2.0 * _safe_float(si.get("CT_Ca", 0.0))
        cations += 2.0 * _safe_float(si.get("CT_Zn", 0.0))
        cations += 2.0 * _safe_float(si.get("CT_Mn", 0.0))
        cations += 2.0 * _safe_float(si.get("CT_Co", 0.0))

        anions += _safe_float(si.get("CT_Cl", 0.0))
        anions += _safe_float(si.get("CT_NO3", 0.0))
        anions += _safe_float(si.get("CT_anion", 0.0))
        anions += 6.0 * _safe_float(si.get("CT_Mo7O24", 0.0))

        return float(cations - anions)

    # ── Root finder ───────────────────────────────────────────────
    def solve_pH(res_fn, pH_lo_: float, pH_hi_: float) -> float:
        a, b = float(pH_lo_), float(pH_hi_)
        fa, fb = float(res_fn(a)), float(res_fn(b))
        if np.isfinite(fa) and np.isfinite(fb) and np.sign(fa) != np.sign(fb):
            return float(brentq(res_fn, a, b, xtol=tol))
        xs = np.linspace(float(pH_min), float(pH_max), int(n_scan))
        fs = np.array([res_fn(float(x)) for x in xs], dtype=float)
        sgn = np.sign(fs)
        idx = np.where(sgn[:-1] * sgn[1:] < 0)[0]
        if len(idx):
            return float(brentq(res_fn, float(xs[idx[0]]),
                                float(xs[idx[0] + 1]), xtol=tol))
        return float(xs[int(np.argmin(np.abs(fs)))])

    # ── Ideal fast path ───────────────────────────────────────────
    is_ideal = getattr(activity_model, "name", "ideal") == "ideal"

    if is_ideal:
        pH_sol = solve_pH(residual, pH_lo, pH_hi)
    else:
        # Davies/activity outer iteration
        I = float(I_init)
        converged = False
        pH_sol = 7.0
        for _ in range(int(max_outer)):
            gamma_H = _gamma(+1, I)
            gamma_OH = _gamma(-1, I)
            res_fn = lambda pH: residual(pH, gamma_H, gamma_OH)
            pH_sol = solve_pH(res_fn, pH_lo, pH_hi)
            sp_tmp = _compute_species_eq(pH_sol, eq_data, Kw, 1.0, 1.0, si)
            I_new = float(ionic_strength_from_speciation(sp_tmp))
            if abs(I_new - I) < float(I_tol):
                I = I_new
                converged = True
                break
            I = float(damping) * I + (1.0 - float(damping)) * I_new

        if not converged:
            import warnings
            warnings.warn(
                f"Activity iteration did not converge after {max_outer} iters",
                RuntimeWarning)

    # ── Compute full species output ───────────────────────────────
    gamma_H = _gamma(+1, 0.0) if is_ideal else _gamma(+1, I)
    gamma_OH = _gamma(-1, 0.0) if is_ideal else _gamma(-1, I)
    sp = _compute_species_eq(pH_sol, eq_data, Kw, gamma_H, gamma_OH, si)

    H = float(sp["H+"])
    pH_conc = -np.log10(max(H, 1e-30))
    sp["pH_conc"] = float(pH_conc)
    sp["aH"] = float(gamma_H * H)
    sp["gamma_H"] = float(gamma_H)
    sp["gamma_OH"] = float(gamma_OH)
    sp["IonicStrength"] = float(ionic_strength_from_speciation(sp))

    if is_ideal:
        sp["pH"] = float(pH_conc)
    else:
        sp["pH"] = float(-np.log10(max(gamma_H * H, 1e-30)))

    sp["logH"] = float(np.log10(max(H, 1e-30)))
    return sp


def _generalised_acid_charge(H, CT, Kas):
    """Compute cation and anion contributions from a polyprotic acid.

    Uses only the provided Kas (which may be a subset — the 'active'
    steps).  Both the denominator and the charge summation use the
    same set of Kas, treating the system as an n_active-protic acid.

    Returns
    -------
    (cations, anions) : tuple of float
        For standard acids, cations=0.  Separated for future extensions.
    """
    n = len(Kas)
    prods = [1.0]
    for Ka in Kas:
        prods.append(prods[-1] * Ka)

    denom = sum(prods[i] * H ** (n - i) for i in range(n + 1))
    if denom <= 0:
        return 0.0, 0.0

    anion_charge = 0.0
    for i in range(1, n + 1):
        alpha_i = prods[i] * H ** (n - i) / denom
        anion_charge += i * CT * alpha_i

    return 0.0, anion_charge


# _CANONICAL_NAMES was deleted in chemistry-unification-3b.
# All recognised acid systems (CO₂, NH₄⁺/NH₃, phosphate, bisulfate)
# now carry ``species_refs`` on their ``EquilibriumDef`` so
# ``_compute_species_eq`` emits by ``Species.id`` directly.
# The only remaining fallback in ``_compute_species_eq`` is the
# deprecated generic ``{name}_HA``/``{name}_A-`` path for string-based
# entries without ``species_refs`` (legacy VFA rows in
# ``EquilibriumSet.bsm2_default()``).  That path is removed in the
# PARTITION_MODEL phase.


def _compute_species_eq(pH, eq_data, Kw, gamma_H, gamma_OH, strong_ions):
    """Compute all species concentrations from the solved pH.

    Uses ALL pKas (not just n_active) for species output, since all
    species exist physically even if some are excluded from the charge
    balance.

    Emission priority (chemistry-unification-3b):

    1. If ``eq_def.species_refs`` is non-empty, emit using the declared
       ``Species.id`` for each ladder position — no synthesised suffixes.
       Example: ``CO2`` Species emits ``out["CO2"]``; ``HCO3-`` Species
       emits ``out["HCO3-"]``.

    2. Otherwise fall back to the legacy paths:

       - Recognised names (entries in :data:`_CANONICAL_NAMES`) emit
         canonical keys (``"CO2aq"``, ``"NH4+"``, etc.).
       - Unrecognised acids (VFAs, custom systems) emit generic
         ``{name}_HA`` / ``{name}_A-`` / ``{name}_BH+`` / ``{name}_B``
         keys via :func:`_species_key`.  These are **deprecated** —
         string-based ``EquilibriumDef`` entries without ``species_refs``
         (e.g. VFA rows in ``EquilibriumSet.bsm2_default()``) still use
         this path. It will be removed in the PARTITION_MODEL phase once
         those entries gain ``species_refs``.
    """
    H = 10.0 ** (-float(pH))
    OH = Kw / (gamma_H * gamma_OH * H) if H > 0 else 0.0

    out = {"H+": float(H), "OH-": float(OH)}

    for eq_def, CT, all_Kas, active_Kas in eq_data:
        refs = eq_def.species_refs  # non-empty → Species-ID emission

        if eq_def.category == "cation_acid":
            # BH⁺ / B system (e.g. NH4+/NH3)
            if not all_Kas:
                continue
            Ka = all_Kas[0]
            denom = Ka + H
            BH = CT * H / denom if denom > 0 else CT
            B = CT - BH
            if refs:
                out[refs[0].id] = float(BH)
                if len(refs) > 1:
                    out[refs[1].id] = float(B)
            else:
                # Deprecated legacy path: string-based entry without
                # species_refs (VFA rows in EquilibriumSet.bsm2_default()).
                out[f"{eq_def.name}_BH+"] = float(BH)
                out[f"{eq_def.name}_B"] = float(B)

        elif eq_def.category in ("acid", "inorganic_acid"):
            # Polyprotic acid: compute all species using all_Kas
            n = len(all_Kas)
            alphas = _polyprotic_alpha(H, all_Kas)

            if refs:
                # Species-ID emission: emit min(n+1, len(refs)) entries.
                for i in range(min(n + 1, len(refs))):
                    out[refs[i].id] = float(CT * alphas[i])
            else:
                # Deprecated legacy path: string-based entry without
                # species_refs (VFA rows in EquilibriumSet.bsm2_default()).
                for i in range(n + 1):
                    key = _species_key(eq_def.name, n, i)
                    out[key] = float(CT * alphas[i])

    # Strong ions (echoed)
    for key, charge_map in [
        ("CT_K", ("K+", 1)), ("CT_Na", ("Na+", 1)),
        ("CT_cation", ("Cation(inert)", 1)),
        ("CT_Cl", ("Cl-", 1)), ("CT_NO3", ("NO3-", 1)),
        ("CT_anion", ("Anion(inert)", 1)),
        ("CT_Mg", ("Mg++", 1)), ("CT_Ca", ("Ca++", 1)),
        ("CT_Zn", ("Zn++", 1)), ("CT_Mn", ("Mn++", 1)),
        ("CT_Co", ("Co++", 1)), ("CT_Mo7O24", ("Mo7O24------", 1)),
    ]:
        val = _safe_float(strong_ions.get(key, 0.0))
        out[charge_map[0]] = float(val)

    return out


def _species_key(name, n_protons, i):
    """Generate output key for species i of an n-protic acid.

    i=0 is the fully protonated form (neutral), i=n is fully deprotonated.
    """
    remaining_H = n_protons - i
    if remaining_H == 0:
        base = "A"
    elif remaining_H == 1:
        base = "HA"
    else:
        base = f"H{remaining_H}A"
    suffix = "-" * i
    return f"{name}_{base}{suffix}"

