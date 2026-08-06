# -*- coding: utf-8 -*-
"""
test_run_standards_suite_with_secondary.py

Unified standards/benchmarks suite + parity plot (predicted vs "standard"/nominal pH).

Primary / NIST/NBS-style cases included:
- NIST carbonate standard (0.025 M NaHCO3 + 0.025 M Na2CO3)  pH(S)=10.012 @ 25°C
- NIST SRM 186 phosphate standards:
    (A) 0.025 m KH2PO4 + 0.025 m Na2HPO4   pH(S)=6.8640 @ 25°C
    (B) 0.008695 m KH2PO4 + 0.03043 m Na2HPO4  pH(S)=7.4157 @ 25°C
- NBS/NIST citrate standard:
    0.05 m potassium dihydrogen citrate (KH2Cit) pH(S)=3.776 @ 25°C
- NBS operational high-pH standard:
    0.01 m NaOH + 0.10 m KCl  pH(S)≈12.45 @ 25°C  (operational)

Secondary / non-certified examples added:
- Acetic acid / sodium acetate (constructed via Henderson-Hasselbalch target pH at 25°C)
- McIlvaine citrate–phosphate buffer examples (constructed from the classic mixing table:
  0.2 M Na2HPO4 + 0.1 M citric acid; total volume 20 mL)

Notes:
- Molality (m) is approximated as molarity (M) here for NIST recipes; this can introduce small offsets.
- Secondary examples are NOT certified standards. They are useful for trend/behavior testing.
"""

import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

from speciation.engine import BisectionChemicalEquilibriumEngine


# -----------------------------
# Helpers to construct cases
# -----------------------------
def acetate_buffer_case(
    *,
    name: str,
    pH_target: float,
    C_total: float = 0.100,  # total acetate = [HAc] + [Ac-] (mol/L)
    pKa: float = 4.76,       # acetic acid pKa at ~25°C (common value)
    T_C: float = 25.0,
    tol: float = 0.10,
) -> dict:
    """
    Construct an acetic acid / acetate buffer with a chosen pH target using HH:
        pH = pKa + log10([Ac-]/[HAc])

    We choose [Ac-] and [HAc] from C_total and pH_target.
    Sodium counterion is set to [Ac-] (as if from sodium acetate).
    """
    r = 10 ** (float(pH_target) - float(pKa))  # [Ac-]/[HAc]
    HAc = float(C_total) / (1.0 + r)
    Ac = float(C_total) - HAc

    return {
        "name": name,
        "T_C": float(T_C),
        "pH_standard": float(pH_target),  # nominal (constructed)
        "tol": float(tol),
        "category": "secondary",
        # Represent acetate as a generic monoprotic acid system:
        "acid_totals": {"ACET": float(C_total)},
        "acid_pKas": {"ACET": [float(pKa)]},
        # Strong ion to balance acetate base supplied as sodium acetate:
        "CT_Na": float(Ac),
        # (Optional) For clarity in debug/printing:
        "_meta": {"HAc": HAc, "Ac-": Ac, "ratio_Ac_to_HAc": r, "pKa": pKa},
    }


def mcilvaine_case_from_table_20mL(
    *,
    pH_label: float,
    V_na2hpo4_mL: float,
    V_citric_mL: float,
    T_C: float = 25.0,
    tol: float = 0.12,
    citric_pKas=(3.13, 4.76, 6.40),
) -> dict:
    """
    Build McIlvaine citrate–phosphate buffer from the classic table:
        Mix Vp mL of 0.2 M Na2HPO4 + Vc mL of 0.1 M citric acid to make 20 mL total.

    Convert to final molarities:
        CT_P = (0.2 * Vp/1000) / 0.02
        CT_CIT = (0.1 * Vc/1000) / 0.02
        CT_Na = 2 * (0.2 * Vp/1000) / 0.02  (since Na2HPO4 has 2 Na per phosphate)
    """
    Vp_L = float(V_na2hpo4_mL) / 1000.0
    Vc_L = float(V_citric_mL) / 1000.0
    Vtot_L = 0.020

    moles_P = 0.2 * Vp_L
    moles_CIT = 0.1 * Vc_L

    CT_P = moles_P / Vtot_L
    CT_CIT = moles_CIT / Vtot_L
    CT_Na = (2.0 * moles_P) / Vtot_L

    return {
        "name": f"McIlvaine citrate–phosphate (table) pH {pH_label:.1f}",
        "T_C": float(T_C),
        "pH_standard": float(pH_label),  # table label (nominal)
        "tol": float(tol),
        "category": "secondary",
        "CT_P": float(CT_P),
        "CT_Na": float(CT_Na),
        "acid_totals": {"CIT": float(CT_CIT)},
        "acid_pKas": {"CIT": [float(x) for x in citric_pKas]},
        "_meta": {
            "Vp_mL_0.2M_Na2HPO4": V_na2hpo4_mL,
            "Vc_mL_0.1M_Citric": V_citric_mL,
            "CT_CIT": CT_CIT,
            "CT_P": CT_P,
            "CT_Na": CT_Na,
        },
    }


# -----------------------------
# Solver wrapper
# -----------------------------
def solve_case(engine: BisectionChemicalEquilibriumEngine, case: dict) -> dict:
    out = engine.solve(
        acid_totals=dict(case.get("acid_totals", {})),
        acid_pKas=dict(case.get("acid_pKas", {})),
        CT_TIC=float(case.get("CT_TIC", 0.0)),
        CT_P=float(case.get("CT_P", 0.0)),
        CT_NH_T=float(case.get("CT_NH_T", 0.0)),
        CT_Na=float(case.get("CT_Na", 0.0)),
        CT_Cl=float(case.get("CT_Cl", 0.0)),
        CT_K=float(case.get("CT_K", 0.0)),
        CT_NO3=float(case.get("CT_NO3", 0.0)),
        CT_SO4=float(case.get("CT_SO4", 0.0)),
        CT_Mg=float(case.get("CT_Mg", 0.0)),
        CT_Ca=float(case.get("CT_Ca", 0.0)),
        CT_Zn=float(case.get("CT_Zn", 0.0)),
        CT_Mn=float(case.get("CT_Mn", 0.0)),
        CT_Co=float(case.get("CT_Co", 0.0)),
        CT_Mo7O24=float(case.get("CT_Mo7O24", 0.0)),
        pH_min=0.0,
        pH_max=14.0,
        n_scan=int(case.get("n_scan", 900)),
        tol=float(case.get("solver_tol", 1e-12)),
    )

    pred_pH = float(out.get("pH", np.nan))
    std_pH = float(case["pH_standard"])
    tol = float(case.get("tol", 0.05))
    ok = abs(pred_pH - std_pH) <= tol

    inputs = {k: case.get(k, 0.0) for k in [
        "CT_TIC","CT_P","CT_NH_T","CT_Na","CT_Cl","CT_K","CT_NO3","CT_SO4","CT_Mg","CT_Ca","CT_Zn","CT_Mn","CT_Co","CT_Mo7O24"
    ]}

    return {
        "name": case["name"],
        "category": case.get("category", "primary"),
        "T_C": float(case.get("T_C", np.nan)),
        "pH_standard": std_pH,
        "pH_pred": pred_pH,
        "abs_err": abs(pred_pH - std_pH),
        "err": pred_pH - std_pH,
        "pass": bool(ok),
        "tol": tol,
        "I": float(out.get("IonicStrength", np.nan)),
        "out": out,
        "inputs": inputs,
        "_meta": case.get("_meta", {}),
    }


def print_case(rec: dict, debug_species: bool = True) -> None:
    print(f"\n[{rec['category']}] {rec['name']}")
    print(f"  T = {rec['T_C']:.2f} °C")
    print("  Inputs:", ", ".join([f"{k}={v:.6g}" for k, v in rec["inputs"].items() if float(v) != 0.0]) or "(all zeros except acids)")
    if rec.get("_meta"):
        # Keep this concise
        pass
    print(f"  Pred:   pH={rec['pH_pred']:.4f}, I={rec['I']:.6g} M")
    print(f"  Ref:    pH={rec['pH_standard']:.4f} (±{rec['tol']:.4f})  ->  {'PASS' if rec['pass'] else 'FAIL'}")
    print(f"  Error:  {rec['err']:+.4f}")

    if debug_species:
        out = rec["out"]
        keys = [
            "H+", "OH-", "aH",
            "CO2aq", "HCO3-", "CO3--",
            "H3PO4", "H2PO4-", "HPO4--", "PO4---",
            "H3Cit", "H2Cit-", "HCit2-", "Cit3-",
            "Na+", "K+", "Cl-",
        ]
        for k in keys:
            if k in out:
                print(f"    {k:8s}: {out[k]}")


if __name__ == "__main__":
    # -----------------------------
    # Engine settings
    # -----------------------------
    LEVEL = 2
    USE_ACTIVITY = False
    ACTIVITY_MODEL = "davies"
    T_C = 25.0

    engine = BisectionChemicalEquilibriumEngine(
        level=LEVEL,
        use_activity=USE_ACTIVITY,
        activity_model=ACTIVITY_MODEL,
        T_C=T_C,
    )

    # -----------------------------
    # Standards / examples in one consolidated section
    # -----------------------------
    CIT_PKA_25C = [3.13, 4.76, 6.40]

    CASES = [
        # --- NIST/NBS-style cases ---
        {
            "name": "NIST Carbonate Standard: 0.025 M NaHCO3 + 0.025 M Na2CO3",
            "category": "primary",
            "T_C": 25.0,
            "pH_standard": 10.012,
            "tol": 0.05,
            "CT_TIC": 0.050,
            "CT_Na": 0.075,
        },
        {
            "name": "NIST Phosphate (SRM 186): 0.025 m KH2PO4 + 0.025 m Na2HPO4",
            "category": "primary",
            "T_C": 25.0,
            "pH_standard": 6.8640,
            "tol": 0.03,
            "CT_P": 0.0500,
            "CT_K": 0.0250,
            "CT_Na": 0.0500,
        },
        {
            "name": "NIST Phosphate (SRM 186): 0.008695 m KH2PO4 + 0.03043 m Na2HPO4",
            "category": "primary",
            "T_C": 25.0,
            "pH_standard": 7.4157,
            "tol": 0.03,
            "CT_P": (0.008695 + 0.03043),
            "CT_K": 0.008695,
            "CT_Na": (2.0 * 0.03043),
        },
        {
            "name": "NBS/NIST Citrate Standard: 0.05 m potassium dihydrogen citrate (KH2Cit)",
            "category": "primary",
            "T_C": 25.0,
            "pH_standard": 3.776,
            "tol": 0.05,
            "CT_K": 0.0500,
            "acid_totals": {"CIT": 0.0500},
            "acid_pKas": {"CIT": CIT_PKA_25C},
        },
        {
            "name": "NBS Operational High-pH: 0.01 m NaOH + 0.10 m KCl",
            "category": "operational",
            "T_C": 25.0,
            "pH_standard": 12.45,
            "tol": 0.15,
            "CT_Na": 0.0100,   # NaOH
            "CT_K": 0.1000,    # KCl
            "CT_Cl": 0.1000,
        },

        # --- Acetate buffer examples (constructed) ---
        acetate_buffer_case(
            name="Acetate buffer example: 0.10 M total acetate targeting pH 4.76 (~pKa)",
            pH_target=4.76, C_total=0.100, pKa=4.76, tol=0.10, T_C=25.0
        ),
        acetate_buffer_case(
            name="Acetate buffer example: 0.10 M total acetate targeting pH 5.20",
            pH_target=5.20, C_total=0.100, pKa=4.76, tol=0.12, T_C=25.0
        ),

        # --- McIlvaine citrate–phosphate examples (from classic 20 mL mixing table) ---
        # Table values: pH -> (mL 0.2M Na2HPO4, mL 0.1M citric acid) for 20 mL total
        mcilvaine_case_from_table_20mL(pH_label=3.0, V_na2hpo4_mL=4.11, V_citric_mL=15.89, tol=0.12, citric_pKas=CIT_PKA_25C),
        mcilvaine_case_from_table_20mL(pH_label=4.0, V_na2hpo4_mL=7.71, V_citric_mL=12.29, tol=0.12, citric_pKas=CIT_PKA_25C),
        mcilvaine_case_from_table_20mL(pH_label=5.0, V_na2hpo4_mL=10.30, V_citric_mL=9.70,  tol=0.12, citric_pKas=CIT_PKA_25C),
        mcilvaine_case_from_table_20mL(pH_label=6.0, V_na2hpo4_mL=12.63, V_citric_mL=7.37,  tol=0.12, citric_pKas=CIT_PKA_25C),
        mcilvaine_case_from_table_20mL(pH_label=7.0, V_na2hpo4_mL=16.47, V_citric_mL=3.53,  tol=0.12, citric_pKas=CIT_PKA_25C),
    ]

    # -----------------------------
    # Run suite
    # -----------------------------
    results = []
    for case in CASES:
        rec = solve_case(engine, case)
        results.append(rec)
        print_case(rec, debug_species=True)

    # Summary table
    print("\nSummary")
    print("{:<18s} {:<70s} {:>10s} {:>10s} {:>10s} {:>7s}".format("Category", "Case", "pH_ref", "pH_pred", "err", "PASS"))
    print("-" * 128)
    for r in results:
        print("{:<18s} {:<70s} {:>10.4f} {:>10.4f} {:>+10.4f} {:>7s}".format(
            r["category"][:18],
            r["name"][:70],
            r["pH_standard"],
            r["pH_pred"],
            r["err"],
            "YES" if r["pass"] else "NO",
        ))

    # -----------------------------
    # Parity plot: predicted vs reference (standards + nominal examples)
    # -----------------------------
    try:
        import matplotlib.pyplot as plt

        x = np.array([r["pH_standard"] for r in results], dtype=float)
        y = np.array([r["pH_pred"] for r in results], dtype=float)
        labels = [r["name"] for r in results]
        cats = [r["category"] for r in results]

        fig, ax = plt.subplots()

        lo = float(np.nanmin([x.min(), y.min()])) - 0.25
        hi = float(np.nanmax([x.max(), y.max()])) + 0.25

        # parity line
        ax.plot([lo, hi], [lo, hi], linestyle="-")

        # Plot by category (no explicit colors; differentiate with markers)
        markers = {"primary": "o", "operational": "s", "secondary": "^"}
        for cat in ["primary", "operational", "secondary"]:
            idx = [i for i, c in enumerate(cats) if c == cat]
            if not idx:
                continue
            ax.scatter(x[idx], y[idx], marker=markers.get(cat, "o"), label=cat)

        # annotate (may be dense; comment out if you prefer)
        for xi, yi, lab in zip(x, y, labels):
            ax.annotate(lab, (xi, yi), textcoords="offset points", xytext=(6, 4), fontsize=8)

        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel("Reference pH (certified or nominal)")
        ax.set_ylabel("Predicted pH (BisectionChemicalEquilibriumEngine)")
        ax.set_title("Parity Plot: Standards + Example Buffers (pH)")
        ax.set_aspect("equal", adjustable="box")
        ax.legend()
        fig.tight_layout()
        plt.show()

    except Exception as e:
        print("\n(Parity plot skipped:", e, ")")
