# -*- coding: utf-8 -*-
"""Element and charge residual of every reaction in the repo's models.

Run by hand from the repo root; pytest does not collect it::

    py tests/performance/balance_audit.py
    py tests/performance/balance_audit.py --only bsm2 adm1 --out balance.json

For each model, every reaction (kinetic and equilibrium) is checked over every element
its species carry and over charge, whatever its ``balance_elements`` says: the residual
of an element is ``sum(coefficient * atoms[element])`` over the reaction's entries, and
of charge ``sum(coefficient * charge)``. A reaction closes when every residual is within
``--atol`` (default 1e-9). The report lists, per model, how many reactions there are and
each one that does not close, with its residuals. It is the measurement behind the
balance table in ``CONSERVATION_LEDGER_CHECKLIST.md`` and behind checkpoint 6.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _path in (REPO, REPO / "models", REPO / "tests" / "standalone", HERE):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))


def _bsm2():
    from vlmodels.adm1.bsm2 import build_bsm2_reactions
    return build_bsm2_reactions(verbose=False).reactions


def _adm1(ethanol=False):
    from vlmodels.adm1.base import build_adm1_reactions
    return build_adm1_reactions(ethanol=ethanol).reactions


def _tank(basis):
    from PyOMES.databases.anaerobic_digestion import AD_BASIC
    from PyOMES.templates.stirred_tank import StirredTankBuilder
    organism_kw = {"balance_basis": basis}
    if "N" in basis:
        organism_kw["n_source_id"] = "NH3"
    cv = (
        StirredTankBuilder()
        .vessel(V_total_L=10.0, T_K=305.15)
        .chemistry(chemistry_db=AD_BASIC)
        .organism("Yeast", **organism_kw)
        .substrate("AceticAcid", mu_max=0.5, Ks=5e-3, yield_gX_gS=0.36)
        .build()
    )
    rs = cv.reaction_system  # the template's is a single reaction, not a system
    return getattr(rs, "reactions", [rs])


def _database(name):
    from PyOMES.databases import anaerobic_digestion, aqueous, bioprocess_basic
    db = {"aqueous": aqueous.AQUEOUS_DEFAULT,
          "bioprocess_basic": bioprocess_basic.BIOPROCESS_BASIC,
          "anaerobic_digestion": anaerobic_digestion.AD_BASIC}[name]
    return db.reactions.reactions if db.reactions is not None else []


def _d2c():
    import ledger_fingerprint as lf
    sim, _, _ = lf.CASES["d2c"]()
    return sim.cvs["main"].reaction_system.reactions


def _user_defined():
    import test_user_defined_model as udm
    return udm._build_reactions().reactions


MODELS = {
    "bsm2": _bsm2,
    "adm1": _adm1,
    "adm1_ethanol": lambda: _adm1(ethanol=True),
    "tank_cho": lambda: _tank("CHO"),
    "tank_chno": lambda: _tank("CHNO"),
    "db_aqueous": lambda: _database("aqueous"),
    "db_bioprocess_basic": lambda: _database("bioprocess_basic"),
    "db_anaerobic_digestion": lambda: _database("anaerobic_digestion"),
    "d2c": _d2c,
    "user_defined": _user_defined,
}


def residuals(reaction):
    """{element: residual, ..., "charge": residual} over every element present."""
    out = {}
    charge = 0.0
    for entry in reaction.stoichiometry:
        sp = entry.species
        for element, count in sp.atoms.items():
            out[element] = out.get(element, 0.0) + entry.coefficient * float(count)
        charge += entry.coefficient * float(getattr(sp, "charge", 0) or 0)
    out["charge"] = charge
    return out


def audit(names, atol):
    report = {}
    for name in names:
        try:
            reactions = list(MODELS[name]())
        except Exception as exc:  # recorded, so one model that cannot build does not stop the rest
            report[name] = {"error": f"{type(exc).__name__}: {exc}"}
            continue
        failing = []
        for rxn in reactions:
            res = residuals(rxn)
            bad = {k: v for k, v in res.items() if abs(v) > atol}
            if bad:
                failing.append({"label": getattr(rxn, "label", "?"),
                                "kind": type(rxn).__name__,
                                "balance_elements": getattr(rxn, "balance_elements", None),
                                "residuals": bad})
        report[name] = {"reactions": len(reactions), "not_closing": failing}
    return report


def _print(report):
    for name, r in report.items():
        if "error" in r:
            print(f"{name}: cannot build ({r['error']})")
            continue
        print(f"{name}: {r['reactions']} reactions, {len(r['not_closing'])} do not close")
        for f in r["not_closing"]:
            res = ", ".join(f"{k} {v:+.4g}" for k, v in sorted(f["residuals"].items()))
            print(f"  {f['label']} ({f['kind']}): {res}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", nargs="+", choices=sorted(MODELS))
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument("--out")
    args = parser.parse_args()
    report = audit(args.only or list(MODELS), args.atol)
    _print(report)
    if args.out:
        Path(args.out).write_text(json.dumps(report, indent=1, default=str), encoding="utf-8")
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
