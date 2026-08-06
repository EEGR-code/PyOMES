# FBA Demo Improvements — Shipped 2026-06-01

## Status: Shipped 2026-06-01

Demo-quality work only — no `src/` or `models/` changes. Closes
item 3 of the 2026-05-29 session agenda (FBA stoichiometry
verification, carried forward from the 2026-05-28 exploration
session that first wrote the FBA demos).

---

## What changed

### `demos/model_api/chemistry/fba/fba_toy.py`

**Documentation fix only.**

The original file carried a `.. warning::` block with a "Future
work" sub-block directing the reader to verify the 7-reaction
network against Orth, Thiele & Palsson (2010) Box 1. Investigation
found that **Box 1 of that paper presents the S-matrix formalism
and steady-state constraint equations** (`S·v = 0`, notation for
stoichiometric coefficients, dimension of null-space) — it is not
a concrete toy reaction list. There is no published Box 1 network
to verify against. The "future work" disclaimer rested on a false
premise.

Changes:
- `.. warning::` → `.. note::`, "Future work" sub-block removed.
- Attribution updated to state clearly that the 7-reaction network
  is an original pedagogical construction; Orth et al 2010 is cited
  for the FBA framework and Mahadevan et al 2002 for the dFBA
  coupling pattern, not for this specific network.
- `main()` print statement corrected: `"Orth-Thiele-Palsson 2010
  network"` → `"original 7-reaction central-metabolism
  construction"`.

---

### `demos/model_api/chemistry/fba/ecoli_core.json`

**Network extended from 14 to 18 reactions to fix anaerobic growth.**

Option B was chosen (fix the existing hand-crafted model) over
Option A (replace with the canonical 95-reaction BiGG
`e_coli_core`). The hand-crafted model is now clearly documented as
an original pedagogical construction; the "UNVERIFIED / FUTURE
WORK" framing is removed.

**Key stoichiometry finding.** The original file carried a note
that adding PFL (pyruvate-formate lyase) would enable anaerobic
growth. This is **incorrect as stated**. In the lumped model, GLYC
produces 2 NADH per glucose; without O₂ the only NADH sink is LDH
(pyr + NADH → lac). The pyr_c and nadh_c steady-state equations
reduce to `2·v_PDH + 3·v_TCA + v_BIOMASS = 0`, forcing
`v_BIOMASS = 0` regardless of whether PFL is present — PFL is
NADH-neutral so it does not break this constraint.

The minimal fix is **PFL + ADH** together:

- **PFL** (`pyr_c → accoa_c + for_e`): NADH-free pyruvate route,
  decouples pyruvate from the LDH NADH-drain.
- **ADH** (`accoa_c + 2·nadh_c → etoh_e`): drains glycolytic NADH
  via acetyl-CoA → ethanol, without competing for pyruvate.

With both reactions present the anaerobic LP can achieve
`v_BIOMASS > 0`. Carbon balance: PFL (3C → 2C + 1C ✓), ADH
(2C → 2C ✓).

New reactions added (4): `PFL`, `ADH`, `EX_for_e`, `EX_etoh_e`.  
New metabolites added (2): `for_e` (formate, C₁H₂O₂), `etoh_e`
(ethanol, C₂H₆O).  
`cv_species_mapping` entries added for both exchange reactions.  
`_comment` block rewritten: documents the 18-reaction scope, lists
known limitations (lumped glycolysis/TCA, simplified biomass
equation, no PPP/anaplerotic/glyoxylate pathways, 18 vs 95
reactions in the canonical BiGG model), directs readers to BiGG
`e_coli_core` for a verified network.

---

### `demos/model_api/chemistry/fba/fba_ecoli_core.py`

- `.. warning::` → `.. note::`, "Future work" sub-block removed;
  attribution updated to match the JSON's reframing.
- `build()`: `"Formate": 0.0` and `"Ethanol": 0.0` added to
  initial `n_mol`.
- `main()` trajectory output: Formate and Ethanol columns added.
- `main()` closing note updated from "growth halts when O₂
  depletes" to describe the aerobic → mixed-acid fermentation shift.

**Smoke-test result:** biomass grows 9× (0.005 → 0.045 mol),
formate and ethanol accumulate in the anaerobic phase, lactate ≈ 0
(LP prefers PFL/ADH over LDH), CO₂ ≈ 0 (PFL produces formate not
CO₂; tiny initial O₂ supply means the aerobic PDH/TCA window is
negligible at 5 d.p.).

---

### README updates

- `demos/model_api/README.md`: table row for `fba_ecoli_core.py`
  updated — "~14-reaction / future replacement" →
  "18-reaction / original pedagogical construction / aerobic +
  anaerobic".
- `demos/README.md`: reading-order description updated —
  "attribution warnings / future replacement" →
  "attribution notes / original hand-crafted networks".

---

## What was not changed

- `fba_toy.py` reaction network: stoichiometry unchanged, as no
  published counterpart exists to verify against.
- `JSONFBASolver` loader class: reusable as-is; no `from_bigg_json()`
  constructor added (Option A deferred).
- No `src/` or `models/` edits.

## References

- Mahadevan, R., Edwards, J. S., & Doyle III, F. J. (2002). Dynamic
  flux balance analysis of diauxic growth in *Escherichia coli*.
  *Biophysical Journal*, 83(3), 1331–1340.
- Orth, J. D., Thiele, I., & Palsson, B. Ø. (2010). What is flux
  balance analysis? *Nature Biotechnology*, 28(3), 245–248.
- Orth, J. D., Fleming, R. M. T., & Palsson, B. Ø. (2011).
  Reconstruction and use of microbial metabolic networks: the core
  *Escherichia coli* metabolic model as an educational guide.
  *EcoSal Plus*, 4(1).
