"""Add NaOH/H2SO4 corrector species, pH 6 setpoint, updated prose."""
import json

path = r'c:\Users\k2473520\PyOMES\demos\model_api\D2Cworkshop\updated_layout\Example2_batch_fermenter.ipynb'
with open(path, 'r', encoding='utf-8') as f:
    nb = json.load(f)


def set_source(cell, text):
    lines = text.split('\n')
    cell['source'] = [l + '\n' for l in lines[:-1]] + ([lines[-1]] if lines[-1] else [])


def find(cid):
    for i, c in enumerate(nb['cells']):
        if c.get('id') == cid:
            return i, c
    raise KeyError(cid)


# ── code-params: PH_SETPOINT 5.0 → 6.0 ──────────────────────────────────────
_, c = find('code-params')
set_source(c, (
    "T_K            = 305.15\n"
    "V_TOTAL_L      = 2.0\n"
    "HEADSPACE_FRAC = 0.20\n"
    "V_GAS = V_TOTAL_L * HEADSPACE_FRAC\n"
    "V_LIQ = V_TOTAL_L * (1.0 - HEADSPACE_FRAC)\n"
    "\n"
    "MU_MAX      = 0.5    # 1/h\n"
    "KS_G_L      = 5e-3  # g/L\n"
    "YIELD       = 0.36   # g biomass / g substrate\n"
    "PH_SETPOINT = 6.0\n"
    "KLA_PER_H   = {\"O2\": 150.0, \"CO2\": 135.0, \"NH3\": 135.0}  # CO2/NH3 at 0.9 × O2 kLa\n"
    "\n"
    "S0_G_L     = 1.2    # g/L initial substrate\n"
    "X0_G_L     = 0.1    # g/L inoculum\n"
    "NH4CL_G_L  = 5.0    # g/L NH4Cl\n"
    "MW_NH4CL   = 53.491  # g/mol\n"
    "KH2PO4_G_L = 30.0   # g/L KH2PO4 (potassium dihydrogen phosphate)\n"
    "MW_KH2PO4  = 136.084  # g/mol\n"
    "\n"
    "TAU_H   = 5.0\n"
    "N_STEPS = 1000\n"
    "\n"
    "print(f\"V_liq = {V_LIQ:.1f} L   V_gas = {V_GAS:.1f} L   T = {T_K:.2f} K\")\n"
    "print(f\"kLa(O2) = {KLA_PER_H['O2']:.0f} /h   vvm = 1   pH setpoint = {PH_SETPOINT}\")\n"
    "print(f\"NH4Cl:  {NH4CL_G_L} g/L  →  {NH4CL_G_L/MW_NH4CL*1000:.1f} mmol/L NH4+ and Cl-\")\n"
    "print(f\"KH2PO4: {KH2PO4_G_L} g/L  →  {KH2PO4_G_L/MW_KH2PO4*1000:.1f} mmol/L total phosphate\")"
))

# ── md-db-config: table row pH setpoint 5.0 → 6.0 ───────────────────────────
_, c = find('md-db-config')
src = ''.join(c['source'])
set_source(c, src.replace(
    '| pH setpoint | 5.0 | Henderson–Hasselbalch near pKₐ |',
    '| pH setpoint | 6.0 | Phosphate buffer range; acetate 94 % dissociated |'
))

# ── md-db-ph: rewrite for pH 6.0 ─────────────────────────────────────────────
_, c = find('md-db-ph')
set_source(c, (
    "### pH setpoint selection\n"
    "\n"
    "Acetic acid has $\\mathrm{p}K_a = 4.756$. At the chosen setpoint pH 6.0,\n"
    "the Henderson–Hasselbalch equation gives the acid/base split:\n"
    "\n"
    "$$\n"
    "\\frac{[\\text{Acetate}^-]}{[\\text{AceticAcid}]} = 10^{\\mathrm{pH} - \\mathrm{p}K_a}\n"
    "= 10^{6.0 - 4.756} = 17.5\n"
    "$$\n"
    "\n"
    "So at pH 6.0, only **5.4 %** of total acetate remains in the undissociated\n"
    "AceticAcid form; **94.6 %** is Acetate⁻. Operating 1.24 units above the\n"
    "$\\mathrm{p}K_a$ is appropriate here because:\n"
    "\n"
    "- Membrane toxicity from the undissociated acid is minimised while a small\n"
    "  undissociated fraction is still available for direct membrane uptake.\n"
    "- The KH₂PO₄ buffer is most effective in the pH 6–8 range\n"
    "  (H₂PO₄⁻/HPO₄²⁻, $\\mathrm{p}K_a = 7.2$); at pH 6.0 phosphate provides\n"
    "  robust chemical buffering that reduces the controller dose required.\n"
    "\n"
    "The liquid is brought to this setpoint before the simulation begins by\n"
    "`ControlVolume.equilibrate_to_pH`, which adds NaOH (if base is needed) or\n"
    "H₂SO₄ (if acid is needed) in the minimum quantity required."
))

# ── md-db-predictions: update prediction 4 basis ─────────────────────────────
_, c = find('md-db-predictions')
src = ''.join(c['source'])
set_source(c, src.replace(
    'Buffer capacity near pKa reduces excursions',
    'Phosphate buffering (H₂PO₄⁻/HPO₄²⁻) reduces controller excursions'
))

# ── md-species-note: mention corrector species ───────────────────────────────
_, c = find('md-species-note')
set_source(c, (
    "### Species, kinetics, and reactions\n"
    "\n"
    "PEKILO biomass assumed CH₁.₆₁O₀.₅₆ (MW = 24.626 — generic fungal formula).\n"
    "NH₄Cl (5 g/L) and KH₂PO₄ (30 g/L) are each defined as an undissociated\n"
    "Species; their complete dissociations (NH₄Cl ⇌ NH₄⁺ + Cl⁻, log K = 3;\n"
    "KH₂PO₄ ⇌ K⁺ + H₂PO₄⁻, log K = 3) are specified as EquilibriumReactions\n"
    "so the initial composition is expressed in terms of the salt rather than\n"
    "its constituent ions. The full phosphate ladder\n"
    "(H₃PO₄ / H₂PO₄⁻ / HPO₄²⁻ / PO₄³⁻, pKₐ 2.15 / 7.20 / 12.35) and the\n"
    "NH₄⁺ ⇌ NH₃ + H⁺ equilibrium (pKₐ = 9.25) are included downstream.\n"
    "NaOH and H₂SO₄ are defined as Species with their dissociation equilibria\n"
    "(NaOH ⇌ Na⁺ + OH⁻; H₂SO₄ ⇌ HSO₄⁻ + H⁺, then HSO₄⁻ ⇌ SO₄²⁻ + H⁺)\n"
    "so that `equilibrate_to_pH` can use either to reach the pH setpoint."
))

# ── code-species: add NAOH, NA_PLUS, H2SO4, HSO4_MINUS, SO4_2MINUS ───────────
_, c = find('code-species')
set_source(c, (
    "ACETIC_ACID   = Species(id=\"AceticAcid\",  atoms={\"C\":2,\"H\":4,\"O\":2},       charge=0)\n"
    "ACETATE_MINUS = Species(id=\"Acetate-\",    atoms={\"C\":2,\"H\":3,\"O\":2},       charge=-1)\n"
    "PEKILO        = Species(id=\"PEKILO\",      atoms={\"C\":1,\"H\":1.61,\"O\":0.56}, charge=0, MW=24.626)  # biomass formula assumed generic fungal — replace when PEKILO-specific data are available\n"
    "NH4CL         = Species(id=\"NH4Cl\",       atoms={\"N\":1,\"H\":4,\"Cl\":1},      charge=0)\n"
    "NH4_PLUS      = Species(id=\"NH4+\",        atoms={\"N\":1,\"H\":4},              charge=+1)\n"
    "CHLORIDE      = Species(id=\"Cl-\",         atoms={\"Cl\":1},                    charge=-1)\n"
    "KH2PO4        = Species(id=\"KH2PO4\",      atoms={\"K\":1,\"H\":2,\"P\":1,\"O\":4}, charge=0)\n"
    "K_PLUS        = Species(id=\"K+\",          atoms={\"K\":1},                     charge=+1)\n"
    "H2PO4_MINUS   = Species(id=\"H2PO4-\",     atoms={\"H\":2,\"P\":1,\"O\":4},        charge=-1)\n"
    "HPO4_2MINUS   = Species(id=\"HPO4--\",     atoms={\"H\":1,\"P\":1,\"O\":4},        charge=-2)\n"
    "PO4_3MINUS    = Species(id=\"PO4---\",     atoms={\"P\":1,\"O\":4},               charge=-3)\n"
    "H3PO4         = Species(id=\"H3PO4\",       atoms={\"H\":3,\"P\":1,\"O\":4},        charge=0)\n"
    "# pH corrector species\n"
    "NAOH          = Species(id=\"NaOH\",        atoms={\"Na\":1,\"O\":1,\"H\":1},       charge=0)\n"
    "NA_PLUS       = Species(id=\"Na+\",         atoms={\"Na\":1},                    charge=+1)\n"
    "H2SO4         = Species(id=\"H2SO4\",       atoms={\"H\":2,\"S\":1,\"O\":4},        charge=0)\n"
    "HSO4_MINUS    = Species(id=\"HSO4-\",       atoms={\"H\":1,\"S\":1,\"O\":4},        charge=-1)\n"
    "SO4_2MINUS    = Species(id=\"SO4--\",       atoms={\"S\":1,\"O\":4},               charge=-2)\n"
    "\n"
    "rxn_growth = ReactionBuilder.monod_aerobic_growth(\n"
    "    substrate=ACETIC_ACID, biomass=PEKILO,\n"
    "    mu_max_per_h=MU_MAX, Ks_gL=KS_G_L, yield_gX_gS=YIELD,\n"
    "    label=\"growth_on_AceticAcid\",\n"
    ")\n"
    "\n"
    "# Read off nu_O2 from the actual stoichiometry (used in validation)\n"
    "nu_O2_act = abs(next(e.coefficient for e in rxn_growth.stoichiometry\n"
    "                     if e.species.id == \"O2\"))\n"
    "print(f\"\\u03bd_O2 from stoichiometry: {nu_O2_act:.4f} mol O2 / mol substrate\")"
))

# ── code-reactions: add corrector equilibria, update ReactionSystem ───────────
_, c = find('code-reactions')
set_source(c, (
    "rxn_water = EquilibriumReaction(\n"
    "    \"H2O,aq <-> H+,aq + OH-,aq\",\n"
    "    log_K=-14.0, dH_J_per_mol=55900.0, T_ref_K=298.15, label=\"eq_water\",\n"
    ")\n"
    "rxn_co2_aq = EquilibriumReaction(\n"
    "    \"CO2,aq + H2O,aq <-> HCO3-,aq + H+,aq\",\n"
    "    log_K=-6.35, dH_J_per_mol=7646.0, T_ref_K=298.15, total_id=\"CO2\", label=\"eq_CO2\",\n"
    ")\n"
    "rxn_co3 = EquilibriumReaction(\n"
    "    \"HCO3-,aq <-> H+,aq + CO3--,aq\", log_K=-10.33, label=\"eq_HCO3\",\n"
    ")\n"
    "rxn_dissoc = EquilibriumReaction(\n"
    "    \"AceticAcid,aq <-> Acetate-,aq + H+,aq\",\n"
    "    species={\"AceticAcid\": ACETIC_ACID, \"Acetate-\": ACETATE_MINUS},\n"
    "    log_K=-4.756, label=\"eq_AceticAcid\",\n"
    ")\n"
    "# Salt dissociations — log K = 3 represents essentially complete dissolution\n"
    "rxn_nh4cl = EquilibriumReaction(\n"
    "    \"NH4Cl,aq <-> NH4+,aq + Cl-,aq\",\n"
    "    species={\"NH4Cl\": NH4CL, \"NH4+\": NH4_PLUS, \"Cl-\": CHLORIDE},\n"
    "    log_K=3.0, label=\"eq_NH4Cl\",\n"
    ")\n"
    "rxn_kh2po4 = EquilibriumReaction(\n"
    "    \"KH2PO4,aq <-> K+,aq + H2PO4-,aq\",\n"
    "    species={\"KH2PO4\": KH2PO4, \"K+\": K_PLUS, \"H2PO4-\": H2PO4_MINUS},\n"
    "    log_K=3.0, label=\"eq_KH2PO4\",\n"
    ")\n"
    "# pH corrector dissociations\n"
    "rxn_naoh = EquilibriumReaction(\n"
    "    \"NaOH,aq <-> Na+,aq + OH-,aq\",\n"
    "    species={\"NaOH\": NAOH, \"Na+\": NA_PLUS},\n"
    "    log_K=3.0, label=\"eq_NaOH\",\n"
    ")\n"
    "# H2SO4 dissociation ladder — first step essentially complete, second pKa 1.99\n"
    "rxn_h2so4 = EquilibriumReaction(\n"
    "    \"H2SO4,aq <-> HSO4-,aq + H+,aq\",\n"
    "    species={\"H2SO4\": H2SO4, \"HSO4-\": HSO4_MINUS},\n"
    "    log_K=3.0, label=\"eq_H2SO4\",\n"
    ")\n"
    "rxn_hso4 = EquilibriumReaction(\n"
    "    \"HSO4-,aq <-> SO4--,aq + H+,aq\",\n"
    "    species={\"HSO4-\": HSO4_MINUS, \"SO4--\": SO4_2MINUS},\n"
    "    log_K=-1.99, dH_J_per_mol=-22600.0, T_ref_K=298.15, label=\"eq_HSO4\",\n"
    ")\n"
    "rxn_nh4 = EquilibriumReaction(\n"
    "    \"NH4+,aq <-> NH3,aq + H+,aq\",\n"
    "    species={\"NH4+\": NH4_PLUS},\n"
    "    log_K=-9.25, dH_J_per_mol=52215.0, T_ref_K=298.15, total_id=\"NH3\", label=\"eq_NH4\",\n"
    ")\n"
    "# Phosphate speciation ladder — pKa values at 25 °C (Martell & Smith, 1977)\n"
    "rxn_h3po4 = EquilibriumReaction(\n"
    "    \"H3PO4,aq <-> H2PO4-,aq + H+,aq\",\n"
    "    species={\"H3PO4\": H3PO4, \"H2PO4-\": H2PO4_MINUS},\n"
    "    log_K=-2.148, dH_J_per_mol=-8000.0, T_ref_K=298.15, label=\"eq_H3PO4\",\n"
    ")\n"
    "rxn_h2po4 = EquilibriumReaction(\n"
    "    \"H2PO4-,aq <-> HPO4--,aq + H+,aq\",\n"
    "    species={\"H2PO4-\": H2PO4_MINUS, \"HPO4--\": HPO4_2MINUS},\n"
    "    log_K=-7.198, dH_J_per_mol=3600.0, T_ref_K=298.15, label=\"eq_H2PO4\",\n"
    ")\n"
    "rxn_hpo4 = EquilibriumReaction(\n"
    "    \"HPO4--,aq <-> PO4---,aq + H+,aq\",\n"
    "    species={\"HPO4--\": HPO4_2MINUS, \"PO4---\": PO4_3MINUS},\n"
    "    log_K=-12.35, dH_J_per_mol=14600.0, T_ref_K=298.15, label=\"eq_HPO4\",\n"
    ")\n"
    "\n"
    "rxn_system = ReactionSystem(\n"
    "    [rxn_growth, rxn_water, rxn_co2_aq, rxn_dissoc, rxn_co3,\n"
    "     rxn_nh4cl, rxn_kh2po4, rxn_naoh, rxn_h2so4, rxn_hso4,\n"
    "     rxn_nh4, rxn_h3po4, rxn_h2po4, rxn_hpo4],\n"
    "    label=\"batch_fermenter_chemistry\",\n"
    ")\n"
    "print(f\"ReactionSystem: {len(rxn_system.kinetic_reactions)} kinetic, \"\n"
    "      f\"{len(rxn_system.single_phase_equilibria)} equilibria\")"
))

# ── code-phases: add corrector species at 0.0 ────────────────────────────────
_, c = find('code-phases')
set_source(c, (
    "n_total_gas = V_GAS / (R_L_ATM_MOL_K * T_K)\n"
    "\n"
    "gas = GasPhase(\n"
    "    n_mol={\"O2\": n_total_gas*0.2095, \"CO2\": n_total_gas*0.0004,\n"
    "           \"N2\": n_total_gas*0.7901, \"NH3\": 0.0},\n"
    "    V_L=V_GAS, T_K=T_K,\n"
    ")\n"
    "\n"
    "def _henry_n(sp):\n"
    "    pm = AD_BASIC.partition_models[sp]\n"
    "    return pm.H_ref * gas.p_atm.get(sp, 0.0) * 101325 / 1000 * V_LIQ\n"
    "\n"
    "liquid = LiquidPhase(\n"
    "    n_mol={\n"
    "        # Organic substrate and biomass\n"
    "        ACETIC_ACID.id:   (S0_G_L    / float(ACETIC_ACID.MW)) * V_LIQ,\n"
    "        ACETATE_MINUS.id: 0.0,\n"
    "        PEKILO.id:        (X0_G_L    / float(PEKILO.MW))      * V_LIQ,\n"
    "        # Salts as undissociated species; ions start at 0 and are\n"
    "        # populated by the dissociation equilibria on the first advance\n"
    "        NH4CL.id:         (NH4CL_G_L  / float(NH4CL.MW))     * V_LIQ,\n"
    "        NH4_PLUS.id:      0.0,\n"
    "        CHLORIDE.id:      0.0,\n"
    "        KH2PO4.id:        (KH2PO4_G_L / float(KH2PO4.MW))    * V_LIQ,\n"
    "        K_PLUS.id:        0.0,\n"
    "        H2PO4_MINUS.id:   0.0,\n"
    "        HPO4_2MINUS.id:   0.0,\n"
    "        H3PO4.id:         0.0,\n"
    "        PO4_3MINUS.id:    0.0,\n"
    "        # pH corrector species — start at zero, dosed by equilibrate_to_pH\n"
    "        NAOH.id:          0.0,\n"
    "        NA_PLUS.id:       0.0,\n"
    "        H2SO4.id:         0.0,\n"
    "        HSO4_MINUS.id:    0.0,\n"
    "        SO4_2MINUS.id:    0.0,\n"
    "        # Dissolved gases at Henry-law equilibrium with initial headspace\n"
    "        \"O2\":    _henry_n(\"O2\"),   \"CO2\":   _henry_n(\"CO2\"),\n"
    "        \"N2\":    _henry_n(\"N2\"),   \"NH3\":   _henry_n(\"NH3\"),\n"
    "        # Carbonate and proton species — derived on first advance\n"
    "        \"HCO3-\": 0.0,              \"CO3--\": 0.0,\n"
    "        \"OH-\":   0.0,              \"H+\":    1e-7 * V_LIQ,\n"
    "    },\n"
    "    V_L=V_LIQ, T_K=T_K,\n"
    ")\n"
    "\n"
    "transfer_models = {\n"
    "    \"O2\":  KineticTransferModel(AD_BASIC.partition_models[\"O2\"],\n"
    "               k_transfer=KLA_PER_H[\"O2\"]),\n"
    "    \"CO2\": KineticTransferModel(AD_BASIC.partition_models[\"CO2\"],\n"
    "               k_transfer=KLA_PER_H[\"CO2\"], transfer_basis=\"molecular\"),\n"
    "    \"NH3\": KineticTransferModel(AD_BASIC.partition_models[\"NH3\"],\n"
    "               k_transfer=KLA_PER_H[\"NH3\"], transfer_basis=\"molecular\"),\n"
    "    \"N2\":  EquilibriumTransferModel(AD_BASIC.partition_models[\"N2\"]),\n"
    "}\n"
    "\n"
    "cv = ControlVolume(\n"
    "    phases={\"gas\": gas, \"liquid\": liquid},\n"
    "    transfer_models=transfer_models,\n"
    "    reaction_system=rxn_system,\n"
    "    label=\"batch_fermenter\",\n"
    ")\n"
    "print(\"CV assembled\")"
))

# ── code-ph-correction: OH- → NaOH ───────────────────────────────────────────
_, c = find('code-ph-correction')
set_source(c, "n_corr = cv.equilibrate_to_pH(\"NaOH\", PH_SETPOINT)")

with open(path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"Done. {len(nb['cells'])} cells.")
