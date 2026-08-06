"""Replace acetic acid substrate with glucose throughout Example2."""
import json

path = r'c:\Users\k2473520\PyOMES\demos\model_api\D2Cworkshop\updated_layout\Example2_batch_fermenter.ipynb'
with open(path, 'r', encoding='utf-8') as f:
    nb = json.load(f)


def set_source(cell, text):
    lines = text.split('\n')
    cell['source'] = [l + '\n' for l in lines[:-1]] + ([lines[-1]] if lines[-1] else [])


def find(cid):
    for c in nb['cells']:
        if c.get('id') == cid:
            return c
    raise KeyError(cid)


# ── md-background ─────────────────────────────────────────────────────────────
set_source(find('md-background'), (
    "## Background\n"
    "\n"
    "### From closed well to sparged fermenter\n"
    "\n"
    "The MTP well (Example 1) is O₂-limited by design: the sealed headspace is a\n"
    "finite reservoir and growth must complete within it. A sparged bench fermenter\n"
    "breaks this constraint by continuously supplying fresh air through a sparger —\n"
    "growth is no longer limited by the headspace inventory but by the rate at which\n"
    "O₂ can transfer from the rising bubbles into the liquid.\n"
    "\n"
    "This introduces three new design questions that Example 1 did not need to answer:\n"
    "\n"
    "1. **kLa selection** — the volumetric mass transfer coefficient is the primary\n"
    "   determinant of oxygen supply. Too low and DO crashes; too high and dissolved\n"
    "   CO₂ is stripped before the speciation engine can buffer it.\n"
    "\n"
    "2. **Sparging rate** — vvm (gas volume per liquid volume per minute) sets both\n"
    "   the molar O₂ supply and the headspace pressure build-up that the vent must\n"
    "   manage. They cannot be chosen independently.\n"
    "\n"
    "3. **pH control** — aerobic respiration on glucose produces CO₂, which dissolves\n"
    "   and acidifies the broth through the carbonate equilibria. A proportional\n"
    "   controller doses acid or base, but the dose rate must remain physically\n"
    "   feasible (concentration and pump speed limits).\n"
    "\n"
    "### Why batch first?\n"
    "\n"
    "Batch operation is the simplest continuous-gas scenario: the liquid inventory\n"
    "is constant, there are no feed or drain flows, and the only steady-state\n"
    "condition is the final endpoint (substrate exhaustion). Understanding batch\n"
    "behaviour first simplifies the analysis of continuous operation, where feed and\n"
    "bleed flows add a second design degree of freedom."
))

# ── md-db-config: update pH setpoint basis ────────────────────────────────────
c = find('md-db-config')
src = ''.join(c['source'])
set_source(c, src.replace(
    '| pH setpoint | 6.0 | Phosphate buffer range; acetate 94 % dissociated |',
    '| pH setpoint | 6.0 | Phosphate buffer range; within optimal fungal growth window |'
))

# ── md-db-otr: update nu_O2 value in prose ────────────────────────────────────
c = find('md-db-otr')
src = ''.join(c['source'])
set_source(c, src.replace(
    'stoichiometric O₂ coefficient $\\nu_{O_2} = 1.015$ from the elemental\nbalance:',
    'stoichiometric O₂ coefficient $\\nu_{O_2} = 3.044$ from the elemental\nbalance for glucose oxidation to PEKILO biomass:'
))

# ── code-db-otr: update MW_S and nu_O2 ───────────────────────────────────────
c = find('code-db-otr')
src = ''.join(c['source'])
src = src.replace(
    "MW_S   = 60.052   # g/mol  AceticAcid",
    "MW_S   = 180.156  # g/mol  Glucose"
)
src = src.replace(
    "nu_O2  = 1.015    # mol O2 / mol substrate (elemental balance)",
    "nu_O2  = 3.044    # mol O2 / mol substrate (elemental balance, glucose → PEKILO biomass)"
)
set_source(c, src)

# ── md-db-ph: rewrite for glucose / pH 6.0 ───────────────────────────────────
set_source(find('md-db-ph'), (
    "### pH setpoint selection\n"
    "\n"
    "We operate at pH 6.0 for two complementary reasons.\n"
    "\n"
    "**Fungal growth window.** PEKILO growth rate is maintained well across\n"
    "pH 5–7; pH 6.0 sits near the centre of this range and avoids the\n"
    "inhibitory extremes at the edges.\n"
    "\n"
    "**Phosphate buffer effectiveness.** KH₂PO₄ (30 g/L) provides the\n"
    "principal chemical buffering via the H₂PO₄⁻ / HPO₄²⁻ couple\n"
    "($\\mathrm{p}K_a = 7.2$). At pH 6.0 the buffer is 1.2 units below\n"
    "this $\\mathrm{p}K_a$, so virtually all dissolved phosphate is in the\n"
    "H₂PO₄⁻ form — this is where the buffer has its largest base-absorbing\n"
    "capacity and resists the downward pH drift from CO₂ production most\n"
    "effectively.\n"
    "\n"
    "CO₂ dissolved at this setpoint exists mainly as CO₂(aq) rather than\n"
    "HCO₃⁻: the CO₂ / HCO₃⁻ equilibrium ($\\mathrm{p}K_a = 6.35$) gives\n"
    "\n"
    "$$\n"
    "\\frac{[\\text{HCO}_3^-]}{[\\text{CO}_2(\\text{aq})]} = 10^{6.0 - 6.35} = 0.45\n"
    "$$\n"
    "\n"
    "so dissolved inorganic carbon is 69 % CO₂(aq) and 31 % HCO₃⁻ at the\n"
    "setpoint, which keeps the bicarbonate buffering demand small.\n"
    "\n"
    "The liquid is brought to this setpoint before the simulation begins by\n"
    "`ControlVolume.equilibrate_to_pH`, which adds NaOH (if base is needed) or\n"
    "H₂SO₄ (if acid is needed) in the minimum quantity required."
))

# ── md-impl-header: remove Henderson-Hasselbalch reference ───────────────────
c = find('md-impl-header')
src = ''.join(c['source'])
set_source(c, src.replace(
    "and pH 5.0 from the Henderson–Hasselbalch analysis.",
    "and pH 6.0 from the design analysis."
))

# ── md-params-note: same fix ──────────────────────────────────────────────────
c = find('md-params-note')
src = ''.join(c['source'])
set_source(c, src.replace(
    "and pH 5.0 from the Henderson–Hasselbalch analysis.",
    "and pH 6.0 from the design analysis."
))

# ── md-species-note ───────────────────────────────────────────────────────────
set_source(find('md-species-note'), (
    "### Species, kinetics, and reactions\n"
    "\n"
    "Glucose (C₆H₁₂O₆, MW = 180.156 g/mol) is the carbon and energy source;\n"
    "PEKILO biomass is modelled as CH₁.₆₁O₀.₅₆ (MW = 24.626 — generic fungal\n"
    "formula). Glucose is not an electrolyte, so no acid–base equilibria are\n"
    "added for the substrate itself.\n"
    "\n"
    "NH₄Cl (5 g/L) and KH₂PO₄ (30 g/L) are each defined as an undissociated\n"
    "Species; their complete dissociations (NH₄Cl ⇌ NH₄⁺ + Cl⁻, log K = 3;\n"
    "KH₂PO₄ ⇌ K⁺ + H₂PO₄⁻, log K = 3) are specified as EquilibriumReactions\n"
    "so the initial composition is expressed in terms of the salt rather than\n"
    "its constituent ions. The full phosphate ladder\n"
    "(H₃PO₄ / H₂PO₄⁻ / HPO₄²⁻ / PO₄³⁻, pKₐ 2.15 / 7.20 / 12.35) and the\n"
    "NH₄⁺ ⇌ NH₃ + H⁺ equilibrium (pKₐ = 9.25) are included downstream.\n"
    "NaOH and H₂SO₄ are defined as Species with their dissociation equilibria\n"
    "(NaOH ⇌ Na⁺ + OH⁻; H₂SO₄ ⇌ HSO₄⁻ + H⁺, then HSO₄⁻ ⇌ SO₄²⁻ + H⁺)\n"
    "so that `equilibrate_to_pH` can use either to reach the pH setpoint."
))

# ── code-species: swap substrate ──────────────────────────────────────────────
set_source(find('code-species'), (
    "GLUCOSE       = Species(id=\"Glucose\",    atoms={\"C\":6,\"H\":12,\"O\":6},         charge=0)\n"
    "PEKILO        = Species(id=\"PEKILO\",      atoms={\"C\":1,\"H\":1.61,\"O\":0.56},   charge=0, MW=24.626)  # biomass formula assumed generic fungal — replace when PEKILO-specific data are available\n"
    "NH4CL         = Species(id=\"NH4Cl\",       atoms={\"N\":1,\"H\":4,\"Cl\":1},        charge=0)\n"
    "NH4_PLUS      = Species(id=\"NH4+\",        atoms={\"N\":1,\"H\":4},                charge=+1)\n"
    "CHLORIDE      = Species(id=\"Cl-\",         atoms={\"Cl\":1},                      charge=-1)\n"
    "KH2PO4        = Species(id=\"KH2PO4\",      atoms={\"K\":1,\"H\":2,\"P\":1,\"O\":4},  charge=0)\n"
    "K_PLUS        = Species(id=\"K+\",          atoms={\"K\":1},                       charge=+1)\n"
    "H2PO4_MINUS   = Species(id=\"H2PO4-\",     atoms={\"H\":2,\"P\":1,\"O\":4},          charge=-1)\n"
    "HPO4_2MINUS   = Species(id=\"HPO4--\",     atoms={\"H\":1,\"P\":1,\"O\":4},          charge=-2)\n"
    "PO4_3MINUS    = Species(id=\"PO4---\",     atoms={\"P\":1,\"O\":4},                 charge=-3)\n"
    "H3PO4         = Species(id=\"H3PO4\",       atoms={\"H\":3,\"P\":1,\"O\":4},          charge=0)\n"
    "# pH corrector species\n"
    "NAOH          = Species(id=\"NaOH\",        atoms={\"Na\":1,\"O\":1,\"H\":1},         charge=0)\n"
    "NA_PLUS       = Species(id=\"Na+\",         atoms={\"Na\":1},                      charge=+1)\n"
    "H2SO4         = Species(id=\"H2SO4\",       atoms={\"H\":2,\"S\":1,\"O\":4},          charge=0)\n"
    "HSO4_MINUS    = Species(id=\"HSO4-\",       atoms={\"H\":1,\"S\":1,\"O\":4},          charge=-1)\n"
    "SO4_2MINUS    = Species(id=\"SO4--\",       atoms={\"S\":1,\"O\":4},                 charge=-2)\n"
    "\n"
    "rxn_growth = ReactionBuilder.monod_aerobic_growth(\n"
    "    substrate=GLUCOSE, biomass=PEKILO,\n"
    "    mu_max_per_h=MU_MAX, Ks_gL=KS_G_L, yield_gX_gS=YIELD,\n"
    "    label=\"growth_on_Glucose\",\n"
    ")\n"
    "\n"
    "# Read off nu_O2 from the actual stoichiometry (used in validation)\n"
    "nu_O2_act = abs(next(e.coefficient for e in rxn_growth.stoichiometry\n"
    "                     if e.species.id == \"O2\"))\n"
    "print(f\"Glucose MW: {float(GLUCOSE.MW):.3f} g/mol\")\n"
    "print(f\"ν_O2 from stoichiometry: {nu_O2_act:.4f} mol O2 / mol substrate\")"
))

# ── code-reactions: drop rxn_dissoc, update ReactionSystem ───────────────────
set_source(find('code-reactions'), (
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
    "    [rxn_growth, rxn_water, rxn_co2_aq, rxn_co3,\n"
    "     rxn_nh4cl, rxn_kh2po4, rxn_naoh, rxn_h2so4, rxn_hso4,\n"
    "     rxn_nh4, rxn_h3po4, rxn_h2po4, rxn_hpo4],\n"
    "    label=\"batch_fermenter_chemistry\",\n"
    ")\n"
    "print(f\"ReactionSystem: {len(rxn_system.kinetic_reactions)} kinetic, \"\n"
    "      f\"{len(rxn_system.single_phase_equilibria)} equilibria\")"
))

# ── md-phases-note: pH 5 → pH 6 ──────────────────────────────────────────────
c = find('md-phases-note')
src = ''.join(c['source'])
set_source(c, src.replace('trace at\n   pH 5 but', 'trace at\n   pH 6 but'))

# ── code-phases: swap substrate species ───────────────────────────────────────
set_source(find('code-phases'), (
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
    "        # Substrate and biomass\n"
    "        GLUCOSE.id:       (S0_G_L    / float(GLUCOSE.MW))    * V_LIQ,\n"
    "        PEKILO.id:        (X0_G_L    / float(PEKILO.MW))     * V_LIQ,\n"
    "        # Salts as undissociated species; ions start at 0 and are\n"
    "        # populated by the dissociation equilibria on the first advance\n"
    "        NH4CL.id:         (NH4CL_G_L  / float(NH4CL.MW))    * V_LIQ,\n"
    "        NH4_PLUS.id:      0.0,\n"
    "        CHLORIDE.id:      0.0,\n"
    "        KH2PO4.id:        (KH2PO4_G_L / float(KH2PO4.MW))   * V_LIQ,\n"
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

# ── code-validation-otr: MW_S and substrate key ───────────────────────────────
set_source(find('code-validation-otr'), (
    "liq     = result.liquid_mol[\"main\"]\n"
    "gas_mol = result.gas_mol[\"main\"]\n"
    "t       = result.t_h\n"
    "MW_S    = float(GLUCOSE.MW)\n"
    "MW_X    = float(PEKILO.MW)\n"
    "\n"
    "# Headspace O2 partial pressure at each timestep\n"
    "n_gas_tot = sum(gas_mol[sp] for sp in gas_mol)\n"
    "P_tot_t   = n_gas_tot * R_L_ATM_MOL_K * T_K / V_GAS  # atm\n"
    "y_O2_t    = gas_mol[\"O2\"] / n_gas_tot\n"
    "P_O2_t    = y_O2_t * P_tot_t\n"
    "\n"
    "# Henry equilibrium concentration (mol/L) at headspace partial pressure\n"
    "pm_O2    = AD_BASIC.partition_models[\"O2\"]\n"
    "C_star_t = pm_O2.H_ref * P_O2_t * 101325 / 1000\n"
    "\n"
    "# Dissolved O2 (mol/L)\n"
    "C_O2_t   = liq[\"O2\"] / V_LIQ\n"
    "\n"
    "# OTR (mol/h)\n"
    "OTR_t = KLA_PER_H[\"O2\"] * (C_star_t - C_O2_t) * V_LIQ\n"
    "\n"
    "# OUR from Monod kinetics (mol/h)\n"
    "S_gL = liq[\"Glucose\"] / V_LIQ * MW_S\n"
    "X_gL = liq[\"PEKILO\"]  / V_LIQ * MW_X\n"
    "mu_t = MU_MAX * S_gL / (KS_G_L + S_gL)\n"
    "OUR_t = nu_O2_act * (mu_t / YIELD * X_gL / MW_S * V_LIQ)\n"
    "\n"
    "# DO percent saturation\n"
    "DO_pct = C_O2_t / np.maximum(C_star_t, 1e-30) * 100\n"
    "\n"
    "print(\"=\" * 50)\n"
    "print(\"OTR / OUR validation\")\n"
    "print(\"=\" * 50)\n"
    "print(f\"\\nPeak OTR  : {OTR_t.max()*1000:.1f} mmol/h\")\n"
    "print(f\"Peak OUR  : {OUR_t.max()*1000:.1f} mmol/h\")\n"
    "print(f\"Min OTR/OUR ratio: {(OTR_t/np.maximum(OUR_t,1e-30)).min():.2f}\")\n"
    "status1 = \"PASS\" if (OTR_t >= OUR_t * 0.9).all() else \"FAIL\"\n"
    "print(f\"OTR \\u2265 OUR throughout: {status1}\")\n"
    "\n"
    "print(f\"\\nMinimum DO : {DO_pct.min():.1f}%  (design basis: > 0%)\")\n"
    "status2 = \"PASS\" if DO_pct.min() > 0.5 else \"FAIL\"\n"
    "print(f\"DO > 0 throughout   : {status2}\")"
))

# ── code-validation-yield: substrate key ─────────────────────────────────────
set_source(find('code-validation-yield'), (
    "C_S0 = liq[\"Glucose\"][0]  / V_LIQ * MW_S\n"
    "C_Sf = liq[\"Glucose\"][-1] / V_LIQ * MW_S\n"
    "C_X0 = liq[\"PEKILO\"][0]   / V_LIQ * MW_X\n"
    "C_Xf = liq[\"PEKILO\"][-1]  / V_LIQ * MW_X\n"
    "delta_S  = C_S0 - C_Sf\n"
    "delta_X  = C_Xf - C_X0\n"
    "Y_obs    = delta_X / delta_S if delta_S > 1e-6 else float(\"nan\")\n"
    "X_pred   = YIELD * S0_G_L + X0_G_L\n"
    "\n"
    "print(\"Prediction 3 — endpoint biomass:\")\n"
    "print(f\"  Predicted X_f : {X_pred:.3f} g/L\")\n"
    "print(f\"  Simulated X_f : {C_Xf:.3f} g/L\")\n"
    "err3 = abs(C_Xf - X_pred) / X_pred * 100\n"
    "print(f\"  Error         : {err3:.1f}%  [{('PASS' if err3 < 5 else 'CHECK')}]\")\n"
    "print(f\"\\n  Observed yield ΔX/ΔS = {Y_obs:.3f}  (parameter Y = {YIELD:.3f})\")\n"
    "err_y = abs(Y_obs - YIELD) / YIELD * 100\n"
    "print(f\"  Yield error          : {err_y:.1f}%  [{('PASS' if err_y < 2 else 'CHECK')}]\")\n"
    "\n"
    "pH     = result.pH[\"main\"]\n"
    "ph_val = pH[np.isfinite(pH)]\n"
    "max_dev = np.abs(ph_val - PH_SETPOINT).max()\n"
    "print(f\"\\nPrediction 4 — pH tracking:\")\n"
    "print(f\"  pH range  : {ph_val.min():.3f} – {ph_val.max():.3f}\")\n"
    "print(f\"  Max deviation from setpoint: {max_dev:.3f}\")\n"
    "print(f\"  [{('PASS' if max_dev < 0.5 else 'CHECK')}]\")"
))

# ── code-plots: substrate label ───────────────────────────────────────────────
c = find('code-plots')
src = ''.join(c['source'])
src = src.replace('C_S = liq["AceticAcid"] / V_LIQ * MW_S', 'C_S = liq["Glucose"] / V_LIQ * MW_S')
src = src.replace('label="Acetic acid (g/L)"', 'label="Glucose (g/L)"')
set_source(c, src)

# ── md-discussion: fix acetic acid reference ──────────────────────────────────
c = find('md-discussion')
src = ''.join(c['source'])
set_source(c, src.replace(
    "roughly 10×\nthe stoichiometric NaOH demand from the acetic acid dissociation shift. If the\n"
    "simulation shows the controller saturating at this limit, it would indicate the\n"
    "setpoint cannot be maintained at the chosen growth rate — a practical constraint\n"
    "that would require either a stronger base solution or a lower pH target.",
    "sufficient to comfortably exceed the NaOH demand from CO₂-driven acidification\n"
    "during peak growth. If the simulation shows the controller saturating at this\n"
    "limit, it would indicate the setpoint cannot be maintained at the chosen growth\n"
    "rate — a practical constraint that would require either a stronger base solution\n"
    "or a higher pH target."
))

with open(path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"Done. {len(nb['cells'])} cells.")
