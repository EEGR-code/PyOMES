{
  "meta": {
    "notes": [
      "Cases are recipe-based (materials added) and then equilibrated with a gas phase.",
      "Compatibility: mirrors keys used in run_test_NIST_buffer_standards_with_recipe.py (name, pH_ref, tol, recipe_mol_L, acid_pKas) plus gas_equilibration.",
      "gas_equilibration uses either pCO2_atm directly OR gas_y (mole fractions) + P_total_atm."
    ]
  },
  "CASES": [
    {
      "category": "gas_equilibrated",
      "name": "High-purity water equilibrated with ambient air (CO2 ~ 400 ppm) at 25C",
      "pH_ref": 5.60,
      "tol": 0.15,
      "T_C": 25.0,
      "recipe_mol_L": {},
      "acid_pKas": {},
      "gas_equilibration": {
        "mode": "open_fixed_pCO2",
        "P_total_atm": 1.0,
        "gas_y": { "CO2": 0.0004, "O2": 0.2095, "N2": 0.7901 },
        "kH_CO2_mol_L_atm": 0.033
      }
    },
    {
      "category": "gas_equilibrated",
      "name": "Physiology-style bicarbonate system: 24 mM NaHCO3 equilibrated with CO2 at 40 mmHg (~5.3%) and 37C",
      "pH_ref": 7.40,
      "tol": 0.15,
      "T_C": 37.0,
      "recipe_mol_L": { "NaHCO3": 0.024 },
      "acid_pKas": {},
      "gas_equilibration": {
        "mode": "open_fixed_pCO2",
        "P_total_atm": 1.0,
        "pCO2_atm": 0.0526,
        "kH_CO2_mol_L_atm": 0.030
      }
    }
  ]
}
