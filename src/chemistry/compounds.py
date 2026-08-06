# -*- coding: utf-8 -*-
"""Standalone chemical compound database.

This module provides molecular weights (g/mol) and atom compositions for all
compounds used by the fermenter simulation, **without** depending on BioSTEAM
or thermosteam.

Usage
-----
>>> from PyOMES.chemistry.compounds import ChemicalRegistry, Chemical
>>> reg = ChemicalRegistry.default()
>>> reg["AceticAcid"].MW
60.052
>>> reg["AceticAcid"].atoms
{'C': 2, 'H': 4, 'O': 2}

Extending
---------
To add a custom organism or compound, use :meth:`ChemicalRegistry.register`::

    reg.register(Chemical("MyOrganism", MW=24.0, atoms={"C": 1, "H": 1.8, "O": 0.5}))
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass(frozen=True)
class Chemical:
    """Lightweight standalone chemical definition."""
    ID: str
    MW: float                          # g/mol
    atoms: Dict[str, float] = field(default_factory=dict)
    phase: str = "l"                   # l, g, s  (informational only)
    CAS: Optional[str] = None

    def __post_init__(self):
        if self.MW <= 0:
            raise ValueError(f"MW must be > 0 for '{self.ID}', got {self.MW}")


class ChemicalRegistry:
    """Registry of :class:`Chemical` objects, keyed by ID.

    Acts like a dict but also exposes a ``.IDs`` property for compatibility
    with code that used ``self.chemicals.IDs`` from bioSTEAM.
    """

    def __init__(self, chemicals: Optional[Dict[str, Chemical]] = None):
        self._chems: Dict[str, Chemical] = dict(chemicals or {})

    # --- dict-like API ---
    def __getitem__(self, key: str) -> Chemical:
        try:
            return self._chems[key]
        except KeyError:
            raise KeyError(f"Chemical '{key}' not found in registry. "
                           f"Available: {sorted(self._chems.keys())}")

    def __contains__(self, key: str) -> bool:
        return key in self._chems

    def __iter__(self):
        return iter(self._chems)

    def __len__(self):
        return len(self._chems)

    @property
    def IDs(self):
        """Return tuple of registered chemical IDs (bioSTEAM compatibility)."""
        return tuple(self._chems.keys())

    def register(self, chem: Chemical) -> None:
        """Add or replace a chemical in the registry."""
        self._chems[chem.ID] = chem

    def get(self, key: str, default=None):
        return self._chems.get(key, default)

    # --- Factory ---
    @classmethod
    def default(cls) -> "ChemicalRegistry":
        """Return a registry pre-populated with common fermentation compounds."""
        reg = cls()
        for c in _DEFAULT_CHEMICALS:
            reg.register(c)
        return reg


# ============================================================================
# Default compound data
# ============================================================================
# Sources: PubChem, NIST, standard reference handbooks.
# Atom compositions are exact stoichiometric counts (float for pseudo-species).

_DEFAULT_CHEMICALS = [
    # --- Solvents ---
    Chemical("Water",        MW=18.015, atoms={"H": 2, "O": 1}, phase="l", CAS="7732-18-5"),
    Chemical("H2O",          MW=18.015, atoms={"H": 2, "O": 1}, phase="l", CAS="7732-18-5"),

    # --- Organic acids (substrates) ---
    Chemical("AceticAcid",   MW=60.052,  atoms={"C": 2, "H": 4, "O": 2},  phase="l", CAS="64-19-7"),
    Chemical("PropionicAcid",MW=74.079,  atoms={"C": 3, "H": 6, "O": 2},  phase="l", CAS="79-09-4"),
    Chemical("ButyricAcid",  MW=88.106,  atoms={"C": 4, "H": 8, "O": 2},  phase="l", CAS="107-92-6"),
    Chemical("CitricAcid",   MW=192.124, atoms={"C": 6, "H": 8, "O": 7},  phase="l", CAS="77-92-9"),

    # --- Gases ---
    Chemical("O2",  MW=31.998, atoms={"O": 2},         phase="g", CAS="7782-44-7"),
    Chemical("CO2", MW=44.009, atoms={"C": 1, "O": 2}, phase="g", CAS="124-38-9"),
    Chemical("N2",  MW=28.014, atoms={"N": 2},         phase="g", CAS="7727-37-9"),
    Chemical("NH3", MW=17.031, atoms={"N": 1, "H": 3}, phase="l", CAS="7664-41-7"),

    # --- Inorganic acids / bases ---
    Chemical("H3PO4", MW=97.994,  atoms={"H": 3, "P": 1, "O": 4},              phase="l", CAS="7664-38-2"),
    Chemical("KOH",   MW=56.106,  atoms={"K": 1, "O": 1, "H": 1},              phase="s", CAS="1310-58-3"),
    Chemical("NaOH",  MW=39.997,  atoms={"Na": 1, "O": 1, "H": 1},             phase="s", CAS="1310-73-2"),
    Chemical("NaHCO3",MW=84.007,  atoms={"Na": 1, "H": 1, "C": 1, "O": 3},     phase="s", CAS="144-55-8"),
    Chemical("HCl",   MW=36.461,  atoms={"H": 1, "Cl": 1},                      phase="l", CAS="7647-01-0"),
    Chemical("H2SO4", MW=98.079,  atoms={"H": 2, "S": 1, "O": 4},              phase="l", CAS="7664-93-9"),

    # --- Salts / nutrients ---
    Chemical("AmmoniumSulfate", MW=132.14, atoms={"N": 2, "H": 8, "S": 1, "O": 4}, phase="s", CAS="7783-20-2"),
    Chemical("(NH4)2SO4",       MW=132.14, atoms={"N": 2, "H": 8, "S": 1, "O": 4}, phase="s", CAS="7783-20-2"),
    Chemical("NH4Cl",           MW=53.491, atoms={"N": 1, "H": 4, "Cl": 1},         phase="s", CAS="12125-02-9"),
    Chemical("KH2PO4",  MW=136.086, atoms={"K": 1, "H": 2, "P": 1, "O": 4},  phase="s", CAS="7778-77-0"),
    Chemical("NaH2PO4", MW=119.977, atoms={"Na": 1, "H": 2, "P": 1, "O": 4}, phase="s", CAS="7558-80-7"),
    Chemical("Na2HPO4", MW=141.959, atoms={"Na": 2, "H": 1, "P": 1, "O": 4}, phase="s", CAS="7558-79-4"),
    Chemical("Na2SO4",  MW=142.04,  atoms={"Na": 2, "S": 1, "O": 4},          phase="s", CAS="7757-82-6"),
    Chemical("NaCl",    MW=58.44,   atoms={"Na": 1, "Cl": 1},                 phase="s", CAS="7647-14-5"),
    Chemical("KCl",     MW=74.55,   atoms={"K": 1, "Cl": 1},                  phase="s", CAS="7447-40-7"),
    Chemical("MgSO4",   MW=120.37,  atoms={"Mg": 1, "S": 1, "O": 4},         phase="s", CAS="10034-99-8"),
    Chemical("CaCl2",   MW=110.98,  atoms={"Ca": 1, "Cl": 2},                 phase="s", CAS="10043-52-4"),
    Chemical("CaSO4",   MW=136.14,  atoms={"Ca": 1, "S": 1, "O": 4},         phase="s", CAS="7778-18-9"),
    Chemical("ZnSO4",   MW=161.47,  atoms={"Zn": 1, "S": 1, "O": 4},         phase="s", CAS="7446-20-0"),
    Chemical("MnCl2",   MW=125.84,  atoms={"Mn": 1, "Cl": 2},                 phase="s", CAS="7773-01-5"),
    Chemical("CoCl2",   MW=129.84,  atoms={"Co": 1, "Cl": 2},                 phase="s", CAS="7646-79-9"),
    Chemical("AmmoniumMolybdate", MW=1163.95, atoms={"N": 6, "H": 24, "Mo": 7, "O": 24}, phase="s"),

    # --- Default biomass pseudo-species ---
    # CHO-only yeast: CH1.61O0.56
    Chemical("Yeast_CHO", MW=24.626, atoms={"C": 1, "H": 1.61, "O": 0.56}, phase="s"),
    # CHNO yeast: CH1.61O0.56N0.16
    Chemical("Yeast",     MW=26.868, atoms={"C": 1, "H": 1.61, "O": 0.56, "N": 0.16}, phase="s"),
]
