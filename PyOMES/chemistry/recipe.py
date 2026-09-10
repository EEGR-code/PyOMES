"""Recipe-style liquid definition utilities.

This module provides a user-friendly way to define an aqueous solution using
"what you weighed out" (e.g., grams of Na2SO4, mL of 1 M HCl) and convert it to
`AqueousTotalsUser` for use with the standalone speciation interface.

Design goals
------------
* No dependency on BioSTEAM/thermosteam for molecular weights or stoichiometry.
* Reuse the unified chemistry registry (ion normalization and salt dissociation).
* Keep the feature set intentionally small and explicit; raise clear errors for
  unknown compounds.

Scope
-----
This is *not* intended to be a full electrolyte thermodynamics database. It is
an input convenience layer that converts recipes into totals (strong ions, TIC,
phosphate, TAN, and optional alkalinity bookkeeping).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict

from .types import AqueousTotalsUser
from .registry import normalize_ion_label, resolve_compound


# --- Tiny units helpers (pure scalars) ---
g: float = 1.0
kg: float = 1000.0
mg: float = 1e-3
L: float = 1.0
mL: float = 1e-3


def _as_float(x) -> float:
    try:
        return float(x)
    except Exception as e:
        raise TypeError(f"Expected a number, got {type(x)!r}") from e


@dataclass
class SolutionRecipe:
    """A simple recipe builder for aqueous totals."""

    volume_L: float = 1.0
    T_C: float = 25.0

    # Accumulated totals in moles (not per-L until conversion)
    _strong_ions_mol: Dict[str, float] = field(default_factory=dict)
    _TIC_mol: float = 0.0
    _P_tot_mol: float = 0.0
    _TAN_mol: float = 0.0
    _TA_eq: float = 0.0
    _adds: list = field(default_factory=list)

    def _add_ion_mol(self, ion_label: str, mol: float) -> None:
        ion = normalize_ion_label(ion_label)
        self._strong_ions_mol[ion] = self._strong_ions_mol.get(ion, 0.0) + float(mol)

    def add_compound_g(self, name: str, mass_g: float) -> "SolutionRecipe":
        """Add a compound by mass (g)."""
        mass_g = _as_float(mass_g)
        rec = resolve_compound(name)
        MW = float(rec["MW_g_mol"])
        if MW <= 0:
            raise ValueError(f"Invalid MW for {name!r}: {MW}")
        n_mol = mass_g / MW
        return self.add_compound_mol(name, n_mol)

    def add_compound_mol(self, name: str, mol: float) -> "SolutionRecipe":
        """Add a compound by amount (mol)."""
        mol = _as_float(mol)
        rec = resolve_compound(name)

        for ion, stoich in (rec.get("ions", {}) or {}).items():
            self._add_ion_mol(ion, mol * float(stoich))

        self._TIC_mol += mol * float(rec.get("TIC_mol_per_mol", 0.0) or 0.0)
        self._P_tot_mol += mol * float(rec.get("P_tot_mol_per_mol", 0.0) or 0.0)
        self._TAN_mol += mol * float(rec.get("TAN_mol_per_mol", 0.0) or 0.0)
        self._TA_eq += mol * float(rec.get("alk_eq_per_mol", 0.0) or 0.0)

        self._adds.append(("compound", name, mol))
        return self

    def add_stock(self, name: str, molarity_mol_L: float, volume_mL: float) -> "SolutionRecipe":
        """Add a stock solution dose (molarity in mol/L, volume in mL)."""
        molarity_mol_L = _as_float(molarity_mol_L)
        volume_mL = _as_float(volume_mL)
        mol = molarity_mol_L * (volume_mL * mL)
        self._adds.append(("stock", name, molarity_mol_L, volume_mL))
        return self.add_compound_mol(name, mol)

    def to_totals_user(self) -> AqueousTotalsUser:
        """Convert accumulated recipe entries to `AqueousTotalsUser` (mol/L)."""
        V = float(self.volume_L)
        if V <= 0:
            raise ValueError("volume_L must be > 0")
        strong_ions_mol_L = {k: v / V for k, v in self._strong_ions_mol.items() if abs(v) > 0}
        return AqueousTotalsUser(
            T_K=float(self.T_C) + 273.15,
            TIC_mol_L=float(self._TIC_mol) / V,
            P_tot_mol_L=float(self._P_tot_mol) / V,
            TAN_mol_L=float(self._TAN_mol) / V,
            TA_eq_L=float(self._TA_eq) / V,
            strong_ions_mol_L=strong_ions_mol_L,
        )

    def to_totals(self) -> AqueousTotalsUser:
        return self.to_totals_user()

    def summary(self) -> str:
        lines = [f"SolutionRecipe(V={self.volume_L:g} L, T={self.T_C:g} °C)"]
        for entry in self._adds:
            lines.append(f"  - {entry}")
        totals = self.to_totals_user()
        lines.append("Totals (mol/L):")
        lines.append(f"  TIC: {totals.TIC_mol_L:.6g}")
        lines.append(f"  Ptot: {totals.P_tot_mol_L:.6g}")
        lines.append(f"  TAN: {totals.TAN_mol_L:.6g}")
        lines.append(f"  TA:  {totals.TA_eq_L:.6g} eq/L")
        if totals.strong_ions_mol_L:
            lines.append("  Strong ions:")
            for k, v in sorted(totals.strong_ions_mol_L.items()):
                lines.append(f"    {k}: {v:.6g}")
        return "\n".join(lines)
