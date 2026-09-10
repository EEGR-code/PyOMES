# -*- coding: utf-8 -*-
"""Lightweight feed composition type for the fermenter framework.

:class:`FeedState` stores species concentrations and provides g/L accessors
via :class:`~PyOMES.chemistry.compounds.ChemicalRegistry`.  It is the
canonical input type for feed-composition data throughout PyOMES.

Quick start
-----------
>>> from PyOMES.stream_adapter import FeedState
>>> feed = FeedState.from_mass_concentrations(
...     {"AceticAcid": 1.0, "Yeast": 0.1},   # g/L
...     volume_L=1.0,
... )
>>> feed.g_L("AceticAcid")
1.0
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional

from PyOMES.chemistry.compounds import ChemicalRegistry


@dataclass
class FeedState:
    """Lightweight representation of a liquid feed (or broth) composition.

    All concentrations are stored internally as **mol/L**.  Mass-based
    concentrations (g/L) are derived on-the-fly from the chemical registry.

    Parameters
    ----------
    concentrations_mol_L : dict
        Species-ID → concentration in mol/L.
    volume_L : float
        Reference liquid volume (litres).  Used for inventory calculations.
    T_K : float
        Temperature in Kelvin.
    registry : ChemicalRegistry | None
        Chemical database used for MW lookups. If ``None``, the default
        registry is created automatically.
    """
    concentrations_mol_L: Dict[str, float] = field(default_factory=dict)
    volume_L: float = 1.0
    T_K: float = 305.15
    registry: Optional[ChemicalRegistry] = field(default=None, repr=False)

    def __post_init__(self):
        if self.registry is None:
            self.registry = ChemicalRegistry.default()

    # ------------------------------------------------------------------
    # Concentration accessors (mirror the helpers on the old unit class)
    # ------------------------------------------------------------------
    def mol_L(self, species_id: str) -> float:
        """Return molar concentration (mol/L) for *species_id*."""
        return float(self.concentrations_mol_L.get(species_id, 0.0))

    def g_L(self, species_id: str) -> float:
        """Return mass concentration (g/L) for *species_id*."""
        mol = self.mol_L(species_id)
        if mol <= 0.0:
            return 0.0
        chem = self.registry.get(species_id)
        if chem is None:
            return 0.0
        return mol * chem.MW

    def has(self, species_id: str) -> bool:
        """True if *species_id* is present with a non-zero concentration."""
        return self.concentrations_mol_L.get(species_id, 0.0) > 0.0

    @property
    def species_ids(self):
        """Return tuple of species IDs present in this feed."""
        return tuple(self.concentrations_mol_L.keys())

    # ------------------------------------------------------------------
    # Convenient constructors
    # ------------------------------------------------------------------
    @classmethod
    def from_mass_concentrations(
        cls,
        mass_g_L: Dict[str, float],
        volume_L: float = 1.0,
        T_K: float = 305.15,
        registry: Optional[ChemicalRegistry] = None,
    ) -> "FeedState":
        """Create a FeedState from mass concentrations (g/L).

        Parameters
        ----------
        mass_g_L : dict
            Species-ID → concentration in g/L.
        """
        if registry is None:
            registry = ChemicalRegistry.default()
        mol_L: Dict[str, float] = {}
        for sid, gL in mass_g_L.items():
            gL = float(gL)
            if gL <= 0.0:
                continue
            chem = registry.get(sid)
            if chem is None:
                raise KeyError(
                    f"Chemical '{sid}' not found in registry. "
                    f"Register it first or use concentrations_mol_L directly."
                )
            mol_L[sid] = gL / chem.MW
        return cls(
            concentrations_mol_L=mol_L,
            volume_L=volume_L,
            T_K=T_K,
            registry=registry,
        )

    @classmethod
    def from_mixed_concentrations(
        cls,
        mass_g_L: Optional[Dict[str, float]] = None,
        molar_mol_L: Optional[Dict[str, float]] = None,
        volume_L: float = 1.0,
        T_K: float = 305.15,
        registry: Optional[ChemicalRegistry] = None,
    ) -> "FeedState":
        """Create a FeedState from a mix of g/L and mol/L inputs.

        This is handy when some species are most naturally specified in g/L
        (e.g. biomass, substrates) and others in mol/L (e.g. salts, acids).
        """
        if registry is None:
            registry = ChemicalRegistry.default()
        mol_L: Dict[str, float] = {}

        # mol/L entries first
        if molar_mol_L:
            for sid, c in molar_mol_L.items():
                c = float(c)
                if c > 0.0:
                    mol_L[sid] = c

        # g/L entries (converted)
        if mass_g_L:
            for sid, gL in mass_g_L.items():
                gL = float(gL)
                if gL <= 0.0:
                    continue
                chem = registry.get(sid)
                if chem is None:
                    raise KeyError(f"Chemical '{sid}' not found in registry.")
                mol_L[sid] = mol_L.get(sid, 0.0) + gL / chem.MW

        return cls(
            concentrations_mol_L=mol_L,
            volume_L=volume_L,
            T_K=T_K,
            registry=registry,
        )

    # ------------------------------------------------------------------
    # BioSTEAM bridge (optional, only if bioSTEAM is installed)
    # ------------------------------------------------------------------
    @classmethod
    def from_biosteam_stream(
        cls,
        stream,
        volume_L: Optional[float] = None,
        registry: Optional[ChemicalRegistry] = None,
    ) -> "FeedState":
        """Convert a BioSTEAM/thermosteam ``Stream`` to a ``FeedState``.

        This is an optional convenience for users migrating from bioSTEAM.
        BioSTEAM is **not** required to use the fermenter; this method only
        works if bioSTEAM is already installed and the stream object is valid.

        Parameters
        ----------
        stream : biosteam.Stream
            A BioSTEAM stream with ``F_vol``, ``imol``, ``imass`` attributes.
        volume_L : float, optional
            Override liquid volume.  If None, derived from the stream.
        registry : ChemicalRegistry, optional
            Chemical database.  If None, the default registry is used.
        """
        if registry is None:
            registry = ChemicalRegistry.default()

        F_vol_m3_hr = float(getattr(stream, "F_vol", 0.0) or 0.0)
        if F_vol_m3_hr <= 0.0:
            return cls(registry=registry, T_K=float(getattr(stream, "T", 305.15)))

        # imol is kmol/hr, F_vol is m³/hr → kmol/m³ = mol/L
        mol_L: Dict[str, float] = {}
        for cid in stream.chemicals.IDs:
            try:
                c = float(stream.imol[cid]) / F_vol_m3_hr  # mol/L
                if c > 0.0:
                    mol_L[cid] = c
            except (TypeError, ValueError, AttributeError, KeyError):
                continue

        if volume_L is None:
            # Use 1 L as default reference volume (concentrations are already per-L)
            volume_L = 1.0

        return cls(
            concentrations_mol_L=mol_L,
            volume_L=volume_L,
            T_K=float(getattr(stream, "T", 305.15)),
            registry=registry,
        )
