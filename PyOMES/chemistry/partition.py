# -*- coding: utf-8 -*-
"""Phase-partition protocols and the ideal-gas multispecies VLE model.

A PartitionModel describes how a species distributes between two phases at
equilibrium.  ``PartitionModel`` is the per-species protocol and
``MultispeciesPartitionModel`` the protocol for models that solve all
species at once; ``MultispeciesVLEPartition`` implements the latter for an
ideal gas.  The single-species implementations (Henry, Raoult and Ksp) live
in :mod:`PyOMES.reactions.equilibrium.interphase`, since each also acts as an
equilibrium constraint.  Future non-linear models (Langmuir, Freundlich)
implement the same protocol.

Capacity parameters are phase-agnostic floats so the protocol is usable for
gas-liquid (volumes in L), solid-liquid (solid mass in kg), etc.

The partition_ratio() return value selects the solver path in KineticGasLiquidLink:
  - float  → analytical exponential step solution (exact for linear models)
  - None   → ODE instantaneous-rate path (required for non-linear models)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Optional, Protocol, runtime_checkable

from ..units import R_L_ATM_PER_MOL_K


@runtime_checkable
class PartitionModel(Protocol):
    """Protocol for phase-partition relationships."""

    def equilibrium_a_moles(
        self,
        n_total: float,
        capacity_a: float,
        capacity_b: float,
        T_K: float,
        alpha: float = 1.0,
    ) -> float:
        """Moles in phase 'a' at equilibrium given n_total across both phases."""
        ...

    def partition_ratio(
        self,
        capacity_a: float,
        capacity_b: float,
        T_K: float,
        alpha: float = 1.0,
    ) -> Optional[float]:
        """Dimensionless partition ratio for the analytical step solution.

        Returns the ratio ``n_a_eq / n_b_eq`` at equilibrium for linear
        models.  Return ``None`` for non-linear models (Langmuir, etc.)
        that require the ODE instantaneous-rate path instead.
        """
        ...


# ── Multispecies partition model protocol and VLE implementation ────────────

@runtime_checkable
class MultispeciesPartitionModel(Protocol):
    """Phase-partition protocol for models where equilibrium of one species
    depends on the concentrations of all other species.

    Use cases:
    - Coupled VLE with non-ideal gas EOS (Peng-Robinson mixing rules)
    - Competitive adsorption (Langmuir multicomponent)
    - NRTL liquid–liquid partitioning

    Unlike ``PartitionModel`` (per-species, independent), this protocol
    solves all species simultaneously in a single call.
    """

    def equilibrium_all_a_moles(
        self,
        n_total_all: Dict[str, float],
        capacity_a: float,
        capacity_b: float,
        T_K: float,
    ) -> Dict[str, float]:
        """Moles of each species in phase 'a' at simultaneous equilibrium.

        Parameters
        ----------
        n_total_all : dict
            ``{species_id: n_total_i}`` — total moles of each species
            summed across both phases.
        capacity_a : float
            Phase 'a' capacity (volume in L for gas-liquid; mass in kg
            for solid-liquid).
        capacity_b : float
            Phase 'b' capacity.
        T_K : float
            Temperature (K).

        Returns
        -------
        dict
            ``{species_id: n_a_i}`` — moles in phase 'a' at equilibrium.
            Species absent from the return dict are treated as n_a = 0.
        """
        ...


def _kH_mol_L_atm_from_ref(H_ref: float, dlnH: float, T_K: float, T_ref: float) -> float:
    """Convert Sander kH (mol m⁻³ Pa⁻¹) to mol/(L·atm) and apply van 't Hoff."""
    H_mol_L_atm = (H_ref / 1000.0) * 101325.0
    return H_mol_L_atm * math.exp(dlnH * (1.0 / T_K - 1.0 / T_ref))


@dataclass(frozen=True)
class MultispeciesVLEPartition:
    """Coupled gas-liquid VLE for multiple volatile species, ideal gas only.

    Each species decouples and the solution is analytical (same formula
    as ``HenryEquilibrium`` without the ``alpha`` correction). There is
    no EOS parameter — non-ideal (Peng-Robinson) coupling is not
    implemented.

    Satisfies ``MultispeciesPartitionModel``:
        ``equilibrium_all_a_moles(n_total_all, V_liq, V_gas, T_K)``

    Parameters
    ----------
    kH_ref : dict
        Henry solubility at T_ref per species (mol m⁻³ Pa⁻¹, Sander
        convention).  Only species with kH_ref entries are solved.
    dlnH : dict
        d(ln kH)/d(1/T) per species (K).  Missing entries default to 0.
    T_ref : float
        Reference temperature for kH_ref (K).  Default 298.15 K.

    Examples
    --------
    >>> from PyOMES.chemistry import MultispeciesVLEPartition
    >>> # H2S and CO2 in a 1 L liquid / 0.1 L headspace system
    >>> H_H2S = 0.10 * 1000 / 101325       # 0.10 mol/L/atm → Sander units
    >>> H_CO2 = 3.4e-4 * 1000 / 101325     # 3.4e-4 mol/L/atm → Sander units
    >>> vle = MultispeciesVLEPartition(
    ...     kH_ref={"H2S": H_H2S, "CO2": H_CO2},
    ...     dlnH={"H2S": 2100.0, "CO2": 2400.0},
    ... )
    >>> n_liq = vle.equilibrium_all_a_moles(
    ...     {"H2S": 0.01, "CO2": 0.02}, V_liq=1.0, V_gas=0.1, T_K=308.15
    ... )
    """

    kH_ref: Dict[str, float]
    dlnH: Dict[str, float] = field(default_factory=dict)
    T_ref: float = 298.15

    def equilibrium_all_a_moles(
        self,
        n_total_all: Dict[str, float],
        capacity_a: float,
        capacity_b: float,
        T_K: float,
    ) -> Dict[str, float]:
        """Solve VLE for all species using the ideal gas law (analytical, decoupled).

        For each species i with Henry constant kH_i(T):
            r_i = kH_i × R × T × V_liq / V_gas
            n_liq_i = r_i × n_total_i / (1 + r_i)

        Species not in ``kH_ref`` are ignored (returned as absent).
        """
        result = {}
        for species_id, n_total in n_total_all.items():
            if species_id not in self.kH_ref:
                continue
            H_ref = self.kH_ref[species_id]
            dlnH  = self.dlnH.get(species_id, 0.0)
            kH = _kH_mol_L_atm_from_ref(H_ref, dlnH, T_K, self.T_ref)
            r = kH * R_L_ATM_PER_MOL_K * T_K * capacity_a / capacity_b
            result[species_id] = r * float(n_total) / (1.0 + r)
        return result
