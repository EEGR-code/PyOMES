# -*- coding: utf-8 -*-
"""Phase state classes for ControlVolume-based simulation.

Each phase tracks its composition internally in **moles** (not concentrations).
Concentrations (mol/L, g/L) are derived properties computed from moles and
volume.  This ensures that fluxes (mol/h) can be applied directly without
unit conversion, and that mass conservation is trivially verifiable by
summing moles across all phases.

Three concrete phases are provided:

* :class:`GasPhase` — ideal-gas headspace (derives P, y_i, p_i from n, V, T)
* :class:`LiquidPhase` — dissolved species (derives mol/L, g/L from n, V)
* :class:`SolidPhase` — undissolved material (placeholder for future extensions)

All three satisfy the :class:`Phase` protocol.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Protocol, runtime_checkable

from ..control.descriptors import MutableScalar as _MS

# Ideal gas constant in L·atm/(mol·K)
R_L_ATM_MOL_K = 0.0820574


def _require_positive(name: str, value: float) -> float:
    """Validate that a physical parameter is strictly positive."""
    value = float(value)
    if value <= 0:
        raise ValueError(f"{name} must be positive, got {value}")
    return value


def _require_non_negative(name: str, value: float) -> float:
    """Validate that a physical parameter is non-negative."""
    value = float(value)
    if value < 0:
        raise ValueError(f"{name} must be non-negative, got {value}")
    return value


# ════════════════════════════════════════════════════════════════════════
#  Water vapour saturation pressure (Raoult's law)
# ════════════════════════════════════════════════════════════════════════
# Water vapour in the headspace is always at the Raoult's law saturation
# pressure:  p_H2O = γ_w × x_w × P_sat(T) ≈ P_sat(T)
# for a dilute aqueous solution (γ_w ≈ 1, x_w ≈ 1).
#
# Two P_sat(T) correlations are provided:
#
# 1. Clausius-Clapeyron — matches BSM2 (Rosen & Jeppsson 2006):
#      P_sat(T) = P_ref × exp(ΔH_vap/R × (1/T_ref − 1/T))
#    with P_ref = 0.0313 bar at 298.15 K, ΔH_vap = 43983 J/mol.
#
# 2. Antoine — NIST standard (slightly more accurate, <1% difference):
#      log₁₀(P_mmHg) = 8.07131 − 1730.63/(233.426 + T_°C)
#
# Both agree to <1% across 20–40°C.  Clausius-Clapeyron is the default
# for BSM2 compatibility.

import math as _math

# BSM2 Clausius-Clapeyron parameters
_CC_P_REF_BAR = 0.0313          # P_sat at T_ref (bar)
_CC_T_REF_K = 298.15            # reference temperature (K)
_CC_DH_VAP_OVER_R = 5290.0     # ΔH_vap / R (K) — gives ΔH_vap = 43983 J/mol
_BAR_TO_ATM = 1.0 / 1.01325


def water_vapour_P_sat_atm(T_K: float, method: str = "clausius_clapeyron") -> float:
    """Saturation vapour pressure of water (atm) via Raoult's law.

    This is the analytical water vapour pressure for a headspace in
    equilibrium with a liquid water surface.  It is a thermodynamic
    property of the system — not a dynamic state variable.

    Parameters
    ----------
    T_K : float
        Temperature (Kelvin).
    method : str
        ``"clausius_clapeyron"`` (default, matches BSM2) or ``"antoine"``.

    Returns
    -------
    float
        Saturation vapour pressure in atm.

    Examples
    --------
    >>> water_vapour_P_sat_atm(308.15)  # 35°C
    0.05494...
    """
    if method == "clausius_clapeyron":
        P_bar = _CC_P_REF_BAR * _math.exp(
            _CC_DH_VAP_OVER_R * (1.0 / _CC_T_REF_K - 1.0 / float(T_K)))
        return P_bar * _BAR_TO_ATM
    elif method == "antoine":
        T_C = float(T_K) - 273.15
        T_C = max(1.0, min(100.0, T_C))
        log10_P_mmHg = 8.07131 - 1730.63 / (233.426 + T_C)
        return 10.0 ** log10_P_mmHg / 760.0
    else:
        raise ValueError(
            f"Unknown method {method!r}. Options: "
            f"'clausius_clapeyron', 'antoine'")


# ════════════════════════════════════════════════════════════════════════
#  Protocol
# ════════════════════════════════════════════════════════════════════════

@runtime_checkable
class Phase(Protocol):
    """Contract that every phase must satisfy."""

    @property
    def T_K(self) -> float: ...

    @property
    def V_L(self) -> float: ...

    @property
    def n_mol(self) -> Dict[str, float]:
        """Moles of each species in this phase."""
        ...

    def total_mol(self) -> Dict[str, float]:
        """Return a copy of the moles dict (for diagnostics / snapshots)."""
        ...

    def apply_flux(
        self,
        flux_mol_per_h: Dict[str, float],
        dt_h: float,
        *,
        clamp: bool = True,
    ) -> None:
        """Add ``flux × dt`` moles to this phase.

        Positive flux = material entering this phase.
        Negative flux = material leaving this phase.

        Parameters
        ----------
        clamp : bool
            If True (default), the result is floored at 0.0 per
            species — the long-standing behavior every caller outside
            the clamp-aware ``StepSolver``s relies on. Pass ``False``
            only when the caller has already run its own
            ``clamp_fn`` (including a deliberate ``clamp_fn=None``)
            over the combined delta dict upstream — floor-clamping
            again here would silently defeat ``clamp_fn=None``'s
            purpose of exposing negative inventory for diagnosis.
        """
        ...

    def snapshot(self) -> "Phase":
        """Return an independent deep copy of this phase."""
        ...


# ════════════════════════════════════════════════════════════════════════
#  Gas phase
# ════════════════════════════════════════════════════════════════════════

class GasPhase:
    """Ideal-gas headspace phase.

    State is stored as moles of each gas species.  Pressure, mole fractions,
    and partial pressures are derived on access via the ideal gas law.

    Parameters
    ----------
    n_mol : dict
        Initial moles of each gas species, e.g. ``{"O2": 0.5, "CO2": 0.01, "N2": 1.8}``.
    V_L : float
        Gas-phase volume (litres).
    T_K : float
        Temperature (Kelvin).
    """

    T_K: float = _MS(float, positive=True)
    V_L: float = _MS(float, positive=True)

    def __init__(
        self,
        n_mol: Dict[str, float],
        V_L: float,
        T_K: float,
        properties: Optional[Dict[str, Any]] = None,
    ):
        self._n_mol = {k: float(v) for k, v in n_mol.items()}
        self.properties: Dict[str, Any] = dict(properties) if properties else {}
        # C6 lifecycle gating: set by Simulation.__init__; None when
        # unowned. Lockable setters consult this via raise_if_running.
        self._context = None
        self.V_L = V_L  # validated + stored by MutableScalar descriptor
        self.T_K = T_K

    @property
    def n_mol(self) -> Dict[str, float]:
        return self._n_mol

    # ── Derived gas-phase quantities ───────────────────────────────────

    @property
    def n_total(self) -> float:
        """Total moles of gas."""
        return sum(self._n_mol.values())

    @property
    def P_atm(self) -> float:
        """Total pressure (atm) from ideal gas law: P = nRT/V."""
        if self._V_L <= 0.0:
            return 0.0
        return self.n_total * R_L_ATM_MOL_K * self._T_K / self._V_L

    @property
    def y(self) -> Dict[str, float]:
        """Mole fractions (dimensionless)."""
        n_tot = self.n_total
        if n_tot <= 0.0:
            return {k: 0.0 for k in self._n_mol}
        return {k: v / n_tot for k, v in self._n_mol.items()}

    @property
    def p_atm(self) -> Dict[str, float]:
        """Partial pressures (atm): p_i = n_i RT / V."""
        if self._V_L <= 0.0:
            return {k: 0.0 for k in self._n_mol}
        factor = R_L_ATM_MOL_K * self._T_K / self._V_L
        return {k: v * factor for k, v in self._n_mol.items()}

    # ── Water vapour (analytical, Raoult's law) ───────────────────────

    @property
    def P_water_sat_atm(self) -> float:
        """Saturation vapour pressure of water (atm) at the gas temperature.

        This is the Raoult's law equilibrium vapour pressure for a
        headspace in contact with a liquid water surface.  It is an
        analytical thermodynamic property, not derived from tracked
        gas-phase water moles.

        Uses the Clausius-Clapeyron equation matching BSM2.
        """
        return water_vapour_P_sat_atm(self._T_K)

    @property
    def P_total_wet_atm(self) -> float:
        """Total headspace pressure including analytical water vapour (atm).

        ``P_total_wet = P_dry_species + P_water_sat(T)``

        This is the physically correct total pressure that a gas outlet
        should use for its driving force, since the headspace is always
        saturated with water vapour.

        :attr:`P_atm` returns only the pressure from tracked species
        (the "dry" pressure).  This property adds the analytical water
        vapour contribution.
        """
        return self.P_atm + self.P_water_sat_atm

    # ── Protocol methods ───────────────────────────────────────────────

    def total_mol(self) -> Dict[str, float]:
        return dict(self._n_mol)

    def apply_flux(
        self,
        flux_mol_per_h: Dict[str, float],
        dt_h: float,
        *,
        clamp: bool = True,
    ) -> None:
        dt = float(dt_h)
        for species, rate in flux_mol_per_h.items():
            current = self._n_mol.get(species, 0.0)
            new_value = current + float(rate) * dt
            self._n_mol[species] = max(0.0, new_value) if clamp else new_value

    def snapshot(self) -> "GasPhase":
        return GasPhase(
            n_mol=dict(self._n_mol),
            V_L=self._V_L,
            T_K=self._T_K,
            properties=copy.deepcopy(self.properties),
        )

    def __repr__(self):
        return (f"GasPhase(n_total={self.n_total:.4g} mol, "
                f"P={self.P_atm:.4g} atm, V={self._V_L:.4g} L, T={self._T_K:.1f} K)")


# ════════════════════════════════════════════════════════════════════════
#  Liquid phase
# ════════════════════════════════════════════════════════════════════════

class LiquidPhase:
    """Liquid phase with dissolved species.

    Internal accounting is in moles.  Concentrations are derived from
    ``n / V``.  Mass concentrations (g/L) require a molecular-weight
    lookup, provided via the ``mw`` dict.

    Parameters
    ----------
    n_mol : dict
        Initial moles of each dissolved species.
    V_L : float
        Liquid volume (litres).
    T_K : float
        Temperature (Kelvin).
    mw : dict or None
        Species → molecular weight (g/mol).  Used for g/L reporting.
        Species not in this dict will report ``None`` for mass concentrations.
    speciation : dict or None
        Last speciation output (pH, ion distribution, etc.).
        Updated externally by the speciation engine.
    properties : dict or None
        Physical property values (viscosity, density, etc.).
        Updated externally by property models.
    """

    T_K: float = _MS(float, positive=True)
    V_L: float = _MS(float, positive=True)

    def __init__(
        self,
        n_mol: Dict[str, float],
        V_L: float,
        T_K: float,
        mw: Optional[Dict[str, float]] = None,
        speciation: Optional[Dict[str, Any]] = None,
        properties: Optional[Dict[str, Any]] = None,
    ):
        self._n_mol = {k: float(v) for k, v in n_mol.items()}
        self.mw = dict(mw) if mw else {}
        self.speciation = dict(speciation) if speciation else {}
        self.properties = dict(properties) if properties else {}
        # C6 lifecycle gating (see GasPhase.__init__ for rationale).
        self._context = None
        self.V_L = V_L
        self.T_K = T_K

    @property
    def n_mol(self) -> Dict[str, float]:
        return self._n_mol

    # ── Derived liquid-phase quantities ────────────────────────────────

    @property
    def concentrations_mol_L(self) -> Dict[str, float]:
        """Molar concentrations (mol/L)."""
        if self._V_L <= 0.0:
            return {k: 0.0 for k in self._n_mol}
        return {k: v / self._V_L for k, v in self._n_mol.items()}

    @property
    def concentrations_g_L(self) -> Dict[str, float]:
        """Mass concentrations (g/L).  Requires ``mw`` for each species."""
        if self._V_L <= 0.0:
            return {}
        out = {}
        for k, n in self._n_mol.items():
            mw_k = self.mw.get(k)
            if mw_k is not None and mw_k > 0.0:
                out[k] = (n * mw_k) / self._V_L
        return out

    @property
    def pH(self) -> float:
        """Current pH derived from ``n_mol["H+"] / V_L``.

        Raises
        ------
        ValueError
            If ``"H+"`` is missing from ``n_mol``, or if either
            ``n_mol["H+"]`` or ``V_L`` is non-positive. A
            speciation-less control volume has no defined pH; callers
            must check ``"H+" in phase.n_mol`` before accessing if
            they want to handle that case gracefully.
        """
        if "H+" not in self._n_mol:
            raise ValueError(
                "LiquidPhase.pH: 'H+' is not present in n_mol — "
                "no speciation engine has populated derived species "
                "on this phase. Check `\"H+\" in phase.n_mol` first "
                "if this CV may not have speciation attached."
            )
        H_mol = self._n_mol["H+"]
        if H_mol <= 0.0 or self._V_L <= 0.0:
            raise ValueError(
                f"LiquidPhase.pH: cannot derive pH from "
                f"n_mol['H+']={H_mol!r}, V_L={self._V_L!r}. "
                f"Both must be strictly positive."
            )
        from math import log10
        return -log10(H_mol / self._V_L)

    @property
    def ionic_strength(self) -> Optional[float]:
        val = self.speciation.get("IonicStrength")
        return float(val) if val is not None else None

    # ── Engine-only privileged writeback ──────────────────────────────

    def _refresh_derived(self, values: Dict[str, float]) -> None:
        """Engine-only path: overwrite n_mol entries with values from
        a speciation solve.

        Intended caller is :meth:`BisectionChemicalEquilibriumEngine.solve` after a
        successful equilibrium solve. ``values`` is a dict of
        ``{species_id: n_mol}`` in absolute mole units (not
        concentration); the engine multiplies by ``V_L`` before
        calling. The method bulk-writes the given keys; species not
        in ``values`` are left untouched.

        The underscore prefix and this docstring are the only access
        controls — Python convention. Don't call this from rate-law
        code, link transfer code, or boundary code; those should use
        :meth:`apply_flux`, which is the public mutation path that
        preserves conservation accounting.
        """
        for species_id, n in values.items():
            self._n_mol[species_id] = float(n)

    # ── Protocol methods ───────────────────────────────────────────────

    def total_mol(self) -> Dict[str, float]:
        return dict(self._n_mol)

    def apply_flux(
        self,
        flux_mol_per_h: Dict[str, float],
        dt_h: float,
        *,
        clamp: bool = True,
    ) -> None:
        dt = float(dt_h)
        for species, rate in flux_mol_per_h.items():
            current = self._n_mol.get(species, 0.0)
            new_value = current + float(rate) * dt
            self._n_mol[species] = max(0.0, new_value) if clamp else new_value

    def snapshot(self) -> "LiquidPhase":
        return LiquidPhase(
            n_mol=dict(self._n_mol),
            V_L=self._V_L,
            T_K=self._T_K,
            mw=dict(self.mw),
            speciation=copy.deepcopy(self.speciation),
            properties=copy.deepcopy(self.properties),
        )

    def __repr__(self):
        n_species = len(self._n_mol)
        # pH access can raise if speciation isn't attached; repr must
        # be safe to call on any phase, so guard the lookup.
        try:
            pH_repr = f"{self.pH:.3f}"
        except ValueError:
            pH_repr = "n/a"
        return (f"LiquidPhase({n_species} species, "
                f"V={self._V_L:.4g} L, T={self._T_K:.1f} K, "
                f"pH={pH_repr})")


# ════════════════════════════════════════════════════════════════════════
#  Solid phase (placeholder for future extensions)
# ════════════════════════════════════════════════════════════════════════

class SolidPhase:
    """Solid phase (undissolved substrate, precipitate, carrier, etc.).

    Parameters
    ----------
    n_mol : dict
        Moles of each solid-phase species.
    V_L : float
        Volume occupied by the solid (litres).  For packed beds or
        suspensions, this is the solid volume, not the bulk volume.
    T_K : float
        Temperature (Kelvin).
    mw : dict or None
        Species → molecular weight (g/mol).
    surface_area_m2 : float
        Total external surface area (m²).  Relevant for dissolution and
        heterogeneous reaction rate calculations.
    particle_diameter_m : float or None
        Representative particle diameter (m).
    porosity : float
        Intra-particle porosity (0–1, dimensionless).
    """

    T_K: float = _MS(float, positive=True)
    V_L: float = _MS(float, non_negative=True)

    def __init__(
        self,
        n_mol: Dict[str, float],
        V_L: float,
        T_K: float,
        mw: Optional[Dict[str, float]] = None,
        surface_area_m2: float = 0.0,
        particle_diameter_m: Optional[float] = None,
        porosity: float = 0.0,
        properties: Optional[Dict[str, Any]] = None,
    ):
        self._n_mol = {k: float(v) for k, v in n_mol.items()}
        self.mw = dict(mw) if mw else {}
        self.surface_area_m2 = _require_non_negative("SolidPhase.surface_area_m2", surface_area_m2)
        self.particle_diameter_m = float(particle_diameter_m) if particle_diameter_m is not None else None
        self.porosity = float(porosity)
        self.properties: Dict[str, Any] = dict(properties) if properties else {}
        # C6 lifecycle gating (see GasPhase.__init__ for rationale).
        self._context = None
        self.V_L = V_L
        self.T_K = T_K

    @property
    def n_mol(self) -> Dict[str, float]:
        return self._n_mol

    @property
    def mass_g(self) -> Dict[str, float]:
        """Mass of each solid species (g)."""
        out = {}
        for k, n in self._n_mol.items():
            mw_k = self.mw.get(k)
            if mw_k is not None and mw_k > 0.0:
                out[k] = n * mw_k
        return out

    def total_mol(self) -> Dict[str, float]:
        return dict(self._n_mol)

    def apply_flux(
        self,
        flux_mol_per_h: Dict[str, float],
        dt_h: float,
        *,
        clamp: bool = True,
    ) -> None:
        dt = float(dt_h)
        for species, rate in flux_mol_per_h.items():
            current = self._n_mol.get(species, 0.0)
            new_value = current + float(rate) * dt
            self._n_mol[species] = max(0.0, new_value) if clamp else new_value

    def snapshot(self) -> "SolidPhase":
        return SolidPhase(
            n_mol=dict(self._n_mol),
            V_L=self._V_L,
            T_K=self._T_K,
            mw=dict(self.mw),
            surface_area_m2=self.surface_area_m2,
            particle_diameter_m=self.particle_diameter_m,
            porosity=self.porosity,
            properties=copy.deepcopy(self.properties),
        )

    def __repr__(self):
        n_species = len(self._n_mol)
        total_g = sum(self.mass_g.values())
        return (f"SolidPhase({n_species} species, "
                f"{total_g:.4g} g, A={self.surface_area_m2:.4g} m²)")
