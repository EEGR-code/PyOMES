# -*- coding: utf-8 -*-
"""Equilibrium set: a declarative collection of acid-base equilibria.

An :class:`EquilibriumSet` formally describes which acid-base systems
participate in a speciation charge balance.  Each system is represented
by an :class:`EquilibriumDef` that carries its thermodynamic parameters
(pKa, dH, reference temperature, correction method) and charge balance
role (category, number of active dissociation steps).

Usage
-----
Start from an empty set and add systems:

>>> eq = EquilibriumSet()
>>> eq.set_water(pKw=14.0, correction="van_t_hoff", dH_J_per_mol=55900.0)
>>> eq.add("CO2", category="inorganic_acid",
...        pKas=(6.35, 10.33), n_active=1,
...        correction="van_t_hoff", dH_J_per_mol=(7646.0, 14900.0))
>>> eq.add("S_ac", category="acid", pKas=(4.76,))

Systems can be removed again by name:

>>> eq.remove("S_ac")

Categories
----------
Each equilibrium falls into one of these categories, which determines
how it enters the charge balance:

``"acid"``
    Weak acid (HA → A⁻ + H⁺).  Anionic species contribute negative
    charge.  Handles monoprotic and polyprotic.  Examples: VFAs, H₂S.

``"cation_acid"``
    Cationic acid (BH⁺ → B + H⁺).  The protonated form is a cation.
    Example: NH₄⁺/NH₃.

``"inorganic_acid"``
    Polyprotic inorganic acid.  Same algebra as ``"acid"`` but
    semantically distinct.  Examples: CO₂/carbonate, phosphate.

``"strong_ion"``
    Fully dissociated, no equilibrium.  Fixed charge.  Included for
    completeness; typically handled via the context dict.
    Examples: Na⁺, Cl⁻.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

from ....units import R_J_PER_MOL_K as _R_J

_VALID_CATEGORIES = ("acid", "cation_acid", "inorganic_acid", "strong_ion")
_VALID_CORRECTIONS = ("none", "van_t_hoff")


# ════════════════════════════════════════════════════════════════════════
#  EquilibriumDef
# ════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class EquilibriumDef:
    """One equilibrium system in the charge balance.

    Parameters
    ----------
    name : str
        Unique identifier (e.g. ``"CO2"``, ``"S_ac"``, ``"phosphate"``).
    category : str
        Charge balance role: ``"acid"``, ``"cation_acid"``,
        ``"inorganic_acid"``, or ``"strong_ion"``.
    pKas : tuple of float
        Reference pKa value(s) at ``T_ref_K``.  Length equals the
        number of dissociation steps defined (``n_protons``).
    n_active : int
        Number of dissociation steps used in the charge balance.
        Must be ``<= len(pKas)``.  Defaults to ``len(pKas)``.

        This allows storing full thermodynamic data while using only
        a subset in the charge balance.  For example, CO₂ can be declared
        diprotic with ``pKas=(6.35, 10.33)`` and ``n_active=1``
        so only the first dissociation (HCO₃⁻) participates in the
        charge balance.  The second pKa is still available for other
        calculations (e.g. gas transfer speciation corrections).
    T_ref_K : float
        Reference temperature (K) for the pKa values.
    correction : str
        Temperature correction method: ``"none"`` or ``"van_t_hoff"``.
    dH_J_per_mol : tuple of float
        Van 't Hoff enthalpy (J/mol) per dissociation step.
        Length must match ``pKas``.
    total_key : str
        Key in the context dict for the total concentration of this
        species.  E.g. ``"CT_TIC"`` for CO₂, ``"CT_S_ac"`` for acetate.
    """
    name: str
    category: str
    pKas: Tuple[float, ...]
    n_active: int
    T_ref_K: float = 298.15
    correction: str = "none"
    dH_J_per_mol: Tuple[float, ...] = ()
    total_key: str = ""
    # Ordered ladder of Species objects (most protonated → least protonated).
    # When non-empty, _compute_species_eq emits using Species.id directly
    # instead of synthesised {name}_HA / {name}_A- keys.
    # Set by EquilibriumSet.add(species_refs=...) or engine._add_acid.
    species_refs: Tuple[Any, ...] = ()

    def __post_init__(self):
        if not isinstance(self.pKas, tuple):
            object.__setattr__(self, "pKas", tuple(self.pKas))
        if not isinstance(self.dH_J_per_mol, tuple):
            object.__setattr__(self, "dH_J_per_mol", tuple(self.dH_J_per_mol))
        if len(self.dH_J_per_mol) == 0:
            object.__setattr__(self, "dH_J_per_mol", (0.0,) * len(self.pKas))
        if not isinstance(self.species_refs, tuple):
            object.__setattr__(self, "species_refs", tuple(self.species_refs))

    @property
    def n_protons(self) -> int:
        """Total number of dissociation steps defined."""
        return len(self.pKas)

    def pKas_at_T(self, T_K: float) -> Tuple[float, ...]:
        """Compute all pKa(s) at temperature T_K."""
        if self.correction == "none":
            return self.pKas
        if abs(T_K - self.T_ref_K) < 0.01:
            return self.pKas
        result = []
        for pKa_ref, dH in zip(self.pKas, self.dH_J_per_mol):
            if abs(dH) < 1e-10:
                result.append(pKa_ref)
            else:
                Ka_ref = 10.0 ** (-pKa_ref)
                Ka_T = Ka_ref * math.exp(
                    -dH / _R_J * (1.0 / T_K - 1.0 / self.T_ref_K))
                result.append(-math.log10(max(Ka_T, 1e-30)))
        return tuple(result)

    def active_pKas_at_T(self, T_K: float) -> Tuple[float, ...]:
        """Compute only the active pKa(s) at temperature T_K.

        Returns the first ``n_active`` entries from :meth:`pKas_at_T`.
        These are the values used in the charge balance.
        """
        return self.pKas_at_T(T_K)[:self.n_active]


# ════════════════════════════════════════════════════════════════════════
#  WaterDef (reused within EquilibriumSet)
# ════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class WaterDef:
    """Water autoionisation constant definition."""
    pKw: float = 14.0
    T_ref_K: float = 298.15
    correction: str = "none"
    dH_J_per_mol: float = 0.0

    def Kw_at_T(self, T_K: float) -> float:
        Kw_ref = 10.0 ** (-self.pKw)
        if self.correction == "none" or abs(self.dH_J_per_mol) < 1e-10:
            return Kw_ref
        if abs(T_K - self.T_ref_K) < 0.01:
            return Kw_ref
        return Kw_ref * math.exp(
            -self.dH_J_per_mol / _R_J * (1.0 / T_K - 1.0 / self.T_ref_K))

    def pKw_at_T(self, T_K: float) -> float:
        return -math.log10(max(self.Kw_at_T(T_K), 1e-30))


# ════════════════════════════════════════════════════════════════════════
#  EquilibriumSet
# ════════════════════════════════════════════════════════════════════════

class EquilibriumSet:
    """Ordered collection of equilibrium definitions for a charge balance.

    Build from scratch or load a preset, then add/remove systems.
    """

    def __init__(self, *, T_ref_K: float = 298.15):
        """Create an empty set with a default reference temperature.

        Parameters
        ----------
        T_ref_K : float
            Default reference temperature (K) used when ``add()`` is
            called without specifying ``T_ref_K``.
        """
        self.T_ref_K = T_ref_K
        self._equilibria: Dict[str, EquilibriumDef] = {}
        self._water: WaterDef = WaterDef(pKw=14.0, T_ref_K=T_ref_K)
        self._insertion_order: List[str] = []

    # ── Water ─────────────────────────────────────────────────────────

    def set_water(self, *, pKw=None, correction="none",
                  dH_J_per_mol=None, T_ref_K=None):
        """Define water autoionisation.

        Returns self for method chaining.
        """
        if pKw is None:
            pKw = 14.0

        t_ref = float(T_ref_K) if T_ref_K is not None else self.T_ref_K

        dH_internal = float(dH_J_per_mol) if dH_J_per_mol is not None else 0.0

        correction = correction.lower().strip()
        if correction not in _VALID_CORRECTIONS:
            raise ValueError(
                f"Unknown correction {correction!r} for water. "
                f"Options: {_VALID_CORRECTIONS}")
        if correction == "van_t_hoff" and abs(dH_internal) < 1e-10:
            raise ValueError(
                "correction='van_t_hoff' for water requires nonzero dH")

        self._water = WaterDef(
            pKw=pKw, T_ref_K=t_ref,
            correction=correction, dH_J_per_mol=dH_internal)
        return self

    @property
    def water(self) -> WaterDef:
        """The water autoionisation definition."""
        return self._water

    # ── Add / remove / query ──────────────────────────────────────────

    def add(self, name, *, category, pKas=None,
            n_active=None, correction="none",
            dH_J_per_mol=None, T_ref_K=None,
            total_key=None, species_refs=()):
        """Add or replace an equilibrium definition.

        Parameters
        ----------
        name : str
            Unique identifier for this equilibrium system.
        category : str
            ``"acid"``, ``"cation_acid"``, ``"inorganic_acid"``,
            or ``"strong_ion"``.
        pKas : tuple of float
            Dissociation constant(s), as pKa.
        n_active : int, optional
            Number of dissociation steps active in the charge balance.
            Defaults to ``len(pKas)``.  Set lower to store inactive
            pKas as metadata (e.g. ``n_active=1`` for monoprotic
            CO₂ while storing both pKa₁ and pKa₂).
        correction : str
            ``"none"`` or ``"van_t_hoff"``.
        dH_J_per_mol : tuple of float
            Van 't Hoff enthalpy (J/mol) per dissociation step.
        T_ref_K : float, optional
            Reference temperature (K).  Defaults to the set's own
            ``T_ref_K``.
        total_key : str, optional
            Context dict key for total concentration.  Defaults to
            ``"CT_{name}"`` (e.g. ``"CT_CO2"``, ``"CT_S_ac"``).

        Returns
        -------
        self

        Examples
        --------
        >>> eq.add("CO2", category="inorganic_acid",
        ...        pKas=(6.35, 10.33), n_active=1,
        ...        correction="van_t_hoff",
        ...        dH_J_per_mol=(7646.0, 14900.0))
        >>> eq.add("H2S", category="acid", pKas=(7.0,),
        ...        correction="van_t_hoff",
        ...        dH_J_per_mol=(20000.0,), T_ref_K=293.15)
        >>> eq.add("S_ac", category="acid", pKas=(4.76,))
        """
        import warnings

        # ── Validate category ─────────────────────────────────────
        category = category.lower().strip()
        if category not in _VALID_CATEGORIES:
            raise ValueError(
                f"Unknown category {category!r} for '{name}'. "
                f"Options: {_VALID_CATEGORIES}")

        # ── pKas ────────────────────────────────────────────────────
        if category == "strong_ion":
            pKas_internal = ()
        else:
            if pKas is None:
                raise ValueError(f"add('{name}'): must provide pKas")
            if not isinstance(pKas, (tuple, list)):
                pKas = (pKas,)
            pKas_internal = tuple(float(p) for p in pKas)

        # ── n_active ──────────────────────────────────────────────
        if n_active is None:
            n_active = len(pKas_internal)
        else:
            n_active = int(n_active)
            if n_active < 0:
                raise ValueError(
                    f"n_active must be >= 0, got {n_active} for '{name}'")
            if n_active > len(pKas_internal):
                raise ValueError(
                    f"n_active ({n_active}) cannot exceed number of pKas "
                    f"({len(pKas_internal)}) for '{name}'")

        # ── Reference temperature ─────────────────────────────────
        t_ref = float(T_ref_K) if T_ref_K is not None else self.T_ref_K

        # ── Enthalpy ────────────────────────────────────────────────
        if category == "strong_ion":
            dH_internal = ()
        elif dH_J_per_mol is not None:
            if not isinstance(dH_J_per_mol, (tuple, list)):
                dH_J_per_mol = (dH_J_per_mol,)
            dH_internal = tuple(float(d) for d in dH_J_per_mol)
        else:
            dH_internal = None

        # ── Validate correction + dH ─────────────────────────────
        correction = correction.lower().strip()
        if correction not in _VALID_CORRECTIONS:
            raise ValueError(
                f"Unknown correction {correction!r} for '{name}'. "
                f"Options: {_VALID_CORRECTIONS}")

        if category != "strong_ion":
            if correction == "van_t_hoff":
                if dH_internal is None:
                    raise ValueError(
                        f"correction='van_t_hoff' for '{name}' requires "
                        f"dH_J_per_mol")
                if len(dH_internal) != len(pKas_internal):
                    raise ValueError(
                        f"dH length ({len(dH_internal)}) must match pKas "
                        f"length ({len(pKas_internal)}) for '{name}'")
            elif correction == "none":
                if dH_internal is not None and any(
                        abs(d) > 1e-10 for d in dH_internal):
                    warnings.warn(
                        f"add('{name}'): dH_J_per_mol provided but "
                        f"correction='none' — values will be ignored.",
                        UserWarning, stacklevel=2)
                dH_internal = (0.0,) * len(pKas_internal)

        # ── Default total_key ─────────────────────────────────────
        if total_key is None:
            total_key = f"CT_{name}"

        # ── Store ─────────────────────────────────────────────────
        is_new = name not in self._equilibria
        self._equilibria[name] = EquilibriumDef(
            name=name, category=category,
            pKas=pKas_internal, n_active=n_active,
            T_ref_K=t_ref, correction=correction,
            dH_J_per_mol=dH_internal if dH_internal is not None else (),
            total_key=total_key,
            species_refs=tuple(species_refs))
        if is_new:
            self._insertion_order.append(name)
        return self

    def remove(self, name: str) -> "EquilibriumSet":
        """Remove an equilibrium by name.

        Raises
        ------
        KeyError
            If ``name`` is not in the set.
        """
        if name not in self._equilibria:
            raise KeyError(
                f"Equilibrium '{name}' not found. "
                f"Available: {self.names}")
        del self._equilibria[name]
        self._insertion_order.remove(name)
        return self

    def has(self, name: str) -> bool:
        """Check if an equilibrium is registered."""
        return name in self._equilibria

    def get(self, name: str) -> EquilibriumDef:
        """Get an equilibrium definition by name.

        Raises
        ------
        KeyError
            If ``name`` is not in the set.
        """
        if name not in self._equilibria:
            raise KeyError(
                f"Equilibrium '{name}' not found. "
                f"Available: {self.names}")
        return self._equilibria[name]

    @property
    def names(self) -> List[str]:
        """List of registered equilibrium names (insertion order)."""
        return list(self._insertion_order)

    def __len__(self) -> int:
        return len(self._equilibria)

    def __contains__(self, name: str) -> bool:
        return name in self._equilibria

    def __iter__(self):
        """Iterate over EquilibriumDef objects in insertion order."""
        for name in self._insertion_order:
            yield self._equilibria[name]

    def __getitem__(self, name: str) -> EquilibriumDef:
        return self.get(name)

    # ── Repr ──────────────────────────────────────────────────────────

    def __repr__(self) -> str:
        entries = []
        for eq_def in self:
            active_str = (f"{eq_def.n_active}/{eq_def.n_protons}"
                          if eq_def.n_active < eq_def.n_protons
                          else str(eq_def.n_protons))
            entries.append(
                f"{eq_def.name}({eq_def.category}, "
                f"pKas={eq_def.pKas}, active={active_str}, "
                f"correction={eq_def.correction!r})")
        water_str = (f"water(pKw={self._water.pKw}, "
                     f"correction={self._water.correction!r})")
        return (f"EquilibriumSet([{water_str}, "
                + ", ".join(entries) + "])")

    def summary(self) -> str:
        """Human-readable summary of all registered equilibria."""
        lines = [f"EquilibriumSet ({len(self)} equilibria, "
                 f"T_ref_default={self.T_ref_K} K):"]
        lines.append(f"  Water: pKw={self._water.pKw} at "
                     f"{self._water.T_ref_K} K, "
                     f"correction={self._water.correction!r}")
        for eq_def in self:
            active_str = (f" (active: {eq_def.n_active}/{eq_def.n_protons})"
                          if eq_def.n_active < eq_def.n_protons else "")
            lines.append(
                f"  {eq_def.name:12s}  {eq_def.category:18s}  "
                f"pKas={eq_def.pKas}{active_str}  "
                f"T_ref={eq_def.T_ref_K} K  "
                f"correction={eq_def.correction!r}  "
                f"key={eq_def.total_key!r}")
        return "\n".join(lines)
