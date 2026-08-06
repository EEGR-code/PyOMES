# -*- coding: utf-8 -*-
"""PHREEQCChemicalEquilibriumEngine — PHREEQC-backed speciation engine.

Wraps `phreeqpython` (v1.6+) to satisfy
:class:`~PyOMES.chemical_equilibrium.protocols.ChemicalEquilibriumEngineProtocol` as a black-box
external engine.  The engine is a drop-in replacement for
:class:`~PyOMES.chemical_equilibrium.nr_engine.NRChemicalEquilibriumEngine` wherever PHREEQC's
thermodynamic database, activity model, or ion-pair library is preferred.

Two translation layers bridge the naming conventions:

- **Input** (`component_map`): PyOMES component ID → PHREEQC element key,
  e.g. ``{"CO2": "C", "NH3": "N(-3)", "H2S": "S(-2)", "SO4": "S(6)"}``.
- **Output** (`species_map`): PHREEQC species name → PyOMES species name,
  applied to every key in ``sol.species`` before writing back to ``phase.n_mol``.
  The built-in default :func:`phreeqc_to_vlsim` converts charge-number notation
  (``Ca+2``) to repeated-symbol notation (``Ca++``) and applies a small
  exceptions dict for species like ``Fe+2`` → ``Fe2+``.

Dependencies
------------
`phreeqpython` is an **optional** dependency (``pip install PyOMES[phreeqc]``).
A clear :exc:`ImportError` with installation instructions is raised at
instantiation when the package is absent.

Usage
-----
::

    from PyOMES.chemical_equilibrium.phreeqc_engine import PHREEQCChemicalEquilibriumEngine

    engine = PHREEQCChemicalEquilibriumEngine(
        {"CO2": 1.0, "NH3": 0.5},          # priming composition, mmol/L
        component_map={"CO2": "C", "NH3": "N(-3)"},
    )
    result = engine.solve(totals={"CO2": 0.010, "NH3": 0.004})  # mol/L
    print(result.pH, result.ionic_strength)
"""
from __future__ import annotations

import re
from typing import Any, Callable, Dict, FrozenSet, Optional

from .protocols import EquilibriumResult

# ---------------------------------------------------------------------------
# Name-translation utilities (no phreeqpython dependency)
# ---------------------------------------------------------------------------

_PHREEQC_SPECIES_EXCEPTIONS: Dict[str, str] = {
    "Fe+2": "Fe2+",
    "Fe+3": "Fe3+",
}

_CHARGE_REGEX = re.compile(r"^(.*?)([+-])(\d+)$")
_OXIDATION_STATE_RE = re.compile(r"\([^)]+\)")


def phreeqc_to_vlsim(phreeqc_name: str) -> str:
    """Convert a PHREEQC species name to the PyOMES naming convention.

    Rules applied in order:

    1. **Exceptions dict**: ``Fe+2`` → ``Fe2+``, ``Fe+3`` → ``Fe3+``.
    2. **Structural regex**: charge-number suffix → repeated symbol:
       ``Ca+2`` → ``Ca++``, ``SO4-2`` → ``SO4--``.
    3. **Pass-through**: singly-charged or un-matched names (``Na+``,
       ``OH-``, ``HCO3-``, ``H+``) are returned unchanged.
    """
    if phreeqc_name in _PHREEQC_SPECIES_EXCEPTIONS:
        return _PHREEQC_SPECIES_EXCEPTIONS[phreeqc_name]
    m = _CHARGE_REGEX.match(phreeqc_name)
    if m:
        base, sign, n_str = m.groups()
        return base + sign * int(n_str)
    return phreeqc_name


def _strip_oxidation_state(key: str) -> str:
    """Remove PHREEQC oxidation-state notation for use with ``sol.change()``.

    Examples: ``"S(-2)"`` → ``"S"``, ``"C(4)"`` → ``"C"``, ``"Ca"`` → ``"Ca"``.
    Required because the PHREEQC CHANGE command accepts only bare element names.
    """
    return _OXIDATION_STATE_RE.sub("", key)


# ---------------------------------------------------------------------------
# PHREEQCChemicalEquilibriumEngine
# ---------------------------------------------------------------------------

class PHREEQCChemicalEquilibriumEngine:
    """PHREEQC-backed black-box speciation engine.

    Satisfies :class:`~PyOMES.chemical_equilibrium.protocols.ChemicalEquilibriumEngineProtocol`.
    Does **not** satisfy Gray/White-box protocols — wrap with
    :class:`~PyOMES.chemical_equilibrium.numerical_gradient.NumericalGradientEquilibriumEngine` to
    obtain ``jacobian_dz_dy()``.

    Parameters
    ----------
    components : dict
        ``{vlsim_component_id: initial_concentration_mmol_L}`` — initial
        composition used for the priming solve at construction.  Every key
        must appear in ``component_map``.
    component_map : dict
        ``{vlsim_component_id: phreeqc_element_key}`` — maps PyOMES component
        IDs to PHREEQC element keys used in ``add_solution_raw`` / ``change()``.
        Use oxidation-state notation where needed (e.g. ``"S(-2)"``).
    species_map : callable or None
        ``(phreeqc_name: str) -> vlsim_name: str``.  Applied to every key in
        ``sol.species`` when building the output dict and the
        algebraic-species cache.  Pass ``None`` to disable translation (PHREEQC
        names are used verbatim).  Default: :func:`phreeqc_to_vlsim`.
    T_C : float
        Operating temperature (°C).  Default 25.0.
    use_warmstart : bool
        Compatibility flag retained for callers that already pass it.
        PyOMES always creates a fresh PHREEQC solution for each solve because
        ``phreeqpython.Solution.change()`` is not equivalent to resetting
        absolute component totals, especially for oxidation-state keys such as
        ``C(4)`` and ``N(-3)``.  Default ``True`` for API compatibility.
    """

    def __init__(
        self,
        components: Dict[str, float],
        *,
        component_map: Dict[str, str],
        species_map: Optional[Callable[[str], str]] = phreeqc_to_vlsim,
        T_C: float = 25.0,
        use_warmstart: bool = True,
    ) -> None:
        try:
            from phreeqpython import PhreeqPython as _PhreeqPython
        except ImportError:
            raise ImportError(
                "PHREEQCChemicalEquilibriumEngine requires phreeqpython>=1.6.  "
                "Install with:\n"
                "    pip install phreeqpython\n"
                "or: pip install 'PyOMES[phreeqc]'"
            ) from None

        self._pp = _PhreeqPython()
        self.component_map: Dict[str, str] = dict(component_map)
        self.species_map: Optional[Callable[[str], str]] = species_map
        self.T_C: float = float(T_C)
        self.use_warmstart: bool = bool(use_warmstart)
        self.n_solve_calls: int = 0
        self._sol = None
        self._priming_components: Dict[str, float] = dict(components)

        # Priming solve — components values are already in mmol/L
        priming_sol = self._fresh_solution(
            {k: float(v) for k, v in components.items()}, T_C
        )

        # Cache algebraic species (conservative superset from the priming point)
        raw_keys = list(priming_sol.species.keys())
        if species_map is not None:
            self._algebraic_species_cache: FrozenSet[str] = frozenset(
                species_map(k) for k in raw_keys
            )
        else:
            self._algebraic_species_cache = frozenset(raw_keys)

        priming_sol.forget()

    # ------------------------------------------------------------------
    # ChemicalEquilibriumEngineProtocol
    # ------------------------------------------------------------------

    def solve(self, **kwargs: Any) -> EquilibriumResult:
        """Solve aqueous equilibrium via PHREEQC.

        Call patterns::

            # Phase-based (primary path from ControlVolume.advance)
            result = engine.solve(phases={"liquid": liq_phase})
            result.apply_to_phases({"liquid": liq_phase})  # explicit commit

            # Direct (testing / scripting)
            result = engine.solve(totals={"CO2": 0.010, "NH3": 0.004}, T_K=298.15)

        Parameters
        ----------
        phases : dict, optional
            ``{"liquid": phase_obj}`` — reads component totals from
            ``phase.n_mol / phase.V_L`` for each key in ``component_map``.
        totals : dict, optional
            ``{component_id: mol_L}`` — direct total concentrations.
        T_K : float, optional
            Temperature override (K).

        Returns
        -------
        EquilibriumResult
            ``pH``, ``logH``, ``ionic_strength`` are the meta fields;
            ``species_mol_L`` carries one entry per PHREEQC species (name
            translated by ``species_map``). ``solve()`` makes no
            side-effecting writes to phases — call
            ``result.apply_to_phases(phases)`` explicitly to commit.
        """
        phases = kwargs.pop("phases", None)
        T_K_override = kwargs.pop("T_K", None)
        T_C = float(T_K_override) - 273.15 if T_K_override is not None else self.T_C

        if phases is not None:
            totals_mol_L = self._read_from_phases(phases)
        else:
            totals_mol_L = dict(kwargs.pop("totals", {}))

        self.n_solve_calls += 1

        # mol/L → mmol/L, keyed by PHREEQC element key
        composition_mmol: Dict[str, float] = {}
        for comp_id, mol_L in totals_mol_L.items():
            phreeqc_key = self.component_map.get(comp_id)
            if phreeqc_key is not None:
                composition_mmol[phreeqc_key] = float(mol_L) * 1000.0

        # Convert PHREEQC keys back: _fresh_solution expects {vlsim_id: mmol_L}.
        # Fresh solves are the correctness baseline.  phreeqpython's
        # Solution.change() path mutates an existing solution and does not
        # reliably match a clean absolute-total solve.
        totals_mmol = {k: v * 1000.0 for k, v in totals_mol_L.items()}
        sol = self._fresh_solution(totals_mmol, T_C)

        out = self._build_output(sol)

        sol.forget()

        species_mol_L = {
            k: float(v) for k, v in out.items()
            if k not in ("pH", "logH", "IonicStrength")
        }
        return EquilibriumResult(
            pH=float(out["pH"]),
            logH=out.get("logH"),
            ionic_strength=out.get("IonicStrength"),
            species_mol_L=species_mol_L,
        )

    def algebraic_species(self) -> FrozenSet[str]:
        """Species IDs written to ``phase.n_mol`` by this engine.

        Returns a conservative superset discovered during the priming solve:
        all species known to PHREEQC at the initial conditions, mapped through
        ``species_map``.  Species with negligible concentration at the priming
        point are included but will have near-zero values in practice.
        """
        return self._algebraic_species_cache

    def reset_cache(self) -> None:
        """Discard any retained PHREEQC ``Solution`` object.

        Current solves always use fresh PHREEQC solutions, so this is normally
        a no-op retained for protocol and API compatibility.
        """
        if self._sol is not None:
            self._sol.forget()
            self._sol = None

    def reset_counters(self) -> None:
        """Reset ``n_solve_calls`` to zero."""
        self.n_solve_calls = 0

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _fresh_solution(self, components_mmol: Dict[str, float], T_C: float):
        """Create a new PHREEQC Solution.

        Parameters
        ----------
        components_mmol : dict
            ``{vlsim_component_id: mmol_L}`` — values in mmol/L.
        T_C : float
            Temperature (°C).
        """
        raw: Dict[str, Any] = {}
        for comp_id, mmol_L in components_mmol.items():
            phreeqc_key = self.component_map.get(comp_id)
            if phreeqc_key is not None and mmol_L > 0.0:
                raw[phreeqc_key] = float(mmol_L)
        raw["temp"] = T_C
        raw["pH"] = "7 charge"
        raw["units"] = "mmol/L"
        return self._pp.add_solution_raw(raw)

    def _build_output(self, sol) -> Dict[str, Any]:
        """Extract ``pH``, ``logH``, ``IonicStrength``, and all species from *sol*."""
        out: Dict[str, Any] = {
            "pH": float(sol.pH),
            "logH": float(-sol.pH),
            "IonicStrength": float(sol.I),
        }
        for phreeqc_name, conc_mol_L in sol.species.items():
            vlsim_name = (
                self.species_map(phreeqc_name)
                if self.species_map is not None
                else phreeqc_name
            )
            out[vlsim_name] = float(conc_mol_L)
        return out

    def _read_from_phases(self, phases) -> Dict[str, float]:
        """Read component totals (mol/L) from a phases dict.

        Reads ``n_mol.get(comp_id, 0.0) / V_L`` for every component in
        ``component_map``.  This is a direct lookup; it does not sum across
        derived species (unlike NRChemicalEquilibriumEngine which uses tableau
        component membership).  The caller is expected to track total
        component amounts under the component ID key in ``n_mol``.
        """
        liq = phases.get("liquid") if hasattr(phases, "get") else None
        if liq is None:
            return {}
        V_L = float(getattr(liq, "V_L", 1.0))
        if V_L <= 0.0:
            return {}
        n_mol = getattr(liq, "n_mol", {})
        return {
            comp_id: float(n_mol.get(comp_id, 0.0)) / V_L
            for comp_id in self.component_map
        }
