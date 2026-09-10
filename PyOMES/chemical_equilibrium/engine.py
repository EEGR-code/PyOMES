# -*- coding: utf-8 -*-
"""Bisection speciation engine — one of three ChemicalEquilibriumEngineProtocol
implementations (the others are :class:`~PyOMES.chemical_equilibrium.nr_engine.NRChemicalEquilibriumEngine`
and :class:`~PyOMES.chemical_equilibrium.phreeqc_engine.PHREEQCChemicalEquilibriumEngine`).

:class:`BisectionChemicalEquilibriumEngine` solves aqueous acid-base equilibria from declared
:class:`~PyOMES.reactions.equilibrium.EquilibriumReaction` instances via 1-D
bisection over the charge balance. It handles single acid-base ladders and
independent ladders that only interact through the charge balance (e.g.
carbonate + ammonia) — not arbitrary cross-component networks or gas-liquid/
solid-liquid folding; use ``NRChemicalEquilibriumEngine`` for those. Chemistry
is fully expressed through the declared reactions; the engine provides no
default chemistry of its own.

Two solve paths are available inside :meth:`BisectionChemicalEquilibriumEngine.solve`:

- **Declared-reaction path** (primary): when the engine carries an
  :class:`~PyOMES.chemistry.equilibria.EquilibriumSet` (built by
  :meth:`~BisectionChemicalEquilibriumEngine.from_reactions`), concentrations are read
  from the phase state and the charge balance is solved via
  :func:`~PyOMES.chemical_equilibrium.acid_base.solve_from_equilibrium_set`.
- **Legacy path** (backward compatibility): when explicit
  ``acid_totals``/``acid_pKas``/``CT_TIC`` keyword arguments are
  supplied without an ``EquilibriumSet``, the solver falls back to
  :func:`~PyOMES.chemical_equilibrium.acid_base.solve_acid_base` directly.

:meth:`BisectionChemicalEquilibriumEngine.solve` returns an immutable
:class:`~PyOMES.chemical_equilibrium.protocols.EquilibriumResult`. It makes no
side-effecting writes to phases; call ``result.apply_to_phases(phases)``
explicitly to commit derived species to ``phase.n_mol``.
"""
from __future__ import annotations

import logging
import warnings

from typing import Any, Dict

import numpy as np

logger = logging.getLogger(__name__)

from .acid_base import solve_acid_base, solve_from_equilibrium_set
from .activity_models import make_activity_model
from .protocols import EquilibriumResult


# Special-case total_key mapping for species whose tracking key deviates
# from the default "CT_{id}" convention used by EquilibriumSet.add().
_TOTAL_KEY_OVERRIDES = {
    "CO2":     "CT_TIC",
    "NH4+":    "CT_NH_T",
    "NH4":     "CT_NH_T",
    "NH3":     "CT_NH_T",
    "TIC":     "CT_TIC",
    "CT_NH_T": "CT_NH_T",
}

# Species written back to phase.n_mol by EquilibriumResult.apply_to_phases()
# for this engine. Fixed since chemistry-unification-3b; deliberately does
# NOT include generic polyprotic-ladder keys (e.g. "{name}_HA") — those were
# never part of the writeback and land in EquilibriumResult.extra instead.
_CANONICAL_WRITEBACK_SPECIES = (
    "H+", "OH-",
    "CO2", "HCO3-", "CO3--",
    "NH3", "NH4+",
    "H3PO4", "H2PO4-", "HPO4--", "PO4---",
    "H2S", "HS-",
    "HSO4-", "SO4--",
    "K+", "Na+", "Cl-",
    "Mg++", "Ca++",
    "Zn++", "Mn++", "Co++", "Mo7O24------",
)

# Meta keys consumed into EquilibriumResult's named fields — excluded from
# both species_mol_L and extra.
_META_KEYS = frozenset({
    "pH", "pH_conc", "logH", "aH", "gamma_H", "gamma_OH",
    "IonicStrength", "charge_residual", "alphas",
})


class BisectionChemicalEquilibriumEngine:
    """Unified speciation engine.

    Chemistry is declared through
    :class:`~PyOMES.reactions.equilibrium.EquilibriumReaction` instances
    bound via :meth:`from_reactions`; the engine solves the charge balance
    for the supplied totals.  No chemistry is hardcoded.

    Construct via :meth:`from_reactions` to bind a pre-built
    :class:`~PyOMES.chemistry.equilibria.EquilibriumSet`.
    """

    def __init__(
        self,
        *,
        use_activity: bool = False,
        activity_model: str = "davies",
        T_C: float = 25.0,
        thermo=None,
        use_warmstart: bool = True,
    ):
        if thermo is not None:
            self._liquid_activity = thermo.liquid_activity
            use_activity = thermo.use_activity
            activity_model = thermo.activity_model
        else:
            self._liquid_activity = None

        self.use_activity = bool(use_activity)
        self.activity_model = str(activity_model)
        self.thermo = thermo
        self.T_C = float(T_C)
        self.use_warmstart = bool(use_warmstart)

        self._logH_last = None
        self._I_last = None
        self.n_solve_calls = 0

        # Gas-liquid/solid-liquid EquilibriumConstraint items supplied to
        # from_reactions() — outside this engine's tableau (this phase
        # doesn't fold them in), but visible for downstream consumers
        # (KineticGasLiquidLink) rather than silently discarded.
        self.cross_phase_constraints: tuple = ()

    # ------------------------------------------------------------------
    # Factory: build an engine from declared equilibrium reactions
    # ------------------------------------------------------------------
    @classmethod
    def from_reactions(
        cls,
        equilibrium_reactions,
        *,
        activity_model: str = "davies",
        use_activity: bool = False,
        T_K: float = 298.15,
        **engine_kwargs,
    ) -> "BisectionChemicalEquilibriumEngine":
        """Build a :class:`BisectionChemicalEquilibriumEngine` from declared equilibrium reactions.

        Each item in ``equilibrium_reactions`` must satisfy
        :class:`~PyOMES.reactions.equilibrium.EquilibriumConstraint`
        (``stoichiometry``/``log_K``/``dH_J_per_mol``/``T_ref_K``) —
        typically an
        :class:`~PyOMES.reactions.equilibrium.EquilibriumReaction`, but
        also :class:`~PyOMES.chemistry.partition.HenryEquilibrium`,
        :class:`~PyOMES.chemistry.partition.KspEquilibrium`, or
        :class:`~PyOMES.chemistry.partition.RaoultEquilibrium`. Each
        item is classified via
        :func:`~PyOMES.reactions.equilibrium.classify_equilibrium_constraint`
        from its stoichiometry's phase tags. Single-phase
        (``"acid_base"``) items are further classified into
        ``"water"``, ``"acid"``, or ``"cation_acid"`` from the charge
        convention, and folded into an
        :class:`~PyOMES.chemistry.equilibria.EquilibriumSet` carrying
        the pKa values and Van 't Hoff temperature parameters. The
        engine is constructed pre-loaded with this set; subsequent
        ``solve(**kwargs)`` calls without an explicit
        ``equilibrium_set`` use this one by default.

        Stoichiometry contract (single-phase / ``"acid_base"`` items):

        - **Water dissociation:** the only reaction whose products
          include both H⁺ and OH⁻. Recognised regardless of the
          species id of the consumed water (defaults to ``"H2O"``).
        - **Acid (``HA ⇌ A⁻ + H⁺``):** one reactant with charge 0
          (or any non-positive charge), one anionic product, and H⁺.
        - **Cation acid (``BH⁺ ⇌ B + H⁺``):** one reactant with
          positive charge, one neutral or less-positive product, and
          H⁺.

        **Gas-liquid / solid-liquid items are not silently dropped.**
        Items that classify as ``"gas_liquid"`` or ``"solid_liquid"``
        (e.g. gas-liquid partition declarations like ``CO2(gas) ⇌
        CO2aq(liquid)``, or a ``HenryEquilibrium``/``KspEquilibrium``/
        ``RaoultEquilibrium``) are outside this engine's tableau in
        this phase (this phase does not fold them in) but are visible
        via ``engine.cross_phase_constraints`` for downstream
        consumers (``KineticGasLiquidLink``) rather than discarded.

        Parameters
        ----------
        equilibrium_reactions : iterable of EquilibriumConstraint
            Equilibrium constraints to incorporate. Order is preserved.
        activity_model, use_activity, T_K : misc
            Forwarded to :class:`BisectionChemicalEquilibriumEngine.__init__`.
        **engine_kwargs : dict
            Additional keyword arguments to forward to
            :class:`BisectionChemicalEquilibriumEngine.__init__`.

        Returns
        -------
        BisectionChemicalEquilibriumEngine
            Engine ready to ``solve(...)``. The constructed
            ``EquilibriumSet`` lives on ``engine._equilibrium_set``
            and is used as the default whenever ``solve()`` does not
            receive an explicit ``equilibrium_set``. Gas-liquid/
            solid-liquid items live on ``engine.cross_phase_constraints``.
        """
        from collections import defaultdict
        from ..chemistry.equilibria import EquilibriumSet
        from ..reactions.equilibrium import (
            EquilibriumConstraint, EquilibriumReaction,
            classify_equilibrium_constraint,
        )

        eq_set = EquilibriumSet(T_ref_K=T_K)
        water_set = False
        # Acid reactions keyed by their resolved name (total_id or acid species id).
        # Groups with more than one reaction share a total_id and form a polyprotic
        # ladder that must be merged into a single EquilibriumDef.
        acid_groups: "dict[str, list]" = defaultdict(list)
        cross_phase_constraints: list = []

        for rxn in equilibrium_reactions:
            if not isinstance(rxn, EquilibriumConstraint):
                raise ValueError(
                    f"BisectionChemicalEquilibriumEngine.from_reactions: item "
                    f"{getattr(rxn, 'label', '?')!r} is "
                    f"{type(rxn).__name__}, expected an "
                    "EquilibriumConstraint-conforming type "
                    "(EquilibriumReaction, HenryEquilibrium, "
                    "KspEquilibrium, RaoultEquilibrium, ...)."
                )
            # Gas-liquid/solid-liquid constraints are partition/
            # precipitation declarations consumed elsewhere
            # (KineticGasLiquidLink, NRChemicalEquilibriumEngine) — not
            # single-phase acid-base equilibria. Expose, don't discard.
            if classify_equilibrium_constraint(rxn) != "acid_base":
                cross_phase_constraints.append(rxn)
                continue
            if not isinstance(rxn, EquilibriumReaction):
                raise ValueError(
                    f"BisectionChemicalEquilibriumEngine.from_reactions: single-phase "
                    f"(acid-base) constraint {getattr(rxn, 'label', '?')!r} "
                    f"is {type(rxn).__name__}; only EquilibriumReaction "
                    "carries acid-base pKa semantics."
                )
            kind = _classify_equilibrium(rxn)
            if kind == "water":
                if water_set:
                    raise ValueError(
                        "BisectionChemicalEquilibriumEngine.from_reactions: more than "
                        "one water-dissociation reaction supplied."
                    )
                _add_water(eq_set, rxn)
                water_set = True
            else:
                reactants = [
                    e for e in rxn.stoichiometry
                    if e.coefficient < 0 and e.species.id not in _SOLVENT_IDS
                ]
                acid_id = reactants[0].species.id
                name = getattr(rxn, "total_id", None) or acid_id
                acid_groups[name].append((kind, rxn))

        for name, group in acid_groups.items():
            if len(group) == 1:
                kind, rxn = group[0]
                _add_acid(eq_set, rxn, category=kind)
            else:
                _add_polyprotic_acid(eq_set, name, group)

        T_C = float(T_K) - 273.15
        engine = cls(
            use_activity=use_activity,
            activity_model=activity_model,
            T_C=T_C,
            **engine_kwargs,  # thermo, use_warmstart, etc.
        )
        engine._equilibrium_set = eq_set
        engine.cross_phase_constraints = tuple(cross_phase_constraints)
        return engine

    def solve(self, **kwargs) -> EquilibriumResult:
        # Inject the engine's bound EquilibriumSet when the caller doesn't
        # supply one.
        if "equilibrium_set" not in kwargs:
            bound = getattr(self, "_equilibrium_set", None)
            if bound is not None:
                kwargs["equilibrium_set"] = bound

        # Read totals + strong ions from `phases` when supplied.
        phases = kwargs.pop("phases", None)
        if phases is not None:
            _populate_totals_from_phases(
                phases, kwargs, kwargs.get("equilibrium_set"),
            )
            _populate_strong_ions_from_phases(phases, kwargs)

        self.n_solve_calls += 1

        T_K_override = kwargs.pop("T_K", None)
        T_K = float(T_K_override) if T_K_override is not None else 273.15 + self.T_C

        if self.use_warmstart:
            if self._logH_last is not None and "logH_guess" not in kwargs:
                kwargs["logH_guess"] = self._logH_last
            if self._I_last is not None and "I_init" not in kwargs:
                kwargs["I_init"] = self._I_last

        am = (
            self._liquid_activity
            if self._liquid_activity is not None
            else make_activity_model(self.use_activity, self.activity_model)
        )
        eq_set = kwargs.pop("equilibrium_set", None)

        if eq_set is not None:
            # ----------------------------------------------------------
            # Declared-reaction path: EquilibriumSet drives the chemistry
            # ----------------------------------------------------------
            concentrations: Dict[str, float] = {}
            acid_totals = kwargs.get("acid_totals", {})
            for eq_def in eq_set:
                if eq_def.total_key in kwargs:
                    concentrations[eq_def.total_key] = float(kwargs[eq_def.total_key])
                elif eq_def.total_key == "CT_TIC":
                    concentrations["CT_TIC"] = float(kwargs.get("CT_TIC", 0.0))
                elif eq_def.total_key == "CT_NH_T":
                    concentrations["CT_NH_T"] = float(kwargs.get("CT_NH_T", 0.0))
                elif eq_def.total_key == "CT_P":
                    concentrations["CT_P"] = float(kwargs.get("CT_P", 0.0))
                elif eq_def.total_key == "CT_SO4":
                    concentrations["CT_SO4"] = float(kwargs.get("CT_SO4", 0.0))
                elif eq_def.name in acid_totals:
                    concentrations[eq_def.total_key] = float(acid_totals[eq_def.name])

            strong_ions = dict(kwargs.get("strong_kwargs", {}))
            for si_key in (
                "CT_K", "CT_Na", "CT_Cl", "CT_NO3",
                "CT_Mg", "CT_Ca", "CT_Zn", "CT_Mn",
                "CT_Co", "CT_Mo7O24",
                "CT_cation", "CT_anion",
            ):
                if si_key in kwargs and si_key not in strong_ions:
                    strong_ions[si_key] = float(kwargs[si_key])

            out = solve_from_equilibrium_set(
                equilibrium_set=eq_set,
                concentrations=concentrations,
                strong_ions=strong_ions,
                T_K=T_K,
                activity_model=am,
                pH_min=kwargs.get("pH_min", -0.5),
                pH_max=kwargs.get("pH_max", 20.0),
                n_scan=kwargs.get("n_scan", 200),
                tol=kwargs.get("tol", 1e-12),
                max_outer=kwargs.get("max_outer", 20),
                I_init=kwargs.get("I_init", 0.01),
                I_tol=kwargs.get("I_tol", 1e-8),
                damping=kwargs.get("damping", 0.3),
                logH_guess=kwargs.get("logH_guess"),
            )
        else:
            # ----------------------------------------------------------
            # Legacy path: explicit acid_totals/acid_pKas/CT_* kwargs
            # ----------------------------------------------------------
            out = solve_acid_base(
                acid_totals=kwargs.get("acid_totals", {}),
                acid_pKas=kwargs.get("acid_pKas", {}),
                CT_TIC=kwargs.get("CT_TIC", 0.0),
                CT_NH_T=kwargs.get("CT_NH_T", 0.0),
                CT_P=kwargs.get("CT_P", 0.0),
                CT_K=kwargs.get("CT_K", 0.0),
                CT_Na=kwargs.get("CT_Na", 0.0),
                CT_Cl=kwargs.get("CT_Cl", 0.0),
                CT_NO3=kwargs.get("CT_NO3", 0.0),
                CT_SO4=kwargs.get("CT_SO4", 0.0),
                CT_Mg=kwargs.get("CT_Mg", 0.0),
                CT_Ca=kwargs.get("CT_Ca", 0.0),
                CT_Zn=kwargs.get("CT_Zn", 0.0),
                CT_Mn=kwargs.get("CT_Mn", 0.0),
                CT_Co=kwargs.get("CT_Co", 0.0),
                CT_Mo7O24=kwargs.get("CT_Mo7O24", 0.0),
                Kw=kwargs.get("Kw", 1e-14),
                pKa1_TIC=kwargs.get("pKa1_TIC", 6.35),
                pKa2_TIC=kwargs.get("pKa2_TIC", 10.33),
                pKa_NH=kwargs.get("pKa_NH", 9.25),
                pKa1_P=kwargs.get("pKa1_P", 2.15),
                pKa2_P=kwargs.get("pKa2_P", 7.20),
                pKa3_P=kwargs.get("pKa3_P", 12.35),
                pKa_HSO4=kwargs.get("pKa_HSO4", 1.99),
                pH_min=kwargs.get("pH_min", -0.5),
                pH_max=kwargs.get("pH_max", 20.0),
                n_scan=kwargs.get("n_scan", 200),
                tol=kwargs.get("tol", 1e-12),
                activity_model=am,
                T_K=T_K,
                max_outer=kwargs.get("max_outer", 20),
                I_init=kwargs.get("I_init", 0.01),
                I_tol=kwargs.get("I_tol", 1e-6),
                damping=kwargs.get("damping", 0.3),
                logH_guess=kwargs.get("logH_guess"),
            )

        if hasattr(out, "as_dict"):
            out = out.as_dict()
        else:
            out = dict(out)

        if self.use_warmstart:
            if "logH" in out and out["logH"] is not None:
                self._logH_last = float(out["logH"])
            if "IonicStrength" in out and out["IonicStrength"] is not None:
                self._I_last = float(out["IonicStrength"])

        # ------------------------------------------------------------------
        # Output normalization
        # ------------------------------------------------------------------
        # The declared-reaction path emits "CO2" (canonical Species.id);
        # the legacy path emits "CO2aq" (from solve_acid_base directly).
        # Only normalise truly old label variants, not "CO2" itself.
        if "CO2aq" not in out:
            if "CO2(aq)" in out:
                out["CO2aq"] = out["CO2(aq)"]
            elif "CO2_aq" in out:
                out["CO2aq"] = out["CO2_aq"]

        if "pH" not in out or out.get("pH") is None:
            H = None
            for k in ("H_molL", "H+"):
                if k in out and out.get(k) is not None:
                    try:
                        H = float(out[k])
                        break
                    except (TypeError, ValueError):
                        pass
            if H is not None and H > 0.0:
                out["pH"] = float(-np.log10(H))

        # ------------------------------------------------------------------
        # alphas channel
        # ------------------------------------------------------------------
        alphas: Dict[str, float] = {}

        def _alpha(mol_conc, total):
            if mol_conc is None or total is None:
                return None
            try:
                t = float(total)
                if t <= 0.0:
                    return None
                v = float(mol_conc) / t
            except (TypeError, ValueError):
                return None
            return min(1.0, max(0.0, v))

        CT_TIC = kwargs.get("CT_TIC")
        if CT_TIC is not None:
            co2_conc = out.get("CO2aq") if "CO2aq" in out else out.get("CO2")
            a = _alpha(co2_conc, CT_TIC)
            if a is not None:
                alphas["CO2aq"] = a

        CT_NH_T = kwargs.get("CT_NH_T")
        if CT_NH_T is not None:
            a = _alpha(out.get("NH3"), CT_NH_T)
            if a is not None:
                alphas["NH3"] = a

        for name, CT in (kwargs.get("acid_totals") or {}).items():
            mol_key = f"{name}_HA"
            if mol_key in out:
                a = _alpha(out.get(mol_key), CT)
                if a is not None:
                    alphas[mol_key] = a

        # ------------------------------------------------------------------
        # Build EquilibriumResult. species_mol_L is restricted to the fixed
        # canonical tuple this engine has always written back (state-
        # unification C3) — generic polyprotic-ladder keys (e.g. "{name}_HA")
        # were never part of that writeback and land in `extra` instead, not
        # species_mol_L, so apply_to_phases() doesn't change writeback
        # behavior.
        # ------------------------------------------------------------------
        species_mol_L: Dict[str, float] = {}
        extra: Dict[str, Any] = {}
        for key, val in out.items():
            if key in _META_KEYS:
                continue
            if key in _CANONICAL_WRITEBACK_SPECIES and val is not None:
                try:
                    species_mol_L[key] = float(val)
                except (TypeError, ValueError):
                    pass
            else:
                extra[key] = val

        pH_val = out.get("pH")
        result = EquilibriumResult(
            pH=float(pH_val) if pH_val is not None else float("nan"),
            pH_conc=out.get("pH_conc"),
            logH=out.get("logH"),
            aH=out.get("aH"),
            gamma_H=out.get("gamma_H"),
            gamma_OH=out.get("gamma_OH"),
            ionic_strength=out.get("IonicStrength"),
            charge_residual=out.get("charge_residual"),
            species_mol_L=species_mol_L,
            alphas=alphas,
            extra=extra,
        )
        return result

    # ------------------------------------------------------------------
    # Convenience helpers
    # ------------------------------------------------------------------
    def get_CO2aq_from_totals(self, **kwargs) -> float:
        """Return dissolved molecular CO2(aq) (mol/L) for the given aqueous totals."""
        result = self.solve(**kwargs)
        # "CO2aq" is not in the canonical writeback tuple (only "CO2" is),
        # so it may land in `extra` rather than `species_mol_L` depending on
        # which solve path produced it — check both.
        for key in ("CO2aq", "CO2"):
            val = result.species_mol_L.get(key)
            if val is None:
                val = result.extra.get(key)
            if val is not None:
                return float(val or 0.0)
        return 0.0

    @property
    def logH_warmstart(self) -> float | None:
        """Current warm-start log[H+] hint, or None if not yet populated."""
        return self._logH_last

    @property
    def I_warmstart(self) -> float | None:
        """Current warm-start ionic strength hint (mol/L), or None if not yet populated."""
        return self._I_last

    def reset_cache(self):
        """Reset warm-start caches (pH and ionic strength)."""
        self._logH_last = None
        self._I_last = None

    def reset_counters(self):
        """Reset solve-call counter(s)."""
        self.n_solve_calls = 0

    def algebraic_species(self) -> frozenset:
        """Return species IDs written to ``phase.n_mol`` by ``_refresh_derived``.

        These form the algebraic state ``z`` for any DAE formulation; the
        complement within ``phase.n_mol`` is the differential state ``y``.

        Returns ``frozenset()`` when no ``EquilibriumSet`` is bound (legacy
        kwargs path — all species are differential in that case).
        """
        eq_set = getattr(self, "_equilibrium_set", None)
        if eq_set is None:
            return frozenset()
        result = {"H+", "OH-"}
        for eq_def in eq_set:
            for sp in eq_def.species_refs:
                result.add(sp.id)
        return frozenset(result)


class ChemicalEquilibriumEngine(BisectionChemicalEquilibriumEngine):
    """Deprecated alias for :class:`BisectionChemicalEquilibriumEngine`.

    The unqualified name read as "the" default/canonical engine, when it is
    actually one of three peer implementations of
    ``ChemicalEquilibriumEngineProtocol`` (the others being
    ``NRChemicalEquilibriumEngine`` and ``PHREEQCChemicalEquilibriumEngine``)
    — specifically the original, simplest one (1-D bisection over the charge
    balance; single/independent acid-base ladders only). Kept for one phase
    as a subclass (not a factory function) so ``isinstance`` checks and
    ``from_reactions()`` (a ``cls``-based classmethod) both keep working
    unchanged; emits ``DeprecationWarning`` on construction.
    """

    def __init__(self, *args: Any, **kwargs: Any):
        warnings.warn(
            "ChemicalEquilibriumEngine is deprecated; use "
            "BisectionChemicalEquilibriumEngine instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(*args, **kwargs)


# ════════════════════════════════════════════════════════════════════════
#  Helpers for BisectionChemicalEquilibriumEngine.from_reactions
# ════════════════════════════════════════════════════════════════════════

_SOLVENT_IDS = ("H2O",)


def _classify_equilibrium(rxn) -> str:
    """Return ``"water"`` | ``"acid"`` | ``"cation_acid"`` for a
    declared equilibrium reaction.
    """
    products_ids = {
        e.species.id for e in rxn.stoichiometry if e.coefficient > 0
    }
    if "H+" in products_ids and "OH-" in products_ids:
        return "water"

    reactants = [
        e for e in rxn.stoichiometry
        if e.coefficient < 0 and e.species.id not in _SOLVENT_IDS
    ]
    if not reactants:
        raise ValueError(
            f"Equilibrium reaction {rxn.label!r} has no non-solvent "
            f"reactant (no entry with coefficient < 0)."
        )
    if len(reactants) > 1:
        raise ValueError(
            f"Equilibrium reaction {rxn.label!r} has multiple "
            "non-solvent reactants — only single-step monoprotic "
            "dissociations are supported by "
            "BisectionChemicalEquilibriumEngine.from_reactions in this phase."
        )
    acid = reactants[0].species
    if int(acid.charge) > 0:
        return "cation_acid"
    return "acid"


def _add_water(eq_set, rxn) -> None:
    """Translate a water-dissociation EquilibriumReaction into ``eq_set.set_water``."""
    pKw = -float(rxn.log_K)
    if rxn.dH_J_per_mol is not None and abs(rxn.dH_J_per_mol) > 1e-12:
        eq_set.set_water(
            pKw=pKw,
            correction="van_t_hoff",
            dH_J_per_mol=float(rxn.dH_J_per_mol),
            T_ref_K=float(rxn.T_ref_K),
        )
    else:
        eq_set.set_water(pKw=pKw, correction="none", T_ref_K=float(rxn.T_ref_K))


_H_PLUS_ID = "H+"


def _add_polyprotic_acid(eq_set, name: str, group: list) -> None:
    """Merge a group of single-step EquilibriumReactions that share a ``total_id``
    into one polyprotic EquilibriumDef.

    Reactions are sorted by ascending pKa (= -log_K descending).  The merged
    species_refs span the full protonation ladder from most-protonated acid to
    least-protonated base, assembled in order of deprotonation.
    """
    # Sort ascending by pKa (most negative log_K = largest pKa last)
    group = sorted(group, key=lambda t: float(t[1].log_K), reverse=True)

    pKas = tuple(-float(rxn.log_K) for _, rxn in group)
    category = group[0][0]
    T_ref_K = float(group[0][1].T_ref_K)
    total_key = _TOTAL_KEY_OVERRIDES.get(name, f"CT_{name}")

    # Build species_refs by walking the chain in pKa order
    seen_ids: set = set()
    species_refs: list = []
    for _, rxn in group:
        reactants = [
            e for e in rxn.stoichiometry
            if e.coefficient < 0 and e.species.id not in _SOLVENT_IDS
        ]
        acid_sp = reactants[0].species
        if acid_sp.id not in seen_ids:
            species_refs.append(acid_sp)
            seen_ids.add(acid_sp.id)
        base_sps = [
            e.species for e in rxn.stoichiometry
            if e.coefficient > 0
            and e.species.id not in _SOLVENT_IDS
            and e.species.id != _H_PLUS_ID
        ]
        base_sps.sort(key=lambda sp: -int(sp.charge))
        for sp in base_sps:
            if sp.id not in seen_ids:
                species_refs.append(sp)
                seen_ids.add(sp.id)

    # dH: use van_t_hoff only when ALL steps carry dH; otherwise none
    dH_vals = [rxn.dH_J_per_mol for _, rxn in group]
    if all(v is not None and abs(float(v)) > 1e-12 for v in dH_vals):
        eq_set.add(name,
                   category=category,
                   pKas=pKas,
                   T_ref_K=T_ref_K,
                   total_key=total_key,
                   species_refs=tuple(species_refs),
                   correction="van_t_hoff",
                   dH_J_per_mol=tuple(float(v) for v in dH_vals))
    else:
        eq_set.add(name,
                   category=category,
                   pKas=pKas,
                   T_ref_K=T_ref_K,
                   total_key=total_key,
                   species_refs=tuple(species_refs),
                   correction="none")


def _add_acid(eq_set, rxn, *, category: str) -> None:
    """Translate an acid EquilibriumReaction into ``eq_set.add(...)``.

    The ``EquilibriumDef`` name is the reaction's ``total_id`` when set,
    otherwise the acid form's id.
    """
    reactants = [
        e for e in rxn.stoichiometry
        if e.coefficient < 0 and e.species.id not in _SOLVENT_IDS
    ]
    acid_species = reactants[0].species
    acid_id = acid_species.id
    name = getattr(rxn, "total_id", None) or acid_id

    base_species = [
        e.species for e in rxn.stoichiometry
        if e.coefficient > 0
        and e.species.id not in _SOLVENT_IDS
        and e.species.id != _H_PLUS_ID
    ]
    base_species.sort(key=lambda sp: -int(sp.charge))
    species_refs = tuple([acid_species] + base_species)

    pKa = -float(rxn.log_K)
    total_key = _TOTAL_KEY_OVERRIDES.get(name, f"CT_{name}")

    add_kwargs = dict(
        category=category,
        pKas=(pKa,),
        T_ref_K=float(rxn.T_ref_K),
        total_key=total_key,
        species_refs=species_refs,
    )
    if rxn.dH_J_per_mol is not None and abs(rxn.dH_J_per_mol) > 1e-12:
        add_kwargs["correction"] = "van_t_hoff"
        add_kwargs["dH_J_per_mol"] = (float(rxn.dH_J_per_mol),)
    else:
        add_kwargs["correction"] = "none"

    eq_set.add(name, **add_kwargs)


# Mapping from canonical strong-ion species IDs to engine-internal CT_* keys.
_STRONG_ION_SPECIES_TO_KEY = {
    "K+":          "CT_K",
    "Na+":         "CT_Na",
    "Cl-":         "CT_Cl",
    "NO3-":        "CT_NO3",
    "Mg++":        "CT_Mg",
    "Ca++":        "CT_Ca",
    "Zn++":        "CT_Zn",
    "Mn++":        "CT_Mn",
    "Co++":        "CT_Co",
    "Mo7O24------": "CT_Mo7O24",
    "S_cat":       "CT_cation",
    "S_an":        "CT_anion",
}


def _populate_strong_ions_from_phases(phases, kwargs, *, liquid_key: str = "liquid") -> None:
    """Read strong-ion concentrations from ``phases[liquid_key].n_mol``."""
    if phases is None:
        return
    liq = phases.get(liquid_key) if hasattr(phases, "get") else None
    if liq is None:
        return
    V_L = float(getattr(liq, "V_L", 1.0))
    if V_L <= 0.0:
        return
    n_mol = getattr(liq, "n_mol", {})
    strong = dict(kwargs.get("strong_kwargs", {}) or {})
    for species_id, ct_key in _STRONG_ION_SPECIES_TO_KEY.items():
        if ct_key in strong:
            continue
        if species_id in n_mol:
            strong[ct_key] = float(n_mol[species_id]) / V_L
    if strong:
        kwargs["strong_kwargs"] = strong


def _populate_totals_from_phases(phases, kwargs, eq_set, *, liquid_key: str = "liquid") -> None:
    """Compute acid totals from ``phases[liquid_key].n_mol`` by summing
    each equilibrium's canonical-ladder species.
    """
    if eq_set is None or phases is None:
        return
    liq = phases.get(liquid_key) if hasattr(phases, "get") else None
    if liq is None:
        return
    V_L = float(getattr(liq, "V_L", 1.0))
    if V_L <= 0.0:
        return
    n_mol = getattr(liq, "n_mol", {})

    acid_totals = dict(kwargs.get("acid_totals", {}) or {})
    for eq_def in eq_set:
        key = eq_def.total_key
        if key in kwargs:
            continue

        if eq_def.species_refs:
            CT_mol = sum(float(n_mol.get(sp.id, 0.0))
                         for sp in eq_def.species_refs)
            if CT_mol <= 0.0:
                acid_id = eq_def.species_refs[0].id
                CT_mol = float(n_mol.get(acid_id, 0.0))
                if CT_mol <= 0.0:
                    CT_mol = float(n_mol.get(eq_def.name, 0.0))
        else:
            CT_mol = float(n_mol.get(eq_def.name, 0.0))
        CT = CT_mol / V_L

        if eq_def.category == "acid" and eq_def.name not in _TOTAL_KEY_OVERRIDES:
            if eq_def.name not in acid_totals:
                acid_totals[eq_def.name] = CT
        else:
            kwargs[key] = CT
    if acid_totals:
        kwargs["acid_totals"] = acid_totals
