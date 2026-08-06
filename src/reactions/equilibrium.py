# -*- coding: utf-8 -*-
"""Equilibrium reaction — algebraic constraint with an equilibrium constant.

An :class:`EquilibriumReaction` bundles a validated stoichiometric
template with an equilibrium constant (``log_K``), optional Van 't Hoff
temperature-dependence parameters, and an optional ``total_id`` mapping
to the species id in a phase's ``n_mol`` dict that carries the
*total* concentration this equilibrium operates on.

This is one of three independent reaction declaration classes —
:class:`~PyOMES.reactions.kinetic.KineticReaction`,
:class:`EquilibriumReaction`, and
:class:`~PyOMES.reactions.blackbox.BlackBoxReactionModel`. There is no
shared base class; shared validation logic lives as free functions in
``_shared.py``.

Equilibrium reactions have no rate; they are routed to
:class:`~PyOMES.chemical_equilibrium.engine.BisectionChemicalEquilibriumEngine` as algebraic
constraints (with ``log_K`` setting their constant), or to
:class:`~PyOMES.core.gas_liquid_link.KineticGasLiquidLink` as partition
declarations (cross-phase reactions). They do not implement
``compute_rates``.

Example
-------
>>> from PyOMES.reactions import EquilibriumReaction, StoichiometryEntry
>>> from PyOMES.chemistry import Species
>>> H2O = Species(id="H2O", atoms={"H": 2, "O": 1})
>>> H_plus = Species(id="H+", atoms={"H": 1}, charge=+1)
>>> OH_minus = Species(id="OH-", atoms={"O": 1, "H": 1}, charge=-1)
>>> rxn = EquilibriumReaction(
...     stoichiometry=[
...         StoichiometryEntry(species=H2O,      phase="liquid", coefficient=-1.0),
...         StoichiometryEntry(species=H_plus,   phase="liquid", coefficient=+1.0),
...         StoichiometryEntry(species=OH_minus, phase="liquid", coefficient=+1.0),
...     ],
...     log_K=-14.0,
...     label="water",
... )
"""

from __future__ import annotations

import math
from typing import (
    Dict, List, Literal, Optional, Protocol, Sequence, Union, runtime_checkable,
)

from ..chemistry.species import Species
from ..units import R_J_PER_MOL_K as _R_J_MOL_K
from .stoichiometry import StoichiometryEntry, _parse_stoichiometry
from ._shared import (
    coerce_and_validate,
    fmt_stoichiometry_string,
    is_cross_phase_from_entries,
    phases_from_entries,
    show_balance_from_entries,
    species_ids_from_entries,
)

_LOG10_E = 1.0 / math.log(10.0)


@runtime_checkable
class EquilibriumConstraint(Protocol):
    """Structural contract for the mass-action equilibrium family.

    Any type exposing these four attributes — a stoichiometry, a
    reference-temperature ``log_K``, an optional Van 't Hoff
    ``dH_J_per_mol``, and the reference temperature they were measured
    at — can be routed through :func:`vant_hoff_log_K` and the
    equilibrium-classification/tableau-building machinery in
    ``PyOMES.chemical_equilibrium``, regardless of whether it also satisfies
    :class:`~PyOMES.chemistry.partition.PartitionModel` (as
    :class:`~PyOMES.chemistry.partition.HenryEquilibrium`,
    :class:`~PyOMES.chemistry.partition.KspEquilibrium`, and
    :class:`~PyOMES.chemistry.partition.RaoultEquilibrium` all do).

    ``log_K``/``dH_J_per_mol`` are plain attributes carrying the
    *reference* mass-action constant (at ``T_ref_K``) — not methods.
    Composition-dependent corrections (``γ_i``, ``φ_i``) are not part
    of this protocol; temperature correction is handled separately by
    :func:`vant_hoff_log_K`.
    """

    stoichiometry: Sequence[StoichiometryEntry]
    log_K: float
    dH_J_per_mol: Optional[float]
    T_ref_K: float


def vant_hoff_log_K(constraint: EquilibriumConstraint, T_K: float) -> float:
    """Van 't Hoff temperature-corrected log10(K) for any EquilibriumConstraint.

    Returns ``constraint.log_K`` unchanged when ``dH_J_per_mol`` is
    ``None``/~0 or when ``T_K`` is at the reference temperature.
    """
    log_K_ref = float(constraint.log_K)
    dH = constraint.dH_J_per_mol
    T_ref_K = float(constraint.T_ref_K)
    if dH is None or abs(dH) < 1e-30:
        return log_K_ref
    if abs(T_K - T_ref_K) < 1e-10:
        return log_K_ref
    # ln K(T) = ln K(T_ref) − (ΔH/R) (1/T − 1/T_ref)
    delta_ln_K = -(float(dH) / _R_J_MOL_K) * (1.0 / float(T_K) - 1.0 / T_ref_K)
    return log_K_ref + delta_ln_K * _LOG10_E


def classify_equilibrium_constraint(
    item: EquilibriumConstraint,
) -> Literal["acid_base", "gas_liquid", "solid_liquid"]:
    """Classify an :class:`EquilibriumConstraint` from its stoichiometry's phase tags.

    Not an ``isinstance(item, EquilibriumReaction)`` check —
    :class:`~PyOMES.chemistry.partition.HenryEquilibrium`,
    :class:`~PyOMES.chemistry.partition.KspEquilibrium`, and
    :class:`~PyOMES.chemistry.partition.RaoultEquilibrium` are siblings,
    not subclasses, of :class:`EquilibriumReaction`. Any
    ``EquilibriumConstraint``-conforming item is classified purely from
    the distinct phases present across its ``stoichiometry``:

    - any entry with ``phase == "solid"`` → ``"solid_liquid"``
      (takes priority — a solid entry always means a precipitation
      equilibrium, regardless of what else is present).
    - more than one distinct phase (and no solid) → ``"gas_liquid"``.
    - a single phase → ``"acid_base"``.

    Raises
    ------
    ValueError
        If ``item.stoichiometry`` is empty (e.g. a
        ``PartitionModel``-only ``HenryEquilibrium``/``RaoultEquilibrium``
        constructed without ``gas_species``/``liquid_species`` — such an
        instance has no reaction row to classify).
    """
    entries = item.stoichiometry
    if not entries:
        raise ValueError(
            f"classify_equilibrium_constraint: {item!r} has an empty "
            "stoichiometry — nothing to classify (e.g. a PartitionModel-"
            "only HenryEquilibrium/RaoultEquilibrium with gas_species/"
            "liquid_species unset)."
        )
    phases = phases_from_entries(entries)
    if "solid" in phases:
        return "solid_liquid"
    if len(phases) > 1:
        return "gas_liquid"
    return "acid_base"


class EquilibriumReaction:
    """A single reaction satisfied as an algebraic equilibrium constraint.

    Parameters
    ----------
    stoichiometry : list of StoichiometryEntry or str
        All participants, either as a list of :class:`StoichiometryEntry`
        objects or as a human-readable string such as
        ``"CO2,aq + H2O,aq <-> HCO3-,aq + H+,aq"``.
        String format: terms separated by `` + ``, arrow ``<->`` separating
        reactants from products, each term optionally prefixed by a
        coefficient and suffixed with ``,phase`` (aq/l → liquid, g → gas,
        s → solid). Coefficients are per unit extent of reaction
        (positive = produced, negative = consumed). Sign convention for
        ``log_K``: products / reactants — for ``HA ⇌ A⁻ + H⁺`` with
        ``Ka = 10^-pKa``, ``log_K = -pKa``.
    species : dict[str, Species], optional
        Caller-supplied species for locally declared IDs not in
        ``common_species``. Only used when *stoichiometry* is a string.
    log_K : float, optional
        Base-10 logarithm of the equilibrium constant. Required for
        single-phase reactions; **optional** for cross-phase
        reactions (e.g. gas-liquid partition declarations like
        ``CO2(gas) ⇌ CO2aq(liquid)``) because the partition constant
        lives on the consuming link (Henry's law on
        :class:`~PyOMES.core.gas_liquid_link.KineticGasLiquidLink`),
        not on the reaction.
    dH_J_per_mol : float, optional
        Van 't Hoff reaction enthalpy (J/mol) at ``T_ref_K``. When
        supplied, :meth:`BisectionChemicalEquilibriumEngine.from_reactions` applies
        Van 't Hoff temperature correction to ``log_K`` for each
        speciation solve.
    T_ref_K : float
        Reference temperature for ``log_K`` and ``dH_J_per_mol``
        (default ``298.15``). Only used when ``dH_J_per_mol`` is
        supplied.
    total_id : str, optional
        Species id in the phase ``n_mol`` dict that carries the
        *total* concentration this equilibrium operates on (e.g.
        BSM2 tracks total ammoniacal N as ``n_mol["NH3"]`` even
        though the reaction is written ``NH4⁺ ⇌ NH3 + H⁺``). When
        unset, :meth:`BisectionChemicalEquilibriumEngine.from_reactions` defaults this
        to the acid form's id.
    balance_elements : sequence of str or None
        Elements to validate.  ``None`` (default) infers them from the
        atoms present across all participants — every element that
        appears anywhere in the stoichiometry is checked.  Pass an
        explicit sequence to restrict the check to a subset.
    balance_atol : float
        Absolute tolerance for each element (default ``1e-10``).
    label : str
        Human-readable label (optional, for diagnostics/logging).

    Raises
    ------
    StoichiometryError
        If the stoichiometry does not close for any element in
        *balance_elements*.
    ValueError
        If ``log_K`` is missing on a single-phase equilibrium reaction.
    """

    def __init__(
        self,
        stoichiometry: Union[Sequence[StoichiometryEntry], str],
        *,
        species: Optional[Dict[str, Species]] = None,
        log_K: Optional[float] = None,
        dH_J_per_mol: Optional[float] = None,
        T_ref_K: float = 298.15,
        total_id: Optional[str] = None,
        balance_elements: Optional[Sequence[str]] = None,
        balance_atol: float = 1e-10,
        label: str = "",
    ):
        if isinstance(stoichiometry, str):
            stoichiometry = _parse_stoichiometry(
                stoichiometry, species, reaction_type="equilibrium"
            )
        entries = coerce_and_validate(
            stoichiometry, balance_elements, balance_atol
        )
        if log_K is None and not is_cross_phase_from_entries(entries):
            raise ValueError(
                "EquilibriumReaction requires log_K for single-phase "
                "reactions. Cross-phase equilibrium reactions "
                "(partition declarations consumed by gas-liquid "
                "links) may omit log_K."
            )

        self.stoichiometry: List[StoichiometryEntry] = entries
        self.log_K: Optional[float] = (
            float(log_K) if log_K is not None else None
        )
        self.dH_J_per_mol: Optional[float] = (
            float(dH_J_per_mol) if dH_J_per_mol is not None else None
        )
        self.T_ref_K: float = float(T_ref_K)
        self.total_id: Optional[str] = (
            str(total_id) if total_id is not None else None
        )
        from ._shared import _infer_elements
        self.balance_elements = (
            tuple(_infer_elements(entries))
            if balance_elements is None
            else tuple(balance_elements)
        )
        self.label = str(label)

    @property
    def species_ids(self) -> List[str]:
        """Return the sorted, unique species IDs in this reaction."""
        return species_ids_from_entries(self.stoichiometry)

    @property
    def phases(self) -> List[str]:
        """Return the sorted, unique phase keys in this reaction."""
        return phases_from_entries(self.stoichiometry)

    @property
    def is_cross_phase(self) -> bool:
        """Whether the stoichiometry spans more than one phase.

        Cross-phase equilibrium reactions are *partition
        declarations* — they tell a
        :class:`~PyOMES.core.gas_liquid_link.KineticGasLiquidLink`
        which liquid molecular form corresponds to a given gas
        species. They are not thermodynamic constraints solved by
        :class:`~PyOMES.chemical_equilibrium.engine.BisectionChemicalEquilibriumEngine`; the
        engine's
        :meth:`~PyOMES.chemical_equilibrium.engine.BisectionChemicalEquilibriumEngine.from_reactions`
        skips them silently.
        """
        return is_cross_phase_from_entries(self.stoichiometry)

    def plot_vant_hoff(
        self,
        *,
        T_range_K: tuple = (273.15, 373.15),
        n_points: int = 200,
        ax=None,
        show_pka: bool = True,
        title: Optional[str] = None,
    ) -> tuple:
        """Plot the Van 't Hoff temperature dependence of this reaction.

        Delegates to :func:`PyOMES.reactions.plots.plot_vant_hoff`.
        Matplotlib is imported lazily — the reactions package does not
        require a display environment.

        Parameters
        ----------
        T_range_K : (float, float)
            Temperature range in Kelvin.  Default 0 – 100 °C.
        n_points : int
            Number of curve points.  Default 200.
        ax : matplotlib Axes, optional
            Axes to draw on; a new figure is created when ``None``.
        show_pka : bool
            Add a secondary pKa axis.  Default ``True``.
        title : str, optional
            Axes title; defaults to ``self.label``.

        Returns
        -------
        fig, ax
            The matplotlib Figure and primary Axes.
        """
        from .plots import plot_vant_hoff
        return plot_vant_hoff(
            self,
            T_range_K=T_range_K,
            n_points=n_points,
            ax=ax,
            show_pka=show_pka,
            title=title,
        )

    def show_stoichiometry(self) -> None:
        """Print the reaction stoichiometry in human-readable string format."""
        print(f"{self.label or type(self).__name__}:")
        print(f"  {fmt_stoichiometry_string(self.stoichiometry, arrow='<->')}")

    def show_balance(self, *, elements: Optional[Sequence[str]] = None) -> None:
        """Print per-element residuals for all elements present in the stoichiometry.

        Parameters
        ----------
        elements : sequence of str, optional
            Elements to check.  Defaults to every element found across the
            atoms of all participants — no need to specify them manually.
        """
        show_balance_from_entries(self.stoichiometry, self.label, elements)

    def __repr__(self) -> str:
        species = ", ".join(
            f"{e.coefficient:+.3g} {e.species.id}" for e in self.stoichiometry
        )
        lbl = f" [{self.label}]" if self.label else ""
        return f"EquilibriumReaction({species}{lbl})"
