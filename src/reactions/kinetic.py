# -*- coding: utf-8 -*-
"""Kinetic reaction — single reaction with a callable rate law.

A :class:`KineticReaction` bundles a validated stoichiometric template
with a callable rate law. The stoichiometry is validated at construction
time for elemental closure, so any ``compute_rates`` call is guaranteed
to produce mass-consistent source terms.

This is one of three independent reaction declaration classes —
:class:`KineticReaction`, :class:`EquilibriumReaction`, and
:class:`~PyOMES.reactions.blackbox.BlackBoxReactionModel`. There is no
shared base class; shared validation logic lives as free functions in
``_shared.py``.

Example
-------
>>> from PyOMES.reactions import KineticReaction, StoichiometryEntry
>>> from PyOMES.chemistry import Species
>>> A = Species(id="A", atoms={"C": 1})
>>> B = Species(id="B", atoms={"C": 1})
>>> rxn = KineticReaction(
...     stoichiometry=[
...         StoichiometryEntry(species=A, phase="liquid", coefficient=-1.0),
...         StoichiometryEntry(species=B, phase="liquid", coefficient=+1.0),
...     ],
...     rate_fn=lambda env: 0.1 * env.S("A") * env.V_L,
... )
>>> result = rxn.compute_rates(env)
"""

from __future__ import annotations

from typing import Callable, Dict, List, Optional, Sequence, Union

from ..chemistry.species import Species
from .stoichiometry import StoichiometryEntry, _parse_stoichiometry
from .environment import ReactionEnvironment
from ._shared import (
    coerce_and_validate,
    fmt_stoichiometry_string,
    is_cross_phase_from_entries,
    phases_from_entries,
    show_balance_from_entries,
    species_ids_from_entries,
)


class KineticReaction:
    """A single reaction integrated via its rate law.

    Parameters
    ----------
    stoichiometry : list of StoichiometryEntry or str
        All participants, either as a list of :class:`StoichiometryEntry`
        objects or as a human-readable string such as
        ``"CO2,g -> CO2,aq"``.
        String format: terms separated by `` + ``, arrow ``->`` separating
        reactants from products, each term optionally prefixed by a
        coefficient and suffixed with ``,phase`` (aq/l → liquid, g → gas,
        s → solid). Coefficients are per unit extent of reaction
        (positive = produced, negative = consumed).
    species : dict[str, Species], optional
        Caller-supplied species for locally declared IDs not in
        ``common_species``. Only used when *stoichiometry* is a string.
    rate_fn : callable
        ``rate_fn(env: ReactionEnvironment) -> float`` returning the
        **extensive** rate of reaction in mol/h.
    balance_elements : sequence of str or None
        Elements to validate.  ``None`` (default) infers them from the
        atoms present across all participants.  Pass an explicit sequence
        to restrict the check to a subset.
    balance_atol : float
        Absolute tolerance for each element (default ``1e-10``).
    label : str
        Human-readable label (optional, for diagnostics/logging).

    Raises
    ------
    StoichiometryError
        If the stoichiometry does not close for any element in
        *balance_elements*.
    """

    def __init__(
        self,
        stoichiometry: Union[Sequence[StoichiometryEntry], str],
        rate_fn: Callable[[ReactionEnvironment], float],
        *,
        species: Optional[Dict[str, Species]] = None,
        balance_elements: Optional[Sequence[str]] = None,
        balance_atol: float = 1e-10,
        label: str = "",
    ):
        if isinstance(stoichiometry, str):
            stoichiometry = _parse_stoichiometry(
                stoichiometry, species, reaction_type="kinetic"
            )
        self.stoichiometry: List[StoichiometryEntry] = coerce_and_validate(
            stoichiometry, balance_elements, balance_atol
        )
        self.rate_fn = rate_fn
        from ._shared import _infer_elements
        self.balance_elements = (
            tuple(_infer_elements(self.stoichiometry))
            if balance_elements is None
            else tuple(balance_elements)
        )
        self.label = str(label)

    def compute_rates(
        self, env: ReactionEnvironment
    ) -> Dict[str, Dict[str, float]]:
        """Evaluate the rate law and return phase-keyed source terms.

        Parameters
        ----------
        env : ReactionEnvironment
            Current local conditions.

        Returns
        -------
        dict of dict
            ``{phase_key: {species_id: mol_per_h}}``.
        """
        r = float(self.rate_fn(env))

        sources: Dict[str, Dict[str, float]] = {}
        for entry in self.stoichiometry:
            phase = entry.phase
            sp_id = entry.species.id
            if phase not in sources:
                sources[phase] = {}
            sources[phase][sp_id] = (
                sources[phase].get(sp_id, 0.0)
                + entry.coefficient * r
            )
        return sources

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

        Kinetic reactions are typically single-phase, but the flag is
        exposed for symmetry with
        :class:`~PyOMES.reactions.equilibrium.EquilibriumReaction`.
        """
        return is_cross_phase_from_entries(self.stoichiometry)

    def show_stoichiometry(self) -> None:
        """Print the reaction stoichiometry in human-readable string format."""
        print(f"{self.label or type(self).__name__}:")
        print(f"  {fmt_stoichiometry_string(self.stoichiometry, arrow='->')}")

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
        return f"KineticReaction({species}{lbl})"
