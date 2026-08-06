# -*- coding: utf-8 -*-
"""Adapter for opaque external metabolic simulators.

:class:`BlackBoxReactionModel` wraps an external model (FBA solver,
genome-scale kinetic model, proprietary simulator, etc.) behind the
:class:`ReactionModel` protocol.  Because the internal stoichiometry
is not inspectable, elemental balance is checked at **runtime** rather
than at construction.

Example
-------
>>> model = BlackBoxReactionModel(
...     external_model=my_fba_solver,
...     flux_mapping={
...         "EX_ac_e": FluxEntry("AceticAcid", "liquid", {"C":2,"H":4,"O":2}, sign=+1),
...         "EX_o2_e": FluxEntry("O2",         "liquid", {"O":2},             sign=-1),
...         ...
...     },
...     on_imbalance="warn",
... )
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Optional, Sequence

from .environment import ReactionEnvironment


class MassBalanceWarning(UserWarning):
    """Issued when a black-box model returns mass-imbalanced fluxes."""
    pass


class MassBalanceError(ValueError):
    """Raised when a black-box model returns mass-imbalanced fluxes
    and ``on_imbalance="raise"``."""
    pass


@dataclass(frozen=True)
class FluxEntry:
    """Mapping from an external model flux ID to a CV species.

    Parameters
    ----------
    species_id : str
        CV / chemistry registry species ID.
    phase : str
        Target phase key (``"liquid"``, ``"gas"``).
    atoms : dict
        Elemental composition (for runtime balance checking).
    sign : float
        Sign convention adapter.  ``+1`` if the external model's
        positive output means production, ``-1`` if it means
        consumption.  Default ``+1``.
    """

    species_id: str
    phase: str
    atoms: Dict[str, float] = field(default_factory=dict)
    sign: float = 1.0


class BlackBoxReactionModel:
    """Adapter for an opaque external metabolic simulator.

    The external model is called each timestep with conditions derived
    from the CV's state.  Its outputs are mapped to phase-keyed source
    terms via :class:`FluxEntry` objects.

    Mass balance is checked at runtime because the stoichiometry is
    not available at construction time.

    Parameters
    ----------
    external_model : object
        Any object with a ``solve(inputs: dict) -> dict`` method.
        The method receives a dict of conditions (built by
        ``_build_model_inputs``) and returns a dict of
        ``{external_flux_id: value}`` where value is in mol/h.
    flux_mapping : dict
        ``{external_flux_id: FluxEntry}`` mapping the external
        model's output IDs to CV species, phases, and atoms.
    balance_elements : sequence of str
        Elements to check at runtime (default ``("C", "H", "O")``).
    balance_atol : float
        Tolerance for elemental residual (default ``1e-8``).
    on_imbalance : str
        Policy when balance check fails: ``"warn"`` (default),
        ``"raise"``, ``"log"``, or ``"ignore"``.
    input_builder : callable or None
        Optional ``f(env, flux_mapping) -> dict`` that translates a
        :class:`ReactionEnvironment` into the external model's input
        format.  If ``None``, a default builder is used.
    """

    def __init__(
        self,
        external_model: Any,
        flux_mapping: Dict[str, FluxEntry],
        *,
        balance_elements: Sequence[str] = ("C", "H", "O"),
        balance_atol: float = 1e-8,
        on_imbalance: str = "warn",
        input_builder: Optional[Callable] = None,
    ):
        self.external_model = external_model
        self.flux_mapping = dict(flux_mapping)
        self.balance_elements = tuple(balance_elements)
        self.balance_atol = float(balance_atol)
        self.on_imbalance = str(on_imbalance).lower().strip()
        self._input_builder = input_builder
        self._imbalance_history: list = []

    @property
    def imbalance_history(self) -> list:
        """Access the list of recorded balance-check failures."""
        return list(self._imbalance_history)

    def compute_rates(
        self, env: ReactionEnvironment
    ) -> Dict[str, Dict[str, float]]:
        """Call the external model and return phase-keyed source terms.

        Parameters
        ----------
        env : ReactionEnvironment
            Current local conditions.

        Returns
        -------
        dict of dict
            ``{phase_key: {species_id: mol_per_h}}``.
        """
        # 1. Build inputs for the external model
        inputs = self._build_model_inputs(env)

        # 2. Call the black box
        raw_fluxes = self.external_model.solve(inputs)

        # 3. Map outputs to phase-keyed source terms
        sources = self._map_fluxes(raw_fluxes)

        # 4. Runtime balance check
        self._check_balance(sources)

        return sources

    def _build_model_inputs(self, env: ReactionEnvironment) -> dict:
        """Translate a ReactionEnvironment into external model inputs.

        Override this method (or provide ``input_builder`` at
        construction) to match a specific external model API.
        """
        if self._input_builder is not None:
            return self._input_builder(env, self.flux_mapping)

        # Default: provide concentrations keyed by CV species ID,
        # plus environmental scalars.
        inputs: Dict[str, Any] = {}
        for _ext_id, entry in self.flux_mapping.items():
            inputs[entry.species_id] = env.concentrations.get(
                entry.species_id, 0.0
            )
        inputs["_temperature_K"] = env.T_K
        inputs["_pH"] = env.pH
        inputs["_volume_L"] = env.V_L
        return inputs

    def _map_fluxes(
        self, raw_fluxes: Dict[str, float]
    ) -> Dict[str, Dict[str, float]]:
        """Map raw external model outputs to phase-keyed source terms."""
        sources: Dict[str, Dict[str, float]] = {}
        for ext_id, value in raw_fluxes.items():
            entry = self.flux_mapping.get(ext_id)
            if entry is None:
                continue  # unmapped output — skip
            rate = float(value) * float(entry.sign)
            phase = entry.phase
            if phase not in sources:
                sources[phase] = {}
            sources[phase][entry.species_id] = (
                sources[phase].get(entry.species_id, 0.0) + rate
            )
        return sources

    def _check_balance(
        self, sources: Dict[str, Dict[str, float]]
    ) -> None:
        """Check elemental balance of the returned source terms."""
        for elem in self.balance_elements:
            total = 0.0
            for phase, species_rates in sources.items():
                for species_id, rate in species_rates.items():
                    # Find atoms for this species from the flux mapping
                    atoms = self._atoms_for(species_id)
                    total += rate * atoms.get(elem, 0.0)

            if abs(total) > self.balance_atol:
                record = {
                    "element": elem,
                    "residual": total,
                    "sources": {
                        ph: dict(sr) for ph, sr in sources.items()
                    },
                }
                self._imbalance_history.append(record)

                msg = (
                    f"Black-box model mass imbalance: element '{elem}' "
                    f"residual = {total:+.6e}"
                )
                if self.on_imbalance == "raise":
                    raise MassBalanceError(msg)
                elif self.on_imbalance == "warn":
                    warnings.warn(msg, MassBalanceWarning, stacklevel=3)
                # "log" and "ignore" — recorded in history, no exception/warning

    def _atoms_for(self, species_id: str) -> Dict[str, float]:
        """Look up elemental composition for a species from the flux mapping."""
        for entry in self.flux_mapping.values():
            if entry.species_id == species_id:
                return dict(entry.atoms)
        return {}
