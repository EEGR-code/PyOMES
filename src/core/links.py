# -*- coding: utf-8 -*-
"""Inter-CV material transport links.

A :class:`CVLink` computes molar fluxes between two ControlVolumes.
The :class:`MultiCVSystem` applies all links before advancing each CV,
so that inter-zone transport is resolved between timesteps.

Concrete implementations:

- :class:`AdvectiveLink` — bulk liquid circulation (carries all species
  in proportion to their concentration in the source phase).
- :class:`DiffusiveLink` — concentration-driven transfer of selected
  species (e.g. gas mixing between headspace zones, or eddy-diffusion
  transport of dissolved species between liquid zones).

Example
-------
>>> link = AdvectiveLink(
...     source_cv_key="sparger", source_phase_key="liquid",
...     sink_cv_key="bulk",     sink_phase_key="liquid",
...     Q_L_per_h=500.0,
... )
>>> flow = link.compute_flow(cvs, dt_h=0.01)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Protocol, Sequence, runtime_checkable


# ════════════════════════════════════════════════════════════════════════
#  Protocol
# ════════════════════════════════════════════════════════════════════════

@runtime_checkable
class CVLink(Protocol):
    """Contract for an inter-CV material transport link.

    A link connects a source phase in one CV to a sink phase in another
    CV (or the same CV, for internal recirculation).

    ``compute_flow`` returns species molar fluxes (mol/h) **from source
    to sink**.  Positive values mean material leaves the source and
    enters the sink.  The :class:`MultiCVSystem` applies the fluxes
    symmetrically: ``source.apply_flux(-flow, dt)`` and
    ``sink.apply_flux(+flow, dt)``.
    """

    @property
    def source_cv_key(self) -> str:
        """Key of the source ControlVolume in the MultiCVSystem."""
        ...

    @property
    def source_phase_key(self) -> str:
        """Phase key within the source CV (e.g. ``"liquid"``)."""
        ...

    @property
    def sink_cv_key(self) -> str:
        """Key of the sink ControlVolume."""
        ...

    @property
    def sink_phase_key(self) -> str:
        """Phase key within the sink CV."""
        ...

    @property
    def label(self) -> str:
        """Human-readable label for this link."""
        ...

    def compute_flow(
        self,
        cvs: Dict[str, Any],
        dt_h: float,
    ) -> Dict[str, float]:
        """Compute species molar fluxes from source to sink.

        Parameters
        ----------
        cvs : dict
            ``{cv_key: ControlVolume}`` — the full system, so the link
            can read concentrations from both source and sink.
        dt_h : float
            Timestep duration (hours).

        Returns
        -------
        dict
            ``{species_id: flux_mol_per_h}`` — positive = source→sink.
        """
        ...


# ════════════════════════════════════════════════════════════════════════
#  Link result (for diagnostics)
# ════════════════════════════════════════════════════════════════════════

@dataclass
class LinkFlowRecord:
    """Record of a single link's flow computation.

    Attributes
    ----------
    link_label : str
        Human-readable label of the link.
    source : str
        ``"source_cv_key.source_phase_key"``.
    sink : str
        ``"sink_cv_key.sink_phase_key"``.
    flow_mol_per_h : dict
        Species molar fluxes applied.
    total_mol_transferred : float
        Sum of ``|flux| * dt`` across all species.
    """

    link_label: str = ""
    source: str = ""
    sink: str = ""
    flow_mol_per_h: Dict[str, float] = field(default_factory=dict)
    total_mol_transferred: float = 0.0


# ════════════════════════════════════════════════════════════════════════
#  AdvectiveLink — bulk liquid (or gas) circulation
# ════════════════════════════════════════════════════════════════════════

@dataclass
class AdvectiveLink:
    """Bulk advective transport between two phases in different CVs.

    All species in the source phase are carried at their local
    concentration, proportional to the volumetric flow rate.  This
    models liquid circulation between zones, or gas recirculation
    between headspace compartments.

    The molar flux for each species is::

        flux_i (mol/h) = C_i (mol/L) × Q (L/h)

    where ``C_i`` is the concentration in the source phase and ``Q`` is
    the volumetric flow rate.

    Parameters
    ----------
    source_cv_key : str
        Key of the source ControlVolume.
    source_phase_key : str
        Phase key in the source CV (typically ``"liquid"``).
    sink_cv_key : str
        Key of the sink ControlVolume.
    sink_phase_key : str
        Phase key in the sink CV.
    Q_L_per_h : float
        Volumetric flow rate (L/h).  Must be ≥ 0.
    species_filter : list or None
        If provided, only these species are transferred.  If ``None``
        (default), all species in the source phase are transferred.
    label : str
        Human-readable label.
    """

    _source_cv_key: str
    _source_phase_key: str
    _sink_cv_key: str
    _sink_phase_key: str
    Q_L_per_h: float
    species_filter: Optional[Sequence[str]] = None
    _label: str = ""

    @property
    def source_cv_key(self) -> str:
        return self._source_cv_key

    @property
    def source_phase_key(self) -> str:
        return self._source_phase_key

    @property
    def sink_cv_key(self) -> str:
        return self._sink_cv_key

    @property
    def sink_phase_key(self) -> str:
        return self._sink_phase_key

    @property
    def target_phase_key(self) -> str:
        """Alias for ``sink_phase_key`` (Phase D matrix-assembly convention)."""
        return self._sink_phase_key

    @property
    def label(self) -> str:
        if self._label:
            return self._label
        return f"advective:{self._source_cv_key}.{self._source_phase_key}->{self._sink_cv_key}.{self._sink_phase_key}"

    def compute_flow(
        self,
        cvs: Dict[str, Any],
        dt_h: float,
    ) -> Dict[str, float]:
        Q = float(self.Q_L_per_h)
        if Q <= 0.0:
            return {}

        src_cv = cvs.get(self._source_cv_key)
        if src_cv is None:
            return {}
        src_phase = src_cv.phases.get(self._source_phase_key)
        if src_phase is None:
            return {}

        V_src = float(src_phase.V_L)
        n_mol = src_phase.n_mol

        flow: Dict[str, float] = {}
        species = self.species_filter if self.species_filter is not None else list(n_mol.keys())

        for sp in species:
            n = float(n_mol.get(sp, 0.0))
            if n <= 0.0:
                continue
            C = n / V_src  # mol/L
            flux = C * Q   # mol/h

            # Safety: don't transfer more than what's available in one timestep
            max_flux = n / max(float(dt_h), 1e-30)
            flux = min(flux, max_flux)

            if flux > 0.0:
                flow[sp] = flux

        return flow

    def __repr__(self):
        return (f"AdvectiveLink({self._source_cv_key}.{self._source_phase_key}"
                f" -> {self._sink_cv_key}.{self._sink_phase_key}, "
                f"Q={self.Q_L_per_h:.4g} L/h)")


# ════════════════════════════════════════════════════════════════════════
#  DiffusiveLink — concentration-driven species transfer
# ════════════════════════════════════════════════════════════════════════

@dataclass
class DiffusiveLink:
    """Concentration-driven transfer of specific species between CVs.

    For each species, the flux is driven by the concentration difference
    between the source and sink phases::

        flux_i (mol/h) = kLa_i (1/h) × V_eff (L) × (C_source_i - C_sink_i)

    where ``kLa_i`` is the species-specific mass transfer coefficient,
    ``V_eff`` is an effective volume for the transfer (defaults to the
    harmonic mean of source and sink volumes), and ``C`` are molar
    concentrations.

    The sign convention is: positive flux = from higher concentration
    to lower.  If sink has higher concentration than source, the flux
    reverses direction (material flows from sink to source through the
    ``-flow`` / ``+flow`` application in the MultiCVSystem).

    Parameters
    ----------
    source_cv_key : str
        Key of one ControlVolume (conventionally the "source").
    source_phase_key : str
        Phase key in the source CV.
    sink_cv_key : str
        Key of the other ControlVolume.
    sink_phase_key : str
        Phase key in the sink CV.
    kLa : dict
        ``{species_id: kLa_value_per_h}`` for each species to transfer.
    V_eff_L : float or None
        Effective volume (L).  If ``None``, uses the harmonic mean of
        the source and sink phase volumes.
    label : str
        Human-readable label.
    """

    _source_cv_key: str
    _source_phase_key: str
    _sink_cv_key: str
    _sink_phase_key: str
    kLa: Dict[str, float] = field(default_factory=dict)
    V_eff_L: Optional[float] = None
    _label: str = ""

    @property
    def source_cv_key(self) -> str:
        return self._source_cv_key

    @property
    def source_phase_key(self) -> str:
        return self._source_phase_key

    @property
    def sink_cv_key(self) -> str:
        return self._sink_cv_key

    @property
    def sink_phase_key(self) -> str:
        return self._sink_phase_key

    @property
    def target_phase_key(self) -> str:
        """Alias for ``sink_phase_key`` (Phase D matrix-assembly convention)."""
        return self._sink_phase_key

    @property
    def label(self) -> str:
        if self._label:
            return self._label
        return f"diffusive:{self._source_cv_key}.{self._source_phase_key}<->{self._sink_cv_key}.{self._sink_phase_key}"

    def compute_flow(
        self,
        cvs: Dict[str, Any],
        dt_h: float,
    ) -> Dict[str, float]:
        src_cv = cvs.get(self._source_cv_key)
        snk_cv = cvs.get(self._sink_cv_key)
        if src_cv is None or snk_cv is None:
            return {}

        src_phase = src_cv.phases.get(self._source_phase_key)
        snk_phase = snk_cv.phases.get(self._sink_phase_key)
        if src_phase is None or snk_phase is None:
            return {}

        V_src = float(src_phase.V_L)
        V_snk = float(snk_phase.V_L)

        # Effective volume for the transfer
        if self.V_eff_L is not None:
            V_eff = max(float(self.V_eff_L), 1e-30)
        else:
            # Harmonic mean
            V_eff = 2.0 * V_src * V_snk / (V_src + V_snk)

        flow: Dict[str, float] = {}
        for sp, kla in self.kLa.items():
            kla = float(kla)
            if kla <= 0.0:
                continue

            C_src = float(src_phase.n_mol.get(sp, 0.0)) / V_src
            C_snk = float(snk_phase.n_mol.get(sp, 0.0)) / V_snk

            dC = C_src - C_snk
            flux = kla * V_eff * dC  # mol/h (positive = src→snk)

            # Safety: don't transfer more than available
            if flux > 0.0:
                max_flux = float(src_phase.n_mol.get(sp, 0.0)) / max(float(dt_h), 1e-30)
                flux = min(flux, max_flux)
            elif flux < 0.0:
                max_flux = float(snk_phase.n_mol.get(sp, 0.0)) / max(float(dt_h), 1e-30)
                flux = max(flux, -max_flux)

            if abs(flux) > 0.0:
                flow[sp] = flux

        return flow

    def __repr__(self):
        species = ", ".join(sorted(self.kLa.keys()))
        return (f"DiffusiveLink({self._source_cv_key}.{self._source_phase_key}"
                f" <-> {self._sink_cv_key}.{self._sink_phase_key}, "
                f"species=[{species}])")
