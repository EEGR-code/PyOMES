# -*- coding: utf-8 -*-
"""TransferModel types for intra-CV phase mass transfer.

Two types cover all kinetic/equilibrium transfer cases:

* :class:`KineticTransferModel` — rate-limited transfer driven by a
  concentration difference.  ``transfer_basis`` selects whether the
  driving force acts on the **total** dissolved pool (default, correct
  for non-speciating gases such as O₂, N₂, CH₄, H₂) or on the
  **molecular** (uncharged) dissolved form (BSM2-style, correct for
  gases whose dissolved form is an acid or base: CO₂, NH₃, H₂S, VFAs).

* :class:`EquilibriumTransferModel` — instantaneous equilibrium
  partitioning (algebraic constraint, no ``k_transfer``).  Appropriate
  for species whose equilibration timescale is fast relative to the
  integration step (e.g. N₂, CH₄ in long-step whole-farm models).

These types are passed via the ``transfer_models`` kwarg of
:class:`~PyOMES.core.control_volume.ControlVolume`::

    cv = ControlVolume(
        phases={"gas": gas_phase, "liquid": liquid_phase},
        transfer_models={
            "O2":  KineticTransferModel(HenryEquilibrium(H_ref=1.3e-5, dlnH=1500.0), k_transfer=150.0),
            "CO2": KineticTransferModel(HenryEquilibrium(H_ref=3.4e-4, dlnH=2400.0), k_transfer=200.0,
                                        transfer_basis="molecular"),
            "N2":  EquilibriumTransferModel(HenryEquilibrium(H_ref=6.4e-6, dlnH=1300.0)),
        },
    )

The ``transfer_basis`` parameter on :class:`KineticTransferModel`:

``"total"`` (default)
    The driving force is ``k_transfer × V × (C*_total − C_total)``.
    ``alpha`` (the molecular fraction from the speciation engine) is
    passed to the partition model and corrects the effective Henry
    constant: ``kH_eff = kH / alpha``.  Correct for non-speciating
    gases; also a valid approximation for speciating gases.

``"molecular"``
    The driving force is ``k_transfer × V × (C*_mol − C_mol)``, where
    ``C_mol = alpha × C_total``.  ``alpha = 1.0`` is passed to the
    partition model (raw Henry constant, no amplification).  The
    relaxation rate is ~alpha× slower than the ``"total"`` formulation,
    matching BSM2 CO₂ kinetics.  Use for any volatile acid/base:
    CO₂, NH₃, H₂S, AceticAcid, Propionate, Butyrate, Valerate, …

See ``docs/design/TRANSFER_MODEL.md`` for the full design rationale.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..chemistry.partition import PartitionModel

_VALID_TRANSFER_BASES = frozenset({"total", "molecular"})


@dataclass
class KineticTransferModel:
    """Rate-limited inter-phase transfer for a single species.

    Parameters
    ----------
    partition_model : PartitionModel
        Thermodynamic model that describes the equilibrium distribution
        between the two phases (e.g. :class:`~PyOMES.chemistry.HenryEquilibrium`).
    k_transfer : float
        Volumetric mass transfer coefficient (h⁻¹).  Equivalent to
        ``kLa`` for gas-liquid systems.
    transfer_basis : str
        ``"total"`` (default) — driving force on the total dissolved
        pool; ``"molecular"`` — driving force on the molecular
        (uncharged) dissolved form only.  See module docstring for
        details.

    Raises
    ------
    ValueError
        If ``transfer_basis`` is not ``"total"`` or ``"molecular"``.
    """

    partition_model: "PartitionModel"
    k_transfer: float
    transfer_basis: str = "total"

    def __post_init__(self) -> None:
        if self.transfer_basis not in _VALID_TRANSFER_BASES:
            raise ValueError(
                f"transfer_basis must be 'total' or 'molecular', "
                f"got {self.transfer_basis!r}"
            )
        self.k_transfer = float(self.k_transfer)


@dataclass
class EquilibriumTransferModel:
    """Instantaneous equilibrium partitioning for a single species.

    No ``k_transfer`` — the species is distributed between phases at
    each step according to the partition model alone.  Appropriate when
    the equilibration timescale is fast relative to the integration
    step.

    Parameters
    ----------
    partition_model : PartitionModel
        Thermodynamic model describing the equilibrium distribution
        (e.g. :class:`~PyOMES.chemistry.HenryEquilibrium`).
    """

    partition_model: "PartitionModel"
