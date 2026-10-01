# -*- coding: utf-8 -*-
"""Configuration dataclasses for fermenter construction.

Each dataclass owns one concern (vessel geometry, gas feed, chemistry,
gas-liquid transfer, organism, substrate, simulation parameters) and
can be validated independently.  Together they fully describe a
fermenter simulation.

All configs are plain dataclasses with ``__post_init__`` validation.
They support round-trip serialisation via ``to_dict()`` and
``from_dict()`` for JSON/YAML workflows and parameter sweeps.
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field, fields, asdict
from enum import Enum
from typing import Any, Dict, List, Optional, Sequence, Union

from PyOMES.chemistry.species import Species
from PyOMES.thermo import ActivityModel, make_activity_model


# ════════════════════════════════════════════════════════════════════════
#  Transfer mode enum
# ════════════════════════════════════════════════════════════════════════

class TransferMode(str, Enum):
    """Gas-liquid transfer mode for a single species."""
    KINETIC = "kinetic"
    EQUILIBRIUM = "equilibrium"
    NONE = "none"


# ════════════════════════════════════════════════════════════════════════
#  VesselConfig
# ════════════════════════════════════════════════════════════════════════

@dataclass
class VesselConfig:
    """Vessel geometry, temperature, and initial gas composition.

    Parameters
    ----------
    V_total_L : float
        Total vessel volume (litres).  Must be > 0.
    headspace_frac : float
        Fraction of total volume that is headspace (0 < frac < 1).
    T_K : float
        Temperature (Kelvin).
    P_init_atm : float
        Initial headspace pressure (atm).
    gas_composition : dict
        ``{species_id: mole fraction}`` of the initial headspace gas;
        normalised by its sum. Empty (the default) means the headspace
        starts with no gas. E.g.
        :data:`~PyOMES.databases.bioprocess_basic.AIR`.
    """

    V_total_L: float = 2.0
    headspace_frac: float = 0.20
    T_K: float = 305.15
    P_init_atm: float = 1.0
    gas_composition: Dict[str, float] = field(default_factory=dict)

    def __post_init__(self):
        if self.V_total_L <= 0:
            raise ValueError(f"V_total_L must be > 0, got {self.V_total_L}")
        if not (0.0 < self.headspace_frac < 1.0):
            raise ValueError(
                f"headspace_frac must be in (0, 1), got {self.headspace_frac}"
            )
        if self.T_K <= 0:
            raise ValueError(f"T_K must be > 0, got {self.T_K}")
        if self.P_init_atm <= 0:
            raise ValueError(f"P_init_atm must be > 0, got {self.P_init_atm}")
        self.gas_composition = dict(self.gas_composition)
        negative = {sp: y for sp, y in self.gas_composition.items() if y < 0}
        if negative:
            raise ValueError(f"gas_composition fractions must be >= 0, got {negative}")

    @property
    def V_headspace_L(self) -> float:
        return self.V_total_L * self.headspace_frac

    @property
    def V_liquid_L(self) -> float:
        return self.V_total_L * (1.0 - self.headspace_frac)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "VesselConfig":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ════════════════════════════════════════════════════════════════════════
#  GasFeedConfig
# ════════════════════════════════════════════════════════════════════════

@dataclass
class GasFeedConfig:
    """Continuous gas feed (sparging) parameters.

    Parameters
    ----------
    vvm_min : float
        Gas volume per liquid volume per minute (L_gas/L_liq/min).
        Set to 0 for no sparging (e.g. well plate).
    composition : dict
        Inlet gas mole fractions by species id, e.g.
        ``{"O2": 0.21, "N2": 0.79}``; normalised internally. No default:
        a feed with ``vvm_min > 0`` must say what it feeds.
    P_inlet_atm : float
        Inlet gas pressure (atm).
    """

    vvm_min: float = 1.0
    composition: Dict[str, float] = field(default_factory=dict)
    P_inlet_atm: float = 1.0

    def __post_init__(self):
        if self.vvm_min < 0:
            raise ValueError(f"vvm_min must be >= 0, got {self.vvm_min}")
        if self.P_inlet_atm <= 0:
            raise ValueError(f"P_inlet_atm must be > 0, got {self.P_inlet_atm}")
        if self.vvm_min > 0 and not self.composition:
            raise ValueError(
                "A gas feed with vvm_min > 0 needs a composition, "
                "e.g. composition={'O2': 0.21, 'N2': 0.79}."
            )
        # Normalise composition
        raw = dict(self.composition)
        y_sum = sum(max(0.0, float(v)) for v in raw.values())
        if y_sum > 0.0:
            self.composition = {k: max(0.0, float(v)) / y_sum for k, v in raw.items()}
        else:
            self.composition = {k: 0.0 for k in raw}

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "GasFeedConfig":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ════════════════════════════════════════════════════════════════════════
#  SpeciesTransferConfig (per-species entry)
# ════════════════════════════════════════════════════════════════════════

@dataclass
class SpeciesTransferConfig:
    """Transfer configuration for a single gas species.

    Parameters
    ----------
    mode : TransferMode or str
        ``"kinetic"``, ``"equilibrium"``, or ``"none"``.
    kLa_per_h : float
        Volumetric mass transfer coefficient (1/h).  Only used when
        ``mode`` is ``"kinetic"``.
    henry_mol_L_atm : float or None
        Henry constant (mol/L/atm) at the vessel temperature.  If
        ``None``, the factory takes the species' partition model from the
        model's ``chemistry_db``.
    """

    mode: TransferMode = TransferMode.EQUILIBRIUM
    kLa_per_h: float = 0.0
    henry_mol_L_atm: Optional[float] = None

    def __post_init__(self):
        if isinstance(self.mode, str):
            self.mode = TransferMode(self.mode.lower().strip())
        if self.kLa_per_h < 0:
            raise ValueError(f"kLa_per_h must be >= 0, got {self.kLa_per_h}")

    def to_dict(self) -> dict:
        d = asdict(self)
        d["mode"] = self.mode.value
        return d


@dataclass
class TransferConfig:
    """Gas-liquid transfer configuration for all species.

    Parameters
    ----------
    species : dict
        ``{species_id: SpeciesTransferConfig}``.  Species not listed
        are not transferred; the default is no transfer.
    """

    species: Dict[str, SpeciesTransferConfig] = field(default_factory=dict)

    def __post_init__(self):
        # Convert any raw dicts to SpeciesTransferConfig
        cleaned = {}
        for sp, cfg in self.species.items():
            if isinstance(cfg, dict):
                cleaned[sp] = SpeciesTransferConfig(**cfg)
            elif isinstance(cfg, SpeciesTransferConfig):
                cleaned[sp] = cfg
            else:
                raise TypeError(
                    f"Expected SpeciesTransferConfig or dict for species {sp!r}, "
                    f"got {type(cfg).__name__}"
                )
        self.species = cleaned

    @classmethod
    def kinetic(
        cls,
        kLa: Dict[str, float],
        equilibrium: Sequence[str] = (),
    ) -> "TransferConfig":
        """Kinetic transfer for the species in *kLa* (``{id: kLa_per_h}``),
        then equilibrium transfer for those in *equilibrium*."""
        species = {
            sp: SpeciesTransferConfig(mode=TransferMode.KINETIC, kLa_per_h=k)
            for sp, k in kLa.items()
        }
        for sp in equilibrium:
            species[sp] = SpeciesTransferConfig(mode=TransferMode.EQUILIBRIUM)
        return cls(species=species)

    @classmethod
    def equilibrium(cls, species_ids: Sequence[str]) -> "TransferConfig":
        """Equilibrium transfer for each of *species_ids*."""
        return cls(species={
            sp: SpeciesTransferConfig(mode=TransferMode.EQUILIBRIUM)
            for sp in species_ids
        })

    def to_dict(self) -> dict:
        return {"species": {sp: cfg.to_dict() for sp, cfg in self.species.items()}}

    @classmethod
    def from_dict(cls, d: dict) -> "TransferConfig":
        species = {}
        for sp, cfg in d.get("species", {}).items():
            if isinstance(cfg, dict):
                species[sp] = SpeciesTransferConfig(**cfg)
            else:
                species[sp] = cfg
        return cls(species=species)


# ════════════════════════════════════════════════════════════════════════
#  ChemistryConfig
# ════════════════════════════════════════════════════════════════════════

@dataclass
class ChemistryConfig:
    """Speciation and aqueous chemistry settings.

    Parameters
    ----------
    activity_model : str or activity model object
        ``"ideal"`` (default), ``"davies"``, ``"sit"``, or a model object such
        as ``SITLiquidModel(epsilon=...)``. Checked at construction by
        :func:`~PyOMES.thermo.make_activity_model` and stored as given.

    Notes
    -----
    :meth:`to_dict` keeps a model object as the object itself, so
    :meth:`from_dict` gets the same model back. A config with a named model
    can be written to JSON or YAML; one holding a model object cannot.
    """

    activity_model: Union[str, ActivityModel] = "ideal"

    def __post_init__(self) -> None:
        make_activity_model(self.activity_model)  # raises on a bad name or object

    def to_dict(self) -> dict:
        # Not asdict(): it would turn a dataclass model (e.g. SITLiquidModel)
        # into a plain dict.
        return {
            f.name: (
                getattr(self, f.name) if f.name == "activity_model"
                else copy.deepcopy(getattr(self, f.name))
            )
            for f in fields(self)
        }

    @classmethod
    def from_dict(cls, d: dict) -> "ChemistryConfig":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ════════════════════════════════════════════════════════════════════════
#  OrganismConfig
# ════════════════════════════════════════════════════════════════════════

def _check_definition(name: str, given: Any, atoms: Any, MW: Any) -> None:
    """Check that *given* (an id or a ``Species``) with *atoms* / *MW* is
    one complete way of naming or defining a species."""
    if isinstance(given, Species):
        if atoms is not None or MW is not None:
            raise ValueError(
                f"{name}={given.id!r} is a Species, which carries its own atoms "
                "and MW; do not pass atoms= or MW= with it."
            )
    elif not isinstance(given, str) or not given:
        raise TypeError(
            f"{name} must be a species id (str) or a Species, got {given!r}"
        )
    if atoms is None and MW is not None:
        raise ValueError(
            f"MW= for {name}={given!r} needs atoms=: a molar mass alone does "
            "not define a species (with atoms=, MW defaults to the value "
            "computed from them)."
        )
    if MW is not None and MW <= 0:
        raise ValueError(f"MW must be > 0, got {MW}")


def _definition_dict(config: Any) -> dict:
    """``asdict`` that keeps a ``Species`` field as the object itself."""
    return {
        f.name: (
            getattr(config, f.name) if isinstance(getattr(config, f.name), Species)
            else copy.deepcopy(getattr(config, f.name))
        )
        for f in fields(config)
    }


@dataclass
class OrganismConfig:
    """Organism identity and composition.

    The organism is named by an id, defined by an id with ``atoms``, or given
    as a :class:`~PyOMES.chemistry.species.Species`. The factory resolves it
    against the model's species: an id alone must be among them; a
    definition is added to them, and one that differs from an existing
    species of the same id raises unless ``overwrite=True``.

    Parameters
    ----------
    organism : str or Species
        The organism's id (e.g. ``"Yeast"``, ``"E_coli"``) or its ``Species``.
    atoms : dict or None
        Elemental composition, e.g. ``{"C": 1, "H": 1.61, "O": 0.56, "N": 0.16}``.
        With an id, defines the organism.
    MW : float or None
        Molecular weight (g/mol). Only with ``atoms``; computed from them
        when not given.
    balance_basis : str
        ``"CHO"`` or ``"CHNO"``.
    n_source_id : str or None
        Id of the nitrogen source, resolved against the model's species.
        Required for ``"CHNO"``.
    o2_id, co2_id, h2o_id : str
        Ids of the oxygen consumed and the carbon dioxide and water produced
        by the growth reaction (default ``"O2"``, ``"CO2"``, ``"H2O"``),
        resolved against the model's species.
    overwrite : bool
        Replace a model species of the same id whose data differ from this
        definition, instead of raising.
    """

    organism: Union[str, Species]
    atoms: Optional[Dict[str, float]] = None
    MW: Optional[float] = None
    balance_basis: str = "CHO"
    n_source_id: Optional[str] = None
    o2_id: str = "O2"
    co2_id: str = "CO2"
    h2o_id: str = "H2O"
    overwrite: bool = False

    def __post_init__(self):
        self.balance_basis = self.balance_basis.upper().strip()
        if self.balance_basis not in ("CHO", "CHNO"):
            raise ValueError(
                f"balance_basis must be 'CHO' or 'CHNO', got {self.balance_basis!r}"
            )
        _check_definition("organism", self.organism, self.atoms, self.MW)
        if self.balance_basis == "CHNO" and self.n_source_id is None:
            raise ValueError(
                "balance_basis='CHNO' needs a nitrogen source: pass "
                "n_source_id, e.g. n_source_id='NH3'."
            )

    def to_dict(self) -> dict:
        return _definition_dict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "OrganismConfig":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ════════════════════════════════════════════════════════════════════════
#  SubstrateConfig
# ════════════════════════════════════════════════════════════════════════

@dataclass
class SubstrateConfig:
    """Substrate identity, kinetic parameters, and yield.

    The substrate is named or defined as the organism is (see
    :class:`OrganismConfig`) and resolved against the model's species by
    the same rules.

    Parameters
    ----------
    substrate : str or Species
        The substrate's id (e.g. ``"AceticAcid"``) or its ``Species``.
    atoms : dict or None
        Elemental composition, e.g. ``{"C": 2, "H": 4, "O": 2}``. With an
        id, defines the substrate.
    MW : float or None
        Molecular weight (g/mol). Only with ``atoms``; computed from them
        when not given.
    mu_max : float
        Maximum specific growth rate (1/h).
    Ks : float
        Monod half-saturation constant (g/L).
    yield_gX_gS : float
        Biomass yield (g biomass / g substrate consumed).  Must be > 0.
    kinetics : GrowthKinetics or None
        Rate law; default Monod with ``mu_max`` and ``Ks``.
    overwrite : bool
        Replace a model species of the same id whose data differ from this
        definition, instead of raising.
    """

    substrate: Union[str, Species]
    atoms: Optional[Dict[str, float]] = None
    MW: Optional[float] = None
    mu_max: float = 0.5
    Ks: float = 5e-3
    yield_gX_gS: float = 0.36
    kinetics: Optional[Any] = None
    overwrite: bool = False

    def __post_init__(self):
        if self.kinetics is None:
            # Validate mu_max and Ks only when using default Monod
            if self.mu_max <= 0:
                raise ValueError(f"mu_max must be > 0, got {self.mu_max}")
            if self.Ks < 0:
                raise ValueError(f"Ks must be >= 0, got {self.Ks}")
        if self.yield_gX_gS <= 0:
            raise ValueError(f"yield_gX_gS must be > 0, got {self.yield_gX_gS}")
        _check_definition("substrate", self.substrate, self.atoms, self.MW)

    def to_dict(self) -> dict:
        return _definition_dict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "SubstrateConfig":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ════════════════════════════════════════════════════════════════════════
#  SimulationConfig
# ════════════════════════════════════════════════════════════════════════

@dataclass
class SimulationConfig:
    """Simulation control parameters.

    Parameters
    ----------
    tau : float
        Batch time (hours).  Must be > 0.
    n_steps : int
        Number of timesteps.  Must be > 0.
    t_lag : float
        Lag phase duration (hours).  Growth is suppressed for ``t < t_lag``.
    """

    tau: float = 5.0
    n_steps: int = 1000
    t_lag: float = 0.0

    def __post_init__(self):
        if self.tau <= 0:
            raise ValueError(f"tau must be > 0, got {self.tau}")
        if self.n_steps <= 0:
            raise ValueError(f"n_steps must be > 0, got {self.n_steps}")
        if self.t_lag < 0:
            raise ValueError(f"t_lag must be >= 0, got {self.t_lag}")

    @property
    def dt_h(self) -> float:
        """Timestep duration (hours)."""
        return self.tau / self.n_steps

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "SimulationConfig":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})
