# -*- coding: utf-8 -*-
"""Configuration dataclasses for fermenter construction.

Each dataclass owns one concern (vessel geometry, gas feed, chemistry,
gas-liquid transfer, organism, substrate, simulation parameters) and
can be validated independently.  Together they fully describe a
fermenter simulation and can be consumed by
``Fermenter.from_config()`` (Stage E).

All configs are plain dataclasses with ``__post_init__`` validation.
They support round-trip serialisation via ``to_dict()`` and
``from_dict()`` for JSON/YAML workflows and parameter sweeps.
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Any, Dict, List, Optional, Sequence


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
    yO2_init : float
        Initial O₂ mole fraction in headspace.
    yCO2_init : float
        Initial CO₂ mole fraction in headspace.
    yN2_init : float or None
        Initial N₂ mole fraction.  If ``None``, computed as
        ``1 − yO2 − yCO2`` (balance).
    """

    V_total_L: float = 2.0
    headspace_frac: float = 0.20
    T_K: float = 305.15
    P_init_atm: float = 1.0
    yO2_init: float = 0.2095
    yCO2_init: float = 0.0004
    yN2_init: Optional[float] = None

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
        if self.yN2_init is None:
            self.yN2_init = max(0.0, 1.0 - self.yO2_init - self.yCO2_init)

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
        Inlet gas mole fractions, e.g. ``{"O2": 0.21, "N2": 0.79}``.
        Normalised internally.
    P_inlet_atm : float
        Inlet gas pressure (atm).
    """

    vvm_min: float = 1.0
    composition: Dict[str, float] = field(
        default_factory=lambda: {"O2": 0.21, "N2": 0.79}
    )
    P_inlet_atm: float = 1.0

    def __post_init__(self):
        if self.vvm_min < 0:
            raise ValueError(f"vvm_min must be >= 0, got {self.vvm_min}")
        if self.P_inlet_atm <= 0:
            raise ValueError(f"P_inlet_atm must be > 0, got {self.P_inlet_atm}")
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
        ``None``, the factory will compute it from the temperature-
        dependent correlation.
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
        are not transferred.
    kLa_CO2_ratio : float
        If CO₂ kLa is not explicitly set, derive it as
        ``kLa_O2 × kLa_CO2_ratio``.  Default 0.9 (diffusivity scaling).
    """

    species: Dict[str, SpeciesTransferConfig] = field(default_factory=dict)
    kLa_CO2_ratio: float = 0.9

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
    def default_kinetic(cls, kLa_O2: float = 150.0, kLa_CO2_ratio: float = 0.9) -> "TransferConfig":
        """Create a default config with kinetic O₂/CO₂ and equilibrium N₂."""
        return cls(
            species={
                "O2": SpeciesTransferConfig(
                    mode=TransferMode.KINETIC, kLa_per_h=kLa_O2,
                ),
                "CO2": SpeciesTransferConfig(
                    mode=TransferMode.KINETIC, kLa_per_h=kLa_O2 * kLa_CO2_ratio,
                ),
                "N2": SpeciesTransferConfig(mode=TransferMode.EQUILIBRIUM),
            },
            kLa_CO2_ratio=kLa_CO2_ratio,
        )

    @classmethod
    def default_equilibrium(cls) -> "TransferConfig":
        """Create a default config with all species at equilibrium."""
        return cls(
            species={
                "O2": SpeciesTransferConfig(mode=TransferMode.EQUILIBRIUM),
                "CO2": SpeciesTransferConfig(mode=TransferMode.EQUILIBRIUM),
                "N2": SpeciesTransferConfig(mode=TransferMode.EQUILIBRIUM),
            },
        )

    def to_dict(self) -> dict:
        return {
            "species": {sp: cfg.to_dict() for sp, cfg in self.species.items()},
            "kLa_CO2_ratio": self.kLa_CO2_ratio,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "TransferConfig":
        species = {}
        for sp, cfg in d.get("species", {}).items():
            if isinstance(cfg, dict):
                species[sp] = SpeciesTransferConfig(**cfg)
            else:
                species[sp] = cfg
        return cls(
            species=species,
            kLa_CO2_ratio=d.get("kLa_CO2_ratio", 0.9),
        )


# ════════════════════════════════════════════════════════════════════════
#  ChemistryConfig
# ════════════════════════════════════════════════════════════════════════

@dataclass
class ChemistryConfig:
    """Speciation and aqueous chemistry settings.

    Parameters
    ----------
    use_activity : bool
        Enable activity coefficient corrections.
    activity_model : str
        Activity model name (``"davies"``, ``"ideal"``).
    acid_pKas : dict
        ``{acid_id: pKa_or_list}`` for weak acid systems.
    """

    use_activity: bool = False
    activity_model: str = "davies"
    acid_pKas: Dict[str, Any] = field(
        default_factory=lambda: {
            "AceticAcid": 4.76,
            "PropionicAcid": 4.87,
            "ButyricAcid": 4.82,
            "CitricAcid": [3.13, 4.76, 6.40],
        }
    )

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "ChemistryConfig":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ════════════════════════════════════════════════════════════════════════
#  OrganismConfig
# ════════════════════════════════════════════════════════════════════════

@dataclass
class OrganismConfig:
    """Organism identity and composition.

    If ``atoms`` and ``MW`` are ``None``, the factory will look them up
    from the :class:`~fermenter.chemistry.compounds.ChemicalRegistry`.

    Parameters
    ----------
    organism_id : str
        Identifier (e.g. ``"Yeast"``, ``"E_coli"``).  Must match either
        a registry entry or have explicit ``atoms``/``MW``.
    atoms : dict or None
        Elemental composition, e.g. ``{"C": 1, "H": 1.61, "O": 0.56, "N": 0.16}``.
    MW : float or None
        Molecular weight (g/mol).
    balance_basis : str
        ``"CHO"`` or ``"CHNO"``.
    n_source_id : str
        Nitrogen source chemical ID (for CHNO mode).
    """

    organism_id: str = "Yeast"
    atoms: Optional[Dict[str, float]] = None
    MW: Optional[float] = None
    balance_basis: str = "CHO"
    n_source_id: str = "NH3"

    def __post_init__(self):
        self.balance_basis = self.balance_basis.upper().strip()
        if self.balance_basis not in ("CHO", "CHNO"):
            raise ValueError(
                f"balance_basis must be 'CHO' or 'CHNO', got {self.balance_basis!r}"
            )
        if self.MW is not None and self.MW <= 0:
            raise ValueError(f"MW must be > 0, got {self.MW}")

    def resolve(self, registry=None) -> "OrganismConfig":
        """Return a copy with atoms/MW filled from the registry if needed.

        Parameters
        ----------
        registry : ChemicalRegistry or None
            If ``None``, uses the default registry.

        Returns
        -------
        OrganismConfig
            A new config with atoms and MW guaranteed non-None.

        Raises
        ------
        KeyError
            If the organism is not in the registry and atoms/MW are not set.
        """
        if self.atoms is not None and self.MW is not None:
            return copy.copy(self)

        if registry is None:
            from PyOMES.chemistry.compounds import ChemicalRegistry
            registry = ChemicalRegistry.default()

        chem = registry[self.organism_id]
        atoms = dict(self.atoms) if self.atoms is not None else dict(
            getattr(chem, "atoms", {}) or {}
        )
        MW = self.MW if self.MW is not None else float(
            getattr(chem, "MW", 0.0) or 0.0
        )
        return OrganismConfig(
            organism_id=self.organism_id,
            atoms=atoms,
            MW=MW,
            balance_basis=self.balance_basis,
            n_source_id=self.n_source_id,
        )

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "OrganismConfig":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ════════════════════════════════════════════════════════════════════════
#  SubstrateConfig
# ════════════════════════════════════════════════════════════════════════

@dataclass
class SubstrateConfig:
    """Substrate identity, kinetic parameters, and yield.

    If ``atoms`` and ``MW`` are ``None``, the factory will look them up
    from the :class:`~fermenter.chemistry.compounds.ChemicalRegistry`.

    Parameters
    ----------
    substrate_id : str
        Chemical identifier (e.g. ``"AceticAcid"``).
    atoms : dict or None
        Elemental composition, e.g. ``{"C": 2, "H": 4, "O": 2}``.
    MW : float or None
        Molecular weight (g/mol).
    mu_max : float
        Maximum specific growth rate (1/h).
    Ks : float
        Monod half-saturation constant (g/L).
    yield_gX_gS : float
        Biomass yield (g biomass / g substrate consumed).  Must be > 0.
    """

    substrate_id: str = "AceticAcid"
    atoms: Optional[Dict[str, float]] = None
    MW: Optional[float] = None
    mu_max: float = 0.5
    Ks: float = 5e-3
    yield_gX_gS: float = 0.36
    kinetics: Optional[Any] = None

    def __post_init__(self):
        if self.kinetics is None:
            # Validate mu_max and Ks only when using default Monod
            if self.mu_max <= 0:
                raise ValueError(f"mu_max must be > 0, got {self.mu_max}")
            if self.Ks < 0:
                raise ValueError(f"Ks must be >= 0, got {self.Ks}")
        if self.yield_gX_gS <= 0:
            raise ValueError(f"yield_gX_gS must be > 0, got {self.yield_gX_gS}")
        if self.MW is not None and self.MW <= 0:
            raise ValueError(f"MW must be > 0, got {self.MW}")

    def resolve(self, registry=None) -> "SubstrateConfig":
        """Return a copy with atoms/MW filled from the registry if needed.

        Parameters
        ----------
        registry : ChemicalRegistry or None
            If ``None``, uses the default registry.

        Returns
        -------
        SubstrateConfig
            A new config with atoms and MW guaranteed non-None.

        Raises
        ------
        KeyError
            If the substrate is not in the registry and atoms/MW are not set.
        """
        if self.atoms is not None and self.MW is not None:
            return copy.copy(self)

        if registry is None:
            from PyOMES.chemistry.compounds import ChemicalRegistry
            registry = ChemicalRegistry.default()

        chem = registry[self.substrate_id]
        atoms = dict(self.atoms) if self.atoms is not None else dict(
            getattr(chem, "atoms", {}) or {}
        )
        MW = self.MW if self.MW is not None else float(
            getattr(chem, "MW", 0.0) or 0.0
        )
        return SubstrateConfig(
            substrate_id=self.substrate_id,
            atoms=atoms,
            MW=MW,
            mu_max=self.mu_max,
            Ks=self.Ks,
            yield_gX_gS=self.yield_gX_gS,
            kinetics=self.kinetics,
        )

    def to_dict(self) -> dict:
        return asdict(self)

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
