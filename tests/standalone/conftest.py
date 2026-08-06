"""Shared fixtures for the standalone test suite.

Every fixture here is bioSTEAM-free.  The suite must pass even when
bioSTEAM is not installed.
"""

import os
import sys
import pytest

# Ensure the src and models directories are on the path
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
_SRC = os.path.join(_ROOT, "src")
_MODELS = os.path.join(_ROOT, "models")
for p in (_SRC, _MODELS):
    if p not in sys.path:
        sys.path.insert(0, p)


# ── Chemical registry ──────────────────────────────────────────────────

@pytest.fixture
def registry():
    """Default ChemicalRegistry with all built-in compounds."""
    from PyOMES.chemistry.compounds import ChemicalRegistry
    return ChemicalRegistry.default()


# ── Feed states ────────────────────────────────────────────────────────

@pytest.fixture
def simple_feed(registry):
    """Minimal feed: acetate + yeast only."""
    from PyOMES.stream_adapter import FeedState
    return FeedState.from_mass_concentrations(
        mass_g_L={"AceticAcid": 1.0, "Yeast": 0.1},
        volume_L=1.0,
        T_K=305.15,
        registry=registry,
    )


@pytest.fixture
def rich_feed(registry):
    """Realistic feed with substrates, N-source, salts, and trace metals."""
    from PyOMES.stream_adapter import FeedState
    return FeedState.from_mixed_concentrations(
        mass_g_L={
            "AceticAcid": 1.0,
            "Yeast": 0.1,
        },
        molar_mol_L={
            "NH3": 0.015,
            "AmmoniumSulfate": 0.015,
            "KH2PO4": 0.007,
            "MgSO4": 0.004,
            "ZnSO4": 1.4e-3,
            "CaCl2": 5.0e-3,
            "MnCl2": 4.0e-3,
            "CoCl2": 7.7e-4,
        },
        volume_L=1.0,
        T_K=305.15,
        registry=registry,
    )


# ── Controllers ────────────────────────────────────────────────────────

@pytest.fixture
def pressure_ctrl():
    from PyOMES import PressureReliefController
    return PressureReliefController(
        P_set_atm=1.10, diameter_m=5e-2, sample_period_s=1,
    )


@pytest.fixture
def ph_ctrl():
    from PyOMES import PHController
    return PHController(
        setpoint=6.6, chemical_id="H3PO4", base_chemical_id="KOH",
        Kp=1.0, Ki=1.0, max_add_molL_hr=10,
    )


