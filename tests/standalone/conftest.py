"""Shared fixtures for the standalone test suite.

Every fixture here is bioSTEAM-free.  The suite must pass even when
bioSTEAM is not installed.
"""

import os
import sys
import pytest

# Ensure the repo root (for PyOMES) and models directory are on the path
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
_MODELS = os.path.join(_ROOT, "models")
for p in (_ROOT, _MODELS):
    if p not in sys.path:
        sys.path.insert(0, p)


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
        setpoint=6.6, acid_dose="H3PO4", base_dose={"K+": 1, "OH-": 1},
        Kp=1.0, Ki=1.0, max_add_molL_hr=10,
    )


