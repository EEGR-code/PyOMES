"""Guard: the gas constant R is written out once, in ``PyOMES/units.py``.

Everything else in ``PyOMES/`` and ``models/`` should import ``R_J_PER_MOL_K`` or
``R_L_ATM_PER_MOL_K`` from :mod:`PyOMES.units`. Separate copies drift apart (the
rounded 0.0820574 sits 4.1e-7 above the CODATA value), and a gas-liquid
calculation that mixes two of them is quietly inconsistent.

The scan reads Python number tokens, so comments, docstrings and strings are
ignored, and any spelling of the value is caught (``8.314``, ``8.31446``,
``0.0820574``, ``8.2057e-2``). A token counts as R when it is within 0.1 % of
8.314462618 J/(mol·K) or 0.082057366 L·atm/(mol·K).

``_KNOWN_COPIES`` lists the duplicates that still exist, as
``{path: {token: how many times}}``. It should only shrink: replace a copy with an
import, then delete its entry. A new literal that is not listed fails the first
test, and a listed copy that has gone fails the second, so the list cannot go
stale.
"""

import tokenize
from pathlib import Path

import pytest

from PyOMES import units

REPO_ROOT = Path(__file__).resolve().parents[2]
SCAN_DIRS = ("PyOMES", "models")
ROOT_FILE = "PyOMES/units.py"

# Fixed reference values for detection, deliberately not read from ``units`` so
# that changing the root cannot change what the scan looks for.
_R_J_PER_MOL_K = 8.314462618
_R_L_ATM_PER_MOL_K = 0.082057366
_REL_TOL = 1e-3

_KNOWN_COPIES = {
    "PyOMES/chemistry/partition.py": {"0.0820574": 1},
    "PyOMES/control/cv_loops.py": {"0.08205736608095958": 1},
    "PyOMES/core/phases.py": {"0.0820574": 1},
    "PyOMES/equilibria/peng_robinson.py": {"0.0820574": 1},
    "PyOMES/reactions/plots.py": {"8.314": 1},
    "models/vlmodels/adm1/base.py": {"8.31446": 1},
    "models/vlmodels/adm1/bsm2.py": {"8.31446": 1},
    "models/vlmodels/adm1/bsm2_direct.py": {"8.31446": 1},
}


def _is_gas_constant(value: float) -> bool:
    return any(
        abs(value / ref - 1.0) < _REL_TOL
        for ref in (_R_J_PER_MOL_K, _R_L_ATM_PER_MOL_K)
    )


def _scan_gas_constant_literals() -> dict:
    """Return ``{path: {token: count}}`` for every R-valued number token."""
    found: dict = {}
    for scan_dir in SCAN_DIRS:
        for path in sorted((REPO_ROOT / scan_dir).rglob("*.py")):
            if "__pycache__" in path.parts:
                continue
            rel = path.relative_to(REPO_ROOT).as_posix()
            if rel == ROOT_FILE:
                continue
            with tokenize.open(path) as fh:
                for tok in tokenize.generate_tokens(fh.readline):
                    if tok.type != tokenize.NUMBER:
                        continue
                    try:
                        value = float(tok.string.replace("_", ""))
                    except ValueError:  # hex, complex
                        continue
                    if _is_gas_constant(value):
                        counts = found.setdefault(rel, {})
                        counts[tok.string] = counts.get(tok.string, 0) + 1
    return found


def test_no_unlisted_gas_constant_literals():
    found = _scan_gas_constant_literals()
    unlisted = []
    for rel, counts in found.items():
        for token, n in counts.items():
            allowed = _KNOWN_COPIES.get(rel, {}).get(token, 0)
            if n > allowed:
                unlisted.append(f"{rel}: {token} (found {n}, listed {allowed})")
    assert not unlisted, (
        "Gas-constant literal outside PyOMES/units.py. Import R_J_PER_MOL_K or "
        "R_L_ATM_PER_MOL_K from PyOMES.units instead:\n  " + "\n  ".join(unlisted)
    )


def test_known_copies_list_has_no_stale_entries():
    found = _scan_gas_constant_literals()
    stale = []
    for rel, tokens in _KNOWN_COPIES.items():
        for token, n in tokens.items():
            actual = found.get(rel, {}).get(token, 0)
            if actual != n:
                stale.append(f"{rel}: {token} (listed {n}, found {actual})")
    assert not stale, (
        "_KNOWN_COPIES no longer matches the source. Remove or update the "
        "entry:\n  " + "\n  ".join(stale)
    )


def test_units_holds_the_reference_values():
    """The root itself: CODATA 2018 for J, and L·atm agreeing to 1e-12."""
    assert units.R_J_PER_MOL_K == 8.31446261815324
    assert units.R_L_ATM_PER_MOL_K == pytest.approx(0.08205736608096, rel=1e-12)
    # The two roots describe the same constant: R[L·atm] = R[J] / (Pa/atm * L/m3).
    assert units.R_L_ATM_PER_MOL_K == pytest.approx(
        units.R_J_PER_MOL_K / (units.PA_PER_ATM / 1000.0), rel=1e-12
    )


def test_scan_detects_a_new_literal(tmp_path, monkeypatch):
    """The scan would notice a fresh copy, so an empty result is meaningful."""
    pkg = tmp_path / "PyOMES"
    pkg.mkdir()
    (pkg / "units.py").write_text("R = 8.31446261815324\n")
    (pkg / "new_module.py").write_text(
        "# R = 8.314 in a comment is ignored\n"
        'DOC = "0.0820574 in a string is ignored"\n'
        "R_GAS = 8.31446\n"
        "R_LA = 8.2057e-2\n"
        "N = 1_000_000\n"
    )
    monkeypatch.setitem(globals(), "REPO_ROOT", tmp_path)
    monkeypatch.setitem(globals(), "SCAN_DIRS", ("PyOMES",))
    assert _scan_gas_constant_literals() == {
        "PyOMES/new_module.py": {"8.31446": 1, "8.2057e-2": 1}
    }
