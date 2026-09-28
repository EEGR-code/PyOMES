"""Guard: the relation exp(−(E/R)(1/T − 1/T_ref)) is written out once.

The van 't Hoff, Clausius–Clapeyron and Arrhenius corrections all share the
exponent −(E/R)·(1/T − 1/T_ref). It is written in exactly one place,
:func:`PyOMES.thermo.temperature_correction.ln_correction`; everything else in
``PyOMES/`` and ``models/`` calls that function or one of the wrappers built on it
(``vant_hoff_K``, ``vant_hoff_log_K``, ``henry_constant``, ``clausius_clapeyron``,
``arrhenius_factor``). Separate copies drift apart: before there was one source,
fifteen copies differed in arithmetic order and in their edge-case rules, so the
same constant could come out differently by path.

The scan walks each module's syntax tree for a subtraction of two reciprocals,
``1/a − 1/b`` (the numerator a literal 1 or 1.0, either operand order, any
denominator: ``T_K``, ``float(T_K)``, ``self.T_ref_K``, an array), which is how
every copy wrote the ``1/T − 1/T_ref`` term. Comments, docstrings and strings are
not code and are ignored.

Tests are deliberately not scanned: ``test_temperature_correction.py`` writes the
formulas out in full as independent references, to pin the kernel's arithmetic bit
for bit.
"""

import ast
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCAN_DIRS = ("PyOMES", "models")
KERNEL_FILE = "PyOMES/thermo/temperature_correction.py"


def _is_reciprocal(node: ast.AST) -> bool:
    return (
        isinstance(node, ast.BinOp)
        and isinstance(node.op, ast.Div)
        and isinstance(node.left, ast.Constant)
        and node.left.value == 1
        and not isinstance(node.left.value, bool)
    )


def _reciprocal_differences(source: str) -> int:
    """Count ``1/a − 1/b`` expressions in *source*."""
    return sum(
        1
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.BinOp)
        and isinstance(node.op, ast.Sub)
        and _is_reciprocal(node.left)
        and _is_reciprocal(node.right)
    )


def _scan() -> dict:
    """Return ``{path: count}`` for every scanned file with the pattern."""
    found = {}
    for scan_dir in SCAN_DIRS:
        for path in sorted((REPO_ROOT / scan_dir).rglob("*.py")):
            if "__pycache__" in path.parts:
                continue
            n = _reciprocal_differences(path.read_text(encoding="utf-8"))
            if n:
                found[path.relative_to(REPO_ROOT).as_posix()] = n
    return found


def test_relation_is_written_only_in_the_kernel():
    found = _scan()
    assert found == {KERNEL_FILE: 1}, (
        "The 1/T − 1/T_ref term appears outside "
        f"{KERNEL_FILE} (or not exactly once inside it). Call ln_correction or one "
        "of its wrappers instead of writing the formula out:\n  "
        + "\n  ".join(f"{p}: {n}" for p, n in found.items())
    )


def test_scan_catches_every_form_it_must():
    """The forms every copy used, so an empty result elsewhere is meaningful."""
    caught = [
        "x = 1.0 / T_K - 1.0 / T_ref",                        # plain
        "x = -dH / R * (1.0 / T_K - 1.0 / self.T_ref_K)",     # inside an exponent
        "x = (1.0 / float(T_K) - 1.0 / float(T_ref_K))",      # float() conversions
        "x = Ea_R * (1.0 / T_ref_K - 1.0 / T_K)",             # reversed order
        "x = k * (1.0 / rxn.T_ref_K - 1.0 / T_arr)",          # attribute and array
        "x = 1 / T - 1 / T0",                                 # integer numerators
        "def f():\n    return math.exp(d * (1.0 / T - 1.0 / Tr))",  # in a function
    ]
    for src in caught:
        assert _reciprocal_differences(src) == 1, src


def test_scan_ignores_what_is_not_the_pattern():
    ignored = [
        '"""The 1.0 / T - 1.0 / T_ref term, in a docstring."""',
        "# 1.0 / T - 1.0 / T_ref in a comment",
        "s = '1.0 / T - 1.0 / T_ref'",
        "x = 1.0 / a + 1.0 / b",       # a sum, not a difference
        "x = 2.0 / a - 1.0 / b",       # numerator is not 1
        "x = a / T - b / T_ref",
        "x = 1.0 / (T - T_ref)",
        "flag = True / 2 - True / 3",  # bool is not the literal 1
    ]
    for src in ignored:
        assert _reciprocal_differences(src) == 0, src
