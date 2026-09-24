"""Guard: ``PyOMES/chemistry/`` depends on nothing in ``PyOMES`` except ``units``.

``chemistry/`` holds species declarations and the phase-partition protocols. It
sits below ``reactions/``, ``core/`` and everything else, so no module in it may
import from them. Unlike ``test_import_graph_acyclic.py``, which counts only
imports that run at import time, this counts *every* import statement in the
package at any depth: inside function and method bodies and inside
``if TYPE_CHECKING:`` blocks too. A deferred import would not break an import
today but would still put the dependency back.

Relative imports are resolved against each file's own dotted module name.
``from PyOMES import name`` and ``from .. import name`` name a sibling
subpackage, so ``name`` is the dependency. Imports inside docstring examples
(``>>>`` lines) are text, not statements, and are not counted.
"""

import ast
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PACKAGE = "PyOMES"
LAYER = "chemistry"
ALLOWED = {"chemistry", "units"}


def _module_name(path: Path, root: Path) -> str:
    parts = list(path.relative_to(root).with_suffix("").parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def _imported_subpackages(source: str, module: str, is_package: bool) -> set:
    """Names of the ``PyOMES`` subpackages/modules any import in *source* reaches.

    ``import PyOMES`` alone reports ``"PyOMES"`` (the root pulls in everything),
    and ``from PyOMES import *`` reports ``"*"``; neither is in ``ALLOWED``.
    """
    found = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            for alias in node.names:
                parts = alias.name.split(".")
                if parts[0] == PACKAGE:
                    found.add(parts[1] if len(parts) > 1 else PACKAGE)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                parts = module.split(".")
                base = parts if is_package else parts[:-1]
                if node.level > 1:
                    base = base[: len(base) - (node.level - 1)]
                target = base + (node.module.split(".") if node.module else [])
            else:
                target = node.module.split(".") if node.module else []
            if not target or target[0] != PACKAGE:
                continue
            if len(target) > 1:
                found.add(target[1])
            else:
                found.update(alias.name for alias in node.names)
    return found


def test_detector_catches_every_import_form():
    # Self-check against synthetic sources, so the empty result in the real
    # test below is meaningful and not an untested no-op.
    def found(src, module="PyOMES.chemistry.mod", is_package=False):
        return _imported_subpackages(src, module, is_package)

    # Module-level, relative and absolute.
    assert found("from ..reactions.equilibrium import X") == {"reactions"}
    assert found("from PyOMES.reactions import X") == {"reactions"}
    assert found("import PyOMES.reactions.equilibrium") == {"reactions"}
    # Deferred: function body and TYPE_CHECKING block.
    assert found("def f():\n    from PyOMES.reactions import X") == {"reactions"}
    assert found("def f():\n    from ..core import X") == {"core"}
    assert found("if TYPE_CHECKING:\n    from ..thermo import X") == {"thermo"}
    # Names that are themselves subpackages.
    assert found("from PyOMES import reactions") == {"reactions"}
    assert found("from .. import reactions, core") == {"reactions", "core"}
    assert found("import PyOMES") == {"PyOMES"}
    assert found("from PyOMES import *") == {"*"}
    # Package __init__ resolves one level higher than a plain module.
    assert found("from ..reactions import X", "PyOMES.chemistry", True) == {"reactions"}
    assert found("from .partition import X", "PyOMES.chemistry", True) == {"chemistry"}
    # Allowed and irrelevant imports.
    assert found("from ..units import R") == {"units"}
    assert found("from .species import S\nfrom . import common_species") == {"chemistry"}
    assert found("import math\nfrom dataclasses import dataclass\nimport numpy") == set()


def test_chemistry_imports_only_units():
    pkg_dir = REPO_ROOT / PACKAGE / LAYER
    files = sorted(p for p in pkg_dir.rglob("*.py") if "__pycache__" not in p.parts)
    assert (pkg_dir / "partition.py") in files, "layer directory not found or empty"

    offenders = {}
    for path in files:
        module = _module_name(path, REPO_ROOT)
        source = path.read_text(encoding="utf-8")
        extra = _imported_subpackages(source, module, path.name == "__init__.py") - ALLOWED
        if extra:
            offenders[path.relative_to(REPO_ROOT).as_posix()] = sorted(extra)

    assert not offenders, (
        f"{PACKAGE}/{LAYER}/ may import only {sorted(ALLOWED)} from {PACKAGE}, at any "
        "depth (function-level and TYPE_CHECKING imports included). Offenders: "
        + "; ".join(f"{f} -> {', '.join(names)}" for f, names in offenders.items())
    )
