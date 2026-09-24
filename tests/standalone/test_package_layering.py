"""Guards on which ``PyOMES`` packages and modules may import which.

- ``PyOMES/chemistry/`` depends on nothing in ``PyOMES`` except ``units``.
  ``chemistry/`` holds species declarations and the phase-partition protocols.
  It sits below ``reactions/``, ``core/`` and everything else, so no module in it
  may import from them.
- Inside ``PyOMES/reactions/``, the ``kinetic/`` and ``equilibrium/`` folders stay
  separate. Neither imports the other, ``blackbox.py`` imports neither, neither
  imports ``reaction_system`` except under ``if TYPE_CHECKING:``, and neither
  imports names from the ``PyOMES.reactions`` package root. They may import the
  shared top-level modules (``stoichiometry``, ``environment``, ``_shared``) and
  their own folder.

Unlike ``test_import_graph_acyclic.py``, which counts only imports that run at
import time, these count *every* import statement at any depth: inside function
and method bodies and inside ``if TYPE_CHECKING:`` blocks too. A deferred import
would not break an import today but would still put the dependency back. The
``reaction_system`` rule is the one place a ``TYPE_CHECKING`` import is allowed.

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

REACTIONS = f"{PACKAGE}.reactions"
KINETIC = f"{REACTIONS}.kinetic"
EQUILIBRIUM = f"{REACTIONS}.equilibrium"
REACTION_SYSTEM = f"{REACTIONS}.reaction_system"


def _module_name(path: Path, root: Path) -> str:
    parts = list(path.relative_to(root).with_suffix("").parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def _is_type_checking(test: ast.expr) -> bool:
    return (isinstance(test, ast.Name) and test.id == "TYPE_CHECKING") or (
        isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING"
    )


def _imports(source: str, module: str, is_package: bool) -> list:
    """Every import statement in *source*, at any depth.

    Returns ``(from_module, name, type_checking_only)`` tuples: ``import a.b``
    gives ``("a.b", None, ...)`` and ``from m import n`` gives ``(m, "n", ...)``,
    with relative ``m`` resolved to an absolute dotted name. Imports in the body
    of an ``if TYPE_CHECKING:`` block are flagged; its ``else`` branch is not.
    """
    found = []

    def visit(node, tc):
        if isinstance(node, ast.Import):
            found.extend((alias.name, None, tc) for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                parts = module.split(".")
                base = parts if is_package else parts[:-1]
                if node.level > 1:
                    base = base[: len(base) - (node.level - 1)]
                target = ".".join(base + ([node.module] if node.module else []))
            else:
                target = node.module or ""
            found.extend((target, alias.name, tc) for alias in node.names)
        if isinstance(node, ast.If) and _is_type_checking(node.test):
            for child in node.body:
                visit(child, True)
            for child in node.orelse:
                visit(child, tc)
            return
        for child in ast.iter_child_nodes(node):
            visit(child, tc)

    visit(ast.parse(source), False)
    return found


def _full_target(from_module: str, name) -> str:
    return f"{from_module}.{name}" if name else from_module


def _imported_subpackages(source: str, module: str, is_package: bool) -> set:
    """Names of the ``PyOMES`` subpackages/modules any import in *source* reaches.

    ``import PyOMES`` alone reports ``"PyOMES"`` (the root pulls in everything),
    and ``from PyOMES import *`` reports ``"*"``; neither is in ``ALLOWED``.
    """
    found = set()
    for from_module, name, _ in _imports(source, module, is_package):
        parts = _full_target(from_module, name).split(".")
        if parts[0] == PACKAGE:
            found.add(parts[1] if len(parts) > 1 else PACKAGE)
    return found


def _within(target: str, prefix: str) -> bool:
    return target == prefix or target.startswith(prefix + ".")


def _reactions_layout_violations(source: str, module: str, is_package: bool) -> list:
    """Imports in a ``reactions/`` module that cross the kinetic/equilibrium seam."""
    if _within(module, KINETIC):
        forbidden = [EQUILIBRIUM]
    elif _within(module, EQUILIBRIUM):
        forbidden = [KINETIC]
    elif module == f"{REACTIONS}.blackbox":
        forbidden = [KINETIC, EQUILIBRIUM]
    else:
        return []
    in_folder = not module.endswith(".blackbox")
    out = []
    for from_module, name, tc in _imports(source, module, is_package):
        target = _full_target(from_module, name)
        stmt = f"from {from_module} import {name}" if name else f"import {from_module}"
        if any(_within(target, f) for f in forbidden):
            out.append(f"{stmt} (crosses the kinetic/equilibrium seam)")
        elif in_folder and _within(target, REACTION_SYSTEM) and not tc:
            out.append(f"{stmt} (reaction_system is allowed only under TYPE_CHECKING)")
        elif in_folder and from_module == REACTIONS and name is not None:
            out.append(f"{stmt} (names from the package root; import the defining module)")
    return out


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


def test_reactions_layout_detector_catches_every_crossing():
    # Self-check for the reactions/ layout rule, as above.
    def bad(src, module, is_package=False):
        return len(_reactions_layout_violations(src, module, is_package))

    builder = "PyOMES.reactions.kinetic.builder"
    plots = "PyOMES.reactions.equilibrium.plots"
    # Crossing the seam, in every form and at any depth.
    assert bad("from PyOMES.reactions.equilibrium.constraint import X", builder) == 1
    assert bad("from ..equilibrium.reaction import X", builder) == 1
    assert bad("from ..equilibrium import reaction", builder) == 1
    assert bad("from .. import equilibrium", builder) == 1
    assert bad("import PyOMES.reactions.equilibrium.reaction", builder) == 1
    assert bad("def f():\n    from ..equilibrium.reaction import X", builder) == 1
    assert bad("if TYPE_CHECKING:\n    from ..equilibrium.reaction import X", builder) == 1
    assert bad("from ..kinetic.reaction import X", plots) == 1
    assert bad("from ..kinetic import rate_laws", "PyOMES.reactions.equilibrium", True) == 1
    assert bad("from .kinetic.reaction import X", "PyOMES.reactions.blackbox") == 1
    assert bad("from .equilibrium.constraint import X", "PyOMES.reactions.blackbox") == 1
    # reaction_system only under TYPE_CHECKING, and not in its else branch.
    assert bad("if TYPE_CHECKING:\n    from PyOMES.reactions.reaction_system import R", plots) == 0
    assert bad("from PyOMES.reactions.reaction_system import R", plots) == 1
    assert bad("def f():\n    from ..reaction_system import R", plots) == 1
    assert bad("if TYPE_CHECKING:\n    pass\nelse:\n    from ..reaction_system import R", plots) == 1
    # Names from the package root.
    assert bad("from PyOMES.reactions import StoichiometryEntry", builder) == 1
    assert bad("from .. import stoichiometry", builder) == 1
    # Allowed: shared top-level modules, the own folder, anything outside reactions/.
    assert bad("from PyOMES.reactions.stoichiometry import S", builder) == 0
    assert bad("from PyOMES.reactions.environment import E\nfrom PyOMES.reactions._shared import f", builder) == 0
    assert bad("from .rate_laws import Monod\nfrom .reaction import K", builder) == 0
    assert bad("from .constraint import vant_hoff_log_K", "PyOMES.reactions.equilibrium.interphase") == 0
    assert bad("from PyOMES.chemistry.species import S\nfrom PyOMES.units import R\nimport math", plots) == 0
    assert bad("from .environment import E", "PyOMES.reactions.blackbox") == 0
    # Files outside the three places the rule covers are not checked.
    assert bad("from .kinetic.reaction import X\nfrom .equilibrium.reaction import Y",
               "PyOMES.reactions.reaction_system") == 0


def test_reactions_kinetic_and_equilibrium_stay_separate():
    reactions_dir = REPO_ROOT / PACKAGE / "reactions"
    files = sorted(
        p for folder in ("kinetic", "equilibrium")
        for p in (reactions_dir / folder).rglob("*.py") if "__pycache__" not in p.parts
    ) + [reactions_dir / "blackbox.py"]
    for expected in ("kinetic/reaction.py", "equilibrium/constraint.py", "blackbox.py"):
        assert reactions_dir / expected in files, f"reactions/{expected} not found"

    offenders = {}
    for path in files:
        module = _module_name(path, REPO_ROOT)
        source = path.read_text(encoding="utf-8")
        problems = _reactions_layout_violations(source, module, path.name == "__init__.py")
        if problems:
            offenders[path.relative_to(REPO_ROOT).as_posix()] = problems

    assert not offenders, (
        f"{PACKAGE}/reactions/kinetic/ and equilibrium/ must not import each other, "
        "blackbox.py must import neither, and neither folder may import reaction_system "
        "outside TYPE_CHECKING or import names from the package root (all import depths "
        "counted). Offenders: "
        + "; ".join(f"{f} -> {', '.join(p)}" for f, p in offenders.items())
    )
