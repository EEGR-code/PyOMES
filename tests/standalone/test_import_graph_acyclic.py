"""Guard: the top-level module import graph of ``PyOMES/`` has no cycles.

Only imports that actually run at import time count: module- and
class-body-level ``import`` / ``from ... import`` statements. Imports inside a
function or method body (deferred until the function is called) and imports
inside ``if TYPE_CHECKING:`` (never executed) are excluded, because neither
can produce a real ``ImportError`` at import time. The one real cycle in the
repo, ``chemistry.equilibria`` <-> ``chemistry.thermo_params``, is deliberately
made of two such lazy imports and is not caught here; see the cleanup plan's
audit.

Relative imports are resolved against each file's own dotted module name.
``from pkg import name`` is resolved to the submodule ``pkg.name`` when that
submodule exists, and to ``pkg`` itself otherwise (covering both "import a
submodule" and "import a name defined in ``pkg/__init__.py``").
"""

import ast
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PACKAGE = "PyOMES"


def _module_name(path: Path, root: Path) -> str:
    parts = list(path.relative_to(root).with_suffix("").parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def _is_type_checking(test: ast.expr) -> bool:
    return (isinstance(test, ast.Name) and test.id == "TYPE_CHECKING") or (
        isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING"
    )


class _TopLevelImportVisitor(ast.NodeVisitor):
    """Collects ``(target_module, names_or_None)`` for import-time imports only."""

    def __init__(self, module: str, is_package: bool):
        self.module = module
        self.is_package = is_package
        self._depth = 0
        self.imports = []

    def visit_FunctionDef(self, node):
        self._depth += 1
        self.generic_visit(node)
        self._depth -= 1

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_If(self, node):
        if _is_type_checking(node.test):
            # The TYPE_CHECKING body never runs; the else branch does.
            for stmt in node.orelse:
                self.visit(stmt)
            return
        self.generic_visit(node)

    def visit_Import(self, node):
        if self._depth:
            return
        for alias in node.names:
            self.imports.append((alias.name, None))

    def visit_ImportFrom(self, node):
        if self._depth:
            return
        if node.level:
            parts = self.module.split(".")
            base = parts if self.is_package else parts[:-1]
            if node.level > 1:
                base = base[: len(base) - (node.level - 1)]
            target = ".".join(base + ([node.module] if node.module else []))
        else:
            target = node.module
        self.imports.append((target, [a.name for a in node.names]))


def _resolve(target: str, names, known: set) -> list:
    if names is None:
        return [target] if target in known else []
    candidates = [f"{target}.{n}" for n in names if f"{target}.{n}" in known]
    if candidates:
        return candidates
    return [target] if target in known else []


def _build_top_level_graph(root: Path, package: str) -> dict:
    pkg_dir = root / package
    files = sorted(p for p in pkg_dir.rglob("*.py") if "__pycache__" not in p.parts)
    known = {_module_name(p, root) for p in files}
    graph = {m: set() for m in known}
    for path in files:
        module = _module_name(path, root)
        visitor = _TopLevelImportVisitor(module, path.name == "__init__.py")
        visitor.visit(ast.parse(path.read_text(encoding="utf-8")))
        for target, names in visitor.imports:
            if target is None or not (target == package or target.startswith(package + ".")):
                continue
            for resolved in _resolve(target, names, known):
                if resolved != module:
                    graph[module].add(resolved)
    return graph


def _find_cycle(graph: dict):
    """Return a list of modules forming a cycle, or None if the graph is acyclic."""
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {node: WHITE for node in graph}
    path = []
    cycle = []

    def visit(node):
        color[node] = GRAY
        path.append(node)
        for nxt in sorted(graph.get(node, ())):
            if nxt not in color:
                continue
            if color[nxt] == GRAY:
                cycle.extend(path[path.index(nxt):] + [nxt])
                return True
            if color[nxt] == WHITE and visit(nxt):
                return True
        path.pop()
        color[node] = BLACK
        return False

    for node in sorted(graph):
        if color[node] == WHITE and visit(node):
            return cycle
    return None


def test_top_level_import_graph_is_acyclic():
    # Self-check: the cycle detector against tiny synthetic graphs, so an
    # empty result below is meaningful and not just an untested no-op.
    assert _find_cycle({"a": {"b"}, "b": {"c"}, "c": {"a"}}) == ["a", "b", "c", "a"]
    assert _find_cycle({"a": {"b"}, "b": {"c"}, "c": set()}) is None

    graph = _build_top_level_graph(REPO_ROOT, PACKAGE)
    cycle = _find_cycle(graph)
    assert cycle is None, (
        "Top-level import cycle in PyOMES/ (function-level and TYPE_CHECKING "
        "imports don't count): " + " -> ".join(cycle)
    )
