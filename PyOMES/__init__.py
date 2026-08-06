"""Local source-tree shim for importing PyOMES without installation."""

from __future__ import annotations

from pathlib import Path

_SRC_ROOT = Path(__file__).resolve().parent.parent / "src"
__path__ = [str(_SRC_ROOT)]
__file__ = str(_SRC_ROOT / "__init__.py")

with open(__file__, "rb") as _source:
    exec(compile(_source.read(), __file__, "exec"), globals())
