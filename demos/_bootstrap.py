"""Bootstrap: add the project models directory to sys.path.

Import this module at the top of any demo script to make the
``vlmodels`` package importable. ``PyOMES`` itself is exposed via the
editable install (``pip install -e .``), which resolves directly to
the ``PyOMES/`` package at the repo root; we leave that path alone
here so the install's resolution wins.

    import _bootstrap  # noqa: F401
    from PyOMES.config import FermenterBuilder
    from vlmodels.fermenter.types import FermenterState
"""
import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parent.parent
_models = str(_repo_root / "models")
if _models not in sys.path:
    sys.path.insert(0, _models)
