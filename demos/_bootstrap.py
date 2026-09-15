"""Bootstrap: add the project models directory to sys.path.

Import this module at the top of any demo script that needs the
``vlmodels`` package (the concrete model implementations under
``models/``, e.g. ``vlmodels.adm1``). ``PyOMES`` itself is exposed via
the editable install (``pip install -e .``), which resolves directly
to the ``PyOMES/`` package at the repo root; we leave that path alone
here so the install's resolution wins. Demos that only need `PyOMES`
(e.g. the `builder/` demos, via `PyOMES.templates.stirred_tank`) don't
need this bootstrap at all.

    import _bootstrap  # noqa: F401
    from vlmodels.adm1.bsm2 import build_bsm2_cv
"""
import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parent.parent
_models = str(_repo_root / "models")
if _models not in sys.path:
    sys.path.insert(0, _models)
