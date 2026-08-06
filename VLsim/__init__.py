"""Compatibility shim for the former VLsim import name."""

from __future__ import annotations

import PyOMES as _pyomes
from PyOMES import *  # noqa: F401,F403

__all__ = getattr(_pyomes, "__all__", [])
__version__ = getattr(_pyomes, "__version__", None)
__path__ = list(getattr(_pyomes, "__path__", []))
