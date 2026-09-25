# -*- coding: utf-8 -*-
"""Resolve an ``activity_model`` argument (a name or a model) to a model object."""

from __future__ import annotations

from typing import Union

from .protocols import ActivityModel
from .ideal import IdealLiquidModel
from .davies import DaviesLiquidModel
from .sit import SITLiquidModel

# Accepted names are the models' own ``name`` labels.
_MODELS_BY_NAME = {
    cls.name: cls for cls in (IdealLiquidModel, DaviesLiquidModel, SITLiquidModel)
}


def make_activity_model(activity_model: Union[str, ActivityModel] = "ideal") -> ActivityModel:
    """Return the activity model for *activity_model*.

    Parameters
    ----------
    activity_model : str or model object
        ``"ideal"``, ``"davies"`` or ``"sit"`` (case-insensitive) builds that
        model with its default settings. Any other object must be an activity
        model already (for example ``SITLiquidModel(epsilon=...)`` or a
        user-written model with a ``gamma(z, I_molL, *, T_K)`` method) and is
        returned unchanged.

    Raises
    ------
    ValueError
        If *activity_model* is a string that is not one of the accepted names.
    TypeError
        If *activity_model* is neither a string nor an object with a callable
        ``gamma`` method.
    """
    if not isinstance(activity_model, str):
        if not callable(getattr(activity_model, "gamma", None)):
            raise TypeError(
                f"activity_model must be a model name or an activity model with a "
                f"gamma() method; got {type(activity_model).__name__}."
            )
        return activity_model
    cls = _MODELS_BY_NAME.get(activity_model.strip().lower())
    if cls is None:
        raise ValueError(
            f"Unknown activity_model {activity_model!r}; "
            f"expected one of {sorted(_MODELS_BY_NAME)} or an activity model object."
        )
    return cls()
