"""Liquid-phase non-ideality: activity-coefficient models and the water
properties they use.

- ``protocols.py``: the :class:`LiquidPhaseModel`, :class:`ActivityModel` and
  :class:`DifferentiableLiquidModel` protocols.
- ``ideal.py``: :class:`IdealLiquidModel` (γ_i = 1).
- ``davies.py``: :class:`DaviesLiquidModel` (Davies equation).
- ``sit.py``: :class:`SITLiquidModel` (Specific Ion Interaction Theory), with
  the ``SIT_EPSILON`` and ``ION_CHARGES`` tables.
- ``water_properties.py``: temperature-dependent water density, dielectric
  constant, Debye–Hückel ``A`` and the mol/L to mol/kg-water conversion.
- ``factory.py``: :func:`make_activity_model`, which builds a model from a
  ``(use_activity, activity_model)`` pair.

Nothing is re-exported here. The public names are exported from
:mod:`PyOMES.thermo`, except ``DifferentiableLiquidModel`` and
``water_kg_per_L``, which are imported from their own modules.
"""
