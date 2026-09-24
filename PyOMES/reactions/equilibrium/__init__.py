# -*- coding: utf-8 -*-
"""Equilibrium reactions: algebraic constraints with an equilibrium constant.

- ``reaction.py``: :class:`EquilibriumReaction`, the general mass-action
  equilibrium, together with the :class:`EquilibriumConstraint` protocol,
  :func:`vant_hoff_log_K` and :func:`classify_equilibrium_constraint`.
- ``interphase.py``: :class:`HenryEquilibrium`, :class:`RaoultEquilibrium`
  and :class:`KspEquilibrium`, named physical-law constraints relating one
  species across two phases, which double as ``PartitionModel`` instances.
- ``plots.py``: Van 't Hoff and speciation plots (matplotlib, imported
  lazily).

Nothing is re-exported here; the public names are exported from
:mod:`PyOMES.reactions`.
"""
