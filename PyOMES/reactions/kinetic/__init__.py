# -*- coding: utf-8 -*-
"""Kinetic reactions: reactions integrated through a rate law.

- ``reaction.py``: :class:`KineticReaction`, a validated stoichiometry
  paired with a callable rate law.
- ``rate_laws.py``: growth rate laws (:class:`Monod`, :class:`Contois`,
  :class:`Andrews` and others) that build rate functions for
  ``KineticReaction``.
- ``builder.py``: :class:`ReactionBuilder`, factories that derive
  ``KineticReaction`` stoichiometry from substrate and biomass formulas.

Nothing is re-exported here; the public names are exported from
:mod:`PyOMES.reactions`.
"""
