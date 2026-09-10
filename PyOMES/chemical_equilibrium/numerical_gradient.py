# -*- coding: utf-8 -*-
"""NumericalGradientEquilibriumEngine — gray-box wrapper via central finite differences.

Wraps any :class:`~PyOMES.chemical_equilibrium.protocols.ChemicalEquilibriumEngineProtocol` and
satisfies :class:`~PyOMES.chemical_equilibrium.protocols.GrayBoxEngineProtocol` by
approximating :math:`\\partial z / \\partial y` numerically.

Per-component perturbation magnitude::

    eps_j = max(eps_abs, eps_rel * |T_j|)

Central difference::

    ∂c_i/∂T_j ≈ [h(T + eps_j·eⱼ) − h(T − eps_j·eⱼ)] / (2·eps_j)

When ``T_j < eps_abs`` (near-zero total), the minus perturbation is clamped at
zero and a one-sided difference is computed automatically.

Requires ``2 × n_components`` inner ``solve()`` calls per Jacobian evaluation.

.. note::
    The wrapped engine must accept ``totals={component_id: mol_L, ...}`` as a
    ``solve()`` kwarg.  :class:`~PyOMES.chemical_equilibrium.nr_engine.NRChemicalEquilibriumEngine`
    supports this natively; :class:`~PyOMES.chemical_equilibrium.engine.BisectionChemicalEquilibriumEngine`
    (pKa-ladder) does not.

.. note::
    This engine cannot satisfy
    :class:`~PyOMES.chemical_equilibrium.protocols.ResidualCapable` or
    :class:`~PyOMES.chemical_equilibrium.protocols.SplitJacobianCapable`; it has no access
    to the algebraic residual ``g``, only to the solution map ``h``.
"""
from __future__ import annotations

from typing import Any, Dict, FrozenSet

import numpy as np

from .protocols import (
    EquilibriumResult,
    GrayBoxEngineProtocol,
    ChemicalEquilibriumEngineProtocol,
    SpeciationJacobian,
)


class NumericalGradientEquilibriumEngine:
    """Gray-box speciation engine wrapper using central finite differences.

    Parameters
    ----------
    engine : ChemicalEquilibriumEngineProtocol
        Inner engine to wrap.  Must accept ``totals=`` in ``solve()``.
    eps_abs : float
        Absolute perturbation floor (mol/L).  Dominates near zero totals.
        Default ``1e-10``.
    eps_rel : float
        Relative perturbation scale (dimensionless).  Dominates in the
        millimolar range.  Default ``1e-5``.
        Set to ``0`` for pure absolute mode; set ``eps_abs=0`` for pure
        relative mode.
    """

    def __init__(
        self,
        engine: ChemicalEquilibriumEngineProtocol,
        *,
        eps_abs: float = 1e-10,
        eps_rel: float = 1e-5,
    ) -> None:
        if not isinstance(engine, ChemicalEquilibriumEngineProtocol):
            raise TypeError(
                f"NumericalGradientEquilibriumEngine requires a ChemicalEquilibriumEngineProtocol; "
                f"got {type(engine).__name__!r}."
            )
        self._engine = engine
        self.eps_abs = float(eps_abs)
        self.eps_rel = float(eps_rel)

    # ------------------------------------------------------------------
    # ChemicalEquilibriumEngineProtocol delegation
    # ------------------------------------------------------------------

    @property
    def n_solve_calls(self) -> int:
        """Total solve calls made by the inner engine (includes perturbation calls)."""
        return self._engine.n_solve_calls

    def solve(self, **kwargs: Any) -> EquilibriumResult:
        return self._engine.solve(**kwargs)

    def algebraic_species(self) -> FrozenSet[str]:
        return self._engine.algebraic_species()

    def reset_cache(self) -> None:
        self._engine.reset_cache()

    def reset_counters(self) -> None:
        self._engine.reset_counters()

    # ------------------------------------------------------------------
    # SolutionJacobianCapable: gray-box addition
    # ------------------------------------------------------------------

    def jacobian_dz_dy(self, **kwargs: Any) -> SpeciationJacobian:
        """Compute :math:`\\partial z / \\partial y` by central finite differences.

        Parameters
        ----------
        totals : dict
            ``{component_id: C_total_mol_L}`` — the point at which to evaluate
            the Jacobian.  Must be provided explicitly; ``phases=`` is not
            supported here because the component-extraction logic is
            engine-specific.
        **kwargs
            Additional kwargs forwarded to the inner ``solve()`` (e.g.
            ``strong_ions=``, ``T_K=``).

        Returns
        -------
        SpeciationJacobian
            ``dz_dy[i, j] = ∂c_i / ∂T_j`` where rows are algebraic species
            (sorted alphabetically) and columns are component totals in the
            order supplied by ``totals``.

        Raises
        ------
        ValueError
            If ``totals=`` is not provided.
        """
        totals = kwargs.pop("totals", None)
        if totals is None:
            if "phases" in kwargs:
                raise ValueError(
                    "NumericalGradientEquilibriumEngine.jacobian_dz_dy does not support "
                    "phases=. Pass totals={component_id: mol_L, ...} explicitly."
                )
            raise ValueError(
                "NumericalGradientEquilibriumEngine.jacobian_dz_dy requires totals= kwarg "
                "(dict mapping component ID → total concentration in mol/L)."
            )

        component_ids = tuple(totals.keys())
        T = np.array([float(totals[c]) for c in component_ids])

        alg_ids = tuple(sorted(self._engine.algebraic_species()))
        n_alg = len(alg_ids)
        n_comp = len(component_ids)
        dz_dy = np.zeros((n_alg, n_comp))

        for j, (comp_id, T_j) in enumerate(zip(component_ids, T)):
            eps_j = max(self.eps_abs, self.eps_rel * abs(T_j))

            T_plus = T_j + eps_j
            T_minus = max(T_j - eps_j, 0.0)
            denom = T_plus - T_minus
            if denom == 0.0:
                continue

            out_plus = self._engine.solve(
                totals={**totals, comp_id: T_plus}, **kwargs
            )
            out_minus = self._engine.solve(
                totals={**totals, comp_id: T_minus}, **kwargs
            )

            for i, sp_id in enumerate(alg_ids):
                v_plus = out_plus.species_mol_L.get(sp_id)
                v_minus = out_minus.species_mol_L.get(sp_id)
                if v_plus is None or v_minus is None:
                    continue
                try:
                    dz_dy[i, j] = (float(v_plus) - float(v_minus)) / denom
                except (TypeError, ValueError):
                    pass

        return SpeciationJacobian(
            dz_dy=dz_dy,
            algebraic_ids=alg_ids,
            component_ids=component_ids,
        )
