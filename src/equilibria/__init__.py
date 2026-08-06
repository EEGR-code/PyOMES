"""Gas-phase / VLE equilibrium utilities (equations of state, Henry partitioning).

The legacy coupled speciation+VLE fermenter bridge (``ProcessCoupledEquilibrator``,
``CoupledEquilibriumFactory``, ``EquilibratorResult``) was removed as dead code
left over from CUFERMENTER_SUNSET (see ``legacy-cleanup-shipped``) — it had no
live callers and depended on the deleted ``Fermenter``/``fermenter_unit.py``
object shape. Use ``PyOMES.core.transfer_models.EquilibriumTransferModel`` for
gas-liquid equilibrium coupling in current code.
"""
