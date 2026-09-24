"""Gas-phase non-ideality: equations of state.

- ``protocols.py``: the :class:`GasEOS` protocol.
- ``ideal.py``: :class:`IdealGasEOS` (Z = 1).
- ``peng_robinson.py``: :class:`PengRobinsonEOS`, with
  :class:`CriticalProperties` and the ``BIOGAS_SPECIES`` and ``BIOGAS_KIJ``
  parameter tables.

Nothing is re-exported here; the public names are exported from
:mod:`PyOMES.thermo`.
"""
