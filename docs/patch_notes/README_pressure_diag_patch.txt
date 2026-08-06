Pressure diagnostics patch

Changes in fermenter_unit.py:
- Added explicit pre- and post-valve pressure logging:
  * P_hs_pre_valve_atm
  * P_hs_post_valve_atm
- Added vent_total_mol and diagnostic stage metadata:
  * diag_stage
  * diag_stage_order
- plot_diagnostics() now sorts diagnostics by time/stage before plotting.
- plot_diagnostics() now plots pre-valve and post-valve pressure when available.

Purpose:
This avoids mixing controller-input pressure with post-vent pressure in the same trace,
which previously made the valve diagnostics hard to interpret.
