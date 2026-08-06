Legacy Pressure Control Removal Patch
======================================
Date: 2026-03-06

Summary
-------
This patch removes the legacy (pre-controller-API) pressure control mechanism
from the Fermenter simulation loop, ensuring that ALL pressure relief is now
handled exclusively through the control system API (PressureReliefController,
InstantPressureReliefController, or SmoothPressureReliefController).

The legacy mechanism was a set of module-level functions (_apply_pressure_control,
_vent_to_pressure_atm, _smooth_vent_fraction) that were called directly from within
the CO2 helper functions and the Stage 1 ODE RHS. The new mechanism delegates all
venting to _apply_valve_control_safely(), which calls the registered controller's
compute() and rates() methods.

Files modified
--------------
1. src/fermenter/fermenter_unit.py  (main changes)
2. scripts/test_run_bioSTEAM.py     (test script cleanup)

Phase 1 — P_relief_atm initialisation
--------------------------------------
- Added _infer_pressure_relief_setpoint() method to CUFermentationSpeciation.
  Searches registered controllers for P_set_atm attribute and falls back to
  P_init_atm + 0.2 if none found.
- P_relief_atm is now set in __init__ from the controller, eliminating the need
  for callers to set F.P_relief_atm = X manually.
- Removed manual F.P_relief_atm = 1.2 from test_run_bioSTEAM.py (previously
  inconsistent with PressureReliefController.P_set_atm=1.10).

Phase 2 — Disable legacy venting inside CO2 helpers
----------------------------------------------------
- Changed both calls to _co2_mass_transfer_step_with_relief() and
  _equilibrate_co2_with_relief[_CO2only]() to always pass
  pressure_control="none".
- The helpers now compute CO2 equilibration / mass transfer without performing
  any headspace venting internally (they return zero for all vent_* fields).
- The module-level functions (_apply_pressure_control etc.) are left in place
  but always short-circuit with mode='none'.

Phase 3 — Unconditional post-CO2 valve control
-----------------------------------------------
- Changed the guard for the post-CO2 _apply_valve_control_safely() block from:
    if pressure_control == 'valve' and control_system is not None
  to:
    if control_system is not None and self._control_has_tag('pressure_relief')
  This matches the pattern already used in the co2_mode='none' branch and makes
  all three CO2 mode paths consistent.
- Updated all diagnostics blocks:
    - Valve-open inference now uses self._control_has_tag('pressure_relief')
      instead of checking pressure_control == 'valve'.
    - Diagnostics 'mode' field set to 'control_api' instead of the legacy
      attribute value.

Phase 4 — Remove legacy pressure control from Stage 1 ODE RHS
--------------------------------------------------------------
- Removed the pc_mode / vent_k / vent_w variable setup that read from
  self.pressure_control.
- Removed two "valve-based pressure relief" blocks inside rhs_aug() that
  called self.control_system.step() within the ODE RHS (problematic due to
  side effects).
- Removed two "smooth pressure relief" blocks inside rhs_aug() that directly
  computed venting rates using self.P_relief_atm.
- Stage 1 pressure relief is now handled solely by the event_relief mechanism
  (scipy solve_ivp discontinuous event detection + solver restarts with
  proportionally scaled-down headspace moles).

What was NOT changed (and why)
------------------------------
- Module-level functions (_vent_to_pressure_atm, _smooth_vent_fraction,
  _apply_pressure_control): Left in place. After Phase 2, they always
  short-circuit with mode='none'. Safe to remove in a future cleanup.
- _control_pressure_floor_atm, _valve_substep_settings: Needed by the new
  controller path. Left as-is.
- event_relief mechanism in Stage 1: This is the correct way to handle
  discontinuous pressure events inside an ODE. Left as-is.
- InstantPressureReliefController, SmoothPressureReliefController in
  control/loops.py: New-style controller implementations. Left as-is.

Testing
-------
Run: python scripts/test_run_bioSTEAM.py

Expected changes vs prior behaviour:
- Pressure relief during CO2 equilibration/kinetic steps is now performed
  AFTER the CO2 step (via the controller API), not INSIDE the helper function.
  This may produce slightly different venting dynamics at the per-timestep
  level, but the overall physics is equivalent.
- P_relief_atm is now 1.10 atm (from PressureReliefController.P_set_atm)
  instead of the manually-set 1.20 atm. This is intentional — it makes the
  event_relief threshold consistent with the controller setpoint.
