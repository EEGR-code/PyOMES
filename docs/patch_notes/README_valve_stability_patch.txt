Valve stability patch

Changes in fermenter_unit.py:
- Added _control_pressure_floor_atm(unit) to derive a conservative explicit-step pressure floor
  from the controller set pressure and actuator back pressure.
- Added _apply_valve_control_safely(...) to apply valve venting with:
  * inventory limiting (cannot vent more gas than exists)
  * pressure-floor limiting (cannot vent below the explicit-step floor)
  * internal substepping for the valve update
- Updated split-solver valve application in the CO2 kinetic and equilibrium branches to use
  the safe helper and to log the resulting pre/post-valve pressures.
- Updated the no-CO2 split branch valve application similarly.

Intent:
- prevent catastrophic over-venting to ~0 atm when a large valve diameter is used with an explicit timestep
- keep diagnostics physically interpretable
