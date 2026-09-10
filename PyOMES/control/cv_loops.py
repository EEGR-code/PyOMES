# -*- coding: utf-8 -*-
"""CV-native concrete controllers.

Pattern 1 controllers: each controller reads a
:class:`~PyOMES.core.snapshot.CVSnapshot` (or
:class:`~PyOMES.core.snapshot.SimulationSnapshot`) and returns a
:class:`~PyOMES.control.actions.ControlAction` from a single
``compute(state, dt_h)`` method. Controller-internal state
(setpoint, integral, sampling-hold buffer) lives on ``self``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Union

from PyOMES.control.actions import ControlAction
from PyOMES.core.snapshot import CVSnapshot, SimulationSnapshot


# ════════════════════════════════════════════════════════════════════════
#  pH controller — Pattern 1 / CV-native
# ════════════════════════════════════════════════════════════════════════


@dataclass
class PHController:
    """pH controller via acid/base dosing — CV-native port.

    Reads ``snapshot.pH``, ``snapshot.V_liq_L``, and the optional
    ``snapshot.sensors["CT_P_mol_L"]`` phosphate cap hint. Returns
    a :class:`ControlAction` with the acid- or base-dosing flux
    written into ``flux_applied["liquid"][chemical_id]`` (mol/h)
    plus the same total in ``dosed_mol`` for audit.

    Behaviour and tuning surface mirror the legacy
    :class:`PyOMES.control.loops.PHController` 1:1, plus the
    Pattern 1 collapse:

    * No ``Commands`` mutation; controller owns its integral and
      mode state on ``self``.
    * No separate ``Actuator`` step; ``compute()`` builds the
      full action in one shot.

    Attributes
    ----------
    setpoint : float
        Target pH.
    Kp, Ki : float
        Proportional and integral gains on the per-volume add rate
        (mol/L/h).
    max_add_molL_hr : float
        Cap on the per-volume add rate (mol/L/h). Historical name
        kept for compatibility.
    CT_P_max : float
        Phosphate cap hint (mol/L). If ``snapshot.sensors["CT_P_mol_L"]``
        is present, the dosing rate is clamped so the phosphate
        total doesn't exceed this.
    deadband : float
        Dead-band around the setpoint. ``|err| <= deadband`` → no
        dosing.
    chemical_id : str
        Acid compound id (must validate against
        :mod:`PyOMES.chemistry.registry`).
    base_chemical_id : Optional[str]
        Base compound id. If ``None``, base dosing is disabled
        (acid-only legacy mode).
    sample_period_h, sample_period_s : Optional[float]
        Optional sampling period. When set, the orchestrator only
        invokes ``compute()`` at sampling instants; between
        samples the held :class:`ControlAction` is reused
        (zero-order hold).
    target_cv_key : Optional[str]
        Optional multi-CV target. When set, the orchestrator
        passes ``sim_snap.cvs[target_cv_key]`` as the state.
    control_tags : tuple of str
        Categorisation tags propagated for downstream tooling.
    enable_diagnostics : bool
        If True, ``self.diag`` accumulates a per-call time series.
    """

    setpoint: float
    Kp: float = 0.5
    Ki: float = 0.0
    max_add_molL_hr: float = 0.2
    CT_P_max: float = 2.0
    deadband: float = 0.0

    chemical_id: str = "H3PO4"
    base_chemical_id: Optional[str] = None

    sample_period_s: Optional[float] = None
    sample_period_h: Optional[float] = None
    target_cv_key: Optional[str] = None
    label: str = "pH"
    control_tags: tuple = field(
        default_factory=lambda: ("pH", "liquid_addition"),
    )
    enable_diagnostics: bool = False

    # Controller-internal state (carried on self for Pattern 1).
    _I_err: float = field(default=0.0, init=False, repr=False)
    _last_mode: str = field(default="", init=False, repr=False)
    diag: Dict[str, List[float]] = field(
        default_factory=dict, init=False, repr=False,
    )

    def __post_init__(self):
        try:
            from ..chemistry.registry import validate_compound_id
        except ImportError:
            return  # registry not available; skip validation
        validate_compound_id(
            self.chemical_id,
            context=f"PHController acid (chemical_id={self.chemical_id!r})",
        )
        if self.base_chemical_id is not None:
            validate_compound_id(
                self.base_chemical_id,
                context=(
                    f"PHController base "
                    f"(base_chemical_id={self.base_chemical_id!r})"
                ),
            )

    def reset(self) -> None:
        """Clear integral state and diagnostics.

        Called by :meth:`Simulation.run` at the top of every run
        (decision 10: always destructive).
        """
        self._I_err = 0.0
        self._last_mode = ""
        if self.enable_diagnostics:
            self._reset_diagnostics()

    def _reset_diagnostics(self) -> None:
        self.diag = {
            "t_h": [],
            "pH": [],
            "setpoint": [],
            "err": [],
            "err_eff": [],
            "mode": [],
            "acid_mol_h": [],
            "base_mol_h": [],
            "add_rate_mol_L_h": [],
            "V_liq_L": [],
            "I_err": [],
        }

    def _diag_append(
        self,
        state: CVSnapshot,
        err: float,
        err_eff: float,
        mode: str,
        acid_mol_h: float = 0.0,
        base_mol_h: float = 0.0,
        add_rate_mol_L_h: Optional[float] = None,
    ) -> None:
        if not self.enable_diagnostics:
            return
        if not self.diag:
            self._reset_diagnostics()
        self.diag["t_h"].append(float(state.t_h))
        pH_val = state.pH if state.pH is not None else float("nan")
        self.diag["pH"].append(float(pH_val))
        self.diag["setpoint"].append(float(self.setpoint))
        self.diag["err"].append(float(err))
        self.diag["err_eff"].append(float(err_eff))
        self.diag["mode"].append(str(mode))
        self.diag["acid_mol_h"].append(float(acid_mol_h))
        self.diag["base_mol_h"].append(float(base_mol_h))
        self.diag["add_rate_mol_L_h"].append(
            float(add_rate_mol_L_h)
            if add_rate_mol_L_h is not None
            else float("nan")
        )
        self.diag["V_liq_L"].append(float(state.V_liq_L))
        self.diag["I_err"].append(float(self._I_err))

    def compute(
        self,
        state: Union[CVSnapshot, SimulationSnapshot],
        dt_h: float,
    ) -> ControlAction:
        """Compute one step's action from the snapshot.

        Returns a :class:`ControlAction` with ``flux_applied``
        populated when dosing fires; empty otherwise. ``dosed_mol``
        is populated with the integrated total over ``dt_h`` for
        audit.
        """
        # Resolve which CV to read (multi-CV controllers can set
        # target_cv_key; otherwise the orchestrator should have
        # passed a single-CV snapshot).
        cv_snap = self._resolve_cv(state)
        cv_key = cv_snap.cv_key
        t_h = cv_snap.t_h
        V_liq = max(float(cv_snap.V_liq_L), 0.0)

        # Empty-action fallback.
        empty = ControlAction(
            controller_label=self.label,
            target_cv_key=cv_key,
            t_h=t_h,
            dt_h=float(dt_h),
        )

        pH = cv_snap.pH
        if pH is None:
            self._diag_append(
                cv_snap, err=0.0, err_eff=0.0, mode="no_pH",
            )
            return empty

        err = float(pH) - float(self.setpoint)
        db = float(self.deadband)
        err_eff = 0.0 if abs(err) <= db else err

        # Too basic → dose acid.
        if err_eff > 0.0:
            if self._last_mode != "acid":
                self._I_err = 0.0
                self._last_mode = "acid"
            self._I_err += err_eff * max(float(dt_h), 0.0)
            add_rate = float(self.Kp) * err_eff + float(self.Ki) * self._I_err
            add_rate = max(0.0, min(float(self.max_add_molL_hr), add_rate))
            mol_h = add_rate * V_liq

            # Optional phosphate cap.
            ct_p = cv_snap.sensors.get("CT_P_mol_L")
            if (
                ct_p is not None
                and float(self.CT_P_max) > 0.0
                and float(dt_h) > 0.0
            ):
                max_add_mol = max(
                    0.0, (float(self.CT_P_max) - float(ct_p)) * V_liq,
                )
                mol_h = min(mol_h, max_add_mol / max(float(dt_h), 1e-30))

            mol_h = max(0.0, float(mol_h))
            self._diag_append(
                cv_snap, err=err, err_eff=err_eff, mode="acid",
                acid_mol_h=mol_h, add_rate_mol_L_h=add_rate,
            )
            return ControlAction(
                controller_label=self.label,
                target_cv_key=cv_key,
                t_h=t_h,
                dt_h=float(dt_h),
                flux_applied={"liquid": {self.chemical_id: mol_h}},
                dosed_mol={self.chemical_id: mol_h * float(dt_h)},
            )

        # Too acidic → dose base (when configured).
        if err_eff < 0.0 and self.base_chemical_id is not None:
            if self._last_mode != "base":
                self._I_err = 0.0
                self._last_mode = "base"
            err_b = -err_eff
            self._I_err += err_b * max(float(dt_h), 0.0)
            add_rate = float(self.Kp) * err_b + float(self.Ki) * self._I_err
            add_rate = max(0.0, min(float(self.max_add_molL_hr), add_rate))
            mol_h = max(0.0, float(add_rate * V_liq))
            self._diag_append(
                cv_snap, err=err, err_eff=err_eff, mode="base",
                base_mol_h=mol_h, add_rate_mol_L_h=add_rate,
            )
            return ControlAction(
                controller_label=self.label,
                target_cv_key=cv_key,
                t_h=t_h,
                dt_h=float(dt_h),
                flux_applied={"liquid": {self.base_chemical_id: mol_h}},
                dosed_mol={self.base_chemical_id: mol_h * float(dt_h)},
            )

        # Dead-band / no-action.
        self._diag_append(
            cv_snap, err=err, err_eff=err_eff, mode="none",
        )
        return empty

    def _resolve_cv(
        self, state: Union[CVSnapshot, SimulationSnapshot],
    ) -> CVSnapshot:
        """Resolve a ``CVSnapshot`` from either a single-CV or
        multi-CV input."""
        if isinstance(state, CVSnapshot):
            return state
        # SimulationSnapshot — pick the targeted CV.
        if self.target_cv_key is not None:
            if self.target_cv_key not in state.cvs:
                raise KeyError(
                    f"PHController(target_cv_key={self.target_cv_key!r}) "
                    f"not found in SimulationSnapshot.cvs: "
                    f"{sorted(state.cvs.keys())}"
                )
            return state.cvs[self.target_cv_key]
        # No target — fall back to the first CV (single-CV ensembles).
        if len(state.cvs) != 1:
            raise RuntimeError(
                f"PHController received a SimulationSnapshot with "
                f"{len(state.cvs)} CVs but has no target_cv_key set. "
                f"Set target_cv_key to disambiguate."
            )
        return next(iter(state.cvs.values()))

    # ── Diagnostics surface (matplotlib helper preserved) ─────────────

    def plot_diagnostics(self, show: bool = True):
        """Plot pH-loop diagnostics — pH/setpoint trace and dosing
        rates. Mirrors the legacy plot_diagnostics method.

        Returns ``(fig1, fig2)`` or ``None`` if diagnostics were
        never enabled. Raises :class:`ImportError` if matplotlib
        is missing.
        """
        if not self.diag:
            return None
        try:
            import matplotlib.pyplot as plt  # type: ignore
        except Exception as e:
            raise ImportError(
                "matplotlib is required for PHController.plot_diagnostics()"
            ) from e
        t = self.diag.get("t_h", [])
        pH = self.diag.get("pH", [])
        sp = self.diag.get("setpoint", [])
        acid = self.diag.get("acid_mol_h", [])
        base = self.diag.get("base_mol_h", [])
        fig1 = plt.figure()
        ax1 = fig1.gca()
        ax1.plot(t, pH, label="pH")
        ax1.plot(t, sp, label="setpoint")
        ax1.set_xlabel("Time (h)")
        ax1.set_ylabel("pH")
        ax1.legend()
        fig2 = plt.figure()
        ax2 = fig2.gca()
        ax2.plot(t, acid, label="acid (mol/h)")
        ax2.plot(t, base, label="base (mol/h)")
        ax2.set_xlabel("Time (h)")
        ax2.set_ylabel("Dosing rate (mol/h)")
        ax2.legend()
        if show:
            plt.show()
        return fig1, fig2


# ════════════════════════════════════════════════════════════════════════
#  Dissolved-oxygen controller (agitation → kLa) — Pattern 1 / CV-native
# ════════════════════════════════════════════════════════════════════════


def _default_power_law_kLa(
    rpm: float, *, kLa_ref: float, rpm_ref: float, exponent: float,
) -> float:
    """Default agitation→kLa correlation: kLa = kLa_ref * (rpm/rpm_ref)^exponent."""
    if rpm_ref <= 0.0 or rpm <= 0.0:
        return 0.0
    return float(kLa_ref) * (float(rpm) / float(rpm_ref)) ** float(exponent)


@dataclass
class DOAgitationController:
    """DO controller via agitation speed — CV-native port.

    PI loop on dissolved-oxygen error → agitation RPM → kLa(O2)
    via an equipment-specific correlation. kLa(CO2) is coupled
    via a fixed diffusivity ratio (default 0.9, ≈√(D_CO2/D_O2)).

    Reads ``snapshot.sensors["DO_mol_L"]``. Returns a
    :class:`ControlAction` with ``params_changed`` populated:

    * ``"kLa.O2"``: float — new kLa_O2 (1/h)
    * ``"kLa.CO2"``: float — new kLa_CO2 (1/h)

    The orchestrator's C9 apply path dispatches these onto the
    target CV's unique :class:`KineticGasLiquidLink` via the
    Pattern B unchecked setter (``link._set_kLa_unchecked``).

    Construct via :meth:`from_pct_saturation` for convenience —
    pct-of-air-saturation setpoint is more familiar than mol/L.

    Anti-windup: when the agitation output is clamped at
    ``rpm_max`` (and error keeps pushing higher) or ``rpm_min``
    (and error keeps pushing lower), the integral term is
    rolled back to prevent wind-up.

    Attributes mirror the legacy :class:`PyOMES.control.loops.DOAgitationController`
    1:1; the controller-internal state (``_I_err``, ``_rpm``,
    ``_kLa_O2``, ``_kLa_CO2``) lives on ``self`` per Pattern 1.
    """

    setpoint_mol_L: float
    Kp: float = 2e5
    Ki: float = 5e4
    rpm_min: float = 100.0
    rpm_max: float = 1200.0
    rpm_initial: float = 400.0

    # Agitation→kLa mapping.
    agitation_to_kLa: Optional[Any] = None
    kLa_ref: float = 150.0
    rpm_ref: float = 400.0
    kLa_exponent: float = 2.0

    # CO2 kLa coupling.
    kLa_CO2_ratio: float = 0.9

    # Tuning.
    deadband_mol_L: float = 0.0

    # Reference 100% saturation DO for diagnostics display.
    DO_sat_ref_mol_L: Optional[float] = None

    # Optional sampling & multi-CV target.
    sample_period_s: Optional[float] = None
    sample_period_h: Optional[float] = None
    target_cv_key: Optional[str] = None
    label: str = "DO_agitation"

    # Tags and flags.
    control_tags: tuple = field(
        default_factory=lambda: ("DO", "agitation"),
    )
    requires_segmented_ode: bool = True

    # Diagnostics.
    enable_diagnostics: bool = False

    # Internal state (carried on self per Pattern 1).
    _I_err: float = field(default=0.0, init=False, repr=False)
    _rpm: float = field(default=0.0, init=False, repr=False)
    _kLa_O2: float = field(default=0.0, init=False, repr=False)
    _kLa_CO2: float = field(default=0.0, init=False, repr=False)
    diag: Dict[str, List[float]] = field(
        default_factory=dict, init=False, repr=False,
    )

    def __post_init__(self):
        self._rpm = float(self.rpm_initial)
        self._kLa_O2 = self._rpm_to_kLa(self._rpm)
        self._kLa_CO2 = self._kLa_O2 * float(self.kLa_CO2_ratio)
        if self.DO_sat_ref_mol_L is None:
            # Standard air-equilibrated water as reference: kH(O2) × p(O2)_air.
            self.DO_sat_ref_mol_L = 1.3e-3 * 0.2095

    @classmethod
    def from_pct_saturation(
        cls,
        pct: float,
        kH_O2_mol_L_atm: float = 1.3e-3,
        pO2_atm: float = 0.2095,
        **kwargs,
    ) -> "DOAgitationController":
        """Construct with the setpoint expressed as a percentage of
        air-saturation DO (e.g. 30.0 for 30%)."""
        DO_sat = float(kH_O2_mol_L_atm) * float(pO2_atm)
        setpoint = DO_sat * float(pct) / 100.0
        return cls(setpoint_mol_L=setpoint, DO_sat_ref_mol_L=DO_sat, **kwargs)

    def _rpm_to_kLa(self, rpm: float) -> float:
        rpm = float(max(0.0, rpm))
        if self.agitation_to_kLa is not None:
            return float(self.agitation_to_kLa(rpm))
        return _default_power_law_kLa(
            rpm,
            kLa_ref=self.kLa_ref,
            rpm_ref=self.rpm_ref,
            exponent=self.kLa_exponent,
        )

    def reset(self) -> None:
        self._I_err = 0.0
        self._rpm = float(self.rpm_initial)
        self._kLa_O2 = self._rpm_to_kLa(self._rpm)
        self._kLa_CO2 = self._kLa_O2 * float(self.kLa_CO2_ratio)
        if self.enable_diagnostics:
            self._reset_diagnostics()

    def _reset_diagnostics(self) -> None:
        self.diag = {
            "t_h": [],
            "DO_mol_L": [],
            "setpoint_mol_L": [],
            "err": [],
            "err_eff": [],
            "rpm": [],
            "kLa_O2": [],
            "kLa_CO2": [],
            "I_err": [],
        }

    def _diag_append(
        self, state: CVSnapshot, err: float, err_eff: float,
    ) -> None:
        if not self.enable_diagnostics:
            return
        if not self.diag:
            self._reset_diagnostics()
        DO = state.sensors.get("DO_mol_L")
        self.diag["t_h"].append(float(state.t_h))
        self.diag["DO_mol_L"].append(
            float(DO) if DO is not None else float("nan")
        )
        self.diag["setpoint_mol_L"].append(float(self.setpoint_mol_L))
        self.diag["err"].append(float(err))
        self.diag["err_eff"].append(float(err_eff))
        self.diag["rpm"].append(float(self._rpm))
        self.diag["kLa_O2"].append(float(self._kLa_O2))
        self.diag["kLa_CO2"].append(float(self._kLa_CO2))
        self.diag["I_err"].append(float(self._I_err))

    def compute(
        self,
        state: Union[CVSnapshot, SimulationSnapshot],
        dt_h: float,
    ) -> ControlAction:
        cv_snap = self._resolve_cv(state)
        empty = ControlAction(
            controller_label=self.label,
            target_cv_key=cv_snap.cv_key,
            t_h=cv_snap.t_h,
            dt_h=float(dt_h),
        )

        DO = cv_snap.sensors.get("DO_mol_L")
        if DO is None:
            return empty

        try:
            DO = float(DO)
        except (TypeError, ValueError):
            return empty

        err = float(self.setpoint_mol_L) - DO
        db = float(self.deadband_mol_L)
        err_eff = 0.0 if abs(err) <= db else err

        # PI → RPM increment, clamped to bounds.
        self._I_err += err_eff * max(float(dt_h), 0.0)
        rpm_raw = float(self.Kp) * err_eff + float(self.Ki) * self._I_err
        rpm = float(max(self.rpm_min, min(self.rpm_max, self._rpm + rpm_raw)))

        # Anti-windup: if the output saturated against the error
        # direction, undo the integral accumulation this step.
        if rpm == self.rpm_max and err_eff > 0.0:
            self._I_err -= err_eff * max(float(dt_h), 0.0)
        elif rpm == self.rpm_min and err_eff < 0.0:
            self._I_err -= err_eff * max(float(dt_h), 0.0)

        self._rpm = rpm
        self._kLa_O2 = self._rpm_to_kLa(rpm)
        self._kLa_CO2 = self._kLa_O2 * float(self.kLa_CO2_ratio)

        self._diag_append(cv_snap, err=err, err_eff=err_eff)

        return ControlAction(
            controller_label=self.label,
            target_cv_key=cv_snap.cv_key,
            t_h=cv_snap.t_h,
            dt_h=float(dt_h),
            params_changed={
                "internal_interfaces[KineticGasLiquidLink].kLa.O2": float(self._kLa_O2),
                "internal_interfaces[KineticGasLiquidLink].kLa.CO2": float(self._kLa_CO2),
            },
        )

    def _resolve_cv(
        self, state: Union[CVSnapshot, SimulationSnapshot],
    ) -> CVSnapshot:
        if isinstance(state, CVSnapshot):
            return state
        if self.target_cv_key is not None:
            if self.target_cv_key not in state.cvs:
                raise KeyError(
                    f"DOAgitationController(target_cv_key="
                    f"{self.target_cv_key!r}) not found in "
                    f"SimulationSnapshot.cvs: {sorted(state.cvs.keys())}"
                )
            return state.cvs[self.target_cv_key]
        if len(state.cvs) != 1:
            raise RuntimeError(
                f"DOAgitationController received a SimulationSnapshot "
                f"with {len(state.cvs)} CVs but has no target_cv_key set. "
                f"Set target_cv_key to disambiguate."
            )
        return next(iter(state.cvs.values()))


# Backward-compatible alias (legacy code uses DOController as the name).
DOController = DOAgitationController


# ════════════════════════════════════════════════════════════════════════
#  Cascade DO controller — Pattern 1 / CV-native
# ════════════════════════════════════════════════════════════════════════


@dataclass
class DOCascadeController:
    """Cascade DO controller — CV-native port.

    Three-tier cascade for dissolved oxygen control. Each step
    inspects the prior tier's saturation state and escalates only
    when needed:

    * **Tier 1 — Agitation (kLa):** PI loop on DO error → RPM →
      kLa(O2) via the same correlation used by
      :class:`DOAgitationController`. kLa(CO2) coupled via
      ``kLa_CO2_ratio``.
    * **Tier 2 — Gas composition:** Activates when agitation
      saturates at ``rpm_max`` (enrich O2 toward ``yO2_max``) or
      ``rpm_min`` (strip O2 toward ``yO2_min``). When agitation
      is within range, ``_yO2`` relaxes toward ``yO2_baseline`` at
      rate ``yO2_relax_rate`` (per hour).
    * **Tier 3 — Gas flow rate (optional):** Activates when both
      Tier 1 AND Tier 2 are saturated (only when
      ``enable_vvm_control=True``).

    The currently-active tier (1, 2, or 3) is exposed on
    ``self._active_tier`` for diagnostics.

    Emits a :class:`ControlAction` with ``params_changed``:

    * ``"kLa.O2"``, ``"kLa.CO2"`` (always)
    * ``"gas_feed.y.O2"``, ``"gas_feed.y.N2"`` (always — even when
      Tier 1 is doing the work, yO2 relaxes toward baseline and the
      action communicates the current value)
    * ``"gas_feed.vvm_min"`` (only when ``enable_vvm_control=True``)

    Behaviour and parameter surface mirror the legacy
    :class:`PyOMES.control.loops.DOCascadeController` 1:1.
    """

    setpoint_mol_L: float

    # Tier 1: agitation
    Kp: float = 2e5
    Ki: float = 5e4
    rpm_min: float = 100.0
    rpm_max: float = 1200.0
    rpm_initial: float = 400.0
    agitation_to_kLa: Optional[Any] = None
    kLa_ref: float = 150.0
    rpm_ref: float = 400.0
    kLa_exponent: float = 2.0
    kLa_CO2_ratio: float = 0.9

    # Tier 2: gas composition
    yO2_baseline: float = 0.21
    yO2_min: float = 0.0
    yO2_max: float = 1.0
    Kp_gas: float = 2.0
    yO2_relax_rate: float = 1.0

    # Tier 3: gas flow rate (optional)
    enable_vvm_control: bool = False
    vvm_baseline: float = 2.0
    vvm_min: float = 0.5
    vvm_max: float = 10.0
    Kp_vvm: float = 20.0

    # Tuning
    deadband_mol_L: float = 0.0

    # Reference saturation for diagnostics
    DO_sat_ref_mol_L: Optional[float] = None

    # Optional sampling & multi-CV target
    sample_period_s: Optional[float] = None
    sample_period_h: Optional[float] = None
    target_cv_key: Optional[str] = None
    label: str = "DO_cascade"

    # Tags and flags
    control_tags: tuple = field(
        default_factory=lambda: ("DO", "agitation", "gas_composition"),
    )
    requires_segmented_ode: bool = True

    # Diagnostics
    enable_diagnostics: bool = False

    # Internal state (carried on self per Pattern 1)
    _I_err: float = field(default=0.0, init=False, repr=False)
    _rpm: float = field(default=0.0, init=False, repr=False)
    _kLa_O2: float = field(default=0.0, init=False, repr=False)
    _kLa_CO2: float = field(default=0.0, init=False, repr=False)
    _yO2: float = field(default=0.21, init=False, repr=False)
    _vvm: float = field(default=2.0, init=False, repr=False)
    _active_tier: int = field(default=1, init=False, repr=False)
    diag: Dict[str, List[float]] = field(
        default_factory=dict, init=False, repr=False,
    )

    def __post_init__(self):
        self._rpm = float(self.rpm_initial)
        self._kLa_O2 = self._rpm_to_kLa(self._rpm)
        self._kLa_CO2 = self._kLa_O2 * float(self.kLa_CO2_ratio)
        self._yO2 = float(self.yO2_baseline)
        self._vvm = float(self.vvm_baseline)
        if self.DO_sat_ref_mol_L is None:
            self.DO_sat_ref_mol_L = 1.3e-3 * 0.2095

    @classmethod
    def from_pct_saturation(
        cls,
        pct: float,
        kH_O2_mol_L_atm: float = 1.3e-3,
        pO2_atm: float = 0.2095,
        **kwargs,
    ) -> "DOCascadeController":
        DO_sat = float(kH_O2_mol_L_atm) * float(pO2_atm)
        setpoint = DO_sat * float(pct) / 100.0
        return cls(setpoint_mol_L=setpoint, DO_sat_ref_mol_L=DO_sat, **kwargs)

    def _rpm_to_kLa(self, rpm: float) -> float:
        rpm = float(max(0.0, rpm))
        if self.agitation_to_kLa is not None:
            return float(self.agitation_to_kLa(rpm))
        return _default_power_law_kLa(
            rpm,
            kLa_ref=self.kLa_ref,
            rpm_ref=self.rpm_ref,
            exponent=self.kLa_exponent,
        )

    def reset(self) -> None:
        self._I_err = 0.0
        self._rpm = float(self.rpm_initial)
        self._kLa_O2 = self._rpm_to_kLa(self._rpm)
        self._kLa_CO2 = self._kLa_O2 * float(self.kLa_CO2_ratio)
        self._yO2 = float(self.yO2_baseline)
        self._vvm = float(self.vvm_baseline)
        self._active_tier = 1
        if self.enable_diagnostics:
            self._reset_diagnostics()

    def _reset_diagnostics(self) -> None:
        self.diag = {
            "t_h": [], "DO_mol_L": [], "setpoint_mol_L": [],
            "err": [], "err_eff": [], "rpm": [],
            "kLa_O2": [], "kLa_CO2": [], "yO2": [], "vvm": [],
            "active_tier": [], "I_err": [],
        }

    def _diag_append(
        self, state: CVSnapshot, err: float, err_eff: float,
    ) -> None:
        if not self.enable_diagnostics:
            return
        if not self.diag:
            self._reset_diagnostics()
        DO = state.sensors.get("DO_mol_L")
        self.diag["t_h"].append(float(state.t_h))
        self.diag["DO_mol_L"].append(
            float(DO) if DO is not None else float("nan")
        )
        self.diag["setpoint_mol_L"].append(float(self.setpoint_mol_L))
        self.diag["err"].append(float(err))
        self.diag["err_eff"].append(float(err_eff))
        self.diag["rpm"].append(float(self._rpm))
        self.diag["kLa_O2"].append(float(self._kLa_O2))
        self.diag["kLa_CO2"].append(float(self._kLa_CO2))
        self.diag["yO2"].append(float(self._yO2))
        self.diag["vvm"].append(float(self._vvm))
        self.diag["active_tier"].append(int(self._active_tier))
        self.diag["I_err"].append(float(self._I_err))

    def compute(
        self,
        state: Union[CVSnapshot, SimulationSnapshot],
        dt_h: float,
    ) -> ControlAction:
        cv_snap = self._resolve_cv(state)
        empty = ControlAction(
            controller_label=self.label,
            target_cv_key=cv_snap.cv_key,
            t_h=cv_snap.t_h,
            dt_h=float(dt_h),
        )

        DO = cv_snap.sensors.get("DO_mol_L")
        if DO is None:
            return empty
        try:
            DO = float(DO)
        except (TypeError, ValueError):
            return empty

        dt = max(float(dt_h), 0.0)
        err = float(self.setpoint_mol_L) - DO
        db = float(self.deadband_mol_L)
        err_eff = 0.0 if abs(err) <= db else err

        # ── Tier 1: PI → RPM → kLa ──
        self._I_err += err_eff * dt
        rpm_raw = float(self.Kp) * err_eff + float(self.Ki) * self._I_err
        rpm = float(max(self.rpm_min, min(self.rpm_max, self._rpm + rpm_raw)))

        agit_saturated_high = (rpm >= self.rpm_max) and (err_eff > 0.0)
        agit_saturated_low = (rpm <= self.rpm_min) and (err_eff < 0.0)
        if agit_saturated_high or agit_saturated_low:
            self._I_err -= err_eff * dt

        self._rpm = rpm
        self._kLa_O2 = self._rpm_to_kLa(rpm)
        self._kLa_CO2 = self._kLa_O2 * float(self.kLa_CO2_ratio)

        # ── Tier 2: Gas composition ──
        DO_sat_ref = (
            float(self.DO_sat_ref_mol_L) if self.DO_sat_ref_mol_L else 1.0
        )
        err_frac = err_eff / DO_sat_ref if DO_sat_ref > 0 else 0.0

        tier = 1
        if agit_saturated_high and dt > 0:
            self._yO2 += float(self.Kp_gas) * err_frac * dt
            tier = 2
        elif agit_saturated_low and dt > 0:
            self._yO2 += float(self.Kp_gas) * err_frac * dt
            tier = 2
        elif dt > 0:
            # Relax toward baseline
            relax = float(self.yO2_relax_rate) * dt
            diff = float(self.yO2_baseline) - self._yO2
            if abs(diff) > 1e-12:
                step = min(abs(diff), abs(relax)) * (
                    1.0 if diff > 0 else -1.0
                )
                self._yO2 += step

        self._yO2 = float(
            max(self.yO2_min, min(self.yO2_max, self._yO2))
        )

        # ── Tier 3: Gas flow rate (optional) ──
        if self.enable_vvm_control and dt > 0:
            comp_saturated_high = (
                (self._yO2 >= self.yO2_max) and (err_eff > 0.0)
            )
            comp_saturated_low = (
                (self._yO2 <= self.yO2_min) and (err_eff < 0.0)
            )
            if comp_saturated_high or comp_saturated_low:
                self._vvm += float(self.Kp_vvm) * err_frac * dt
                tier = 3
            else:
                vvm_diff = float(self.vvm_baseline) - self._vvm
                if abs(vvm_diff) > 1e-6:
                    vvm_relax = float(self.yO2_relax_rate) * dt
                    self._vvm += min(abs(vvm_diff), abs(vvm_relax)) * (
                        1.0 if vvm_diff > 0 else -1.0
                    )
            self._vvm = float(
                max(self.vvm_min, min(self.vvm_max, self._vvm))
            )

        self._active_tier = tier
        self._diag_append(cv_snap, err=err, err_eff=err_eff)

        # Build the action's params_changed.
        params: Dict[str, Any] = {
            "internal_interfaces[KineticGasLiquidLink].kLa.O2": float(self._kLa_O2),
            "internal_interfaces[KineticGasLiquidLink].kLa.CO2": float(self._kLa_CO2),
            "boundaries[GasFeed].y.O2": float(self._yO2),
            "boundaries[GasFeed].y.N2": float(max(0.0, 1.0 - self._yO2)),
        }
        if self.enable_vvm_control:
            params["boundaries[GasFeed].vvm_min"] = float(self._vvm)

        return ControlAction(
            controller_label=self.label,
            target_cv_key=cv_snap.cv_key,
            t_h=cv_snap.t_h,
            dt_h=float(dt_h),
            params_changed=params,
        )

    def _resolve_cv(
        self, state: Union[CVSnapshot, SimulationSnapshot],
    ) -> CVSnapshot:
        if isinstance(state, CVSnapshot):
            return state
        if self.target_cv_key is not None:
            if self.target_cv_key not in state.cvs:
                raise KeyError(
                    f"DOCascadeController(target_cv_key="
                    f"{self.target_cv_key!r}) not found in "
                    f"SimulationSnapshot.cvs: {sorted(state.cvs.keys())}"
                )
            return state.cvs[self.target_cv_key]
        if len(state.cvs) != 1:
            raise RuntimeError(
                f"DOCascadeController received a SimulationSnapshot "
                f"with {len(state.cvs)} CVs but has no target_cv_key set. "
                f"Set target_cv_key to disambiguate."
            )
        return next(iter(state.cvs.values()))


# ════════════════════════════════════════════════════════════════════════
#  Pressure-relief controllers — Pattern 1 / CV-native
# ════════════════════════════════════════════════════════════════════════

# R in (L·atm)/(mol·K) for ideal gas.
_R_L_ATM_PER_MOL_K = 0.08205736608095958


def _smooth_vent_fraction(
    excess_atm: float,
    dt_h: float,
    k_vent_per_h: float,
    smooth_width_atm: float,
) -> float:
    """Smooth venting fraction (matches legacy
    :func:`PyOMES.control.loops._smooth_vent_fraction` 1:1)."""
    import math
    excess = float(excess_atm)
    dt = float(max(0.0, dt_h))
    k = float(max(0.0, k_vent_per_h))
    w = float(max(1e-12, smooth_width_atm))
    # Smooth ReLU: 0.5*(x + sqrt(x^2 + w^2)) ~ max(0, x)
    relu = 0.5 * (excess + math.sqrt(excess * excess + w * w))
    # First-order ramp: 1 - exp(-k * relu_normalised * dt)
    if dt <= 0.0 or k <= 0.0 or relu <= 0.0:
        return 0.0
    return float(1.0 - math.exp(-k * (relu / w) * dt))


def _resolve_cv_snapshot(
    state: Union[CVSnapshot, SimulationSnapshot],
    target_cv_key: Optional[str],
    class_name: str,
) -> CVSnapshot:
    """Shared CV resolution helper for the pressure-relief controllers
    (duplicates the controller-local _resolve_cv methods, factored out
    so the three relief variants don't each carry the same boilerplate)."""
    if isinstance(state, CVSnapshot):
        return state
    if target_cv_key is not None:
        if target_cv_key not in state.cvs:
            raise KeyError(
                f"{class_name}(target_cv_key={target_cv_key!r}) not "
                f"found in SimulationSnapshot.cvs: "
                f"{sorted(state.cvs.keys())}"
            )
        return state.cvs[target_cv_key]
    if len(state.cvs) != 1:
        raise RuntimeError(
            f"{class_name} received a SimulationSnapshot with "
            f"{len(state.cvs)} CVs but has no target_cv_key set. "
            f"Set target_cv_key to disambiguate."
        )
    return next(iter(state.cvs.values()))


@dataclass
class InstantPressureReliefController:
    """Idealised instant-relief pressure controller — CV-native port.

    When ``snapshot.P_gas_atm > P_set_atm``, vents enough gas
    (proportionally across species via mole-fractions) to bring the
    pressure down to ``P_set_atm`` over ``dt_h``. Emits
    ``flux_applied = {"gas": {ID: -mol/h}}`` per species.

    This is the "idealised" sibling — no valve model, no choked-flow
    nozzle equations. Use :class:`PressureReliefController` or
    :class:`SmoothPressureReliefController` for physical valve
    modelling.
    """

    P_set_atm: float
    sample_period_s: Optional[float] = None
    sample_period_h: Optional[float] = None
    target_cv_key: Optional[str] = None
    label: str = "pressure_relief_instant"
    control_tags: tuple = field(
        default_factory=lambda: ("pressure_relief", "gas_vent"),
    )

    def reset(self) -> None:
        # No internal state to clear.
        pass

    def compute(
        self,
        state: Union[CVSnapshot, SimulationSnapshot],
        dt_h: float,
    ) -> ControlAction:
        cv = _resolve_cv_snapshot(
            state, self.target_cv_key,
            class_name="InstantPressureReliefController",
        )
        empty = ControlAction(
            controller_label=self.label,
            target_cv_key=cv.cv_key,
            t_h=cv.t_h,
            dt_h=float(dt_h),
        )
        if float(dt_h) <= 0.0:
            return empty
        if cv.P_gas_atm <= float(self.P_set_atm):
            return empty
        if cv.V_gas_L <= 0.0 or cv.T_K <= 0.0:
            return empty

        n_tot = sum(max(0.0, v) for v in cv.n_gas_mol.values())
        if n_tot <= 0.0:
            return empty

        # Ideal-gas target: n_target = P_set * V / (R * T)
        n_target = (
            float(self.P_set_atm) * float(cv.V_gas_L)
            / (_R_L_ATM_PER_MOL_K * float(cv.T_K))
        )
        vent_mol = max(0.0, n_tot - max(0.0, n_target))
        if vent_mol <= 0.0:
            return empty

        ndot_total = vent_mol / float(dt_h)
        flux_per_species = {
            ID: -float(y) * float(ndot_total)
            for ID, y in cv.y_gas.items()
            if float(y) > 0.0
        }
        vented_per_species = {
            ID: float(y) * float(vent_mol)
            for ID, y in cv.y_gas.items()
            if float(y) > 0.0
        }
        return ControlAction(
            controller_label=self.label,
            target_cv_key=cv.cv_key,
            t_h=cv.t_h,
            dt_h=float(dt_h),
            flux_applied={"gas": flux_per_species},
            vented_mol=vented_per_species,
        )


@dataclass
class SmoothPressureReliefController:
    """Differentiable pressure-relief controller — CV-native port.

    Venting fraction ramps smoothly above ``P_set_atm`` via a
    smooth-ReLU + first-order rate response, ensuring the
    pressure-vs-vent curve is differentiable (useful for
    adjoint / sensitivity work). Mirrors the legacy
    :class:`PyOMES.control.loops.SmoothPressureReliefController` 1:1.
    """

    P_set_atm: float
    k_vent_per_h: float = 500.0
    smooth_width_atm: float = 0.01
    sample_period_s: Optional[float] = None
    sample_period_h: Optional[float] = None
    target_cv_key: Optional[str] = None
    label: str = "pressure_relief_smooth"
    control_tags: tuple = field(
        default_factory=lambda: ("pressure_relief", "gas_vent"),
    )

    def reset(self) -> None:
        pass

    def compute(
        self,
        state: Union[CVSnapshot, SimulationSnapshot],
        dt_h: float,
    ) -> ControlAction:
        cv = _resolve_cv_snapshot(
            state, self.target_cv_key,
            class_name="SmoothPressureReliefController",
        )
        empty = ControlAction(
            controller_label=self.label,
            target_cv_key=cv.cv_key,
            t_h=cv.t_h,
            dt_h=float(dt_h),
        )
        if float(dt_h) <= 0.0:
            return empty

        n_tot = sum(max(0.0, v) for v in cv.n_gas_mol.values())
        if n_tot <= 0.0:
            return empty

        frac = _smooth_vent_fraction(
            excess_atm=cv.P_gas_atm - float(self.P_set_atm),
            dt_h=float(dt_h),
            k_vent_per_h=float(self.k_vent_per_h),
            smooth_width_atm=float(self.smooth_width_atm),
        )
        if frac <= 0.0:
            return empty

        vent_mol = frac * n_tot
        ndot_total = vent_mol / float(dt_h)
        flux_per_species = {
            ID: -float(y) * float(ndot_total)
            for ID, y in cv.y_gas.items()
            if float(y) > 0.0
        }
        vented_per_species = {
            ID: float(y) * float(vent_mol)
            for ID, y in cv.y_gas.items()
            if float(y) > 0.0
        }
        return ControlAction(
            controller_label=self.label,
            target_cv_key=cv.cv_key,
            t_h=cv.t_h,
            dt_h=float(dt_h),
            flux_applied={"gas": flux_per_species},
            vented_mol=vented_per_species,
        )


# ════════════════════════════════════════════════════════════════════════
#  Nozzle physics helpers (isentropic orifice flow)
# ════════════════════════════════════════════════════════════════════════

_R_UNIV = 8.31446261815324  # J/mol/K
_DEFAULT_GAMMA: Dict[str, float] = {
    "O2": 1.40, "N2": 1.40, "CO2": 1.30, "H2O": 1.33, "Air": 1.40,
}
_DEFAULT_MW_KG: Dict[str, float] = {
    "O2": 31.998e-3, "N2": 28.014e-3, "CO2": 44.009e-3, "H2O": 18.015e-3,
}
_MW_AIR_KG = 28.97e-3


def _clamp(v: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, v))


def _nozzle_mw(chemicals: Any, ID: str) -> float:
    if ID in _DEFAULT_MW_KG:
        return _DEFAULT_MW_KG[ID]
    try:
        if chemicals is not None and ID in chemicals:
            return float(chemicals[ID].MW) / 1000.0
    except (TypeError, ValueError, AttributeError):
        pass
    return _MW_AIR_KG


def _nozzle_gamma(ID: str) -> float:
    return float(_DEFAULT_GAMMA.get(ID, 1.35))


def _choked_pressure_ratio(gamma: float) -> float:
    g = float(gamma)
    return float((2.0 / (g + 1.0)) ** (g / (g - 1.0)))


def mixture_gamma_and_mw(chemicals: Any, y: Dict[str, float]) -> tuple:
    """Compute mixture heat-capacity ratio and mean MW from mole fractions."""
    cp_mix = 0.0
    mw_mix = 0.0
    for ID, yi in y.items():
        yi = float(yi)
        if yi <= 0.0:
            continue
        g = _nozzle_gamma(ID)
        cp_i = g / (g - 1.0) * _R_UNIV
        cp_mix += yi * cp_i
        mw_mix += yi * _nozzle_mw(chemicals, ID)
    mw_mix = max(mw_mix, 1e-9)
    cp_mix = max(cp_mix, _R_UNIV * 1.0001)
    gamma_mix = cp_mix / (cp_mix - _R_UNIV)
    gamma_mix = _clamp(gamma_mix, 1.05, 1.67)
    return float(gamma_mix), float(mw_mix)


def mass_flow_kg_s(
    P0_Pa: float,
    Pb_Pa: float,
    T_K: float,
    gamma: float,
    mw_kg_per_mol: float,
    Cd: float,
    A_m2: float,
) -> float:
    """Isentropic nozzle mass-flow rate (kg/s); supports choked and subsonic regimes."""
    import math
    P0 = max(float(P0_Pa), 1e-3)
    Pb = max(float(Pb_Pa), 1e-3)
    T = max(float(T_K), 1e-9)
    g = _clamp(float(gamma), 1.05, 1.67)
    Cd = max(float(Cd), 0.0)
    A = max(float(A_m2), 0.0)
    if Cd == 0.0 or A == 0.0:
        return 0.0
    R_spec = _R_UNIV / max(float(mw_kg_per_mol), 1e-12)
    pr = _clamp(Pb / P0, 0.0, 1.0)
    pr_crit = _choked_pressure_ratio(g)
    if pr <= pr_crit:
        term = (
            math.sqrt(g / (R_spec * T))
            * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
        )
        return float(Cd * A * P0 * term)
    inside = max(pr ** (2.0 / g) - pr ** ((g + 1.0) / g), 0.0)
    term = math.sqrt((2.0 * g / (R_spec * T * (g - 1.0))) * inside)
    return float(Cd * A * P0 * term)


@dataclass
class PressureReliefController:
    """Physical valve-model pressure-relief controller — CV-native port.

    Isentropic-nozzle pressure-relief controller. Supports smooth-sigmoid
    and on/off valve-opening laws (``controller_kind="smooth"`` or
    ``"onoff"``). Nozzle physics computed via :func:`mixture_gamma_and_mw`
    and :func:`mass_flow_kg_s` (module-level helpers in this file).

    The valve opens proportionally to pressure overshoot;
    species-specific molar venting follows the gas-mixture's
    average gamma and molecular weight.
    """

    P_set_atm: float
    diameter_m: float
    controller_kind: str = "smooth"
    width_atm: float = 0.01
    hysteresis_atm: float = 0.0
    Cd: float = 0.75
    P_back_atm: float = 1.0
    max_open: float = 1.0
    sample_period_s: Optional[float] = None
    sample_period_h: Optional[float] = None
    target_cv_key: Optional[str] = None
    label: str = "pressure_relief_valve"
    control_tags: tuple = field(
        default_factory=lambda: ("pressure_relief", "gas_vent"),
    )

    # Internal state — on/off latch (matches legacy
    # OnOffPressureController._is_open).
    _is_open: bool = field(default=False, init=False, repr=False)
    _last_valve_open: float = field(default=0.0, init=False, repr=False)

    def reset(self) -> None:
        self._is_open = False
        self._last_valve_open = 0.0

    def _valve_open_smooth(self, P: float) -> float:
        """Sigmoid (logistic) valve-opening law."""
        import math
        w = max(float(self.width_atm), 1e-12)
        x = max(-60.0, min(60.0, (P - float(self.P_set_atm)) / w))
        u = 1.0 / (1.0 + math.exp(-x))
        return max(0.0, min(float(self.max_open), float(u) * float(self.max_open)))

    def _valve_open_onoff(self, P: float) -> float:
        """Bang-bang with hysteresis (latched)."""
        hi = float(self.P_set_atm)
        lo = hi - float(self.hysteresis_atm)
        if self._is_open:
            if P <= lo:
                self._is_open = False
        else:
            if P >= hi:
                self._is_open = True
        return float(self.max_open) if self._is_open else 0.0

    def compute(
        self,
        state: Union[CVSnapshot, SimulationSnapshot],
        dt_h: float,
    ) -> ControlAction:
        import math

        cv = _resolve_cv_snapshot(
            state, self.target_cv_key,
            class_name="PressureReliefController",
        )
        empty = ControlAction(
            controller_label=self.label,
            target_cv_key=cv.cv_key,
            t_h=cv.t_h,
            dt_h=float(dt_h),
        )

        # Update valve_open state.
        kind = str(self.controller_kind or "smooth").lower().strip()
        if kind == "onoff":
            valve_open = self._valve_open_onoff(cv.P_gas_atm)
        else:
            valve_open = self._valve_open_smooth(cv.P_gas_atm)
        self._last_valve_open = valve_open

        if valve_open <= 0.0:
            return empty
        if cv.P_gas_atm <= float(self.P_back_atm):
            return empty

        n_tot = sum(max(0.0, v) for v in cv.n_gas_mol.values())
        if n_tot <= 0.0:
            return empty

        gamma_mix, mw_mix = mixture_gamma_and_mw(None, cv.y_gas)
        diameter = max(0.0, float(self.diameter_m))
        area_m2 = math.pi * diameter * diameter / 4.0
        P0_Pa = cv.P_gas_atm * 101325.0
        Pb_Pa = float(self.P_back_atm) * 101325.0
        mdot_kg_s = mass_flow_kg_s(
            P0_Pa=P0_Pa,
            Pb_Pa=Pb_Pa,
            T_K=cv.T_K,
            gamma=gamma_mix,
            mw_kg_per_mol=mw_mix,
            Cd=float(self.Cd),
            A_m2=area_m2 * float(valve_open),
        )
        if mdot_kg_s <= 0.0:
            return empty
        # Convert kg/s → mol/h.
        ndot_total = (mdot_kg_s / max(mw_mix, 1e-12)) * 3600.0

        flux_per_species = {
            ID: -float(y) * float(ndot_total)
            for ID, y in cv.y_gas.items()
            if float(y) > 0.0
        }
        vented_per_species = {
            ID: float(y) * float(ndot_total) * float(dt_h)
            for ID, y in cv.y_gas.items()
            if float(y) > 0.0
        }
        return ControlAction(
            controller_label=self.label,
            target_cv_key=cv.cv_key,
            t_h=cv.t_h,
            dt_h=float(dt_h),
            flux_applied={"gas": flux_per_species},
            vented_mol=vented_per_species,
        )
