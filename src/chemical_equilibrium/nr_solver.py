# -*- coding: utf-8 -*-
"""Newton-Raphson speciation solver.

Solves the equilibrium speciation system defined by an :class:`NRTableau`
for a given set of component totals and strong-ion concentrations.  The
primary entry point is :func:`solve_nr`.

Mathematical summary
--------------------
Unknowns: ``x = [x_H+, x_m1, x_m2, ...]`` where ``x_k = log10(a_k)``
(log-activity of master species ``k``).

Every species concentration is expressed analytically:
  ``c_j = 10^(log_K'_j + Σ_k  ν_jk · x_k) / γ_j``

Residuals (one per master):
  * Mass balance  (for each non-H⁺ component i):
    ``R_i = Σ_j  ν_{j,i} · c_j  −  C_{i,total}``
  * Charge balance (closes H⁺):
    ``R_charge = Σ_j  z_j · c_j  +  Σ_s  z_s · C_{s,strong}``

Analytic Jacobian (``∂c_j/∂x_k = ln(10) · ν_{jk} · c_j``):
  ``J[i,k] = ln(10) · Σ_j  ν_{j,i} · ν_{j,k} · c_j``
  ``J[charge,k] = ln(10) · Σ_j  z_j · ν_{j,k} · c_j``

Activity outer loop
-------------------
Davies (or other) activity coefficients are updated in an outer fixed-point
iteration on ionic strength, identical to the existing engine.  The NR loop
runs to convergence at each ionic-strength iterate.

Warmstart
---------
Consecutive speciation calls are nearly identical (consecutive timesteps).
Pass a :class:`NRSolverCache` object shared across calls; the cache stores
the converged ``x`` vector and ionic strength from the previous call so the
NR loop typically needs only 2–4 iterations.
"""
from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .nr_tableau import NRTableau, SecondaryEntry
from ..core.phases import R_L_ATM_MOL_K as _R_L_ATM_MOL_K

logger = logging.getLogger(__name__)

_LN10 = float(np.log(10.0))


# ─────────────────────────────────────────────────────────────────────────────
#  Warmstart cache
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class NRSolverCache:
    """Persistent state for warmstarting consecutive solve calls."""
    log_x: Optional[np.ndarray] = None   # last converged log10-activities
    I_last: Optional[float] = None        # last converged ionic strength
    jacobian: Optional[np.ndarray] = None  # last converged NR Jacobian (Phase 4)


# ─────────────────────────────────────────────────────────────────────────────
#  Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _gamma_safe(activity_model, z: int, I: float, T_K: float) -> float:
    """Call activity_model.gamma(z, I) with optional T_K, returning float."""
    I = max(0.0, float(I))
    try:
        return float(activity_model.gamma(z, I, T_K=float(T_K)))
    except TypeError:
        pass
    try:
        return float(activity_model.gamma(z, I, float(T_K)))
    except TypeError:
        pass
    return float(activity_model.gamma(z, I))


def _compute_concentrations(
    tableau: NRTableau,
    x: np.ndarray,
    gammas: Dict[str, float],
    max_log_activity: float = 50.0,
) -> Dict[str, float]:
    """Return concentration (mol/L) for every species in the tableau.

    Parameters
    ----------
    x : ndarray, shape (m,)
        Log10-activities of masters in ``tableau.masters`` order.
    gammas : dict
        ``{species_id: gamma}`` — activity coefficients.
        Missing species default to ``gamma = 1``.
    max_log_activity : float
        Log10-activity values are clamped to ``[-max_log_activity,
        +max_log_activity]`` before exponentiation to prevent floating-point
        overflow.  Default ``50.0`` (concentrations below 1e-50 or above 1e50
        are treated as the respective limit).

    Returns
    -------
    dict mapping species_id → concentration (mol/L) for liquid-phase
        species and masters; **partial pressure (atm)** for any gas-phase
        secondary (``sec.phase == "gas"``, folded by CP1 of
        ``LAYER1_GAP_CLOSURE``) — the log-linear formula's activity-space
        consistency (``a_liquid = γ_liquid·C_liquid``) already yields the
        correct Henry's-law partial pressure directly (``p_gas =
        a_liquid/kH``, verified against ``HenryEquilibrium``'s own
        ``p_i = γ_i·C_i/kH`` convention — see CP2 design notes), so no
        additional γ correction is applied here beyond what the existing
        per-species ``gammas`` division already does for every secondary.
        Callers accumulating mass balance (:func:`_residual_and_jacobian`)
        must convert this to an equivalent liquid-normalized concentration
        before summing it against liquid-phase concentrations — see that
        function's ``gas_scale`` parameter.
    """
    c: Dict[str, float] = {}
    lo, hi = -float(max_log_activity), float(max_log_activity)

    # Masters: log(a_k) = x[k], c_k = 10^x[k] / gamma_k
    for idx, m_id in enumerate(tableau.masters):
        log_a = float(np.clip(x[idx], lo, hi))
        g = gammas.get(m_id, 1.0)
        c[m_id] = 10.0 ** log_a / max(g, 1e-30)

    # Secondaries: log(a_j) = log_K' + Σ ν_jk x_k
    # Keyed by sec.c_key, not sec.species_id: a gas-phase secondary's
    # species_id commonly collides with its liquid-phase parent (e.g. both
    # "CO2") — see SecondaryEntry.c_key's docstring for why that matters.
    for sec in tableau.secondaries:
        log_a = sec.log_K_prime + sum(
            nu_k * x[tableau.masters.index(k_id)]
            for k_id, nu_k in sec.nu.items()
        )
        log_a = float(np.clip(log_a, lo, hi))
        g = gammas.get(sec.c_key, 1.0)
        c[sec.c_key] = 10.0 ** log_a / max(g, 1e-30)

    return c


def _residual_and_jacobian(
    tableau: NRTableau,
    x: np.ndarray,
    c: Dict[str, float],
    totals: Dict[str, float],
    strong_charge: float,
    pin_specs: Optional[List[Tuple[int, int, float]]] = None,
    gas_scale: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute residual vector R and Jacobian J.

    Parameters
    ----------
    c : dict
        Concentrations from :func:`_compute_concentrations` — mol/L for
        liquid species, partial pressure (atm) for gas-phase secondaries.
    totals : dict
        ``{master_id: C_total}`` for each non-H⁺ master. For a component
        with folded gas-liquid rows (CP1/CP2 of ``LAYER1_GAP_CLOSURE``),
        this is total moles *across both phases* divided by ``V_liq_L``
        (matching the existing ``PartitionModel.equilibrium_a_moles``
        convention where ``n_total`` already spans both phases) — not
        liquid-only, as it was before this phase.
    strong_charge : float
        Net charge from strong ions (Σ z_s · C_s_strong), mol/L.
    pin_specs : list of (row_idx, x_col_idx, pin_val), optional
        Pin constraints for absent components (CT = 0).  For each entry,
        the mass-balance row is replaced by ``R[row] = x[col] − pin_val``
        with a unit Jacobian, forcing the master x-value to ``pin_val``.
    gas_scale : float
        Conversion factor from a gas-phase secondary's partial pressure
        (atm) to an equivalent liquid-normalized concentration (mol/L),
        so it can be summed into the same mass-balance row as ordinary
        liquid concentrations: ``n_gas/V_liq = p_gas·V_gas/(R·T·V_liq)``,
        i.e. ``gas_scale = V_gas_L/(R_L_ATM_MOL_K·T_K·V_liq_L)``. Applied
        only to secondaries with ``phase == "gas"`` — liquid species and
        masters are unaffected (``gas_scale`` has no effect when the
        tableau has no folded gas rows; default ``0.0`` matches that
        case). The Jacobian entry scales identically since the conversion
        is linear in ``c_j``.

    Returns
    -------
    R : ndarray, shape (m,)
    J : ndarray, shape (m, m)
    """
    m = len(tableau.masters)   # number of NR unknowns
    R = np.zeros(m)
    J = np.zeros((m, m))

    # Rows that are replaced by pin equations — skip normal accumulation
    pinned_rows: set = set()
    if pin_specs:
        for row_idx, _, _ in pin_specs:
            pinned_rows.add(row_idx)

    # Indices: 0 = H+ (charge balance), 1..m-1 = components (mass balances)
    # BUT we arrange: R[0..m-2] = mass balances, R[m-1] = charge balance.
    # This keeps the charge balance as the last row, closing H+.
    # Component i corresponds to masters[i+1] (masters[0] = H+).
    # x[0] = log(a_H+), x[1] = log(a_master1), ...

    # Build look-up: master_id → column index
    m_idx = {m_id: k for k, m_id in enumerate(tableau.masters)}

    # Accumulate over ALL species (masters + secondaries)
    def _accum(sp_id: str, charge: int, nu: Dict[str, float], log_K_prime: float):
        c_j = c[sp_id]
        if c_j <= 0.0:
            return

        # Mass balance contribution for each non-H+ component
        for row_idx, comp in enumerate(tableau.components):
            if row_idx in pinned_rows:
                continue   # row is a pin equation; skip normal accumulation
            nu_ji = nu.get(comp.master_id, 0.0)
            if nu_ji == 0.0:
                continue
            R[row_idx] += nu_ji * c_j
            # Jacobian: ∂R[row]/∂x[col] = ln10 * nu_ji * nu_jcol * c_j
            for col_master_id, nu_jk in nu.items():
                col = m_idx.get(col_master_id)
                if col is None:
                    continue
                J[row_idx, col] += _LN10 * nu_ji * nu_jk * c_j

        # Charge balance contribution (last row)
        charge_row = m - 1
        if charge != 0:
            R[charge_row] += charge * c_j
            for col_master_id, nu_jk in nu.items():
                col = m_idx.get(col_master_id)
                if col is None:
                    continue
                J[charge_row, col] += _LN10 * charge * nu_jk * c_j

    # Masters
    for idx, m_id in enumerate(tableau.masters):
        charge = tableau.master_charges[m_id]
        nu = {m_id: 1.0}
        _accum(m_id, charge, nu, 0.0)

    # Secondaries — mass-balance row uses element_stoichiometry so that
    # cross-component species contribute to every component whose master
    # appears in their formula.  The Jacobian still iterates the full nu
    # (including H+) so d(c_j)/d(x_H+) is correctly propagated.
    for sec in tableau.secondaries:
        c_j = c.get(sec.c_key, 0.0)
        if c_j <= 0.0:
            continue
        if sec.phase == "gas":
            # Convert partial pressure (atm) to an equivalent liquid-
            # normalized concentration (mol/L) before accumulating —
            # see the gas_scale parameter doc above.
            c_j = c_j * gas_scale
            if c_j <= 0.0:
                continue
        # Mass-balance rows
        for row_idx, comp in enumerate(tableau.components):
            if row_idx in pinned_rows:
                continue
            coeff = sec.element_stoichiometry.get(comp.master_id, 0.0)
            if coeff == 0.0:
                continue
            R[row_idx] += coeff * c_j
            for col_master_id, nu_jk in sec.nu.items():
                col = m_idx.get(col_master_id)
                if col is None:
                    continue
                J[row_idx, col] += _LN10 * coeff * nu_jk * c_j
        # Charge-balance row
        charge_row = m - 1
        if sec.charge != 0:
            R[charge_row] += sec.charge * c_j
            for col_master_id, nu_jk in sec.nu.items():
                col = m_idx.get(col_master_id)
                if col is None:
                    continue
                J[charge_row, col] += _LN10 * sec.charge * nu_jk * c_j

    # Subtract totals from mass-balance rows (skip pinned rows)
    for row_idx, comp in enumerate(tableau.components):
        if row_idx in pinned_rows:
            continue
        CT = totals.get(comp.master_id, 0.0)
        R[row_idx] -= CT

    # Add strong-ion charge to charge-balance row
    R[m - 1] += strong_charge

    # Apply pin constraints: replace absent-component mass balance rows
    if pin_specs:
        for row_idx, x_col, pin_val in pin_specs:
            R[row_idx] = float(x[x_col]) - pin_val
            J[row_idx, :] = 0.0
            J[row_idx, x_col] = 1.0

    return R, J


def _ionic_strength(
    c: Dict[str, float],
    charges: Dict[str, int],
    strong_ions: Dict[str, float],
) -> float:
    """Compute ionic strength I = 0.5 Σ z² c from speciation + strong ions."""
    I = 0.0
    for sp_id, conc in c.items():
        z = charges.get(sp_id, 0)
        if z != 0:
            I += z * z * conc
    # Strong ions: {CT_Na: mol/L, CT_Cl: mol/L, ...}  charges known from sign
    _STRONG_CHARGES = {
        "CT_K": +1, "CT_Na": +1, "CT_cation": +1,
        "CT_Cl": -1, "CT_NO3": -1, "CT_anion": -1,
        "CT_Mg": +2, "CT_Ca": +2, "CT_Zn": +2, "CT_Mn": +2, "CT_Cu": +2, "CT_Co": +2,
        "CT_Fe2": +2,
        "CT_Mo7O24": -6, "CT_MoO4": -2,
    }
    for key, val in strong_ions.items():
        z = _STRONG_CHARGES.get(key, 0)
        if z != 0:
            I += z * z * float(val)
    return 0.5 * I


# ─────────────────────────────────────────────────────────────────────────────
#  Public solver
# ─────────────────────────────────────────────────────────────────────────────

def solve_nr(
    tableau: NRTableau,
    totals: Dict[str, float],
    strong_ions: Dict[str, float],
    *,
    T_K: float = 298.15,
    activity_model,
    tol: float = 1e-10,
    max_inner: int = 50,
    max_outer: int = 20,
    max_backtrack: int = 10,
    cache: Optional[NRSolverCache] = None,
    I_init: float = 0.01,
    I_tol: float = 1e-8,
    damping: float = 0.3,
    max_log_activity: float = 50.0,
    min_component_total: float = 1e-20,
    retain_jacobian: bool = False,
    V_liq_L: Optional[float] = None,
    V_gas_L: Optional[float] = None,
) -> Dict[str, Any]:
    """Solve the NR speciation system.

    Parameters
    ----------
    tableau : NRTableau
        Pre-built log-linear tableau (temperature-corrected log_K values).
    totals : dict
        ``{master_id: C_total}`` (mol/L) for each non-H⁺ master.
        Missing keys are treated as zero.
    strong_ions : dict
        ``{CT_Na: mol/L, CT_Cl: mol/L, ...}`` — fully dissociated ions.
        These enter the charge balance but have no mass-balance equation.
    T_K : float
        Operating temperature (K), used for activity coefficients.
    activity_model : ActivityModel
        Provides ``gamma(z, I)`` or ``gamma(z, I, T_K=...)``.
    tol : float
        Convergence threshold: ``||R||_∞ < tol``.
    max_inner / max_outer : int
        Newton iteration and activity-outer loop limits.
    max_backtrack : int
        Maximum backtracking halvings per Newton step.
    cache : NRSolverCache, optional
        Mutable warmstart cache shared across calls.
    I_init : float
        Initial ionic strength if ``cache.I_last`` is unset.
    I_tol : float
        Outer-loop convergence threshold on |ΔI|.
    damping : float
        Damping factor for ionic-strength update (0 < damping < 1).
    max_log_activity : float
        Log10-activity values are clamped to ``[-max_log_activity,
        +max_log_activity]`` before exponentiation to prevent floating-point
        overflow.  Default ``50.0``.
    min_component_total : float
        Components whose total concentration is below this threshold (mol/L)
        are treated as absent.  Their master log-activity is pinned to
        ``-max_log_activity`` throughout the NR solve, avoiding divergence
        toward ``−∞``.  Default ``1e-20``.
    V_liq_L, V_gas_L : float, optional
        Liquid- and gas-phase volumes (litres). Required (both, > 0) when
        *tableau* carries any gas-phase secondary (folded by CP1 of
        ``LAYER1_GAP_CLOSURE``) — used to convert that secondary's partial
        pressure into an equivalent liquid-normalized concentration for
        the mass-balance row it's attached to (CP2). Ignored when the
        tableau has no folded gas rows.

    Returns
    -------
    dict
        ``{species_id: concentration, "pH": ..., "logH": ...,
           "IonicStrength": ..., "gamma_H": ..., "gamma_OH": ...,
           "pH_conc": ..., "aH": ..., "charge_residual": ...}``.
        For a folded gas-phase secondary, its entry is a partial pressure
        (atm), not a concentration — see :func:`_compute_concentrations`.

    Raises
    ------
    ValueError
        If *tableau* carries any gas-phase secondary but ``V_liq_L``/
        ``V_gas_L`` are not both supplied (> 0) — folding a gas-liquid row
        without volumes would silently drop its mass-balance contribution
        (``gas_scale`` would be ambiguous), which this refuses rather than
        doing silently.
    """
    m = len(tableau.masters)   # number of NR unknowns
    pin_val = -float(max_log_activity)

    gas_species = sorted({sec.species_id for sec in tableau.secondaries if sec.phase == "gas"})
    if gas_species:
        if not (V_liq_L is not None and V_gas_L is not None
                and float(V_liq_L) > 0.0 and float(V_gas_L) > 0.0):
            raise ValueError(
                f"solve_nr: tableau folds gas-liquid secondaries {gas_species} "
                "(LAYER1_GAP_CLOSURE CP1/CP2), which requires V_liq_L and "
                "V_gas_L (both > 0) to convert gas partial pressure into an "
                "equivalent mass-balance contribution. Pass them explicitly, "
                "or via phases={'liquid':..., 'gas':...} on "
                "NRChemicalEquilibriumEngine.solve()."
            )
        gas_scale = float(V_gas_L) / (_R_L_ATM_MOL_K * float(T_K) * float(V_liq_L))
    else:
        gas_scale = 0.0

    # ── Charge dictionary for all species ────────────────────────────
    # Keyed by sec.c_key (not sec.species_id) for the same reason as
    # _compute_concentrations's c dict — see SecondaryEntry.c_key.
    all_charges: Dict[str, int] = dict(tableau.master_charges)
    for sec in tableau.secondaries:
        all_charges[sec.c_key] = sec.charge

    # ── Strong-ion net charge (constant across NR iterations) ─────────
    _STRONG_CHARGES = {
        "CT_K": +1, "CT_Na": +1, "CT_cation": +1,
        "CT_Cl": -1, "CT_NO3": -1, "CT_anion": -1,
        "CT_Mg": +2, "CT_Ca": +2, "CT_Zn": +2, "CT_Mn": +2, "CT_Cu": +2, "CT_Co": +2,
        "CT_Fe2": +2,
        "CT_Mo7O24": -6, "CT_MoO4": -2,
    }
    strong_charge = sum(
        _STRONG_CHARGES.get(key, 0) * float(val)
        for key, val in strong_ions.items()
    )

    # ── Initial x vector ──────────────────────────────────────────────
    if cache is not None and cache.log_x is not None and len(cache.log_x) == m:
        x = cache.log_x.copy()
    else:
        x = _initial_x(tableau, totals)

    # ── Pin absent components ─────────────────────────────────────────
    # Components with CT < min_component_total cannot be solved by the normal
    # mass-balance equation (it would require x → -∞).  Replace the
    # mass-balance row with a pin equation x[col] = pin_val, keeping the
    # NR system square and well-conditioned.
    #
    # Also handle the reverse transition: if a warmstart cache carried x[col]
    # at (or near) pin_val for a component that is now active, the Jacobian
    # contributions from that component are ≈ 0 and convergence would stall.
    # Reset those components to the heuristic initial guess.
    pin_specs: List[Tuple[int, int, float]] = []
    for row_idx, comp in enumerate(tableau.components):
        CT = float(totals.get(comp.master_id, 0.0))
        x_col = tableau.masters.index(comp.master_id)
        if CT < float(min_component_total):
            pin_specs.append((row_idx, x_col, pin_val))
            x[x_col] = pin_val
        elif x[x_col] <= pin_val + 1.0:
            # Was pinned in previous call; reset to sensible starting point
            x[x_col] = np.log10(max(CT, 1e-10))

    # ── Initial ionic strength ────────────────────────────────────────
    I = (cache.I_last if cache is not None and cache.I_last is not None
         else float(I_init))

    is_ideal = getattr(activity_model, "name", "ideal") == "ideal"

    # ── Outer activity loop ───────────────────────────────────────────
    outer_converged = False
    for _outer in range(int(max_outer)):
        # Build gamma dict for all species at current I
        gammas = _build_gammas(tableau, all_charges, I, T_K, activity_model)

        # ── Inner NR loop ─────────────────────────────────────────────
        inner_converged = False
        for _inner in range(int(max_inner)):
            c = _compute_concentrations(tableau, x, gammas,
                                        max_log_activity=max_log_activity)
            R, J = _residual_and_jacobian(tableau, x, c, totals, strong_charge,
                                          pin_specs=pin_specs, gas_scale=gas_scale)

            norm_R = float(np.max(np.abs(R)))
            if norm_R < tol:
                inner_converged = True
                break

            # Solve J Δx = -R
            try:
                dx = np.linalg.solve(J, -R)
            except np.linalg.LinAlgError:
                # Singular Jacobian — fall back to gradient-descent step
                dx = -R * 0.01
                logger.warning("NR Jacobian singular at iter %d/%d", _outer, _inner)

            # Backtracking line search
            alpha = 1.0
            norm_R_sq = float(np.dot(R, R))
            for _bt in range(int(max_backtrack)):
                x_trial = x + alpha * dx
                c_trial = _compute_concentrations(tableau, x_trial, gammas,
                                                  max_log_activity=max_log_activity)
                R_trial, _ = _residual_and_jacobian(
                    tableau, x_trial, c_trial, totals, strong_charge,
                    pin_specs=pin_specs, gas_scale=gas_scale,
                )
                if float(np.dot(R_trial, R_trial)) < norm_R_sq:
                    break
                alpha *= 0.5

            x = x + alpha * dx

        if not inner_converged:
            logger.warning(
                "NR inner loop did not converge (outer=%d, ||R||=%.2e)",
                _outer, float(np.max(np.abs(R))),
            )

        if is_ideal:
            outer_converged = True
            break

        # Update ionic strength
        I_new = _ionic_strength(c, all_charges, strong_ions)
        if abs(I_new - I) < float(I_tol):
            I = I_new
            outer_converged = True
            break
        I = float(damping) * I + (1.0 - float(damping)) * I_new

    if not outer_converged and not is_ideal:
        warnings.warn(
            f"NR speciation: activity outer loop did not converge "
            f"(max_outer={max_outer}, |ΔI| last={abs(I_new - I):.2e}).",
            RuntimeWarning,
        )

    # ── Final species evaluation ───────────────────────────────────────
    gammas = _build_gammas(tableau, all_charges, I, T_K, activity_model)
    c = _compute_concentrations(tableau, x, gammas, max_log_activity=max_log_activity)

    # Final residual (and Jacobian when caller requests retention)
    R_final, J_final = _residual_and_jacobian(tableau, x, c, totals, strong_charge,
                                               pin_specs=pin_specs, gas_scale=gas_scale)
    charge_residual = float(R_final[-1])   # charge balance is the last row

    # ── Build output dict ─────────────────────────────────────────────
    # Concentrations for all species
    out: Dict[str, Any] = dict(c)

    # pH and related
    c_H = float(c.get("H+", 1e-7))
    gamma_H = float(gammas.get("H+", 1.0))
    gamma_OH = float(gammas.get("OH-", 1.0))
    a_H = gamma_H * c_H
    pH_conc = -np.log10(max(c_H, 1e-30))
    pH = -np.log10(max(a_H, 1e-30))

    out["pH"] = float(pH)
    out["pH_conc"] = float(pH_conc)
    out["logH"] = float(np.log10(max(c_H, 1e-30)))
    out["aH"] = float(a_H)
    out["gamma_H"] = float(gamma_H)
    out["gamma_OH"] = float(gamma_OH)
    out["IonicStrength"] = float(I)
    out["charge_residual"] = float(charge_residual)

    # Update warmstart cache (and optional Jacobian retention)
    if cache is not None:
        cache.log_x = x.copy()
        cache.I_last = float(I)
        cache.jacobian = J_final.copy()

    # Pack white-box state into output for retain_jacobian=True callers
    if retain_jacobian:
        out["_jacobian_matrix"] = J_final.copy()
        out["_log_activities"] = x.copy()
        out["_ionic_strength_final"] = float(I)

    return out


# ─────────────────────────────────────────────────────────────────────────────
#  Private helpers
# ─────────────────────────────────────────────────────────────────────────────

def _initial_x(tableau: NRTableau, totals: Dict[str, float]) -> np.ndarray:
    """Make a reasonable initial guess for the NR unknown vector."""
    m = len(tableau.masters)
    x = np.zeros(m)
    # H+: start at pH 7 (log10([H+]) = -7)
    x[0] = -7.0
    # Each non-H+ master: log10 of its total concentration (or a small default)
    for idx, m_id in enumerate(tableau.masters[1:], start=1):
        CT = totals.get(m_id, 0.0)
        x[idx] = np.log10(max(CT, 1e-10))
    return x


def _build_gammas(
    tableau: NRTableau,
    all_charges: Dict[str, int],
    I: float,
    T_K: float,
    activity_model,
) -> Dict[str, float]:
    """Compute activity coefficients for all species in the tableau."""
    gammas: Dict[str, float] = {}
    for sp_id, z in all_charges.items():
        gammas[sp_id] = _gamma_safe(activity_model, z, I, T_K)
    return gammas
