# -*- coding: utf-8 -*-
"""Visualization utilities for reaction objects.

Matplotlib is imported lazily inside each function so that importing
:mod:`PyOMES.reactions` never requires a display environment or a matplotlib
installation.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, List, Optional, Sequence, Tuple

if TYPE_CHECKING:
    from .equilibrium import EquilibriumReaction
    from .reaction_system import ReactionSystem

_R_GAS = 8.314  # J / (mol · K)

# Species excluded from ladder-member detection — they are the solvent/proton
# channel and never form the "pool" being speciated.
_SOLVENT_IDS = frozenset({"H+", "OH-", "H2O"})


def plot_vant_hoff(
    rxn: "EquilibriumReaction",
    *,
    T_range_K: Tuple[float, float] = (273.15, 373.15),
    n_points: int = 200,
    ax=None,
    show_pka: bool = True,
    title: Optional[str] = None,
) -> "tuple":
    """Plot the Van 't Hoff temperature dependence of an EquilibriumReaction.

    Uses the Van 't Hoff equation:

    .. math::

        \\log_{10} K(T) = \\log_{10} K(T_{\\mathrm{ref}})
            + \\frac{\\Delta H^\\circ}{R \\ln 10}
              \\left( \\frac{1}{T_{\\mathrm{ref}}} - \\frac{1}{T} \\right)

    When ``dH_J_per_mol`` is ``None`` on the reaction, the equilibrium
    constant is temperature-independent and a flat line is drawn with an
    annotation noting that no Van 't Hoff data is available.

    Parameters
    ----------
    rxn : EquilibriumReaction
        The reaction to plot.  Must have a ``log_K`` value.
    T_range_K : (float, float)
        Temperature range to plot in Kelvin.  Default ``(273.15, 373.15)``
        (0 – 100 °C).
    n_points : int
        Number of temperature points in the curve.  Default 200.
    ax : matplotlib Axes, optional
        Axes to draw on.  A new figure and axes are created when
        ``None``.
    show_pka : bool
        When ``True`` (default), add a secondary y-axis showing
        pKa = –log₁₀ K.
    title : str, optional
        Axes title.  Defaults to the reaction label.

    Returns
    -------
    fig, ax
        The matplotlib Figure and primary Axes objects.

    Raises
    ------
    ValueError
        If ``rxn.log_K`` is ``None`` (cross-phase partition declarations
        have no equilibrium constant).
    ImportError
        If matplotlib is not installed.
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        raise ImportError(
            "plot_vant_hoff requires matplotlib and numpy. "
            "Install them with:  pip install matplotlib numpy"
        ) from exc

    if rxn.log_K is None:
        raise ValueError(
            f"EquilibriumReaction '{rxn.label}' has no log_K "
            "(cross-phase partition declarations omit the equilibrium constant). "
            "plot_vant_hoff requires a numeric log_K."
        )

    T_arr = np.linspace(T_range_K[0], T_range_K[1], n_points)
    T_C   = T_arr - 273.15

    if rxn.dH_J_per_mol is not None:
        log_K_arr = rxn.log_K + (rxn.dH_J_per_mol / (_R_GAS * math.log(10))) * (
            1.0 / rxn.T_ref_K - 1.0 / T_arr
        )
        temperature_dependent = True
    else:
        log_K_arr = np.full_like(T_arr, rxn.log_K)
        temperature_dependent = False

    # ── figure setup ─────────────────────────────────────────────────────────
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4))
    else:
        fig = ax.get_figure()

    ax.plot(T_C, log_K_arr, color="steelblue", linewidth=2)

    # Mark the reference point
    T_ref_C = rxn.T_ref_K - 273.15
    ax.scatter(
        [T_ref_C], [rxn.log_K],
        color="steelblue", s=60, zorder=5,
        label=f"$T_{{\\mathrm{{ref}}}}$ = {T_ref_C:.1f} °C, "
              f"$\\log_{{10}} K$ = {rxn.log_K:.3f}",
    )

    if not temperature_dependent:
        ax.annotate(
            "No Van 't Hoff data (dH not specified)\n— constant log K assumed",
            xy=(0.5, 0.92), xycoords="axes fraction",
            ha="center", fontsize=9, color="gray",
            bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="gray"),
        )

    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("$\\log_{10} K$")
    ax.set_title(title if title is not None else rxn.label)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # ── optional pKa secondary axis ───────────────────────────────────────────
    if show_pka:
        ax2 = ax.twinx()
        ax2.plot(T_C, -log_K_arr, alpha=0)  # invisible — just sets the scale
        ax2.set_ylabel("pKa  ($-\\log_{10} K$)")
        # Sync limits so both axes stay consistent
        y1_lo, y1_hi = ax.get_ylim()
        ax2.set_ylim(-y1_hi, -y1_lo)

    fig.tight_layout()
    return fig, ax


# ── Speciation ladder helpers ─────────────────────────────────────────────────

def _build_ladder_from_system(
    rxn_system: "ReactionSystem",
    anchor_id: str,
) -> "tuple[List[str], List[float], list]":
    """Return (chain_ids, log_Ks, chain_rxns) for the ladder anchored at anchor_id.

    Scans ``rxn_system.single_phase_equilibria`` for reactions of the form::

        acid_species [+ H2O] ⇌ base_species + H⁺

    (i.e. exactly one non-solvent reactant, one non-solvent product, H⁺
    produced).  Builds a directed graph acid→base and walks it to return an
    ordered chain starting from the most-protonated form.

    Raises
    ------
    ValueError
        If ``anchor_id`` is not found in any ladder reaction, if the
        ladder is not a simple linear chain, or if fewer than two species
        are connected.
    """
    edges: list = []  # (acid_id, base_id, log_K, rxn)
    for rxn in rxn_system.single_phase_equilibria:
        if rxn.log_K is None:
            continue
        non_solvent_reactants = [
            e for e in rxn.stoichiometry
            if e.coefficient < 0 and e.species.id not in _SOLVENT_IDS
        ]
        non_solvent_products = [
            e for e in rxn.stoichiometry
            if e.coefficient > 0 and e.species.id not in _SOLVENT_IDS
        ]
        h_produced = any(
            e.species.id == "H+" and e.coefficient > 0
            for e in rxn.stoichiometry
        )
        if (
            len(non_solvent_reactants) == 1
            and len(non_solvent_products) == 1
            and h_produced
        ):
            edges.append((
                non_solvent_reactants[0].species.id,
                non_solvent_products[0].species.id,
                rxn.log_K,
                rxn,
            ))

    # Build adjacency maps
    forward: dict = {}   # acid_id -> (base_id, log_K, rxn)
    backward: dict = {}  # base_id -> (acid_id, log_K, rxn)
    all_ids: set = set()
    for acid_id, base_id, log_K, rxn in edges:
        forward[acid_id] = (base_id, log_K, rxn)
        backward[base_id] = (acid_id, log_K, rxn)
        all_ids.add(acid_id)
        all_ids.add(base_id)

    if anchor_id not in all_ids:
        raise ValueError(
            f"Anchor species '{anchor_id}' not found in any acid-base ladder "
            f"reaction within the ReactionSystem. "
            f"Available species: {sorted(all_ids)}"
        )

    # Walk backward to find the root (most protonated species)
    root = anchor_id
    seen: set = {root}
    while root in backward:
        prev = backward[root][0]
        if prev in seen:
            raise ValueError(
                f"Cycle detected in speciation ladder near '{root}'."
            )
        seen.add(prev)
        root = prev

    # Walk forward from root to build the ordered chain
    chain_ids: List[str] = [root]
    chain_rxns: list = []
    log_Ks: List[float] = []
    current = root
    while current in forward:
        base_id, log_K, rxn = forward[current]
        chain_ids.append(base_id)
        log_Ks.append(log_K)
        chain_rxns.append(rxn)
        current = base_id

    if len(chain_ids) < 2:
        raise ValueError(
            f"Speciation ladder anchored at '{anchor_id}' has only one "
            "species — at least two are needed to plot a distribution."
        )

    return chain_ids, log_Ks, chain_rxns


def plot_speciation(
    rxn_system: "ReactionSystem",
    anchor_id: str,
    *,
    pH_range: Tuple[float, float] = (0.0, 14.0),
    T_K: Optional[float] = 298.15,
    n_points: int = 500,
    ax=None,
    title: Optional[str] = None,
) -> "tuple":
    """Plot molar-fraction vs pH for the speciation ladder anchored at *anchor_id*.

    Builds the acid-base chain automatically by scanning
    ``rxn_system.single_phase_equilibria`` for reactions of the form
    ``acid ⇌ base + H⁺``.  The chain is ordered from the most-protonated
    form (highest H⁺ count) to the least.

    Alpha fractions are computed analytically using the standard polyprotic
    formula.  When ``T_K`` is supplied and a reaction carries ``dH_J_per_mol``,
    the Van 't Hoff correction is applied before computing the fractions.

    Parameters
    ----------
    rxn_system : ReactionSystem
        The system whose ``single_phase_equilibria`` are searched for
        the ladder.
    anchor_id : str
        Species ID of any member of the ladder (e.g. ``"CO2"``,
        ``"HCO3-"``, ``"AceticAcid"``).  The function walks to the
        most-protonated end automatically.
    pH_range : (float, float)
        pH axis limits.  Default ``(0, 14)``.
    T_K : float or None
        Temperature for Van 't Hoff correction.  ``None`` uses stored
        ``log_K`` values without correction.  Default ``298.15`` K.
    n_points : int
        Number of pH points.  Default 500.
    ax : matplotlib Axes, optional
        Axes to draw on.  A new figure is created when ``None``.
    title : str, optional
        Axes title.  Defaults to ``"Speciation: <root> ladder"``.

    Returns
    -------
    fig, ax
        The matplotlib Figure and primary Axes objects.

    Raises
    ------
    ValueError
        If ``anchor_id`` is not found in the system, the ladder is
        non-linear, or fewer than two species are connected.
    ImportError
        If matplotlib or numpy are not installed.
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        raise ImportError(
            "plot_speciation requires matplotlib and numpy. "
            "Install them with:  pip install matplotlib numpy"
        ) from exc

    chain_ids, log_Ks, chain_rxns = _build_ladder_from_system(rxn_system, anchor_id)

    # Apply optional Van 't Hoff correction to get temperature-specific pKa values
    pKas: List[float] = []
    for log_K, rxn in zip(log_Ks, chain_rxns):
        if T_K is not None and rxn.dH_J_per_mol is not None:
            log_K_T = log_K + (rxn.dH_J_per_mol / (_R_GAS * math.log(10))) * (
                1.0 / rxn.T_ref_K - 1.0 / T_K
            )
        else:
            log_K_T = log_K
        pKas.append(-log_K_T)  # pKa = -log10(Ka)

    # ── Analytical alpha fractions ─────────────────────────────────────────────
    # For n species with Ka_1 … Ka_(n-1):
    #   beta_0 = 1,  beta_k = Ka_1 * … * Ka_k
    #   D(h)   = sum_k  beta_k * h^(n-1-k)
    #   alpha_k = beta_k * h^(n-1-k) / D(h)
    pH = np.linspace(pH_range[0], pH_range[1], n_points)
    h = 10.0 ** (-pH)

    n = len(chain_ids) - 1   # number of dissociation steps
    Kas = [10.0 ** (-pKa) for pKa in pKas]
    betas = [1.0]
    prod = 1.0
    for Ka in Kas:
        prod *= Ka
        betas.append(prod)

    D = np.zeros(n_points)
    for k, beta in enumerate(betas):
        D += beta * h ** (n - k)

    alphas = [beta * h ** (n - k) / D for k, beta in enumerate(betas)]

    # ── Plot ──────────────────────────────────────────────────────────────────
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = ax.get_figure()

    colors = plt.cm.tab10.colors
    for k, (sp_id, alpha) in enumerate(zip(chain_ids, alphas)):
        ax.plot(pH, alpha, label=sp_id, color=colors[k % len(colors)], linewidth=2)

    # Vertical dashed lines at each pKa, annotated above the plot
    for i, pKa in enumerate(pKas):
        if pH_range[0] <= pKa <= pH_range[1]:
            ax.axvline(pKa, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
            ax.text(
                pKa + 0.2, 0.97,
                f"pKa{i + 1} = {pKa:.2f}",
                ha="left", va="top", fontsize=8, color="gray",
                transform=ax.get_xaxis_transform(),
            )

    ax.set_xlabel("pH")
    ax.set_ylabel("Molar fraction α")
    ax.set_xlim(pH_range)
    ax.set_ylim(0.0, 1.05)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    _title = title if title is not None else f"Speciation: {chain_ids[0]} ladder"
    if T_K is not None:
        _title += f"  (T = {T_K - 273.15:.1f} °C)"
    ax.set_title(_title)

    fig.tight_layout()
    return fig, ax
