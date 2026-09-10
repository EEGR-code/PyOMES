# -*- coding: utf-8 -*-
"""Convenience builders for constructing validated KineticReaction objects.

:class:`ReactionBuilder` provides static methods that derive complete
stoichiometric coefficients from high-level specifications (substrate
formula, biomass formula, yield, balance mode) and produce validated
:class:`~PyOMES.reactions.kinetic.KineticReaction` objects.

High-level convenience method (Monod kinetics + stoichiometry in one call):

>>> rxn = ReactionBuilder.monod_aerobic_growth(
...     substrate    = ACETIC_ACID,   # Species object
...     biomass      = YEAST,         # Species object
...     mu_max_per_h = 0.5,
...     Ks_gL        = 5e-3,
...     yield_gX_gS  = 0.36,
... )

Lower-level method for custom rate functions:

>>> rxn = ReactionBuilder.aerobic_growth(
...     substrate_id    = "AceticAcid",
...     substrate_atoms = {"C": 2, "H": 4, "O": 2},
...     MW_substrate    = 60.052,
...     biomass_id      = "Yeast",
...     biomass_atoms   = {"C": 1, "H": 1.61, "O": 0.56},
...     MW_biomass      = 24.626,
...     yield_gX_gS     = 0.36,
...     rate_fn         = my_custom_rate,
... )
"""

from __future__ import annotations

from typing import Callable, Dict, Mapping, Optional, Sequence

from .stoichiometry import StoichiometryEntry
from .kinetic import KineticReaction
from .environment import ReactionEnvironment
from ..chemistry.species import Species
from ..chemistry.common_species import CO2 as _CO2, H2O as _H2O, NH3 as _NH3


class ReactionBuilder:
    """Static factory methods for constructing validated reactions."""

    @staticmethod
    def aerobic_growth(
        substrate_id: str,
        substrate_atoms: Mapping[str, float],
        MW_substrate: float,
        biomass_id: str,
        biomass_atoms: Mapping[str, float],
        MW_biomass: float,
        yield_gX_gS: float,
        rate_fn: Callable[[ReactionEnvironment], float],
        *,
        balance: str = "CHO",
        n_source_id: str = "NH3",
        n_source_atoms: Mapping[str, float] | None = None,
        label: str = "",
        species_overrides: Mapping[str, "Species"] | None = None,
    ) -> KineticReaction:
        """Build an aerobic growth reaction from a yield and elemental formulas.

        Derives the stoichiometric coefficients for O₂ (consumed),
        CO₂ (produced), and H₂O (produced) from elemental balance.
        In CHNO mode, also derives the nitrogen source demand.

        The reaction basis is **1 mol of substrate consumed** (coefficient
        = −1).  The rate function should therefore return the substrate
        consumption rate in mol/h (extensive).

        Parameters
        ----------
        substrate_id : str
            Species ID for the substrate (e.g. ``"AceticAcid"``).
        substrate_atoms : Mapping[str, float]
            Elemental composition of the substrate,
            e.g. ``{"C": 2, "H": 4, "O": 2}``.  Any read-only mapping
            is accepted, including ``Species.atoms`` directly.
        MW_substrate : float
            Molecular weight of the substrate (g/mol).
        biomass_id : str
            Species ID for the organism (e.g. ``"Yeast"``).
        biomass_atoms : Mapping[str, float]
            Elemental composition of the biomass *per formula unit*,
            e.g. ``{"C": 1, "H": 1.61, "O": 0.56}`` for CH₁.₆₁O₀.₅₆.
            ``Species.atoms`` can be passed directly.
        MW_biomass : float
            Molecular weight of one biomass formula unit (g/mol).
        yield_gX_gS : float
            Yield coefficient in g biomass per g substrate consumed.
        rate_fn : callable
            ``rate_fn(env) -> float`` returning the extensive substrate
            consumption rate (mol substrate / h).
        balance : str
            ``"CHO"`` (default) or ``"CHNO"``.  When ``"CHNO"``, the
            nitrogen source is included in the stoichiometry.
        n_source_id : str
            Species ID for the nitrogen source (default ``"NH3"``).
            Only used when ``balance`` includes ``"N"``.
        n_source_atoms : dict or None
            Elemental composition of the nitrogen source.
            Defaults to ``{"N": 1, "H": 3}`` for NH₃.
        label : str
            Human-readable label for the reaction.

        Returns
        -------
        KineticReaction
            A fully validated ``KineticReaction`` with stoichiometry
            closed for the requested elements.

        Raises
        ------
        StoichiometryError
            If the derived stoichiometry does not close (should not happen
            when formulas and yields are self-consistent, but serves as a
            safety net against transcription errors).
        ValueError
            If inputs are non-physical (zero MW, negative yield, etc.).
        """
        if MW_substrate <= 0.0:
            raise ValueError(f"MW_substrate must be > 0, got {MW_substrate}")
        if MW_biomass <= 0.0:
            raise ValueError(f"MW_biomass must be > 0, got {MW_biomass}")
        if yield_gX_gS < 0.0:
            raise ValueError(f"yield_gX_gS must be >= 0, got {yield_gX_gS}")

        # Normalise balance mode
        balance_upper = balance.upper().strip()
        balance_elements = list(balance_upper)  # ["C","H","O"] or ["C","H","N","O"]

        # mol biomass produced per mol substrate consumed
        Y_mol = yield_gX_gS * MW_substrate / MW_biomass

        # Substrate atoms
        Cs = float(substrate_atoms.get("C", 0.0))
        Hs = float(substrate_atoms.get("H", 0.0))
        Os = float(substrate_atoms.get("O", 0.0))

        # Biomass atoms (per formula unit)
        Cx = float(biomass_atoms.get("C", 0.0))
        Hx = float(biomass_atoms.get("H", 0.0))
        Ox = float(biomass_atoms.get("O", 0.0))

        # Per 1 mol substrate consumed (nS = 1):
        # Carbon balance:  Cs = Y_mol * Cx + nCO2
        nCO2 = Cs - Y_mol * Cx

        if "N" in balance_upper:
            # CHNO mode
            Ns = float(substrate_atoms.get("N", 0.0))
            Nx = float(biomass_atoms.get("N", 0.0))

            if n_source_atoms is None:
                n_source_atoms = {"N": 1, "H": 3}
            H_nso = float(n_source_atoms.get("H", 0.0))
            N_nso = float(n_source_atoms.get("N", 0.0))
            O_nso = float(n_source_atoms.get("O", 0.0))

            # Nitrogen balance: Ns + N_nso * nNH3 = Y_mol * Nx
            if N_nso > 0.0:
                nNsource = (Y_mol * Nx - Ns) / N_nso
            else:
                nNsource = 0.0

            # Hydrogen balance: Hs + H_nso * nNsource = Y_mol * Hx + 2 * nH2O
            nH2O = (Hs + H_nso * nNsource - Y_mol * Hx) / 2.0

            # Oxygen balance: Os + O_nso * nNsource + 2*nO2 = Y_mol * Ox + 2*nCO2 + nH2O
            nO2 = (Y_mol * Ox + 2.0 * nCO2 + nH2O - Os - O_nso * nNsource) / 2.0

            _ov = species_overrides or {}
            sub_sp  = _ov.get(substrate_id) or Species(id=substrate_id, atoms=dict(substrate_atoms), MW=MW_substrate)
            o2_sp   = Species(id="O2", atoms={"O": 2})
            nsrc_sp = _ov.get(n_source_id) or Species(id=n_source_id, atoms=dict(n_source_atoms))
            bio_sp  = _ov.get(biomass_id) or Species(id=biomass_id, atoms=dict(biomass_atoms), MW=MW_biomass)
            co2_sp  = _CO2
            h2o_sp  = _H2O
            entries = [
                StoichiometryEntry(species=sub_sp, phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=o2_sp,  phase="liquid", coefficient=-nO2),
                StoichiometryEntry(species=nsrc_sp, phase="liquid", coefficient=-nNsource),
                StoichiometryEntry(species=bio_sp,  phase="liquid", coefficient=Y_mol),
                StoichiometryEntry(species=co2_sp,  phase="liquid", coefficient=nCO2),
                StoichiometryEntry(species=h2o_sp,  phase="liquid", coefficient=nH2O),
            ]
        else:
            # CHO mode
            # Hydrogen balance: Hs = Y_mol * Hx + 2 * nH2O
            nH2O = (Hs - Y_mol * Hx) / 2.0

            # Oxygen balance: Os + 2*nO2 = Y_mol * Ox + 2*nCO2 + nH2O
            nO2 = (Y_mol * Ox + 2.0 * nCO2 + nH2O - Os) / 2.0

            _ov = species_overrides or {}
            sub_sp = _ov.get(substrate_id) or Species(id=substrate_id, atoms=dict(substrate_atoms), MW=MW_substrate)
            o2_sp  = Species(id="O2", atoms={"O": 2})
            bio_sp = _ov.get(biomass_id) or Species(id=biomass_id, atoms=dict(biomass_atoms), MW=MW_biomass)
            co2_sp = _CO2
            h2o_sp = _H2O
            entries = [
                StoichiometryEntry(species=sub_sp, phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=o2_sp,  phase="liquid", coefficient=-nO2),
                StoichiometryEntry(species=bio_sp, phase="liquid", coefficient=Y_mol),
                StoichiometryEntry(species=co2_sp, phase="liquid", coefficient=nCO2),
                StoichiometryEntry(species=h2o_sp, phase="liquid", coefficient=nH2O),
            ]

        if not label:
            label = f"aerobic_{substrate_id}_{biomass_id}"

        return KineticReaction(
            stoichiometry=entries,
            rate_fn=rate_fn,
            balance_elements=balance_elements,
            label=label,
        )

    @staticmethod
    def monod_aerobic_growth(
        substrate: "Species",
        biomass: "Species",
        mu_max_per_h: float,
        Ks_gL: float,
        yield_gX_gS: float,
        *,
        Ko2_gL: Optional[float] = None,
        balance: str = "CHO",
        label: str = "",
    ) -> KineticReaction:
        """Aerobic growth with Monod substrate kinetics.

        Convenience wrapper around :meth:`aerobic_growth` that takes
        :class:`~PyOMES.chemistry.species.Species` objects directly and
        builds the Monod rate closure internally, so stoichiometry,
        elemental balance, and kinetics are all derived from the same
        minimal set of parameters.

        Rate law (extensive, mol substrate consumed / h):

        .. math::

            r = \\frac{\\mu_{\\max} \\cdot S_{g/L}}{K_S + S_{g/L}}
                \\cdot \\frac{C_{O_2}}{K_{O_2} + C_{O_2}}
                \\cdot \\frac{X_{g/L}}{Y} \\cdot \\frac{V_L}{\\text{MW}_S}

        The O₂ Monod term is omitted when ``Ko2_gL`` is ``None``
        (backward-compatible default).

        Parameters
        ----------
        substrate : Species
            Substrate species.  ``id``, ``atoms``, and ``MW`` are read
            directly — no need to pass them separately.
        biomass : Species
            Biomass species.
        mu_max_per_h : float
            Maximum specific growth rate (h⁻¹).
        Ks_gL : float
            Substrate half-saturation constant (g/L).
        yield_gX_gS : float
            Yield coefficient (g biomass / g substrate consumed).
        Ko2_gL : float or None
            O₂ half-saturation constant (g/L).  When provided, multiplies
            µ by ``O2_gL / (Ko2_gL + O2_gL)`` where ``O2_gL`` is derived
            from ``env.concentrations["O2"] * 32.0``.  Pass ``None``
            (default) to omit the O₂ Monod term entirely.
        balance : str
            Element set for stoichiometric closure: ``"CHO"`` (default)
            or ``"CHNO"``.
        label : str
            Human-readable label for the reaction.  Defaults to
            ``"monod_growth_{substrate.id}"``.

        Returns
        -------
        KineticReaction

        Example
        -------
        >>> rxn = ReactionBuilder.monod_aerobic_growth(
        ...     substrate    = ACETIC_ACID,
        ...     biomass      = YEAST,
        ...     mu_max_per_h = 0.5,
        ...     Ks_gL        = 5e-3,
        ...     yield_gX_gS  = 0.36,
        ...     Ko2_gL       = 0.2e-3,
        ...     label        = "growth_on_AceticAcid",
        ... )
        """
        MW_S = float(substrate.MW)
        MW_X = float(biomass.MW)
        _MW_O2 = 32.0

        def _rate_fn(env: ReactionEnvironment) -> float:
            S_gL = env.concentrations.get(substrate.id, 0.0) * MW_S
            X_gL = env.concentrations.get(biomass.id, 0.0) * MW_X
            if X_gL <= 1e-30 or S_gL <= 0.0:
                return 0.0
            mu = mu_max_per_h * S_gL / (Ks_gL + S_gL)
            if Ko2_gL is not None:
                O2_gL = env.concentrations.get("O2", 0.0) * _MW_O2
                mu *= O2_gL / (Ko2_gL + O2_gL)
            return mu / yield_gX_gS * X_gL / MW_S * env.V_L

        return ReactionBuilder.aerobic_growth(
            substrate_id    = substrate.id,
            substrate_atoms = substrate.atoms,
            MW_substrate    = MW_S,
            biomass_id      = biomass.id,
            biomass_atoms   = biomass.atoms,
            MW_biomass      = MW_X,
            yield_gX_gS     = yield_gX_gS,
            rate_fn         = _rate_fn,
            balance         = balance,
            label           = label or f"monod_growth_{substrate.id}",
            species_overrides={substrate.id: substrate, biomass.id: biomass},
        )

    @staticmethod
    def from_coefficients(
        entries: Sequence[StoichiometryEntry],
        rate_fn: Callable[[ReactionEnvironment], float],
        *,
        balance_elements: Sequence[str] = ("C", "H", "O"),
        balance_atol: float = 1e-10,
        label: str = "",
    ) -> KineticReaction:
        """Build a kinetic reaction from explicit stoichiometric coefficients.

        This is the most general factory: the user provides all entries
        with pre-computed coefficients. Elemental balance is validated
        at construction.

        Parameters
        ----------
        entries : sequence of StoichiometryEntry
            All reaction participants with coefficients.
        rate_fn : callable
            Rate law: ``rate_fn(env) -> float`` (mol/h).
        balance_elements : sequence of str
            Elements to validate.
        balance_atol : float
            Tolerance for elemental residual.
        label : str
            Human-readable label.

        Returns
        -------
        KineticReaction
        """
        return KineticReaction(
            stoichiometry=entries,
            rate_fn=rate_fn,
            balance_elements=balance_elements,
            balance_atol=balance_atol,
            label=label,
        )
