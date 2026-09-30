# -*- coding: utf-8 -*-
"""Convenience builders for constructing validated KineticReaction objects.

:class:`ReactionBuilder` provides static methods that derive complete
stoichiometric coefficients from high-level specifications (substrate
formula, biomass formula, yield, balance mode) and produce validated
:class:`~PyOMES.reactions.kinetic.reaction.KineticReaction` objects.

Every participant is a :class:`~PyOMES.chemistry.species.Species` the
caller passes in: the substrate, the biomass, O₂, CO₂ and H₂O, and, for a
CHNO balance, the nitrogen source.

High-level convenience method (Monod kinetics + stoichiometry in one call):

>>> rxn = ReactionBuilder.monod_aerobic_growth(
...     substrate    = ACETIC_ACID,   # Species objects throughout
...     biomass      = YEAST,
...     mu_max_per_h = 0.5,
...     Ks_gL        = 5e-3,
...     yield_gX_gS  = 0.36,
...     o2=O2, co2=CO2, h2o=H2O,
... )

Lower-level method for custom rate functions:

>>> rxn = ReactionBuilder.aerobic_growth(
...     substrate   = ACETIC_ACID,
...     biomass     = YEAST,
...     o2=O2, co2=CO2, h2o=H2O,
...     yield_gX_gS = 0.36,
...     rate_fn     = my_custom_rate,
... )
"""

from __future__ import annotations

from typing import Callable, Optional, Sequence

from PyOMES.reactions.stoichiometry import StoichiometryEntry
from .reaction import KineticReaction
from PyOMES.reactions.environment import ReactionEnvironment
from .rate_laws import Monod, DualSubstrateMonod
from PyOMES.chemistry.species import Species


class ReactionBuilder:
    """Static factory methods for constructing validated reactions."""

    @staticmethod
    def aerobic_growth(
        substrate: Species,
        biomass: Species,
        *,
        o2: Species,
        co2: Species,
        h2o: Species,
        yield_gX_gS: float,
        rate_fn: Callable[[ReactionEnvironment], float],
        balance: str = "CHO",
        n_source: Optional[Species] = None,
        label: str = "",
    ) -> KineticReaction:
        """Build an aerobic growth reaction from a yield and elemental formulas.

        Derives the stoichiometric coefficients for O₂ (consumed),
        CO₂ (produced), and H₂O (produced) from elemental balance.
        In CHNO mode, also derives the nitrogen source demand.

        The reaction basis is **1 mol of substrate consumed** (coefficient
        = −1).  The rate function should therefore return the substrate
        consumption rate in mol/h (extensive).

        Every participant is the :class:`~PyOMES.chemistry.species.Species`
        passed in; ids, atoms and molar masses are read from them. The
        balance treats ``o2``, ``co2`` and ``h2o`` as O₂, CO₂ and H₂O; the
        reaction's own element-balance check rejects species whose atoms
        differ.

        Parameters
        ----------
        substrate : Species
            The substrate consumed (e.g. acetic acid).
        biomass : Species
            The biomass formed, one formula unit (e.g. CH₁.₆₁O₀.₅₆).
        o2, co2, h2o : Species
            Oxygen consumed, and carbon dioxide and water produced.
        yield_gX_gS : float
            Yield coefficient in g biomass per g substrate consumed.
        rate_fn : callable
            ``rate_fn(env) -> float`` returning the extensive substrate
            consumption rate (mol substrate / h).
        balance : str
            ``"CHO"`` (default) or ``"CHNO"``.  When ``"CHNO"``, the
            nitrogen source is included in the stoichiometry.
        n_source : Species or None
            The nitrogen source consumed. Required when ``balance``
            includes ``"N"``; ignored otherwise.
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
            If inputs are non-physical (zero MW, negative yield, etc.), or
            ``balance`` includes ``"N"`` and no ``n_source`` is given.
        """
        MW_substrate = float(substrate.MW)
        MW_biomass = float(biomass.MW)
        substrate_atoms = substrate.atoms
        biomass_atoms = biomass.atoms
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

            if n_source is None:
                raise ValueError(
                    f"balance={balance!r} includes N, so a nitrogen source is "
                    "needed: pass n_source=<Species>, e.g. NH3."
                )
            n_source_atoms = n_source.atoms
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

            entries = [
                StoichiometryEntry(species=substrate, phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=o2,        phase="liquid", coefficient=-nO2),
                StoichiometryEntry(species=n_source,  phase="liquid", coefficient=-nNsource),
                StoichiometryEntry(species=biomass,   phase="liquid", coefficient=Y_mol),
                StoichiometryEntry(species=co2,       phase="liquid", coefficient=nCO2),
                StoichiometryEntry(species=h2o,       phase="liquid", coefficient=nH2O),
            ]
        else:
            # CHO mode
            # Hydrogen balance: Hs = Y_mol * Hx + 2 * nH2O
            nH2O = (Hs - Y_mol * Hx) / 2.0

            # Oxygen balance: Os + 2*nO2 = Y_mol * Ox + 2*nCO2 + nH2O
            nO2 = (Y_mol * Ox + 2.0 * nCO2 + nH2O - Os) / 2.0

            entries = [
                StoichiometryEntry(species=substrate, phase="liquid", coefficient=-1.0),
                StoichiometryEntry(species=o2,        phase="liquid", coefficient=-nO2),
                StoichiometryEntry(species=biomass,   phase="liquid", coefficient=Y_mol),
                StoichiometryEntry(species=co2,       phase="liquid", coefficient=nCO2),
                StoichiometryEntry(species=h2o,       phase="liquid", coefficient=nH2O),
            ]

        if not label:
            label = f"aerobic_{substrate.id}_{biomass.id}"

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
        o2: Species,
        co2: Species,
        h2o: Species,
        Ko2_gL: Optional[float] = None,
        balance: str = "CHO",
        n_source: Optional[Species] = None,
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
        o2, co2, h2o : Species
            Oxygen consumed, and carbon dioxide and water produced; passed
            to :meth:`aerobic_growth`.
        Ko2_gL : float or None
            O₂ half-saturation constant (g/L).  When provided, multiplies
            µ by ``O2_gL / (Ko2_gL + O2_gL)`` where ``O2_gL`` is derived
            from ``env.concentrations[o2.id] * o2.MW``.  Pass ``None``
            (default) to omit the O₂ Monod term entirely.  ``Ko2_gL=0.0``
            with zero O₂ present returns a rate of 0.0 rather than raising
            (the 0/0 case is guarded).
        balance : str
            Element set for stoichiometric closure: ``"CHO"`` (default)
            or ``"CHNO"``.
        n_source : Species or None
            Nitrogen source; required when ``balance`` includes ``"N"``.
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
        ...     o2=O2, co2=CO2, h2o=H2O,
        ...     Ko2_gL       = 0.2e-3,
        ...     label        = "growth_on_AceticAcid",
        ... )
        """
        MW_S = float(substrate.MW)
        MW_X = float(biomass.MW)

        if Ko2_gL is None:
            kin = Monod(mu_max=mu_max_per_h, Ks=Ks_gL)
        else:
            kin = DualSubstrateMonod(
                mu_max=mu_max_per_h, Ks=Ks_gL,
                secondary_id=o2.id, Ko=Ko2_gL,
                secondary_in_mol_L=False, secondary_MW=float(o2.MW),
            )
        rate_fn = kin.make_rate_fn(
            organism_id=biomass.id,
            substrate_id=substrate.id,
            MW_organism=MW_X,
            MW_substrate=MW_S,
            yield_gX_gS=yield_gX_gS,
        )

        return ReactionBuilder.aerobic_growth(
            substrate, biomass,
            o2=o2, co2=co2, h2o=h2o,
            yield_gX_gS = yield_gX_gS,
            rate_fn     = rate_fn,
            balance     = balance,
            n_source    = n_source,
            label       = label or f"monod_growth_{substrate.id}",
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
