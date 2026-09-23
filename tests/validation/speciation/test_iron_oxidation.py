# -*- coding: utf-8 -*-
"""Coverage for 07_iron_oxidation.ipynb / 08_iron_oxidation_and_precipitation.ipynb.

Both notebooks previously reported Fe2+ conversion of 0.2% against a stated
target of >70% (P1's own [CHECK] in its saved output). Root cause: the
Singer & Stumm (1970) rate law is calibrated against p(O2) in atm, but
`ControlVolume._build_reaction_environment` only ever populates a kinetic
`rate_fn`'s `env.concentrations` from the *liquid* phase's `n_mol` -- gas-
phase state is never merged in (confirmed by reading both call sites,
control_volume.py:700-748 and :924-1003). The old notebook's K_SS=5e11 was
that literature k applied directly to aqueous [O2] in mol/L with no
unit conversion, understating the rate by roughly the Henry's-law factor
(kH_O2 ~ 1.3e-3 mol/(L.atm), i.e. ~760x) -- regardless of pH, since the
[OH-]**2 term alone can't close a gap that large. See
docs/dev/implementation/upcoming/REACTION_ENVIRONMENT_PHASE_EXPOSURE.md for
the follow-up idea (exposing gas-phase partial pressure to rate_fn directly,
so this conversion wouldn't need to happen by hand) logged as a result of
this fix.

The corrected model:
- K_SS is derived from a verified literature constant (k = 1.33e12
  M^-2 atm^-1 s^-1, Singer & Stumm 1970 Eq. 22, cited via USGS's PHREEQC
  docs Example 9) via K_SS = k / kH_O2, applied against aqueous [O2].
- O2 is supplied by a real GasPhase (1000 L air headspace) via an
  EquilibriumTransferModel, not a fixed depletable pool -- same
  "pseudo-unlimited reservoir" idiom as
  docs/tutorials/ArXiv_preprint/_generate_notebooks.py's kinetic CO2 demo.
- The batch run uses SimultaneousAdaptiveSolver(method="BDF",
  use_engine_jacobian=True) rather than the default explicit solver: this
  reaction is stiff on the timescale that matters (t1/2 ~ 11 min inside a
  2h run). Switching to BDF also eliminates the ConservationWarnings a
  fixed-step explicit solver accumulates from many large per-step
  Newton-Raphson speciation re-solves once the reaction is this fast --
  confirmed a genuine fix, not just a suppressed report (Fe mass-balance
  drift dropped ~10x). It does trip a *different*, confirmed-false-positive
  AccuracyWarning for implicit solvers -- see TestBatchOxidation's
  docstring and docs/dev/implementation/upcoming/
  SCIPY_REJECTION_CHECK_SOLVER_AWARENESS.md.
- O2(aq) starts pre-equilibrated with the headspace, not at 0 -- a second,
  independent bug caught by noticing an O2(aq) spike right after t=0 in the
  notebook's plot: starting at 0 and letting the EquilibriumTransferModel
  jump it to equilibrium on the solver's very first internal step is an
  artificial discontinuity, not a physical initial condition, and it alone
  roughly quadrupled the number of BDF steps an otherwise-identical run
  needed (~70 vs ~15 accepted steps).

`TestSingerStummRateLaw` is a fast, no-Simulation regression guard on the
rate law's pH-direction (higher pH -> faster oxidation, since r ~ [OH-]**2)
and 4th-order structure -- the relationship whose absence let the original
bug through unnoticed. `TestBatchOxidation` regression-tests the corrected,
literature-derived K_SS against the notebooks' own P1/P3/P4/P5 predictions.
"""
from __future__ import annotations

import pytest

from PyOMES.chemistry import Species
from PyOMES.reactions import HenryEquilibrium
from PyOMES.chemistry.common_species import (
    H2O, H_plus, OH_minus,
    H3PO4, H2PO4_minus, HPO4_2minus, PO4_3minus,
    K_plus,
)
from PyOMES.reactions import (
    KineticReaction, EquilibriumReaction, ReactionSystem, StoichiometryEntry,
)
from PyOMES.reactions.environment import ReactionEnvironment
from PyOMES.core import (
    LiquidPhase, GasPhase, ControlVolume, Simulation, EquilibriumTransferModel,
)
from PyOMES.units import R_L_ATM_PER_MOL_K
from PyOMES.core.solvers import SimultaneousAdaptiveSolver

# Singer & Stumm (1970) Eq. 22, p(O2) basis (M^-2 atm^-1 s^-1 -> per hour).
K_SS_LITERATURE_ATM = 1.33e12 * 3600
O2_HENRY = HenryEquilibrium(H_ref=1.3e-5, dlnH=1500.0)  # matches bioprocess_basic.py
K_SS = K_SS_LITERATURE_ATM / O2_HENRY._kH_mol_L_atm(298.15)

Fe2_plus = Species(id="Fe2+", atoms={"Fe": 1}, charge=+2)
Fe3_plus = Species(id="Fe3+", atoms={"Fe": 1}, charge=+3)
FeOH_2p = Species(id="FeOH2+", atoms={"Fe": 1, "O": 1, "H": 1}, charge=+2)
FeOH2_p = Species(id="Fe(OH)2+", atoms={"Fe": 1, "O": 2, "H": 2}, charge=+1)
O2_aq = Species(id="O2", atoms={"O": 2}, charge=0)


def _e(sp, coeff):
    return StoichiometryEntry(species=sp, phase="liquid", coefficient=coeff)


def _rate_fn(env):
    return K_SS * env.S("Fe2+") * env.S("O2") * env.S("OH-") ** 2 * env.V_L


def _make_reaction_system():
    rxn_water = EquilibriumReaction(
        stoichiometry=[_e(H2O, -1), _e(H_plus, +1), _e(OH_minus, +1)],
        log_K=-14.0, label="water",
    )
    rxn_p1 = EquilibriumReaction(
        stoichiometry=[_e(H3PO4, -1), _e(H2PO4_minus, +1), _e(H_plus, +1)],
        log_K=-2.15, total_id="H3PO4", label="P_pKa1",
    )
    rxn_p2 = EquilibriumReaction(
        stoichiometry=[_e(H2PO4_minus, -1), _e(HPO4_2minus, +1), _e(H_plus, +1)],
        log_K=-7.20, total_id="H3PO4", label="P_pKa2",
    )
    rxn_p3 = EquilibriumReaction(
        stoichiometry=[_e(HPO4_2minus, -1), _e(PO4_3minus, +1), _e(H_plus, +1)],
        log_K=-12.35, total_id="H3PO4", label="P_pKa3",
    )
    rxn_fe3_h1 = EquilibriumReaction(
        stoichiometry=[_e(Fe3_plus, -1), _e(H2O, -1), _e(FeOH_2p, +1), _e(H_plus, +1)],
        log_K=-2.19, total_id="Fe3+", label="Fe3_h1",
    )
    rxn_fe3_h2 = EquilibriumReaction(
        stoichiometry=[_e(FeOH_2p, -1), _e(H2O, -1), _e(FeOH2_p, +1), _e(H_plus, +1)],
        log_K=-3.48, total_id="Fe3+", label="Fe3_h2",
    )
    rxn_oxidation = KineticReaction(
        stoichiometry=[
            _e(Fe2_plus, -4), _e(O2_aq, -1), _e(H2O, -2),
            _e(Fe3_plus, +4), _e(OH_minus, +4),
        ],
        rate_fn=_rate_fn,
        label="Fe2_O2_oxidation",
    )
    return ReactionSystem(
        [rxn_oxidation, rxn_water, rxn_p1, rxn_p2, rxn_p3, rxn_fe3_h1, rxn_fe3_h2],
        label="Fe_O2_phosphate_medium_test",
        solver="newton_raphson",
    )


def _make_cv(pH_target, V_L=1.0, V_gas=1000.0, T_K=298.15):
    """Phosphate-buffered liquid (0.88 mM Fe2+, the notebooks' actual --
    not earlier-draft-claimed -- loading) + a large air headspace supplying
    O2 via Henry's-law equilibrium, equilibrated to pH_target.

    O2 starts pre-equilibrated with the headspace (not 0), matching the
    notebooks' own fix: a medium open to air has always been in contact
    with it, so starting at 0 mM produces an artificial discontinuity on
    the very first internal solver step (confirmed: this alone took an
    otherwise-identical BDF run from ~15 accepted steps to ~70)."""
    n_gas_total = (1.0 * V_gas) / (R_L_ATM_PER_MOL_K * T_K)
    gas = GasPhase(
        n_mol={"O2": n_gas_total * 0.2095, "N2": n_gas_total * 0.7905},
        V_L=V_gas, T_K=T_K,
    )
    n_o2_init = O2_HENRY._kH_mol_L_atm(T_K) * 0.2095 * V_L
    liquid = LiquidPhase(
        n_mol={
            K_plus.id: 0.2204,  # counter-ion for H2PO4- -- charge balance
            H2PO4_minus.id: 0.2204, H3PO4.id: 0.0, HPO4_2minus.id: 0.0, PO4_3minus.id: 0.0,
            "Fe2+": 8.78e-4, "Fe3+": 0.0, "FeOH2+": 0.0, "Fe(OH)2+": 0.0,
            H_plus.id: 2e-5, OH_minus.id: 0.0, H2O.id: 55.5,
            O2_aq.id: n_o2_init,
        },
        V_L=V_L, T_K=T_K,
    )
    cv = ControlVolume(
        phases={"gas": gas, "liquid": liquid},
        transfer_models={"O2": EquilibriumTransferModel(O2_HENRY)},
        reaction_system=_make_reaction_system(),
        label="test_batch_Fe_O2_oxidation",
    )
    cv.equilibrate_to_pH("KOH", pH_target)
    return cv


class TestSingerStummRateLaw:
    """Rate law direction check -- the relationship whose absence let the
    original bug (a rate constant tuned as if raising pH would slow the
    reaction) through unnoticed. r ~ [OH-]**2, so pH 6.5 must be *faster*
    than the medium's natural pH 4.68, not slower."""

    def test_rate_increases_with_pH(self):
        env_low = ReactionEnvironment(
            V_L=1.0, concentrations={"Fe2+": 8.78e-4, "O2": 2.76e-4, "OH-": 10 ** (4.68 - 14)},
        )
        env_high = ReactionEnvironment(
            V_L=1.0, concentrations={"Fe2+": 8.78e-4, "O2": 2.76e-4, "OH-": 10 ** (6.5 - 14)},
        )
        rate_low = _rate_fn(env_low)
        rate_high = _rate_fn(env_high)
        assert rate_high > rate_low
        # [OH-] scales 10**(6.5-4.68) ~= 66x; rate ~ [OH-]**2 -> ~4380x
        assert rate_high / rate_low == pytest.approx(10 ** (2 * 1.82), rel=0.05)

    def test_rate_law_is_fourth_order_overall(self):
        """1 (Fe2+) + 1 (O2) + 2 (OH-) = 4th order -- doubling any single
        concentration input should scale the rate by exactly that factor
        (Fe2+/O2, linear) or its square (OH-, quadratic)."""
        base = ReactionEnvironment(
            V_L=1.0, concentrations={"Fe2+": 1e-3, "O2": 1e-3, "OH-": 1e-8},
        )
        doubled_fe2 = ReactionEnvironment(
            V_L=1.0, concentrations={"Fe2+": 2e-3, "O2": 1e-3, "OH-": 1e-8},
        )
        doubled_oh = ReactionEnvironment(
            V_L=1.0, concentrations={"Fe2+": 1e-3, "O2": 1e-3, "OH-": 2e-8},
        )
        assert _rate_fn(doubled_fe2) / _rate_fn(base) == pytest.approx(2.0)
        assert _rate_fn(doubled_oh) / _rate_fn(base) == pytest.approx(4.0)

    def test_K_SS_matches_henrys_law_conversion_of_literature_constant(self):
        """K_SS = k_literature / kH_O2 -- not an independently chosen value."""
        kH_O2 = O2_HENRY._kH_mol_L_atm(298.15)
        assert K_SS * kH_O2 == pytest.approx(K_SS_LITERATURE_ATM, rel=1e-9)


@pytest.mark.filterwarnings("ignore::PyOMES.monitoring.accuracy.AccuracyWarning")
class TestBatchOxidation:
    """Regression coverage for the corrected, literature-derived K_SS against
    the notebooks' own P1 (>70% Fe2+ conversion), P3 (Fe(OH)2+ dominant), P4
    (Fe mass balance), and P5 (gas reservoir barely moves) predictions, at
    the notebooks' actual pH 6.5 operating point.

    BDF's AccuracyWarning (check_scipy_rejections) is a confirmed false
    positive for implicit solvers -- see the module docstring and
    docs/dev/implementation/upcoming/SCIPY_REJECTION_CHECK_SOLVER_AWARENESS.md."""

    def _run(self, tau_h=2.0, n_steps=200):
        # BDF (implicit, adaptive-step) rather than the default explicit
        # solver -- this reaction is stiff on the timescale that matters
        # (t1/2 ~ 11 min inside a 2h run). use_engine_jacobian=True since
        # this engine satisfies GrayBoxEngineProtocol. Matches the notebooks
        # (07/08_iron_oxidation*.ipynb) exactly.
        cv = _make_cv(pH_target=6.5)
        bdf_solver = SimultaneousAdaptiveSolver(method="BDF", use_engine_jacobian=True)
        sim = Simulation(cvs={"main": cv}, label="test_Fe_O2_oxidation", solver=bdf_solver)
        return cv, sim.run(tau_h=tau_h, n_steps=n_steps)

    def test_fe2_conversion_exceeds_70_percent_in_2h(self):
        _, result = self._run()
        liq = result.liquid_mol["main"]
        conv = 1 - liq["Fe2+"][-1] / liq["Fe2+"][0]
        assert conv > 0.70

    def test_fe_III_pool_dominated_by_second_hydrolysis_product(self):
        """At pH 6.5 (pKh2=3.48 << 6.5), Fe(OH)2+ dominates -- not FeOH2+,
        which is what the notebooks' pre-fix narrative (written for the
        medium's natural pH 4.68) claimed."""
        _, result = self._run()
        liq = result.liquid_mol["main"]
        fe3_total = liq["Fe3+"][-1] + liq["FeOH2+"][-1] + liq["Fe(OH)2+"][-1]
        assert fe3_total > 0
        assert liq["Fe(OH)2+"][-1] / fe3_total > 0.9
        assert liq["Fe(OH)2+"][-1] > liq["FeOH2+"][-1]

    def test_total_iron_conserved(self):
        _, result = self._run()
        liq = result.liquid_mol["main"]
        fe_total = liq["Fe2+"] + liq["Fe3+"] + liq["FeOH2+"] + liq["Fe(OH)2+"]
        drift = abs(fe_total - fe_total[0]).max() / fe_total[0]
        assert drift < 1e-3

    def test_gas_reservoir_stays_effectively_unlimited(self):
        """Confirms the 1000 L : 1 L headspace ratio is large enough that
        the O2 reservoir barely moves -- the load-bearing assumption behind
        treating it as 'pseudo-unlimited' rather than verifying it."""
        _, result = self._run()
        gas = result.gas_mol["main"]
        drawdown = 1 - gas["O2"][-1] / gas["O2"][0]
        assert drawdown < 0.01

    def test_pH_stays_near_operating_point(self):
        _, result = self._run()
        pH = result.pH["main"]
        pH_valid = pH[pH == pH]  # drop NaN
        assert abs(pH_valid[-1] - 6.5) < 0.5

    def test_O2_starts_pre_equilibrated_not_at_zero(self):
        """Regression guard: O2 must start at its Henry's-law equilibrium
        value, not 0 -- starting at 0 produces an artificial discontinuity
        on the solver's first internal step (a real bug caught by noticing
        an O2(aq) spike in the notebook's plot right after t=0)."""
        cv = _make_cv(pH_target=6.5)
        o2_initial_mM = cv.phases["liquid"].n_mol["O2"] / cv.phases["liquid"].V_L * 1e3
        expected_mM = O2_HENRY._kH_mol_L_atm(298.15) * 0.2095 * 1e3
        assert o2_initial_mM == pytest.approx(expected_mM, rel=1e-6)
        assert o2_initial_mM > 0.1  # not the old bug's 0.0
