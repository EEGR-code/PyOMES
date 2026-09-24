# -*- coding: utf-8 -*-
"""Tests for PyOMES.reactions.rate_laws (checkpoint 12, decisions D5/D6).

This module had zero test coverage before this phase (checkpoint-1
audit). Covers:

- Fixed-point values for all eight laws' ``mu()``.
- ``Andrews`` and ``ContoisAndrews`` pinned against a frozen Haldane-form
  reference, since they are the reference for a later composable-
  inhibition redesign (see the cleanup plan's D5/D6 note).
- ``make_rate_fn`` extensive-rate (mol/h) outputs for ``Monod`` and
  ``DualSubstrateMonod``.
- A frozen-reference fingerprint proving
  ``ReactionBuilder.monod_aerobic_growth`` (now delegating to ``Monod``
  / ``DualSubstrateMonod``) is bit-identical to its pre-refactor inline
  closures, on both the plain-Monod and the ``Ko2_gL`` (O2) paths, and
  the one deliberate edge-case change: ``Ko2_gL=0.0`` with zero O2
  present now returns 0.0 instead of raising ``ZeroDivisionError``.
"""
from __future__ import annotations

import math
import random

import pytest

from PyOMES.reactions import (
    Monod, Contois, Andrews, ContoisAndrews, Tessier, Moser, Blackman,
    DualSubstrateMonod,
)


class _Env:
    """Minimal stand-in for ReactionEnvironment: concentrations + V_L."""

    def __init__(self, concentrations, V_L):
        self.concentrations = concentrations
        self.V_L = V_L


# ════════════════════════════════════════════════════════════════════════
#  Fixed-point values for all eight laws
# ════════════════════════════════════════════════════════════════════════

class TestMonod:
    def test_half_saturation_point(self):
        # S == Ks -> mu == mu_max / 2
        kin = Monod(mu_max=0.5, Ks=0.01)
        assert kin.mu(0.01, 1.0) == pytest.approx(0.25, rel=1e-12)

    def test_fixed_point(self):
        kin = Monod(mu_max=0.5, Ks=0.01)
        assert kin.mu(0.02, 1.0) == pytest.approx(1.0 / 3.0, rel=1e-12)

    def test_zero_substrate(self):
        assert Monod(mu_max=0.5, Ks=0.01).mu(0.0, 1.0) == 0.0

    def test_label(self):
        assert Monod().label == "Monod"


class TestContois:
    def test_fixed_point(self):
        kin = Contois(mu_max=0.5, Ks=0.1)
        assert kin.mu(0.02, 1.0) == pytest.approx(0.08333333333333333, rel=1e-12)

    def test_zero_biomass_returns_zero(self):
        # Density-dependent: no biomass -> S/X undefined, guarded to 0.
        assert Contois(mu_max=0.5, Ks=0.1).mu(0.02, 0.0) == 0.0

    def test_label(self):
        assert Contois().label == "Contois"


class TestTessier:
    def test_fixed_point(self):
        kin = Tessier(mu_max=0.5, Ks=0.01)
        expected = 0.5 * (1.0 - math.exp(-0.02 / 0.01))
        assert kin.mu(0.02, 1.0) == pytest.approx(expected, rel=1e-12)

    def test_large_S_approaches_mu_max(self):
        kin = Tessier(mu_max=0.5, Ks=0.01)
        assert kin.mu(5.0, 1.0) == pytest.approx(0.5, rel=1e-6)

    def test_label(self):
        assert Tessier().label == "Tessier"


class TestMoser:
    def test_reduces_to_monod_when_n_is_1(self):
        S_gL, X_gL = 0.02, 1.0
        moser = Moser(mu_max=0.5, Ks=0.01, n=1.0)
        monod = Monod(mu_max=0.5, Ks=0.01)
        assert moser.mu(S_gL, X_gL) == pytest.approx(monod.mu(S_gL, X_gL), rel=1e-12)

    def test_fixed_point_n2(self):
        kin = Moser(mu_max=0.5, Ks=0.0001, n=2.0)
        assert kin.mu(0.02, 1.0) == pytest.approx(0.4, rel=1e-12)

    def test_label(self):
        assert Moser().label == "Moser"


class TestBlackman:
    def test_half_saturation_point(self):
        # S == Ks -> mu == mu_max / 2 (same convention as Monod).
        kin = Blackman(mu_max=0.5, Ks=0.01)
        assert kin.mu(0.01, 1.0) == pytest.approx(0.25, rel=1e-12)

    def test_saturates_above_2Ks(self):
        kin = Blackman(mu_max=0.5, Ks=0.01)
        assert kin.mu(0.02, 1.0) == pytest.approx(0.5, rel=1e-12)
        assert kin.mu(5.0, 1.0) == pytest.approx(0.5, rel=1e-12)

    def test_label(self):
        assert Blackman().label == "Blackman"


# ════════════════════════════════════════════════════════════════════════
#  Andrews / ContoisAndrews — pinned as the Haldane-form reference for a
#  later composable-inhibition redesign (see the cleanup plan's D5/D6 note)
# ════════════════════════════════════════════════════════════════════════

def _haldane_mu(mu_max: float, Ks: float, Ki: float, r: float) -> float:
    """Frozen reference: the Haldane form both Andrews and ContoisAndrews use.

    μ = μ_max × r / (Ks + r + r²/Ki), where r is S (Andrews) or S/X
    (ContoisAndrews). Kept here, independent of rate_laws.py, as the
    fixed target for the later composable-inhibition design (an
    optional Ki on the base law) to reproduce exactly.
    """
    denom = Ks + r + r * r / Ki
    if denom <= 0.0:
        return 0.0
    return mu_max * r / denom


class TestAndrewsHaldaneForm:
    PARAMS = [
        (0.5, 0.01, 1.0),
        (0.8, 0.02, 50.0),
        (0.3, 0.005, 0.5),
    ]
    S_VALUES = [0.0, 1e-4, 0.01, 0.02, 0.5, 5.0, 50.0]

    def test_matches_haldane_reference(self):
        for mu_max, Ks, Ki in self.PARAMS:
            kin = Andrews(mu_max=mu_max, Ks=Ks, Ki=Ki)
            for S in self.S_VALUES:
                assert kin.mu(S, 1.0) == pytest.approx(
                    _haldane_mu(mu_max, Ks, Ki, S), rel=1e-12
                )

    def test_fixed_point(self):
        kin = Andrews(mu_max=0.5, Ks=0.01, Ki=1.0)
        assert kin.mu(0.02, 1.0) == pytest.approx(0.32894736842105265, rel=1e-12)

    def test_peaks_below_mu_max_for_finite_Ki(self):
        # Substrate inhibition: mu at the analytical peak S* = sqrt(Ks*Ki)
        # is strictly less than mu_max.
        mu_max, Ks, Ki = 0.5, 0.01, 1.0
        kin = Andrews(mu_max=mu_max, Ks=Ks, Ki=Ki)
        S_peak = math.sqrt(Ks * Ki)
        assert kin.mu(S_peak, 1.0) < mu_max

    def test_independent_of_biomass(self):
        # Andrews is substrate-only: X_gL does not appear in the formula.
        kin = Andrews(mu_max=0.5, Ks=0.01, Ki=1.0)
        assert kin.mu(0.02, 1.0) == kin.mu(0.02, 100.0)

    def test_label(self):
        assert Andrews().label == "Andrews"


class TestContoisAndrewsHaldaneForm:
    PARAMS = [
        (0.5, 0.1, 1.0),
        (0.8, 0.05, 20.0),
    ]
    RATIOS = [0.0, 1e-4, 0.02, 0.1, 1.0, 10.0]

    def test_matches_haldane_reference_via_ratio(self):
        for mu_max, Ks, Ki in self.PARAMS:
            kin = ContoisAndrews(mu_max=mu_max, Ks=Ks, Ki=Ki)
            for r in self.RATIOS:
                X_gL = 1.0
                S_gL = r * X_gL
                assert kin.mu(S_gL, X_gL) == pytest.approx(
                    _haldane_mu(mu_max, Ks, Ki, r), rel=1e-12
                )

    def test_fixed_point(self):
        kin = ContoisAndrews(mu_max=0.5, Ks=0.1, Ki=1.0)
        assert kin.mu(0.02, 1.0) == pytest.approx(0.08305647840531562, rel=1e-12)

    def test_zero_biomass_returns_zero(self):
        assert ContoisAndrews(mu_max=0.5, Ks=0.1, Ki=1.0).mu(0.02, 0.0) == 0.0

    def test_scale_invariant_in_S_and_X(self):
        # Density-dependent: only the ratio S/X matters.
        kin = ContoisAndrews(mu_max=0.5, Ks=0.1, Ki=1.0)
        assert kin.mu(0.02, 1.0) == pytest.approx(kin.mu(0.2, 10.0), rel=1e-12)

    def test_label(self):
        assert ContoisAndrews().label == "ContoisAndrews"


# ════════════════════════════════════════════════════════════════════════
#  DualSubstrateMonod
# ════════════════════════════════════════════════════════════════════════

class TestDualSubstrateMonod:
    def test_make_rate_fn_extensive_rate(self):
        kin = DualSubstrateMonod(mu_max=0.5, Ks=0.01, secondary_id="O2",
                                  Ko=1e-4, secondary_in_mol_L=True)
        rate_fn = kin.make_rate_fn(organism_id="X", substrate_id="S",
                                    MW_organism=24.6, MW_substrate=180.0,
                                    yield_gX_gS=0.5)
        env = _Env({"S": 0.001, "X": 0.01, "O2": 2e-4}, 10.0)
        assert rate_fn(env) == pytest.approx(0.00863157894736842, rel=1e-12)

    def test_zero_Ko_with_zero_secondary_returns_zero(self):
        """The D5 edge case: Ko=0 and zero secondary substrate present."""
        kin = DualSubstrateMonod(mu_max=0.5, Ks=0.01, secondary_id="O2",
                                  Ko=0.0, secondary_in_mol_L=True)
        rate_fn = kin.make_rate_fn(organism_id="X", substrate_id="S",
                                    MW_organism=24.6, MW_substrate=180.0,
                                    yield_gX_gS=0.5)
        env = _Env({"S": 0.001, "X": 0.01, "O2": 0.0}, 10.0)
        assert rate_fn(env) == 0.0

    def test_secondary_in_gL_conversion(self):
        # secondary_in_mol_L=False: secondary read in mol/L then * MW.
        kin_gl = DualSubstrateMonod(mu_max=0.5, Ks=0.01, secondary_id="O2",
                                     Ko=1e-4, secondary_in_mol_L=False,
                                     secondary_MW=32.0)
        rate_fn = kin_gl.make_rate_fn(organism_id="X", substrate_id="S",
                                       MW_organism=24.6, MW_substrate=180.0,
                                       yield_gX_gS=0.5)
        # 2e-4 mol/L O2 -> converted to g/L via *32.0 -> different Ko regime
        # than the mol/L-direct case above; just check it runs and is finite.
        env = _Env({"S": 0.001, "X": 0.01, "O2": 2e-4}, 10.0)
        result = rate_fn(env)
        assert result > 0.0

    def test_label_includes_secondary_id(self):
        assert DualSubstrateMonod(secondary_id="O2").label == "DualSubstrateMonod(+O2)"


# ════════════════════════════════════════════════════════════════════════
#  make_rate_fn — extensive rate (mol/h) for the shared base class
# ════════════════════════════════════════════════════════════════════════

class TestMakeRateFnOutputs:
    def test_monod_extensive_rate(self):
        kin = Monod(mu_max=0.5, Ks=0.01)
        rate_fn = kin.make_rate_fn(organism_id="X", substrate_id="S",
                                    MW_organism=24.6, MW_substrate=180.0,
                                    yield_gX_gS=0.5)
        env = _Env({"S": 0.001, "X": 0.01}, 10.0)
        assert rate_fn(env) == pytest.approx(0.012947368421052633, rel=1e-12)

    def test_zero_biomass_gives_zero_rate(self):
        kin = Monod(mu_max=0.5, Ks=0.01)
        rate_fn = kin.make_rate_fn(organism_id="X", substrate_id="S",
                                    MW_organism=24.6, MW_substrate=180.0,
                                    yield_gX_gS=0.5)
        env = _Env({"S": 0.001, "X": 0.0}, 10.0)
        assert rate_fn(env) == 0.0

    def test_zero_substrate_gives_zero_rate(self):
        kin = Monod(mu_max=0.5, Ks=0.01)
        rate_fn = kin.make_rate_fn(organism_id="X", substrate_id="S",
                                    MW_organism=24.6, MW_substrate=180.0,
                                    yield_gX_gS=0.5)
        env = _Env({"S": 0.0, "X": 0.01}, 10.0)
        assert rate_fn(env) == 0.0


# ════════════════════════════════════════════════════════════════════════
#  Fingerprint: ReactionBuilder.monod_aerobic_growth now delegates to
#  Monod / DualSubstrateMonod. Frozen reference reproduces the pre-move
#  inline closures exactly (checkpoint-1 audit: bit-identical on 200,000
#  log-uniform physical points and a 144-point edge grid).
# ════════════════════════════════════════════════════════════════════════

def _reference_monod_rate_fn(mu_max_per_h, Ks_gL, yield_gX_gS, MW_S, MW_X,
                              sub_id, bio_id, Ko2_gL=None, MW_O2=32.0):
    """Frozen copy of the pre-checkpoint-12 inline closure in
    ReactionBuilder.monod_aerobic_growth (builder.py, formerly line 298)."""

    def rate_fn(env):
        S_gL = env.concentrations.get(sub_id, 0.0) * MW_S
        X_gL = env.concentrations.get(bio_id, 0.0) * MW_X
        if X_gL <= 1e-30 or S_gL <= 0.0:
            return 0.0
        mu = mu_max_per_h * S_gL / (Ks_gL + S_gL)
        if Ko2_gL is not None:
            O2_gL = env.concentrations.get("O2", 0.0) * MW_O2
            mu *= O2_gL / (Ko2_gL + O2_gL)  # ZeroDivisionError at Ko2_gL=O2_gL=0
        return mu / yield_gX_gS * X_gL / MW_S * env.V_L

    return rate_fn


class TestMonodAerobicGrowthFingerprint:
    """Grid comparison of ReactionBuilder.monod_aerobic_growth's current
    (delegating) rate_fn against the frozen pre-refactor reference."""

    MU_MAX, KS_GL, YIELD, MW_S, MW_X = 0.6, 0.01, 0.45, 180.156, 24.6

    @staticmethod
    def _grid(n=500, seed=12345):
        rng = random.Random(seed)
        out = []
        for _ in range(n):
            S = 10 ** rng.uniform(-6, 0)
            X = 10 ** rng.uniform(-6, 0)
            O2 = 10 ** rng.uniform(-8, -2)
            V = 10 ** rng.uniform(-1, 3)
            out.append((S, X, O2, V))
        return out

    def _current_rate_fn(self, Ko2_gL):
        from PyOMES.reactions import ReactionBuilder
        from PyOMES.chemistry.species import Species
        sub = Species(id="Glucose", atoms={"C": 6, "H": 12, "O": 6}, charge=0, MW=self.MW_S)
        bio = Species(id="Ecoli", atoms={"C": 1, "H": 1.8, "O": 0.5, "N": 0.2}, charge=0, MW=self.MW_X)
        rxn = ReactionBuilder.monod_aerobic_growth(
            substrate=sub, biomass=bio, mu_max_per_h=self.MU_MAX, Ks_gL=self.KS_GL,
            yield_gX_gS=self.YIELD, Ko2_gL=Ko2_gL,
        )
        return rxn.rate_fn, sub.id, bio.id

    def test_plain_monod_path_bit_identical(self):
        rate_fn, sub_id, bio_id = self._current_rate_fn(Ko2_gL=None)
        ref_fn = _reference_monod_rate_fn(
            self.MU_MAX, self.KS_GL, self.YIELD, self.MW_S, self.MW_X,
            sub_id, bio_id, Ko2_gL=None,
        )
        for S, X, O2, V in self._grid():
            env = _Env({sub_id: S, bio_id: X, "O2": O2}, V)
            assert rate_fn(env) == ref_fn(env)

    def test_o2_path_bit_identical(self):
        Ko2_gL = 0.2e-3
        rate_fn, sub_id, bio_id = self._current_rate_fn(Ko2_gL=Ko2_gL)
        ref_fn = _reference_monod_rate_fn(
            self.MU_MAX, self.KS_GL, self.YIELD, self.MW_S, self.MW_X,
            sub_id, bio_id, Ko2_gL=Ko2_gL,
        )
        for S, X, O2, V in self._grid():
            env = _Env({sub_id: S, bio_id: X, "O2": O2}, V)
            assert rate_fn(env) == ref_fn(env)

    def test_edge_case_zero_Ko2_zero_O2(self):
        """D5: the one intentional behaviour change. The old inline closure
        (see _reference_monod_rate_fn) raises ZeroDivisionError here; the
        current DualSubstrateMonod-backed implementation returns 0.0."""
        rate_fn, sub_id, bio_id = self._current_rate_fn(Ko2_gL=0.0)
        env = _Env({sub_id: 0.01, bio_id: 0.01, "O2": 0.0}, 1.0)
        assert rate_fn(env) == 0.0

        ref_fn = _reference_monod_rate_fn(
            self.MU_MAX, self.KS_GL, self.YIELD, self.MW_S, self.MW_X,
            sub_id, bio_id, Ko2_gL=0.0,
        )
        with pytest.raises(ZeroDivisionError):
            ref_fn(env)
