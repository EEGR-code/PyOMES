"""Unit tests for the HPLC column chromatography model.

Tests construction, properties, inlet functions, species types,
pH-dependent retention, mass conservation, and simulation results.
"""

import pytest
import warnings
import numpy as np
from scipy.integrate import trapezoid

from vlmodels.hplc.column import (
    HPLCColumn,
    LangmuirSpecies,
    PartitionSpecies,
    HPLCResult,
    pulse_inlet,
    step_inlet,
    injection_inlet,
    _f_protonated,
)
from PyOMES.numerics.spatial import (
    available_advection_schemes,
    available_dispersion_schemes,
)

# Minimal species list for tests that don't focus on species behaviour
_MINIMAL_SP = [PartitionSpecies("test_sp", K_D=0.5, kf=30)]


# ═══════════════════════════════════════════════════════════════════════
#  _f_protonated helper
# ═══════════════════════════════════════════════════════════════════════

class TestFProtonated:
    def test_at_pKa_is_half(self):
        assert abs(_f_protonated(4.76, 4.76) - 0.5) < 1e-10

    def test_low_pH_approaches_one(self):
        assert _f_protonated(1.0, 4.76) > 0.999

    def test_high_pH_approaches_zero(self):
        assert _f_protonated(10.0, 4.76) < 0.0001

    def test_vectorised(self):
        pH = np.array([2.0, 4.76, 8.0])
        result = _f_protonated(pH, 4.76)
        assert result.shape == (3,)
        assert result[0] > 0.99
        assert abs(result[1] - 0.5) < 1e-10
        assert result[2] < 0.01

    def test_lactic_acid_at_pH_2(self):
        """Lactic acid (pKa 3.86) at pH 2: ~98.6% protonated."""
        f = _f_protonated(2.0, 3.86)
        assert 0.98 < f < 0.995


# ═══════════════════════════════════════════════════════════════════════
#  Species dataclasses
# ═══════════════════════════════════════════════════════════════════════

class TestSpecies:
    def test_langmuir_defaults(self):
        sp = LangmuirSpecies("glucose")
        assert sp.name == "glucose"
        assert sp.q_max > 0
        assert sp.K_L > 0
        assert sp.kf > 0

    def test_partition_defaults(self):
        sp = PartitionSpecies("ethanol")
        assert sp.name == "ethanol"
        assert sp.K_D == 0.5
        assert sp.pKa is None

    def test_partition_with_pKa(self):
        sp = PartitionSpecies("lactic acid", K_D=0.674, pKa=3.86)
        assert sp.pKa == 3.86

    def test_langmuir_custom(self):
        sp = LangmuirSpecies("sucrose", q_max=0.01, K_L=50.0, kf=20.0, MW=342.3)
        assert sp.q_max == 0.01
        assert sp.K_L == 50.0
        assert sp.MW == 342.3


# ═══════════════════════════════════════════════════════════════════════
#  HPLCColumn construction
# ═══════════════════════════════════════════════════════════════════════

class TestColumnConstruction:
    def test_default_construction(self):
        sp = [PartitionSpecies("test", K_D=0.5)]
        col = HPLCColumn(species=sp)
        assert col.n_cells == 100
        assert col.length_cm == 25.0
        assert col.void_fraction == 0.40

    def test_with_species(self):
        sp = [LangmuirSpecies("A"), PartitionSpecies("B")]
        col = HPLCColumn(species=sp)
        assert len(col.species) == 2
        assert col._n_sp == 2

    def test_n_cells_and_target_plates_exclusive(self):
        """Cannot specify both n_cells and target_plates."""
        with pytest.raises(ValueError, match="n_cells or target_plates"):
            HPLCColumn(species=_MINIMAL_SP, n_cells=100, target_plates=5000)

    def test_target_plates_sets_n_cells(self):
        col = HPLCColumn(species=_MINIMAL_SP, target_plates=500, advection_scheme="upwind")
        # Upwind: N_plates ≈ N_cells, but rounding may differ by a few
        assert abs(col.n_cells - 500) < 10

    def test_invalid_advection_scheme_raises(self):
        with pytest.raises(ValueError, match="Unknown advection"):
            HPLCColumn(species=_MINIMAL_SP, advection_scheme="nonexistent")

    def test_invalid_dispersion_scheme_raises(self):
        with pytest.raises(ValueError, match="Unknown dispersion"):
            HPLCColumn(species=_MINIMAL_SP, dispersion_scheme="nonexistent")

    def test_advection_scheme_stored(self):
        col = HPLCColumn(species=_MINIMAL_SP, advection_scheme="tvd_vanleer")
        assert col.advection_scheme == "tvd_vanleer"

    def test_dispersion_scheme_stored(self):
        col = HPLCColumn(species=_MINIMAL_SP, dispersion_scheme="central_4th")
        assert col.dispersion_scheme == "central_4th"

    def test_mixed_species(self):
        """Can mix Langmuir and Partition species in the same column."""
        sp = [LangmuirSpecies("A", q_max=0.01, K_L=10.0),
              PartitionSpecies("B", K_D=0.3)]
        col = HPLCColumn(species=sp, n_cells=20)
        assert col._n_sp == 2

    def test_mobile_phase_pH(self):
        sp = [PartitionSpecies("acid", K_D=0.5, pKa=4.0)]
        col = HPLCColumn(species=sp, mobile_phase_pH=2.1, K_D_reference_pH=2.1)
        assert col.mobile_phase_pH == 2.1
        assert col.K_D_reference_pH == 2.1

    def test_K_D_neutral_back_calculation(self):
        """K_D_neutral should be K_D / f_HA(reference_pH, pKa)."""
        sp = [PartitionSpecies("lactic", K_D=0.674, pKa=3.86)]
        col = HPLCColumn(species=sp, mobile_phase_pH=2.1, K_D_reference_pH=2.1)
        f_HA = _f_protonated(2.1, 3.86)
        expected = 0.674 / f_HA
        assert abs(col._K_D_neutral["lactic"] - expected) < 1e-10


# ═══════════════════════════════════════════════════════════════════════
#  Column properties
# ═══════════════════════════════════════════════════════════════════════

class TestColumnProperties:
    def test_column_volume(self):
        col = HPLCColumn(species=_MINIMAL_SP, length_cm=30.0, diameter_cm=0.78)
        expected = np.pi / 4 * 0.78**2 * 30.0
        assert abs(col.column_volume_mL - expected) < 1e-6

    def test_dead_time(self):
        col = HPLCColumn(species=_MINIMAL_SP, length_cm=30.0, diameter_cm=0.78,
                         void_fraction=0.40, flow_rate_mL_min=1.0)
        V_void = col.column_volume_mL * 0.40
        expected = V_void / 1.0
        assert abs(col.dead_time_min - expected) < 1e-6

    def test_effective_plates_upwind(self):
        """Upwind: N_eff ≈ N_cells (numerical diffusion dominates)."""
        col = HPLCColumn(species=_MINIMAL_SP, n_cells=200, advection_scheme="upwind")
        assert abs(col.effective_plates - 200) < 5

    def test_effective_plates_tvd_higher(self):
        """TVD should give more effective plates than upwind at same n_cells."""
        col_up = HPLCColumn(species=_MINIMAL_SP, n_cells=200, advection_scheme="upwind")
        col_tvd = HPLCColumn(species=_MINIMAL_SP, n_cells=200, advection_scheme="tvd_vanleer")
        assert col_tvd.effective_plates > col_up.effective_plates * 10

    def test_peclet_number(self):
        col = HPLCColumn(species=_MINIMAL_SP, length_cm=25.0, flow_rate_mL_min=1.0,
                         void_fraction=0.40, D_ax_default=0.001)
        assert col.peclet_number > 0

    def test_summary_contains_key_info(self):
        sp = [PartitionSpecies("glucose", K_D=0.24)]
        col = HPLCColumn(species=sp, advection_scheme="tvd_vanleer",
                         mobile_phase_pH=2.1)
        s = col.summary()
        assert "tvd_vanleer" in s
        assert "glucose" in s
        assert "Mobile pH" in s
        assert "2.10" in s


# ═══════════════════════════════════════════════════════════════════════
#  Inlet functions
# ═══════════════════════════════════════════════════════════════════════

class TestPulseInlet:
    def test_nonzero_during_pulse(self):
        inlet = pulse_inlet({"A": 0.01}, t_start_min=1.0, duration_min=0.5)
        result = inlet(1.25)
        assert result["A"] == 0.01

    def test_zero_outside_pulse(self):
        inlet = pulse_inlet({"A": 0.01}, t_start_min=1.0, duration_min=0.5)
        assert inlet(0.5)["A"] == 0.0
        assert inlet(2.0)["A"] == 0.0

    def test_boundary_inclusive(self):
        inlet = pulse_inlet({"A": 0.01}, t_start_min=1.0, duration_min=0.5)
        assert inlet(1.0)["A"] == 0.01   # start
        assert inlet(1.5)["A"] == 0.01   # end


class TestStepInlet:
    def test_zero_before_start(self):
        inlet = step_inlet({"A": 0.01}, t_start_min=1.0)
        assert inlet(0.5)["A"] == 0.0

    def test_nonzero_after_start(self):
        inlet = step_inlet({"A": 0.01}, t_start_min=1.0)
        assert inlet(1.0)["A"] == 0.01
        assert inlet(100.0)["A"] == 0.01


class TestInjectionInlet:
    def test_mobile_phase_before_injection(self):
        inlet = injection_inlet(
            sample_conc={"glucose": 0.025},
            mobile_phase_conc={"sulfate": 0.001},
            t_start_min=1.0, injection_volume_uL=20.0,
            flow_rate_mL_min=0.6, tau_mix_s=0.5,
        )
        result = inlet(0.5)
        assert result["sulfate"] == 0.001
        assert result["glucose"] == 0.0

    def test_mobile_phase_long_after_injection(self):
        inlet = injection_inlet(
            sample_conc={"glucose": 0.025},
            mobile_phase_conc={"sulfate": 0.001},
            t_start_min=1.0, injection_volume_uL=20.0,
            flow_rate_mL_min=0.6, tau_mix_s=0.5,
        )
        result = inlet(100.0)
        assert abs(result["sulfate"] - 0.001) < 1e-10
        assert abs(result["glucose"]) < 1e-10

    def test_sample_during_injection(self):
        inlet = injection_inlet(
            sample_conc={"glucose": 0.025, "sulfate": 0.0},
            mobile_phase_conc={"glucose": 0.0, "sulfate": 0.001},
            t_start_min=0.5, injection_volume_uL=20.0,
            flow_rate_mL_min=0.6, tau_mix_s=0.5,
        )
        # Mid-injection: glucose should be rising, sulfate dropping
        dur = 20e-3 / 0.6  # injection duration in min
        mid = 0.5 + dur / 2
        result = inlet(mid)
        assert result["glucose"] > 0
        assert result["sulfate"] < 0.001

    def test_tau_zero_is_rectangular(self):
        """tau_mix_s=0 should give a perfect rectangular pulse."""
        inlet = injection_inlet(
            sample_conc={"A": 0.01},
            mobile_phase_conc={"A": 0.0},
            t_start_min=1.0, injection_volume_uL=20.0,
            flow_rate_mL_min=0.6, tau_mix_s=0.0,
        )
        dur = 20e-3 / 0.6
        assert inlet(0.9)["A"] == 0.0
        assert inlet(1.0 + dur / 2)["A"] == 0.01
        assert inlet(1.0 + dur + 0.01)["A"] == 0.0

    def test_mass_conservation(self):
        """Integrated mass should match injection_volume × concentration."""
        C_sample = 0.025  # mol/L
        V_inj_uL = 20.0
        Q = 0.6  # mL/min
        inlet = injection_inlet(
            sample_conc={"A": C_sample},
            mobile_phase_conc={"A": 0.0},
            t_start_min=0.5, injection_volume_uL=V_inj_uL,
            flow_rate_mL_min=Q, tau_mix_s=0.5,
        )
        t = np.linspace(0, 5, 10000)
        C = np.array([inlet(ti)["A"] for ti in t])
        mass_injected = trapezoid(C, t) * Q  # mol/L × min × mL/min = mol × 1e-3
        expected = C_sample * V_inj_uL * 1e-3  # mol/L × mL × 1e-3 L/mL
        assert abs(mass_injected - expected) / expected < 0.01  # <1% error


# ═══════════════════════════════════════════════════════════════════════
#  Single-species simulations
# ═══════════════════════════════════════════════════════════════════════

class TestSingleSpeciesPartition:
    """Test a single partition species simulation."""

    @pytest.fixture
    def column_and_result(self):
        sp = [PartitionSpecies("glucose", K_D=0.243, kf=30)]
        col = HPLCColumn(
            length_cm=30.0, diameter_cm=0.78,
            void_fraction=0.30, bulk_density_g_mL=0.40,
            flow_rate_mL_min=0.6, n_cells=50,
            species=sp, intraparticle_porosity=0.50,
            advection_scheme="tvd_vanleer",
        )
        inlet = pulse_inlet({"glucose": 0.025}, t_start_min=0.5, duration_min=0.033)
        result = col.simulate(t_end_min=20, inlet_fn=inlet, n_output=300)
        return col, result

    def test_returns_hplc_result(self, column_and_result):
        _, result = column_and_result
        assert isinstance(result, HPLCResult)

    def test_solver_success(self, column_and_result):
        _, result = column_and_result
        assert result.solver_info["success"]

    def test_peak_detected(self, column_and_result):
        _, result = column_and_result
        C = result.C_outlet["glucose"]
        assert C.max() > 0

    def test_peak_position_reasonable(self, column_and_result):
        col, result = column_and_result
        C = result.C_outlet["glucose"]
        t_R = result.t_min[np.argmax(C)]
        # Expected: t_R ≈ (ε + (1-ε)×ε_p×K_D) × V_col / Q
        V_col = col.column_volume_mL
        eps, eps_p, K_D, Q = 0.30, 0.50, 0.243, 0.6
        t_R_expected = (eps + (1 - eps) * eps_p * K_D) * V_col / Q
        assert abs(t_R - t_R_expected) / t_R_expected < 0.10  # <10% error

    def test_outlet_concentrations_nonnegative(self, column_and_result):
        _, result = column_and_result
        assert np.all(result.C_outlet["glucose"] >= 0)

    def test_spatial_profiles_available(self, column_and_result):
        col, result = column_and_result
        assert result.C_all["glucose"].shape == (300, col.n_cells)
        assert result.q_all["glucose"].shape == (300, col.n_cells)


class TestSingleSpeciesLangmuir:
    """Test a single Langmuir species simulation."""

    def test_peak_elutes(self):
        sp = [LangmuirSpecies("sucrose", q_max=0.005, K_L=24.4, kf=25)]
        col = HPLCColumn(
            length_cm=30.0, diameter_cm=0.78,
            void_fraction=0.40, bulk_density_g_mL=0.60,
            flow_rate_mL_min=0.6, n_cells=50,
            species=sp, advection_scheme="tvd_vanleer",
        )
        inlet = pulse_inlet({"sucrose": 0.01}, t_start_min=0.5, duration_min=0.033)
        result = col.simulate(t_end_min=15, inlet_fn=inlet, n_output=200)
        assert result.C_outlet["sucrose"].max() > 0


# ═══════════════════════════════════════════════════════════════════════
#  pH-dependent retention
# ═══════════════════════════════════════════════════════════════════════

class TestPHDependentRetention:
    """Test that pH correction shifts acid retention times."""

    def _run_at_pH(self, pH_op, pH_ref=2.1):
        sp = [
            PartitionSpecies("glucose", K_D=0.243, kf=30, pKa=None),
            PartitionSpecies("lactic", K_D=0.674, kf=40, pKa=3.86),
        ]
        col = HPLCColumn(
            length_cm=30.0, diameter_cm=0.78,
            void_fraction=0.30, bulk_density_g_mL=0.40,
            flow_rate_mL_min=0.6, n_cells=50,
            species=sp, intraparticle_porosity=0.50,
            advection_scheme="tvd_vanleer",
            mobile_phase_pH=pH_op, K_D_reference_pH=pH_ref,
        )
        inlet = pulse_inlet(
            {"glucose": 0.025, "lactic": 0.020},
            t_start_min=0.5, duration_min=0.033,
        )
        return col.simulate(t_end_min=25, inlet_fn=inlet, n_output=300)

    def test_glucose_unchanged_with_pH(self):
        """Glucose (no pKa) should elute at the same time at any pH."""
        r1 = self._run_at_pH(2.1)
        r2 = self._run_at_pH(5.0)
        t1 = r1.t_min[np.argmax(r1.C_outlet["glucose"])]
        t2 = r2.t_min[np.argmax(r2.C_outlet["glucose"])]
        assert abs(t1 - t2) < 0.2  # <0.2 min difference

    def test_acid_elutes_earlier_at_higher_pH(self):
        """Lactic acid should elute earlier at pH 5 than pH 2.1."""
        r1 = self._run_at_pH(2.1)
        r2 = self._run_at_pH(5.0)
        t1 = r1.t_min[np.argmax(r1.C_outlet["lactic"])]
        t2 = r2.t_min[np.argmax(r2.C_outlet["lactic"])]
        assert t2 < t1 - 1.0  # at least 1 minute earlier

    def test_reference_pH_reproduces_original(self):
        """At pH == reference_pH, retention should match no-pH case."""
        sp = [PartitionSpecies("lactic", K_D=0.674, kf=40, pKa=3.86)]
        col_no_pH = HPLCColumn(
            length_cm=30.0, diameter_cm=0.78,
            void_fraction=0.30, bulk_density_g_mL=0.40,
            flow_rate_mL_min=0.6, n_cells=50,
            species=sp, intraparticle_porosity=0.50,
            advection_scheme="tvd_vanleer",
        )
        col_pH = HPLCColumn(
            length_cm=30.0, diameter_cm=0.78,
            void_fraction=0.30, bulk_density_g_mL=0.40,
            flow_rate_mL_min=0.6, n_cells=50,
            species=sp, intraparticle_porosity=0.50,
            advection_scheme="tvd_vanleer",
            mobile_phase_pH=2.1, K_D_reference_pH=2.1,
        )
        inlet = pulse_inlet({"lactic": 0.020}, t_start_min=0.5, duration_min=0.033)
        r1 = col_no_pH.simulate(t_end_min=25, inlet_fn=inlet, n_output=300)
        r2 = col_pH.simulate(t_end_min=25, inlet_fn=inlet, n_output=300)
        t1 = r1.t_min[np.argmax(r1.C_outlet["lactic"])]
        t2 = r2.t_min[np.argmax(r2.C_outlet["lactic"])]
        assert abs(t1 - t2) < 0.2


# ═══════════════════════════════════════════════════════════════════════
#  Summary and repr
# ═══════════════════════════════════════════════════════════════════════

class TestSummary:
    def test_summary_species_listed(self):
        sp = [PartitionSpecies("ethanol", K_D=1.77),
              LangmuirSpecies("glucose", q_max=0.005, K_L=24.4)]
        col = HPLCColumn(species=sp, n_cells=50)
        s = col.summary()
        assert "ethanol" in s
        assert "glucose" in s
        assert "Langmuir" in s
        assert "Partition" in s

    def test_summary_pH_info(self):
        sp = [PartitionSpecies("acid", K_D=0.5, pKa=4.0)]
        col = HPLCColumn(species=sp, mobile_phase_pH=3.0, K_D_reference_pH=2.1)
        s = col.summary()
        assert "Mobile pH" in s
        assert "K_D_neutral" in s


# ═══════════════════════════════════════════════════════════════════════
#  Input validation (Phase A)
# ═══════════════════════════════════════════════════════════════════════

class TestInputValidation:
    """Tests that invalid physical parameters are rejected at construction."""

    def test_zero_flow_rate_raises(self):
        with pytest.raises(ValueError, match="flow_rate_mL_min.*positive"):
            HPLCColumn(species=_MINIMAL_SP, flow_rate_mL_min=0)

    def test_negative_flow_rate_raises(self):
        with pytest.raises(ValueError, match="flow_rate_mL_min.*positive"):
            HPLCColumn(species=_MINIMAL_SP, flow_rate_mL_min=-1)

    def test_zero_void_fraction_raises(self):
        with pytest.raises(ValueError, match="void_fraction"):
            HPLCColumn(species=_MINIMAL_SP, void_fraction=0)

    def test_void_fraction_one_raises(self):
        with pytest.raises(ValueError, match="void_fraction"):
            HPLCColumn(species=_MINIMAL_SP, void_fraction=1.0)

    def test_negative_void_fraction_raises(self):
        with pytest.raises(ValueError, match="void_fraction"):
            HPLCColumn(species=_MINIMAL_SP, void_fraction=-0.1)

    def test_zero_length_raises(self):
        with pytest.raises(ValueError, match="length_cm.*positive"):
            HPLCColumn(species=_MINIMAL_SP, length_cm=0)

    def test_zero_diameter_raises(self):
        with pytest.raises(ValueError, match="diameter_cm.*positive"):
            HPLCColumn(species=_MINIMAL_SP, diameter_cm=0)

    def test_zero_temperature_raises(self):
        with pytest.raises(ValueError, match="T_K.*positive"):
            HPLCColumn(species=_MINIMAL_SP, T_K=0)

    def test_empty_species_raises(self):
        with pytest.raises(ValueError, match="species"):
            HPLCColumn(species=[])

    def test_negative_K_D_raises(self):
        with pytest.raises(ValueError, match="K_D.*non-negative"):
            HPLCColumn(species=[PartitionSpecies("x", K_D=-0.5)])

    def test_negative_q_max_raises(self):
        with pytest.raises(ValueError, match="q_max.*positive"):
            HPLCColumn(species=[LangmuirSpecies("x", q_max=-1)])

    def test_negative_K_L_raises(self):
        with pytest.raises(ValueError, match="K_L.*non-negative"):
            HPLCColumn(species=[LangmuirSpecies("x", q_max=0.01, K_L=-1)])

    def test_negative_kf_raises(self):
        with pytest.raises(ValueError, match="kf.*non-negative"):
            HPLCColumn(species=[PartitionSpecies("x", K_D=0.5, kf=-1)])

    def test_duplicate_species_names_raises(self):
        sp = [PartitionSpecies("A", K_D=0.2), PartitionSpecies("A", K_D=0.8)]
        with pytest.raises(ValueError, match="Duplicate species"):
            HPLCColumn(species=sp)

    def test_negative_D_ax_raises(self):
        with pytest.raises(ValueError, match="D_ax_default.*non-negative"):
            HPLCColumn(species=_MINIMAL_SP, D_ax_default=-0.01)

    def test_negative_bulk_density_raises(self):
        with pytest.raises(ValueError, match="bulk_density_g_mL.*non-negative"):
            HPLCColumn(species=_MINIMAL_SP, bulk_density_g_mL=-0.1)


class TestInjectionInletValidation:
    """Tests that invalid injection parameters are rejected."""

    def test_zero_flow_rate_raises(self):
        with pytest.raises(ValueError, match="flow_rate_mL_min.*positive"):
            injection_inlet(sample_conc={"A": 0.01}, flow_rate_mL_min=0)

    def test_negative_tau_raises(self):
        with pytest.raises(ValueError, match="tau_mix_s.*non-negative"):
            injection_inlet(sample_conc={"A": 0.01}, tau_mix_s=-1)

    def test_zero_volume_raises(self):
        with pytest.raises(ValueError, match="injection_volume_uL.*positive"):
            injection_inlet(sample_conc={"A": 0.01}, injection_volume_uL=0)

    def test_negative_volume_raises(self):
        with pytest.raises(ValueError, match="injection_volume_uL.*positive"):
            injection_inlet(sample_conc={"A": 0.01}, injection_volume_uL=-10)


# ═══════════════════════════════════════════════════════════════════════
#  Ambiguous configuration warnings (Phase C)
# ═══════════════════════════════════════════════════════════════════════

class TestConfigurationWarnings:
    """Tests that ambiguous-but-valid configurations emit warnings."""

    def test_mobile_phase_pH_without_reference_warns(self):
        sp = [PartitionSpecies("acid", K_D=0.5, pKa=3.86)]
        with pytest.warns(UserWarning, match="K_D_reference_pH"):
            HPLCColumn(species=sp, mobile_phase_pH=5.0)

    def test_mobile_phase_pH_without_reference_no_pKa_no_warn(self):
        """No warning if no species have pKa (pH correction is irrelevant)."""
        sp = [PartitionSpecies("glucose", K_D=0.5, pKa=None)]
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            HPLCColumn(species=sp, mobile_phase_pH=5.0)

    def test_mobile_phase_pH_with_reference_no_warn(self):
        """No warning when K_D_reference_pH is explicitly set."""
        sp = [PartitionSpecies("acid", K_D=0.5, pKa=3.86)]
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            HPLCColumn(species=sp, mobile_phase_pH=5.0, K_D_reference_pH=2.1)

    def test_high_void_fraction_warns(self):
        sp = [PartitionSpecies("test", K_D=0.5)]
        with pytest.warns(UserWarning, match="void_fraction.*unusually high"):
            HPLCColumn(species=sp, void_fraction=0.96)

    def test_normal_void_fraction_no_warn(self):
        sp = [PartitionSpecies("test", K_D=0.5)]
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            HPLCColumn(species=sp, void_fraction=0.40)

    def test_C_init_unknown_species_warns(self):
        sp = [PartitionSpecies("glucose", K_D=0.2, kf=30)]
        col = HPLCColumn(species=sp, n_cells=20)
        inlet = pulse_inlet({"glucose": 0.01}, t_start_min=0.5, duration_min=0.1)
        with pytest.warns(UserWarning, match="C_init.*not in the column"):
            col.simulate(t_end_min=5, inlet_fn=inlet, n_output=20,
                        C_init={"fructose": 0.01})

    def test_q_init_unknown_species_warns(self):
        sp = [PartitionSpecies("glucose", K_D=0.2, kf=30)]
        col = HPLCColumn(species=sp, n_cells=20)
        inlet = pulse_inlet({"glucose": 0.01}, t_start_min=0.5, duration_min=0.1)
        with pytest.warns(UserWarning, match="q_init.*not in the column"):
            col.simulate(t_end_min=5, inlet_fn=inlet, n_output=20,
                        q_init={"fructose": 0.001})


# ═══════════════════════════════════════════════════════════════════════
#  Post-construction immutability (G.1)
# ═══════════════════════════════════════════════════════════════════════

class TestColumnFrozen:
    """Tests that the column and its species are immutable after construction."""

    def test_species_dataclass_frozen(self):
        sp = PartitionSpecies("glucose", K_D=0.243)
        with pytest.raises(AttributeError):
            sp.K_D = 0.5

    def test_langmuir_dataclass_frozen(self):
        sp = LangmuirSpecies("sucrose", q_max=0.01)
        with pytest.raises(AttributeError):
            sp.q_max = 0.05

    def test_species_is_tuple(self):
        sp = [PartitionSpecies("A", K_D=0.2)]
        col = HPLCColumn(species=sp, n_cells=20)
        assert isinstance(col.species, tuple)

    def test_frozen_flow_rate(self):
        col = HPLCColumn(species=_MINIMAL_SP)
        with pytest.raises(AttributeError, match="Cannot modify"):
            col.flow_rate_mL_min = 1.0

    def test_frozen_void_fraction(self):
        col = HPLCColumn(species=_MINIMAL_SP)
        with pytest.raises(AttributeError, match="Cannot modify"):
            col.void_fraction = 0.5

    def test_frozen_length(self):
        col = HPLCColumn(species=_MINIMAL_SP)
        with pytest.raises(AttributeError, match="Cannot modify"):
            col.length_cm = 50.0

    def test_frozen_species(self):
        col = HPLCColumn(species=_MINIMAL_SP)
        with pytest.raises(AttributeError, match="Cannot modify"):
            col.species = (PartitionSpecies("new", K_D=0.1),)

    def test_non_frozen_label(self):
        col = HPLCColumn(species=_MINIMAL_SP)
        col.label = "renamed"
        assert col.label == "renamed"


# ═══════════════════════════════════════════════════════════════════════
#  Recycle mode
# ═══════════════════════════════════════════════════════════════════════

class TestRecycleMode:
    """Tests for closed-loop recycle chromatography."""

    @pytest.fixture
    def col_and_inlet(self):
        sp = [
            PartitionSpecies("glucose", K_D=0.243, kf=30, pKa=None),
            PartitionSpecies("lactic", K_D=0.674, kf=40, pKa=3.86),
        ]
        col = HPLCColumn(
            length_cm=30.0, diameter_cm=0.78,
            void_fraction=0.30, bulk_density_g_mL=0.40,
            flow_rate_mL_min=0.6, n_cells=50,
            species=sp, intraparticle_porosity=0.50,
            advection_scheme="tvd_vanleer",
            mobile_phase_pH=2.1, K_D_reference_pH=2.1,
        )
        inlet = pulse_inlet(
            {"glucose": 0.025, "lactic": 0.020},
            t_start_min=0.5, duration_min=0.033,
        )
        return col, inlet

    def test_recycle_produces_repeated_peaks(self, col_and_inlet):
        """Recycle at dead time should produce multiple peaks per species."""
        col, inlet = col_and_inlet
        t_dead = col.dead_time_min
        r = col.simulate(t_end_min=60, inlet_fn=inlet, n_output=800,
                         recycle_after_min=t_dead)
        assert r.solver_info["success"]
        # Glucose should have multiple peaks (recycles every ~9 min)
        C = r.C_outlet["glucose"]
        threshold = C.max() * 0.05
        n_peaks = sum(1 for i in range(1, len(C)-1)
                      if C[i] > C[i-1] and C[i] > C[i+1] and C[i] > threshold)
        assert n_peaks >= 3

    def test_recycle_peak_spacing_matches_retention(self, col_and_inlet):
        """Peak spacing in recycle mode should match predicted t_R."""
        col, inlet = col_and_inlet
        t_dead = col.dead_time_min
        r = col.simulate(t_end_min=60, inlet_fn=inlet, n_output=800,
                         recycle_after_min=t_dead)
        C = r.C_outlet["glucose"]
        threshold = C.max() * 0.05
        peak_times = [r.t_min[i] for i in range(1, len(C)-1)
                      if C[i] > C[i-1] and C[i] > C[i+1] and C[i] > threshold]
        assert len(peak_times) >= 3
        spacings = [peak_times[i+1] - peak_times[i]
                    for i in range(len(peak_times)-1)]
        t_R = col.predicted_retention_min("glucose")
        mean_spacing = np.mean(spacings)
        assert abs(mean_spacing - t_R) / t_R < 0.05  # <5% error

    def test_no_recycle_gives_single_peaks(self, col_and_inlet):
        """Without recycle, each species should have exactly one peak."""
        col, inlet = col_and_inlet
        r = col.simulate(t_end_min=30, inlet_fn=inlet, n_output=400)
        C = r.C_outlet["glucose"]
        threshold = C.max() * 0.05
        n_peaks = sum(1 for i in range(1, len(C)-1)
                      if C[i] > C[i-1] and C[i] > C[i+1] and C[i] > threshold)
        assert n_peaks == 1

    def test_recycle_none_is_single_pass(self, col_and_inlet):
        """recycle_after_min=None should behave identically to no recycle."""
        col, inlet = col_and_inlet
        r1 = col.simulate(t_end_min=25, inlet_fn=inlet, n_output=300,
                          recycle_after_min=None)
        r2 = col.simulate(t_end_min=25, inlet_fn=inlet, n_output=300)
        np.testing.assert_allclose(
            r1.C_outlet["glucose"], r2.C_outlet["glucose"], atol=1e-10)

    def test_recycle_after_last_elution(self, col_and_inlet):
        """Switch after last elution: only slow species recycles."""
        col, inlet = col_and_inlet
        t_last = col.predicted_last_elution_min()
        r = col.simulate(t_end_min=60, inlet_fn=inlet, n_output=800,
                         recycle_after_min=t_last)
        # Glucose elutes before switch → single peak
        C_glu = r.C_outlet["glucose"]
        thresh_glu = C_glu.max() * 0.05
        n_glu = sum(1 for i in range(1, len(C_glu)-1)
                    if C_glu[i] > C_glu[i-1] and C_glu[i] > C_glu[i+1]
                    and C_glu[i] > thresh_glu)
        assert n_glu == 1
        # Lactic acid is slower — may still be partially in column
        # (at least should not crash)
        assert r.solver_info["success"]


class TestPredictedRetention:
    """Tests for predicted_retention_min and predicted_last_elution_min."""

    def test_glucose_prediction(self):
        sp = [PartitionSpecies("glucose", K_D=0.243, kf=30)]
        col = HPLCColumn(
            length_cm=30.0, diameter_cm=0.78,
            void_fraction=0.30, bulk_density_g_mL=0.40,
            flow_rate_mL_min=0.6, n_cells=50,
            species=sp, intraparticle_porosity=0.50,
            advection_scheme="tvd_vanleer",
        )
        t_R = col.predicted_retention_min("glucose")
        assert 8.0 < t_R < 11.0  # ~9.2 min expected

    def test_unknown_species_raises(self):
        sp = [PartitionSpecies("glucose", K_D=0.243, kf=30)]
        col = HPLCColumn(species=sp, n_cells=20)
        with pytest.raises(ValueError, match="not in column"):
            col.predicted_retention_min("fructose")

    def test_last_elution_is_max(self):
        sp = [
            PartitionSpecies("glucose", K_D=0.243, kf=30),
            PartitionSpecies("ethanol", K_D=1.774, kf=25),
        ]
        col = HPLCColumn(
            length_cm=30.0, diameter_cm=0.78,
            void_fraction=0.30, bulk_density_g_mL=0.40,
            flow_rate_mL_min=0.6, n_cells=50,
            species=sp, intraparticle_porosity=0.50,
        )
        assert col.predicted_last_elution_min() == col.predicted_retention_min("ethanol")

    def test_langmuir_prediction(self):
        sp = [LangmuirSpecies("sucrose", q_max=0.005, K_L=24.4, kf=25)]
        col = HPLCColumn(
            length_cm=30.0, diameter_cm=0.78,
            void_fraction=0.40, bulk_density_g_mL=0.60,
            flow_rate_mL_min=0.6, n_cells=50,
            species=sp,
        )
        t_R = col.predicted_retention_min("sucrose")
        assert t_R > col.dead_time_min  # retained beyond void


# ═══════════════════════════════════════════════════════════════════════
#  Peak shaving / divert functions
# ═══════════════════════════════════════════════════════════════════════

class TestDivertFunctions:
    """Tests for peak-shaving divert functions in recycle mode."""

    @pytest.fixture
    def col_3sp(self):
        sp = [
            PartitionSpecies("sucrose", K_D=0.040, kf=25, pKa=None),
            PartitionSpecies("glucose", K_D=0.243, kf=30, pKa=None),
            PartitionSpecies("fructose", K_D=0.339, kf=30, pKa=None),
        ]
        col = HPLCColumn(
            length_cm=30.0, diameter_cm=0.78,
            void_fraction=0.30, bulk_density_g_mL=0.40,
            flow_rate_mL_min=0.6, n_cells=50,
            species=sp, intraparticle_porosity=0.50,
            advection_scheme="tvd_vanleer",
            mobile_phase_pH=2.1, K_D_reference_pH=2.1,
        )
        return col

    def test_divert_windows_zeros_during_window(self):
        from vlmodels.hplc.column import make_divert_windows
        divert = make_divert_windows([(5.0, 10.0), (20.0, 25.0)])
        C = np.array([0.1, 0.2, 0.3])
        # Inside window → zeros
        result = divert(7.0, C)
        np.testing.assert_array_equal(result, np.zeros(3))
        # Outside window → pass through
        result = divert(15.0, C)
        np.testing.assert_array_equal(result, C)

    def test_divert_threshold(self):
        from vlmodels.hplc.column import make_divert_threshold
        divert = make_divert_threshold(0.01)
        C_high = np.array([0.005, 0.006, 0.001])  # sum = 0.012 > 0.01
        C_low = np.array([0.003, 0.002, 0.001])   # sum = 0.006 < 0.01
        np.testing.assert_array_equal(divert(1.0, C_high), np.zeros(3))
        np.testing.assert_array_equal(divert(1.0, C_low), C_low)

    def test_peak_shave_schedule_generates_windows(self, col_3sp):
        from vlmodels.hplc.column import make_peak_shave_schedule
        divert = make_peak_shave_schedule(col_3sp, ["sucrose"], n_passes=5)
        # At sucrose peak time (pass 1) → should divert
        t_R = col_3sp.predicted_retention_min("sucrose")
        C = np.array([0.1, 0.2, 0.3])
        result = divert(t_R, C)
        np.testing.assert_array_equal(result, np.zeros(3))
        # Far from any sucrose peak → should pass through
        result = divert(t_R * 1.5, C)
        np.testing.assert_array_equal(result, C)

    def test_peak_shave_removes_sucrose(self, col_3sp):
        """With peak shaving, sucrose signal should diminish across passes."""
        from vlmodels.hplc.column import make_peak_shave_schedule, pulse_inlet
        col = col_3sp
        inlet = pulse_inlet(
            {"sucrose": 0.015, "glucose": 0.025, "fructose": 0.025},
            t_start_min=0.5, duration_min=0.033,
        )
        divert = make_peak_shave_schedule(col, ["sucrose"], n_passes=10)
        r = col.simulate(t_end_min=50, inlet_fn=inlet, n_output=800,
                         recycle_after_min=col.dead_time_min,
                         divert_fn=divert)
        assert r.solver_info["success"]
        # Sucrose should be mostly gone in the second half
        C_suc = r.C_outlet["sucrose"]
        midpoint = len(C_suc) // 2
        first_half_max = C_suc[:midpoint].max()
        second_half_max = C_suc[midpoint:].max()
        assert second_half_max < first_half_max * 0.7

    def test_divert_none_is_full_recycle(self, col_3sp):
        """divert_fn=None should give identical results to no divert."""
        col = col_3sp
        inlet = pulse_inlet(
            {"sucrose": 0.015, "glucose": 0.025, "fructose": 0.025},
            t_start_min=0.5, duration_min=0.033,
        )
        r1 = col.simulate(t_end_min=30, inlet_fn=inlet, n_output=300,
                          recycle_after_min=col.dead_time_min, divert_fn=None)
        r2 = col.simulate(t_end_min=30, inlet_fn=inlet, n_output=300,
                          recycle_after_min=col.dead_time_min)
        np.testing.assert_allclose(
            r1.C_outlet["glucose"], r2.C_outlet["glucose"], atol=1e-10)
