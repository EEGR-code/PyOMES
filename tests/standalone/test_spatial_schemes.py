"""Unit tests for the shared numerics spatial schemes module.

Tests the advection and dispersion scheme registry, individual scheme
behaviour, mass conservation, and numerical diffusion properties.
"""

import pytest
import numpy as np

from PyOMES.numerics.spatial import (
    available_advection_schemes,
    available_dispersion_schemes,
    get_advection_scheme,
    get_dispersion_scheme,
    register_advection_scheme,
    register_dispersion_scheme,
)


# ═══════════════════════════════════════════════════════════════════════
#  Registry tests
# ═══════════════════════════════════════════════════════════════════════

class TestRegistry:
    def test_available_advection_schemes(self):
        schemes = available_advection_schemes()
        assert "upwind" in schemes
        assert "tvd_vanleer" in schemes

    def test_available_dispersion_schemes(self):
        schemes = available_dispersion_schemes()
        assert "central_2nd" in schemes
        assert "central_4th" in schemes

    def test_get_advection_returns_callable(self):
        fn = get_advection_scheme("upwind", u=1.0, dz=0.1, N=10, n_sp=1)
        assert callable(fn)

    def test_get_dispersion_returns_callable(self):
        D_ax = np.array([0.01])
        fn = get_dispersion_scheme("central_2nd", D_ax=D_ax, dz=0.1, N=10, n_sp=1)
        assert callable(fn)

    def test_unknown_advection_raises(self):
        with pytest.raises(ValueError, match="Unknown advection"):
            get_advection_scheme("nonexistent", u=1.0, dz=0.1, N=10, n_sp=1)

    def test_unknown_dispersion_raises(self):
        with pytest.raises(ValueError, match="Unknown dispersion"):
            get_dispersion_scheme("nonexistent", D_ax=np.array([0.01]),
                                 dz=0.1, N=10, n_sp=1)

    def test_case_insensitive(self):
        fn1 = get_advection_scheme("Upwind", u=1.0, dz=0.1, N=10, n_sp=1)
        fn2 = get_advection_scheme("UPWIND", u=1.0, dz=0.1, N=10, n_sp=1)
        assert callable(fn1) and callable(fn2)


# ═══════════════════════════════════════════════════════════════════════
#  Advection scheme tests
# ═══════════════════════════════════════════════════════════════════════

class TestUpwindAdvection:
    """Tests for the first-order upwind advection scheme."""

    def test_returns_correct_shape(self):
        N, n_sp = 20, 2
        advect = get_advection_scheme("upwind", u=1.0, dz=0.1, N=N, n_sp=n_sp)
        C = np.zeros((N, n_sp))
        C_in = np.zeros((1, n_sp))
        result = advect(C, C_in)
        assert result.shape == (N, n_sp)

    def test_zero_concentration_gives_zero_flux(self):
        N, n_sp = 10, 1
        advect = get_advection_scheme("upwind", u=1.0, dz=0.1, N=N, n_sp=n_sp)
        C = np.zeros((N, n_sp))
        C_in = np.zeros((1, n_sp))
        result = advect(C, C_in)
        assert np.allclose(result, 0.0)

    def test_uniform_concentration_gives_zero_flux(self):
        """Uniform field should produce zero advective derivative."""
        N, n_sp = 20, 1
        advect = get_advection_scheme("upwind", u=1.0, dz=0.1, N=N, n_sp=n_sp)
        C = np.ones((N, n_sp)) * 0.5
        C_in = np.ones((1, n_sp)) * 0.5  # inlet matches
        result = advect(C, C_in)
        assert np.allclose(result, 0.0, atol=1e-14)

    def test_inlet_pulse_enters_first_cell(self):
        """A nonzero inlet with empty column should produce positive dC/dt in cell 0."""
        N, n_sp = 10, 1
        advect = get_advection_scheme("upwind", u=1.0, dz=0.1, N=N, n_sp=n_sp)
        C = np.zeros((N, n_sp))
        C_in = np.ones((1, n_sp))
        result = advect(C, C_in)
        assert result[0, 0] > 0  # material entering

    def test_mass_conservation_step(self):
        """Sum of dC/dt × dz should equal net influx - outflux."""
        N, n_sp = 50, 1
        u, dz = 2.0, 0.1
        advect = get_advection_scheme("upwind", u=u, dz=dz, N=N, n_sp=n_sp)
        # Gaussian profile in column
        z = np.linspace(0, (N - 1) * dz, N)
        C = np.exp(-((z - 2.0) ** 2) / 0.1).reshape(N, 1)
        C_in = np.zeros((1, n_sp))
        dCdt = advect(C, C_in)
        # Net change in total mass = inlet flux - outlet flux
        # = u * (C_in - C[-1]) / dz ... integrated
        total_dCdt = np.sum(dCdt[:, 0]) * dz
        inlet_flux = u * C_in[0, 0]
        outlet_flux = u * C[-1, 0]
        expected = inlet_flux - outlet_flux
        assert abs(total_dCdt - expected) < 1e-10


class TestTVDAdvection:
    """Tests for the TVD van Leer advection scheme."""

    def test_returns_correct_shape(self):
        N, n_sp = 20, 2
        advect = get_advection_scheme("tvd_vanleer", u=1.0, dz=0.1, N=N, n_sp=n_sp)
        C = np.zeros((N, n_sp))
        C_in = np.zeros((1, n_sp))
        result = advect(C, C_in)
        assert result.shape == (N, n_sp)

    def test_uniform_concentration_gives_zero_flux(self):
        N, n_sp = 20, 1
        advect = get_advection_scheme("tvd_vanleer", u=1.0, dz=0.1, N=N, n_sp=n_sp)
        C = np.ones((N, n_sp)) * 0.5
        C_in = np.ones((1, n_sp)) * 0.5
        result = advect(C, C_in)
        assert np.allclose(result, 0.0, atol=1e-14)

    def test_mass_conservation_step(self):
        N, n_sp = 50, 1
        u, dz = 2.0, 0.1
        advect = get_advection_scheme("tvd_vanleer", u=u, dz=dz, N=N, n_sp=n_sp)
        z = np.linspace(0, (N - 1) * dz, N)
        C = np.exp(-((z - 2.0) ** 2) / 0.1).reshape(N, 1)
        C_in = np.zeros((1, n_sp))
        dCdt = advect(C, C_in)
        total_dCdt = np.sum(dCdt[:, 0]) * dz
        inlet_flux = u * C_in[0, 0]
        outlet_flux = u * C[-1, 0]
        expected = inlet_flux - outlet_flux
        assert abs(total_dCdt - expected) < 1e-10

    def test_less_diffusive_than_upwind(self):
        """TVD should produce a narrower peak than upwind after transport."""
        N, n_sp = 200, 1
        u, dz, dt = 1.0, 0.1, 0.02  # CFL = u×dt/dz = 0.2, stable
        z = np.linspace(0, (N - 1) * dz, N)

        # Initial narrow Gaussian
        C0 = np.exp(-((z - 5.0) ** 2) / 0.05).reshape(N, 1)
        initial_std = np.sqrt(np.average((z - 5.0)**2, weights=C0[:, 0]))

        widths = {}
        for scheme_name in ("upwind", "tvd_vanleer"):
            advect = get_advection_scheme(scheme_name, u=u, dz=dz, N=N, n_sp=n_sp)
            C = C0.copy()
            C_in = np.zeros((1, n_sp))
            for _ in range(200):
                C = C + dt * advect(C, C_in)
                C = np.maximum(C, 0)

            # Measure peak width (std dev of the mass distribution)
            total_mass = np.sum(C[:, 0])
            if total_mass > 1e-10:
                mean_z = np.sum(z * C[:, 0]) / total_mass
                var_z = np.sum((z - mean_z)**2 * C[:, 0]) / total_mass
                widths[scheme_name] = np.sqrt(var_z)
            else:
                widths[scheme_name] = float("inf")

        # TVD should produce a narrower (less diffused) peak
        assert widths["tvd_vanleer"] < widths["upwind"]


# ═══════════════════════════════════════════════════════════════════════
#  Dispersion scheme tests
# ═══════════════════════════════════════════════════════════════════════

class TestDispersion:
    """Tests for central difference dispersion schemes."""

    def test_central_2nd_shape(self):
        N, n_sp = 20, 2
        D_ax = np.array([0.01, 0.02])
        disp = get_dispersion_scheme("central_2nd", D_ax=D_ax, dz=0.1, N=N, n_sp=n_sp)
        C = np.zeros((N, n_sp))
        C_in = np.zeros((1, n_sp))
        result = disp(C, C_in)
        assert result.shape == (N, n_sp)

    def test_central_4th_shape(self):
        N, n_sp = 20, 2
        D_ax = np.array([0.01, 0.02])
        disp = get_dispersion_scheme("central_4th", D_ax=D_ax, dz=0.1, N=N, n_sp=n_sp)
        C = np.zeros((N, n_sp))
        C_in = np.zeros((1, n_sp))
        result = disp(C, C_in)
        assert result.shape == (N, n_sp)

    def test_uniform_gives_zero(self):
        """Uniform concentration → zero Laplacian."""
        for scheme in ("central_2nd", "central_4th"):
            N, n_sp = 30, 1
            D_ax = np.array([0.05])
            disp = get_dispersion_scheme(scheme, D_ax=D_ax, dz=0.1, N=N, n_sp=n_sp)
            C = np.ones((N, n_sp)) * 1.5
            C_in = np.ones((1, n_sp)) * 1.5  # match uniform field
            result = disp(C, C_in)
            assert np.allclose(result, 0.0, atol=1e-12), f"{scheme} failed"

    def test_delta_spreads_symmetrically(self):
        """A delta function should produce symmetric spreading."""
        N, n_sp = 51, 1
        D_ax = np.array([0.1])
        dz = 0.1
        disp = get_dispersion_scheme("central_2nd", D_ax=D_ax, dz=dz, N=N, n_sp=n_sp)
        C = np.zeros((N, n_sp))
        C[25, 0] = 1.0  # delta at centre
        C_in = np.zeros((1, n_sp))
        dCdt = disp(C, C_in)
        # Centre should decrease (diffusing away)
        assert dCdt[25, 0] < 0
        # Neighbours should increase (receiving material)
        assert dCdt[24, 0] > 0
        assert dCdt[26, 0] > 0
        # Symmetric
        assert abs(dCdt[24, 0] - dCdt[26, 0]) < 1e-14

    def test_mass_conservation(self):
        """Dispersion with zero-flux BCs should conserve total mass."""
        for scheme in ("central_2nd", "central_4th"):
            N, n_sp = 50, 1
            D_ax = np.array([0.05])
            dz = 0.1
            disp = get_dispersion_scheme(scheme, D_ax=D_ax, dz=dz, N=N, n_sp=n_sp)
            # Gaussian profile centred away from boundaries
            z = np.linspace(0, (N - 1) * dz, N)
            C = np.exp(-((z - 2.5) ** 2) / 0.2).reshape(N, 1)
            C_in = np.zeros((1, n_sp))
            dCdt = disp(C, C_in)
            # Total mass change should be ~zero (zero-flux BCs)
            total_change = np.sum(dCdt[:, 0]) * dz
            assert abs(total_change) < 1e-10, f"{scheme}: mass change {total_change}"
