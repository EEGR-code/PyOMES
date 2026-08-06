"""Tests for the opt-in multi-component component-basis pathway."""

import pytest

from PyOMES.chemical_equilibrium.multicomponent import (
    ComponentBasis, ComponentEquilibriumReaction, ComponentSpecies,
    MultiComponentTableau,
    ComponentReactionNetwork,
    component_totals,
    IdealMassActionSolver, MultiComponentProblem,
    mass_action_residuals,
)


def test_complex_keeps_each_independent_component() -> None:
    basis = ComponentBasis(("Fe2", "phosphate", "citrate", "Ca"))
    fe_citrate = ComponentSpecies("FeCit-", -1, {"Fe2": 1, "citrate": 1})
    ca_phosphate = ComponentSpecies("CaHPO4", 0, {"Ca": 1, "phosphate": 1})

    assert basis.vector(fe_citrate) == (1.0, 0.0, 1.0, 0.0)
    assert basis.vector(ca_phosphate) == (0.0, 1.0, 0.0, 1.0)


def test_unknown_component_is_not_silently_discarded() -> None:
    basis = ComponentBasis(("Fe2", "phosphate"))
    unknown = ComponentSpecies("FeCit-", -1, {"Fe2": 1, "citrate": 1})

    with pytest.raises(ValueError, match="citrate"):
        basis.vector(unknown)


def test_tableau_accepts_balanced_fe_citrate_complexation() -> None:
    basis = ComponentBasis(("Fe2", "citrate"))
    fe = ComponentSpecies("Fe2+", +2, {"Fe2": 1})
    citrate = ComponentSpecies("Cit---", -3, {"citrate": 1})
    complexed = ComponentSpecies("FeCit-", -1, {"Fe2": 1, "citrate": 1})
    reaction = ComponentEquilibriumReaction(
        "Fe2_citrate", {fe: -1, citrate: -1, complexed: +1}, log_k=6.10,
    )

    tableau = MultiComponentTableau.build(basis, (fe, citrate, complexed), (reaction,))
    assert tableau.component_matrix == ((1.0, 0.0), (0.0, 1.0), (1.0, 1.0))


def test_tableau_rejects_component_unbalanced_reaction() -> None:
    basis = ComponentBasis(("Fe2", "citrate"))
    fe = ComponentSpecies("Fe2+", +2, {"Fe2": 1})
    citrate = ComponentSpecies("Cit---", -3, {"citrate": 1})
    invalid = ComponentEquilibriumReaction("bad", {fe: -1, citrate: +1}, log_k=1.0)

    with pytest.raises(ValueError, match="not component-balanced"):
        MultiComponentTableau.build(basis, (fe, citrate), (invalid,))


def test_network_rejects_duplicate_equilibrium_constraint() -> None:
    fe = ComponentSpecies("Fe2+", +2, {"Fe2": 1})
    citrate = ComponentSpecies("Cit---", -3, {"citrate": 1})
    complexed = ComponentSpecies("FeCit-", -1, {"Fe2": 1, "citrate": 1})
    reaction = ComponentEquilibriumReaction("Fe2_citrate", {fe: -1, citrate: -1, complexed: 1}, 6.1)

    with pytest.raises(ValueError, match="linearly dependent"):
        ComponentReactionNetwork.build((fe, citrate, complexed), (reaction, reaction))


def test_inventory_sums_free_and_complexed_material() -> None:
    basis = ComponentBasis(("Fe2", "citrate"))
    fe = ComponentSpecies("Fe2+", +2, {"Fe2": 1})
    complexed = ComponentSpecies("FeCit-", -1, {"Fe2": 1, "citrate": 1})

    assert component_totals(basis, {fe: 0.002, complexed: 0.0003}) == {
        "Fe2": pytest.approx(0.0023), "citrate": pytest.approx(0.0003),
    }


def test_mass_action_residual_is_zero_at_equilibrium() -> None:
    fe = ComponentSpecies("Fe2+", +2, {"Fe2": 1})
    citrate = ComponentSpecies("Cit---", -3, {"citrate": 1})
    complexed = ComponentSpecies("FeCit-", -1, {"Fe2": 1, "citrate": 1})
    reaction = ComponentEquilibriumReaction("Fe2_citrate", {fe: -1, citrate: -1, complexed: 1}, 6.0)
    network = ComponentReactionNetwork.build((fe, citrate, complexed), (reaction,))

    # log[FeCit-] - log[Fe2+] - log[Cit---] = log K = 6
    assert mass_action_residuals(network, [-3.0, -3.0, 0.0]) == pytest.approx([0.0])


def test_ideal_solver_solves_fe_citrate_association() -> None:
    basis = ComponentBasis(("Fe2", "citrate"))
    fe = ComponentSpecies("Fe2+", +2, {"Fe2": 1})
    citrate = ComponentSpecies("Cit---", -3, {"citrate": 1})
    complexed = ComponentSpecies("FeCit-", -1, {"Fe2": 1, "citrate": 1})
    reaction = ComponentEquilibriumReaction("Fe2_citrate", {fe: -1, citrate: -1, complexed: 1}, 6.0)
    tableau = MultiComponentTableau.build(basis, (fe, citrate, complexed), (reaction,))
    network = ComponentReactionNetwork.build((fe, citrate, complexed), (reaction,))
    problem = MultiComponentProblem(
        basis, tableau, network, {"Fe2": 1e-3, "citrate": 1e-3}, strong_charge_mol_L=1e-3,
    )
    result = IdealMassActionSolver().solve(problem)

    assert result["FeCit-"] > 9e-4
    assert result["Fe2+"] == pytest.approx(result["Cit---"], rel=1e-7)
