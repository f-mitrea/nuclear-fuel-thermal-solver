import numpy as np
import pytest
import sympy as sp

from variable_conductivity import (
    K_DEFAULT,
    assemble_pellet_jacobian,
    build_conductivity_functions,
    compute_pellet_residual,
    convergence_study,
    solve_variable_conductivity,
)

T_WALL = 800.0
Q_VOL = 400.0e6
R_PELLET = 4.0e-3
TOL = 1.0e-8
MAX_ITER = 20

# k(T) = 1/(A + B T)
A_COEF = 0.0375
B_COEF = 2.165e-4


def kirchhoff_exact(radius):
    arg = B_COEF * Q_VOL * (R_PELLET**2 - radius**2) / 4
    return ((A_COEF + B_COEF * T_WALL) * np.exp(arg) - A_COEF) / B_COEF


def test_matches_kirchhoff_closed_form():
    radius, T, dr = solve_variable_conductivity(
        sp.sympify(K_DEFAULT), T_WALL, Q_VOL, R_PELLET, 100, TOL, MAX_ITER
    )

    assert T[0] == pytest.approx(1202.8913, abs=1e-3)
    assert np.max(np.abs(T - kirchhoff_exact(radius))) < 2e-2


def test_constant_k_gives_parabola():
    k = 2.3
    radius, T, dr = solve_variable_conductivity(
        sp.sympify(k), T_WALL, Q_VOL, R_PELLET, 100, TOL, MAX_ITER
    )
    T_exact = T_WALL + Q_VOL * (R_PELLET**2 - radius**2) / (4 * k)

    assert np.max(np.abs(T - T_exact)) < 1e-9


def test_grid_convergence_is_second_order():
    _, _, orders = convergence_study(
        sp.sympify(K_DEFAULT),
        T_WALL,
        Q_VOL,
        R_PELLET,
        TOL,
        MAX_ITER,
        [25, 50, 100, 200],
    )

    assert np.all(np.abs(orders - 2.0) < 0.05)


def test_jacobian_matches_finite_differences():
    k, dk, d2k, _, _ = build_conductivity_functions(sp.sympify(K_DEFAULT))
    n = 20
    dr = R_PELLET / (n - 1)
    radius = np.arange(n) * dr
    rng = np.random.default_rng(1)
    T = kirchhoff_exact(radius) + rng.uniform(-5, 5, n)
    v = rng.uniform(-1, 1, n)
    eps = 1e-4

    F_plus = compute_pellet_residual(T + eps * v, radius, dr, Q_VOL, T_WALL, k, dk)
    F_minus = compute_pellet_residual(T - eps * v, radius, dr, Q_VOL, T_WALL, k, dk)
    fd = (F_plus - F_minus) / (2 * eps)

    J = assemble_pellet_jacobian(T, radius, dr, k, dk, d2k)

    assert np.allclose(J @ v, fd, rtol=1e-6, atol=1e-3)


def test_too_few_nodes_raises():
    with pytest.raises(ValueError):
        solve_variable_conductivity(
            sp.sympify(K_DEFAULT), T_WALL, Q_VOL, R_PELLET, 2, TOL, MAX_ITER
        )


def test_non_positive_radius_raises():
    with pytest.raises(ValueError):
        solve_variable_conductivity(
            sp.sympify(K_DEFAULT), T_WALL, Q_VOL, 0.0, 100, TOL, MAX_ITER
        )


def test_newton_not_converged_raises():
    with pytest.raises(RuntimeError):
        solve_variable_conductivity(
            sp.sympify(K_DEFAULT), T_WALL, Q_VOL, R_PELLET, 100, TOL, 1
        )
