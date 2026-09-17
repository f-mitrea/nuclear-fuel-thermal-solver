from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from linear_solvers import solve_tridiagonal_system

K_DEFAULT = "1/(0.0375 + 2.165e-4*T)"
NODES_DEFAULT = [25, 50, 100, 200, 400]
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"


def compute_pellet_residual(T, radius, dr, q_vol, T_wall, k_UO2, dk_dT):

    N = len(T)
    F = np.zeros(N)

    F[0] = 4 * k_UO2(T[0]) * (T[1] - T[0]) / dr**2 + q_vol

    for i in range(1, N - 1):
        dT_dr = (T[i + 1] - T[i - 1]) / (2 * dr)
        d2T_dr2 = (T[i + 1] - 2 * T[i] + T[i - 1]) / dr**2
        F[i] = (
            k_UO2(T[i]) * d2T_dr2
            + k_UO2(T[i]) / radius[i] * dT_dr
            + dk_dT(T[i]) * dT_dr**2
            + q_vol
        )

    F[N - 1] = T[N - 1] - T_wall

    return F


def assemble_pellet_jacobian(T, radius, dr, k_UO2, dk_dT, d2k_dT2):

    N = len(T)
    J = np.zeros((N, N))

    J[0, 0] = (4 * dk_dT(T[0]) * (T[1] - T[0]) - 4 * k_UO2(T[0])) / dr**2
    J[0, 1] = 4 * k_UO2(T[0]) / dr**2
    J[N - 1, N - 1] = 1

    for i in range(1, N - 1):
        dT_dr = (T[i + 1] - T[i - 1]) / (2 * dr)
        d2T_dr2 = (T[i + 1] - 2 * T[i] + T[i - 1]) / dr**2

        J[i, i - 1] = (
            k_UO2(T[i]) / dr**2
            - k_UO2(T[i]) / (2 * dr * radius[i])
            - dk_dT(T[i]) * dT_dr / dr
        )
        J[i, i + 1] = (
            k_UO2(T[i]) / dr**2
            + k_UO2(T[i]) / (2 * dr * radius[i])
            + dk_dT(T[i]) * dT_dr / dr
        )
        J[i, i] = (
            -2 * k_UO2(T[i]) / dr**2
            + dk_dT(T[i]) * d2T_dr2
            + dk_dT(T[i]) / radius[i] * dT_dr
            + d2k_dT2(T[i]) * dT_dr**2
        )

    return J


def solve_pellet_newton(
    T_old, radius, dr, tol, max_iter, k_UO2, dk_dT, d2k_dT2, q_vol, T_wall
):

    n_iter = 0
    err = 1.0
    T_new = T_old

    while np.max(np.abs(err)) > tol and n_iter < max_iter:
        J = assemble_pellet_jacobian(T_old, radius, dr, k_UO2, dk_dT, d2k_dT2)
        F = compute_pellet_residual(T_old, radius, dr, q_vol, T_wall, k_UO2, dk_dT)
        b = -F

        delta_T = solve_tridiagonal_system(J, b)

        T_new = T_old + delta_T
        err = T_new - T_old
        T_old = T_new
        n_iter = n_iter + 1
    final_err = np.max(np.abs(err))
    if final_err > tol:
        raise RuntimeError(
            f"Newton did not converge: {n_iter} iterations, "
            f"last increment {final_err:.3e}, tolerance {tol:.3e}"
        )

    return T_new


def build_conductivity_functions(k_expr):

    T_sym = sp.Symbol("T")

    dk_expr = sp.diff(k_expr, T_sym)
    d2k_expr = sp.diff(k_expr, T_sym, 2)

    k_UO2 = sp.lambdify(T_sym, k_expr, "numpy")
    dk_dT = sp.lambdify(T_sym, dk_expr, "numpy")
    d2k_dT2 = sp.lambdify(T_sym, d2k_expr, "numpy")

    return k_UO2, dk_dT, d2k_dT2, dk_expr, d2k_expr


def kirchhoff_error(k_expr, T, radius, T_wall, q_vol, R_pellet):
    T_sym = sp.Symbol("T")

    K_expr = sp.integrate(k_expr, T_sym)
    K_func = sp.lambdify(T_sym, K_expr, "numpy")

    theta_num = K_func(T) - K_func(T_wall)
    theta_exact = q_vol * (R_pellet**2 - radius**2) / 4

    error = np.max(np.abs(theta_num - theta_exact)) / np.max(np.abs(theta_exact))

    return error


def solve_variable_conductivity(
    k_expr, T_wall, q_vol, R_pellet, n_nodes, tol, max_iter
):
    if n_nodes < 3:
        raise ValueError(f"n_nodes must be at least 3, got {n_nodes}")

    if R_pellet <= 0:
        raise ValueError(f"R_pellet must be positive, got {R_pellet}")

    k_UO2, dk_dT, d2k_dT2, _, _ = build_conductivity_functions(k_expr)

    dr = R_pellet / (n_nodes - 1)

    radius = np.zeros(n_nodes)
    for i in range(n_nodes):
        radius[i] = i * dr

    T = np.zeros(n_nodes)
    for i in range(n_nodes):
        T[i] = T_wall + q_vol * (R_pellet**2 - radius[i] ** 2) / (4 * k_UO2(T_wall))

    T = solve_pellet_newton(
        T, radius, dr, tol, max_iter, k_UO2, dk_dT, d2k_dT2, q_vol, T_wall
    )

    return radius, T, dr


def convergence_study(k_expr, T_wall, q_vol, R_pellet, tol, max_iter, node_counts):
    dr_values = []
    errors = []

    for n_nodes in node_counts:
        radius, T, dr = solve_variable_conductivity(
            k_expr, T_wall, q_vol, R_pellet, n_nodes, tol, max_iter
        )

        errors.append(kirchhoff_error(k_expr, T, radius, T_wall, q_vol, R_pellet))
        dr_values.append(dr)

    dr_values = np.array(dr_values)
    errors = np.array(errors)

    observed_orders = np.log(errors[:-1] / errors[1:]) / np.log(
        dr_values[:-1] / dr_values[1:]
    )

    return dr_values, errors, observed_orders


def main(
    k_expr=K_DEFAULT,
    T_wall=800.0,
    q_vol=400.0e6,
    R_pellet=4.0e-3,
    n_nodes=100,
    tol=1.0e-8,
    max_iter=20,
    node_counts=NODES_DEFAULT,
):

    k_expr = sp.sympify(k_expr)

    _, _, _, dk_expr, d2k_expr = build_conductivity_functions(k_expr)

    print("Conductivity law and derivatives obtained by symbolic differentiation")
    print(f"k(T)    = {k_expr}")
    print(f"dk/dT   = {sp.simplify(dk_expr)}")
    print(f"d2k/dT2 = {sp.simplify(d2k_expr)}")

    radius, T, _ = solve_variable_conductivity(
        k_expr, T_wall, q_vol, R_pellet, n_nodes, tol, max_iter
    )

    print(f"\nCentreline temperature : {T[0]:.4f} K")
    print(f"Wall temperature : {T[-1]:.4f} K")

    dr_values, errors, observed_orders = convergence_study(
        k_expr, T_wall, q_vol, R_pellet, tol, max_iter, node_counts
    )

    print("\nGrid convergence, error measured on the conduction integral")
    print("   nodes     dr [m]      rel. error      order")

    print(f"   {node_counts[0]:4d}   {dr_values[0]:.3e}   {errors[0]:.3e}")

    for i in range(1, len(node_counts)):
        print(
            f"   {node_counts[i]:4d}   {dr_values[i]:.3e}   {errors[i]:.3e}   {observed_orders[i-1]:6.2f}"
        )

    FIGURES_DIR.mkdir(exist_ok=True)

    plt.figure()
    plt.plot(radius * 1e3, T, "-o", markersize=3)
    plt.xlabel("Radius [mm]")
    plt.ylabel("Temperature [K]")
    plt.title("Temperature profile, temperature-dependent conductivity")
    plt.grid(True)
    plt.savefig(FIGURES_DIR / "variable_k_profile.png", dpi=150)

    plt.figure()
    plt.loglog(dr_values, errors, "o-", label="Kirchhoff error")
    plt.xlabel("Grid spacing dr [m]")
    plt.ylabel("Relative error on the conduction integral")
    plt.legend()
    plt.grid(True, which="both")
    plt.savefig(FIGURES_DIR / "variable_k_convergence.png", dpi=150)

    plt.show()


if __name__ == "__main__":
    main()
