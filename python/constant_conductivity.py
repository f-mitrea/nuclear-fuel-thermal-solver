from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from linear_solvers import solve_tridiagonal_system

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"


def assemble_constant_conductivity_system(
    p_coeff, q_coeff, r_coeff, dr, n_nodes, T_wall, q_vol, k_UO2
):

    A = np.zeros((n_nodes, n_nodes))
    b = np.zeros(n_nodes)
    gamma = q_vol * dr**2 / k_UO2

    for i in range(n_nodes):
        if i > 0 and i < n_nodes - 1:
            b[i] = dr**2 * r_coeff[i]
            A[i, i] = 2 + dr**2 * q_coeff[i]
            A[i, i + 1] = -1 + dr / 2 * p_coeff[i]
            A[i, i - 1] = -1 - dr / 2 * p_coeff[i]
        elif i == 0:
            b[i] = gamma
            A[i, i] = 4
            A[i, i + 1] = -4
        elif i == n_nodes - 1:
            b[i] = dr**2 * r_coeff[i] + (1 - dr / 2 * p_coeff[i]) * T_wall
            A[i, i] = 2 + dr**2 * q_coeff[i]
            A[i, i - 1] = -1 - dr / 2 * p_coeff[i]
    return A, b


def solve_constant_conductivity(
    k_UO2=2.3, T_wall=800.0, q_vol=400.0e6, R_pellet=4.0e-3, n_nodes=100
):

    dr = R_pellet / n_nodes

    radius = np.zeros(n_nodes)
    for i in range(n_nodes):
        radius[i] = i * dr

    p_coeff = np.zeros(n_nodes)
    p_coeff[1:] = -1.0 / radius[1:]

    r_coeff = q_vol / k_UO2 * np.ones(n_nodes)
    q_coeff = np.zeros(n_nodes)

    A, b = assemble_constant_conductivity_system(
        p_coeff, q_coeff, r_coeff, dr, n_nodes, T_wall, q_vol, k_UO2
    )

    T_num = solve_tridiagonal_system(A, b)

    radius = np.append(radius, R_pellet)
    T_num = np.append(T_num, T_wall)

    T_exact = T_wall + q_vol * (R_pellet**2 - radius**2) / (4 * k_UO2)

    return radius, T_num, T_exact


def main(k_UO2=2.3, T_wall=800.0, q_vol=400.0e6, R_pellet=4.0e-3, n_nodes=100):

    radius, T_num, T_exact = solve_constant_conductivity(
        k_UO2, T_wall, q_vol, R_pellet, n_nodes
    )
    err = np.abs(T_num - T_exact)
    FIGURES_DIR.mkdir(exist_ok=True)

    print(f"Center temperature, numerical : {T_num[0]:.4f} K")
    print(f"Center temperature, analytical : {T_exact[0]:.4f} K")
    print(f"Maximum absolute error : {np.max(err):.6e} K")

    plt.figure()
    plt.plot(radius * 1e3, T_num, "o", label="Numerical")
    plt.plot(radius * 1e3, T_exact, "-", label="Analytical")
    plt.xlabel("Radius [mm]")
    plt.ylabel("Temperature [K]")
    plt.legend()
    plt.grid(True)
    plt.savefig(FIGURES_DIR / "constant_k_profile.png", dpi=150)

    plt.figure()
    plt.plot(radius * 1e3, err, "-")
    plt.xlabel("Radius [mm]")
    plt.ylabel("|T_num - T_exact| [K]")
    plt.grid(True)
    plt.savefig(FIGURES_DIR / "constant_k_error.png", dpi=150)

    plt.show()


if __name__ == "__main__":
    main()
