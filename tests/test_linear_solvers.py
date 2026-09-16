import numpy as np
from linear_solvers import solve_tridiagonal_system


def test_thomas_matches_numpy():
    rng = np.random.default_rng(0)
    n = 50
    lower = rng.uniform(-1, 1, n - 1)
    upper = rng.uniform(-1, 1, n - 1)
    diag = 4 + rng.uniform(0, 1, n)
    A = np.diag(diag) + np.diag(lower, -1) + np.diag(upper, 1)
    b = rng.uniform(-1, 1, n)

    x = solve_tridiagonal_system(A, b)

    assert np.allclose(x, np.linalg.solve(A, b), rtol=0, atol=1e-12)
