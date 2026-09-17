import numpy as np


def solve_tridiagonal_system(A, B):
    N = A.shape[0]

    X = np.zeros(N)
    Y = np.zeros(N)
    u = np.zeros(N)
    v = np.zeros(N)
    alphas = np.zeros(N)

    for i in range(N - 1):
        v[i] = A[i, i + 1]

    u[0] = A[0, 0]

    for j in range(1, N):
        alphas[j] = A[j, j - 1] / u[j - 1]
        u[j] = A[j, j] - alphas[j] * v[j - 1]

    X[0] = B[0]

    for k in range(1, N):
        X[k] = B[k] - alphas[k] * X[k - 1]

    Y[N - 1] = X[N - 1] / u[N - 1]

    for l in range(N - 2, -1, -1):
        Y[l] = (X[l] - v[l] * Y[l + 1]) / u[l]

    return Y


if __name__ == "__main__":
    A = np.array(
        [
            [4.0, -1.0, 0.0, 0.0],
            [-1.0, 4.0, -1.0, 0.0],
            [0.0, -1.0, 4.0, -1.0],
            [0.0, 0.0, -1.0, 4.0],
        ]
    )

    B = np.array([5.0, 5.0, 5.0, 5.0])

    Y_thomas = solve_tridiagonal_system(A, B)
    Y_numpy = np.linalg.solve(A, B)

    print(f"Thomas algorithm  : {Y_thomas}")
    print(f"NumPy reference   : {Y_numpy}")
    print(f"Max abs deviation : {np.max(np.abs(Y_thomas - Y_numpy)):.3e}")
