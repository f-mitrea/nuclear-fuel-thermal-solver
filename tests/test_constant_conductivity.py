import numpy as np
import pytest
from constant_conductivity import solve_constant_conductivity


def test_centreline_temperature():
    radius, T_num, T_exact = solve_constant_conductivity()

    assert T_num[0] == pytest.approx(860.8696, abs=1e-4)


@pytest.mark.parametrize("n_nodes", [10, 100, 400])
def test_matches_parabola_at_any_grid(n_nodes):
    radius, T_num, T_exact = solve_constant_conductivity(n_nodes=n_nodes)

    assert np.max(np.abs(T_num - T_exact)) < 1e-9
