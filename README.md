# nuclear-fuel-thermal-solver

This is a small solver for the temperature inside a nuclear fuel pellet. It's 1D (only the radius matters), steady state, and uses finite differences.

I started it for three reasons. I'd just passed a numerical methods exam and found the subject more interesting than I expected, so I wanted to use those methods on an actual physical problem instead of textbook exercises. I also want to specialise in nuclear engineering, and a fuel pellet is about the simplest nuclear problem that still has some real physics in it. And I wanted to get better at Python, which I knew much less well than MATLAB.

That's also why there's a MATLAB folder. I wrote the first version in MATLAB because that's where I was comfortable after the exam, got it working, and then ported it to Python. The Python code is the one I test and keep updated. The MATLAB files are my original prototype and I've left them exactly as they were.

## The problem

A UO2 pellet generates heat everywhere inside it (q''', in W/m^3) and that heat has to leave through the outer surface. At steady state the temperature T(r) satisfies

```
(1/r) d/dr ( r k(T) dT/dr ) + q''' = 0
```

The centre is a symmetry point, so dT/dr = 0 there. The surface at r = R is held at a fixed T_wall.

I solved it twice.

**Constant conductivity.** If k is a number, the equation is linear and the exact answer is a parabola:

```
T(r) = T_wall + q''' (R^2 - r^2) / (4k)
```

**Conductivity that depends on temperature.** For real UO2, k drops as the fuel gets hotter. I used

```
k(T) = 1 / (0.0375 + 2.165e-4 T)      [W/m/K, T in K]
```

which is the phonon term of the correlation recommended by Harding and Martin (J. Nucl. Mater. 166 (1989) 223-226). Their full formula has a second term for high temperatures, but at the temperatures in this pellet it's below 1e-6 W/m/K against a k of about 5, so I left it out.

With k(T) the equation becomes nonlinear. What made this case nice to work on is that it still has an exact solution. If you use the Kirchhoff transform K(T) = integral of k dT, the equation turns back into the linear one, K(T) - K(T_wall) = q''' (R^2 - r^2)/4, and for this particular k(T) you can invert it by hand:

```
T(r) = [ (A + B T_wall) exp( B q''' (R^2 - r^2) / 4 ) - A ] / B      A = 0.0375, B = 2.165e-4
```

So I can check the numerical solution against an exact one in both cases, which was the whole point.

The default numbers are R = 4 mm, q''' = 150 MW/m^3 and T_wall = 600 K.

## How the code works

I discretise the equation with second-order central differences on a uniform grid. The centre needs special handling because of the 1/r term. Taking the limit as r goes to 0 gives 4 k (T1 - T0)/dr^2 + q''' = 0, and that's the equation I use for the first node.

In the constant-k case the result is a tridiagonal linear system. I solve it with the Thomas algorithm, which I wrote myself in `linear_solvers.py`.

In the variable-k case I use Newton's method. I wrote the residual and the Jacobian by hand. The derivatives of k(T) are computed by SymPy, so if you want a different conductivity law you only have to change the expression. Newton stops when the temperature change between two iterations is below the tolerance, and if it doesn't get there within the maximum number of iterations it raises an error instead of returning a wrong answer quietly.

The folders:

```
python/
  linear_solvers.py            Thomas algorithm
  constant_conductivity.py     constant k, compared with the parabola
  variable_conductivity.py     variable k, Newton, convergence study
tests/                         12 pytest tests
figures/                       the plots, created by the scripts
matlab_prototype/              my first MATLAB version
```

## Running it

I've only tried this on Windows with Python 3.12.

```
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python -m pytest
```

For me that's 12 passed in about 15 seconds. The tests check the centre temperatures against the exact values, the Newton result against the Kirchhoff solution, that the convergence order is really 2, that my hand-written Jacobian agrees with one computed by finite differences, and that wrong inputs (too few nodes, zero radius, Newton not converging) give an error.

## Results

To get the numbers and the plots:

```
python python\constant_conductivity.py
python python\variable_conductivity.py
```

Each script prints its results and saves the figures into `figures/`.

With constant k (2.3 W/m/K) the centre temperature is 860.8696 K, the same as the parabola, and the difference stays below 1e-9 K on 10, 100 and 400 nodes. That looks too good, and there's a reason for it: that a central difference is exact on a parabola, so there's no discretisation error at all, just rounding. It tells me the matrix is assembled correctly, but it doesn't say anything about how the error shrinks with the grid.

![Constant k, difference from the exact solution](figures/constant_k_error.png)

The variable-k case is where that shows up. The centre temperature is 707.2566 K. That's lower than in the constant-k case because with this law k is between 5.2 and 6.0 W/m/K inside the pellet, more than twice the 2.3 I'd assumed. Here's what happens to the error as I refine the grid. I measure it on the Kirchhoff variable, as a relative error:

```
nodes   dr [m]      rel. error   order
  25    1.667e-04   1.734e-04
  50    8.163e-05   4.163e-05    2.00
 100    4.040e-05   1.020e-05    2.00
 200    2.010e-05   2.525e-06    2.00
 400    1.003e-05   6.281e-07    2.00
```

Each time dr halves, the error gets four times smaller, which is exactly what second order means.

![Convergence with variable k](figures/variable_k_convergence.png)

![Temperature profile with variable k](figures/variable_k_profile.png)

## What's still wrong with it

There are a few things I know aren't right yet.

My Thomas function takes the whole N x N matrix even though it only uses three diagonals. It works, but it wastes memory for no reason, and passing three vectors would be the proper way.

`n_nodes` means different things in the two solvers. With constant k, dr = R/n_nodes and the wall point gets added at the end, so you have n_nodes + 1 points. With variable k, dr = R/(n_nodes - 1) and you have exactly n_nodes.

Newton only checks how much the temperature changed, not whether the residual is actually small. Most of the time those go together, but not always.

The convergence error is computed on the Kirchhoff variable and not on T directly. The order is the same because one is a monotone function of the other, but it would be cleaner to check T.

The plotting code sits in the same files as the solvers, so you can't use a solver without also importing matplotlib.

The tests take 15 seconds, which feels slow for 12 tests. I haven't looked into why yet.

On the physics side it's very simplified. It's steady state only, the heat generation is uniform, there's no gap or cladding (the pellet surface temperature is just given), and k doesn't depend on burnup or porosity. Also, Harding and Martin give their correlation for 773 to 3120 K, and my pellet sits between 600 and 707 K, so I'm using it a bit below its range.

The MATLAB prototype has a couple of typos in file and folder names (`costant`, `coductivity`) and hard-coded paths. I'm leaving it as it was on purpose.

If I work on this again, I'll first fix Thomas and the `n_nodes` mismatch, because they're small changes and the tests I already have would tell me if I broke something. After that I'd add the residual check to Newton and move the plots into their own script. The physics comes last, starting with the gap and the cladding, which I'll need anyway if I couple this with reactor kinetics.

## Author

Fabio Mitrea, third-year Energy Engineering student at Sapienza University of Rome.
