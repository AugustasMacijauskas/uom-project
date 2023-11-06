# AUTOGENERATED! DO NOT EDIT! File to edit: ../nbs/00_core.ipynb.

# %% auto 0
__all__ = ['setup_poisson_problem', 'solve_sparse_linear_system', 'solve_sparse_linear_system_pytorch', 'calculate_force',
           'calculate_velocity']

# %% ../nbs/00_core.ipynb 4
import numpy as np
from scipy import sparse

import torch

# %% ../nbs/00_core.ipynb 6
def setup_poisson_problem(N):
    """
    Setup a 2D Poisson problem. Only square domains are supported.
    """

    nx = ny = N + 1

    # Initialize grid
    x_grid, y_grid = np.meshgrid(
        np.linspace(0, 1, nx), np.linspace(0, 1, ny), indexing="ij"
    )

    # Initialize the value of the vorticity on the grid
    w = 2 * np.pi ** 2 * np.sin(np.pi * x_grid) * np.sin(np.pi * y_grid)

    exact_solution = np.sin(np.pi * x_grid) * np.sin(np.pi * y_grid)

    return w, exact_solution, nx, ny, x_grid, y_grid

# %% ../nbs/00_core.ipynb 7
SPARSE_ALGORITHM_DICT = {
    "base": sparse.linalg.spsolve,
    "cg": sparse.linalg.cg,
    "bicgstab": sparse.linalg.bicgstab,
}

# %% ../nbs/00_core.ipynb 8
def solve_sparse_linear_system(A, b, algorithm="base", **kwargs):
    """
    Solve a sparse linear system using the specified algorithm.
    """

    solver = SPARSE_ALGORITHM_DICT[algorithm]

    if algorithm == "base" and len(b.shape) == 1:
        b = sparse.csr_matrix(b[:, None])

    soln = solver(A, b, **kwargs)

    return soln if algorithm == "base" else soln[0]


# %% ../nbs/00_core.ipynb 9
def lu_solve(A, b):
    LU, pivots = torch.linalg.lu_factor(A)
    return torch.linalg.lu_solve(LU, pivots, b)


TORCH_ALGORITHM_DICT = {
    "base": torch.linalg.solve,
    "lu_solve": lu_solve,
    "lstsq": torch.linalg.lstsq,
}


# %% ../nbs/00_core.ipynb 10
def solve_sparse_linear_system_pytorch(A, b, algorithm="base", **kwargs):
    """
    Solve a sparse linear system using the specified algorithm.
    """

    solver = TORCH_ALGORITHM_DICT[algorithm]

    out = solver(A, b, **kwargs)

    if algorithm == "lstsq": out = out.solution

    return out


# %% ../nbs/00_core.ipynb 11
def calculate_force(psi):
    h = 1 / (psi.shape[0] - 1)
    ny = psi.shape[1]
    
    return np.sum(
        2 * psi[1:-1, ny - 1] -
        5 * psi[1:-1, ny - 2] +
        4 * psi[1:-1, ny - 3] -
        psi[1:-1, ny - 4]
    ) / h

# %% ../nbs/00_core.ipynb 12
def calculate_velocity(psi):
    h = 1 / (psi.shape[0] - 1)

    nx = psi.shape[0]
    
    u = np.zeros_like(psi)
    v = np.zeros_like(psi)
    
    u[1:-1, 1:-1] = (psi[1:-1, 2:] - psi[1:-1, :-2]) / (2 * h)
    v[1:-1, 1:-1] = -(psi[2:, 1:-1] - psi[:-2, 1:-1]) / (2 * h)
    
    u[1:-1, -1] = np.sin(np.pi * np.arange(1, nx - 1) * h) ** 2
    
    return u, v 
