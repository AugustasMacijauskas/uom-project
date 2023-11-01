# AUTOGENERATED! DO NOT EDIT! File to edit: ../nbs/01_poisson_solvers.ipynb.

# %% auto 0
__all__ = ['setup_poisson_problem', 'poisson_gauss_seidel_with_sor_solver', 'poisson_non_iterative_solver',
           'poisson_newton_solver', 'poisson_newton_alternative_solver']

# %% ../nbs/01_poisson_solvers.ipynb 5
from functools import partial

import numpy as np
import scipy
from scipy import sparse


# %% ../nbs/01_poisson_solvers.ipynb 8
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

# %% ../nbs/01_poisson_solvers.ipynb 10
def poisson_gauss_seidel_with_sor_solver(
    w, r=None, verbose=False, log_middle_values=False,
):
    N = w.shape[0] - 1
    nx = ny = N + 1

    h = 1 / N

    # Handle the SOR parameter
    if r is None:
        r = 2 / (1 + np.pi / N) # optimal value

    psi = np.zeros((nx, ny)) # streamfunction

    middle_values = []

    for iteration in range(1, 4 * nx + 1):
        for i in range(1, nx - 1):
            for j in range(1, ny - 1):
                psi[i, j] = (1 - r) * psi[i, j] + r * (
                    psi[i - 1, j] + psi[i + 1, j] +
                    psi[i, j - 1] + psi[i, j + 1] +
                    w[i, j] * h ** 2
                ) / 4
        
        if verbose:
            print(f"{iteration=}; psi(0.5, 0.5) = {psi[nx // 2, ny // 2]}")
        
        if log_middle_values:
            middle_values.append(psi[nx // 2, ny // 2])
    
    return psi, iteration, middle_values


# %% ../nbs/01_poisson_solvers.ipynb 25
def construct_laplacian_kernel_matrix(N, h=None):
    '''
    Construct the matrix that defines the linear system in sparse form
    '''

    if h is None:
        h = 1 / (N + 1)

    index_matrix = np.arange(N ** 2).reshape((N, N))
    all_indices = np.arange(N ** 2)
    
    rows = [np.arange(N ** 2)]
    cols = [np.arange(N ** 2)]
    values = [[-4 / h ** 2] * N ** 2]

    # Create grid
    x_grid, y_grid = np.meshgrid(
        np.arange(0, N), np.arange(0, N), indexing="ij"
    )

    # Masks for validity
    valid_above = x_grid > 0
    valid_below = x_grid < N - 1
    valid_left = y_grid > 0
    valid_right = y_grid < N - 1

    # Get the indices
    rows.append(all_indices[valid_above.flatten()])
    cols.append(index_matrix[(
        (x_grid - 1)[valid_above].flatten(), y_grid[valid_above].flatten()
    )].flatten())
    values.append([1 / h ** 2] * (N ** 2 - N))

    rows.append(all_indices[valid_below.flatten()])
    cols.append(index_matrix[(
        (x_grid + 1)[valid_below].flatten(), y_grid[valid_below].flatten()
    )].flatten())
    values.append([1 / h ** 2] * (N ** 2 - N))

    rows.append(all_indices[valid_left.flatten()])
    cols.append(index_matrix[(
        x_grid[valid_left].flatten(), (y_grid - 1)[valid_left].flatten()
    )].flatten())
    values.append([1 / h ** 2] * (N ** 2 - N))

    rows.append(all_indices[valid_right.flatten()])
    cols.append(index_matrix[(
        x_grid[valid_right].flatten(), (y_grid + 1)[valid_right].flatten()
    )].flatten())
    values.append([1 / h ** 2] * (N ** 2 - N))

    # Flatten the lists
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    values = np.concatenate(values)

    # Create the sparse matrix from the above information
    return sparse.csr_matrix((values, (rows, cols)), shape=(N ** 2, N ** 2))


# %% ../nbs/01_poisson_solvers.ipynb 26
SPARSE_ALGORITHM_DICT = {
    "base": sparse.linalg.spsolve,
    "cg": sparse.linalg.cg,
    "bicgstab": sparse.linalg.bicgstab,
}

def solve_sparse_linear_system(A, b, algorithm="base"):
    """
    Solve a sparse linear system using the specified algorithm.
    """

    solver = SPARSE_ALGORITHM_DICT[algorithm]

    if algorithm == "base": b = sparse.csr_matrix(b[:, None])

    soln = solver(A, b)

    return soln if algorithm == "base" else soln[0]


# %% ../nbs/01_poisson_solvers.ipynb 27
def poisson_non_iterative_solver(w, algorithm="bicgstab"):
    N = w.shape[0] - 2

    # Construct the matrix that defines the linear system
    kernel_matrix = construct_laplacian_kernel_matrix(N)
    
    # Cast vorticity to the required form
    w = -w[1:-1, 1:-1].flatten()

    # Solve the linear system
    psi = solve_sparse_linear_system(A=kernel_matrix, b=w, algorithm=algorithm)
    
    psi = psi.reshape(N, N)

    psi = np.pad(psi, (1, 1), mode="constant")
    
    return psi


# %% ../nbs/01_poisson_solvers.ipynb 34
def get_jacobian(N):
    # Set h=1 and note that we actually compute h ** 2 * f_current below.
    # This is to done have better scaling.
    return construct_laplacian_kernel_matrix(N, h=1)
    

def f(x, w):
    # this N is different from before
    N = int(np.sqrt(x.shape[0]))
    h = 1 / (N + 1)
    
    x = x.reshape(N, N)
    x = np.pad(x, (1, 1), mode="constant", constant_values=0)
    
    f = -4 * x[1:-1, 1:-1] + x[2:, 1:-1] + x[:-2, 1:-1] + x[1:-1, 2:] + x[1:-1, :-2]
    f = f + h ** 2 * w[1:-1, 1:-1] # note we compute h ** 2 * f_current
    
    return f.reshape((-1, 1))


def newton_iterator(
    w, f, get_jacobian,
    algorithm="bicgstab", TOL=1e-8, max_iter=10, quiet=True
):
    '''
        - w: vorticity
        - f: evaluates the function given x, w, h
        - get_jacobian: evaluates the Jacobian given N
    '''

    N = w.shape[0] - 1
    
    n_iter = 0 # number of iterations

    # Initialization
    x = np.zeros(((N - 1) ** 2, 1))
    
    # Initialization
    f_current = f(x, w)
    
    # Check if the initial guess is a solution
    f_norm = scipy.linalg.norm(f_current)
    if f_norm <= TOL:
        if not quiet:
            print(f"n_iter={n_iter}")

        return x, n_iter
    
    
    while n_iter < max_iter:
        n_iter += 1
        # Set h=1 to have better scaling and note that we actually compute h ** 2 * f_current
        jacobian = sparse.csr_matrix(get_jacobian(N=N-1)) # Sparsify the jacobian
        
        # dx = scipy.sparse.linalg.spsolve(jacobian, -f_current).reshape((-1, 1))
        dx = solve_sparse_linear_system(
            A=jacobian, b=-f_current, algorithm=algorithm
        ).reshape((-1, 1))
        x_next = x + dx
        
        f_current = f(x_next, w)
        
        f_norm = scipy.linalg.norm(f_current)
        if not quiet:
            print(f"iter={n_iter}; |residual|={f_norm}; |dx|={scipy.linalg.norm(dx)}")
        if f_norm <= TOL:
            break
        
        x = x_next

    if not quiet:
        print(f"n_iter={n_iter}")
    
    return x_next, n_iter


# %% ../nbs/01_poisson_solvers.ipynb 35
def poisson_newton_solver(
    w, algorithm="bicgstab", TOL=1e-8, max_iter=10, quiet=True
):
    solution, _ = newton_iterator(
        w=w, f=f, get_jacobian=get_jacobian,
        algorithm=algorithm,
        TOL=TOL, max_iter=max_iter, quiet=quiet
    )
    
    N = w.shape[0] - 1
    psi = solution.reshape(N - 1, N - 1)

    psi = np.pad(psi, (1, 1), mode="constant", constant_values=0)
    
    return psi


# %% ../nbs/01_poisson_solvers.ipynb 36
def alternative_jacobian(x):
    N = int(np.sqrt(x.shape[0]))
    return get_jacobian(N)


def alternative_f(x, w):
    return f(x, w).flatten()


# %% ../nbs/01_poisson_solvers.ipynb 37
def poisson_newton_alternative_solver(w, **kwargs):
    N = w.shape[0] - 1

    solution = scipy.optimize.root(
        fun=partial(alternative_f, w=w), x0=np.zeros(((N - 1) ** 2, )),
        # jac=alternative_jacobian,
        **kwargs
    )
    
    psi = solution.x.reshape(N - 1, N - 1)

    psi = np.pad(psi, (1, 1), mode="constant", constant_values=0)
    
    return psi, solution

