# AUTOGENERATED! DO NOT EDIT! File to edit: ../nbs/03_streamfunction_vorticity_newton.ipynb.

# %% auto 0
__all__ = ['get_standard_basis_vector', 'get_jacobian', 'f', 'reconstruct_w', 'newton_iterator', 'newton_solver']

# %% ../nbs/03_streamfunction_vorticity_newton.ipynb 6
from . import core, poisson_solvers, streamfunction_vorticity

import numpy as np
import scipy
from scipy import sparse

# %% ../nbs/03_streamfunction_vorticity_newton.ipynb 9
# Given the vorticity, solve the Poisson eqn. to find the streamfunction
def get_standard_basis_vector(size, i):
    vec = np.zeros((size, 1))
    vec[i] = 1.0
    
    return vec


def get_jacobian(f, x, Re, kernel_matrix):
    N = int(np.sqrt(x.shape[0] // 2 + 1))
    h = 1 / N

    f_evaluated = f(x=x, Re=Re, kernel_matrix=kernel_matrix)
    
    # Jacobian is sparse
    return sparse.csr_matrix(
        np.hstack([(
            f(
                x=x + h*get_standard_basis_vector(size=x.shape[0], i=i), Re=Re,
                kernel_matrix=kernel_matrix
            ) -
            f_evaluated
        ) for i in range(x.shape[0])])
    )


def f(x, Re, U_wall_top, kernel_matrix):
    N = int(np.sqrt(x.shape[0] // 2 + 1))
    h = 1 / N

    psi = x[:(N-1)**2]
    w_left   = x[(N-1)**2 + 0*(N-1) : (N-1)**2 + 1*(N-1)][:, 0]
    w_right  = x[(N-1)**2 + 1*(N-1) : (N-1)**2 + 2*(N-1)][:, 0]
    w_bottom = x[(N-1)**2 + 2*(N-1) : (N-1)**2 + 3*(N-1)][:, 0]
    w_top    = x[(N-1)**2 + 3*(N-1) : (N-1)**2 + 4*(N-1)][:, 0]
    w_middle = x[(N-1)**2 + 4*(N-1) :]

    # Calculate the equations coming from the Poisson equation
    f_poisson = kernel_matrix @ psi
    f_poisson = f_poisson + h ** 2 * w_middle

    psi = psi.reshape(N-1, N-1)

    # Calculate contributions coming from the vorticity transport equation
    w_middle = w_middle.reshape(N-1, N-1)
    
    # Calculate the sides first
    # y = 0, U_wall = 0
    f_w_bottom = h ** 2 * (w_middle[:, 0] + 3 * w_bottom) + 8 * psi[:, 0]
    # y = 1, U_wall is known here
    f_w_top = h ** 2 * (w_middle[:, -1] + 3 * w_top) + 8 * (
        h * U_wall_top + psi[:, -1]
    )
    # x = 0
    f_w_left = h ** 2 * (w_middle[0, :] + 3 * w_left) + 8 * psi[0, :]
    # x = 1
    f_w_right = h ** 2 * (w_middle[-1, :] + 3 * w_right) + 8 * psi[-1, :]

    f_w_middle = -4 * w_middle
    f_w_middle[:-1, :] += w_middle[1:, :]
    f_w_middle[-1:, :] += w_right
    f_w_middle[1:, :] += w_middle[:-1, :]
    f_w_middle[:1, :] += w_left
    f_w_middle[:, :-1] += w_middle[:, 1:]
    f_w_middle[:, -1] += w_top
    f_w_middle[:, 1:] += w_middle[:, :-1]
    f_w_middle[:, 0] += w_bottom

    f_w_middle[1:-1, 1:-1] += Re * (
        (psi[2:, 1:-1] - psi[:-2, 1:-1]) * (w_middle[1:-1, 2:] - w_middle[1:-1, :-2]) -
        (psi[1:-1, 2:] - psi[1:-1, :-2]) * (w_middle[2:, 1:-1] - w_middle[:-2, 1:-1])
    ) / 4
    f_w_middle[:1, 1:-1] += Re * (
        psi[1, 1:-1] * (w_middle[0, 2:] - w_middle[0, :-2]) -
        (psi[0, 2:] - psi[0, :-2]) * (w_middle[1, 1:-1] - w_left[1:-1])
    ) / 4
    f_w_middle[-1:, 1:-1] -= Re * (
        psi[-2, 1:-1] * (w_middle[-1, 2:] - w_middle[-1, :-2]) +
        (psi[-1, 2:] - psi[-1, :-2]) * (w_right[1:-1] - w_middle[-2, 1:-1])
    ) / 4
    f_w_middle[1:-1, 0] += Re * (
        (psi[2:, 0] - psi[:-2, 0]) * (w_middle[1:-1, 1] - w_bottom[1:-1]) -
        psi[1:-1, 1] * (w_middle[2:, 0] - w_middle[:-2, 0])
    ) / 4
    f_w_middle[1:-1, -1] += Re * (
        (psi[2:, -1] - psi[:-2, -1]) * (w_top[1:-1] - w_middle[1:-1, -2]) +
        psi[1:-1, -2] * (w_middle[2:, -1] - w_middle[:-2, -1])
    ) / 4
    f_w_middle[0, 0] += Re * (
        psi[1, 0] * (w_middle[0, 1] - w_bottom[0]) -
        psi[0, 1] * (w_middle[1, 0] - w_left[0])
    ) / 4
    f_w_middle[-1, 0] -= Re * (
        psi[-2, 0] * (w_middle[-1, 1] - w_bottom[-1]) +
        psi[-1, 1] * (w_right[0] - w_middle[-2, 0])
    ) / 4
    f_w_middle[0, -1] += Re * (
        psi[1, -1] * (w_top[0] - w_middle[0, -2]) +
        psi[0, -2] * (w_middle[1, -1] - w_left[-1])
    ) / 4
    f_w_middle[-1, -1] -= Re * (
        psi[-2, -1] * (w_top[-1] - w_middle[-1, -2]) -
        psi[-1, -2] * (w_right[-1] - w_middle[-2, -1])
    ) / 4

    return np.concatenate([
        f_poisson[:, 0], f_w_left, f_w_right, f_w_bottom, f_w_top, f_w_middle.flatten()
    ], axis=0).reshape(-1, 1)


def reconstruct_w(w_tmp, N):
    w = np.zeros((N+1, N+1))
    
    w[:1, 1:-1] = w_tmp[0*(N-1):1*(N-1)].T
    w[-1:, 1:-1] = w_tmp[1*(N-1):2*(N-1)].T
    w[1:-1, :1] = w_tmp[2*(N-1):3*(N-1)]
    w[1:-1, -1:] = w_tmp[3*(N-1):4*(N-1)]
    w[1:-1, 1:-1] = w_tmp[4*(N-1):].reshape((N - 1, N - 1))
    
    return w


def newton_iterator(
    f, get_jacobian, N, Re,
    algorithm="base", TOL=1e-8, max_iter=10, quiet=True
):
    '''
        - f: evaluates the function given x, Re
        - get_jacobian: evaluates the Jacobian given N, h
        - N: number of grid points
        - h: grid size
        - Re: Reynolds number
    '''

    h = 1 / N

    n_iter = 0 # number of iterations

    # Initialization
    # Size (N - 1) ** 2         + (N + 1) ** 2    - 4
    # Size (for streamfunction) + (for vorticity) - (corners of vorticity)
    x = np.zeros(((N - 1) ** 2 + (N + 1) ** 2 - 4, 1))
    A = poisson_solvers.construct_laplacian_kernel_matrix(N=N-1, h=1)
    f_current = f(x=x, Re=Re, kernel_matrix=A)
    
    # Check if the initial guess is a solution
    f_norm = scipy.linalg.norm(f_current)
    if f_norm <= TOL:
        if not quiet:
            print(f"n_iter={n_iter}")

        return x, n_iter
    
    while n_iter < max_iter:
        n_iter += 1
        jacobian = get_jacobian(f=f, x=x, Re=Re, kernel_matrix=A)

        kwargs = {}
        if algorithm == "bicgstab": kwargs["tol"] = 1e-6
        dx = core.solve_sparse_linear_system(
            A=jacobian, b=-h * f_current, algorithm=algorithm, **kwargs
        ).reshape((-1, 1))
        x_next = x + dx
        
        f_current = f(x=x_next, Re=Re, kernel_matrix=A)
        
        f_norm = scipy.linalg.norm(f_current)
        if not quiet:
            print(f"iter={n_iter}; residual={f_norm}; dx={scipy.linalg.norm(dx)}")
        if f_norm <= TOL:
            break
        
        x = x_next
        
    if not quiet:
        print(f"n_iter={n_iter}")
    
    return x_next, n_iter


# %% ../nbs/03_streamfunction_vorticity_newton.ipynb 10
def newton_solver(
    f, get_jacobian, N, Re,
    algorithm="base", TOL=1e-8, max_iter=10, quiet=True
):
    assert algorithm != "cg", (
        "Linear system is not positive definite, cg algorithm cannot be used"
    )

    # Initialization
    # Size (N - 1) ** 2 (for streamfunction) + (N + 1) ** 2 (for vorticity) - 4 (corners of vorticity)
    x0 = np.zeros(((N - 1) ** 2 + (N + 1) ** 2 - 4, 1))

    solution, n_iter = newton_iterator(
        f=f, get_jacobian=get_jacobian, N=N, Re=Re,
        algorithm=algorithm,
        TOL=TOL, max_iter=max_iter, quiet=quiet
    )
    
    psi, w = solution[:(N - 1) ** 2], solution[(N - 1) ** 2:]
    
    # Get final psi
    psi = psi.reshape(N - 1, N - 1)
    psi = np.pad(psi, (1, 1), mode="constant", constant_values=0)
    
    # Get final w
    w = reconstruct_w(w_tmp=w, N=N)
    w = w.reshape(N + 1, N + 1)
    
    return w, psi, n_iter

