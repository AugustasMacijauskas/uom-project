# AUTOGENERATED! DO NOT EDIT! File to edit: ../nbs/05_streamfunction_vorticity_pytorch.ipynb.

# %% auto 0
__all__ = ['newton_solver_pytorch']

# %% ../nbs/05_streamfunction_vorticity_pytorch.ipynb 4
from . import core

import numpy as np

import torch
import torch.nn.functional as F

# %% ../nbs/05_streamfunction_vorticity_pytorch.ipynb 8
# Given the vorticity, solve the Poisson eqn. to find the streamfunction
def f_pytorch(x, Re, U_wall_top):
    N = int(np.sqrt(x.shape[0] // 2 + 1))
    h = 1 / N

    psi = x[:(N-1)**2].reshape(N-1, N-1)
    w_left   = x[(N-1)**2 + 0*(N-1) : (N-1)**2 + 1*(N-1)][:, 0]
    w_right  = x[(N-1)**2 + 1*(N-1) : (N-1)**2 + 2*(N-1)][:, 0]
    w_bottom = x[(N-1)**2 + 2*(N-1) : (N-1)**2 + 3*(N-1)][:, 0]
    w_top    = x[(N-1)**2 + 3*(N-1) : (N-1)**2 + 4*(N-1)][:, 0]
    w_middle = x[(N-1)**2 + 4*(N-1) :].reshape(N-1, N-1)

    # Calculate the equations coming from the Poisson equation
    f_poisson = -4 * psi + h ** 2 * w_middle
    f_poisson[:-1, :] += psi[1:, :]
    f_poisson[1:, :] += psi[:-1, :]
    f_poisson[:, :-1] += psi[:, 1:]
    f_poisson[:, 1:] += psi[:, :-1]

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

    return torch.hstack([
        f_poisson.flatten(), f_w_left, f_w_right, f_w_bottom, f_w_top,
        f_w_middle.flatten()
    ])[:, None]


def reconstruct_w_pytorch(w_tmp, N, device):
    w = torch.zeros((N+1, N+1), device=device)

    w[:1, 1:-1] = w_tmp[0*(N-1):1*(N-1)].T
    w[-1:, 1:-1] = w_tmp[1*(N-1):2*(N-1)].T
    w[1:-1, :1] = w_tmp[2*(N-1):3*(N-1)]
    w[1:-1, -1:] = w_tmp[3*(N-1):4*(N-1)]
    w[1:-1, 1:-1] = w_tmp[4*(N-1):].reshape((N - 1, N - 1))

    return w


def newton_iterator_pytorch(
    f, get_jacobian, N, Re, U_wall_top, device,
    algorithm="base", TOL=1e-8, max_iter=10, quiet=True
):
    '''
        - f: evaluates the function given x, Re
        - get_jacobian: evaluates the Jacobian given N, h
        - N: number of grid points
        - h: grid size
        - Re: Reynolds number
    '''

    n_iter = 0 # number of iterations

    # Initialization
    # Size (N - 1) ** 2         + (N + 1) ** 2    - 4
    # Size (for streamfunction) + (for vorticity) - (corners of vorticity)
    x = torch.zeros(((N - 1) ** 2 + (N + 1) ** 2 - 4, 1), device=device).double()
    f_current = f(x=x, Re=Re, U_wall_top=U_wall_top)

    # Check if the initial guess is a solution
    f_norm = torch.linalg.norm(f_current)
    if f_norm <= TOL:
        if not quiet:
            print(f"n_iter={n_iter}")

        return x, n_iter

    while n_iter < max_iter:
        n_iter += 1
        jacobian = get_jacobian(x, Re, U_wall_top).squeeze()

        kwargs = {}
        if algorithm == "lstsq" and device.type == "cpu":
            kwargs["driver"] = "gels"
        dx = core.solve_sparse_linear_system_pytorch(
            A=jacobian, b=-f_current, algorithm=algorithm, **kwargs
        )
        x_next = x + dx

        f_current = f(x=x_next, Re=Re, U_wall_top=U_wall_top)

        f_norm = torch.linalg.norm(f_current)
        if not quiet:
            print(f"iter={n_iter}; residual={f_norm}; dx={torch.linalg.norm(dx)}")
        if f_norm <= TOL:
            break

        x = x_next

    if not quiet:
        print(f"n_iter={n_iter}")

    return x_next, n_iter


# %% ../nbs/05_streamfunction_vorticity_pytorch.ipynb 9
def newton_solver_pytorch(
    f, get_jacobian, N, Re, U_wall_top, device,
    algorithm="base", TOL=1e-8, max_iter=10, quiet=True
):

    solution, n_iter = newton_iterator_pytorch(
        f=f, get_jacobian=get_jacobian, N=N, Re=Re, U_wall_top=U_wall_top,
        algorithm=algorithm, device=device,
        TOL=TOL, max_iter=max_iter, quiet=quiet
    )

    psi, w = solution[:(N - 1) ** 2], solution[(N - 1) ** 2:]

    # Get final psi
    psi = psi.reshape(N - 1, N - 1)
    psi = F.pad(psi, (1, 1, 1, 1), mode="constant", value=0)

    # Get final w
    w = reconstruct_w_pytorch(w_tmp=w, N=N, device=device)
    w = w.reshape(N + 1, N + 1)

    return w, psi, n_iter

