# AUTOGENERATED! DO NOT EDIT! File to edit: ../nbs/06_streamfunction_vorticity_jax.ipynb.

# %% auto 0
__all__ = ['newton_solver_jax']

# %% ../nbs/06_streamfunction_vorticity_jax.ipynb 4
from . import core

from functools import partial

import numpy as np

import jax
import jax.numpy as jnp
jax.config.update("jax_enable_x64", True) # enable JAX to use double precision


# %% ../nbs/06_streamfunction_vorticity_jax.ipynb 7
@jax.jit
def f_jax(x, Re, U_wall_top):
    N = int(np.sqrt(x.shape[0] // 2 + 1))
    h = 1 / N

    psi = x[:(N-1)**2].reshape(N-1, N-1)
    w_left   = x[(N-1)**2 + 0*(N-1) : (N-1)**2 + 1*(N-1)]
    w_right  = x[(N-1)**2 + 1*(N-1) : (N-1)**2 + 2*(N-1)]
    w_bottom = x[(N-1)**2 + 2*(N-1) : (N-1)**2 + 3*(N-1)]
    w_top    = x[(N-1)**2 + 3*(N-1) : (N-1)**2 + 4*(N-1)]
    w_middle = x[(N-1)**2 + 4*(N-1) :].reshape(N-1, N-1)

    # Calculate the equations coming from the Poisson equation
    f_poisson = -4 * psi
    f_poisson = f_poisson.at[:-1, :].add(psi[1:, :])
    f_poisson = f_poisson.at[1:, :].add(psi[:-1, :])
    f_poisson = f_poisson.at[:, :-1].add(psi[:, 1:])
    f_poisson = f_poisson.at[:, 1:].add(psi[:, :-1])
    f_poisson = f_poisson + h ** 2 * w_middle

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
    f_w_middle = f_w_middle.at[:-1, :].add(w_middle[1:, :])
    f_w_middle = f_w_middle.at[-1:, :].add(w_right)
    f_w_middle = f_w_middle.at[1:, :].add(w_middle[:-1, :])
    f_w_middle = f_w_middle.at[:1, :].add(w_left)
    f_w_middle = f_w_middle.at[:, :-1].add(w_middle[:, 1:])
    f_w_middle = f_w_middle.at[:, -1].add(w_top)
    f_w_middle = f_w_middle.at[:, 1:].add(w_middle[:, :-1])
    f_w_middle = f_w_middle.at[:, 0].add(w_bottom)

    f_w_middle = f_w_middle.at[1:-1, 1:-1].add(Re * (
        (psi[2:, 1:-1] - psi[:-2, 1:-1]) * (w_middle[1:-1, 2:] - w_middle[1:-1, :-2]) -
        (psi[1:-1, 2:] - psi[1:-1, :-2]) * (w_middle[2:, 1:-1] - w_middle[:-2, 1:-1])
    ) / 4)
    f_w_middle = f_w_middle.at[:1, 1:-1].add(Re * (
        psi[1, 1:-1] * (w_middle[0, 2:] - w_middle[0, :-2]) -
        (psi[0, 2:] - psi[0, :-2]) * (w_middle[1, 1:-1] - w_left[1:-1])
    ) / 4)
    f_w_middle = f_w_middle.at[-1:, 1:-1].add(-Re * (
        psi[-2, 1:-1] * (w_middle[-1, 2:] - w_middle[-1, :-2]) +
        (psi[-1, 2:] - psi[-1, :-2]) * (w_right[1:-1] - w_middle[-2, 1:-1])
    ) / 4)
    f_w_middle = f_w_middle.at[1:-1, 0].add(Re * (
        (psi[2:, 0] - psi[:-2, 0]) * (w_middle[1:-1, 1] - w_bottom[1:-1]) -
        psi[1:-1, 1] * (w_middle[2:, 0] - w_middle[:-2, 0])
    ) / 4)
    f_w_middle = f_w_middle.at[1:-1, -1].add(Re * (
        (psi[2:, -1] - psi[:-2, -1]) * (w_top[1:-1] - w_middle[1:-1, -2]) +
        psi[1:-1, -2] * (w_middle[2:, -1] - w_middle[:-2, -1])
    ) / 4)
    f_w_middle = f_w_middle.at[0, 0].add(Re * (
        psi[1, 0] * (w_middle[0, 1] - w_bottom[0]) -
        psi[0, 1] * (w_middle[1, 0] - w_left[0])
    ) / 4)
    f_w_middle = f_w_middle.at[-1, 0].add(-Re * (
        psi[-2, 0] * (w_middle[-1, 1] - w_bottom[-1]) +
        psi[-1, 1] * (w_right[0] - w_middle[-2, 0])
    ) / 4)
    f_w_middle = f_w_middle.at[0, -1].add(Re * (
        psi[1, -1] * (w_top[0] - w_middle[0, -2]) +
        psi[0, -2] * (w_middle[1, -1] - w_left[-1])
    ) / 4)
    f_w_middle = f_w_middle.at[-1, -1].add(-Re * (
        psi[-2, -1] * (w_top[-1] - w_middle[-1, -2]) -
        psi[-1, -2] * (w_right[-1] - w_middle[-2, -1])
    ) / 4)

    return jnp.concatenate([
        f_poisson.flatten(), f_w_left, f_w_right, f_w_bottom, f_w_top,
        f_w_middle.flatten(),
    ], axis=0)


# %% ../nbs/06_streamfunction_vorticity_jax.ipynb 16
def reconstruct_w_jax(w_tmp, N):
    w = jnp.zeros((N+1, N+1))

    w = w.at[0, 1:-1].set(w_tmp[0*(N-1):1*(N-1)])
    w = w.at[-1, 1:-1].set(w_tmp[1*(N-1):2*(N-1)])
    w = w.at[1:-1, 0].set(w_tmp[2*(N-1):3*(N-1)])
    w = w.at[1:-1, -1].set(w_tmp[3*(N-1):4*(N-1)])
    w = w.at[1:-1, 1:-1].set(w_tmp[4*(N-1):].reshape((N - 1, N - 1)))

    return w


@partial(jax.jit, static_argnums=(0, 5,))
def newton_iteration_step(
    f, x, Re, U_wall_top, f_current, algorithm="lineax_lu",
):
    jacobian = jax.jacfwd(f)(x, Re, U_wall_top)

    dx = core.solve_sparse_linear_system_jax(
        A=jacobian, b=-f_current, algorithm=algorithm,
    )
    x = x + dx

    f_current = f(x=x, Re=Re, U_wall_top=U_wall_top)
    f_norm = jnp.linalg.norm(f_current)

    return x, dx, f_current, f_norm


def newton_iterator_jax(
    f, N, Re, U_wall_top,
    algorithm="lineax_lu", TOL=1e-8, max_iter=10, quiet=True
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
    x = jnp.zeros(((N - 1) ** 2 + (N + 1) ** 2 - 4, ), dtype=jnp.float64)
    f_current = f(x=x, Re=Re, U_wall_top=U_wall_top)

    # Check if the initial guess is a solution
    f_norm = jnp.linalg.norm(f_current)
    if f_norm <= TOL:
        if not quiet:
            print(f"n_iter={n_iter}")

        return x, n_iter

    # Iterate
    while n_iter < max_iter:
        n_iter += 1

        x, dx, f_current, f_norm = newton_iteration_step(
            f=f, x=x, Re=Re, U_wall_top=U_wall_top,
            f_current=f_current,
            algorithm=algorithm,
        )

        # Check for convergence
        if not quiet:
            print(f"iter={n_iter}; residual={f_norm}; dx={jnp.linalg.norm(dx)}")
        if f_norm <= TOL:
            break

    if not quiet:
        print(f"n_iter={n_iter}")

    return x, n_iter


# %% ../nbs/06_streamfunction_vorticity_jax.ipynb 17
def newton_solver_jax(
    f, N, Re, U_wall_top,
    algorithm="lineax_lu", TOL=1e-8, max_iter=10, quiet=True
):

    solution, n_iter = newton_iterator_jax(
        f=f, N=N, Re=Re, U_wall_top=U_wall_top,
        algorithm=algorithm,
        TOL=TOL, max_iter=max_iter, quiet=quiet
    )

    psi, w = solution[:(N - 1) ** 2], solution[(N - 1) ** 2:]

    # Get final psi
    psi = psi.reshape(N - 1, N - 1)
    psi = jnp.pad(psi, (1, 1), mode="constant", constant_values=0)

    # Get final w
    w = reconstruct_w_jax(w_tmp=w, N=N)
    w = w.reshape(N + 1, N + 1)

    return w, psi, n_iter

