# AUTOGENERATED! DO NOT EDIT! File to edit: ../nbs/02_streamfunction_vorticity_iterative_solver/02_streamfunction_vorticity_iterative_solver.ipynb.

# %% auto 0
__all__ = ['streamfunction_vorticity_iterative_solver']

# %% ../nbs/02_streamfunction_vorticity_iterative_solver/02_streamfunction_vorticity_iterative_solver.ipynb 6
from . import core, poisson_solvers

import numpy as np

# %% ../nbs/02_streamfunction_vorticity_iterative_solver/02_streamfunction_vorticity_iterative_solver.ipynb 9
def update_vorticity_bcs(w, psi, U_wall_top, nx, ny, h):
    # First calculate BCs to 1st order

    # y = 0: U_wall = 0
    w[1:nx - 1, 0] = -(
        (psi[1:nx - 1, 1] - psi[1:nx - 1, 0]) / h - 0
    ) / (0.5 * h)
    # y = 1: U_wall_top is given here
    w[1:nx - 1, ny - 1] = -(
        U_wall_top -
        (psi[1:nx - 1, ny - 1] - psi[1:nx - 1, ny - 2]) / h
    ) / (0.5 * h)
    
    # x = 0: U_wall = 0
    w[0, 1:ny - 1] = (
        -(psi[1, 1:ny - 1] - psi[0, 1:ny - 1]) / h - 0
    ) / (0.5 * h)
    # x = 1: U_wall = 0
    w[nx - 1, 1:ny - 1] = (
        0 - (-(psi[nx - 1, 1:ny - 1] - psi[nx - 2, 1:ny - 1]) / h)
    ) / (0.5 * h)

    
    # Compute the 2nd order correction to the BCs
    w[1:nx - 1, 0] = (4 * w[1:nx - 1, 0] - w[1:nx - 1, 1]) / 3 # y = 0
    w[1:nx - 1, ny - 1] = (
        4 * w[1:nx - 1, ny - 1] - w[1:nx - 1, ny - 2]
    ) / 3 # y = 1
    
    w[0, 1:ny - 1] = (4 * w[0, 1:ny - 1] - w[1, 1:ny - 1]) / 3 # x = 0
    w[nx - 1, 1:ny - 1] = (
        4 * w[nx - 1, 1:ny - 1] - w[nx - 2, 1:ny - 1]
    ) / 3 # x = 1

    return w


def get_dw_dt(dw_dt, w, psi, Re, nx, ny, h):
    dw_dt[1:nx - 1, 1:ny - 1] = (
        -(
            (psi[1:nx - 1, 2:ny] - psi[1:nx - 1, 0:ny - 2]) *
            (w[2:nx, 1:ny - 1] - w[0:nx - 2, 1:ny - 1])
        ) / (4 * h ** 2) + 
        (
            (psi[2:nx, 1:ny - 1] - psi[0:nx - 2, 1:ny - 1]) *
            (w[1:nx - 1, 2:ny] - w[1:nx - 1, 0:ny - 2])
        ) / (4 * h ** 2) + 
        (
            w[2:nx, 1:ny - 1] + w[0:nx - 2, 1:ny - 1] +
            w[1:nx - 1, 2:ny] + w[1:nx - 1, 0:ny - 2] -
            4 * w[1:nx - 1 , 1:ny - 1]
        ) / (Re * h ** 2)
    )

    return dw_dt


# %% ../nbs/02_streamfunction_vorticity_iterative_solver/02_streamfunction_vorticity_iterative_solver.ipynb 10
# Incorporates both the Poisson problem and time-stepping for the vorticity
def streamfunction_vorticity_iterative_solver(
    N, Re, tfinal, U_wall_top, dt=None,
    algorithm="base", print_every=0.0
):
    nx = ny = N + 1 # i.e. N = 10
    h = 1 / N # h = dx (=dy) = 1 / N; x-size is 1, y-size is ny/nx
    
    w = np.zeros(shape=(nx, ny)) # vorticity
    psi = np.zeros(shape=(nx, ny)) # streamfunction
    dw_dt = np.zeros(shape=(nx, ny)) # dw/dt    
    
    if dt is None:
        dt = np.round(0.2 * Re * h ** 2, 8) # 0.8 * marginal value
    
    omega_mid = []
    
    t = 0
    while t < tfinal:
        # Start by solving the Poisson problem
        psi = poisson_solvers.poisson_non_iterative_solver(
            w, algorithm=algorithm
        )

        # Now time-step vorticity
        w = update_vorticity_bcs(w, psi, U_wall_top, nx, ny, h)

        # We can now find dw/dt
        dw_dt = get_dw_dt(dw_dt, w, psi, Re, nx, ny, h)

        # Finally, update w
        w[1:nx-1, 1:ny-1] = w[1:nx-1, 1:ny-1] + dw_dt[1:nx-1, 1:ny-1] * dt
        
        # One pass done, increment time
        t = np.round(t + dt, 8)
        
        # Print
        if print_every > 0:
            if 0 <= np.round(t % print_every, 8) < dt:
                print(f"t={t:.5f}; w(0.5, 0.5)={w[nx // 2, ny // 2]}")
        
        omega_mid.append(w[nx // 2, ny // 2])
    
    return w, psi, np.array(omega_mid)

