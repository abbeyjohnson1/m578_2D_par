# Abbey Johnson
# February 22, 2019

# import numpy and math modules
import numpy as np
import math

# flux subroutine
def FLUX(F, Mp1, D, U, x, b, time):
    # left boundary condition (Dirichlet)
    U[0] = 1.
    # right boundary condition (Dirichlet) (caution: time-dependent boundary condition)
    # right boundary condition, impose exact solution at x = b
    U[Mp1] = math.erfc(b/(2*math.sqrt(D*time)))
    # left boundary flux
    # F[1] = -(D*(U[1]-U[0])/(x[1]-x[0]))
    # flux at internal faces
    for ii in range(1, (Mp1+1)):
        F[ii] = -D * (U[ii]-U[ii-1]) / (x[ii]-x[ii-1])
    # right boundary flux
    # F[Mp1] = -(D*(U[Mp1]-U[M]))/(x[Mp1]-x[M])
    return (F)

# PDE subroutine
def PDE(Mp1, dt, dx, F, U, b, D, time):
    # internal concentration
    for jj in range(1, Mp1):
        U[jj] = U[jj]+(dt/dx)*(F[jj]-F[jj+1])
    # update boundary concentrations after internal concentrations have been updated
    # left boundary concentration
    # U[0] = U[1]+(F[1]/D)*(x[1]-x[0])
    # left Dirichlet condition
    U[0] = 1.
    # right boundary concentration
    # U[Mp1] = U[M]-(F[Mp1]/D)*(x[Mp1]-x[M])
    # right Dirichlet condition (caution: time-dependent) (exact solution at x = b)
    U[Mp1] = math.erfc(b/(2*math.sqrt(D*time)))
    return (U)
