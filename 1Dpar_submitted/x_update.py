# Abbey L. Johnson

# imports
import numpy as np
import math

# flux subroutine
def FLUX(nWRs, Me, F, Mp1, D, U, x):
    # flux at internal faces
    for ii in range(1, (Mp1 + 1)):
        F[ii] = -D * (U[ii] - U[ii-1]) / (x[ii] - x[ii-1])
    return(F)

# PDE subroutine
def PDE(nWRs, Me, Mp1, dt, dx, F, U, D, time, x):
    # internal concentration
    for jj in range(1, Mp1):
        U[jj] = U[jj] + (dt / dx) * (F[jj] - F[jj+1])
    # update boundary concentrations after internal concentrations updated
    # left boundary concentration (Dirichlet BC)
    if (Me == 1):
        U[0] = 1.0
    # right boundary concentration (Dirichlet BC, exact soln. at x = b)
    if (Me == nWRs):
        U[Mp1] = math.erfc(x[Mp1] / (2 * math.sqrt(D * time)))
    return(U)

def COMPARISON(x, U, Mp1, time, D, u_exact):
    # ERR = 0
    for q in range(0, (Mp1+1)):
        # compute exact solution
        u_exact[q] = math.erfc(x[q] / (2*math.sqrt(D * time)))
        # compute error at x[i]
        # ERRi = abs(U[q] - u_exact[q])
        # compute max error
        # ERR = max(ERRi, ERR)
    # return(ERR, u_exact)
    return(u_exact)
