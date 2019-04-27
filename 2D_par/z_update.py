# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

# import modules
import numpy as np
import math

# flux subroutine

# flux subroutine
def FLUX(r, Mr1, z, Mz1, U, time, Fr, D, Fz):
    # boundary conditions
    for i in range(1, Mr1):
        U[i][0] = 0.0
        U[i][Mz1] = 0.0
    for j in range(1, Mz1):
        U[0][j] = 0.0
        U[Mr1][j] = np.exp(-time) * np.log(2.) * np.sin(z[j])
    # radial fluxes
    for i in range(1, Mr1 + 1):
        for j in range(1, Mz1):
            Fr[i][j] = -D * (U[i][j] - U[i-1][j]) / (r[i] - r[i-1])
    # axial fluxes
    for i in range(1, Mr1):
        for j in range(1, Mz1 + 1):
            Fz[i][j] = -D * (U[i][j] - U[i][j-1]) / (z[j] - z[j-1])
    return(U, Fr, Fz)

###############################################################################

# PDE subroutine
def PDE(Mr1, Mz1, Ar, Fr, Az, Fz, U, dt, dV):
    for i in range(1, Mr1):
        for j in range(1, Mz1):
            # radial fluxes and areas
            radial = Ar[i][j] * Fr[i][j] - Ar[i+1][j] * Fr[i+1][j]
            # axial fluxes and areas
            axial = Az[i][j] * Fz[i][j] - Az[i][j+1] * Fz[i][j+1]
            # update (internal) temperatures
            U[i][j] = U[i][j] + (dt / dV[i][j]) * (radial + axial)
    return(U)


###############################################################################

# subroutine comparing approximation with exact solution
def COMPARISON(Mr, Mz, u_exact, time, r, z, U):
    # comp.write('# i j U(i,j) u_exact(i,j) error(i,j) at time = %f \n' % time)
    ERR = 0.0
    for i in range(0, Mr + 2):
        for j in range(0, Mz + 2):
            # compute exact solution
            u_exact[i][j] = np.exp(-time) * np.log(r[i]) * np.sin(z[j])
            # compute error at r[i], z[j]
            ERRij = abs(U[i][j] - u_exact[i][j])
            # print to comparison output file
            # comp.write('%i %i %e %e %e \n' % (i, j, U[i][j], u_exact[i][j], ERRij))
            # compute max error
            ERR = max(ERRij, ERR)
    # comp.write('\n')
    # comp.write('Maximum error at time %f is %e \n' % (time, ERR))
    # comp.write('\n')
    return(ERR)

###############################################################################
