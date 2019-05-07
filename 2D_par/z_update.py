# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

import numpy as np

# >>> flux subroutine >>>
def FLUX(glob_Mr1, glob_Mr2, loc_Mz1, loc_Mz2, Me, nWRs, U, time, r, z, D):
    # initialize arrays locally for each worker
    Fr = np.zeros((glob_Mr2, loc_Mz2), dtype = np.float64)
    Fz = np.zeros((glob_Mr2, loc_Mz2), dtype = np.float64)
    # boundary conditions (only applicable to first and last worker)
    if (Me == 1):
        for i in range(1, glob_Mr1):
            U[i][0] = 0.0
    if (Me == nWRs):
        for i in range(1, glob_Mr1):
            U[i][loc_Mz1] = 0.0
    # boundary conditions (applicable to each worker)
    for j in range(1, loc_Mz1):
        U[0][j] = 0.0
        U[glob_Mr1][j] = np.exp(-time) * np.log(2.0) * np.sin(z[j])
    # radial fluxes
    for i in range(1, glob_Mr2):
        for j in range(1, loc_Mz1):
            Fr[i][j] = -D * (U[i][j] - U[i-1][j]) / (r[i] - r[i-1])
    # axial fluxes
    for i in range(1, glob_Mr1):
        for j in range(1, loc_Mz2):
            Fz[i][j] = -D * (U[i][j] - U[i][j-1]) / (z[j] - z[j-1])
    return(U, Fr, Fz)
# <<< end flux subroutine <<<

# >>> PDE subroutine >>>
def PDE(glob_Mr1, loc_Mz1, Ar, Fr, Az, Fz, U, dt, dV):
    for i in range(1, glob_Mr1):
        for j in range(1, loc_Mz1):
            # radial fluxes and areas
            radial = Ar[i][j] * Fr[i][j] - Ar[i+1][j] * Fr[i+1][j]
            # axial fluxes and areas
            axial = Az[i][j] * Fz[i][j] - Az[i][j+1] * Fz[i][j+1]
            # update (internal) temperatures
            U[i][j] = U[i][j] + (dt / dV[i][j]) * (radial + axial)
    return(U)
# <<< end PDE subroutine <<<
