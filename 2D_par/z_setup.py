# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

import numpy as np

# >>> initialization subroutine >>>
def INIT(Mr2, Mz2, r, z):
    U = np.zeros((Mr2, Mz2), dtype = np.float64)
    # impose initial condition, u(r,z,0)
    for i in range(0, Mr2):
        for j in range(0, Mz2):
            # note: natural log, not log base 10
            U[i][j] = np.log(r[i]) * np.sin(z[j])
    return(U)
# <<< end initialization subroutine <<<

# >>> global mesh subroutine >>>
def glob_MESH(glob_Mr1, glob_Mr2, glob_Mz1, glob_Mz2, Rin, dr, Rout, dz, Z):
    glob_r = np.zeros(glob_Mr2, dtype = np.float64)
    glob_z = np.zeros(glob_Mz2, dtype = np.float64)
    # mesh in r-direction
    glob_r[0] = Rin
    glob_r[1] = Rin + (dr / 2)
    for i in range(2, glob_Mr1):
        glob_r[i] = glob_r[1] + (i - 1) * dr
    glob_r[glob_Mr1] = Rout
    # mesh in z-direction
    glob_z[0] = 0.0
    glob_z[1] = 0.0 + (dz / 2)
    for j in range(2, glob_Mz1):
        glob_z[j] = glob_z[1] + (j - 1) * dz
    glob_z[glob_Mz1] = Z
    return(glob_r, glob_z)
# <<< end global mesh subroutine <<<

# >>> mesh subroutine >>>
def MESH(glob_Mr1, glob_Mr2, loc_Mz, loc_Mz1, loc_Mz2, Rin, Rout, dr, dz, Me, nWRs, Z):
    # 2D mesh, control volumes, radial areas, and axial areas
    # initialize arrays locally for each worker
    r = np.zeros(glob_Mr2, dtype = np.float64)
    z = np.zeros(loc_Mz2, dtype = np.float64)
    Ar = np.zeros((glob_Mr2, loc_Mz2), dtype = np.float64)
    Az = np.zeros((glob_Mr2, loc_Mz2), dtype = np.float64)
    dV = np.zeros((glob_Mr2, loc_Mz2), dtype = np.float64)
    # mesh in r-direction
    r[0] = Rin
    r[1] = Rin + (dr / 2)
    for i in range(2, glob_Mr1):
        r[i] = r[1] + (i - 1) * dr
    r[glob_Mr1] = Rout
    # mesh in z-direction (changes depending on worker)
    loc_z0 = loc_Mz * (Me - 1) * dz - (dz / 2)
    z[0] = loc_z0
    z[1] = loc_z0 + dz
    if (Me == 1):
        z[0] = 0.0
        z[1] = 0.0 + (dz / 2)
    for j in range(2, loc_Mz2):
        z[j] = z[j-1] + dz
    if (Me == nWRs):
        z[loc_Mz1] = Z
    # areas of radial faces
    for i in range(1, glob_Mr2):
        for j in range(1, loc_Mz1):
            Ar[i][j] = 2 * np.pi * dz * (r[i] - dr / 2)
    # areas of axial faces
    for i in range(1, glob_Mr1):
        for j in range(1, loc_Mz2):
            Az[i][j] = 2 * np.pi * dr * r[i]
    # control volumes
    for i in range(1, glob_Mr1):
        for j in range(1, loc_Mz1):
            dV[i][j] = 2 * np.pi * dr * dz * r[i]
    return(r, z, Ar, Az, dV)
# <<< end mesh subroutine <<<
