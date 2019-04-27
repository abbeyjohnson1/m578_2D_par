# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

# import modules
import numpy as np

# mesh subroutine

# 2D mesh, control volumes, radial areas, and axial areas
def MESH(r, Rin, dr, Mr1, Rout, z, dz, Mz1, Z, Mr, Ar, Mz, Az, dV):
    # mesh in r-direction
    # left boundary
    r[0] = Rin
    # internal nodes
    r[1] = Rin + dr / 2
    for i in range(2, Mr1):
        r[i] = r[1] + (i - 1) * dr
    # right boundary
    r[Mr1] = Rout
    # mesh in z-direction
    # left boundary
    z[0] = 0.0
    # internal nodes
    z[1] = 0.0 + dz / 2
    for j in range(2, Mz1):
        z[j] = z[1] + (j - 1) * dz
    # right boundary
    z[Mz1] = Z
    # areas of radial faces
    for i in range(1, Mr + 2):
        for j in range(1, Mz1):
            Ar[i][j] = 2 * np.pi * dz * (r[i] - dr / 2)
    # areas of axial faces
    for i in range(1, Mr1):
        for j in range(1, Mz + 2):
            Az[i][j] = 2 * np.pi * dr * r[i]
    # control volumes
    for i in range(1, Mr1):
        for j in range(1, Mz1):
            dV[i][j] = 2 * np.pi * dr * dz * r[i]
    return(r, z, Ar, Az, dV)

###############################################################################

# initialization subroutine

# initialization subroutine (for initial temperature profile/initial condition)
def INIT(Mr, Mz, U, r, z):
    # impose initial condition, u(r,z,0)
    for i in range(0, Mr + 2):
        for j in range(0, Mz + 2):
            # note: natural log, not log base 10
            U[i][j] = np.log(r[i]) * np.sin(z[j])
    return(U)

###############################################################################
