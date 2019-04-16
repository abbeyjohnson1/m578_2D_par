# Abbey L. Johnson

# import numpy module
import numpy as np

# mesh subroutine
def MESH(x, glob_a, dx, M, Mp1, glob_b, Me, nWRs):
    # a-value local to each worker
    loc_a = M * (Me - 1) * dx - (dx / 2)
    # left boundary of mesh
    x[0] = loc_a
    x[1] = loc_a + dx
    # left boundary of mesh (different if first worker)
    # NOTE TO SELF: would maybe be more efficient if I made this an if-else loop?
    if (Me == 1):
        x[0] = glob_a
        x[1] = glob_a + (dx / 2)
    # internal nodes of mesh
    for i in range(2, Mp1+1):
        x[i] = x[i-1] + dx
    # right boundary of mesh (only necessary to impose for last worker)
    if (Me == nWRs):
        x[Mp1] = glob_b
    return(x)

# initialization subroutine
def INIT(Mp1, Me):
    # initial concentration profile/initial time condition
    U = np.zeros(Mp1 + 1)
    # note: only for worker 1
    if (Me == 1):
        U[0] = 1.0
    return(U)
