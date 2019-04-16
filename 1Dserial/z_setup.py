# Abbey Johnson
# Febrary 22, 2019

# import numpy module
import numpy as np

# mesh subroutine
# mesh of M uniform control volumes in the interval [a, b] up to time t_end
def MESH(x, a, dx, Mp1, b):
    # left boundary of mesh
    x[0] = a
    # internal nodes of mesh
    x[1] = a + dx/2
    for i in range(2, Mp1):
        x[i] = x[1]+(i-1)*dx
    # right boundary of mesh
    x[Mp1] = b
    return (x)

# initialization subroutine (for the initial concentration profile / initial condition)
def INIT(U, Mp1):
    for j in range(0, (Mp1+1)):
        U[j] = 0.0
    return (U)
