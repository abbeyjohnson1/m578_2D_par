# Abbey L. Johnosn
# February 22, 2019

# import numpy module
import numpy as np

# import math module
import math

# input subroutine: read from data file: MM, tend, factor, dtout, D, a, b
def INPUT(filename):
    data = np.loadtxt(filename)
    return (data)

# output subroutine: print data to file
def OUTPUT(filename, x, U):
    # print concentration profile,  x(i), U(i) into output file
    np.savetxt(filename, np.transpose([x, U]))

def COMPARISON(x, U, Mp1, time, D, u_exact):
    ERR = 0
    for q in range(0, (Mp1+1)):
        # compute exact solution
        u_exact[q] = math.erfc(x[q]/(2*math.sqrt(D*time)))
        # compute error at x[i]
        ERRi = abs(U[q]-u_exact[q])
        # compute max error
        ERR = max(ERRi, ERR)
    return(ERR)
