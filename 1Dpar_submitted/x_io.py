# Abbey L. Johnosn

# import numpy module
import numpy as np

# input subroutine: read parameters from data file
def INPUT(filename):
    data = np.loadtxt(filename)
    return (data)

# output subroutine: print data to file
def OUTPUT(filename, U, u_exact, ERR):
    # print concentration profile
    np.savetxt(filename, np.transpose([U, u_exact, ERR]))
