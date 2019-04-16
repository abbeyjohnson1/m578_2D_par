# Abbey L. Johnosn

# import numpy module
import numpy as np

# input subroutine: read parameters from data file
def INPUT(filename):
    data = np.loadtxt(filename)
    return (data)

# output subroutine: print data to file
def OUTPUT(filename, x, U):
    # print concentration profile,  x(i), U(i) into output file
    np.savetxt(filename, np.transpose([x, U]))
