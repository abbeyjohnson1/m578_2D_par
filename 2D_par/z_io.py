# Abbey L. Johnson

# import modules
import numpy as np

# input subroutine: read parameters from data file
def INPUT(filename):
    data = np.loadtxt(filename)
    return(data)


################################################################################

# output subroutine: print stuff??? to file
def OUTPUT():
    1

# output subroutine
# def OUTPUT(time, r, Mr, z, Mz, U):
#    # print temperature profile to output file
#    prof.write('# temperature profile of U(r,z,t) at t = %f \n' % time)
#    prof.write('# i j r(i) z(j) U(r,z,t) \n')
#    # note: maybe should print to file in matrix format?
#    for i in range(0, Mr + 2):
#        for j in range (0, Mz + 2):
#            prof.write('%i %i %f %f %e \n' % (i, j, r[i], z[j], U[i][j]))
#    prof.write('\n')

################################################################################
