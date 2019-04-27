# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

# import modules
import numpy as np



################################################################################
# do additional output, one output subroutine per file profile, error, run info etc

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
# output to outf
    # print parameter values etc. to output file
    print('Parameters for this run: ', file = outf)
    print('MMr = %f ' % MMr, file = outf)
    print('MMz = %f ' % MMz, file = outf)
    print('Rin = %f ' % Rin, file = outf)
    print('Rout = %f' % Rout, file = outf)
    print('t_end = %f ' % tend, file = outf)
    print('factor = %f ' % factor, file = outf)
    print('dtout = %f ' % dtout, file = outf)
    print('D = %f \n' % D, file = outf)
