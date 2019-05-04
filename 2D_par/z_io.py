# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

import numpy as np

# >>> input subroutine >>>
def INPUT(filename):
    # read parameters from initialization data file
    with open(filename) as init_file:
        for line in init_file:
            line = line.partition('#')[0]
            input = line.split()
    init_file.close()
    return(input)
# <<< end input subroutine <<<

# >>> concentration profile output >>>
def OUTPUT(prof, time, U, glob_r, glob_z):
    # record mesh corresponding to temperature profiles
    if (time == 0.0):
        prof.write('# mesh in z-direction \n')
        np.savetxt(prof, glob_z)
        prof.write('# mesh in r-direction \n')
        np.savetxt(prof, glob_r)
    # print temperature profile
    prof.write('# temperature profile of U(r,z,t) at t = %f \n' % time)
    np.savetxt(prof, U, delimeter = '\t')
    prof.write('\n')
# <<< end concentration profile output <<<

# >>> comparison subroutine >>>
def COMPARISON(comp, time, U, glob_Mr2, glob_Mz2, glob_r, glob_z):
    # compare numerical approximation with exact solution
    comp.write('# i j U(i,j) u_exact(i,j) error(i,j) at time = %f \n' % time)
    ERR = 0.0
    u_exact = np.zeros((glob_Mr2, glob_Mz2), dtype = np.float64)
    for i in range(0, glob_Mr2):
        for j in range(0, glob_Mz2):
            # compute exact solution
            u_exact[i][j] = np.exp(-time) * np.log(glob_r[i]) * np.sin(glob_z[j])
            # compute error at r[i], z[j]
            ERRij = abs(U[i][j] - u_exact[i][j])
            # print to comparison output file
            comp.write('%i %i %e %e %e \n' % (i, j, U[i][j], u_exact[i][j], ERRij))
            # compute max error
            ERR = max(ERRij, ERR)
    comp.write('\n')
    comp.write('Maximum error at time %f is %e \n' % (time, ERR))
    comp.write('\n')
    return(ERR)
# <<< end comparison subroutine <<<

# >>> runtime output subroutine >>>
def OUTPUT_runtime(outf, outputs):
    # print summary of runtime information to output file
    time = np.float64(outputs[0])
    nsteps = np.int_(outputs[1])
    print('Master done: exiting at t = %f after %i steps. \n' % (time, nsteps), file = outf)
    ERR = np.float64(outputs[2])
    print('Maximum error is %e at time %f \n' % (ERR, time), file = outf)
    MMr = np.int_(outputs[3])
    print('Parameters for this run: ', file = outf)
    print('MMr = %f ' % MMr, file = outf)
    MMz = np.int_(outputs[4])
    print('MMz = %f ' % MMz, file = outf)
    Rin = np.float64(outputs[5])
    print('Rin = %f ' % Rin, file = outf)
    Rout = np.float64(outputs[6])
    print('Rout = %f' % Rout, file = outf)
    tend = np.float64(outputs[7])
    print('t_end = %f ' % tend, file = outf)
    factor = np.float64(outputs[8])
    print('factor = %f ' % factor, file = outf)
    dtout = np.float64(outputs[9])
    print('dtout = %f ' % dtout, file = outf)
    D = np.float64(outputs[10])
    print('D = %f \n' % D, file = outf)
    nWRs = np.int_(outputs[11])
    print('Number of workers  = %i \n' % nWRs, file = outf)
# <<< end runtime output subroutine <<<
