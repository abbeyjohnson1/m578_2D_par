# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

# import modules, MPI, and subroutines
import math
import numpy as np
from mpi4py import MPI
from z_messaging import RECV_output_MPI
from z_io import INPUT, OUTPUT
from z_setup import MESH, INIT
from z_update import COMPARISON

# >>> master >>>

# NEED ARGUMENTS
def MASTER(comm):

    # check to see if z_mainMR file is running
    print('z_mainMR is being run')

    # read parameters from initialization data file
    with open("./Init.txt") as init_file:
        for line in init_file:
            line = line.partition('#')[0]
            input = line.split()

    # assign parameter values
    # NOTE TO SELF: make sure to always use even M divisible by nWRs!!!

    # nodes in r-direction
    MMr = np.float64(input[0])
    # nodes in z-direction
    MMz = np.float64(input[1])
    # inner radius of cylinder
    Rin = np.float64(input[2])
    # outer radius of cylinder
    Rout = np.float64(input[3])
    # end time
    tend = np.float64(input[4])
    # factor for stability
    factor = np.float64(input[5])
    # time-step when to output
    dtout = np.float64(input[6])
    # diffusivity
    D = np.float64(input[7])


# should these files be opened and closed here?
    # open files for outputting info
    outf = open('o_out.txt', 'w')
    prof = open('o_prof.txt', 'w')
    comp = open('o_compare.txt', 'w')

    # compute other parameters

    # height of cylindar
    Z = np.pi

    dr = np.float64(1.0 / MMr)

    glob_Mr = int( (Rout - Rin) * MMr )
    glob_Mr1 = int(Mr + 1)
    glob_Mr2 = int(Mr + 2)

    dz = np.float64(1.0 / MMz)

    glob_Mz = int( (Z - 0) * MMz )
    glob_Mz1 = int(Mz + 1)
    glob_Mz2 = int(Mz + 2)

    dtEXPL = np.float64(1 / (2 * D * (1 / (dr * dr) + 1 / (dz * dz) ) ) )

    # time-step (for stability of explicit scheme)
    dt = np.float64(factor * dtEXPL)

    # declare variables and dimensions of arrays
# not sure that this needs to be here

    U = np.zeros((Mr + 2, Mz + 2), dtype = np.float64)
    u_exact = np.zeros((Mr + 2, Mz + 2), dtype = np.float64)

    # maximum number of iterations
    MaxSteps = int(tend / dt) + 2

    # when to print outputs
    tout = dtout

    # initialize time-stepping variables
    nsteps = 0
    time = 0.0

#??? call mesh subroutine

#??? call initialization subroutine

#??? call output routine to record concentration profile at time = 0

    # pack integers in iparms array
    iparms = np.array([ , , , ], dtype = np.int)

    # pack reals in parms array
    parms = np.array([ , , , ], dtype = np.float64)

    # send parms arrays to everyone
    comm.bcast(iparms, root = 0)
    comm.bcast(parms, root = 0)

    # time-stepping loop:
    # workers do computation and send output to master every dtout
    # master recieves from each worker its secton of output arrays(s)

    # begin time-stepping
    for nsteps in range(1, MaxSteps):

        # synchronize everyone
        comm.Barrier()

        # update time = time + dt
        time = nsteps * dt

        # check if ready to output
        if time >= tout:

# need arguments
            # recieve output from workers when time = dtout
            U = RECV_output_MPI(comm)

# global or local Mr and Mz
            # call comparison subroutine
            ERR = COMPARISON(Mr, Mz, u_exact, time, r, z, U)

# call output subroutine

            # update tout
            tout = tout + dtout

    # end time-stepping

    # print run-time information to screen
    print('MASTER DONE: exiting at t = %f after %i steps.' % (time, nsteps))
    print('MASTER DONE: maximum error = {}' .format(maxERR))

    # print information to outf file
    print('Parameters for this run: ', file = outf)
    print('MMr = %f ' % MMr, file = outf)
    print('MMz = %f ' % MMz, file = outf)
    print('Rin = %f ' % Rin, file = outf)
    print('Rout = %f' % Rout, file = outf)
    print('t_end = %f ' % tend, file = outf)
    print('factor = %f ' % factor, file = outf)
    print('dtout = %f ' % dtout, file = outf)
    print('D = %f \n' % D, file = outf)

    # print run-time information to output file
    print('MASTER DONE: exiting at t = %f after %i steps. \n' % (time, nsteps), file = outf)
#    print('Maximum error is %e at time %f \n' % (ERR, time), file = outf)

    # close files
    init_file.close()
    prof.close()
    comp.close()
    outf.close()

# <<< end master <<<
