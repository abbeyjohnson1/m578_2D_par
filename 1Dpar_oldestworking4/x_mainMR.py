# Abbey L. Johnson

# necessary imports
import math
import numpy as np
from mpi4py import MPI
from x_messaging import RECV_output_MPI
from x_io import INPUT, OUTPUT
from x_setup import MESH, INIT
from x_update import COMPARISON

def MASTER(nWRs, mster, comm, myID):

    # check to see if mainMR.py file is running
    # print("z_mainMR is being run")

    # input parameter valus from file
    MM, tend, factor, dtout, D, glob_a, glob_b = INPUT('./Init.txt')

    # print parameter values to parms_out.txt
    with open('parms_out.txt', 'w') as outf:
        print('Parameters for this run: ', file = outf)
        print('MM = %f ' % MM, file = outf)
        print('t_end = %f ' % tend, file = outf)
        print('factor = %f ' % factor, file = outf)
        print('dtout = %f ' % dtout, file = outf)
        print('D = %f ' % D, file = outf)
        print('a = %f ' % glob_a, file = outf)
        print('b = %f ' % glob_b, file = outf)

# not sure about this output file part???
    # create output file for concentration profile
    filename = 'o_prof.txt'

    # compute parameters

    # delta x
    dx = np.float64(1.0 / MM)

    # M = int((b-a)/dx)
    glob_M = int((glob_b - glob_a) * MM)

    Mp1 = int(glob_M + 1)

    dtEXPL = np.float64((dx * dx)/(2 * D))

    # time-step (for stability of the explicit scheme)
    dt = np.float64(factor * dtEXPL)

    # maximum number of iterations
    MaxSteps = int(tend / dt) + 1

    # when to print outputs
    tout = dtout

    # initialize time-stepping variables
    nsteps = 0
    time = 0.0

# this part may not be necessary???
    # call initialization subroutine
    U = np.zeros((glob_M + 2), dtype = np.float64)
    #U = INIT(Mp1)

# don't think this part is necessary???
    # initialize other arrays
    #F = np.zeros((glob_M + 2), dtype = np.float64)
    #u_exact = np.zeros((glob_M + 2), dtype = np.float64)

# not sure how to record output at time 0???
    # call output routine to record the concentration profile at time t = 0
    # OUTPUT(filename, x, U)

    # pack integers in iparms array
    iparms = np.array([glob_M, Mp1, MaxSteps], dtype=np.int)

    # pack reals in parms array
    parms = np.array([dx, dtEXPL, dt, tend, factor, dtout, D, glob_a, glob_b, tout], dtype = np.float64)

    # send iparms to everyone
    comm.bcast(iparms, root = 0)

    # send parms to everyone
    comm.bcast(parms, root = 0)

    # workers do computation and send output to master every dtout
    # master recieves from each worker its secton of output arrays(s)

# LATER: insert STS scheme to improve speedup!?!?!

# not sure this is necessary????
    comm.Barrier()

    # begin time-stepping
    for nsteps in range(1, MaxSteps + 1):

        # synchronize everyone
        comm.Barrier()

        # update time = time + dt
        time = nsteps * dt

        if time >= tout:

# not sure if this should be global M or local M in the argument
            # recieve output from workers when time=dtout
            RECV_output_MPI(nWRs, glob_M, U, comm)

# not sure about this part???
            # call compare subroutine
            # ERR = COMPARISON(x, U, Mp1, time, D, u_exact)

# not sure abou this part????
            # call output subroutine
            # OUTPUT(filename, x, U)

# check to see if i am importing too many subroutines???

            # update tout
            tout = tout + dtout

    # end time-stepping

    # print run-time information to screen
    # print(' MASTER DONE: exiting at t = %f after %i steps.' % (time, nsteps))
