# Abbey L. Johnson

# necessary imports
import math
import numpy as np
from mpi4py import MPI
from x_messaging import RECV_output_MPI, RECV_output_exact_MPI
from x_io import INPUT, OUTPUT
from x_setup import MESH, INIT
from x_update import COMPARISON

def MASTER(nWRs, mster, comm, myID):

    # check to see if mainMR.py file is running
    # print("z_mainMR is being run")

    # input parameter valus from file
    MM, tend, factor, dtout, D, glob_a, glob_b = INPUT('./Init.txt')

    # output file
    outf = open('o_out.txt', 'w')

    # print parameter values etc to output file
    print('Parameters for this run: ', file = outf)
    print('MM = %f ' % MM, file = outf)
    print('t_end = %f ' % tend, file = outf)
    print('factor = %f ' % factor, file = outf)
    print('dtout = %f ' % dtout, file = outf)
    print('D = %f ' % D, file = outf)
    print('a = %f ' % glob_a, file = outf)
    print('b = %f ' % glob_b, file = outf)
    print(' ', file = outf)
    print('Number of workers = %i' % nWRs, file = outf)
    print(' ', file = outf)

    # compute parameters

    # delta x
    dx = np.float64(1.0 / MM)

    # make sure to always use an even M divisible by nWRs!!!
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

    # maybe call output routine to record the concentration profile at time t = 0

    # pack integers in iparms array
    iparms = np.array([glob_M, MaxSteps, nsteps], dtype=np.int)

    # pack reals in parms array
    parms = np.array([dx, dt, dtout, D, glob_a, glob_b, tout, time], dtype = np.float64)

    # send iparms to everyone
    comm.bcast(iparms, root = 0)

    # send parms to everyone
    comm.bcast(parms, root = 0)

    # time-stepping loop:
    # workers do computation and send output to master every dtout
    # master recieves from each worker its secton of output arrays(s)

    # NOTE TO SELF: INSERT STS SCHEME TO IMPROVE SPEEDUP!!!

    # begin time-stepping
    for nsteps in range(1, MaxSteps + 1):

        # synchronize everyone
        comm.Barrier()

        # update time = time + dt
        time = nsteps * dt

        if time >= tout:

            # recieve output from workers when time = dtout
            U = RECV_output_MPI(nWRs, glob_M, comm)

            # receive exact solution from worker when time = dtout
            u_exact = RECV_output_exact_MPI(nWRs, glob_M, comm)

            # NOTE TO SELF: rearrange to call compare subroutine
            # instead of RECV_output_exact_MPI subroutine
            # ERR, u_exact = COMPARISON(x, U, Mp1, time, D, u_exact)

            # compute error
            ERR = abs(U - u_exact)

            # compute max error
            maxERR = np.amax(ERR)

            # NOTE TO SELF: put mesh into output file
            # call output subroutine and print info to file
            print('# U, u_exact, error at time = {}' .format(time), file = outf)
            OUTPUT(outf, U, u_exact, ERR)
            print(' ', file = outf)
            print('Maximum error = {} at time = {}' .format(maxERR, time), file = outf)
            print(' ', file = outf)

            # update tout
            tout = tout + dtout

    # end time-stepping

    # print run-time information to screen
    print('MASTER DONE: exiting at t = %f after %i steps.' % (time, nsteps))
    print('MASTER DONE: maximum error = {}' .format(maxERR))

    outf.close()
