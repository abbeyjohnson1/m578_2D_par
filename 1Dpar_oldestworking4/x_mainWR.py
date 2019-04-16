# Abbey L. Johnson

# necessary imports
import math
import numpy as np
from mpi4py import MPI
from x_messaging import EXCHANGE_bry_MPI, SEND_output_MPI
from x_setup import INIT, MESH
from x_update import PDE, FLUX

def WORKER(nWRs, Me, comm):

    iparms = comm.bcast(None, root = 0)
    parms = comm.bcast(None, root = 0)

    comm.Barrier()

    # check to see if worker.py file is running
    # print("z_mainWR number {} is being run" .format(Me))

# check to see if I can get rid of any of these

    # unpack initial iparms array that was broadcasted

    glob_M, Mp1, MaxSteps = iparms

# check to see if I can get rid of any of these

    # unpack parms array that was broadcasted

    dx, dtEXPL, dt, tend, factor, dtout, D, glob_a, glob_b, tout = parms

    # declare variables local to workers

    # local M (number of nodes/control volumes)
    M = int(glob_M / nWRs)
    Mp1 = M + 1

    NodeUP = Me + 1
    NodeDN = Me - 1

    # declare mesh local to worker
    x = np.zeros((M + 2), dtype = np.float64)
    x = MESH(x, glob_a, dx, M, Mp1, glob_b, Me, nWRs)
    # print('Me = {}, Mesh = {}' .format(Me, x))

    # initialize U array in the worker
    U = INIT(Mp1, Me)

# not sure if this is necessary????
    # initialize other arrays locally for each worker
    F = np.zeros((M + 2), dtype = np.float64)
# not sure if this is necessary????
    u_exact = np.zeros((M + 2), dtype = np.float64)

    # timestepping loop:
    # workers do computation and send output to master every dtout
    # master recieves from each worker its section of output array(s)

# i could also just broadcast these from parms and iparms in master???
    # initialize time-stepping variables
    nsteps = 0
    time = 0.0

    # LATER: IMPLEMENT STS SCHEME TO IMPROVE SPEEDUP

    # begin time-stepping
    for nsteps in range(1, MaxSteps + 1):

        # synchronize everyone
        comm.Barrier()

        # exchange "boundary" values with neighbors
        EXCHANGE_bry_MPI(nWRs, Me, NodeUP, NodeDN, M, U, comm)
        # print("Me = {}, U = {}" .format(Me, U))

        # update time = time + dt
        time = nsteps * dt

# check to see if I can get rid of any arguments
        # call flux subroutine
        F = FLUX(nWRs, Me, F, Mp1, D, U, x, time, glob_b)

# check to see if I can get rid of any arguments
        # call PDE subroutine
        U = PDE(nWRs, Me, Mp1, dt, dx, F, U, D, time, x, glob_b)
        # print("Me = {}, U = {}" .format(Me, U))

        # send output to master every dtout
        if time >= tout:

            SEND_output_MPI(Me, NodeUP, NodeDN, M, U, comm, nWRs)

            # update tout
            tout = tout + dtout

    # end time-stepping

# print run-time information to screen
# print('WORKER DONE: exiting at t = %f after %i steps.' % (time, nsteps))
