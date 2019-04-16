# Abbey L. Johnson

# necessary imports
import math
import numpy as np
from mpi4py import MPI
from x_messaging import EXCHANGE_bry_MPI, SEND_output_MPI, SEND_output_exact_MPI
from x_setup import INIT, MESH
from x_update import PDE, FLUX, COMPARISON

def WORKER(nWRs, Me, comm):

    # check to see if worker.py file is running
    # print("z_mainWR number {} is being run" .format(Me))

    # unpack iparms array that was broadcasted
    iparms = comm.bcast(None, root = 0)
    glob_M, MaxSteps, nsteps = iparms

    # unpack parms array that was broadcasted
    parms = comm.bcast(None, root = 0)
    dx, dt, dtout, D, glob_a, glob_b, tout, time = parms

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

    # initialize arrays locally for each worker
    F = np.zeros((M + 2), dtype = np.float64)
    u_exact = np.zeros((M + 2), dtype = np.float64)

    # timestepping loop:
    # workers do computation and send output to master every dtout
    # master recieves from each worker its section of output array(s)

    # NOTE TO SELF: IMPLEMENT STS SCHEME TO IMPROVE SPEEDUP

    # begin time-stepping
    for nsteps in range(1, MaxSteps + 1):

        # synchronize everyone
        comm.Barrier()

        # exchange "boundary" values with neighbors
        EXCHANGE_bry_MPI(nWRs, Me, NodeUP, NodeDN, M, U, comm)
        # print("Me = {}, U = {}" .format(Me, U))

        # update time = time + dt
        time = nsteps * dt

        # call flux subroutine
        F = FLUX(nWRs, Me, F, Mp1, D, U, x)

        # call PDE subroutine
        U = PDE(nWRs, Me, Mp1, dt, dx, F, U, D, time, x)
        # print("Me = {}, U = {}" .format(Me, U))

        # send output to master every dtout
        if time >= tout:

            u_exact = COMPARISON(x, U, Mp1, time, D, u_exact)

            SEND_output_MPI(Me, NodeUP, NodeDN, M, U, comm, nWRs)

            SEND_output_exact_MPI(Me, NodeUP, NodeDN, M, u_exact, comm, nWRs)

            # update tout
            tout = tout + dtout

    # end time-stepping

# print run-time information to screen
# print('WORKER DONE: exiting at t = %f after %i steps.' % (time, nsteps))
