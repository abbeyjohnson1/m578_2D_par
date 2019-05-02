# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

# >>> worker >>>

# imports
import math
import numpy as np
from mpi4py import MPI
from z_messaging import EXCHANGE_bry_MPI, SEND_output_MPI
from z_setup import INIT, MESH
from x_update import PDE, FLUX

# NEED ARGUMENTS
def WORKER(comm):

    # check to see if worker file is running
    print('z_mainWR number {} is being run' .format(ME))

    # upack iparms arrray that was broadcasted
    iparms = comm.bcast(None, root = 0)
    ?, ?, ? = iparms

    # upack parms array that was broadcasted
    parms = comm.bcast(None, root = 0)
    ?, ?, ? = parms

    # declare variables local to workers

    # local Mz (number of nodes/control volumes)
    loc_Mz = int(glob_Mz / nWRs)
    loc_Mz1 = int(loc_Mz + 1)
    loc_Mz2 = int(loc_Mz + 2)

    NodeUP = Me + 1
    NodeDN = Me - 1

    # declare mesh, which is local to worker
    r = np.zeros(glob_Mr + 2, dtype = np.float64)
    z = np.zeros(loc_Mz + 2, dtype = np.float64)
    Ar = np.zeros((glob_Mr + 2, loc_Mz + 2), dtype = np.float64)
    Az = np.zeros((glob_Mr + 2, loc_Mz + 2), dtype = np.float64)
    dV = np.zeros((glob_Mr + 2, loc_Mz + 2), dtype = np.float64)
# need arguments
    r, z, Ar, Az, dV = MESH()
    print('Me = {}, r mesh = {}' .format(Me, r))
    print('Me = {}, z mesh = {}' .format(Me, z))

    # initialize U
# need ARGUMENTS
    U = INIT()

    # initialize arrays locally for each worker
    Fr = np.zeros((glob_Mr + 2, loc_Mz + 2), dtype = np.float64)
    Fz = np.zeros((glob_Mr + 2, loc_Mz + 2), dtype = np.float64)

    # time-stepping loop:
    # workers do computation and send output to master every dtout
    # master receives from each worker its section of output arrays

    # begin time-stepping
    for nsteps in range(1, MaxSteps):

        # synchronize everyone
        comm.Barrier()

        # exchange "boundary" values with neighbros
        EXCHANGE_bry_MPI()

        # update time = time + dt
        time = nsteps * dt

# need arguments here
        # call flux subroutine
        FLUX()

# need arguments here
        # call PDE subroutine
        PDE()

        # send output to master every dtout
        if time >= tout:

            SEND_output_MPI()

            # update tout
            tout = tout + dtout

        # end time-stepping

# print run-time information to screen
print('WORKER DONE: exiting at t = %f after %i steps.' % (time, nsteps))

# <<< end worker <<<
