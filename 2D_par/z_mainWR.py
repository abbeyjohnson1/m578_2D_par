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
def WORKER(comm, nWRs, Me):

    # check to see if worker file is running
    print('z_mainWR number {} is being run' .format(Me))

    # upack iparms arrray that was broadcasted
    iparms = comm.bcast(None, root = 0)
    glob_Mr, glob_Mr1, glob_Mr2, glob_Mz, MaxSteps, nsteps= iparms

    # upack parms array that was broadcasted
    parms = comm.bcast(None, root = 0)
    dr, dz, dt, tout, dtout, D, Rin, Rout, Z, time = parms

    # declare variables local to workers

    # local Mz (number of nodes/control volumes)
    loc_Mz = int(glob_Mz / nWRs)
    loc_Mz1 = int(loc_Mz + 1)
    loc_Mz2 = int(loc_Mz + 2)

    NodeUP = Me + 1
    NodeDN = Me - 1

# need arguments
    r, z, Ar, Az, dV = MESH()
    print('Me = {}, r mesh = {}' .format(Me, r))
    print('Me = {}, z mesh = {}' .format(Me, z))

    # initialize U
    U = INIT(glob_Mr2, loc_Mz2, r, z)

    # time-stepping loop:
    # workers do computation and send output to master every dtout
    # master receives from each worker its section of output arrays

    # begin time-stepping
    for nsteps in range(1, MaxSteps):

        # synchronize everyone
        comm.Barrier()

        # exchange "boundary" values with neighboring workers
        EXCHANGE_bry_MPI(comm, nWRs, Me, NodeUP, NodeDN, glob_Mr, loc_Mz, U)

        # update time = time + dt
        time = nsteps * dt

        # call flux subroutine
        U, Fr, Fz = FLUX(glob_Mr1, glob_Mr2, loc_Mz1, loc_Mz2,
        Me, nWRs, U, time, r, z, D)

        # call PDE subroutine
        U = PDE(glob_Mr1, loc_Mz1, Ar, Fr, Az, Fz, U, dt, dV)

        # send output to master every dtout
        if time >= tout:

            # send output from worker to master
            SEND_output_MPI(comm, Me, U)

            # update tout
            tout = tout + dtout

        # end time-stepping

# print run-time information to screen
print('WORKER DONE: exiting at t = %f after %i steps.' % (time, nsteps))

# <<< end worker <<<
