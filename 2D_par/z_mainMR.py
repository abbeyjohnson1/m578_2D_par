# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

# imports
import numpy as np
from mpi4py import MPI
from z_messaging import RECV_output_MPI
from z_io import INPUT, OUTPUT, OUTPUT_runtime, COMPARISON
from z_setup import glob_MESH, INIT

# >>> master >>>
def MASTER(comm, nWRs, myID):
    # check to see if z_mainMR file is running
    # print('z_mainMR is being run')
    # read input data from file
    input = INPUT('./Init.txt')
    # assign parameter values, always use even Mz divisible by nWRs!!!
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
    # open files for outputting info
    outf = open('o_out.txt', 'w')
    prof = open('o_prof.txt', 'w')
    comp = open('o_compare.txt', 'w')
    # compute other parameters
    # height of cylindar
    Z = np.pi
    # dr
    dr = np.float64(1.0 / MMr)
    # global Mr values
    glob_Mr = int( (Rout - Rin) * MMr )
    glob_Mr1 = int(glob_Mr + 1)
    glob_Mr2 = int(glob_Mr + 2)
    # dz
    dz = np.float64(1.0 / MMz)
    # global Mz values
    glob_Mz = int( (Z - 0.0) / dz )
    glob_Mz1 = int(glob_Mz + 1)
    glob_Mz2 = int(glob_Mz + 2)
    # time-step dt (for stability of explicit scheme)
    dtEXPL = np.float64(1 / (2 * D * (1 / (dr * dr) + 1 / (dz * dz) ) ) )
    dt = np.float64(factor * dtEXPL)
    # maximum number of iterations
    MaxSteps = int(tend / dt) + 2
    # when to print outputs
    tout = dtout
    # initialize time-stepping variables
    nsteps = 0
    time = 0.0
    # pack integers in iparms array
    iparms = np.array([glob_Mr, glob_Mr1, glob_Mr2, glob_Mz, MaxSteps, nsteps], dtype = np.int)
    # pack reals in parms array
    parms = np.array([dr, dz, dt, tout, dtout, D, Rin, Rout, Z, time, tend], dtype = np.float64)
    # send parms arrays to everyone
    comm.bcast(iparms, root = 0)
    comm.bcast(parms, root = 0)
    # compute global mesh
    glob_r, glob_z = glob_MESH(glob_Mr1, glob_Mr2, glob_Mz1, glob_Mz2, Rin, dr, Rout, dz, Z)
    # record initial information at time = 0
    # U = INIT(glob_Mr2, glob_Mz2, glob_r, glob_z)
    # record concentration profile at time = 0
    # OUTPUT(prof, time, U, glob_r, glob_z)
    # time-stepping loop:
    # workers do computation and send output to master every dtout
    # master recieves from each worker its secton of output arrays(s)
    for nsteps in range(1, MaxSteps):
        # synchronize everyone
        comm.Barrier()
        # update time = time + dt
        time = nsteps * dt
        # check if ready to output
        if time >= tout:
            # recieve output from workers when time = dtout
            U = RECV_output_MPI(comm, nWRs, glob_Mr2, glob_Mz)
            # call comparison subroutine
            ERR = COMPARISON(comp, time, U, glob_Mr2, glob_Mz2, glob_r, glob_z)
            # record concentration profile
            # OUTPUT(prof, time, U, glob_r, glob_z)
            # print error at dtout to file
            print('Maximum error = {} at time = {} \n' .format(ERR, time), file = outf)
            # update tout
            tout = tout + dtout
    # end time-stepping
    # print run-time information to screen
    # print('MASTER DONE: exiting at t = %f after %i steps.' % (time, nsteps))
    print('MASTER DONE: maximum error at end time = {}' .format(ERR))
    # print summary of information to output file
    outputs = np.array([time, nsteps, ERR, MMr, MMz, Rin, Rout, tend, factor, dtout, D, nWRs])
    OUTPUT_runtime(outf, outputs)
    # close files
    prof.close()
    comp.close()
    outf.close()
# <<< end master <<<
