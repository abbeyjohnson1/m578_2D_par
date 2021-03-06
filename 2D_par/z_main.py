# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

# problem describes axially symmetric heat conduction in a hollow cylinder:
# Rin < r < Rout and 0 < z < Z
# imposed temperature conditions on all boundaries, starting with given temp.

# note: this problem has exact solution

# note: code is parallelized along the z-direction

# MAIN PROGRAM
# terminal prompt: mpirun -n nPROCS python z_main.py

# imports
from mpi4py import MPI
from z_mainMR import MASTER
from z_mainWR import WORKER

# >>> start mpi >>>
# defines default communicator, which contains all processes
comm = MPI.COMM_WORLD
# nPROC = number of processes = nWRs + 1
nPROC = comm.Get_size()
# returns process ID of current process; number between 0 and nPROC - 1
# myID = rank of process
myID = comm.Get_rank()
# master rank
mster = 0
# number of workers
nWRs = nPROC - 1
# print('>>> main >>> running on {} workers, myID = {}' .format(nWRs, myID))
# start 0, 1, ... , nWRs tasks
# master
if (myID == mster):
    # start CPU timer on master
    tt0 = MPI.Wtime()
    # call master
    MASTER(comm, nWRs, mster)
    # end timer
    tt1 = MPI.Wtime()
    # calculate run time
    tt = tt1 - tt0
    print('>>> main >>> master timing = {} seconds on {} workers' .format(tt, nWRs))
# workers
else:
    # call worker, now MPI is running
    WORKER(comm, nWRs, myID)
    # print('bye-bye from worker number {}' .format(myID))
# <<< end mpi <<<
