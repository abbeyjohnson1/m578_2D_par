# Abbey L. Johnson

# MAIN PROGRAM

# import MPI
from mpi4py import MPI
# import MASTER subroutine
from x_mainMR import MASTER
# import WORKER subroutine
from x_mainWR import WORKER

#>>>mpi>>>

# startup
# MPI.Init() is automatically called when MPI module is imported from mpi4py package
# MPI.Init()

# defines the default communicator, which contains all processes
comm = MPI.COMM_WORLD

# nPROC specified at mpirun/mipexec
# nPROC = number of processes = nWRs+1 to use in this run
nPROC = comm.Get_size()

# returns the process ID of current process; a number between 0 and nPROC-1
# myID = rank of a process
myID = comm.Get_rank()

# master rank
mster = 0

# number of workers
nWRs = nPROC - 1

print(">>>main>>> running on {} workers, myID = {}" .format(nWRs, myID))

# start 0,1,...,nWRs tasks

if (myID == mster):
    # start CPU timer on master
    tt0 = MPI.Wtime()
    # call MASTER
    MASTER(nWRs, mster, comm, myID)
    # end timer
    tt1 = MPI.Wtime()
    # calculate run time
    tt = tt1 - tt0
    print(">>main>> Master timing = {} seconds on {} workers" .format(tt, nWRs))
else:
    # call worker, now MPI is running
    WORKER(nWRs, myID, comm)
    print("Bye-bye from worker: {}" .format(myID))

# termination
# not necessary to call MPI.Finalize() to ensure MPI finalization in Python
# MPI.Finalize()

#<<<mpi<<<
