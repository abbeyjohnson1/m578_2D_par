# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

# imports
from mpi4py import MPI
import numpy as np

# >>> receive output subroutine >>>
def RECV_output_MPI(comm, nWRs, glob_Mr2, glob_Mz):
    # only master calls this subroutine
    U = np.zeros((glob_Mr2, glob_Mz + 2))
    # recalculate local Mz
    loc_Mz = int(glob_Mz / nWRs)
    for k in range(1, nWRs + 1):
        Jme = int((k - 1) * loc_Mz + 1)
        Jend = int(Jme + loc_Mz)
        msgtag = 1000 + k
        msg = comm.recv(source = MPI.ANY_SOURCE, tag = msgtag)
        # append msg to U
        U[ : , Jme : Jend] = msg[ : , 1 : -1]
        # if receiving output from first worker
        if (k == 1):
            # save first column of first worker
            U[ : , 0] = msg[ : , 0]
        # if receiving output from last worker
        if (k == nWRs):
            # save the last column of last worker
            U[ : , -1] = msg[ : , -1]
    return(U)
# <<< end receive output subroutine <<<

# >>> send output subroutine >>>
def SEND_output_MPI(comm, Me, U):
    # every worker sends values to master for output
    mster = 0
    msgtag = 1000 + Me
    # call mpi_send
    comm.send(U, dest = mster, tag = msgtag)
# <<< end send output subroutine <<<

# >>> subroutine for workers exchanging "boundary" values >>>
def EXCHANGE_bry_MPI(comm, nWRs, Me, NodeUP, NodeDN, glob_Mr, loc_Mz, U):
    msgUP = 10
    msgDN = 20
    # send second column to neighbor down
    if (Me != 1):
        msg = U[ : , 1]
        comm.send(msg, dest = NodeDN, tag = msgDN)
    #  recieve second column from neighbor up; save as upper boundary
    if (Me != nWRs):
        msg = comm.recv(None, source = NodeUP, tag = msgDN)
        U[ : , -1] = msg
    # send the second-to-last column to neighbor up
    if (Me != nWRs):
        msg = U[ : , -2]
        comm.send(msg, dest = NodeUP, tag = msgUP)
    # receive second-to-last column from neighbor down; save as lower boundary
    if (Me != 1):
        msg = comm.recv(None, source = NodeDN, tag = msgUP)
        U[ : , 0] = msg
# >>> end exchange boundary subroutine >>>
