# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

from mpi4py import mpi4py
import numpy as np

def RECV_output_MPI(comm, U, ):
    # only master calls this subroutine

    # size of U(i,:) array
    I2 = Mr + 2

    # size of U(:,j) array
    J2 = Mz + 2

    for i in range(1, nWRs + 1):
        Jme = (i-1) * Mz + 1
        msgtag = 1000 + i
        msg = comm.recv(source = i, tag = msgtag)
        if (i == 1):
            U = np.append(U, )
        if (i == nWRs):
            U = np.append(U, )
        else:
            U = np.append(U, )
    return(U)

def SEND_output_MPI(comm, Me, NodeUP, NodeDN, Mr, Mz, U):
    # every worker does this

    # everybody sends values to Master for output
    mster = 0

    msgtag = 1000 + Me

    # call mpi_send
    comm.send(U, dest = mster, tag = msgtag)


def EXCHANGE_bry_MPI(comm, nWRs, Me, NodeUP, NodeDN, Mr, Mz, U):

    # test by setting all u values for each worker to myID

    Jup = Mz
    Jup1 = Jup + 1
    I2 = Mr +2

    # send bottom row to neighbor down
    if (Me != 1):
        msg = U[:, I2]
        comm.send(msg, dest = NodeDN, tag = msgDN)

    #  recieve bottom row from neighbor up and save as upper boundary
    if (Me != nWRs):
        msg = comm.recv(None, source = NodeUP, tag = msgDN)
        U[:, Jup1] = msg

    # send the top row to neighbor up
    if (Me != nWRs):
        msg = U[:, Jup]
        comm.send(msg, dest = NodeUP, tag = msgUP)

    # receive top row from neighbor down and save as lower boundary
    if (Me != 1):
        msg = comm.recv(None, source = NodeDN, tag = msgUP)
        U[:, 0] = msg
