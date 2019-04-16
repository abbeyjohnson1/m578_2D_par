# Abbey L. Johnson

# NOTE TO SELF: the send and receive subroutines need improvement
# NOTE TO SELF: ideally I want just 1 recv function and 1 send subroutine

# 3 subroutines for message passing

# import MPI
from mpi4py import MPI
import numpy as np

# only for the master
def RECV_output_MPI(nWRs, glob_M, comm):
    # receive values from everybody for output
    # NOTE TO SELF: need to make this more generalized!!!
    M = int(glob_M / nWRs)
    msg = np.zeros(nWRs)
    U = np.array([1.0])
    for i in range(1, nWRs + 1):
        msgtag = 1000 + i
        msg = comm.recv(source = MPI.ANY_SOURCE, tag = msgtag)
        # print('OUTPUT', msg)
        U = np.append(U, msg[1:M+1])
        if (i == nWRs):
            U = np.append(U, msg[M+1])
    return(U)

# only for the master (to recieve u_exact)
def RECV_output_exact_MPI(nWRs, glob_M, comm):
    # receive exact U values from everybody for output
    # NOTE TO SELF: need to make this more generalized!!!
    M = int(glob_M / nWRs)
    msg = np.zeros(nWRs)
    u_exact = np.array([1.0])
    for j in range(1, nWRs + 1):
        msgtag = 3000 + j
        msg = comm.recv(source = MPI.ANY_SOURCE, tag = msgtag)
        u_exact = np.append(u_exact, msg[1:M+1])
        if (j == nWRs):
            u_exact = np.append(u_exact, msg[M+1])
    return(u_exact)

# for the workers
def SEND_output_MPI(Me, NodeUP, NodeDN, M, U, comm, nWRs):
    msgtag = 1000 + Me
    comm.send(U, dest = 0, tag = msgtag)

# for the workers (to send u_exact)
def SEND_output_exact_MPI(Me, NodeUP, NodeDN, M, u_exact, comm, nWRs):
    msgtag = 3000 + Me
    comm.send(u_exact, dest = 0, tag = msgtag)

# for the workers exchanging their "boundary" values
def EXCHANGE_bry_MPI(nWRs, Me, NodeUP, NodeDN, M, U, comm):

    Jup = M
    Jup1 = Jup + 1

    msgUP = 10
    msgDN = 20

    # send bottom row to neighbor down
    if (Me != 1):
        msg = U[1]
        comm.send(msg, dest = NodeDN, tag = msgDN)

    #  recieve bottom row from neighbor up and save as upper boundary
    if (Me != nWRs):
        msg = comm.recv(None, source = NodeUP, tag = msgDN)
        U[Jup1] = msg

    # send the top row to neighbor up
    if (Me != nWRs):
        msg = U[Jup]
        comm.send(msg, dest = NodeUP, tag = msgUP)

    # receive top row from neighbor down and save as lower boundary
    if (Me != 1):
        msg = comm.recv(None, source = NodeDN, tag = msgUP)
        U[0] = msg
