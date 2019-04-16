# Abbey L. Johnson

# 3 subroutines for message passing

# import MPI
from mpi4py import MPI

# only for the master
def RECV_output_MPI(nWRs, M, U, comm):
    # receive values from everybody for output
    for i in range(1, nWRs + 1):
        # Jme = (i-1)* M + 1
        msgtag = 1000 + i
        msg = comm.recv(source = MPI.ANY_SOURCE, tag = msgtag)
        print('OUTPUT', msg)

# for the workers
def SEND_output_MPI(Me, NodeUP, NodeDN, M, U, comm, nWRs):
    msgtag = 1000 + Me
    comm.send(U, dest = 0, tag = msgtag)

# for the workers exchanging their "boundary" values
def EXCHANGE_bry_MPI(nWRs, Me, NodeUP, NodeDN, M, U, comm):

    Jup = M
    Jup1 = Jup + 1

    msgUP = 10
    msgDN = 20

    # I2 = M + 2

    # send bottom row to neighbor down
    if (Me != 1):
        msg = U[1]
        comm.send(msg, dest = NodeDN, tag = msgDN)

    #  recieve bottom row from neighbor up and save as upper boundary
    if (Me != nWRs):
        msg = comm.recv(None, source = NodeUP, tag = msgDN)
        #U[Jup1] = msg
        U[-1] = msg

    # send the top row to neighbor up
    if (Me != nWRs):
        # msg = U[Jup]
        msg = U[-2]
        comm.send(msg, dest = NodeUP, tag = msgUP)

    # receive top row from neighbor down and save as lower boundary
    if (Me != 1):
        msg = comm.recv(None, source = NodeDN, tag = msgUP)
        U[0] = msg
