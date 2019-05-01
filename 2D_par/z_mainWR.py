# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

# >>> workers >>>

# imports
import math
import numpy as np
from mpi4py import MPI
from z_messaging import EXCHANGE_bry_MPI, SEND_output_MPI
from z_setup import INIT, MESH
from x_update import PDE, FLUX

# NEED ARGUMENTS
def WORKER(comm):
