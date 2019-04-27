# Abbey L. Johnson
# Axially Symmetric Heat Transfer
# Parallelized Explicit Finite Volume Code

# >>> workers >>>

# import modules
import numpy as np
import math

# import subroutines from other files
from z_io import INPUT, OUTPUT
from z_setup import INIT, MESH
from z_update import FLUX, PDE, COMPARISON

# Abbey L. Johnson
# Lab 6: Axially Symmetric Heat Transfer
# Explicit Finite Volume Code

# problem describes axially symmetric heat conduction in a hollow cylinder:
# Rin < r < Rout and 0 < z < Z
# imposed temperature conditions on all boundaries, starting with given temp.
# note: this problem has exact solution

###############################################################################

# import modules
import numpy as np

# check to make sure file is being run directly
if __name__ == "__main__":
   print ("Lab 9 serial is being run directly \n")
else:
   print ("Lab 9 serial is being imported \n")

# open files for outputting info
prof = open('o_prof.txt', 'w')
comp = open('o_compare.txt', 'w')
outf = open('o_out.txt', 'w')

###############################################################################

# read parameters from init data file
with open("./Init.txt") as init_file:
    for line in init_file:
        line = line.partition('#')[0]
        contents = line.split()

# assign parameter values
# nodes in r-direction
MMr = np.float64(contents[0])
# nodes in z-direction
MMz = np.float64(contents[1])
# inner radius of cylinder
Rin = np.float64(contents[2])
# outer radius of cylinder
Rout = np.float64(contents[3])
# end time
tend = np.float64(contents[4])
# factor for stability
factor = np.float64(contents[5])
# time-step when to output
dtout = np.float64(contents[6])
# diffusivity
D = np.float64(contents[7])

# set pi to full precision
Pi = 4 * np.arctan(1.0)
# height of cylinder
Z = Pi
# note: maybe should read Z in from init file?
# note: Z = pi could create problems, so code fluxes for non-uniform grid

dr = np.float64(1.0 / MMr)
dz = np.float64(1.0 / MMz)

# Mr = int( (Rout - Rin) / dr )
Mr = int( (Rout - Rin) * MMr )
Mr1 = int(Mr + 1)

# Mz = int( (Z - 0.0) / dz )
Mz = int( (Z - 0.0) * MMz )
Mz1 = int(Mz + 1)

dtEXPL = np.float64(1 / (2 * D * (1/(dr * dr) + 1/(dz * dz))))

# time-step (for stability of the explicit scheme)
dt = np.float64(factor * dtEXPL)

# declare variables and dimensions of arrays
r = np.zeros(Mr + 2, dtype = np.float64)
z = np.zeros(Mz + 2, dtype = np.float64)
Ar = np.zeros((Mr + 2, Mz + 2), dtype = np.float64)
Az = np.zeros((Mr + 2, Mz + 2), dtype = np.float64)
Fr = np.zeros((Mr + 2, Mz + 2), dtype = np.float64)
Fz = np.zeros((Mr + 2, Mz + 2), dtype = np.float64)
dV = np.zeros((Mr + 2, Mz + 2), dtype = np.float64)
U = np.zeros((Mr + 2, Mz + 2), dtype = np.float64)
u_exact = np.zeros((Mr + 2, Mz + 2), dtype = np.float64)

# call mesh subroutine
r, z, Ar, Az, dV = MESH(r, Rin, dr, Mr1, Rout, z, dz, Mz1, Z, Mr, Ar, Mz, Az, dV)
# never forget to check the mesh
# print('r = {}' .format(r))
# print('z = {}' .format(z))

# initialize
nsteps = 0
time = 0.0
tout = dtout
MaxSteps = int(tend / dt) + 1

# call initialization subroutine
U = INIT(Mr, Mz, U, r, z)

# call output routine to record the concentration profile at time t = 0
# OUTPUT(time, r, Mr, z, Mz, U)

###############################################################################

# note: should insert STS scheme for speedup!!!

# begin time-stepping
for nsteps in range(1, MaxSteps + 1):
    # update time = time + dt
    # BCs are time dependent, so update time before calling FLUX
    time = nsteps * dt
    # call flux subroutine
    U, Fr, Fz = FLUX(r, Mr1, z, Mz1, U, time, Fr, D, Fz)
    # call PDE subroutine
    U = PDE(Mr1, Mz1, Ar, Fr, Az, Fz, U, dt, dV)
    # check if ready to output
    if time >= tout:
        # call compare subroutine
        ERR = COMPARISON(Mr, Mz, u_exact, time, r, z, U)
        # call output subroutine
        # OUTPUT(time, r, Mr, z, Mz, U)
        # update tout
        tout = tout + dtout

# print run-time information to screen
print('DONE: exiting at t = %f after %i steps. \n' % (time, nsteps))
print('Maximum error is %e at time %f \n' % (ERR, time))

# end time-stepping

###############################################################################

# print final outputs

print('Parameters for this run: ', file = outf)
print('MMr = %f ' % MMr, file = outf)
print('MMz = %f ' % MMz, file = outf)
print('Rin = %f ' % Rin, file = outf)
print('Rout = %f' % Rout, file = outf)
print('t_end = %f ' % tend, file = outf)
print('factor = %f ' % factor, file = outf)
print('dtout = %f ' % dtout, file = outf)
print('D = %f \n' % D, file = outf)


# print run-time information to output file
print('DONE: exiting at t = %f after %i steps. \n' % (time, nsteps), file = outf)

# print maximum error at t_end for i=0:M+1
print('Maximum error is %e at time %f \n' % (ERR, time), file = outf)

# close files
init_file.close()
prof.close()
comp.close()
outf.close()

################################################################################
