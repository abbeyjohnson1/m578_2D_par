# Abbey Johnson
# February 22, 2019

# import math and mumpy modules
import math
import numpy as np

# import subroutines from other files
from z_io import INPUT, OUTPUT, COMPARISON
from z_setup import MESH, INIT
from z_update import FLUX, PDE

# check to see if Lab $ is being run directly
if __name__ == "__main__":
    print ("Lab4 1D_serial is being run directly")
else:
    print ("Lab4 is being imported")

# create output file
filename = 'data_out.txt'

# BEGIN TIME-STEPPING SCHEME

# input parameter valus from file
MM, tend, factor, dtout, D, a, b = INPUT('./Init.txt')

dx = np.float64(1.0/MM)

# M = int((b-a)/dx)
M = int((b-a)*MM)

Mp1 = int(M+1)

dtEXPL = np.float64((dx*dx)/(2*D))

# time-step (for stability of the explicit scheme)
dt = np.float64(factor*dtEXPL)

# declare variables and dimensions of arrays
x = np.zeros((M+2), dtype=np.float64)
U = np.zeros((M+2), dtype=np.float64)
F = np.zeros((M+2), dtype=np.float64)
u_exact = np.zeros((M+2), dtype=np.float64)

# call mesh subroutine
x = MESH(x, a, dx, Mp1, b)

# initialize
nsteps = 0
time = 0.0
tout = dtout
MaxSteps = int(tend/dt)+1

# call initialization subroutine
U = INIT(U, Mp1)

# call output routine to record the concentration profile at time t = 0
OUTPUT(filename, x, U)

# begin time-stepping
for nsteps in range(1, (MaxSteps+1)):
    # update time = time + dt
    time = nsteps * dt
    # call flux subroutine
    F = FLUX(F, Mp1, D, U, x, b, time)
    # call PDE subroutine
    U = PDE(Mp1, dt, dx, F, U, b, D, time)
    print(U)
    if time >= tout:
        # call compare subroutine
        ERR = COMPARISON(x, U, Mp1, time, D, u_exact)
        # call output subroutine
        OUTPUT(filename, x, U)
        # update tout
        tout = tout + dtout
    # end if-statement
# end time-stepping for-loop

# print run-time information to screen
print('DONE: exiting at t = %f after %i steps.' % (time, nsteps))

# print maximum error at t_end for i=0:M+1
print('Maximum error is %e \n' % ERR)

print('Parameters for this run: ')
print('MM = %f ' % MM)
print('t_end = %f ' % tend)
print('factor = %f ' % factor)
print('dtout = %f ' % dtout)
print('D = %f ' % D)
print('a = %f ' % a)
print('b = %f ' % b)

# END TIME-STEPPING SCHEME
