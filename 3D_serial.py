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
   print ("Lab 6 is being run directly \n")
else:
   print ("Lab 6 is being imported \n")

# open files for outputting info
prof = open('o_prof.txt', 'w')
comp = open('o_compare.txt', 'w')
outf = open('o_out.txt', 'w')

###############################################################################

# begin subroutines

# 2D mesh, control volumes, radial areas, and axial areas
def MESH(r, Rin, dr, Mr1, Rout, z, dz, Mz1, Z, Mr, Ar, Mz, Az, dV):
    # mesh in r-direction
    # left boundary
    r[0] = Rin
    # internal nodes
    r[1] = Rin + dr / 2
    for i in range(2, Mr1):
        r[i] = r[1] + (i - 1) * dr
    # right boundary
    r[Mr1] = Rout
    # mesh in z-direction
    # left boundary
    z[0] = 0.0
    # internal nodes
    z[1] = 0.0 + dz / 2
    for j in range(2, Mz1):
        z[j] = z[1] + (j - 1) * dz
    # right boundary
    z[Mz1] = Z
    # areas of radial faces
    for i in range(1, Mr + 2):
        for j in range(1, Mz1):
            Ar[i][j] = 2 * np.pi * dz * (r[i] - dr / 2)
    # areas of axial faces
    for i in range(1, Mr1):
        for j in range(1, Mz + 2):
            Az[i][j] = 2 * np.pi * dr * r[i]
    # control volumes
    for i in range(1, Mr1):
        for j in range(1, Mz1):
            dV[i][j] = 2 * np.pi * dr * dz * r[i]
    return(r, z, Ar, Az, dV)

# initialization subroutine (for initial temperature profile/initial condition)
def INIT(Mr, Mz, U, r, z):
    # impose initial condition, u(r,z,0)
    for i in range(0, Mr + 2):
        for j in range(0, Mz + 2):
            # note: natural log, not log base 10
            U[i][j] = np.log(r[i]) * np.sin(z[j])
    return(U)

# output subroutine
# def OUTPUT(time, r, Mr, z, Mz, U):
#    # print temperature profile to output file
#    prof.write('# temperature profile of U(r,z,t) at t = %f \n' % time)
#    prof.write('# i j r(i) z(j) U(r,z,t) \n')
#    # note: maybe should print to file in matrix format?
#    for i in range(0, Mr + 2):
#        for j in range (0, Mz + 2):
#            prof.write('%i %i %f %f %e \n' % (i, j, r[i], z[j], U[i][j]))
#    prof.write('\n')

# flux subroutine
def FLUX(r, Mr1, z, Mz1, U, time, Fr, D, Fz):
    # boundary conditions
    for i in range(1, Mr1):
        U[i][0] = 0.0
        U[i][Mz1] = 0.0
    for j in range(1, Mz1):
        U[0][j] = 0.0
        U[Mr1][j] = np.exp(-time) * np.log(2.) * np.sin(z[j])
    # radial fluxes
    for i in range(1, Mr1 + 1):
        for j in range(1, Mz1):
            Fr[i][j] = -D * (U[i][j] - U[i-1][j]) / (r[i] - r[i-1])
    # axial fluxes
    for i in range(1, Mr1):
        for j in range(1, Mz1 + 1):
            Fz[i][j] = -D * (U[i][j] - U[i][j-1]) / (z[j] - z[j-1])
    return(U, Fr, Fz)

# PDE subroutine
def PDE(Mr1, Mz1, Ar, Fr, Az, Fz, U, dt, dV):
    for i in range(1, Mr1):
        for j in range(1, Mz1):
            # radial fluxes and areas
            radial = Ar[i][j] * Fr[i][j] - Ar[i+1][j] * Fr[i+1][j]
            # axial fluxes and areas
            axial = Az[i][j] * Fz[i][j] - Az[i][j+1] * Fz[i][j+1]
            # update (internal) temperatures
            U[i][j] = U[i][j] + (dt / dV[i][j]) * (radial + axial)
    return(U)

# subroutine comparing approximation with exact solution
def COMPARISON(Mr, Mz, u_exact, time, r, z, U):
    # comp.write('# i j U(i,j) u_exact(i,j) error(i,j) at time = %f \n' % time)
    ERR = 0.0
    for i in range(0, Mr + 2):
        for j in range(0, Mz + 2):
            # compute exact solution
            u_exact[i][j] = np.exp(-time) * np.log(r[i]) * np.sin(z[j])
            # compute error at r[i], z[j]
            ERRij = abs(U[i][j] - u_exact[i][j])
            # print to comparison output file
            # comp.write('%i %i %e %e %e \n' % (i, j, U[i][j], u_exact[i][j], ERRij))
            # compute max error
            ERR = max(ERRij, ERR)
    # comp.write('\n')
    # comp.write('Maximum error at time %f is %e \n' % (time, ERR))
    # comp.write('\n')
    return(ERR)

# end subroutines

###############################################################################

# read parameters from init data file
with open("./Lab6Init.txt") as init_file:
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
