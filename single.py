# -*- coding: utf-8 -*-
"""
This is a program that calculates the production and degradation of a diffusible component (U) and its membrane association (M) at the both ends of the cell.
"""

#----Package---
import numpy as np
from scipy import sparse
import time

def find_dx_dt(x0, xL, Nx, D):
    # calculate dx
    dx = (xL - x0)/(Nx - 1)
    # estimate time step
    dt = 0.1*(dx**2/(2*D)) 
    return dx, dt

def f(m):
    delta = 0.2 
    Km =30
    n=6
    mm = (delta + (1-delta)*(m**n/(Km**n+ m**n)))
    return mm

#-----------------------------------------
Nx = 10 # GRID POINTS on space interval
x0 = 0
xL = 2 #um
t0 = 0
tFinal=3600 #seconds
D = 0.32  #um**2/s # difusion 
dx, dt = find_dx_dt(x0, xL, Nx, D)  # GRID POINTS on time interval
Nt = int(tFinal/dt)
tF = dt*Nt
xspan = np.linspace(x0, xL, Nx)
tspan = np.linspace(t0, tF, Nt+1)

uss = 60#1000
mss = 3000#5000

ks=0.0024
kd=0.0005
kp = kd*uss
kb = ks*mss/uss

# diffusion matrix A
s = D*dt/dx**2
a0 = 1 - 2*s
main_diag_a0 = a0*np.ones((1,Nx))
off_diag_a0 = s*np.ones((1, Nx-1))
a = main_diag_a0.shape[1]
diagonalsA = [main_diag_a0, off_diag_a0, off_diag_a0]
A = sparse.diags(diagonalsA, [0,-1,1], shape=(a,a)).toarray()
A[0,0] = 1-s
A[-1,-1] = 1-s

#%%
U = np.zeros((Nx, Nt+1))
M = np.zeros((Nx, Nt+1))
# ---initial condition
U[0,0] = 0
M[0,0] = 0
start_time = time.time()

# constant production for u at the first grid
u2 = np.zeros((Nx,1)).ravel() 
u2[0] = dt*kp 

for k in range(0, Nt):
    #---dm/dt---
    m1 = kb*U[0, k]*f(M[0, k])-ks*M[0, k]
    m2 = kb*U[-1, k]*f(M[-1, k])-ks*M[-1, k]
        
    M[0, k+1] = M[0, k] + dt*m1
    M[-1, k+1] = M[-1, k] + dt*m2    
        #---du/dt            
    #degrdation for u
    u3 = -dt*kd*U[:, k]
    # association to membrane
    u4 = np.zeros((Nx,1)).ravel()
    u4[0] = -m1*dt
    u4[-1] = -m2*dt
    # diffusion    
    U[:,k+1] = A.dot(U[:, k]) + u2 + u3 + u4
        
    #print(str(k))

# change sec into min, allocating a smaller array U_n for making plots
p_num = 60/dt
print_n = np.arange(0, Nt, p_num)
tspan_n = np.zeros(len(print_n))
U_n = np.zeros([Nx, len(print_n)])
M_n =  np.zeros([Nx, len(print_n)])

for i, j  in enumerate(print_n):
    tspan_n[i] = tspan[int(j)]
    U_n[:, i] = U[:,int(j)]
    M_n[:, i] = M[:,int(j)]

end_time = time.time()
print('Simulation time: {0} sec'.format(end_time - start_time))

ratio=M_n[0,-1]/M_n[-1,-1]
print('kd=', kd, 'ks=', ks, 'ratio=', ratio)
data = [tspan_n, xspan, U_n,  M_n]
np.save('single_data',data)


