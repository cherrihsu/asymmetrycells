# -*- coding: utf-8 -*-
"""
This is a program that simulates the production and degradation of a diffusible component (U) and its membrane association (M) at the both ends of the cell.  In particular we scan for the degradation rates for U and dissociation rates for M
"""
import numpy as np
import multiprocessing as mp
import time  
from functools import partial
from scipy import sparse

#def saveList(myList,filename):
#    # the filename should mention the extension 'npy'
#    np.save(filename, myList)
#    print("Saved successfully!")

def find_dx_dt(x0, xL, Nx, D):
    # calculate dx
    dx = (xL - x0)/(Nx - 1)
    # estimate time step
    dt = 0.1*(dx**2/(2*D)) 
    return dx, dt

def f(m):
    delta = 0.2 
    Km =30
    mm = (delta + (1-delta)*(m**6/(Km**6+ m**6)))
    return mm


def main(para):
    #-----------------------------------------
    D = 0.32 #um**2/s -> um**2/min # difuusion 
    Nx = 10 # GRID POINTS on space interval
    x0 = 0
    xL = 2 #um
    t0 = 0
    dx, dt = find_dx_dt(x0, xL, Nx, D)
    tFinal=3600
    Nt = int(tFinal/dt)
    tF = dt*Nt

    uss = 60#1000
    mss = 3000#5000

    kd = para[0]
    ks = para[1]
    kp = kd*uss
    kb = ks*mss/uss
    
    U = np.zeros((Nx, Nt+1))
    M = np.zeros((Nx, Nt+1))
    U[0,0] = 0
    M[0,0] = 0

    s = D*dt/(dx**2)
    a0 = 1 - 2*s
    main_diag_a0 = a0*np.ones((1,Nx))
    off_diag_a0 = s*np.ones((1, Nx-1))
    a = main_diag_a0.shape[1]
    diagonalsA = [main_diag_a0, off_diag_a0, off_diag_a0]
    A = sparse.diags(diagonalsA, [0,-1,1], shape=(a,a)).toarray()
    A[0,0] = 1-s
    A[-1,-1] = 1-s
    
    # production
    u2 = np.zeros((Nx,1)).ravel() 
    u2[0] = dt*kp 
    
    #---Euler method
    for k in range(0, Nt):
        #---dm/dt---
        m1 = kb*U[0, k]*f(M[0, k])-ks*M[0, k]
        m2 = kb*U[-1, k]*f(M[-1, k])-ks*M[-1, k]
        
        M[0, k+1] = M[0, k] + dt*m1
        M[-1, k+1] = M[-1, k] + dt*m2    
        #---du/dt    
        
        #degrdation
        u3 = -dt*kd*U[:, k]
        # dimer
        u4 = np.zeros((Nx,1)).ravel()
        u4[0] = -m1*dt
        u4[-1] = -m2*dt
        # diffusion
    
        U[:,k+1] = A.dot(U[:, k]) + u2 + u3 + u4
        
    #----polar diffenrence
    #diff_M1 = [(M[Nx-1,-1]- M[0,-1])/M[Nx-1,-1]]
    #U_mean = np.mean(U[:,-1])
   
    return [M[0,-1], M[-1,-1], U[0,-1], U[-1,-1]]


def async_multicore(main, paras, size, pool_n):
    pool = mp.Pool(processes = pool_n)    # Open multiprocessing pool
    result = []
    para_save = np.array([[]]).reshape(0, size)
    #do computation
    for para in paras:
        res = pool.apply_async(main, args=(para,))
        result.append(res)
        #para_reshape = para.reshape(1, size)
        #para_save = np.concatenate((para_save, para_reshape), axis=0)
        
    pool.close()
    pool.join()
    
    return result 
  
if __name__ == "__main__" :  
    
    name = "scan_k_data.npy"
    pool_n = 2 # number of processes in parallel
    para_size = 50 # number of grids to scan
    kd = np.logspace(-6, -2, para_size)
    ks = np.logspace(-4, -0, para_size)
    para = np.meshgrid(kd, ks)
    tmp = np.concatenate((np.array([para[0].flatten()]),
                          np.array([para[1].flatten()])), axis=0)
    paras = np.rot90(tmp)
    
    print('Opening {0} CPUs for simulation...'.format(pool_n))
    print('Preparing parameter sets...')
    print('Preparing simulation space...')
    print('Loading all simulations...')
    holder = partial(main,)
    start_time = time.time()    
    result = async_multicore(holder, paras, 2, pool_n)
    end_time = time.time()
    print('Finished all simulations with {0} sec'.format(end_time - start_time))
    print('Saving data...')
    data = [p.get() for p in result]
    data2 = [data, [paras]]   
    np.save(name,data2)
    print('Finishing program')