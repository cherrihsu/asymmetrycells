# -*- coding: utf-8 -*-
"""
This program reads the npy datafile generated from scan_k.py and makes plots.
"""
import numpy as np
import matplotlib.pyplot as plt

DIR = "./"
data=np.load(DIR + "scan_k_data.npy", allow_pickle=True)


#%% data processing
D = data[0]
kskd = data[1][0]
kd = kskd[:, 0]
ks = kskd[:, 1]

MU = np.array(D)
M_H = MU[:,0]
M_T = MU[:,1]
U_H = MU[:, 2]
U_T = MU[:, 3] 
U = U_H/U_T
M = M_H/M_T
UM = (U_H+M_H)/(U_T + M_T) 
# scanning ks and kd in the same region
# square matrix
m_size = int(len(D)**0.5)

U = np.reshape(U, (m_size,m_size))
M = np.reshape(M, (m_size,m_size))
MH = np.reshape(M_H, (m_size,m_size))
MT = np.reshape(M_T, (m_size,m_size))
UM = np.reshape(UM, (m_size,m_size))
T = np.reshape(kd, (m_size,m_size))
X = np.reshape(ks, (m_size,m_size))

#%%
# plot U
fig, ax = plt.subplots()
U_min, U_max = -np.abs(U).max(), np.abs(U).max()
c = ax.pcolormesh(T, X, U, cmap = 'Blues')
fig.set_figheight(6)
fig.set_figwidth(8)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#ax.set_title('reaction -defusssion', size = 14), vmin=U_min, vmax=U_max
# set the limits of the plot to the limits of the data
ax.axis([T.min(), T.max(),X.min(), X.max()], size = 18)
ax.set_xlabel('k$_d$ (1/s)', size = 18)
ax.set_ylabel('k$_s$ (1/s)', size = 18)
ax.set_xscale("log")
ax.set_yscale("log")
fig.colorbar(c, ax=ax)
name = 'U ratio result'
filename = DIR + name + '.pdf'
fig.savefig(filename) # ,dpi= 1000)    


# plot M
fig, ax = plt.subplots()
U_min, U_max = -np.abs(M).max(), np.abs(M).max()
c = ax.pcolormesh(T, X, M, cmap = 'Blues')
fig.set_figheight(6)
fig.set_figwidth(8)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#ax.set_title('reaction -defusssion', size = 14), vmin=U_min, vmax=U_max
# set the limits of the plot to the limits of the data
ax.axis([T.min(), T.max(),X.min(), X.max()], size = 18)
ax.set_xlabel('k$_d$ (1/s)', size = 18)
ax.set_ylabel('k$_s$ (1/s)', size = 18)
ax.set_xscale("log")
ax.set_yscale("log")
fig.colorbar(c, ax=ax)
name = 'M ratio result'
filename = DIR + name + '.pdf'
fig.savefig(filename) # ,dpi= 1000)  
#plot U+M
fig, ax = plt.subplots()
U_min, U_max = -np.abs(UM).max(), np.abs(UM).max()
c = ax.pcolormesh(T, X, UM, cmap = 'Blues')
fig.set_figheight(6)
fig.set_figwidth(8)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#ax.set_title('reaction -defusssion', size = 14), vmin=U_min, vmax=U_max
# set the limits of the plot to the limits of the data
ax.axis([T.min(), T.max(),X.min(), X.max()])        
ax.set_xlabel('k$_d$ (1/s)', size = 18)
ax.set_ylabel('k$_s$ (1/s)', size = 18)
ax.set_xscale("log")
ax.set_yscale("log")
fig.colorbar(c, ax=ax) 
name = 'U + M ratio result'
filename = DIR + name + '.pdf'
fig.savefig(filename) # ,dpi= 1000) 

