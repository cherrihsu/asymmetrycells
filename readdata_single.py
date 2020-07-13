# -*- coding: utf-8 -*-
"""
This program reads the npy datafile generated from single.py and makes plots.
"""
import numpy as np
import matplotlib.pyplot as plt

DIR = "./"
data=np.load(DIR + "single_data.npy", allow_pickle=True)

#%%
tspan_n = data[0]  # s
xspan = data[1] 
M = data[3]
U = data[2] 
tspan_min = tspan_n / 60 # min
#nm = 0.06022 # 1 nm = 0.06022 copies/grid
#%% Series drawing
# U and M finial state
name = 'U and M finial state'
filename = DIR + name + '.pdf'
plt.figure(figsize=(8.5,6), linewidth = 1.5)
plt.plot(xspan, M[:,-1])
plt.plot(xspan, U[:, -1])
plt.plot(xspan, ( M[:,-1] + U[:, -1]))
plt.xlabel('Length ($\mu$m)', fontsize = 20)
plt.ylabel('Concentration (nM)', fontsize = 20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(['m', 'u', 'total'], fontsize = 18)
plt.savefig(filename) # dpi= 1500)

#%% total
name = 'total'
Total = (M[:,-1] + U[:, -1])
Total_n = (M[-1,-1] + U[-1, -1])

filename = DIR + name + '.pdf'
plt.figure(figsize=(8.5,6), linewidth = 1.5)
plt.plot(xspan, (Total/Total_n))
plt.xlabel('Length ($\mu$m)', fontsize = 20)
plt.ylabel('Concentration (arb. unit)', fontsize = 20)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.savefig(filename) # dpi= 1500)


#%%
#　U finial state 
name = 'U finial state'
filename = DIR + name + '.pdf'
plt.figure(figsize=(8.5,6), linewidth = 1.5)
plt.plot(xspan, U[:, -1])
plt.xlabel('length ($\mu$m)', fontsize = 20)
plt.ylabel('Concentration (nM)',  fontsize = 20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig(filename) # dpi= 1000)

#　M finial state 
name = 'M finial state' 
filename = DIR + name + '.pdf'
plt.figure(figsize=(8.5,6), linewidth = 1.5)
plt.plot(xspan, M[:, -1])
plt.xlabel('Length ($\mu$m)', fontsize = 20)
plt.ylabel('Concentration (nM)', fontsize = 20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig(filename) # dpi= 1000)

#%%
# U time series
name = 'U time series' 
filename = DIR + name + '.pdf'
plt.figure(figsize=(8.5,6), linewidth = 1.5)
plt.plot(tspan_min, U[0, :])
plt.plot(tspan_min, U[-1,:])
plt.legend(['head', 'tail'], fontsize = 18)
plt.xlabel('Time (Minutes)', fontsize = 20)
plt.ylabel('concentration (nM)', fontsize = 20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig(filename) # dpi= 1000)

#%%
# M time series
name = 'M time series' 
filename = DIR + name + '.pdf'
plt.figure(figsize=(8.5,6), linewidth = 1.5)
plt.plot(tspan_min, M[0, :])
plt.plot(tspan_min, M[-1,:])
plt.xlabel('Time (Minutes)', fontsize = 20)
plt.ylabel('concentration (nM)', fontsize = 20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(['head', 'tail'], fontsize = 18)
plt.savefig(filename) # dpi= 1000)

#%% heatmap
x0 = 0
xL = 2
Nx = 10
xspan_1 = np.linspace(x0, xL, Nx+1)
tspan_min = tspan_n/60
T, X = np.meshgrid(tspan_min, xspan_1)

# U heat map
fig, ax = plt.subplots()
fig.set_figheight(6)
fig.set_figwidth(8.5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
U_min, U_max = -np.abs(U).max(), np.abs(U).max()
c = ax.pcolormesh(T, X, U)
#ax.set_title('reaction -defusssion', size = 14), vmin=U_min, vmax=U_max
#set the limits of the plot to the limits of the data
ax.axis([T.min(), T.max(),X.min(), X.max()])
ax.set_xlabel('Time (Minutes)', size = 20)
ax.set_ylabel('Space ($\mu$m)', size = 20)
fig.colorbar(c, ax=ax)
name = 'heat_map_U'
filename = DIR + 'heat_map_U.pdf'
fig.savefig(filename) # dpi= 1000)


# heat_map_M
fig, ax = plt.subplots()
fig.set_figheight(6)
fig.set_figwidth(8.5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
M_min, M_max = -np.abs(M).max(), np.abs(M).max()
c = ax.pcolormesh(T, X, M)
#ax.set_title('reaction -defusssion', size = 14), vmin=U_min, vmax=U_max
# set the limits of the plot to the limits of the data
ax.axis([ T.min(), T.max(),X.min(), X.max()])
ax.set_xlabel('Time (Minutes)', size = 20)
ax.set_ylabel('Space ($\mu$m)', size = 20)
fig.colorbar(c, ax=ax)
name = 'heat_map_M'
filename = DIR + name + '.pdf'
fig.savefig(filename) # dpi= 1000)    

 
