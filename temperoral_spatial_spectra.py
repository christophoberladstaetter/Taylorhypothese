# (c) UIBK, Christoph Oberladstätter, Bachelorarbeit
# Erstellt 27.05.2022

#This script is used to compare the spatial and temporal spectras of the energy

# This script is for loading all the simulation data from the UIBK
# workstation and then visualize the simulation data
# (C) Christoph Oberladstätter, 2022, UIBK
# Hintergrundenergie abziehen also die Hintergrundgeschwindigkeit

import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# from scipy.fft import fft, ifft
import math
import psutil
import statistics
from functools import partial
# read in the griddata file as well as the counter

simdata= "simdata4"
grid= np.loadtxt(simdata+"/info.dat",delimiter = '  ') #delimiter are two blanks
xpixel = grid[0]
ypixel = grid[1]    
xpixel = int(xpixel)             # number of xpixel ---remove two again
ypixel = int(ypixel)            # number of ypixel
counter= int(grid[2])             # .dat counter
reynold = float(grid[3])          # reynoldsnumber
hypervisc = float(grid[4])        # hypervis
dt = float(grid[5])              # time step
itstp = float(grid[6])           # steps in between output

obs = np.loadtxt(simdata+"/obstacle.dat", delimiter = '  ').T[2] #transpose and 2nd 
obs = obs.reshape((xpixel,ypixel))
cmap = 'seismic'

# inverted obstacle for first subplot
obs_mask = copy.copy(obs)
# inverted obstacle for secons subplot
obs_mask_r = copy.copy(obs)

# get a copy of the gray color map
cmap_obs = copy.copy(plt.cm.get_cmap('gray_r'))
# set how the colormap handles 'bad' values
cmap_obs.set_bad(alpha=0)
# make all zero values 'bad' values
obs_mask[obs_mask<0.5] = np.nan

# get a copy of the gray color map
cmap_obs_r = copy.copy(plt.cm.get_cmap('gray'))
# set how the colormap handles 'bad' values
cmap_obs_r.set_bad(alpha=0)

# make all zero values 'bad' values
obs_mask_r[obs_mask_r>0.5] = np.nan

# Define position where to carry out FFT and over how many pixels
ypos=109
diff=15

# Define the wavevector k
k = np.linspace(1,xpixel,xpixel) #Define the k-Vector


# Calculate the fft in the x-direction
E_sum_spatial = np.zeros(xpixel)

for i in range(600,601):
    # get data
    E = np.loadtxt(simdata+"/e_mat"+str(i)+".dat").T[0]
    # Check if dimension fit bcs otherwise plot doesnt work 
    E = E.reshape((xpixel,ypixel))
    # Define the Fouriertransform in x-Axis
   
    for j in range(ypos,ypos+diff):
        for i in range(0,xpixel):
            E_sum_spatial[i]=+E[i,j]

    # Perform Fouriertransform and extract the Amplitude
    E_f_sum_spatial = np.fft.fft(E_sum_spatial)
    E_f_amp_spatial = 2/(len(E_sum_spatial)*(diff+1)) * np.abs( np.sqrt(E_f_sum_spatial.real*E_f_sum_spatial.real+E_f_sum_spatial.imag*E_f_sum_spatial.imag) )

# Calculate the fft in the x-direction
t_length=500
x_pos =50
y_pos=109
t_start= 300
E_sum_temp = np.zeros(t_length)
for i in range(t_start,t_start+t_length):
    # get data
    E = np.loadtxt(simdata+"/e_mat"+str(i)+".dat").T[0]
    # Check if dimension fit bcs otherwise plot doesnt work 
    E = E.reshape((xpixel,ypixel))
    # Obtain the value
    E_sum_temp[i-t_start]=E[x_pos,y_pos]
  
    

  
# Perform Fouriertransform and extract the Amplitude
E_f_sum_temp = np.fft.fft(E_sum_temp)
E_f_amp_temp = 2/(len(E_f_sum_temp)*(diff+1)) * np.abs( np.sqrt(E_f_sum_temp.real*E_f_sum_temp.real+E_f_sum_temp.imag*E_f_sum_temp.imag) )
f=np.linspace(0,t_length)/t_length*max(k)




fig = plt.figure()
#clear axis to plot graph in new chart
line2 = plt.loglog(f, E_f_amp_temp[:len(f)],marker="o", markersize=1, markeredgecolor="red", markerfacecolor="green", label='temporal')
# take just half the indices bcs the fourier is symmetrical about the nyquist frequency
line1 = plt.loglog(k[:xpixel//2], E_f_amp_spatial[:xpixel//2],marker="o", markersize=1, markeredgecolor="red", markerfacecolor="green", label='spatial')
plt.text(3,0.1,"Reynoldsnumber:" + str(reynold) +
             "\n Hyperviskosität:" + str(hypervisc))
    
# Make labels
plt.xlabel("log(k)", fontsize=20)
plt.ylabel("log(E)", fontsize=20)
plt.legend(loc= 'lower left')
plt.show()

