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

simdata= "simdata8"
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
ypos = 100
diff = 70
# As well define the x-position so that we can skip effects from the boundary
xpos1 = 50
xpos2 = xpixel-xpos1
delta = xpos2-xpos1
# Define the wavevector k
k = np.linspace(1,delta, delta) #Define the k-Vector
# Define the Fouriertransform in x-Axis
E_sum = np.zeros(delta)
E_sum_spatial = np.zeros(delta)
E_f_tot_spa= np.zeros(delta)
arrdiff=2
# the indices of Energy array used for averaging
for h in range(1200,1200+arrdiff):
    # get data
    E = np.loadtxt(simdata+"/e_mat"+str(h)+".dat").T[0]
    # Check if dimension fit bcs otherwise plot doesnt work 
    E = E.reshape((xpixel,ypixel))
    for j in range(ypos,ypos+diff):
        for i in range(xpos1,xpos2):
            E_sum[i-xpos1]+=E[i,j]

    #   Perform Fouriertransform and extract the Amplitude
    E_f_sum= np.fft.fft(E_sum)
    E_f_sum_amp = np.abs( np.sqrt(E_f_sum.real*E_f_sum.real+E_f_sum.imag*E_f_sum.imag) )
    # total Ef
    E_f_tot_spa +=E_f_sum_amp
E_f_tot_spa= E_f_tot_spa*2/(arrdiff*delta*diff)
# Calculate the fft in the x-direction
t_length=int(100)
y_pos=180
t_start= 1200
E_sum_temp = np.zeros((delta,t_length))
E_sum = np.zeros((delta,t_length))
for s in range(t_start,t_start+t_length):
    # get data
    E = np.loadtxt(simdata+"/e_mat"+str(s)+".dat").T[0]
    # Check if dimension fit bcs otherwise plot doesnt work 
    E = E.reshape((xpixel,ypixel))
    # Obtain the value
    for i in range(xpos1,xpos2):
        E_sum[i-xpos1,s-t_start]+=E[i,ypos]
    E_sum_temp+=E_sum
  
    

  
# Perform Fouriertransform and extract the Amplitude
E_f_amp_temp_ges = np.zeros(t_length)
for i in range(delta):
    E_f_sum_temp = np.fft.fft(E_sum_temp[i])
    E_f_amp_temp = 2/(len(E_f_sum_temp)*delta/diff*t_length) * np.abs( np.sqrt(E_f_sum_temp.real*E_f_sum_temp.real+E_f_sum_temp.imag*E_f_sum_temp.imag) )
    E_f_amp_temp_ges += E_f_amp_temp 


E_f_amp_temp_ges= E_f_amp_temp_ges/(delta*t_length)
f = np.linspace(0,t_length)/t_length*delta/2





fig = plt.figure()
#clear axis to plot graph in new chart
line2 = plt.loglog(f, E_f_amp_temp_ges[:len(f)],marker="o", markersize=1, markeredgecolor="red", markerfacecolor="green", label='temporal')
# take just half the indices bcs the fourier is symmetrical about the nyquist frequency
line1 = plt.loglog(k[:int(delta//2)], E_f_tot_spa[:int(delta//2)],marker="o", markersize=1, markeredgecolor="red", markerfacecolor="green", label='spatial')
plt.text(3,0.1,"Reynoldsnumber:" + str(reynold) +
             "\n Hyperviskosität:" + str(hypervisc))
    
# Make labels
plt.xlabel("$k$", fontsize=15)
plt.ylabel("$E$", fontsize=15)

plt.legend(loc= 'lower left')
plt.show()

