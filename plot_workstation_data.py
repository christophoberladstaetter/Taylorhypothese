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
#
simdata= "simdata5"
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
ypos=100
diff=20

# Define the wavevector k
k = np.linspace(1,xpixel,xpixel) #Define the k-Vector
k = np.linspace(1,xpixel,xpixel) #Define the k-Vector
 # get data
E = np.loadtxt(simdata+"/e_mat0.dat").T[0]
# Check if dimension fit bcs otherwise plot doesnt work 
E = E.reshape((xpixel,ypixel))
# Define the Fouriertransform in x-Axis
E_sum = np.zeros(xpixel)

for j in range(ypos,ypos+diff):
  for i in range(0,xpixel):
    E_sum[i]=+E[i,j]

# Perform Fouriertransform and extract the Amplitude
E_f_sum = np.fft.fft(E_sum)
E_f_amp = 2/(len(E_sum)*(diff+1)) * np.abs( np.sqrt(E_f_sum.real*E_f_sum.real+E_f_sum.imag*E_f_sum.imag) )
file = simdata+"/w2d0.dat"

# Check if dimension fit 
data = np.loadtxt(file, delimiter= '  ').T[2]
data = data.reshape((xpixel,ypixel))
fig,(ax1, ax2) = plt.subplots(2,1, figsize=(7,15))
im = ax2.imshow(data, cmap=cmap)
def my_function(i, E_f_amp_sum):

    # get data
    file = simdata+"/w2d"+str(i)+".dat"
    E = np.loadtxt(simdata+"/e_mat"+str(i)+".dat").T[0]
    simtime = dt*itstp*i
    # Check if dimension fit
    data = np.loadtxt(file, delimiter= '  ').T[2]
    data = data.reshape((xpixel,ypixel))
    E = E.reshape((xpixel,ypixel))

    # Define the Fouriertransform in x-Axis
    E_sum = np.zeros(xpixel)
    for j in range(ypos,ypos+diff):
        for i in range(0,xpixel):
            E_sum[i]=+E[i,j]
    
    # Subtract the mean energy bcs I dont know???
    #meanE = statistics.mean(E_sum)  
    #for i in range(0,xpixel):
    #  E_sum[i]=E_sum[i]-meanE

    # Perform Fouriertransform
    E_f_sum = np.fft.fft(E_sum)
    E_f_amp = 2/(len(E_sum)*(diff+1)) * np.abs( np.sqrt(E_f_sum.real*E_f_sum.real+E_f_sum.imag*E_f_sum.imag) )
    
    # Update the sum--> not working for now
    E_f_amp_sum= E_f_amp_sum+E_f_amp

    #clear axis to plot graph in new chart
    ax1.cla()
    ax2.cla()

    # take just half the indices bcs the fourier is symmetrical about the nyquist frequency
    ax1.loglog(k[:xpixel//2], E_f_amp_sum[:xpixel//2],marker="o", markersize=1, markeredgecolor="red", markerfacecolor="green")
    ax1.text(3,0.1,"Reynoldsnumber:" + str(reynold) +
             "\n Hyperviskosität:" + str(hypervisc)+
             "\n Simulationszeit:" + str(simtime))
    
    # Make labels
    ax1.set_xlabel("log(k)", fontsize=20)
    ax1.set_ylabel("log(E)", fontsize=20)
   
    
    im = ax2.imshow(data, cmap=cmap)
    # Plot a title
    ax2.set_title("Vorticity")

    # mark where the FFT is calculated
    ax2.axvline(x=ypos, color = 'grey', alpha=0.7)
    ax2.axvline(x=ypos+diff, color = 'grey', alpha=0.7)


E_f_amp_sum = E_f_amp 
# animate
ani = FuncAnimation(fig, func=partial(my_function,E_f_amp_sum=E_f_amp_sum), interval=100)
fig.colorbar( im, ax=ax2) 
plt.show()



