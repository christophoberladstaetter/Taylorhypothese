# (c) Christoph Oberladstätter, UIBK
# This file is intended to be used for making easy plots
# for the thesis and the presentation


import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# from scipy.fft import fft, ifft
import math
import psutil
import statistics
from functools import partial
import matplotlib.animation as animation
from IPython import display
# read in the griddata file as well as the counter
#
simdata= "simdata6"
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
cmap_obs.set_bad(alpha=0.8)
# make all zero values 'bad' values
obs_mask[obs_mask<0.5] = np.nan

# get a copy of the gray color map
cmap_obs_r = copy.copy(plt.cm.get_cmap('gray'))
# set how the colormap handles 'bad' values
cmap_obs_r.set_bad(alpha=0.8)

# make all zero values 'bad' values
obs_mask_r[obs_mask_r>0.5] = np.nan

# Define position where to carry out FFT and over how many pixels
ypos = 230
diff = 5
# As well define the x-position so that we can skip effects from the boundary
xpos1 = 80
xpos2 = xpixel-xpos1
delta = xpos2-xpos1
# Define the wavevector k
k = np.linspace(1,delta, delta) #Define the k-Vector

 # get data
E = np.loadtxt(simdata+"/e_mat0.dat").T[0]
# Check if dimension fit bcs otherwise plot doesnt work 
E = E.reshape((xpixel,ypixel))
# Define the Fouriertransform in x-Axis
E_sum = np.zeros(delta)

for j in range(ypos,ypos+diff):
  for i in range(xpos1,xpos2):
    E_sum[i-xpos1]=+E[i,j]

# Perform Fouriertransform and extract the Amplitude
E_f_sum = np.fft.fft(E_sum)
E_f_amp = 2/(len(E_sum)*(diff+1)) * np.abs( np.sqrt(E_f_sum.real*E_f_sum.real+E_f_sum.imag*E_f_sum.imag) )
file = simdata+"/w2d0.dat"

# Check if dimension fit 
data = np.loadtxt(file, delimiter= '  ').T[2]
data = data.reshape((xpixel,ypixel))

# Define the kind of figure 
fig,(ax2, ax3) = plt.subplots(1,2, figsize=(7,7))
im = ax2.imshow(data, cmap=cmap)

file = open(simdata+"/w2d"+str(500)+".dat",'r')
simtime = dt*itstp*i

# Check if dimension fit
data = np.loadtxt(file, delimiter= '  ').T[2]
data = data.reshape((xpixel,ypixel))


#clear axis to plot graph in new chart
ax2.cla()

#Plot the voritcity
im = ax2.imshow(data, cmap=cmap)
im2 = ax3.imshow(obs, cmap=cmap_obs) 
ax2.tick_params(axis='both', which='major', labelsize=15)
ax2.tick_params(axis='both', which='minor', labelsize=15)
# Plot a title
ax2.set_title("Vorticity")
ax3.set_title("Obstacle")


fig.colorbar( im, ax=ax2) 
plt.show()

# saving to m4 using ffmpeg writer
#writervideo = animation.FFMpegWriter(fps=30)
#ani.save('try.mp4', writer=writervideo)
#plt.close()