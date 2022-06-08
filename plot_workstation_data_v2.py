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
import matplotlib.animation as animation
from IPython import display
# read in the griddata file as well as the counter
#
simdata= "simdata1"
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
ypos = 160
diff = 5
# As well define the x-position so that we can skip effects from the boundary
xpos1 = 50
xpos2 = xpixel-xpos1
delta = xpos2-xpos1
# Define the wavevector k
k = np.linspace(1,delta, delta) #Define the k-Vector

file = simdata+"/w2d"+str(0)+".dat"



# Check if dimension fit
data = np.loadtxt(file, delimiter= '  ').T[2]
data = data.reshape((xpixel,ypixel))


# Define the kind of figure 
fig,ax2 = plt.subplots(1,1, figsize=(14,7))
#Plot the voritcity
im = ax2.imshow(data, cmap=cmap)
def my_function(i):
    try:
        file = open(simdata+"/w2d"+str(i)+".dat",'r')
    except IOError:
        return False
    # get data
    file = simdata+"/w2d"+str(i)+".dat"
   
    simtime = dt*itstp*i
    ax2.cla()
    # Check if dimension fit
    data = np.loadtxt(file, delimiter= '  ').T[2]
    data = data.reshape((xpixel,ypixel))
    ax2.text(110,30,"Reynoldsnumber:" + str(reynold) +
             "\n Hyperviskosität:" + str(hypervisc)+
             "\n Simulationszeit:" + "{:.3f}".format(simtime))
    
    #Plot the voritcity
    im = ax2.imshow(data, cmap=cmap)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax2.tick_params(axis='both', which='minor', labelsize=15)
    # Plot a title
    ax2.set_title("Vorticity")



# animate and save the animation as a video
fig.colorbar(im, ax=ax2) 

#ani = FuncAnimation(fig, my_function, interval=2)
ani = FuncAnimation(fig, my_function, interval=5, blit=False,repeat=True,save_count = 1467)
plt.show()
plt.close()
# saving to m4 using ffmpeg writer
plt.close()
writervideo = animation.FFMpegWriter(fps=30)
ani.save('praesentation_vid_v5.mp4', writer=writervideo)
