# This is a script which continously plots when simulating
# For that just run the script in the terminal
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# from scipy.fft import fft, ifft
import math
import psutil

#* decide what to plot here ---------------------------
# velocity
#file = 'v2d.dat'

# vorticity-omega
file = 'w2d.dat'

# flow potential
# file = 'w2d.dat'
# ----------------------------------------------------


# read in the griddata file
grid= np.loadtxt("grid_counter.dat",delimiter = '  ') #delimiter are two blanks
xpixel = grid[0]
ypixel = grid[1]
xpixel = int(xpixel)
ypixel = int(ypixel)
counter= int(grid[2])




if file == 'w2d.dat': cmap = 'seismic'
else: cmap = 'Reds'

obs = np.loadtxt("obstacle.dat", delimiter = '  ').T[2] #transpose and 2nd 
obs = obs.reshape((xpixel,ypixel))
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
k = np.linspace(1,xpixel,xpixel) #Define the k-Vector
 # get data
E = np.loadtxt("energie_matrix.dat")
# Check if dimension fit bcs otherwise plot doesnt work 

E = E.reshape((xpixel,ypixel))
# Define the Fouriertransform in x-Axis
E_sum = np.zeros(xpixel)
# Average over a bigger y-Interval
diff=5
for j in range(100,100+diff):
  for i in range(0,xpixel):
    E_sum[i]=+E[i,j]

# Perform Fouriertransform
E_f_sum = np.fft.rfft(E_sum)
E_f_amp = 2/(len(E_sum)*(diff+1)) * np.abs( np.sqrt(E_f_sum.real*E_f_sum.real+E_f_sum.imag*E_f_sum.imag) )


# function to update data
def my_function(i):
  # get data
  E = np.loadtxt("energie_matrix.dat")
  data = np.loadtxt(file, delimiter= '  ').T[2]
  data = data.reshape((xpixel,ypixel))
  if(len(E)==xpixel*ypixel):
    E = E.reshape((xpixel,ypixel))

    # Define the Fouriertransform in x-Axis
    E_sum = np.zeros(xpixel)
    # Average over a bigger y-Interval
    for j in range(100,100+diff):
      for i in range(0,xpixel):
        E_sum[i]=+E[i,j]

    # Perform Fouriertransform
    E_f_sum = np.fft.rfft(E_sum)
    E_f_amp = 2/(len(E_sum)*(diff+1)*xpixel*ypixel) * np.abs( np.sqrt(E_f_sum.real*E_f_sum.real+E_f_sum.imag*E_f_sum.imag) )
    k = np.linspace(1,xpixel,xpixel) #Define the k-Vector

    #clear axis to plot graph in new chart
    ax1.cla()
    ax2.cla()


    #plot data1
    # take just half the indices bcs the fourier is symmetrical about the nyquist frequency
    ax1.plot(np.log10(k[0:int(xpixel/2)+1]), np.log10(E_f_amp[0:int(xpixel/2)+1]),marker="o", markersize=1, markeredgecolor="red", markerfacecolor="green")
    #ax.scatter(len(E_f_amp)-1, np.log(E_f_amp[1]))


  
  ax2.imshow(data, cmap=cmap)
# start with zeros arrays


# define and adjust figure
fig,(ax1, ax2) = plt.subplots(2,1)
ax1.scatter(np.log10(k[0:int(xpixel/2)+1]), np.log10(E_f_amp[0:int(xpixel/2)+1]))
#ax1.set_facecolor('#DEDEDE')
#ax2.set_facecolor('#DEDEDE')

# animate
ani = FuncAnimation(fig, my_function, interval=1000)
plt.show()

