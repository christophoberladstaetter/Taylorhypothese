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
from scipy.optimize import curve_fit
from scipy.signal import blackman
from scipy.signal import hanning
import scipy as sci


# read in the griddata file as well as the counter
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
ypos = 80
diff = 4
# As well define the x-position so that we can skip effects from the boundary
xpos1 = 40
xpos2 = xpixel-xpos1
delta = xpos2-xpos1
fitdata_spa=20

# Define the wavevector k
k = np.linspace(1,delta, delta) #Define the k-Vector

# Define the Fouriertransform in x-Axis
E_sum = np.zeros(delta)
E_sum_spatial = np.zeros(delta)
E_f_tot_spa= np.zeros(delta)
arrdiff=1
t_start= 200

# the indices of Energy array used for averaging
#for h in range(t_start,t_start+arrdiff):
    # get data
E = np.loadtxt(simdata+"/e_mat"+str(t_start)+".dat").T[0]
# Check if dimension fit bcs otherwise plot doesnt work 
E = E.reshape((xpixel,ypixel))
for j in range(ypos,ypos+diff):
    for i in range(xpos1,xpos2):
        E_sum[i-xpos1]=+E[i,j]

# Create hanning window
w = sci.signal.windows.hann(int(delta))
#   Perform Fouriertransform and extract the Amplitude
E_res = sci.fft.fft(E_sum*w)
#E_f_sum= np.fft.fft(E_sum)
# Use the windowed fft
E_f_sum=E_res
E_f_sum_amp =  np.sqrt(E_f_sum.real*E_f_sum.real+E_f_sum.imag*E_f_sum.imag) 
# total Ef
E_f_tot_spa +=E_f_sum_amp
E_f_tot_spa= E_f_tot_spa*2/(arrdiff*delta)


# objective function
def objective(x,a,b):
    return a+b*x
# fit curve
popt, pcov = curve_fit(objective,np.log10(k[0:fitdata_spa]),np.log10(E_f_tot_spa[0:fitdata_spa]), p0=[0.477,-1.6], maxfev=55000)
# Extract the fitting parameters and get errors
fit_a_spa,fit_b_spa= popt
E_f_spa_fit = pow(k[0:fitdata_spa],fit_b_spa)* pow(fit_a_spa,10)
spa_err = np.sqrt(np.diag(pcov))




# Calculate the fft in the x-direction
t_length=1000
x_pos =110
x_diff =3
y_diff=3
y_pos=ypos+diff
t_length_fit=30


E_sum_temp = np.zeros(t_length)
for i in range(t_start,t_start+t_length):
    # get data
    E = np.loadtxt(simdata+"/e_mat"+str(i)+".dat").T[0]
    # Check if dimension fit bcs otherwise plot doesnt work 
    E = E.reshape((xpixel,ypixel))
    # Obtain the value
    for l in range(x_pos, x_pos+x_diff):
        for j in range(y_pos, y_pos+y_diff):
            E_sum_temp[i-t_start]+=E[l,j]
  
    

# Create hanning window
w2 = sci.signal.windows.hann(t_length)
#   Perform Fouriertransform and extract the Amplitude
E_f_sum_temp = sci.fft.fft(E_sum_temp*w2)
# Perform Fouriertransform and extract the Amplitude
#E_f_sum_temp = np.fft.fft(E_sum_temp)
E_f_amp_temp = 2/(len(E_f_sum_temp)*x_diff*y_diff/diff*t_length) * np.abs( np.sqrt(E_f_sum_temp.real*E_f_sum_temp.real+E_f_sum_temp.imag*E_f_sum_temp.imag) )

# Adjust the frequency to the wavenumber
f = np.linspace(1,t_length,t_length)/t_length*delta

popt, _ = curve_fit(objective,np.log10(f[1:t_length_fit]),np.log10(E_f_amp_temp[1:t_length_fit]), p0=[0.5,-1.6], maxfev=55000)
# Extract the fitting parameters
fit_a_temp,fit_b_temp= popt
E_f_temp_fit = pow(f[1:t_length_fit],fit_b_temp)*pow(fit_a_temp,10)
temp_err = np.sqrt(np.diag(pcov))

# Check if dimension fit
# get data
file = simdata+"/w2d"+str(5)+".dat"
data = np.loadtxt(file, delimiter= '  ').T[2]
data = data.reshape((xpixel,ypixel))
E = E.reshape((xpixel,ypixel))



# Define the kind of figure 
fig,(ax1, ax2) = plt.subplots(1,2, figsize=(14,7))
# take just half the indices bcs the fourier is symmetrical about the nyquist frequency
line1 = ax1.loglog(f[0:len(f)], E_f_amp_temp[0:len(f)],marker="o", markersize=1,color='red', markeredgecolor="blue", markerfacecolor="green", label='temporal')
line2 = ax1.loglog(k[:int(delta//2)], E_f_tot_spa[:int(delta//2)],marker="o",color='black', markersize=1, markeredgecolor="red", markerfacecolor="green", label='spatial')
line5 = ax1.loglog(k[0:fitdata_spa], E_f_spa_fit[0:fitdata_spa]*pow(np.e,np.log(E_f_tot_spa[0])-np.log(E_f_spa_fit[0])+1),color='black',linestyle='dashed', label='spatial-slope')
line6 = ax1.loglog(f[1:t_length_fit], E_f_temp_fit*pow(np.e,np.log(E_f_amp_temp[1])-np.log(E_f_temp_fit[0])),color='red',linestyle='dashed', label='temporal-slope')

ax1.text(5,3,"Reynoldsnumber:" + str(reynold) +
            "\n Hyperviskosität:" + str(hypervisc))
ax1.text(1,0.001,"Fitmodel k * x (+ d ): \nslope spatial:" +"{:.3f}".format(fit_b_spa)+ "+/-" +"{:.3f}".format(spa_err[1]) +"\nslope temporal:"+"{:.3f}".format(fit_b_temp) +"+/-"+"{:.3f}".format(temp_err[1]))
# Make labels
ax1.set_xlabel("$k$ and $\omega$", fontsize=20)
ax1.set_ylabel("$E$", fontsize=20)
ax1.grid()
ax1.tick_params(axis='both', which='major', labelsize=15)
ax1.tick_params(axis='both', which='minor', labelsize=15)
ax1.set_xlim([1, 400])
ax1.set_ylim([0.0000000000000001, 10])
#Plot the voritcity
im = ax2.imshow(data, cmap=cmap)
ax2.tick_params(axis='both', which='major', labelsize=15)
ax2.tick_params(axis='both', which='minor', labelsize=15)
# Set limits and plot a grid, axis ticks
#ax1.set_xlim([0.9, 100])
#ax1.set_ylim([0.00000001, 10])
# Plot a title
ax2.set_title("Vorticity")

# mark where the FFT is calculated
line4 = ax2.plot([y_pos,y_pos+y_diff,y_pos+y_diff,y_pos,y_pos], [x_pos,x_pos,x_pos+x_diff,x_pos+x_diff,x_pos ], label= 'temporal', color= 'red')
line3 = ax2.plot([ypos,ypos+diff,ypos+diff,ypos,ypos],[xpos1,xpos1,xpos2,xpos2,xpos1], label='spatial', color='black')

ax1.legend(loc='lower right')
plt.show()
