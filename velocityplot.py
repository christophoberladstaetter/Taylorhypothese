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


x=np.linspace(0,100, num=2000)
y= 3+0.1*np.sin(310*x)+0.22*np.sin(8*x)+0.21*np.sin(3*x)+0.21*np.sin(0.01*x)

fig, ax = plt.subplots()

ax.plot(x,y)
plt.rcParams.update({'font.size': 22})
ax.set(xlabel='$t$ / (s)', ylabel='$v$ / (m/s)')
ax.grid()
ax.set_ylim([0,4])

ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(axis='both', which='minor', labelsize=15)
plt.show()