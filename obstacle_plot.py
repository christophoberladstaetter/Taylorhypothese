import numpy as np
import matplotlib.pyplot as plt

obs = np.loadtxt("obstacle.dat", delimiter = '  ').T[2] # T is for transposing the matrix and [2] gives back 2nd array
obs = obs.reshape((64,256))

plt.imshow(obs)
plt.colorbar()
plt.tight_layout()
plt.show()
