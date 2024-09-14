import numpy as np
import os
import matplotlib.pyplot as plt

# Load the data
time = np.load(os.path.join('np_store', 'batch_times.npy'))
std = np.load(os.path.join('np_store', 'std_list.npy'))

# Flatten the arrays (only if necessary)
time = np.array(time).flatten()
std = np.array(std).flatten()

# Assuming you want to plot all points in the arrays:
x = std  # x-values for plotting
y = time  # y-values for plotting

# Plotting
plt.figure(1)
plt.suptitle(" ")
plt.xlabel("Standard Deviation")
plt.xscale('linear')
plt.ylabel("Batch Times (s)")
plt.scatter(x, y, s=0.9, alpha=1, color='b', label="Data Points")
plt.legend()
plt.show()
