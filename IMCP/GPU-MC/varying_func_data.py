import numpy as np
import os
import matplotlib.pyplot as plt

# Load the data
avg_calc_int1 = np.array(np.load(os.path.join('np_store_GPU', 'avg_calc_int_GPU.npy')))
std_list1 = np.array(np.load(os.path.join('np_store_GPU', 'std_list_GPU.npy')))
history_count_list1 = np.array(np.load(os.path.join('np_store_GPU', 'history_count_list_GPU.npy')))
batch_times1 = np.array(np.load(os.path.join('np_store_GPU', 'batch_times_GPU.npy')))

np.savetxt("calc_int.csv", avg_calc_int1, delimiter = ",")
np.savetxt("stdv_list.csv", std_list1, delimiter = ",")
np.savetxt("batch_times.csv", batch_times1, delimiter = ",")