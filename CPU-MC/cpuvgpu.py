import numpy as np
import os
import matplotlib.pyplot as plt
import statsmodels.api as sm


# Load the data
avg_calc_int1 = np.array(np.load(os.path.join('CPU-MC/np_store_CPU', 'avg_calc_int_CPU.npy')))
std_list1 = np.array(np.load(os.path.join('CPU-MC/np_store_CPU', 'std_list_CPU.npy')))
history_count_list1 = np.array(np.load(os.path.join('CPU-MC/np_store_CPU', 'history_count_list_CPU.npy')))
batch_times1 = np.array(np.load(os.path.join('CPU-MC/np_store_CPU', 'batch_times_CPU.npy')))

avg_calc_int2 = np.array(np.load(os.path.join('GPU-MC/np_store_GPU', 'avg_calc_int_GPU.npy')))
std_list2 = np.array(np.load(os.path.join('GPU-MC/np_store_GPU', 'std_list_GPU.npy')))
history_count_list2 = np.array(np.load(os.path.join('GPU-MC/np_store_GPU', 'history_count_list_GPU.npy')))
batch_times2 = np.array(np.load(os.path.join('GPU-MC/np_store_GPU', 'batch_times_GPU.npy')))



x1 = history_count_list1
y1 = batch_times1
x2 = history_count_list2
y2 = batch_times2

X1 = sm.add_constant(x1)
model = sm.OLS(y1, X1).fit()
x1_fit = np.linspace(x1.min(), x1.max(), 500)
X1_fit = sm.add_constant(x1_fit)
y1_fit = model.predict(X1_fit)
X2 = sm.add_constant(x2)
model = sm.OLS(y2, X2).fit()
x2_fit = np.linspace(x2.min(), x2.max(), 500)
X2_fit = sm.add_constant(x2_fit)
y2_fit = model.predict(X2_fit)
plt.plot(history_count_list1, batch_times1, label=f"CPU")
plt.plot(history_count_list2,batch_times2, label=f"GPU")
plt.plot(x1_fit, y1_fit, 'r-', alpha=0.3, label=f"Fit CPU : {model.params[1]:.3e}*x + {model.params[0]:.3e}")
plt.plot(x2_fit, y2_fit, 'r-', alpha=0.3, label=f"Fit GPU: {model.params[1]:.3e}*x + {model.params[0]:.3e}")
plt.ylim(0, max(y1))
plt.xlabel('History Count')
plt.ylabel('Batch Time (s)')
plt.legend()
plt.show()