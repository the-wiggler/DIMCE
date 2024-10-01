import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv("results.csv")

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Plot variance
axes[0, 0].scatter(data['history_count'], data['variance'], s=0.4, color='b')
axes[0, 0].set_title("History Count v. Variance", fontsize=10)
axes[0, 0].set_xlabel("History Count")
axes[0, 0].set_ylabel("Variance")
axes[0, 0].set_xscale('log')
axes[0, 0].set_yscale('log')

# Plot standard deviation
axes[0, 1].scatter(data['history_count'], data['stdv'], s=0.4, color='r')
axes[0, 1].set_title("History Count v. Standard Deviation", fontsize=10)
axes[0, 1].set_xlabel("History Count")
axes[0, 1].set_ylabel("Standard Deviation")
axes[0, 1].set_xscale('log')
axes[0, 1].set_yscale('log')

# Plot relative error
axes[0, 2].scatter(data['history_count'], data['rel_err'], s=1, color='r')
axes[0, 2].set_title("History Count v. Relative Error", fontsize=10)
axes[0, 2].set_xlabel("History Count")
axes[0, 2].set_ylabel("Relative Error")
axes[0, 2].set_xscale('log')
axes[0, 2].set_yscale('log')

# Plot batch time
axes[1, 0].scatter(data['history_count'], data['batch_time'], s=1, color='g')
axes[1, 0].set_title("History Count v. Batch Time", fontsize=10)
axes[1, 0].set_xlabel("History Count")
axes[1, 0].set_ylabel("Batch Time")
axes[1, 0].set_xscale('log')
axes[1, 0].set_yscale('log')

# Plot batch time vs relative error
axes[1, 1].scatter(data['batch_time'], data['rel_err'], s=1, color='g')
axes[1, 1].set_title("Batch Time v. Relative Error", fontsize=10)
axes[1, 1].set_xlabel("Batch Time")
axes[1, 1].set_ylabel("Relative Error")
axes[1, 1].set_xscale('log')
axes[1, 1].set_yscale('log')

# Plot calculated integral vs histories
axes[1, 2].scatter(data['history_count'], data['IMC_val'], s=1, color='g')
axes[1, 2].set_title("Histories v. Calculated Integral", fontsize=10)
axes[1, 2].set_xlabel("Histories")
axes[1, 2].set_ylabel("IMC Value")

plt.tight_layout()

plt.show()
