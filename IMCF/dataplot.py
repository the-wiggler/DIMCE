import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv("results.csv")
# data2 = pd.read_csv("resultsSERIAL.csv")
print(data.columns)
print(data.head())


# batch
# history
# calc_int
# stddev
# batch_time

# Plot stddev
plt.figure(1)
plt.scatter(data['history_count'], data['variance'], s=1, color='r')
plt.title(r"Integration of $\int_{0}^{1} \frac{\ln(1+x)}{x}\, dx$ ~ History Count v. Variance of Integration", fontsize=14)
plt.xlabel("History Count")
plt.ylabel("Variance")

plt.figure(3)
plt.scatter(data['history_count'], data['stdv'], s=1, color='r')
plt.title(r"Integration of $\int_{0}^{1} \frac{\ln(1+x)}{x}\, dx$ ~ History Count v. Standard Deviation of Integration", fontsize=14)
plt.xlabel("History Count")
plt.ylabel("Standard Deviation")



# # Plot calc_int
# plt.figure(2)
# plt.scatter(data['history'], data['calc_int'], s=1, color='r')
# plt.title(r"Integration of $\int_{0}^{1} \frac{\ln(1+x)}{x}\, dx$ ~ History Count v. Calculated Integration", fontsize=14)
# plt.xlabel("History Count")
# plt.ylabel("Calculated Integration")

# # # Plot batch_time
# plt.figure(3)
# plt.scatter(data['history'], data['batch_time'], s=1, color='r')
# plt.scatter(data2['history'], data2['batch_time'], s=1, color='b')
# plt.title(r"Integration of $\int_{0}^{1} \frac{\ln(1+x)}{x}\, dx$ ~ History Count v. Batch Time", fontsize=14)
# plt.gca().set_title('*Using OpenMP Parallelization', loc='right', pad=20)
# plt.xlabel("History Count")
# plt.ylabel("Batch Time")

# plt.figure(4)
# plt.scatter(data['history'], data['percent_err'], s=1, color='r')
# plt.title(r"Integration of $\int_{0}^{1} \frac{\ln(1+x)}{x}\, dx$ ~ History Count v. Percent Error from Known Value", fontsize=14)
# plt.gca().set_title('*Using OpenMP Parallelization', loc='right', pad=20)
# plt.xlabel("History Count")
# plt.ylabel("Percent Error")

# data = pd.read_csv("varbatchOMP.csv")
# data2 = pd.read_csv("varbatchSERIAL.csv")

# plt.figure()
# plt.plot(data['batches'], data[' compute_time'], color='r', marker='o', label='OMP Parallelization Method')
# plt.plot(data2['batches'], data2[' compute_time'], color='b', marker='o', label='Serial Computation')
# plt.title(r"Integration Time of $\int_{0}^{1} \frac{\ln(1+x)}{x}\, dx$ With Varying Batches", fontsize=14)
# plt.gca().set_title('Histories Per Batch: 10000000', loc='right', pad=20)
# plt.xlabel("Batches")
# plt.ylabel("Compute Time")
# plt.legend()  # Add a legend to differentiate lines
# plt.grid(True)  

plt.show()