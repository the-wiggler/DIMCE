import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("results.csv")
# data2 = pd.read_csv("results2.csv")

# batch
# history
# calc_int
# stddev
# batch_time

# Plot stddev
plt.figure()
plt.scatter(data['history'], data['stddev'], s=1, color='r')
plt.title("Stddev vs History")
plt.xlabel("History")
plt.ylabel("Stddev")

# Plot calc_int
plt.figure()
plt.scatter(data['history'], data['calc_int'], s=1, color='r')
plt.title("Calc Int vs History")
plt.xlabel("History")
plt.ylabel("Calc Int")

# # Plot batch_time
plt.figure()
plt.scatter(data['history'], data['batch_time'], s=1, color='r')
# plt.scatter(data2['history2'], data2['batch_time2'], s=1, color='b')
plt.title("Batch Time vs History")
plt.xlabel("History")
plt.ylabel("Batch Time")

plt.show()