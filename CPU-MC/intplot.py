import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


class PlotArrays:
    def __init__(self, calc_int, avg_calc_int, std_list, history_count_list, batch_times):
        self.calc_int = calc_int
        self.avg_calc_int = avg_calc_int
        self.std_list = std_list
        self.history_count_list = history_count_list
        self.batch_times = batch_times

    def rungraphs(self):
        self.intvhistory(self.history_count_list, self.avg_calc_int)
        self.stdvhistory(self.history_count_list, self.std_list)
        self.timevhistory(self.history_count_list, self.batch_times)
        plt.show()

    def intvhistory(self, history_count_list, avg_calc_int):
        x = np.array(history_count_list).flatten()
        y = np.array(avg_calc_int).flatten()
        x = np.array(history_count_list).flatten()
        y = np.array(avg_calc_int).flatten()
        # curve fit function
        def func(x, a, b):
            return a * x + b
        popt, pcov = curve_fit(func, x, y)
        x_fit = np.linspace(min(x), max(x), 500)
        y_fit = func(x_fit, *popt)
        a, b = popt

        plt.figure(1)
        plt.suptitle("Integral Calculation Value v. History Count")
        plt.title(r'$\int_{0}^{5} f(x) = \sin{(x)}\, dx$', fontsize=10)
        plt.xlabel("History Count")
        plt.ylabel("Integral Calculation Value")
        plt.xscale('linear')
        plt.scatter(x, y, s=0.05, alpha=1, color='r')
        plt.plot(x_fit, y_fit, 'r-', label=f"Fit: {a:.3e}*x + {b:.3e}")
        plt.legend()

    def stdvhistory(self, history_count_list, std_list):
        x = np.array(history_count_list).flatten()
        y = np.array(std_list).flatten()
        # Curve Fit
        def func(x, a, b):
            return a / x + b
        popt, pcov = curve_fit(func, x, y)
        x_fit = np.linspace(min(x), max(x), 500)
        y_fit = func(x_fit, *popt)
        a, b = popt
        # plot
        plt.figure(2)
        plt.suptitle("Standard Deviation v. History Count")
        plt.xlabel("History Count")
        plt.xscale('linear')
        plt.ylabel("Standard Deviation")
        plt.scatter(x, y, s=0.9, alpha=1, color='b', label='Data')
        plt.plot(x_fit, y_fit, color='red', label=fr'Fitted curve: $y = \frac{{{a:.3f}}}{{x}}$ + {{{b:.3f}}}')
        plt.legend()

    def timevhistory(self, history_count_list, batch_times):
        x = np.array(history_count_list).flatten()
        y = np.array(batch_times).flatten()
        # curve fit function
        def func(x, a, b):
            return a * x + b
        popt, pcov = curve_fit(func, x, y)
        x_fit = np.linspace(min(x), max(x), 500)
        y_fit = func(x_fit, *popt)
        a, b = popt
        # plot
        plt.figure(3)
        plt.suptitle("Batch Times v. History Count")
        plt.title(r'$\int_{0}^{5} f(x) = {f}, dx$', fontsize=10)
        plt.xlabel("History Count")
        plt.xscale('linear')
        plt.ylabel("Batch Times (s)")
        plt.scatter(x, y, s=0.9, alpha=1, color='b', label="Data Points")
        plt.plot(x_fit, y_fit, 'r-', label=f"Fit: {a:.3e}*x + {b:.3e}")
        plt.legend()

    # def inprogress(self):
        # batch size (int_per_batch) v. accuracy (# of batches constant)
        # x = int_per_batch
        # y = avg_calc_int
        # plt.figure(1)
        # plt.suptitle("Integral Calculation Value v. History Count")
        # plt.title(r'$\int_{0}^{5} f(x) = \sin{(x)}\, dx$', fontsize=10)
        # plt.xlabel("History Count")
        # plt.ylabel("Integral Calculation Value")
        # plt.xscale('linear')
        # plt.scatter(x , y, s=0.05, alpha=1, color='r')

        # chunk_size v. execution time (all else constant)
        # different functions v. execution time


plot_arrays = PlotArrays(
    calc_int=np.load(os.path.join('np_store', 'calc_int.npy')),
    avg_calc_int=np.load(os.path.join('np_store', 'avg_calc_int.npy')),
    std_list=np.load(os.path.join('np_store', 'std_list.npy')),
    history_count_list=np.load(os.path.join('np_store', 'history_count_list.npy')),
    batch_times=np.load(os.path.join('np_store', 'batch_times.npy')),
)

