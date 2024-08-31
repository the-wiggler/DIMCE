import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm


class PlotArrays:
    def __init__(self, avg_calc_int, std_list, history_count_list, batch_times):
        self.data = pd.DataFrame({
            'history_count': np.array(history_count_list),
            'avg_calc_int': np.array(avg_calc_int),
            'std': np.array(std_list),
            'batch_times': np.array(batch_times)
        })

    def rungraphs(self):
        self.intvhistory()
        self.stdvhistory()
        self.timevhistory()
        plt.show()

    def intvhistory(self):
        plt.figure(1)
        plt.suptitle("Integral Calculation Value v. History Count")
        plt.title(r'$\int_{0}^{5} f(x) = \sin{(x)}\, dx$', fontsize=10)
        plt.xlabel("History Count")
        plt.ylabel("Integral Calculation Value")
        plt.xscale('linear')

        x = self.data['history_count']
        y = self.data['avg_calc_int']

        plt.scatter(x, y, s=0.05, alpha=1, color='r')

        X = sm.add_constant(x)
        model = sm.OLS(y, X).fit()
        x_fit = np.linspace(x.min(), x.max(), 500)
        X_fit = sm.add_constant(x_fit)
        y_fit = model.predict(X_fit)

        plt.plot(x_fit, y_fit, 'r-', label=f"Fit: {model.params[1]:.3e}*x + {model.params[0]:.3e}")
        plt.ylim(0, y.max())
        plt.legend()

    def stdvhistory(self):
        plt.figure(2)
        plt.suptitle("Standard Deviation v. History Count")
        plt.xlabel("History Count")
        plt.xscale('linear')
        plt.ylabel("Standard Deviation")

        x = self.data['history_count']
        y = self.data['std']

        plt.scatter(x, y, s=0.9, alpha=1, color='b', label='Data')

        def inverse_x(x, a, b):
            return a / x ** 2 + b

        model = sm.OLS(y, sm.add_constant(1 / x)).fit()
        x_fit = np.linspace(x.min(), x.max(), 500)
        y_fit = inverse_x(x_fit, model.params[1], model.params[0])

        plt.plot(x_fit, y_fit, color='red', label=fr'Fitted curve: $y = \frac{{{model.params[1]:.3f}}}{{x^2}}$ + {{model.params[0]:.3f}}')
        plt.ylim(0, y.max())
        plt.legend()

    def timevhistory(self):
        plt.figure(3)
        plt.suptitle("Batch Times v. History Count")
        plt.title(r'$\int_{0}^{5} f(x) = {f}, dx$', fontsize=10)
        plt.xlabel("History Count")
        plt.xscale('linear')
        plt.ylabel("Batch Times (s)")

        x = self.data['history_count']
        y = self.data['batch_times']

        plt.scatter(x, y, s=0.9, alpha=1, color='b', label="Data Points")

        X = sm.add_constant(x)
        model = sm.OLS(y, X).fit()
        x_fit = np.linspace(x.min(), x.max(), 500)
        X_fit = sm.add_constant(x_fit)
        y_fit = model.predict(X_fit)

        plt.plot(x_fit, y_fit, 'r-', label=f"Fit: {model.params[1]:.3e}*x + {model.params[0]:.3e}")
        plt.ylim(0, y.max())
        plt.legend()


# Usage
plot_arrays = PlotArrays(
    avg_calc_int=np.load(os.path.join('np_store_GPU', 'avg_calc_int_GPU.npy')),
    std_list=np.load(os.path.join('np_store_GPU', 'std_list_GPU.npy')),
    history_count_list=np.load(os.path.join('np_store_GPU', 'history_count_list_GPU.npy')),
    batch_times=np.load(os.path.join('np_store_GPU', 'batch_times_GPU.npy')),
)