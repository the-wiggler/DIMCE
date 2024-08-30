import numpy as np
from tqdm import tqdm
import time
import os

class IntegralMC:
    def __init__(self, batches, int_per_batch, chunk_size, a, b, histories, f, hist_factor, repetition_factor):
        self.batches = batches
        self.int_per_batch = int_per_batch
        self.chunk_size = chunk_size
        self.a = a
        self.b = b
        self.histories = histories
        self.f = f
        self.hist_factor = hist_factor
        self.repetition_factor = repetition_factor

    def calc(self):
        start = time.time()
        calc_int, history_count_list, batch_times = self.int_loop(self.batches, self.histories, self.int_per_batch,
                                                                  self.a, self.b, self.chunk_size, self.f,
                                                                  self.hist_factor, self.repetition_factor)
        self.int_print(self.batches, history_count_list, calc_int)
        avg_calc_int = self.avg_matrix(calc_int)
        std_list = self.std_calc_int(self.batches, calc_int, self.repetition_factor)
        self.sv_matrix(calc_int, avg_calc_int, std_list, history_count_list, batch_times)
        print(f'EXECUTION TIME: {time.time() - start} SECONDS')

    def int_estimate(self, n, a, b, chunk_size, f):
        sum_function_output = 0
        total_points = 0
        num_chunks = (n + chunk_size - 1) // chunk_size
        with tqdm(total=num_chunks) as pbar:
            while n > 0:
                current_chunk = min(n, chunk_size)
                x = np.random.uniform(a, b, current_chunk)
                sum_function_output += np.sum(f(x))
                total_points += current_chunk
                n -= current_chunk
                pbar.update(1)
            estimate_output = (b - a) * (sum_function_output / total_points)
            return estimate_output

    def int_loop(self, batches, histories, int_per_batch, a, b, chunk_size, f, hist_factor, repetition_factor):
        calc_int = []
        history_count_list = []
        batch_times = []
        for _ in range(repetition_factor):
            delta_histories = histories
            for i in range(batches):
                batch_start_time = time.time()
                num = []
                for i in range(int_per_batch):
                    num.append(self.int_estimate(delta_histories, a, b, chunk_size, f))
                history_count_list.append([delta_histories])
                calc_int.append(num)
                delta_histories += hist_factor
                batch_time = time.time() - batch_start_time
                batch_times_y = [batch_time]
                batch_times.append(batch_times_y)
        batch_times = np.array(batch_times)
        calc_int = np.array(calc_int)
        history_count_list = np.array(history_count_list)
        print(history_count_list)
        return calc_int, history_count_list, batch_times

    def int_print(self, batches, int_per_batch, calc_int):
        print(f"INTEGRAL OUTPUT ARRAY: [{batches} BATCHES WITH {int_per_batch} ESTIMATES]")
        calc_int = np.array(calc_int)
        print(calc_int)

    def avg_matrix(self, calc_int):
        print("INTEGRAL BATCH AVERAGE VALUES")
        avg_calc_int = []
        for i in range(len(calc_int)):
            avg_calc_int_y = [np.average(calc_int[i])]
            avg_calc_int.append(avg_calc_int_y)
        avg_calc_int = np.array(avg_calc_int)
        print(avg_calc_int)
        return avg_calc_int

    def std_calc_int(self, batches, calc_int, repetition_factor):
        std_list = []
        for l in range(repetition_factor):
            for i in range(batches):
                std_y = [np.std(calc_int[i])]
                std_list.append(std_y)
        std_list = np.array(std_list)
        print("STANDARD DEVIATION FOR INTEGRAL OUTPUT ARRAY")
        print(std_list)
        return std_list

    def sv_matrix(self, calc_int, avg_calc_int, std_list, history_count_list, batch_times):
        np.save(os.path.join('np_store', 'calc_int.npy'), calc_int)
        np.save(os.path.join('np_store', 'avg_calc_int.npy'), avg_calc_int)
        np.save(os.path.join('np_store', 'std_list.npy'), std_list)
        np.save(os.path.join('np_store', 'history_count_list.npy'), history_count_list)
        np.save(os.path.join('np_store', 'batch_times.npy'), batch_times)
