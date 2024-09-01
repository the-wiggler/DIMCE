import cupy as cp
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
                x = cp.random.uniform(a, b, current_chunk)
                sum_function_output += cp.sum(f(x))
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
            for _ in range(batches):
                batch_start_time = time.time()
                num = [self.int_estimate(delta_histories, a, b, chunk_size, f) for _ in range(int_per_batch)]
                history_count_list.append(delta_histories)
                calc_int.extend(num)
                delta_histories += hist_factor
                # int_per_batch += 1
                batch_time = time.time() - batch_start_time
                batch_times.append(batch_time)
        print(calc_int)
        return cp.array(calc_int), cp.array(history_count_list), cp.array(batch_times)

    def avg_matrix(self, calc_int):
        print("INTEGRAL BATCH AVERAGE VALUES")
        avg_calc_int = cp.mean(calc_int.reshape(-1, self.int_per_batch), axis=1)
        print(avg_calc_int)
        return avg_calc_int

    def std_calc_int(self, batches, calc_int, repetition_factor):
        print("STANDARD DEVIATION FOR INTEGRAL OUTPUT ARRAY")
        std_list = cp.std(calc_int.reshape(-1, self.int_per_batch), axis=1)
        print(std_list)
        return std_list

    def sv_matrix(self, calc_int, avg_calc_int, std_list, history_count_list, batch_times):
        np.save(os.path.join('np_store_GPU', 'avg_calc_int_GPU.npy'), avg_calc_int.get())
        np.save(os.path.join('np_store_GPU', 'std_list_GPU.npy'), std_list.get())
        np.save(os.path.join('np_store_GPU', 'history_count_list_GPU.npy'), history_count_list.get())
        np.save(os.path.join('np_store_GPU', 'batch_times_GPU.npy'), batch_times.get())
