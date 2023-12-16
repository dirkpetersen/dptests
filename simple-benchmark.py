#! /usr/bin/env python3

# a simple benchmark that spends fifty/fifty in CPU and file io tasks
# So far in AWS: 
# m7a (Epyc): 520, m7i (Xeon): 380, m7g (Graviton3): 340, m4 (Xeon): 153

import time, random, hashlib, os

class BigBadClass:
    def run_benchmark(self, duration=10):
        start_time = time.time()
        elapsed = 0
        total_iterations = 0
        temp_file_path = '/dev/shm/simple_benchmark_test.txt'

        while elapsed < duration:
            # Perform a CPU-intensive task
            self.perform_intensive_calculation()

            # Perform a file I/O operation in memory
            self.perform_file_io_operation(temp_file_path)

            total_iterations += 1
            elapsed = time.time() - start_time

        # Clean up
        if os.path.exists(temp_file_path):
            os.remove(temp_file_path)

        print(f"Benchmark completed in {elapsed:.2f} seconds.")
        print(f"Total iterations: {total_iterations}")

    def perform_intensive_calculation(self):
        # Simulate a CPU-bound task such as calculating hash values
        for _ in range(10000):
            hashlib.sha256(str(random.random()).encode()).hexdigest()

    def perform_file_io_operation(self, file_path):
        # Write and read from a memory-based file
        with open(file_path, 'w') as file:
            for _ in range(15000):
                file.write(str(random.random()) + '\n')

        with open(file_path, 'r') as file:
            for _ in file:
                pass

# Example usage
big_bad = BigBadClass()
big_bad.run_benchmark()