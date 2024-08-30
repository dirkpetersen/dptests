#!/usr/bin/env python3

import os
import sys
import ctypes
import ctypes.util
import time
import numpy as np

def check_cuda_availability():
    try:
        cuda_path = ctypes.util.find_library('cudart')
        if cuda_path is None:
            return False, "CUDA runtime library (libcudart.so) not found. CUDA may not be installed or not in the library path."

        ctypes.CDLL(cuda_path)
        return True, "CUDA runtime library found."
    except Exception as e:
        return False, f"Error loading CUDA runtime library: {str(e)}"

cuda_available, cuda_message = check_cuda_availability()

if not cuda_available:
    print(cuda_message)
    print('CUDA is not available. GPU functions will be disabled.')
    sys.exit(1)

try:
    import pycuda.autoinit
    import pycuda.driver as drv
    from pycuda import gpuarray
    from pycuda.elementwise import ElementwiseKernel
except ImportError:
    print('PyCUDA not found. Please run: python3 -m pip install --upgrade --no-cache pycuda')
    sys.exit(1)

def test_gpu():
    cuda_visible_devices = os.environ.get('CUDA_VISIBLE_DEVICES')
    print(f"CUDA_VISIBLE_DEVICES: {cuda_visible_devices}")

    print(f"PyCUDA version: {'.'.join(map(str, drv.get_version()))}")
    print(f"CUDA available: {drv.get_version() != (0, 0, 0)}")

    print(f"\nNumber of CUDA devices: {drv.Device.count()}")

    for i in range(drv.Device.count()):
        gpu = drv.Device(i)
        print(f"\nCUDA Device {i}:")
        print(f"  Name: {gpu.name()}")
        print(f"  Compute Capability: {gpu.compute_capability()}")
        print(f"  Total Memory: {gpu.total_memory() // (1024**2)} MB")

    current_device = pycuda.autoinit.device
    print(f"\nCurrent CUDA device:")
    print(f"  Name: {current_device.name()}")
    print(f"  Compute Capability: {current_device.compute_capability()}")
    print(f"  Total Memory: {current_device.total_memory() // (1024**2)} MB")

    if cuda_visible_devices is not None:
        specified_gpu_ids = [int(x) for x in cuda_visible_devices.split(',')]
        print(f"Device ID(s) specified by CUDA_VISIBLE_DEVICES: {specified_gpu_ids}")

        current_device_in_visible = any(
            current_device.name() == drv.Device(i).name() and
            current_device.compute_capability() == drv.Device(i).compute_capability()
            for i in specified_gpu_ids
        )

        if current_device_in_visible:
            print("VERIFIED: PyCUDA is using a GPU specified by CUDA_VISIBLE_DEVICES")
        else:
            print("WARNING: PyCUDA might not be using a GPU specified by CUDA_VISIBLE_DEVICES")
            print("This might indicate a problem with CUDA_VISIBLE_DEVICES setting or its recognition.")

    print("\nStarting GPU test...")
    try:
        # Prepare data
        n = 10_000_000
        a = np.ones(n, dtype=np.float32)
        b = np.full(n, 2, dtype=np.float32)

        # Transfer data to GPU
        a_gpu = gpuarray.to_gpu(a)
        b_gpu = gpuarray.to_gpu(b)

        # Define the elementwise addition operation
        vector_add = ElementwiseKernel(
            "float *a, float *b, float *c",
            "c[i] = a[i] + b[i]",
            "vector_add"
        )

        # Perform vector addition
        start_time = time.time()
        result_gpu = gpuarray.empty_like(a_gpu)
        vector_add(a_gpu, b_gpu, result_gpu)
        gpu_time = time.time() - start_time

        # Transfer result back to CPU
        result = result_gpu.get()

        print(f"GPU computation time: {gpu_time:.4f} seconds")
        print(f"First 10 elements of input 'a': {a[:10]}")
        print(f"First 10 elements of input 'b': {b[:10]}")
        print(f"First 10 elements of result: {result[:10]}")
        print(f"Last 10 elements of result: {result[-10:]}")

        # Verify result
        expected_result = a + b
        if np.allclose(result, expected_result):
            print("GPU test successful: Vector addition completed correctly on GPU.")
        else:
            print("GPU test failed: Results do not match expected values.")
            print(f"Max difference: {np.max(np.abs(result - expected_result))}")
            print(f"Mean difference: {np.mean(np.abs(result - expected_result))}")

    except Exception as e:
        print(f"GPU test failed: {str(e)}")

if __name__ == "__main__":
    print("Starting GPU test...")
    test_gpu()

