#!/usr/bin/env python3

import os
import sys
import ctypes
import ctypes.util
import random
import array
import math
import time

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
    from numba import cuda
    import pynvml
    numba_cuda_available = True
except ImportError:
    print('Numba CUDA or pynvml not found. Please run: python3 -m pip install --upgrade --no-cache numba pynvml')
    sys.exit(1)

def decode_if_bytes(value):
    return value.decode('utf-8') if isinstance(value, bytes) else value

def get_gpu_info():
    try:
        pynvml.nvmlInit()
        num_gpus = pynvml.nvmlDeviceGetCount()
        gpu_info = []
        for i in range(num_gpus):
            handle = pynvml.nvmlDeviceGetHandleByIndex(i)
            name = decode_if_bytes(pynvml.nvmlDeviceGetName(handle))
            uuid = decode_if_bytes(pynvml.nvmlDeviceGetUUID(handle))
            memory = pynvml.nvmlDeviceGetMemoryInfo(handle)
            gpu_info.append((i, name, uuid, memory.free))
        pynvml.nvmlShutdown()
        return gpu_info
    except Exception as e:
        print(f"Error getting GPU info: {str(e)}")
        return []

@cuda.jit
def vector_add(a, b, result):
    i = cuda.grid(1)
    if i < result.shape[0]:
        result[i] = a[i] + b[i]

def test_gpu():
    cuda_visible_devices = os.environ.get('CUDA_VISIBLE_DEVICES')
    print(f"CUDA_VISIBLE_DEVICES: {cuda_visible_devices}")

    print(f"Numba CUDA version: {cuda.runtime.get_version()}")
    print(f"CUDA available: {cuda.is_available()}")

    print("\nSystem GPU Information:")
    gpu_info = get_gpu_info()
    if gpu_info:
        for gpu in gpu_info:
            print(f"  GPU {gpu[0]}: {gpu[1]}, UUID: ....{gpu[2][-4:]}, Free Memory: {gpu[3] / (1024**2):.0f} MB")
    else:
        print("  No GPU information available.")

    print(f"\nNumber of CUDA devices visible to Numba: {len(cuda.gpus)}")

    for i, gpu in enumerate(cuda.gpus):
        print(f"\nNumba CUDA Device {i}:")
        print(f"  Name: {gpu.name}")
        print(f"  Compute Capability: {gpu.compute_capability}")

        pynvml.nvmlInit()
        handle = pynvml.nvmlDeviceGetHandleByIndex(i)
        uuid = decode_if_bytes(pynvml.nvmlDeviceGetUUID(handle))
        pynvml.nvmlShutdown()
        print(f"  UUID: ....{uuid[-4:]}")

    current_device = cuda.get_current_device()
    print(f"\nCurrent CUDA device index (from Numba's perspective): {current_device.id}")

    pynvml.nvmlInit()
    handle = pynvml.nvmlDeviceGetHandleByIndex(current_device.id)
    actual_gpu_id = pynvml.nvmlDeviceGetIndex(handle)
    actual_gpu_uuid = decode_if_bytes(pynvml.nvmlDeviceGetUUID(handle))
    pynvml.nvmlShutdown()

    print(f"Actual Device ID of GPU Numba is using: {actual_gpu_id}")
    print(f"UUID of GPU Numba is using: ....{actual_gpu_uuid[-4:]}")

    if cuda_visible_devices is not None:
        specified_gpu_ids = [int(x) for x in cuda_visible_devices.split(',')]
        print(f"Device ID(s) specified by CUDA_VISIBLE_DEVICES: {specified_gpu_ids}")

        if actual_gpu_id in specified_gpu_ids:
            print("VERIFIED: Numba is using a GPU specified by CUDA_VISIBLE_DEVICES")
            print(f"Numba's device {current_device.id} corresponds to actual GPU {actual_gpu_id}")
        else:
            print("WARNING: Numba is NOT using a GPU specified by CUDA_VISIBLE_DEVICES")
            print("This might indicate a problem with CUDA_VISIBLE_DEVICES setting or its recognition.")

    print("\nStarting GPU test...")
    try:
        # Prepare data
        n = 1_000_000
        a = array.array('f', (random.random() for _ in range(n)))
        b = array.array('f', (random.random() for _ in range(n)))
        result = array.array('f', (0.0 for _ in range(n)))

        # Copy data to device
        d_a = cuda.to_device(a)
        d_b = cuda.to_device(b)
        d_result = cuda.to_device(result)

        # Set up grid and block sizes
        threads_per_block = 256
        blocks_per_grid = (n + threads_per_block - 1) // threads_per_block

        # Launch kernel
        vector_add[blocks_per_grid, threads_per_block](d_a, d_b, d_result)

        # Copy result back to host
        d_result.copy_to_host(result)

        # Verify result
        for i in range(n):
            if not math.isclose(result[i], a[i] + b[i], rel_tol=1e-5):
                raise ValueError(f"Mismatch at index {i}: {result[i]} != {a[i] + b[i]}")
        print("GPU test successful: Vector addition completed correctly on GPU.")

    except Exception as e:
        print(f"GPU test failed: {str(e)}")

if __name__ == "__main__":
    print("Starting GPU test...")
    test_gpu()
