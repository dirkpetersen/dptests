#!/usr/bin/env python3

import os
import sys
try:
    import pycuda.autoinit
    import pycuda.driver as drv
    from pycuda import gpuarray
    import pynvml
except ImportError:
    print('please run first: python3 -m pip install --upgrade --no-cache pycuda pynvml')
    sys.exit(1)

def get_free_gpu_memory(device):
    return device.total_memory() - pycuda.driver.mem_get_info()[1]

def decode_if_bytes(value):
    return value.decode('utf-8') if isinstance(value, bytes) else value

def get_gpu_info():
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

def test_gpu():
    cuda_visible_devices = os.environ.get('CUDA_VISIBLE_DEVICES')
    print(f"CUDA_VISIBLE_DEVICES: {cuda_visible_devices}")
    print(f"PyCUDA version: {'.'.join(map(str, drv.get_version()))}")
    print(f"CUDA available: {drv.get_version() != (0, 0, 0)}")

    print("\nSystem GPU Information:")
    for gpu in get_gpu_info():
        print(f"  GPU {gpu[0]}: {gpu[1]}, UUID: ....{gpu[2][-4:]}, Free Memory: {gpu[3] / (1024**2):.0f} MB")

    if drv.get_version() != (0, 0, 0):
        print(f"\nCUDA version: {'.'.join(map(str, drv.get_version()))}")
        print(f"Number of CUDA devices visible to PyCUDA: {drv.Device.count()}")

        for i in range(drv.Device.count()):
            device = drv.Device(i)
            print(f"\nPyCUDA CUDA Device {i}:")
            print(f"  Name: {device.name()}")
            print(f"  Capability: {device.compute_capability()}")

            pynvml.nvmlInit()
            handle = pynvml.nvmlDeviceGetHandleByIndex(i)
            uuid = decode_if_bytes(pynvml.nvmlDeviceGetUUID(handle))
            pynvml.nvmlShutdown()
            print(f"  UUID: ....{uuid[-4:]}")

        current_device = pycuda.autoinit.device
        current_device_index = current_device.get_attribute(drv.device_attribute.PCI_DEVICE_ID)
        print(f"\nCurrent CUDA device index: {current_device_index}")

        pynvml.nvmlInit()
        handle = pynvml.nvmlDeviceGetHandleByIndex(current_device_index)
        actual_gpu_uuid = decode_if_bytes(pynvml.nvmlDeviceGetUUID(handle))
        pynvml.nvmlShutdown()

        print(f"UUID of GPU PyCUDA is using: ....{actual_gpu_uuid[-4:]}")

        if cuda_visible_devices is not None:
            specified_gpu_ids = [int(x) for x in cuda_visible_devices.split(',')]
            print(f"Device ID(s) specified by CUDA_VISIBLE_DEVICES: {specified_gpu_ids}")

            if current_device_index in specified_gpu_ids:
                print("VERIFIED: PyCUDA is using a GPU specified by CUDA_VISIBLE_DEVICES")
                print(f"PyCUDA's device {current_device_index} corresponds to actual GPU {current_device_index}")
            else:
                print("WARNING: PyCUDA is NOT using a GPU specified by CUDA_VISIBLE_DEVICES")
                print("This might indicate a problem with CUDA_VISIBLE_DEVICES setting or its recognition.")

        free_memory = get_free_gpu_memory(current_device)
        print(f"\nFree GPU memory before test: {free_memory / 1e9:.2f} GB")

        try:
            x = gpuarray.to_gpu(drv.pagelocked_empty((100, 100), dtype=float))
            x.fill(0.5)
            y = gpuarray.empty_like(x)
            drv.Context.synchronize()
            print("GPU test successful: GPUArray operations completed on GPU.")

            free_memory_after = get_free_gpu_memory(current_device)
            print(f"Free GPU memory after test: {free_memory_after / 1e9:.2f} GB")

        except Exception as e:
            print(f"GPU test failed: {str(e)}")
    else:
        print("No CUDA GPUs are available for PyCUDA.")

if __name__ == "__main__":
    print("Starting GPU test...")
    test_gpu()
