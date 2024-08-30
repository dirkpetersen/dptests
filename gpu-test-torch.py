#! /usr/bin/env python3

import os, sys
try :
    import torch
    import pynvml
except ImportError:
    print('please run first: python3 -m pip install --upgrade --no-cache torch pynvml')
    sys.exit(1)

def get_free_gpu_memory(device):
    return torch.cuda.get_device_properties(device).total_memory - torch.cuda.memory_allocated(device)

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
    print(f"PyTorch version: {torch.__version__}")
    print(f"CUDA available: {torch.cuda.is_available()}")

    print("\nSystem GPU Information:")
    for gpu in get_gpu_info():
        print(f"  GPU {gpu[0]}: {gpu[1]}, UUID: ....{gpu[2][-4:]}, Free Memory: {gpu[3] / (1024**2):.0f} MB")

    if torch.cuda.is_available():
        print(f"\nCUDA version: {torch.version.cuda}")
        print(f"Number of CUDA devices visible to PyTorch: {torch.cuda.device_count()}")

        for i in range(torch.cuda.device_count()):
            print(f"\nPyTorch CUDA Device {i}:")
            print(f"  Name: {torch.cuda.get_device_name(i)}")
            print(f"  Capability: {torch.cuda.get_device_capability(i)}")

            pynvml.nvmlInit()
            handle = pynvml.nvmlDeviceGetHandleByIndex(i)
            uuid = decode_if_bytes(pynvml.nvmlDeviceGetUUID(handle))
            pynvml.nvmlShutdown()
            print(f"  UUID: ....{uuid[-4:]}")

        current_device = torch.cuda.current_device()
        print(f"\nCurrent CUDA device index (from PyTorch's perspective): {current_device}")

        pynvml.nvmlInit()
        handle = pynvml.nvmlDeviceGetHandleByIndex(current_device)
        actual_gpu_id = pynvml.nvmlDeviceGetIndex(handle)
        actual_gpu_uuid = decode_if_bytes(pynvml.nvmlDeviceGetUUID(handle))
        pynvml.nvmlShutdown()

        print(f"Actual Device ID of GPU PyTorch is using: {actual_gpu_id}")
        print(f"UUID of GPU PyTorch is using: ....{actual_gpu_uuid[-4:]}")

        if cuda_visible_devices is not None:
            specified_gpu_ids = [int(x) for x in cuda_visible_devices.split(',')]
            print(f"Device ID(s) specified by CUDA_VISIBLE_DEVICES: {specified_gpu_ids}")

            if actual_gpu_id in specified_gpu_ids:
                print("VERIFIED: PyTorch is using a GPU specified by CUDA_VISIBLE_DEVICES")
                print(f"PyTorch's device {current_device} corresponds to actual GPU {actual_gpu_id}")
            else:
                print("WARNING: PyTorch is NOT using a GPU specified by CUDA_VISIBLE_DEVICES")
                print("This might indicate a problem with CUDA_VISIBLE_DEVICES setting or its recognition.")

        free_memory = get_free_gpu_memory(current_device)
        print(f"\nFree GPU memory before test: {free_memory / 1e9:.2f} GB")

        try:
            x = torch.rand(100, 100).to(f"cuda:{current_device}")
            y = torch.matmul(x, x.t())
            print("GPU test successful: Tensor operations completed on GPU.")

            free_memory_after = get_free_gpu_memory(current_device)
            print(f"Free GPU memory after test: {free_memory_after / 1e9:.2f} GB")

        except Exception as e:
            print(f"GPU test failed: {str(e)}")
    else:
        print("No CUDA GPUs are available for PyTorch.")

if __name__ == "__main__":
    print("Starting GPU test...")
    test_gpu()