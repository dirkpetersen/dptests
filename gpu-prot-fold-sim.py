#! /usr/bin/env python3

import cupy as cp
#import numpy as np
import threading
#from queue import Queue
import math

def check_gpu_availability():
    """Check if CUDA GPU is available and return device count"""
    try:
        device_count = cp.cuda.runtime.getDeviceCount()
        if device_count == 0:
            raise RuntimeError("No CUDA devices found")
        
        # Try to get properties of first device
        with cp.cuda.Device(0):
            mem_info = cp.cuda.runtime.memGetInfo()
            free_memory = mem_info[0]
            total_memory = mem_info[1]
            print(f"Found {device_count} CUDA device(s)")
            print(f"GPU 0: {free_memory/1e9:.1f}GB free of {total_memory/1e9:.1f}GB total")
        
        return device_count
    except Exception as e:
        print(f"Error checking GPU: {str(e)}")
        return 0

# System Configuration
AVAILABLE_GPUS = check_gpu_availability()
if AVAILABLE_GPUS == 0:
    raise RuntimeError("No available CUDA devices")

# Get actual GPU memory
with cp.cuda.Device(0):
    mem_info = cp.cuda.runtime.memGetInfo()
    GPU_MEMORY_PER_DEVICE = mem_info[1]  # Total memory
    
# Very conservative memory usage - only use 50% of available memory
MEMORY_SAFETY_MARGIN = 0.9
BYTES_PER_ATOM = 4 * (3 * 2 + 1 + 1)  # positions, new_positions, sequence, energies

# Calculate maximum atoms while keeping memory usage reasonable
ATOMS_PER_GPU = int((GPU_MEMORY_PER_DEVICE * MEMORY_SAFETY_MARGIN) / BYTES_PER_ATOM)

print(f"Maximum atoms per GPU: {ATOMS_PER_GPU:,}")
MAX_TOTAL_ATOMS = ATOMS_PER_GPU * AVAILABLE_GPUS
print(f"Maximum total atoms: {MAX_TOTAL_ATOMS:,}")

class GPUMemoryTracker:
    @staticmethod
    def get_memory_info(gpu_id):
        try:
            with cp.cuda.Device(gpu_id):
                mem_info = cp.cuda.runtime.memGetInfo()
                return mem_info[0], mem_info[1]
        except Exception as e:
            print(f"Error getting memory info for GPU {gpu_id}: {str(e)}")
            return 0, 0
    
    @staticmethod
    def print_memory_usage(gpu_id):
        free, total = GPUMemoryTracker.get_memory_info(gpu_id)
        if total > 0:
            used = total - free
            print(f"GPU {gpu_id}: Using {used/1e9:.2f}GB of {total/1e9:.2f}GB "
                  f"({used/total*100:.1f}%)")

class MemoryEfficientProteinFolding:
    def __init__(self, sequence_length=None):
        self.n_gpus = AVAILABLE_GPUS
        
        # Use a very conservative sequence length if none provided
        if sequence_length is None:
            self.sequence_length = min(1000, MAX_TOTAL_ATOMS)
        else:
            self.sequence_length = min(sequence_length, MAX_TOTAL_ATOMS)
        
        self.atoms_per_gpu = max(1, self.sequence_length // self.n_gpus)
        print(f"Atoms per GPU: {self.atoms_per_gpu:,}")
        
        self.gpu_data = {}
        self._initialize_gpu_memory()
    
    def _initialize_gpu_memory(self):
        """Initialize GPU memory with conservative sizes"""
        try:
            with cp.cuda.Device(0):  # Only use first GPU for now
                # Small test allocation to verify memory
                test_array = cp.zeros((10, 3), dtype=cp.float32)
                del test_array
                
                # Initialize actual data with small sizes
                self.gpu_data[0] = {
                    'positions': cp.random.uniform(
                        -1, 1, 
                        (min(1000, self.atoms_per_gpu), 3),
                        dtype=cp.float32
                    ),
                    'new_positions': cp.zeros(
                        (min(1000, self.atoms_per_gpu), 3),
                        dtype=cp.float32
                    ),
                    'sequence': cp.array(
                        [i % 2 for i in range(min(1000, self.atoms_per_gpu))],
                        dtype=cp.int32
                    ),
                    'energies': cp.zeros(min(1000, self.atoms_per_gpu), dtype=cp.float32)
                }
                
                GPUMemoryTracker.print_memory_usage(0)
        except Exception as e:
            print(f"Error initializing GPU memory: {str(e)}")
            raise

    def _calculate_energy(self, positions, sequence):
        """Calculate energy using simple operations"""
        try:
            # Use smaller chunks for calculation
            chunk_size = 100
            total_energy = 0.0
            
            for i in range(0, len(positions), chunk_size):
                chunk_pos = positions[i:i+chunk_size]
                chunk_seq = sequence[i:i+chunk_size]
                
                # Simple distance calculation
                diff = chunk_pos[:, None, :] - chunk_pos[None, :, :]
                distances = cp.sqrt(cp.sum(diff * diff, axis=2) + 1e-10)
                
                # Simple energy calculation
                energy = cp.sum(1.0 / distances)
                total_energy += float(energy)
            
            return total_energy
        except Exception as e:
            print(f"Error in energy calculation: {str(e)}")
            return 0.0

    def _simulation_step(self, gpu_id, step, temperature):
        """Perform one simulation step"""
        try:
            with cp.cuda.Device(gpu_id):
                positions = self.gpu_data[gpu_id]['positions']
                sequence = self.gpu_data[gpu_id]['sequence']
                
                # Small random movement
                displacement = cp.random.uniform(-0.1, 0.1, positions.shape)
                new_positions = positions + displacement
                
                if step % 100 == 0:
                    # Calculate energies occasionally
                    old_energy = self._calculate_energy(positions, sequence)
                    new_energy = self._calculate_energy(new_positions, sequence)
                    
                    if new_energy < old_energy:
                        self.gpu_data[gpu_id]['positions'] = new_positions
        except Exception as e:
            print(f"Error in simulation step: {str(e)}")

    def run_simulation(self, n_steps=1000):
        """Run a short simulation"""
        try:
            for step in range(n_steps):
                if step % 100 == 0:
                    print(f"Step {step}/{n_steps}")
                self._simulation_step(0, step, 1.0)
                
                if step % 100 == 0:
                    GPUMemoryTracker.print_memory_usage(0)
        except Exception as e:
            print(f"Error in simulation: {str(e)}")

def main():
    try:
        print("Initializing simulation...")
        simulator = MemoryEfficientProteinFolding(sequence_length=100000)  # Use more atoms
        
        print("Starting simulation...")
        simulator.run_simulation(n_steps=100)  # Run for just 100 steps
        
    except Exception as e:
        print(f"Error during simulation: {str(e)}")
    finally:
        print("\nFinal memory usage:")
        GPUMemoryTracker.print_memory_usage(0)

if __name__ == "__main__":
    main()
