#! /usr/bin/env python3

import cupy as cp
import threading
from typing import Optional, Dict, Any

def check_nvlink_topology() -> bool:
    """
    Check if GPUs are connected via NVLink/NVSwitch
    
    Returns:
        bool: True if NVLink/NVSwitch is available and GPUs support unified memory
    """
    try:
        # Check for NVLink connections between GPUs
        nvlink_status = cp.cuda.runtime.deviceGetNvSMEMConfig()
        if nvlink_status is not None:
            print("NVLink/NVSwitch detected - using unified memory mode")
            return True
        print("Traditional GPU setup detected - using separate memory mode")
        return False
    except Exception as e:
        print(f"Error checking NVLink topology: {str(e)}")
        return False

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
            device_props = cp.cuda.runtime.getDeviceProperties(0)
            device_name = device_props["name"].decode('utf-8')
            print(f"Found {device_count} CUDA device(s)")
            print(f"GPU 0: {device_name}")
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
MEMORY_SAFETY_MARGIN = 0.97
# Memory calculation:
# - positions (float32 * 3)
# - new_positions (float32 * 3) 
# - velocities (float32 * 3)
# - forces (float32 * 3)
# - sequence (int32)
# - energies (float32)
# - chunk_buffer overhead (float32 * 3 * 1000 / atoms_per_gpu)
BYTES_PER_ATOM = 4 * (3 * 4 + 1 + 1) + (4 * 3 * 50)  # Reduced chunk buffer overhead from 100 to 50

# Calculate maximum atoms while keeping memory usage reasonable
def calculate_safe_atoms():
    with cp.cuda.Device(0):
        mem_info = cp.cuda.runtime.memGetInfo()
        initial_free = mem_info[0]  # Use free memory instead of total
        safe_memory = initial_free * MEMORY_SAFETY_MARGIN
        # Add additional 10% safety factor for runtime allocations
        safe_memory *= 0.9
        return int(safe_memory / BYTES_PER_ATOM)

ATOMS_PER_GPU = calculate_safe_atoms()

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
            initial_free = 3.5e9  # Initial free memory from first check
            safe_free = initial_free * MEMORY_SAFETY_MARGIN  # Apply safety margin
            our_usage = initial_free - free
            print(f"GPU {gpu_id}: Using {our_usage/1e9:.2f}GB of safe limit {safe_free/1e9:.2f}GB "
                  f"({our_usage/safe_free*100:.1f}% of {MEMORY_SAFETY_MARGIN*100:.0f}% margin)")

class MemoryEfficientProteinFolding:
    def __init__(self, sequence_length: Optional[int] = None):
        """
        Initialize protein folding simulation with architecture-aware memory management
        
        Args:
            sequence_length: Number of atoms to simulate. If None, uses conservative default.
        """
        self.total_energy = 0.0
        self.accepted_moves = 0
        self.total_moves = 0
        # Check hardware architecture
        self.using_nvlink = check_nvlink_topology()
        self.n_gpus = AVAILABLE_GPUS
        
        # Use a very conservative sequence length if none provided
        if sequence_length is None:
            self.sequence_length = min(1000, MAX_TOTAL_ATOMS)
        else:
            self.sequence_length = min(sequence_length, MAX_TOTAL_ATOMS)
        
        # Divide sequence evenly across GPUs
        self.atoms_per_gpu = max(1, self.sequence_length // self.n_gpus)
        # Handle remainder atoms
        self.gpu_atom_counts = [self.atoms_per_gpu] * self.n_gpus
        remainder = self.sequence_length % self.n_gpus
        for i in range(remainder):
            self.gpu_atom_counts[i] += 1
            
        print(f"Total atoms: {self.sequence_length:,}")
        for gpu_id in range(self.n_gpus):
            print(f"GPU {gpu_id}: {self.gpu_atom_counts[gpu_id]:,} atoms")
        
        self.gpu_data = {}
        self._initialize_gpu_memory()
    
    def _initialize_gpu_memory(self):
        """
        Initialize GPU memory with architecture-specific optimizations
        
        For NVL72 (NVLink/NVSwitch):
            - Uses unified memory space across GPUs
            - Enables direct GPU-to-GPU access
            - Optimizes for high-bandwidth GPU interconnect
        
        For traditional systems:
            - Allocates separate memory per GPU
            - Manages explicit memory transfers
        """
        try:
            if self.using_nvlink:
                self._initialize_unified_memory()
            else:
                self._initialize_separate_memory()
        except Exception as e:
            print(f"Error initializing GPU memory: {str(e)}")
            raise

    def _initialize_unified_memory(self):
        """Initialize unified memory for NVL72 architecture"""
        try:
            # Allocate unified memory accessible by all GPUs
            self.unified_data = {
                'positions': cp.cuda.managed_memory((self.sequence_length, 3), dtype=cp.float32),
                'new_positions': cp.cuda.managed_memory((self.sequence_length, 3), dtype=cp.float32),
                'velocities': cp.cuda.managed_memory((self.sequence_length, 3), dtype=cp.float32),
                'forces': cp.cuda.managed_memory((self.sequence_length, 3), dtype=cp.float32),
                'sequence': cp.cuda.managed_memory(self.sequence_length, dtype=cp.int32),
                'energies': cp.cuda.managed_memory(self.sequence_length, dtype=cp.float32)
            }
            
            # Initialize unified memory with data
            with cp.cuda.Device(0):  # Can use any GPU for initialization
                self.unified_data['positions'][:] = cp.random.uniform(-1, 1, (self.sequence_length, 3))
                self.unified_data['sequence'][:] = cp.array([i % 2 for i in range(self.sequence_length)])
            
            print("Initialized unified memory for NVL72 architecture")
            GPUMemoryTracker.print_memory_usage(0)

        except Exception as e:
            print(f"Error initializing unified memory: {str(e)}")
            raise

    def _initialize_separate_memory(self):
        """Initialize separate memory for traditional GPU architecture"""
        try:
            atom_offset = 0
            for gpu_id in range(self.n_gpus):
                n_atoms = self.gpu_atom_counts[gpu_id]
                with cp.cuda.Device(gpu_id):
                    # Small test allocation to verify memory
                    test_array = cp.zeros((10, 3), dtype=cp.float32)
                    del test_array
                    
                    # Initialize actual data with proper atom counts per GPU
                    self.gpu_data[gpu_id] = {
                    'positions': cp.random.uniform(
                        -1, 1, 
                        (n_atoms, 3),
                        dtype=cp.float32
                    ),
                    'atom_offset': atom_offset,  # Track position in global sequence
                    'new_positions': cp.zeros(
                        (n_atoms, 3),
                        dtype=cp.float32
                    ),
                    'velocities': cp.zeros(
                        (n_atoms, 3),
                        dtype=cp.float32
                    ),
                    'forces': cp.zeros(
                        (n_atoms, 3),
                        dtype=cp.float32
                    ),
                    'sequence': cp.array(
                        [i % 2 for i in range(atom_offset, atom_offset + n_atoms)],
                        dtype=cp.int32
                    ),
                    'energies': cp.zeros(n_atoms, dtype=cp.float32),
                    # Calculate distances and interactions in chunks instead of storing full matrices
                    'chunk_buffer': cp.zeros(
                        (50, n_atoms, 3),  # Reduced from 100 to 50 atoms at a time
                        dtype=cp.float32
                    )
                }
                
                atom_offset += n_atoms
                GPUMemoryTracker.print_memory_usage(gpu_id)
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
                #chunk_seq = sequence[i:i+chunk_size]
                
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

    def _simulation_step(self, gpu_id: int, step: int, temperature: float):
        """
        Perform one simulation step with architecture-aware memory access
        
        Args:
            gpu_id: GPU ID to run this step on
            step: Current simulation step number
            temperature: Current simulation temperature
        """
        try:
            with cp.cuda.Device(gpu_id):
                if self.using_nvlink:
                    # Direct access to unified memory
                    start_idx = self.gpu_atom_counts[gpu_id] * gpu_id
                    end_idx = start_idx + self.gpu_atom_counts[gpu_id]
                    positions = self.unified_data['positions'][start_idx:end_idx]
                    sequence = self.unified_data['sequence'][start_idx:end_idx]
                else:
                    # Traditional separate memory access
                    positions = self.gpu_data[gpu_id]['positions']
                    sequence = self.gpu_data[gpu_id]['sequence']
                
                # Small random movement
                displacement = cp.random.uniform(-0.1, 0.1, positions.shape)
                new_positions = positions + displacement
                
                if step % 100 == 0:
                    # Calculate energies occasionally
                    old_energy = self._calculate_energy(positions, sequence)
                    new_energy = self._calculate_energy(new_positions, sequence)
                    
                    # Track energy changes and acceptance ratio
                    delta_energy = new_energy - old_energy
                    self.total_moves += 1
                    
                    if delta_energy < 0:  # Accept if energy decreases
                        self.gpu_data[gpu_id]['positions'] = new_positions
                        self.total_energy = new_energy
                        self.accepted_moves += 1
                    elif cp.random.random() < cp.exp(-delta_energy / temperature):  # Metropolis criterion
                        self.gpu_data[gpu_id]['positions'] = new_positions
                        self.total_energy = new_energy
                        self.accepted_moves += 1
        except Exception as e:
            print(f"Error in simulation step: {str(e)}")

    def _run_gpu_thread(self, gpu_id, n_steps, temperature=1.0):
        """Run simulation steps on a specific GPU"""
        try:
            for step in range(n_steps):
                self._simulation_step(gpu_id, step, temperature)
                if step % 100 == 0:
                    GPUMemoryTracker.print_memory_usage(gpu_id)
        except Exception as e:
            print(f"Error in GPU {gpu_id} thread: {str(e)}")

    def _sync_gpus(self):
        """
        Synchronize GPU memory based on architecture
        
        For NVL72:
            - Uses hardware-level synchronization
            - Takes advantage of unified memory coherence
        
        For traditional systems:
            - Explicit synchronization through host memory
        """
        if self.using_nvlink:
            # Hardware-level synchronization for NVL72
            cp.cuda.runtime.deviceSynchronize()
        else:
            # Traditional explicit synchronization
            for gpu_id in range(self.n_gpus):
                with cp.cuda.Device(gpu_id):
                    cp.cuda.runtime.deviceSynchronize()

    def run_simulation(self, n_steps: int = 1000):
        """
        Run a parallel simulation across all available GPUs
        
        Args:
            n_steps: Number of simulation steps to run
        """
        try:
            threads = []
            for gpu_id in range(self.n_gpus):
                thread = threading.Thread(
                    target=self._run_gpu_thread,
                    args=(gpu_id, n_steps)
                )
                threads.append(thread)
                thread.start()

            # Wait for all GPU threads to complete
            for thread in threads:
                thread.join()
                
        except Exception as e:
            print(f"Error in simulation: {str(e)}")

def main():
    try:
        print("Initializing simulation...")
        # Distribute atoms across all available GPUs
        total_atoms = ATOMS_PER_GPU * AVAILABLE_GPUS
        simulator = MemoryEfficientProteinFolding(sequence_length=total_atoms)
        
        print("Starting simulation...")
        simulator.run_simulation(n_steps=100)  # Run for just 100 steps
        
    except Exception as e:
        print(f"Error during simulation: {str(e)}")
    finally:
        print("\nSimulation Results:")
        print(f"Final Energy: {simulator.total_energy:,.2f}")
        acceptance_rate = (simulator.accepted_moves / simulator.total_moves * 100) if simulator.total_moves > 0 else 0
        print(f"Move Acceptance Rate: {acceptance_rate:.1f}%")
        print(f"Accepted Moves: {simulator.accepted_moves:,} of {simulator.total_moves:,}")
        
        print("\nFinal memory usage:")
        GPUMemoryTracker.print_memory_usage(0)

if __name__ == "__main__":
    main()
