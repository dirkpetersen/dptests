#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <curand.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <thread>
#include <iomanip>
#include <cmath>

#define CHECK_CUDA(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << " - " << cudaGetErrorString(error) << std::endl; \
            exit(1); \
        } \
    } while(0)

#define CHECK_CUBLAS(call) \
    do { \
        cublasStatus_t status = call; \
        if (status != CUBLAS_STATUS_SUCCESS) { \
            std::cerr << "cuBLAS error at " << __FILE__ << ":" << __LINE__ << std::endl; \
            exit(1); \
        } \
    } while(0)

__global__ void memoryBandwidthKernel(float* data, size_t size, int iterations) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    
    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = idx; i < size; i += stride) {
            data[i] = data[i] * 1.001f + 0.001f;
        }
    }
}

__global__ void computeIntensiveKernel(float* data, size_t size, int iterations) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    
    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = idx; i < size; i += stride) {
            float val = data[i];
            val = sinf(val) + cosf(val);
            val = expf(val * 0.1f);
            val = logf(val + 1.0f);
            val = sqrtf(val);
            data[i] = val;
        }
    }
}

class GPUBenchmark {
private:
    int numGPUs;
    std::vector<int> gpuIds;
    std::vector<cudaStream_t> streams;
    std::vector<cublasHandle_t> cublasHandles;
    std::vector<float*> deviceMemory;
    std::vector<size_t> memoryPerGPU;
    
public:
    GPUBenchmark(int num_gpus) : numGPUs(num_gpus) {
        int totalGPUs;
        CHECK_CUDA(cudaGetDeviceCount(&totalGPUs));
        
        if (numGPUs > totalGPUs) {
            std::cerr << "Requested " << numGPUs << " GPUs but only " << totalGPUs << " available" << std::endl;
            numGPUs = totalGPUs;
        }
        
        gpuIds.resize(numGPUs);
        streams.resize(numGPUs);
        cublasHandles.resize(numGPUs);
        deviceMemory.resize(numGPUs);
        memoryPerGPU.resize(numGPUs);
        
        for (int i = 0; i < numGPUs; i++) {
            gpuIds[i] = i;
            CHECK_CUDA(cudaSetDevice(i));
            
            size_t free, total;
            CHECK_CUDA(cudaMemGetInfo(&free, &total));
            memoryPerGPU[i] = free * 0.8; // Use 80% of available memory
            
            CHECK_CUDA(cudaMalloc(&deviceMemory[i], memoryPerGPU[i]));
            CHECK_CUDA(cudaStreamCreate(&streams[i]));
            CHECK_CUBLAS(cublasCreate(&cublasHandles[i]));
            CHECK_CUBLAS(cublasSetStream(cublasHandles[i], streams[i]));
            
            std::cout << "GPU " << i << ": Allocated " << memoryPerGPU[i] / (1024*1024*1024) << " GB" << std::endl;
        }
    }
    
    ~GPUBenchmark() {
        for (int i = 0; i < numGPUs; i++) {
            CHECK_CUDA(cudaSetDevice(i));
            if (deviceMemory[i]) CHECK_CUDA(cudaFree(deviceMemory[i]));
            if (streams[i]) CHECK_CUDA(cudaStreamDestroy(streams[i]));
            if (cublasHandles[i]) CHECK_CUBLAS(cublasDestroy(cublasHandles[i]));
        }
    }
    
    void initializeData() {
        std::cout << "Initializing data on all GPUs..." << std::endl;
        
        for (int i = 0; i < numGPUs; i++) {
            CHECK_CUDA(cudaSetDevice(i));
            
            curandGenerator_t gen;
            curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
            curandSetPseudoRandomGeneratorSeed(gen, time(NULL) + i);
            curandSetStream(gen, streams[i]);
            
            size_t numElements = memoryPerGPU[i] / sizeof(float);
            curandGenerateUniform(gen, deviceMemory[i], numElements);
            
            curandDestroyGenerator(gen);
        }
        
        for (int i = 0; i < numGPUs; i++) {
            CHECK_CUDA(cudaSetDevice(i));
            CHECK_CUDA(cudaStreamSynchronize(streams[i]));
        }
    }
    
    void runMemoryBandwidthTest(int durationMs) {
        std::cout << "Running memory bandwidth test..." << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        auto end = start + std::chrono::milliseconds(durationMs);
        
        int iterations = 0;
        while (std::chrono::high_resolution_clock::now() < end) {
            for (int i = 0; i < numGPUs; i++) {
                CHECK_CUDA(cudaSetDevice(i));
                
                size_t numElements = memoryPerGPU[i] / sizeof(float);
                int numBlocks = std::min(65535, (int)((numElements + 255) / 256));
                int threadsPerBlock = 256;
                
                memoryBandwidthKernel<<<numBlocks, threadsPerBlock, 0, streams[i]>>>(
                    deviceMemory[i], numElements, 10);
            }
            iterations++;
        }
        
        for (int i = 0; i < numGPUs; i++) {
            CHECK_CUDA(cudaSetDevice(i));
            CHECK_CUDA(cudaStreamSynchronize(streams[i]));
        }
        
        auto actualEnd = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(actualEnd - start);
        
        size_t totalBytes = 0;
        for (int i = 0; i < numGPUs; i++) {
            totalBytes += memoryPerGPU[i];
        }
        
        double bandwidth = (double)(totalBytes * iterations * 2) / (1024.0 * 1024.0 * 1024.0) / (duration.count() / 1000.0);
        std::cout << "Memory bandwidth: " << std::fixed << std::setprecision(2) << bandwidth << " GB/s" << std::endl;
    }
    
    void runComputeTest(int durationMs) {
        std::cout << "Running compute-intensive test..." << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        auto end = start + std::chrono::milliseconds(durationMs);
        
        int iterations = 0;
        while (std::chrono::high_resolution_clock::now() < end) {
            for (int i = 0; i < numGPUs; i++) {
                CHECK_CUDA(cudaSetDevice(i));
                
                size_t numElements = memoryPerGPU[i] / sizeof(float);
                int numBlocks = std::min(65535, (int)((numElements + 255) / 256));
                int threadsPerBlock = 256;
                
                computeIntensiveKernel<<<numBlocks, threadsPerBlock, 0, streams[i]>>>(
                    deviceMemory[i], numElements, 100);
            }
            iterations++;
        }
        
        for (int i = 0; i < numGPUs; i++) {
            CHECK_CUDA(cudaSetDevice(i));
            CHECK_CUDA(cudaStreamSynchronize(streams[i]));
        }
        
        auto actualEnd = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(actualEnd - start);
        
        size_t totalElements = 0;
        for (int i = 0; i < numGPUs; i++) {
            totalElements += memoryPerGPU[i] / sizeof(float);
        }
        
        double flops = (double)(totalElements * iterations * 100 * 6) / (duration.count() / 1000.0) / 1e12;
        std::cout << "Compute performance: " << std::fixed << std::setprecision(2) << flops << " TFLOPS" << std::endl;
    }
    
    void runMatrixMultiplication(int durationMs) {
        std::cout << "Running matrix multiplication test..." << std::endl;
        
        std::vector<float*> matrixA(numGPUs);
        std::vector<float*> matrixB(numGPUs);
        std::vector<float*> matrixC(numGPUs);
        
        for (int i = 0; i < numGPUs; i++) {
            CHECK_CUDA(cudaSetDevice(i));
            
            size_t availableElements = memoryPerGPU[i] / sizeof(float);
            int matrixSize = (int)sqrt(availableElements / 3);
            matrixSize = (matrixSize / 32) * 32; // Align to 32
            
            size_t matrixBytes = matrixSize * matrixSize * sizeof(float);
            
            CHECK_CUDA(cudaMalloc(&matrixA[i], matrixBytes));
            CHECK_CUDA(cudaMalloc(&matrixB[i], matrixBytes));
            CHECK_CUDA(cudaMalloc(&matrixC[i], matrixBytes));
            
            curandGenerator_t gen;
            curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
            curandSetPseudoRandomGeneratorSeed(gen, time(NULL) + i);
            curandSetStream(gen, streams[i]);
            
            curandGenerateUniform(gen, matrixA[i], matrixSize * matrixSize);
            curandGenerateUniform(gen, matrixB[i], matrixSize * matrixSize);
            
            curandDestroyGenerator(gen);
            
            std::cout << "GPU " << i << ": Matrix size " << matrixSize << "x" << matrixSize << std::endl;
        }
        
        auto start = std::chrono::high_resolution_clock::now();
        auto end = start + std::chrono::milliseconds(durationMs);
        
        int iterations = 0;
        while (std::chrono::high_resolution_clock::now() < end) {
            for (int i = 0; i < numGPUs; i++) {
                CHECK_CUDA(cudaSetDevice(i));
                
                size_t availableElements = memoryPerGPU[i] / sizeof(float);
                int matrixSize = (int)sqrt(availableElements / 3);
                matrixSize = (matrixSize / 32) * 32;
                
                const float alpha = 1.0f, beta = 0.0f;
                CHECK_CUBLAS(cublasSgemm(cublasHandles[i],
                    CUBLAS_OP_N, CUBLAS_OP_N,
                    matrixSize, matrixSize, matrixSize,
                    &alpha,
                    matrixA[i], matrixSize,
                    matrixB[i], matrixSize,
                    &beta,
                    matrixC[i], matrixSize));
            }
            iterations++;
        }
        
        for (int i = 0; i < numGPUs; i++) {
            CHECK_CUDA(cudaSetDevice(i));
            CHECK_CUDA(cudaStreamSynchronize(streams[i]));
        }
        
        auto actualEnd = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(actualEnd - start);
        
        double totalFlops = 0;
        for (int i = 0; i < numGPUs; i++) {
            size_t availableElements = memoryPerGPU[i] / sizeof(float);
            int matrixSize = (int)sqrt(availableElements / 3);
            matrixSize = (matrixSize / 32) * 32;
            totalFlops += 2.0 * matrixSize * matrixSize * matrixSize;
        }
        
        double tflops = totalFlops * iterations / (duration.count() / 1000.0) / 1e12;
        std::cout << "Matrix multiplication: " << std::fixed << std::setprecision(2) << tflops << " TFLOPS" << std::endl;
        
        for (int i = 0; i < numGPUs; i++) {
            CHECK_CUDA(cudaSetDevice(i));
            CHECK_CUDA(cudaFree(matrixA[i]));
            CHECK_CUDA(cudaFree(matrixB[i]));
            CHECK_CUDA(cudaFree(matrixC[i]));
        }
    }
    
    void printGPUInfo() {
        std::cout << "\n=== GPU Information ===" << std::endl;
        for (int i = 0; i < numGPUs; i++) {
            CHECK_CUDA(cudaSetDevice(i));
            
            cudaDeviceProp prop;
            CHECK_CUDA(cudaGetDeviceProperties(&prop, i));
            
            std::cout << "GPU " << i << ": " << prop.name << std::endl;
            std::cout << "  Compute Capability: " << prop.major << "." << prop.minor << std::endl;
            std::cout << "  Memory: " << prop.totalGlobalMem / (1024*1024*1024) << " GB" << std::endl;
            std::cout << "  SMs: " << prop.multiProcessorCount << std::endl;
            std::cout << "  Max Threads per SM: " << prop.maxThreadsPerMultiProcessor << std::endl;
        }
        std::cout << std::endl;
    }
};

int main(int argc, char* argv[]) {
    int numGPUs = 1;
    
    if (argc > 1) {
        numGPUs = std::atoi(argv[1]);
        if (numGPUs <= 0) {
            std::cerr << "Invalid number of GPUs: " << numGPUs << std::endl;
            return 1;
        }
    }
    
    std::cout << "=== Multi-GPU CUDA Benchmark ===" << std::endl;
    std::cout << "Requested GPUs: " << numGPUs << std::endl;
    
    try {
        GPUBenchmark benchmark(numGPUs);
        benchmark.printGPUInfo();
        benchmark.initializeData();
        
        const int testDuration = 15000; // 15 seconds per test
        
        benchmark.runMemoryBandwidthTest(testDuration);
        benchmark.runComputeTest(testDuration);
        benchmark.runMatrixMultiplication(testDuration);
        
        std::cout << "\n=== Benchmark Complete ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}