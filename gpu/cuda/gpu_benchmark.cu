#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cufft.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <memory>
#include <iomanip>

#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << " - " << cudaGetErrorString(error) << std::endl; \
            exit(1); \
        } \
    } while(0)

#define CUBLAS_CHECK(call) \
    do { \
        cublasStatus_t status = call; \
        if (status != CUBLAS_STATUS_SUCCESS) { \
            std::cerr << "cuBLAS error at " << __FILE__ << ":" << __LINE__ << std::endl; \
            exit(1); \
        } \
    } while(0)

class GPUBenchmark {
private:
    int num_gpus;
    std::vector<cudaStream_t> streams;
    std::vector<cublasHandle_t> cublas_handles;
    std::vector<float*> d_matrices_a, d_matrices_b, d_matrices_c;
    std::vector<float*> d_memory_test;
    std::vector<cufftComplex*> d_fft_data;
    std::vector<cufftHandle> fft_plans;
    
    static const int MATRIX_SIZE = 8192;
    static const int FFT_SIZE = 1048576;
    static const size_t MEMORY_SIZE = 16ULL * 1024 * 1024 * 1024; // 16GB per GPU
    
public:
    GPUBenchmark(int gpus) : num_gpus(gpus) {
        int available_gpus;
        cudaGetDeviceCount(&available_gpus);
        
        if (num_gpus > available_gpus) {
            std::cerr << "Requested " << num_gpus << " GPUs but only " << available_gpus << " available" << std::endl;
            exit(1);
        }
        
        std::cout << "Initializing benchmark on " << num_gpus << " GPUs..." << std::endl;
        
        streams.resize(num_gpus);
        cublas_handles.resize(num_gpus);
        d_matrices_a.resize(num_gpus);
        d_matrices_b.resize(num_gpus);
        d_matrices_c.resize(num_gpus);
        d_memory_test.resize(num_gpus);
        d_fft_data.resize(num_gpus);
        fft_plans.resize(num_gpus);
        
        for (int i = 0; i < num_gpus; i++) {
            CUDA_CHECK(cudaSetDevice(i));
            
            // Create streams
            CUDA_CHECK(cudaStreamCreate(&streams[i]));
            
            // Create cuBLAS handles
            CUBLAS_CHECK(cublasCreate(&cublas_handles[i]));
            CUBLAS_CHECK(cublasSetStream(cublas_handles[i], streams[i]));
            
            // Allocate matrices for GEMM
            size_t matrix_bytes = MATRIX_SIZE * MATRIX_SIZE * sizeof(float);
            CUDA_CHECK(cudaMalloc(&d_matrices_a[i], matrix_bytes));
            CUDA_CHECK(cudaMalloc(&d_matrices_b[i], matrix_bytes));
            CUDA_CHECK(cudaMalloc(&d_matrices_c[i], matrix_bytes));
            
            // Allocate memory test buffer
            CUDA_CHECK(cudaMalloc(&d_memory_test[i], MEMORY_SIZE));
            
            // Allocate FFT data
            CUDA_CHECK(cudaMalloc(&d_fft_data[i], FFT_SIZE * sizeof(cufftComplex)));
            
            // Create FFT plan
            cufftPlan1d(&fft_plans[i], FFT_SIZE, CUFFT_C2C, 1);
            cufftSetStream(fft_plans[i], streams[i]);
            
            // Initialize matrices with random data
            initializeMatrices(i);
            
            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, i);
            std::cout << "GPU " << i << ": " << prop.name << " (" 
                      << prop.totalGlobalMem / (1024*1024*1024) << " GB)" << std::endl;
        }
    }
    
    ~GPUBenchmark() {
        for (int i = 0; i < num_gpus; i++) {
            CUDA_CHECK(cudaSetDevice(i));
            
            cudaStreamDestroy(streams[i]);
            cublasDestroy(cublas_handles[i]);
            cufftDestroy(fft_plans[i]);
            
            cudaFree(d_matrices_a[i]);
            cudaFree(d_matrices_b[i]);
            cudaFree(d_matrices_c[i]);
            cudaFree(d_memory_test[i]);
            cudaFree(d_fft_data[i]);
        }
    }
    
    void initializeMatrices(int gpu_id) {
        CUDA_CHECK(cudaSetDevice(gpu_id));
        
        size_t matrix_elements = MATRIX_SIZE * MATRIX_SIZE;
        size_t matrix_bytes = matrix_elements * sizeof(float);
        
        std::vector<float> h_matrix(matrix_elements);
        for (size_t i = 0; i < matrix_elements; i++) {
            h_matrix[i] = static_cast<float>(rand()) / RAND_MAX;
        }
        
        CUDA_CHECK(cudaMemcpyAsync(d_matrices_a[gpu_id], h_matrix.data(), matrix_bytes, 
                                  cudaMemcpyHostToDevice, streams[gpu_id]));
        CUDA_CHECK(cudaMemcpyAsync(d_matrices_b[gpu_id], h_matrix.data(), matrix_bytes, 
                                  cudaMemcpyHostToDevice, streams[gpu_id]));
        
        // Initialize FFT data
        std::vector<cufftComplex> h_fft_data(FFT_SIZE);
        for (int i = 0; i < FFT_SIZE; i++) {
            h_fft_data[i].x = static_cast<float>(rand()) / RAND_MAX;
            h_fft_data[i].y = static_cast<float>(rand()) / RAND_MAX;
        }
        CUDA_CHECK(cudaMemcpyAsync(d_fft_data[gpu_id], h_fft_data.data(), 
                                  FFT_SIZE * sizeof(cufftComplex), 
                                  cudaMemcpyHostToDevice, streams[gpu_id]));
    }
    
    double runMatrixMultiply(int iterations = 100) {
        std::cout << "Running matrix multiplication benchmark..." << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        const float alpha = 1.0f, beta = 0.0f;
        
        for (int iter = 0; iter < iterations; iter++) {
            for (int i = 0; i < num_gpus; i++) {
                CUDA_CHECK(cudaSetDevice(i));
                CUBLAS_CHECK(cublasSgemm(cublas_handles[i], CUBLAS_OP_N, CUBLAS_OP_N,
                                       MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE,
                                       &alpha,
                                       d_matrices_a[i], MATRIX_SIZE,
                                       d_matrices_b[i], MATRIX_SIZE,
                                       &beta,
                                       d_matrices_c[i], MATRIX_SIZE));
            }
        }
        
        for (int i = 0; i < num_gpus; i++) {
            CUDA_CHECK(cudaSetDevice(i));
            CUDA_CHECK(cudaStreamSynchronize(streams[i]));
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(end - start).count();
        
        double flops = static_cast<double>(num_gpus) * iterations * 2.0 * 
                      MATRIX_SIZE * MATRIX_SIZE * MATRIX_SIZE;
        double tflops = flops / elapsed / 1e12;
        
        std::cout << "Matrix Multiply: " << std::fixed << std::setprecision(2) 
                  << tflops << " TFLOPS" << std::endl;
        
        return elapsed;
    }
    
    double runFFTBenchmark(int iterations = 1000) {
        std::cout << "Running FFT benchmark..." << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        for (int iter = 0; iter < iterations; iter++) {
            for (int i = 0; i < num_gpus; i++) {
                CUDA_CHECK(cudaSetDevice(i));
                cufftExecC2C(fft_plans[i], d_fft_data[i], d_fft_data[i], CUFFT_FORWARD);
            }
        }
        
        for (int i = 0; i < num_gpus; i++) {
            CUDA_CHECK(cudaSetDevice(i));
            CUDA_CHECK(cudaStreamSynchronize(streams[i]));
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(end - start).count();
        
        std::cout << "FFT Performance: " << std::fixed << std::setprecision(2)
                  << (num_gpus * iterations) / elapsed << " FFTs/sec" << std::endl;
        
        return elapsed;
    }
    
    double runMemoryBandwidthTest() {
        std::cout << "Running memory bandwidth test..." << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        const int iterations = 10;
        for (int iter = 0; iter < iterations; iter++) {
            for (int i = 0; i < num_gpus; i++) {
                CUDA_CHECK(cudaSetDevice(i));
                CUDA_CHECK(cudaMemsetAsync(d_memory_test[i], iter % 256, MEMORY_SIZE, streams[i]));
            }
        }
        
        for (int i = 0; i < num_gpus; i++) {
            CUDA_CHECK(cudaSetDevice(i));
            CUDA_CHECK(cudaStreamSynchronize(streams[i]));
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(end - start).count();
        
        double bytes_processed = static_cast<double>(num_gpus) * iterations * MEMORY_SIZE;
        double bandwidth_gbps = bytes_processed / elapsed / 1e9;
        
        std::cout << "Memory Bandwidth: " << std::fixed << std::setprecision(2)
                  << bandwidth_gbps << " GB/s" << std::endl;
        
        return elapsed;
    }
    
    void runComprehensiveBenchmark(int duration_seconds = 60) {
        std::cout << "\n=== GPU Benchmark Started ===" << std::endl;
        std::cout << "Target duration: " << duration_seconds << " seconds" << std::endl;
        std::cout << "Using " << num_gpus << " GPU(s)" << std::endl;
        
        auto benchmark_start = std::chrono::high_resolution_clock::now();
        double total_elapsed = 0;
        
        while (total_elapsed < duration_seconds) {
            auto cycle_start = std::chrono::high_resolution_clock::now();
            
            runMatrixMultiply(20);
            runFFTBenchmark(100);
            runMemoryBandwidthTest();
            
            auto cycle_end = std::chrono::high_resolution_clock::now();
            total_elapsed = std::chrono::duration<double>(cycle_end - benchmark_start).count();
            
            std::cout << "Elapsed: " << std::fixed << std::setprecision(1) 
                      << total_elapsed << "s" << std::endl;
        }
        
        std::cout << "\n=== Benchmark Complete ===" << std::endl;
        std::cout << "Total runtime: " << std::fixed << std::setprecision(2) 
                  << total_elapsed << " seconds" << std::endl;
    }
};

int main(int argc, char* argv[]) {
    int num_gpus = 1;
    int duration = 60;
    
    if (argc > 1) {
        num_gpus = std::atoi(argv[1]);
        if (num_gpus <= 0) {
            std::cerr << "Number of GPUs must be positive" << std::endl;
            return 1;
        }
    }
    
    if (argc > 2) {
        duration = std::atoi(argv[2]);
        if (duration <= 0) {
            std::cerr << "Duration must be positive" << std::endl;
            return 1;
        }
    }
    
    std::cout << "GPU Comprehensive Benchmark" << std::endl;
    std::cout << "Requested GPUs: " << num_gpus << std::endl;
    std::cout << "Duration: " << duration << " seconds" << std::endl;
    
    try {
        GPUBenchmark benchmark(num_gpus);
        benchmark.runComprehensiveBenchmark(duration);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}