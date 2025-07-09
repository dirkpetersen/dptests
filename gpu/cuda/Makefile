NVCC = nvcc
CUDA_PATH = /usr/local/cuda
CFLAGS = -O3 -std=c++14
NVCCFLAGS = -O3 -std=c++14 -arch=sm_90 -lcublas -lcufft
INCLUDES = -I$(CUDA_PATH)/include
LIBS = -L$(CUDA_PATH)/lib64 -lcudart -lcublas -lcufft

TARGET = gpu_benchmark
SOURCE = gpu_benchmark.cu

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) -o $(TARGET) $(SOURCE) $(LIBS)

clean:
	rm -f $(TARGET)

test: $(TARGET)
	./$(TARGET) 1 10

benchmark-1gpu: $(TARGET)
	./$(TARGET) 1 60

benchmark-2gpu: $(TARGET)
	./$(TARGET) 2 60

benchmark-4gpu: $(TARGET)
	./$(TARGET) 4 60

benchmark-8gpu: $(TARGET)
	./$(TARGET) 8 60

.PHONY: all clean test benchmark-1gpu benchmark-2gpu benchmark-4gpu benchmark-8gpu