NVCC = nvcc
CUDA_FLAGS = -O3 -std=c++14 -arch=sm_87 -Xcompiler -fPIC
LIBS = -lcublas -lcurand

TARGET = gpu_benchmark
SOURCE = gpu_benchmark.cu

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(NVCC) $(CUDA_FLAGS) -o $(TARGET) $(SOURCE) $(LIBS)

clean:
	rm -f $(TARGET)

test: $(TARGET)
	./$(TARGET) 1
	./$(TARGET) 2
	./$(TARGET) 4

.PHONY: all clean test