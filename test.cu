#include <stdio.h>
#include <cuda.h>

__global__ void init(float *a) {
	int idx = threadIdx.x + (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x;
	a[idx] = idx;
}

__global__ void kernel(float *a) {
	int idx = threadIdx.x + (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x;
	a[idx] = sin(a[idx]);
}

int main() {
	const int N = 64*1024*1024;
	float *a;
	cudaMalloc(&a, N*sizeof(float));
	printf("Error: %s\n", cudaGetErrorString(cudaGetLastError()));
	cudaEvent_t start, end;
	cudaEventCreate(&start);
	cudaEventCreate(&end);
	printf("Error: %s\n", cudaGetErrorString(cudaGetLastError()));
	float time = 0.f;
	for (int bs = 16; bs <= 512; bs += 16) {
		dim3 grid(N/bs/1024, 1024);
		cudaEventRecord(start, 0);
		init<<< grid, bs >>>(a);
		cudaEventRecord(end, 0);
		cudaEventSynchronize(end);
		printf("Error: %s\n", cudaGetErrorString(cudaGetLastError()));
		cudaEventElapsedTime(&time, start, end);
		printf("[%10d] init : % 10.2fms\n", bs, time);
		cudaEventRecord(start, 0);
		kernel<<< grid, bs >>>(a);
		cudaEventRecord(end, 0);
		cudaEventSynchronize(end);
		cudaEventElapsedTime(&time, start, end);
		printf("[%10d] kern : % 10.2fms\n", bs, time);
	}
}
