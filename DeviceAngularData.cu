#include "AngularData.h"

#include "common.cuh"

#include <stdio.h>

__global__ void copy_omega(REAL *dst, REAL *src, idx slm) {
	idx aslm = blockDim.x;
	#pragma unroll
	for (int s = 0; s < 3; s++) {
		if ((threadIdx.x < slm) && (blockIdx.x < slm))
			dst[s*aslm*aslm + aslm * blockIdx.x + threadIdx.x] = src[s*slm*slm + slm * blockIdx.x + threadIdx.x];
		else
			dst[s*aslm*aslm + aslm * blockIdx.x + threadIdx.x] = 0;
	}
}

__global__ void copy_omega_pos(idx *dst, idx *src, idx slm) {
	idx aslm = blockDim.x;
	#pragma unroll
	for (int s = 0; s < 3; s++) {
		if ((threadIdx.x < slm) && (blockIdx.x < slm))
			dst[s*aslm*aslm + aslm * blockIdx.x + threadIdx.x] = src[s*slm*slm + slm * blockIdx.x + threadIdx.x];
		else
			dst[s*aslm*aslm + aslm * blockIdx.x + threadIdx.x] = 0;
	}
}

__global__ void copy_On(REAL *dst, REAL *src, idx slm) {
	idx aslm = blockDim.x;
	if ((threadIdx.x < slm) && (blockIdx.x < slm))
		dst[aslm * blockIdx.x + threadIdx.x] = src[slm * blockIdx.x + threadIdx.x];
	else
		dst[aslm * blockIdx.x + threadIdx.x] = 0;
}

DeviceAngularData::DeviceAngularData(const AngularData &host) {
	slm = host.slm;
	aslm = align_power(host.slm, COALESCED_NUM(REAL));
	cudaMalloc((void **)&omega, 3*aslm*aslm*sizeof(REAL));
	cudaMalloc((void **)&omega_pos, 3*aslm*aslm*sizeof(idx));
	cudaMalloc((void **)&Ox, aslm*aslm*sizeof(REAL));
	cudaMalloc((void **)&Oy, aslm*aslm*sizeof(REAL));
	cudaMalloc((void **)&Oz, aslm*aslm*sizeof(REAL));

	printf("DeviceAngularData:\n");
	printf("\tomega     = %p\n", omega);
	printf("\tomega_pos = %p\n", omega_pos);
	printf("\tOx        = %p\n", Ox);
	printf("\tOy        = %p\n", Oy);
	printf("\tOz        = %p\n", Oz);

	void *tmp;
	dim3 block;
	cudaMalloc((void **)&tmp, 3*aslm*aslm*sizeof(REAL));

	cudaMemcpy(tmp, host.omega, 3*slm*slm*sizeof(REAL), cudaMemcpyHostToDevice);
	copy_omega<<<aslm,aslm>>>(omega, (REAL *)tmp, slm);

	cudaMemcpy(tmp, host.omega_pos, 3*slm*slm*sizeof(idx), cudaMemcpyHostToDevice);
	copy_omega_pos<<<aslm,aslm>>>(omega_pos, (idx *)tmp, slm);

	cudaMemcpy(tmp, host.Ox, slm*slm*sizeof(REAL), cudaMemcpyHostToDevice);
	copy_On<<<aslm,aslm>>>(Ox, (REAL *)tmp, slm);

	cudaMemcpy(tmp, host.Oy, slm*slm*sizeof(REAL), cudaMemcpyHostToDevice);
	copy_On<<<aslm,aslm>>>(Oy, (REAL *)tmp, slm);

	cudaMemcpy(tmp, host.Oz, slm*slm*sizeof(REAL), cudaMemcpyHostToDevice);
	copy_On<<<aslm,aslm>>>(Oz, (REAL *)tmp, slm);

	cudaFree(tmp);
}

DeviceAngularData::~DeviceAngularData() {
	cudaFree(omega);
	cudaFree(omega_pos);
	cudaFree(Ox);
	cudaFree(Oy);
	cudaFree(Oz);
}
