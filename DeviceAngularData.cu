#include "AngularData.h"

#include "common.cuh"

#include <stdio.h>

#ifndef _
#define _(x) do { \
	if ((x) != cudaSuccess) { \
	fprintf(stderr, "File %s line %d, %s failed with error `%s'\n", __FILE__, __LINE__, #x, cudaGetErrorString(cudaGetLastError())); \
	fflush(stderr); } \
} while (0)
#endif

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
	_(cudaMalloc((void **)&omega, 3*aslm*aslm*sizeof(REAL)));
	_(cudaMalloc((void **)&omega_pos, 3*aslm*aslm*sizeof(idx)));
	_(cudaMalloc((void **)&Ox, aslm*aslm*sizeof(REAL)));
	_(cudaMalloc((void **)&Oy, aslm*aslm*sizeof(REAL)));
	_(cudaMalloc((void **)&Oz, aslm*aslm*sizeof(REAL)));

	printf("DeviceAngularData:\n");
	printf("\tomega     = %p\n", omega);
	printf("\tomega_pos = %p\n", omega_pos);
	printf("\tOx        = %p\n", Ox);
	printf("\tOy        = %p\n", Oy);
	printf("\tOz        = %p\n", Oz);

	void *tmp;
	dim3 block;
	_(cudaMalloc((void **)&tmp, 3*aslm*aslm*sizeof(REAL)));

	_(cudaMemcpy(tmp, host.omega, 3*slm*slm*sizeof(REAL), cudaMemcpyHostToDevice));
	copy_omega<<<aslm,aslm>>>(omega, (REAL *)tmp, slm);

	_(cudaMemcpy(tmp, host.omega_pos, 3*slm*slm*sizeof(idx), cudaMemcpyHostToDevice));
	copy_omega_pos<<<aslm,aslm>>>(omega_pos, (idx *)tmp, slm);

	_(cudaMemcpy(tmp, host.Ox, slm*slm*sizeof(REAL), cudaMemcpyHostToDevice));
	copy_On<<<aslm,aslm>>>(Ox, (REAL *)tmp, slm);

	_(cudaMemcpy(tmp, host.Oy, slm*slm*sizeof(REAL), cudaMemcpyHostToDevice));
	copy_On<<<aslm,aslm>>>(Oy, (REAL *)tmp, slm);

	_(cudaMemcpy(tmp, host.Oz, slm*slm*sizeof(REAL), cudaMemcpyHostToDevice));
	copy_On<<<aslm,aslm>>>(Oz, (REAL *)tmp, slm);

	_(cudaFree(tmp));
}

DeviceAngularData::~DeviceAngularData() {
	_(cudaFree(omega));
	_(cudaFree(omega_pos));
	_(cudaFree(Ox));
	_(cudaFree(Oy));
	_(cudaFree(Oz));
}
