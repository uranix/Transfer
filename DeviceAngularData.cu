#include "AngularData.h"

#include "common.cuh"

__global__ copy_omega(REAL *dst, REAL *src, idx slm) {
	idx aslm = blockDim.x;
	#pragma unroll
	for (int s = 0; s < 3; s++) {
		if ((threadIdx.x < slm) && (blockIdx.x < slm))
			dst[s*aslm*aslm + aslm * blockIdx.x + threadIdx.x] = src[s*slm*slm + slm * blockIdx.x + threadIdx.x];
		else
			dst[s*aslm*aslm + aslm * blockIdx.x + threadIdx.x] = 0;
	}
}

DeviceAngularData::DeviceAngularData(const AngularData &host) {
	slm = align_power(host.slm, COALESCED_NUM(REAL));
	cudaMalloc(omega, 3*slm*slm*sizeof(REAL));
	cudaMalloc(omega_pos, 3*slm*slm*sizeof(idx));
	cudaMalloc(Ox, slm*slm*sizeof(REAL));
	cudaMalloc(Oy, slm*slm*sizeof(REAL));
	cudaMalloc(Oz, slm*slm*sizeof(REAL));

	void *tmp;
	dim3 block;
	cudaMalloc(tmp, 3*slm*slm*sizeof(REAL));

	cudaMemcpy(tmp, host.omega, cudaMemcpyHostToDevice);
	copy_omega<<<slm,slm>>>(omega, (REAL *)tmp, host.slm);

	cudaMemcpy(tmp, host.omega_pos, cudaMemcpyHostToDevice);
	copy_omega_pos<<<slm,slm>>>(omega, (idx *)tmp, host.slm);

	cudaMemcpy(tmp, host.Ox, cudaMemcpyHostToDevice);
	copy_On<<<slm,slm>>>(Ox, (REAL *)tmp, host.slm);

	cudaMemcpy(tmp, host.Oy, cudaMemcpyHostToDevice);
	copy_On<<<slm,slm>>>(Oy, (REAL *)tmp, host.slm);

	cudaMemcpy(tmp, host.Oz, cudaMemcpyHostToDevice);
	copy_On<<<slm,slm>>>(Oz, (REAL *)tmp, host.slm);

	cudaFree(tmp);
}

DeviceAngularData::~DeviceAngularData() {
	cudaFree(omega);
	cudaFree(omega_pos);
	cudaFree(Ox);
	cudaFree(Oy);
	cudaFree(Oz);
}
