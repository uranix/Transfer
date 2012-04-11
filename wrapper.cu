#include "kernels.h"
#include "kernels.cu"

void *deviceAlloc(size_t size) {
	void *ret;
	cudaMalloc(&ret, size);
	return ret;
}

void deviceFree(void *mem) {
	cudaFree(mem);
}

void computeRhs(const DeviceMeshDataRaw meshdata, const DeviceAngularDataRaw angdata, REAL *f, REAL *r) {
	dim3 grid(meshdata.nPlow, meshdata.nPhigh);
	dim3 block(angdata.aslm);
	volumePart<<<grid, block>>>(meshdata, angdata, f, r);
	surfacePart<<<grid, block>>>(meshdata, angdata, f, r);
}

