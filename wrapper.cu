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

void copyToDev(void *dst, void *src, size_t sz) {
	cudaMemcpy(dst, src, sz, cudaMemcpyHostToDevice);
}

void copyToHost(void *dst, void *src, size_t sz) {
	cudaMemcpy(dst, src, sz, cudaMemcpyDeviceToHost);
}

void computeRhs(const DeviceMeshData *meshdata, const DeviceAngularData *angdata, REAL *f, REAL *Af, REAL *b) {
	dim3 grid(meshdata->nPlow, meshdata->nPhigh);
	dim3 block(angdata->aslm);
	volumePart<<<grid, block>>>(*reinterpret_cast<const DeviceMeshDataRaw *>(meshdata), *reinterpret_cast<const DeviceAngularDataRaw *>(angdata), f, Af);
	surfacePart<<<grid, block>>>(*reinterpret_cast<const DeviceMeshDataRaw *>(meshdata), *reinterpret_cast<const DeviceAngularDataRaw *>(angdata), f, Af);
	rightHandSide<<<grid, block>>>(*reinterpret_cast<const DeviceMeshDataRaw *>(meshdata), *reinterpret_cast<const DeviceAngularDataRaw *>(angdata), b);
}

