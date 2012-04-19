#include "kernels.h"
#include "kernels.cu"
#include <stdio.h>

#ifndef _
#define _(x) do { \
	fprintf(stderr, "File %s line %d, %s is going to fail ...\n", __FILE__, __LINE__, #x); \
	fflush(stderr); \
	(x); \
	fprintf(stderr, "File %s line %d, %s failed with error `%s'\n", __FILE__, __LINE__, #x, cudaGetErrorString(cudaGetLastError())); \
	fflush(stderr); \
} while (0)
#endif

void *deviceAlloc(size_t size) {
	void *ret;
	_(cudaMalloc(&ret, size));
	return ret;
}

void deviceFree(void *mem) {
	_(cudaFree(mem));
}

void copyToDev(void *dst, void *src, size_t sz) {
	_(cudaMemcpy(dst, src, sz, cudaMemcpyHostToDevice));
}

void copyToHost(void *dst, void *src, size_t sz) {
	_(cudaMemcpy(dst, src, sz, cudaMemcpyDeviceToHost));
}

void computeRhs(const DeviceMeshData *meshdata, const DeviceAngularData *angdata, REAL *f, REAL *Af, REAL *b) {
	dim3 grid(meshdata->nPlow, meshdata->nPhigh);
	dim3 block(angdata->slm); /* no need of extra threads in block */
	volumePart<<<grid, block>>>(*reinterpret_cast<const DeviceMeshDataRaw *>(meshdata), *reinterpret_cast<const DeviceAngularDataRaw *>(angdata), f, Af);
	_(cudaDeviceSynchronize());
	surfacePart<<<grid, block>>>(*reinterpret_cast<const DeviceMeshDataRaw *>(meshdata), *reinterpret_cast<const DeviceAngularDataRaw *>(angdata), f, Af);
	_(cudaDeviceSynchronize());
	rightHandSide<<<grid, block>>>(*reinterpret_cast<const DeviceMeshDataRaw *>(meshdata), *reinterpret_cast<const DeviceAngularDataRaw *>(angdata), b);
	_(cudaDeviceSynchronize());
}

