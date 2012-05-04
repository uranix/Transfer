#include "CudaContext.h"
#include <stdio.h>
#include <math.h>

#ifndef _
#define _(x) do { \
	if ((x) != cudaSuccess) { \
	fprintf(stderr, "%s:%d: %s failed with error `%s'\n", __FILE__, __LINE__, \
		#x, cudaGetErrorString(cudaGetLastError())); \
	fflush(stderr); } \
} while (0)
#endif

#include "kernels.cu"

void CudaContext::setDevice(int dev) {
	_(cudaSetDevice(dev));
}

void *CudaContext::deviceAlloc(size_t size) {
	void *ret;
	_(cudaMalloc(&ret, size));
	_(cudaMemset(ret, 0, size));
	return ret;
}

void CudaContext::deviceFree(void *mem) {
	_(cudaFree(mem));
}

void CudaContext::copyToDev(void *dst, void *src, size_t sz) {
	_(cudaMemcpy(dst, src, sz, cudaMemcpyHostToDevice));
}

void CudaContext::copyToHost(void *dst, void *src, size_t sz) {
	_(cudaMemcpy(dst, src, sz, cudaMemcpyDeviceToHost));
}

REAL *CudaContext::getHostSmall(void *src) {
	REAL *p = (REAL *)malloc(angdata->slm * meshdata->nP * sizeof(REAL)), 
		 *q = (REAL *)malloc(angdata->aslm * meshdata->nP * sizeof(REAL));
	_(cudaMemcpy(q, src, angdata->aslm * meshdata->nP * sizeof(REAL), cudaMemcpyDeviceToHost));
	for (idx i = 0; i < meshdata->nP; i++)
		for (idx j = 0; j < angdata->aslm; j++)
			if (j < angdata->slm)
				p[i * angdata->slm + j] = 
					q[i * angdata->aslm + j];
	
	free(q);
	return p;
}

void CudaContext::computeRhs(REAL *b) {
	dim3 grid(meshdata->nPlow, meshdata->nPhigh);
	dim3 block(angdata->aslm); 
	rightHandSide<<<grid, block>>>(
			*reinterpret_cast<const DeviceMeshDataRaw *>(meshdata), 
			*reinterpret_cast<const DeviceAngularDataRaw *>(angdata), b);
	_(/*rightHandSide*/cudaDeviceSynchronize());
}

void CudaContext::computeLhs(REAL *f, REAL *Af) {
	dim3 grid(meshdata->nPlow, meshdata->nPhigh);
	dim3 block(angdata->aslm); 
	volumePart<<<grid, block>>>(
			*reinterpret_cast<const DeviceMeshDataRaw *>(meshdata), 
			*reinterpret_cast<const DeviceAngularDataRaw *>(angdata), f, Af);
	_(/*volumePart*/cudaDeviceSynchronize());
	surfacePart<<<grid, block>>>(
			*reinterpret_cast<const DeviceMeshDataRaw *>(meshdata), 
			*reinterpret_cast<const DeviceAngularDataRaw *>(angdata), f, Af);
	_(/*surfacePart*/cudaDeviceSynchronize());
}

/* x += wy*y */
void CudaContext::addProd(REAL *x, const REAL *y, const REAL wy) {
	dim3 grid(meshdata->nPlow, meshdata->nPhigh);
	dim3 block(angdata->aslm); 
	addProdKern<<<grid, block>>>(meshdata->nP, angdata->aslm, x, y, wy);
	_(/*addProd*/cudaDeviceSynchronize());
}

/* x = wx*x + y */
void CudaContext::mulAdd(REAL *x, const REAL wx, const REAL *y) {
	dim3 grid(meshdata->nPlow, meshdata->nPhigh);
	dim3 block(angdata->aslm); 
	mulAddKern<<<grid, block>>>(meshdata->nP, angdata->aslm, x, wx, y);
	_(/*mulAdd*/cudaDeviceSynchronize());
}

void CudaContext::mulAddProd(REAL *x, const REAL wx, const REAL *y, const REAL wy) {
	dim3 grid(meshdata->nPlow, meshdata->nPhigh);
	dim3 block(angdata->aslm); 
	mulAddProdKern<<<grid, block>>>(meshdata->nP, angdata->aslm, x, wx, y, wy);
	_(/*mullAddProd*/cudaDeviceSynchronize());
}

REAL CudaContext::norm(const REAL *x) {
	dim3 grid(1, 1);
	dim3 block(ASLM_MAX);
	normKern<<<grid, block>>>(meshdata->nP, angdata->aslm, angdata->slm, x, red);
	REAL hred;
	copyToHost(&hred, red, sizeof(REAL));
	hred /= meshdata->nP * angdata->slm;
	return sqrt(hred);
}	

REAL CudaContext::dot(const REAL *x, const REAL *y) {
	dim3 grid(1, 1);
	dim3 block(ASLM_MAX);
	dotKern<<<grid, block>>>(meshdata->nP, angdata->aslm, angdata->slm, x, y, red);
	REAL hred;
	copyToHost(&hred, red, sizeof(REAL));
	return hred;
}	
