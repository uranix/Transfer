#include "CudaContext.h"
#include <stdio.h>
#include <math.h>

#ifndef _
#define _(x) do { \
	if ((x) != cudaSuccess) { \
	fprintf(stderr, "File %s line %d, %s failed with error `%s'\n", __FILE__, __LINE__, #x, cudaGetErrorString(cudaGetLastError())); \
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
	for (idx i = 0; i < meshdata->nP; i++) {
		for (idx j = 0; j < angdata->aslm; j++)
			if (j < angdata->slm)
				p[i * angdata->slm + j] = 
					q[i * angdata->aslm + j];
	}
	free(q);
	return p;
}

void CudaContext::computeRhs(REAL *b) {
	dim3 grid(meshdata->nPlow, meshdata->nPhigh);
	dim3 block(angdata->slm); /* no need of extra threads in block */
	rightHandSide<<<grid, block>>>(*reinterpret_cast<const DeviceMeshDataRaw *>(meshdata), *reinterpret_cast<const DeviceAngularDataRaw *>(angdata), b);
	_(/*rightHandSide*/cudaDeviceSynchronize());
}

void CudaContext::computeLhs(REAL *f, REAL *Af) {
	dim3 grid(meshdata->nPlow, meshdata->nPhigh);
	dim3 block(angdata->slm); /* no need of extra threads in block */
	volumePart<<<grid, block>>>(*reinterpret_cast<const DeviceMeshDataRaw *>(meshdata), *reinterpret_cast<const DeviceAngularDataRaw *>(angdata), f, Af);
	_(/*volumePart*/cudaDeviceSynchronize());
	surfacePart<<<grid, block>>>(*reinterpret_cast<const DeviceMeshDataRaw *>(meshdata), *reinterpret_cast<const DeviceAngularDataRaw *>(angdata), f, Af);
	_(/*surfacePart*/cudaDeviceSynchronize());
}

__global__ void addProdKern(idx nP, idx aslm, REAL *x, const REAL *y, const REAL wy) {
	int lm = threadIdx.x;
	int vertex = blockIdx.x + blockIdx.y * gridDim.x;

	int i = vertex * aslm + lm;

	if (vertex < nP)
		x[i] += wy * y[i];
}

__global__ void mulAddKern(idx nP, idx aslm, REAL *x, const REAL wx, const REAL *y) {
	int lm = threadIdx.x;
	int vertex = blockIdx.x + blockIdx.y * gridDim.x;

	int i = vertex * aslm + lm;

	if (vertex < nP)
		x[i] = wx*x[i]+y[i];
}

__global__ void mulAddProdKern(idx nP, idx aslm, REAL *x, const REAL wx, const REAL *y, const REAL wy) {
	int lm = threadIdx.x;
	int vertex = blockIdx.x + blockIdx.y * gridDim.x;

	int i = vertex * aslm + lm;

	if (vertex < nP)
		x[i] = wx*x[i]+wy*y[i];
}

/* x += wy*y */
void CudaContext::addProd(REAL *x, const REAL *y, const REAL wy) {
	dim3 grid(meshdata->nPlow, meshdata->nPhigh);
	dim3 block(angdata->slm); /* no need of extra threads in block */
	addProdKern<<<grid, block>>>(meshdata->nP, angdata->aslm, x, y, wy);
	_(/*addProd*/cudaDeviceSynchronize());
}

/* x = wx*x + y */
void CudaContext::mulAdd(REAL *x, const REAL wx, const REAL *y) {
	dim3 grid(meshdata->nPlow, meshdata->nPhigh);
	dim3 block(angdata->slm); /* no need of extra threads in block */
	mulAddKern<<<grid, block>>>(meshdata->nP, angdata->aslm, x, wx, y);
	_(/*mulAdd*/cudaDeviceSynchronize());
}

void CudaContext::mulAddProd(REAL *x, const REAL wx, const REAL *y, const REAL wy) {
	dim3 grid(meshdata->nPlow, meshdata->nPhigh);
	dim3 block(angdata->slm); /* no need of extra threads in block */
	mulAddProdKern<<<grid, block>>>(meshdata->nP, angdata->aslm, x, wx, y, wy);
	_(/*mullAddProd*/cudaDeviceSynchronize());
}

__global__ void normKern(idx nP, idx aslm, idx slm, const REAL *x, REAL *res) {
	__shared__ REAL reduce[ASLM_MAX];
	idx lm = threadIdx.x;

	reduce[lm] = 0;
	if (lm < slm) {
		for (idx i = lm, j = nP*aslm; i < j; i += aslm) {
			REAL q = x[i];
			reduce[lm] += q*q;
		}
	} 
	__syncthreads();
#pragma unroll
	for (idx s = ASLM_MAX >> 1; s > 0; s>>=1) {
		if (lm < s)
			reduce[lm] += reduce[lm + s];
		__syncthreads();
	}
	if (lm == 0)
		res[0] = reduce[0];
}

__global__ void dotKern(idx nP, idx aslm, idx slm, const REAL *x, const REAL *y, REAL *res) {
	__shared__ REAL reduce[ASLM_MAX];
	idx lm = threadIdx.x;

	reduce[lm] = 0;
	if (lm < slm) {
		for (idx i = lm, j = nP*aslm; i < j; i += aslm) {
			reduce[lm] += x[i]*y[i];
		}
	} 
	__syncthreads();
#pragma unroll
	for (idx s = ASLM_MAX >> 1; s > 0; s>>=1) {
		if (lm < s)
			reduce[lm] += reduce[lm + s];
		__syncthreads();
	}
	if (lm == 0)
		res[0] = reduce[0];
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
