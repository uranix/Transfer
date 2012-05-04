#include "CudaContext.h"
#include <stdio.h>
#include <math.h>

#ifndef _
#define _(x) do { \
	CUresult __res = (x); \
	if (__res != CUDA_SUCCESS) { \
	fprintf(stderr, "%s:%d: %s failed with error `%d'\n", __FILE__, __LINE__, \
		#x, __res); \
	fflush(stderr); } \
} while (0)
#endif

#include "common.cuh"

/*
#include "kernels.cu"
*/

#define LOAD_KERNEL(x) \
	_(cuModuleGetFunction(&this->x, _mod, #x))

CudaContext::CudaContext(const int dev, const MeshData &dmd, const AngularData &dad) {
	CUdevice handle;
	_(cuInit(0));
	_(cuDeviceGet(&handle, dev));
	_(cuCtxCreate(&this->_ctx, CU_CTX_SCHED_AUTO, handle));
	int major, minor;
	_(cuDeviceComputeCapability(&major, &minor, handle));
	char fn[1024];
	sprintf(fn, "kernels.sm_%d%d.cubin", major, minor);
	printf("Loading %s module\n", fn);
	_(cuModuleLoad(&_mod, fn));
	LOAD_KERNEL(addProdKern);
	LOAD_KERNEL(mulAddKern);
	LOAD_KERNEL(mulAddProdKern);
	LOAD_KERNEL(normKern);
	LOAD_KERNEL(dotKern);
	LOAD_KERNEL(rightHandSide);
	LOAD_KERNEL(volumePart);
	LOAD_KERNEL(surfacePart);
	
	meshdata = new DeviceMeshData(this, dmd);
	angdata = new DeviceAngularData(this, dad);
	red = (REAL *)deviceAlloc(sizeof(REAL));
}

CudaContext::~CudaContext(){
	deviceFree(red);
	delete meshdata;
	delete angdata;
	_(cuModuleUnload(_mod));
	_(cuCtxDetach(this->_ctx));

}

void *CudaContext::deviceAlloc(size_t size) const {
	CUdeviceptr ret;
	_(cuMemAlloc(&ret, size));
	_(cuMemsetD32(ret, 0, size >> 2));
	return (void *)ret;
}

idx CudaContext::N() const {
	return meshdata->nP * angdata->aslm;
}

REAL *CudaContext::allocVector() const {
	return (REAL *)deviceAlloc(N() * sizeof(REAL));
}

void CudaContext::deviceFree(void *mem) const {
	_(cuMemFree((CUdeviceptr)mem));
}

void CudaContext::copyToDev(void *dst, void *src, size_t sz) const {
	_(cuMemcpyHtoD((CUdeviceptr)dst, src, sz));
}

void CudaContext::copyToHost(void *dst, void *src, size_t sz) const {
	_(cuMemcpyDtoH(dst, (CUdeviceptr)src, sz));
}

REAL *CudaContext::getHostSmall(void *src) const {
	REAL *p = (REAL *)malloc(angdata->slm * meshdata->nP * sizeof(REAL)), 
		 *q = (REAL *)malloc(angdata->aslm * meshdata->nP * sizeof(REAL));
	copyToHost(q, src, N() * sizeof(REAL));
	for (idx i = 0; i < meshdata->nP; i++)
		for (idx j = 0; j < angdata->aslm; j++)
			if (j < angdata->slm)
				p[i * angdata->slm + j] = q[i * angdata->aslm + j];
	
	free(q);
	return p;
}
void CudaContext::computeRhs(REAL *b) const {
	const void *params[3] = {	meshdata, angdata, &b };
	_(cuLaunchKernel(rightHandSide, 
		meshdata->nPlow, meshdata->nPhigh, 1,
		angdata->aslm, 1, 1,
		0, 0,
		const_cast<void **>(params), 0));
}

void CudaContext::computeLhs(REAL *f, REAL *Af) const {
	const void *params[4] = { meshdata, angdata, &f, &Af };
	_(cuLaunchKernel(volumePart, 
		meshdata->nPlow, meshdata->nPhigh, 1,
		angdata->aslm, 1, 1,
		0, 0,
		const_cast<void **>(params), 0));
	_(cuLaunchKernel(surfacePart, 
		meshdata->nPlow, meshdata->nPhigh, 1,
		angdata->aslm, 1, 1,
		0, 0,
		const_cast<void **>(params), 0));
}

/* x += wy*y */
void CudaContext::addProd(REAL *x, const REAL *y, const REAL wy) const {
	const void *params[5] = { &meshdata->nP, &angdata->aslm, &x, &y, &wy };
	_(cuLaunchKernel(addProdKern, 
		meshdata->nPlow, meshdata->nPhigh, 1,
		angdata->aslm, 1, 1,
		0, 0,
		const_cast<void **>(params), 0));
}

/* x = wx*x + y */
void CudaContext::mulAdd(REAL *x, const REAL wx, const REAL *y) const {
	const void *params[5] = { &meshdata->nP, &angdata->aslm, &x, &wx, &y };
	_(cuLaunchKernel(mulAddKern, 
		meshdata->nPlow, meshdata->nPhigh, 1,
		angdata->aslm, 1, 1,
		0, 0,
		const_cast<void **>(params), 0));
}

void CudaContext::mulAddProd(REAL *x, const REAL wx, const REAL *y, const REAL wy) const {
	const void *params[6] = { &meshdata->nP, &angdata->aslm, &x, &wx, &y, &wy };
	_(cuLaunchKernel(mulAddProdKern, 
		meshdata->nPlow, meshdata->nPhigh, 1,
		angdata->aslm, 1, 1,
		0, 0,
		const_cast<void **>(params), 0));
}

REAL CudaContext::norm(const REAL *x) const {
	const void *params[5] = { &meshdata->nP, &angdata->aslm, &angdata->slm, &x, &red };
	_(cuLaunchKernel(normKern, 
		1, 1, 1,
		ASLM_MAX, 1, 1,
		0, 0,
		const_cast<void **>(params), 0));
	REAL hred;
	copyToHost(&hred, red, sizeof(REAL));
	hred /= meshdata->nP * angdata->slm;
	return sqrt(hred);
}	

REAL CudaContext::dot(const REAL *x, const REAL *y) const {
	const void *params[6] = { &meshdata->nP, &angdata->aslm, &angdata->slm, &x, &y, &red };
	_(cuLaunchKernel(dotKern, 
		1, 1, 1,
		ASLM_MAX, 1, 1,
		0, 0,
		const_cast<void **>(params), 0));
	REAL hred;
	copyToHost(&hred, red, sizeof(REAL));
	return hred;
}
