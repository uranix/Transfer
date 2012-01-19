#include "kernels.h"

#include "kernels.cu"

#ifdef __cplusplus
extern "C" {
#endif

void *deviceAlloc(size_t size) {
	void *ret;
	cudaMalloc(&ret, size);
	return ret;
}

void deviceFree(void *mem) {
	cudaFree(mem);
}

void computeRhs(MeshData *meshdata, AngularData *angdata, REAL *f, REAL *r) {
	dim3 grid(meshdata->nPlow, meshdata->nPhigh);
	dim3 block(align_power(angdata->slm, COALESCED_NUM(REAL)));
	volumePart<<<grid, block>>>(meshdata->nPlow * meshdata->nPhigh, meshdata->tetstart, meshdata->tetidx, meshdata->tetpos, meshdata->mesh, 
		angdata->slm, angdata->omega, angdata->omega_pos, f, r);
	surfacePart<<<grid, block>>>(meshdata->nPlow * meshdata->nPhigh, meshdata->facestart, meshdata->faceidx, meshdata->facepos, meshdata->bnd, 
		angdata->slm, angdata->Ox, angdata->Oy, angdata->Oz, f, r);
}

#ifdef __cplusplus
}
#endif
