#ifndef __CUDACONTEXT_H__
#define __CUDACONTEXT_H__

#include <cuda.h>

#include "MeshData.h"
#include "AngularData.h"

struct CudaContext {
private:
	CUcontext _ctx;
	CUmodule _mod;
	CUfunction addProdKern;
	CUfunction mulAddKern;
	CUfunction mulAddProdKern;
	CUfunction normKern;
	CUfunction dotKern;
	CUfunction rightHandSide;
	CUfunction volumePart;
	CUfunction surfacePart;
public:	
	REAL *red;
	const DeviceMeshData *meshdata;
	const DeviceAngularData *angdata;
	
	CudaContext(const int dev, const MeshData &dmd, const AngularData &dad);
	~CudaContext();

	void *deviceAlloc(size_t size) const;
	REAL *allocVector() const;
	idx N() const;
	void deviceFree(void *mem) const;
	void copyToDev(void *dst, void *src, size_t sz) const;
	void copyToHost(void *dst, void *src, size_t sz) const;
	REAL *getHostSmall(void *src) const;
	void computeLhs(REAL *f, REAL *Af) const;
	void computeRhs(REAL *b) const;
	void addProd(REAL *x, const REAL *y, const REAL wy) const;
	void mulAdd(REAL *x, const REAL wx, const REAL *y) const;
	void mulAddProd(REAL *x, const REAL wx, const REAL *y, const REAL wy) const;
	REAL norm(const REAL *x) const;
	REAL dot(const REAL *x, const REAL *y) const;
private:
	CudaContext(const CudaContext &);
};

#endif
