#ifndef __KERNELS_H__

#include "MeshData.h"
#include "AngularData.h"

struct CudaContext {
	REAL *red;
	static void setDevice(int dev);
	const DeviceMeshData *meshdata;
	const DeviceAngularData *angdata;
	CudaContext (int dev, const DeviceMeshData *dmd, const DeviceAngularData *dad) : meshdata(dmd), angdata(dad) {
		red = (REAL *)deviceAlloc(sizeof(REAL));
	}
	void *deviceAlloc(size_t size);
	void deviceFree(void *mem);
	void copyToDev(void *dst, void *src, size_t sz);
	void copyToHost(void *dst, void *src, size_t sz);
	void computeLhs(REAL *f, REAL *Af);
	void computeRhs(REAL *b);
	void addProd(REAL *x, const REAL *y, const REAL wy);
	void mulAdd(REAL *x, const REAL wx, const REAL *y);
	void mulAddProd(REAL *x, const REAL wx, const REAL *y, const REAL wy);
	REAL norm(const REAL *x);
	REAL dot(const REAL *x, const REAL *y);
};

#endif
