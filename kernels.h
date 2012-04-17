#ifndef __KERNELS_H__

#include "MeshData.h"
#include "AngularData.h"

void *deviceAlloc(size_t size);
void deviceFree(void *mem);
void copyToDev(void *dst, void *src, size_t sz);
void copyToHost(void *dst, void *src, size_t sz);
void computeRhs(const DeviceMeshData *meshdata, const DeviceAngularData *angdata, REAL *f, REAL *Af, REAL *b);

#endif
