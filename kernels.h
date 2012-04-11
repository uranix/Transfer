#ifndef __KERNELS_H__

#include "MeshData.h"
#include "AngularData.h"

void *deviceAlloc(size_t size);
void deviceFree(void *mem);
void computeRhs(const DeviceMeshData *meshdata, const DeviceAngularData *angdata, REAL *f, REAL *r);

#endif
