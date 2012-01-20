#ifndef __KERNELS_H__

#include <string.h>

#include "MeshData.h"
#include "AngularData.h"

void *deviceAlloc(size_t size);
void deviceFree(void *mem);
void computeRhs(MeshData *meshdata, AngularData *angdata, REAL *f, REAL *r);

#endif
