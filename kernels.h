#ifndef __KERNELS_H__

#include <string.h>

#include "MeshData.h"
#include "AngularData.h"

#ifdef __cplusplus
extern "C" {
#endif

void *deviceAlloc(size_t size);
void deviceFree(void *mem);
void computeRhs(struct MeshData *meshdata, struct AngularData *angdata, REAL *f, REAL *r);

#ifdef __cplusplus
}
#endif

#endif
