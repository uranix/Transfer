#include "aft.h"
REAL orient2d (REAL *pa, REAL *pb, REAL *pc);
REAL orient3d (REAL *pa, REAL *pb, REAL *pc, REAL *pd);
void detinit();

REAL det2i3(REAL *vertex, int v1, int v2, int v3);
int idet2i3(REAL *vertex, int v1, int v2, int v3);
REAL det3i4(REAL *vertex, int v1, int v2, int v3, int v4);
int idet3i4(REAL *vertex, int v1, int v2, int v3, int v4);
