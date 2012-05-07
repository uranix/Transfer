#ifndef __EIGENVALUE_H__
#define __EIGENVALUE_H__

#include "CudaContext.h"

double eigest(CudaContext *ctx, REAL *v, REAL *v0, REAL *w, double *extra);

#endif
