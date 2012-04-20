#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "AngularData.h"
#include "MeshData.h"
#include "kernels.h"

int main() {
	AngularData ad(1); /* maxk = 1, maxl = 2*/
	MeshData md("mesh.vol");

	DeviceAngularData dad(ad);
	DeviceMeshData dmd(md);

	CudaContext *ctx = new CudaContext (7, &dmd, &dad);

	REAL *f = (REAL *)ctx->deviceAlloc(dad.aslm * dmd.nP * sizeof(REAL));
	REAL *Af = (REAL *)ctx->deviceAlloc(dad.aslm * dmd.nP * sizeof(REAL));
	REAL *b = (REAL *)ctx->deviceAlloc(dad.aslm * dmd.nP * sizeof(REAL));

	REAL *_f = new REAL[dad.aslm * dmd.nP];
	REAL *_Af = new REAL[dad.aslm * dmd.nP];
	REAL *_b = new REAL[dad.aslm * dmd.nP];

	for (int i = 0; i < dad.aslm * dmd.nP; i++)
		_f[i] = 1;

	REAL tau = 0.01;

	ctx->copyToDev(f, _f, dad.aslm * dmd.nP * sizeof(REAL));
	for (int it = 0; it < 1; it++) {
		printf("%d iteration \n", it);
		ctx->computeRhs(f, Af, b);
		ctx->addProd(b, Af, -1);
		ctx->addProd(f, b, tau);
	}
	ctx->copyToHost(_f, f, dad.aslm * dmd.nP * sizeof(REAL));
	ctx->copyToHost(_b, b, dad.aslm * dmd.nP * sizeof(REAL));
	ctx->copyToHost(_Af, Af, dad.aslm * dmd.nP * sizeof(REAL));

	return 0;
}
